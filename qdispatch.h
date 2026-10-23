/* qdispatch.h: generic argument dispatch for element-wise operations on
 * Perl data.  Uses only the Perl API (no PDL), so it can live in a
 * PDL-free XS module.
 *
 * Every operation of the API maps fixed-shape inputs to a fixed-shape
 * output, element by element.  An operation is described by a qd_op:
 * the shapes of its one or two inputs and of its output, which input the
 * result takes its form from, which inputs may also be the output (in
 * place), and a kernel that does the arithmetic for one element.
 * qd_dispatch() does the rest for every operation:
 *
 *   1. classify each input: a single value, an array of values, columns,
 *      or an object (ndarray).  The input's rank says how deeply a single
 *      value is nested, so the depth of the first-element chain decides;
 *   2. hand ndarrays and other objects to the PDL side, through a function
 *      pointer it publishes in $QBench::QD_PDL;
 *   3. check lengths, choose the result's form, prepare the output (new,
 *      into a caller's container, or in place) and check in place;
 *   4. loop: read each multi-valued input, run the kernel, write the
 *      result.  Single-valued inputs are read once, before the loop;
 *   5. truncate the output to the result's length.
 */
#ifndef QDISPATCH_H
#define QDISPATCH_H

#include <string.h>

/* forced inlining for the per-element helpers: left to itself, the compiler
   made them separate functions, costing about 120 instructions per element */
#if defined(__GNUC__)
#  define QD_ALWAYS_INLINE static inline __attribute__((always_inline))
#else
#  define QD_ALWAYS_INLINE static inline
#endif

struct qd_op;
struct qd_operand;

/* the element loop of one operation; see QD_DEFINE_RUN */
typedef void (*qd_run_fn)(pTHX_ const struct qd_op *op, const struct qd_operand *a,
                          const struct qd_operand *b, const struct qd_operand *o,
                          SSize_t n, double *ia, double *ib);

typedef struct {
    int         rank;       /* 0 scalar, 1 vector, 2 3x3 matrix */
    int         ncomp;      /* number of doubles in one value: 1, 2, 3, 4 or 9 */
    const char *keys[4];    /* column keys (rank 1); keys[0] NULL: no columns form */
} qd_shape;

typedef struct qd_op {
    const char *name;
    const char *argname[2];
    int         n_in;       /* 1 or 2 inputs */
    qd_shape    in[2];
    qd_shape    out;
    int         primary[2]; /* the result takes the form of in[primary[0]] if it is
                               multi-valued, otherwise of in[primary[1]]; -1: none */
    unsigned    inplace;    /* bit i set: in[i] may also be the output */
    qd_run_fn   run;        /* its element loop, made by QD_DEFINE_RUN */
    int         pdl_op;     /* index into the PDL side's table */
} qd_op;

enum { QD_SINGLE = 1, QD_ARRAY, QD_COLUMNS, QD_OBJECT };

/* number of columns of a shape: only vectors (rank 1, at most 4 components)
   have a columns form; 0 otherwise */
#define QD_NCOLS(s) ((s)->rank == 1 && (s)->ncomp <= 4 ? (s)->ncomp : 0)

typedef struct qd_operand {
    const qd_shape *shape;
    int     form;
    SSize_t n;              /* number of values; 1 for a single value */
    SV     *sv;             /* the argument */
    AV     *av;             /* SINGLE (rank > 0): the value; ARRAY: the outer array */
    AV     *col[4];         /* COLUMNS */
    int     fresh;          /* output only: a container just created here, so values
                               are stored as new scalars, with no lookups */
} qd_operand;

/* the PDL side's entry point, published in $QBench::QD_PDL */
typedef SV *(*qd_pdl_fn)(pTHX_ const qd_op *op, SV **in, SV *out);

/* --- classifying an input ------------------------------------------------ */

static void
qd_classify(pTHX_ const qd_op *op, const char *what, const qd_shape *s, SV *sv, qd_operand *o)
{
    SV *rv;
    o->shape = s; o->sv = sv; o->n = 1; o->av = NULL;
    if (!SvROK(sv)) {
        if (s->rank == 0 && SvOK(sv)) { o->form = QD_SINGLE; return; }
        croak("%s: %s is not a recognised data form", op->name, what);
    }
    rv = SvRV(sv);
    if (SvOBJECT(rv)) { o->form = QD_OBJECT; return; }

    if (SvTYPE(rv) == SVt_PVHV) {
        int k;
        if (s->rank != 1 || !s->keys[0]) croak("%s: %s has no columns form", op->name, what);
        for (k = 0; k < QD_NCOLS(s); k++) {
            SV **e = hv_fetch((HV *)rv, s->keys[k], strlen(s->keys[k]), 0);
            if (!e || !SvROK(*e) || SvTYPE(SvRV(*e)) != SVt_PVAV)
                croak("%s: %s columns need a '%s' array reference", op->name, what, s->keys[k]);
            o->col[k] = (AV *)SvRV(*e);
        }
        o->n = av_count(o->col[0]);
        for (k = 1; k < QD_NCOLS(s); k++)
            if (av_count(o->col[k]) != (Size_t)o->n)
                croak("%s: %s columns have different lengths", op->name, what);
        o->form = QD_COLUMNS;
        return;
    }

    if (SvTYPE(rv) != SVt_PVAV) croak("%s: %s is not a recognised data form", op->name, what);
    o->av = (AV *)rv;
    if (av_count(o->av) == 0) { o->form = QD_ARRAY; o->n = 0; return; }   /* empty array of values */
    {
        /* depth of the first-element chain, up to rank + 1 */
        AV *a = o->av;
        int depth = 1;
        while (depth <= s->rank) {
            SV **e;
            if (av_count(a) == 0) break;
            e = av_fetch(a, 0, 0);
            if (!e || !SvROK(*e) || SvTYPE(SvRV(*e)) != SVt_PVAV) break;
            a = (AV *)SvRV(*e);
            depth++;
        }
        if (depth == s->rank)     { o->form = QD_SINGLE; return; }
        if (depth == s->rank + 1) { o->form = QD_ARRAY; o->n = av_count(o->av); return; }
        croak("%s: %s is not nested as a single value or an array of values", op->name, what);
    }
}

/* --- reading -------------------------------------------------------------- */

QD_ALWAYS_INLINE void
qd_get(pTHX_ const qd_op *op, AV *av, double *d, int n, const char *what)
{
    int i;
    if (av_count(av) < (Size_t)n) croak("%s: %s must have %d elements", op->name, what, n);
    for (i = 0; i < n; i++) {
        SV **e = av_fetch(av, i, 0);
        d[i] = e ? SvNV(*e) : 0;
    }
}

QD_ALWAYS_INLINE AV *
qd_elem_av(pTHX_ const qd_op *op, AV *av, SSize_t i, const char *what)
{
    SV **e = av_fetch(av, i, 0);
    if (!e || !SvROK(*e) || SvTYPE(SvRV(*e)) != SVt_PVAV)
        croak("%s: %s element %ld is not an array reference", op->name, what, (long)i);
    return (AV *)SvRV(*e);
}

/* one value of rank >= 1 held in array v */
QD_ALWAYS_INLINE void
qd_get_value(pTHX_ const qd_op *op, const qd_shape *s, AV *v, double *d, const char *what)
{
    int r;
    if (s->rank == 1) { qd_get(aTHX_ op, v, d, s->ncomp, what); return; }
    if (av_count(v) < 3) croak("%s: %s must have 3 rows", op->name, what);
    for (r = 0; r < 3; r++) qd_get(aTHX_ op, qd_elem_av(aTHX_ op, v, r, what), d + 3 * r, 3, what);
}

/* value i of input o */
QD_ALWAYS_INLINE void
qd_read(pTHX_ const qd_op *op, const qd_operand *o, SSize_t i, double *d, const char *what)
{
    switch (o->form) {
    case QD_SINGLE:
        if (o->shape->rank == 0) d[0] = SvNV(o->sv);
        else qd_get_value(aTHX_ op, o->shape, o->av, d, what);
        break;
    case QD_ARRAY:
        if (o->shape->rank == 0) { SV **e = av_fetch(o->av, i, 0); d[0] = e ? SvNV(*e) : 0; }
        else qd_get_value(aTHX_ op, o->shape, qd_elem_av(aTHX_ op, o->av, i, what), d, what);
        break;
    case QD_COLUMNS: {
        int k;
        for (k = 0; k < QD_NCOLS(o->shape); k++) {
            SV **e = av_fetch(o->col[k], i, 0);
            d[k] = e ? SvNV(*e) : 0;
        }
        break;
    }
    }
}


/* --- writing -------------------------------------------------------------- */

QD_ALWAYS_INLINE void
qd_set(pTHX_ AV *av, SSize_t i, double val)
{
    SV **e = av_fetch(av, i, 1);
    if (!e) croak("cannot store element %ld", (long)i);
    sv_setnv(*e, val);
    SvSETMAGIC(*e);
}

/* the array in element i of av, created if it is not an array ref */
QD_ALWAYS_INLINE AV *
qd_slot(pTHX_ AV *av, SSize_t i)
{
    AV *v;
    SV **e = av_fetch(av, i, 1);
    if (!e) croak("cannot store element %ld", (long)i);
    if (SvROK(*e) && SvTYPE(SvRV(*e)) == SVt_PVAV) return (AV *)SvRV(*e);
    v = newAV();
    sv_setsv(*e, sv_2mortal(newRV_noinc((SV *)v)));
    return v;
}

/* one value of rank >= 1 into array v, reusing its scalars; resized exactly */
QD_ALWAYS_INLINE void
qd_put_value(pTHX_ const qd_shape *s, AV *v, const double *d)
{
    int k, r;
    if (s->rank == 1) {
        av_extend(v, s->ncomp - 1);
        for (k = 0; k < s->ncomp; k++) qd_set(aTHX_ v, k, d[k]);
        av_fill(v, s->ncomp - 1);
        return;
    }
    av_extend(v, 2);
    for (r = 0; r < 3; r++) {
        AV *row = qd_slot(aTHX_ v, r);
        av_extend(row, 2);
        for (k = 0; k < 3; k++) qd_set(aTHX_ row, k, d[3 * r + k]);
        av_fill(row, 2);
    }
    av_fill(v, 2);
}

/* fill empty array v with one value of rank >= 1, as new scalars */
QD_ALWAYS_INLINE void
qd_fill_new(pTHX_ const qd_shape *s, AV *v, const double *d)
{
    int k, r;
    if (s->rank == 1) {
        av_extend(v, s->ncomp - 1);
        for (k = 0; k < s->ncomp; k++) av_store(v, k, newSVnv(d[k]));
        return;
    }
    av_extend(v, 2);
    for (r = 0; r < 3; r++) {
        AV *row = newAV();
        av_extend(row, 2);
        for (k = 0; k < 3; k++) av_store(row, k, newSVnv(d[3 * r + k]));
        av_store(v, r, newRV_noinc((SV *)row));
    }
}

/* value i of output o */
QD_ALWAYS_INLINE void
qd_write(pTHX_ const qd_operand *o, SSize_t i, const double *d)
{
    if (o->fresh) {                                    /* a container created by qd_dispatch */
        switch (o->form) {
        case QD_SINGLE: qd_fill_new(aTHX_ o->shape, o->av, d); return;
        case QD_ARRAY: {
            AV *v = newAV();
            qd_fill_new(aTHX_ o->shape, v, d);
            av_store(o->av, i, newRV_noinc((SV *)v));
            return;
        }
        case QD_COLUMNS: {
            int k;
            for (k = 0; k < QD_NCOLS(o->shape); k++) av_store(o->col[k], i, newSVnv(d[k]));
            return;
        }
        }
    }
    switch (o->form) {
    case QD_SINGLE: qd_put_value(aTHX_ o->shape, o->av, d); break;
    case QD_ARRAY:  qd_put_value(aTHX_ o->shape, qd_slot(aTHX_ o->av, i), d); break;
    case QD_COLUMNS: {
        int k;
        for (k = 0; k < QD_NCOLS(o->shape); k++) qd_set(aTHX_ o->col[k], i, d[k]);
        break;
    }
    }
}

/* --- the element loop ------------------------------------------------------ */

/* QD_DEFINE_RUN(fn, kernel, SA, SB, SO) defines fn, the element loop for one
   operation: its kernel, kernel(a, b, r), and the shapes of its inputs and
   output (SA, SB, SO: qd_shape initialisers, as in its descriptor), are
   compile-time constants here, so the compiler inlines the kernel, unrolls
   the loops over components and drops the tests on rank.  The descriptor
   must use the same shapes.  (One shared loop, calling the kernel through
   a pointer and reading the shapes from the descriptor, cost about 60-80
   instructions more per element than a hand-written loop.)

   Single-valued inputs are read before the loop (into ia, ib).  With the
   helpers above forced inline, the loop makes no calls except into the Perl
   API.  Compiling a further copy for each combination of forms was tried:
   it saved about 10 instructions per element (about 2%), cost 37-41 more
   per call for single values, and added 66 KB of code, so it was dropped. */
#define QD_DEFINE_RUN(fn, kernel, SA, SB, SO)                                        \
static void                                                                          \
fn(pTHX_ const qd_op *op, const qd_operand *a, const qd_operand *b,                  \
   const qd_operand *o, SSize_t n, double *ia, double *ib)                           \
{                                                                                    \
    /* the shapes as constants, and local copies of the operands: nothing  */       \
    /* outside this function can change them, so the compiler keeps them  */       \
    /* in registers across the Perl API calls */                                    \
    static const qd_shape sa = SA, sb = SB, so = SO;                                 \
    qd_operand la = *a, lb, lo = *o;                                                 \
    double r[9];                                                                     \
    SSize_t j;                                                                       \
    la.shape = &sa; lo.shape = &so;                                                  \
    if (b) { lb = *b; lb.shape = &sb; }                                              \
    for (j = 0; j < n; j++) {                                                        \
        if (la.form != QD_SINGLE) qd_read(aTHX_ op, &la, j, ia, op->argname[0]);     \
        if (b && lb.form != QD_SINGLE) qd_read(aTHX_ op, &lb, j, ib, op->argname[1]); \
        kernel(ia, ib, r);                                                           \
        qd_write(aTHX_ &lo, j, r);                                                   \
    }                                                                                \
}

/* --- the dispatcher ------------------------------------------------------- */

/* args: the inputs, then the optional output; returns a new reference to the result */
static SV *
qd_dispatch(pTHX_ const qd_op *op, SV **args, I32 items)
{
    qd_operand in[2], out;
    double ia[9], ib[9];
    SV *outsv;
    int i, k, have_out, rform = QD_SINGLE, any_multi = 0;
    SSize_t n = 1;

    if (items < op->n_in || items > op->n_in + 1)
        croak("%s: expected %d or %d arguments, got %d", op->name, op->n_in, op->n_in + 1, (int)items);
    outsv = items > op->n_in ? args[op->n_in] : &PL_sv_undef;
    have_out = SvOK(outsv);

    for (i = 0; i < op->n_in; i++) qd_classify(aTHX_ op, op->argname[i], &op->in[i], args[i], &in[i]);

    /* ndarrays and other objects: the PDL side */
    if (in[0].form == QD_OBJECT || (op->n_in > 1 && in[1].form == QD_OBJECT)
        || (have_out && SvROK(outsv) && SvOBJECT(SvRV(outsv)))) {
        static qd_pdl_fn fn = NULL;
        if (!fn) {
            SV *hook = get_sv("QBench::QD_PDL", 0);
            if (!hook || !SvIOK(hook)) croak("%s: ndarrays need the PDL part of the module", op->name);
            fn = INT2PTR(qd_pdl_fn, SvIV(hook));
        }
        return fn(aTHX_ op, args, outsv);
    }

    /* the length of the result, and its form */
    for (i = 0; i < op->n_in; i++) {
        if (in[i].form == QD_SINGLE) continue;
        if (any_multi && in[i].n != n)
            croak("%s: multi-valued arguments have different lengths (%ld and %ld)", op->name, (long)n, (long)in[i].n);
        n = in[i].n;
        any_multi = 1;
    }
    for (i = 0; i < 2; i++) {
        int p = op->primary[i];
        if (p >= 0 && p < op->n_in && in[p].form != QD_SINGLE) { rform = in[p].form; break; }
    }
    if (rform == QD_COLUMNS && !op->out.keys[0]) croak("%s: the result has no columns form", op->name);

    /* the output: new, or a caller's container (possibly an input: in place) */
    out.shape = &op->out; out.form = rform; out.n = n; out.sv = outsv; out.fresh = !have_out;
    if (!have_out) {
        if (rform == QD_COLUMNS) {
            HV *hv = newHV();
            outsv = sv_2mortal(newRV_noinc((SV *)hv));
            for (k = 0; k < QD_NCOLS(&op->out); k++) {
                out.col[k] = newAV();
                (void)hv_store(hv, op->out.keys[k], strlen(op->out.keys[k]), newRV_noinc((SV *)out.col[k]), 0);
            }
        }
        else {
            out.av = newAV();
            outsv = sv_2mortal(newRV_noinc((SV *)out.av));
        }
    }
    else {
        SV *rv;
        if (!SvROK(outsv)) croak("%s: the output must be a reference", op->name);
        rv = SvRV(outsv);
        if (rform == QD_COLUMNS) {
            if (SvTYPE(rv) != SVt_PVHV) croak("%s: the output must be a hash reference, for columns", op->name);
            for (k = 0; k < QD_NCOLS(&op->out); k++) {
                SV **e = hv_fetch((HV *)rv, op->out.keys[k], strlen(op->out.keys[k]), 1);
                if (!e) croak("%s: cannot store column '%s'", op->name, op->out.keys[k]);
                if (!(SvROK(*e) && SvTYPE(SvRV(*e)) == SVt_PVAV)) {
                    if (SvOK(*e)) croak("%s: output column '%s' must be an array reference", op->name, op->out.keys[k]);
                    sv_setsv(*e, sv_2mortal(newRV_noinc((SV *)newAV())));
                }
                out.col[k] = (AV *)SvRV(*e);
            }
        }
        else {
            if (SvTYPE(rv) != SVt_PVAV) croak("%s: the output must be an array reference", op->name);
            out.av = (AV *)rv;
        }
        /* in place: the output is an input's container, or has its column arrays */
        for (i = 0; i < op->n_in; i++) {
            int same = SvROK(args[i]) && SvRV(args[i]) == rv;
            if (!same && in[i].form == QD_COLUMNS && rform == QD_COLUMNS && in[i].shape->ncomp == op->out.ncomp) {
                same = 1;
                for (k = 0; k < QD_NCOLS(&op->out); k++) if (in[i].col[k] != out.col[k]) { same = 0; break; }
            }
            if (same && !(op->inplace & (1u << i)))
                croak("%s: %s cannot also be the output", op->name, op->argname[i]);
        }
    }

    if (n && rform == QD_ARRAY) av_extend(out.av, n - 1);
    if (n && rform == QD_COLUMNS) for (k = 0; k < QD_NCOLS(&op->out); k++) av_extend(out.col[k], n - 1);

    /* single-valued inputs are read once */
    if (in[0].form == QD_SINGLE) qd_read(aTHX_ op, &in[0], 0, ia, op->argname[0]);
    if (op->n_in > 1 && in[1].form == QD_SINGLE) qd_read(aTHX_ op, &in[1], 0, ib, op->argname[1]);

    op->run(aTHX_ op, &in[0], op->n_in > 1 ? &in[1] : NULL, &out, n, ia, ib);

    if (rform == QD_ARRAY) av_fill(out.av, n - 1);
    if (rform == QD_COLUMNS) for (k = 0; k < QD_NCOLS(&op->out); k++) av_fill(out.col[k], n - 1);
    return SvREFCNT_inc(outsv);
}

#endif
