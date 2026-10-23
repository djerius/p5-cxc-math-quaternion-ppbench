/* qrot.h: rotation code shared by the PDL::PP kernels and the plain XS
 * functions in qbench.pd.  Plain C, no Perl or PDL API, so the same file
 * can back a PDL-free XS module.
 *
 * q is a scalar-last unit quaternion (x, y, z, w), Hamilton convention;
 * v' = q v q*  computed as  v' = v + 2w(u x v) + 2 u x (u x v),  u = (x, y, z).
 */
#ifndef QROT_H
#define QROT_H

#include <math.h>

/* canonical sign: q and -q are the same rotation; keep w >= 0 */
static inline void
qrot_canon(double *q)
{
    if (q[3] < 0) { q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3]; }
}

/* Hamilton product r = a b: apply b first, then a */
static inline void
qrot_mult(const double *a, const double *b, double *r)
{
    double x =  a[3]*b[0] - a[2]*b[1] + a[1]*b[2] + a[0]*b[3];
    double y =  a[2]*b[0] + a[3]*b[1] - a[0]*b[2] + a[1]*b[3];
    double z = -a[1]*b[0] + a[0]*b[1] + a[3]*b[2] + a[2]*b[3];
    double w = -a[0]*b[0] - a[1]*b[1] - a[2]*b[2] + a[3]*b[3];
    r[0] = x; r[1] = y; r[2] = z; r[3] = w;
    qrot_canon(r);
}

/* fractional rotation: same axis as q, t times its angle (shorter arc) */
static inline void
qrot_pow(const double *q, double t, double *r)
{
    double c[4] = { q[0], q[1], q[2], q[3] }, half, s, k;
    qrot_canon(c);
    half = acos(c[3] > 1 ? 1 : c[3]);
    s = sin(half);
    if (s < 1e-15) { r[0] = r[1] = r[2] = 0; r[3] = 1; return; }
    k = sin(t * half) / s;
    r[0] = c[0] * k; r[1] = c[1] * k; r[2] = c[2] * k; r[3] = cos(t * half);
    qrot_canon(r);
}

/* rotation matrix, row-major: m[3*row + col]; m . v equals qrot_vec(q, v) */
static inline void
qrot_to_matrix(const double *q, double *m)
{
    double x = q[0], y = q[1], z = q[2], w = q[3];
    m[0] = 1 - 2*(y*y + z*z); m[1] = 2*(x*y - z*w);     m[2] = 2*(x*z + y*w);
    m[3] = 2*(x*y + z*w);     m[4] = 1 - 2*(x*x + z*z); m[5] = 2*(y*z - x*w);
    m[6] = 2*(x*z - y*w);     m[7] = 2*(y*z + x*w);     m[8] = 1 - 2*(x*x + y*y);
}

static inline void
qrot_vec(const double *q, const double *v, double *r)
{
    double x = q[0], y = q[1], z = q[2], w = q[3];
    double a = v[0], b = v[1], c = v[2];
    double tx = 2*(y*c - z*b), ty = 2*(z*a - x*c), tz = 2*(x*b - y*a);
    r[0] = a + w*tx + (y*tz - z*ty);
    r[1] = b + w*ty + (z*tx - x*tz);
    r[2] = c + w*tz + (x*ty - y*tx);
}

#endif
