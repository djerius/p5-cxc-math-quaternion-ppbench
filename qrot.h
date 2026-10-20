/* qrot.h: rotation code shared by the PDL::PP kernels and the plain XS
 * functions in qbench.pd.  Plain C, no Perl or PDL API, so the same file
 * can back a PDL-free XS module.
 *
 * q is a scalar-last unit quaternion (x, y, z, w), Hamilton convention;
 * v' = q v q*  computed as  v' = v + 2w(u x v) + 2 u x (u x v),  u = (x, y, z).
 */
#ifndef QROT_H
#define QROT_H

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
