#include <stdio.h>
#include <math.h>
#include "lamb.h"

/* Cross product: c = a x b */
void cross(double a[3], double b[3], double c[3]) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

double dot(double a[3], double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double mag(double a[3]) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

void print_vec(const char *label, double v[3]) {
    printf("  %s = (%+.15e, %+.15e, %+.15e)\n", label, v[0], v[1], v[2]);
}

int main(void) {
    double mu = 1.0;

    /* 90-degree transfer: r1 along +x, r2 along +y, both at unit distance */
    double r1[3] = {1.0, 0.0, 0.0};
    double r2[3] = {0.0, 1.0, 0.0};

    /* Time of flight: quarter period of circular orbit with r=1, mu=1 */
    /* Period = 2*pi*sqrt(r^3/mu) = 2*pi, quarter = pi/2 */
    double tof = M_PI / 2.0;
    int nrev = 0;

    double v1_pro[3], v2_pro[3];
    double v1_ret[3], v2_ret[3];

    printf("=== Lambert Test: 90-degree transfer, r=1, mu=1 ===\n\n");
    printf("Input:\n");
    print_vec("r1", r1);
    print_vec("r2", r2);
    printf("  mu   = %.15e\n", mu);
    printf("  |r1| = %.15e\n", mag(r1));
    printf("  |r2| = %.15e\n", mag(r2));
    printf("  tof  = %.15e  (pi/2)\n", tof);
    printf("  nrev = %d\n\n", nrev);

    /* --- Prograde: +tof --- */
    printf("=== PROGRADE (tdelt = +%.15e) ===\n", tof);
    int code_pro = lambert(mu, r1, r2, nrev, +tof, v1_pro, v2_pro);
    printf("  return code = %d\n", code_pro);
    print_vec("v1_pro", v1_pro);
    print_vec("v2_pro", v2_pro);

    double h_pro[3];
    cross(r1, v1_pro, h_pro);
    print_vec("h_pro (r1 x v1)", h_pro);
    printf("  |h_pro| = %.15e\n", mag(h_pro));

    /* Energy at r1 */
    double speed_pro_1 = mag(v1_pro);
    double energy_pro = 0.5 * speed_pro_1 * speed_pro_1 - mu / mag(r1);
    printf("  |v1_pro| = %.15e\n", speed_pro_1);
    printf("  energy_pro (at r1) = %.15e\n", energy_pro);

    /* Energy at r2 */
    double speed_pro_2 = mag(v2_pro);
    double energy_pro_r2 = 0.5 * speed_pro_2 * speed_pro_2 - mu / mag(r2);
    printf("  |v2_pro| = %.15e\n", speed_pro_2);
    printf("  energy_pro (at r2) = %.15e\n", energy_pro_r2);

    /* h from r2 x v2 as consistency check */
    double h_pro2[3];
    cross(r2, v2_pro, h_pro2);
    print_vec("h_pro (r2 x v2)", h_pro2);
    printf("  |h_pro (r2xv2)| = %.15e\n", mag(h_pro2));

    printf("\n");

    /* --- Retrograde: -tof --- */
    printf("=== RETROGRADE (tdelt = -%.15e) ===\n", tof);
    int code_ret = lambert(mu, r1, r2, nrev, -tof, v1_ret, v2_ret);
    printf("  return code = %d\n", code_ret);
    print_vec("v1_ret", v1_ret);
    print_vec("v2_ret", v2_ret);

    double h_ret[3];
    cross(r1, v1_ret, h_ret);
    print_vec("h_ret (r1 x v1)", h_ret);
    printf("  |h_ret| = %.15e\n", mag(h_ret));

    /* Energy at r1 */
    double speed_ret_1 = mag(v1_ret);
    double energy_ret = 0.5 * speed_ret_1 * speed_ret_1 - mu / mag(r1);
    printf("  |v1_ret| = %.15e\n", speed_ret_1);
    printf("  energy_ret (at r1) = %.15e\n", energy_ret);

    /* Energy at r2 */
    double speed_ret_2 = mag(v2_ret);
    double energy_ret_r2 = 0.5 * speed_ret_2 * speed_ret_2 - mu / mag(r2);
    printf("  |v2_ret| = %.15e\n", speed_ret_2);
    printf("  energy_ret (at r2) = %.15e\n", energy_ret_r2);

    /* h from r2 x v2 */
    double h_ret2[3];
    cross(r2, v2_ret, h_ret2);
    print_vec("h_ret (r2 x v2)", h_ret2);
    printf("  |h_ret (r2xv2)| = %.15e\n", mag(h_ret2));

    printf("\n");

    /* --- Comparison --- */
    printf("=== COMPARISON ===\n");
    printf("  h_pro z-component = %+.15e\n", h_pro[2]);
    printf("  h_ret z-component = %+.15e\n", h_ret[2]);
    printf("  h_pro direction   = %s\n", h_pro[2] > 0 ? "+z (prograde in xy-plane)" : "-z (retrograde in xy-plane)");
    printf("  h_ret direction   = %s\n", h_ret[2] > 0 ? "+z (prograde in xy-plane)" : "-z (retrograde in xy-plane)");

    return 0;
}
