#ifndef POLY_2B_H2O_CL_V2X_H
#define POLY_2B_H2O_CL_V2X_H

namespace h2o_cl {

//
// this is the polynomial used by x2b_h2o_cl_v2<4> (including gradients)
//

struct poly_2b_h2o_cl_v2x {
    static const unsigned degree = 4;
    static const unsigned n_vars = 8;
    static const unsigned size = 105;

    static double eval(const double a[105],
                       const double x[8]);

    static double eval(const double a[105],
                       const double x[8],
                             double g[8]);
};

} // namespace h2o_cl

#endif // POLY_2B_H2O_CL_V2X_H
