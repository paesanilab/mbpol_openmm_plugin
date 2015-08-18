#ifndef POLY_3B_H2O_CL_V2X_H
#define POLY_3B_H2O_CL_V2X_H

namespace h2o_cl {

//
// this is the polynomial used by x3b_h2o_cl_v2 (including gradients)
//

struct poly_3b_h2o_cl_v2x {
    static const unsigned n_vars = 21;
    static const unsigned size = 924;

    static double eval(const double a[924],
                       const double x[21]);

    static double eval(const double a[924],
                       const double x[21],
                             double g[21]);
};

} // namespace h2o_cl

#endif // POLY_3B_H2O_CL_V2X_H
