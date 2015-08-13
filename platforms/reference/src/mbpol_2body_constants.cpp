#include "mbpol_2body_constants.h"
#include <algorithm>
#include <cmath>
#include "openmm/reference/RealVec.h"

using OpenMM::RealVec;
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_exp(const double& r0, const double& k,
                       RealVec& O1, RealVec& O2)
{

    g = O1 - O2;

    const double r = std::sqrt(g.dot(g));

    const double exp1 = std::exp(k*(r0 - r));
    const double gg = - k*exp1/r;

    g *= gg;

    return exp1;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_coul(const double& r0, const double& k,
        RealVec& O1, RealVec& O2)
{
    g = O1 - O2;

    const double rsq = g.dot(g);
    const double r = std::sqrt(rsq);

    const double exp1 = std::exp(k*(r0 - r));
    const double rinv = 1.0/r;
    const double val = exp1*rinv;

    const double gg = - (k + rinv)*val*rinv;

    g *= gg;

    return val;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void variable::grads(const double& gg, RealVec& O1, RealVec& O2) const
{
        RealVec d = g*gg;

        O1 += d;
        O2 -= d;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::setup(const RealVec& O, const RealVec& H1, const RealVec& H2,
                    const double& in_plane_g, const double& out_of_plane_g,
                    RealVec& x1, RealVec& x2)
{
    oh1 = H1 - O;
    oh2 = H2 - O;

    const RealVec v = oh1.cross(oh2);
    const RealVec in_plane = O + (oh1 + oh2) * 0.5*in_plane_g;
    const RealVec out_of_plane = v * out_of_plane_g;

    x1 = in_plane + out_of_plane;
    x2 = in_plane - out_of_plane;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::grads(RealVec& g1, RealVec& g2,
                    const double& in_plane_g, const double& out_of_plane_g,
                    RealVec& O, RealVec& H1, RealVec& H2 ) const
{
    const RealVec gm = g1-g2;

    const RealVec t1 = oh2.cross(gm);

    const RealVec t2 = oh1.cross(gm);

    const RealVec gsum = g1 + g2;
    const RealVec in_plane = gsum*0.5*in_plane_g;

    const RealVec gh1 = in_plane + t1*out_of_plane_g;
    const RealVec gh2 = in_plane - t2*out_of_plane_g;

    O += gsum - (gh1 + gh2); // O
    H1 += gh1; // H1
    H2 += gh2; // H2
}

double f_switch(const double& r, double& g)
{
    if (r > r2f) {
        g = 0.0;
        return 0.0;
    } else if (r > r2i) {
        const double t1 = M_PI/(r2f - r2i);
        const double x = (r - r2i)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}
double f_switch_chloride(const double& r, double& g)
{
    if (r > r2f_chloride) {
        g = 0.0;
        return 0.0;
    } else if (r > r2i_chloride) {
        const double t1 = M_PI/(r2f_chloride - r2i_chloride);
        const double x = (r - r2i_chloride)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}
