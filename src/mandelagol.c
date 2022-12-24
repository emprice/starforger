#include <math.h>
#include <stdlib.h>
#include <emscripten.h>
#include "gsl/gsl_sf_ellint.h"

#define PRECISION   (GSL_PREC_SINGLE)

/**
 * @brief Convenience function for squaring the input value
 * @param[in] x Value
 * @return Square of value, @f$ x^2 @f$
 */
inline double sqr(double x)
{
    return x * x;
}

/**
 * @brief Evaluates the Mandel & Agol (2008) quadratic limb darkened
 * transit model at a single value of the input parameters
 *
 * @param[in] p Relative size of the planet (planet radius divided by
 * stellar radius)
 * @param[in] zin Relative distance from the center of the planet to
 * the center of the star (semimajor axis divided by stellar radius)
 * @param[in] gamma1 Linear limb darkening coefficient
 * @param[in] gamma2 Quadratic limb darkening coefficient
 *
 * @return Relative flux from the star
 *
 * @attention Mandel & Agol (2002) uses the sign convention of
 * Gradshteyn & Ryzhik (1994) for the elliptic integrals, but GSL
 * uses a different convention for the parameter @f$ n @f$.
 * When evaluating @f$ \Pi @f$, make the substitution @f$ n \rightarrow -n @f$
 * from what is given in Mandel & Agol (2002). See section 8.111 in
 * http://fisica.ciens.ucv.ve/~svincenz/TISPISGIMR.pdf, for example.
 */
EMSCRIPTEN_KEEPALIVE
double mandelagol(double p, double zin, double gamma1, double gamma2)
{
    double c2 = gamma1 + 2 * gamma2, c4 = -gamma2;
    double c0 = 1 - c2 - c4;
    double Omega = c0 / 4 + c2 / 6 + c4 / 8;
    double lambdad = 0, lambdae = 0, etad = 0;

    double z = fabs(zin);
    double p2 = sqr(p), z2 = sqr(z);
    double a = sqr(z - p), b = sqr(z + p);
    double k = sqrt((1 - a) / (4 * p * z)), q = p2 - z2;
    double q2 = sqr(q);

    double kappa1 = acos((1 - p2 + z2) / (2 * z));
    double kappa0 = acos((p2 + z2 - 1) / (2 * p * z));
    double eta2 = 0.5 * p2 * (p2 + 2 * z2);

    if ((fabs(1 - p) < z) && (z <= 1 + p))
    {
        lambdae = M_1_PI * (p2 * kappa0 + kappa1 -
            sqrt(0.25 * (4 * z2 - sqr(1 + z2 - p2))));
    }
    else if (z <= 1 - p)
    {
        lambdae = p2;
    }
    else if ((p > 1) && (z <= p - 1))
    {
        lambdae = 1;
    }

    if (((0 < p) && (1 + p <= z)) || (p == 0))
    {
        // case I

        lambdad = etad = 0;
    }
    else if (((0.5 + fabs(p - 0.5) < z) && (z < 1 + p)) ||
             ((0.5 < p) && (fabs(1 - p) <= z) && (z < p)))
    {
        // case II and VIII

        double Kk = gsl_sf_ellint_Kcomp(k, PRECISION);
        double Ek = gsl_sf_ellint_Ecomp(k, PRECISION);
        double Pk = gsl_sf_ellint_Pcomp(k, (1 - a) / a, PRECISION);

        lambdad = (((1 - b) * (2 * b + a - 3) - 3 * q * (b - 2)) * Kk +
            4 * p * z * (z2 + 7 * p2 - 4) * Ek - 3 * (q / a) * Pk) /
            (9 * M_PI * sqrt(p * z));
        etad = 0.5 * M_1_PI * (kappa1 + 2 * eta2 * kappa0 +
            -0.25 * (1 + 5 * p2 + z2) * sqrt((1 - a) * (b - 1)));
    }
    else if (((p < 0.5) && (p < z) && (z < 1 - p)) ||
             ((p < 1) && (0 < z) && (z < 0.5 - fabs(p - 0.5))))
    {
        // case III and IX

        double Kk = gsl_sf_ellint_Kcomp(1. / k, PRECISION);
        double Ek = gsl_sf_ellint_Ecomp(1. / k, PRECISION);
        double Pk = gsl_sf_ellint_Pcomp(1. / k, (b - a) / a, PRECISION);

        lambdad = 2 * ((1 - 5 * z2 + p2 + q2) * Kk +
            (1 - a) * (z2 + 7 * p2 - 4) * Ek +
            -3 * (q / a) * Pk) / (9 * M_PI * sqrt(1 - a));
        etad = eta2;
    }
    else if ((p < 0.5) && (z == 1 - p))
    {
        // case IV

        lambdad = (2 * M_1_PI / 3) * acos(1 - 2 * p) -
            (4 * M_1_PI / 9) * (3 + 2 * p - 8 * p2);
        etad = eta2;
    }
    else if ((p < 0.5) && (z == p))
    {
        // case V

        double Kk = gsl_sf_ellint_Kcomp(2 * p, PRECISION);
        double Ek = gsl_sf_ellint_Ecomp(2 * p, PRECISION);

        lambdad = (1. / 3) + (2 * M_1_PI / 9) *
            (4 * (2 * p2 - 1) * Ek + (1 - 4 * p2) * Kk);
        etad = eta2;
    }
    else if ((p == 0.5) && (z == 0.5))
    {
        // case VI

        lambdad = (1. / 3) - (4 * M_1_PI / 9);
        etad = 3. / 32;
    }
    else if ((0.5 < p) && (z == p))
    {
        // case VII

        double Kk = gsl_sf_ellint_Kcomp(0.5 / p, PRECISION);
        double Ek = gsl_sf_ellint_Ecomp(0.5 / p, PRECISION);

        lambdad = (1. / 3) + (16 * p * M_1_PI / 9) * (2 * p2 - 1) * Ek -
            ((1 - 4 * p2) * (3 - 8 * p2) / (9 * M_PI * p)) * Kk;
        etad = 0.5 * M_1_PI * (kappa1 + 2 * eta2 * kappa0 +
            -0.25 * (1 + 5 * p2 + z2) * sqrt((1 - a) * (b - 1)));
    }
    else if ((p < 1) && (z == 0))
    {
        // case X

        lambdad = -1 * (2. / 3) * pow(1 - p2, 1.5);
        etad = eta2;
    }
    else if ((1 < p) && (0 <= z) && (z < p - 1))
    {
        // case XI

        lambdae = 1;
        etad = 0.5;
    }

    return 1. - (0.25 / Omega) * ((1 - c2) * lambdae +
        c2 * (lambdad + (2. / 3.) * (p > z)) - c4 * etad);
}

EMSCRIPTEN_KEEPALIVE
void mandelagolArray(double p, const double *z,
    double gamma1, double gamma2, double *flux, size_t n)
{
    for (size_t i = 0; i < n; ++i)
    {
        flux[i] = mandelagol(p, z[i], gamma1, gamma2);
    }
}

/* vim: set ft=c.doxygen: */
