/**
 * @file blackbody.c
 *
 * This mapping from XYZ tristimulus values to the RGB colorspace is heavily
 * inspired by and based on https://astronomy.stackexchange.com/a/40014/48526,
 * but using the analytic approximations to the color matching functions found
 * at https://en.wikipedia.org/wiki/CIE_1931_color_space#Analytical_approximation
 * and the exact conversion matrix from XYZ to RGB found at
 * https://en.wikipedia.org/wiki/CIE_1931_color_space#Construction_of_the_CIE_XYZ_color_space_from_the_Wright%E2%80%93Guild_data
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <emscripten.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_integration.h"

/// Planck constant @f$ h @f$ in erg s
#define PLANCK      (6.626e-27)

/// Speed of light in vacuum @f$ c @f$ in cm/s
#define LIGHT       (2.998e10)

/// Boltzmann constant @f$ k @f$ in erg/K
#define BOLTZ       (1.381e-16)

/// Centimeters per nanometer conversion factor
#define CENTI       (1e-7)

/// Absolute tolerance for integrals over wavelength
#define ABSTOL      (1e-4)

/// Relative tolerance for integrals over wavelength
#define RELTOL      (1e-3)

/// Maximum number of integration samples
#define LIMIT       (1000)

/// Lower wavelength integration limit in nm
#define XMIN        (300)

/// Upper wavelength integration limit in nm
#define XMAX        (800)

/// Exact conversion matrix data for mapping RGB to XYZ, given by CIE standards
static const double Adata[] = {
    0.49000, 0.31000, 0.20000,
    0.17697, 0.81240, 0.01063,
    0.00000, 0.01000, 0.99000
};

/**
 * @brief Computes the intensity of the blackbody spectrum @f$ B_\lambda @f$
 *
 * @param[in] lam_cm Wavelength in cm
 * @param[in] Tbb Temperature in K
 *
 * @return Intensity in cgs units
 */
EMSCRIPTEN_KEEPALIVE
static double blackbody(double lam_cm, double Tbb)
{
    return ((2 * PLANCK * LIGHT * LIGHT) / gsl_pow_5(lam_cm)) /
       (exp((PLANCK * LIGHT) / (lam_cm * BOLTZ * Tbb)) - 1);
}

/**
 * @brief Helper to evaluate a skew Gaussian function
 *
 * @param[in] x Point where the density is computed
 * @param[in] mu Center of the distribution
 * @param[in] s1 Standard deviation for the left side
 * @param[in] s2 Standard deviation for the right side
 *
 * @return Density
 */
EMSCRIPTEN_KEEPALIVE
static double gaussian(double x, double mu, double s1, double s2)
{
    double dx = x - mu;
    return (x < mu) ? exp(-0.5 * gsl_pow_2(dx / s1)) :
                      exp(-0.5 * gsl_pow_2(dx / s2));
}

/**
 * @brief Computes the value of @f$ B_\lambda \bar{x} @f$
 *
 * @param[in] lam_nm Wavelength in nm
 * @param[in] params Fixed parameters
 *
 * @return Integrand
 */
EMSCRIPTEN_KEEPALIVE
static double integrandX(double lam_nm, void *params)
{
    double *Tbb = (double *) params;
    double lam_cm = lam_nm * CENTI;

    double xbar = 1.056 * gaussian(lam_nm, 599.8, 37.9, 31.0) +
                  0.362 * gaussian(lam_nm, 442.0, 16.0, 26.7) +
                 -0.065 * gaussian(lam_nm, 501.1, 20.4, 26.2);
    return xbar * blackbody(lam_cm, *Tbb);
}

/**
 * @brief Computes the value of @f$ B_\lambda \bar{y} @f$
 *
 * @param[in] lam_nm Wavelength in nm
 * @param[in] params Fixed parameters
 *
 * @return Integrand
 */
EMSCRIPTEN_KEEPALIVE
static double integrandY(double lam_nm, void *params)
{
    double *Tbb = (double *) params;
    double lam_cm = lam_nm * CENTI;

    double ybar = 0.821 * gaussian(lam_nm, 568.8, 46.9, 40.5) +
                  0.286 * gaussian(lam_nm, 530.9, 16.3, 31.1);
    return ybar * blackbody(lam_cm, *Tbb);
}

/**
 * @brief Computes the value of @f$ B_\lambda \bar{z} @f$
 *
 * @param[in] lam_nm Wavelength in nm
 * @param[in] params Fixed parameters
 *
 * @return Integrand
 */
EMSCRIPTEN_KEEPALIVE
static double integrandZ(double lam_nm, void *params)
{
    double *Tbb = (double *) params;
    double lam_cm = lam_nm * CENTI;

    double zbar = 1.217 * gaussian(lam_nm, 437.0, 11.8, 36.0) +
                  0.681 * gaussian(lam_nm, 459.0, 26.0, 13.8);
    return zbar * blackbody(lam_cm, *Tbb);
}

/**
 * @brief Computes an approximate RGB triple for the given blackbody temperature
 *
 * @param[in] Tbb Blackbody temperature in K
 * @param[out] rgb Array of (red, green, blue) values
 */
EMSCRIPTEN_KEEPALIVE
void blackbodyToRgb(double Tbb, double *rgb)
{
    double xyz[3], LUdata[9], err;

    gsl_function function;
    function.params = &Tbb;

    gsl_integration_workspace *wk = gsl_integration_workspace_alloc(LIMIT);

    // calculate X = integral of B_lambda * x_bar
    function.function = &integrandX;
    gsl_integration_qag(&function, XMIN, XMAX, ABSTOL, RELTOL,
        LIMIT, GSL_INTEG_GAUSS51, wk, xyz + 0, &err);

    // calculate Y = integral of B_lambda * y_bar
    function.function = &integrandY;
    gsl_integration_qag(&function, XMIN, XMAX, ABSTOL, RELTOL,
        LIMIT, GSL_INTEG_GAUSS51, wk, xyz + 1, &err);

    // calculate Z = integral of B_lambda * z_bar
    function.function = &integrandZ;
    gsl_integration_qag(&function, XMIN, XMAX, ABSTOL, RELTOL,
        LIMIT, GSL_INTEG_GAUSS51, wk, xyz + 2, &err);

    gsl_integration_workspace_free(wk);

    // normalization fix
    xyz[0] /= xyz[1];
    xyz[2] /= xyz[1];
    xyz[1]  = 1;

    gsl_vector_view xview = gsl_vector_view_array(rgb, 3);
    gsl_vector_const_view bview = gsl_vector_const_view_array(xyz, 3);

    memcpy(LUdata, Adata, 9 * sizeof(double));
    gsl_matrix_view Aview = gsl_matrix_view_array(LUdata, 3, 3);

    // solve [X, Y, Z] = A [R, G, B] and save in output array
    int signum;
    gsl_permutation *p = gsl_permutation_alloc(3);
    gsl_linalg_LU_decomp(&Aview.matrix, p, &signum);
    gsl_linalg_LU_solve(&Aview.matrix, p, &bview.vector, &xview.vector);
    gsl_permutation_free(p);
}

/* vim: set ft=c.doxygen: */
