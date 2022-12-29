/**
 * @file orbit.c
 *
 * Tools for computing the orbital displacement of a secondary relative to
 * a primary body
 */

#define _GNU_SOURCE

#include <math.h>
#include <stdlib.h>
#include <emscripten.h>
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv2.h"
#include "gsl/gsl_integration.h"

/// Initial step size
#define H0          (1e-4)

/// Absolute solution tolerance
#define ABSTOL      (1e-4)

/// Relative solution tolerance
#define RELTOL      (1e-3)

/// Maximum number of rootfinding iterations
#define MAXIT       (1000)

/// Maximum number of integration samples
#define LIMIT       (1000)

/// Orbital parameters used in the ODE function
typedef struct
{
    double a;       ///< Semimajor axis
    double b;       ///< Semiminor axis
    double e;       ///< Eccentricity of orbit
    double omega;   ///< Argument of periastron
    double Omega;   ///< Orbital frequency
    double slr;     ///< Semi-latus rectum of ellipse
    double yval;    ///< Target projection value
}
orbit_t;

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
 * @brief Computes the orbital radius at a given angle
 *
 * @param[in] theta Angle around the ellipse (theta = 0 points to observer)
 * @param[in] p Orbit parameters
 */
inline double radius(double theta, const orbit_t *p)
{
    return p->slr / (1 + p->e * cos(theta - p->omega));
}

/**
 * @brief Computes the projection of the planet's position on the
 * plane of the "sky"
 *
 * @param[in] theta Angle around the ellipse (theta = 0 points to observer)
 * @param[in] params Orbit parameters
 */
EMSCRIPTEN_KEEPALIVE
static double projection_value(double theta, void *params)
{
    orbit_t *p = (orbit_t *) params;

    double costh = cos(theta);
    double costh1 = cos(theta - p->omega);
    return p->slr * costh / (1 + p->e * costh1) - p->yval;
}

/**
 * @brief Computes the differential area of the orbital ellipse
 *
 * @param[in] theta Angle around the ellipse (theta = 0 points to observer)
 * @param[in] params Orbit parameters
 *
 * @return @f$ \delta A = r\!\left(\theta\right)^2 @f$
 */
EMSCRIPTEN_KEEPALIVE
static double area_integrand(double theta, void *params)
{
    orbit_t *p = (orbit_t *) params;
    double r = radius(theta, p);
    return sqr(r);
}

/**
 * @brief Evaluates the time derivative @f$ f\!\left(t, y\right) @f$
 *
 * @param[in] t Current time
 * @param[in] y Current function value @f$ \theta\!\left(t\right) @f$
 * @param[out] ydot Current time derivative @f$ \dot{\theta}\!\left(t\right) @f$
 * @param[in] params Input parameters
 *
 * @return Success code
 */
EMSCRIPTEN_KEEPALIVE
static int kepler_law_function(double t, const double y[], double ydot[], void *params)
{
    orbit_t *p = (orbit_t *) params;
    double theta = y[0];
    double r = p->slr / (1 + p->e * cos(theta - p->omega));

    ydot[0] = p->Omega * (p->a * p->b / sqr(r));
    return 0;
}

/**
 * @brief Evaluates the partial derivatives of @f$ f\!\left(t, y\right) @f$
 *
 * @param[in] t Current time
 * @param[in] y Current function value @f$ \theta\!\left(t\right) @f$
 * @param[out] dfdy Jacobian matrix @f$ \partial f / \partial \theta @f$
 * @param[out] dfdt Derivative vector @f$ \partial f / \partial t @f$
 * @param[in] params Input parameters
 *
 * @return Success code
 */
EMSCRIPTEN_KEEPALIVE
static int kepler_law_jacobian(double t, const double y[],
    double dfdy[], double dfdt[], void *params)
{
    orbit_t *p = (orbit_t *) params;
    double theta = y[0];

    double sinth, costh;
    sincos(theta - p->omega, &sinth, &costh);

    dfdt[0] = 0;
    dfdy[0] = -2 * p->e * p->Omega * (p->a * p->b / sqr(p->slr)) *
        (1 + p->e * costh) * sinth;
    return 0;
}

/**
 * @brief Solves for the angle @f$ \theta @f$ for which the planet is
 * visible with projected distance @f$ y @f$
 *
 * @param[in] solver Root solver objet
 * @param[in] y Projected distance on the plane of the "sky"
 * @param[in] p Orbit parameters
 *
 * @return Angle @f$ \theta @f$
 */
EMSCRIPTEN_KEEPALIVE
static double solve_angle_from_projection(gsl_root_fsolver *solver,
    double y, orbit_t *p)
{
    double root;
    size_t niter = 0;

    gsl_function fn;
    fn.function = &projection_value;
    fn.params = p;

    // find the angle between the limiting values; using an algorithm
    // with derivatives (like newton's method) can't guarantee that the
    // root we find will be in the right range
    p->yval = y;
    double thmin = (y > 0) ? 0 : M_PI_2;
    double thmax = (y > 0) ? M_PI_2 : M_PI;
    gsl_root_fsolver_set(solver, &fn, thmin, thmax);

    do
    {
        // iterate on finding the root
        gsl_root_fsolver_iterate(solver);
        root = gsl_root_fsolver_root(solver);
        double resid = projection_value(root, p);

        // check for convergence
        if (fabs(resid) < ABSTOL) break;
    }
    while (niter++ < MAXIT);

    return root;
}

/**
 * @brief Computes the displacement vector (radius and angle) of an
 * elliptic orbit as a function of time
 *
 * @param[in] M Mass of the primary body
 * @param[in] m Mass of the secondary body
 * @param[in] a Semimajor axis
 * @param[in] e Eccentricity
 * @param[in] omega Argument of periastron
 * @param[in] ymin Leftmost requested projected distance
 * @param[in] ymax Rightmost requested projected distance
 * @param[out] t Times at which the orbit is sampled
 * @param[out] r Radius from the center of mass
 * @param[out] theta Angle from the center of mass
 * @param[out] x Cartesian @f$ x @f$ coordinate along orbit
 * @param[out] y Cartesian @f$ y @f$ coordinate along orbit
 * @param[in] n Number of time, radius, and angle samples
 */
EMSCRIPTEN_KEEPALIVE
void orbitArray(double M, double m, double a, double e, double omega,
    double ymin, double ymax, double *t, double *r, double *theta,
    double *x, double *y, size_t n)
{
    orbit_t p;

    // save the inputs in a blob to send to the solver
    p = (orbit_t) { .e = e, .omega = omega, .a = a,
        .b = a * sqrt(1 - e * e), .slr = a * (1 - e * e),
        .Omega = sqrt((M + m) / gsl_pow_3(a)) };

    // solve for limiting angle values
    gsl_root_fsolver *solver =
        gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    double thmax = solve_angle_from_projection(solver, ymin, &p);
    double thmin = solve_angle_from_projection(solver, ymax, &p);
    gsl_root_fsolver_free(solver);

    gsl_integration_workspace *wk = gsl_integration_workspace_alloc(LIMIT);

    gsl_function integrand;
    integrand.function = area_integrand;
    integrand.params = &p;

    // find the amount of time elapsed between thmin and thmax
    double ttot, err;
    gsl_integration_qag(&integrand, thmin, thmax, ABSTOL,
        RELTOL, LIMIT, GSL_INTEG_GAUSS31, wk, &ttot, &err);
    ttot /= p.a * p.b * p.Omega;

    gsl_integration_workspace_free(wk);

    // initial conditions
    r[0] = radius(thmin, &p);
    theta[0] = thmin;

    for (size_t i = 0; i < n; ++i)
    {
        // calculate linear timesteps
        t[i] = i * ttot / (n - 1);
    }

    // definition of the ODE system
    gsl_odeiv2_system system;
    system.function = &kepler_law_function;
    system.jacobian = &kepler_law_jacobian;
    system.dimension = 1;
    system.params = &p;

    // allocate the driver for an implicit integration
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&system,
        gsl_odeiv2_step_msbdf, H0, ABSTOL, RELTOL);

    // set current values from the initial conditions
    double tcur = t[0], th = theta[0];

    for (size_t i = 1; i < n; ++i)
    {
        // step to the next time
        gsl_odeiv2_driver_apply(driver, &tcur, t[i], &th);

        // save the results to the output arrays
        r[i] = radius(th, &p);
        theta[i] = th;
    }

    // free allocated resources
    gsl_odeiv2_driver_free(driver);

    for (size_t i = 0; i < n; ++i)
    {
        double sinth, costh;
        sincos(theta[i], &sinth, &costh);

        // calculate cartesian coords
        x[i] = r[i] * costh;
        y[i] = r[i] * sinth;
    }
}

/* vim: set ft=c.doxygen: */
