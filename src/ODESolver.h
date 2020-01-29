#ifndef _ODESOLVER_HEADER
#define _ODESOLVER_HEADER
#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

//===================================================
// 
// This is a wrapper around the GSL library to easily
// solve ODEs and return the data in whatever format
// you want.
//
// Supplying an x-array the solution will be stored at
// each of the points in the array (must be monotonic)
//
// Example use: 
// Solve the system dy0/dx=y1, dy1/dx=-y0 with
// y0 = 2.0 and y1 = -2.0 on the interval [0,1]
// and store the solution at x = 0.0, 0.25, 0.5 and 1.0
//---------------------------------------------------
// ODEFunction dydx = [&](double x, const double *y, double *dydx){
//   dydx[0] =  y[1];
//   dydx[1] = -y[0];
//   return GSL_SUCCESS;
// };
// Vector yini{2.0, -2.0}
// Vector x_array{0.0, 0.25, 0.5, 1.0};
// ODESolver ode;
// ode.solve(dydx, x_array, yini);
// auto solution = ode.get_data();
//---------------------------------------------------
//
// Choices of steppers (fiducial one set below):
// gsl_odeiv2_step_rk2;
// gsl_odeiv2_step_rk4;
// gsl_odeiv2_step_rkf45;
// gsl_odeiv2_step_rkck;
// gsl_odeiv2_step_rk8pd;
// gsl_odeiv2_step_rk2imp;
// gsl_odeiv2_step_rk4imp;
// gsl_odeiv2_step_bsimp;
// gsl_odeiv2_step_rk1imp;
// gsl_odeiv2_step_msadams;
// gsl_odeiv2_step_msbdf;
//
//===================================================

// The fiducial stepper if not provided by the user
#define _FIDUCIAL_STEPPER gsl_odeiv2_step_rk4

using ODEFunctionPointer         = int (*)(double, const double[], double [], void*);
using ODEFunctionPointerJacobian = int (*)(double, const double[], double [], double[], void*);
using ODEFunction                = std::function<int(double, const double *, double *)>;
using ODEFunctionJacobian        = std::function<int(double, const double *, double *, double *)>;
using Vector                     = std::vector<double>;
using Vector2D                   = std::vector<Vector>;

extern ODEFunctionJacobian *no_jacobian_ptr;

class ODESolver{
  private:

    const gsl_odeiv2_step_type * stepper = _FIDUCIAL_STEPPER;

    // Fiducial accuracy parameters
    double hstart    = 1e-3;
    double abserr    = 1e-7;
    double relerr    = 1e-7;

    // Show information about the solution as we integrate
#ifdef _FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE
    bool verbose     = true;
#else
    bool verbose     = false;
#endif

    int nequations   = 1;
    int num_x_points = 0;

    std::vector< Vector > data;
    std::vector< Vector > derivative_data;

  public:

    ODESolver() = default;
    ODESolver(double hstart, double abserr, double relerr);
    ODESolver(const ODESolver & rhs) = delete;
    ODESolver& operator=(const ODESolver & rhs) = delete;

    void solve(
        ODEFunctionPointer ode_equation, 
        void *parameters, 
        Vector& xarr, 
        Vector& yinitial, 
        const gsl_odeiv2_step_type * stepper = _FIDUCIAL_STEPPER,
        ODEFunctionPointerJacobian jacobian = nullptr);
    void solve(
        ODEFunction & ode_equation, 
        Vector& xarr, 
        Vector& yinitial, 
        const gsl_odeiv2_step_type * stepper = _FIDUCIAL_STEPPER,
        ODEFunctionJacobian & jacobian = *no_jacobian_ptr
        );

    void set_verbose(bool onoff);
    void set_accuracy(const double hstart, const double abserr, const double relerr);

    // Get the solution at the end point
    Vector get_final_data() const;
    double get_final_data_by_component(int icomponent) const;

    // Get all the data z_ij = ( y_i(xarr_j) )
    Vector2D get_data() const;

    // Get all the data transposed z_ji = ( y_i(xarr_j) )
    Vector2D get_data_transpose() const;

    // Get the data for a particular component y_i(xarr)
    Vector get_data_by_component(int icomponent) const;

    // Get the data at a particular x-index y(xarr_i)
    Vector get_data_by_xindex(int ix) const;

    // Get all the data for the derivatives ( dy_i/dx(xarr) )_i=1^nequations
    Vector2D get_derivative_data() const;

    // Get the data dy_i/dx(xarr) for the derivatives for a particular component
    Vector get_derivative_data_by_component(int icomponent) const;
};

#endif
