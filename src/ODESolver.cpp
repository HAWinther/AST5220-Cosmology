#include "ODESolver.h"

const char *ODE_ERROR = "Error in the ODE routines\n";

//==================================================================
// Wrappers for conversions between std::function and the raw 
// function pointer required by GSL
//==================================================================

ODEFunctionJacobian no_jacobian = [](double t, const double *y, double *dfdy, double *dfdt){
  return GSL_SUCCESS;
};
ODEFunctionJacobian *no_jacobian_ptr = &no_jacobian;

struct ParamWrapper{
  ODEFunction & f;
  ODEFunctionJacobian & J;
  ParamWrapper(ODEFunction & f,  ODEFunctionJacobian & J) : f(f), J(J){}
};

int ode_wrapper(double x, const double *y, double *dydx, void *param){
  struct ParamWrapper *p = (struct ParamWrapper *) param;
  ODEFunction & deriv = p->f;
  return deriv(x, y, dydx);
}

int jacobian_wrapper(double t, const double *y, double *dfdy, double *dfdt, void *param){
  struct ParamWrapper *p = (struct ParamWrapper *) param;
  ODEFunctionJacobian & jacobian = p->J;
  return jacobian(t, y, dfdy, dfdt);
}

//==================================================================
// Constructors
//==================================================================

ODESolver::ODESolver(double hstart, double abserr, double relerr) : hstart(hstart), abserr(abserr), relerr(relerr) {}

void ODESolver::solve(
    ODEFunction & ode_equation, 
    Vector& xarr, 
    Vector& yinitial, 
    const gsl_odeiv2_step_type *stepper,
    ODEFunctionJacobian & jacobian){
  struct ParamWrapper equations(ode_equation, jacobian);
  solve(ode_wrapper, &equations, xarr, yinitial, stepper, jacobian_wrapper);
}

//==================================================================
// Class methods
//==================================================================

void ODESolver::solve(
    ODEFunctionPointer ode_equation, 
    void *parameters, 
    Vector& xarr, 
    Vector& yinitial, 
    const gsl_odeiv2_step_type *stepper,
    ODEFunctionPointerJacobian jacobian){

  // Store the number of equations and size of x-array
  nequations = yinitial.size();
  num_x_points = xarr.size();
  if(num_x_points < 2) {
    std::cout <<"Error ODESolver The xarray needs to have atleast size 2\n";
    throw ODE_ERROR;
  }
  if(nequations < 1) {
    std::cout << "Error ODESolver The yarray needs to have atleast size 1\n";
    throw ODE_ERROR;
  }

  // Set up the ODE system
  gsl_odeiv2_system ode_system = {ode_equation, jacobian, size_t(yinitial.size()), parameters};
  gsl_odeiv2_driver* ode_driver = gsl_odeiv2_driver_alloc_y_new (&ode_system, stepper, hstart, abserr, relerr);

  // Initialize with the initial condition
  double x = xarr[0];
  Vector y(yinitial);
  Vector dydx(nequations);
  ode_equation(x, y.data(), dydx.data(), parameters);

  // Allocate memory for the the results: data[i][j] = y_j(x_i)
  data            = std::vector< Vector >(xarr.size(), yinitial);
  derivative_data = std::vector< Vector >(xarr.size(), dydx);
    
  if(verbose){
    std::cout << "ODESolver step " << std::setw(5) << 0 << " / " << num_x_points-1 << " x: [" << std::setw(10) << x << "] ";
    std::cout << "y: [";
    for(auto &ycomp : y) std::cout << " " << std::setw(10) << ycomp << " ";
    std::cout << "]" << std::endl;
  }

  // Solve it step by step...
  for(int i = 1; i < num_x_points; i++){
    const double xnew = xarr[i];

    if(verbose){
      std::cout << "ODESolver step " << std::setw(5) << i << " / " << num_x_points-1 << " x: [" << std::setw(10) << xnew << "] ";
    }

    // Integrate from x -> xnew
    const int status = gsl_odeiv2_driver_apply(ode_driver, &x, xnew, y.data() );
    if(status != GSL_SUCCESS) {
      std::cout << "Error ODESolver GSL returned non-successful\n";
      throw ODE_ERROR;
    }

    if(verbose){
      std::cout << "y: [";
      for(auto &ycomp : y) std::cout << " " << std::setw(10) << ycomp << " ";
      std::cout << "]" << std::endl;
    }

    // Store the derivative
    ode_equation(x, y.data(), dydx.data(), parameters);

    // Copy over result
    data[i] = y;
    derivative_data[i] = dydx;
  }
}

void ODESolver::set_verbose(bool onoff){
  verbose = onoff;
}

Vector2D ODESolver::get_data() const{
  return data;
}

Vector2D ODESolver::get_data_transpose() const{
  auto data_transpose = std::vector< Vector >(nequations, Vector(num_x_points));
  for(int ix = 0; ix < num_x_points; ix++){
    for(int icomponent = 0; icomponent < nequations; icomponent++){
      data_transpose[icomponent][ix] = data[ix][icomponent];
    }
  }
  return data_transpose;
}

Vector ODESolver::get_final_data() const{
  Vector res(nequations);
  for(int icomponent = 0; icomponent < nequations; icomponent++)
    res[icomponent] = data[num_x_points-1][icomponent];
  return res;
}

double ODESolver::get_final_data_by_component(int icomponent) const{
  return data[num_x_points-1][icomponent];
}

Vector ODESolver::get_data_by_component(int icomponent) const{
  Vector res(num_x_points);
  for(int ix = 0; ix < num_x_points; ix++)
    res[ix] = data[ix][icomponent];
  return res;
}

Vector ODESolver::get_data_by_xindex(int ix) const{
  Vector res(nequations);
  for(int icomponent = 0; icomponent < nequations; icomponent++)
    res[icomponent] = data[ix][icomponent];
  return res;
}

Vector2D ODESolver::get_derivative_data() const{
  return derivative_data;
}

Vector ODESolver::get_derivative_data_by_component(int icomponent) const{
  Vector res(num_x_points);
  for(int ix = 0; ix < num_x_points; ix++)
    res[ix] = derivative_data[ix][icomponent];
  return res;
}

void ODESolver::set_accuracy(const double h, const double a, const double r){
  hstart = h;
  abserr = a;
  relerr = r;
}

