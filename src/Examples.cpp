#include "Utils.h"
#include <iomanip>
  
//================================================
// Example for how to create and use a 1D spline
//================================================

void example_make_spline_basic(){

  const double xmin = 0.0;
  const double xmax = 1.0;
  const int    npts = 10;
  
  Vector x_array = Utils::linspace(xmin, xmax, npts);
  Vector y_array = exp(x_array);

  Spline f_spline;
  f_spline.create(x_array, y_array, "Function y = exp(x)");

  std::cout << "e^log(2) = " << f_spline( log(2) ) << "\n";
}


void example_make_spline(){
  
  //=================================================================
  //=================================================================

  // The range and number of points
  const double xmin = 0.0;
  const double xmax = 2.0*M_PI;
  const int    npts = 100;

  // A test function to generate some data with
  // and its derivatives (to test the spline derivative routines)
  auto func      = [&](double x){ return  sin(x); };
  auto dfuncdx   = [&](double x){ return  cos(x); };
  auto ddfuncddx = [&](double x){ return -sin(x); };

  // Make an x-array with npts points between xmin and xmax and a y=y(x) array
  Vector x_array = Utils::linspace(xmin, xmax, npts);
  Vector y_array(npts);

  // Fill the y-array = f(x_array)
  // Easier way: std::transform(x_array.begin(), x_array.end(), y_array.begin(), func);
  for(size_t i = 0; i < x_array.size(); i++) 
    y_array[i] = func(x_array[i]);

  // Make the spline
  Spline f_spline;
  f_spline.create(x_array, y_array, "Test Spline");

  // Check that it gives OK results
  std::cout << "Example getting function values from a spline\n";
  std::cout << "# x   f_spline(x)   f_exact(x)\n";
  for(int i = 0; i < 4; i++){
    double x = (i+1.0)/4.0;
    std::cout << std::setw(8) << x           << " ";
    std::cout << std::setw(8) << f_spline(x) << " ";
    std::cout << std::setw(8) << func(x)     << "\n";
  }
  
  //=================================================================
  //=================================================================

  // If you evaluate it out of bounds it will return the closest value
  std::cout << "Spline gives that f(-1.0) = " << f_spline(-1.0) << " which is f(0) = " << f_spline(0.0) << "\n"; 

  // To get a warning when this happens either set _FIDUCIAL_SPLINE_WARNING = true in Spline.h
  // or turn it on for any given spline by running (very useful for debugging)
  // If you don't want warnings by default then set _FIDUCIAL_SPLINE_WARNING = false in Spline.h
  f_spline.set_out_of_bounds_warning(true);

  // Test
  std::cout << "Example getting a warning if its out of bounds\n";
  std::cout << "This next line should give a warning: \n";
  f_spline(xmin-1.0);
  
  //=================================================================
  //=================================================================

  // We can also evaluate derivatives from the spline
  std::cout << "Example getting derivative from 1D spline\n";
  std::cout << "# x   dfdx_spline(x)   dfdx_exact(x)\n";
  for(int i = 0; i < 4; i++){
    double x = (i+1.0)/4.0;
    std::cout << std::setw(8) << x                   << " ";
    std::cout << std::setw(8) << f_spline.deriv_x(x) << " ";
    std::cout << std::setw(8) << dfuncdx(x)          << "\n";
  }
  
  // We can also evaluate second derivatives from the spline
  std::cout << "Example getting second derivative from 1D spline\n";
  std::cout << "# x   ddfddx_spline(x)   ddfddx_exact(x)\n";
  for(int i = 0; i < 4; i++){
    double x = (i+1.0)/4.0;
    std::cout << std::setw(8) << x                   << " ";
    std::cout << std::setw(8) << f_spline.deriv_xx(x) << " ";
    std::cout << std::setw(8) << ddfuncddx(x)   << "\n";
  }
  
  //=================================================================
  //=================================================================
}

//================================================
// Example for how to create and use a 2D spline
//================================================

void example_make_2D_spline(){
 
  // The range and number of points
  const double xmin = 0.0;
  const double xmax = 2.0*M_PI;
  const int    nx   = 100;
  const double ymin = 0.0;
  const double ymax = 2.0*M_PI;
  const int    ny   = 101;

  // A test function to generate some data with
  // and its derivatives (to test the spline derivative routines)
  auto function = [&](double x, double y){
    return x + y + x*y;
  };

  // Make an x-array, an y-array and a z = z(x,y) array
  Vector x_array = Utils::linspace(xmin, xmax, nx);
  Vector y_array = Utils::linspace(ymin, ymax, ny);
  Vector z_array = Vector(nx*ny);
  
  // Fill the z array
  for(int ix = 0; ix < nx; ix++)
    for(int iy = 0; iy < ny; iy++)
      z_array[ix + nx * iy] = function(x_array[ix], y_array[iy]);

  // Make a spline
  Spline2D z_spline(x_array, y_array, z_array, "Test 2D spline");

  // Test that it gives good results
  std::cout << "Example getting function values from 2D spline\n";
  std::cout << "# x   y  z_spline(x,y)   z_exact(x,y)\n";
  for(int ix = 0; ix < 4; ix++){
    for(int iy = 0; iy < 4; iy++){
      double x = (ix+1.0)/4.0;
      double y = (iy+1.0)/4.0;
      std::cout << std::setw(8) << x             << " ";
      std::cout << std::setw(8) << y             << " ";
      std::cout << std::setw(8) << z_spline(x,y) << " ";
      std::cout << std::setw(8) << function(x,y) << "\n";
    }
  }
  
  //=================================================================
  //=================================================================

  // We have the same functionality as for the 1D spline
  // and we can compute up to second order derivatives using
  // z_spline.deriv_x(x,y)
  // z_spline.deriv_y(x,y)
  // z_spline.deriv_xx(x,y)
  // z_spline.deriv_yy(x,y)
  // z_spline.deriv_xy(x,y)
  
  // Test that it gives good results
  std::cout << "Example getting derivatives from 2D spline\n";
  std::cout << "# x   y  ddzdxdy_spline(x,y)   ddzdxdy_exact(x,y)\n";
  for(int ix = 0; ix < 4; ix++){
    for(int iy = 0; iy < 4; iy++){
      double x = (ix+1.0)/4.0;
      double y = (iy+1.0)/4.0;
      std::cout << std::setw(8) << x             << " ";
      std::cout << std::setw(8) << y             << " ";
      std::cout << std::setw(8) << z_spline.deriv_xy(x,y) << " ";
      std::cout << std::setw(8) << 1.0           << "\n";
    }
  }
  
  //=================================================================
  //=================================================================
  
  // Alternative way of making a spline: from a 2D vector z[ix][iy]
  Vector2D zz_array(nx, Vector(ny));
  for(int ix = 0; ix < nx; ix++)
    for(int iy = 0; iy < ny; iy++){
      zz_array[ix][iy] = function(x_array[ix], y_array[iy]);
    }
  Spline2D zz_spline(x_array, y_array, zz_array, "Test 2D spline");

  //=================================================================
  //=================================================================

}

//================================================
// Example for how to solve an ODE
//================================================

void example_single_ODE(){

  //=================================================================
  //=================================================================

  // The range [xmin, xmax] to solve the ODE over and an array
  // of x-values for which we store the solution on
  const double xmin = 0.0;
  const double xmax = 2.0*M_PI;
  const int    npts = 10;
  Vector x_array = Utils::linspace(xmin, xmax, npts);

  // The ODE y' = 2*x (solution is x^2)
  ODEFunction dydx = [&](double x, const double *y, double *dydx){
    dydx[0] = 2.0*x;
    return GSL_SUCCESS;
  };

  // Set the IC (with the given IC the solution is y(x) = x^2
  double yini = 0.0;
  Vector y_ic{yini};

  // Solve the ODE
  ODESolver ode;
  ode.solve(dydx, x_array, y_ic);

  // Get the solution (we only have one component so index 0 holds the solution)
  auto y_array = ode.get_data_by_component(0);

  // Output the result and compare to analytical result
  std::cout << "Example solving an ODE and comparing to analytical result\n";
  for(int i = 0; i < x_array.size(); i++){
    std::cout << " xi          = " << std::setw(10) << x_array[i];
    std::cout << " y(xi)       = " << std::setw(10) << y_array[i];
    std::cout << " y_exact(xi) = " << std::setw(10) << x_array[i]*x_array[i] << "\n";
  }
}


void example_solve_coupled_ODE(){

  //=================================================================
  //=================================================================

  // The range [xmin, xmax] to solve the ODE over and an array
  // of x-values for which we store the solution on
  const double xmin = 0.0;
  const double xmax = 2.0*M_PI;
  const int    npts = 10;
  Vector x_array = Utils::linspace(xmin, xmax, npts);

  // The ODE y0' = y1 ;  y1' = -y0 
  ODEFunction dydx = [&](double x, const double *y, double *dydx){
    dydx[0] =  y[1];
    dydx[1] = -y[0];
    return GSL_SUCCESS;
  };

  // The analytical solution to the ODE above 
  // y0(x) = y(0)sin(x) + y'(0)cos(x)
  // and y1(x) = y0'(x)
  auto analytical_solution = [&](Vector & yic, double x){
    Vector y(2);
    y[0] =  yic[0] * cos(x) + yic[1] * sin(x);
    y[1] = -yic[0] * sin(x) + yic[1] * cos(x);
    return y;
  };

  // Set the IC (with the given IC the solution is y0(x) = sin(x)
  // and y1(x) = cos(x) )
  double y0_ini = 0.0;
  double y1_ini = 1.0;
  Vector y_ic{y0_ini, y1_ini};

  // Solve the ODE
  ODESolver ode;
  ode.solve(dydx, x_array, y_ic);

  // Get the data: this is a Vector2D with data[i] = {y0(xi), y1(xi)} 
  auto all_data = ode.get_data();

  // Output the result and compare to analytical result
  std::cout << "Example solving an ODE and comparing to analytical result\n";
  for(int i = 0; i < all_data.size(); i++){
    auto y_exact = analytical_solution(y_ic, x_array[i]);
    std::cout << " xi     = " << std::setw(10) << x_array[i];
    std::cout << " y0(xi) = " << std::setw(10) << all_data[i][0]  << " y0_exact(xi) = " << std::setw(10) << y_exact[0];
    std::cout << " y1(xi) = " << std::setw(10) << all_data[i][1]  << " y1_exact(xi) = " << std::setw(10) << y_exact[1] << "\n";
  }

  //=================================================================
  //=================================================================

  // Other ways of getting the data: just the final point - this is a Vector with {y0(xend), y1(xend)} 
  auto final_data = ode.get_final_data();
  std::cout << "Example extracting data from the ODE solver\n";
  std::cout << " y0(xmax) = " << final_data[0]  << " y0_exact(xmax) = " << sin(xmax);
  std::cout << " y1(xmax) = " << final_data[1]  << " y1_exact(xmax) = " << cos(xmax) << "\n";

  // We can also get the data for the derivatives dydx[]
  auto derivative_data = ode.get_derivative_data();

  //...and we can get the data at x_array[ix] by calling
  auto data_x_point_2 = ode.get_data_by_xindex(2);
  std::cout << "y(" << x_array[2] << ") = " << data_x_point_2[0] << " y_exact = " << sin(x_array[2]) << "\n";

  //=================================================================
  //=================================================================

  // If you want higher accuracy then call this function before solving
  // Here hstart is the first step-size to try, abserr is the absolute error limit
  // and relerr is the relative error limit, e.g. (1e-3, 1e-10, 1e-10)
  // Lower values of abserr and relerr will lead to a longer time to solve it
  double hstart = 1e-3, abserr = 1e-10, relerr = 1e-10;
  ode.set_accuracy(hstart, abserr, relerr);

  //=================================================================
  //=================================================================

  // If you want the solver to print the solution as it integrates along 
  // run ode.set_verbose(true); before solving
  std::cout << "Example running the ODE solver with verbose on\n";
  ode.set_verbose(true);
  ode.solve(dydx, x_array, y_ic);

  //=================================================================
  //=================================================================

  // We can also choose the ODE method, e.g. to use Runge-Kutta 2
  // instead of the fiducial RK 4 method. See ODESolver.h
  // for a list of availiable methods, but be aware that some of them
  // might also require the Jacobian of the system to be provided!
  ode.solve(dydx, x_array, y_ic, gsl_odeiv2_step_rk2);
}

void other_stuff(){
  
  //=================================================================
  //=================================================================
  
  // Do not use [new] to allocate memory. Its a bad idea.
  // You will forget to delete it and get memory leaks
  double *arr = new double[100];
  //...
  //...
  delete[] arr;
  
  // Instead when making an array use a standard container
  // like Vector == std::vector<double>
  Vector arr1(100);
  
  //=================================================================
  //=================================================================

  // To make a lineary spaced array with n points between xstart and xend 
  // you can use Utils::linspace
  double x_min = 0.0;
  double x_max = 1.0;
  int n = 10;
  auto arr2 = Utils::linspace(x_min, x_max, n);
  for(auto &x: arr2)
    std::cout << "x = " << x << "\n";
  
  //=================================================================
  //=================================================================

  // If you want to evaluate a standard mathematical function on
  // all elements in a Vector you can just do this (macro in Utils)
  auto arr3 = sin(arr2);
  for(auto &x: arr3)
    std::cout << "sin(x) = " << x << "\n";
 
  //=================================================================
  //=================================================================

  // You can make local functions (lambdas) just as you would make a double
  // or similar (remember ; at the end)
  auto func = [&](double x){
    return x*x;
  };
  std::cout << "5^2 = " << func(5) << "\n";
 
  //=================================================================
  //=================================================================
  
}

int main(int argc, char **argv){
  std::cout << "\n==============================\n\n";
  example_make_spline_basic();
  
  std::cout << "\n==============================\n\n";
  example_make_spline();
  
  std::cout << "\n==============================\n\n";
  example_single_ODE();
  
  std::cout << "\n==============================\n\n";
  example_solve_coupled_ODE();
  
  std::cout << "\n==============================\n\n";
  example_make_2D_spline();

  std::cout << "\n==============================\n\n";
  other_stuff();
}
