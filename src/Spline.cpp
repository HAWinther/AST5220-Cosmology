#include "Spline.h"

const char *SPLINE_ERROR = "Error in spline routines\n";

//====================================================
// Constructors
//====================================================

Spline::Spline(
    double *x, 
    double *y, 
    const int nx, 
    std::string splinename, 
    const gsl_interp_type *interpoltype){
  create(x, y, nx, splinename, interpoltype);
}

Spline::Spline(
    Vector &x, 
    Vector &y, 
    std::string splinename, 
    const gsl_interp_type *interpoltype){
  create(x, y, splinename, interpoltype);
}
    
Spline::Spline(std::string name) : name(name) {}

//====================================================
// Assignment constructor
//====================================================
Spline& Spline::operator=(const Spline& rhs){
  // We just create the spline from scratch instead of copying all the data
  if(rhs.spline)
    create(rhs.spline->x, rhs.spline->y, rhs.size_x, rhs.name, rhs.interpoltype_used);
  else
    free();
  return *this;
}

//====================================================
// Copy constructor
//====================================================
Spline::Spline(const Spline& rhs){
  // We just create the spline from scratch instead of copying all the data
  if(rhs.spline)
    create(rhs.spline->x, rhs.spline->y, rhs.size_x, rhs.name, rhs.interpoltype_used);
  else
    free();
}
    
//====================================================
// Create a GSL spline
//====================================================
void Spline::create(
    double *x, 
    double *y, 
    const int nx, 
    std::string splinename, 
    const gsl_interp_type *interpoltype){
 
  // Clean up if we already have a spline allocated
  if(spline) free();

  // Set class variables
  xmin   = x[0];
  xmax   = x[nx-1];
  dx_min = (x[1]-x[0])/2.0;
  dx_max = (x[nx-1]-x[nx-2])/2.0;
  size_x = nx;
  name = splinename;
  interpoltype_used = interpoltype;

  // Make the spline
  spline = gsl_spline_alloc(interpoltype, nx);
  gsl_spline_init(spline, x, y, nx);

  // Make accelerators (one per thread if OpenMP)
#ifdef _USEOPENMP
  int nthreads = 1;
#pragma omp parallel
  {
    int id = omp_get_thread_num();
    if(id == 0) nthreads = omp_get_num_threads();
  }
  xaccs = std::vector<gsl_interp_accel*>(nthreads);
  for(auto & xa: xaccs)
    xa = gsl_interp_accel_alloc();
#else
  xacc = gsl_interp_accel_alloc();
#endif
}

void Spline::create(
    Vector &x, 
    Vector &y, 
    std::string splinename, 
    const gsl_interp_type *interpoltype){
  if(x.size() != y.size()){
    std::cout << "Error Spline::create [" << splinename << "]: x and y array must have the same number of elements\n";
    throw SPLINE_ERROR;
  } 
  create(x.data(), y.data(), int(x.size()), splinename, interpoltype);
}

//====================================================
// Evaluate the function
// NB: using closest points if out-of-bounds!
//====================================================

double Spline::operator()(double x) const{
  return eval(x);
}

double Spline::eval(const double x) const{
  if(!spline) {
    std::cout << "Error Spline::eval [" << name << "] Spline has not been created!\n";
    throw SPLINE_ERROR;
  }

  // If out of bounds show a warning and set x to boundary value
  out_of_bounds_check(x);
  double xx = x;
  if(x < xmin) xx = xmin;
  if(x > xmax) xx = xmax;

  // Return f, f' or f'' depending on value of deriv
#ifdef _USEOPENMP
  gsl_interp_accel *xacc_thread = xaccs[omp_get_thread_num()];
#else
  gsl_interp_accel *xacc_thread = xacc;
#endif
  return gsl_spline_eval(spline, xx, xacc_thread);
}

double Spline::eval_deriv(const double x, const int deriv) const{
  if(!spline) {
    std::cout << "Error Spline::eval_deriv Spline [" << name << "] has not been created!\n";
    throw SPLINE_ERROR;
  }
  if(deriv < 0 || deriv > 2) {
    std::cout << "Error Spline::eval_deriv Spline [" << name << "] Wrong value for deriv ";
    std::cout << "It's only possible to eval 0=f, 1=f' and 2=f''!\n";
    throw SPLINE_ERROR;
  }

  // If out of bounds show a warning and set x to boundary value
  out_of_bounds_check(x);
  double xx = x;
  if(x < xmin) xx = xmin;
  if(x > xmax) xx = xmax;

  // Return f, f' or f'' depending on value of deriv
#ifdef _USEOPENMP
  gsl_interp_accel *xacc_thread = xaccs[omp_get_thread_num()];
#else
  gsl_interp_accel *xacc_thread = xacc;
#endif
  if(deriv == 0)
    return gsl_spline_eval(spline, xx, xacc_thread);
  else if(deriv == 1)
    return gsl_spline_eval_deriv(spline, xx, xacc_thread);
  else if(deriv == 2)
    return gsl_spline_eval_deriv2(spline, xx, xacc_thread);
  return 0.0;
}

//====================================================
// Free up memory
//====================================================
void Spline::free(){
  if(spline){

    // Reset class variables
    xmin = xmax = 0.0;
    size_x = 0;
    dx_min = dx_max = 0.0;
    
    // Free the spline
    gsl_spline_free(spline);
    spline = nullptr;

    // Free accelerators
#ifdef _USEOPENMP
    for(auto & xa: xaccs)
      gsl_interp_accel_free(xa);
    std::vector<gsl_interp_accel *>().swap(xaccs);
#else
    gsl_interp_accel_free(xacc);
    xacc = nullptr;
#endif
  }
}

//====================================================
// Rest of the class methods
//====================================================
void Spline::out_of_bounds_check(const double x) const{
  if(out_of_bounds_warning){
    if(x < xmin - dx_min || x > xmax + dx_max){
      std::cout << "Warning Spline[" << name << "] ";
      std::cout << "x = " << x << " is out of bounds (" << xmin << "," << xmax << ")\n";
    }
  }
}
double Spline::deriv_x(const double x) const{
  return eval_deriv(x, 1);
}
double Spline::deriv_xx(const double x) const{
  return eval_deriv(x, 2);
}
std::pair<double,double> Spline::get_xrange() const{
  return {xmin, xmax};
}
std::string Spline::get_name() const{
  return name;
}
void Spline::set_out_of_bounds_warning(bool v){
  out_of_bounds_warning = v;
}

//====================================================
// Destructor
//====================================================
Spline::~Spline(){
  free();
}

//====================================================
// Constructors
//====================================================

Spline2D::Spline2D(
    double *x,
    double *y, 
    double *z, 
    const int nx, 
    const int ny, 
    std::string splinename, 
    const gsl_interp2d_type *interpoltype) {
  create(x, y, z, nx, ny, splinename, interpoltype);
}

Spline2D::Spline2D(
    Vector &x, 
    Vector &y, 
    Vector &z, 
    std::string splinename, 
    const gsl_interp2d_type *interpoltype){
  create(x, y, z, splinename, interpoltype);
}

Spline2D::Spline2D(
    Vector &x, 
    Vector &y, 
    Vector2D &z, 
    std::string splinename, 
    const gsl_interp2d_type *interpoltype){
  create(x, y, z, splinename, interpoltype);
}

Spline2D::Spline2D(std::string name) : name(name) {}

//====================================================
// Assignment constructor
//====================================================
Spline2D& Spline2D::operator=(const Spline2D& rhs){
  // We just create the spline from scratch instead of copying all the data
  if(rhs.spline)
    create(rhs.spline->xarr, rhs.spline->yarr, rhs.spline->zarr, 
        rhs.size_x, rhs.size_y,rhs.name, rhs.interpoltype_used);
  else
    free();
  return *this;
}

//====================================================
// Copy constructor
//====================================================
Spline2D::Spline2D(const Spline2D& rhs){
  // We just create the spline from scratch instead of copying all the data
  if(rhs.spline)
    create(rhs.spline->xarr, rhs.spline->yarr, rhs.spline->zarr, 
        rhs.size_x, rhs.size_y,rhs.name, rhs.interpoltype_used);
  else
    free();
}

// For more easy evaluation
double Spline2D::operator()(double x, double y) const{
  return eval(x, y);
}

//====================================================
// Create a GSL 2D spline
//====================================================
void Spline2D::create(
    double *x, 
    double *y, 
    double *z, 
    const int nx, 
    const int ny, 
    std::string splinename, 
    const gsl_interp2d_type *interpoltype){
  
  // Clean up if we already have a spline
  if(spline) free();

  // Set class variables
  xmin   = x[0];
  xmax   = x[nx-1];
  ymin   = y[0];
  ymax   = y[ny-1];
  dx_min = (x[1]-x[0])/2.0;
  dx_max = (x[nx-1]-x[nx-2])/2.0;
  dy_min = (y[1]-y[0])/2.0;
  dy_max = (y[ny-1]-y[ny-2])/2.0;
  size_x = nx;
  size_y = ny;
  name = splinename;
  interpoltype_used = interpoltype;

  // Create spline
  spline = gsl_spline2d_alloc(interpoltype, nx, ny);
  gsl_spline2d_init(spline, x, y, z, nx, ny);

  // Create accelerators (one per thread if OpenMP)
#ifdef _USEOPENMP
  int nthreads = 1;
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    if(id == 0) nthreads = omp_get_num_threads();
  }
  xaccs = std::vector<gsl_interp_accel*>(nthreads);
  yaccs = std::vector<gsl_interp_accel*>(nthreads);
  for(auto & xa: xaccs)
    xa = gsl_interp_accel_alloc();
  for(auto & ya: yaccs)
    ya = gsl_interp_accel_alloc();
#else
  xacc   = gsl_interp_accel_alloc();
  yacc   = gsl_interp_accel_alloc();
#endif
}

void Spline2D::create(
    Vector& x, 
    Vector& y, 
    Vector& z, 
    std::string splinename, 
    const gsl_interp2d_type *interpoltype){
  if(x.size()*y.size() != z.size()){ 
    std::cout << "Error Spline2D: z array have wrong number of elements nx*ny != nz\n";
    throw SPLINE_ERROR;
  }
  create(x.data(), y.data(), z.data(), x.size(), y.size(), splinename, interpoltype);
}

void Spline2D::create(
    Vector& x, 
    Vector& y, 
    Vector2D& z, 
    std::string splinename, 
    const gsl_interp2d_type *interpoltype){
  int nz_x = z.size();
  int nz_y = nz_x == 0 ? 0 : z[0].size();
  if(nz_x * nz_y != x.size() * y.size()){
    std::cout << "Error Spline2D: z array have wrong dimensions\n";
    throw SPLINE_ERROR;
  }

  for(int i = 0; i < nz_x; i++){
    if(z[i].size() != nz_y){
      std::cout << "Error Spline2D: z array have wrong dimensions\n";
      throw SPLINE_ERROR;
    }
  }
  Vector f(nz_x*nz_y);
  for(int iy = 0; iy < nz_y; iy++){
    for(int ix = 0; ix < nz_x; ix++){
      f[ix + nz_x * iy] = z[ix][iy];
    }
  }
  create(x.data(), y.data(), f.data(), x.size(), y.size(), splinename, interpoltype);
}

//====================================================
// Lookup a value from a GSL 2D spline
// Use closest points for out-of-bounds
//====================================================
double Spline2D::eval(const double x, const double y) const{
  if(!spline) {
    std::cout << "Error Spline2D::eval Spline has not been created!\n";
    exit(1);
  }

  // If out of bounds show a warning and set x,y to boundary value
  out_of_bounds_check(x,y);
  double xx = x, yy = y;
  if(x < xmin) xx = xmin; 
  if(x > xmax) xx = xmax; 
  if(y < ymin) yy = ymin; 
  if(y > ymax) yy = ymax;

#ifdef _USEOPENMP
  gsl_interp_accel *xacc_thread = xaccs[omp_get_thread_num()];
  gsl_interp_accel *yacc_thread = yaccs[omp_get_thread_num()];
#else
  gsl_interp_accel *xacc_thread = xacc;
  gsl_interp_accel *yacc_thread = yacc;
#endif
  return gsl_spline2d_eval(spline, xx, yy, xacc_thread, yacc_thread);
}

double Spline2D::eval_deriv(
    const double x, 
    const double y, 
    const int derivx, 
    const int derivy) const{
  // Map (dx,dy) => n = derivx + 3*derivy 
  // which gives 0 = f, 1 = f_x, 2 = f_xx, 3 = f_y, 4 = f_xy and (f_xyy), 6 = f_yy
  const int n = derivx + 3*derivy;

  if(!spline) {
    std::cout << "Error Spline2D::eval_deriv Spline [" << name << "] has not been created!\n";
    throw SPLINE_ERROR;
  }
  if(n < 0 || n > 7) {
    std::cout << "Error Spline2D::eval_deriv Spline [" << name << "] n is out of bounds. ";
    std::cout << "It's only possible to eval (0,0)=f, (1,0)=f_x, (0,1)=f_y, (2,0)=f_xx, (0,2)=f_yy and (1,1)=f_xy!\n";
    throw SPLINE_ERROR;
  }

  // If out of bounds show a warning and set x,y to boundary value
  out_of_bounds_check(x,y);
  double xx = x, yy = y;
  if(x < xmin) xx = xmin; 
  if(x > xmax) xx = xmax; 
  if(y < ymin) yy = ymin; 
  if(y > ymax) yy = ymax;

#ifdef _USEOPENMP
  gsl_interp_accel *xacc_thread = xaccs[omp_get_thread_num()];
  gsl_interp_accel *yacc_thread = yaccs[omp_get_thread_num()];
#else
  gsl_interp_accel *xacc_thread = xacc;
  gsl_interp_accel *yacc_thread = yacc;
#endif

  return derivfunc[n](spline, xx, yy, xacc_thread, yacc_thread); 
}

//====================================================
// Free up memory
//====================================================
void Spline2D::free(){
  if(spline){

    // Reset class variables
    xmin = xmax = 0.0;
    ymin = ymax = 0.0;
    size_x = size_y = 0;
    dx_min = dx_max = 0.0;
    dy_min = dy_max = 0.0;

    // Free the spline data
    gsl_spline2d_free(spline);
    spline = nullptr;

    // Free accelerators
#ifdef _USEOPENMP
    for(auto & xa: xaccs)
      gsl_interp_accel_free(xa);
    std::vector<gsl_interp_accel *>().swap(xaccs);
    for(auto & ya: yaccs)
      gsl_interp_accel_free(ya);
    std::vector<gsl_interp_accel *>().swap(yaccs);
#else
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
    xacc = nullptr;
    yacc = nullptr;
#endif
  }
}

//====================================================
// Destructor
//====================================================
Spline2D::~Spline2D(){
  free();
}

//====================================================
// Rest of the class methods
//====================================================

void Spline2D::out_of_bounds_check(const double x, const double y) const{
  if(out_of_bounds_warning){
    if(x < xmin - dx_min || x > xmax + dx_max){
      std::cout << "Warning Spline2D[" << name << "] ";
      std::cout << "x = " << x << " is out of bounds (" << xmin << "," << xmax << ")\n";
    }
    if(y < ymin - dy_min || y > ymax + dy_max){
      std::cout << "Warning Spline2D[" << name << "] ";
      std::cout <<"y = " << y << " is out of bounds (" << ymin << "," << ymax << ")\n";
    }
  }
}
double Spline2D::deriv_x(const double x, const double y) const{
  return eval_deriv(x,y,1,0);
}
double Spline2D::deriv_xx(const double x, const double y) const{
  return eval_deriv(x,y,2,0);
}
double Spline2D::deriv_y(const double x, const double y) const{
  return eval_deriv(x,y,0,1);
}
double Spline2D::deriv_yy(const double x, const double y) const{
  return eval_deriv(x,y,0,2);
}
double Spline2D::deriv_xy(const double x, const double y) const{
  return eval_deriv(x,y,1,1);
}
std::pair<double,double> Spline2D::get_xrange() const{
  return {xmin, xmax};
}
std::pair<double,double> Spline2D::get_yrange() const{
  return {ymin, ymax};
}
std::string Spline2D::get_name() const{
  return name;
}
void Spline2D::set_out_of_bounds_warning(bool v){
  out_of_bounds_warning = v;
}

//====================================================
//====================================================
// Testing of Spline
//====================================================
//====================================================

double test_function(double x){
  return x*x*x;
}
double test_function_deriv_x(double x){
  return 3*x*x;
}
double test_function_deriv_xx(double x){
  return 6*x;
}
void test_Spline(){
  std::cout << "Running tests for Spline...\n";
  const double epsilon = 1e-10;
  const double xmin    = 0.0;
  const double xmax    = 2.0;
  const int    nx      = 100;

  // Create test-data
  Vector x(nx), y(nx);
  for(int ix = 0; ix < nx; ix++)
    x[ix] = xmin + (xmax-xmin) * ix / double(nx-1);
  for(int ix = 0; ix < nx; ix++)
    y[ix] = test_function(x[ix]);

  // Some random points to check it at (we dont check close to other end as the condition y''=0 is not satisfied there)
  Vector points(1000);
  for(auto &p : points)
    p = xmin + (xmax-xmin) * ( rand() % 10000 )/20000.0;

  // The function that does the testing
  auto run_test = [&](Spline &f){
    std::cout << "Running test for spline: " << f.get_name() << "\n";
    auto xrange = f.get_xrange();
    assert(fabs(xmin - xrange.first)  < epsilon);
    assert(fabs(xmax - xrange.second) < epsilon);
    for(auto & point: points){
      double xx    = point;
      double yy    = f(xx);
      double yy_x  = f.deriv_x (xx);
      double yy_xx = f.deriv_xx(xx);
      double error          = fabs(yy    - test_function(xx));
      double error_deriv_x  = fabs(yy_x  - test_function_deriv_x (xx));
      double error_deriv_xx = fabs(yy_xx - test_function_deriv_xx(xx));
      assert(error          < epsilon);
      assert(error_deriv_x  < epsilon);
      assert(error_deriv_xx < epsilon);
    }
  };

  // Test that creating a spline works
  Spline f(x,y,"test_spline");
  run_test(f);

  // Check operator= works properly
  Spline g;
  g = f;
  run_test(g);

  // Check copy constructor works properly
  Spline h(g);
  run_test(h);

  // Evaluate out of bounds and check that we get the same as the end-points
  std::cout << "Test out of bounds:\n";
  assert( fabs(h(xmin-1.0) - h(xmin)) < epsilon );
  assert( fabs(h(xmax+1.0) - h(xmax)) < epsilon );

  // Try to free a spline and then use it
  try {
    std::cout << "Trying to use a spline that does not exist to see if an error is thrown:\n";
    h.free();
    h(0.0);
  } catch(const char * error) {
    std::cout << "ErrorMessage = " << error;
  }
  
  // Create a spline using new, do a lookup and delete it
  Spline *s = new Spline(x,y,"pointerspline");
  assert( fabs(s->eval(xmin) - test_function(xmin)) < epsilon );
  delete s;
  
  // Try to create a spline with wrong size
  auto xnew = x;
  xnew.push_back(xmax+1.0);
  try {
    std::cout << "Trying to create a spline from a wrongly sized vector:\n";
    Spline q(xnew,y);
  } catch(const char * error) {
    std::cout << "ErrorMessage = " << error;
  }

  std::cout << "End testing!\n";
}

//====================================================
//====================================================
// Testing of Spline2D
//====================================================
//====================================================

double test_function_2D(double x, double y){
  return x*x*x + y*y*y + x*y;
}
double test_function_2D_deriv_x(double x, double y){
  return 3*x*x + y;
}
double test_function_2D_deriv_y(double x, double y){
  return 3*y*y + x;
}
double test_function_2D_deriv_xx(double x, double y){
  return 6*x;
}
double test_function_2D_deriv_yy(double x, double y){
  return 6*y;
}
double test_function_2D_deriv_xy(double x, double y){
  return 1.0;
}

void test_Spline2D(){
  std::cout << "Running tests for Spline2D...\n";
  const double epsilon = 1e-10;
  const double xmin    = 0.0;
  const double xmax    = 2.0;
  const double ymin    = 0.0;
  const double ymax    = 2.0;
  const int    nx      = 100;
  const int    ny      = 100;
  
  // Create test-data
  Vector x(nx), y(ny), z(nx*ny);
  for(int ix = 0; ix < nx; ix++)
    x[ix] = xmin + (xmax-xmin) * ix / double(nx-1);
  for(int iy = 0; iy < ny; iy++)
    y[iy] = ymin + (ymax-ymin) * iy / double(ny-1);
  for(int ix = 0; ix < nx; ix++)
    for(int iy = 0; iy < ny; iy++)
      z[ix + nx*iy] = test_function_2D(x[ix], y[iy]);
  
  // Some random points to check it at (we dont check close to other end as the condition y''=0 is not satisfied there)
  std::vector<std::pair<double,double>> points(1000);
  for(auto &p : points){
    p.first  = xmin + (xmax-xmin) * ( rand() % 10000 )/20000.0;
    p.second = ymin + (ymax-ymin) * ( rand() % 10000 )/20000.0;
  }

  // The function that does the testing
  auto run_test = [&](Spline2D &f){
    std::cout << "Running test for spline: " << f.get_name() << "\n";
    auto xrange = f.get_xrange();
    auto yrange = f.get_yrange();
    assert(fabs(xmin - xrange.first)  < epsilon);
    assert(fabs(xmax - xrange.second) < epsilon);
    assert(fabs(ymin - yrange.first)  < epsilon);
    assert(fabs(ymax - yrange.second) < epsilon);
    for(auto & point: points){
      double xx = point.first;
      double yy = point.second;
      double zz    = f(xx,yy);
      double zz_x  = f.deriv_x (xx,yy);
      double zz_y  = f.deriv_y (xx,yy);
      double zz_xx = f.deriv_xx(xx,yy);
      double zz_yy = f.deriv_yy(xx,yy);
      double zz_xy = f.deriv_xy(xx,yy);

      double error          = fabs(zz    - test_function_2D(xx,yy));
      double error_deriv_x  = fabs(zz_x  - test_function_2D_deriv_x (xx, yy));
      double error_deriv_y  = fabs(zz_y  - test_function_2D_deriv_y (xx, yy));
      double error_deriv_xx = fabs(zz_xx - test_function_2D_deriv_xx(xx, yy));
      double error_deriv_yy = fabs(zz_yy - test_function_2D_deriv_yy(xx, yy));
      double error_deriv_xy = fabs(zz_xy - test_function_2D_deriv_xy(xx, yy));
      
      assert(error          < epsilon);
      assert(error_deriv_x  < epsilon);
      assert(error_deriv_y  < epsilon);
      assert(error_deriv_xx < epsilon);
      assert(error_deriv_yy < epsilon);
      assert(error_deriv_xy < epsilon);
    }
  };

  // Test that creating a spline works
  Spline2D f(x,y,z,"test_spline");
  run_test(f);

  // Check operator= works properly
  Spline2D g;
  g = f;
  run_test(g);

  // Check copy constructor works properly
  Spline2D h(g);
  run_test(h);
  
  // Evaluate out of bounds and check that we get the same as the end-points
  std::cout << "Test out of bounds:\n";
  assert( fabs(h(xmin    ,ymin-1.0) - h(xmin,ymin)) < epsilon );
  assert( fabs(h(xmin-1.0,ymin    ) - h(xmin,ymin)) < epsilon );
  assert( fabs(h(xmin-1.0,ymin-1.0) - h(xmin,ymin)) < epsilon );
  assert( fabs(h(xmax    ,ymax+1.0) - h(xmax,ymax)) < epsilon );
  assert( fabs(h(xmax+1.0,ymax    ) - h(xmax,ymax)) < epsilon );
  assert( fabs(h(xmax+1.0,ymax+1.0) - h(xmax,ymax)) < epsilon );

  // Try to free a spline and then use it
  try {
    std::cout << "Trying to use a spline that does not exist to see if an error is thrown:\n";
    h.free();
    h(0.0,0.0);
  } catch(const char * error) {
    std::cout << "ErrorMessage = " << error;
  }
  
  // Create a spline using new, do a lookup and delete it
  Spline2D *s = new Spline2D(x,y,z,"pointerspline");
  assert( fabs(s->eval(xmin,ymin) - test_function_2D(xmin,ymin)) < epsilon );
  delete s;

  // Try to create a spline with wrong size
  auto xnew = x;
  xnew.push_back(xmax+1.0);
  try {
    std::cout << "Trying to create a spline from a wrongly sized vector:\n";
    Spline2D q(xnew,y,z);
  } catch(const char * error) {
    std::cout << "ErrorMessage = " << error;
  }
  
  std::cout << "End testing!\n";
}

