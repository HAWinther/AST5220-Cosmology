#ifndef _GSLSPLINE_HEADER
#define _GSLSPLINE_HEADER
#ifdef _USEOPENMP
#include <omp.h>
#endif
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

//====================================================
// 
// This is a wrapper class for easy use of GSL splines
// It is thread safe. The fiducial splines are
// cubic splines for 1D and bicubic for 2D assuming 
// natural (y''=0) boundary conditions at the end-points
//
// The fiducial choice allows to get the function value
// and up to second order derivatives. NB: if you change
// the fiducial choice to a lower order spline, 
// e.g to gsl_interp_linear, then the second derivative 
// deriv_xx etc. will not work
//
// If evaluating out of bounds it will return the closest
// value. To get a warning when this happens call
// set_out_of_bounds_warning(true) before using the spline
// or set the define _FIDUCIAL_SPLINE_WARNING to be true
//
//====================================================

#ifndef _FIDUCIAL_INTERPOL_TYPE
#define _FIDUCIAL_INTERPOL_TYPE    gsl_interp_cspline
#endif
#ifndef _FIDUCIAL_INTERPOL_TYPE_2D
#define _FIDUCIAL_INTERPOL_TYPE_2D gsl_interp2d_bicubic
#endif
#ifndef _FIDUCIAL_SPLINE_WARNING
#ifdef _SPLINE_WARNINGS_ON
#define _FIDUCIAL_SPLINE_WARNING   true
#else
#define _FIDUCIAL_SPLINE_WARNING   false
#endif
#endif

using Vector   = std::vector<double>;
using Vector2D = std::vector<Vector>;

//==============================================================================
//==============================================================================

class Spline{
  private:

    // GSL spline
    gsl_spline *spline = nullptr;

    // If we use threads then we must have a unique accelerator per thread
#ifdef _USEOPENMP
    std::vector<gsl_interp_accel *> xaccs;
#else
    gsl_interp_accel *xacc = nullptr;
#endif

    const gsl_interp_type *interpoltype_used = _FIDUCIAL_INTERPOL_TYPE;

    // Info about the spline
    int size_x       = 0;
    double xmin      = 0.0;
    double xmax      = 0.0;
    double dx_min    = 0.0;
    double dx_max    = 0.0;
    std::string name = "NoName";
  
    // Print warnings if out of bounds if wanted
    bool out_of_bounds_warning = _FIDUCIAL_SPLINE_WARNING;
    void out_of_bounds_check(const double x) const;

  public:
    
    // Constructors and destructor
    Spline() = default;
    Spline(std::string name);
    Spline(
        double *x, 
        double *y, 
        const int nx, 
        std::string splinename = "NoName", 
        const gsl_interp_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE);
    Spline(
        Vector& x, 
        Vector &y, 
        std::string splinename = "NoName", 
        const gsl_interp_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE);
    ~Spline();

    operator bool() const{
      return (spline != nullptr);
    }

    // Copy and assignment operator
    Spline(const Spline& rhs);
    Spline& operator=(const Spline& rhs);
    
    // Methods that create the spline
    void create(
        double *x, 
        double *y, 
        const int nx, 
        std::string splinename = "NoName", 
        const gsl_interp_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE);
    void create(
        Vector &x, 
        Vector &y, 
        std::string splinename = "NoName", 
        const gsl_interp_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE);
    
    // Methods for spline lookup of function and its derivatives
    double operator()(const double x) const;
    double eval(const double x) const;
    double eval_deriv(const double x, const int deriv) const;
    double deriv_x(const double x) const;
    double deriv_xx(const double x) const;
  
    // Some useful info
    std::pair<double,double> get_xrange() const;
    std::string get_name() const;
    void set_out_of_bounds_warning(bool v);
    
    // Clean up
    void free();
};

// Test the spline routines
void test_Spline();

//==============================================================================
//==============================================================================

class Spline2D{
  private:

    // GSL spline
    gsl_spline2d *spline   = nullptr;
    
    // If we use threads then we must have a unique accelerator per thread
#ifdef _USEOPENMP
    std::vector<gsl_interp_accel *> xaccs;
    std::vector<gsl_interp_accel *> yaccs;
#else
    gsl_interp_accel *xacc = nullptr;
    gsl_interp_accel *yacc = nullptr;
#endif
    const gsl_interp2d_type *interpoltype_used = _FIDUCIAL_INTERPOL_TYPE_2D;
    
    // Info about the spline
    int size_x       = 0;
    int size_y       = 0;
    double xmin      = 0.0;
    double xmax      = 0.0;
    double ymin      = 0.0;
    double ymax      = 0.0;
    double dx_min    = 0.0;
    double dx_max    = 0.0;
    double dy_min    = 0.0;
    double dy_max    = 0.0;
    std::string name = "NoName";
    
    // Print warnings if out of bounds if wanted
    bool out_of_bounds_warning = _FIDUCIAL_SPLINE_WARNING;
    void out_of_bounds_check(const double x, const double y) const;
  
    // A list of all the (up to second order) derivative functions in GSL,
    // We map Dx^nx Dy^ny f => nx + 3*ny in the list for use in eval_deriv
    typedef double (*evalfunc)(
        const gsl_spline2d *, 
        const double, 
        const double, 
        gsl_interp_accel *, gsl_interp_accel *);
    std::vector<evalfunc> derivfunc{
        gsl_spline2d_eval,
        gsl_spline2d_eval_deriv_x,
        gsl_spline2d_eval_deriv_xx,
        gsl_spline2d_eval_deriv_y,
        gsl_spline2d_eval_deriv_xy, 
        nullptr, 
        gsl_spline2d_eval_deriv_yy};

  public:

    // Constructors and destructor
    Spline2D() = default;
    Spline2D(std::string name);
    Spline2D(
        double *x, 
        double *y, 
        double *z, 
        const int nx, 
        const int ny, 
        std::string splinename = "NoName", 
        const gsl_interp2d_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE_2D);
    Spline2D(
        Vector& x, 
        Vector& y, 
        Vector& z, 
        std::string splinename = "NoName", 
        const gsl_interp2d_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE_2D);
    Spline2D(
        Vector& x, 
        Vector& y, 
        Vector2D& z, 
        std::string splinename = "NoName", 
        const gsl_interp2d_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE_2D);
    ~Spline2D();
    
    // Copy and assignment operator
    Spline2D(const Spline2D& rhs);
    Spline2D& operator=(const Spline2D& rhs);
    
    // Methods for creating the spline
    void create(
        double *x, 
        double *y, 
        double *z, 
        const int nx, 
        const int ny, 
        std::string splinename = "NoName", 
        const gsl_interp2d_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE_2D);
    void create(
        Vector& x, 
        Vector& y, 
        Vector& z, 
        std::string splinename = "NoName", 
        const gsl_interp2d_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE_2D);
    void create(
        Vector& x, 
        Vector& y, 
        Vector2D& z, 
        std::string splinename = "NoName", 
        const gsl_interp2d_type *interpoltype = _FIDUCIAL_INTERPOL_TYPE_2D);

    // Methods for spline lookup of function and its derivatives
    double operator()(const double x, const double y) const;
    double eval(const double x, const double y) const;
    double eval_deriv(const double x, const double y, const int derivx, const int derivy) const;
    double deriv_x(const double x, const double y) const;
    double deriv_xx(const double x, const double y) const;
    double deriv_xy(const double x, const double y) const;
    double deriv_y(const double x, const double y) const;
    double deriv_yy(const double x, const double y) const;
    
    // Some useful info
    std::pair<double,double> get_xrange() const;
    std::pair<double,double> get_yrange() const;
    std::string get_name() const;
    void set_out_of_bounds_warning(bool v);

    // Clean up
    void free();
};

// Test the spline routines
void test_Spline2D();

//==============================================================================
//==============================================================================

#endif
