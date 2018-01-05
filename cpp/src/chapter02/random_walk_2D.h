#ifndef CPP_RANDOM_WALK_2D_H
#define CPP_RANDOM_WALK_2D_H

/**
 * @file
 *
 * @author Hendrix Demers <hendrix.demers@mail.mcgill.ca>
 * @copyright 2018
 */

// Forwards declarations
// C system headers
// C++ system header
// Library headers
// Project headers
// Project private headers

struct Point
{
    double x;
    double y;
    double cx;
    double cy;
};

int set_up_screen();
void initialize();
Point initialize_coordinates();
Point new_coordinates(const double theta_rad, const Point current_point, const double step);
void plot_xy();
Point reset_coordinates(const double xn, const double yn, const double CA, const double CB);
void how_far(const Point current_point, const double step);


#endif //CPP_RANDOM_WALK_2D_H
