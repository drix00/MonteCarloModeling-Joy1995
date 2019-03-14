/**
 * @file
 *
 * This stimulates a simple random walk with equal-length steps.
 *
 * @author Hendrix Demers <hendrix.demers@mail.mcgill.ca>
 * @copyright 2018
 */

//   Copyright 2017 Hendrix Demers
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

// C system headers
// C++ system header
#include <iostream>
#include <cmath>
#include <cstdlib>
// Library headers
// Precompiled header
#pragma hdrstop
// Current declaration header file of this implementation file.
#include "random_walk_2D.h"
// Project headers
// Project private headers

// Global and constant variables/functions.

const double two_pi = 6.28318;

/**
 * Gets the required input data to run the simulation.
 */
int set_up_screen()
{
    std::cout << "Random walk simulation" << std::endl;
    std::cout << "How many steps?" << std::endl;
    int number_steps = 0;
    std::cin >> number_steps;

    return number_steps;
}

/**
 * Identify which graphics cards is in use and initialize it.
 */
void initialize()
{

}

/**
 * Set up the starting values of all the parameters.
 */
Point initialize_coordinates()
{
    Point point;
    point.x = 0;
    point.y = 0;
    point.cx = 1.0;
    point.cy = 0.0;

//    randomize();

    return point;
}

/**
 * Compute the new coordinates xn, yn, given x, y, cx, cy and theta.
 *
 * @return new coordinates xn, yn.
 */
Point new_coordinates(const double theta_rad, const Point current_point, const double step)
{
    const double V1 = std::cos(theta_rad);
    const double V2 = std::sin(theta_rad);

    ;
    // New direction cosines.
    const double CA = current_point.cx * V1 - current_point.cy * V2;
    const double CB = current_point.cy * V1 + current_point.cx * V2;

    // New coordinates
    const double xn = current_point.x + step * CA;
    const double yn = current_point.y + step * CB;

    const Point new_point = reset_coordinates(xn, yn, CA, CB);

    return new_point;
}

/**
 * Plots the step on the screen. Since all screen coordinates are integers, this conversion is made first.
 * The real X, Y coordinates are separated from the plotting coordinates, which put X=Y=0 at the point hstart, vstart.
 */
void plot_xy()
{

}

/**
 * Shifts coordinate reference X, Y to new coordinates XN, Yn and resets the direction cosines.
 */
Point reset_coordinates(const double xn, const double yn, const double CA, const double CB)
{
    Point new_point;

    new_point.x = xn;
    new_point.y = yn;
    new_point.cx = CA;
    new_point.cy = CB;

    return new_point;
}

/**
 * Computes the distance traveled in the walk.
 */
void how_far(const Point current_point, const double step)
{
    double distance = std::sqrt((current_point.x - 0.0)*(current_point.x - 0.0) + (current_point.y - 0.0)*(current_point.y - 0.0));
    distance /= step;

    std::cout << "The walker traveled " << distance << " steps from origin." << std::endl;
}


int main() {
    int number_steps = set_up_screen();

    initialize();

    const double step = 1.0;

    Point point = initialize_coordinates();

    for (int step_id=0; step_id < number_steps; ++step_id)
    {
        const double theta_rad = two_pi * static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);

        point = new_coordinates(theta_rad, point, step);

        plot_xy();

    }

    how_far(point, step);

    return 0;
}