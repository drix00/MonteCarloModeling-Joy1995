# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:50:35 2015

@author: hdemers
"""
import logging
import math
import random

import matplotlib.pyplot as plt

TWO_PI = 2.0*math.pi
step = 1.0
x_start = 0.0
y_start = 0.0


def random_walk(number_steps):
    """
    This stimulates a simple random walk with equal-length steps.
    """
    logging.debug("random_walk")
    xs = []
    ys = []
    
    x, y, CX, CY = initialize_coordinates()
    xs.append(x)
    ys.append(y)
    
    # Loop start here.
    for step_id in range(number_steps):
        theta_rad = TWO_PI*random.random()
        
        x, y, CX, CY = new_coord(x, y, CX, CY, theta_rad)
        xs.append(x)
        ys.append(y)
        
    distance = how_far(x, y)
    
    return distance, xs, ys


def initialize_coordinates():
    """
    Set up the starting values of all the parameters.
    """
    x = 0.0
    y = 0.0
    CX = 1.0
    CY = 0.0
    
    random.seed(random.randint(-100000, 100000))
    
    return x, y, CX, CY


def new_coord(x, y, CX, CY, theta_rad):
    """
    Compute the  new coordinates xn, yn, given x, y, CX, CY and theta.
    """
    
    V1 = math.cos(theta_rad)
    V2 = math.sin(theta_rad)
    
    # New direction cosines.
    CA = CX*V1 - CY*V2
    CB = CY*V1 + CX*V2
    
    xn = x + step*CA
    yn = y + step*CB
    
    return xn, yn, CA, CB


def how_far(x, y):
    """
    Computes the distance traveled in the walk.
    """
    distance = math.sqrt((x - x_start)*(x - x_start) + (y - y_start)*(y - y_start))
    distance = distance/step
    
    return distance


def run_figure2_4():
    maximum_number_steps = 800
    range_step = 100
    
    number_steps_list = range(range_step, maximum_number_steps+range_step, range_step)
    number_repetitions = 5
    
    means = []
    distances = {}
    for number_steps in number_steps_list:
        total = 0.0
        for iteration in range(number_repetitions):
            logging.debug("Number of steps: %i for iteration %i", number_steps, iteration+1)
            distance, _xs, _ys = random_walk(number_steps)
            distances.setdefault(iteration, []).append(distance)
            total += distance
        means.append(total/float(number_repetitions))
        
    plt.figure()
    
    for iteration in range(number_repetitions):
        plt.plot(number_steps_list, distances[iteration], '.k')
    
    plt.plot(number_steps_list, means, 'o')
    
    number_steps_list_theory = list(range(0, maximum_number_steps, 1))
    distance_theory = [math.sqrt(number_steps) for number_steps in number_steps_list_theory]
    plt.plot(number_steps_list_theory, distance_theory, '-')
    
    plt.xlabel('step')
    plt.ylabel('distance (steps)')
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)


def profile_figure2_4():
    maximum_number_steps = 800
    range_step = 100

    number_steps_list = range(range_step, maximum_number_steps + range_step, range_step)
    number_repetitions = 5

    means = []
    distances = {}
    for number_steps in number_steps_list:
        total = 0.0
        for iteration in range(number_repetitions):
            logging.debug("Number of steps: %i for iteration %i", number_steps, iteration + 1)
            distance, _xs, _ys = random_walk(number_steps)
            distances.setdefault(iteration, []).append(distance)
            total += distance
        means.append(total / float(number_repetitions))


def run_figure2_4_better():
    maximum_number_steps = 10000
    range_step = 100
    
    number_steps_list = range(range_step, maximum_number_steps+range_step, range_step)
    number_repetitions = 100
    
    means = []
    distances = {}
    for number_steps in number_steps_list:
        total = 0.0
        for iteration in range(number_repetitions):
            logging.debug("Number of steps: %i for iteration %i", number_steps, iteration+1)
            distance, _xs, _ys = random_walk(number_steps)
            distances.setdefault(iteration, []).append(distance)
            total += distance
        means.append(total/float(number_repetitions))
        
    plt.figure()
    
    for iteration in range(number_repetitions):
        plt.plot(number_steps_list, distances[iteration], '.k')
    
    plt.plot(number_steps_list, means, 'o')
    
    number_steps_list_theory = list(range(0, maximum_number_steps, 1))
    distance_theory = [math.sqrt(number_steps) for number_steps in number_steps_list_theory]
    plt.plot(number_steps_list_theory, distance_theory, '-')
    
    plt.xlabel('step')
    plt.ylabel('distance (steps)')
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)


def run_trajectories():
    number_steps = 200
    number_repetitions = 5
    for iteration in range(number_repetitions):
        distance, xs, ys = random_walk(number_steps)
        logging.info("Distance: %s", distance)
                
        plt.figure()
        plt.title("Iteration %i" % (iteration+1))
        
        plt.plot(xs, ys)
        
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.gca().set_aspect('equal', 'datalim')
        plt.grid(True)  


if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    
    # run_trajectories()
    # run_figure2_4()
    profile_figure2_4()
    # run_figure2_4_better()
    
    # plt.show()
