# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 18:11:49 2015

@author: hdemers
"""
import logging
import math
import random

import matplotlib.pyplot as plt

two_pi = 2.0*math.pi
# Cutoff energy in keV in bulk case
e_min = 0.5

def plurial_scatter(inc_energy, at_num, at_wht, density, tilt_deg, traj_num, debug=False):
    """
    THis program performs a Monte Carlo trajectory simulation using a screened Rutherford cross section and 
    plural scattering approximation.
    """    
    logging.info("Plural scattering Monte Carlo simulation")
    
    mn_ion_pot = compute_mean_ionization_potential(at_num)
    
    # 57.4 convert degrees to radians.
    s_tilt = math.sin(tilt_deg/57.4)
    c_tilt = math.cos(tilt_deg/57.4)
    
    p0, rf = rutherford_factor(at_num, inc_energy)
    
    # Now get the range and step length in microns for these conditions.
    m_t_step, step = compute_bethe_range(inc_energy, mn_ion_pot, at_num, at_wht, density)
    E, e_smooth = profile(inc_energy, m_t_step, mn_ion_pot, at_num, at_wht)
    
    if debug:
        random.seed(-1)
    else:
        random.seed(random.randint(-1000000, 1000000))    
        
    # Monte Carlo loop.
    bk_sct, num = init_counters()
    
    while num < traj_num:
        x, y, z, cx, cy, cz = reset_coordinates(s_tilt, c_tilt)
        for k in range(0, 49):
            sp, cp, ga = p_scatter(rf, E, k)
            xn, yn, zn, ca, cb, cc = new_coord(step, cp, sp, ga, x, y, z, cx, cy, cz)
            
            if zn <= 0.0:
                num, bk_sct = back_scatter(num, bk_sct)
                break
            else:
                cx = ca
                cy = cb
                cz = cc
                
                x = xn
                y = yn
                z = zn
        else:
            num = num + 1
            
        show_traj_num(num, traj_num, debug)
        
    # End of the Monte Carlo loop.
    bse_coefficient, bse_coefficient_error = show_BS_coeff(bk_sct, traj_num, debug)
    return bse_coefficient, bse_coefficient_error
        
def rutherford_factor(at_num, inc_energy):
    """
    Computes the Rutherford scattering factor for this incident energy.
    """    
    # Duncumb's function.
    p0 = 0.394 * math.pow(at_num, 0.4) / inc_energy
    
    # Scattering constants b/2 in Eq. (4.3).
    rf = 0.0072 * at_num / p0
    
    return p0, rf
        
def compute_bethe_range(inc_energy, mn_ion_pot, at_num, at_wht, density):
    """
    This calulates the range in microns assuming the modified Bethe stopping power Eq. (3.21).
    """
    fs = 0.0
    
    # a Simpson's rule integration.
    # 20 equal steps.
    for m in range(1, 21):
        energy = (m - 1.0) * inc_energy / 20.0
        f = 1.0 / stop_pwr(energy, mn_ion_pot, at_num, at_wht)
        
        l = 2
        if m % 2 == 0:
            l = 4
        if m == 1 or m == 21:
            l = 1
            
        fs = fs + l * f
        
    # Now use this to find the range and step length for these conditions.
    # in g/cm2
    bethe_range = fs * inc_energy / 60.0
    m_t_step = bethe_range / 50.0
    
    # in microns
    bethe_range = bethe_range * 10000.0 / density
    logging.info("Range is %.f microns", bethe_range)
    
    step = bethe_range / 50.0
    
    return m_t_step, step
    
def profile(inc_energy, m_t_step, mn_ion_pot, at_num, at_wht):
    """
    Compute 50-step energy profile for electrons beam.
    """
    E = []
    E.append(inc_energy)
    
    for m in range(2, 51):
        A1 = m_t_step * stop_pwr(E[-1], mn_ion_pot, at_num, at_wht)
        A2 = m_t_step * stop_pwr(E[-1] - A1/2.0, mn_ion_pot, at_num, at_wht)
        A3 = m_t_step * stop_pwr(E[-1] - A2/2.0, mn_ion_pot, at_num, at_wht)
        A4 = m_t_step * stop_pwr(E[-1] - A3, mn_ion_pot, at_num, at_wht)
        
        E.append(E[-1] - (A1 + 2.0 * A2 + 2.0 * A3 + A4) / 6.0)
            
    E.append(0.0)
    assert len(E) == 51
    
    e_smooth = []
    # A little smoothing of the profile.    
    for m in range(1, 49):
        e_smooth.append((E[m] + E[m+1])/2.0)
        
    assert len(E) == 51
    return E, e_smooth
    
def compute_mean_ionization_potential(at_num):
    """
    Calculate the mean ionization potential J using the Berger-Selzer analytical fit in units of keV.
    """
    mm_ion_pot = (9.76 * at_num + (58.5 / math.pow(at_num, 0.19)))*0.001
    return mm_ion_pot
    
def stop_pwr(energy, mn_ion_pot, at_num, at_wht):
    """
    Calculates the stopping power using the modified Bethe expression of Eq. (3.21) in units of keV/g/cm2.
    """
    # Just in case
    if energy < 0.05:
        logging.warning("energy %f < 0.05", energy)
        energy = 0.05
        
    factor = math.log(1.166*(energy + 0.85*mn_ion_pot) / mn_ion_pot)
    stop_power = factor * 78500.0 * at_num / (at_wht * energy)
    
    return stop_power
    
def init_counters():
    bk_sct = 0
    num = 0
    
    return bk_sct, num
    
def reset_coordinates(s_tilt, c_tilt):
    """
    Reinitialize all the elctron variables at start of each new trajectory.
    """
    x = 0.0
    y = 0.0
    z = 0.0
    
    cx = 0.0
    cy = s_tilt
    cz = c_tilt
    
    return x, y, z, cx, cy, cz
    
def p_scatter(rf, E, k):
    """
    Calculates scattering angles using the plural scattering model with the small angle correction as given in Eq. (4.10).
    """
    # First call the random number generator function.
    nu = math.sqrt(random.random())
    nu = 1.0/nu - 1.0
    an = nu * rf / E[k]
    
    # and use this to find the scattering angles.
    sp = (an + an) / (1.0 + (an * an))
    cp = (1.0 - (an*an)) / (1.0 + (an * an))
    
    # and the azimuthal scattering angle.
    ga = two_pi * random.random()
    
    return sp, cp, ga
    
def new_coord(step, cp, sp, ga, x, y, z, cx, cy, cz):
    """
    Gets xn, yn, zn from x, y, z and scattering angles using Eqs. (3.12) to (3.15).
    """
    # The coordinate rotation angles are.
    if cz == 0:
        cz = 0.000001
        
    an_m = (-cx / cz)
    an_n = 1.0 / math.sqrt(1.0 + (an_m * an_m))
    
    # Save computation time by getting all the trancendentals first.
    v1 = an_n * sp
    v2 = an_m * an_n * sp
    v3 = math.cos(ga)
    v4 = math.sin(ga)

    # Find the new direction cosines.
    ca = (cx * cp) + (v1 * v3) + (cy * v2 * v4)
    cb = (cy * cp) + (v4 * (cz * v1 - cx * v2))
    cc = (cz * cp) + (v2 * v3) - (cy * v1 * v4)

    # and get the new coordinates.
    xn = x + step * ca
    yn = y + step * cb
    zn = z + step * cc
    
    return xn, yn, zn, ca, cb, cc
    
def back_scatter(num, bk_sct):
    """
    Handles case of backscattered electrons.
    """
    num = num + 1
    bk_sct = bk_sct + 1
    
    return num, bk_sct
    
def show_traj_num(num, traj_num, debug):
    """
    Updates thermometer display for % of trajectories done.
    """
    
    percent_done = num / float(traj_num) * 100.0
    if debug and percent_done % 5 ==0:
        print("%3i %%" % percent_done)
        
def show_BS_coeff(bk_sct, traj_num, debug):
    """
    Displays BS coefficient on thermometer scale.
    """
    bse_coefficient = bk_sct / float(traj_num)
    bse_coefficient_error = math.sqrt(bse_coefficient / traj_num * (1.0 - bse_coefficient))
    
    if debug:
        print("BSE: %f +- %f" % (bse_coefficient, 3.0*bse_coefficient_error))
    
    return bse_coefficient, bse_coefficient_error
    
def run_carbon():
    inc_energy = 10.0
    at_num = 6
    at_wht = 12.011
    density = 2.62
    bse_book = 0.0625724721707
    
    tilt_deg = 0.0
    traj_num = 1000
    
    bse_coefficient, bse_coefficient_error = plurial_scatter(inc_energy, at_num, at_wht, density, tilt_deg, traj_num, debug=True)
    
    print("Z = %3i" % at_num)
    print("BSE book = %f" % bse_book)
    print("BSE MC   = %f +- %f" % (bse_coefficient, 3.0*bse_coefficient_error))
    print("Diff BSE = %f" % abs(bse_book - bse_coefficient))
    
def run_figure_6_1():
    print("Start figure 6.1 simulation")
    
    inc_energy = 10.0
    tilt_deg = 0.0
    traj_num = 100000
    
    elements = []
    elements.append((6, 12.011, 2.62, 0.0625724721707))
    elements.append((13, 26.98, 2.7, 0.148075863868))
    elements.append((14, 28.09, 2.33, 0.170487882653))
    elements.append((26, 55.85, 7.86, 0.275218141234))
    elements.append((29, 63.55, 8.96, 0.29910424397))
    elements.append((47, 107.87, 10.5, 0.421003159787))
    elements.append((79, 196.97, 19.3, 0.520543686224))
    
    bse_results = []
    for element in elements:
        at_num, at_wht, density, bse_book = element
        
        bse_coefficient, bse_coefficient_error = plurial_scatter(inc_energy, at_num, at_wht, density, tilt_deg, traj_num)
        
        print("Z = %3i" % at_num)
        print("BSE book = %f" % bse_book)
        print("BSE MC   = %f +- %f" % (bse_coefficient, 3.0*bse_coefficient_error))
        print("Diff BSE = %f" % abs(bse_book - bse_coefficient))
        bse_results.append(bse_coefficient)
        
    line = ""
    for bse in bse_results:
        line += "%f " % (bse)
    print(line)
    
if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    
    #run_carbon()
    run_figure_6_1()
    
    plt.show()
    