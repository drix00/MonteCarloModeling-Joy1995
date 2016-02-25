# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 18:11:49 2015

@author: hdemers
"""
import logging
import math
import random
import csv

import matplotlib.pyplot as plt

two_pi = 2.0*math.pi
# Cutoff energy in keV in bulk case
e_min = 0.5

def single_scatter(inc_energy, at_num, at_wht, density, traj_num, thickness=None, debug=False):
    """
    A single scattering Monte Carlo simulation which uses the screened Rutherford cross section.
    """
    mn_ion_pot = compute_mean_ionization_potential(at_num)

    if thickness is None:
        thick = estimate_range(inc_energy, density)
        thin = False
    else:
        thick = thickness
        thin = True


    al_a, er, lam_a, sg_a = get_constants(at_num, inc_energy, at_wht, density)

    if debug:
        random.seed(-1)
    else:
        random.seed(random.randint(-1000000, 1000000))

    # The Monte Carlo loop.
    bk_sct, num = zero_counters()

    while num < traj_num:
        s_en, x, y, z, cx, cy, cz = reset_coordinates(inc_energy)

        # Allow initial entrance of electron.
        step = -lambda_mfp(s_en, al_a, sg_a, lam_a) * math.log(random.random())
        zn = step

        if zn > thick:
            num = straight_through(num)
            continue
        else:
            y = 0.0
            z = zn
            s_en = reset_next_step(step, s_en, density, mn_ion_pot, at_num, at_wht)

        # Start the single scattering loop.

        while s_en > e_min:
            step = -lambda_mfp(s_en, al_a, sg_a, lam_a) * math.log(random.random())
            cp, sp, ga = s_scatter(s_en, al_a)
            xn, yn, zn, ca, cb, cc = new_coord(step, cp, sp, ga, x, y, z, cx, cy, cz)

            if zn <= 0.0:
                num, bk_sct = back_scatter(num, bk_sct)
                break
            elif zn > thick:
                num, yn = transmit_electron(num, thick, y, z, cb, cc)
                break
            else:
                cx = ca
                cy = cb
                cz = cc

                x = xn
                y = yn
                z = zn
                s_en = reset_next_step(step, s_en, density, mn_ion_pot, at_num, at_wht)
        else:
            num = num + 1
            show_traj_num(num, traj_num, debug)

    # End of the Monte Carlo loop.
    bse_coefficient, bse_coefficient_error = show_BS_coeff(bk_sct, traj_num, debug)
    return bse_coefficient, bse_coefficient_error

def compute_mean_ionization_potential(at_num):
    """
    Calculate J the mean ionization potential mn_ion_pot using the Berger-Selzer analytical fit.
    """
    mm_ion_pot = (9.76 * at_num + (58.5 / math.pow(at_num, 0.19)))*0.001
    return mm_ion_pot

def estimate_range(inc_energy, density):
    """
    Estimate the beam range.
    """
    thick = 700.0 * math.pow(inc_energy, 1.66) / density
    if thick < 1000.0:
        thick = 1000.0

    return thick

def stop_pwr(energy, mn_ion_pot, at_num, at_wht):
    """
    This computes the stopping power in keV/g/cm2 using the modified Bethe experssion of Eq. (3.21).
    """
    # Just in case
    if energy < 0.05:
        logging.warning("energy %f < 0.05", energy)
        energy = 0.05

    factor = math.log(1.166*(energy + 0.85*mn_ion_pot) / mn_ion_pot)
    stop_power = factor * 78500.0 * at_num / (at_wht * energy)

    return stop_power

def lambda_mfp(energy, al_a, sg_a, lam_a):
    """
    Compute elastic MGP for single scattering model.
    """
    al = al_a / energy
    ak = al * (1.0 + al)

    # Giving the sg scross section in cm2 as
    sg = sg_a / (energy * energy * ak)
    # and lambda in angstroms is
    lambda_mfp_A = lam_a/sg

    return lambda_mfp_A

def get_constants(at_num, inc_energy, at_wht, density):
    """
    Computes some constants needed by the program.
    """
    al_a = math.pow(at_num, 0.67)*3.43e-3

    # Relativistically correct the beam energy for up to 500 keV.
    er = (inc_energy + 511.0) / (inc_energy + 1022.0)
    er = er * er
    # lambda in cm.
    lam_a = at_wht / (density*6.0e23)
    # Put into angstroms
    lam_a = lam_a * 1.0e8
    sg_a = at_num * at_num * 12.56 * 5.21e-21*er

    return al_a, er, lam_a, sg_a

def reset_coordinates(inc_energy):
    """
    Reset coordinates at start of each trajectory.
    """
    s_en = inc_energy
    x = 0.0
    y = 0.0
    z = 0.0

    cx = 0.0
    cy = 0.0
    cz = 1.0

    return s_en, x, y, z, cx, cy, cz

def zero_counters():
    bk_sct = 0
    num = 0

    return bk_sct, num

def s_scatter(energy, al_a):
    """
    Calculates scattering angle using screened Rutherford cross section.
    """
    al = al_a / energy
    R1 = random.random()
    cp = 1.0 - ((2.0 * al *R1) / (1.0 + al - R1))
    sp = math.sqrt(1.0 - cp*cp)
    # And get the azimuthal scattering angle.
    ga = two_pi * random.random()

    return cp, sp, ga

def new_coord(step, cp, sp, ga, x, y, z, cx, cy, cz):
    """
    Get xn, yn, zn from x, y, z and scattering angles.
    """
    # Find the transformation angles.
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

def straight_through(num):
    """
    Handles caes where initial entry exceeds thickness.
    """
    num = num + 1

    return num

def back_scatter(num, bk_sct):
    """
    Handles case of backscattered electrons.
    """
    num = num + 1
    bk_sct = bk_sct + 1

    return num, bk_sct

def transmit_electron(num, thick, y, z, cb, cc):
    """
    Handles case of transmitted electron.
    """
    num = num + 1
    # Length of path from z to bottom face.
    ll = (thick - z) / cc
    # hence the exit y-coordinate
    yn = y + ll * cb

    return num, yn

def reset_next_step(step, s_en, density, mn_ion_pot, at_num, at_wht):
    """
    Resets variables for next trajectory step.
    """
    # Find the energy loss on thos step.
    del_E = step * stop_pwr(s_en, mn_ion_pot, at_num, at_wht) * density * 1.0e-8
    # So the current energy is.
    s_en = s_en - del_E

    return s_en

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

    traj_num = 10000

    bse_coefficient, bse_coefficient_error = single_scatter(inc_energy, at_num, at_wht, density, traj_num, debug=True)

    print("Z = %3i" % at_num)
    print("BSE book = %f" % bse_book)
    print("BSE MC   = %f +- %f" % (bse_coefficient, 3.0*bse_coefficient_error))
    print("Diff BSE = %f" % abs(bse_book - bse_coefficient))

def run_figure_6_1():
    print("Start figure 6.1 simulation")

    inc_energy = 10.0
    traj_num = 10000

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

        bse_coefficient, bse_coefficient_error = single_scatter(inc_energy, at_num, at_wht, density, traj_num)

        print("Z = %3i" % at_num)
        print("BSE book = %f" % bse_book)
        print("BSE MC   = %f +- %f" % (bse_coefficient, 3.0*bse_coefficient_error))
        print("Diff BSE = %f" % abs(bse_book - bse_coefficient))
        bse_results.append(bse_coefficient)

    line = ""
    for bse in bse_results:
        line += "%f " % (bse)
    print(line)

def run_figure_6_5():
    print("Start figure 6.5 simulation")

    traj_num = 100000

    elements = []
    elements.append((6, 12.011, 2.62, 0.0625724721707))
    elements.append((13, 26.98, 2.7, 0.148075863868))
    elements.append((14, 28.09, 2.33, 0.170487882653))
    elements.append((26, 55.85, 7.86, 0.275218141234))
    elements.append((29, 63.55, 8.96, 0.29910424397))
    elements.append((47, 107.87, 10.5, 0.421003159787))
    elements.append((79, 196.97, 19.3, 0.520543686224))

    with open("figure_6.5_HD_100ke.csv", 'w', newline='\n') as bse_results_file:
        bse_writer = csv.writer(bse_results_file)
        row = ["energy (kV)"]
        for element in elements:
            row.append(element[0])
        bse_writer.writerow(row)

        for energy_keV in [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60]:
            row = [energy_keV]
            for element in elements:
                at_num, at_wht, density, bse_book = element

                bse_coefficient, bse_coefficient_error = single_scatter(energy_keV, at_num, at_wht, density, traj_num)

                print("%3i %2i %.4f %.4f" % (energy_keV, at_num, bse_coefficient, bse_coefficient_error))
                row.append(bse_coefficient)

            bse_writer.writerow(row)
            bse_results_file.flush()

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)

    #run_carbon()
    #run_figure_6_1()
    run_figure_6_5()

    plt.show()
