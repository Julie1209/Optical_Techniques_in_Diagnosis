# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:26:39 2023

@author: cdpss
"""

#!/usr/bin/env python
# coding: utf-8

import numpy as np
import module_singleLayer
import matplotlib.pyplot as plt
from tqdm import tqdm


# 1/cm
mua = 10
mus = 90

g = 0

n_air=1
n_tissue=1.5

absorbed_prob = mua/(mua+mus)
scattered_prob = 1 - absorbed_prob

# photon disappearing threshold
disappear_threshold = 0.0001

# photon survived probability (playing roulette)
photon_survived_prob = 1/20

# cm
medium_thickness = 0.2

# photon number in one run
photon_num = 10000

# running times
run_times = 5

# Define and initialize roulette_win_times
roulette_win_times = 0  

# In[13]:


for run_index in range(run_times):

    ####### run simulation

    A = 0
    R = 0

    for photon_index in tqdm(range(photon_num), desc="Run {} ".format(run_index+1)):

        new_photon = module_singleLayer.photon(variable_weight=True)

        new_photon.travel(mua, mus)
        
        transmitted_percentage, reflected_percentage= new_photon.get_boundary_split_ratio(n1=n_air,n2=n_tissue)
        R=R+new_photon.weight*reflected_percentage
        new_photon.weight=new_photon.weight*transmitted_percentage
        
        while True:
            if new_photon.position[2] > 0:        
                A = A + new_photon.weight*absorbed_prob
                new_photon.weight = new_photon.weight*scattered_prob

                if new_photon.weight > disappear_threshold:
                    new_photon.scattering(g)
                    new_photon.travel(mua, mus)
                else:
                    new_photon.play_roulette(photon_survived_prob)
                    if new_photon.weight == 0:
                        break
                    else:
                        roulette_win_times += 1
                        new_photon.scattering(g)
                        new_photon.travel(mua, mus)
            else:
                transmitted_percentage, reflected_percentage= new_photon.get_boundary_split_ratio(n1=n_tissue,n2=n_air)
                print(transmitted_percentage)
                    
                R = R + new_photon.weight*transmitted_percentage
                new_photon.weight = new_photon.weight*reflected_percentage
                new_photon.reflect_updatePosDir()

    # print("roulette_win_times:", roulette_win_times)
    print("roulette_win_times:",roulette_win_times)
    print("R:", R/photon_num)
    print("A:", A/photon_num)
    print("R+A:", (R + A)/photon_num)
    
    # plot the graph
    plt.pie([R, A], labels=["R ", "A"], autopct='%1.3f%%')
    title = "Isotropic scattering(g={}), photon_num={}, run_{}".format(g, photon_num, run_index+1)
    plt.title(title)
    plt.savefig("{}.png".format(title), bbox_inches = "tight")
    plt.show()





