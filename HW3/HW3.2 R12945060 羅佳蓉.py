#!/usr/bin/env python
# coding: utf-8

import numpy as np
import module
import matplotlib.pyplot as plt
from tqdm import tqdm

####### parameter setting

# 1/cm
mua = 10
mus = 90

# anisotropic factor
g = 0.75

# absorbed and scattered probability
absorbed_prob = mua/(mua+mus)
scattered_prob = 1 - absorbed_prob

# photon disappearing threshold
disappear_threshold = 0.001

# photon survived probability (playing roulette)
photon_survived_prob = 1/20

# cm
medium_thickness = 0.02

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
    T = 0
    R = 0

    for photon_index in tqdm(range(photon_num), desc="Run {} ".format(run_index+1)):

        # photon initializaion
        new_photon = module.photon(variable_weight=True)

        # go into the tissue
        new_photon.travel(mua, mus)

        # travel in the tissue
        while True:
            if 0 <= new_photon.position[2] < medium_thickness:        
                A = A + new_photon.weight*absorbed_prob
                new_photon.weight = new_photon.weight*scattered_prob

                if new_photon.weight > disappear_threshold:
                    new_photon.scattering(g)
                    new_photon.travel(mua, mus)
                else:
                    new_photon.play_roulette(photon_survived_prob)
                    if new_photon.weight == 0:
                        roulette_win_times += 1
                        break

            # exit tissue and update R, T
            elif new_photon.position[2] < 0:
                R += new_photon.weight
                break
            else:
                T += new_photon.weight
                break

    # print("roulette_win_times:", roulette_win_times)
    print("R:", R/photon_num)
    print("T:", T/photon_num)
    print("A:", A/photon_num)
    print("R+T+A:", (R + T + A)/photon_num)
    
    # plot the graph
    plt.pie([R, T, A], labels=["R ", "T", "A"], autopct='%1.3f%%')
    title = "Anisotropic scattering(g={}), photon_num={}, run_{}".format(g, photon_num, run_index+1)
    plt.title(title)
    plt.savefig("{}.png".format(title), bbox_inches = "tight")
    plt.show()





