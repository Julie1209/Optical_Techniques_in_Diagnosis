#!/usr/bin/env python
# coding: utf-8

# In[13]:


import numpy as np
import module
import matplotlib.pyplot as plt
from tqdm import tqdm


# In[23]:


####### parameter setting

# 1/cm
mua = 10
mus = 90

# absorbed probability
absorbed_prob = mua/(mua+mus)

# cm
medium_thickness = 0.02

# photon number in one run
photon_num = 10000

# running times
run_times = 5


# In[35]:


for run_index in range(run_times):

    ####### run simulation

    absorbed_photons = []
    transmitted_photons = []
    reflected_photons = []

    photon_tracks = {}

    for photon_index in tqdm(range(photon_num), desc="Run {} ".format(run_index+1)):

        # create list for recording each photon position
        photon_y = np.array([])
        photon_z = np.array([])

        # photon initializaion
        new_photon = module.photon()
        photon_y = np.append(photon_y, new_photon.position[1])
        photon_z = np.append(photon_z, new_photon.position[2])

        # go into the tissue
        new_photon.travel(mua, mus)
        photon_y = np.append(photon_y, new_photon.position[1])
        photon_z = np.append(photon_z, new_photon.position[2])

        # travel in the tissue
        while True:
            if 0 <= new_photon.position[2] < medium_thickness:        
                rnd = np.random.uniform()
                if rnd < absorbed_prob:
                    absorbed_photons.append(photon_index)
                    break
                else:
                    new_photon.scattering()
                    new_photon.travel(mua, mus)
                    photon_y = np.append(photon_y, new_photon.position[1])
                    photon_z = np.append(photon_z, new_photon.position[2])

            # exit tissue and update R, T
            elif new_photon.position[2] < 0:
                reflected_photons.append(photon_index)
                break
            else:
                transmitted_photons.append(photon_index)
                break

        photon_tracks[photon_index] = [photon_y, photon_z]

    print("reflectance: {},  theoretical value: 0.3616".format(len(reflected_photons)/photon_num))
    print("transmittance: {},  theoretical value: 0.3565".format(len(transmitted_photons)/photon_num))
    print("absorbance: {}".format(len(absorbed_photons)/photon_num))
    print("reflectance + transmittance + absorbance:", len(reflected_photons+transmitted_photons+absorbed_photons)/photon_num)
    
    # plot the graph
    plt.pie([len(reflected_photons), len(transmitted_photons), len(absorbed_photons)], 
            labels=["R ", " T", "A"], autopct='%1.2f%%')
    title = "Isotropic scattering, photon_num={}, run_{}".format(photon_num, run_index+1)
    plt.title(title)
    plt.tight_layout()
    plt.savefig("{}.png".format(title), bbox_inches = "tight")
    plt.show()




