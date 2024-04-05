#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import module_singleLayer
import matplotlib.pyplot as plt
from tqdm import tqdm

####### parameter setting

# refractive index
n_air = 1
n_tissue = 1.37

# 1/cm
mua = 6
mus = 414

# anisotropic factor
g = 0.91

# absorbed and scattered probability
absorbed_prob = mua/(mua+mus)
scattered_prob = 1 - absorbed_prob

# photon disappearing threshold
disappear_threshold = 0.1

# photon survived probability (playing roulette)
photon_survived_prob = 1/20

# cm (tissue thickness is 1.5mm in question)
medium_thickness = 0.15

# grid structure
grid_width = 0.3    # cm
grid_depth = 0.15   # cm
delta_r = 0.01      # cm
delta_z = 0.01      # cm

# source type ("Pencil", "Uniform", "Gaussian")
source = "Gaussian"
beam_power = "7.85e-3" # W
distribution_r = 0.05 # cm


# photon number in one run
photon_num = 10000

# running times
run_times = 1


get_ipython().run_line_magic('matplotlib', 'inline')

R_all = np.empty(run_times)
T_all = np.empty(run_times)

for run_index in range(run_times):

    # initialize total absorption, total transmittance, total reflectance, absorption_matrix, roulette_win_times
    A = 0
    T = 0
    R = 0
    absorption_matrix = np.zeros(shape=(int(grid_width/delta_r), int(grid_depth/delta_z)))
    roulette_win_times = 0
    
    # start to run simulation !!
    for photon_index in tqdm(range(photon_num), desc="Run {} ".format(run_index+1)):

        # photon initializaion
        new_photon = module_singleLayer.photon(variable_weight=True, source_type=source, source_radius=distribution_r)

        # launch the photon (start from z = 0)
        new_photon.travel(mua, mus)
        
        # experience the first reflection
        transmitted_percentage, reflected_percentage = new_photon.get_boundary_split_ratio(n1=n_air, n2=n_tissue)
        R = R + new_photon.weight*reflected_percentage
        new_photon.weight = new_photon.weight*transmitted_percentage

        # go into the tissue
        photon_tissue_interaction_index = 0
        while True:
            
            # determine if photon is in tissue
            if 0 <= new_photon.position[2] < medium_thickness: # True --> simulation keep going (0大於等於還是只要大於就好呢？ 還是光子一開始的z直接設-0.0000000001??)
                
                # update W, A, absorption_matrix
#                 if photon_tissue_interaction_index != 0:
                photon_distance_to_original_axis = np.sqrt(new_photon.position[0]**2 + new_photon.position[1]**2)
                if photon_distance_to_original_axis < grid_width and new_photon.position[2] < grid_depth:
                    index_r = int(photon_distance_to_original_axis//delta_r)
                    index_z = int(new_photon.position[2]//delta_z)
                    absorption_matrix[index_r][index_z] = absorption_matrix[index_r][index_z] + new_photon.weight*absorbed_prob

                A = A + new_photon.weight*absorbed_prob
                new_photon.weight = new_photon.weight*scattered_prob
                photon_tissue_interaction_index += 1

                # check if weight is too small ?
                if new_photon.weight > disappear_threshold:
                    new_photon.scattering(g)
                    new_photon.travel(mua, mus)
                else:
                    new_photon.play_roulette(photon_survived_prob)
                    if new_photon.weight == 0:                        
                        break # -> next photon
                    else:
                        roulette_win_times += 1
                        new_photon.scattering(g)
                        new_photon.travel(mua, mus)                        

            else: # False --> photon leave tissue, do reflection (different from hw_3)
                
                transmitted_percentage, reflected_percentage = new_photon.get_boundary_split_ratio(n1=n_tissue, n2=n_air)
                
                if new_photon.position[2] < 0: # update R
                
                    # collect photon weight back to surface from tissue
                    R = R + new_photon.weight*transmitted_percentage
                    new_photon.weight = new_photon.weight*reflected_percentage

                    # reflect from tissue upper boundary to inner tissue (update position and direction)
                    new_photon.reflect_updatePosDir()
                    
                else:                          # update T
                    
                    # collect photon weight transmitting to lower air layer
                    T = T + new_photon.weight*transmitted_percentage
                    new_photon.weight = new_photon.weight*reflected_percentage
                    
                    # reflect from tissue lower boundary to inner tissue (update position and direction)
                    new_photon.transmit_updatePosDir(medium_thickness)

    print("roulette_win_times:", roulette_win_times)
    print("R: {:5f}, T: {:5f}, A: {:5f}".format(R/photon_num, T/photon_num, A/photon_num))
    print("R+T+A:", (R+T+A)/photon_num)
    print("Weight in absorption matrix:", absorption_matrix.sum()/photon_num)
    
    # plot the graph (R, T, A pie chart)
    plt.pie([R, T, A], labels=["R (theoretical value=22%)", "T (theoretical value=1.45%)", "A"], autopct='%1.3f%%')
    title = "Anisotropic scattering(g={}), photon_num={}, run_{}".format(g, photon_num, run_index+1)
    plt.title(title)
#     plt.savefig("hw_4\{}.png".format(title), bbox_inches = "tight")
    plt.show()
    
    # calculte and plot fluence rate
    module_singleLayer.calculate_fluenceRate_andPlotHist3D(grid_width=grid_width, grid_depth=grid_depth, 
                                         delta_r=delta_r, delta_z=delta_z, 
                                         absorption_matrix=absorption_matrix, photon_num=photon_num, 
                                         mua=mua, source_type=source, beam_power=beam_power
                                        )
    
    # store the R, T in the current time
    R_all[run_index] = R/photon_num
    T_all[run_index] = T/photon_num
    
print("Summary:")
print("R mean = {}, std = {}".format(R_all.mean(), R_all.std(ddof=1)))
print("T mean = {}, std = {}".format(T_all.mean(), T_all.std(ddof=1)))


#get_ipython().run_line_magic('matplotlib', 'notebook')
module_singleLayer.calculate_fluenceRate_andPlotHist3D(
    grid_width=grid_width, grid_depth=grid_depth, 
                                         delta_r=delta_r, delta_z=delta_z, 
                                         absorption_matrix=absorption_matrix, photon_num=photon_num, 
                                         mua=mua, source_type=source, beam_power=beam_power
                                        )



