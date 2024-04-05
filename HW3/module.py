import numpy as np
import matplotlib.pyplot as plt

class photon:
    def __init__(self, variable_weight=False):
        self.position = np.array([0, 0, 0], dtype=float)
        self.cx = 0
        self.cy = 0
        self.cz = 1
        if variable_weight is True:
            self.weight = 1

    def _get_pathlength(self, mua, mus):
        pathlength = -(np.log(np.random.uniform()) / (mua+mus))
        return pathlength
    
    def scattering(self, g=0):
        phi = 2*np.pi*np.random.uniform()
        if g is 0:
            # isotropic scattering
            theta = np.arccos(2*np.random.uniform() - 1)
        else:
            # Henyey-Greenstein phase function
            theta = np.arccos(1/(2*g)*(1+g**2-((1-g**2)/(1-g+2*g*np.random.uniform()))**2))
        
        if abs(self.cz) > 0.99999:
            self.cx = np.sin(theta)*np.cos(phi)
            self.cy = np.sin(theta)*np.sin(phi)
            self.cz = np.cos(theta)*self.cz/abs(self.cz)
        else:
            self.cx = (np.sin(theta)/np.sqrt(1-self.cz**2)) * (self.cx*self.cz*np.cos(phi) - self.cy*np.sin(phi)) + self.cx*np.cos(theta)
            self.cy = (np.sin(theta)/np.sqrt(1-self.cz**2)) * (self.cy*self.cz*np.cos(phi) + self.cx*np.sin(phi)) + self.cy*np.cos(theta)
            self.cz = -np.sin(theta)*np.cos(phi)*np.sqrt(1-self.cz**2) + self.cz*np.cos(theta)
            
    def travel(self, mua, mus):
        travel_distance = self._get_pathlength(mua, mus)
        self.position[0] += self.cx*travel_distance
        self.position[1] += self.cy*travel_distance
        self.position[2] += self.cz*travel_distance
        
    def play_roulette(self, photon_survived_prob):
        if np.random.uniform() < photon_survived_prob:
            self.weight = 1/photon_survived_prob * self.weight
        else:
            self.weight = 0

# outer function
def plot_photon_path(photon_tracks, photon_index=0):
    fig, ax = plt.subplots()
#     for photon_index in photon_tracks.keys():        
    ax.plot(photon_tracks[photon_index][0], photon_tracks[photon_index][1], marker="o", color="black")
    ax.quiver(photon_tracks[photon_index][0][:-1], photon_tracks[photon_index][1][:-1], 
              photon_tracks[photon_index][0][1:]-photon_tracks[photon_index][0][:-1], 
              photon_tracks[photon_index][1][1:]-photon_tracks[photon_index][1][:-1], 
              scale_units='xy', angles='xy', scale=1.2, headwidth=5, headlength=7, width=0.005)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_facecolor("#FFDEC9")
    ax.axhline(y=0, label="surface", color="#B81313")
    ax.axhline(y=0.02, label="boundary", linestyle="--", color="#B81313")        
    ax.legend()
    plt.show()
#         print("final photon position:", new_photon.position[1:3])
#         print("scattering times:", len(photon_y)-2)
        


