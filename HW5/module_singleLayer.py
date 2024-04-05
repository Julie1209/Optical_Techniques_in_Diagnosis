import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class photon:
    def __init__(self, variable_weight=False, source_type="Pencil", source_radius=0):
        if source_type == "Pencil":
            self.position = np.array([0, 0, 0], dtype=float)
        if source_type == "Uniform":
            r = source_radius * np.sqrt(np.random.uniform())
            self.position = np.array([r, 0, 0], dtype=float)
        if source_type == "Gaussian":
            r = source_radius * np.sqrt(-np.log(np.random.uniform())/2)
            self.position = np.array([r, 0, 0], dtype=float)
        self.cx = 0
        self.cy = 0
        self.cz = 1
        phi = None
        theta = None
        if variable_weight is True:
            self.weight = 1

    def _get_pathlength(self, mua, mus):
        pathlength = -(np.log(np.random.uniform()) / (mua+mus))
        return pathlength
    
    def scattering(self, g=0):
        # phi: 0~360 degrees
        # theta: 0~180 degrees
        
        # get phi and theta
        phi = 2*np.pi*np.random.uniform()
        if g is 0:
            # isotropic scattering
            theta = np.arccos(2*np.random.uniform() - 1)
        else:
            # Henyey-Greenstein phase function
            theta = np.arccos(1/(2*g)*(1+g**2-((1-g**2)/(1-g+2*g*np.random.uniform()))**2))
        
        # update direction with rotation matrix (use phi & theta calculated above)
        if abs(self.cz) > 0.99999:
            self.cxprime = np.sin(theta)*np.cos(phi)
            self.cyprime = np.sin(theta)*np.sin(phi)
            self.czprime = np.cos(theta)*self.cz/abs(self.cz)
        else:
            self.cxprime = (np.sin(theta)/np.sqrt(1-self.cz**2)) * (self.cx*self.cz*np.cos(phi) - self.cy*np.sin(phi)) + self.cx*np.cos(theta)
            self.cyprime = (np.sin(theta)/np.sqrt(1-self.cz**2)) * (self.cy*self.cz*np.cos(phi) + self.cx*np.sin(phi)) + self.cy*np.cos(theta)
            self.czprime = -np.sin(theta)*np.cos(phi)*np.sqrt(1-self.cz**2) + self.cz*np.cos(theta)
        self.cx = self.cxprime
        self.cy = self.cyprime
        self.cz = self.czprime
            
    def travel(self, mua, mus):
        travel_distance = self._get_pathlength(mua, mus)
        self.position[0] += self.cx*travel_distance
        self.position[1] += self.cy*travel_distance
        self.position[2] += self.cz*travel_distance
        
    def _calculate_Rs_Rp(self, incident_angle, n1, n2):
        refractive_angle = np.arcsin(n1/n2*np.sin(incident_angle))
        Rs = ((n1*np.cos(incident_angle)-n2*np.cos(refractive_angle))/(n1*np.cos(incident_angle)+n2*np.cos(refractive_angle)))**2
        Rp = ((n1*np.cos(refractive_angle)-n2*np.cos(incident_angle))/(n1*np.cos(refractive_angle)+n2*np.cos(incident_angle)))**2
        return Rs, Rp                
    
    def get_boundary_split_ratio(self, n1, n2):
        incident_angle = np.arccos(abs(self.cz))
        # judge if total internal reflection may occur
        if n1 < n2:  # total internal reflection never occur, do not calculate criticle angle
            Rs, Rp = self._calculate_Rs_Rp(incident_angle, n1, n2)
            r = (Rs+Rp)/2
            t = 1-r
        else:        # total internal reflection may occur, calculate criticle angle
            critical_angle = np.arcsin(n2/n1)        
            if incident_angle >= critical_angle:
                r = 1
                t = 0
            else:
                Rs, Rp = self._calculate_Rs_Rp(incident_angle, n1, n2)
                r = (Rs+Rp)/2
                t = 1-r
        return t, r
    
    def reflect_updatePosDir(self):
        self.position[2] = -self.position[2]
        self.cz = -self.cz
        
    def transmit_updatePosDir(self, medium_thickness):
        d = self.position[2] - medium_thickness
        self.position[2] = self.position[2] - 2*d
        self.cz = -self.cz
    
    def play_roulette(self, photon_survived_prob):
        if np.random.uniform() < photon_survived_prob:
            self.weight = 1/photon_survived_prob * self.weight
        else:
            self.weight = 0

# outer function
def plot_photon_path(photon_tracks, photon_index=0, arrowwidth=0.005):
    fig, ax = plt.subplots()
#     for photon_index in photon_tracks.keys():        
    ax.plot(photon_tracks[photon_index][0], photon_tracks[photon_index][1], marker="o", color="black")
    ax.quiver(photon_tracks[photon_index][0][:-1], photon_tracks[photon_index][1][:-1], 
              photon_tracks[photon_index][0][1:]-photon_tracks[photon_index][0][:-1], 
              photon_tracks[photon_index][1][1:]-photon_tracks[photon_index][1][:-1], 
              scale_units='xy', angles='xy', scale=1.2, headwidth=5, headlength=7, width=arrowwidth)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_facecolor("#FFDEC9")
    ax.axhline(y=0, label="surface", color="#B81313")
    ax.axhline(y=0.02, label="boundary", linestyle="--", color="#B81313")        
    ax.legend()
    plt.show()
#         print("final photon position:", new_photon.position[1:3])
#         print("scattering times:", len(photon_y)-2)

def plot_photon_path(photon_tracks, photon_index=0, arrowwidth=0.005):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
#     for photon_index in photon_tracks.keys():        
    ax.plot(photon_tracks[photon_index][0], photon_tracks[photon_index][1], marker="o", color="black")
    ax.quiver(photon_tracks[photon_index][0][:-1], photon_tracks[photon_index][1][:-1], 
              photon_tracks[photon_index][0][1:]-photon_tracks[photon_index][0][:-1], 
              photon_tracks[photon_index][1][1:]-photon_tracks[photon_index][1][:-1], 
              scale_units='xy', angles='xy', scale=1.2, headwidth=5, headlength=7, width=arrowwidth)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_facecolor("#FFDEC9")
    ax.axhline(y=0, label="surface", color="#B81313")
    ax.axhline(y=0.02, label="boundary", linestyle="--", color="#B81313")        
    ax.legend()
    plt.show()
    
def calculate_fluenceRate_andPlot(grid_width, grid_depth, delta_r, delta_z, absorption_matrix, photon_num, mua, source_type):
    # calculte fluence rate
    absorption_matrix_V = np.empty_like(absorption_matrix)
    for i_r in range(int(grid_width/delta_r)):
        V_i_r = (2*i_r+1) * np.pi * delta_z * (delta_r**2)  # cm^3
        for i_z in range(int(grid_depth/delta_z)):
            absorption_matrix_V[i_r][i_z] = absorption_matrix[i_r][i_z]/(V_i_r*photon_num)
    fluence_rate = absorption_matrix_V/mua
    
    # plot fluence rate 3D plot
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection='3d')
    R_grid, Z_grid = np.meshgrid(np.arange(grid_width/delta_r), np.arange(grid_depth/delta_z))
    ax.plot_surface(R_grid.T, Z_grid.T, fluence_rate)
    ax.set_xlabel('R', labelpad=10); 
    ax.set_ylabel('Z', labelpad=10); 
    ax.set_zlabel('Fluence_rate', labelpad=10);
    ax.set_title("Fluence_rate_distribution, Source_type = {}".format(source_type))
    plt.show()
    
def calculate_fluenceRate_andPlotHist3D(grid_width, grid_depth, delta_r, delta_z, absorption_matrix, photon_num, mua, source_type, beam_power="1", plotAbsorptionMatrix=False):
    # plot absorption_matrix
    if plotAbsorptionMatrix:
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, projection='3d')
        R_grid, Z_grid = np.meshgrid(np.arange(grid_width/delta_r), np.arange(grid_depth/delta_z))

        x_pos = R_grid.T.ravel()*delta_r
        y_pos = Z_grid.T.ravel()*delta_z
        z_pos = 0

        dx = dy = delta_r/2 * np.ones_like(z_pos)
        dz = absorption_matrix.ravel()

        ax.bar3d(x_pos, y_pos, z_pos, dx, dy, dz, shade=True)
        ax.set_xlabel('Width [cm]', labelpad=10); 
        ax.set_ylabel('Depth [cm]', labelpad=10); 
        ax.set_zlabel('# of Photons Absorbed [-]', labelpad=10);
        ax.set_title("Absorption_distribution of scattered photons, photon_num = {}, Source_type = {}".format(photon_num, source_type))
        plt.show()
    
    # calculte fluence rate for the given beam power (default = 1W)
    absorption_matrix_V = np.empty_like(absorption_matrix)
    for i_r in range(int(grid_width/delta_r)):
        V_i_r = (2*i_r+1) * np.pi * delta_z * (delta_r**2)  # cm^3
        for i_z in range(int(grid_depth/delta_z)):
            absorption_matrix_V[i_r][i_z] = absorption_matrix[i_r][i_z]/(V_i_r*photon_num)
    fluence_rate = absorption_matrix_V/mua * float(beam_power)
    
    # plot fluence rate 3D plot
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection='3d')
    R_grid, Z_grid = np.meshgrid(np.arange(grid_width/delta_r), np.arange(grid_depth/delta_z))
    
    x_pos = R_grid.T.ravel()*delta_r
    y_pos = Z_grid.T.ravel()*delta_z
    z_pos = 0
    
    dx = dy = delta_r/2 * np.ones_like(z_pos)
    dz = fluence_rate.ravel()
    
    ax.bar3d(x_pos, y_pos, z_pos, dx, dy, dz, shade=True)
    ax.set_xlabel('Width [cm]', labelpad=10); 
    ax.set_ylabel('Depth [cm]', labelpad=10); 
    ax.set_zlabel('Fluence_rate [J/cm^2]', labelpad=10);
    ax.set_title("Fluence_rate_distribution, Incident beam power = {}W, Source_type = {}".format(beam_power, source_type))
    plt.show()
    
def absorptionMatrix_PlotHist3D(grid_width, grid_depth, delta_r, delta_z, absorption_matrix, photon_num, source_type):
    # plot absorption matrix 3D plot
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection='3d')
    R_grid, Z_grid = np.meshgrid(np.arange(grid_width/delta_r), np.arange(grid_depth/delta_z))
    
    x_pos = R_grid.T.ravel()*delta_r
    y_pos = Z_grid.T.ravel()*delta_z
    z_pos = 0
    
    dx = dy = delta_r/2 * np.ones_like(z_pos)
    dz = absorption_matrix.ravel()
    
    ax.bar3d(x_pos, y_pos, z_pos, dx, dy, dz, shade=True)
    ax.set_xlabel('Width [cm]', labelpad=10); 
    ax.set_ylabel('Depth [cm]', labelpad=10); 
    ax.set_zlabel('# of Photons Absorbed [-]', labelpad=10);
    ax.set_title("Absorption_distribution of scattered photons, Source_type = {}".format(source_type))
    plt.show()
    
        


