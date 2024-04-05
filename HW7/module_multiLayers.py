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
        self.position_pseudoStart = np.empty(3)
        self.travel_distance = None
        self.current_layer = 0
        self.next_layer = None
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
        # record current position to the pseudoStart
        self.position_pseudoStart = self.position.copy()
        
        # photon start to travel
        self.travel_distance = self._get_pathlength(mua, mus)
        self.position[0] += self.cx*self.travel_distance
        self.position[1] += self.cy*self.travel_distance
        self.position[2] += self.cz*self.travel_distance
        
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
    
    def reflect_updatePosDir(self, tissue_model):
        if self.cz < 0:  # upward to downward
            incident_pathlength = (tissue_model["Depth"][self.current_layer][0] - self.position_pseudoStart[2]) / self.cz
            # update position
            self.position[2] = -self.position[2] + 2*tissue_model["Depth"][self.current_layer][0]
        else:            # downward to upward
            incident_pathlength = (tissue_model["Depth"][self.current_layer][1] - self.position_pseudoStart[2]) / self.cz
            # update position
            self.position[2] = -self.position[2] + 2*tissue_model["Depth"][self.current_layer][1]
        # update pseudoStart point
        self.position_pseudoStart[0] = self.position_pseudoStart[0] + self.cx*incident_pathlength
        self.position_pseudoStart[1] = self.position_pseudoStart[1] + self.cy*incident_pathlength
        self.position_pseudoStart[2] = self.position_pseudoStart[2] + self.cz*incident_pathlength
        # update direction
        self.cz = -self.cz
        
    def transmit_updatePosDir(self, tissue_model):
        # calculate S1 & new self.travel_distance (new self.travel_distance is equivalent to S2)
        if self.cz < 0:  # upward
            S1 = (tissue_model["Depth"][self.current_layer][0] - self.position_pseudoStart[2]) / self.cz
        else:            # downward
            S1 = (tissue_model["Depth"][self.current_layer][1] - self.position_pseudoStart[2]) / self.cz
        self.travel_distance = (tissue_model["Media"][self.current_layer]["mua"]+tissue_model["Media"][self.current_layer]["mus"])/(tissue_model["Media"][self.next_layer]["mua"]+tissue_model["Media"][self.next_layer]["mus"])*(self.travel_distance - S1)
        # update position
        self.position[0] = self.position_pseudoStart[0] + self.cx*(S1+self.travel_distance)
        self.position[1] = self.position_pseudoStart[1] + self.cy*(S1+self.travel_distance)
        self.position[2] = self.position_pseudoStart[2] + self.cz*(S1+self.travel_distance)
        # update current layer
        self.current_layer = self.next_layer
        # update pseudoStart point
        self.position_pseudoStart[0] = self.position_pseudoStart[0] + self.cx*S1
        self.position_pseudoStart[1] = self.position_pseudoStart[1] + self.cy*S1
        self.position_pseudoStart[2] = self.position_pseudoStart[2] + self.cz*S1
    
    def refract_updatePosDir(self, tissue_model):
        # calculate S1 & new self.travel_distance (new self.travel_distance is equivalent to S2)
        if self.cz < 0:  # upward
            S1 = (tissue_model["Depth"][self.current_layer][0] - self.position_pseudoStart[2]) / self.cz
        else:            # downward
            S1 = (tissue_model["Depth"][self.current_layer][1] - self.position_pseudoStart[2]) / self.cz
        self.travel_distance = (tissue_model["Media"][self.current_layer]["mua"]+tissue_model["Media"][self.current_layer]["mus"])/(tissue_model["Media"][self.next_layer]["mua"]+tissue_model["Media"][self.next_layer]["mus"])*(self.travel_distance - S1)
        # move photon into boundary
        self.position[0] = self.position_pseudoStart[0] + self.cx*S1
        self.position[1] = self.position_pseudoStart[1] + self.cy*S1
        self.position[2] = self.position_pseudoStart[2] + self.cz*S1
        # update pseudoStart point
        self.position_pseudoStart = self.position.copy()
        # update direction (due to refraction)
        incident_angle = np.arccos(abs(self.cz))
        n1 = tissue_model["Media"][new_photon.current_layer]["n"]
        n2 = tissue_model["Media"][new_photon.next_layer]["n"]
        refractive_angle = np.arcsin(n1/n2*np.sin(incident_angle))
        self.cx = n1/n2 * self.cx
        self.cy = n1/n2 * self.cy
        self.cz = np.cos(refractive_angle) * np.sign(self.cz)
        # move photon to final position
        self.position[0] = self.position[0] + self.cx*self.travel_distance
        self.position[1] = self.position[1] + self.cy*self.travel_distance
        self.position[2] = self.position[2] + self.cz*self.travel_distance
        # update current layer
        self.current_layer = self.next_layer        
    
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
    
def calculate_fluenceRate_andPlotHist3D(tissue_model, absorption_matrix, plotAbsorptionMatrix=False):
    grid_width = tissue_model["Session"]["grid_width"]
    grid_depth = tissue_model["Session"]["grid_depth"]
    delta_r = tissue_model["Session"]["delta_r"]
    delta_z = tissue_model["Session"]["delta_z"]
    photon_num = tissue_model["Session"]["photon_num"]
    source_type = tissue_model["Session"]["source"]
    beam_power = tissue_model["Session"]["beam_power"]
    
    # plot absorption_matrix (per volume -- 1/cm^3)
    if plotAbsorptionMatrix:
        # transform to per volume
        absorption_matrix_V = np.empty_like(absorption_matrix)
        for i_r in range(int(grid_width/delta_r)):
            V_i_r = (2*i_r+1) * np.pi * delta_z * (delta_r**2)  # cm^3
            for i_z in range(int(grid_depth/delta_z)):
                absorption_matrix_V[i_r][i_z] = absorption_matrix[i_r][i_z]/V_i_r
        absorption_matrix_log = np.ma.log10(absorption_matrix_V).filled(0)
        
        # plot bar3D
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111, projection='3d')
        R_grid, Z_grid = np.meshgrid(np.arange(grid_width/delta_r), np.arange(grid_depth/delta_z))

        x_pos = R_grid.T.ravel()*delta_r
        y_pos = Z_grid.T.ravel()*delta_z
        z_pos = 0

        dx = delta_r/2 * np.ones_like(z_pos)
        dy = delta_z/2 * np.ones_like(z_pos)
        dz = absorption_matrix_log.ravel()

        ax.bar3d(x_pos, y_pos, z_pos, dx, dy, dz, shade=True)
        ax.set_xlabel('Width [cm]', labelpad=10); 
        ax.set_ylabel('Depth [cm]', labelpad=10); 
        ax.set_zlabel('# of Photons Absorbed (in log scale)', labelpad=10);
        ax.set_title("Absorption_distribution of scattered photons, photon_num = {}, Source_type = {}".format(photon_num, source_type))
        plt.show()
        
        # plot 2D (with r & z)
        R = np.arange(0, grid_width, delta_r)
        Z = np.arange(0, grid_depth, delta_z)
        plt.plot(R, absorption_matrix_V.sum(axis=1), '-o')
        plt.xlabel("r (\u0394r = {}cm)".format(delta_r))
        plt.ylabel("# of Photons Absorbed [$1/cm^3$]")
        plt.title("Absorption_distribution of scattered photons (r dimension), photon_num = {}, Source_type = {}".format(photon_num, source_type), pad = 15)
        plt.show()
        plt.plot(Z, absorption_matrix_V.sum(axis=0), '-o')
        plt.xlabel("z (\u0394z = {}cm)".format(delta_z))
        plt.ylabel("# of Photons Absorbed [$1/cm^3$]")
        plt.title("Absorption_distribution of scattered photons (z dimension), photon_num = {}, Source_type = {}".format(photon_num, source_type), pad = 15)
        plt.show()
        
    # calculte fluence rate for the given beam power (default = 1W)
    absorption_matrix_VN = np.empty_like(absorption_matrix)
    for i_r in range(int(grid_width/delta_r)):
        V_i_r = (2*i_r+1) * np.pi * delta_z * (delta_r**2)  # cm^3
        for i_z in range(int(grid_depth/delta_z)):
            absorption_matrix_VN[i_r][i_z] = absorption_matrix[i_r][i_z]/(V_i_r*photon_num)
    fluence_rate = np.empty_like(absorption_matrix_VN)
    for tissue_layer in range(1, tissue_model["Session"]["num_layer"]-1):
        start = int(tissue_model["Depth"][tissue_layer][0]/delta_z)
        end = int(tissue_model["Depth"][tissue_layer][1]/delta_z)
        fluence_rate[:, start:end] = absorption_matrix_VN[:, start:end]/tissue_model["Media"][tissue_layer]["mua"]
    fluence_rate = fluence_rate * float(beam_power)
    
    # plot fluence rate 3D plot
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, projection='3d')
    R_grid, Z_grid = np.meshgrid(np.arange(grid_width/delta_r), np.arange(grid_depth/delta_z))
    
    x_pos = R_grid.T.ravel()*delta_r
    y_pos = Z_grid.T.ravel()*delta_z
    z_pos = 0
    
    dx = delta_r/2 * np.ones_like(z_pos)
    dy = delta_z/2 * np.ones_like(z_pos)
    dz = fluence_rate.ravel()
    
    ax.bar3d(x_pos, y_pos, z_pos, dx, dy, dz, shade=True)
    ax.set_xlabel('Width [cm]', labelpad=10)
    ax.set_ylabel('Depth [cm]', labelpad=10)
    ax.set_zlabel('Fluence_rate [1/cm^2]', labelpad=10)
    ax.set_title("Fluence_rate of scattered photons, Incident beam power = {}W, Source_type = {}".format(beam_power, source_type))
    plt.show()
    
    # plot 2D (with r & z)
    R = np.arange(0, grid_width, delta_r)
    Z = np.arange(0, grid_depth, delta_z)
    plt.plot(R, fluence_rate.sum(axis=1), '-o')
    plt.xlabel("r (\u0394r = {}cm)".format(delta_r))
    plt.ylabel("fluence rate [$W/cm^2$]")
    plt.title("Fluence_rate of scattered photons (r dimension), Source_type = {}".format(source_type), pad = 15)
#     plt.title("Fluence_rate of scattered photons, Incident beam power = {}W, Source_type = {}".format(beam_power, source_type))
    plt.show()
    plt.plot(Z, fluence_rate.sum(axis=0), '-o')
    plt.xlabel("z (\u0394z = {}cm)".format(delta_z))
    plt.ylabel("fluence rate [$W/cm^2$]")
    plt.title("Fluence_rate of scattered photons (z dimension), Source_type = {}".format(source_type), pad = 15)
#     plt.title("Fluence_rate of scattered photons, Incident beam power = {}W, Source_type = {}".format(beam_power, source_type))
    plt.show()


