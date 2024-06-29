"""------------------------------------------------------------------------------------------------
   -------------------------- Author: Lewis Russell : SUPAASP Class Student -----------------------
   ---------------------------------------------------------------------------------------------"""

import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams['figure.dpi'] = 300

def read_data(directory,type):
    data = {}
    for filename in os.listdir(directory):
        if filename.endswith(type,5,12):
            filepath = os.path.join(directory, filename)
            print(f"Reading data from {filepath}")
            with open(filepath, 'r') as file:
                lines = file.readlines()
                num_particles = int(lines[0].split(":")[1].strip())
                num_timesteps = int(lines[1].split(":")[1].strip())
                box_lengths = list(map(float, lines[2].split(":")[1].split(',')))
                timesteps = []
                for i in range(6, len(lines), num_particles + 2):
                    time = float(lines[i].split(":")[1].strip())
                    iteration_step = int(lines[i+1].split(":")[1].strip())
                    coordinates = np.array([list(map(float, line.split())) 
                                            for line in lines[i+2:i+2+num_particles]])
                    timesteps.append({'time': time, 'iteration_step': iteration_step, 
                                      'coordinates': coordinates})
                #print(timesteps['time'])
                data[filename] = {'num_particles': num_particles, 'num_timesteps': num_timesteps, 
                                  'box_lengths': box_lengths, 'timesteps': timesteps}
    return data

def plot_data_2D(data):
    for filename, filedata in data.items():
        time = np.array([timestep['time'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])

        fig, ((ax1)) = plt.subplots(1, 1, \
                            constrained_layout=True, \
                            sharex=True)
        for k in range(num_particles):
            ax1.plot(time, (coordinates[:,k,0]-coordinates[0,k,0]), label="X_{}(m)".format(k+1))
            ax1.plot(time, (coordinates[:,k,1]-coordinates[0,k,1]), label="Y_{}(m)".format(k+1))
            ax1.plot(time, (coordinates[:,k,2]-coordinates[0,k,2]), label="Z_{}(m)".format(k+1))
        
        ax1.set_ylabel("Coordinate (m)")
        ax1.set_title("N={}".format(num_particles))
        ax1.grid()
        #ax1.legend()
        plt.xlabel("t(s)")
        #plt.savefig(f'post_processing//trajectories_2D_{filename}.png')  # Save the plot to a file
        plt.savefig(f'trajectories_2D_{filename}.png')  # Save the plot to a file
        plt.close()

def plot_data_3D(data):

    for filename, filedata in data.items():
        fig = plt.figure(figsize=plt.figaspect(0.9),constrained_layout=True)
        tick_fontsize = 6
        titles_size = 8
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        #ax1
        ax1 = fig.add_subplot(2,2,1, projection='3d', proj_type='persp')
        ax1.set_xticklabels([])
        for i in range(num_particles):
            #ax1.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
            #        marker='.', linewidth=1)
            ax1.plot(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                    marker=',', linewidth=0.1)
        elev, azim = 0, 0
        ax1.view_init(elev=elev,azim=azim)
        ax1.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
            fontsize=titles_size, y=0.85)
        ax1.xaxis.set_tick_params(labelsize=tick_fontsize)
        ax1.yaxis.set_tick_params(labelsize=tick_fontsize)
        ax1.zaxis.set_tick_params(labelsize=tick_fontsize)
        #ax2
        ax2 = fig.add_subplot(2,2,2, projection='3d', proj_type='persp')
        ax2.set_yticklabels([])
        for i in range(num_particles):
            #ax2.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
            #          marker='.', linewidth=1)
            ax2.plot(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker=',', linewidth=0.1)
        elev, azim = 0, -90
        ax2.view_init(elev=elev,azim=azim)
        ax2.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
        fontsize=titles_size, y=0.85)
        ax2.xaxis.set_tick_params(labelsize=tick_fontsize)
        ax2.yaxis.set_tick_params(labelsize=tick_fontsize)
        ax2.zaxis.set_tick_params(labelsize=tick_fontsize)
        #ax3
        ax3 = fig.add_subplot(2,2,3, projection='3d', proj_type='persp')
        for i in range(num_particles):
            #ax3.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
            #          marker='.', linewidth=1)
            ax3.plot(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker=',', linewidth=0.1)
        elev, azim = 30, 60
        ax3.view_init(elev=elev,azim=azim)
        ax3.view_init(elev=elev,azim=azim)
        ax3.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
            fontsize=titles_size, y=1)
        ax3.xaxis.set_tick_params(labelsize=tick_fontsize)
        ax3.yaxis.set_tick_params(labelsize=tick_fontsize)
        ax3.zaxis.set_tick_params(labelsize=tick_fontsize)
        #ax4
        ax4 = fig.add_subplot(2,2,4, projection='3d', proj_type='persp')
        for i in range(num_particles):
            #ax4.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
            #          marker='.', linewidth=1)
            ax4.plot(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker=',', linewidth=0.2)
        elev, azim = 30, -60
        ax4.view_init(elev=elev,azim=azim)
        ax4.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
            fontsize=titles_size, y=1)
        ax4.xaxis.set_tick_params(labelsize=tick_fontsize)
        ax4.yaxis.set_tick_params(labelsize=tick_fontsize)
        ax4.zaxis.set_tick_params(labelsize=tick_fontsize)
        #plot_end
        fig.suptitle("Particle Trajectories, Perspective FoV Visualisation", \
                 fontsize=12, y=0.93)
        #plt.savefig(f'post_processing//trajectories_3D_{filename}.png')  # Save the plot to a file
        plt.savefig(f'trajectories_3D_{filename}.png')  # Save the plot to a file
        plt.close()

def plot_MSD(data):

    for filename, filedata in data.items():
        time = np.array([timestep['iteration_step'] for timestep in filedata['timesteps']])
        MSDtime = np.array([timestep['time'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
    
        MSD_theory = 6*(1/25)*(MSDtime)
        msd = []
        for i in range(len(time)):
            hold = []
            for k in range(num_particles):
                #msd_pre = np.linalg.norm((r[k][3*i:3*(i+1)]) - (r[k][0:3]))**2
                msd_pre = np.linalg.norm((coordinates[i,k,0:3] - coordinates[0,k,0:3]))**2
                #coordinates[:,i,0]
                hold.append(msd_pre)
            msd.append(np.sum(hold)/num_particles)
        meanT = np.mean(MSD_theory)
        meanA = np.mean(msd)
        print('\n #############--------------------#############', \
        "\n For N={} particles, time averaged: \n".format(num_particles), \
        '\n Theoretical MSD =', meanT, \
        '\n Actual MSD =', meanA, \
        '\n #############--------------------#############')
        #plt.scatter(t,msd,marker=".", s = 1)
        #plt.scatter(t,MSD_theory,marker=".", s = 1)
        plt.plot(time/1000,msd)
        plt.plot(time/1000,MSD_theory)
        plt.xlabel(r"Iteration Step ($x10^{-3}$)")
        plt.ylabel("MSD")
        plt.legend(
            ["Actual MSD (Mean = {:.2f})".format(meanA),
             "Theoretical MSD (Mean = {:.2f})".format(meanT)], 
            loc='upper right')
        plt.grid()
        plt.title("MSD for N={} particles".format(num_particles))
        #plt.savefig(f'post_processing//MSD_{filename}.png')  # Save the plot to a file
        plt.savefig(f'MSD_{filename}.png')  # Save the plot to a file
        plt.close()
        
def plot_MSD_multiple(data):
    fig, ax = plt.subplots()  # Initialize the figure and axes outside of the loop
    meanA_list = []  # Create an empty list to store meanA values

    for filename, filedata in data.items():
        time = np.array([timestep['iteration_step'] for timestep in filedata['timesteps']])
        MSDtime = np.array([timestep['time'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
    
        MSD_theory = 6*(1/25)*(MSDtime)
        msd = []
        for i in range(len(time)):
            hold = []
            for k in range(num_particles):
                msd_pre = np.linalg.norm((coordinates[i,k,0:3] - coordinates[0,k,0:3]))**2
                hold.append(msd_pre)
            msd.append(np.sum(hold)/num_particles)
        meanT = np.mean(MSD_theory)
        meanA = np.mean(msd)
        meanA_list.append(meanA)  # Append meanA value to the list
        print('\n #############--------------------#############', \
        "\n For N={} particles, time averaged: \n".format(num_particles), \
        '\n Theoretical MSD =', meanT, \
        '\n Actual MSD =', meanA, \
        '\n #############--------------------#############')
        ax.plot(time/1000,msd, label="Observed MSD (Mean = {:.2f})".format(meanA))  # Plot on the same axes
    
    ax.plot(time/1000,MSD_theory, label="Theoretical MSD (Mean = {:.2f})".format(meanT))  # Plot on the same axes

    ax.set_xlabel(r"Iteration Step ($x10^{-3}$)")
    ax.set_ylabel("MSD")
    ax.legend(
        loc='upper left')  # Include meanA values in the legend
    ax.grid()
    ax.set_title("MSD for N={} particles".format(num_particles))
    plt.savefig(f'MSD_multiple_{len(data)}_runs.png')  # Save the plot to a file
    plt.close()
    
def plot_MSD_multipleLJ(data):
    fig, ax = plt.subplots()  # Initialize the figure and axes outside of the loop
    meanA_list = []  # Create an empty list to store meanA values
    LJ=[0,4,8,10]
    LJi=0
    for filename, filedata in data.items():
        time = np.array([timestep['iteration_step'] for timestep in filedata['timesteps']])
        MSDtime = np.array([timestep['time'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
    
        MSD_theory = 6*(1/25)*(MSDtime)
        msd = []
        for i in range(len(time)):
            hold = []
            for k in range(num_particles):
                msd_pre = np.linalg.norm((coordinates[i,k,0:3] - coordinates[0,k,0:3]))**2
                hold.append(msd_pre)
            msd.append(np.sum(hold)/num_particles)
        meanT = np.mean(MSD_theory)
        meanA = np.mean(msd)
        meanA_list.append(meanA)  # Append meanA value to the list
        print('\n #############--------------------#############', \
        "\n For N={} particles, time averaged: \n".format(num_particles), \
        '\n Theoretical MSD =', meanT, \
        '\n Actual MSD =', meanA, \
        '\n #############--------------------#############')
        ax.plot(time/1000,msd, 
                label="$\epsilon,\sigma={}$ (Observed MSD = {:.2f})".format(LJ[LJi],meanA))  # Plot on the same axes
        LJi+=1
    
    ax.plot(time/1000,MSD_theory, 
            label="Theoretical MSD (no $E_p$) = {:.2f}".format(meanT))  # Plot on the same axes

    ax.set_xlabel(r"Iteration Step ($x10^{-3}$)")
    ax.set_ylabel("MSD")
    ax.legend(
        loc='upper left')  # Include meanA values in the legend
    ax.grid()
    ax.set_title("MSD for N={} particles, Various LJ Scales".format(num_particles))
    plt.savefig(f'MSD_multiple_LJ_{len(data)}_runs.png')  # Save the plot to a file
    plt.close()

def plot_kinetic_energy(data):
    for filename, filedata in data.items():
        time = np.array([timestep['iteration_step'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        velocities = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
        
        kinetic_energy_theory = 1.5 * num_particles * np.ones(len(time))
        kinetic_energy = []
        for i in range(len(time)):
            hold = []
            for k in range(num_particles):
                kinetic_energy_pre = 0
                for d in range(3):
                    kinetic_energy_pre += 0.5 * (velocities[i,k,d])**2
                hold.append(kinetic_energy_pre)
            kinetic_energy.append(np.sum(hold))
        
        mean_theory = np.mean(kinetic_energy_theory)
        mean_actual = np.mean(kinetic_energy)
        
        print('\n #############--------------------#############', \
        "\n For N={} particles, time averaged: \n".format(num_particles), \
        '\n Theoretical Mean Kinetic Energy =', mean_theory, \
        '\n Actual Mean Kinetic Energy =', mean_actual, \
        '\n #############--------------------#############')
        
        plt.plot(time/1000, kinetic_energy)
        plt.plot(time/1000, kinetic_energy_theory)
        plt.xlabel(r"Iteration Step ($x10^{-3}$)")
        plt.ylabel(r"E$_{k}$")
        plt.legend(
            ["Observed Kinetic Energy (Mean = {:.2f})".format(mean_actual),
             "Theoretical Kinetic Energy (Mean = {:.2f})".format(mean_theory)], 
            loc='upper right')
        plt.grid()
        plt.title("N={}".format(num_particles))
        #plt.savefig(f'post_processing//kinetic_energy_{filename}.png')  # Save the plot to a file
        plt.savefig(f'kinetic_energy_{filename}.png')  # Save the plot to a file
        plt.close()
        
def plot_kinetic_energy_running_avg(data):
    for filename, filedata in data.items():
        time = np.array([timestep['iteration_step'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        velocities = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
        
        kinetic_energy_theory = 1.5 * num_particles * np.ones(len(time))
        kinetic_energy = []
        running_avg = []
        cumulative_sum = 0
        
        for i in range(len(time)):
            hold = []
            for k in range(num_particles):
                kinetic_energy_pre = 0
                for d in range(3):
                    kinetic_energy_pre += 0.5 * (velocities[i,k,d])**2
                hold.append(kinetic_energy_pre)
            kinetic_energy.append(np.sum(hold))
            cumulative_sum += kinetic_energy[i]
            running_avg.append(cumulative_sum / (i+1))
        
        mean_theory = np.mean(kinetic_energy_theory)
        mean_actual = np.mean(kinetic_energy)
        
        print('\n #############--------------------#############', \
        "\n For N={} particles, time averaged: \n".format(num_particles), \
        '\n Theoretical Mean Kinetic Energy =', mean_theory, \
        '\n Actual Mean Kinetic Energy =', mean_actual, \
        '\n #############--------------------#############')
        
        plt.scatter(time/1000, running_avg, marker='x', s=2)
        plt.plot(time/1000, kinetic_energy_theory, color='red')
        plt.xlabel(r"Iteration Step ($x10^{-3}$)")
        plt.ylabel(r"E$_{k}$")
        #plt.ylim(0.9 * np.min(running_avg), 1.1 * np.max(running_avg))
        plt.legend(
            ["Running Average Kinetic Energy",
             "Theoretical Kinetic Energy (Mean = {:.2f})".format(mean_theory)], 
            loc='upper right')
        plt.grid()
        plt.title("N={}".format(num_particles))
        #plt.savefig(f'post_processing//kinetic_energy_{filename}.png')  # Save the plot to a file
        plt.savefig(f'kinetic_energy_run_avg_{filename}.png')  # Save the plot to a file
        plt.close()
        
def Kokkos_Perf_Plot():
    Natoms = np.array([10,50,100,150,200])
    Vanilla = np.array([1.76,38.88,153.16,342.55,607.42])
    KK1 = np.array([5.80,67.23,236.64,508.51,885.24])
    KK2 = np.array([4.00,47.69,171.67,374.17,654.34])
    KK4 = np.array([3.02,29.25,99.82,220.26,383.66])
    KK8 = np.array([2.57,17.21,56.97,118.88,206.45])
    # 1st plot
    fig, ax1 = plt.subplots(1, 1, \
                            constrained_layout=True, \
                            sharex=True)
    ax1.plot(Natoms, Vanilla/Vanilla, label="Vanilla C++", linestyle='--')
    ax1.plot(Natoms, Vanilla/KK1, label="KOKKOS, 1 OpenMP Thread")
    ax1.plot(Natoms, Vanilla/KK2, label="KOKKOS, 2 OpenMP Threads")
    ax1.plot(Natoms, Vanilla/KK4, label="KOKKOS, 4 OpenMP Threads")
    ax1.plot(Natoms, Vanilla/KK8, label="KOKKOS, 8 OpenMP Threads")
    #ax1.scatter(Natoms, Vanilla/Vanilla, marker='x', color='blue')
    ax1.scatter(Natoms, Vanilla/KK1, marker='x', color='orange')
    ax1.scatter(Natoms, Vanilla/KK2, marker='x', color='green')
    ax1.scatter(Natoms, Vanilla/KK4, marker='x', color='red')
    ax1.scatter(Natoms, Vanilla/KK8, marker='x', color='purple')
    ax1.legend(loc='upper right',bbox_to_anchor=(0.8,1))
    ax1.set_ylabel("Runtime Ratio (Vanilla/Runtime'$x$')")
    ax1.set_xlabel("Number of Atoms")
    ax1.set_xticks(Natoms)
    ax1.grid()
    ax1.legend()
    plt.ylim(0,4)
    plt.title("Ratio Comparison of Vanilla C++ and Kokkos Runtimes")
    plt.savefig('Performance_Ratio.png')  # Save the plot to a file
    plt.close()
    # 2nd plot
    fig, ax1 = plt.subplots(1, 1, \
                            constrained_layout=True, \
                            sharex=True)
    ax1.set_ylabel("Total Runtime (seconds)")
    ax1.set_xlabel("Number of Atoms")
    ax1.set_xticks(Natoms)
    ax1.grid()
    ax1.plot(Natoms, Vanilla, label="Vanilla C++", linestyle=':')
    ax1.plot(Natoms, KK1, label="KOKKOS, 1 OpenMP Thread")
    ax1.plot(Natoms, KK2, label="KOKKOS, 2 OpenMP Threads")
    ax1.plot(Natoms, KK4, label="KOKKOS, 4 OpenMP Threads")
    ax1.plot(Natoms, KK8, label="KOKKOS, 8 OpenMP Threads")
    ax1.scatter(Natoms, Vanilla, marker='x', color='blue')
    ax1.scatter(Natoms, KK1, marker='x', color='orange')
    ax1.scatter(Natoms, KK2, marker='x', color='green')
    ax1.scatter(Natoms, KK4, marker='x', color='red')
    ax1.scatter(Natoms, KK8, marker='x', color='purple')
    ax1.set_ylabel("Total Runtime (seconds)")
    #plt.ylim(0,5)
    ax1.legend(loc='upper left')
    plt.title("Time Comparison of Vanilla C++ and Kokkos Runtimes")
    plt.savefig('Performance_Raw.png')  # Save the plot to a file
    plt.close()

# Usage:
# VS sets path from git -> "../src_written_data" not required
# If not in VSCode, just comment on/off *_data as needed
#
#pos_data = read_data("src_written_data","pos.txt")
#vel_data = read_data("src_written_data","vel.txt")
pos_data = read_data("../src_written_data","pos.txt")
vel_data = read_data("../src_written_data","vel.txt")
#plot_data_2D(pos_data)
#plot_data_3D(pos_data)
plot_MSD(pos_data)
#plot_MSD_multiple(pos_data)
#plot_MSD_multipleLJ(pos_data)
plot_kinetic_energy(vel_data)
#plot_kinetic_energy_running_avg(vel_data)
#plot_data_2D(vel_data)
Kokkos_Perf_Plot()