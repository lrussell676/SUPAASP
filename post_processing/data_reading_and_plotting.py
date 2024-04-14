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

'''def plot_data(data):
    for filename, filedata in data.items():
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        #for timestep in filedata['timesteps']:
            #time = timestep['time']
            #coordinates = timestep['coordinates']
            ##ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], label=f'Time: {time}')
        print(coordinates.shape)
        print(num_particles)
        for i in range(num_particles):
            ax.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker=',', linewidth=1)
            ax.plot(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker=',', linewidth=0.1)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'Particle Trajectories')
        ##ax.legend()
        plt.savefig(f'post_processing//trajectories_{filename}.png')  # Save the plot to a file

# Usage:
# VS sets path from git -> "../src_written_data" not required
pos_data = read_data("src_written_data","pos.txt")
plot_data(pos_data)'''

"""def plot_data(data):
    for filename, filedata in data.items():
        for timestep in filedata['timesteps']:
            time = timestep['time']
            coordinates = timestep['coordinates']
            fig = plt.figure(figsize=(10, 8))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], label=f'Time: {time}')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            ax.set_title(f'Coordinates at Time: {time} for {filename}')
            ax.legend()
            plt.show()"""

def plot_data(data):

    for filename, filedata in data.items():
        fig = plt.figure(figsize=plt.figaspect(0.9),constrained_layout=True)
        tick_fontsize = 6
        titles_size = 8
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        print(coordinates)
        print(num_particles)
        print(coordinates.shape)
        #ax1
        ax1 = fig.add_subplot(2,2,1, projection='3d', proj_type='persp')
        ax1.set_xticklabels([])
        for i in range(num_particles):
            ax1.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                    marker='.', linewidth=1)
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
            ax2.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker='.', linewidth=1)
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
            ax3.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker='.', linewidth=1)
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
            ax4.scatter(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker='.', linewidth=1)
            ax4.plot(coordinates[:,i,0], coordinates[:,i,1], coordinates[:,i,2],
                      marker=',', linewidth=0.1)
        elev, azim = 30, -60
        ax4.view_init(elev=elev,azim=azim)
        ax4.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
            fontsize=titles_size, y=1)
        ax4.xaxis.set_tick_params(labelsize=tick_fontsize)
        ax4.yaxis.set_tick_params(labelsize=tick_fontsize)
        ax4.zaxis.set_tick_params(labelsize=tick_fontsize)
        #plot_end
        fig.suptitle("Particle Trajectories", \
                 fontsize=12, y=0.93)
        plt.savefig(f'post_processing//trajectories_{filename}.png')  # Save the plot to a file

    """fig = plt.figure(figsize=plt.figaspect(0.9),constrained_layout=True)
    tick_fontsize = 6
    titles_size = 8
    #ax1
    ax1 = fig.add_subplot(2,2,1, projection='3d', proj_type='persp')
    ax1.set_xticklabels([])
    elev, azim = 0, 0
    ax1.view_init(elev=elev,azim=azim)
    ax1.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
        fontsize=titles_size, y=0.85)
    ax1.xaxis.set_tick_params(labelsize=tick_fontsize)
    ax1.yaxis.set_tick_params(labelsize=tick_fontsize)
    ax1.zaxis.set_tick_params(labelsize=tick_fontsize)
    for k in range(len(data[:,0])):
            ax1.plot(x[k][:steps], y[k][:steps], z[k][:steps], \
                      marker=',', linewidth=0.01)
    #ax2
    ax2 = fig.add_subplot(2,2,2, projection='3d', proj_type='persp')
    ax2.set_yticklabels([])
    elev, azim = 0, -90
    ax2.view_init(elev=elev,azim=azim)
    ax2.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
        fontsize=titles_size, y=0.85)
    ax2.xaxis.set_tick_params(labelsize=tick_fontsize)
    ax2.yaxis.set_tick_params(labelsize=tick_fontsize)
    ax2.zaxis.set_tick_params(labelsize=tick_fontsize)
    for k in range(len(data[:,0])):
            ax2.plot(x[k][:steps], y[k][:steps], z[k][:steps], \
                      marker=',', linewidth=0.01)
    #ax3
    ax3 = fig.add_subplot(2,2,3, projection='3d', proj_type='persp')
    ax3.view_init(elev=30,azim=60)
    elev, azim = 30, 60
    ax3.view_init(elev=elev,azim=azim)
    ax3.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
        fontsize=titles_size, y=1)
    ax3.xaxis.set_tick_params(labelsize=tick_fontsize)
    ax3.yaxis.set_tick_params(labelsize=tick_fontsize)
    ax3.zaxis.set_tick_params(labelsize=tick_fontsize)
    for k in range(len(data[:,0])):
            ax3.plot(x[k][:steps], y[k][:steps], z[k][:steps], \
                      marker=',', linewidth=0.01)
    #ax4
    ax4 = fig.add_subplot(2,2,4, projection='3d', proj_type='persp')
    elev, azim = 30, -60
    ax4.view_init(elev=elev,azim=azim)
    ax4.set_title('elev={:.0f}, azim={:.0f}'.format(elev,azim), \
        fontsize=titles_size, y=1)
    ax4.xaxis.set_tick_params(labelsize=tick_fontsize)
    ax4.yaxis.set_tick_params(labelsize=tick_fontsize)
    ax4.zaxis.set_tick_params(labelsize=tick_fontsize)
    for k in range(len(data[:,0])):
            ax4.plot(x[k][:steps], y[k][:steps], z[k][:steps], \
                      marker=',', linewidth=0.01)
    fig.suptitle("Particle Trajectories", \
                 fontsize=12, y=0.93)
    plt.show()"""

def MSD_plot(data):

    for filename, filedata in data.items():
        time = np.array([timestep['time'] for timestep in filedata['timesteps']])
        num_particles = filedata['num_particles']
        coordinates = np.array([timestep['coordinates'] for timestep in filedata['timesteps']])
    
    print(time.shape)
    MSD_theory = 6*(1/25)*(time)
    msd = []
    for i in range(len(time)):
        hold = []
        for k in range(num_particles):
            #msd_pre = np.linalg.norm((r[k][3*i:3*(i+1)]) - (r[k][0:3]))**2
            msd_pre = np.linalg.norm((coordinates[i,k,0:3] - coordinates[0,k,0:3]))**2
            #coordinates[:,i,0]
            hold.append(msd_pre)
        msd.append(np.sum(hold)/num_particles)

    print('\n #############--------------------#############', \
      "\n For N={} particles, time averaged: \n".format(num_particles), \
      '\n Theoretical MSD =', np.mean(MSD_theory), \
      '\n Actual MSD =', np.mean(msd), \
      '\n #############--------------------#############')
#   plt.scatter(t,msd,marker=".", s = 1)
    #plt.scatter(t,MSD_theory,marker=".", s = 1)
    print(time)
    plt.plot(time,msd)
    plt.plot(time,MSD_theory)
    plt.xlabel("Time (s)")
    plt.ylabel("MSD")
    plt.legend(["Actual MSD","Theoretical MSD"], loc='upper right')
    plt.grid()
    plt.title("MSD for N={} particles".format(num_particles))
    plt.savefig(f'post_processing//MSD_{filename}.png')  # Save the plot to a file

# Usage:
# VS sets path from git -> "../src_written_data" not required
pos_data = read_data("src_written_data","pos.txt")
#plot_data(pos_data)

MSD_plot(pos_data)