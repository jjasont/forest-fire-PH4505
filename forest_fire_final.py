# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 15:13:03 2018

@author: Jason (U1440158A) & Victor Getty
For the completion of PH4505 Course Project of the topic Forest-Fire Simulation
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(color_codes=True)
from pylab import rcParams
rcParams['figure.figsize'] = 15, 7.5

#%%
FIRE = -1
EMPTY = 0
TREE = 1

def neighbour_entry(i, j, forest_size, nn = 'N', pbc = False):
    if nn == 'N':
        neighbor = np.array([[(i-1), j],[i, (j-1)],[(i+1), j],[i, (j+1)]])
    elif nn == 'M':
        neighbor = np.array([[(i-1), j],[i, (j-1)],[(i+1), j],[i, (j+1)], [(i-1), (j-1)],[(i+1), (j-1)],[(i+1), (j-1)],[(i-1), (j+1)]])
    if pbc:
        neighbor = np.mod(neighbor, forest_size)
    return neighbor

def find_radius(cluster_list):
    x_var,y_var = np.var(cluster_list,axis=0)
    radius = np.sqrt(x_var+y_var)
    return radius

def flatten_clusterdata(cluster_data):
    lists = sorted(cluster_data.items()) # sorted by key, return a list of tuples

    cluster_size, y_temp = zip(*lists) # unpack a list of pairs into two tuples
    cum_radius, number_cluster = zip(*y_temp)
    radius_cluster = np.array(cum_radius)/np.array(number_cluster)
    
    # =========== Primitive Method ==========
    #    cluster_size = np.zeros(len(cluster_data))
    #    number_cluster = np.zeros(len(cluster_data))
    #    radius_cluster = np.zeros(len(cluster_data))
    #    i = 0
    #    for key, value in cluster_data.items():
    #        cluster_size[i] = key
    #        number_cluster[i] = value[1]
    #        radius_cluster[i] = value[0]/value[1]
    #        i += 1
    # =======================================
    return cluster_size, number_cluster, radius_cluster

def cluster_distribution(forest_grid, pbc = False):
    forest_size = np.shape(forest_grid)
    cluster_index = np.full(forest_size, np.inf)
    
    forest_index = {}
    index = 0
    TREE = 1
    for i in range(0,forest_size[0]):
        for j in range(0, forest_size[1]):
            if forest_grid[i,j] == 1:
                up_left_down_right = neighbour_entry(i, j, forest_size)


                # Determine Valid Neighborhood
                if not pbc: # Only consider index range from 0 to L-1
                    pos_idx = np.all(up_left_down_right >= 0, axis = 1) #neighbor with positive index
                    in_bound_idx = (up_left_down_right[:, 0] <= forest_size[0]-1) & (up_left_down_right[:, 1] <= forest_size[1]-1) #neighbor within size of lattice
                    idx_true = pos_idx & in_bound_idx
                    nbor = up_left_down_right[idx_true]
                else:
                    nbor = up_left_down_right


                # Determined cluster indexation of neighbor
                clust_nbor = cluster_index[(nbor[:,0],nbor[:,1])]
                min_index = min(clust_nbor)
                
                if np.isinf(min_index):
                    cluster_index[i,j] = index
                    min_index = index
                    index += 1
                else:
                    cluster_index[i,j] = min_index
                
                if forest_index.get(min_index) is None:
                    forest_index[min_index] = [[i,j]]
                else:
                    value_temp = forest_index[min_index]
                    value_temp.append([i,j])
                    forest_index[min_index] = value_temp
                
                for coor, index_nbor in zip(nbor,clust_nbor): #loop for the existing neighbor and the value of cluster index at existing neighbor
                    if (forest_grid[coor[0], coor[1]] == TREE) & (index_nbor > min_index):
                        cluster_index[coor[0], coor[1]] = min_index
                        if forest_index.get(index_nbor) != None:
                            value_temp = forest_index[min_index]
                            for member in forest_index[index_nbor]:
                                cluster_index[member[0],member[1]] = min_index
                                value_temp.append([member[0],member[1]])
                            forest_index[min_index] = value_temp
                            forest_index.pop(index_nbor, None)
    #cluster_size_dict = {}
    cluster_radius_dict = {}
    
    for key,value in forest_index.items():
        cluster_size = len(value)
        radius = find_radius(value)
        if cluster_radius_dict.get(cluster_size) is None:
            #cluster_size_dict[cluster_size] = 1
            cluster_radius_dict[cluster_size] = np.array([radius,1])
        else:
            value_radius = cluster_radius_dict[cluster_size]
            value_radius += np.array([radius,1])
            cluster_radius_dict[cluster_size] = value_radius
    return cluster_radius_dict


def cell_update(forest_grid, forest_size, p, f, equil_iter, pbc = False, immune = 0):
    #Forest_grid_temp is used to create the forest grid at next time step.
    forest_grid_temp = forest_grid

    
    num_tree = 0
    num_empty = 0
    num_fire = 0
    for i in range(forest_size[0]):
        for j in range(forest_size[1]):
            if forest_grid_temp[i,j] == FIRE:
                num_fire += 1
                forest_grid[i,j] = EMPTY
            elif forest_grid_temp[i,j] == EMPTY:
                num_empty +=1
                if np.random.random() <= p:
                    forest_grid[i,j] = TREE
            elif forest_grid_temp[i,j] == TREE:
                num_tree +=1
                neighbor = neighbour_entry(i, j, forest_size, nn = 'N', pbc = pbc)
                if not pbc:
                    pos_idx = np.all(neighbor >= 0, axis = 1) #neighbor with positive index
                    in_bound_idx = (neighbor[:, 0] <= forest_size[0]-1) & (neighbor[:, 1] <= forest_size[1]-1) #neighbor within size of lattice
                    idx_true = pos_idx & in_bound_idx
                    nbor = neighbor[idx_true]
                else:
                    nbor = neighbor
                
                for item in nbor:
                    if forest_grid[item[0],item[1]] == FIRE:
                        ##
                        # Snippet of immune to fire code
                        if immune == 0: #next to fire directly burnt
                            forest_grid[i,j] = FIRE
                            break
                        else:
                            if np.random.random() > immune: #resistance/immune of tree
                                forest_grid[i,j] = FIRE
                                break
                        ##
                        
                    else:
                        if np.random.random() <= f:
                            forest_grid[i,j] = FIRE
    
    ####
    #cluster measurement
    if equil_iter >= 0:
        cluster_radius_dict = cluster_distribution(forest_grid, pbc = pbc)
        return num_tree, num_empty, num_fire, forest_grid, cluster_radius_dict
    else:
        return num_tree, num_empty, num_fire, forest_grid
    ####
    
    
    
    

def forest_fire_main(forest_size,p,pbf,iteration, pbc = False, sweep = False, immune = 0, warm_up = 200):
#Here we implement the forest fire model using the periodic boundary cond.
# N- it defines the size of 2-dimensional grid
# p- Tree will be grown at an empty site if random number generated is less 
#than p 
# pbf is the ratio of p and f where f is defined below
# f- A green tree will start burning if either of its neighbours is already
#burning or a random number generated is less than the lightening parameter f 
#Iter- it defines the number of time steps we want to evolve the system.
#If gridbit is 1, the inital grid will be the one given as Forest_grid,
#otherwise it is generated with 50-50 prob to grow a tree or to have empty
#site

#0-tree
#1-empty site
#2-fire

#This code is specifically written to calculate various distributions at
#every step and then average over all of them.
#This is important since in the SOC state, the distributions change in every step. 

#It_steps store is a vector which stores the correct time for measuring the
#various distributions (associated index is "me")
    
#    FIRE = -1
#    EMPTY = 0
#    TREE = 1
    #pbc = False
    import os
    if immune == 0:
        directory = 'simulation_result\\forest-fire_' + str(forest_size[0]) + 'x' + str(forest_size[1]) + '_' + str(p) + '_' + str(pbf) + '_' + str(iteration)
    else:
        directory = 'simulation_result\\forest-fire_' + str(forest_size[0]) + 'x' + str(forest_size[1]) + '_' + str(p) + '_' + str(pbf) + '_' + str(iteration) + '_immune' + str(immune)
    if not os.path.exists(directory):
        os.makedirs(directory)
    #warm_up_step = 199
    #only calculate the radius,cluster distribution after iteration > warm_up_step
    #neighbourhood = ((-1,0), (0,-1), (0, 1), (1,0))
    f = p/pbf
    #To accomodate periodic boundary condition
    forest_grid = np.random.choice([0,1],size=(forest_size))
    #K = 0
    #number = 0
    #radius = 0
    #me = 0
    trees = np.zeros(iteration)
    empty = np.zeros(iteration)
    fires = np.zeros(iteration)
    cluster_data = {}
    for iterate in range(iteration):
        
        #Forest_grid_temp is used to create the forest grid at next time step.
        #forest_grid_temp = forest_grid
        #grid = np.random.uniform(size=(N,N))
        
        #At every step, a matrix of random numbers of size N x N is generated.
        #It is used in 2 scenarios: either a green tree (lightening
        #probability) or an empty site (probability to grow tree)    
        
        #index = 0
        #index_grid = np.zeros((N,N))
#        num_tree = 0
#        num_empty = 0
#        num_fire = 0
        equil_iter = iterate-warm_up
        if (not sweep) & np.any(np.array(forest_size) > 200):
            if (equil_iter >= 0) & ((equil_iter)%100 == 0):
                print('Iteration {0} completed'.format(iterate+1))
            elif (equil_iter == -1):
                print('Equilibration completed')
                
        if equil_iter >= 0:
            num_tree, num_empty, num_fire, forest_grid, cluster_radius_dict = cell_update(forest_grid, forest_size, p, f, equil_iter, pbc = pbc, immune = immune)
            for key,value in cluster_radius_dict.items():
                value_temp = cluster_radius_dict[key]
                if cluster_data.get(key) == None:
                    cluster_data[key] = value_temp
                else:
                    value_cluster_data = cluster_data[key]
                    cluster_data[key] = value_cluster_data + value_temp
        else:
            #print('no_clust')
            num_tree, num_empty, num_fire, forest_grid = cell_update(forest_grid, forest_size, p, f, equil_iter, pbc = pbc, immune = immune)
        trees[iterate] = num_tree
        empty[iterate] = num_empty
        fires[iterate] = num_fire
        #dict3 = defaultdict(list)
        
        #cluster_data value consist of 2 value.
        #array[0] is the cumulative sum of radius (w/ given size) across iteration and configuration
        #array[1] is the cumulative sum of number of cluster (w/ given size) across iteration and configuration
    
    cluster_size, number_cluster, radius_cluster = flatten_clusterdata(cluster_data)        
    #plt.plot([trees, empty, fires])
    #forest-fire_rowsxcols_p_pbf_iter
    file_name = directory+'/forest-fire_' + str(forest_size[0]) + 'x' + str(forest_size[1]) + '_' + str(p) + '_' + str(pbf) + '_' + str(iteration) + '_immune' + str(immune)
    
    #file_name_number = file_name + 'number.png'
#    fig1, ax1 = plt.subplot()
#    ax1.plot(trees, 'g', label = 'Tree')
#    ax1.plot(empty, 'k', label = 'Empty')
#    ax1.plot(fires, 'r', label = 'Fire')
#    ax1.set_xlabel('Time')
#    ax1.legend()
#    ax1.savefig(file_name + '_number.png', bbox_inches='tight')
#    ax1.show()
    
    fig1, (ax1) = plt.subplots(nrows = 1, ncols = 1, figsize = [15,7.5])
    ax1.plot(trees, 'g', label = 'Tree')
    ax1.plot(empty, 'k', label = 'Empty')
    ax1.plot(fires, 'r', label = 'Fire')
    ax1.set_xlabel('Time')
    ax1.legend()
    fig1.savefig(file_name + '_number.png', bbox_inches='tight')
    if not sweep:
        fig1.show()
    
    fig2, (ax2) = plt.subplots(nrows = 1, ncols = 1, figsize = [15,7.5])
    ax2.plot(trees/np.size(forest_grid), 'g', label = 'Tree Density')
    ax2.plot(empty/np.size(forest_grid), 'k', label = 'Empty Density')
    ax2.plot(fires/np.size(forest_grid), 'r', label = 'Fire Density')
    ax2.set_xlabel('Time')
    ax2.legend()
    fig2.savefig(file_name + '_density.png', bbox_inches='tight')
    fig2.show()
    
    
    fig3, (ax3) = plt.subplots(nrows = 1, ncols = 1, figsize = [15,7.5])
    ax3.plot(cluster_size, number_cluster)
    ax3.set_title('Cluster Size vs Number of Cluster')
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('Cluster Size, s')
    ax3.set_ylabel('Number of Cluster, N(s)')
    fig3.savefig(file_name + '_numcluster.png', bbox_inches='tight')
    fig3.show()
    
    fig4, (ax4) = plt.subplots(nrows = 1, ncols = 1, figsize = [15,7.5])
    ax4.plot(cluster_size[1:], radius_cluster[1:])
    ax4.set_title('Cluster Size vs Radius of Cluster')
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('Cluster Size, s')
    ax4.set_ylabel('Radius of Cluster, R(s)')
    fig4.savefig(file_name + '_Rcluster.png', bbox_inches='tight')
    fig4.show()
    
    if not sweep:
        fig1.show()
        fig2.show()
        fig3.show()
        fig4.show()
    if sweep:
        plt.close(fig1)
        plt.close(fig2)
        plt.close(fig3)
        plt.close(fig4)
    
    np.savez(file_name, p = p, pbf = pbf, forest_size = forest_size, iteration = iteration, pbc = pbc, warm_up = warm_up, trees = trees, empty = empty, fires = fires, cluster_size = cluster_size, number_cluster = number_cluster, radius_cluster = radius_cluster)
    return [trees, empty, fires, cluster_size, number_cluster, radius_cluster]
                    
#                    for dx,dy in neighbourhood:
#                        if forest_grid[(i+dx)%N,(j+dy)%N] == FIRE:
#                            forest_grid = FIRE
#                            break
#                        
#                        
#                        
#                    else:
#                        if np.random.random() <= f:
#                            forest_grid[i,j] = FIRE
                    #if iteration > warm_up_step:
                                    
#%% VARY SIZE
#% fix p and fix pbf
p = 0.05
pbf = 100
forest_size = [[10,10],[25,25],[50,50],[100,100]]
iteration = 1000
result = []
sweep = True
#fig1 = plt.figure(figsize = [15,7.5])
for size, num in zip(forest_size, range(1,len(forest_size)+1)):
    result_temp = forest_fire_main(size,p,pbf,iteration, pbc = False, sweep = sweep)
    print('The {0}-th forest size completed'.format(num))
    result.append(result_temp)
    if num == 1:
        fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
    
    ax1.plot(result_temp[3],result_temp[4], label = 'Forest Size = ' + str(size[0]) + 'x' + str(size[1]))
    ax1.legend()
    ax2.plot(result_temp[3][1:],result_temp[5][1:], label = 'Forest Size = ' + str(size[0]) + 'x' + str(size[1]))
    ax2.legend()

ax1.set_xlabel('Cluster Size, s')
ax1.set_ylabel('Number of Cluster, N(s)')
ax1.set_title('Number of Cluster vs Cluster Size, p = ' + str(p) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xlabel('Cluster Size, s')
ax2.set_ylabel('Radius of Cluster, R(s)')
ax2.set_title('Radius of Cluster vs Cluster Size, p = ' + str(p) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax2.set_xscale('log')
ax2.set_yscale('log')

fig.savefig('forest-size vary, fix p, pbf.png', bbox_inches='tight')
fig.show()

np.savez('forest-size vary, fix p, pbf', p = p, pbf = pbf, forest_size = forest_size, iteration = iteration, result = result, data = ['trees', 'empty', 'fires', 'cluster_size', 'number_cluster', 'radius_cluster'])

#%% VARY PBF
#% fix p and size
p = 0.05
pbf = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000]
forest_size = [100,100]
iteration = 1000
result = []
sweep = True
#fig1 = plt.figure(figsize = [15,7.5])
for pbf_trial, num in zip(pbf, range(1,len(pbf)+1)):
    result_temp = forest_fire_main(forest_size,p, pbf_trial, iteration, pbc = False, sweep = sweep)
    print('The {0}-th pbf completed'.format(num))
    result.append(result_temp)
    if num == 1:
        fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
    
    ax1.plot(result_temp[3],result_temp[4], label = 'pbf = ' + str (pbf_trial))
    ax1.legend()
    ax2.plot(result_temp[3][1:],result_temp[5][1:], label = 'pbf = ' + str (pbf_trial))
    ax2.legend()


ax1.set_xlabel('Cluster Size, s')
ax1.set_ylabel('Number of Cluster, N(s)')
ax1.set_title('Number of Cluster vs Cluster Size, p = ' + str(p)  + ', Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', iteration = ' + str (iteration))
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xlabel('Cluster Size, s')
ax2.set_ylabel('Radius of Cluster, R(s)')
ax2.set_title('Radius of Cluster vs Cluster Size, p = ' + str(p)  + ', Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', iteration = ' + str (iteration))
ax2.set_xscale('log')
ax2.set_yscale('log')

fig.savefig('fix forest-size, fix p, vary pbf.png', bbox_inches='tight')
fig.show()

np.savez('fix forest-size, fix p, vary pbf', p = p, pbf = pbf, forest_size = forest_size, iteration = iteration, result = result, data = ['trees', 'empty', 'fires', 'cluster_size', 'number_cluster', 'radius_cluster'])

#%% VARY SIZE AND PBF
#% Fix p
p = 0.5
pbf = [1, 5, 10, 50, 100, 500, 1000, 5000, 10000]
#pbf = [1]
forest_size = [[10,10],[25,25],[50,50],[100,100]]
#forest_size = [[10,10]]
iteration = 1000
result = []
num = 1
sweep = True
#fig1 = plt.figure(figsize = [15,7.5])
for pbf_trial in pbf:
    for size in forest_size:
        result_temp = forest_fire_main(size,p,pbf_trial,iteration, pbc = False, sweep = sweep)
        print('The {0}-th forest size completed'.format(num))
        result.append(result_temp)
        if num == 1:
            fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
        
        ax1.plot(result_temp[3],result_temp[4], label = 'Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', pbf = ' + str (pbf_trial))
        #ax1.legend()
        ax2.plot(result_temp[3][1:],result_temp[5][1:], label = 'Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', pbf = ' + str (pbf_trial))
        #ax2.legend()
        num += 1

ax1.set_xlabel('Cluster Size, s')
ax1.set_ylabel('Number of Cluster, N(s)')
ax1.set_title('Number of Cluster vs Cluster Size, p = ' + str(p)  + ', iteration = ' + str (iteration))
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xlabel('Cluster Size, s')
ax2.set_ylabel('Radius of Cluster, R(s)')
ax2.set_title('Radius of Cluster vs Cluster Size, p = ' + str(p) + ', iteration = ' + str (iteration))
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend(loc=7, bbox_to_anchor=(1.45, 0.5))

fig.savefig('forest-size vary, pbf vary, fix p.png', bbox_inches='tight')
fig.show()

np.savez('forest-size vary, pbf vary, fix p', p = p, pbf = pbf, forest_size = forest_size, iteration = iteration, result = result, data = ['trees', 'empty', 'fires', 'cluster_size', 'number_cluster', 'radius_cluster'])
#%%
#%% Immunity
#% Fix p
p = 0.05
pbf = 10000
immune_prob = np.linspace(0,1,21)#[1, 5, 10, 50, 100, 500, 1000, 5000, 10000]
#pbf = [1]
forest_size = [100,100]
#forest_size = [[10,10]]
iteration = 1000
result = []
num = 1
sweep = True
#fig1 = plt.figure(figsize = [15,7.5])
for immune in immune_prob:
    result_temp = forest_fire_main(size,p,pbf_trial,iteration, pbc = False, sweep = sweep, immune = immune)
    
    result.append(result_temp)
    if num == 1:
        fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
    
    ax1.plot(result_temp[3],result_temp[4], label = 'Probability Immune = ' + str(immune) )
    #ax1.legend()
    ax2.plot(result_temp[3][1:],result_temp[5][1:], label = 'Probability Immune = ' + str(immune))
    print('The {0}-th forest size completed'.format(num))
    #ax2.legend()
    num += 1

ax1.set_xlabel('Cluster Size, s')
ax1.set_ylabel('Number of Cluster, N(s)')
ax1.set_title('Number of Cluster vs Cluster Size, p = ' + str(p)  + ', Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', pbf = ' + str (pbf_trial) + ', iteration = ' + str (iteration))
ax1.set_xscale('log')
ax1.set_yscale('log')


ax2.set_xlabel('Cluster Size, s')
ax2.set_ylabel('Radius of Cluster, R(s)')
ax2.set_title('Radius of Cluster vs Cluster Size, p = ' + str(p) + ', Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', pbf = ' + str (pbf_trial) + ', iteration = ' + str (iteration))
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend(loc=7, bbox_to_anchor=(1.45, 0.5))

fig.savefig('vary_immune1.png', bbox_inches='tight')
fig.show()

np.savez('vary_immune1', p = p, pbf = pbf, immune = immune_prob, forest_size = forest_size, iteration = iteration, result = result, data = ['trees', 'empty', 'fires', 'cluster_size', 'number_cluster', 'radius_cluster'])

# ==============================================================================
#% Average Density for different immunity rate

for i in range(0,len(immune_prob)):
    if i == 0:
        fig2, (ax1) = plt.subplots(nrows = 1, ncols = 1, figsize = [20,10])
        density_immune = []
        density_immune_err = []
    density_immune.append(np.mean(result[i][0][200:]/1e4))
    density_immune_err.append(np.std(result[i][0][200:]/1e4))


ax1.errorbar(immune_prob,density_immune, yerr = density_immune_err, marker = 's')
ax1.set_xlabel('Immune Rate/Probability')
ax1.set_ylabel('Density of Tree, rho_{tree}')
ax1.set_title('Density of Tree vs Immune Probability, p = ' + str(p)  + ', Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
fig2.savefig('DensityVsImmune.png.png', bbox_inches='tight')
fig2.show()

# ==============================================================================
#% Average Number of Cluster for different immunity rate

for i in range(0,len(immune_prob)):
    if i == 0:
        fig0, (ax3) = plt.subplots(nrows = 1, ncols = 1, figsize = [20,10])
        ncluster_immune = []
        ncluster_immune_err = []
    ncluster_immune.append(np.mean(np.array(result[i][4][1:])))
    ncluster_immune_err.append(np.std(np.array(result[i][4][1:])))

ax3.errorbar(immune_prob, ncluster_immune, yerr = ncluster_immune_err, marker = 's')
ax3.set_xlabel('Immune Rate/Probability')
ax3.set_ylabel('Average Number of Cluster')
ax3.set_title('Average Number of Cluster vs Immune Probability, p = ' + str(p)  + ', Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax3.set_yscale('log')
fig0.savefig('AverageNumClusterVsImmune.png', bbox_inches='tight')
fig0.show()

#%% Average Cluster size for different immunity rate

for i in range(0,len(immune_prob)):
    if i == 0:
        fig1, (ax3) = plt.subplots(nrows = 1, ncols = 1, figsize = [20,10])
        density_immune = []
        density_immune_err = []
        density_max =[]
        density_min = []
    norm = np.sum(np.array(result[i][4][1:]))
    clust_sizetemp = np.array(result[i][3][1:])
    density_max.append(np.max(clust_sizetemp))
    density_min.append(np.min(clust_sizetemp))
    mean_temp = np.sum(np.array(result[i][4][1:])*clust_sizetemp/norm)
    density_immune.append(mean_temp)
    density_immune_err.append(np.sqrt(np.mean((np.repeat(clust_sizetemp, np.array(result[i][4][1:], dtype = int)) - mean_temp)**2)))

ax3.plot(immune_prob, density_max, marker = 's', label = 'Max')
ax3.plot(immune_prob, density_min, marker = 's', label = 'Min')
ax3.errorbar(immune_prob, density_immune, yerr = density_immune_err, marker = 's', label = 'Mean +/- SD')
ax3.set_xlabel('Immune Rate/Probability')
ax3.set_ylabel('Cluster Size')
ax3.set_title('Cluster Size vs Immune Probability, p = ' + str(p)  + ', Forest Size = ' + str(size[0]) + 'x' + str(size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax3.legend()
#ax3.set_yscale('log')
ax3.set_yscale('log')
fig1.savefig('Cluster Size vs immune2.png', bbox_inches='tight')
fig1.show()

#%% Critical Exponent Calculation
#p = files5['p']
#pbf = files5['pbf']
#forest_size = files5['forest_size']
#iteration = files5['iteration']
#tree_density = files5['trees']/(forest_size[0] * forest_size[1])
#cluster_size = files5['cluster_size']
#number_cluster = files5['number_cluster']
#radius_cluster = files5['radius_cluster']

tree_density = result[0]/(forest_size[0] * forest_size[1])
cluster_size = result[3]
number_cluster = result[4]
radius_cluster = result[5]

from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a * x**(-b) + c
popt, pcov = curve_fit(func, cluster_size[0:500], number_cluster[0:500])
tau1 = popt[1]
tau1_err = np.sqrt(np.diag(pcov))[1]

def func1(x, a, b, c = 0):
    return a * x**(1/b) + c
popt1, pcov1 = curve_fit(func1, cluster_size[1:1000], radius_cluster[1:1000])
mu1 = popt1[1]
mu1_err = np.sqrt(np.diag(pcov1))[1]

#nu1 = 1/mu1
#nu1_err = mu1_err/nu1**2
s = pbf * (1-(tree_density[200:]))/(tree_density[200:])
s_average1 = np.mean(s)
s_err1 = np.std(s)
Aprime = np.log(s_average1)
Aprime_err = s_err1/Aprime
nu3 = np.log(func1(s_average1, popt1[0], popt1[1], 0))/np.log(1e4) #modify func1 to x**b
nu3_err = nu3 *np.sqrt( (s_err1/s_average1)**2 + (np.sqrt(np.diag(pcov1))[1]/popt[1])**2 + ((np.diag(pcov1))[0] / popt1[0])**2) / np.log(1e4)
#lambda = np.log(np.max(cluster_size))/np.log(1e4)
lambda_ = np.mean([0.96253898060552301,0.9824823900211469,0.9864409807777933])
#Out[559]: 0.97715411713482103

lambda_err = np.std([0.96253898060552301,0.9824823900211469,0.9864409807777933])
print('Critical Exponent tau = {0:.5f} +/- {1:.5f}'.format(tau1, tau1_err))
print('Critical Exponent mu = {0:.5f} +/- {1:.5f}'.format(mu1, mu1_err))
print('Critical Exponent nu = {0:.5f} +/- {1:.5f}'.format(nu3,nu3_err))
print('Critical Exponent lambda = {0:.5f} +/- {1:.5f} (From system size 100, 250, 500)'.format(lambda_,lambda_err))


#%% Critical Exponent Calculation
#% Use 250x250 data
#d = 2
#start_from = -200
#mu = np.log(np.array(result_temp1[3][1:]))/np.log(np.array(result_temp1[5][1:]))
#mu_avg = np.mean(mu[start_from:])
#mu_std = np.std(mu[start_from:])
#nu = np.log(np.array(result_temp1[5][1:]))/np.log(10000) #10000 is pbf/zeta
#nu_avg = np.mean(nu[start_from:])
#nu_std = np.std(nu[start_from:])
#cov_mu_nu = np.cov(mu[start_from:],nu[start_from:])[0,1]
#lambda_ = mu_avg*nu_avg
#lambda_std = lambda_ * np.sqrt((mu_std/mu_avg)**2 + (nu_std/nu_avg)**2 + 2*cov_mu_nu/(mu_avg*nu_avg))
#tau_avg = d/mu_avg + 1
#tau_std = abs((d/mu_avg)*(-1)*mu_std/tau_avg**2)
#
#
#from scipy.optimize import curve_fit
#def func(x, a, b, c):
#    return a * x**(-b) + c
#popt, pcov = curve_fit(func, result_temp1[3][0:500], result_temp1[4][0:500])
#tau1 = popt[1]
#tau1_err = np.sqrt(np.diag(pcov1))[1]
#
#def func1(x, a, b, c = 0):
#    return a * x**(1/b) + c
#popt1, pcov1 = curve_fit(func1, result_temp1[3][1:1000], result_temp1[5][1:1000])
#mu1 = popt1[1]
#mu1_err = np.sqrt(np.diag(pcov1))[1]
#
##nu1 = 1/mu1
##nu1_err = mu1_err/nu1**2
#s = 1e4 * (1-(result_temp1[0][200:]/500**2))/(result_temp1[0][200:]/500**2)
#s_average1 = np.mean(s)
#s_err1 = np.std(s)
#Aprime = np.log(s_average1)
#Aprime_err = s_err1/Aprime
#nu3 = np.log(func1(s_average, popt1[0], popt1[1], 0))/np.log(1e4) #modify func1 to x**b
#nu3_err = nu3 *np.sqrt( (s_err1/s_average1)**2 + (np.sqrt(np.diag(pcov1))[1]/popt[1])**2 + ((np.diag(pcov1))[0] / popt1[0])**2) / np.log(1e4)
##lambda = np.log(np.max(result_temp[3]))/np.log(1e4)
#lambda_ = np.mean([0.96253898060552301,0.9824823900211469,0.9864409807777933])
##Out[559]: 0.97715411713482103
#
#lambda_err = np.std([0.96253898060552301,0.9824823900211469,0.9864409807777933])
##Out[560]: 0.010460059656368073