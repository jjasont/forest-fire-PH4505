# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 15:18:10 2019

@author: Jason Tanuwijaya
"""
import os
import numpy as np

FIRE = -1
EMPTY = 0
TREE = 1

import matplotlib.pyplot as plt
#import seaborn as sns
#
#sns.set(color_codes=True)
from pylab import rcParams
rcParams['figure.figsize'] = 15, 7.5

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

def neighbour_entry(i, j, forest_size, nn = 'N', pbc = False):
    """
        Retrieve the neighbouring element of indicated location
        
        Keyword arguments:
          i -- the row i of the forest grid, start from 0 (int)
          j -- the column j of the forest grid, start from 0 (int)
          forest_size -- the size of forest grid system (np.size())
        
          nn -- nearest neighbour used for implementation (default 'N')
            'N' : Von Neumann neighbour. 4-adjacent cells around the centre
            'M' : Moore neighbour. All 8 cells surrounding the centre
            
        Return:
          neighbour -- np.array listing down all the neighbouring indices
    """
    list_neighbour = []
    if nn == 'N':
        list_neighbour = [[(i-1), j],[i, (j-1)],[(i+1), j],[i, (j+1)]]
    elif nn == 'M':
        list_neighbour = [(i-1), j], [i, (j-1)], [(i+1), j], [i, (j+1)],  \
            [(i-1), (j-1)],[(i+1), (j-1)],[(i+1), (j-1)],[(i-1), (j+1)]
    
    neighbour = np.array(list_neighbour)
    if pbc:
        neighbour = np.mod(neighbour, forest_size)
    else:
        pos_idx = np.all(neighbour >= 0, axis = 1) # neighbor with positive index
        row_in_system = (neighbour[:, 0] <= forest_size[0]-1) # row index in system
        col_in_system = (neighbour[:, 1] <= forest_size[1]-1) # col index in system
        
        # neighbor within size of lattice
        in_bound_idx =  row_in_system & col_in_system
        idx_true = pos_idx & in_bound_idx
        
        neighbour = neighbour[idx_true]
        
    return neighbour

def find_radius(cluster_list):
    x_var,y_var = np.var(cluster_list,axis=0)
    radius = np.sqrt(x_var+y_var)
    return radius

def cluster_distribution(forest_grid, pbc = False):
    forest_size = np.shape(forest_grid)
    cluster_index = np.full(forest_size, np.inf)
    
    forest_index = {}
    index = 0
    TREE = 1
    for i in range(0,forest_size[0]):
        for j in range(0, forest_size[1]):
            if forest_grid[i,j] == 1:
                nbor = neighbour_entry(i, j, forest_size)

                # Determined cluster indexation of neighbor
                clust_nbor = cluster_index[(nbor[:,0],nbor[:,1])]
                min_index = min(clust_nbor)
                
                if np.isinf(min_index):
                    cluster_index[i,j] = index
                    min_index = index
                    index += 1
                else:
                    cluster_index[i,j] = min_index
                
                if forest_index.get(min_index) == None:
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

def cell_update(forest_grid, p, f, cluster_result=False, pbc = False, immune = 0):
    """
        Update the forest_grid according to the rule of forest-fire
        
        Key arguments:
         forest_grid -- The previous state of the forest
         p -- Parameter that control tree growth. 
              Tree will be grown at an empty site if random number 
              generated is less than p. Range from 0 to 1
         f -- Parameter that control lightning to struck tree
              A green tree will start burning if either of its neighbours 
              is already burning or a random number generated is less than 
              the lightening parameter f.
         cluster_result -- Flag for cluster_result request. Use accordingly,
                           cluster calculation is expensive (default False)
         pbc -- The flag parameter to use periodic boundary condition. Periodic
                boundary condition useful to approximate the system on a large
                (infinite) scale from small system
                 (default False)
         immune -- Parameter that control tree resistance. This custom
                   parameter introduced on the case where given neighbouring
                   tree are dying due fire/infection yet the tree observed are
                   not affected due any resistance factor presence such as
                   wet tree/mutated tree/etc. Range from 0 to 1 (default 0.0)
         
    """
    #Forest_grid_temp is used to create the forest grid at next time step.
    forest_grid_temp = forest_grid
    
    row_size = np.size(forest_grid, 0)
    col_size = np.size(forest_grid, 1)
    forest_size = np.shape(forest_grid)

    
    num_tree = 0
    num_empty = 0
    num_fire = 0
    
    for i in range(row_size):
        for j in range(col_size):
            if forest_grid_temp[i,j] == FIRE:
                num_fire += 1
                forest_grid[i,j] = EMPTY
            elif forest_grid_temp[i,j] == EMPTY:
                num_empty +=1
                if np.random.random() <= p:
                    forest_grid[i,j] = TREE
            elif forest_grid_temp[i,j] == TREE:
                num_tree +=1
                nbor = neighbour_entry(i, j, forest_size, nn = 'N', pbc = pbc)
                
                for item in nbor:
                    if forest_grid[item[0],item[1]] == FIRE:
                        # Snippet of immune to fire code
                        if immune == 0: #next to fire directly burnt
                            forest_grid[i,j] = FIRE
                            break
                        else:
                            if np.random.random() > immune: #resistance/immune of tree
                                forest_grid[i,j] = FIRE
                                break
                        
                    else:
                        if np.random.random() <= f:
                            forest_grid[i,j] = FIRE
    
    ####
    #cluster measurement
    if cluster_result:
        cluster_radius_dict = cluster_distribution(forest_grid, pbc = pbc)
        return num_tree, num_empty, num_fire, forest_grid, cluster_radius_dict
    else:
        return num_tree, num_empty, num_fire, forest_grid

def main(forest_size, p, pbf, iteration, \
    pbc = False, sweep = False, immune = 0, warm_up = 200):
    """
       Main forest function
       
       Keyword arguments:
         forest_size -- The size of the system of format [row_size, col_size]
         p -- Parameter that control tree growth. 
              Tree will be grown at an empty site if random number 
              generated is less than p. Range from 0 to 1
         pbf -- The ratio of parameter p and f.
                f is implicit parameter controlling the lightning parameter
                i.e. tree died.
                A green tree will start burning if either of its neighbours 
                is already burning or a random number generated is less than 
                the lightening parameter f.
         iteration -- The total number of time steps to observe the system
         pbc -- The flag parameter to use periodic boundary condition. Periodic
                boundary condition useful to approximate the system on a large
                (infinite) scale from small system
                 (default False)
         sweep -- The flag parameter to indicate whether function used to run
                  on queued multiple system size. When it is true, less verbose
                  output displayed (default False)
         immune -- Parameter that control tree resistance. This custom
                   parameter introduced on the case where given neighbouring
                   tree are dying due fire/infection yet the tree observed are
                   not affected due any resistance factor presence such as
                   wet tree/mutated tree/etc. Range from 0 to 1 (default 0.0)
         warm_up -- Number of iteration for the system considered to achieve
                    equilibrium state. Equilibriation is important as large
                    system tend to have different characteristic on the initial
                    stage of the system due to random initialization
                    (default 200)
    """
    
    # Initalize folder to store result
    folder_name = 'simulation_result\\forest-fire_{col_size}x{row_size}_{p}_{pbf}_{iteration}'
    directory = folder_name.format(col_size = forest_size[0], \
                                   row_size = forest_size[1], p = p, \
                                   pbf = pbf, iteration=iteration)
    if immune < 0:
        return "Invalid value for immune. immune range from 0 to 1"
    elif immune > 0:
        directory = directory + '_immune' + str(immune)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    
    # Initialize the f, lightining parameter value
    f = p/pbf
    #To accomodate periodic boundary condition
    forest_grid = np.random.choice([0,1],size=(forest_size))
    
    # Store the number of each state per iteration and cluster_data
    trees = np.zeros(iteration)
    empty = np.zeros(iteration)
    fires = np.zeros(iteration)
    cluster_data = {}
    
    # Iterate according to the number of iteration inputted
    for iterate in range(iteration+warm_up):
        
        # Iteration for Equlibrium
        equil_iter = iterate-warm_up
        
        # Only do this when not sweeping across different system size
        # and when the system size us large enough
        if (not sweep) & np.any(np.array(forest_size) > 200):
            if iterate < warm_up:
                print("Equilibration - Iteration {0} completed.".format(iterate+1))
            elif iterate == warm_up:
                print("Equlibration Completed")

                
        # Iteration post-equlibirum
        n_iterate_post_warmup = iterate-warm_up
        if iterate > warm_up:
            num_tree, num_empty, num_fire, forest_grid, cluster_radius_dict = cell_update(forest_grid, p, f, cluster_result=True, pbc = pbc, immune = immune)
            for key, value in cluster_radius_dict.items():
                value_temp = cluster_radius_dict[key]
                if cluster_data.get(key) is None:
                    cluster_data[key] = value_temp
                else:
                    value_cluster_data = cluster_data[key]
                    cluster_data[key] = value_cluster_data + value_temp
        else:
            #print('no_clust')
            num_tree, num_empty, num_fire, forest_grid = cell_update(forest_grid, p, f, pbc = pbc, immune = immune)
        
        
        trees[n_iterate_post_warmup] = num_tree
        empty[n_iterate_post_warmup] = num_empty
        fires[n_iterate_post_warmup] = num_fire
        
    cluster_size, number_cluster, radius_cluster = flatten_clusterdata(cluster_data)        