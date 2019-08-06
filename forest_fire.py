import numpy as np

from itertools import product
from collections import Counter
from tqdm import tqdm

class ForestFire:

    def __init__(self, grid_row_size, 
                 grid_col_size, tree_growth_param,
                 pbf, iteration_step,
                 periodic_boundary_condition=False,
                 immunity_rate=0,
                 warm_up_iteraton_step=200,
                 nn_type = 'N',
                 cluster_result=False):
        
        assert grid_row_size > 0, "Invalid value for grid_row_size. Only positive value"
        assert grid_col_size > 0, "Invalid value for grid_col_size. Only positive value"
        assert immunity_rate >= 0, "Invalid value for immunity_rate. Value range from 0 to 1"
        assert tree_growth_param >= 0, "Invalid value for tree_growth_param. Value range from 0 to 1"
        assert pbf > 0, "Invalid value for pbf. Only positive value"
        assert tree_growth_param <= pbf, "Invalid value for pbf. pbf should be larger than equal of tree_growth_param"

        # Supplied Parameter
        self.grid_row_size = grid_row_size
        self.grid_col_size = grid_col_size
        self.tree_growth_param = tree_growth_param
        self.pbf = pbf
        self.iteration_step = iteration_step
        self.periodic_boundary_condition = periodic_boundary_condition
        self.immunity_rate = immunity_rate
        self.warm_up_iteraton_step = warm_up_iteraton_step
        self.nn_type = nn_type
        self.cluster_result = cluster_result

        # Derived Parameter
        self.forest_size = [grid_row_size, grid_col_size]
        self.fire_prob_param = tree_growth_param/pbf
        self.forest_grid = np.random.choice([0,1],size=([grid_row_size, grid_col_size]))
        self.is_warm_up = 0 if warm_up_iteraton_step else 1

        # Constant Parameter
        self.FIRE = -1
        self.EMPTY = 0
        self.TREE = 1

        print("Success")

    def cell_update(self, generate_cluster_result=False):
        forest_grid_temp = self.forest_grid
        forest_grid_new = self.forest_grid  # to be updated

        row_size = self.grid_row_size
        col_size = self.grid_col_size

        num_tree = 0
        num_empty = 0
        num_fire = 0

        flatten_grid = [cell_val for row in forest_grid_temp for cell_val in row]
        cell_val_count = Counter(flatten_grid)

        # Count the number of cases
        num_tree = cell_val_count.get(self.TREE, 0)
        num_empty = cell_val_count.get(self.EMPTY, 0)
        num_fire = cell_val_count.get(self.FIRE, 0)

        # Work through each cell for update
        for i, (row, col) in tqdm(enumerate(product(range(row_size), range(col_size))), desc='Updating cell', ascii='True'):
            if forest_grid_temp[row, col] == self.FIRE:
                forest_grid_new[row, col] = self.EMPTY
            elif forest_grid_temp[row, col] == self.EMPTY:
                if np.random.random() <= self.tree_growth_param:
                    forest_grid_new[row, col] = self.TREE
            elif forest_grid_temp[row, col] == self.TREE:

                nbors = self.neighbour_entry(row, col)

                for nbor in nbors:
                    nbor_row = nbor[0]
                    nbor_col = nbor[1]
                    if self.forest_grid[nbor_row, nbor_col] == self.FIRE:
                        # Snippet of immune to fire code
                        if self.immunity_rate == 0: #next to fire directly burnt
                            forest_grid_new[row, col] = self.FIRE
                            break
                        else:
                            if np.random.random() > self.immunity_rate: #resistance/immune of tree
                                forest_grid_new[row, col] = self.FIRE
                                break
                        
                    else:
                        if np.random.random() <= self.fire_prob_param:
                            forest_grid_new[row, col] = self.FIRE

        # Update forest_grid config
        self.forest_grid = forest_grid_new

        if self.cluster_result:
            cluster_radius_dict = self.cluster_distribution()
            return num_tree, num_empty, num_fire, cluster_radius_dict
        else:
            return num_tree, num_empty, num_fire, {}


        

    def neighbour_entry(self, row_index, col_index):
        list_neighbour = []
        if self.nn_type == 'N':
            list_neighbour = [[(row_index-1), col_index],
                              [row_index, (col_index-1)],
                              [(row_index+1), col_index],
                              [row_index, (col_index+1)]]
        elif self.nn_type == 'M':
            list_neighbour = [[(row_index-1), col_index], 
                              [row_index, (col_index-1)], 
                              [(row_index+1), col_index], 
                              [row_index, (col_index+1)],
                              [(row_index-1), (col_index-1)],
                              [(row_index+1), (col_index-1)],
                              [(row_index+1), (col_index-1)],
                              [(row_index-1), (col_index+1)]]

        neighbour = np.array(list_neighbour)

        if self.periodic_boundary_condition:
            neighbour = np.mod(neighbour, self.forest_size)

           # Proposed Neighboring 
           # row_edge_neighbour_cell = (neighbour[:,0] >= forest_size[0])
           # neighbour[row_edge_neighbour_cell, 0] = np.mod(neighbour[row_edge_neighbour_cell,0], forest_size[0])
           
           # col_edge_neighbour_cell = (neighbour[:,1] >= forest_size[1])
           # neighbour[col_edge_neighbour_cell, 1] = np.mod(neighbour[col_edge_neighbour_cell,1], forest_size[1])
        else:
            pos_idx = np.all(neighbour >= 0, axis = 1) # neighbor with positive index
            row_in_system = (neighbour[:, 0] <= self.grid_row_size-1) # row index in system
            col_in_system = (neighbour[:, 1] <= self.grid_col_size-1) # col index in system
            
            # neighbor within size of lattice
            in_bound_idx =  row_in_system & col_in_system
            idx_true = pos_idx & in_bound_idx
            
            neighbour = neighbour[idx_true]
            
        return neighbour

    def cluster_distribution(self):
        """
        Build cluster index and retrieve the distribution of the clusters
        
        Note: 
        - Not fully implementing pbc due to vagueness of cluster definition
          in pbc realm
          Example:
              [1, 1, 0, 1]
              [1, 0, 0, 0]
              [0, 1, 1, 0]
              [1, 0, 0, 1]
              
              with Von Neumann Nearest Neighbour
              non-PBC case: 5 cluster, use row_index and col_index directly 
                             to calculate cluster radius
              PBC case: 2 cluster, to calculate radius for cluster on the edge, 
                          need to shift until no cluster located across 
                          periodic boundary i.e.
                          
             [1, 1, 0, 1]     [1, 1, 1, 0]     [1, 1, 0, 0]
             [1, 0, 0, 0]     [0, 1, 0, 0]     [1, 1, 1, 0]
             [0, 1, 1, 0]  >  [0, 0, 1, 1]  >  [0, 1, 0, 0]
             [1, 0, 0, 1]     [1, 1, 0, 0]     [0, 0, 1, 1]
             
             Hence, indexing cluster in PBC case a little intricate
                          
              
        """
        print("Processing cluster")
        
        # Initialize cluster_index array
        row_size = self.grid_row_size
        col_size = self.grid_col_size
        forest_size = self.forest_size
        cluster_index = np.full(forest_size, np.inf)
        
        # Initialize parameter used for cluster indexing
        # First cluster found assigned indexation=0
        # forest_index will list down the cell coordinate belonging to the cluster
        # forest_index = {cluster_id : [[a,b],[a,b+1],[a,b-1],[a-1,b]] ...}
        cluster_indices_coordinates = {}
        indexation = 0
        
        for i, (row, col) in tqdm(enumerate(product(range(row_size), range(col_size))), desc='Calculting cluster cell', ascii='True'):
            if self.forest_grid[row, col] == self.TREE:
                nbors = self.neighbour_entry(row, col)

                # Determined cluster indexation of neighbor
                clust_nbors = cluster_index[(nbors[:,0],nbors[:,1])]
                min_index = min(clust_nbors) # retrieve minimum index
                
                if np.isinf(min_index): # assign index if minimum is infinity
                    cluster_index[row, col] = indexation
                    min_index = indexation
                    indexation += 1
                else:
                    cluster_index[row, col] = min_index
                
                # Store Cluster Index and Corresponding Coordinate
                if cluster_indices_coordinates.get(min_index) is None:  # New Index, record first coordinate
                    cluster_indices_coordinates[min_index] = [[row, col]]
                else:  # Existing index, add new coordinate
                    cluster_indices_coordinates[min_index] += [[row, col]]
                
                # Check for the cluster index of neighbouring cell 
                for nbor_coor, clust_index_nbor in zip(nbors, clust_nbors):
                    # Assign smallest cluster index if there exist neighbouring cell of smaller index cluster
                    if (self.forest_grid[nbor_coor[0], nbor_coor[1]] == self.TREE) & (clust_index_nbor > min_index):
                        cluster_index[nbor_coor[0], nbor_coor[1]] = min_index

                        if cluster_indices_coordinates.get(clust_index_nbor) is not None:  # check if there is a set of coordinates with the larger index
                            for member in cluster_indices_coordinates[clust_index_nbor]:
                                cluster_index[member[0],member[1]] = min_index
                            cluster_indices_coordinates[min_index] += cluster_indices_coordinates[clust_index_nbor]  # combine larger index cluster to the smaller one
                            cluster_indices_coordinates.pop(clust_index_nbor, None)  # remove larger cluster index

        cluster_radius_dict = {}
        
        # cluster_radius_dict has the following format
        # cluster_radius_dict = {cluster_size: 
        #                       array([sum of radius for all the cluster this size,
        #                              occurence]) ...}

        for cluster_index, cluster_coordinates in cluster_indices_coordinates.items():
            cluster_size = len(cluster_coordinates)
            radius = self.find_cluster_radius(cluster_coordinates)
            radius_sqr = radius**2
            
            if cluster_radius_dict.get(cluster_size) is None:
                #cluster_size_dict[cluster_size] = 1
                cluster_radius_dict[cluster_size] = np.array([radius,1, radius_sqr])
            else:
                value_radius = cluster_radius_dict[cluster_size]
                value_radius += np.array([radius,1, radius_sqr])
                cluster_radius_dict[cluster_size] = value_radius
        return cluster_radius_dict

    def find_cluster_radius(self, cluster_list):
        x_var,y_var = np.var(cluster_list,axis=0)
        radius = np.sqrt(x_var+y_var)
        return radius


if __name__ == '__main__':
    import sys
    a = ForestFire(500, 500, 0.5, 10, 1000, cluster_result=True)
    print(a.forest_grid)
    print(sys.getsizeof(a.forest_grid))
    num_tree, num_empty, num_fire, cluster = a.cell_update()
    print(a.forest_grid)
    print(num_tree, num_empty, num_fire, cluster)
    print(sys.getsizeof(a.forest_grid))