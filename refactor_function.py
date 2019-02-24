# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 14:38:08 2019

@author: Jason Tanuwijaya
"""

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