# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 15:18:10 2019

@author: Jason Tanuwijaya
"""
import os


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
                    equilibrium state (default 200)
    """
    folder_name = 'simulation_result\\forest-fire_{col_size}x{row_size}_{p}_{pbf}_{iteration}'
    directory = folder_name.format(col_size = forest_size[0], \
                                   row_size = forest_size[1], p = p, \
                                   pbf = pbf, iteration=iteration)
    if immune != 0:
        directory = directory + '_immune' + str(immune)
    if not os.path.exists(directory):
        os.makedirs(directory)