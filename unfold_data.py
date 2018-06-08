# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 15:13:03 2018

@author: Jason Tanuwijaya
"""
#%% Load packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes = True)

#%% Section to unfold data of fix p and pbf/zeta and varying system size.
files1 = np.load('forest-size vary, fix p, pbf2.npz')
# File Explanation
# *pbfN.npz
# N = none. p = 0.5, zeta = 10
# N = 1. p = 0.5, zeta = 100
# N = 2. p = 0.05, zeta = 10
# N = 3. p = 0.05, zeta = 100

#file_name1 = files1.files

p = files1['p']
pbf = files1['pbf']
forest_size = files1['forest_size']
iteration = files1['iteration']
result = files1['result']

fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
for num in range(0,len(forest_size)):
    ax1.plot(result[num][3],result[num][4], label = 'Forest Size = ' + str(forest_size[num][0]) + 'x' + str(forest_size[num][1]))
    ax1.legend()
    ax2.plot(result[num][3][1:],result[num][5][1:], label = 'Forest Size = ' + str(forest_size[num][0]) + 'x' + str(forest_size[num][1]))
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

#fig.savefig('forest-size vary, fix p, pbf3.png', bbox_inches='tight')
fig.show()
#%% Section to unfold data of fix p and system size and varying zeta/pbf.
files2 = np.load('fix forest-size, fix p, vary pbf.npz')
#file_name2 = files2.files


p = files2['p']
pbf = files2['pbf']
forest_size = files2['forest_size']
iteration = files2['iteration']
result = files2['result']

fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
for num in range(0,len(pbf)):
    ax1.plot(result[num][3],result[num][4], label = 'pbf = ' + str (pbf[num]))
    ax1.legend()
    ax2.plot(result[num][3][1:],result[num][5][1:], label = 'pbf = ' + str (pbf[num]))
    ax2.legend()
ax1.set_xlabel('Cluster Size, s')
ax1.set_ylabel('Number of Cluster, N(s)')
ax1.set_title('Number of Cluster vs Cluster Size, p = ' + str(p) + ', Forest Size = ' + str(forest_size[0]) + 'x' + str(forest_size[1]) + ', iteration = ' + str (iteration))
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xlabel('Cluster Size, s')
ax2.set_ylabel('Radius of Cluster, R(s)')
ax2.set_title('Radius of Cluster vs Cluster Size, p = ' + str(p) + ', Forest Size = ' + str(forest_size[0]) + 'x' + str(forest_size[1]) + ', iteration = ' + str (iteration))
ax2.set_xscale('log')
ax2.set_yscale('log')

#fig.savefig('forest-size vary, fix p, pbf3.png', bbox_inches='tight')
fig.show()

#%% Section to unfold data of fix p and varyizing system size and zeta/pbf.
files3 = np.load('forest-size vary, pbf vary, fix p1.npz')
# File Explanation
# *pN.npz
# N = none. p = 0.5
# N = 1. p = 0.05

#file_name3 = files2.files


p = files3['p']
pbf = files3['pbf']
forest_size = files3['forest_size']
iteration = files3['iteration']
result = files3['result']

fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
num = 0
for num_pbf in range(0,len(pbf)):
    for num_forest in range(0,len(forest_size)):
        result_temp = result[num]
        ax1.plot(result_temp[3],result_temp[4], label = 'Forest Size = ' + str(forest_size[num_forest][0]) + 'x' + str(forest_size[num_forest][1]) + ', pbf = ' + str (pbf[num_pbf]))
        ax2.plot(result_temp[3][1:],result_temp[5][1:], label = 'Forest Size = ' + str(forest_size[num_forest][0]) + 'x' + str(forest_size[num_forest][1]) + ', pbf = ' + str (pbf[num_pbf]))
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

#fig.savefig('forest-size vary, pbf vary, fix p.png', bbox_inches='tight')
fig.show()
#%% Section to unfold data of fix p, system size and zeta/pbf and varying immunity.

files4 = np.load('vary_immune1.npz')
# File Explanation
# *pN.npz
# N = none. p = 0.5
# N = 1. p = 0.05

#file_name3 = files2.files


p = files4['p']
pbf = files4['pbf']
forest_size = files4['forest_size']
iteration = files4['iteration']
immune = files4['immune']
result = files4['result']

fig, (ax1, ax2) = plt.subplots(nrows = 1, ncols = 2, figsize = [20,10])
for immune_num in range(0,len(immune)):
    ax1.plot(result[immune_num][3],result[immune_num][4], label = 'Probability Immune = ' + str(immune[immune_num]) )
    ax2.plot(result[immune_num][3][1:],result[immune_num][5][1:], label = 'Probability Immune = ' + str(immune[immune_num]))

ax1.set_xlabel('Cluster Size, s')
ax1.set_ylabel('Number of Cluster, N(s)')
ax1.set_title('Number of Cluster vs Cluster Size, p = ' + str(p)  + ', Forest Size = ' + str(forest_size[0]) + 'x' + str(forest_size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax1.set_xscale('log')
ax1.set_yscale('log')


ax2.set_xlabel('Cluster Size, s')
ax2.set_ylabel('Radius of Cluster, R(s)')
ax2.set_title('Radius of Cluster vs Cluster Size, p = ' + str(p) + ', Forest Size = ' + str(forest_size[0]) + 'x' + str(forest_size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend(loc=7, bbox_to_anchor=(1.45, 0.5))

#fig.savefig('vary_immune1.png', bbox_inches='tight')
fig.show()

# ======================================================================
#% Average Number of Cluster for Different Immunity Rate

fig0, (ax3) = plt.subplots(nrows = 1, ncols = 1, figsize = [20,10])
ncluster_immune = []
ncluster_immune_err = []
for i in range(0,len(immune)):
    ncluster_immune.append(np.mean(np.array(result[i][4][1:])))
    ncluster_immune_err.append(np.std(np.array(result[i][4][1:])))
    #, label = 'Probability Immune = ' + str(immune_prob[i]) )
    #ax1.legend()
    #ax2.plot(result[i][3][1:],result[i][5][1:], '.-', label = 'Probability Immune = ' + str(immune_prob[i]))
    #print('The {0}-th forest size completed'.format(num))
    #ax2.legend()
    num += 1

ax3.errorbar(immune, ncluster_immune, yerr = ncluster_immune_err, marker = 's')
#ax3.errorbar(immune_prob, density_immune, yerr = density_immune_err, marker = 's')
ax3.set_xlabel('Immune Rate/Probability')
ax3.set_ylabel('Average Number of Cluster')
#ax3.set_ylabel('Average Cluster Size')
ax3.set_title('Average Number of Cluster vs Immune Probability, p = ' + str(p)  + ', Forest Size = ' + str(forest_size[0]) + 'x' + str(forest_size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax3.set_yscale('log')
#ax3.set_yscale('log')
fig0.savefig('Average Number of Cluster vs immune2.png', bbox_inches='tight')
fig0.show()



# ============================================================
#% Average Cluster Size for Different Immunity Rate

fig1, (ax3) = plt.subplots(nrows = 1, ncols = 1, figsize = [20,10])
clust_size_immune = []
clust_size_immune_err = []
clust_size_max =[]
clust_size_min = []
for i in range(0,len(immune)):
    norm = np.sum(np.array(result[i][4][1:]))
    clust_sizetemp = np.array(result[i][3][1:])
    clust_size_max.append(np.max(clust_sizetemp))
    clust_size_min.append(np.min(clust_sizetemp))
    mean_temp = np.sum(np.array(result[i][4][1:])*clust_sizetemp/norm)
    clust_size_immune.append(mean_temp)
    clust_size_immune_err.append(np.sqrt(np.mean((np.repeat(clust_sizetemp, np.array(result[i][4][1:], dtype = int)) - mean_temp)**2)))

ax3.plot(immune, clust_size_max, marker = 's', label = 'Max')
ax3.plot(immune, clust_size_min, marker = 's', label = 'Min')
ax3.errorbar(immune, clust_size_immune, yerr = clust_size_immune_err, marker = 's', label = 'Mean +/- SD')
ax3.set_xlabel('Immune Rate/Probability')
ax3.set_ylabel('Cluster Size')
ax3.set_title('Cluster Size vs Immune Probability, p = ' + str(p)  + ', Forest Size = ' + str(forest_size[0]) + 'x' + str(forest_size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
ax3.legend()
ax3.set_yscale('log')
#fig1.savefig('Cluster Size vs immune2.png', bbox_inches='tight')
fig1.show()

# ============================================================
#% Average Forest Density for Different Immunity Rate

for i in range(0,len(immune)):
    if i == 0:
        fig2, (ax1) = plt.subplots(nrows = 1, ncols = 1, figsize = [20,10])
        density_immune = []
        density_immune_err = []
    density_immune.append(np.mean(result[i][0][200:]/1e4))
    density_immune_err.append(np.std(result[i][0][200:]/1e4))

ax1.errorbar(immune,density_immune, yerr = density_immune_err, marker = 's')
ax1.set_xlabel('Immune Rate/Probability')
ax1.set_ylabel('Density of Tree, rho_{tree}')
ax1.set_title('Density of Tree vs Immune Probability, p = ' + str(p)  + ', Forest Size = ' + str(forest_size[0]) + 'x' + str(forest_size[1]) + ', pbf = ' + str (pbf) + ', iteration = ' + str (iteration))
#fig2.savefig('density vs immune1.png', bbox_inches='tight')
fig2.show()