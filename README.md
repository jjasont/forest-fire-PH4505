# forest-fire-PH4505
## PH4505 Computational Physics Project - Forest Fire Simulation

## Files Information
### `forest_fire_final.py`
- Python script to run the forest-fire simulation
- Also, run various parameter change of the system (system size, p, p/f)
- Calculate the critical exponent

### `unfold_data.py`
- Python script to 'unfold' the saved data and plotted the figure in the folder of relevant_file
- Calculate the critical exponent

### DensityVs_fp.png
- Plot of tree density vs f/p for a system size of 50 by 50 and various choice of f and p

### forest-size vary, fix p, pbf.png; forest-size vary, fix p, pbf1.png;forest-size vary, fix p, pbf2.png;forest-size vary, fix p, pbf3.png
- Plot of number of cluster given size and radius of cluster for different system size
while setting p and p/f to be the same across comparison
- __pbf.png is for p/f = 10, p = 0.5
- __pbf1.png is for p/f = 100, p = 0.5
- __pbf2.png is for p/f = 10, p = 0.05
- __pbf3.png is for p/f = 100, p = 0.05
- Corresponding simulation data saved in the similar naming

### fix forest-size, fix p, vary pbf.png; 
- Plot of number of cluster given size and radius of cluster for different p/f value
while setting p and system size to be the same across comparison
- Corresponding simulation data saved in the similar naming

### forest-size vary, pbf vary, fix p.png; forest-size vary, pbf vary, fix p1.png
- Plot of number of cluster given size and radius of cluster for different system size
and p/f value while setting p to be the same across comparison
- __pbf.png is for p = 0.5
- __pbf1.png is for p = 0.05

### vary_immune1.npz
- Simulation data of forest fire given various value of additional immunity probability (from 0 [basic forest fire] to 1 [all tree immune from neighborhood fire])
- The summarized simulated data plotted and saved in the following figure
AverageNumClusterVsImmune.png; (Average number of cluster as a function of immunity rate)
ClusterSizeVsImmune.png; (Average cluster size as a function of immunity rate)
DensityVsImmune.png; (Average tree density as a function of immunity rate)
- Data that involves immune rate included in simulation_result/immune_variation/immune###.
- immune###, the first number is number before decimal and the rest is after decimal

### simulation_result folder
- Each folder named according to the parameter of the simulation
- 'forest-fire_#Rowsx#Cols_p_pbf_iteration'
   * #Rows, #Cols : The number of rows and columns in the system
   * p		  : The value of tree growth probability
   * pbf	  : The p/f value (or zeta in Jason's report)
   * iteration	  : The number of iteration for the simulation result
- Each folder contains their corresponding simulation data and figure