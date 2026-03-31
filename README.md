# simulating-blockchains
We simulate consensus mining blockchains with propagation delay.



# **README.md draft**
The Abelian Sandpile Model is the common name for the original Bak-Tang-Wiesenfeld Model which was the first discovered example of a system to display Self-Organized Criticality. The Abelian Sandpile Model attempts to predict and model avalanches within a sandpile. This model was first introduced in a 1987 paper by Bak, Tang, Wisenfeld. Since then many extension of the model have been built and much research has gone into understanding the statistics of the model.

To learn more about The Abelian Sandpile Model and Self-Organized Criticality please see the oringial 1987 scientific paper, *Self-Organized Criticality* by Bak, Tang, Weisenfeld.

[1987 Bak, Tang, Wiesenfeld.](http://www.chialvo.net/Curso/Cordoba2005/ClasesPowerpoints/Presentacion9/PapersClase9/soc2.pdf)
___

# **About this code**
The following files contain code to simulate The Abelian Sandpile Model subject to various conditions. As well as code to simulate The General Wildfire Model which is an extension of The Abelian Wildfire Model.

- sandpile.py
- wildfire.py
- Functions.py
- InitialConditionComparison.py

This code is intended for scientific purpose use and hence there are various parameters available to tune for both the Sandpile Model and the General Wildfire Model. Output for both models consists of various data which is explained in detail below.
___

# **Pip Install Packages**
The following needs to be installed for the code to run.

- numpy
- scipy
- matplotlib

___

# **Explanation of Files and Output**

**1. <u> sandpile.py** </u>\
Run this file to simulate the Sandpile Model.
You will be asked for user input for the following questions which determine certain parameters of the Model.

* How many stationary state avalanche occurances do you need between 100 and 1000?:
* Choose length of sandpile grid between 10 and 50:
* Choose initial condition.
Randomized Configuration (0) or Level Configuration (1)? Enter 0 or 1:
* Choose level configuration between 0,1,2,3:

Output data:
- Four plots containing the distributions of sandpile hights 0, 1, 2, 3 for the given input.
- Plot containing the mean number of sandgrains per lattice site over time.
- Plot containing the distribution of the simulated times between avalanches against a geoemetric distribtuion of rate 0.41
- Four loglog plots containing frequency of avalanche metrics: Topples, Loss, Area, Length.


**2. <u> sandpile_InitialConditions.py** </u>\
Run this file to view how convergence to stationary state is independent of initial configruations.
We chose intial configurations of;
1. Lattice Uniformly Randomized on {0, 1, 2, 3}
2. Lattice Level Set for Level set chosen from {0, 1, 2, 3}

Output data:
- Plot showing mean number of sandgrains per site over time for each initial configuration.

**3. <u> wildfire.py** </u>\
Run this file to simulate the General Wildfire Model.\
There are various parameters to this model which are coded in at the top and are available for adjustment as commented.
* Length of grid dimension.
* Initial Tree density.
* Wind Direction or no Wind.
* Pg = Probability of site with no tree growing a tree.
* Pf = Probability of tree next to fire catching fire.
* Pe = Probability of tree on fire burning out and becoming an empty site.
* f = Probability of a tree catching fire by lightning strike.
* Number of Iterations

Please feel free to adjust the parameters to your liking.\
Some parameter combinations will make the fire burn out very fast, others will make the fire consume the grid quickly.

In order to simulate a stationary wildfire we recommend using the following probabilies as a starting point for adjustment.
* Pg = 0.2
* Pf = 0.9
* Pe = 0.5
* f = 0.005

Output data:
- Time series plot containing the proportion of grid with tree on fire, tree not on fire and empty sites.
- Snap shot of the fire after the specified iteraitons.
- Live animation of the wildfire.\
 If a wind direction was chosen then you will be able to view the affect of the wind on the model. The fire tends to drift towards the direction of the wind. This is a global behaviour over the whole grid even though it is only locally defined in our code.

**4. <u> Functions.py** </u>\
No need to run this file!\
This file contains various functions which have been used in the above files.
