This is DeepCarb version 0.147

DeepCarb is a simple carbon cycle box model consisting of three boxes
(atmosphere, surface ocean, deep ocean) and two infinite reservoirs
(terrestrial silicate rock and deep ocean CaCO3). It is designed to
assess the impact of carbon cycle perturbations on atmospheric CO2 and
seawater carbonate chemistry, particularly over i) long timescales and
ii) in a Monte Carlo fashion with random permutations of key model
parameterisations.

By default, the model can be run in one of two ways. Setting the variable
nMC = 1 will result in a simulation run only once using the default set
of model parameterisations. Setting nMC = n will result in n model
simulations in which key parameterisations are randomly varied over wide
ranges. These parameterisations are i) the relationship between the rate 
of carbon drawdown from silicate weathering and temperature, ii) the
slope of the relationship between CaCO3 dissolution/preservation to the 
saturation state of the deep ocean box, and iii) the slope of the
relationship between CaCO3 production in the surface ocean and the
saturation state of the surface ocean box.

To maintain flexibility, this script is not a function. Before running
the model, at least four parameters need to be manually changed in the
script. These are:
  1) nMC - the number of simulations to run, as described above.
  2) parfor/for - use the parallel computing toolbox to speed the model
     up when running multiple simulations. To do so, manually change
     for to parfor on line 152.
  3) ageEnd - the time in years at which the simulation should stop. Note
     that rather than preallocating all matrices, the model loads a set of
     stable simulations which includes the necessary variables. These
     have a length of 15 Myr. If longer simulations are required, the
     matrices need to be modified in the accompanying data files from 
     which the initial stable conditions are loaded.
  4) The model needs to be perturbed in some way. Examples of how this
     can be done are given below andthroughout the script. For example, 
     atmospheric CO2 could be modified or prescribed (atCO2in), weathering 
     rates could be modified (w12), seawater [Ca] could be changed either 
     directly (oceanCa), or via the Ca/total cation ratio of the
     weathering product (weatherCa).

In addition, the script produces a plot of the results. Edit linlog to 
control whether the data are shown on a linear or logarithmic scale 
(1 = linear, 0  = log).

The results are contained in the cell array dataOut which has dimensions 
n*1. The format of each cell is an ageEnd*15 matrix, with the following 
columns:
  1) Total C in the atmosphere box (Gt)
  2) Total C in the surface ocean box (Gt)
  3) TAlk in the surface ocean box (umol/kg)
  4) DIC in the surface ocean box (umol/kg)
  5) pH in the surface ocean box (total scale)
  6) Omega calcite in the surface ocean box at 10.3 mM [Ca]
  7) TAlk in the deep ocean box (umol/kg)
  8) DIC in the deep ocean box (umol/kg)
  9) pH in the deep ocean box (total)
  10) Omega calcite in the deep ocean box at 10.3 mM [Ca]
  11) Seawater [Ca] (mM)
  12) C diffusion between the atmosphere and surface ocean box (Gt/year)
      Negative values represent net diffusion to the atmosphere
  13) Net CaCO3 burial/dissolution rate (moles/year)
      Negative values represent net burial
  14) Carbon drawdown as a result of silicate weathering (Gt/year)
  15) Change in CaCO3 export from the surface ocean box (moles/year)


This script depends on the following functions which are packaged here
but are not my work. Please give appropriate credit to the original
authors.
- A stripped-down version of co2sys ('co2syslite') is used, which runs
approximately 1000 times faster than co2sys. Please credit (e.g.) van
Heuven et al. (2011) CO2SYS v 1.1, MATLAB program developed for CO2 
system calculations.
- The function parforwait is required to provide a waitbar within a
parfor loop. See:
https://uk.mathworks.com/matlabcentral/fileexchange/71083-waitbar-for-parfor
- The parallel computing toolbox is required if parfor is used.
- The random number generator (normrnd) requires the statistics and
machine learning toolbox


------------------------------- EXAMPLES -------------------------------

Where these example recommend changing the values of nMC, linlog, ageEnd,
and the use or otherwise of the parallel computing toolbox, these can be
found on lines 135, 140, 136, and 156 respectively. In the latter case,
'for' simply needs to be replaced with 'parfor'.

1. Simulate the anthropogenic C release and run the model once using the
default settings. Set nMC = 1, linlog = 0 (to view the long and short
term effects, ageEnd = 1e6, and use a for loop rather than parfor.
Uncomment lines 274-275 to load the anthropogenic C release data, and
comment line 280 specifying no atmospheric CO2 perturbation.

2. Simulate the PETM and run the model with 100 random permutations of
the key model parameterisations. Set nMC = 100, linlog = 1, ageEnd =
2.5e5, and use the parallel computing toolbox (change the for loop on
line 156 to a parfor loop). Uncomment line 263 to perturb the model with
one possible PETM C-release scenario, and comment line 280 specifying no
atmospheric CO2 perturbation.

3. Prescribe an atmospheric CO2 doubling over a range of doubling rates.
Set nMC = 50, linlog = 0, ageEnd = 1e7, and use the parallel computing
toolbox. Uncomment the code on lines 339-347, which will vary the
doubling rate from ~10 to 1e6 years. Note that setting nMC>1 will result
in random selections of the key model parameterisations described above.
To maintain consistency between the doubling rate experiments, remove
this feature by uncommenting line 235.

4. Decrease seawater [Ca] starting from Eocene boundary conditions and
run a full Monte Carlo simulation. Set nMC = 1000, linlog = 1, ageEnd =
15e6, and use the parallel computing toolbox. Uncomment line 293 to drive
[Ca] down by 10 mM over 8 Myr, starting at model year 1 million, and then
let the model equilibrate. Comment out line 294 which specifies no
seawater [Ca] perturbation. Comment out lines 162-163 which load the 
modern stable conditions, and uncomment lines 165-166 to load Eocene
boundary conditions.

5. Drive CO2 down by increasing the weathering rate independently of
climate, and run a full Monte Carlo simulation. Set nMC = 1000, linlog =
1, ageEnd = 15e6, and use the parallel computing toolbox. Uncomment lines
383-389 to increase the silicate weathering rate by 5% over 4 Myr,
starting at model year 1 million, and then hold it at that point. Comment
out lines 162-163 which load the modern stable conditions, and uncomment
lines 165-166 to load Eocene boundary conditions.

------------------------------------------------------------------------
