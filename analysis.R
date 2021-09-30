

###############################################################
###   SPATIALLY EXPLICIT SIMULATION IN DYNAMIC LANDSCAPES   ###
###############################################################
#
#   Scripts required for full analysis:       'job_submission.slurm'; 
#                                             'batch_submission.slurm';
#                                             '2_VICFireFunctions_SPARTAN.R'
#
#   Data inputs required to run analysis:     'occ_speciesGroup.tif';
#                                             'det.method1_speciesGroup.tif';
#                                             'det.method2_speciesGroup.tif';
#                                             'det.method3_speciesGroup.tif';
#                                             'det.method4_speciesGroup.tif';
#                                             'species.det.csv';
#                                             'species.covar.csv';
#                                             'Z_speciesGroup_5k.CAZ_DEA.rank.compressed.tif';
#                                             'burnt.tif';
#                                             'remote.tif';
#                                             'vic.road.tif';
#                                             'exisitingSitesSPDF.shp';
#                                             'repeatVisits.csv.tif';
#                                             'method.cost.array.csv.tif';

###########################################################################################
###   STEP 1: Define command line arguments and indices. Set up the parameter           ###
###           space for running jobs on Spartan                                         ###
###########################################################################################

# Clear workspace
rm(list=ls())

# Create command line for interacting with job_submission.slurm script.
command_args<- commandArgs(trailingOnly = TRUE)
parameter_index<- as.numeric(command_args[1])

## Possible ID options for each parameter. 
max.budget_options = c(5000000,  9000000, 20000000)
all.existing.sites_options = c(TRUE)
group_options = c("Frogs", "Birds", "Mammals", "Reptiles")
repeat.scenario_options = c(1, 2)
method.scenario_options = c(1, 2)
survey.year_options = c(1, 2) 
site.ratio_options = c(0.5, 0.8)

# Define grid of all possible parameter combinations. MAKE SURE that the length of the testGrid corresponse to the sequence length in 'batch_submission.slurm' before proceeding.
testGrid<-expand.grid(max.budget_options,
                      group_options,    
                      repeat.scenario_options,
                      method.scenario_options,
                      survey.year_options,
                      site.ratio_options,
                      all.existing.sites_options)


# Create indices for interacting with job_submission.slurm script.
max.budget<-testGrid[parameter_index, 1] 
group<- testGrid[parameter_index, 2]
repeat.scenario<- testGrid[parameter_index, 3]
method.scenario<- testGrid[parameter_index, 4]
survey.year<- testGrid[parameter_index, 5]
site.ratio<- testGrid[parameter_index, 6]
all.existing.sites<- testGrid[parameter_index, 7] 

colnames(testGrid)<- c("max.budget_options", "group_options", "repeat.scenario_options", "method.scenario_options", "survey.year_options", "site.ratio_options", "all.existing.sites_options")

# Check for package dependices and install/load in packages as required. -> Remove when in package <-
.libPaths("/home/asmart1/R/lib/3.6")
lib = .libPaths()[1]
repo<- "https://cran.ms.unimelb.edu.au/" # Set mirror to download packages.
pkgs = c("unmarked", "raster", "rgdal", "doParallel", "foreach", "abind", "sf", "stringr", "dplyr", "parallel", "future", "doFuture") # Define packages.
new.pkgs <- pkgs[!(pkgs %in% installed.packages(lib=lib)[,"Package"])]
if(length(new.pkgs)) install.packages(new.pkgs, lib=lib, repos=repo)
inst = lapply(pkgs, library, character.only = TRUE)

# Load in functions for the analysis.
source("src/functions.R")
  
###########################################################################################
###   STEP 2: Load in occupancy and detectability layer for each species group          ###
###########################################################################################   
  
# Simulations require raster layers of occupancy and detectability for each species. These raster layers need to be built beforehand and 
# should be loaded into R raster stacks. Each species can be detected with up to 3 detection methods. Any combination of 3 methods is allowed.
# For example, species 1 might be detected with 2 methods. These two methods could be method 1 and 2, method 1 and 3 or method 2 and 3.
# Importantly, if you're evaulating power for only one or two methods, you still must create a raster stack for the third unused method.
# All raster stacks must have the same dimensions and cell resolution. The package will check all raster layers and provide warning if there is a mismatch. 
# Cells not considered for monitoriong should be assigned an NA value

group = group
occ<- stack(sprintf("spotR/layers/occ_%s.tif", group)) # Stack species within a group. If loading in raw HDMs make sure to divide by 100 (want a value range of 0 - 1).
occ<-occ/100 
det.method1<- stack(paste(sprintf("spotR/layers/det.method1_%s.tif", group))) # Load in species detectability layers for method 1 combined into a raster stack. Note, if you only assess one species this should still be a stack.
det.method2<- stack(paste(sprintf("spotR/layers/det.method2_%s.tif", group))) 
det.method3<- stack(paste(sprintf("spotR/layers/det.method3_%s.tif", group))) 
det.method4<- stack(paste(sprintf("spotR/layers/det.method4_%s.tif", group))) 

#########################################################################################
###   STEP 3: Load files specifying which methods apply to each species and a         ###
###           statistical model relating occupancy and detectability to covariates    ###
#########################################################################################

# Simulations require a list that specifies which detection methods are relevant to each species. 
# A '1' in the species.list means that detection method is relevant to a species, a '0' means that it is not relevant. 
# If you type 'species.list', you will see what is required. The number of rows depends on the number of species in your raster stacks. 
# The occupancy and detectability raster stacks loaded above can remain static over time or be updated in response to fire events at monitoring sites.
# The species.cov file specifies the statistical model between occupancy/detectability and environmental covariates.
# Note, if you only want to model a static landscape (i.e. no fire), you still need to provide a species.cov file - just fill it with zero's.

if (method.scenario == 1) {species.detection<- read.csv("spotR/inputs/species.det.csv", stringsAsFactors = FALSE)} else{ # Load in species detecion methods.
  species.detection<- read.csv("spotR/inputs/species.det_2.csv", stringsAsFactors = FALSE)} 

species.detection.all<- species.detection
species.detection<- subset(species.detection, species.detection[,1]==group) # Subset by target group.
species.detection<-species.detection[,-1] # Remote group col.
species.list<- species.detection
species.list[,2:5]<- data.frame(lapply(species.list[2:5], function(x) as.integer(x!=""))) # Binary method used or not.

# Load in covariates of detection.
species.cov<- read.csv("spotR/inputs/species.covar.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species. This wont be needed for the fire modelling.
names<- colnames(species.cov)
covariats<- as.data.frame(matrix(0, nrow=length(species.detection[,2]), ncol=length(names)-1))
species.cov<- cbind(species.detection[,2], covariats)
colnames(species.cov)<- names
species.cov[,1]<- as.character(species.cov[,1])

#########################################################################################
###   STEP 4: Load in additional raster layers for simulations                        ###
#########################################################################################

#Three additional raster layers are required for simulations - a parks layer, a remote layer and a vegetation layer
#The parks layer is needed if you want to estimate power within smaller level management units (hereafter referred to as parks).
#Each park should be defined by an integer value starting at 1. For instance, if there are 2 parks, all cells in the first park should be given a value of 1, and all cells in the second park should be given a value of 2 and so forth
#Type 'parks' before laoding in your own layer to see what this raster should look like
#If you only want to estimate power at a landscape level, you'll still have to load in a parks layer. In this case, just load in a raster layer with all cells equal to 1
#A remote layer is needed to divide monitoring sites between different regions of the landscape
#The remote layer should contain 1s and 0s only, with 1s defining remote areas, 0s non-remote areas (i.e burnt and un-burnt).
#You might not want to divide your monitoring sites between remote and non-remote areas. In this case, just make all your cells in the remote layer equal to 1
#In the example data, the veg layer is needed to simulate fire at monitoring sites. Type 'veg' to see what it looks like
#If you don't intend on modelling fire, the program still requires a veg layer. Once again, just make all values equal to 1.

parks<- raster("spotR/layers/burnt.tif") # Load in layer of burnt and unburnt areas.
zonationOutput<- raster(sprintf("zonation/Zout/Z_%s_5k.CAZ_DEA.rank.compressed.tif", group)) # Load in Zonation output.
zonationOutput[zonationOutput < 0.9]<- NA #Mask out bottom 90% of Zonation layer.
remote<- raster("spotR/layers/remote.tif") # Load in remote layer - this is a map of the edge of burnt areas. We might load in the top 10% of the landscape here determined with Zonation.
vic.road<- raster("spotR/layers/vic.road.tif") # Load in time to destination layer for cost function.
veg<- remote # Fire severity classes for site stratification.
init.occ.burnt <- 0.01 # Set initial occupancy in burnt cells.


#########################################################################################
###   STEP 4: (optional) Decide whether to model fire and load in fire history layer  ###
#########################################################################################

# If you model fire at monitoring sites, you have to load in two additional raster stacks - one stack of the environmental covariates in the study region and one stack of the fire history for each cell
# The stack covar is a raster stack of environmental covariates. It is used to update occupancy and detectability layers for fire pronce species given the statistical models 
# In the example dataset, the layers correspond to the following covariates - #1=dem, 2=creek, 3-moist, 4=Tfire, 5=Ffire, 6=roughness, 7=temp
# The fire stack defines the fire history over a pre-determined time period, with 1s represting fire events, 0s representing no fires in a cell for that year
# The example fire raster stack has fire history for the proceding 15 years. This is used to simulate further fire events for each new year of simulations

 model.fire<- FALSE # Is fire modelled at sites. Set to TRUE to model fire, FALSE otherise.
if (model.fire == TRUE) {
   covar<- stack("~/Kakadu/Paper2/Mapsv2/layers2.tif") # Load in covariate raster stack.
   fire<- stack("~/Kakadu/Paper2/Mapsv2/fire.tif") # Load in fire history raster stack.
 }

###########################################################################################
###   STEP 5: Define suite of input parameter for the select sites and cost functions   ###
###########################################################################################

n.species<- nlayers(occ) # Calculates the number of species.
n.park<- unique(parks) # Calculates the number of parks.

###########################################################################################
###   STEP 6: Define where and how to select monitoring sites                           ###
###########################################################################################

# A three monitoring scenarios can be evaulated. Power can be assessed through:
# 1) Randomly paired sites  - sites are randomly selected across the landscape, but paired across burnt and unburnt zones.
# 2) Zonation paired sites - sites are randomly selected across the top 10% of the occupancy cells defined via the Zonation input layer. Sites are paired across burnt and unburnt areas.
# 3) Pre-selected sites - sites are pre-defined via a shape file.

existing.sites<- FALSE # If you want to simulate monitoring at pre-selected sites, set to TRUE. Otherwise, set to FALSE.
sites<- shapefile("spotR/inputs/exisitingSitesShapefile/exisitingSitesSPDF.shp") # If you specifying TRUE above, load in a shapefile of existing sites (as a spatial points dataframe)
sites<- unique(sites@coords) # Extract site co-ordinates.
include.site<- extract(zonationOutput, sites[,1:2]) # Crop to sites within stuyd area.
sites<- sites[!is.na(include.site),] # Remove NA sites.

all.existing.sites<- all.existing.sites #.If you want to monitor via '(2) Randomly paired sites' set to TRUE. If FALSE sites will be selected in random pairs (option 1).

if (existing.sites == TRUE & all.existing.sites == TRUE) {n.sites<- nrow(sites)} else { #If you're monitoring all pre-selected sites, the number of sites is equal to the number of rows in the XY coordinate list.
  n.sites<- 100} #If you're monitoring a subset of the pre-selected sites, enter the MAXIMUM number of sites you wish to survey. Note this can be greater than the number of pre-selected sites. 

init.sites <- n.sites # Set global site max parameter.
R<- site.ratio # This defines the ratio of burnt to unburnt sites. For example, R = 0.4 means that 40% of sites will be in burnt areas and 60% will be in non-burnt areas. If you don't have burnt or non-burnt areas, make all cells in your remote layer equal to 1 and make R = 1. 
new.site.selection<- "pairs" # If you don't pre-select sites, specify how sites will be selected. Currently pairs is the only option.

###########################################################################################
###   STEP 7: Define when to monitor                                                    ###
###########################################################################################

# Users can specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to each site for each detection method.

Tmax<- 10 # Define length of monitoring program in years.
if (survey.year == 1){s.years<- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)} else {s.years<- c(1, 3, 5, 7, 10) } # Reduce the years in which bats surveys occur as they are expensive.
if (repeat.scenario == 1){n.method<- read.csv("spotR/inputs/repeatVisits.csv")} else {n.method<- read.csv("spotR/inputs/repeatVisits_2.csv")} # Set the scenario of repeat visits to a site. The second option has a reduce number of site repeats.
colnames(n.method)<- c("Group", "Method1", "Method2", "Method3", "Method4")
n.method<- as.numeric(n.method[n.method$Group==group, 2:5]) # Set the number of repeats for a given species group.

###########################################################################################
###   STEP 7: Define the power analysis parameters                                      ###
###########################################################################################

# Users must now decide whether calculate power within smaller level management units, how species are expected to decline over time, the Type I error rate,
# whether to conduct a one or two tailed test, the magnitude of the reduction in occupancy over the monitoring program and the number of simulations.

park.level<- FALSE # Set to TRUE to estimate power at a park level, FALSE to estimate power just at a landscape-level.
decline<- "random" # Random is the only option - it means we simulate a constant decline in occupancy across space.
trend<- "increasing" # Set to 'decreasing' for a decreasing trend in occupancy and 'increasing' for an increasing trend.
model.variation<- FALSE # Set if the model takes in variation.
variation<- c(0.273, 0.376, 0.523, 0.073, 0.154, 0.035, 0.523, 0.229, 0.073, 0.229, 0.035, 0.376) # Set model variation 
alpha<- 0.1 # Set significance level. Choose from 0.01, 0.05 or 0.1.
two.tailed<- FALSE # Decide on a one-tailed or two-tailed signficance test. Set to TRUE for a two tailed test, FALSE for a one-tailed test.
effect.size<- c(0.1, 0.3, 0.5, 0.7, 0.9) # Decide on the proportional reduction in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes.
nsims<- 100

###########################################################################################
###   STEP 7: Set up monitoring sites                                                   ###
###########################################################################################

# This section calls the select.sites function to return monitoring sites for the first simulation given the information provided above.
# Occupancy and detectability values are then extracted from the sites for further processing. If all pre-selected sites are to be surveyed this function is not called again.

plotter<- FALSE # Set to TRUE to plot monitoring sites at the start of each simuation.
xy.sites<- select.sites(sites, n.sites, R, all.existing.sites, existing.sites, new.site.selection, plotter, occ, zonationOutput) # Define monitoring sites.
xy.sites2<- SpatialPoints(xy.sites[,1:2], proj4string = CRS("+proj=lcc +lat_1=-36 +lat_2=-38 +lat_0=-37 +lon_0=145 +x_0=2500000 +y_0=2500000 +ellps=GRS80 +units=m +no_defs ")) # Turn to spatial to allocated site 'zone'.
park.ID<- parks[cellFromXY(parks, xy.sites2)] # Extract park values at monitoring sites.
park.ID[is.na(park.ID)]<- 3
xy.sites<- cbind(xy.sites[,1:3], park.ID) # Combined XY coordinates of sites with park values.

###########################################################################################
###   STEP 8: Compute cost to set up sites, and remove sites based on budget            ###
###########################################################################################

# This section calls upon cost.to.survey function to return the cost to monitor a site given the species present, travel time and the detection methods used to find each species (taking into account survey duration and repeats),
# Up to four detection methods are assigned from species.detection (Step 3), with costs defined via detection.costs.
# Site locations are taken from xy.sites (Step 10), and travel time to each site is extracted from vic.road within the cost.to.survey.site function.

overhead.cost<- FALSE # Set whether to include set up costs in budgeting.

# Define proportion of monitoring budget for each species group Plants get 40% of budget (correspondance with DELWP).
max.budget<- (max.budget* (n.species/45))

# Load in method cost table.
method.cost.array<- read.csv("spotR/inputs/method.cost.array.csv")

# Run cost function, outputs cost table header and plot of site locations.
project.costs <- cost.to.survey.site(zonationOutput, occ, sites, xy.sites, species.detection, species.detection.all, existing.sites,  n.sites, s.years, n.method, group, vic.road, method.cost.array, overhead.cost, max.budget, plotter)

# Remove excess sites over budget.
program.cost <- remove.excess.site(new.site.selection, project.costs, max.budget, all.existing.sites, det.method1.time, det.method2.time, det.method3.time, det.method4.time, occ.time, occ, plotter, xy.sites)
program.cost <- sum(program.cost$program.cost) # Total cost of monitoring program 

###########################################################################################
###   STEP 9: Define which sites are management, and the benefit of management             ###
###########################################################################################
man.effect.size <- 0
prop.man <- 0.5
managed <- rep(0, nrow(xy.sites))
managed[sample(nrow(xy.sites), round(nrow(xy.sites)*prop.man))] <- 1
xy.sites <- cbind(xy.sites, managed)

###########################################################################################
###   STEP 10:  Perform a pre-simulation check of layers and inputs                     ###
###########################################################################################

# Check all raster layers have the same dimensions and then run the power analysis. 
# A warning message will be displayed if there is a mismatch between rasters.
# This types of analysis can take a long time on a desktop. The time required to run a simulation will increase with the number of species, the number of cells in each raster, the number of sites, whether siets are randomly selected each simulation, whether fire is modelled, whether power is estimated at a park level etc
# Simulations can be sped up by running them in parallel (or via computing clusters). Excessively large datasets should be run on a computer with at least 10 cores .  

check.inputs(occ) # This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK .

###########################################################################################
###   STEP 11:  Define inputs for analysis on the HPC computing resource                ###
###########################################################################################

n.cores<- 5  # Define the number of cores per node. This needs to reflects the number in 'job_submission.slurm', with one additional core to run background tasks. 
start.time<- proc.time() # Record start time for run comparison.
registerDoFuture() # Register the cores via the 'Future' and 'doFuture' R package.
plan(multiprocess, workers = n.cores) # This function detects the OSX and computer specification and automatically set's up the inputs for parralel computing. Set to plan(sequential) to run in sequence locally.
###########################################################################################
###   STEP 12:  Run the simulation and save out objects for plotting and analysis       ###
###########################################################################################

# Run the power analysis. Set %dopar% to run in parallel, %do% otherwise. Note, you cannot track progress in parallel. If you run in sequence first, you can see how long it takes to do each step of the simulation.
pwr<- foreach(i = 1:n.cores,.combine = '+',.packages = c("raster", "rgdal")) %dopar% {
  
  set.seed(1234 + i) # Set unique seed for each node of the analysis.
  
  .libPaths("/home/asmart1/R/lib/3.6") # Point to library on HPC cluster.
  library("unmarked") # Load in node specific dependacies.
  library("rgdal")
  
  sapply(effect.size, run.power, 
         nsims = nsims, 
         alpha = alpha, 
         decline = decline, 
         Tmax = Tmax, 
         s.years = s.years, 
         trend = trend,
         sites = sites, 
         model.fire = model.fire, 
         species.list = species.list, 
         park.level = park.level, 
         xy.sites = xy.sites, 
         R = R,
         two.tailed = two.tailed, 
         plotter = plotter, 
         n.sites = n.sites, 
         n.species = n.species, 
         n.park = n.park, 
         all.existing.sites = all.existing.sites,
         existing.sites = existing.sites, 
         new.site.selection = new.site.selection,
         n.method = n.method,
         occ = occ,
         det.method1 = det.method1,
         det.method2 = det.method2,
         det.method3 = det.method3,
         det.method4 = det.method4,
         veg = veg,
         remote = remote,
         parks = parks,
         variation = variation,
         model.variation = model.variation,
         zonationOutput = zonationOutput,
         species.detection = species.detection,
         species.detection.all = species.detection.all,
         group = group, 
         vic.road = vic.road, 
         method.cost.array = method.cost.array, 
         overhead.cost = overhead.cost, 
         max.budget = max.budget,
         project.costs = project.costs,
         init.sites = init.sites,
         prop.man = prop.man,
         man.effect.size = man.effect.size,
         init.occ.burnt = init.occ.burnt,
         program.cost = program.cost)
}

time.elapsed<- proc.time() - start.time # Record end time for run comparison.
time.elapsed #  Print elapsed time.
plan(sequential) # Close nodes post simulation.

############################################################################################
###   STEP 13:  Save outputs and plot results                                            ###
############################################################################################

Results<- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores)

# Produce table of power values per species.
pwrtble<-matrix(NA, nrow = 5, ncol = nrow(species.list))
for (i in 1:5){
  pwrtble[i,]<- Results[[1]][,,i][1,] 
}
colnames(pwrtble)<- species.list$Species
row.names(pwrtble)<- c("0.1", "0.3", "0.5", "0.7", "0.9")
powerTable<- t(pwrtble)

# Write table to local.
write.table(powerTable, file = sprintf("out/tabs/powerTable_%04d_group=%s_zonation=%s_site.ratio=%s_repeat.scenario=%s_method.scenario=%s_survey.year=%s.txt", parameter_index, group, all.existing.sites, site.ratio, repeat.scenario, method.scenario, survey.year), sep = ",", row.names = T, col.names = T )

# Save out model object.
save(pwr, file =  sprintf("out/output_%04d_group=%s_zonation=%s_site.ratio=%s_repeat.scenario=%s_method.scenario=%s_survey.year=%s.RData", parameter_index, group, all.existing.sites, site.ratio, repeat.scenario, method.scenario, survey.year))

# Save out diagnostic plot.
pdf(file = sprintf("out/figs/power_efs_%04d_group=%s_zonation=%s_site.ratio=%s_repeat.scenario=%s_method.scenario=%s_survey.year=%s.pdf", parameter_index, group, all.existing.sites, site.ratio, repeat.scenario, method.scenario, survey.year),  height = 8.5, width = 8.5)
Results<- plot.results(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores)
dev.off()

