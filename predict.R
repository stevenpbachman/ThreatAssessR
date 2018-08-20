#####################################################################################
### this will get variables for the smaller species i.e. <100,000 occurrences
### TO DO - collect bigkeys and make a table so that these can be run separately via the download function

# Load this first - START ####################################################################

TDWG_to_IUCN_version1_UTF_8 <- read.delim("TDWG_to_IUCN_version3_UTF-8.txt", encoding="UTF-8")
raster.tdwg = raster("rasters/tdwg3.tiff")
raster.eco = raster("rasters/eco2017.tif")
tdwg_raster = read_delim("tdwg_raster.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
IUCN_taxonomy = read.delim("IUCN_taxonomy.csv", encoding="UTF-8", sep = ",")
TDWG_realms = read.delim("TDWG_realms.csv", encoding="UTF-8", sep = ",")

# set path for results - default setting is outputs folder in working directory
path = paste0(getwd(),"/outputs/")

# Load this first - END ####################################################################

#########################################################################################################
# 1. ok, now load in a list of species to assess
# the list should have a column with a binomial and another column with IPNI ID

# Read in test species as an example
example_species = read.csv("outputs/example_species.csv")

redlist2018 = read.csv("C:/Users/sb42kg/OneDrive - The Royal Botanic Gardens, Kew/02_Publications/spb/Machine learning - big data red list/05_R/03_Model_Predict/RedLeast_RedList_join_run_on_2018-08-17.csv")



# get the family from Plants of the world online (POWO)
#family = "Arecaceae"
#powo.species = get.POWO(family)

# define full_name (binomial) and ID_list (IPNI ID)
full_name = example_species[,4] # work through this sequence and save, might need to split into batches of 500 - 1000
ID_list = example_species[,7]

# user details for IUCN point file and SIS csv outputs. These are required fields, but will be the same for all assessments
Firstname = 'Steve'
Lastname = 'Bachman'
Email = 's.bachman@kew.org'
Initials = 'S.P.'
Affiliation = 'Royal Botanic Gardens, Kew'
# point compiler needs to be specific format: last name, initial and affiliation
pointcompiler = paste0(Lastname,", ",substr(Firstname,1,1),". (", Affiliation,")")
pointcitation = Affiliation  
credittype = 'Assessor'

# Run ThreatAssessR
t = proc.time()
LC = ThreatAssess(full_name = full_name, ID_list = ID_list, LC.points = FALSE, SIS.files = FALSE, path = path)
proc.time()- t





