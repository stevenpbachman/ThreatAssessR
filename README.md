# ThreatAssessR
Tool to predict threat status of plants and generate full Red List assessments for species predicted as Least Concern.

# How to run:

Open 01b_FUNCTIONS_Get_RedList_variables and source on save. THis loads all the relevant libraries and functions

# Still to do:
To do:
* investigate outlier cleaning with zizka location tools in R
*Report GBIF download DOI – not sure how to get this if you only use occ_search()
*Check the dataset ID and use that to get the DOI e.g. iNaturalist.org (2018). iNaturalist Research-grade Observations. Occurrence dataset https://doi.org/10.15468/ab3s5x accessed via GBIF.org on 2018-08-01.
library(rgbif)

`#get the dataset`
`test = datasets(uuid = '7bac0ff4-f762-11e1-a439-00145eb45e9a')`

`# get the DOI`
`test$data$doi`

`# use this to build the references for the species assessment`
`#Dep. Biology, Univ. Autónoma de Madrid (2018). Universidad Autónoma de Madrid, Biología, Acalypha. Occurrence dataset`
`#https://doi.org/10.15468/yhsqkx accessed via GBIF.org on 2018-08-01.`

*	Check FAM variable for modelling
	Either add new family code and see what happens
	Need mapping file that converts family to integer – you need to convert back after the modelling stage
	How many families are missing? Get from POWO? 
*	Try crosstalk for HTML widget
*	Add unique ID to gbif points data?
*	Sort out problem with small species numbers – returns error
*	Add estimate for how long it will take based on number of collections – need to plot this on timing test (sample at 10, 100, 1000, 10,000 ?). Random list from POWO – or Red List dataset?
*	Add message to script to estimate how long it will take to get the data, based on the size of the points per species e.g. This is a large list, or a list with of species with large numbers of points. It may take XX time to complete.
*	Add HTML report options to output
* Leaflet for interactive maps, link points?
*	Use SIS csv files to generate report – you need to join on species name ‘i’ from the loop
*	Also add data from CSV files to the repot
*	Get SIS data first independent of whether user asks for it so that it can be used for the report. Then, later in script, if user asks for SIS, then save the SIS files already generated
*	Add lime local model output to explain predicted result
*	As above, but in PDF format!
*	Clean up old code no longer needed e.g. functions– save somewhere else
*	Make code elegant! How – check in with Baz and Justin. Tidyverse – piping etc.
*	Add options to query by species or by region, or by csv upload or by family? Genus?
*	This is mostly done, but needs to be clearly written up. Markdown?
*	Check countries – is this working and reporting correctly – check commas, and make sure that original TDWG code is returned (not just the error message)
*	Think about Ebbe Nielsen prize – 5th September deadline
*	Try to model habitat
*	C:\Users\sb42kg\OneDrive - The Royal Botanic Gardens, Kew\02_Publications\spb\Machine learning - big data red list\01_Data\Ecoregion vs IUCN classification
*	Add plant growth form – POWO
*	Comment throughout code – functions
*	Trial with large datasets – palms?
*	Share with Justin and Baz, ask for help turning into a package
*	Add link to relevant specialist groups – add contact details – email address for review – or save this for when it is in SIS?
*	Species priority list per country? Or TDWG region

Already fixed:

*~~	Add model to assign rating – with probability (Not threatened (LC), Possibly threatened (CR< EN< VU< NT)~~
*~~	Check merge (join) – could be the wring kind of join causing multiplication and memory loss~~
*~~	get version on github~~
*~~	Investigate raster memory issues:~~
*~~	https://cran.r-project.org/web/packages/raster/vignettes/functions.pdf~~
*~~	Remove LC labels with thresholds, but keep code.~~


