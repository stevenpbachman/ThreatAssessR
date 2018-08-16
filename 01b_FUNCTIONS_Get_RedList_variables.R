#functions.R
#Contains all of the functions needed to perform the actual analysis. 
#source()'ing this file should have no side effects other than loading up the 
#function definitions. This means that you can modify this file and reload it 
#without having to go back an repeat steps 1 & 2 which can take a long time to run for large data sets.

### RedLeast - gather parameters for plant species using GBIF and POWO (distribution, EOO, AOO etc.)
### Compare results with predefined model to determine if species are Least Concern according to IUCN criteria 
### For LC species you can generate IUCN SIS point files and SIS csv files for Red List submission of LC species
### Steve Bachman & Justin Moat & Robin Freeman

# other options 
options(stringsAsFactors = FALSE)

# install libraries for Red List download
#install.packages("rredlist")

# load libraries for Red List download
library(rredlist)

# 0. RED LIST
# function to get the top level species info based on each ID - WORKING ###
rl.apply = function(x, key){
  sp <- rl_search(id=x,key=rlkey)
  sp = sp$result
  return(sp)
}

# get IDs for all plants on the Red List ###
rl_plant_IDs = function(rlkey){
  
  # get IDs (all taxa on Red List) from pages of 10,000 and keep going until page has no more species  
  i=1
  testcount = 1
  rl = rl_sp(0, key = rlkey)
  rl.all = rl$result
  
  while ((testcount)>'0'){
    test = rl_sp(i, key = rlkey)
    testcount = test$count
    rl.all = rbind(rl.all,test$result)
    i = i+1
  }
  # filter on plants 'PLANTAE'
  plants <- rl.all[ which(rl.all$kingdom_name=='PLANTAE'),] 
  
  # filter out the first column 'taxon ID'
  plants.id = plants[,1]
  
  # run through all IDs using lapply
  apply.rl = lapply(plants.id,rl.apply)
  rl_df = as.data.frame(do.call(rbind,apply.rl))
  
  #get RL version number
  #version = rl_version(key = rlkey)
  
  #get date
  #today = Sys.Date()
  
  #file path
  #filename = paste0("PLANTAE_RL_version_",version,"_run_on_",today)

  return(rl_df)
  
}

# install libraries for RedLeast 
#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("rgbif")
#install.packages("rgdal")
#install.packages("sp")
#install.packages("rCAT")
#install.packages("data.table")
#install.packages("raster")
#install.packages("jsonlite")

# load relevant libraries for Red Least
library(plyr)
library(tidyverse)
library(rgbif)
library(rgdal)
library(sp)
library(rCAT)
library(data.table)
library(raster)
library(jsonlite)
library(devtools)
library(pryr)

############### FUNCTIONS ##############

# get species name from gbif by key
add_names = function(BigKeys){
  
  tax_name = name_usage(key=BigKeys)
  full_name = tax_name$data$species
}

# error catcher - updates results.table when there are errors
update.errors = function(result.table, error.species, error.text) {
  if (nrow(error.species) < 1)
  {
    result.table = result.table
  }
  else{
    as.data.frame(error.species)
    error.species$eWarning = error.text
    result = full_join(result.table,
                       error.species,
                       by = "POWO_ID",
                       copy = FALSE)
    tidytable = result %>% unite(Warning, Warning, eWarning, sep = ";")
    tidytable$Warning = gsub('NA', "", tidytable$Warning)
    tidytable$Warning = gsub(';;', "", tidytable$Warning)
    return(tidytable)
  }
}

# 00  - get names from POWO API
# 0.1 - get species names and ID - user supplies family. POWO may fix page size issue later
get.POWO = function(family, pagesize = 24) {
  #get list first
  full_url =  paste(
    "http://powo.science.kew.org/api/1/search?page=0&page.size=",
    pagesize,
    "&f=accepted_names%2Cspecies_f&q=family%3A",
    family,
    sep = ""
  )
  raw.data <- readLines(full_url, warn = "F", encoding = "UTF-8")
  rd <- fromJSON(raw.data)
  pages = rd$totalPages
  last = pages - 1
  list = c(1:last)
  all.results = rd$results
  accepted = all.results$accepted
  author = all.results$author
  name = all.results$name
  url = all.results$url
  results = cbind(accepted, author, name, url, family = family)
  
  for (l in list) {
    spe_url =  paste(
      "http://powo.science.kew.org/api/1/search?page=",
      l,
      "&page.size=",
      pagesize,
      "&f=accepted_names%2Cspecies_f&q=family%3A",
      family,
      sep = ""
    )
    raw.data <- readLines(spe_url, warn = "F", encoding = "UTF-8")
    rd <- fromJSON(raw.data)
    all.results = rd$results
    accepted = all.results$accepted
    author = all.results$author
    name = all.results$name
    url = all.results$url
    all.results = cbind(accepted, author, name, url, family = family)
    results = rbind(results, all.results)
  }
  fam.results = as.data.frame(results)
  fam.results = fam.results %>% unite(full_name, name, author, sep = " ")
  fam.results = fam.results %>% separate(url, into = c("base_url", "POWO_ID"), sep = 31)
  
  return(fam.results)
}

# 0.2 - get names by checkin binomial (note - doesn;t check authors)
# check if name (binomial) is accepted on POWO. ID is the binomial as a list
check.accepted.POWO = function(binomial) {
  
  #binomial = "Poa annua"
  
  # use ID full name to search API  
  full_url =  paste("http://plantsoftheworldonline.org/api/1/search?q=species:", binomial, sep = "")
  
  # get raw json data
  raw.data <- readLines(full_url, warn = "F", encoding = "UTF-8")
  
  # organise
  rd <- fromJSON(raw.data)
  
  if (length(rd$results) == 0) {
    # add new column to results and call it warning
    
    accepted = "NA"          
    author = ""  
    kingdom = ""  
    name = ""  
    rank  = "" 
    #synonymOf = ""
    base_url = ""    
    IPNI_ID = ""
    search_name = binomial  
    results = data.frame(accepted,author, kingdom, name, rank, base_url, IPNI_ID,search_name)
    
  } else {
    
    # make data frame
    results = as.data.frame(rd$results)
    
    # add original search term
    results$search_name = binomial
    
    # filter on 
    if (nrow(results) >1) {
      
      accepted = "FALSE"          
      author = ""  
      kingdom = ""  
      name = ""  
      rank  = "" 
      #synonymOf = ""
      base_url = ""    
      IPNI_ID = ""
      search_name = binomial  
      results = data.frame(accepted,author, kingdom, name, rank, base_url, IPNI_ID,search_name)
      
    } else {
      
      # PROBLEM HERE - var with no author - replace with blank field?
      if (!"author" %in% colnames(results)) {
        results$author = NA
      }
      
      # only include these fields - you don't want synonym of
      results = subset(results, select=c(accepted, author, kingdom,name,rank, url, search_name))
      
      # split url to get url and IPNI_ID
      results = results %>% separate(url, into = c("base_url", "IPNI_ID"), sep = "names:")
      
      # take out any results where it matched on something that wasn't species 
      results = subset(results, rank == "Species") 
      
      # check if 
      if (nrow(results) < 1) {
        # add new column to results and call it warning
        
        accepted = "NA"          
        author = ""  
        kingdom = ""  
        name = ""  
        rank  = "" 
        #synonymOf = ""
        base_url = ""    
        IPNI_ID = ""
        search_name = binomial  
        results = data.frame(accepted,author, kingdom, name, rank, base_url, IPNI_ID,search_name)
        
      }
      # make data frame
      results = as.data.frame(results)
    }
    
    return(results)
    
  }
}

# 1 - get gbif name key using name_suggest and pick 1st suggestion if >1 or add blank row if no match, then get GBIF taxonomy
# 1.1
LC.gbif.key = function (full_name) {
  gbif.key = name_backbone(name = full_name, rank = 'species',kingdom = 'Plantae', strict = TRUE, verbose = FALSE) ### NOTE: try name_backbone instead, so we can make filter on kingdom = PLANTAE
  
  # add blank row when there is no value returned
  if(gbif.key$matchType=="NONE"){
    key = "NA"
    canonicalName = "NA"
    gbif.key = data.frame(key,canonicalName, stringsAsFactors = FALSE)
  }  else {
  # now reformat output of 
    key = gbif.key$usageKey
    canonicalName = gbif.key$canonicalName
    gbif.key = data.frame(key, canonicalName,stringsAsFactors = FALSE)
  }
  # modify field headings to keep clean
  colnames(gbif.key)[which(names(gbif.key) == "key")] = "GBIF_SuggestedKey"
  colnames(gbif.key)[which(names(gbif.key) == "canonicalName")] = "GBIF_SuggestedName"
  gbif.key = gbif.key[, 1:2]
  gbif.key = data.frame(gbif.key, stringsAsFactors = FALSE)
  return(gbif.key)
}
# 1.2 - results.keys
LC.get.keys = function(full_name) {
  apply.gbif.key = lapply(full_name, LC.gbif.key)
  keys = do.call(rbind, apply.gbif.key)
  #keys.results = cbind(keys,LC.table[1],ID_list[2]) %>% as.data.frame()
  # update
  return(keys)
}
# 1.3 - get higher taxonomy for taxonomy.csv
LC.gbif.tax = function (result.table.keys) {
  if ((result.table.keys) == "NA") {
    # add new column to results and call it warning
    kin = ""
    phy = ""
    ord = ""
    fam = ""
    GBIF_SuggestedGen = ""
    GBIF_SuggestedSpe = ""
    GBIF_Suggestedauth = ""
    GBIF_SuggestedNameStatus = ""
    GBIF_SuggestedKey = ""
    GBIF_AcceptedKey = ""
    Warning = "No name match in GBIF"
    taxlist = data.frame(
      GBIF_SuggestedKey,
      kin,
      phy,
      ord,
      fam,
      GBIF_SuggestedGen,
      GBIF_SuggestedSpe,
      GBIF_Suggestedauth,
      GBIF_SuggestedNameStatus,
      GBIF_AcceptedKey,
      Warning
    )
  } else {
    # filter out NA
    # filter on key only - suggested key
    gbif_nam_search = name_usage(result.table.keys, rank = 'family')
    kin = gbif_nam_search$data$kingdom
    phy = gbif_nam_search$data$phylum
    ord = gbif_nam_search$data$order
    fam = gbif_nam_search$data$family
    GBIF_SuggestedGen = gbif_nam_search$data$genus
    GBIF_SuggestedSpe = gbif_nam_search$data$species
    GBIF_Suggestedauth = gbif_nam_search$data$authorship
    GBIF_SuggestedNameStatus = gbif_nam_search$data$taxonomicStatus
    GBIF_SuggestedKey = as.character(gbif_nam_search$data$key)
    GBIF_AcceptedKey = as.character(gbif_nam_search$data$speciesKey)
    Warning = ""
    taxlist = data.frame(
      GBIF_SuggestedKey,
      kin,
      phy,
      ord,
      fam,
      GBIF_SuggestedGen,
      GBIF_SuggestedSpe,
      GBIF_Suggestedauth,
      GBIF_SuggestedNameStatus,
      GBIF_AcceptedKey,
      Warning
    ) #gbifkey
  }
  search.df = c(
    "PROPARTE_SYNONYM",
    "DOUBTFUL",
    "HETEROTYPIC_SYNONYM",
    "HOMOTYPIC_SYNONYM",
    "MISAPPLIED",
    "SYNONYM"
  )
  taxlist$Warning[taxlist$GBIF_SuggestedNameStatus %in% search.df] = "GBIF matched name not treated as accepted"
  return(taxlist = as.data.frame(taxlist))
}
# 1.4 - combine functions into one query to get keys and tax
LC.get.keys.tax = function(full_name, ID_list) {
  result.keys = LC.get.keys(full_name) # 1 - lapply on gbif key and merge back to get results table which will be appended to later
  result.table.keys = cbind(ID_list, full_name, result.keys)
  result.tax = lapply(result.table.keys[, 3], LC.gbif.tax)
  result.tax = do.call(bind_rows, result.tax)
  result.table.tax = cbind(ID_list, full_name, result.tax)
  colnames(result.table.tax)[which(names(result.table.tax) == "ID_list")] = "POWO_ID"
  return(result.table.tax)
}

# 2 - get tdwg areas # and convert to country csv
# 2.0
LC.tdwg = function(ID) {
  #ID = "665137-1"
  full_url =  paste("http://plantsoftheworldonline.org/api/1/taxon/urn:lsid:ipni.org:names:", ID, sep = "")
  raw.data <- readLines(full_url, warn = "F", encoding = "UTF-8")
  rd <- fromJSON(raw.data)
  if (length(rd$distributions) == 0) {
    nat_dis = data.frame(
      tdwgCode = "NA",
      featureId = "NA",
      tdwgLevel = "NA",
      establishment = "NA",
      LEVEL3_NAM = "NA",
      POWO_ID = ID,
      TDWG_Name = "NA",
      ID = "",
      countryoccurrencename = "NA",
      countryoccurrencelookup = "NA",
      area = "NA",
      forest_loss = "NA",
      hfoot = "NA",
      w_forest_loss = "NA",
      w_hfoot = "NA")
  }
  
  else{
    nat_dis = rd$distributions
    nat_dis["POWO_ID"] = ID
    colnames(nat_dis)[which(names(nat_dis) == "name")] = "LEVEL3_NAM"
    colnames(TDWG_to_IUCN_version1_UTF_8)[which(names(TDWG_to_IUCN_version1_UTF_8) == "Level.3.code")] =
      "tdwgCode"
    nat_dis = merge(nat_dis, TDWG_to_IUCN_version1_UTF_8, by = "tdwgCode")
    
    #replace Panamá with Panama
    nat_dis$LEVEL3_NAM = gsub("Panamá", "Panama", nat_dis$LEVEL3_NAM)
    return(nat_dis)
  }
}

LC.countries = function(unique.sp, result.tdwg) {
  #unique.sp = unique.sp[1]
  result.tdwg.reduced = subset(result.tdwg, subset = result.tdwg$POWO_ID ==
                                 unique.sp)
  colnames(TDWG_to_IUCN_version1_UTF_8)[which(names(TDWG_to_IUCN_version1_UTF_8) == "Level.3.code")] =
    "tdwgCode"
  country.merge = merge(result.tdwg.reduced, TDWG_to_IUCN_version1_UTF_8, by = "tdwgCode")
  country.merge = distinct(country.merge, countryoccurrencename, .keep_all = TRUE)
  powoid = country.merge$POWO_ID
  countries = data.frame(
    internal_taxon_id = powoid,
    CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceLookup = country.merge$countryoccurrencelookup,
    CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName = country.merge$countryoccurrencename,
    CountryOccurrence.CountryOccurrenceSubfield.presence = 'Extant',
    CountryOccurrence.CountryOccurrenceSubfield.origin = country.merge$establishment,
    CountryOccurrence.CountryOccurrenceSubfield.seasonaility = 'Resident'
  )
  return(countries)
}

# 2.1 - collect errors i.e. no TDWG data and update results table
LC.tdwg.apply = function(ID_list) {
  apply.countries = lapply(ID_list, LC.tdwg)
  tdwg = do.call(rbind, apply.countries)
  return(tdwg) # add tdwg, to outputs
}
# 2.2
LC.tdwg.errors = function(result.tdwg) {
  errors = subset(result.tdwg, (LEVEL3_NAM == "NA"))
  errors = errors[, 6, drop = FALSE]
  return(errors)
}

# 3 - use GBIF keys to get points

# 3.0 - count occurrences from GBIF key
LC_count_occs = function(gbifkey){
  count = occ_count(taxonKey=gbifkey, georeferenced=TRUE)
  #count = cbind(count,gbifkey)
  return(count)
}
# this function get the data from gbif
download_bigkeys = function(dkeys){
  dkeys = dkeys[1]
  occ_download_get(dkeys, overwrite = TRUE) %>% occ_download_import()
}
# clean points
big_clean = function(big_points){
  
  
  if (!"recordNumber" %in% colnames(big_points)) {
    big_points$recordNumber = NA
    as.character(big_points$recordNumber)
  }
  
  if (!"decimalLongitude" %in% colnames(big_points)) {
    big_points$decimalLongitude = NA
    as.character(big_points$decimalLongitude)
  }
  
  if (!"decimalLatitude" %in% colnames(big_points)) {
    big_points$decimalLatitude = NA
    as.character(big_points$decimalLatitude)
  }
  
  if (!"year" %in% colnames(big_points)) {
    big_points$year = NA
    as.character(big_points$year)
  }
  
  if (!"datasetKey" %in% colnames(big_points)) {
    big_points$datasetKey = NA
  }
  
  if (!"basisOfRecord" %in% colnames(big_points)) {
    big_points$basisOfRecord = NA
  }
  
  if (!"catalogNumber" %in% colnames(big_points)) {
    big_points$catalogNumber = NA
  }
  
  if (!"recordedBy" %in% colnames(big_points)) {
    big_points$recordedBy = NA
  }
  
  if (!"issues" %in% colnames(big_points)) {
    big_points$issues = NA
  }
  
  if (!"institutionCode" %in% colnames(big_points)) {
    big_points$institutionCode = NA
  }
  
  if (!"country" %in% colnames(big_points)) {
    big_points$country = NA
  }
  
  if (!"familyKey" %in% colnames(big_points)) {
    big_points$familyKey = NA
    as.character(big_points$familyKey)
  }
  
  if (!"scientificName" %in% colnames(big_points)) {
    big_points$scientificName = NA
  }
  
  
  big_points_clean = subset(big_points,    select = c(
    'basisOfRecord',
    'datasetKey',
    'taxonKey',
    'familyKey',
    'scientificName',
    'decimalLongitude',
    'decimalLatitude',
    'year',
    'issues',
    'country',
    'recordNumber',
    'catalogNumber',
    'recordedBy',
    'institutionCode'
  ))
}
# this function kicks off the downloads for all of the big keys
# as we can only have three at a time, it waits 10 mins (a little less than max) and then starts on the next in the list
get_bigkeys = function(big_gbifkeys){
  #taxonKey = as.data.frame(big_gbifkeys$gbifkey)
  download_text = paste0("occ_download(","'taxonKey = ",big_gbifkeys,"','hasGeospatialIssue = FALSE','hasCoordinate = TRUE')")
  eval(parse(text = download_text))
  Sys.sleep(900) # change to 10 minutes i.e. 600
}

gbif_queue <- function(...) {
  reqs <- lazyeval::lazy_dots(...)
  results <- list()
  groups <- split(reqs, ceiling(seq_along(reqs)/3))
  
  for (i in seq_along(groups)) {
    cat("running group of three: ", i)
    res <- lapply(groups[[i]], function(w) {
      tmp <- tryCatch(lazyeval::lazy_eval(w), error = function(e) e)
      if (inherits(tmp, "error")) {
        "http request error"
      } else {
        tmp
      }
    })
    
    # filter out errors
    res_noerrors <- Filter(function(x) inherits(x, "occ_download"), res)
    still_running <- TRUE
    while (still_running) {
      metas <- lapply(res_noerrors, occ_download_meta)
      status <- vapply(metas, "[[", "", "status", USE.NAMES = FALSE)
      still_running <- !all(tolower(status) %in% c('succeeded', 'killed'))
      Sys.sleep(2)
    }
    results[[i]] <- res
  }
  
  results <- unlist(results, recursive = FALSE)
  
  return(results)
}

# 3.1 - get points from gbif using clean keys
LC.gbif.points = function(gbifkey, result.tdwg) {
  #gbifkey = result.table
  #gbifkey = result.table[,3]
  
  if (gbifkey == "") {
    basisOfRecord = as.character("NA")
    datasetKey = as.character("NA")
    taxonKey = as.character("NA")
    familyKey = as.integer("-999")
    scientificName = as.character("NA")
    decimalLongitude = as.numeric("-999")
    decimalLatitude = as.numeric("-999")
    year = as.integer("-999")
    issues = as.character("NA")
    country = as.character("NA")
    recordNumber = as.character("NA")
    catalogNumber = as.character("NA")
    recordedBy = as.character("NA")
    institutionCode = as.character("NA")
    
    res = data.frame(
      datasetKey,
      basisOfRecord,
      taxonKey,
      familyKey,
      scientificName,
      decimalLongitude,
      decimalLatitude,
      year,
      issues,
      country,
      recordNumber,
      catalogNumber,
      recordedBy,
      institutionCode,
      stringsAsFactors = FALSE
    )
  } else    {
    #res = occ_search(taxonKey = gbifkey, hasGeospatialIssue = FALSE , hasCoordinate = TRUE, fields=c('taxonKey', 'familyKey', 'scientificName','basisOfRecord','issues', 'decimalLatitude','hasCoordinate', 'decimalLongitude','year','recordNumber','recordedBy','catalogNumber', 'datasetKey', 'country','institutionCode'),limit=200000)
    
    res = occ_data(
      taxonKey = gbifkey,
      hasGeospatialIssue = FALSE,
      hasCoordinate = TRUE,
      #geometry = 'POLYGON((-105.253164768219 25.917780721001716,-79.237539768219 25.917780721001716,-79.237539768219 5.746349439479163,-105.253164768219 5.746349439479163,-105.253164768219 25.917780721001716))',
      limit = 200000
    )
    
    res$data$taxonKey = gbifkey
    res = res$data
    res = as.data.frame(res)
    
    
    if (!"recordNumber" %in% colnames(res)) {
      res$recordNumber = NA
      as.character(res$recordNumber)
    }
    
    if (!"decimalLongitude" %in% colnames(res)) {
      res$decimalLongitude = NA
      as.character(res$decimalLongitude)
    }
    
    if (!"decimalLatitude" %in% colnames(res)) {
      res$decimalLatitude = NA
      as.character(res$decimalLatitude)
    }
    
    if (!"year" %in% colnames(res)) {
      res$year = NA
      as.character(res$year)
    }
    
    if (!"datasetKey" %in% colnames(res)) {
      res$datasetKey = NA
    }
    
    if (!"basisOfRecord" %in% colnames(res)) {
      res$basisOfRecord = NA
    }
    
    if (!"catalogNumber" %in% colnames(res)) {
      res$catalogNumber = NA
    }
    
    if (!"recordedBy" %in% colnames(res)) {
      res$recordedBy = NA
    }
    
    if (!"issues" %in% colnames(res)) {
      res$issues = NA
    }
    
    if (!"institutionCode" %in% colnames(res)) {
      res$institutionCode = NA
    }
    
    if (!"country" %in% colnames(res)) {
      res$country = NA
    }
    
    if (!"familyKey" %in% colnames(res)) {
      res$familyKey = NA
      as.character(res$familyKey)
    }
    
    if (!"scientificName" %in% colnames(res)) {
      res$scientificName = NA
    }
  }
  
  res = subset(
    res,
    select = c(
      'basisOfRecord',
      'datasetKey',
      'taxonKey',
      'familyKey',
      'scientificName',
      'decimalLongitude',
      'decimalLatitude',
      'year',
      'issues',
      'country',
      'recordNumber',
      'catalogNumber',
      'recordedBy',
      'institutionCode'
    )
  )
  return(res)
}
# 3.2 - clean the points so they are ready for TDWG raster
clean.points = function(points) {
  points.clean = points[!is.na(points$decimalLongitude),]
  if (nrow(points.clean) < 1) {
    
  }  else {
    long = points.clean$decimalLongitude
    lat = points.clean$decimalLatitude
    coords = cbind(long, lat)
    coords <- data.matrix(coords)
    sp = SpatialPoints(coords)
    return(sp)
  }
}
# 3.3 - make copy of clean points as data.frame to merge back to later
clean.points.df = function(points) {
  points.clean.df = points[!is.na(points$decimalLongitude),]
  if (nrow(points) < 1) {
    
  }  else {
    return(points.clean.df)
  }
} # check if this needed

# 4 - get TDWG values and clip points so that only native are being used - check native
# 4.1  extract tdwg raster values
LC.extract.tdwg = function(points.clean.df, raster.tdwg) {
  if (nrow(points.clean.df) < 1) {
    
  }  else {
    decimalLongitude = points.clean.df$decimalLongitude
    decimalLatitude = points.clean.df$decimalLatitude
    coords = cbind(decimalLongitude, decimalLatitude)
    coords <- data.matrix(coords)
    coords = unique(coords)
    sp = SpatialPoints(coords)
    
    
    # extract values from raster for each point lat long
    extract.out = raster::extract(raster.tdwg, sp, df = TRUE)
    extract.out = cbind(extract.out, coords)
    
    # change column names for raster out and points so that you can do merge - also subset to get rid of NA
    colnames(extract.out)[which(names(extract.out) == "tdwg3")] = "VALUE"
    merge = merge(extract.out, tdwg_raster, by = "VALUE", all.x = TRUE)
    tdwg.merge.out = subset(merge, LEVEL3_NAM != "NA") # this gets rid of TDWG areas that had no points in
    #tdwg.merge.out = tdwg.merge.out[order(tdwg.merge.out$ID),]
    # final bind back to points
    point.clean.bind = merge(
      points.clean.df,
      tdwg.merge.out[, c(3, 4, 5)],
      by = c("decimalLongitude", "decimalLatitude"),
      all = TRUE
    )
    return(point.clean.bind)
  }
}
# 4.2 merge tdwg.extract with result.tdwg to only get native points in native distributions
native = function(tdwg.extract, powo.list, result.tdwg) {
  if (is.null(tdwg.extract)) {
    
  }  else {
    l = powo.list[1,]
    powo.max = nrow(powo.list)
    powo.list = powo.list[2:powo.max,]
    
    points.sub = subset(tdwg.extract, subset = (tdwg.extract["POWO_ID"] ==
                                                  l))
    #points.sub = tdwg.extract[ which==powo.list,]
    tdwg.sub = subset(result.tdwg, subset = (result.tdwg["POWO_ID"] == l))
    vars = "LEVEL3_NAM"
    tdwg.sub = tdwg.sub[vars]
    NativeWarning = ""
    first.merge = merge(points.sub, tdwg.sub, by = "Level.3.code")
    first.merge = cbind(first.merge, NativeWarning = NativeWarning)
    
    if (nrow(first.merge) < 1) {
      #add blank row with empty results
      NativeWarning = "No GBIF points in native range"
      first.merge = cbind(tdwg.extract, NativeWarning = NativeWarning)
    }
    
    for (l in powo.list) {
      #l = "77110432-1"
      points.sub = subset(tdwg.extract, subset = (tdwg.extract["POWO_ID"] ==
                                                    l))
      #points.sub = tdwg.extract[ which==powo.list,]
      tdwg.sub = subset(result.tdwg, subset = (result.tdwg["POWO_ID"] == l))
      vars = "LEVEL3_NAM"
      tdwg.sub = tdwg.sub[vars]
      list.merge = merge(points.sub, tdwg.sub, by = "LEVEL3_NAM")
      NativeWarning = ""
      list.merge = cbind(list.merge, NativeWarning = NativeWarning)
      list.merge = rbind(list.merge, first.merge)
      return(list.merge)
    }
  }
}

# 4.3 query ecoregions to get ecoregion count
LC.extract.eco = function(unique.sp, native_only, raster.eco)
  if (nrow(native_only) < 1) {
    
  }  else {
    unique.sp = unique.sp[1]
    species.points = subset(native_only, subset = native_only$POWO_ID ==
                              unique.sp)
    decimalLongitude =  species.points$decimalLongitude
    decimalLatitude =  species.points$decimalLatitude
    coords = cbind(decimalLongitude, decimalLatitude)
    coords <- data.matrix(coords)
    coords = unique(coords)
    sp = SpatialPoints(coords)
    
    # extract values from raster for each point lat long
    extract.out = raster::extract(raster.eco, sp, df = TRUE)
    extract.out = cbind(extract.out, coords)
    
    ecocount = length(unique(extract.out$eco2017))
    eco.res = data.frame(EcoregionCount = ecocount, POWO_ID = unique.sp)
    return(eco.res)
  }

# 5 calculate EOO and AOO
LC.eoo.aoo = function(unique.sp, native_only) {
  #get the points and project
  sp = subset(native_only, subset = native_only$POWO_ID == unique.sp)
  sp.ddlong = sp$decimalLongitude
  sp.ddlat = sp$decimalLatitude
  mypointsll = data.frame(lat = sp.ddlat, long = sp.ddlong)
  centreofpoints <- trueCOGll(mypointsll)
  mypointsxy <- simProjWiz(mypointsll, centreofpoints)
  
  #EOO and AOO calculation
  EOOm2 <- EOOarea(mypointsxy)
  EOOkm2 <- EOOm2 / 1000000
  EOOkm2abs = abs(EOOkm2)
  rec_count = nrow(mypointsxy)
  cellsizem <- 10000
  AOOnocells <- AOOsimp (mypointsxy, cellsizem)
  
  # extra field to find if any recent collection
  coll.last10yrs =
    
    eoo.aoo.res = data.frame(
      EOO = EOOkm2abs,
      AOO = AOOnocells,
      RecordCount = rec_count,
      POWO_ID = unique.sp
    )
  
  #eoo.aoo.res = do.call(rbind, res)
  return(eoo.aoo.res)
}

#6 get IUCN point file
LC.point.file = function(unique.sp, native_only, result.table) {
  #unique.sp = unique.sp[1]
  #native_only = native_only[,1]
  pointsLC = subset(native_only, subset = native_only$POWO_ID == unique.sp)
  GetYrCompiled = substring(Sys.Date(), 1, 4)
  GetSpatialRef = "WGS84"
  gbif_dataset = 'http://www.gbif.org/dataset/'
  res.tab = subset(result.table, subset = result.table$POWO_ID == unique.sp)
  Binomial = res.tab$GBIF_SuggestedSpe
  respath = paste0(path, Binomial, ".csv")
  iucn_point_file = data.frame(
    Binomial = Binomial,
    Presence = '1',
    Origin = '1',
    Seasonal = '1',
    Compiler = pointcompiler,
    Year = GetYrCompiled,
    Dec_Lat = pointsLC$decimalLatitude,
    Dec_Long = pointsLC$decimalLongitude,
    SpatialRef = GetSpatialRef,
    Event_Year = pointsLC$year,
    Citation = pointcitation,
    Source = paste0(gbif_dataset, pointsLC$datasetKey),
    BasisOfRec = pointsLC$basisOfRecord,
    CatalogNo = pointsLC$catalogNumber,
    recordedBy = pointsLC$recordedBy,
    recordNumber = pointsLC$recordNumber
  )
  FinalResult = write.table(
    iucn_point_file,
    respath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
}

#7 generate the SIS files
LC.credits = function(LC.results) {
  #credits = data.frame(LC.results)[1,]
  powoid = LC.results[["POWO_ID"]]
  credits = data.frame(
    internal_taxon_id = powoid,
    credit_type = credittype,
    firstName = Firstname,
    lastName = Lastname,
    initials = Initials,
    Order = "1",
    email = Email,
    affiliation = Affiliation,
    user_id = "1"
  )
  return(credits)
}
LC.taxonomy = function(IUCN_taxonomy, LC.results) {
  #LC.results = result.table
  
  # filter to get those where GBIF and IUCN families matched and use IUCN higher taxonomy
  # leave those that don't match as blank and let user decide
  tax = full_join(IUCN_taxonomy, LC.results, by = "fam", copy = FALSE)
  tax_matched = subset(
    tax,
    select = c(
      Kingdom,
      Phylum,
      Class ,
      Order,
      fam,
      IUCN_fam,
      POWO_ID,
      GBIF_SuggestedGen,
      GBIF_SuggestedSpe,
      GBIF_Suggestedauth
    )
  )
  tax_good = tax_matched[complete.cases(tax_matched),]
  tax_unmatched = subset(tax_matched, (!is.na(tax_matched[, 7])))
  tax_unmatched = subset(tax_unmatched, (is.na(tax_unmatched[, 6])))
  tax_comb = rbind(tax_good, tax_unmatched)
  tax_comb$internal_taxon_id = tax_comb$POWO_ID
  tax_comb$kingdom = tax_comb$Kingdom
  tax_comb$phylum = tax_comb$Phylum
  tax_comb$classname = tax_comb$Class
  tax_comb$ordername = tax_comb$Order
  tax_comb$family = tax_comb$IUCN_fam
  tax_comb$genus = tax_comb$GBIF_SuggestedGen
  tax_comb$species = tax_comb$GBIF_SuggestedSpe
  tax_comb$taxonomicAuthority = tax_comb$GBIF_Suggestedauth
  tax_comb = tax_comb %>% separate(species, into = c("gen", "species"), sep = " ") %>% as.data.frame()
  tax_comb = subset(tax_comb, select = -c(1:10, 18))
  return(tax_comb)
}
LC.biorealms = function(TDWG_realms, result.tdwg) {
  #add realms using TDWG_realms table - first change column name to get match with result.tdwg
  colnames(TDWG_realms) [which(names(TDWG_realms) == "LEVEL3_COD")] = "tdwgCode"
  #merge result.tdwg and TDWG_realms so we have realms linked to POWO ID
  realms = merge(TDWG_realms, result.tdwg, by.x = "tdwgCode", by.y = "tdwgCode")
  # summarise and collapse to get single column with all unique realms - nice!
  biorealm_summary = realms %>% group_by(POWO_ID) %>% summarise(newREALM = paste(unique(REALM), collapse =
                                                                                   ","))
  make_biorealms = data.frame(POWO_ID = biorealm_summary$POWO_ID,
                              biogeographicrealm.realm = (gsub("," , "\\|", biorealm_summary$newREALM)))
  return(make_biorealms)
  
}
LC.assessments = function(LC.results, biorealms) {
  #result.table = as.data.frame(result.table)[1,]
  #assessLC = subset(result.table, subset= result.table$POWO_ID==unique.sp)
  #join biorealms to assessments and update
  LC.results = merge(LC.results, biorealms, by = "POWO_ID")
  biogeographicrealm.realm = LC.results[["biogeographicrealm.realm"]]
  Rating = 'LC'   #Rating
  Rationale = paste(
    LC.results[["full_name"]],
    " has a wide distribution range across ",
    LC.results$TDWGCount,
    " countries and is known from ",
    LC.results$RecordCount,
    " occurrence records. Due to the wide range and lack of any major threats to this species at present the rating of Least Concern is the most appropriate category.",
    sep = "",
    collapse = NULL
  )
  PopTrend = 'Stable'
  System = 'Terrestrial'
  mapstatus = 'Incomplete'
  sysdate = Sys.Date()
  day = substr(sysdate, 9, 10)
  month = substr(sysdate, 6, 7)
  year = substr(sysdate, 1, 4)
  date = paste0(day, "/", month, "/", year)
  version = "3.1"
  language = "English"
  internal_taxon_id = LC.results[["POWO_ID"]]
  range_nar = ""
  pop_nar = ""
  hab_nar = ""
  threat_nar = ""
  redlistcriteria.ismanual = "TRUE"
  
  assessments = data.frame(
    internal_taxon_id = internal_taxon_id,
    RedListRationale.value = Rationale,
    mapstatus.status = mapstatus,
    RedListAssessmentDate.value = date,
    RedListCriteria.critVersion = version,
    RedListCriteria.manualCategory = Rating,
    PopulationTrend.value = PopTrend,
    System.value = System,
    language.value = language,
    rangedocumentation.narrative = range_nar,
    populationdocumentation.narrative = pop_nar,
    habitatdocumentation.narrative = hab_nar,
    threatsdocumentation.value = threat_nar,
    redlistcriteria.ismanual = redlistcriteria.ismanual,
    biogeographicrealm.realm = biogeographicrealm.realm
  )
  
  return(assessments)
}
LC.habitats = function(LC.results) {
  powoid = LC.results[["POWO_ID"]]
  habitats = data.frame(
    internal_taxon_id = powoid,
    GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsLookup = "18",
    GeneralHabitats.GeneralHabitatsSubfield.GeneralHabitatsName = "Unknown",
    GeneralHabitats.GeneralHabitatsSubfield.majorImportance = "",
    GeneralHabitats.GeneralHabitatsSubfield.season = "",
    GeneralHabitats.GeneralHabitatsSubfield.suitability = "Unknown"
  )
  return(habitats)
}
LC.plantspecific = function(LC.results) {
  powoid = LC.results[["POWO_ID"]]
  plantspecific = data.frame(
    internal_taxon_id = powoid,
    PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsLookup = "",
    PlantGrowthForms.PlantGrowthFormsSubfield.PlantGrowthFormsName = ""
  )
  return(plantspecific)
}
LC.allfields = function(LC.results) {
  powoid = LC.results[["POWO_ID"]]
  trendderiv = "Suspected"
  eoores = LC.results[["EOO"]]
  eoojust = "Minimum Convex Polygon (MCP) based on attached point file."
  nothreats.nothreats = "TRUE"
  threatsunknown.value = "FALSE"
  allfields = data.frame(
    internal_taxon_id = powoid,
    CurrentTrendDataDerivation.value = trendderiv,
    EOO.range = eoores,
    EOO.justification = eoojust,
    nothreats.nothreats = nothreats.nothreats,
    threatsunknown.value = threatsunknown.value
  )
  return(allfields)
}
LC.countries = function(unique.sp, result.tdwg) {
  #unique.sp = unique.sp[1]
  result.tdwg.reduced = subset(result.tdwg, subset = result.tdwg$POWO_ID ==
                                 unique.sp)
  colnames(TDWG_to_IUCN_version1_UTF_8)[which(names(TDWG_to_IUCN_version1_UTF_8) == "Level.3.code")] =
    "tdwgCode"
  #country.merge = merge(result.tdwg.reduced,TDWG_to_IUCN_version1_UTF_8, by = "tdwgCode")
  country.merge = distinct(result.tdwg.reduced, countryoccurrencename, .keep_all = TRUE)
  powoid = country.merge$POWO_ID
  countries = data.frame(
    internal_taxon_id = powoid,
    CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceLookup = country.merge$countryoccurrencelookup,
    CountryOccurrence.CountryOccurrenceSubfield.CountryOccurrenceName = country.merge$countryoccurrencename,
    CountryOccurrence.CountryOccurrenceSubfield.presence = 'Extant',
    CountryOccurrence.CountryOccurrenceSubfield.origin =
      country.merge$establishment,
    CountryOccurrence.CountryOccurrenceSubfield.seasonaility = 'Resident'
  )
  return(countries)
}
LC.references = function(LC.results){
  #powoid = trees_7079[["ID"]]
  #powoid = powoid[1:10]
  powoid = LC.results[["POWO_ID"]]
  
  type = c("electronic source","electronic source","electronic source" )
  author = c("Board of Trustees, RBG Kew","Moat, J.", "Chamberlain, S")
  year = c("2018","2017", "2017")
  title = c("Plants of the World Online Portal", "rCAT: Conservation Assessment Tools. R package version 0.1.5.","rgbif: Interface to the Global 'Biodiversity' Information Facility API. R package version
            0.9.9.")
  place_published = c("Richmond, UK", "", "")
  url = c("http://powo.science.kew.org/","https://CRAN.R-project.org/package=rCAT", "https://CRAN.R-project.org/package=rgbif")
  reference_type = c("Assessment","Assessment", "Assessment")
  references = data.frame(type, author, year, title, place_published, url, reference_type)
  references$internal_taxon_id = powoid[1] 
  
  # make loop to create table with refs below and add powoid to each table, then bind together
  list = powoid[-1]
  
  for (l in list){
    type = c("electronic source","electronic source","electronic source" )
    author = c("Board of Trustees, RBG Kew","Moat, J.", "Chamberlain, S")
    year = c("2018","2017", "2017")
    title = c("Plants of the World Online Portal", 
              "rCAT: Conservation Assessment Tools. R package version 0.1.5.",
              "rgbif: Interface to the Global 'Biodiversity' Information Facility API. R package version
              0.9.9.")
    place_published = c("Richmond, UK", "", "")
    url = c("http://powo.science.kew.org/","https://CRAN.R-project.org/package=rCAT", "https://CRAN.R-project.org/package=rgbif")
    reference_type = c("Assessment","Assessment", "Assessment")
    all.references = data.frame(type, author, year, title, place_published, url, reference_type)
    all.references$internal_taxon_id = l 
    references = rbind(references,all.references)
  }
  
  return(references)
  
}
LC.SIS.all = function(LC.results, unique.sp,result.tdwg,IUCN_taxonomy,biorealms) {
  sis.credits = LC.credits(LC.results)
  creditspath = paste0(path, "credits", ".csv")
  save.SIScredits = write.table(
    sis.credits,
    creditspath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
  
  sis.taxonomy = LC.taxonomy(LC.results, IUCN_taxonomy)
  taxonomypath = paste0(path, "taxonomy", ".csv")
  save.SIStaxonomy = write.table(
    sis.taxonomy,
    taxonomypath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
  
  sis.assessments = LC.assessments(LC.results, biorealms)
  assessmentspath = paste0(path, "assessments", ".csv")
  save.SISassessments = write.table(
    sis.assessments,
    assessmentspath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
  
  sis.habitats = LC.habitats(LC.results)
  habitatspath = paste0(path, "habitats", ".csv")
  save.SIShabitats = write.table(
    sis.habitats,
    habitatspath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
  
  sis.plantspecific = LC.plantspecific(LC.results)
  plantspecificpath = paste0(path, "plantspecific", ".csv")
  save.SISplantspecific = write.table(
    sis.plantspecific,
    plantspecificpath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
  
  sis.allfields = LC.allfields(LC.results)
  allfieldspath = paste0(path, "allfields", ".csv")
  save.SISallfields = write.table(
    sis.allfields,
    allfieldspath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
  
  sis.countries = lapply(unique.sp, LC.countries, result.tdwg)
  sis.countries = do.call(rbind, sis.countries)
  countriespath = paste0(path, "countries", ".csv")
  save.SIScountries = write.table(
    sis.countries,
    countriespath,
    row.names = FALSE,
    na = "",
    sep = ","
  )
  
  sis.references = LC.references(LC.results)
  referencespath = paste0(path, "references", ".csv")
  save.SISreferences = write.table(sis.references,referencespath, row.names = FALSE,na = "",sep = ",")
}

#8 Add result status: LC or 'Possibly threatened' based on predefined thresholds
#This is where you add the LC thresholds and decide what is LC e.g. >5 ecoregions, EOO > 40,000 etc.

LC.results = function(result.table) {
  #result.table = as.data.frame(result.table)
  dplyr::mutate(
    result.table,
    TDWGResult = ifelse(TDWGCount >= 4 & !is.na(TDWGCount), "LC", "Possibly threatened"),
    EcoregionResult = ifelse(EcoregionCount >= 5 & !is.na(EcoregionCount), "LC", "Possibly threatened"),
    EOOResult = ifelse(EOO > 30000 & !is.na(EOO), "LC", "Possibly threatened"),
    AOOResult = ifelse(AOO > 30 & !is.na(AOO), "LC", "Possibly threatened"),
    RecordResult = ifelse(RecordCount > 40 & !is.na(RecordCount), "LC", "Possibly threatened")
  )
}

LC.results.errors = function(result.table.errors){
  #add the extra rows here
  result.table.errors$TDWGResult = ""
  result.table.errors$EcoregionResult = ""
  result.table.errors$EOOResult = ""
  result.table.errors$AOOResult = ""
  result.table.errors$RecordResult = ""
  return(result.table.errors)
}

#8.1 make empty columns to ensure rbind works for good error and good species frames << NOT WORKING
error.tab = function(key.errors) {
  key.errors$TDWGCount = ""
  key.errors$EcoregionCount = ""
  key.errors$EOO = ""
  key.errors$AOO = ""
  key.errors$RecordCount = ""
  result.table.1 = key.errors
  rm(result.table)
}

#8.2 something extra - LC.variables - only if we have native points
# now I have area of TDWG, forest loss, weighted forest loss, human footp and weighted human footp
# I need to get sum of area and forest loss, human footp per species range e.g. one or more TDWG regions
# I need to divide sum of forest loss/hfoot by sum of area
# these are then normalised by area for each species and should then update the final result table
LC.variables = function(result.tdwg.native, unique.sp) {
  #normalise tdwg variables by area (tdwg L3 regions are different sizes)
  
  # make sure relevant fields are numbers - it won't summarise on characters
  result.tdwg.native = transform(result.tdwg.native, tdwgLevel = as.integer(tdwgLevel),
                                 ID = as.integer(ID),
                                 area = as.numeric(area),
                                 forest_loss = as.integer(forest_loss),
                                 hfoot = as.numeric(hfoot),
                                 w_forest_loss = as.numeric(w_forest_loss),
                                 w_hfoot = as.numeric(w_hfoot)
  )
  
  extras = result.tdwg.native %>% group_by(POWO_ID) %>% summarise(
    sum_w_forest_loss = sum(w_forest_loss, na.rm = TRUE),
    sum_w_hfoot = sum(w_hfoot, na.rm = TRUE),
    sum_area = sum(area, na.rm = TRUE),
    norm_forest_loss = sum_w_forest_loss / sum_area,  
    norm_hfoot = sum_w_hfoot /  sum_area)
  extras = subset(extras, select = c(POWO_ID, norm_forest_loss, norm_hfoot))
}

#########################################################################
#9 These lines bring it all together ####################################
#########################################################################

RedLeast = function(full_name,ID_list,path,LC.points = FALSE,SIS.files = FALSE, BigKeys) {

  # 1 - Run the full_name against GBIF query to see if you get a name match. Return warning value in table if no match
  keys.tax = LC.get.keys.tax(full_name, ID_list)
  key.errors = subset(keys.tax, (nchar(keys.tax$Warning) > 3))
  key.matches = subset(keys.tax, (nchar(keys.tax$Warning) < 3))
  
  # Run ID (POWO IDs) against POWO api to get TDWG regions
  result.tdwg = LC.tdwg.apply(ID_list) 
  # Find the errors where there is no TDWG info in POWO
  tdwg.errors = LC.tdwg.errors(result.tdwg)
  # Use update.errors function to add error value to warning column
  result.table = update.errors(keys.tax, tdwg.errors, error.text = "No TDWG distribution data")
  #check for empty TDWG results and make table to join to main results at the end
  error.species = subset(result.table, (nchar(result.table$Warning) > 3))
  good.species = subset(result.table, (nchar(result.table$Warning) < 3))
  # print to check progress

  if (nrow(error.species) >= 1) {
    error.species$TDWGCount = ""
    error.species$EcoregionCount = ""
    error.species$EOO = ""
    error.species$AOO = ""
    error.species$RecordCount = ""
    # add/remove next two lines below when you do/don't want to calcualte forest loss, hfoot etc. 
    error.species$norm_forest_loss = ""
    error.species$norm_hfoot = ""
    result.table.errors = error.species
    rm(result.table)
  }
  
  # 2 - are there GBIF point errors i.e. no georeferneced occurrence points
  if (nrow(good.species) >= 1) {
    result.table = good.species
    gbifkey = result.table[, 3]
    
    # count occurrences for each 'good species' and combine with gbif key
    count_occs = lapply(gbifkey,LC_count_occs)
    count_occs_tibb = tibble(count_occs) %>% unlist()
    gbifkey_df = as.data.frame(gbifkey)
    all_count = cbind(count_occs_tibb,gbifkey_df)
    
    # Filter out the big keys
    small_gbifkeys = subset(all_count, count_occs <=200000) %>% as.data.frame() # change to whatever is the threshold e.g. 200,000
    small_gbifkeys = small_gbifkeys[,2]
    big_keys = subset(all_count, count_occs >200000) %>% as.data.frame() # change to whatever is the threshold e.g. 200,000?
    big_keys = big_keys[,2]
    print(paste0(big_keys))
    
    # not get the points
    apply.points = lapply(small_gbifkeys, LC.gbif.points) #if GBIF keys exists....
    points = do.call(rbind.fill, apply.points)
  
    # add POWO_ID to points because POWO ID is the unique value. join points.with.tdwg to result.table.tdwg
    powo.ids = cbind(POWO_ID = result.table$POWO_ID, taxonKey = result.table$GBIF_SuggestedKey) %>% as.data.frame()
    points.with.powoid = merge(points, powo.ids, by = "taxonKey")
    smr = ddply(points.with.powoid, "POWO_ID",summarise,n = length(decimalLatitude),n_na = sum(is.na(decimalLatitude)), prop_missing = n_na / n)
    gbif.point.errors = subset(smr, prop_missing == 1)
    gbif.point.errors = gbif.point.errors[1]
    result.table = update.errors(result.table, gbif.point.errors, error.text = "No GBIF points")

    # another break point to filter out error species that have no points
    noPoints.species = subset(result.table, (nchar(result.table$Warning) > 3))
    withPoints.species = subset(result.table, (nchar(result.table$Warning) < 3))

    if (nrow(noPoints.species) >= 1) {
      noPoints.species$TDWGCount = ""
      noPoints.species$EcoregionCount = ""
      noPoints.species$EOO = ""
      noPoints.species$AOO = ""
      noPoints.species$RecordCount = ""
      # add/remove next two lines below when you do/don't want to calcualte forest loss, hfoot etc.
      noPoints.species$norm_forest_loss = ""
      noPoints.species$norm_hfoot = ""
      rm(result.table)
    }
    
    # check if noPoints has any rows - if not, do nothing
    if (nrow(noPoints.species) >=1) {
      
      # if noPoints has rows, then check if error table exists 
      if (exists("result.table.errors")) {
        
        # if error table exists, then rbind noPoints
        result.table.errors = rbind(result.table.errors,noPoints.species)
      } else {
        
        # otherwise make error table from points
        result.table.errors = noPoints.species
      }
    }
    
    
    #3 - Check species points against species' native TDWG range
    
    if (nrow(withPoints.species) >= 1) {
      #now I only need clean points i.e. use merge on POWO ID to remove taxa with no points 
      points.with.powoid = semi_join(points.with.powoid, withPoints.species, by = "POWO_ID", copy = FALSE)    
      points.clean.df = clean.points.df(points.with.powoid)
      rm(points)
      rm(points.with.powoid)
      tdwg.extract = LC.extract.tdwg(points.clean.df, raster.tdwg) 
      rm(points.clean.df)
      powo.list = unique(tdwg.extract[15]) %>% as.data.frame()
      result.tdwg.native = subset(result.tdwg, establishment == "Native")
      colnames(result.tdwg.native)[which(names(result.tdwg.native) == "TDWG_Name")] = "Native_NAM"
      #colnames(result.tdwg)[which(names(result.tdwg) == "LEVEL3_NAM")] = "Native_NAM"
      countTDWG = ddply(result.tdwg.native,"POWO_ID", summarise,TDWGCount = length(tdwgCode))
      result.table = merge(withPoints.species, countTDWG, by = "POWO_ID")
  
      #4 - run final metrics on clean point data from native range
      # find the native only point
      native_only = inner_join(tdwg.extract,result.tdwg.native[, c("POWO_ID", "Native_NAM")], by = c("POWO_ID" = "POWO_ID", "LEVEL3_NAM" = "Native_NAM"))
      
      # now find the non-native points and compare with native points. 
      tdwg.extract.unmatched = anti_join(tdwg.extract,result.tdwg.native, by = c("POWO_ID" = "POWO_ID", "LEVEL3_NAM" = "Native_NAM"))
      
      #If prop. non-native is 1 then mark as error 
      sum.tdwg.extract = ddply(tdwg.extract, "POWO_ID", summarise,n = length(decimalLatitude))
      sum.unmatched = ddply(tdwg.extract.unmatched, "POWO_ID", summarise, n_non_native = length(decimalLatitude))
      sum.merged = merge(sum.tdwg.extract,sum.unmatched, all = TRUE)
      sum.merged$prop_non_native = (sum.merged$n_non_native / sum.merged$n)
      gbif.native.errors = subset(sum.merged, prop_non_native == 1)
      result.table = update.errors(result.table, gbif.native.errors[1], error.text = "No GBIF points in range")
      
      # get 'points out of range' errors and add them to result.table.errors, then leave clean list for final analysis
      PointsOutofRange.species = subset(result.table, (nchar(result.table$Warning) > 3))
      PointsInRange.species = subset(result.table, (nchar(result.table$Warning) < 3))
      
      # need to update result.table.errors here with points out of range
      if (nrow(PointsOutofRange.species) >= 1){
        
        PointsOutofRange.species$EcoregionCount = ""
        PointsOutofRange.species$EOO = ""
        PointsOutofRange.species$AOO = ""
        PointsOutofRange.species$RecordCount = ""
        # add/remove next two lines below when you do/don't want to calcualte forest loss, hfoot etc.
        PointsOutofRange.species$norm_forest_loss = ""
        PointsOutofRange.species$norm_hfoot = ""
        rm(result.table)
      }
      
      # check if PointsOutofRange.species has any rows - if not, do nothing
      if (nrow(PointsOutofRange.species) >=1) {
        
        # if noPoints has rows, then check if error table exists 
        if (exists("result.table.errors")) {
          
          # if error table exists, then rbind noPoints
          result.table.errors = rbind(result.table.errors,PointsOutofRange.species)
        } else {
          
          # otherwise make error table from points
          result.table.errors = PointsOutofRange.species
        }
      }
      
      # Get list of unique species by POWO_IDs and run ecoregions function to get count of ecoregions per species
      if (nrow(PointsInRange.species) >= 1) {
        
        native.points = semi_join(native_only, PointsInRange.species, by = "POWO_ID", copy = FALSE)
        unique.sp = unique(native.points$POWO_ID)
        eco.res = lapply(unique.sp, LC.extract.eco, native.points, raster.eco)
        eco.res = do.call(rbind, eco.res)
        # Join the ecoregion count to the result table
        result.table = full_join(PointsInRange.species, eco.res, by = "POWO_ID", copy = FALSE)
        
        # Run the native only points against the EOO and AOO scripts using rCAT package
        eoo.aoo.res = lapply(unique.sp, LC.eoo.aoo, native_only)
        eoo.aoo.res = do.call(rbind, eoo.aoo.res)
        # Join the EOO and AOO results to the result table
        result.table = full_join(result.table, eoo.aoo.res, by = "POWO_ID", copy = FALSE)
        
        # Calcualte biorealms for each species
        biorealms = LC.biorealms(TDWG_realms, result.tdwg)
        
        # now add the extra variables - forest loss and hfoot
        normalise.vars = LC.variables(result.tdwg.native,unique.sp)
        unique.sp.vars = as.data.frame(unique.sp)
        colnames(unique.sp.vars)[which(names(unique.sp.vars) == "unique.sp")] = "POWO_ID"
        var.update = merge(unique.sp.vars,normalise.vars, by = "POWO_ID")
        result.table = full_join(result.table, var.update, by = "POWO_ID", copy = FALSE)
        
        #if user declares LC.points = TRUE then save the point csv's
        #note that it doesn't mean species is LC - it will save all those that have point data
        if (LC.points) {
          lapply(unique.sp, LC.point.file, native_only, result.table)
        }
        
        #if user declares SIS.files = TRUE then save the SIS csv's
        if (SIS.files) {
          LC.SIS.all(result.table,unique.sp,result.tdwg,IUCN_taxonomy,biorealms)
        }
        
        
      }
    }
  }
  
  #now save the results file - checks for either error or good species only and saves
  #Or - if both error and good species exist, rbind these and then save
  
  # make below into a function e.g.save.results = function()
  respath = paste0(path, "Results", ".csv")
  if (!exists("result.table.errors")) {
    #result.table = LC.results(result.table)
    save.results = write.table(result.table, respath,row.names = FALSE,na = "", sep = ",")
  } else if (!exists("result.table")) {
    #result.table.errors = LC.results(result.table.errors)
    save.results = write.table(result.table.errors,respath,row.names = FALSE,na = "",sep = ",")
  } else {
    #result.table = LC.results(result.table)
    #result.table.errors = LC.results.errors(result.table.errors)
    result.table = rbind(result.table, result.table.errors)
    save.results = write.table(result.table,respath,row.names = FALSE,na = "",sep = ",")
  }

  print(paste0("Results written to disk"))
  time_now = Sys.time()
  print(time_now)
  
  
  }


