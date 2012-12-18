### Andrew Trister
### Sage Bionetworks
### Seattle, WA
### 20121217

require(synapseClient)



setwd("~/GBM-Age/")

#login to Synapse
synapseLogin()

## start by loading all of the data
source("./src/populateData.R")

#now start looking at the differences in Age using elastic net in TCGA
source("./src/enet-Age.R")