#######################################################################################################################
#Code from Andy Houseman to process IDAT files into FunNorm-normalized signals
#######################################################################################################################

#Necessary packages
  library(minfi) #1.12.0
  library(IlluminaHumanMethylation450kmanifest) #0.4.0
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #0.2.1

#Location of the IDAT files
setwd("/Users/kevinjohnson/Documents/NDRI_Breast_5hmC/Preprocessing_Data/IDAT_Files")
  chips <- list.files(pattern="^[[:digit:]]*") #Find numeric files/Chips 
  nChips <- length(chips) #How many unique microarrays

#Locate associated IDAT files, create list.
  arrayList <- list()
  for(i in 1:nChips){
    ff <- grep("Grn.idat$",list.files(path=chips[i],pattern=chips[i]),value=TRUE)
    arrayList[[i]] <- paste(chips[i],gsub("[_]Grn.idat","",ff),sep="/")
  }
  names(arrayList) <- chips

#Manifest file that also contains som indication of oxBS and BS treatment
setwd("/Users/kevinjohnson/Documents/NDRI_Breast_5hmC/I.Preprocessing_Data/IDAT_Files")
  meta <- read.delim("metaData-breast-150318.txt", sep="\t", head=TRUE, 
  stringsAsFactors=FALSE)
  meta$ArrayID <- with(meta, paste(Beadchip,Position,sep="_"))
  meta$Label.ID.[meta$OX_BS==1] <- gsub("B$","",meta$ID[meta$OX_BS==1]) #Samples with 'B' in name were oxBS treated samples
  meta$Label.ID.[meta$OX_BS==0] <- gsub("A$","",meta$ID[meta$OX_BS==0]) #'A' in sample name indicated BS treated sample

  with(meta, all(Label.ID.==gsub("[A|B]$","",ID))) # Check consistency of IDs

#Merge the oxBS and BS assays for the same sample
  mrg <- merge(
    data.frame(LabelID=meta$Label.ID.[meta$OX_BS==0],
       ixBS=as.numeric(rownames(meta)[meta$OX_BS==0]),
       stringsAsFactors=FALSE),
    data.frame(LabelID=meta$Label.ID.[meta$OX_BS==1],
       ixOxBS=as.numeric(rownames(meta)[meta$OX_BS==1]),
       stringsAsFactors=FALSE)
)

#Identify the indices of a paired sample (i.e. B10A in row 1, matches B10B in row 2)
  meta$Pair[mrg$ixOxBS] <- mrg$ixBS 
  meta$Pair[mrg$ixBS] <- mrg$ixOxBS 
#'ArrayID' represents the location of the sample in which well and on what slide 
  rownames(meta) <- meta$ArrayID

#Select those arrays with BS only treatment
arraysBS <- intersect(unlist(arrayList),
  with(meta, paste(Beadchip,ArrayID,sep="/")[OX_BS==0]))

#Select those arrays with oxBS treatment
arraysOxBS <- intersect(unlist(arrayList),
  with(meta, paste(Beadchip,ArrayID,sep="/")[OX_BS==1]))

#####################################
setwd("/Users/kevinjohnson/Documents/NormalBreast_Glioma/Preprocessing_Data/IDAT_Files")
# Read IDAT file using the minfi functions
  datListBS <- read.450k(arraysBS,verbose=TRUE)
  datListOxBS <- read.450k(arraysOxBS,verbose=TRUE)
######################################

#Andy rewrote the FunNorm function so that he could pull out the functional-normalization
#signals (as opposed to betas and M-values). Later Andy will use these combined signals to
#estimate 5mC and 5hmC using a maximum-likelihood approach
preprocessFunnormRedGreen <- function (rgSet, nPCs = 2, sex = NULL, verbose = TRUE)
{
    minfi:::.isRG(rgSet)
    rgSet <- updateObject(rgSet)
    if (verbose)
        cat("[preprocessFunnorm] Mapping to genome\n")
    gmSet <- mapToGenome(rgSet)
    subverbose <- max(as.integer(verbose) - 1L, 0)
    if (verbose)
        cat("[preprocessFunnorm] Quantile extraction\n")
    extractedData <- minfi:::.extractFromRGSet450k(rgSet)
    if (is.null(sex)) {
        gmSet <- addSex(gmSet, getSex(gmSet, cutoff = -3))
        sex <- rep(1L, length(gmSet$predictedSex))
        sex[gmSet$predictedSex == "F"] <- 2L
    }
    rm(rgSet)
    if (verbose)
        cat("[preprocessFunnorm] Normalization\n")
    CN <- getCN(gmSet)
    minfi:::.normalizeFunnorm450k(object = gmSet, extractedData = extractedData,
        sex = sex, nPCs = nPCs, verbose = subverbose)
}
  rgBS <-  preprocessFunnormRedGreen(datListBS)
  rgOxBS <-  preprocessFunnormRedGreen(datListOxBS)
  
#'Methylated' channel signal
  methBS <- assay(rgBS,"Meth")
  methOxBS <- assay(rgOxBS,"Meth")

#'Unmethylated' channel signal
  unmethBS <- assay(rgBS,"Unmeth")
  unmethOxBS <- assay(rgOxBS,"Unmeth")

#Save output. Next step usually requires running a job on Discovery,
setwd("/Users/kevinjohnson/Documents/NormalBreast_Glioma/Preprocessing_Data")
save(list=c("methBS","methOxBS","unmethBS","unmethOxBS"), file="OxyBreastSignals-FunNorm.RData")
