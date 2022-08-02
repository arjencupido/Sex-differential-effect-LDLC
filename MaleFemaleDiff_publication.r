# Custom functions used to analyze the data for: https://www.ahajournals.org/doi/10.1161/JAHA.121.024248
# For the preparation of the dataset for analysis, see https://github.com/arjencupido/UKB-pipeline-Utrecht
############################

# Load packages ---------
library(dplyr)
library(sjmisc)
library(ggplot2)
library(meta)
library(tableone)
library(ggpubr)
library(jtools)
library(broom)
library(survival)
library(survminer)
library(readxl)
library(reshape2)
library(forestplot)
library(gridExtra)
library(rcompanion)
library(DescTools)

# Extra functions ----------------------

# The 'not in' function
`%notin%` <- Negate(`%in%`)


# Function for aligning UKB variants so that the score will be calculated per increasing weight (e.g. increasing LDL-C).
UKB_WeighVarPos <- function(d, UKB_variants_list, weights_list ){
  # Takes a weights list, which is a dataframe with RSID, effect allele, alternative allele, beta, and prefereably effect allele frequency
  # Takes a dataframe with all the variants from the UK Biobank, which can be extracted using Bgen (the --list option)
  # Takes the dataframe with all the individual data on genetic variants.

  varchanges <- list()

  for (i in 1:length(weights_list$rsid)){

    v <- weights_list$rsid[i] # extract RSID

    weight <- subset( weights_list, rsid == paste0(v) ) # Extract external GWAS data for RSID

    weight$POS_EXP <- ifelse( weight$BETA_EXP > 0, toupper( paste0( weight$REF_EXP )), toupper( paste0( weight$ALT_EXP ))) # We will align on the allele with an increasing effect. Here we decide which allele is the positive effect allele.

    UKBvar <- subset( UKB_variants_list, rsid == paste0(v) ) # Extract UKB data for RSID

    varchanges[[i]] <- data.frame("status" = NA,"RSID" = v, # dataframe to check alignment afterwards
                                  "REF_EXP" = weight$REF_EXP, "ALT_EXP" =  weight$ALT_EXP, "POS_EXP" = weight$POS_EXP,
                                  "REF_UKB" =  UKBvar$alternative_alleles, "ALT_UKB" =  UKBvar$first_allele, "NumberOfAlleles" = UKBvar$number_of_alleles,
                                  "EAFUKB" = mean(d[, paste0(v)], na.rm=T)/2, "EAF_RECODED"=NA)
    xx <- data.frame()

    if( as.character( weight$POS_EXP ) != as.character( UKBvar$alternative_alleles )){ # Recoding of variant if necessary. The alternative allele is the counted allele in the UK Biobank (e.g. 0,1,2 = 0, 1 or 2 alternative alleles)

      xx <- d[, paste0(v)]
      xx[xx == 2] <- 99
      xx[xx == 0] <- 2
      xx[xx == 99] <- 0

      varchanges[[i]][1] <- "changed"

    } else {

      xx <- d[, paste0(v)]

      varchanges[[i]][1] <- "same"
    }

    #Create allele with Weight

    d[, paste0("r_", v)] <- xx # unweighted allele variable

    d[, paste0("w_" ,v)] <- xx * abs(as.numeric(weight$BETA_EXP)) # Weighted allele variable

    varchanges[[i]][10] <- mean(d[, paste0("r_", v)], na.rm=T)/2 # Calculate positive effect allele frequency in UKB

  }
  varchanges <- do.call( rbind, varchanges ) # Create dataframe with metadata from the recoding above.

  return( list( 'data' = d,
                "variant_changes" = varchanges ))
}


# Function for constructing the following variables: unweighted, weighted, weighted dichotomized. Takes de master dataframe, a vector with RSIDS and a name for the genetic instrument.
UKB_ConstructScores <- function(data, rsids, ScoreName){
  # Takes a dataframe
  # Takes a list with RSIDS for the genetic score
  # Takes a name for the score.

  d <- data

  r_rsids <- paste("r", rsids, sep = "_")

  w_rsids <- paste('w', rsids, sep = "_")

  d[,paste0("u_",ScoreName)] <- select(d, r_rsids) %>% rowSums(na.rm=T) # Unweighted

  d[,paste0("w_",ScoreName)] <- select(d, all_of(w_rsids)) %>% rowSums(na.rm=T) # Weighted

  d[,paste0("dich_",ScoreName)] <- dicho(d[,paste0("w_",ScoreName)], dich.by = "median") # Dichotomized weighted score

  d[,paste0("udich_",ScoreName)] <- dicho(d[,paste0("u_",ScoreName)], dich.by = "median") # Dichotomized unweighted score

  return(d)
}

# Analysis function
UKB_MODELS_BINOMIAL <- function(datafr, linear_outcome, logistic_outcome, dependent_variable){ # Gives you linear and logistic regression models for certain defined outcomes and for a certain independent variable. Corrected for age and sex automatically.
  # Takes a dataframe
  # Takes a vector with linear outcomes and a vector with logistic outcomes
  # Takes a string with your dependent variables

  require(dplyr)
  require(psychometric)

  datafr <- data.frame(datafr)

  # Place to save all results
  reslist <- list()

  # Define dependent variables
  variables <- c(paste0(dependent_variable))

  # Linear regression part
  for (i in 1:length(linear_outcome)){

    # define outcome
    outcome <- linear_outcome[i]

    # Define formula
    f <- as.formula(paste(outcome, # Formula for in lm call
                          paste(variables, collapse = " + "),
                          sep = " ~ "))

    # Perform regression
    fit1 <- lm(data=datafr, formula = f)

    # Clean up results
    tidied <- tidy(fit1)
    glanced <- glance(fit1)

    total <- list("glanced" = glanced, "tidied" = tidied, "CI_r2" = CI.Rsq(glanced$r.squared, glanced$nobs, length(fit1$coefficients)))
    reslist[[paste0(linear_outcome[i])]] <- total

  }

  for (i in 1:length(logistic_outcome)){

    # define outcome
    outcome <- logistic_outcome[i]

    # Define formula
    f <- as.formula(paste(outcome, # Formula for in lm call
                          paste(variables, collapse = " + "),
                          sep = " ~ "))

    # Perform regression
    fit1 <- glm(data = datafr, formula = f, family = "quasibinomial")


    # cLean up results
    tidied <- tidy(fit1)
    tidied$es_exp <- exp(tidied$estimate)

    # additional statistics
    glanced <- glance(fit1)


    total <- list("glanced" = glanced, "tidied" = tidied)

    reslist[[paste0(logistic_outcome[i])]] <- total
  }
  return(reslist)
}

# EXAMPLE script ---------------
# Load  data

load(file = "<FILE>") # Assume that data has already been filtered for discordant self-reported and genetic sex, and for ethnicity, as well as for kinship (see ukbtools R package)

# Read details on all variants extracted from the .list file from BGEN
UKBann <- list()

for (i in 1:22) { # 22 chromosome list files, with SNPS from each chromosome

  setwd('<WHERE BGEN .LIST FILES ARE') # Where the list files are

  x<- read.delim(file = paste0("CHR", i, ".list"), header = T, sep = "\t", comment.char="#", stringsAsFactors = F, colClasses = "character")

  UKBann[[i]] <- x

}

datasetvariants <- do.call(rbind, UKBann) # Generate dataframe from lists with variants extracted, and their allelic annotation

# Load genetic scores constructed using 2 sample MR
SCORE_LDLC <- read.csv(file = "")%>%
  rename(rsid = SNP,
         BETA_EXP = beta,
         REF_EXP = EA,
         ALT_EXP = NEA)

# Prepare dataset for analysis: -------------

# Subset all rsids from the dataset to a vector

rsids <-  grep(colnames(df), pattern = "^rs", value = T)

# Round all imputed rsids to the nearest 1 number integer, so that we can flip them. (e.g. from homozygous for 1 allele to homozygous for the other allele).

for(i in 1:length(rsids)){

  df[,paste0(rsids[i])]  <- round(df[,paste0(rsids[i])])

}

# Calculate positive genetic scores (e.g. 1 unit increase in LDLC GRS)

res <- UKB_WeighVarPos(d=df,UKB_variants_list = datasetvariants, weights_list = SCORE_LDLC) # LDLC

df <- res$data

write.csv(res$variant_changes, "LDLC_variant_changes.csv") # Reference metadata for changes

# Stratify and create scores ------------
df <- UKB_ConstructScores(data = df, rsids = SCORE_LDLC$rsid, ScoreName = "LDLC")


# Regression analyses -------------

linear_outcome <-c("Linear_outcome1", "Logistic_outcome2")

logistic_outcome <- c("Logistic_outcome1", "Logistic_outcome2")

Res_wLDLC <- UKB_MODELS_BINOMIAL(datafr=df, linear_outcome = linear_outcome, logistic_outcome = logistic_outcome, dependent_variable = "w_LDLC+SEX+AGE") # WEighted LDL score, age and sex
save(Res_wLDLC, file="Res_wLDLC.rda")

Res_wLDLCint <- UKB_MODELS_BINOMIAL(datafr=df, linear_outcome = linear_outcome, logistic_outcome = logistic_outcome, dependent_variable = "w_LDLC*SEX+AGE")
save(Res_wLDLCint, file="Res_wLDLCint.rda")

# Leave one out -----------------

LOUres<- list()

for(i in 1:nrow(SCORE_LDLC)){

  out = SCORE_LDLC$rsid[i]

  SCORE_LDLC2 <- SCORE_LDLC[-i,]

  df <- UKB_ConstructScores(data = df, rsids = SCORE_LDLC2$rsid, ScoreName = "LDLC")

  o <- tidy(glm(data=df, formula = CVD ~ CURRENTAGE + SEX * w_LDLC, family = 'quasibinomial')) %>% subset(term != c("(Intercept)", "CURRENTAGE", "SEXmale"))

  o$rsid = out

  LOUres[[out]] <- o
}

LOUres <- do.call(rbind, LOUres)
write.csv(LOUres, "leave_one_out.csv")
