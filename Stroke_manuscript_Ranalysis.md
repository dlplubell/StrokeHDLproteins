---
title: "Statistical Analysis to accompany: "
subtitle: "High-density lipoprotein carries markers that track with recovery from stroke"
author: "Plubell D.L., et al., Circulation Research 2020"
date: "1/12/2020"
output:
    html_document:
      keep_md: true
      code_folding: hide
      toc: true
      toc_float: true
      theme: spacelab
      highlight: haddock
---
__Updated June 5, 2020 by DLP__


# Outline 

* __Abstract__
* __Demographic information__
    + _Table 1._ Non-stroke v Stroke demographics
* __HDL protein differential abundance testing__
    + Shapiro-Wilk test for normality of each protein
    + _Table 2._ Significant changes
    + _Figure 2._ Significantly changing HDL protein abundance distributions
    + _Figure 3._ Paired significantly changing HDL protein abundances
* __Correlation testing__
    + _Figure 4._ Correlation plot of all patient measures
* __Stroke recovery linear regression analysis__
    + _Table 3._ HDL protein logFC from 24h to 96 h post stroke with significant correlation with stroke recovery
* __Cholesterol efflux capacity analysis__
    + _Table 4._ Cholesterol efflux analysis
    + _Figure 5._ Cholesterol efflux distribution plot 
* __Supplemental figures/tables__
    + _Supplemental Figure 1._ Effect of normalization on N15 APOA1 spike-in variance distribution
    + _Supplemental Table 4._ Correlation results matrices
    + _Supplemental Figure 4._ Plots of HDL protein logFC correlation with stroke recovery
    + _Supplemental Table 5._ Linear regression statistics for HDL protein logFC and stroke recovery
    + _Supplemental Figure 5._ Plots of plasma lipid level correlation with stroke recovery 
    + _Supplemental Table 6._ Linear regression statistics for plasma lipid levels and stroke recovery
    + _Supplemental Table 7._ Linear regression statistics for plasma lipid levels and HDL protein logFC
    + _Supplemental Figure 6._ CEC distributions of stroke tPA status
    + _Supplemental Table 8._ T-test statistics for plasma lipid levels by stroke tPA status
* __Additional info for reviewers__
  
<br>

# 1. Abstract

<br>
Prospective cohort studies question the value of HDL-C for stroke risk prediction. We investigated the relationship between long-term functional recovery and changes in HDL protein composition and function (cholesterol efflux capacity, or CEC) in patients after acute ischemic stroke at two time points (24 h, 35 patients; 96 h, 20 patients) and in 35 control subjects. When compared to control subject after adjustments for sex and HDL-C levels, thirteen proteins some of which participate in acute phase response and platelet activation ( APMAP, GPLD1, APOE, IHH, ITIH4, SAA2, APOA4, CLU, SELL, ANTRX2, PON1, SERPINA1, and APOF) were significantly (adj. p<0.05) altered in stroke HDL at 96h. The first eight of these proteins were also significantly altered at 24h.  Consistent with inflammatory remodeling, CEC was reduced by 32% (P<0.001) at both time points. Baseline stroke severity adjusted regression model showed that changes within 96 hours post stroke in APOF, APOL1, APMAP, APOC4, APOM, PCYOX1, PON1, APOE, and PPBP correlate with stroke recovery scores (R2=0.37-0.75, adjusted p<0.05).  APOF (R2=0.74), and APOL1 (R2=0.60) continued to significantly correlate with recovery scores after accounting for tPA treatment. We conclude that changes in HDL associated proteins during early acute phase of stroke associate with recovery. 
<br>
<br>
This document contains the statistical analysis performed on general patient demographics, patient plasma lipid measures, patient cholesterol efflux capacity, stroke severity and recovery scores, and patient HDL proteomics data. Patient demographics and plasma measures are documented in Supplemental Table 1, and proteomic assay data are documented in Supplemental Table 1, and more thoroughly in Supplemental Table 3. All proteomics data is also deposited in Panorama Public with the ProteomeExchange identifier PXD015001. 
<br>

### Workspace set-up & Dataframe import 


```r
#Packages used in this analysis:
library(ggplot2, quietly = TRUE)
library(RColorBrewer, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(grid, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(pheatmap, quietly = TRUE)
library(psych, quietly = TRUE)
library(Hmisc, quietly = TRUE)
library(corrplot, quietly = TRUE)
library(mclust, quietly = TRUE)
library(purrr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(plyr, quietly = TRUE)
library(Gmisc, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(limma, quietly = TRUE)
library(scales, quietly=TRUE)
library(stargazer, quietly = TRUE)
library(knitr, quietly = TRUE)
library(kableExtra, quietly = TRUE)
library(statmod, quietly = TRUE)
library(xlsx, quietly=TRUE)
```



```r
#First read in the 'metadata' stored in supplemental Table 1, sheet 2
#Since the read.xlsx has issues with determining the class of each column they are specified
m <- as.data.frame(read.xlsx2("SupplementalTable_1.xlsx", 2, header=TRUE, 
                  colClasses =  c("Sample_no."="numeric", "Sex"="character",
                                        "Age_blood"="numeric","Group"="character", 
                                        "Stroke"="numeric", "Pair"="character", 
                                        "Pair_ID"="character", "Ischemic_Type"="character", 
                                        "tPa"="character", "Thrombectomy"="character", 
                                        "TICIB."="character", "Hour"="numeric", 
                                        "Cholesterol"="numeric", "Triglycerides"="numeric", 
                                        "HDLC"="numeric", "LDLC"="numeric", "Statin"="character",
                                        "NIHSS_baseline"="numeric", "NIHSS_3mo"="numeric", 
                                        "mRS"="numeric", "CEC"="numeric")))

#Now read in the protein assay data, stored in supplemental Table 1, sheet 3
p <- read.xlsx2("SupplementalTable_1.xlsx", 3, header=TRUE, 
                colClasses = c("numeric", "character", "character", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", "numeric"))

#for data wrangling I choose to combine everything into one large table
all <- as.data.frame(cbind(m, p))

#This removes a redundant column that appeared in both tables
#It also removes a sample that is a 96 hour sample with no paired 24 hour
#The 24 hour sample that corresponds was not analyzed by mass spec due to sample handling issues
all <- all[-1,-22]
#3 proteins analyzed by mass spec assay had to high of variance in system suitability replicates
#For specifics; see supplemental Figure 1.
#Therefore, I also remove these so they don't have a chance of being included in any analysis.
all$HP <- NULL
all$GC <- NULL
all$SAA4 <- NULL
```

<br>  


# 2. Demographic information

## Table 1: Non-stroke v Stroke sample demographics {.tabset}

Because the some of the demographic statistics are only applicable to certain groups, multiple tables are generated here and then compiled/merged into Table 1 in excel for the final manuscript table (for ease of formatting).


```r
#Creating subsets of the metadata for the 
meta_all <- as.data.frame(all[1:23])

#metadata dataframe for the 24h stroke samples
meta24 <- select(filter(meta_all, Group=="S24"), c(Sex, Hour, Age_blood, Group, 
                              Cholesterol, Triglycerides, HDLC, LDLC, Statin, 
                              NIHSS_baseline, NIHSS_3mo, mRS, CEC, Ischemic_Type, 
                              tPa, Thrombectomy, TICI2B)) 

#metadata dataframe for the non-stroke 'control' samples
metaC <- select(filter(meta_all, Group=="C"), c(Sex, Hour, Age_blood, Group, 
                               Cholesterol, Triglycerides, HDLC, LDLC, Statin, 
                               NIHSS_baseline, NIHSS_3mo, mRS, CEC, Ischemic_Type, 
                               tPa, Thrombectomy, TICI2B))

#metadata dataframe for the 96h stroke samples
meta96 <- select(filter(meta_all, Group=="S96"), c(Sex, Hour, Age_blood, Group, 
                              Cholesterol, Triglycerides, HDLC, LDLC, Statin, 
                              NIHSS_baseline, NIHSS_3mo, mRS, CEC, Ischemic_Type, 
                              tPa, Thrombectomy, TICI2B))

#metadata dataframe for both the 24h and 96h stroke samples
meta_stroke <- rbind(meta24, meta96)

#metadata dataframe for both the 24h and control samples
meta <- rbind(meta24, metaC)
```


### Non-stroke control vs. 24 hr statistics


```r
#Function to write out statistics: number of samples in each category & (percent of total) -or- mean & (sd)
getTable1Stats <- function(x, digits = 2, ...){
  getDescriptionStatsBy(x = x, 
                        by = meta$Group,
                        digits = digits,
                        continuous_fn = describeMean,
                        header_count = TRUE,
                        #statistics = TRUE,
                        ...)
  
}

#generating lists for each category of interest
t1 <- list()

t1[["Sex"]] <- getTable1Stats(meta$Sex)
t1[["Statin"]] <- getTable1Stats(meta$Statin)
t1[["Age"]] <- getTable1Stats(meta$Age_blood)
t1[["Cholesterol"]] <- getTable1Stats(meta$Cholesterol)
t1[["HDLC"]] <- getTable1Stats(meta$HDLC)
t1[["LDLC"]] <- getTable1Stats(meta$LDLC)
t1[["Triglycerides"]] <- getTable1Stats(meta$Triglycerides)

#merging statistics into a output table
mergeDesc(t1, 
          htmlTable_args = list(css.rgroup=""))
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>C<br />
 No. 35</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>S24<br />
 No. 35</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>S96<br />
 No. 0</th>
</tr>
</thead>
<tbody> 
<tr><td colspan='4' style=''>Sex</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;FEMALE</td>
<td style='text-align: center;'>26 (74.29%)</td>
<td style='text-align: center;'>14 (40.00%)</td>
<td style='text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;MALE</td>
<td style='text-align: center;'>8 (22.86%)</td>
<td style='text-align: center;'>21 (60.00%)</td>
<td style='text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='text-align: center;'>1 (2.86%)</td>
<td style='text-align: center;'>0 (0.00%)</td>
<td style='text-align: center;'>-</td>
</tr> 
<tr><td colspan='4' style=''>Statin</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='text-align: center;'>4 (11.43%)</td>
<td style='text-align: center;'>0 (0.00%)</td>
<td style='text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;No</td>
<td style='text-align: center;'>25 (71.43%)</td>
<td style='text-align: center;'>25 (71.43%)</td>
<td style='text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Yes</td>
<td style='text-align: center;'>6 (17.14%)</td>
<td style='text-align: center;'>10 (28.57%)</td>
<td style='text-align: center;'>-</td>
</tr> 
<tr><td colspan='4' style=''>Age</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Mean (SD)</td>
<td style='text-align: center;'>63.93 (&plusmn;3.42)</td>
<td style='text-align: center;'>68.20 (&plusmn;10.51)</td>
<td style='text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Missing</td>
<td style='text-align: center;'>4 (11.43%)</td>
<td style='text-align: center;'>0 (0%)</td>
<td style='text-align: center;'>-</td>
</tr>
<tr>
<td style='text-align: left;'>Cholesterol</td>
<td style='text-align: center;'>208.34 (&plusmn;43.98)</td>
<td style='text-align: center;'>148.83 (&plusmn;40.54)</td>
<td style='text-align: center;'></td>
</tr>
<tr>
<td style='text-align: left;'>HDLC</td>
<td style='text-align: center;'>63.86 (&plusmn;20.92)</td>
<td style='text-align: center;'>43.09 (&plusmn;15.99)</td>
<td style='text-align: center;'></td>
</tr>
<tr>
<td style='text-align: left;'>LDLC</td>
<td style='text-align: center;'>169.26 (&plusmn;52.22)</td>
<td style='text-align: center;'>128.43 (&plusmn;49.47)</td>
<td style='text-align: center;'></td>
</tr>
<tr>
<td style='border-bottom: 2px solid grey; text-align: left;'>Triglycerides</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>123.89 (&plusmn;86.49)</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>113.46 (&plusmn;80.64)</td>
<td style='border-bottom: 2px solid grey; text-align: center;'></td>
</tr>
</tbody>
</table>

### Stroke 24 hr vs. 96 hr stroke related statistics


```r
#Function to write out statistics: number of samples in each category & (percent of total) -or- mean & (sd)
#This time for just stroke samples
getTable1Stats <- function(x, digits = 2, ...){
  getDescriptionStatsBy(x = x, 
                        by = meta_stroke$Group,
                        digits = digits,
                        continuous_fn = describeMedian,
                        header_count = TRUE,
                        #statistics = TRUE,
                        ...)
}

#generating lists for each category of interest
t2 <- list()

t2[["Ischemic Type"]] <- getTable1Stats(meta_stroke$Ischemic_Type)
t2[["tPA"]] <- getTable1Stats(meta_stroke$tPa)
t2[["Thrombectomy"]] <- getTable1Stats(meta_stroke$Thrombectomy)
t2[["TICI2B"]] <- getTable1Stats(meta_stroke$TICI2B)
t2[["NIHSS-baseline"]] <- getTable1Stats(meta_stroke$NIHSS_baseline)
t2[["NIHSS-3mo"]] <- getTable1Stats(meta_stroke$NIHSS_3mo)
t2[["mRS"]] <- getTable1Stats(meta_stroke$mRS)

#merging statistics into a output table
mergeDesc(t2, 
          htmlTable_args = list(css.rgroup=""))
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>C<br />
 No. 0</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>S24<br />
 No. 35</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>S96<br />
 No. 20</th>
</tr>
</thead>
<tbody> 
<tr><td colspan='4' style=''>Ischemic Type</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;CE</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>8 (22.86%)</td>
<td style='text-align: center;'>4 (20.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;LV</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>22 (62.86%)</td>
<td style='text-align: center;'>14 (70.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>0 (0.00%)</td>
<td style='text-align: center;'>0 (0.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;SV</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>5 (14.29%)</td>
<td style='text-align: center;'>2 (10.00%)</td>
</tr> 
<tr><td colspan='4' style=''>tPA</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>0 (0.00%)</td>
<td style='text-align: center;'>0 (0.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;No</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>19 (54.29%)</td>
<td style='text-align: center;'>11 (55.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Yes</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>16 (45.71%)</td>
<td style='text-align: center;'>9 (45.00%)</td>
</tr> 
<tr><td colspan='4' style=''>Thrombectomy</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>0 (0.00%)</td>
<td style='text-align: center;'>0 (0.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;No</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>20 (57.14%)</td>
<td style='text-align: center;'>12 (60.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Yes</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>15 (42.86%)</td>
<td style='text-align: center;'>8 (40.00%)</td>
</tr> 
<tr><td colspan='4' style=''>TICI2B</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;NA</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>20 (57.14%)</td>
<td style='text-align: center;'>12 (60.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;No</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>5 (14.29%)</td>
<td style='text-align: center;'>5 (25.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Yes</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>10 (28.57%)</td>
<td style='text-align: center;'>3 (15.00%)</td>
</tr>
<tr>
<td style='text-align: left;'>NIHSS-baseline</td>
<td style='text-align: center;'></td>
<td style='text-align: center;'>15.00 (9.00 - 20.00)</td>
<td style='text-align: center;'>17.00 (12.25 - 20.00)</td>
</tr> 
<tr><td colspan='4' style=''>NIHSS-3mo</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Median (IQR)</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>6.00 (1.00 - 13.00)</td>
<td style='text-align: center;'> ( - )</td>
</tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Missing</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>2 (5.71%)</td>
<td style='text-align: center;'>20 (100.00%)</td>
</tr> 
<tr><td colspan='4' style=''>mRS</td></tr>
<tr>
<td style='text-align: left;'>&nbsp;&nbsp;Median (IQR)</td>
<td style='text-align: center;'>-</td>
<td style='text-align: center;'>3.00 (2.00 - 4.00)</td>
<td style='text-align: center;'> ( - )</td>
</tr>
<tr>
<td style='border-bottom: 2px solid grey; text-align: left;'>&nbsp;&nbsp;Missing</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>-</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>0 (0%)</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>20 (100.00%)</td>
</tr>
</tbody>
</table>

### All 3 groups


```r
#Function to write out statistics: number of samples in each category 
# & (percent of total) -or- mean & (sd)

#This time for statistics applicable to all 3 groups: specifically for lipid measures
getTable1Stats <- function(x, digits = 2, ...){
  getDescriptionStatsBy(x = x, 
                        by = meta_all$Group,
                        digits = digits,
                        continuous_fn = describeMean,
                        header_count = TRUE,
                        #statistics = TRUE,
                        ...)
  
}

#generating lists for each category of interest
t3 <- list()

t3[["Cholesterol"]] <- getTable1Stats(meta_all$Cholesterol)
t3[["HDLC"]] <- getTable1Stats(meta_all$HDLC)
t3[["LDLC"]] <- getTable1Stats(meta_all$LDLC)
t3[["Triglycerides"]] <- getTable1Stats(meta_all$Triglycerides)

#merging statistics into a output table
mergeDesc(t3, 
          htmlTable_args = list(css.rgroup=""))
```

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>C<br />
 No. 35</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>S24<br />
 No. 35</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>S96<br />
 No. 20</th>
</tr>
</thead>
<tbody>
<tr>
<td style='text-align: left;'>Cholesterol</td>
<td style='text-align: center;'>208.34 (&plusmn;43.98)</td>
<td style='text-align: center;'>148.83 (&plusmn;40.54)</td>
<td style='text-align: center;'>140.65 (&plusmn;40.83)</td>
</tr>
<tr>
<td style='text-align: left;'>HDLC</td>
<td style='text-align: center;'>63.86 (&plusmn;20.92)</td>
<td style='text-align: center;'>43.09 (&plusmn;15.99)</td>
<td style='text-align: center;'>38.80 (&plusmn;9.30)</td>
</tr>
<tr>
<td style='text-align: left;'>LDLC</td>
<td style='text-align: center;'>169.26 (&plusmn;52.22)</td>
<td style='text-align: center;'>128.43 (&plusmn;49.47)</td>
<td style='text-align: center;'>124.84 (&plusmn;55.22)</td>
</tr>
<tr>
<td style='border-bottom: 2px solid grey; text-align: left;'>Triglycerides</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>123.89 (&plusmn;86.49)</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>113.46 (&plusmn;80.64)</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>114.95 (&plusmn;75.76)</td>
</tr>
</tbody>
</table>


# 3. HDL protein differential abundance testing


```r
#Since there are 3 groups, and 2 of the groups are paired, 
#this makes data a bit more complicated to handle, so I just generated 
#new dataframes instead of having to wrangle the data with other methods. 

## Paired stroke sample dataframe generation
PS24 <- all %>% filter(Pair=="Yes" & Group=="S24") %>% droplevels()
PS24rownames <- PS24[,1]
rownames(PS24) <- PS24rownames

PS96 <- all %>% filter(Pair=="Yes" & Group=="S96") %>% droplevels()
PS96rownames <- PS96[,1]
rownames(PS96) <- PS96rownames

paired <- rbind(PS24, PS96)

#The next 3 sets are all the 24h vs all the 96h measures, 
#all the control vs all the 24h measures, and all the control 
# vs all the 96 h measures.
C <- all %>% filter(Group=="C") %>% droplevels()
Crownames <- C[,1]
rownames(C) <- Crownames

S24 <- all %>% filter(Group=="S24") %>% droplevels()
S24rownames <- S24[,1]
rownames(S24) <- S24rownames

#dataframes for significance testing
CvS24 <- rbind(C,S24)
CvS96 <- rbind(C,PS96)
all <- rbind(C,S24,PS96)


#log2 transforming the paired stroke dataframe
logP <- as.data.frame(cbind(paired[,1:23], log2(paired[,24:ncol(paired)])))

#log2 transforming the control & stroke dataframe
logAll <- as.data.frame(cbind(all[,1:23], log2(all[,24:ncol(all)])))

#log2 transforming the control & stroke dataframe
logCvS24 <- as.data.frame(cbind(CvS24[,1:23], log2(CvS24[,24:ncol(CvS24)])))

#log2 transforming the control & stroke dataframe
logCvS96 <- as.data.frame(cbind(CvS96[,1:23], log2(CvS96[,24:ncol(CvS96)])))
```
<br>
<br>


## Shapiro-Wilk test for normality of each protein

A caveat to targeted mass spectrometry measurements is that the data isn't necessarily normally distributed - in part due to fundamentally how the measurements are made. However, within the proteomics community it is standard to use parametric tests for differential abundance testing. In this assay we observe several cases in which there are outlier low abundance measurments below the otherwise fairly normal distribution. This is due to integrating the same transition ions across all samples, including ones lacking of signal. In that case we are likely integrating background noise. Another cause of these outlier measurements could be due to them falling below the limit of quantification or limit of detection. In that case we are integrating the actual signal, but it is lower than expected for a true linear scale. 

<br>  
__Shapiro Wilk statistics for non-stroke control and 24 h post stroke samples__

```r
#Shapiro-Wilk test
#perform the test for each protein measure
p_swnorm <- lapply(logCvS24[, 24:length(logCvS24)], 
                   function(x) shapiro.test(x))
#retrieve the W statistic
results_w <- ldply(p_swnorm, function(x) x$statistic)
#retrieve the p-value
results_p <- ldply(p_swnorm, function(x) x$p.value)

#Report the mean and median p-value, the number and 
#percent with p-values greater than 0.05:
mean_pC24 <- format(mean(results_p$V1), digits=3, format="f")
median_pC24 <- format(median(results_p$V1), digits=2, format="f")
count_C24 <- with(results_p, sum(V1>0.05))
count_tC24 <- length(results_p$V1)
pct_normC24 <- format((count_C24/count_tC24*100), 
                      digits=3, format="f")

#Making a table of the statistics
results_swnorm <- data.frame(cbind("protein" = results_w$.id, 
                                   "W" = results_w$W, 
                                   "Pvalue" = results_p$V1))
kable(results_swnorm) %>%
  kable_styling() %>%
  scroll_box(width = "75%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:75%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> protein </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> W </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:left;"> 0.972916707058013 </td>
   <td style="text-align:left;"> 0.133213707321529 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:left;"> 0.850998646299485 </td>
   <td style="text-align:left;"> 7.20893480131045e-07 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:left;"> 0.935325757522222 </td>
   <td style="text-align:left;"> 0.00133132452005188 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:left;"> 0.968807268942514 </td>
   <td style="text-align:left;"> 0.0780177614688162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:left;"> 0.970876723864567 </td>
   <td style="text-align:left;"> 0.102175968800515 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:left;"> 0.986633591994339 </td>
   <td style="text-align:left;"> 0.754438162438065 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:left;"> 0.966861793121548 </td>
   <td style="text-align:left;"> 0.0605587049879609 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:left;"> 0.878748650952478 </td>
   <td style="text-align:left;"> 6.18471908299661e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:left;"> 0.975551405489702 </td>
   <td style="text-align:left;"> 0.187064049284491 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:left;"> 0.978181468611017 </td>
   <td style="text-align:left;"> 0.260674565194396 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:left;"> 0.98108244092458 </td>
   <td style="text-align:left;"> 0.370316247905908 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:left;"> 0.983055292095996 </td>
   <td style="text-align:left;"> 0.463542070367344 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:left;"> 0.974707805204013 </td>
   <td style="text-align:left;"> 0.167889885291914 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:left;"> 0.908991700425532 </td>
   <td style="text-align:left;"> 8.94663404061348e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:left;"> 0.980089463682507 </td>
   <td style="text-align:left;"> 0.329148277379397 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:left;"> 0.967959550511187 </td>
   <td style="text-align:left;"> 0.0698589503975613 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:left;"> 0.881243067503217 </td>
   <td style="text-align:left;"> 7.60003254224941e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:left;"> 0.993093618705489 </td>
   <td style="text-align:left;"> 0.968573463845914 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:left;"> 0.989880898377397 </td>
   <td style="text-align:left;"> 0.848393987787419 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:left;"> 0.970298332186325 </td>
   <td style="text-align:left;"> 0.0947575551367428 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:left;"> 0.953739023140204 </td>
   <td style="text-align:left;"> 0.0114355859957422 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:left;"> 0.952059647667699 </td>
   <td style="text-align:left;"> 0.00930820804575936 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:left;"> 0.952760309878095 </td>
   <td style="text-align:left;"> 0.0108997857304548 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:left;"> 0.989913504133544 </td>
   <td style="text-align:left;"> 0.850036073390094 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:left;"> 0.95268120820667 </td>
   <td style="text-align:left;"> 0.0100427940148696 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:left;"> 0.917645644716561 </td>
   <td style="text-align:left;"> 0.000231324115819709 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:left;"> 0.825334655916189 </td>
   <td style="text-align:left;"> 1.20397926228046e-07 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:left;"> 0.883536601747945 </td>
   <td style="text-align:left;"> 9.20459626007136e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:left;"> 0.899356235480224 </td>
   <td style="text-align:left;"> 3.65616811799498e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:left;"> 0.93379876957411 </td>
   <td style="text-align:left;"> 0.00112543571708825 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:left;"> 0.978244955701556 </td>
   <td style="text-align:left;"> 0.262736950904108 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:left;"> 0.974814373525079 </td>
   <td style="text-align:left;"> 0.170205061759006 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:left;"> 0.99022206835699 </td>
   <td style="text-align:left;"> 0.865236511357912 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:left;"> 0.694357614558252 </td>
   <td style="text-align:left;"> 8.76703025795617e-11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:left;"> 0.9892816275195 </td>
   <td style="text-align:left;"> 0.817122683724782 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:left;"> 0.964670825272416 </td>
   <td style="text-align:left;"> 0.053904162082101 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:left;"> 0.981337845022924 </td>
   <td style="text-align:left;"> 0.391332926963427 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:left;"> 0.982366809236247 </td>
   <td style="text-align:left;"> 0.429312975758955 </td>
  </tr>
</tbody>
</table></div>
<br>  
For the control and 24 h post stroke protein measures __60.5__ percent had p-values above 0.05; the mean p-value was __0.231__ and median p-value was __0.086__. 4 of the 8 proteins found to be significantly different below are normally distributed.

<br> 

__Shapiro Wilk statistics for non-stroke control and 96 h post stroke samples__


```r
#Shapiro-Wilk test
p_swnorm <- lapply(logCvS96[, 24:length(logCvS96)], 
                   function(x) shapiro.test(x))
results_w <- ldply(p_swnorm, function(x) x$statistic)
results_p <- ldply(p_swnorm, function(x) x$p.value)

#Report the mean and median p-value:
mean_pC96 <-format(mean(results_p$V1), digits=3, format="f")
median_pC96 <-format(median(results_p$V1), digits=2, format="f")
count_C96 <- with(results_p, sum(V1>0.05))
count_tC96 <- length(results_p$V1)
pct_normC96 <- format((count_C96/count_tC96*100), 
                      digits=3, format="f")

#Making a table of the statistics
results_swnorm <- data.frame(cbind("protein" = results_w$.id, 
                                   "W" = results_w$W, 
                                   "Pvalue" = results_p$V1))
kable(results_swnorm) %>%
  kable_styling() %>%
  scroll_box(width = "75%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:75%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> protein </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> W </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:left;"> 0.984569033978003 </td>
   <td style="text-align:left;"> 0.700515710229085 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:left;"> 0.871162761777608 </td>
   <td style="text-align:left;"> 2.84499807010518e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:left;"> 0.93671351307853 </td>
   <td style="text-align:left;"> 0.00619621123103363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:left;"> 0.927106673932437 </td>
   <td style="text-align:left;"> 0.00252479404643427 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:left;"> 0.972724793226624 </td>
   <td style="text-align:left;"> 0.243314936794797 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:left;"> 0.958085649210766 </td>
   <td style="text-align:left;"> 0.0739063838904789 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:left;"> 0.959592435332756 </td>
   <td style="text-align:left;"> 0.0619525845983186 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:left;"> 0.86368113016257 </td>
   <td style="text-align:left;"> 1.68896163332746e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:left;"> 0.982913932682701 </td>
   <td style="text-align:left;"> 0.62072221785666 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:left;"> 0.981518073624968 </td>
   <td style="text-align:left;"> 0.555260493891661 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:left;"> 0.961330687256626 </td>
   <td style="text-align:left;"> 0.074329019460728 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:left;"> 0.961942123589978 </td>
   <td style="text-align:left;"> 0.0792540877486649 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:left;"> 0.965623627914927 </td>
   <td style="text-align:left;"> 0.116646508331159 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:left;"> 0.887365349950404 </td>
   <td style="text-align:left;"> 9.30858666368004e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:left;"> 0.984419975519016 </td>
   <td style="text-align:left;"> 0.69331456484689 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:left;"> 0.975032829468435 </td>
   <td style="text-align:left;"> 0.306406531222257 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:left;"> 0.823234368960334 </td>
   <td style="text-align:left;"> 1.27812028164465e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:left;"> 0.986072595420624 </td>
   <td style="text-align:left;"> 0.772088716124038 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:left;"> 0.98459415894536 </td>
   <td style="text-align:left;"> 0.701728730219006 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:left;"> 0.9890985235947 </td>
   <td style="text-align:left;"> 0.898150407006441 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:left;"> 0.958160040054358 </td>
   <td style="text-align:left;"> 0.0533403825228751 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:left;"> 0.955062780596469 </td>
   <td style="text-align:left;"> 0.0386603198961599 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:left;"> 0.977701834853351 </td>
   <td style="text-align:left;"> 0.40822626697194 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:left;"> 0.968838333279901 </td>
   <td style="text-align:left;"> 0.16320729244716 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:left;"> 0.914715259525487 </td>
   <td style="text-align:left;"> 0.000843152168996652 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:left;"> 0.960008912913546 </td>
   <td style="text-align:left;"> 0.0647131502391475 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:left;"> 0.760527327694844 </td>
   <td style="text-align:left;"> 4.28442460929897e-08 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:left;"> 0.852223677580451 </td>
   <td style="text-align:left;"> 7.82121161608235e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:left;"> 0.880961327104868 </td>
   <td style="text-align:left;"> 5.77062895861446e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:left;"> 0.94776037086204 </td>
   <td style="text-align:left;"> 0.0183371995372958 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:left;"> 0.984287981376001 </td>
   <td style="text-align:left;"> 0.6869328981246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:left;"> 0.964039450481194 </td>
   <td style="text-align:left;"> 0.0987791476894501 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:left;"> 0.977416291459016 </td>
   <td style="text-align:left;"> 0.385568189874355 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:left;"> 0.6561175382853 </td>
   <td style="text-align:left;"> 4.37289426496688e-10 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:left;"> 0.983148600814531 </td>
   <td style="text-align:left;"> 0.631955541742189 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:left;"> 0.965448914989642 </td>
   <td style="text-align:left;"> 0.127703761482916 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:left;"> 0.945152322014669 </td>
   <td style="text-align:left;"> 0.0153750163536789 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:left;"> 0.987523323617749 </td>
   <td style="text-align:left;"> 0.836862065642754 </td>
  </tr>
</tbody>
</table></div>

<br>  
For the control and 96 h post stroke protein measures __63.2__ percent had p-values above 0.05; the mean p-value was __0.248__ and median p-value was __0.077__. 7 of the 12 proteins found to be significantly different below are normally distributed.

<br>

__Shapiro Wilk statistics for paired 24 h and 96 h post stroke samples__


```r
#Shapiro-Wilk test
p_swnorm <- lapply(logP[, 24:length(logP)], function(x) shapiro.test(x))
results_w <- ldply(p_swnorm, function(x) x$statistic)
results_p <- ldply(p_swnorm, function(x) x$p.value)

#Report the mean and median p-value:
mean_pP <- format(mean(results_p$V1), digits=3, format="f")
median_pP <- format(median(results_p$V1), digits=2, format="f")
count_P <- with(results_p, sum(V1>0.05))
count_tP <- length(results_p$V1)
pct_normP <- format((count_P/count_tP*100), 
                    digits=3, format="f")

#Making a table of the statistics
results_swnorm <- data.frame(cbind("protein" = results_w$.id, 
                                   "W" = results_w$W, 
                                   "Pvalue" = results_p$V1))
kable(results_swnorm) %>%
  kable_styling() %>%
  scroll_box(width = "75%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:75%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> protein </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> W </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:left;"> 0.985655513950479 </td>
   <td style="text-align:left;"> 0.88376669169845 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:left;"> 0.949061016965588 </td>
   <td style="text-align:left;"> 0.0703965238725716 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:left;"> 0.982412442338416 </td>
   <td style="text-align:left;"> 0.777891065788042 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:left;"> 0.96655309054376 </td>
   <td style="text-align:left;"> 0.278409430775947 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:left;"> 0.959042522622838 </td>
   <td style="text-align:left;"> 0.155320204347786 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:left;"> 0.940227246166859 </td>
   <td style="text-align:left;"> 0.0759866958040728 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:left;"> 0.968803404779125 </td>
   <td style="text-align:left;"> 0.329607125325641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:left;"> 0.976725030224524 </td>
   <td style="text-align:left;"> 0.569638934504142 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:left;"> 0.968927408519029 </td>
   <td style="text-align:left;"> 0.332646645258803 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:left;"> 0.976097830529137 </td>
   <td style="text-align:left;"> 0.547559021580809 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:left;"> 0.957266586049406 </td>
   <td style="text-align:left;"> 0.134968937619974 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:left;"> 0.98079990785719 </td>
   <td style="text-align:left;"> 0.719205168560376 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:left;"> 0.984753422457826 </td>
   <td style="text-align:left;"> 0.856715992878622 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:left;"> 0.976028923939916 </td>
   <td style="text-align:left;"> 0.545158543686624 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:left;"> 0.981214248622732 </td>
   <td style="text-align:left;"> 0.734471449922815 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:left;"> 0.97123883036428 </td>
   <td style="text-align:left;"> 0.39360169618946 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:left;"> 0.975304793085344 </td>
   <td style="text-align:left;"> 0.520260704319982 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:left;"> 0.97374234269091 </td>
   <td style="text-align:left;"> 0.468796289112221 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:left;"> 0.978875404955038 </td>
   <td style="text-align:left;"> 0.647826880807906 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:left;"> 0.985385394309606 </td>
   <td style="text-align:left;"> 0.875907501723648 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:left;"> 0.885772968681251 </td>
   <td style="text-align:left;"> 0.000756180407918371 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:left;"> 0.974321924173145 </td>
   <td style="text-align:left;"> 0.487505257011207 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:left;"> 0.96743073088793 </td>
   <td style="text-align:left;"> 0.297496085411599 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:left;"> 0.975527191264103 </td>
   <td style="text-align:left;"> 0.527841772832523 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:left;"> 0.929742693542839 </td>
   <td style="text-align:left;"> 0.0157975660341636 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:left;"> 0.97324920002216 </td>
   <td style="text-align:left;"> 0.469560997793684 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:left;"> 0.957805133731255 </td>
   <td style="text-align:left;"> 0.140848283987468 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:left;"> 0.970892131760639 </td>
   <td style="text-align:left;"> 0.383935845603567 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:left;"> 0.972378698303703 </td>
   <td style="text-align:left;"> 0.426676164534604 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:left;"> 0.812850330118147 </td>
   <td style="text-align:left;"> 1.26950046085226e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:left;"> 0.921514171825052 </td>
   <td style="text-align:left;"> 0.00859613598806194 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:left;"> 0.97124661299873 </td>
   <td style="text-align:left;"> 0.393820792204484 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:left;"> 0.978280053178659 </td>
   <td style="text-align:left;"> 0.625883973267704 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:left;"> 0.978203202747782 </td>
   <td style="text-align:left;"> 0.623064232675991 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:left;"> 0.966044304244947 </td>
   <td style="text-align:left;"> 0.267844802456448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:left;"> 0.874501048733442 </td>
   <td style="text-align:left;"> 0.00062049961633631 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:left;"> 0.93071060539253 </td>
   <td style="text-align:left;"> 0.0169902946526589 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:left;"> 0.974311722550804 </td>
   <td style="text-align:left;"> 0.487171947719173 </td>
  </tr>
</tbody>
</table></div>

<br>  
For the 24 h and 96 h post stroke protein measures __84.2__ percent had p-values above 0.05; the mean p-value was __0.397__ and median p-value was __0.41__. Both proteins found to be different by paired t-test analysis are normally distributed (APOF and ANTXR2).

<br>
<br>

## Table 2. Significant HDL protein abundance changes

### Paired t test for 24h - 96 h paired samples


```r
#Since a paired t.test will not accomodate missing values, several proteins will be removed from this analysis
non_missing <- as.data.frame(select(logP, -"LPA", -"ITIH4", -"SAA1"))

p_t.test <- lapply(non_missing[, 24:length(non_missing)], 
                   function(x) t.test(x ~ non_missing$Group, paired = TRUE))
results_p_t.test <- ldply(p_t.test, function(x) x$p.value)

results_p_t.test$adj.P.Val <- p.adjust(results_p_t.test$V1, 
                                       method="BH")
results_p_t.test <-rename(results_p_t.test, c(".id"= "Protein", 
                                              "V1"="P.Value"))

count_pair <- with(results_p_t.test, sum(adj.P.Val<0.05))

#write.csv(results_p_t.test, file="paired_ttest.csv")

kable(results_p_t.test, 
      caption="Paired t-test statistics for 24 h and 96 h post stroke samples") %>%
  kable_styling() %>%
  scroll_box(width = "75%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:75%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Paired t-test statistics for 24 h and 96 h post stroke samples</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Protein </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> P.Value </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> adj.P.Val </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.0232465 </td>
   <td style="text-align:right;"> 0.1627257 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.3853037 </td>
   <td style="text-align:right;"> 0.6714833 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.2534772 </td>
   <td style="text-align:right;"> 0.5849641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 0.0005329 </td>
   <td style="text-align:right;"> 0.0149537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 0.5340493 </td>
   <td style="text-align:right;"> 0.7499932 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.9566294 </td>
   <td style="text-align:right;"> 0.9626844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 0.3369339 </td>
   <td style="text-align:right;"> 0.6551493 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 0.6299003 </td>
   <td style="text-align:right;"> 0.8479427 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.8505786 </td>
   <td style="text-align:right;"> 0.9540000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 0.4837164 </td>
   <td style="text-align:right;"> 0.7370082 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.4843197 </td>
   <td style="text-align:right;"> 0.7370082 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.3341036 </td>
   <td style="text-align:right;"> 0.6551493 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.5357094 </td>
   <td style="text-align:right;"> 0.7499932 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.0008545 </td>
   <td style="text-align:right;"> 0.0149537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 0.1962744 </td>
   <td style="text-align:right;"> 0.5569666 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.9626844 </td>
   <td style="text-align:right;"> 0.9626844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.1863001 </td>
   <td style="text-align:right;"> 0.5569666 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.8722285 </td>
   <td style="text-align:right;"> 0.9540000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0296958 </td>
   <td style="text-align:right;"> 0.1732254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.0221582 </td>
   <td style="text-align:right;"> 0.1627257 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.8170512 </td>
   <td style="text-align:right;"> 0.9540000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.2674122 </td>
   <td style="text-align:right;"> 0.5849641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 0.8613524 </td>
   <td style="text-align:right;"> 0.9540000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 0.1788193 </td>
   <td style="text-align:right;"> 0.5569666 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0625167 </td>
   <td style="text-align:right;"> 0.2735107 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.6616746 </td>
   <td style="text-align:right;"> 0.8577263 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 0.2563235 </td>
   <td style="text-align:right;"> 0.5849641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 0.9499223 </td>
   <td style="text-align:right;"> 0.9626844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.0771167 </td>
   <td style="text-align:right;"> 0.2998981 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 0.2068733 </td>
   <td style="text-align:right;"> 0.5569666 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.0348365 </td>
   <td style="text-align:right;"> 0.1741825 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.3667914 </td>
   <td style="text-align:right;"> 0.6714833 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.4028900 </td>
   <td style="text-align:right;"> 0.6714833 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.0027379 </td>
   <td style="text-align:right;"> 0.0319420 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.7591141 </td>
   <td style="text-align:right;"> 0.9488926 </td>
  </tr>
</tbody>
</table></div>
<br>  
From paired t.test analysis 3 proteins are significantly different between the two time points.
<br>  
<br>

### Extra: Paired Wilcoxon ranked sum test for 24h - 96 h paired samples

At least one of the significantly changing proteins (SAA2) is non-normally distributed by the Shapiro-Wilk test, so the paired measurements are also assessed by wilcoxon signed rank test.


```r
wilcox <- lapply(non_missing[, 24:length(non_missing)], 
                   function(x) pairwise.wilcox.test(x, non_missing$Group, 
                                          p.adjust ="BH", paired = TRUE))

results_wilcox <- ldply(wilcox, function(x) x$p.value)
results_wilcox$adj.P.Val <- p.adjust(results_wilcox$S24, method="BH")
names(results_wilcox)[1] <- "protein"
names(results_wilcox)[2] <- "p.value"

count_pair <- with(results_wilcox, sum(adj.P.Val<0.05))

kable(results_wilcox, 
      caption="Wilcox Paired t-test statistics for 24 h and 96 h post stroke samples") %>%
  kable_styling() %>%
  scroll_box(width = "75%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:75%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Wilcox Paired t-test statistics for 24 h and 96 h post stroke samples</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> protein </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.value </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> adj.P.Val </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.0362339 </td>
   <td style="text-align:right;"> 0.1585233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.6215134 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.2773552 </td>
   <td style="text-align:right;"> 0.5710254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 0.0007076 </td>
   <td style="text-align:right;"> 0.0211620 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 0.7561665 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.9854355 </td>
   <td style="text-align:right;"> 0.9854355 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 0.2773552 </td>
   <td style="text-align:right;"> 0.5710254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 0.6742229 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.6476555 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 0.2454872 </td>
   <td style="text-align:right;"> 0.5710254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.7561665 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.1230927 </td>
   <td style="text-align:right;"> 0.3590202 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.3117943 </td>
   <td style="text-align:right;"> 0.5743579 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.0012093 </td>
   <td style="text-align:right;"> 0.0211620 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 0.2610989 </td>
   <td style="text-align:right;"> 0.5710254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.8694878 </td>
   <td style="text-align:right;"> 0.9221840 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.2942524 </td>
   <td style="text-align:right;"> 0.5721574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.7841263 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0327682 </td>
   <td style="text-align:right;"> 0.1585233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.0214844 </td>
   <td style="text-align:right;"> 0.1585233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.8123550 </td>
   <td style="text-align:right;"> 0.8885133 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.4523754 </td>
   <td style="text-align:right;"> 0.7196882 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 0.7285061 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 0.0758514 </td>
   <td style="text-align:right;"> 0.2654800 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0582581 </td>
   <td style="text-align:right;"> 0.2265591 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.7841263 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 0.2305126 </td>
   <td style="text-align:right;"> 0.5710254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 0.9563293 </td>
   <td style="text-align:right;"> 0.9844567 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.0362339 </td>
   <td style="text-align:right;"> 0.1585233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 0.1230927 </td>
   <td style="text-align:right;"> 0.3590202 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.0239506 </td>
   <td style="text-align:right;"> 0.1585233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.4523754 </td>
   <td style="text-align:right;"> 0.7196882 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.4090977 </td>
   <td style="text-align:right;"> 0.7159209 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.0031528 </td>
   <td style="text-align:right;"> 0.0367832 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.7561665 </td>
   <td style="text-align:right;"> 0.8853039 </td>
  </tr>
</tbody>
</table></div>
<br>  
Again, we find the same 3 proteins to be significantly different (adj. p-value <0.05) by paired Wilcoxon.

<br>
<br>

### Linear model for each protein

#### Model: ~0 + Group + Sex + HDLC {.tabset}

__Linear model statistics for 24h vs non-stroke control protein differences__




```r
all1 <- all
rownames(all1) <- all1$ID

#Removing sample that has NA indicated for sex
all1 <- all1[-9,]
all_mat <- as.matrix(log2(all1[, 24:length(all1)]))
rownames(all_mat) <- all1$ID
tmat <- t(all_mat)

meta <- as.data.frame(all1[, 1:23])
rownames(meta) <- all1$ID

#Model for finding differences due to group with sex and HDLC covariates
design <- model.matrix(~0 + Group + Sex + HDLC, data = meta)
fit <- lmFit(tmat, design = design)

#GroupS24 - GroupS96 contrast not made by linear model due to their paired nature.
contrasts <- makeContrasts(GroupS24 - GroupC, GroupS96 - GroupC, 
                           levels = design)
fit1 <- contrasts.fit(fit, contrasts)
fit1 <- eBayes(fit1,robust = TRUE,trend = TRUE)
out1 <- lapply(1:ncol(contrasts), function(x) topTable(fit1, coef = x, 
                                                       n= Inf, 
                                                       adjust = "BH", 
                                                       sort.by = "p")[,-c(3,6)])

#Number of significant differences
num <- sapply(out1,function(x) sum(x$adj.P<0.05))

# C vs S24 results
S24vC <- out1[[1]]
count_24vC <- with(S24vC, sum(adj.P.Val<0.05))
#write.csv(S24vC, file="S24vC.csv")

# C vs S96 results
S96vC <- out1[[2]]
count_96vC <- with(S96vC, sum(adj.P.Val<0.05))
#write.csv(S96vC, file="S96vC.csv")

#Table of the C v 24 h post stroke results
kable(S24vC, caption="Model: ~0 + Group + Sex + HDLC") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Model: ~0 + Group + Sex + HDLC</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> logFC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AveExpr </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> P.Value </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> adj.P.Val </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 1.6105734 </td>
   <td style="text-align:right;"> 21.52132 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 1.4695627 </td>
   <td style="text-align:right;"> 20.78969 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.7574868 </td>
   <td style="text-align:right;"> 28.56741 </td>
   <td style="text-align:right;"> 0.0000277 </td>
   <td style="text-align:right;"> 0.0003511 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.5885446 </td>
   <td style="text-align:right;"> 23.13573 </td>
   <td style="text-align:right;"> 0.0001068 </td>
   <td style="text-align:right;"> 0.0010148 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -1.1210225 </td>
   <td style="text-align:right;"> 19.79977 </td>
   <td style="text-align:right;"> 0.0003313 </td>
   <td style="text-align:right;"> 0.0025176 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 2.4818640 </td>
   <td style="text-align:right;"> 22.92835 </td>
   <td style="text-align:right;"> 0.0015677 </td>
   <td style="text-align:right;"> 0.0099288 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -0.5554768 </td>
   <td style="text-align:right;"> 23.61617 </td>
   <td style="text-align:right;"> 0.0051057 </td>
   <td style="text-align:right;"> 0.0277169 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.3581448 </td>
   <td style="text-align:right;"> 24.67599 </td>
   <td style="text-align:right;"> 0.0087013 </td>
   <td style="text-align:right;"> 0.0413312 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.6404597 </td>
   <td style="text-align:right;"> 20.50256 </td>
   <td style="text-align:right;"> 0.0122721 </td>
   <td style="text-align:right;"> 0.0518156 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 0.3195771 </td>
   <td style="text-align:right;"> 26.52012 </td>
   <td style="text-align:right;"> 0.0408526 </td>
   <td style="text-align:right;"> 0.1552399 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 0.2862072 </td>
   <td style="text-align:right;"> 24.35341 </td>
   <td style="text-align:right;"> 0.0560477 </td>
   <td style="text-align:right;"> 0.1936194 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.4327166 </td>
   <td style="text-align:right;"> 21.77589 </td>
   <td style="text-align:right;"> 0.0712996 </td>
   <td style="text-align:right;"> 0.2257821 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 1.5999117 </td>
   <td style="text-align:right;"> 25.05205 </td>
   <td style="text-align:right;"> 0.0887907 </td>
   <td style="text-align:right;"> 0.2532865 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> -0.2924746 </td>
   <td style="text-align:right;"> 26.82370 </td>
   <td style="text-align:right;"> 0.0959379 </td>
   <td style="text-align:right;"> 0.2532865 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.1537845 </td>
   <td style="text-align:right;"> 29.39118 </td>
   <td style="text-align:right;"> 0.1053507 </td>
   <td style="text-align:right;"> 0.2532865 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.2501451 </td>
   <td style="text-align:right;"> 21.78068 </td>
   <td style="text-align:right;"> 0.1066470 </td>
   <td style="text-align:right;"> 0.2532865 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -0.2130811 </td>
   <td style="text-align:right;"> 23.30947 </td>
   <td style="text-align:right;"> 0.1724332 </td>
   <td style="text-align:right;"> 0.3854390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> -0.5020297 </td>
   <td style="text-align:right;"> 22.70342 </td>
   <td style="text-align:right;"> 0.1910258 </td>
   <td style="text-align:right;"> 0.3898493 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -0.3235295 </td>
   <td style="text-align:right;"> 21.36931 </td>
   <td style="text-align:right;"> 0.1949247 </td>
   <td style="text-align:right;"> 0.3898493 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.0891988 </td>
   <td style="text-align:right;"> 31.85634 </td>
   <td style="text-align:right;"> 0.2235050 </td>
   <td style="text-align:right;"> 0.4246595 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> -0.1656677 </td>
   <td style="text-align:right;"> 22.88818 </td>
   <td style="text-align:right;"> 0.2610687 </td>
   <td style="text-align:right;"> 0.4724100 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.2132191 </td>
   <td style="text-align:right;"> 24.90489 </td>
   <td style="text-align:right;"> 0.2743829 </td>
   <td style="text-align:right;"> 0.4739341 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.2081944 </td>
   <td style="text-align:right;"> 24.87372 </td>
   <td style="text-align:right;"> 0.3498308 </td>
   <td style="text-align:right;"> 0.5779813 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.0965410 </td>
   <td style="text-align:right;"> 27.91060 </td>
   <td style="text-align:right;"> 0.3750950 </td>
   <td style="text-align:right;"> 0.5899769 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.1811179 </td>
   <td style="text-align:right;"> 23.87678 </td>
   <td style="text-align:right;"> 0.3881427 </td>
   <td style="text-align:right;"> 0.5899769 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> -0.2360763 </td>
   <td style="text-align:right;"> 23.03815 </td>
   <td style="text-align:right;"> 0.5486301 </td>
   <td style="text-align:right;"> 0.7649104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.1068963 </td>
   <td style="text-align:right;"> 23.03265 </td>
   <td style="text-align:right;"> 0.5698335 </td>
   <td style="text-align:right;"> 0.7649104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -0.0824796 </td>
   <td style="text-align:right;"> 21.96019 </td>
   <td style="text-align:right;"> 0.5841704 </td>
   <td style="text-align:right;"> 0.7649104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> -0.0705512 </td>
   <td style="text-align:right;"> 23.00858 </td>
   <td style="text-align:right;"> 0.6014076 </td>
   <td style="text-align:right;"> 0.7649104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.1128546 </td>
   <td style="text-align:right;"> 25.43956 </td>
   <td style="text-align:right;"> 0.6062686 </td>
   <td style="text-align:right;"> 0.7649104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -0.0615707 </td>
   <td style="text-align:right;"> 25.95109 </td>
   <td style="text-align:right;"> 0.6240059 </td>
   <td style="text-align:right;"> 0.7649104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 0.0584719 </td>
   <td style="text-align:right;"> 30.92011 </td>
   <td style="text-align:right;"> 0.7000791 </td>
   <td style="text-align:right;"> 0.8313440 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.0149509 </td>
   <td style="text-align:right;"> 29.08623 </td>
   <td style="text-align:right;"> 0.9152979 </td>
   <td style="text-align:right;"> 0.9692513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.0105296 </td>
   <td style="text-align:right;"> 23.82664 </td>
   <td style="text-align:right;"> 0.9346263 </td>
   <td style="text-align:right;"> 0.9692513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.0092695 </td>
   <td style="text-align:right;"> 22.86559 </td>
   <td style="text-align:right;"> 0.9440157 </td>
   <td style="text-align:right;"> 0.9692513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.0159686 </td>
   <td style="text-align:right;"> 22.17590 </td>
   <td style="text-align:right;"> 0.9557054 </td>
   <td style="text-align:right;"> 0.9692513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.0105744 </td>
   <td style="text-align:right;"> 26.74184 </td>
   <td style="text-align:right;"> 0.9610081 </td>
   <td style="text-align:right;"> 0.9692513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.0327929 </td>
   <td style="text-align:right;"> 21.34046 </td>
   <td style="text-align:right;"> 0.9692513 </td>
   <td style="text-align:right;"> 0.9692513 </td>
  </tr>
</tbody>
</table></div>

For this comparison 8 proteins are significantly different with adjusted P-values<0.05
<br>  
<br>  

__Linear model statistics for 96h vs non-stroke control protein differences__


```r
#Table of the C v 96 h post stroke results
kable(S96vC, caption="Model: ~0 + Group + Sex + HDLC") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Model: ~0 + Group + Sex + HDLC</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> logFC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AveExpr </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> P.Value </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> adj.P.Val </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 1.8441569 </td>
   <td style="text-align:right;"> 21.52132 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0000001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 1.4722726 </td>
   <td style="text-align:right;"> 20.78969 </td>
   <td style="text-align:right;"> 0.0000002 </td>
   <td style="text-align:right;"> 0.0000039 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 4.1546373 </td>
   <td style="text-align:right;"> 22.92835 </td>
   <td style="text-align:right;"> 0.0000111 </td>
   <td style="text-align:right;"> 0.0001406 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.8133407 </td>
   <td style="text-align:right;"> 21.78068 </td>
   <td style="text-align:right;"> 0.0000187 </td>
   <td style="text-align:right;"> 0.0001772 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.6795287 </td>
   <td style="text-align:right;"> 23.13573 </td>
   <td style="text-align:right;"> 0.0001296 </td>
   <td style="text-align:right;"> 0.0009847 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.7305689 </td>
   <td style="text-align:right;"> 28.56741 </td>
   <td style="text-align:right;"> 0.0004495 </td>
   <td style="text-align:right;"> 0.0028470 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.7023154 </td>
   <td style="text-align:right;"> 24.90489 </td>
   <td style="text-align:right;"> 0.0026366 </td>
   <td style="text-align:right;"> 0.0126489 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -0.6997240 </td>
   <td style="text-align:right;"> 23.61617 </td>
   <td style="text-align:right;"> 0.0026629 </td>
   <td style="text-align:right;"> 0.0126489 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -1.0304440 </td>
   <td style="text-align:right;"> 19.79977 </td>
   <td style="text-align:right;"> 0.0042404 </td>
   <td style="text-align:right;"> 0.0179038 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.4022956 </td>
   <td style="text-align:right;"> 25.95109 </td>
   <td style="text-align:right;"> 0.0073113 </td>
   <td style="text-align:right;"> 0.0277830 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.4211697 </td>
   <td style="text-align:right;"> 24.67599 </td>
   <td style="text-align:right;"> 0.0083838 </td>
   <td style="text-align:right;"> 0.0289622 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.2835779 </td>
   <td style="text-align:right;"> 29.39118 </td>
   <td style="text-align:right;"> 0.0115905 </td>
   <td style="text-align:right;"> 0.0367032 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.7108245 </td>
   <td style="text-align:right;"> 20.50256 </td>
   <td style="text-align:right;"> 0.0173352 </td>
   <td style="text-align:right;"> 0.0506722 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 2.5882092 </td>
   <td style="text-align:right;"> 25.05205 </td>
   <td style="text-align:right;"> 0.0195461 </td>
   <td style="text-align:right;"> 0.0530538 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.1870728 </td>
   <td style="text-align:right;"> 27.91060 </td>
   <td style="text-align:right;"> 0.1434371 </td>
   <td style="text-align:right;"> 0.3523334 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.3563333 </td>
   <td style="text-align:right;"> 23.87678 </td>
   <td style="text-align:right;"> 0.1483509 </td>
   <td style="text-align:right;"> 0.3523334 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.1165764 </td>
   <td style="text-align:right;"> 31.85634 </td>
   <td style="text-align:right;"> 0.1744125 </td>
   <td style="text-align:right;"> 0.3898633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 0.2156339 </td>
   <td style="text-align:right;"> 24.35341 </td>
   <td style="text-align:right;"> 0.2157742 </td>
   <td style="text-align:right;"> 0.4555233 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 0.1985329 </td>
   <td style="text-align:right;"> 26.52012 </td>
   <td style="text-align:right;"> 0.2734324 </td>
   <td style="text-align:right;"> 0.5468648 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -0.2989347 </td>
   <td style="text-align:right;"> 21.36931 </td>
   <td style="text-align:right;"> 0.3051884 </td>
   <td style="text-align:right;"> 0.5690547 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.1558148 </td>
   <td style="text-align:right;"> 22.86559 </td>
   <td style="text-align:right;"> 0.3144776 </td>
   <td style="text-align:right;"> 0.5690547 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> -0.1737999 </td>
   <td style="text-align:right;"> 26.82370 </td>
   <td style="text-align:right;"> 0.3950248 </td>
   <td style="text-align:right;"> 0.6348584 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.2341405 </td>
   <td style="text-align:right;"> 21.77589 </td>
   <td style="text-align:right;"> 0.3999422 </td>
   <td style="text-align:right;"> 0.6348584 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.2138963 </td>
   <td style="text-align:right;"> 24.87372 </td>
   <td style="text-align:right;"> 0.4114440 </td>
   <td style="text-align:right;"> 0.6348584 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> -0.3629883 </td>
   <td style="text-align:right;"> 22.70342 </td>
   <td style="text-align:right;"> 0.4176700 </td>
   <td style="text-align:right;"> 0.6348584 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.1926856 </td>
   <td style="text-align:right;"> 26.74184 </td>
   <td style="text-align:right;"> 0.4472140 </td>
   <td style="text-align:right;"> 0.6536204 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.1153719 </td>
   <td style="text-align:right;"> 23.00858 </td>
   <td style="text-align:right;"> 0.4656942 </td>
   <td style="text-align:right;"> 0.6554215 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.1532364 </td>
   <td style="text-align:right;"> 23.03265 </td>
   <td style="text-align:right;"> 0.4864937 </td>
   <td style="text-align:right;"> 0.6602414 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.0885091 </td>
   <td style="text-align:right;"> 23.82664 </td>
   <td style="text-align:right;"> 0.5560835 </td>
   <td style="text-align:right;"> 0.7286611 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.5439231 </td>
   <td style="text-align:right;"> 21.34046 </td>
   <td style="text-align:right;"> 0.5852301 </td>
   <td style="text-align:right;"> 0.7412914 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.1190615 </td>
   <td style="text-align:right;"> 25.43956 </td>
   <td style="text-align:right;"> 0.6420964 </td>
   <td style="text-align:right;"> 0.7870860 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0535073 </td>
   <td style="text-align:right;"> 22.88818 </td>
   <td style="text-align:right;"> 0.7555999 </td>
   <td style="text-align:right;"> 0.8972749 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -0.0306988 </td>
   <td style="text-align:right;"> 29.08623 </td>
   <td style="text-align:right;"> 0.8519550 </td>
   <td style="text-align:right;"> 0.9810391 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -0.0201577 </td>
   <td style="text-align:right;"> 21.96019 </td>
   <td style="text-align:right;"> 0.9089146 </td>
   <td style="text-align:right;"> 0.9917259 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0083315 </td>
   <td style="text-align:right;"> 23.30947 </td>
   <td style="text-align:right;"> 0.9634383 </td>
   <td style="text-align:right;"> 0.9917259 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> -0.0083650 </td>
   <td style="text-align:right;"> 23.03815 </td>
   <td style="text-align:right;"> 0.9854931 </td>
   <td style="text-align:right;"> 0.9917259 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.0054109 </td>
   <td style="text-align:right;"> 22.17590 </td>
   <td style="text-align:right;"> 0.9871668 </td>
   <td style="text-align:right;"> 0.9917259 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 0.0018411 </td>
   <td style="text-align:right;"> 30.92011 </td>
   <td style="text-align:right;"> 0.9917259 </td>
   <td style="text-align:right;"> 0.9917259 </td>
  </tr>
</tbody>
</table></div>

For this comparison 12 proteins are significantly different with adjusted P-values<0.05
<br>
<br>


## Figure 2. Significantly changing HDL protein abundance distributions

The log2 relative intensity is plotted by group for the 8 proteins found to be significantly changing both between the 24 h and 96 h post stroke and non-stroke control.


```r
allViolinplot <- function(a, b){
  cols <- c("C"="#202f52", "S24" = "goldenrod1", "S96" = "firebrick3")
  ggplot(a, aes(x = factor(Group), y = b, fill = Group)) +
    geom_violin(alpha=0.15, width=0.8) +
    geom_boxplot(width=0.4, alpha = 0.8) +
    geom_jitter(fill = "black", size = 2, shape = 21, 
                position = position_jitter(0.1), alpha=0.5) +
    scale_fill_manual(values=cols) +
    scale_x_discrete(limits=c("C","S24","S96")) +
    labs(y="log2 relative intensity") + 
    theme_bw() +
    theme(legend.position = "none", axis.title.x = element_blank())
}

#Violin plots for the 8 significantly changing proteins from non-stroke to 24 h post-stroke
va <- allViolinplot(logAll, logAll$APMAP) +
  ggtitle("APMAP")

vd <- allViolinplot(logAll, logAll$GPLD1) + 
  ggtitle("GPLD1") 

vb <- allViolinplot(logAll, logAll$APOE) + 
  ggtitle("APOE")

vc <- allViolinplot(logAll, logAll$IHH) + 
  ggtitle("IHH")

vg <- allViolinplot(logAll, logAll$ITIH4) +
  ggtitle("ITIH4")

ve <- allViolinplot(logAll, logAll$SAA2) +
  ggtitle("SAA2")

vf <- allViolinplot(logAll, logAll$APOA4) +
  ggtitle("APOA4") 

vh <- allViolinplot(logAll, logAll$CLU) +
  ggtitle("CLU") 


grid.arrange(va,vb,vc,vd,ve,vf,vg,vh, ncol=4)
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

<br>
<br>

## Figure 3. Paired significantly changing HDL protein abundances

The log2 relative intensity is plotted by group for the proteins found to be significantly changing both between the 24 h and 96 h post stroke paired samples. Black lines connect the paired samples between the two time point measurements.


```r
connected_boxplot <- function(a,b,c,d){ggplot(a, aes(x=b, y = c)) +
    geom_boxplot(aes(fill=b, group=b), alpha=0.8, 
                 width=0.5, show.legend = FALSE) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    geom_line(alpha=0.5, size=0.6, aes(group = d)) +
    geom_point(size=2, alpha=0.6, aes(group=b), 
               position = position_dodge(width=0.75)) +
    scale_fill_manual(values = c("goldenrod1","firebrick3")) +
    labs(y="log2 relative intensity") + 
    theme_bw() +
    theme(axis.title.x = element_blank())
}

a <- connected_boxplot(logP, logP$Group, logP$ANTXR2, logP$Pair_ID) + 
  ggtitle("ANTXR2")

b <- connected_boxplot(logP, logP$Group, logP$APOF, logP$Pair_ID) + 
  ggtitle("APOF")

c <- connected_boxplot(logP, logP$Group, logP$SAA2, logP$Pair_ID) + 
  ggtitle("SAA2") + scale_y_continuous(limits = c(16.5, 28.8))

grid.arrange(a,b,c, ncol=3)
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


<br>

# 4. Correlation testing for post stroke changes

__Log2FC values calculated for all protein measures in paired 24 h and 96 h samples__

A new dataframe is created containing the log2FC of proteins for all the paired samples. These will be used for correlation and linear regression analysis. 

```r
#Separating out the two time points of paired data
P24 <- logP %>% filter(Group=="S24")
P24rownames <- P24[,1]
rownames(P24) <- P24rownames
#renaming the measures to include which time point they belong to
colnames(P24) <- paste("S24", colnames(P24), sep = "_")

P96 <- logP %>% filter(Group=="S96")
P96rownames <- P96[,1]
rownames(P96) <- P96rownames
colnames(P96) <- paste("S96", colnames(P96), sep = "_")

#combining the renamed subsets
l_paired <- cbind(P24, P96)

#Actually calculating the log2FC of the pairs for each protein 
#(there's probably a better/more elegant way to do this, 
# but I never went back to figure it out... this works tho lol)
lfc <- data.frame(l_paired$S24_Pair_ID)
lfc$SERPINA1 <- log2(l_paired$S96_SERPINA1/l_paired$S24_SERPINA1)
lfc$ALB <- log2(l_paired$S96_ALB/l_paired$S24_ALB)
lfc$AMBP <- log2(l_paired$S96_AMBP/l_paired$S24_AMBP)
lfc$ANTXR2 <- log2(l_paired$S96_ANTXR2/l_paired$S24_ANTXR2)
lfc$APMAP <- log2(l_paired$S96_APMAP/l_paired$S24_APMAP)
lfc$LPA <- log2(l_paired$S96_LPA/l_paired$S24_LPA)
lfc$APOA1 <- log2(l_paired$S96_APOA1/l_paired$S24_APOA1)
lfc$APOA2 <- log2(l_paired$S96_APOA2/l_paired$S24_APOA2)
lfc$APOA4 <- log2(l_paired$S96_APOA4/l_paired$S24_APOA4)
lfc$APOB <- log2(l_paired$S96_APOB/l_paired$S24_APOB)
lfc$APOC3 <- log2(l_paired$S96_APOC3/l_paired$S24_APOC3)
lfc$APOC4 <- log2(l_paired$S96_APOC4/l_paired$S24_APOC4)
lfc$APOD <- log2(l_paired$S96_APOD/l_paired$S24_APOD)
lfc$APOE <- log2(l_paired$S96_APOE/l_paired$S24_APOE)
lfc$APOF <- log2(l_paired$S96_APOF/l_paired$S24_APOF)
lfc$APOL1 <- log2(l_paired$S96_APOL1/l_paired$S24_APOL1)
lfc$APOM <- log2(l_paired$S96_APOM/l_paired$S24_APOM)
lfc$CAMP <- log2(l_paired$S96_CAMP/l_paired$S24_CAMP)
lfc$CLU <- log2(l_paired$S96_CLU/l_paired$S24_CLU)
lfc$C3 <- log2(l_paired$S96_C3/l_paired$S24_C3)
lfc$PPBP <- log2(l_paired$S96_PPBP/l_paired$S24_PPBP)
lfc$AHSG <- log2(l_paired$S96_AHSG/l_paired$S24_AHSG)
lfc$FGA <- log2(l_paired$S96_FGA/l_paired$S24_FGA)
lfc$HPR <- log2(l_paired$S96_HPR/l_paired$S24_HPR)
lfc$IHH <- log2(l_paired$S96_IHH/l_paired$S24_IHH)
lfc$ITIH4<- log2(l_paired$S96_ITIH4/l_paired$S24_ITIH4)
lfc$LCAT <- log2(l_paired$S96_LCAT/l_paired$S24_LCAT)
lfc$SELL <- log2(l_paired$S96_SELL/l_paired$S24_SELL)
lfc$PCYOX1<- log2(l_paired$S96_PCYOX1/l_paired$S24_PCYOX1)
lfc$GPLD1 <- log2(l_paired$S96_GPLD1/l_paired$S24_GPLD1)
lfc$PF4 <- log2(l_paired$S96_PF4/l_paired$S24_PF4)
lfc$PLTP <- log2(l_paired$S96_PLTP/l_paired$S24_PLTP)
lfc$PON1 <- log2(l_paired$S96_PON1/l_paired$S24_PON1)
lfc$PON3 <- log2(l_paired$S96_PON3/l_paired$S24_PON3)
lfc$RBP4 <- log2(l_paired$S96_RBP4/l_paired$S24_RBP4)
lfc$SAA1 <- log2(l_paired$S96_SAA1/l_paired$S24_SAA1)
lfc$SAA2 <- log2(l_paired$S96_SAA2/l_paired$S24_SAA2)
lfc$TTR <- log2(l_paired$S96_TTR/l_paired$S24_TTR)

meta_p <- logP[1:20,1:23]
lfc_df <- cbind(meta_p, lfc)

#generating a table for manual inspection if so inclined
kable(lfc_df) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Sample_no. </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Sex </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Age_blood </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Group </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Stroke </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Pair </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Pair_ID </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Ischemic_Type </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> tPa </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Thrombectomy </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> TICI2B </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Hour </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Cholesterol </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Triglycerides </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> HDLC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LDLC </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Statin </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> NIHSS_baseline </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> NIHSS_3mo </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> mRS </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CEC </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> MS_RunID </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> MS_Batch </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> l_paired.S24_Pair_ID </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SERPINA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ALB </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AMBP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ANTXR2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APMAP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LPA </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOB </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOC3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOC4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOD </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOE </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOF </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOL1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOM </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CAMP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CLU </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> C3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PPBP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AHSG </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> FGA </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> HPR </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> IHH </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ITIH4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LCAT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SELL </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PCYOX1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> GPLD1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PF4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PLTP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PON1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PON3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> RBP4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SAA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SAA2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> TTR </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 36 </td>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P1 </td>
   <td style="text-align:left;"> SV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 223 </td>
   <td style="text-align:right;"> 507 </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 304.4 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 11.520752 </td>
   <td style="text-align:left;"> R59 </td>
   <td style="text-align:left;"> B5 </td>
   <td style="text-align:left;"> P1 </td>
   <td style="text-align:right;"> -0.0162519 </td>
   <td style="text-align:right;"> -0.0211727 </td>
   <td style="text-align:right;"> 0.0164035 </td>
   <td style="text-align:right;"> -0.0458165 </td>
   <td style="text-align:right;"> 0.0778548 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> -0.0028026 </td>
   <td style="text-align:right;"> -0.0007504 </td>
   <td style="text-align:right;"> -0.0130053 </td>
   <td style="text-align:right;"> -0.0517495 </td>
   <td style="text-align:right;"> 0.0030251 </td>
   <td style="text-align:right;"> 0.0395652 </td>
   <td style="text-align:right;"> 0.0277009 </td>
   <td style="text-align:right;"> -0.0077149 </td>
   <td style="text-align:right;"> 0.0423464 </td>
   <td style="text-align:right;"> -0.0125286 </td>
   <td style="text-align:right;"> 0.0086376 </td>
   <td style="text-align:right;"> -0.0202327 </td>
   <td style="text-align:right;"> 0.0100416 </td>
   <td style="text-align:right;"> 0.0146075 </td>
   <td style="text-align:right;"> 0.0142019 </td>
   <td style="text-align:right;"> -0.0568123 </td>
   <td style="text-align:right;"> -0.0340262 </td>
   <td style="text-align:right;"> 0.0860573 </td>
   <td style="text-align:right;"> -0.0132653 </td>
   <td style="text-align:right;"> 0.0493746 </td>
   <td style="text-align:right;"> 0.0238093 </td>
   <td style="text-align:right;"> 0.0378632 </td>
   <td style="text-align:right;"> -0.0028870 </td>
   <td style="text-align:right;"> 0.0992712 </td>
   <td style="text-align:right;"> -0.0022280 </td>
   <td style="text-align:right;"> -0.0077927 </td>
   <td style="text-align:right;"> 0.0178773 </td>
   <td style="text-align:right;"> 0.0129932 </td>
   <td style="text-align:right;"> -0.0140763 </td>
   <td style="text-align:right;"> 0.0295023 </td>
   <td style="text-align:right;"> 0.0328591 </td>
   <td style="text-align:right;"> -0.0167884 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 64 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P10 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 183 </td>
   <td style="text-align:right;"> 122 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 176.4 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 10.743049 </td>
   <td style="text-align:left;"> R8 </td>
   <td style="text-align:left;"> B1 </td>
   <td style="text-align:left;"> P10 </td>
   <td style="text-align:right;"> -0.0165236 </td>
   <td style="text-align:right;"> -0.0320940 </td>
   <td style="text-align:right;"> -0.0270564 </td>
   <td style="text-align:right;"> -0.0298956 </td>
   <td style="text-align:right;"> 0.0618072 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> -0.0251449 </td>
   <td style="text-align:right;"> -0.0493964 </td>
   <td style="text-align:right;"> -0.0376468 </td>
   <td style="text-align:right;"> 0.1053799 </td>
   <td style="text-align:right;"> 0.0031132 </td>
   <td style="text-align:right;"> -0.0161203 </td>
   <td style="text-align:right;"> 0.0016648 </td>
   <td style="text-align:right;"> -0.0013314 </td>
   <td style="text-align:right;"> 0.0182685 </td>
   <td style="text-align:right;"> 0.0077024 </td>
   <td style="text-align:right;"> -0.0004795 </td>
   <td style="text-align:right;"> -0.0188520 </td>
   <td style="text-align:right;"> -0.0045082 </td>
   <td style="text-align:right;"> 0.0157462 </td>
   <td style="text-align:right;"> 0.0179558 </td>
   <td style="text-align:right;"> -0.0662339 </td>
   <td style="text-align:right;"> -0.0035634 </td>
   <td style="text-align:right;"> 0.0254926 </td>
   <td style="text-align:right;"> -0.0077791 </td>
   <td style="text-align:right;"> 0.0332751 </td>
   <td style="text-align:right;"> -0.0341188 </td>
   <td style="text-align:right;"> 0.0243510 </td>
   <td style="text-align:right;"> -0.0016234 </td>
   <td style="text-align:right;"> 0.0200795 </td>
   <td style="text-align:right;"> 0.0379325 </td>
   <td style="text-align:right;"> 0.0084318 </td>
   <td style="text-align:right;"> 0.0186039 </td>
   <td style="text-align:right;"> 0.0144359 </td>
   <td style="text-align:right;"> -0.0403740 </td>
   <td style="text-align:right;"> 0.0681345 </td>
   <td style="text-align:right;"> 0.0414906 </td>
   <td style="text-align:right;"> -0.0596644 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 82 </td>
   <td style="text-align:right;"> 82 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 66 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P11 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 166 </td>
   <td style="text-align:right;"> 125 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 152.0 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 14.879076 </td>
   <td style="text-align:left;"> R72 </td>
   <td style="text-align:left;"> B6 </td>
   <td style="text-align:left;"> P11 </td>
   <td style="text-align:right;"> -0.0175976 </td>
   <td style="text-align:right;"> -0.0370019 </td>
   <td style="text-align:right;"> -0.0353152 </td>
   <td style="text-align:right;"> -0.0274479 </td>
   <td style="text-align:right;"> 0.0245790 </td>
   <td style="text-align:right;"> 0.1099388 </td>
   <td style="text-align:right;"> 0.0254271 </td>
   <td style="text-align:right;"> 0.0623196 </td>
   <td style="text-align:right;"> 0.0010516 </td>
   <td style="text-align:right;"> -0.0034507 </td>
   <td style="text-align:right;"> -0.0022647 </td>
   <td style="text-align:right;"> -0.0279623 </td>
   <td style="text-align:right;"> 0.0129528 </td>
   <td style="text-align:right;"> -0.0015492 </td>
   <td style="text-align:right;"> 0.0285717 </td>
   <td style="text-align:right;"> 0.0099083 </td>
   <td style="text-align:right;"> 0.0001410 </td>
   <td style="text-align:right;"> 0.0150989 </td>
   <td style="text-align:right;"> 0.0163502 </td>
   <td style="text-align:right;"> -0.0205908 </td>
   <td style="text-align:right;"> -0.0495944 </td>
   <td style="text-align:right;"> -0.0780781 </td>
   <td style="text-align:right;"> -0.0638416 </td>
   <td style="text-align:right;"> -0.0566164 </td>
   <td style="text-align:right;"> -0.0085645 </td>
   <td style="text-align:right;"> -0.1326236 </td>
   <td style="text-align:right;"> 0.0211290 </td>
   <td style="text-align:right;"> -0.0399961 </td>
   <td style="text-align:right;"> 0.0058161 </td>
   <td style="text-align:right;"> 0.0297280 </td>
   <td style="text-align:right;"> -0.0527539 </td>
   <td style="text-align:right;"> 0.0235933 </td>
   <td style="text-align:right;"> 0.0081989 </td>
   <td style="text-align:right;"> 0.0288973 </td>
   <td style="text-align:right;"> -0.0507394 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> -0.0613104 </td>
   <td style="text-align:right;"> -0.0418688 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 51 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 77 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P12 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 97 </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 74.2 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8.785186 </td>
   <td style="text-align:left;"> R49 </td>
   <td style="text-align:left;"> B4 </td>
   <td style="text-align:left;"> P12 </td>
   <td style="text-align:right;"> -0.0362716 </td>
   <td style="text-align:right;"> -0.0667583 </td>
   <td style="text-align:right;"> -0.0661943 </td>
   <td style="text-align:right;"> -0.0934201 </td>
   <td style="text-align:right;"> -0.0095972 </td>
   <td style="text-align:right;"> -0.0666377 </td>
   <td style="text-align:right;"> -0.0098154 </td>
   <td style="text-align:right;"> -0.0514316 </td>
   <td style="text-align:right;"> -0.0495608 </td>
   <td style="text-align:right;"> -0.0288176 </td>
   <td style="text-align:right;"> -0.0452815 </td>
   <td style="text-align:right;"> -0.0223248 </td>
   <td style="text-align:right;"> -0.0420016 </td>
   <td style="text-align:right;"> -0.0117322 </td>
   <td style="text-align:right;"> 0.0463630 </td>
   <td style="text-align:right;"> 0.0047743 </td>
   <td style="text-align:right;"> -0.0088821 </td>
   <td style="text-align:right;"> -0.0061413 </td>
   <td style="text-align:right;"> -0.0401639 </td>
   <td style="text-align:right;"> -0.0449103 </td>
   <td style="text-align:right;"> 0.0747248 </td>
   <td style="text-align:right;"> -0.1593122 </td>
   <td style="text-align:right;"> -0.0249712 </td>
   <td style="text-align:right;"> -0.1084357 </td>
   <td style="text-align:right;"> -0.0356900 </td>
   <td style="text-align:right;"> 0.0305746 </td>
   <td style="text-align:right;"> -0.0596022 </td>
   <td style="text-align:right;"> 0.0520443 </td>
   <td style="text-align:right;"> -0.0103006 </td>
   <td style="text-align:right;"> -0.0119375 </td>
   <td style="text-align:right;"> 0.0698738 </td>
   <td style="text-align:right;"> -0.0389862 </td>
   <td style="text-align:right;"> 0.0013368 </td>
   <td style="text-align:right;"> -0.0221461 </td>
   <td style="text-align:right;"> -0.0789172 </td>
   <td style="text-align:right;"> 0.1512222 </td>
   <td style="text-align:right;"> 0.1727510 </td>
   <td style="text-align:right;"> -0.0838937 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P13 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 169 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 152.2 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 12.364192 </td>
   <td style="text-align:left;"> R68 </td>
   <td style="text-align:left;"> B5 </td>
   <td style="text-align:left;"> P13 </td>
   <td style="text-align:right;"> -0.0167266 </td>
   <td style="text-align:right;"> -0.0731250 </td>
   <td style="text-align:right;"> -0.0031603 </td>
   <td style="text-align:right;"> -0.0953670 </td>
   <td style="text-align:right;"> -0.0200160 </td>
   <td style="text-align:right;"> -0.0061674 </td>
   <td style="text-align:right;"> -0.0285667 </td>
   <td style="text-align:right;"> -0.0578663 </td>
   <td style="text-align:right;"> -0.0219179 </td>
   <td style="text-align:right;"> -0.0014634 </td>
   <td style="text-align:right;"> 0.0102624 </td>
   <td style="text-align:right;"> -0.0192775 </td>
   <td style="text-align:right;"> -0.0155491 </td>
   <td style="text-align:right;"> -0.0375008 </td>
   <td style="text-align:right;"> 0.0560340 </td>
   <td style="text-align:right;"> 0.0038776 </td>
   <td style="text-align:right;"> -0.0026817 </td>
   <td style="text-align:right;"> 0.0116526 </td>
   <td style="text-align:right;"> -0.0807943 </td>
   <td style="text-align:right;"> -0.0181448 </td>
   <td style="text-align:right;"> 0.0076121 </td>
   <td style="text-align:right;"> -0.1063057 </td>
   <td style="text-align:right;"> -0.1598408 </td>
   <td style="text-align:right;"> -0.0456496 </td>
   <td style="text-align:right;"> 0.0070478 </td>
   <td style="text-align:right;"> 0.0511329 </td>
   <td style="text-align:right;"> 0.0129934 </td>
   <td style="text-align:right;"> -0.0747218 </td>
   <td style="text-align:right;"> -0.0218000 </td>
   <td style="text-align:right;"> -0.1002865 </td>
   <td style="text-align:right;"> -0.0068125 </td>
   <td style="text-align:right;"> -0.0646607 </td>
   <td style="text-align:right;"> -0.0035575 </td>
   <td style="text-align:right;"> -0.0150056 </td>
   <td style="text-align:right;"> -0.0220377 </td>
   <td style="text-align:right;"> 0.2696641 </td>
   <td style="text-align:right;"> 0.3899098 </td>
   <td style="text-align:right;"> -0.0582382 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 35 </td>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P14 </td>
   <td style="text-align:left;"> CE </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 175 </td>
   <td style="text-align:right;"> 79 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 145.8 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8.794060 </td>
   <td style="text-align:left;"> R32 </td>
   <td style="text-align:left;"> B3 </td>
   <td style="text-align:left;"> P14 </td>
   <td style="text-align:right;"> 0.0102192 </td>
   <td style="text-align:right;"> -0.0497208 </td>
   <td style="text-align:right;"> -0.0306116 </td>
   <td style="text-align:right;"> -0.0453621 </td>
   <td style="text-align:right;"> -0.0264574 </td>
   <td style="text-align:right;"> 0.1524945 </td>
   <td style="text-align:right;"> -0.0032718 </td>
   <td style="text-align:right;"> -0.0282794 </td>
   <td style="text-align:right;"> -0.0615768 </td>
   <td style="text-align:right;"> 0.0155322 </td>
   <td style="text-align:right;"> -0.0194542 </td>
   <td style="text-align:right;"> 0.0383385 </td>
   <td style="text-align:right;"> 0.0119545 </td>
   <td style="text-align:right;"> -0.0131756 </td>
   <td style="text-align:right;"> 0.0414550 </td>
   <td style="text-align:right;"> -0.0073366 </td>
   <td style="text-align:right;"> -0.0325885 </td>
   <td style="text-align:right;"> 0.0205607 </td>
   <td style="text-align:right;"> -0.0222011 </td>
   <td style="text-align:right;"> 0.0062703 </td>
   <td style="text-align:right;"> 0.0589253 </td>
   <td style="text-align:right;"> -0.0968503 </td>
   <td style="text-align:right;"> -0.1246225 </td>
   <td style="text-align:right;"> 0.0520609 </td>
   <td style="text-align:right;"> -0.0185971 </td>
   <td style="text-align:right;"> 0.0061414 </td>
   <td style="text-align:right;"> -0.0327731 </td>
   <td style="text-align:right;"> -0.0544129 </td>
   <td style="text-align:right;"> -0.0228320 </td>
   <td style="text-align:right;"> -0.0403451 </td>
   <td style="text-align:right;"> 0.0480985 </td>
   <td style="text-align:right;"> -0.0180940 </td>
   <td style="text-align:right;"> 0.0061140 </td>
   <td style="text-align:right;"> 0.0068250 </td>
   <td style="text-align:right;"> -0.0396628 </td>
   <td style="text-align:right;"> 0.1653982 </td>
   <td style="text-align:right;"> 0.2196745 </td>
   <td style="text-align:right;"> -0.0349443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 32 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 69 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P15 </td>
   <td style="text-align:left;"> CE </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 56 </td>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 90.2 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 7.110238 </td>
   <td style="text-align:left;"> R63 </td>
   <td style="text-align:left;"> B5 </td>
   <td style="text-align:left;"> P15 </td>
   <td style="text-align:right;"> 0.0718284 </td>
   <td style="text-align:right;"> 0.0784698 </td>
   <td style="text-align:right;"> 0.0347289 </td>
   <td style="text-align:right;"> -0.0222841 </td>
   <td style="text-align:right;"> -0.0537437 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> 0.0372006 </td>
   <td style="text-align:right;"> 0.0605878 </td>
   <td style="text-align:right;"> 0.0053343 </td>
   <td style="text-align:right;"> -0.0448119 </td>
   <td style="text-align:right;"> -0.0175719 </td>
   <td style="text-align:right;"> -0.0141607 </td>
   <td style="text-align:right;"> 0.0054286 </td>
   <td style="text-align:right;"> -0.0079815 </td>
   <td style="text-align:right;"> 0.0086964 </td>
   <td style="text-align:right;"> -0.0225217 </td>
   <td style="text-align:right;"> 0.0049645 </td>
   <td style="text-align:right;"> 0.0987590 </td>
   <td style="text-align:right;"> 0.0072823 </td>
   <td style="text-align:right;"> 0.0315582 </td>
   <td style="text-align:right;"> 0.0274109 </td>
   <td style="text-align:right;"> 0.0569269 </td>
   <td style="text-align:right;"> -0.0310407 </td>
   <td style="text-align:right;"> -0.0710213 </td>
   <td style="text-align:right;"> -0.0266047 </td>
   <td style="text-align:right;"> -0.2144181 </td>
   <td style="text-align:right;"> 0.0253529 </td>
   <td style="text-align:right;"> -0.0096523 </td>
   <td style="text-align:right;"> -0.0315126 </td>
   <td style="text-align:right;"> -0.0238004 </td>
   <td style="text-align:right;"> 0.0077375 </td>
   <td style="text-align:right;"> -0.0120095 </td>
   <td style="text-align:right;"> -0.0119641 </td>
   <td style="text-align:right;"> 0.0024258 </td>
   <td style="text-align:right;"> 0.0304099 </td>
   <td style="text-align:right;"> 0.1158435 </td>
   <td style="text-align:right;"> 0.1967526 </td>
   <td style="text-align:right;"> 0.0942212 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 80 </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P16 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 206 </td>
   <td style="text-align:right;"> 206 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 216.2 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 10.917753 </td>
   <td style="text-align:left;"> R82 </td>
   <td style="text-align:left;"> B6 </td>
   <td style="text-align:left;"> P16 </td>
   <td style="text-align:right;"> 0.0441420 </td>
   <td style="text-align:right;"> 0.0280163 </td>
   <td style="text-align:right;"> 0.0149106 </td>
   <td style="text-align:right;"> -0.0069982 </td>
   <td style="text-align:right;"> -0.0405769 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> 0.0155979 </td>
   <td style="text-align:right;"> 0.0277761 </td>
   <td style="text-align:right;"> -0.0188757 </td>
   <td style="text-align:right;"> -0.0639185 </td>
   <td style="text-align:right;"> -0.0204340 </td>
   <td style="text-align:right;"> -0.0438743 </td>
   <td style="text-align:right;"> 0.0100484 </td>
   <td style="text-align:right;"> 0.0042136 </td>
   <td style="text-align:right;"> 0.0073324 </td>
   <td style="text-align:right;"> 0.0118713 </td>
   <td style="text-align:right;"> 0.0118553 </td>
   <td style="text-align:right;"> 0.0013059 </td>
   <td style="text-align:right;"> -0.0038940 </td>
   <td style="text-align:right;"> -0.0033242 </td>
   <td style="text-align:right;"> -0.0670321 </td>
   <td style="text-align:right;"> 0.0337936 </td>
   <td style="text-align:right;"> -0.0252425 </td>
   <td style="text-align:right;"> 0.0187371 </td>
   <td style="text-align:right;"> 0.0115773 </td>
   <td style="text-align:right;"> -0.1071406 </td>
   <td style="text-align:right;"> 0.0451326 </td>
   <td style="text-align:right;"> 0.0140349 </td>
   <td style="text-align:right;"> 0.0005731 </td>
   <td style="text-align:right;"> 0.0439461 </td>
   <td style="text-align:right;"> -0.0991108 </td>
   <td style="text-align:right;"> 0.0316630 </td>
   <td style="text-align:right;"> 0.0143936 </td>
   <td style="text-align:right;"> 0.0351317 </td>
   <td style="text-align:right;"> 0.0552261 </td>
   <td style="text-align:right;"> -0.1024582 </td>
   <td style="text-align:right;"> -0.0331058 </td>
   <td style="text-align:right;"> 0.0341269 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 55 </td>
   <td style="text-align:right;"> 55 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P17 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 90 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 47.0 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 8.225830 </td>
   <td style="text-align:left;"> R95 </td>
   <td style="text-align:left;"> B7 </td>
   <td style="text-align:left;"> P17 </td>
   <td style="text-align:right;"> 0.0758322 </td>
   <td style="text-align:right;"> 0.1207182 </td>
   <td style="text-align:right;"> 0.0716522 </td>
   <td style="text-align:right;"> -0.0274765 </td>
   <td style="text-align:right;"> 0.1006750 </td>
   <td style="text-align:right;"> 0.1814324 </td>
   <td style="text-align:right;"> 0.0100928 </td>
   <td style="text-align:right;"> 0.0148463 </td>
   <td style="text-align:right;"> -0.0582991 </td>
   <td style="text-align:right;"> 0.0921999 </td>
   <td style="text-align:right;"> 0.0442151 </td>
   <td style="text-align:right;"> 0.0924935 </td>
   <td style="text-align:right;"> 0.0013346 </td>
   <td style="text-align:right;"> 0.0568322 </td>
   <td style="text-align:right;"> 0.0700041 </td>
   <td style="text-align:right;"> 0.0739595 </td>
   <td style="text-align:right;"> 0.0656045 </td>
   <td style="text-align:right;"> 0.0236497 </td>
   <td style="text-align:right;"> 0.0388477 </td>
   <td style="text-align:right;"> 0.1084633 </td>
   <td style="text-align:right;"> 0.1601678 </td>
   <td style="text-align:right;"> 0.1706840 </td>
   <td style="text-align:right;"> 0.1018648 </td>
   <td style="text-align:right;"> 0.0582445 </td>
   <td style="text-align:right;"> 0.0100356 </td>
   <td style="text-align:right;"> 0.0491863 </td>
   <td style="text-align:right;"> 0.1095967 </td>
   <td style="text-align:right;"> 0.0684682 </td>
   <td style="text-align:right;"> 0.0295924 </td>
   <td style="text-align:right;"> 0.0314961 </td>
   <td style="text-align:right;"> 0.1512434 </td>
   <td style="text-align:right;"> 0.0055718 </td>
   <td style="text-align:right;"> 0.0386431 </td>
   <td style="text-align:right;"> 0.0446510 </td>
   <td style="text-align:right;"> 0.0936712 </td>
   <td style="text-align:right;"> -0.1007456 </td>
   <td style="text-align:right;"> -0.1228410 </td>
   <td style="text-align:right;"> 0.0616486 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 28 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P18 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 104 </td>
   <td style="text-align:right;"> 93 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 95.6 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 5.737339 </td>
   <td style="text-align:left;"> R35 </td>
   <td style="text-align:left;"> B3 </td>
   <td style="text-align:left;"> P18 </td>
   <td style="text-align:right;"> -0.0309279 </td>
   <td style="text-align:right;"> -0.0523985 </td>
   <td style="text-align:right;"> -0.0380789 </td>
   <td style="text-align:right;"> 0.0139750 </td>
   <td style="text-align:right;"> -0.0199593 </td>
   <td style="text-align:right;"> -0.0428676 </td>
   <td style="text-align:right;"> 0.0055553 </td>
   <td style="text-align:right;"> 0.0069617 </td>
   <td style="text-align:right;"> -0.0202215 </td>
   <td style="text-align:right;"> -0.0320383 </td>
   <td style="text-align:right;"> -0.0053369 </td>
   <td style="text-align:right;"> 0.0051857 </td>
   <td style="text-align:right;"> -0.0078150 </td>
   <td style="text-align:right;"> -0.0083247 </td>
   <td style="text-align:right;"> -0.0004000 </td>
   <td style="text-align:right;"> -0.0341965 </td>
   <td style="text-align:right;"> -0.0318527 </td>
   <td style="text-align:right;"> -0.0022945 </td>
   <td style="text-align:right;"> -0.0143749 </td>
   <td style="text-align:right;"> -0.0023729 </td>
   <td style="text-align:right;"> 0.0449815 </td>
   <td style="text-align:right;"> -0.0250035 </td>
   <td style="text-align:right;"> -0.0017748 </td>
   <td style="text-align:right;"> -0.0221846 </td>
   <td style="text-align:right;"> -0.0123760 </td>
   <td style="text-align:right;"> -0.0091763 </td>
   <td style="text-align:right;"> 0.0020679 </td>
   <td style="text-align:right;"> -0.0216303 </td>
   <td style="text-align:right;"> -0.0410220 </td>
   <td style="text-align:right;"> -0.0430002 </td>
   <td style="text-align:right;"> 0.0534480 </td>
   <td style="text-align:right;"> -0.0251716 </td>
   <td style="text-align:right;"> 0.0036321 </td>
   <td style="text-align:right;"> -0.0044785 </td>
   <td style="text-align:right;"> 0.0027456 </td>
   <td style="text-align:right;"> 0.4059867 </td>
   <td style="text-align:right;"> -0.0361021 </td>
   <td style="text-align:right;"> -0.0248089 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 52 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 64 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P19 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 121 </td>
   <td style="text-align:right;"> 120 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> 113.0 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 12.166660 </td>
   <td style="text-align:left;"> R38 </td>
   <td style="text-align:left;"> B3 </td>
   <td style="text-align:left;"> P19 </td>
   <td style="text-align:right;"> 0.0401062 </td>
   <td style="text-align:right;"> -0.0032745 </td>
   <td style="text-align:right;"> 0.0195321 </td>
   <td style="text-align:right;"> -0.0344011 </td>
   <td style="text-align:right;"> 0.0214249 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> 0.0032976 </td>
   <td style="text-align:right;"> -0.0101202 </td>
   <td style="text-align:right;"> -0.0782195 </td>
   <td style="text-align:right;"> -0.0070690 </td>
   <td style="text-align:right;"> -0.0462535 </td>
   <td style="text-align:right;"> -0.0129224 </td>
   <td style="text-align:right;"> 0.0089649 </td>
   <td style="text-align:right;"> -0.0318802 </td>
   <td style="text-align:right;"> 0.0149120 </td>
   <td style="text-align:right;"> -0.0064885 </td>
   <td style="text-align:right;"> -0.0251528 </td>
   <td style="text-align:right;"> -0.0284002 </td>
   <td style="text-align:right;"> -0.0245683 </td>
   <td style="text-align:right;"> 0.0452847 </td>
   <td style="text-align:right;"> 0.0709418 </td>
   <td style="text-align:right;"> 0.0104546 </td>
   <td style="text-align:right;"> -0.1163005 </td>
   <td style="text-align:right;"> -0.0537345 </td>
   <td style="text-align:right;"> -0.0129811 </td>
   <td style="text-align:right;"> 0.0189701 </td>
   <td style="text-align:right;"> -0.0190848 </td>
   <td style="text-align:right;"> 0.0256610 </td>
   <td style="text-align:right;"> -0.0283027 </td>
   <td style="text-align:right;"> -0.0234999 </td>
   <td style="text-align:right;"> 0.0694956 </td>
   <td style="text-align:right;"> -0.0115214 </td>
   <td style="text-align:right;"> -0.0253080 </td>
   <td style="text-align:right;"> -0.0058524 </td>
   <td style="text-align:right;"> -0.0012955 </td>
   <td style="text-align:right;"> 0.1284102 </td>
   <td style="text-align:right;"> 0.2040084 </td>
   <td style="text-align:right;"> 0.0201984 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 59 </td>
   <td style="text-align:right;"> 59 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 61 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P2 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 163 </td>
   <td style="text-align:right;"> 98 </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 114.6 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 13.458247 </td>
   <td style="text-align:left;"> R94 </td>
   <td style="text-align:left;"> B7 </td>
   <td style="text-align:left;"> P2 </td>
   <td style="text-align:right;"> 0.0026147 </td>
   <td style="text-align:right;"> 0.0215960 </td>
   <td style="text-align:right;"> 0.0414324 </td>
   <td style="text-align:right;"> 0.0620152 </td>
   <td style="text-align:right;"> -0.0189674 </td>
   <td style="text-align:right;"> -0.0860071 </td>
   <td style="text-align:right;"> -0.0096215 </td>
   <td style="text-align:right;"> -0.0199660 </td>
   <td style="text-align:right;"> 0.0275890 </td>
   <td style="text-align:right;"> 0.0313606 </td>
   <td style="text-align:right;"> 0.0754457 </td>
   <td style="text-align:right;"> -0.0140414 </td>
   <td style="text-align:right;"> 0.0083043 </td>
   <td style="text-align:right;"> 0.0077431 </td>
   <td style="text-align:right;"> -0.0432057 </td>
   <td style="text-align:right;"> 0.0403149 </td>
   <td style="text-align:right;"> 0.0113847 </td>
   <td style="text-align:right;"> 0.0218154 </td>
   <td style="text-align:right;"> 0.0177423 </td>
   <td style="text-align:right;"> 0.0110380 </td>
   <td style="text-align:right;"> 0.0020467 </td>
   <td style="text-align:right;"> 0.0202878 </td>
   <td style="text-align:right;"> 0.0515429 </td>
   <td style="text-align:right;"> 0.0507594 </td>
   <td style="text-align:right;"> 0.0615510 </td>
   <td style="text-align:right;"> -0.0102744 </td>
   <td style="text-align:right;"> 0.0752324 </td>
   <td style="text-align:right;"> -0.0261372 </td>
   <td style="text-align:right;"> 0.0232921 </td>
   <td style="text-align:right;"> 0.0129318 </td>
   <td style="text-align:right;"> 0.0252866 </td>
   <td style="text-align:right;"> 0.0446081 </td>
   <td style="text-align:right;"> 0.0383404 </td>
   <td style="text-align:right;"> 0.0326664 </td>
   <td style="text-align:right;"> 0.0870033 </td>
   <td style="text-align:right;"> -0.0121376 </td>
   <td style="text-align:right;"> -0.0136071 </td>
   <td style="text-align:right;"> 0.0024867 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 71 </td>
   <td style="text-align:right;"> 71 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 82 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P20 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 133 </td>
   <td style="text-align:right;"> 79 </td>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 98.8 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 13.429480 </td>
   <td style="text-align:left;"> R44 </td>
   <td style="text-align:left;"> B4 </td>
   <td style="text-align:left;"> P20 </td>
   <td style="text-align:right;"> 0.0000700 </td>
   <td style="text-align:right;"> 0.0174107 </td>
   <td style="text-align:right;"> 0.0335286 </td>
   <td style="text-align:right;"> -0.0046480 </td>
   <td style="text-align:right;"> -0.0647983 </td>
   <td style="text-align:right;"> 0.0321625 </td>
   <td style="text-align:right;"> 0.0068160 </td>
   <td style="text-align:right;"> -0.0085280 </td>
   <td style="text-align:right;"> 0.0016257 </td>
   <td style="text-align:right;"> 0.0361619 </td>
   <td style="text-align:right;"> 0.0054260 </td>
   <td style="text-align:right;"> -0.0545610 </td>
   <td style="text-align:right;"> -0.0074218 </td>
   <td style="text-align:right;"> -0.0015882 </td>
   <td style="text-align:right;"> 0.0074391 </td>
   <td style="text-align:right;"> 0.0053098 </td>
   <td style="text-align:right;"> 0.0134257 </td>
   <td style="text-align:right;"> 0.0088376 </td>
   <td style="text-align:right;"> -0.0100281 </td>
   <td style="text-align:right;"> -0.0106701 </td>
   <td style="text-align:right;"> -0.0871984 </td>
   <td style="text-align:right;"> 0.0414038 </td>
   <td style="text-align:right;"> -0.1291456 </td>
   <td style="text-align:right;"> -0.0225656 </td>
   <td style="text-align:right;"> 0.0082959 </td>
   <td style="text-align:right;"> 0.0669792 </td>
   <td style="text-align:right;"> 0.0223451 </td>
   <td style="text-align:right;"> -0.0099909 </td>
   <td style="text-align:right;"> -0.0156262 </td>
   <td style="text-align:right;"> -0.0162554 </td>
   <td style="text-align:right;"> -0.1090733 </td>
   <td style="text-align:right;"> -0.0001758 </td>
   <td style="text-align:right;"> -0.0052864 </td>
   <td style="text-align:right;"> -0.0011826 </td>
   <td style="text-align:right;"> 0.0412160 </td>
   <td style="text-align:right;"> -0.0326447 </td>
   <td style="text-align:right;"> -0.0664542 </td>
   <td style="text-align:right;"> 0.0228862 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 19 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 86 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P3 </td>
   <td style="text-align:left;"> CE </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 177 </td>
   <td style="text-align:right;"> 187 </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 172.4 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 9.061798 </td>
   <td style="text-align:left;"> R71 </td>
   <td style="text-align:left;"> B6 </td>
   <td style="text-align:left;"> P3 </td>
   <td style="text-align:right;"> 0.1105161 </td>
   <td style="text-align:right;"> 0.1037401 </td>
   <td style="text-align:right;"> 0.0112510 </td>
   <td style="text-align:right;"> -0.1043498 </td>
   <td style="text-align:right;"> 0.0644335 </td>
   <td style="text-align:right;"> -0.0346104 </td>
   <td style="text-align:right;"> -0.0234158 </td>
   <td style="text-align:right;"> -0.0331406 </td>
   <td style="text-align:right;"> 0.0544529 </td>
   <td style="text-align:right;"> -0.0612691 </td>
   <td style="text-align:right;"> -0.0160422 </td>
   <td style="text-align:right;"> 0.0302441 </td>
   <td style="text-align:right;"> 0.0091083 </td>
   <td style="text-align:right;"> 0.0075498 </td>
   <td style="text-align:right;"> 0.0123897 </td>
   <td style="text-align:right;"> 0.0046280 </td>
   <td style="text-align:right;"> 0.0057676 </td>
   <td style="text-align:right;"> -0.0024316 </td>
   <td style="text-align:right;"> -0.0015854 </td>
   <td style="text-align:right;"> 0.0390607 </td>
   <td style="text-align:right;"> 0.0369132 </td>
   <td style="text-align:right;"> 0.0721454 </td>
   <td style="text-align:right;"> 0.0728438 </td>
   <td style="text-align:right;"> 0.0105995 </td>
   <td style="text-align:right;"> -0.0309105 </td>
   <td style="text-align:right;"> -0.0658547 </td>
   <td style="text-align:right;"> 0.0399683 </td>
   <td style="text-align:right;"> 0.0719022 </td>
   <td style="text-align:right;"> -0.0045058 </td>
   <td style="text-align:right;"> 0.0116384 </td>
   <td style="text-align:right;"> 0.0159141 </td>
   <td style="text-align:right;"> -0.0089844 </td>
   <td style="text-align:right;"> 0.0047651 </td>
   <td style="text-align:right;"> -0.0188908 </td>
   <td style="text-align:right;"> 0.0306140 </td>
   <td style="text-align:right;"> 0.1112861 </td>
   <td style="text-align:right;"> 0.1348294 </td>
   <td style="text-align:right;"> 0.0339351 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 46 </td>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 70 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P4 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 147 </td>
   <td style="text-align:right;"> 198 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 133.6 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 11.370783 </td>
   <td style="text-align:left;"> R55 </td>
   <td style="text-align:left;"> B4 </td>
   <td style="text-align:left;"> P4 </td>
   <td style="text-align:right;"> 0.0607694 </td>
   <td style="text-align:right;"> 0.0259916 </td>
   <td style="text-align:right;"> 0.0150697 </td>
   <td style="text-align:right;"> -0.0569904 </td>
   <td style="text-align:right;"> 0.0224823 </td>
   <td style="text-align:right;"> 0.0116501 </td>
   <td style="text-align:right;"> -0.0036973 </td>
   <td style="text-align:right;"> -0.0102339 </td>
   <td style="text-align:right;"> 0.0226718 </td>
   <td style="text-align:right;"> 0.0104704 </td>
   <td style="text-align:right;"> -0.0050058 </td>
   <td style="text-align:right;"> 0.0593446 </td>
   <td style="text-align:right;"> -0.0007326 </td>
   <td style="text-align:right;"> 0.0092664 </td>
   <td style="text-align:right;"> 0.0208592 </td>
   <td style="text-align:right;"> 0.0082779 </td>
   <td style="text-align:right;"> -0.0152816 </td>
   <td style="text-align:right;"> 0.0012244 </td>
   <td style="text-align:right;"> 0.0825707 </td>
   <td style="text-align:right;"> 0.0427025 </td>
   <td style="text-align:right;"> 0.0050074 </td>
   <td style="text-align:right;"> 0.0583401 </td>
   <td style="text-align:right;"> 0.0485702 </td>
   <td style="text-align:right;"> 0.0344152 </td>
   <td style="text-align:right;"> -0.0362583 </td>
   <td style="text-align:right;"> 0.0119906 </td>
   <td style="text-align:right;"> 0.0067419 </td>
   <td style="text-align:right;"> -0.0208984 </td>
   <td style="text-align:right;"> -0.0194212 </td>
   <td style="text-align:right;"> -0.0042400 </td>
   <td style="text-align:right;"> 0.0193225 </td>
   <td style="text-align:right;"> -0.0115689 </td>
   <td style="text-align:right;"> 0.0043266 </td>
   <td style="text-align:right;"> -0.0245011 </td>
   <td style="text-align:right;"> -0.0041734 </td>
   <td style="text-align:right;"> 0.1171797 </td>
   <td style="text-align:right;"> 0.1661969 </td>
   <td style="text-align:right;"> 0.0200708 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 75 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P5 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 170 </td>
   <td style="text-align:right;"> 62 </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 145.4 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 7.108044 </td>
   <td style="text-align:left;"> R79 </td>
   <td style="text-align:left;"> B6 </td>
   <td style="text-align:left;"> P5 </td>
   <td style="text-align:right;"> 0.0808519 </td>
   <td style="text-align:right;"> 0.0968676 </td>
   <td style="text-align:right;"> 0.0584943 </td>
   <td style="text-align:right;"> -0.0531321 </td>
   <td style="text-align:right;"> -0.0002564 </td>
   <td style="text-align:right;"> 0.0338647 </td>
   <td style="text-align:right;"> 0.0182500 </td>
   <td style="text-align:right;"> 0.0390181 </td>
   <td style="text-align:right;"> 0.0610191 </td>
   <td style="text-align:right;"> -0.0651618 </td>
   <td style="text-align:right;"> 0.0170171 </td>
   <td style="text-align:right;"> 0.0640835 </td>
   <td style="text-align:right;"> -0.0072691 </td>
   <td style="text-align:right;"> 0.0087556 </td>
   <td style="text-align:right;"> 0.0292337 </td>
   <td style="text-align:right;"> 0.0228625 </td>
   <td style="text-align:right;"> -0.0095346 </td>
   <td style="text-align:right;"> 0.0585014 </td>
   <td style="text-align:right;"> 0.0165114 </td>
   <td style="text-align:right;"> 0.0165065 </td>
   <td style="text-align:right;"> 0.0310057 </td>
   <td style="text-align:right;"> 0.1138115 </td>
   <td style="text-align:right;"> 0.0244538 </td>
   <td style="text-align:right;"> 0.0148904 </td>
   <td style="text-align:right;"> -0.0039673 </td>
   <td style="text-align:right;"> -0.0059240 </td>
   <td style="text-align:right;"> 0.0804353 </td>
   <td style="text-align:right;"> -0.0259703 </td>
   <td style="text-align:right;"> 0.0076896 </td>
   <td style="text-align:right;"> 0.0008810 </td>
   <td style="text-align:right;"> 0.0156787 </td>
   <td style="text-align:right;"> -0.0218403 </td>
   <td style="text-align:right;"> 0.0219584 </td>
   <td style="text-align:right;"> 0.0185511 </td>
   <td style="text-align:right;"> 0.0472187 </td>
   <td style="text-align:right;"> 0.0944576 </td>
   <td style="text-align:right;"> 0.1358721 </td>
   <td style="text-align:right;"> 0.0607302 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 65 </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P6 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 109 </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 86.0 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 6.066922 </td>
   <td style="text-align:left;"> R69 </td>
   <td style="text-align:left;"> B5 </td>
   <td style="text-align:left;"> P6 </td>
   <td style="text-align:right;"> 0.0524508 </td>
   <td style="text-align:right;"> -0.0086915 </td>
   <td style="text-align:right;"> 0.0202029 </td>
   <td style="text-align:right;"> -0.0140439 </td>
   <td style="text-align:right;"> -0.0023883 </td>
   <td style="text-align:right;"> 0.0268196 </td>
   <td style="text-align:right;"> -0.0053101 </td>
   <td style="text-align:right;"> -0.0463372 </td>
   <td style="text-align:right;"> 0.0482856 </td>
   <td style="text-align:right;"> 0.0518124 </td>
   <td style="text-align:right;"> 0.0058412 </td>
   <td style="text-align:right;"> 0.0383240 </td>
   <td style="text-align:right;"> 0.0365546 </td>
   <td style="text-align:right;"> 0.0179778 </td>
   <td style="text-align:right;"> 0.0407431 </td>
   <td style="text-align:right;"> 0.0338091 </td>
   <td style="text-align:right;"> 0.0215280 </td>
   <td style="text-align:right;"> 0.0072130 </td>
   <td style="text-align:right;"> 0.0235945 </td>
   <td style="text-align:right;"> 0.0397855 </td>
   <td style="text-align:right;"> 0.0644094 </td>
   <td style="text-align:right;"> 0.0028046 </td>
   <td style="text-align:right;"> 0.0490920 </td>
   <td style="text-align:right;"> 0.0083515 </td>
   <td style="text-align:right;"> 0.0219151 </td>
   <td style="text-align:right;"> 0.0626677 </td>
   <td style="text-align:right;"> 0.0398245 </td>
   <td style="text-align:right;"> 0.0571733 </td>
   <td style="text-align:right;"> 0.0426975 </td>
   <td style="text-align:right;"> 0.0249095 </td>
   <td style="text-align:right;"> 0.0580147 </td>
   <td style="text-align:right;"> 0.0067026 </td>
   <td style="text-align:right;"> 0.0121223 </td>
   <td style="text-align:right;"> 0.0224465 </td>
   <td style="text-align:right;"> -0.0003232 </td>
   <td style="text-align:right;"> 0.0868791 </td>
   <td style="text-align:right;"> 0.1555783 </td>
   <td style="text-align:right;"> 0.0235610 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 53 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:left;"> MALE </td>
   <td style="text-align:right;"> 58 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P7 </td>
   <td style="text-align:left;"> CE </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 147 </td>
   <td style="text-align:right;"> 95 </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 126.0 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 10.486528 </td>
   <td style="text-align:left;"> R40 </td>
   <td style="text-align:left;"> B3 </td>
   <td style="text-align:left;"> P7 </td>
   <td style="text-align:right;"> 0.0462551 </td>
   <td style="text-align:right;"> 0.0013951 </td>
   <td style="text-align:right;"> 0.0014867 </td>
   <td style="text-align:right;"> -0.0380749 </td>
   <td style="text-align:right;"> -0.0183718 </td>
   <td style="text-align:right;"> 0.0310088 </td>
   <td style="text-align:right;"> -0.0177393 </td>
   <td style="text-align:right;"> -0.0156533 </td>
   <td style="text-align:right;"> -0.0120880 </td>
   <td style="text-align:right;"> -0.0132625 </td>
   <td style="text-align:right;"> -0.0222493 </td>
   <td style="text-align:right;"> -0.0238345 </td>
   <td style="text-align:right;"> 0.0032602 </td>
   <td style="text-align:right;"> 0.0000087 </td>
   <td style="text-align:right;"> 0.0192607 </td>
   <td style="text-align:right;"> 0.0076103 </td>
   <td style="text-align:right;"> -0.0369030 </td>
   <td style="text-align:right;"> -0.0038508 </td>
   <td style="text-align:right;"> -0.0267115 </td>
   <td style="text-align:right;"> 0.0235154 </td>
   <td style="text-align:right;"> 0.0432524 </td>
   <td style="text-align:right;"> -0.0030262 </td>
   <td style="text-align:right;"> 0.0332587 </td>
   <td style="text-align:right;"> -0.0133997 </td>
   <td style="text-align:right;"> -0.0305442 </td>
   <td style="text-align:right;"> -0.0309935 </td>
   <td style="text-align:right;"> -0.0088184 </td>
   <td style="text-align:right;"> -0.0195503 </td>
   <td style="text-align:right;"> -0.0131975 </td>
   <td style="text-align:right;"> -0.0036659 </td>
   <td style="text-align:right;"> 0.0235982 </td>
   <td style="text-align:right;"> -0.0272245 </td>
   <td style="text-align:right;"> 0.0063365 </td>
   <td style="text-align:right;"> -0.0251069 </td>
   <td style="text-align:right;"> 0.0058346 </td>
   <td style="text-align:right;"> 0.3337175 </td>
   <td style="text-align:right;"> 0.0909391 </td>
   <td style="text-align:right;"> -0.0209686 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 41 </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 74 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P8 </td>
   <td style="text-align:left;"> LV </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 133 </td>
   <td style="text-align:right;"> 100 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 105.0 </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 8.736720 </td>
   <td style="text-align:left;"> R13 </td>
   <td style="text-align:left;"> B1 </td>
   <td style="text-align:left;"> P8 </td>
   <td style="text-align:right;"> -0.0220618 </td>
   <td style="text-align:right;"> 0.0845143 </td>
   <td style="text-align:right;"> 0.0243661 </td>
   <td style="text-align:right;"> -0.0367829 </td>
   <td style="text-align:right;"> 0.0101494 </td>
   <td style="text-align:right;"> -0.0113541 </td>
   <td style="text-align:right;"> -0.0107468 </td>
   <td style="text-align:right;"> -0.0425370 </td>
   <td style="text-align:right;"> -0.0000665 </td>
   <td style="text-align:right;"> -0.0211170 </td>
   <td style="text-align:right;"> -0.0326015 </td>
   <td style="text-align:right;"> -0.0091217 </td>
   <td style="text-align:right;"> 0.0053190 </td>
   <td style="text-align:right;"> -0.0256080 </td>
   <td style="text-align:right;"> -0.0127474 </td>
   <td style="text-align:right;"> -0.0158888 </td>
   <td style="text-align:right;"> -0.0056185 </td>
   <td style="text-align:right;"> -0.0144659 </td>
   <td style="text-align:right;"> -0.0288512 </td>
   <td style="text-align:right;"> 0.0338279 </td>
   <td style="text-align:right;"> 0.0512231 </td>
   <td style="text-align:right;"> 0.0953703 </td>
   <td style="text-align:right;"> -0.0272878 </td>
   <td style="text-align:right;"> 0.0761863 </td>
   <td style="text-align:right;"> -0.0308693 </td>
   <td style="text-align:right;"> -0.0559776 </td>
   <td style="text-align:right;"> 0.0008878 </td>
   <td style="text-align:right;"> 0.0171222 </td>
   <td style="text-align:right;"> -0.0233297 </td>
   <td style="text-align:right;"> -0.0263163 </td>
   <td style="text-align:right;"> 0.0584457 </td>
   <td style="text-align:right;"> -0.0123530 </td>
   <td style="text-align:right;"> -0.0056063 </td>
   <td style="text-align:right;"> -0.0142087 </td>
   <td style="text-align:right;"> 0.0151034 </td>
   <td style="text-align:right;"> 0.0659889 </td>
   <td style="text-align:right;"> 0.1034610 </td>
   <td style="text-align:right;"> 0.0106633 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 84 </td>
   <td style="text-align:right;"> 84 </td>
   <td style="text-align:left;"> FEMALE </td>
   <td style="text-align:right;"> 78 </td>
   <td style="text-align:left;"> S24 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:left;"> P9 </td>
   <td style="text-align:left;"> SV </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> No </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 145 </td>
   <td style="text-align:right;"> 90 </td>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 120.0 </td>
   <td style="text-align:left;"> Yes </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 8.579152 </td>
   <td style="text-align:left;"> R58 </td>
   <td style="text-align:left;"> B5 </td>
   <td style="text-align:left;"> P9 </td>
   <td style="text-align:right;"> 0.0464827 </td>
   <td style="text-align:right;"> -0.0069118 </td>
   <td style="text-align:right;"> 0.0157845 </td>
   <td style="text-align:right;"> -0.0383151 </td>
   <td style="text-align:right;"> 0.0023770 </td>
   <td style="text-align:right;"> 0.1128462 </td>
   <td style="text-align:right;"> 0.0135962 </td>
   <td style="text-align:right;"> 0.0051771 </td>
   <td style="text-align:right;"> 0.0631374 </td>
   <td style="text-align:right;"> 0.0178653 </td>
   <td style="text-align:right;"> -0.0196197 </td>
   <td style="text-align:right;"> 0.0071547 </td>
   <td style="text-align:right;"> 0.0078471 </td>
   <td style="text-align:right;"> -0.0218513 </td>
   <td style="text-align:right;"> 0.0479994 </td>
   <td style="text-align:right;"> 0.0062478 </td>
   <td style="text-align:right;"> 0.0211598 </td>
   <td style="text-align:right;"> 0.0045854 </td>
   <td style="text-align:right;"> 0.0188868 </td>
   <td style="text-align:right;"> 0.0080881 </td>
   <td style="text-align:right;"> 0.1090423 </td>
   <td style="text-align:right;"> -0.0148785 </td>
   <td style="text-align:right;"> 0.0106496 </td>
   <td style="text-align:right;"> -0.0037056 </td>
   <td style="text-align:right;"> 0.0013501 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> 0.0074275 </td>
   <td style="text-align:right;"> 0.0199769 </td>
   <td style="text-align:right;"> 0.0148768 </td>
   <td style="text-align:right;"> 0.0213667 </td>
   <td style="text-align:right;"> 0.1113460 </td>
   <td style="text-align:right;"> -0.0019817 </td>
   <td style="text-align:right;"> 0.0102895 </td>
   <td style="text-align:right;"> 0.0000768 </td>
   <td style="text-align:right;"> 0.0109977 </td>
   <td style="text-align:right;"> NaN </td>
   <td style="text-align:right;"> 0.1761154 </td>
   <td style="text-align:right;"> 0.0167255 </td>
  </tr>
</tbody>
</table></div>

<br>  
<br>  

### Shapiro-Wilk test for normality of each protein log2FC


```r
#Shapiro-Wilks test
p_swnorm <- lapply(lfc_df[, 25:length(lfc_df)], function(x) shapiro.test(x))
results_w <- ldply(p_swnorm, function(x) x$statistic)
results_p <- ldply(p_swnorm, function(x) x$p.value)

#Report the mean and median p-value, the number and 
# percent with p-values greater than 0.05:
mean_FC <- format(mean(results_p$V1), digits=3, format="f")
median_FC <- format(median(results_p$V1), digits=3, format="f")
count_FC <- with(results_p, sum(V1>0.05))
count_tFC <- length(results_p$V1)
pct_normFC <- format((count_FC/count_tFC*100), 
                     digits=3, format="f")

#Making a table of the statistics
results_swnorm <- data.frame(cbind("protein" = results_w$.id, 
                                   "W" = results_w$W, 
                                   "Pvalue" = results_p$V1))
kable(results_swnorm) %>%
  kable_styling() %>%
  scroll_box(width = "75%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:75%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> protein </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> W </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:left;"> 0.928396300106443 </td>
   <td style="text-align:left;"> 0.143817832458436 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:left;"> 0.93877423131274 </td>
   <td style="text-align:left;"> 0.22730073133825 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:left;"> 0.958254829678667 </td>
   <td style="text-align:left;"> 0.50960265060482 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:left;"> 0.938767129452875 </td>
   <td style="text-align:left;"> 0.22723022772362 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:left;"> 0.951629095454959 </td>
   <td style="text-align:left;"> 0.392497878998906 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:left;"> 0.945099083579324 </td>
   <td style="text-align:left;"> 0.450813166286681 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:left;"> 0.981399844106021 </td>
   <td style="text-align:left;"> 0.950864055017285 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:left;"> 0.946467589718173 </td>
   <td style="text-align:left;"> 0.316629077694802 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:left;"> 0.962242480033381 </td>
   <td style="text-align:left;"> 0.589558533798349 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:left;"> 0.944126804552666 </td>
   <td style="text-align:left;"> 0.286582328034959 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:left;"> 0.910261889728923 </td>
   <td style="text-align:left;"> 0.0644719114954409 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:left;"> 0.936535758808858 </td>
   <td style="text-align:left;"> 0.206070146455531 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:left;"> 0.908131070533888 </td>
   <td style="text-align:left;"> 0.0587241491033675 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:left;"> 0.91496860875723 </td>
   <td style="text-align:left;"> 0.0793135153717448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:left;"> 0.968196672510003 </td>
   <td style="text-align:left;"> 0.716429486954957 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:left;"> 0.910908326189805 </td>
   <td style="text-align:left;"> 0.0663280399759093 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:left;"> 0.928009886040205 </td>
   <td style="text-align:left;"> 0.141373420337559 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:left;"> 0.834535956233788 </td>
   <td style="text-align:left;"> 0.00296729214154869 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:left;"> 0.958686375989927 </td>
   <td style="text-align:left;"> 0.51795221833363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:left;"> 0.938355716346468 </td>
   <td style="text-align:left;"> 0.223180777634757 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:left;"> 0.960558465012244 </td>
   <td style="text-align:left;"> 0.555065982303779 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:left;"> 0.993708371236164 </td>
   <td style="text-align:left;"> 0.999937929809321 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:left;"> 0.951915723187482 </td>
   <td style="text-align:left;"> 0.397109109563531 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:left;"> 0.980127266919931 </td>
   <td style="text-align:left;"> 0.935734248310555 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:left;"> 0.903394873408244 </td>
   <td style="text-align:left;"> 0.0477688733079339 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:left;"> 0.875564784240334 </td>
   <td style="text-align:left;"> 0.0179884328598313 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:left;"> 0.976959161440231 </td>
   <td style="text-align:left;"> 0.889110646826791 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:left;"> 0.967742377478776 </td>
   <td style="text-align:left;"> 0.706670128272622 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:left;"> 0.970090070386053 </td>
   <td style="text-align:left;"> 0.756775761142612 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:left;"> 0.957097240780783 </td>
   <td style="text-align:left;"> 0.487614292540579 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:left;"> 0.94363352075278 </td>
   <td style="text-align:left;"> 0.280582697280262 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:left;"> 0.969386245835529 </td>
   <td style="text-align:left;"> 0.741861326535586 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:left;"> 0.963257770929561 </td>
   <td style="text-align:left;"> 0.61077636528134 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:left;"> 0.955082912071624 </td>
   <td style="text-align:left;"> 0.450848122593373 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:left;"> 0.983302008653532 </td>
   <td style="text-align:left;"> 0.969307097573842 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:left;"> 0.944168797465841 </td>
   <td style="text-align:left;"> 0.340949637576333 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:left;"> 0.956660007864635 </td>
   <td style="text-align:left;"> 0.479469089073726 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:left;"> 0.978541543780803 </td>
   <td style="text-align:left;"> 0.913905113760977 </td>
  </tr>
</tbody>
</table></div>

<br>  
For the control and 24 h post stroke protein measures __92.1__ percent had p-values above 0.05; the mean p-value was __0.441__ and median p-value was __0.424__. In regression analysis APOF had the strongest signal, and is normally distributed.

<br> 

### Shapiro-Wilk test for normality of each plasma lipid measure


```r
#Shapiro-Wilks test
plipid <- data.frame("Cholesterol"=lfc_df$Cholesterol, 
                     "Triglycerides"=lfc_df$Triglycerides, 
                     "HDLC"=lfc_df$HDLC, 
                     "LDLC"=lfc_df$LDLC, 
                     "CEC"=lfc_df$CEC)
p_swnorm <- lapply(plipid, function(x) shapiro.test(x))
results_w <- ldply(p_swnorm, function(x) x$statistic)
results_p <- ldply(p_swnorm, function(x) x$p.value)

#Making a table of the statistics
results_swnorm <- data.frame(cbind("protein" = results_w$.id, 
                                   "W" = results_w$W, 
                                   "Pvalue" = results_p$V1))
kable(results_swnorm) %>%
  kable_styling(full_width = FALSE, position = "left")
```

<table class="table" style="width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:left;"> protein </th>
   <th style="text-align:left;"> W </th>
   <th style="text-align:left;"> Pvalue </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Cholesterol </td>
   <td style="text-align:left;"> 0.976258522663825 </td>
   <td style="text-align:left;"> 0.877268941485731 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Triglycerides </td>
   <td style="text-align:left;"> 0.631648804337249 </td>
   <td style="text-align:left;"> 6.37929901629401e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HDLC </td>
   <td style="text-align:left;"> 0.972986873137154 </td>
   <td style="text-align:left;"> 0.816276285560426 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LDLC </td>
   <td style="text-align:left;"> 0.895439304757207 </td>
   <td style="text-align:left;"> 0.0338974795223445 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CEC </td>
   <td style="text-align:left;"> 0.970696413674107 </td>
   <td style="text-align:left;"> 0.769512927970785 </td>
  </tr>
</tbody>
</table>

Both triglycerides and LDLC measures are non-normally distributed.
<br>  
<br>  

## Figure 4. Correlation plot of all patient measures

### Figure 4a. Correlation plot of measures



```r
#subset dataframe for clustering
lfc_clus <- lfc_df

#Not all data columns are useful in a correlation or clustering analysis; e.g. identifiers
lfc_clus <- select(lfc_clus, -"Pair", -"Pair_ID", 
                   -"Stroke", -"Hour", -"Group", 
                   -"l_paired.S24_Pair_ID", -"Sample_no.",
                   -"MS_RunID", -"MS_Batch")

lfc_matrix <- data.matrix(lfc_clus[,1:ncol(lfc_clus)])

# Calculating correlation values
# r-value: Pearson + Benjamini-Hochberg
lFC.pears.BH<-corr.test(lfc_matrix, use = "complete",
                        method="pearson",adjust="BH",alpha=.05)
lFC.pval.BH_pears<-lFC.pears.BH$p
lFC.rval.BH_pears<-lFC.pears.BH$r

#From corrplot package
corrplot(lFC.rval.BH_pears, is.corr=FALSE, order="hclust", 
         tl.col = "black", 
         method = "square", 
         tl.cex = 0.75)
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

### Figure 4b. Hierarchical clustering of measures

Hierarchical clustering with the complete method

```r
#generate pearson correlation matrix
cor.pear<-rcorr(as.matrix(lfc_matrix), type="pearson")
#hierarchical clustering
hclust.x<-hclust((as.dist(1-abs(cor.pear$r))),method="complete")
par(cex=1.2)
plot(hclust.x,ylab = "1-|r|",xlab="",main="")
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

<br>
<br>

# 5. Stroke recovery linear regression analysis

In this section we are examining whether any of the protein changes following stroke show a relationship with patient's stroke recovery - as assessed at 3 mos post stroke with NIHSS score. 

## Table 3. HDL protein logFC from 24h to 96 h post stroke with significant correlation with stroke recovery

<br>  

__Model: NIHSS_3mo ~ NIHSS_baseline + Protein Log2FC__ This model includes the baseline NIHSS score assessed on patient presentation with Stroke.


```r
#This function returns all the statistics for the linear models as a dataframe - used for all downstream lm results
lm_out <- function(fit) {
    if (class(fit) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(fit)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    i <- summary(fit)$coef[1]
    se <- summary(fit)$coef[1,2]
    r2 <- summary(fit)$adj.r.squared
    f <- summary(fit)$fstatistic[1]
    df <- summary(fit)$fstatistic[3]
    out <- cbind(i, se, r2, p, f, df)
    return(out)
    out_df <- rbind(out_df, df)
}

#Initializing an output dataframe to contain all the linear model results
out_df <- data.frame()

#Loop to apply this model to each protein logFC measure
for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + i, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

#renaming & formatting the results dataframe to save
lm_NIH3_baseline.Prot <- out_df
rn <- data.frame(colnames(lfc_df[,25:ncol(lfc_df)]))
rownames(lm_NIH3_baseline.Prot) <- rn[,1]

#multiple testing correction since we are testing 
# the model for all protein log2FC measures
lm_NIH3_baseline.Prot$p.adj <- p.adjust(lm_NIH3_baseline.Prot$p, method = "BH")

count_sig_t4 <- with(lm_NIH3_baseline.Prot, sum(p.adj<0.05))

#Results table
kable(lm_NIH3_baseline.Prot) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> -2.9599877 </td>
   <td style="text-align:right;"> 5.582361 </td>
   <td style="text-align:right;"> 0.2169573 </td>
   <td style="text-align:right;"> 0.0624753 </td>
   <td style="text-align:right;"> 3.355091 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> -2.2255191 </td>
   <td style="text-align:right;"> 5.827706 </td>
   <td style="text-align:right;"> 0.1656413 </td>
   <td style="text-align:right;"> 0.1005689 </td>
   <td style="text-align:right;"> 2.687464 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> -1.9353613 </td>
   <td style="text-align:right;"> 5.661177 </td>
   <td style="text-align:right;"> 0.2057388 </td>
   <td style="text-align:right;"> 0.0695092 </td>
   <td style="text-align:right;"> 3.201769 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -3.6100914 </td>
   <td style="text-align:right;"> 5.871157 </td>
   <td style="text-align:right;"> 0.1799770 </td>
   <td style="text-align:right;"> 0.0883106 </td>
   <td style="text-align:right;"> 2.865563 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 1.8049581 </td>
   <td style="text-align:right;"> 5.126595 </td>
   <td style="text-align:right;"> 0.4132690 </td>
   <td style="text-align:right;"> 0.0071715 </td>
   <td style="text-align:right;"> 6.987048 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0408158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -8.4015070 </td>
   <td style="text-align:right;"> 5.224247 </td>
   <td style="text-align:right;"> 0.5989938 </td>
   <td style="text-align:right;"> 0.0041672 </td>
   <td style="text-align:right;"> 9.962364 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.0408158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> -3.9448568 </td>
   <td style="text-align:right;"> 6.113375 </td>
   <td style="text-align:right;"> 0.1738829 </td>
   <td style="text-align:right;"> 0.0933533 </td>
   <td style="text-align:right;"> 2.789099 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -4.1637643 </td>
   <td style="text-align:right;"> 5.999362 </td>
   <td style="text-align:right;"> 0.1870780 </td>
   <td style="text-align:right;"> 0.0827340 </td>
   <td style="text-align:right;"> 2.956108 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -2.7798924 </td>
   <td style="text-align:right;"> 5.648797 </td>
   <td style="text-align:right;"> 0.1966923 </td>
   <td style="text-align:right;"> 0.0756714 </td>
   <td style="text-align:right;"> 3.081250 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> -1.0406449 </td>
   <td style="text-align:right;"> 5.689091 </td>
   <td style="text-align:right;"> 0.2291305 </td>
   <td style="text-align:right;"> 0.0555487 </td>
   <td style="text-align:right;"> 3.526509 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1055426 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 2.3249182 </td>
   <td style="text-align:right;"> 5.830804 </td>
   <td style="text-align:right;"> 0.3171786 </td>
   <td style="text-align:right;"> 0.0223677 </td>
   <td style="text-align:right;"> 4.948350 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0708310 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 6.1513882 </td>
   <td style="text-align:right;"> 5.933812 </td>
   <td style="text-align:right;"> 0.4109442 </td>
   <td style="text-align:right;"> 0.0073874 </td>
   <td style="text-align:right;"> 6.929871 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0408158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -2.9330815 </td>
   <td style="text-align:right;"> 5.778249 </td>
   <td style="text-align:right;"> 0.1672054 </td>
   <td style="text-align:right;"> 0.0991635 </td>
   <td style="text-align:right;"> 2.706598 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 1.2207997 </td>
   <td style="text-align:right;"> 5.228905 </td>
   <td style="text-align:right;"> 0.3801257 </td>
   <td style="text-align:right;"> 0.0108292 </td>
   <td style="text-align:right;"> 6.212459 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0457231 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -6.5098570 </td>
   <td style="text-align:right;"> 3.321309 </td>
   <td style="text-align:right;"> 0.7336031 </td>
   <td style="text-align:right;"> 0.0000192 </td>
   <td style="text-align:right;"> 24.407274 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0007304 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 0.7782944 </td>
   <td style="text-align:right;"> 4.075043 </td>
   <td style="text-align:right;"> 0.5994117 </td>
   <td style="text-align:right;"> 0.0004098 </td>
   <td style="text-align:right;"> 13.718795 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0077859 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -1.3761293 </td>
   <td style="text-align:right;"> 4.911646 </td>
   <td style="text-align:right;"> 0.3989523 </td>
   <td style="text-align:right;"> 0.0085928 </td>
   <td style="text-align:right;"> 6.641972 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0408158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> -2.6856977 </td>
   <td style="text-align:right;"> 5.795587 </td>
   <td style="text-align:right;"> 0.1551256 </td>
   <td style="text-align:right;"> 0.1104737 </td>
   <td style="text-align:right;"> 2.560666 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> -2.3614623 </td>
   <td style="text-align:right;"> 5.721537 </td>
   <td style="text-align:right;"> 0.1797408 </td>
   <td style="text-align:right;"> 0.0885015 </td>
   <td style="text-align:right;"> 2.862579 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> -0.2979278 </td>
   <td style="text-align:right;"> 5.579077 </td>
   <td style="text-align:right;"> 0.2726919 </td>
   <td style="text-align:right;"> 0.0359092 </td>
   <td style="text-align:right;"> 4.186932 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0813626 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.8616384 </td>
   <td style="text-align:right;"> 5.272122 </td>
   <td style="text-align:right;"> 0.3628831 </td>
   <td style="text-align:right;"> 0.0133033 </td>
   <td style="text-align:right;"> 5.841351 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0505527 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -2.3231699 </td>
   <td style="text-align:right;"> 5.918351 </td>
   <td style="text-align:right;"> 0.1575165 </td>
   <td style="text-align:right;"> 0.1081504 </td>
   <td style="text-align:right;"> 2.589218 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> -0.5633727 </td>
   <td style="text-align:right;"> 6.183559 </td>
   <td style="text-align:right;"> 0.1918381 </td>
   <td style="text-align:right;"> 0.0791690 </td>
   <td style="text-align:right;"> 3.017694 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -2.1079974 </td>
   <td style="text-align:right;"> 5.951323 </td>
   <td style="text-align:right;"> 0.1614365 </td>
   <td style="text-align:right;"> 0.1044329 </td>
   <td style="text-align:right;"> 2.636382 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.7424064 </td>
   <td style="text-align:right;"> 5.676802 </td>
   <td style="text-align:right;"> 0.2411964 </td>
   <td style="text-align:right;"> 0.0493501 </td>
   <td style="text-align:right;"> 3.701845 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0987002 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -0.5900380 </td>
   <td style="text-align:right;"> 7.252770 </td>
   <td style="text-align:right;"> 0.2196481 </td>
   <td style="text-align:right;"> 0.0691974 </td>
   <td style="text-align:right;"> 3.251785 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -1.8546436 </td>
   <td style="text-align:right;"> 5.267034 </td>
   <td style="text-align:right;"> 0.3059877 </td>
   <td style="text-align:right;"> 0.0252680 </td>
   <td style="text-align:right;"> 4.747621 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0738604 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> -1.2058696 </td>
   <td style="text-align:right;"> 5.588165 </td>
   <td style="text-align:right;"> 0.2430234 </td>
   <td style="text-align:right;"> 0.0484658 </td>
   <td style="text-align:right;"> 3.728881 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0987002 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 1.2628836 </td>
   <td style="text-align:right;"> 4.861633 </td>
   <td style="text-align:right;"> 0.4518628 </td>
   <td style="text-align:right;"> 0.0043051 </td>
   <td style="text-align:right;"> 8.007066 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0408158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -4.1962602 </td>
   <td style="text-align:right;"> 5.375333 </td>
   <td style="text-align:right;"> 0.2925652 </td>
   <td style="text-align:right;"> 0.0291720 </td>
   <td style="text-align:right;"> 4.515242 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0791810 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.9214570 </td>
   <td style="text-align:right;"> 5.769297 </td>
   <td style="text-align:right;"> 0.2817232 </td>
   <td style="text-align:right;"> 0.0326969 </td>
   <td style="text-align:right;"> 4.333878 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0813626 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -3.2108399 </td>
   <td style="text-align:right;"> 5.814795 </td>
   <td style="text-align:right;"> 0.1720637 </td>
   <td style="text-align:right;"> 0.0949062 </td>
   <td style="text-align:right;"> 2.766490 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> -1.7636244 </td>
   <td style="text-align:right;"> 4.891178 </td>
   <td style="text-align:right;"> 0.4005318 </td>
   <td style="text-align:right;"> 0.0084249 </td>
   <td style="text-align:right;"> 6.679233 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0408158 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -2.6908737 </td>
   <td style="text-align:right;"> 5.071612 </td>
   <td style="text-align:right;"> 0.3520499 </td>
   <td style="text-align:right;"> 0.0150966 </td>
   <td style="text-align:right;"> 5.618295 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0521519 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> -2.5099197 </td>
   <td style="text-align:right;"> 5.703423 </td>
   <td style="text-align:right;"> 0.1815143 </td>
   <td style="text-align:right;"> 0.0870765 </td>
   <td style="text-align:right;"> 2.885031 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> -6.1525867 </td>
   <td style="text-align:right;"> 7.135180 </td>
   <td style="text-align:right;"> 0.3069299 </td>
   <td style="text-align:right;"> 0.0363991 </td>
   <td style="text-align:right;"> 4.321416 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0813626 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> -2.5403770 </td>
   <td style="text-align:right;"> 5.626489 </td>
   <td style="text-align:right;"> 0.2028285 </td>
   <td style="text-align:right;"> 0.0714423 </td>
   <td style="text-align:right;"> 3.162700 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -2.5704560 </td>
   <td style="text-align:right;"> 5.796046 </td>
   <td style="text-align:right;"> 0.1567161 </td>
   <td style="text-align:right;"> 0.1089234 </td>
   <td style="text-align:right;"> 2.579642 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1104737 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_NIH3_baseline.Prot, file="lm_3mos_bs_Prot.csv")
```

Only __9__ significant correlations are reported in Table 4. All protein statistics are reported in Supplemental Table 5, along with the other linear models examined.
<br>
<br>  


# 6. Cholesterol efflux capacity analysis

## Table 4: Cholesterol Efflux analysis

Cholesterol efflux capacity (CEC) was measured by J774A.1 macrophage based assay, with results recorded as a percentage. The relationship of stroke status and cholesterol efflux capacity are assessed by linear model, with correction for HDL-cholesterol measures and sex, age, and LDL-cholesterol.


```r
fit1 = lm(CEC~Stroke+HDLC,
          data=CvS24)

fit2 = lm(CEC~Stroke+Sex+Age_blood+HDLC+LDLC,
          data=CvS24)

#html table with both CEC model statistics
stargazer(fit1, fit2,
          title="Linear Regression Results",
          report="vcsp*",ci=TRUE,
          keep.stat=c("n","rsq","adj.rsq","aic","ser","f"),
          single.row = FALSE,type="html")
```


<table style="text-align:center"><caption><strong>Linear Regression Results</strong></caption>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"></td><td colspan="2"><em>Dependent variable:</em></td></tr>
<tr><td></td><td colspan="2" style="border-bottom: 1px solid black"></td></tr>
<tr><td style="text-align:left"></td><td colspan="2">CEC</td></tr>
<tr><td style="text-align:left"></td><td>(1)</td><td>(2)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Stroke</td><td>-3.435</td><td>-2.450</td></tr>
<tr><td style="text-align:left"></td><td>(-5.330, -1.540)</td><td>(-4.344, -0.555)</td></tr>
<tr><td style="text-align:left"></td><td>p = 0.001<sup>***</sup></td><td>p = 0.014<sup>**</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">SexMALE</td><td></td><td>1.709</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(0.015, 3.404)</td></tr>
<tr><td style="text-align:left"></td><td></td><td>p = 0.053<sup>*</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">SexNA</td><td></td><td>-9.036</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-15.134, -2.939)</td></tr>
<tr><td style="text-align:left"></td><td></td><td>p = 0.006<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Age_blood</td><td></td><td>0.002</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(-0.092, 0.097)</td></tr>
<tr><td style="text-align:left"></td><td></td><td>p = 0.963</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">HDLC</td><td>0.076</td><td>0.088</td></tr>
<tr><td style="text-align:left"></td><td>(0.031, 0.121)</td><td>(0.043, 0.132)</td></tr>
<tr><td style="text-align:left"></td><td>p = 0.002<sup>***</sup></td><td>p = 0.0003<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">LDLC</td><td></td><td>0.030</td></tr>
<tr><td style="text-align:left"></td><td></td><td>(0.015, 0.045)</td></tr>
<tr><td style="text-align:left"></td><td></td><td>p = 0.0003<sup>***</sup></td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td style="text-align:left">Constant</td><td>10.717</td><td>4.185</td></tr>
<tr><td style="text-align:left"></td><td>(7.619, 13.815)</td><td>(-3.677, 12.047)</td></tr>
<tr><td style="text-align:left"></td><td>p = 0.000<sup>***</sup></td><td>p = 0.302</td></tr>
<tr><td style="text-align:left"></td><td></td><td></td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left">Observations</td><td>70</td><td>66</td></tr>
<tr><td style="text-align:left">R<sup>2</sup></td><td>0.409</td><td>0.556</td></tr>
<tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td>0.391</td><td>0.511</td></tr>
<tr><td style="text-align:left">Residual Std. Error</td><td>3.521 (df = 67)</td><td>3.001 (df = 59)</td></tr>
<tr><td style="text-align:left">F Statistic</td><td>23.172<sup>***</sup> (df = 2; 67)</td><td>12.317<sup>***</sup> (df = 6; 59)</td></tr>
<tr><td colspan="3" style="border-bottom: 1px solid black"></td></tr><tr><td style="text-align:left"><em>Note:</em></td><td colspan="2" style="text-align:right"><sup>*</sup>p<0.1; <sup>**</sup>p<0.05; <sup>***</sup>p<0.01</td></tr>
</table>

<br>  
<br>

## Figure 5. Cholesterol efflux distribution plot 

Plotting the percent CEC distributions for each group.


```r
#this plot function is found in section 3 for generation of Figure 2
allViolinplot(all, all$CEC) + labs(y="cholesterol efflux capacity (%)")
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

<br>
<br>

# 7. Supplemental material


## Supplemental Figure 1. Effect of normalization on N15 APOA1 spike-in variance distribution

For this figure I used a Skyline generated .csv report with peptide integrated areas, background signal, sample total ion current (TIC), TIC mean normalization factors and values, and heavy APOA1 peptide mean normalization factors and values.

<br>   

In this dataset two different normalization strategies were assessed. Normalization of peptide intensity to total ion current was determined to be a better method for normalization, with it reducing the variance across system suitability replicates. These replicates are also called rotor control samples due to the nature of HDL isolation being limited to rotor batches for density gradient ultra-centrifugation.


```r
#This is a modified skyline .csv report file. 
#These rotor control samples are only analyzed here to demonstrate that the- 
#mass spec data normalization improved data quality.
RCdf <- as.data.frame(read.csv("RotorControl_PRM.csv"))

#making a new identifier that incorporates both the peptide sequence and whether it was light or heavy
RCdf$PepID =paste(RCdf$Peptide, RCdf$Isotope, sep="_")

RCdf_CV <- select(RCdf, 
                  'PepID', 'Background.subtracted', 
                  'TIC_normalized_area', 
                  'APOA1_normalized_area')
RCdf_CV <- RCdf_CV %>% group_by(PepID) %>% 
                          summarise_each(funs(mean, sd))

#just a simple CV function
CV_fun <- function(mean, sd){
  (sd/mean)*100
}

#Calculating different CVs
RCdf_CV$Background.subtracted_cv <- CV_fun(RCdf_CV$Background.subtracted_mean, 
                                           RCdf_CV$Background.subtracted_sd)
RCdf_CV$TIC_normalized_area_cv <- CV_fun(RCdf_CV$TIC_normalized_area_mean, 
                                         RCdf_CV$TIC_normalized_area_sd)

#Mean and median CVs across all peptides
Pre_mean = mean(RCdf_CV$Background.subtracted_cv)
Pre_med = median(RCdf_CV$Background.subtracted_cv)
prem <- format(median(RCdf_CV$Background.subtracted_cv), digits=3, format="f")
TIC_mean = mean(RCdf_CV$TIC_normalized_area_cv)
TIC_med = median(RCdf_CV$TIC_normalized_area_cv)
TICm <- format(median(RCdf_CV$TIC_normalized_area_cv), digits=3, format="f")

#density plots for CV distribution
ggplot(RCdf_CV) +
  geom_density(aes(x=Background.subtracted_cv), 
               colour='black', fill='orange', alpha=0.4) +
  geom_density(aes(x=TIC_normalized_area_cv), 
               colour='black', fill='dodgerblue3', alpha=0.4) +
  geom_vline(xintercept=TIC_med, color='dodgerblue3', 
             linetype='longdash', size=1) +
  geom_vline(xintercept=Pre_med, color='orange', 
             linetype='dashed', size=1) +
  geom_text(aes(x=TIC_med, label="Post median\n", y=0.008), 
            colour="dodgerblue3", angle=90, 
            text=element_text(size=10)) +
  geom_text(aes(x=Pre_med, label="Pre median\n", y=0.008), 
            colour="goldenrod1", angle=90, 
            text=element_text(size=10)) +
  labs(x="Peptide %CV", y="Density") + 
  ggtitle("Rotor control peptide %CV - Pre/Post-TIC normalization") +
  theme_bw() +
  scale_x_continuous(limits = c(0,55))
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

The median %CV prior to normalization was 21.5, reduced to 16.1 following total ion current normalization. Typically for targeted peptide assays CVs < 20% are preferred. 


## Supplemental Table 4. Correlation results matrices

The p values and r values corresponding to the large correlation plot in Figure 4 are compiled and written to .csv. Both the following tables were combined into a single excel workbook.

__Pearson correlation r statistic matrix__

```r
#write.csv(lFC.pval.BH_pears, "lFC_pvalBH_pears.csv", row.names=FALSE)
#write.csv(lFC.rval.BH_pears, "lFC_rvalBH_pears.csv", row.names=FALSE)

kable(lFC.pval.BH_pears) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Sex </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Age_blood </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ischemic_Type </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> tPa </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Thrombectomy </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> TICI2B </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Cholesterol </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Triglycerides </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> HDLC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LDLC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Statin </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> NIHSS_baseline </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> NIHSS_3mo </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> mRS </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CEC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SERPINA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ALB </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AMBP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ANTXR2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APMAP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LPA </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOB </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOC3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOC4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOD </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOE </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOF </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOL1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOM </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CAMP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CLU </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> C3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PPBP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AHSG </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> FGA </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> HPR </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> IHH </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ITIH4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LCAT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SELL </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PCYOX1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> GPLD1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PF4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PLTP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PON1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PON3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> RBP4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SAA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SAA2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> TTR </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Sex </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.4772483 </td>
   <td style="text-align:right;"> 0.9145204 </td>
   <td style="text-align:right;"> 0.9145204 </td>
   <td style="text-align:right;"> 0.9727465 </td>
   <td style="text-align:right;"> 0.6977462 </td>
   <td style="text-align:right;"> 0.6212312 </td>
   <td style="text-align:right;"> 0.6281754 </td>
   <td style="text-align:right;"> 0.5154820 </td>
   <td style="text-align:right;"> 0.4772483 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.3040382 </td>
   <td style="text-align:right;"> 0.6914639 </td>
   <td style="text-align:right;"> 0.2449623 </td>
   <td style="text-align:right;"> 0.3040562 </td>
   <td style="text-align:right;"> 0.0679049 </td>
   <td style="text-align:right;"> 0.1061220 </td>
   <td style="text-align:right;"> 0.7364382 </td>
   <td style="text-align:right;"> 0.1955318 </td>
   <td style="text-align:right;"> 0.6043755 </td>
   <td style="text-align:right;"> 0.4217927 </td>
   <td style="text-align:right;"> 0.6591875 </td>
   <td style="text-align:right;"> 0.4007046 </td>
   <td style="text-align:right;"> 0.7448671 </td>
   <td style="text-align:right;"> 0.3995072 </td>
   <td style="text-align:right;"> 0.2336319 </td>
   <td style="text-align:right;"> 0.3307208 </td>
   <td style="text-align:right;"> 0.3058270 </td>
   <td style="text-align:right;"> 0.8915819 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.1821335 </td>
   <td style="text-align:right;"> 0.4631920 </td>
   <td style="text-align:right;"> 0.5253909 </td>
   <td style="text-align:right;"> 0.1481069 </td>
   <td style="text-align:right;"> 0.3105036 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.2493874 </td>
   <td style="text-align:right;"> 0.2238112 </td>
   <td style="text-align:right;"> 0.4221956 </td>
   <td style="text-align:right;"> 0.9139271 </td>
   <td style="text-align:right;"> 0.0997341 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.3315460 </td>
   <td style="text-align:right;"> 0.2635613 </td>
   <td style="text-align:right;"> 0.1952102 </td>
   <td style="text-align:right;"> 0.1147992 </td>
   <td style="text-align:right;"> 0.1914335 </td>
   <td style="text-align:right;"> 0.2160365 </td>
   <td style="text-align:right;"> 0.6794488 </td>
   <td style="text-align:right;"> 0.0967334 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age_blood </td>
   <td style="text-align:right;"> 0.6843132 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.8964113 </td>
   <td style="text-align:right;"> 0.1219558 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.7769474 </td>
   <td style="text-align:right;"> 0.9894492 </td>
   <td style="text-align:right;"> 0.8789393 </td>
   <td style="text-align:right;"> 0.2023322 </td>
   <td style="text-align:right;"> 0.7865122 </td>
   <td style="text-align:right;"> 0.5639717 </td>
   <td style="text-align:right;"> 0.6496004 </td>
   <td style="text-align:right;"> 0.9875093 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.6607688 </td>
   <td style="text-align:right;"> 0.9070175 </td>
   <td style="text-align:right;"> 0.3342780 </td>
   <td style="text-align:right;"> 0.6913284 </td>
   <td style="text-align:right;"> 0.8924880 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.7547639 </td>
   <td style="text-align:right;"> 0.3135264 </td>
   <td style="text-align:right;"> 0.6634940 </td>
   <td style="text-align:right;"> 0.9471271 </td>
   <td style="text-align:right;"> 0.9894492 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.9776100 </td>
   <td style="text-align:right;"> 0.6043755 </td>
   <td style="text-align:right;"> 0.7801209 </td>
   <td style="text-align:right;"> 0.8301618 </td>
   <td style="text-align:right;"> 0.8164696 </td>
   <td style="text-align:right;"> 0.7575202 </td>
   <td style="text-align:right;"> 0.7485384 </td>
   <td style="text-align:right;"> 0.6573562 </td>
   <td style="text-align:right;"> 0.9925121 </td>
   <td style="text-align:right;"> 0.7936401 </td>
   <td style="text-align:right;"> 0.6033733 </td>
   <td style="text-align:right;"> 0.9484162 </td>
   <td style="text-align:right;"> 0.7994115 </td>
   <td style="text-align:right;"> 0.6253230 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.9989164 </td>
   <td style="text-align:right;"> 0.6561511 </td>
   <td style="text-align:right;"> 0.9684830 </td>
   <td style="text-align:right;"> 0.4738537 </td>
   <td style="text-align:right;"> 0.7786512 </td>
   <td style="text-align:right;"> 0.4645421 </td>
   <td style="text-align:right;"> 0.9480859 </td>
   <td style="text-align:right;"> 0.9160070 </td>
   <td style="text-align:right;"> 0.8061277 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.5588577 </td>
   <td style="text-align:right;"> 0.6062824 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ischemic_Type </td>
   <td style="text-align:right;"> 0.2320438 </td>
   <td style="text-align:right;"> 0.7829887 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.9456724 </td>
   <td style="text-align:right;"> 0.9456724 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.3603448 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.9544699 </td>
   <td style="text-align:right;"> 0.4373454 </td>
   <td style="text-align:right;"> 0.6943273 </td>
   <td style="text-align:right;"> 0.3234836 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.8150002 </td>
   <td style="text-align:right;"> 0.8941359 </td>
   <td style="text-align:right;"> 0.8915819 </td>
   <td style="text-align:right;"> 0.6659844 </td>
   <td style="text-align:right;"> 0.6252458 </td>
   <td style="text-align:right;"> 0.9853451 </td>
   <td style="text-align:right;"> 0.6659844 </td>
   <td style="text-align:right;"> 0.4104424 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.9192539 </td>
   <td style="text-align:right;"> 0.4873797 </td>
   <td style="text-align:right;"> 0.9570606 </td>
   <td style="text-align:right;"> 0.5253909 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.6548282 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.7840744 </td>
   <td style="text-align:right;"> 0.2592034 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.9544699 </td>
   <td style="text-align:right;"> 0.9044865 </td>
   <td style="text-align:right;"> 0.6008058 </td>
   <td style="text-align:right;"> 0.7952425 </td>
   <td style="text-align:right;"> 0.7994115 </td>
   <td style="text-align:right;"> 0.5347501 </td>
   <td style="text-align:right;"> 0.4828866 </td>
   <td style="text-align:right;"> 0.4408199 </td>
   <td style="text-align:right;"> 0.4701513 </td>
   <td style="text-align:right;"> 0.7273135 </td>
   <td style="text-align:right;"> 0.9151110 </td>
   <td style="text-align:right;"> 0.9875093 </td>
   <td style="text-align:right;"> 0.8941359 </td>
   <td style="text-align:right;"> 0.9480859 </td>
   <td style="text-align:right;"> 0.7663131 </td>
   <td style="text-align:right;"> 0.6607688 </td>
   <td style="text-align:right;"> 0.4772483 </td>
   <td style="text-align:right;"> 0.7943741 </td>
   <td style="text-align:right;"> 0.6164267 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tPa </td>
   <td style="text-align:right;"> 0.8227557 </td>
   <td style="text-align:right;"> 0.0126695 </td>
   <td style="text-align:right;"> 0.8889862 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6281754 </td>
   <td style="text-align:right;"> 0.2743478 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.5107596 </td>
   <td style="text-align:right;"> 0.2418292 </td>
   <td style="text-align:right;"> 0.3386730 </td>
   <td style="text-align:right;"> 0.9456724 </td>
   <td style="text-align:right;"> 0.7576908 </td>
   <td style="text-align:right;"> 0.4286568 </td>
   <td style="text-align:right;"> 0.3097152 </td>
   <td style="text-align:right;"> 0.6253230 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.4595219 </td>
   <td style="text-align:right;"> 0.6008058 </td>
   <td style="text-align:right;"> 0.9411736 </td>
   <td style="text-align:right;"> 0.7051462 </td>
   <td style="text-align:right;"> 0.7848482 </td>
   <td style="text-align:right;"> 0.4772483 </td>
   <td style="text-align:right;"> 0.6474421 </td>
   <td style="text-align:right;"> 0.5665047 </td>
   <td style="text-align:right;"> 0.7357528 </td>
   <td style="text-align:right;"> 0.8741703 </td>
   <td style="text-align:right;"> 0.8941359 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.7952425 </td>
   <td style="text-align:right;"> 0.3394980 </td>
   <td style="text-align:right;"> 0.6803704 </td>
   <td style="text-align:right;"> 0.7313016 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.9998409 </td>
   <td style="text-align:right;"> 0.3078751 </td>
   <td style="text-align:right;"> 0.5342046 </td>
   <td style="text-align:right;"> 0.6470452 </td>
   <td style="text-align:right;"> 0.2700655 </td>
   <td style="text-align:right;"> 0.6607688 </td>
   <td style="text-align:right;"> 0.6166331 </td>
   <td style="text-align:right;"> 0.9913579 </td>
   <td style="text-align:right;"> 0.5243703 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.9954592 </td>
   <td style="text-align:right;"> 0.3699431 </td>
   <td style="text-align:right;"> 0.6050853 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.9145204 </td>
   <td style="text-align:right;"> 0.8220028 </td>
   <td style="text-align:right;"> 0.4392507 </td>
   <td style="text-align:right;"> 0.9831447 </td>
   <td style="text-align:right;"> 0.4595219 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thrombectomy </td>
   <td style="text-align:right;"> 0.8227557 </td>
   <td style="text-align:right;"> 0.9250160 </td>
   <td style="text-align:right;"> 0.8889862 </td>
   <td style="text-align:right;"> 0.3784767 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0018003 </td>
   <td style="text-align:right;"> 0.6144991 </td>
   <td style="text-align:right;"> 0.7028201 </td>
   <td style="text-align:right;"> 0.3528014 </td>
   <td style="text-align:right;"> 0.7357528 </td>
   <td style="text-align:right;"> 0.9456724 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.7400180 </td>
   <td style="text-align:right;"> 0.7631012 </td>
   <td style="text-align:right;"> 0.6196941 </td>
   <td style="text-align:right;"> 0.9544699 </td>
   <td style="text-align:right;"> 0.7547639 </td>
   <td style="text-align:right;"> 0.7631012 </td>
   <td style="text-align:right;"> 0.5342046 </td>
   <td style="text-align:right;"> 0.3209503 </td>
   <td style="text-align:right;"> 0.3514749 </td>
   <td style="text-align:right;"> 0.6252458 </td>
   <td style="text-align:right;"> 0.7130048 </td>
   <td style="text-align:right;"> 0.4613459 </td>
   <td style="text-align:right;"> 0.5243703 </td>
   <td style="text-align:right;"> 0.8078867 </td>
   <td style="text-align:right;"> 0.4796402 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.9198142 </td>
   <td style="text-align:right;"> 0.6668826 </td>
   <td style="text-align:right;"> 0.8581986 </td>
   <td style="text-align:right;"> 0.7448671 </td>
   <td style="text-align:right;"> 0.9380902 </td>
   <td style="text-align:right;"> 0.9112983 </td>
   <td style="text-align:right;"> 0.4093059 </td>
   <td style="text-align:right;"> 0.6573562 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.8016897 </td>
   <td style="text-align:right;"> 0.1914335 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.8960037 </td>
   <td style="text-align:right;"> 0.9121943 </td>
   <td style="text-align:right;"> 0.7181992 </td>
   <td style="text-align:right;"> 0.8164696 </td>
   <td style="text-align:right;"> 0.6005902 </td>
   <td style="text-align:right;"> 0.5861017 </td>
   <td style="text-align:right;"> 0.9192539 </td>
   <td style="text-align:right;"> 0.9352347 </td>
   <td style="text-align:right;"> 0.9539938 </td>
   <td style="text-align:right;"> 0.9161951 </td>
   <td style="text-align:right;"> 0.6606778 </td>
   <td style="text-align:right;"> 0.6144991 </td>
   <td style="text-align:right;"> 0.9471271 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TICI2B </td>
   <td style="text-align:right;"> 0.9484984 </td>
   <td style="text-align:right;"> 0.5822534 </td>
   <td style="text-align:right;"> 0.6269262 </td>
   <td style="text-align:right;"> 0.0816021 </td>
   <td style="text-align:right;"> 0.0000065 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4840927 </td>
   <td style="text-align:right;"> 0.3812442 </td>
   <td style="text-align:right;"> 0.2578932 </td>
   <td style="text-align:right;"> 0.5946005 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.8676437 </td>
   <td style="text-align:right;"> 0.9219944 </td>
   <td style="text-align:right;"> 0.4833604 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.9875093 </td>
   <td style="text-align:right;"> 0.7853940 </td>
   <td style="text-align:right;"> 0.9192539 </td>
   <td style="text-align:right;"> 0.6666071 </td>
   <td style="text-align:right;"> 0.5620014 </td>
   <td style="text-align:right;"> 0.4978452 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.6496004 </td>
   <td style="text-align:right;"> 0.7692291 </td>
   <td style="text-align:right;"> 0.8536794 </td>
   <td style="text-align:right;"> 0.5618687 </td>
   <td style="text-align:right;"> 0.7200398 </td>
   <td style="text-align:right;"> 0.7865122 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.8619808 </td>
   <td style="text-align:right;"> 0.9380902 </td>
   <td style="text-align:right;"> 0.8046712 </td>
   <td style="text-align:right;"> 0.7786512 </td>
   <td style="text-align:right;"> 0.5479784 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.8061277 </td>
   <td style="text-align:right;"> 0.1289304 </td>
   <td style="text-align:right;"> 0.4833604 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.6679530 </td>
   <td style="text-align:right;"> 0.6285719 </td>
   <td style="text-align:right;"> 0.6913284 </td>
   <td style="text-align:right;"> 0.7126080 </td>
   <td style="text-align:right;"> 0.9456724 </td>
   <td style="text-align:right;"> 0.8520815 </td>
   <td style="text-align:right;"> 0.8448507 </td>
   <td style="text-align:right;"> 0.9129462 </td>
   <td style="text-align:right;"> 0.7313016 </td>
   <td style="text-align:right;"> 0.6043755 </td>
   <td style="text-align:right;"> 0.9629801 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cholesterol </td>
   <td style="text-align:right;"> 0.4698901 </td>
   <td style="text-align:right;"> 0.9814642 </td>
   <td style="text-align:right;"> 0.1363445 </td>
   <td style="text-align:right;"> 0.0512695 </td>
   <td style="text-align:right;"> 0.3627854 </td>
   <td style="text-align:right;"> 0.2406412 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.7308369 </td>
   <td style="text-align:right;"> 0.9897566 </td>
   <td style="text-align:right;"> 0.0001508 </td>
   <td style="text-align:right;"> 0.7348579 </td>
   <td style="text-align:right;"> 0.7994115 </td>
   <td style="text-align:right;"> 0.6212312 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.4185564 </td>
   <td style="text-align:right;"> 0.8301618 </td>
   <td style="text-align:right;"> 0.8915819 </td>
   <td style="text-align:right;"> 0.9151110 </td>
   <td style="text-align:right;"> 0.4757301 </td>
   <td style="text-align:right;"> 0.4039211 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.9800795 </td>
   <td style="text-align:right;"> 0.6841514 </td>
   <td style="text-align:right;"> 0.4373688 </td>
   <td style="text-align:right;"> 0.9335152 </td>
   <td style="text-align:right;"> 0.9450891 </td>
   <td style="text-align:right;"> 0.8915819 </td>
   <td style="text-align:right;"> 0.3011440 </td>
   <td style="text-align:right;"> 0.9126519 </td>
   <td style="text-align:right;"> 0.6120971 </td>
   <td style="text-align:right;"> 0.2643314 </td>
   <td style="text-align:right;"> 0.4161712 </td>
   <td style="text-align:right;"> 0.7181992 </td>
   <td style="text-align:right;"> 0.6607688 </td>
   <td style="text-align:right;"> 0.3307208 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.2621030 </td>
   <td style="text-align:right;"> 0.7313016 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.6573562 </td>
   <td style="text-align:right;"> 0.8446452 </td>
   <td style="text-align:right;"> 0.0142002 </td>
   <td style="text-align:right;"> 0.5090339 </td>
   <td style="text-align:right;"> 0.2803367 </td>
   <td style="text-align:right;"> 0.2926542 </td>
   <td style="text-align:right;"> 0.4208246 </td>
   <td style="text-align:right;"> 0.6187062 </td>
   <td style="text-align:right;"> 0.6380263 </td>
   <td style="text-align:right;"> 0.8082713 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.9311527 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Triglycerides </td>
   <td style="text-align:right;"> 0.3705748 </td>
   <td style="text-align:right;"> 0.7496925 </td>
   <td style="text-align:right;"> 0.9255713 </td>
   <td style="text-align:right;"> 0.2616856 </td>
   <td style="text-align:right;"> 0.4743271 </td>
   <td style="text-align:right;"> 0.1501920 </td>
   <td style="text-align:right;"> 0.5175660 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.5291199 </td>
   <td style="text-align:right;"> 0.5036602 </td>
   <td style="text-align:right;"> 0.2897794 </td>
   <td style="text-align:right;"> 0.8520815 </td>
   <td style="text-align:right;"> 0.2449623 </td>
   <td style="text-align:right;"> 0.4701513 </td>
   <td style="text-align:right;"> 0.6903275 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.9139271 </td>
   <td style="text-align:right;"> 0.9112983 </td>
   <td style="text-align:right;"> 0.9636034 </td>
   <td style="text-align:right;"> 0.9146012 </td>
   <td style="text-align:right;"> 0.6140420 </td>
   <td style="text-align:right;"> 0.7400494 </td>
   <td style="text-align:right;"> 0.8646012 </td>
   <td style="text-align:right;"> 0.4312296 </td>
   <td style="text-align:right;"> 0.9151110 </td>
   <td style="text-align:right;"> 0.6806411 </td>
   <td style="text-align:right;"> 0.9121943 </td>
   <td style="text-align:right;"> 0.6668826 </td>
   <td style="text-align:right;"> 0.8480176 </td>
   <td style="text-align:right;"> 0.3986824 </td>
   <td style="text-align:right;"> 0.5709317 </td>
   <td style="text-align:right;"> 0.4733762 </td>
   <td style="text-align:right;"> 0.4470100 </td>
   <td style="text-align:right;"> 0.2752137 </td>
   <td style="text-align:right;"> 0.9069815 </td>
   <td style="text-align:right;"> 0.4828866 </td>
   <td style="text-align:right;"> 0.9311527 </td>
   <td style="text-align:right;"> 0.7881789 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.3758418 </td>
   <td style="text-align:right;"> 0.7066730 </td>
   <td style="text-align:right;"> 0.6506594 </td>
   <td style="text-align:right;"> 0.6967520 </td>
   <td style="text-align:right;"> 0.6066109 </td>
   <td style="text-align:right;"> 0.9954592 </td>
   <td style="text-align:right;"> 0.7042441 </td>
   <td style="text-align:right;"> 0.8588721 </td>
   <td style="text-align:right;"> 0.4735637 </td>
   <td style="text-align:right;"> 0.2622488 </td>
   <td style="text-align:right;"> 0.7313016 </td>
   <td style="text-align:right;"> 0.7590379 </td>
   <td style="text-align:right;"> 0.7212019 </td>
   <td style="text-align:right;"> 0.9380902 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HDLC </td>
   <td style="text-align:right;"> 0.3790003 </td>
   <td style="text-align:right;"> 0.0423289 </td>
   <td style="text-align:right;"> 0.9107748 </td>
   <td style="text-align:right;"> 0.0617735 </td>
   <td style="text-align:right;"> 0.1318525 </td>
   <td style="text-align:right;"> 0.0688713 </td>
   <td style="text-align:right;"> 0.9840106 </td>
   <td style="text-align:right;"> 0.2787671 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.7588192 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.8066668 </td>
   <td style="text-align:right;"> 0.6741025 </td>
   <td style="text-align:right;"> 0.4911550 </td>
   <td style="text-align:right;"> 0.3307208 </td>
   <td style="text-align:right;"> 0.1885416 </td>
   <td style="text-align:right;"> 0.2165099 </td>
   <td style="text-align:right;"> 0.8061277 </td>
   <td style="text-align:right;"> 0.4833604 </td>
   <td style="text-align:right;"> 0.2323097 </td>
   <td style="text-align:right;"> 0.7066730 </td>
   <td style="text-align:right;"> 0.8964113 </td>
   <td style="text-align:right;"> 0.9347131 </td>
   <td style="text-align:right;"> 0.1692739 </td>
   <td style="text-align:right;"> 0.7531198 </td>
   <td style="text-align:right;"> 0.5767151 </td>
   <td style="text-align:right;"> 0.4143197 </td>
   <td style="text-align:right;"> 0.3304656 </td>
   <td style="text-align:right;"> 0.9826077 </td>
   <td style="text-align:right;"> 0.3577065 </td>
   <td style="text-align:right;"> 0.3461724 </td>
   <td style="text-align:right;"> 0.9724767 </td>
   <td style="text-align:right;"> 0.1755108 </td>
   <td style="text-align:right;"> 0.1759068 </td>
   <td style="text-align:right;"> 0.9649381 </td>
   <td style="text-align:right;"> 0.2126416 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.1460815 </td>
   <td style="text-align:right;"> 0.9798903 </td>
   <td style="text-align:right;"> 0.8552529 </td>
   <td style="text-align:right;"> 0.6281754 </td>
   <td style="text-align:right;"> 0.5742916 </td>
   <td style="text-align:right;"> 0.5618687 </td>
   <td style="text-align:right;"> 0.2611891 </td>
   <td style="text-align:right;"> 0.9839432 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.7481592 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.3986824 </td>
   <td style="text-align:right;"> 0.0608425 </td>
   <td style="text-align:right;"> 0.5690653 </td>
   <td style="text-align:right;"> 0.2102899 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LDLC </td>
   <td style="text-align:right;"> 0.2667189 </td>
   <td style="text-align:right;"> 0.6010140 </td>
   <td style="text-align:right;"> 0.1986780 </td>
   <td style="text-align:right;"> 0.1242421 </td>
   <td style="text-align:right;"> 0.5290364 </td>
   <td style="text-align:right;"> 0.3361348 </td>
   <td style="text-align:right;"> 0.0000002 </td>
   <td style="text-align:right;"> 0.2565816 </td>
   <td style="text-align:right;"> 0.5592689 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.8941359 </td>
   <td style="text-align:right;"> 0.8207976 </td>
   <td style="text-align:right;"> 0.3985795 </td>
   <td style="text-align:right;"> 0.9314912 </td>
   <td style="text-align:right;"> 0.5270490 </td>
   <td style="text-align:right;"> 0.9544699 </td>
   <td style="text-align:right;"> 0.6573562 </td>
   <td style="text-align:right;"> 0.8792609 </td>
   <td style="text-align:right;"> 0.4772483 </td>
   <td style="text-align:right;"> 0.3078751 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.5915043 </td>
   <td style="text-align:right;"> 0.9544699 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.2609659 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.8789393 </td>
   <td style="text-align:right;"> 0.9471271 </td>
   <td style="text-align:right;"> 0.1898696 </td>
   <td style="text-align:right;"> 0.7846594 </td>
   <td style="text-align:right;"> 0.3603448 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.6152790 </td>
   <td style="text-align:right;"> 0.7037227 </td>
   <td style="text-align:right;"> 0.4833604 </td>
   <td style="text-align:right;"> 0.2803367 </td>
   <td style="text-align:right;"> 0.6599151 </td>
   <td style="text-align:right;"> 0.3009524 </td>
   <td style="text-align:right;"> 0.9139271 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.5690653 </td>
   <td style="text-align:right;"> 0.6711661 </td>
   <td style="text-align:right;"> 0.0030454 </td>
   <td style="text-align:right;"> 0.3357679 </td>
   <td style="text-align:right;"> 0.1914335 </td>
   <td style="text-align:right;"> 0.2824454 </td>
   <td style="text-align:right;"> 0.2650631 </td>
   <td style="text-align:right;"> 0.4370291 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.6043755 </td>
   <td style="text-align:right;"> 0.3632971 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.7576908 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Statin </td>
   <td style="text-align:right;"> 0.2320438 </td>
   <td style="text-align:right;"> 0.3098161 </td>
   <td style="text-align:right;"> 0.4655722 </td>
   <td style="text-align:right;"> 0.8889862 </td>
   <td style="text-align:right;"> 0.8889862 </td>
   <td style="text-align:right;"> 0.6269262 </td>
   <td style="text-align:right;"> 0.5258127 </td>
   <td style="text-align:right;"> 0.0897938 </td>
   <td style="text-align:right;"> 0.6906961 </td>
   <td style="text-align:right;"> 0.7721493 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.9411736 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.9348425 </td>
   <td style="text-align:right;"> 0.7631012 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.4101880 </td>
   <td style="text-align:right;"> 0.2743478 </td>
   <td style="text-align:right;"> 0.8588721 </td>
   <td style="text-align:right;"> 0.4286568 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.7183016 </td>
   <td style="text-align:right;"> 0.9161951 </td>
   <td style="text-align:right;"> 0.8016897 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.9544699 </td>
   <td style="text-align:right;"> 0.2743478 </td>
   <td style="text-align:right;"> 0.9798903 </td>
   <td style="text-align:right;"> 0.8907941 </td>
   <td style="text-align:right;"> 0.9358488 </td>
   <td style="text-align:right;"> 0.8177333 </td>
   <td style="text-align:right;"> 0.5773918 </td>
   <td style="text-align:right;"> 0.5716226 </td>
   <td style="text-align:right;"> 0.6782650 </td>
   <td style="text-align:right;"> 0.9831447 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.7575202 </td>
   <td style="text-align:right;"> 0.4569638 </td>
   <td style="text-align:right;"> 0.1814966 </td>
   <td style="text-align:right;"> 0.9160070 </td>
   <td style="text-align:right;"> 0.3577065 </td>
   <td style="text-align:right;"> 0.7308369 </td>
   <td style="text-align:right;"> 0.8177333 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.9097531 </td>
   <td style="text-align:right;"> 0.7994115 </td>
   <td style="text-align:right;"> 0.7785247 </td>
   <td style="text-align:right;"> 0.2964747 </td>
   <td style="text-align:right;"> 0.2652301 </td>
   <td style="text-align:right;"> 0.9660226 </td>
   <td style="text-align:right;"> 0.7183016 </td>
   <td style="text-align:right;"> 0.5448245 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NIHSS_baseline </td>
   <td style="text-align:right;"> 0.0069233 </td>
   <td style="text-align:right;"> 0.4063538 </td>
   <td style="text-align:right;"> 0.1103319 </td>
   <td style="text-align:right;"> 0.5569962 </td>
   <td style="text-align:right;"> 0.4793281 </td>
   <td style="text-align:right;"> 0.7359541 </td>
   <td style="text-align:right;"> 0.6186389 </td>
   <td style="text-align:right;"> 0.7124039 </td>
   <td style="text-align:right;"> 0.6911869 </td>
   <td style="text-align:right;"> 0.6587824 </td>
   <td style="text-align:right;"> 0.8803867 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.2643314 </td>
   <td style="text-align:right;"> 0.1427922 </td>
   <td style="text-align:right;"> 0.4878268 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.5844535 </td>
   <td style="text-align:right;"> 0.8177333 </td>
   <td style="text-align:right;"> 0.1755108 </td>
   <td style="text-align:right;"> 0.8451777 </td>
   <td style="text-align:right;"> 0.7183016 </td>
   <td style="text-align:right;"> 0.8948936 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.9161951 </td>
   <td style="text-align:right;"> 0.3842754 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.6461032 </td>
   <td style="text-align:right;"> 0.5476553 </td>
   <td style="text-align:right;"> 0.3976805 </td>
   <td style="text-align:right;"> 0.3811707 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.6008058 </td>
   <td style="text-align:right;"> 0.4061416 </td>
   <td style="text-align:right;"> 0.1906611 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.3658783 </td>
   <td style="text-align:right;"> 0.6253230 </td>
   <td style="text-align:right;"> 0.4772483 </td>
   <td style="text-align:right;"> 0.8591992 </td>
   <td style="text-align:right;"> 0.2609659 </td>
   <td style="text-align:right;"> 0.5550425 </td>
   <td style="text-align:right;"> 0.2126416 </td>
   <td style="text-align:right;"> 0.7181992 </td>
   <td style="text-align:right;"> 0.1692498 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.2131754 </td>
   <td style="text-align:right;"> 0.7677309 </td>
   <td style="text-align:right;"> 0.7348579 </td>
   <td style="text-align:right;"> 0.7326545 </td>
   <td style="text-align:right;"> 0.5167072 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NIHSS_3mo </td>
   <td style="text-align:right;"> 0.0976348 </td>
   <td style="text-align:right;"> 0.9776273 </td>
   <td style="text-align:right;"> 0.9288234 </td>
   <td style="text-align:right;"> 0.1909980 </td>
   <td style="text-align:right;"> 0.5332640 </td>
   <td style="text-align:right;"> 0.8463882 </td>
   <td style="text-align:right;"> 0.3704343 </td>
   <td style="text-align:right;"> 0.0632815 </td>
   <td style="text-align:right;"> 0.6363185 </td>
   <td style="text-align:right;"> 0.1633587 </td>
   <td style="text-align:right;"> 0.6942951 </td>
   <td style="text-align:right;"> 0.0513727 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0474417 </td>
   <td style="text-align:right;"> 0.6033733 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.4933012 </td>
   <td style="text-align:right;"> 0.4144646 </td>
   <td style="text-align:right;"> 0.7869062 </td>
   <td style="text-align:right;"> 0.0858812 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.6599151 </td>
   <td style="text-align:right;"> 0.7028201 </td>
   <td style="text-align:right;"> 0.6470452 </td>
   <td style="text-align:right;"> 0.2336319 </td>
   <td style="text-align:right;"> 0.1683476 </td>
   <td style="text-align:right;"> 0.1260899 </td>
   <td style="text-align:right;"> 0.8483870 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.0101833 </td>
   <td style="text-align:right;"> 0.0063197 </td>
   <td style="text-align:right;"> 0.0999526 </td>
   <td style="text-align:right;"> 0.2630073 </td>
   <td style="text-align:right;"> 0.6689929 </td>
   <td style="text-align:right;"> 0.2057286 </td>
   <td style="text-align:right;"> 0.0820271 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.3820342 </td>
   <td style="text-align:right;"> 0.7923621 </td>
   <td style="text-align:right;"> 0.3315460 </td>
   <td style="text-align:right;"> 0.3842754 </td>
   <td style="text-align:right;"> 0.1769500 </td>
   <td style="text-align:right;"> 0.3307208 </td>
   <td style="text-align:right;"> 0.0679049 </td>
   <td style="text-align:right;"> 0.3430160 </td>
   <td style="text-align:right;"> 0.1433045 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.0065936 </td>
   <td style="text-align:right;"> 0.0677543 </td>
   <td style="text-align:right;"> 0.4519715 </td>
   <td style="text-align:right;"> 0.3368211 </td>
   <td style="text-align:right;"> 0.8211133 </td>
   <td style="text-align:right;"> 0.5388404 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mRS </td>
   <td style="text-align:right;"> 0.4616450 </td>
   <td style="text-align:right;"> 0.3404071 </td>
   <td style="text-align:right;"> 0.6488064 </td>
   <td style="text-align:right;"> 0.1020397 </td>
   <td style="text-align:right;"> 0.5653357 </td>
   <td style="text-align:right;"> 0.2387610 </td>
   <td style="text-align:right;"> 0.9278214 </td>
   <td style="text-align:right;"> 0.2224775 </td>
   <td style="text-align:right;"> 0.4432053 </td>
   <td style="text-align:right;"> 0.8598380 </td>
   <td style="text-align:right;"> 0.8649668 </td>
   <td style="text-align:right;"> 0.0745959 </td>
   <td style="text-align:right;"> 0.0014115 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4480027 </td>
   <td style="text-align:right;"> 0.1980116 </td>
   <td style="text-align:right;"> 0.9727347 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.6958174 </td>
   <td style="text-align:right;"> 0.4856913 </td>
   <td style="text-align:right;"> 0.4370291 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.7066730 </td>
   <td style="text-align:right;"> 0.9660226 </td>
   <td style="text-align:right;"> 0.7400494 </td>
   <td style="text-align:right;"> 0.3307208 </td>
   <td style="text-align:right;"> 0.2156370 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.2715956 </td>
   <td style="text-align:right;"> 0.0346104 </td>
   <td style="text-align:right;"> 0.1417060 </td>
   <td style="text-align:right;"> 0.6352294 </td>
   <td style="text-align:right;"> 0.1705204 </td>
   <td style="text-align:right;"> 0.7771138 </td>
   <td style="text-align:right;"> 0.6940736 </td>
   <td style="text-align:right;"> 0.2630073 </td>
   <td style="text-align:right;"> 0.9180665 </td>
   <td style="text-align:right;"> 0.4738537 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.4061416 </td>
   <td style="text-align:right;"> 0.4217927 </td>
   <td style="text-align:right;"> 0.4312296 </td>
   <td style="text-align:right;"> 0.7786512 </td>
   <td style="text-align:right;"> 0.1265585 </td>
   <td style="text-align:right;"> 0.5431755 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.9471271 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.2367100 </td>
   <td style="text-align:right;"> 0.9192539 </td>
   <td style="text-align:right;"> 0.9921956 </td>
   <td style="text-align:right;"> 0.7846697 </td>
   <td style="text-align:right;"> 0.8608818 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CEC </td>
   <td style="text-align:right;"> 0.0637547 </td>
   <td style="text-align:right;"> 0.4244123 </td>
   <td style="text-align:right;"> 0.7720437 </td>
   <td style="text-align:right;"> 0.3761921 </td>
   <td style="text-align:right;"> 0.3687585 </td>
   <td style="text-align:right;"> 0.4785434 </td>
   <td style="text-align:right;"> 0.1825489 </td>
   <td style="text-align:right;"> 0.4593834 </td>
   <td style="text-align:right;"> 0.2462904 </td>
   <td style="text-align:right;"> 0.2769111 </td>
   <td style="text-align:right;"> 0.5649209 </td>
   <td style="text-align:right;"> 0.0185485 </td>
   <td style="text-align:right;"> 0.3489757 </td>
   <td style="text-align:right;"> 0.2065429 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.8423521 </td>
   <td style="text-align:right;"> 0.9139271 </td>
   <td style="text-align:right;"> 0.4927776 </td>
   <td style="text-align:right;"> 0.5432158 </td>
   <td style="text-align:right;"> 0.9798903 </td>
   <td style="text-align:right;"> 0.4101880 </td>
   <td style="text-align:right;"> 0.6596618 </td>
   <td style="text-align:right;"> 0.8220028 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.9544699 </td>
   <td style="text-align:right;"> 0.2927537 </td>
   <td style="text-align:right;"> 0.5786988 </td>
   <td style="text-align:right;"> 0.5773964 </td>
   <td style="text-align:right;"> 0.9897076 </td>
   <td style="text-align:right;"> 0.9353484 </td>
   <td style="text-align:right;"> 0.9727465 </td>
   <td style="text-align:right;"> 0.7651509 </td>
   <td style="text-align:right;"> 0.7357528 </td>
   <td style="text-align:right;"> 0.6731019 </td>
   <td style="text-align:right;"> 0.1471424 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.2367100 </td>
   <td style="text-align:right;"> 0.7655537 </td>
   <td style="text-align:right;"> 0.8742524 </td>
   <td style="text-align:right;"> 0.6591875 </td>
   <td style="text-align:right;"> 0.7357528 </td>
   <td style="text-align:right;"> 0.3985795 </td>
   <td style="text-align:right;"> 0.5121105 </td>
   <td style="text-align:right;"> 0.4631920 </td>
   <td style="text-align:right;"> 0.1260899 </td>
   <td style="text-align:right;"> 0.6462071 </td>
   <td style="text-align:right;"> 0.3304656 </td>
   <td style="text-align:right;"> 0.3152488 </td>
   <td style="text-align:right;"> 0.9411736 </td>
   <td style="text-align:right;"> 0.8137934 </td>
   <td style="text-align:right;"> 0.7313016 </td>
   <td style="text-align:right;"> 0.7655537 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.0979688 </td>
   <td style="text-align:right;"> 0.7970959 </td>
   <td style="text-align:right;"> 0.7661661 </td>
   <td style="text-align:right;"> 0.6763713 </td>
   <td style="text-align:right;"> 0.9126089 </td>
   <td style="text-align:right;"> 0.9781932 </td>
   <td style="text-align:right;"> 0.6689338 </td>
   <td style="text-align:right;"> 0.6816170 </td>
   <td style="text-align:right;"> 0.1163923 </td>
   <td style="text-align:right;"> 0.9095487 </td>
   <td style="text-align:right;"> 0.7773774 </td>
   <td style="text-align:right;"> 0.2439134 </td>
   <td style="text-align:right;"> 0.0507778 </td>
   <td style="text-align:right;"> 0.0412404 </td>
   <td style="text-align:right;"> 0.6825393 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.1456649 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.8608818 </td>
   <td style="text-align:right;"> 0.2603510 </td>
   <td style="text-align:right;"> 0.2023322 </td>
   <td style="text-align:right;"> 0.3632971 </td>
   <td style="text-align:right;"> 0.1940901 </td>
   <td style="text-align:right;"> 0.3603448 </td>
   <td style="text-align:right;"> 0.5347501 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.3078751 </td>
   <td style="text-align:right;"> 0.0646294 </td>
   <td style="text-align:right;"> 0.5140630 </td>
   <td style="text-align:right;"> 0.0820271 </td>
   <td style="text-align:right;"> 0.4733762 </td>
   <td style="text-align:right;"> 0.1699505 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.6575859 </td>
   <td style="text-align:right;"> 0.1334059 </td>
   <td style="text-align:right;"> 0.1437541 </td>
   <td style="text-align:right;"> 0.3105896 </td>
   <td style="text-align:right;"> 0.6285719 </td>
   <td style="text-align:right;"> 0.8591992 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.7588192 </td>
   <td style="text-align:right;"> 0.1219558 </td>
   <td style="text-align:right;"> 0.1207609 </td>
   <td style="text-align:right;"> 0.7308369 </td>
   <td style="text-align:right;"> 0.2786305 </td>
   <td style="text-align:right;"> 0.0679049 </td>
   <td style="text-align:right;"> 0.2425228 </td>
   <td style="text-align:right;"> 0.1771804 </td>
   <td style="text-align:right;"> 0.4143197 </td>
   <td style="text-align:right;"> 0.7489393 </td>
   <td style="text-align:right;"> 0.0646294 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.0029982 </td>
   <td style="text-align:right;"> 0.1206801 </td>
   <td style="text-align:right;"> 0.4315850 </td>
   <td style="text-align:right;"> 0.2140878 </td>
   <td style="text-align:right;"> 0.5520843 </td>
   <td style="text-align:right;"> 0.5990196 </td>
   <td style="text-align:right;"> 0.7659889 </td>
   <td style="text-align:right;"> 0.8180645 </td>
   <td style="text-align:right;"> 0.0357107 </td>
   <td style="text-align:right;"> 0.4148639 </td>
   <td style="text-align:right;"> 0.4877489 </td>
   <td style="text-align:right;"> 0.3407500 </td>
   <td style="text-align:right;"> 0.2484405 </td>
   <td style="text-align:right;"> 0.9466163 </td>
   <td style="text-align:right;"> 0.6986999 </td>
   <td style="text-align:right;"> 0.0196616 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0065936 </td>
   <td style="text-align:right;"> 0.6596618 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.3820598 </td>
   <td style="text-align:right;"> 0.2131754 </td>
   <td style="text-align:right;"> 0.1788926 </td>
   <td style="text-align:right;"> 0.5359895 </td>
   <td style="text-align:right;"> 0.7257756 </td>
   <td style="text-align:right;"> 0.3105036 </td>
   <td style="text-align:right;"> 0.2609659 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.1692498 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.2449623 </td>
   <td style="text-align:right;"> 0.4370291 </td>
   <td style="text-align:right;"> 0.2165099 </td>
   <td style="text-align:right;"> 0.0679049 </td>
   <td style="text-align:right;"> 0.6573562 </td>
   <td style="text-align:right;"> 0.0002086 </td>
   <td style="text-align:right;"> 0.2105899 </td>
   <td style="text-align:right;"> 0.1219558 </td>
   <td style="text-align:right;"> 0.8425236 </td>
   <td style="text-align:right;"> 0.7146947 </td>
   <td style="text-align:right;"> 0.0684889 </td>
   <td style="text-align:right;"> 0.4367475 </td>
   <td style="text-align:right;"> 0.3461724 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.6668826 </td>
   <td style="text-align:right;"> 0.1762757 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.2715956 </td>
   <td style="text-align:right;"> 0.0213974 </td>
   <td style="text-align:right;"> 0.1219558 </td>
   <td style="text-align:right;"> 0.2603510 </td>
   <td style="text-align:right;"> 0.0059469 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.0093184 </td>
   <td style="text-align:right;"> 0.4610528 </td>
   <td style="text-align:right;"> 0.3738770 </td>
   <td style="text-align:right;"> 0.3449650 </td>
   <td style="text-align:right;"> 0.5654037 </td>
   <td style="text-align:right;"> 0.8410863 </td>
   <td style="text-align:right;"> 0.8307721 </td>
   <td style="text-align:right;"> 0.8057431 </td>
   <td style="text-align:right;"> 0.0494209 </td>
   <td style="text-align:right;"> 0.7510088 </td>
   <td style="text-align:right;"> 0.1744341 </td>
   <td style="text-align:right;"> 0.3273456 </td>
   <td style="text-align:right;"> 0.1783581 </td>
   <td style="text-align:right;"> 0.6911644 </td>
   <td style="text-align:right;"> 0.8194535 </td>
   <td style="text-align:right;"> 0.0049635 </td>
   <td style="text-align:right;"> 0.0000868 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6428056 </td>
   <td style="text-align:right;"> 0.3078751 </td>
   <td style="text-align:right;"> 0.2813276 </td>
   <td style="text-align:right;"> 0.3315460 </td>
   <td style="text-align:right;"> 0.2429087 </td>
   <td style="text-align:right;"> 0.3443578 </td>
   <td style="text-align:right;"> 0.4334720 </td>
   <td style="text-align:right;"> 0.0704253 </td>
   <td style="text-align:right;"> 0.3223092 </td>
   <td style="text-align:right;"> 0.3998941 </td>
   <td style="text-align:right;"> 0.1814699 </td>
   <td style="text-align:right;"> 0.9069815 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.1433045 </td>
   <td style="text-align:right;"> 0.2609659 </td>
   <td style="text-align:right;"> 0.2776812 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.9471271 </td>
   <td style="text-align:right;"> 0.0030454 </td>
   <td style="text-align:right;"> 0.4447703 </td>
   <td style="text-align:right;"> 0.1769500 </td>
   <td style="text-align:right;"> 0.2924356 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.0031048 </td>
   <td style="text-align:right;"> 0.7273135 </td>
   <td style="text-align:right;"> 0.2131754 </td>
   <td style="text-align:right;"> 0.3111872 </td>
   <td style="text-align:right;"> 0.9925673 </td>
   <td style="text-align:right;"> 0.2219724 </td>
   <td style="text-align:right;"> 0.2131754 </td>
   <td style="text-align:right;"> 0.1773547 </td>
   <td style="text-align:right;"> 0.0030454 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.4470100 </td>
   <td style="text-align:right;"> 0.0020105 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 0.5301500 </td>
   <td style="text-align:right;"> 0.7681355 </td>
   <td style="text-align:right;"> 0.9739043 </td>
   <td style="text-align:right;"> 0.8796582 </td>
   <td style="text-align:right;"> 0.2829293 </td>
   <td style="text-align:right;"> 0.4324722 </td>
   <td style="text-align:right;"> 0.2288890 </td>
   <td style="text-align:right;"> 0.9326450 </td>
   <td style="text-align:right;"> 0.6343886 </td>
   <td style="text-align:right;"> 0.2314371 </td>
   <td style="text-align:right;"> 0.0816274 </td>
   <td style="text-align:right;"> 0.6550982 </td>
   <td style="text-align:right;"> 0.6024572 </td>
   <td style="text-align:right;"> 0.4675812 </td>
   <td style="text-align:right;"> 0.2478192 </td>
   <td style="text-align:right;"> 0.7271890 </td>
   <td style="text-align:right;"> 0.4212644 </td>
   <td style="text-align:right;"> 0.3974386 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.9231215 </td>
   <td style="text-align:right;"> 0.7547639 </td>
   <td style="text-align:right;"> 0.2469587 </td>
   <td style="text-align:right;"> 0.3660693 </td>
   <td style="text-align:right;"> 0.7663131 </td>
   <td style="text-align:right;"> 0.6031351 </td>
   <td style="text-align:right;"> 0.6304238 </td>
   <td style="text-align:right;"> 0.9350584 </td>
   <td style="text-align:right;"> 0.2414288 </td>
   <td style="text-align:right;"> 0.4614243 </td>
   <td style="text-align:right;"> 0.3040382 </td>
   <td style="text-align:right;"> 0.9145204 </td>
   <td style="text-align:right;"> 0.8422695 </td>
   <td style="text-align:right;"> 0.9126519 </td>
   <td style="text-align:right;"> 0.5238411 </td>
   <td style="text-align:right;"> 0.4965329 </td>
   <td style="text-align:right;"> 0.9122387 </td>
   <td style="text-align:right;"> 0.3461724 </td>
   <td style="text-align:right;"> 0.6428056 </td>
   <td style="text-align:right;"> 0.4738537 </td>
   <td style="text-align:right;"> 0.5233423 </td>
   <td style="text-align:right;"> 0.9465305 </td>
   <td style="text-align:right;"> 0.5620014 </td>
   <td style="text-align:right;"> 0.7247234 </td>
   <td style="text-align:right;"> 0.9120676 </td>
   <td style="text-align:right;"> 0.4796402 </td>
   <td style="text-align:right;"> 0.9353484 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.8164696 </td>
   <td style="text-align:right;"> 0.4367475 </td>
   <td style="text-align:right;"> 0.2603510 </td>
   <td style="text-align:right;"> 0.9254848 </td>
   <td style="text-align:right;"> 0.0679049 </td>
   <td style="text-align:right;"> 0.3011440 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 0.0405488 </td>
   <td style="text-align:right;"> 0.9307688 </td>
   <td style="text-align:right;"> 0.4315181 </td>
   <td style="text-align:right;"> 0.4817009 </td>
   <td style="text-align:right;"> 0.1083032 </td>
   <td style="text-align:right;"> 0.3078967 </td>
   <td style="text-align:right;"> 0.1679585 </td>
   <td style="text-align:right;"> 0.8256632 </td>
   <td style="text-align:right;"> 0.2397312 </td>
   <td style="text-align:right;"> 0.1012100 </td>
   <td style="text-align:right;"> 0.7229983 </td>
   <td style="text-align:right;"> 0.0298542 </td>
   <td style="text-align:right;"> 0.0051728 </td>
   <td style="text-align:right;"> 0.2417883 </td>
   <td style="text-align:right;"> 0.2913182 </td>
   <td style="text-align:right;"> 0.0705307 </td>
   <td style="text-align:right;"> 0.0253669 </td>
   <td style="text-align:right;"> 0.1009104 </td>
   <td style="text-align:right;"> 0.8480927 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.3279466 </td>
   <td style="text-align:right;"> 0.6794488 </td>
   <td style="text-align:right;"> 0.6043755 </td>
   <td style="text-align:right;"> 0.7869062 </td>
   <td style="text-align:right;"> 0.3700763 </td>
   <td style="text-align:right;"> 0.3211473 </td>
   <td style="text-align:right;"> 0.0646294 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.0997341 </td>
   <td style="text-align:right;"> 0.3357679 </td>
   <td style="text-align:right;"> 0.1036481 </td>
   <td style="text-align:right;"> 0.1584373 </td>
   <td style="text-align:right;"> 0.7840744 </td>
   <td style="text-align:right;"> 0.2609659 </td>
   <td style="text-align:right;"> 0.0256314 </td>
   <td style="text-align:right;"> 0.0432958 </td>
   <td style="text-align:right;"> 0.1955318 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.3189744 </td>
   <td style="text-align:right;"> 0.9288310 </td>
   <td style="text-align:right;"> 0.9897566 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.2307214 </td>
   <td style="text-align:right;"> 0.2650631 </td>
   <td style="text-align:right;"> 0.2735961 </td>
   <td style="text-align:right;"> 0.0333065 </td>
   <td style="text-align:right;"> 0.5253909 </td>
   <td style="text-align:right;"> 0.0608425 </td>
   <td style="text-align:right;"> 0.2681869 </td>
   <td style="text-align:right;"> 0.3009524 </td>
   <td style="text-align:right;"> 0.3699431 </td>
   <td style="text-align:right;"> 0.6019492 </td>
   <td style="text-align:right;"> 0.3195884 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.3514214 </td>
   <td style="text-align:right;"> 0.5514848 </td>
   <td style="text-align:right;"> 0.1748401 </td>
   <td style="text-align:right;"> 0.5976020 </td>
   <td style="text-align:right;"> 0.1311017 </td>
   <td style="text-align:right;"> 0.2518129 </td>
   <td style="text-align:right;"> 0.6814788 </td>
   <td style="text-align:right;"> 0.3613847 </td>
   <td style="text-align:right;"> 0.0576560 </td>
   <td style="text-align:right;"> 0.7770668 </td>
   <td style="text-align:right;"> 0.1908595 </td>
   <td style="text-align:right;"> 0.7034970 </td>
   <td style="text-align:right;"> 0.0140229 </td>
   <td style="text-align:right;"> 0.1976939 </td>
   <td style="text-align:right;"> 0.9617636 </td>
   <td style="text-align:right;"> 0.0424340 </td>
   <td style="text-align:right;"> 0.1513822 </td>
   <td style="text-align:right;"> 0.0861540 </td>
   <td style="text-align:right;"> 0.5508743 </td>
   <td style="text-align:right;"> 0.1122091 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.5107596 </td>
   <td style="text-align:right;"> 0.5021129 </td>
   <td style="text-align:right;"> 0.5476553 </td>
   <td style="text-align:right;"> 0.1363127 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.1814966 </td>
   <td style="text-align:right;"> 0.3315460 </td>
   <td style="text-align:right;"> 0.1814966 </td>
   <td style="text-align:right;"> 0.3011440 </td>
   <td style="text-align:right;"> 0.1755108 </td>
   <td style="text-align:right;"> 0.3554490 </td>
   <td style="text-align:right;"> 0.3040382 </td>
   <td style="text-align:right;"> 0.5275394 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.4126249 </td>
   <td style="text-align:right;"> 0.4088785 </td>
   <td style="text-align:right;"> 0.8192778 </td>
   <td style="text-align:right;"> 0.1471424 </td>
   <td style="text-align:right;"> 0.5021129 </td>
   <td style="text-align:right;"> 0.7313016 </td>
   <td style="text-align:right;"> 0.2824454 </td>
   <td style="text-align:right;"> 0.9147856 </td>
   <td style="text-align:right;"> 0.3974393 </td>
   <td style="text-align:right;"> 0.5140630 </td>
   <td style="text-align:right;"> 0.5448245 </td>
   <td style="text-align:right;"> 0.3078751 </td>
   <td style="text-align:right;"> 0.1265585 </td>
   <td style="text-align:right;"> 0.1004852 </td>
   <td style="text-align:right;"> 0.2725257 </td>
   <td style="text-align:right;"> 0.2824454 </td>
   <td style="text-align:right;"> 0.6285719 </td>
   <td style="text-align:right;"> 0.3105036 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.1865985 </td>
   <td style="text-align:right;"> 0.1046657 </td>
   <td style="text-align:right;"> 0.3903164 </td>
   <td style="text-align:right;"> 0.2312491 </td>
   <td style="text-align:right;"> 0.3737866 </td>
   <td style="text-align:right;"> 0.4795460 </td>
   <td style="text-align:right;"> 0.4856482 </td>
   <td style="text-align:right;"> 0.5340864 </td>
   <td style="text-align:right;"> 0.4902608 </td>
   <td style="text-align:right;"> 0.3326675 </td>
   <td style="text-align:right;"> 0.6934318 </td>
   <td style="text-align:right;"> 0.5030196 </td>
   <td style="text-align:right;"> 0.4223840 </td>
   <td style="text-align:right;"> 0.6231846 </td>
   <td style="text-align:right;"> 0.1742243 </td>
   <td style="text-align:right;"> 0.1392024 </td>
   <td style="text-align:right;"> 0.0466621 </td>
   <td style="text-align:right;"> 0.1180220 </td>
   <td style="text-align:right;"> 0.0645175 </td>
   <td style="text-align:right;"> 0.4486926 </td>
   <td style="text-align:right;"> 0.2620516 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0184227 </td>
   <td style="text-align:right;"> 0.6649028 </td>
   <td style="text-align:right;"> 0.9151110 </td>
   <td style="text-align:right;"> 0.3223092 </td>
   <td style="text-align:right;"> 0.2699089 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.6008058 </td>
   <td style="text-align:right;"> 0.5719593 </td>
   <td style="text-align:right;"> 0.1814966 </td>
   <td style="text-align:right;"> 0.1814966 </td>
   <td style="text-align:right;"> 0.5160666 </td>
   <td style="text-align:right;"> 0.9112983 </td>
   <td style="text-align:right;"> 0.1814699 </td>
   <td style="text-align:right;"> 0.4353874 </td>
   <td style="text-align:right;"> 0.5714886 </td>
   <td style="text-align:right;"> 0.6539417 </td>
   <td style="text-align:right;"> 0.9112983 </td>
   <td style="text-align:right;"> 0.2102899 </td>
   <td style="text-align:right;"> 0.6253230 </td>
   <td style="text-align:right;"> 0.6056486 </td>
   <td style="text-align:right;"> 0.2160365 </td>
   <td style="text-align:right;"> 0.9121943 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.2057286 </td>
   <td style="text-align:right;"> 0.1456649 </td>
   <td style="text-align:right;"> 0.1814966 </td>
   <td style="text-align:right;"> 0.3315460 </td>
   <td style="text-align:right;"> 0.1181558 </td>
   <td style="text-align:right;"> 0.0997341 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 0.4200048 </td>
   <td style="text-align:right;"> 0.4285266 </td>
   <td style="text-align:right;"> 0.8425382 </td>
   <td style="text-align:right;"> 0.4040640 </td>
   <td style="text-align:right;"> 0.4956884 </td>
   <td style="text-align:right;"> 0.4786917 </td>
   <td style="text-align:right;"> 0.9630099 </td>
   <td style="text-align:right;"> 0.7315857 </td>
   <td style="text-align:right;"> 0.7832215 </td>
   <td style="text-align:right;"> 0.9129110 </td>
   <td style="text-align:right;"> 0.5021932 </td>
   <td style="text-align:right;"> 0.7734531 </td>
   <td style="text-align:right;"> 0.4742799 </td>
   <td style="text-align:right;"> 0.4898110 </td>
   <td style="text-align:right;"> 0.4210406 </td>
   <td style="text-align:right;"> 0.0395786 </td>
   <td style="text-align:right;"> 0.0314166 </td>
   <td style="text-align:right;"> 0.0624018 </td>
   <td style="text-align:right;"> 0.1410615 </td>
   <td style="text-align:right;"> 0.3526246 </td>
   <td style="text-align:right;"> 0.2550754 </td>
   <td style="text-align:right;"> 0.0003342 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6607688 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.2752137 </td>
   <td style="text-align:right;"> 0.9777522 </td>
   <td style="text-align:right;"> 0.1943267 </td>
   <td style="text-align:right;"> 0.8769057 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.8555500 </td>
   <td style="text-align:right;"> 0.1460815 </td>
   <td style="text-align:right;"> 0.2371021 </td>
   <td style="text-align:right;"> 0.4061416 </td>
   <td style="text-align:right;"> 0.9238855 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.3367028 </td>
   <td style="text-align:right;"> 0.5562687 </td>
   <td style="text-align:right;"> 0.8480176 </td>
   <td style="text-align:right;"> 0.8046712 </td>
   <td style="text-align:right;"> 0.1471424 </td>
   <td style="text-align:right;"> 0.9684830 </td>
   <td style="text-align:right;"> 0.8536794 </td>
   <td style="text-align:right;"> 0.3985795 </td>
   <td style="text-align:right;"> 0.9465305 </td>
   <td style="text-align:right;"> 0.4714586 </td>
   <td style="text-align:right;"> 0.1680427 </td>
   <td style="text-align:right;"> 0.3310143 </td>
   <td style="text-align:right;"> 0.0974954 </td>
   <td style="text-align:right;"> 0.7769474 </td>
   <td style="text-align:right;"> 0.1700233 </td>
   <td style="text-align:right;"> 0.1062989 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 0.1660394 </td>
   <td style="text-align:right;"> 0.8959466 </td>
   <td style="text-align:right;"> 0.2433362 </td>
   <td style="text-align:right;"> 0.3116187 </td>
   <td style="text-align:right;"> 0.2152724 </td>
   <td style="text-align:right;"> 0.4063237 </td>
   <td style="text-align:right;"> 0.4532876 </td>
   <td style="text-align:right;"> 0.1927703 </td>
   <td style="text-align:right;"> 0.8641687 </td>
   <td style="text-align:right;"> 0.3407496 </td>
   <td style="text-align:right;"> 0.8357455 </td>
   <td style="text-align:right;"> 0.3414801 </td>
   <td style="text-align:right;"> 0.4032959 </td>
   <td style="text-align:right;"> 0.9370619 </td>
   <td style="text-align:right;"> 0.6611455 </td>
   <td style="text-align:right;"> 0.1370299 </td>
   <td style="text-align:right;"> 0.2854980 </td>
   <td style="text-align:right;"> 0.1274474 </td>
   <td style="text-align:right;"> 0.5713807 </td>
   <td style="text-align:right;"> 0.6021506 </td>
   <td style="text-align:right;"> 0.2956862 </td>
   <td style="text-align:right;"> 0.4299190 </td>
   <td style="text-align:right;"> 0.4257642 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6187062 </td>
   <td style="text-align:right;"> 0.7547639 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.4663638 </td>
   <td style="text-align:right;"> 0.9167815 </td>
   <td style="text-align:right;"> 0.4872298 </td>
   <td style="text-align:right;"> 0.9456724 </td>
   <td style="text-align:right;"> 0.9506640 </td>
   <td style="text-align:right;"> 0.5979156 </td>
   <td style="text-align:right;"> 0.4217927 </td>
   <td style="text-align:right;"> 0.9139271 </td>
   <td style="text-align:right;"> 0.4156137 </td>
   <td style="text-align:right;"> 0.4245170 </td>
   <td style="text-align:right;"> 0.6285719 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.9354082 </td>
   <td style="text-align:right;"> 0.4833604 </td>
   <td style="text-align:right;"> 0.9636034 </td>
   <td style="text-align:right;"> 0.5789362 </td>
   <td style="text-align:right;"> 0.6253230 </td>
   <td style="text-align:right;"> 0.4683706 </td>
   <td style="text-align:right;"> 0.6668826 </td>
   <td style="text-align:right;"> 0.9756764 </td>
   <td style="text-align:right;"> 0.9506640 </td>
   <td style="text-align:right;"> 0.6573562 </td>
   <td style="text-align:right;"> 0.9112983 </td>
   <td style="text-align:right;"> 0.8866011 </td>
   <td style="text-align:right;"> 0.2535174 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.5400016 </td>
   <td style="text-align:right;"> 0.9815509 </td>
   <td style="text-align:right;"> 0.9160834 </td>
   <td style="text-align:right;"> 0.5287623 </td>
   <td style="text-align:right;"> 0.2736011 </td>
   <td style="text-align:right;"> 0.5755263 </td>
   <td style="text-align:right;"> 0.1990060 </td>
   <td style="text-align:right;"> 0.8281988 </td>
   <td style="text-align:right;"> 0.0278847 </td>
   <td style="text-align:right;"> 0.0714244 </td>
   <td style="text-align:right;"> 0.6289017 </td>
   <td style="text-align:right;"> 0.8346811 </td>
   <td style="text-align:right;"> 0.0584928 </td>
   <td style="text-align:right;"> 0.5343608 </td>
   <td style="text-align:right;"> 0.6237143 </td>
   <td style="text-align:right;"> 0.2844230 </td>
   <td style="text-align:right;"> 0.5108871 </td>
   <td style="text-align:right;"> 0.1944018 </td>
   <td style="text-align:right;"> 0.3479626 </td>
   <td style="text-align:right;"> 0.1434113 </td>
   <td style="text-align:right;"> 0.0168165 </td>
   <td style="text-align:right;"> 0.8304268 </td>
   <td style="text-align:right;"> 0.7767746 </td>
   <td style="text-align:right;"> 0.3677216 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.2535174 </td>
   <td style="text-align:right;"> 0.5905921 </td>
   <td style="text-align:right;"> 0.3326865 </td>
   <td style="text-align:right;"> 0.1683476 </td>
   <td style="text-align:right;"> 0.3219903 </td>
   <td style="text-align:right;"> 0.1334059 </td>
   <td style="text-align:right;"> 0.0704253 </td>
   <td style="text-align:right;"> 0.9358488 </td>
   <td style="text-align:right;"> 0.4827031 </td>
   <td style="text-align:right;"> 0.1700233 </td>
   <td style="text-align:right;"> 0.6364345 </td>
   <td style="text-align:right;"> 0.6143312 </td>
   <td style="text-align:right;"> 0.7785247 </td>
   <td style="text-align:right;"> 0.5347501 </td>
   <td style="text-align:right;"> 0.2371021 </td>
   <td style="text-align:right;"> 0.1503813 </td>
   <td style="text-align:right;"> 0.4158008 </td>
   <td style="text-align:right;"> 0.3603448 </td>
   <td style="text-align:right;"> 0.2210602 </td>
   <td style="text-align:right;"> 0.4335766 </td>
   <td style="text-align:right;"> 0.6958174 </td>
   <td style="text-align:right;"> 0.1877684 </td>
   <td style="text-align:right;"> 0.4208246 </td>
   <td style="text-align:right;"> 0.2247968 </td>
   <td style="text-align:right;"> 0.4183279 </td>
   <td style="text-align:right;"> 0.1914335 </td>
   <td style="text-align:right;"> 0.4891131 </td>
   <td style="text-align:right;"> 0.5110666 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 0.1649634 </td>
   <td style="text-align:right;"> 0.6271903 </td>
   <td style="text-align:right;"> 0.2746369 </td>
   <td style="text-align:right;"> 0.7428544 </td>
   <td style="text-align:right;"> 0.6378670 </td>
   <td style="text-align:right;"> 0.7155182 </td>
   <td style="text-align:right;"> 0.8623838 </td>
   <td style="text-align:right;"> 0.4504679 </td>
   <td style="text-align:right;"> 0.5487172 </td>
   <td style="text-align:right;"> 0.6244945 </td>
   <td style="text-align:right;"> 0.1613970 </td>
   <td style="text-align:right;"> 0.1528178 </td>
   <td style="text-align:right;"> 0.0272435 </td>
   <td style="text-align:right;"> 0.1157467 </td>
   <td style="text-align:right;"> 0.9120882 </td>
   <td style="text-align:right;"> 0.0253689 </td>
   <td style="text-align:right;"> 0.1028376 </td>
   <td style="text-align:right;"> 0.0033804 </td>
   <td style="text-align:right;"> 0.3833782 </td>
   <td style="text-align:right;"> 0.1086028 </td>
   <td style="text-align:right;"> 0.0528945 </td>
   <td style="text-align:right;"> 0.1096974 </td>
   <td style="text-align:right;"> 0.0528066 </td>
   <td style="text-align:right;"> 0.5521060 </td>
   <td style="text-align:right;"> 0.0673348 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.2135610 </td>
   <td style="text-align:right;"> 0.6125213 </td>
   <td style="text-align:right;"> 0.1427922 </td>
   <td style="text-align:right;"> 0.3417962 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.1080730 </td>
   <td style="text-align:right;"> 0.1460815 </td>
   <td style="text-align:right;"> 0.4061551 </td>
   <td style="text-align:right;"> 0.1802483 </td>
   <td style="text-align:right;"> 0.8109993 </td>
   <td style="text-align:right;"> 0.1814966 </td>
   <td style="text-align:right;"> 0.6428159 </td>
   <td style="text-align:right;"> 0.5620014 </td>
   <td style="text-align:right;"> 0.0608425 </td>
   <td style="text-align:right;"> 0.2650631 </td>
   <td style="text-align:right;"> 0.0020737 </td>
   <td style="text-align:right;"> 0.9154700 </td>
   <td style="text-align:right;"> 0.2449623 </td>
   <td style="text-align:right;"> 0.6668826 </td>
   <td style="text-align:right;"> 0.8950023 </td>
   <td style="text-align:right;"> 0.5342046 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.0679049 </td>
   <td style="text-align:right;"> 0.0514570 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.5373617 </td>
   <td style="text-align:right;"> 0.1073639 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.0581718 </td>
   <td style="text-align:right;"> 0.9563268 </td>
   <td style="text-align:right;"> 0.7773601 </td>
   <td style="text-align:right;"> 0.7714980 </td>
   <td style="text-align:right;"> 0.2337279 </td>
   <td style="text-align:right;"> 0.3066221 </td>
   <td style="text-align:right;"> 0.8847351 </td>
   <td style="text-align:right;"> 0.8095890 </td>
   <td style="text-align:right;"> 0.3205833 </td>
   <td style="text-align:right;"> 0.7500962 </td>
   <td style="text-align:right;"> 0.9115853 </td>
   <td style="text-align:right;"> 0.0059798 </td>
   <td style="text-align:right;"> 0.0144937 </td>
   <td style="text-align:right;"> 0.0485105 </td>
   <td style="text-align:right;"> 0.0913528 </td>
   <td style="text-align:right;"> 0.0043688 </td>
   <td style="text-align:right;"> 0.0719645 </td>
   <td style="text-align:right;"> 0.1095980 </td>
   <td style="text-align:right;"> 0.8658451 </td>
   <td style="text-align:right;"> 0.0024551 </td>
   <td style="text-align:right;"> 0.0337178 </td>
   <td style="text-align:right;"> 0.0777604 </td>
   <td style="text-align:right;"> 0.0822690 </td>
   <td style="text-align:right;"> 0.6804960 </td>
   <td style="text-align:right;"> 0.3317259 </td>
   <td style="text-align:right;"> 0.0475785 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4480027 </td>
   <td style="text-align:right;"> 0.1004852 </td>
   <td style="text-align:right;"> 0.3603448 </td>
   <td style="text-align:right;"> 0.1858813 </td>
   <td style="text-align:right;"> 0.4480027 </td>
   <td style="text-align:right;"> 0.2126416 </td>
   <td style="text-align:right;"> 0.0911555 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.1692739 </td>
   <td style="text-align:right;"> 0.2609412 </td>
   <td style="text-align:right;"> 0.1524893 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.7663131 </td>
   <td style="text-align:right;"> 0.9506640 </td>
   <td style="text-align:right;"> 0.1683476 </td>
   <td style="text-align:right;"> 0.6668826 </td>
   <td style="text-align:right;"> 0.2487154 </td>
   <td style="text-align:right;"> 0.2803367 </td>
   <td style="text-align:right;"> 0.1524320 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.0256314 </td>
   <td style="text-align:right;"> 0.1265585 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.4927776 </td>
   <td style="text-align:right;"> 0.8395371 </td>
   <td style="text-align:right;"> 0.1802483 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.1161997 </td>
   <td style="text-align:right;"> 0.3501741 </td>
   <td style="text-align:right;"> 0.4110497 </td>
   <td style="text-align:right;"> 0.6951055 </td>
   <td style="text-align:right;"> 0.6835708 </td>
   <td style="text-align:right;"> 0.5052819 </td>
   <td style="text-align:right;"> 0.7662600 </td>
   <td style="text-align:right;"> 0.4345469 </td>
   <td style="text-align:right;"> 0.1776334 </td>
   <td style="text-align:right;"> 0.8956390 </td>
   <td style="text-align:right;"> 0.0816177 </td>
   <td style="text-align:right;"> 0.4008841 </td>
   <td style="text-align:right;"> 0.7080153 </td>
   <td style="text-align:right;"> 0.6919169 </td>
   <td style="text-align:right;"> 0.3229459 </td>
   <td style="text-align:right;"> 0.1011542 </td>
   <td style="text-align:right;"> 0.3886925 </td>
   <td style="text-align:right;"> 0.1654134 </td>
   <td style="text-align:right;"> 0.0614960 </td>
   <td style="text-align:right;"> 0.6765504 </td>
   <td style="text-align:right;"> 0.1171872 </td>
   <td style="text-align:right;"> 0.7791662 </td>
   <td style="text-align:right;"> 0.9571754 </td>
   <td style="text-align:right;"> 0.2196445 </td>
   <td style="text-align:right;"> 0.1192650 </td>
   <td style="text-align:right;"> 0.3600452 </td>
   <td style="text-align:right;"> 0.2061473 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.5393244 </td>
   <td style="text-align:right;"> 0.9070175 </td>
   <td style="text-align:right;"> 0.6940736 </td>
   <td style="text-align:right;"> 0.7881789 </td>
   <td style="text-align:right;"> 0.9139271 </td>
   <td style="text-align:right;"> 0.4326261 </td>
   <td style="text-align:right;"> 0.1920253 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.4508706 </td>
   <td style="text-align:right;"> 0.6314647 </td>
   <td style="text-align:right;"> 0.1107600 </td>
   <td style="text-align:right;"> 0.4217927 </td>
   <td style="text-align:right;"> 0.9849491 </td>
   <td style="text-align:right;"> 0.5090339 </td>
   <td style="text-align:right;"> 0.8676437 </td>
   <td style="text-align:right;"> 0.3342780 </td>
   <td style="text-align:right;"> 0.4796402 </td>
   <td style="text-align:right;"> 0.8376058 </td>
   <td style="text-align:right;"> 0.1576294 </td>
   <td style="text-align:right;"> 0.6875233 </td>
   <td style="text-align:right;"> 0.3603448 </td>
   <td style="text-align:right;"> 0.5270490 </td>
   <td style="text-align:right;"> 0.7848482 </td>
   <td style="text-align:right;"> 0.9069815 </td>
   <td style="text-align:right;"> 0.2743478 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.0987612 </td>
   <td style="text-align:right;"> 0.5899027 </td>
   <td style="text-align:right;"> 0.6914667 </td>
   <td style="text-align:right;"> 0.6140334 </td>
   <td style="text-align:right;"> 0.8437193 </td>
   <td style="text-align:right;"> 0.6006790 </td>
   <td style="text-align:right;"> 0.0959299 </td>
   <td style="text-align:right;"> 0.7070916 </td>
   <td style="text-align:right;"> 0.1145571 </td>
   <td style="text-align:right;"> 0.0372023 </td>
   <td style="text-align:right;"> 0.9603161 </td>
   <td style="text-align:right;"> 0.2954302 </td>
   <td style="text-align:right;"> 0.0069940 </td>
   <td style="text-align:right;"> 0.0792950 </td>
   <td style="text-align:right;"> 0.3218000 </td>
   <td style="text-align:right;"> 0.0024857 </td>
   <td style="text-align:right;"> 0.0276351 </td>
   <td style="text-align:right;"> 0.0327824 </td>
   <td style="text-align:right;"> 0.2156439 </td>
   <td style="text-align:right;"> 0.0080778 </td>
   <td style="text-align:right;"> 0.0337024 </td>
   <td style="text-align:right;"> 0.0251788 </td>
   <td style="text-align:right;"> 0.0399089 </td>
   <td style="text-align:right;"> 0.8369456 </td>
   <td style="text-align:right;"> 0.0271228 </td>
   <td style="text-align:right;"> 0.0184426 </td>
   <td style="text-align:right;"> 0.0086510 </td>
   <td style="text-align:right;"> 0.2884486 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4208246 </td>
   <td style="text-align:right;"> 0.0256314 </td>
   <td style="text-align:right;"> 0.1214987 </td>
   <td style="text-align:right;"> 0.4701513 </td>
   <td style="text-align:right;"> 0.0683329 </td>
   <td style="text-align:right;"> 0.0511929 </td>
   <td style="text-align:right;"> 0.2803367 </td>
   <td style="text-align:right;"> 0.1391036 </td>
   <td style="text-align:right;"> 0.0646294 </td>
   <td style="text-align:right;"> 0.5021129 </td>
   <td style="text-align:right;"> 0.4833604 </td>
   <td style="text-align:right;"> 0.5479784 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.1576294 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.0142002 </td>
   <td style="text-align:right;"> 0.3304656 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.0068929 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.2165099 </td>
   <td style="text-align:right;"> 0.1457117 </td>
   <td style="text-align:right;"> 0.0931793 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.7667087 </td>
   <td style="text-align:right;"> 0.6693104 </td>
   <td style="text-align:right;"> 0.9270465 </td>
   <td style="text-align:right;"> 0.1249096 </td>
   <td style="text-align:right;"> 0.4351387 </td>
   <td style="text-align:right;"> 0.9281977 </td>
   <td style="text-align:right;"> 0.8113985 </td>
   <td style="text-align:right;"> 0.1640442 </td>
   <td style="text-align:right;"> 0.9669202 </td>
   <td style="text-align:right;"> 0.5956123 </td>
   <td style="text-align:right;"> 0.7634454 </td>
   <td style="text-align:right;"> 0.1587259 </td>
   <td style="text-align:right;"> 0.0001605 </td>
   <td style="text-align:right;"> 0.0009042 </td>
   <td style="text-align:right;"> 0.9825254 </td>
   <td style="text-align:right;"> 0.2653975 </td>
   <td style="text-align:right;"> 0.7793682 </td>
   <td style="text-align:right;"> 0.7957479 </td>
   <td style="text-align:right;"> 0.0977423 </td>
   <td style="text-align:right;"> 0.1220225 </td>
   <td style="text-align:right;"> 0.0953739 </td>
   <td style="text-align:right;"> 0.7776878 </td>
   <td style="text-align:right;"> 0.7464517 </td>
   <td style="text-align:right;"> 0.2429078 </td>
   <td style="text-align:right;"> 0.1091215 </td>
   <td style="text-align:right;"> 0.1260033 </td>
   <td style="text-align:right;"> 0.1369703 </td>
   <td style="text-align:right;"> 0.7969910 </td>
   <td style="text-align:right;"> 0.1844543 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.2752137 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.9826077 </td>
   <td style="text-align:right;"> 0.7186166 </td>
   <td style="text-align:right;"> 0.2650631 </td>
   <td style="text-align:right;"> 0.7663131 </td>
   <td style="text-align:right;"> 0.9069815 </td>
   <td style="text-align:right;"> 0.7840744 </td>
   <td style="text-align:right;"> 0.4161712 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.6043755 </td>
   <td style="text-align:right;"> 0.7448671 </td>
   <td style="text-align:right;"> 0.2144243 </td>
   <td style="text-align:right;"> 0.8676437 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.8082713 </td>
   <td style="text-align:right;"> 0.1943267 </td>
   <td style="text-align:right;"> 0.3279466 </td>
   <td style="text-align:right;"> 0.9727465 </td>
   <td style="text-align:right;"> 0.6581255 </td>
   <td style="text-align:right;"> 0.6470452 </td>
   <td style="text-align:right;"> 0.9139271 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 0.0530345 </td>
   <td style="text-align:right;"> 0.6517537 </td>
   <td style="text-align:right;"> 0.5945992 </td>
   <td style="text-align:right;"> 0.4497950 </td>
   <td style="text-align:right;"> 0.7211858 </td>
   <td style="text-align:right;"> 0.7287429 </td>
   <td style="text-align:right;"> 0.3593516 </td>
   <td style="text-align:right;"> 0.3152968 </td>
   <td style="text-align:right;"> 0.1346598 </td>
   <td style="text-align:right;"> 0.1371965 </td>
   <td style="text-align:right;"> 0.8706157 </td>
   <td style="text-align:right;"> 0.1496155 </td>
   <td style="text-align:right;"> 0.0000734 </td>
   <td style="text-align:right;"> 0.0180781 </td>
   <td style="text-align:right;"> 0.8681500 </td>
   <td style="text-align:right;"> 0.0041668 </td>
   <td style="text-align:right;"> 0.0516028 </td>
   <td style="text-align:right;"> 0.0138844 </td>
   <td style="text-align:right;"> 0.8218312 </td>
   <td style="text-align:right;"> 0.0090260 </td>
   <td style="text-align:right;"> 0.0299311 </td>
   <td style="text-align:right;"> 0.3454996 </td>
   <td style="text-align:right;"> 0.3907602 </td>
   <td style="text-align:right;"> 0.8895284 </td>
   <td style="text-align:right;"> 0.0161083 </td>
   <td style="text-align:right;"> 0.0100282 </td>
   <td style="text-align:right;"> 0.0348022 </td>
   <td style="text-align:right;"> 0.4644586 </td>
   <td style="text-align:right;"> 0.0005578 </td>
   <td style="text-align:right;"> 0.0067761 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0256314 </td>
   <td style="text-align:right;"> 0.3304656 </td>
   <td style="text-align:right;"> 0.3300747 </td>
   <td style="text-align:right;"> 0.1156032 </td>
   <td style="text-align:right;"> 0.3026619 </td>
   <td style="text-align:right;"> 0.2803367 </td>
   <td style="text-align:right;"> 0.2559459 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.2715956 </td>
   <td style="text-align:right;"> 0.2650631 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.0065936 </td>
   <td style="text-align:right;"> 0.1357212 </td>
   <td style="text-align:right;"> 0.4217927 </td>
   <td style="text-align:right;"> 0.3603448 </td>
   <td style="text-align:right;"> 0.0316554 </td>
   <td style="text-align:right;"> 0.0997341 </td>
   <td style="text-align:right;"> 0.2135456 </td>
   <td style="text-align:right;"> 0.1219558 </td>
   <td style="text-align:right;"> 0.6565642 </td>
   <td style="text-align:right;"> 0.2057286 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.0339683 </td>
   <td style="text-align:right;"> 0.5557714 </td>
   <td style="text-align:right;"> 0.0695556 </td>
   <td style="text-align:right;"> 0.5206781 </td>
   <td style="text-align:right;"> 0.5387020 </td>
   <td style="text-align:right;"> 0.8747172 </td>
   <td style="text-align:right;"> 0.0746189 </td>
   <td style="text-align:right;"> 0.2260388 </td>
   <td style="text-align:right;"> 0.1288080 </td>
   <td style="text-align:right;"> 0.0134680 </td>
   <td style="text-align:right;"> 0.6538947 </td>
   <td style="text-align:right;"> 0.1623025 </td>
   <td style="text-align:right;"> 0.0083755 </td>
   <td style="text-align:right;"> 0.3872225 </td>
   <td style="text-align:right;"> 0.9494514 </td>
   <td style="text-align:right;"> 0.2259998 </td>
   <td style="text-align:right;"> 0.0638182 </td>
   <td style="text-align:right;"> 0.0188149 </td>
   <td style="text-align:right;"> 0.6980201 </td>
   <td style="text-align:right;"> 0.0235701 </td>
   <td style="text-align:right;"> 0.1333337 </td>
   <td style="text-align:right;"> 0.3171095 </td>
   <td style="text-align:right;"> 0.7183392 </td>
   <td style="text-align:right;"> 0.9037517 </td>
   <td style="text-align:right;"> 0.0034242 </td>
   <td style="text-align:right;"> 0.0103524 </td>
   <td style="text-align:right;"> 0.2067705 </td>
   <td style="text-align:right;"> 0.6045755 </td>
   <td style="text-align:right;"> 0.0122557 </td>
   <td style="text-align:right;"> 0.0823888 </td>
   <td style="text-align:right;"> 0.0005293 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6920357 </td>
   <td style="text-align:right;"> 0.5342046 </td>
   <td style="text-align:right;"> 0.1802483 </td>
   <td style="text-align:right;"> 0.4797683 </td>
   <td style="text-align:right;"> 0.2592034 </td>
   <td style="text-align:right;"> 0.5433303 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.1391036 </td>
   <td style="text-align:right;"> 0.1683476 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.1337855 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.3078751 </td>
   <td style="text-align:right;"> 0.5620014 </td>
   <td style="text-align:right;"> 0.2735961 </td>
   <td style="text-align:right;"> 0.1909447 </td>
   <td style="text-align:right;"> 0.0921777 </td>
   <td style="text-align:right;"> 0.1683476 </td>
   <td style="text-align:right;"> 0.0439991 </td>
   <td style="text-align:right;"> 0.4026721 </td>
   <td style="text-align:right;"> 0.2288391 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.2172522 </td>
   <td style="text-align:right;"> 0.5437496 </td>
   <td style="text-align:right;"> 0.9240419 </td>
   <td style="text-align:right;"> 0.3903264 </td>
   <td style="text-align:right;"> 0.8742944 </td>
   <td style="text-align:right;"> 0.6324085 </td>
   <td style="text-align:right;"> 0.1809046 </td>
   <td style="text-align:right;"> 0.2052155 </td>
   <td style="text-align:right;"> 0.9456595 </td>
   <td style="text-align:right;"> 0.3638987 </td>
   <td style="text-align:right;"> 0.3213784 </td>
   <td style="text-align:right;"> 0.1615490 </td>
   <td style="text-align:right;"> 0.0734817 </td>
   <td style="text-align:right;"> 0.0286660 </td>
   <td style="text-align:right;"> 0.5680329 </td>
   <td style="text-align:right;"> 0.0281195 </td>
   <td style="text-align:right;"> 0.1982171 </td>
   <td style="text-align:right;"> 0.0713569 </td>
   <td style="text-align:right;"> 0.8122336 </td>
   <td style="text-align:right;"> 0.5943812 </td>
   <td style="text-align:right;"> 0.0975042 </td>
   <td style="text-align:right;"> 0.0334829 </td>
   <td style="text-align:right;"> 0.0201244 </td>
   <td style="text-align:right;"> 0.3411303 </td>
   <td style="text-align:right;"> 0.8706517 </td>
   <td style="text-align:right;"> 0.0201684 </td>
   <td style="text-align:right;"> 0.0461775 </td>
   <td style="text-align:right;"> 0.8186859 </td>
   <td style="text-align:right;"> 0.2234754 </td>
   <td style="text-align:right;"> 0.1605849 </td>
   <td style="text-align:right;"> 0.1148716 </td>
   <td style="text-align:right;"> 0.4625290 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6688718 </td>
   <td style="text-align:right;"> 0.7357528 </td>
   <td style="text-align:right;"> 0.9146012 </td>
   <td style="text-align:right;"> 0.5121105 </td>
   <td style="text-align:right;"> 0.9151110 </td>
   <td style="text-align:right;"> 0.7347897 </td>
   <td style="text-align:right;"> 0.3386730 </td>
   <td style="text-align:right;"> 0.7055080 </td>
   <td style="text-align:right;"> 0.1460815 </td>
   <td style="text-align:right;"> 0.7448671 </td>
   <td style="text-align:right;"> 0.4161712 </td>
   <td style="text-align:right;"> 0.7785247 </td>
   <td style="text-align:right;"> 0.9727465 </td>
   <td style="text-align:right;"> 0.9084417 </td>
   <td style="text-align:right;"> 0.1655484 </td>
   <td style="text-align:right;"> 0.1471424 </td>
   <td style="text-align:right;"> 0.3658783 </td>
   <td style="text-align:right;"> 0.5844535 </td>
   <td style="text-align:right;"> 0.9776100 </td>
   <td style="text-align:right;"> 0.2640119 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.2752774 </td>
   <td style="text-align:right;"> 0.4153787 </td>
   <td style="text-align:right;"> 0.4863892 </td>
   <td style="text-align:right;"> 0.3405155 </td>
   <td style="text-align:right;"> 0.8035884 </td>
   <td style="text-align:right;"> 0.5882264 </td>
   <td style="text-align:right;"> 0.5011244 </td>
   <td style="text-align:right;"> 0.0824842 </td>
   <td style="text-align:right;"> 0.0298981 </td>
   <td style="text-align:right;"> 0.4754469 </td>
   <td style="text-align:right;"> 0.3165080 </td>
   <td style="text-align:right;"> 0.3457467 </td>
   <td style="text-align:right;"> 0.4383894 </td>
   <td style="text-align:right;"> 0.5831174 </td>
   <td style="text-align:right;"> 0.5271920 </td>
   <td style="text-align:right;"> 0.0046074 </td>
   <td style="text-align:right;"> 0.0494925 </td>
   <td style="text-align:right;"> 0.0834253 </td>
   <td style="text-align:right;"> 0.2725647 </td>
   <td style="text-align:right;"> 0.0717840 </td>
   <td style="text-align:right;"> 0.2775516 </td>
   <td style="text-align:right;"> 0.0333481 </td>
   <td style="text-align:right;"> 0.0602219 </td>
   <td style="text-align:right;"> 0.1867152 </td>
   <td style="text-align:right;"> 0.2364474 </td>
   <td style="text-align:right;"> 0.1709506 </td>
   <td style="text-align:right;"> 0.0062182 </td>
   <td style="text-align:right;"> 0.1937085 </td>
   <td style="text-align:right;"> 0.0031241 </td>
   <td style="text-align:right;"> 0.9663449 </td>
   <td style="text-align:right;"> 0.1135380 </td>
   <td style="text-align:right;"> 0.2829966 </td>
   <td style="text-align:right;"> 0.4378247 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.1437541 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.2978460 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.8834358 </td>
   <td style="text-align:right;"> 0.3011440 </td>
   <td style="text-align:right;"> 0.4772483 </td>
   <td style="text-align:right;"> 0.3782618 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.7051462 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.2449623 </td>
   <td style="text-align:right;"> 0.4701513 </td>
   <td style="text-align:right;"> 0.3307208 </td>
   <td style="text-align:right;"> 0.3342780 </td>
   <td style="text-align:right;"> 0.3991050 </td>
   <td style="text-align:right;"> 0.0886043 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0212810 </td>
   <td style="text-align:right;"> 0.9889108 </td>
   <td style="text-align:right;"> 0.9101427 </td>
   <td style="text-align:right;"> 0.9998409 </td>
   <td style="text-align:right;"> 0.1731679 </td>
   <td style="text-align:right;"> 0.2968288 </td>
   <td style="text-align:right;"> 0.4258075 </td>
   <td style="text-align:right;"> 0.7956044 </td>
   <td style="text-align:right;"> 0.0301263 </td>
   <td style="text-align:right;"> 0.2384351 </td>
   <td style="text-align:right;"> 0.4469264 </td>
   <td style="text-align:right;"> 0.1698578 </td>
   <td style="text-align:right;"> 0.0435942 </td>
   <td style="text-align:right;"> 0.4648984 </td>
   <td style="text-align:right;"> 0.4420590 </td>
   <td style="text-align:right;"> 0.0059611 </td>
   <td style="text-align:right;"> 0.0030308 </td>
   <td style="text-align:right;"> 0.0047962 </td>
   <td style="text-align:right;"> 0.2504284 </td>
   <td style="text-align:right;"> 0.0005580 </td>
   <td style="text-align:right;"> 0.0138824 </td>
   <td style="text-align:right;"> 0.2673959 </td>
   <td style="text-align:right;"> 0.1704776 </td>
   <td style="text-align:right;"> 0.8190754 </td>
   <td style="text-align:right;"> 0.0282837 </td>
   <td style="text-align:right;"> 0.0319594 </td>
   <td style="text-align:right;"> 0.0047349 </td>
   <td style="text-align:right;"> 0.0386002 </td>
   <td style="text-align:right;"> 0.0015975 </td>
   <td style="text-align:right;"> 0.5037617 </td>
   <td style="text-align:right;"> 0.0113254 </td>
   <td style="text-align:right;"> 0.0321779 </td>
   <td style="text-align:right;"> 0.5291227 </td>
   <td style="text-align:right;"> 0.0136077 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.0365837 </td>
   <td style="text-align:right;"> 0.0999526 </td>
   <td style="text-align:right;"> 0.0646294 </td>
   <td style="text-align:right;"> 0.6614673 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.0886043 </td>
   <td style="text-align:right;"> 0.3798126 </td>
   <td style="text-align:right;"> 0.2102899 </td>
   <td style="text-align:right;"> 0.1618261 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.1252776 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.1683476 </td>
   <td style="text-align:right;"> 0.0820271 </td>
   <td style="text-align:right;"> 0.2824454 </td>
   <td style="text-align:right;"> 0.3341755 </td>
   <td style="text-align:right;"> 0.0679049 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.1027681 </td>
   <td style="text-align:right;"> 0.6099165 </td>
   <td style="text-align:right;"> 0.7909334 </td>
   <td style="text-align:right;"> 0.1006361 </td>
   <td style="text-align:right;"> 0.4164781 </td>
   <td style="text-align:right;"> 0.6810747 </td>
   <td style="text-align:right;"> 0.1164003 </td>
   <td style="text-align:right;"> 0.2369801 </td>
   <td style="text-align:right;"> 0.9348276 </td>
   <td style="text-align:right;"> 0.0856471 </td>
   <td style="text-align:right;"> 0.9694980 </td>
   <td style="text-align:right;"> 0.0374958 </td>
   <td style="text-align:right;"> 0.0041655 </td>
   <td style="text-align:right;"> 0.0734100 </td>
   <td style="text-align:right;"> 0.0210011 </td>
   <td style="text-align:right;"> 0.4175528 </td>
   <td style="text-align:right;"> 0.4163539 </td>
   <td style="text-align:right;"> 0.8960122 </td>
   <td style="text-align:right;"> 0.8102904 </td>
   <td style="text-align:right;"> 0.0012254 </td>
   <td style="text-align:right;"> 0.1763370 </td>
   <td style="text-align:right;"> 0.8046084 </td>
   <td style="text-align:right;"> 0.8494651 </td>
   <td style="text-align:right;"> 0.1794559 </td>
   <td style="text-align:right;"> 0.3903940 </td>
   <td style="text-align:right;"> 0.6432672 </td>
   <td style="text-align:right;"> 0.0278584 </td>
   <td style="text-align:right;"> 0.6839354 </td>
   <td style="text-align:right;"> 0.0851298 </td>
   <td style="text-align:right;"> 0.0754377 </td>
   <td style="text-align:right;"> 0.0966410 </td>
   <td style="text-align:right;"> 0.2346617 </td>
   <td style="text-align:right;"> 0.8251158 </td>
   <td style="text-align:right;"> 0.6277371 </td>
   <td style="text-align:right;"> 0.0532986 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.7994115 </td>
   <td style="text-align:right;"> 0.1802483 </td>
   <td style="text-align:right;"> 0.6428159 </td>
   <td style="text-align:right;"> 0.9798903 </td>
   <td style="text-align:right;"> 0.8390914 </td>
   <td style="text-align:right;"> 0.6573562 </td>
   <td style="text-align:right;"> 0.2203442 </td>
   <td style="text-align:right;"> 0.3315460 </td>
   <td style="text-align:right;"> 0.4217927 </td>
   <td style="text-align:right;"> 0.0000030 </td>
   <td style="text-align:right;"> 0.8082713 </td>
   <td style="text-align:right;"> 0.0999526 </td>
   <td style="text-align:right;"> 0.3070777 </td>
   <td style="text-align:right;"> 0.8520815 </td>
   <td style="text-align:right;"> 0.8672454 </td>
   <td style="text-align:right;"> 0.8082713 </td>
   <td style="text-align:right;"> 0.9112983 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.0098348 </td>
   <td style="text-align:right;"> 0.3489275 </td>
   <td style="text-align:right;"> 0.3451777 </td>
   <td style="text-align:right;"> 0.2825090 </td>
   <td style="text-align:right;"> 0.6236386 </td>
   <td style="text-align:right;"> 0.6922061 </td>
   <td style="text-align:right;"> 0.6780016 </td>
   <td style="text-align:right;"> 0.8588498 </td>
   <td style="text-align:right;"> 0.0462935 </td>
   <td style="text-align:right;"> 0.4222197 </td>
   <td style="text-align:right;"> 0.3404631 </td>
   <td style="text-align:right;"> 0.3906473 </td>
   <td style="text-align:right;"> 0.3916520 </td>
   <td style="text-align:right;"> 0.8387850 </td>
   <td style="text-align:right;"> 0.6951464 </td>
   <td style="text-align:right;"> 0.0161675 </td>
   <td style="text-align:right;"> 0.0000005 </td>
   <td style="text-align:right;"> 0.0000251 </td>
   <td style="text-align:right;"> 0.1288539 </td>
   <td style="text-align:right;"> 0.0405821 </td>
   <td style="text-align:right;"> 0.1726904 </td>
   <td style="text-align:right;"> 0.0327910 </td>
   <td style="text-align:right;"> 0.0135445 </td>
   <td style="text-align:right;"> 0.1885373 </td>
   <td style="text-align:right;"> 0.3620007 </td>
   <td style="text-align:right;"> 0.0336480 </td>
   <td style="text-align:right;"> 0.0710108 </td>
   <td style="text-align:right;"> 0.2084213 </td>
   <td style="text-align:right;"> 0.0174980 </td>
   <td style="text-align:right;"> 0.5722324 </td>
   <td style="text-align:right;"> 0.0853577 </td>
   <td style="text-align:right;"> 0.0695974 </td>
   <td style="text-align:right;"> 0.2638595 </td>
   <td style="text-align:right;"> 0.0190691 </td>
   <td style="text-align:right;"> 0.0009823 </td>
   <td style="text-align:right;"> 0.6194984 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.0999526 </td>
   <td style="text-align:right;"> 0.6903275 </td>
   <td style="text-align:right;"> 0.8016883 </td>
   <td style="text-align:right;"> 0.0316554 </td>
   <td style="text-align:right;"> 0.5479784 </td>
   <td style="text-align:right;"> 0.4144646 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.7848482 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.2131754 </td>
   <td style="text-align:right;"> 0.2603510 </td>
   <td style="text-align:right;"> 0.0020316 </td>
   <td style="text-align:right;"> 0.1929722 </td>
   <td style="text-align:right;"> 0.1705204 </td>
   <td style="text-align:right;"> 0.0004749 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.0655140 </td>
   <td style="text-align:right;"> 0.8995500 </td>
   <td style="text-align:right;"> 0.6136475 </td>
   <td style="text-align:right;"> 0.4029134 </td>
   <td style="text-align:right;"> 0.6288212 </td>
   <td style="text-align:right;"> 0.6353082 </td>
   <td style="text-align:right;"> 0.0726584 </td>
   <td style="text-align:right;"> 0.6042235 </td>
   <td style="text-align:right;"> 0.4866210 </td>
   <td style="text-align:right;"> 0.0950031 </td>
   <td style="text-align:right;"> 0.5556538 </td>
   <td style="text-align:right;"> 0.1405648 </td>
   <td style="text-align:right;"> 0.1510948 </td>
   <td style="text-align:right;"> 0.2276385 </td>
   <td style="text-align:right;"> 0.0594434 </td>
   <td style="text-align:right;"> 0.0190907 </td>
   <td style="text-align:right;"> 0.0453884 </td>
   <td style="text-align:right;"> 0.2034114 </td>
   <td style="text-align:right;"> 0.3971937 </td>
   <td style="text-align:right;"> 0.0058144 </td>
   <td style="text-align:right;"> 0.6569681 </td>
   <td style="text-align:right;"> 0.1958927 </td>
   <td style="text-align:right;"> 0.1229038 </td>
   <td style="text-align:right;"> 0.3812233 </td>
   <td style="text-align:right;"> 0.5855562 </td>
   <td style="text-align:right;"> 0.3980818 </td>
   <td style="text-align:right;"> 0.0223533 </td>
   <td style="text-align:right;"> 0.3844695 </td>
   <td style="text-align:right;"> 0.0023539 </td>
   <td style="text-align:right;"> 0.7951242 </td>
   <td style="text-align:right;"> 0.0681656 </td>
   <td style="text-align:right;"> 0.2917739 </td>
   <td style="text-align:right;"> 0.8294382 </td>
   <td style="text-align:right;"> 0.0071676 </td>
   <td style="text-align:right;"> 0.0084865 </td>
   <td style="text-align:right;"> 0.0320618 </td>
   <td style="text-align:right;"> 0.0368348 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.5989343 </td>
   <td style="text-align:right;"> 0.9146012 </td>
   <td style="text-align:right;"> 0.7642200 </td>
   <td style="text-align:right;"> 0.3135264 </td>
   <td style="text-align:right;"> 0.1427922 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.0333065 </td>
   <td style="text-align:right;"> 0.1519730 </td>
   <td style="text-align:right;"> 0.2715956 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.5038062 </td>
   <td style="text-align:right;"> 0.3774174 </td>
   <td style="text-align:right;"> 0.7575202 </td>
   <td style="text-align:right;"> 0.3315460 </td>
   <td style="text-align:right;"> 0.2323097 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 0.0545723 </td>
   <td style="text-align:right;"> 0.6200275 </td>
   <td style="text-align:right;"> 0.6185655 </td>
   <td style="text-align:right;"> 0.0780015 </td>
   <td style="text-align:right;"> 0.0383423 </td>
   <td style="text-align:right;"> 0.0153444 </td>
   <td style="text-align:right;"> 0.5202712 </td>
   <td style="text-align:right;"> 0.4865656 </td>
   <td style="text-align:right;"> 0.0202479 </td>
   <td style="text-align:right;"> 0.8182540 </td>
   <td style="text-align:right;"> 0.2122328 </td>
   <td style="text-align:right;"> 0.3752755 </td>
   <td style="text-align:right;"> 0.6083593 </td>
   <td style="text-align:right;"> 0.6250946 </td>
   <td style="text-align:right;"> 0.5690323 </td>
   <td style="text-align:right;"> 0.1032293 </td>
   <td style="text-align:right;"> 0.0125855 </td>
   <td style="text-align:right;"> 0.0306902 </td>
   <td style="text-align:right;"> 0.2273580 </td>
   <td style="text-align:right;"> 0.1071735 </td>
   <td style="text-align:right;"> 0.0208881 </td>
   <td style="text-align:right;"> 0.3160191 </td>
   <td style="text-align:right;"> 0.3023550 </td>
   <td style="text-align:right;"> 0.6848991 </td>
   <td style="text-align:right;"> 0.2844498 </td>
   <td style="text-align:right;"> 0.3079922 </td>
   <td style="text-align:right;"> 0.0362604 </td>
   <td style="text-align:right;"> 0.0106902 </td>
   <td style="text-align:right;"> 0.2554290 </td>
   <td style="text-align:right;"> 0.5944471 </td>
   <td style="text-align:right;"> 0.4857441 </td>
   <td style="text-align:right;"> 0.4873566 </td>
   <td style="text-align:right;"> 0.5246975 </td>
   <td style="text-align:right;"> 0.0935902 </td>
   <td style="text-align:right;"> 0.0023113 </td>
   <td style="text-align:right;"> 0.3983779 </td>
   <td style="text-align:right;"> 0.0083147 </td>
   <td style="text-align:right;"> 0.3424965 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.9380902 </td>
   <td style="text-align:right;"> 0.4974820 </td>
   <td style="text-align:right;"> 0.3357679 </td>
   <td style="text-align:right;"> 0.9629801 </td>
   <td style="text-align:right;"> 0.7588192 </td>
   <td style="text-align:right;"> 0.5931859 </td>
   <td style="text-align:right;"> 0.5934721 </td>
   <td style="text-align:right;"> 0.1870723 </td>
   <td style="text-align:right;"> 0.4833604 </td>
   <td style="text-align:right;"> 0.3774174 </td>
   <td style="text-align:right;"> 0.2201507 </td>
   <td style="text-align:right;"> 0.4095298 </td>
   <td style="text-align:right;"> 0.6780173 </td>
   <td style="text-align:right;"> 0.1237702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 0.1871999 </td>
   <td style="text-align:right;"> 0.3754998 </td>
   <td style="text-align:right;"> 0.2842474 </td>
   <td style="text-align:right;"> 0.4254902 </td>
   <td style="text-align:right;"> 0.6838472 </td>
   <td style="text-align:right;"> 0.2389744 </td>
   <td style="text-align:right;"> 0.6828410 </td>
   <td style="text-align:right;"> 0.1459182 </td>
   <td style="text-align:right;"> 0.9614231 </td>
   <td style="text-align:right;"> 0.4838598 </td>
   <td style="text-align:right;"> 0.0332284 </td>
   <td style="text-align:right;"> 0.2309434 </td>
   <td style="text-align:right;"> 0.1180242 </td>
   <td style="text-align:right;"> 0.1703009 </td>
   <td style="text-align:right;"> 0.7435586 </td>
   <td style="text-align:right;"> 0.3813780 </td>
   <td style="text-align:right;"> 0.6994536 </td>
   <td style="text-align:right;"> 0.0908291 </td>
   <td style="text-align:right;"> 0.2719253 </td>
   <td style="text-align:right;"> 0.8553603 </td>
   <td style="text-align:right;"> 0.2548294 </td>
   <td style="text-align:right;"> 0.4100186 </td>
   <td style="text-align:right;"> 0.7067705 </td>
   <td style="text-align:right;"> 0.4834797 </td>
   <td style="text-align:right;"> 0.0600685 </td>
   <td style="text-align:right;"> 0.0021181 </td>
   <td style="text-align:right;"> 0.5717516 </td>
   <td style="text-align:right;"> 0.1855377 </td>
   <td style="text-align:right;"> 0.2399264 </td>
   <td style="text-align:right;"> 0.1803442 </td>
   <td style="text-align:right;"> 0.0791538 </td>
   <td style="text-align:right;"> 0.0173237 </td>
   <td style="text-align:right;"> 0.1243603 </td>
   <td style="text-align:right;"> 0.9301114 </td>
   <td style="text-align:right;"> 0.4267376 </td>
   <td style="text-align:right;"> 0.9621129 </td>
   <td style="text-align:right;"> 0.4590060 </td>
   <td style="text-align:right;"> 0.8237186 </td>
   <td style="text-align:right;"> 0.8754422 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0931793 </td>
   <td style="text-align:right;"> 0.1260899 </td>
   <td style="text-align:right;"> 0.8177333 </td>
   <td style="text-align:right;"> 0.1805635 </td>
   <td style="text-align:right;"> 0.8819753 </td>
   <td style="text-align:right;"> 0.8773829 </td>
   <td style="text-align:right;"> 0.6043755 </td>
   <td style="text-align:right;"> 0.4701513 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.2978460 </td>
   <td style="text-align:right;"> 0.4986542 </td>
   <td style="text-align:right;"> 0.7946947 </td>
   <td style="text-align:right;"> 0.4048325 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 0.8181650 </td>
   <td style="text-align:right;"> 0.9240932 </td>
   <td style="text-align:right;"> 0.2372382 </td>
   <td style="text-align:right;"> 0.3655945 </td>
   <td style="text-align:right;"> 0.7815649 </td>
   <td style="text-align:right;"> 0.3409691 </td>
   <td style="text-align:right;"> 0.4162819 </td>
   <td style="text-align:right;"> 0.4901373 </td>
   <td style="text-align:right;"> 0.7174691 </td>
   <td style="text-align:right;"> 0.3136507 </td>
   <td style="text-align:right;"> 0.8335797 </td>
   <td style="text-align:right;"> 0.7244657 </td>
   <td style="text-align:right;"> 0.1525589 </td>
   <td style="text-align:right;"> 0.1852904 </td>
   <td style="text-align:right;"> 0.4199008 </td>
   <td style="text-align:right;"> 0.7245207 </td>
   <td style="text-align:right;"> 0.4973819 </td>
   <td style="text-align:right;"> 0.6862532 </td>
   <td style="text-align:right;"> 0.8916896 </td>
   <td style="text-align:right;"> 0.9836134 </td>
   <td style="text-align:right;"> 0.5211452 </td>
   <td style="text-align:right;"> 0.8049824 </td>
   <td style="text-align:right;"> 0.6322832 </td>
   <td style="text-align:right;"> 0.8688843 </td>
   <td style="text-align:right;"> 0.0217169 </td>
   <td style="text-align:right;"> 0.0754089 </td>
   <td style="text-align:right;"> 0.9028121 </td>
   <td style="text-align:right;"> 0.9727981 </td>
   <td style="text-align:right;"> 0.2964810 </td>
   <td style="text-align:right;"> 0.0257093 </td>
   <td style="text-align:right;"> 0.0757667 </td>
   <td style="text-align:right;"> 0.0265994 </td>
   <td style="text-align:right;"> 0.4879166 </td>
   <td style="text-align:right;"> 0.7558569 </td>
   <td style="text-align:right;"> 0.9251306 </td>
   <td style="text-align:right;"> 0.6880793 </td>
   <td style="text-align:right;"> 0.6266496 </td>
   <td style="text-align:right;"> 0.5667872 </td>
   <td style="text-align:right;"> 0.2512681 </td>
   <td style="text-align:right;"> 0.0068187 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6285719 </td>
   <td style="text-align:right;"> 0.7038012 </td>
   <td style="text-align:right;"> 0.3070777 </td>
   <td style="text-align:right;"> 0.8956092 </td>
   <td style="text-align:right;"> 0.7506619 </td>
   <td style="text-align:right;"> 0.8141569 </td>
   <td style="text-align:right;"> 0.7631012 </td>
   <td style="text-align:right;"> 0.4154550 </td>
   <td style="text-align:right;"> 0.9151110 </td>
   <td style="text-align:right;"> 0.4061416 </td>
   <td style="text-align:right;"> 0.9831447 </td>
   <td style="text-align:right;"> 0.9456724 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0081785 </td>
   <td style="text-align:right;"> 0.9981915 </td>
   <td style="text-align:right;"> 0.2012160 </td>
   <td style="text-align:right;"> 0.9863220 </td>
   <td style="text-align:right;"> 0.8095127 </td>
   <td style="text-align:right;"> 0.6835525 </td>
   <td style="text-align:right;"> 0.7018278 </td>
   <td style="text-align:right;"> 0.4074884 </td>
   <td style="text-align:right;"> 0.3792757 </td>
   <td style="text-align:right;"> 0.4403006 </td>
   <td style="text-align:right;"> 0.1347240 </td>
   <td style="text-align:right;"> 0.0718280 </td>
   <td style="text-align:right;"> 0.0306305 </td>
   <td style="text-align:right;"> 0.1925719 </td>
   <td style="text-align:right;"> 0.5278706 </td>
   <td style="text-align:right;"> 0.0067523 </td>
   <td style="text-align:right;"> 0.0031809 </td>
   <td style="text-align:right;"> 0.0000293 </td>
   <td style="text-align:right;"> 0.3078956 </td>
   <td style="text-align:right;"> 0.0368569 </td>
   <td style="text-align:right;"> 0.0869405 </td>
   <td style="text-align:right;"> 0.0450853 </td>
   <td style="text-align:right;"> 0.0210356 </td>
   <td style="text-align:right;"> 0.2394247 </td>
   <td style="text-align:right;"> 0.1798384 </td>
   <td style="text-align:right;"> 0.0000120 </td>
   <td style="text-align:right;"> 0.0265490 </td>
   <td style="text-align:right;"> 0.2604273 </td>
   <td style="text-align:right;"> 0.0069622 </td>
   <td style="text-align:right;"> 0.3517595 </td>
   <td style="text-align:right;"> 0.0049700 </td>
   <td style="text-align:right;"> 0.0065159 </td>
   <td style="text-align:right;"> 0.0201546 </td>
   <td style="text-align:right;"> 0.0959378 </td>
   <td style="text-align:right;"> 0.0054834 </td>
   <td style="text-align:right;"> 0.4157933 </td>
   <td style="text-align:right;"> 0.0007351 </td>
   <td style="text-align:right;"> 0.1048880 </td>
   <td style="text-align:right;"> 0.1220753 </td>
   <td style="text-align:right;"> 0.0145154 </td>
   <td style="text-align:right;"> 0.3817959 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.6412448 </td>
   <td style="text-align:right;"> 0.1441993 </td>
   <td style="text-align:right;"> 0.3365619 </td>
   <td style="text-align:right;"> 0.7308369 </td>
   <td style="text-align:right;"> 0.3279466 </td>
   <td style="text-align:right;"> 0.0670841 </td>
   <td style="text-align:right;"> 0.0477306 </td>
   <td style="text-align:right;"> 0.0030454 </td>
   <td style="text-align:right;"> 0.2494228 </td>
   <td style="text-align:right;"> 0.3812442 </td>
   <td style="text-align:right;"> 0.0101833 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.0515077 </td>
   <td style="text-align:right;"> 0.4123562 </td>
   <td style="text-align:right;"> 0.2225188 </td>
   <td style="text-align:right;"> 0.2734599 </td>
   <td style="text-align:right;"> 0.5005644 </td>
   <td style="text-align:right;"> 0.4367385 </td>
   <td style="text-align:right;"> 0.0002473 </td>
   <td style="text-align:right;"> 0.4687149 </td>
   <td style="text-align:right;"> 0.3188194 </td>
   <td style="text-align:right;"> 0.0000265 </td>
   <td style="text-align:right;"> 0.5176319 </td>
   <td style="text-align:right;"> 0.3012858 </td>
   <td style="text-align:right;"> 0.1157176 </td>
   <td style="text-align:right;"> 0.5873951 </td>
   <td style="text-align:right;"> 0.1634234 </td>
   <td style="text-align:right;"> 0.5594778 </td>
   <td style="text-align:right;"> 0.1969628 </td>
   <td style="text-align:right;"> 0.5130252 </td>
   <td style="text-align:right;"> 0.5096205 </td>
   <td style="text-align:right;"> 0.0567595 </td>
   <td style="text-align:right;"> 0.8264936 </td>
   <td style="text-align:right;"> 0.3758716 </td>
   <td style="text-align:right;"> 0.9410731 </td>
   <td style="text-align:right;"> 0.9328352 </td>
   <td style="text-align:right;"> 0.1375482 </td>
   <td style="text-align:right;"> 0.8317623 </td>
   <td style="text-align:right;"> 0.4345476 </td>
   <td style="text-align:right;"> 0.7366786 </td>
   <td style="text-align:right;"> 0.0232946 </td>
   <td style="text-align:right;"> 0.5389639 </td>
   <td style="text-align:right;"> 0.0506113 </td>
   <td style="text-align:right;"> 0.0163106 </td>
   <td style="text-align:right;"> 0.5394704 </td>
   <td style="text-align:right;"> 0.2310766 </td>
   <td style="text-align:right;"> 0.1485624 </td>
   <td style="text-align:right;"> 0.0518081 </td>
   <td style="text-align:right;"> 0.2970536 </td>
   <td style="text-align:right;"> 0.0185091 </td>
   <td style="text-align:right;"> 0.9308341 </td>
   <td style="text-align:right;"> 0.6551361 </td>
   <td style="text-align:right;"> 0.4795272 </td>
   <td style="text-align:right;"> 0.3955429 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.1246680 </td>
   <td style="text-align:right;"> 0.0608425 </td>
   <td style="text-align:right;"> 0.2205568 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.3554490 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.6875233 </td>
   <td style="text-align:right;"> 0.2323097 </td>
   <td style="text-align:right;"> 0.2715956 </td>
   <td style="text-align:right;"> 0.5931859 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 0.0071002 </td>
   <td style="text-align:right;"> 0.9409952 </td>
   <td style="text-align:right;"> 0.5126330 </td>
   <td style="text-align:right;"> 0.3916507 </td>
   <td style="text-align:right;"> 0.6517238 </td>
   <td style="text-align:right;"> 0.3809028 </td>
   <td style="text-align:right;"> 0.2603617 </td>
   <td style="text-align:right;"> 0.3556906 </td>
   <td style="text-align:right;"> 0.3062493 </td>
   <td style="text-align:right;"> 0.1216053 </td>
   <td style="text-align:right;"> 0.6549169 </td>
   <td style="text-align:right;"> 0.0459886 </td>
   <td style="text-align:right;"> 0.0030552 </td>
   <td style="text-align:right;"> 0.0149309 </td>
   <td style="text-align:right;"> 0.2635114 </td>
   <td style="text-align:right;"> 0.0125551 </td>
   <td style="text-align:right;"> 0.1288726 </td>
   <td style="text-align:right;"> 0.0471832 </td>
   <td style="text-align:right;"> 0.8074909 </td>
   <td style="text-align:right;"> 0.0757724 </td>
   <td style="text-align:right;"> 0.1583412 </td>
   <td style="text-align:right;"> 0.3542473 </td>
   <td style="text-align:right;"> 0.7155296 </td>
   <td style="text-align:right;"> 0.3234984 </td>
   <td style="text-align:right;"> 0.0535806 </td>
   <td style="text-align:right;"> 0.0634311 </td>
   <td style="text-align:right;"> 0.0651569 </td>
   <td style="text-align:right;"> 0.1208058 </td>
   <td style="text-align:right;"> 0.0047048 </td>
   <td style="text-align:right;"> 0.0479265 </td>
   <td style="text-align:right;"> 0.0000909 </td>
   <td style="text-align:right;"> 0.0045118 </td>
   <td style="text-align:right;"> 0.1808878 </td>
   <td style="text-align:right;"> 0.1476813 </td>
   <td style="text-align:right;"> 0.0450030 </td>
   <td style="text-align:right;"> 0.1185919 </td>
   <td style="text-align:right;"> 0.1780881 </td>
   <td style="text-align:right;"> 0.0526555 </td>
   <td style="text-align:right;"> 0.5589966 </td>
   <td style="text-align:right;"> 0.0323652 </td>
   <td style="text-align:right;"> 0.0996108 </td>
   <td style="text-align:right;"> 0.0192545 </td>
   <td style="text-align:right;"> 0.0133086 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0858852 </td>
   <td style="text-align:right;"> 0.4447703 </td>
   <td style="text-align:right;"> 0.2371021 </td>
   <td style="text-align:right;"> 0.0974954 </td>
   <td style="text-align:right;"> 0.0843294 </td>
   <td style="text-align:right;"> 0.3979588 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.7786512 </td>
   <td style="text-align:right;"> 0.2205568 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 0.0531625 </td>
   <td style="text-align:right;"> 0.2276423 </td>
   <td style="text-align:right;"> 0.8295014 </td>
   <td style="text-align:right;"> 0.9934747 </td>
   <td style="text-align:right;"> 0.3438793 </td>
   <td style="text-align:right;"> 0.4609782 </td>
   <td style="text-align:right;"> 0.0856219 </td>
   <td style="text-align:right;"> 0.9940144 </td>
   <td style="text-align:right;"> 0.0722156 </td>
   <td style="text-align:right;"> 0.0381118 </td>
   <td style="text-align:right;"> 0.6860117 </td>
   <td style="text-align:right;"> 0.5013844 </td>
   <td style="text-align:right;"> 0.1267018 </td>
   <td style="text-align:right;"> 0.2909024 </td>
   <td style="text-align:right;"> 0.2174784 </td>
   <td style="text-align:right;"> 0.0120060 </td>
   <td style="text-align:right;"> 0.0256288 </td>
   <td style="text-align:right;"> 0.1036538 </td>
   <td style="text-align:right;"> 0.2342510 </td>
   <td style="text-align:right;"> 0.0806096 </td>
   <td style="text-align:right;"> 0.2656116 </td>
   <td style="text-align:right;"> 0.0488288 </td>
   <td style="text-align:right;"> 0.1629307 </td>
   <td style="text-align:right;"> 0.3749954 </td>
   <td style="text-align:right;"> 0.1947633 </td>
   <td style="text-align:right;"> 0.4351214 </td>
   <td style="text-align:right;"> 0.0851067 </td>
   <td style="text-align:right;"> 0.2339980 </td>
   <td style="text-align:right;"> 0.0002404 </td>
   <td style="text-align:right;"> 0.7364357 </td>
   <td style="text-align:right;"> 0.0166451 </td>
   <td style="text-align:right;"> 0.1008599 </td>
   <td style="text-align:right;"> 0.5858709 </td>
   <td style="text-align:right;"> 0.0048566 </td>
   <td style="text-align:right;"> 0.0241917 </td>
   <td style="text-align:right;"> 0.1861494 </td>
   <td style="text-align:right;"> 0.0369053 </td>
   <td style="text-align:right;"> 0.0008218 </td>
   <td style="text-align:right;"> 0.3344698 </td>
   <td style="text-align:right;"> 0.7539673 </td>
   <td style="text-align:right;"> 0.7805709 </td>
   <td style="text-align:right;"> 0.1226082 </td>
   <td style="text-align:right;"> 0.0021063 </td>
   <td style="text-align:right;"> 0.0052354 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4541924 </td>
   <td style="text-align:right;"> 0.0378701 </td>
   <td style="text-align:right;"> 0.1417060 </td>
   <td style="text-align:right;"> 0.2803367 </td>
   <td style="text-align:right;"> 0.3304656 </td>
   <td style="text-align:right;"> 0.2150127 </td>
   <td style="text-align:right;"> 0.1890541 </td>
   <td style="text-align:right;"> 0.1391036 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.1178798 </td>
   <td style="text-align:right;"> 0.5879779 </td>
   <td style="text-align:right;"> 0.9769096 </td>
   <td style="text-align:right;"> 0.1430429 </td>
   <td style="text-align:right;"> 0.3287784 </td>
   <td style="text-align:right;"> 0.4948954 </td>
   <td style="text-align:right;"> 0.0911093 </td>
   <td style="text-align:right;"> 0.4803987 </td>
   <td style="text-align:right;"> 0.9710905 </td>
   <td style="text-align:right;"> 0.0871344 </td>
   <td style="text-align:right;"> 0.8008204 </td>
   <td style="text-align:right;"> 0.0276016 </td>
   <td style="text-align:right;"> 0.0188230 </td>
   <td style="text-align:right;"> 0.1612327 </td>
   <td style="text-align:right;"> 0.0145488 </td>
   <td style="text-align:right;"> 0.5167118 </td>
   <td style="text-align:right;"> 0.4355547 </td>
   <td style="text-align:right;"> 0.9896861 </td>
   <td style="text-align:right;"> 0.8675262 </td>
   <td style="text-align:right;"> 0.0008094 </td>
   <td style="text-align:right;"> 0.2933670 </td>
   <td style="text-align:right;"> 0.8084769 </td>
   <td style="text-align:right;"> 0.8922665 </td>
   <td style="text-align:right;"> 0.2209295 </td>
   <td style="text-align:right;"> 0.4673313 </td>
   <td style="text-align:right;"> 0.7793924 </td>
   <td style="text-align:right;"> 0.0222343 </td>
   <td style="text-align:right;"> 0.6852374 </td>
   <td style="text-align:right;"> 0.1146440 </td>
   <td style="text-align:right;"> 0.1604346 </td>
   <td style="text-align:right;"> 0.1864630 </td>
   <td style="text-align:right;"> 0.3083259 </td>
   <td style="text-align:right;"> 0.9488166 </td>
   <td style="text-align:right;"> 0.4820376 </td>
   <td style="text-align:right;"> 0.0522369 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.5977461 </td>
   <td style="text-align:right;"> 0.0220570 </td>
   <td style="text-align:right;"> 0.3350663 </td>
   <td style="text-align:right;"> 0.7474946 </td>
   <td style="text-align:right;"> 0.5463817 </td>
   <td style="text-align:right;"> 0.5170640 </td>
   <td style="text-align:right;"> 0.0522824 </td>
   <td style="text-align:right;"> 0.2036648 </td>
   <td style="text-align:right;"> 0.2106161 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.7943741 </td>
   <td style="text-align:right;"> 0.1471424 </td>
   <td style="text-align:right;"> 0.4061416 </td>
   <td style="text-align:right;"> 0.9126519 </td>
   <td style="text-align:right;"> 0.9129462 </td>
   <td style="text-align:right;"> 0.8148945 </td>
   <td style="text-align:right;"> 0.9146012 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 0.0738278 </td>
   <td style="text-align:right;"> 0.2184494 </td>
   <td style="text-align:right;"> 0.7708085 </td>
   <td style="text-align:right;"> 0.3534787 </td>
   <td style="text-align:right;"> 0.8423484 </td>
   <td style="text-align:right;"> 0.8900850 </td>
   <td style="text-align:right;"> 0.1843571 </td>
   <td style="text-align:right;"> 0.7226379 </td>
   <td style="text-align:right;"> 0.0069113 </td>
   <td style="text-align:right;"> 0.0757873 </td>
   <td style="text-align:right;"> 0.6201531 </td>
   <td style="text-align:right;"> 0.7764427 </td>
   <td style="text-align:right;"> 0.4861053 </td>
   <td style="text-align:right;"> 0.8940001 </td>
   <td style="text-align:right;"> 0.4014175 </td>
   <td style="text-align:right;"> 0.0839127 </td>
   <td style="text-align:right;"> 0.0303174 </td>
   <td style="text-align:right;"> 0.0539628 </td>
   <td style="text-align:right;"> 0.0095782 </td>
   <td style="text-align:right;"> 0.2751066 </td>
   <td style="text-align:right;"> 0.1004713 </td>
   <td style="text-align:right;"> 0.0251791 </td>
   <td style="text-align:right;"> 0.2244389 </td>
   <td style="text-align:right;"> 0.4350723 </td>
   <td style="text-align:right;"> 0.0354280 </td>
   <td style="text-align:right;"> 0.2822244 </td>
   <td style="text-align:right;"> 0.1609974 </td>
   <td style="text-align:right;"> 0.0233356 </td>
   <td style="text-align:right;"> 0.0059644 </td>
   <td style="text-align:right;"> 0.6405169 </td>
   <td style="text-align:right;"> 0.1365540 </td>
   <td style="text-align:right;"> 0.0804851 </td>
   <td style="text-align:right;"> 0.7990068 </td>
   <td style="text-align:right;"> 0.0059257 </td>
   <td style="text-align:right;"> 0.0141824 </td>
   <td style="text-align:right;"> 0.6404042 </td>
   <td style="text-align:right;"> 0.0136892 </td>
   <td style="text-align:right;"> 0.0788292 </td>
   <td style="text-align:right;"> 0.0351609 </td>
   <td style="text-align:right;"> 0.3526255 </td>
   <td style="text-align:right;"> 0.6469534 </td>
   <td style="text-align:right;"> 0.1123617 </td>
   <td style="text-align:right;"> 0.0363531 </td>
   <td style="text-align:right;"> 0.0599884 </td>
   <td style="text-align:right;"> 0.0010443 </td>
   <td style="text-align:right;"> 0.6116335 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.3812442 </td>
   <td style="text-align:right;"> 0.1927086 </td>
   <td style="text-align:right;"> 0.1929830 </td>
   <td style="text-align:right;"> 0.1265585 </td>
   <td style="text-align:right;"> 0.1062989 </td>
   <td style="text-align:right;"> 0.0704253 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.0402320 </td>
   <td style="text-align:right;"> 0.8985488 </td>
   <td style="text-align:right;"> 0.8980177 </td>
   <td style="text-align:right;"> 0.4767392 </td>
   <td style="text-align:right;"> 0.8666870 </td>
   <td style="text-align:right;"> 0.7127354 </td>
   <td style="text-align:right;"> 0.3675530 </td>
   <td style="text-align:right;"> 0.2264721 </td>
   <td style="text-align:right;"> 0.5429312 </td>
   <td style="text-align:right;"> 0.1979896 </td>
   <td style="text-align:right;"> 0.5848226 </td>
   <td style="text-align:right;"> 0.0508294 </td>
   <td style="text-align:right;"> 0.0000850 </td>
   <td style="text-align:right;"> 0.0057073 </td>
   <td style="text-align:right;"> 0.1147013 </td>
   <td style="text-align:right;"> 0.0029444 </td>
   <td style="text-align:right;"> 0.0361168 </td>
   <td style="text-align:right;"> 0.0471377 </td>
   <td style="text-align:right;"> 0.6512195 </td>
   <td style="text-align:right;"> 0.0019923 </td>
   <td style="text-align:right;"> 0.0149703 </td>
   <td style="text-align:right;"> 0.0435018 </td>
   <td style="text-align:right;"> 0.0262186 </td>
   <td style="text-align:right;"> 0.9530192 </td>
   <td style="text-align:right;"> 0.1841555 </td>
   <td style="text-align:right;"> 0.0096074 </td>
   <td style="text-align:right;"> 0.0005545 </td>
   <td style="text-align:right;"> 0.4564189 </td>
   <td style="text-align:right;"> 0.0001000 </td>
   <td style="text-align:right;"> 0.0398931 </td>
   <td style="text-align:right;"> 0.0007246 </td>
   <td style="text-align:right;"> 0.0376901 </td>
   <td style="text-align:right;"> 0.0248756 </td>
   <td style="text-align:right;"> 0.0635082 </td>
   <td style="text-align:right;"> 0.0046703 </td>
   <td style="text-align:right;"> 0.0084656 </td>
   <td style="text-align:right;"> 0.0469672 </td>
   <td style="text-align:right;"> 0.0098583 </td>
   <td style="text-align:right;"> 0.2383699 </td>
   <td style="text-align:right;"> 0.2232273 </td>
   <td style="text-align:right;"> 0.5645504 </td>
   <td style="text-align:right;"> 0.0026288 </td>
   <td style="text-align:right;"> 0.1333579 </td>
   <td style="text-align:right;"> 0.0076443 </td>
   <td style="text-align:right;"> 0.0180989 </td>
   <td style="text-align:right;"> 0.0210322 </td>
   <td style="text-align:right;"> 0.1505057 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.0337976 </td>
   <td style="text-align:right;"> 0.1683476 </td>
   <td style="text-align:right;"> 0.4064472 </td>
   <td style="text-align:right;"> 0.4126249 </td>
   <td style="text-align:right;"> 0.1683476 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.0111633 </td>
   <td style="text-align:right;"> 0.8335098 </td>
   <td style="text-align:right;"> 0.5719938 </td>
   <td style="text-align:right;"> 0.8229356 </td>
   <td style="text-align:right;"> 0.9076095 </td>
   <td style="text-align:right;"> 0.7026117 </td>
   <td style="text-align:right;"> 0.3930946 </td>
   <td style="text-align:right;"> 0.0728892 </td>
   <td style="text-align:right;"> 0.4795858 </td>
   <td style="text-align:right;"> 0.1600978 </td>
   <td style="text-align:right;"> 0.0927290 </td>
   <td style="text-align:right;"> 0.0470797 </td>
   <td style="text-align:right;"> 0.0027043 </td>
   <td style="text-align:right;"> 0.0596070 </td>
   <td style="text-align:right;"> 0.1056930 </td>
   <td style="text-align:right;"> 0.0621267 </td>
   <td style="text-align:right;"> 0.0794289 </td>
   <td style="text-align:right;"> 0.0310178 </td>
   <td style="text-align:right;"> 0.1971386 </td>
   <td style="text-align:right;"> 0.0770697 </td>
   <td style="text-align:right;"> 0.0086776 </td>
   <td style="text-align:right;"> 0.0195598 </td>
   <td style="text-align:right;"> 0.1167438 </td>
   <td style="text-align:right;"> 0.9036474 </td>
   <td style="text-align:right;"> 0.0549757 </td>
   <td style="text-align:right;"> 0.0028435 </td>
   <td style="text-align:right;"> 0.0149371 </td>
   <td style="text-align:right;"> 0.1373461 </td>
   <td style="text-align:right;"> 0.0056150 </td>
   <td style="text-align:right;"> 0.1125680 </td>
   <td style="text-align:right;"> 0.0080623 </td>
   <td style="text-align:right;"> 0.0063548 </td>
   <td style="text-align:right;"> 0.0207254 </td>
   <td style="text-align:right;"> 0.2234544 </td>
   <td style="text-align:right;"> 0.0268852 </td>
   <td style="text-align:right;"> 0.0996080 </td>
   <td style="text-align:right;"> 0.0705490 </td>
   <td style="text-align:right;"> 0.2570216 </td>
   <td style="text-align:right;"> 0.1470778 </td>
   <td style="text-align:right;"> 0.0048502 </td>
   <td style="text-align:right;"> 0.1790858 </td>
   <td style="text-align:right;"> 0.0014548 </td>
   <td style="text-align:right;"> 0.1621513 </td>
   <td style="text-align:right;"> 0.0050182 </td>
   <td style="text-align:right;"> 0.0853011 </td>
   <td style="text-align:right;"> 0.1702644 </td>
   <td style="text-align:right;"> 0.0388773 </td>
   <td style="text-align:right;"> 0.0008584 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.1391036 </td>
   <td style="text-align:right;"> 0.2065979 </td>
   <td style="text-align:right;"> 0.3811707 </td>
   <td style="text-align:right;"> 0.1290883 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.0382633 </td>
   <td style="text-align:right;"> 0.6347952 </td>
   <td style="text-align:right;"> 0.4250435 </td>
   <td style="text-align:right;"> 0.6615392 </td>
   <td style="text-align:right;"> 0.8351745 </td>
   <td style="text-align:right;"> 0.8142314 </td>
   <td style="text-align:right;"> 0.6403618 </td>
   <td style="text-align:right;"> 0.5194040 </td>
   <td style="text-align:right;"> 0.1637730 </td>
   <td style="text-align:right;"> 0.3522408 </td>
   <td style="text-align:right;"> 0.0760275 </td>
   <td style="text-align:right;"> 0.5738482 </td>
   <td style="text-align:right;"> 0.2092582 </td>
   <td style="text-align:right;"> 0.8413087 </td>
   <td style="text-align:right;"> 0.8797371 </td>
   <td style="text-align:right;"> 0.0308587 </td>
   <td style="text-align:right;"> 0.0004037 </td>
   <td style="text-align:right;"> 0.0000206 </td>
   <td style="text-align:right;"> 0.0703512 </td>
   <td style="text-align:right;"> 0.0949319 </td>
   <td style="text-align:right;"> 0.0798987 </td>
   <td style="text-align:right;"> 0.0334629 </td>
   <td style="text-align:right;"> 0.0077119 </td>
   <td style="text-align:right;"> 0.4169298 </td>
   <td style="text-align:right;"> 0.1821457 </td>
   <td style="text-align:right;"> 0.0016430 </td>
   <td style="text-align:right;"> 0.1617653 </td>
   <td style="text-align:right;"> 0.2767678 </td>
   <td style="text-align:right;"> 0.0136636 </td>
   <td style="text-align:right;"> 0.9478809 </td>
   <td style="text-align:right;"> 0.0474201 </td>
   <td style="text-align:right;"> 0.0270672 </td>
   <td style="text-align:right;"> 0.1407224 </td>
   <td style="text-align:right;"> 0.1161856 </td>
   <td style="text-align:right;"> 0.0040532 </td>
   <td style="text-align:right;"> 0.7129535 </td>
   <td style="text-align:right;"> 0.0000103 </td>
   <td style="text-align:right;"> 0.1468724 </td>
   <td style="text-align:right;"> 0.0509872 </td>
   <td style="text-align:right;"> 0.0934868 </td>
   <td style="text-align:right;"> 0.8299436 </td>
   <td style="text-align:right;"> 0.0000259 </td>
   <td style="text-align:right;"> 0.4565195 </td>
   <td style="text-align:right;"> 0.1610039 </td>
   <td style="text-align:right;"> 0.1144425 </td>
   <td style="text-align:right;"> 0.8126443 </td>
   <td style="text-align:right;"> 0.0392128 </td>
   <td style="text-align:right;"> 0.0266730 </td>
   <td style="text-align:right;"> 0.0175003 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.2535174 </td>
   <td style="text-align:right;"> 0.1210124 </td>
   <td style="text-align:right;"> 0.0057176 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 0.0489139 </td>
   <td style="text-align:right;"> 0.0099016 </td>
   <td style="text-align:right;"> 0.2319851 </td>
   <td style="text-align:right;"> 0.2001810 </td>
   <td style="text-align:right;"> 0.4233516 </td>
   <td style="text-align:right;"> 0.5198264 </td>
   <td style="text-align:right;"> 0.4772562 </td>
   <td style="text-align:right;"> 0.5601898 </td>
   <td style="text-align:right;"> 0.0021193 </td>
   <td style="text-align:right;"> 0.1389621 </td>
   <td style="text-align:right;"> 0.9372802 </td>
   <td style="text-align:right;"> 0.5254362 </td>
   <td style="text-align:right;"> 0.1231914 </td>
   <td style="text-align:right;"> 0.9878754 </td>
   <td style="text-align:right;"> 0.6460740 </td>
   <td style="text-align:right;"> 0.1776945 </td>
   <td style="text-align:right;"> 0.0124987 </td>
   <td style="text-align:right;"> 0.0139386 </td>
   <td style="text-align:right;"> 0.8516072 </td>
   <td style="text-align:right;"> 0.1430912 </td>
   <td style="text-align:right;"> 0.0871223 </td>
   <td style="text-align:right;"> 0.1184871 </td>
   <td style="text-align:right;"> 0.5824287 </td>
   <td style="text-align:right;"> 0.8044073 </td>
   <td style="text-align:right;"> 0.0382771 </td>
   <td style="text-align:right;"> 0.1608633 </td>
   <td style="text-align:right;"> 0.2475659 </td>
   <td style="text-align:right;"> 0.5980339 </td>
   <td style="text-align:right;"> 0.0493070 </td>
   <td style="text-align:right;"> 0.4183730 </td>
   <td style="text-align:right;"> 0.0127443 </td>
   <td style="text-align:right;"> 0.0012772 </td>
   <td style="text-align:right;"> 0.3274297 </td>
   <td style="text-align:right;"> 0.1207375 </td>
   <td style="text-align:right;"> 0.0873162 </td>
   <td style="text-align:right;"> 0.7344524 </td>
   <td style="text-align:right;"> 0.0390706 </td>
   <td style="text-align:right;"> 0.5556465 </td>
   <td style="text-align:right;"> 0.1735598 </td>
   <td style="text-align:right;"> 0.2525839 </td>
   <td style="text-align:right;"> 0.1706502 </td>
   <td style="text-align:right;"> 0.0657043 </td>
   <td style="text-align:right;"> 0.0576438 </td>
   <td style="text-align:right;"> 0.0364811 </td>
   <td style="text-align:right;"> 0.0482140 </td>
   <td style="text-align:right;"> 0.8141171 </td>
   <td style="text-align:right;"> 0.0149030 </td>
   <td style="text-align:right;"> 0.1713685 </td>
   <td style="text-align:right;"> 0.0439283 </td>
   <td style="text-align:right;"> 0.0671396 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.4620460 </td>
   <td style="text-align:right;"> 0.1471424 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.4483558 </td>
   <td style="text-align:right;"> 0.3041679 </td>
   <td style="text-align:right;"> 0.6110633 </td>
   <td style="text-align:right;"> 0.9683173 </td>
   <td style="text-align:right;"> 0.3629915 </td>
   <td style="text-align:right;"> 0.3516921 </td>
   <td style="text-align:right;"> 0.0256432 </td>
   <td style="text-align:right;"> 0.5066208 </td>
   <td style="text-align:right;"> 0.3138532 </td>
   <td style="text-align:right;"> 0.0099562 </td>
   <td style="text-align:right;"> 0.5026706 </td>
   <td style="text-align:right;"> 0.5226410 </td>
   <td style="text-align:right;"> 0.6596317 </td>
   <td style="text-align:right;"> 0.5961896 </td>
   <td style="text-align:right;"> 0.5194035 </td>
   <td style="text-align:right;"> 0.5445843 </td>
   <td style="text-align:right;"> 0.0706613 </td>
   <td style="text-align:right;"> 0.2053391 </td>
   <td style="text-align:right;"> 0.0029275 </td>
   <td style="text-align:right;"> 0.3468416 </td>
   <td style="text-align:right;"> 0.3810692 </td>
   <td style="text-align:right;"> 0.0116612 </td>
   <td style="text-align:right;"> 0.0283783 </td>
   <td style="text-align:right;"> 0.7592085 </td>
   <td style="text-align:right;"> 0.2449115 </td>
   <td style="text-align:right;"> 0.2866189 </td>
   <td style="text-align:right;"> 0.6939879 </td>
   <td style="text-align:right;"> 0.7952283 </td>
   <td style="text-align:right;"> 0.0197736 </td>
   <td style="text-align:right;"> 0.4033468 </td>
   <td style="text-align:right;"> 0.4130923 </td>
   <td style="text-align:right;"> 0.1671469 </td>
   <td style="text-align:right;"> 0.9557300 </td>
   <td style="text-align:right;"> 0.1645077 </td>
   <td style="text-align:right;"> 0.1200413 </td>
   <td style="text-align:right;"> 0.6393775 </td>
   <td style="text-align:right;"> 0.0287088 </td>
   <td style="text-align:right;"> 0.1186155 </td>
   <td style="text-align:right;"> 0.4462712 </td>
   <td style="text-align:right;"> 0.6124570 </td>
   <td style="text-align:right;"> 0.9695890 </td>
   <td style="text-align:right;"> 0.1504693 </td>
   <td style="text-align:right;"> 0.0791183 </td>
   <td style="text-align:right;"> 0.5879800 </td>
   <td style="text-align:right;"> 0.0367272 </td>
   <td style="text-align:right;"> 0.6481309 </td>
   <td style="text-align:right;"> 0.0100052 </td>
   <td style="text-align:right;"> 0.1763687 </td>
   <td style="text-align:right;"> 0.1496468 </td>
   <td style="text-align:right;"> 0.0121188 </td>
   <td style="text-align:right;"> 0.2162697 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 0.2336319 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.0075112 </td>
   <td style="text-align:right;"> 0.3550580 </td>
   <td style="text-align:right;"> 0.3650248 </td>
   <td style="text-align:right;"> 0.2139143 </td>
   <td style="text-align:right;"> 0.8962654 </td>
   <td style="text-align:right;"> 0.9269097 </td>
   <td style="text-align:right;"> 0.8586984 </td>
   <td style="text-align:right;"> 0.8754601 </td>
   <td style="text-align:right;"> 0.0451711 </td>
   <td style="text-align:right;"> 0.5566780 </td>
   <td style="text-align:right;"> 0.2930132 </td>
   <td style="text-align:right;"> 0.2681028 </td>
   <td style="text-align:right;"> 0.2877986 </td>
   <td style="text-align:right;"> 0.7267904 </td>
   <td style="text-align:right;"> 0.5694431 </td>
   <td style="text-align:right;"> 0.0024575 </td>
   <td style="text-align:right;"> 0.0000647 </td>
   <td style="text-align:right;"> 0.0000088 </td>
   <td style="text-align:right;"> 0.0958931 </td>
   <td style="text-align:right;"> 0.1076118 </td>
   <td style="text-align:right;"> 0.1029754 </td>
   <td style="text-align:right;"> 0.0081189 </td>
   <td style="text-align:right;"> 0.0098831 </td>
   <td style="text-align:right;"> 0.0672805 </td>
   <td style="text-align:right;"> 0.2625799 </td>
   <td style="text-align:right;"> 0.0102066 </td>
   <td style="text-align:right;"> 0.0320454 </td>
   <td style="text-align:right;"> 0.0810760 </td>
   <td style="text-align:right;"> 0.0070951 </td>
   <td style="text-align:right;"> 0.8197488 </td>
   <td style="text-align:right;"> 0.0433763 </td>
   <td style="text-align:right;"> 0.0561303 </td>
   <td style="text-align:right;"> 0.0741456 </td>
   <td style="text-align:right;"> 0.0056932 </td>
   <td style="text-align:right;"> 0.0027897 </td>
   <td style="text-align:right;"> 0.8061485 </td>
   <td style="text-align:right;"> 0.0000014 </td>
   <td style="text-align:right;"> 0.0574636 </td>
   <td style="text-align:right;"> 0.0130237 </td>
   <td style="text-align:right;"> 0.1686313 </td>
   <td style="text-align:right;"> 0.8864334 </td>
   <td style="text-align:right;"> 0.0001626 </td>
   <td style="text-align:right;"> 0.3344742 </td>
   <td style="text-align:right;"> 0.0531278 </td>
   <td style="text-align:right;"> 0.0175646 </td>
   <td style="text-align:right;"> 0.8248759 </td>
   <td style="text-align:right;"> 0.0033824 </td>
   <td style="text-align:right;"> 0.0267755 </td>
   <td style="text-align:right;"> 0.0154569 </td>
   <td style="text-align:right;"> 0.0000581 </td>
   <td style="text-align:right;"> 0.0208204 </td>
   <td style="text-align:right;"> 0.0583405 </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
</tbody>
</table></div>
<br>  
<br>  
__Pearson correlation adjusted P-value matrix__

```r
kable(lFC.rval.BH_pears) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Sex </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Age_blood </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ischemic_Type </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> tPa </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Thrombectomy </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> TICI2B </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Cholesterol </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Triglycerides </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> HDLC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LDLC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Statin </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> NIHSS_baseline </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> NIHSS_3mo </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> mRS </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CEC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SERPINA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ALB </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AMBP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ANTXR2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APMAP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LPA </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOA4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOB </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOC3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOC4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOD </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOE </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOF </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOL1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> APOM </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CAMP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CLU </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> C3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PPBP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AHSG </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> FGA </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> HPR </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> IHH </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ITIH4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> LCAT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SELL </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PCYOX1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> GPLD1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PF4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PLTP </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PON1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PON3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> RBP4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SAA1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SAA2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> TTR </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Sex </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.1249087 </td>
   <td style="text-align:right;"> -0.3563483 </td>
   <td style="text-align:right;"> -0.0690066 </td>
   <td style="text-align:right;"> -0.0690066 </td>
   <td style="text-align:right;"> -0.0199205 </td>
   <td style="text-align:right;"> 0.2201304 </td>
   <td style="text-align:right;"> 0.2709475 </td>
   <td style="text-align:right;"> -0.2663847 </td>
   <td style="text-align:right;"> 0.3326801 </td>
   <td style="text-align:right;"> 0.3563483 </td>
   <td style="text-align:right;"> -0.7066620 </td>
   <td style="text-align:right;"> -0.4790877 </td>
   <td style="text-align:right;"> -0.2241327 </td>
   <td style="text-align:right;"> 0.5278085 </td>
   <td style="text-align:right;"> -0.4786708 </td>
   <td style="text-align:right;"> -0.7524823 </td>
   <td style="text-align:right;"> -0.6881281 </td>
   <td style="text-align:right;"> -0.1918124 </td>
   <td style="text-align:right;"> -0.5732583 </td>
   <td style="text-align:right;"> -0.2815303 </td>
   <td style="text-align:right;"> -0.3909040 </td>
   <td style="text-align:right;"> -0.2448976 </td>
   <td style="text-align:right;"> -0.4082899 </td>
   <td style="text-align:right;"> -0.1873223 </td>
   <td style="text-align:right;"> -0.4092369 </td>
   <td style="text-align:right;"> -0.5374950 </td>
   <td style="text-align:right;"> -0.4572435 </td>
   <td style="text-align:right;"> -0.4776855 </td>
   <td style="text-align:right;"> -0.0913140 </td>
   <td style="text-align:right;"> -0.5470097 </td>
   <td style="text-align:right;"> -0.5895328 </td>
   <td style="text-align:right;"> -0.3670954 </td>
   <td style="text-align:right;"> -0.3271235 </td>
   <td style="text-align:right;"> -0.6289787 </td>
   <td style="text-align:right;"> -0.4727843 </td>
   <td style="text-align:right;"> -0.6846201 </td>
   <td style="text-align:right;"> -0.5248808 </td>
   <td style="text-align:right;"> -0.5440953 </td>
   <td style="text-align:right;"> -0.3904141 </td>
   <td style="text-align:right;"> 0.0708207 </td>
   <td style="text-align:right;"> -0.6964293 </td>
   <td style="text-align:right;"> -0.5499636 </td>
   <td style="text-align:right;"> -0.7051383 </td>
   <td style="text-align:right;"> -0.5467648 </td>
   <td style="text-align:right;"> -0.4553864 </td>
   <td style="text-align:right;"> -0.5117469 </td>
   <td style="text-align:right;"> -0.5739959 </td>
   <td style="text-align:right;"> -0.6761977 </td>
   <td style="text-align:right;"> -0.5786772 </td>
   <td style="text-align:right;"> 0.5551287 </td>
   <td style="text-align:right;"> 0.2306563 </td>
   <td style="text-align:right;"> -0.7017066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Age_blood </td>
   <td style="text-align:right;"> -0.1249087 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.0847963 </td>
   <td style="text-align:right;"> 0.6675219 </td>
   <td style="text-align:right;"> 0.0290227 </td>
   <td style="text-align:right;"> 0.1684445 </td>
   <td style="text-align:right;"> 0.0071658 </td>
   <td style="text-align:right;"> -0.0981643 </td>
   <td style="text-align:right;"> 0.5691887 </td>
   <td style="text-align:right;"> -0.1602412 </td>
   <td style="text-align:right;"> 0.3056611 </td>
   <td style="text-align:right;"> -0.2519262 </td>
   <td style="text-align:right;"> -0.0086494 </td>
   <td style="text-align:right;"> -0.2877577 </td>
   <td style="text-align:right;"> 0.2426529 </td>
   <td style="text-align:right;"> 0.0791748 </td>
   <td style="text-align:right;"> 0.4523301 </td>
   <td style="text-align:right;"> 0.2244215 </td>
   <td style="text-align:right;"> -0.0907414 </td>
   <td style="text-align:right;"> 0.0267910 </td>
   <td style="text-align:right;"> 0.1821326 </td>
   <td style="text-align:right;"> 0.4705084 </td>
   <td style="text-align:right;"> 0.2405681 </td>
   <td style="text-align:right;"> -0.0403214 </td>
   <td style="text-align:right;"> -0.0071323 </td>
   <td style="text-align:right;"> -0.1489583 </td>
   <td style="text-align:right;"> 0.0168897 </td>
   <td style="text-align:right;"> -0.2822302 </td>
   <td style="text-align:right;"> 0.1650874 </td>
   <td style="text-align:right;"> -0.1311561 </td>
   <td style="text-align:right;"> 0.1385272 </td>
   <td style="text-align:right;"> 0.1802071 </td>
   <td style="text-align:right;"> 0.1856233 </td>
   <td style="text-align:right;"> 0.2472663 </td>
   <td style="text-align:right;"> -0.0042868 </td>
   <td style="text-align:right;"> -0.1563835 </td>
   <td style="text-align:right;"> 0.2829309 </td>
   <td style="text-align:right;"> -0.0389185 </td>
   <td style="text-align:right;"> 0.1520278 </td>
   <td style="text-align:right;"> -0.2682738 </td>
   <td style="text-align:right;"> -0.0293808 </td>
   <td style="text-align:right;"> -0.0006991 </td>
   <td style="text-align:right;"> 0.2488211 </td>
   <td style="text-align:right;"> 0.0228269 </td>
   <td style="text-align:right;"> 0.3595014 </td>
   <td style="text-align:right;"> -0.1659305 </td>
   <td style="text-align:right;"> 0.3662092 </td>
   <td style="text-align:right;"> 0.0393083 </td>
   <td style="text-align:right;"> 0.0647647 </td>
   <td style="text-align:right;"> 0.1457133 </td>
   <td style="text-align:right;"> -0.6841765 </td>
   <td style="text-align:right;"> -0.3090736 </td>
   <td style="text-align:right;"> 0.2794976 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ischemic_Type </td>
   <td style="text-align:right;"> -0.3563483 </td>
   <td style="text-align:right;"> 0.0847963 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.0430331 </td>
   <td style="text-align:right;"> -0.0430331 </td>
   <td style="text-align:right;"> -0.1490712 </td>
   <td style="text-align:right;"> -0.4360514 </td>
   <td style="text-align:right;"> 0.0288072 </td>
   <td style="text-align:right;"> -0.0345530 </td>
   <td style="text-align:right;"> -0.3812499 </td>
   <td style="text-align:right;"> 0.2222222 </td>
   <td style="text-align:right;"> 0.4638749 </td>
   <td style="text-align:right;"> 0.0275456 </td>
   <td style="text-align:right;"> -0.1397713 </td>
   <td style="text-align:right;"> -0.0891744 </td>
   <td style="text-align:right;"> -0.0915318 </td>
   <td style="text-align:right;"> 0.2390248 </td>
   <td style="text-align:right;"> 0.2691527 </td>
   <td style="text-align:right;"> 0.0100892 </td>
   <td style="text-align:right;"> 0.2390585 </td>
   <td style="text-align:right;"> -0.4006884 </td>
   <td style="text-align:right;"> 0.2603396 </td>
   <td style="text-align:right;"> 0.0612116 </td>
   <td style="text-align:right;"> 0.3484210 </td>
   <td style="text-align:right;"> 0.0324905 </td>
   <td style="text-align:right;"> 0.3275359 </td>
   <td style="text-align:right;"> 0.0870459 </td>
   <td style="text-align:right;"> -0.2494951 </td>
   <td style="text-align:right;"> 0.1219456 </td>
   <td style="text-align:right;"> -0.0282349 </td>
   <td style="text-align:right;"> 0.1630347 </td>
   <td style="text-align:right;"> 0.5183584 </td>
   <td style="text-align:right;"> 0.0294007 </td>
   <td style="text-align:right;"> 0.2122196 </td>
   <td style="text-align:right;"> 0.0347987 </td>
   <td style="text-align:right;"> -0.0816276 </td>
   <td style="text-align:right;"> 0.2850472 </td>
   <td style="text-align:right;"> 0.1547731 </td>
   <td style="text-align:right;"> -0.1526560 </td>
   <td style="text-align:right;"> 0.3214070 </td>
   <td style="text-align:right;"> 0.3526735 </td>
   <td style="text-align:right;"> 0.3792686 </td>
   <td style="text-align:right;"> 0.3632190 </td>
   <td style="text-align:right;"> 0.1998870 </td>
   <td style="text-align:right;"> 0.0663445 </td>
   <td style="text-align:right;"> -0.0089270 </td>
   <td style="text-align:right;"> 0.0896694 </td>
   <td style="text-align:right;"> 0.0395150 </td>
   <td style="text-align:right;"> 0.1729750 </td>
   <td style="text-align:right;"> 0.2423324 </td>
   <td style="text-align:right;"> -0.3563901 </td>
   <td style="text-align:right;"> -0.1558881 </td>
   <td style="text-align:right;"> 0.2739832 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> tPa </td>
   <td style="text-align:right;"> -0.0690066 </td>
   <td style="text-align:right;"> 0.6675219 </td>
   <td style="text-align:right;"> -0.0430331 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2666667 </td>
   <td style="text-align:right;"> 0.5003702 </td>
   <td style="text-align:right;"> 0.5504302 </td>
   <td style="text-align:right;"> 0.3359970 </td>
   <td style="text-align:right;"> 0.5311746 </td>
   <td style="text-align:right;"> 0.4485102 </td>
   <td style="text-align:right;"> 0.0430331 </td>
   <td style="text-align:right;"> -0.1796580 </td>
   <td style="text-align:right;"> -0.3873435 </td>
   <td style="text-align:right;"> -0.4736655 </td>
   <td style="text-align:right;"> 0.2678994 </td>
   <td style="text-align:right;"> 0.1282102 </td>
   <td style="text-align:right;"> 0.3694524 </td>
   <td style="text-align:right;"> 0.2851676 </td>
   <td style="text-align:right;"> 0.0466718 </td>
   <td style="text-align:right;"> -0.2144546 </td>
   <td style="text-align:right;"> 0.1617256 </td>
   <td style="text-align:right;"> 0.3569149 </td>
   <td style="text-align:right;"> 0.2531168 </td>
   <td style="text-align:right;"> 0.3045794 </td>
   <td style="text-align:right;"> -0.1924477 </td>
   <td style="text-align:right;"> -0.1009288 </td>
   <td style="text-align:right;"> 0.0893930 </td>
   <td style="text-align:right;"> 0.1204422 </td>
   <td style="text-align:right;"> -0.1546067 </td>
   <td style="text-align:right;"> -0.4478025 </td>
   <td style="text-align:right;"> -0.2299454 </td>
   <td style="text-align:right;"> -0.1961637 </td>
   <td style="text-align:right;"> 0.2603343 </td>
   <td style="text-align:right;"> 0.2876959 </td>
   <td style="text-align:right;"> 0.0000615 </td>
   <td style="text-align:right;"> -0.4753756 </td>
   <td style="text-align:right;"> 0.3225066 </td>
   <td style="text-align:right;"> -0.2537164 </td>
   <td style="text-align:right;"> 0.5055409 </td>
   <td style="text-align:right;"> -0.2421058 </td>
   <td style="text-align:right;"> -0.2736705 </td>
   <td style="text-align:right;"> -0.0052876 </td>
   <td style="text-align:right;"> -0.3282950 </td>
   <td style="text-align:right;"> -0.2596329 </td>
   <td style="text-align:right;"> -0.0025225 </td>
   <td style="text-align:right;"> -0.4294680 </td>
   <td style="text-align:right;"> 0.2803789 </td>
   <td style="text-align:right;"> -0.2168310 </td>
   <td style="text-align:right;"> -0.0689355 </td>
   <td style="text-align:right;"> 0.1344105 </td>
   <td style="text-align:right;"> -0.3800747 </td>
   <td style="text-align:right;"> 0.0122501 </td>
   <td style="text-align:right;"> 0.3695822 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Thrombectomy </td>
   <td style="text-align:right;"> -0.0690066 </td>
   <td style="text-align:right;"> 0.0290227 </td>
   <td style="text-align:right;"> -0.0430331 </td>
   <td style="text-align:right;"> 0.2666667 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.9237604 </td>
   <td style="text-align:right;"> 0.2752151 </td>
   <td style="text-align:right;"> 0.2179904 </td>
   <td style="text-align:right;"> 0.4405867 </td>
   <td style="text-align:right;"> 0.1923222 </td>
   <td style="text-align:right;"> 0.0430331 </td>
   <td style="text-align:right;"> 0.2155896 </td>
   <td style="text-align:right;"> 0.1903892 </td>
   <td style="text-align:right;"> -0.1759329 </td>
   <td style="text-align:right;"> 0.2719383 </td>
   <td style="text-align:right;"> 0.0338402 </td>
   <td style="text-align:right;"> 0.1818629 </td>
   <td style="text-align:right;"> 0.1759026 </td>
   <td style="text-align:right;"> -0.3222404 </td>
   <td style="text-align:right;"> 0.4662228 </td>
   <td style="text-align:right;"> 0.4413547 </td>
   <td style="text-align:right;"> -0.2692017 </td>
   <td style="text-align:right;"> -0.2078159 </td>
   <td style="text-align:right;"> -0.3685675 </td>
   <td style="text-align:right;"> 0.3282038 </td>
   <td style="text-align:right;"> 0.1444066 </td>
   <td style="text-align:right;"> 0.3551515 </td>
   <td style="text-align:right;"> 0.1252168 </td>
   <td style="text-align:right;"> -0.0607473 </td>
   <td style="text-align:right;"> 0.2372384 </td>
   <td style="text-align:right;"> 0.1097366 </td>
   <td style="text-align:right;"> 0.1879125 </td>
   <td style="text-align:right;"> -0.0487667 </td>
   <td style="text-align:right;"> 0.0765954 </td>
   <td style="text-align:right;"> 0.4021136 </td>
   <td style="text-align:right;"> 0.2467022 </td>
   <td style="text-align:right;"> 0.1504787 </td>
   <td style="text-align:right;"> -0.1482612 </td>
   <td style="text-align:right;"> 0.5784861 </td>
   <td style="text-align:right;"> -0.1251021 </td>
   <td style="text-align:right;"> -0.0853650 </td>
   <td style="text-align:right;"> 0.0742457 </td>
   <td style="text-align:right;"> -0.2055219 </td>
   <td style="text-align:right;"> -0.1385398 </td>
   <td style="text-align:right;"> -0.2857829 </td>
   <td style="text-align:right;"> 0.2944537 </td>
   <td style="text-align:right;"> -0.0612862 </td>
   <td style="text-align:right;"> 0.0517411 </td>
   <td style="text-align:right;"> 0.0357834 </td>
   <td style="text-align:right;"> 0.0641090 </td>
   <td style="text-align:right;"> -0.2431921 </td>
   <td style="text-align:right;"> 0.2751016 </td>
   <td style="text-align:right;"> 0.0401972 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TICI2B </td>
   <td style="text-align:right;"> -0.0199205 </td>
   <td style="text-align:right;"> 0.1684445 </td>
   <td style="text-align:right;"> -0.1490712 </td>
   <td style="text-align:right;"> 0.5003702 </td>
   <td style="text-align:right;"> 0.9237604 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3502924 </td>
   <td style="text-align:right;"> 0.4226609 </td>
   <td style="text-align:right;"> 0.5194436 </td>
   <td style="text-align:right;"> 0.2902029 </td>
   <td style="text-align:right;"> 0.1490712 </td>
   <td style="text-align:right;"> 0.1037256 </td>
   <td style="text-align:right;"> -0.0596986 </td>
   <td style="text-align:right;"> -0.3516054 </td>
   <td style="text-align:right;"> 0.2159655 </td>
   <td style="text-align:right;"> -0.0084306 </td>
   <td style="text-align:right;"> 0.1611085 </td>
   <td style="text-align:right;"> 0.0617825 </td>
   <td style="text-align:right;"> -0.2385781 </td>
   <td style="text-align:right;"> 0.3068168 </td>
   <td style="text-align:right;"> 0.3426146 </td>
   <td style="text-align:right;"> -0.2154853 </td>
   <td style="text-align:right;"> -0.2158945 </td>
   <td style="text-align:right;"> -0.2519418 </td>
   <td style="text-align:right;"> 0.1714115 </td>
   <td style="text-align:right;"> -0.1120529 </td>
   <td style="text-align:right;"> 0.3075865 </td>
   <td style="text-align:right;"> 0.2033120 </td>
   <td style="text-align:right;"> -0.1603868 </td>
   <td style="text-align:right;"> -0.0277883 </td>
   <td style="text-align:right;"> -0.1066563 </td>
   <td style="text-align:right;"> -0.0486015 </td>
   <td style="text-align:right;"> -0.1467302 </td>
   <td style="text-align:right;"> 0.1658216 </td>
   <td style="text-align:right;"> 0.3135620 </td>
   <td style="text-align:right;"> 0.1262534 </td>
   <td style="text-align:right;"> 0.1216399 </td>
   <td style="text-align:right;"> -0.1454950 </td>
   <td style="text-align:right;"> 0.6538657 </td>
   <td style="text-align:right;"> -0.3514561 </td>
   <td style="text-align:right;"> -0.2874374 </td>
   <td style="text-align:right;"> -0.1252244 </td>
   <td style="text-align:right;"> -0.2364365 </td>
   <td style="text-align:right;"> -0.2653618 </td>
   <td style="text-align:right;"> -0.2244579 </td>
   <td style="text-align:right;"> 0.2081899 </td>
   <td style="text-align:right;"> 0.0426049 </td>
   <td style="text-align:right;"> -0.1131923 </td>
   <td style="text-align:right;"> -0.1173486 </td>
   <td style="text-align:right;"> -0.0723769 </td>
   <td style="text-align:right;"> -0.1965567 </td>
   <td style="text-align:right;"> 0.2813786 </td>
   <td style="text-align:right;"> 0.0282879 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cholesterol </td>
   <td style="text-align:right;"> 0.2201304 </td>
   <td style="text-align:right;"> 0.0071658 </td>
   <td style="text-align:right;"> -0.4360514 </td>
   <td style="text-align:right;"> 0.5504302 </td>
   <td style="text-align:right;"> 0.2752151 </td>
   <td style="text-align:right;"> 0.3502924 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.1976010 </td>
   <td style="text-align:right;"> -0.0061813 </td>
   <td style="text-align:right;"> 0.9593977 </td>
   <td style="text-align:right;"> -0.1938006 </td>
   <td style="text-align:right;"> -0.1526244 </td>
   <td style="text-align:right;"> -0.2710240 </td>
   <td style="text-align:right;"> -0.0279343 </td>
   <td style="text-align:right;"> 0.3942284 </td>
   <td style="text-align:right;"> 0.1313135 </td>
   <td style="text-align:right;"> -0.0916029 </td>
   <td style="text-align:right;"> 0.0658435 </td>
   <td style="text-align:right;"> -0.3586046 </td>
   <td style="text-align:right;"> -0.4066107 </td>
   <td style="text-align:right;"> 0.1260855 </td>
   <td style="text-align:right;"> -0.2125722 </td>
   <td style="text-align:right;"> 0.0143034 </td>
   <td style="text-align:right;"> 0.2282246 </td>
   <td style="text-align:right;"> -0.3809930 </td>
   <td style="text-align:right;"> -0.0534254 </td>
   <td style="text-align:right;"> -0.0446908 </td>
   <td style="text-align:right;"> 0.0914941 </td>
   <td style="text-align:right;"> -0.4812312 </td>
   <td style="text-align:right;"> -0.0734986 </td>
   <td style="text-align:right;"> -0.2771120 </td>
   <td style="text-align:right;"> -0.5105520 </td>
   <td style="text-align:right;"> 0.3955915 </td>
   <td style="text-align:right;"> -0.2052591 </td>
   <td style="text-align:right;"> -0.2419448 </td>
   <td style="text-align:right;"> -0.4570208 </td>
   <td style="text-align:right;"> -0.1275314 </td>
   <td style="text-align:right;"> -0.5135293 </td>
   <td style="text-align:right;"> 0.1963514 </td>
   <td style="text-align:right;"> -0.1255198 </td>
   <td style="text-align:right;"> -0.2468028 </td>
   <td style="text-align:right;"> -0.1176711 </td>
   <td style="text-align:right;"> -0.8482983 </td>
   <td style="text-align:right;"> -0.3368757 </td>
   <td style="text-align:right;"> -0.4947810 </td>
   <td style="text-align:right;"> -0.4874364 </td>
   <td style="text-align:right;"> -0.3927382 </td>
   <td style="text-align:right;"> -0.2725973 </td>
   <td style="text-align:right;"> -0.2588695 </td>
   <td style="text-align:right;"> -0.1433470 </td>
   <td style="text-align:right;"> 0.2165828 </td>
   <td style="text-align:right;"> 0.6138286 </td>
   <td style="text-align:right;"> -0.0548690 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Triglycerides </td>
   <td style="text-align:right;"> 0.2709475 </td>
   <td style="text-align:right;"> -0.0981643 </td>
   <td style="text-align:right;"> 0.0288072 </td>
   <td style="text-align:right;"> 0.3359970 </td>
   <td style="text-align:right;"> 0.2179904 </td>
   <td style="text-align:right;"> 0.4226609 </td>
   <td style="text-align:right;"> 0.1976010 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3248869 </td>
   <td style="text-align:right;"> 0.3393993 </td>
   <td style="text-align:right;"> 0.4891688 </td>
   <td style="text-align:right;"> -0.1133281 </td>
   <td style="text-align:right;"> -0.5286057 </td>
   <td style="text-align:right;"> -0.3632492 </td>
   <td style="text-align:right;"> 0.2252365 </td>
   <td style="text-align:right;"> 0.1260281 </td>
   <td style="text-align:right;"> -0.0708605 </td>
   <td style="text-align:right;"> -0.0757404 </td>
   <td style="text-align:right;"> 0.0260635 </td>
   <td style="text-align:right;"> -0.0678586 </td>
   <td style="text-align:right;"> -0.2759877 </td>
   <td style="text-align:right;"> -0.1900140 </td>
   <td style="text-align:right;"> -0.1055000 </td>
   <td style="text-align:right;"> 0.3859239 </td>
   <td style="text-align:right;"> -0.0668582 </td>
   <td style="text-align:right;"> -0.2296133 </td>
   <td style="text-align:right;"> 0.0742155 </td>
   <td style="text-align:right;"> 0.2375354 </td>
   <td style="text-align:right;"> -0.1155071 </td>
   <td style="text-align:right;"> -0.4100490 </td>
   <td style="text-align:right;"> -0.3023832 </td>
   <td style="text-align:right;"> -0.3606594 </td>
   <td style="text-align:right;"> -0.3761779 </td>
   <td style="text-align:right;"> 0.4991277 </td>
   <td style="text-align:right;"> 0.0797681 </td>
   <td style="text-align:right;"> -0.3528549 </td>
   <td style="text-align:right;"> 0.0548097 </td>
   <td style="text-align:right;"> 0.1588479 </td>
   <td style="text-align:right;"> 0.2121357 </td>
   <td style="text-align:right;"> -0.4267040 </td>
   <td style="text-align:right;"> -0.2104401 </td>
   <td style="text-align:right;"> -0.2513375 </td>
   <td style="text-align:right;"> -0.2206988 </td>
   <td style="text-align:right;"> -0.2791451 </td>
   <td style="text-align:right;"> 0.0023138 </td>
   <td style="text-align:right;"> -0.2150772 </td>
   <td style="text-align:right;"> 0.1091440 </td>
   <td style="text-align:right;"> -0.3603460 </td>
   <td style="text-align:right;"> -0.5131760 </td>
   <td style="text-align:right;"> -0.1967517 </td>
   <td style="text-align:right;"> 0.1782288 </td>
   <td style="text-align:right;"> 0.2026866 </td>
   <td style="text-align:right;"> 0.0483113 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HDLC </td>
   <td style="text-align:right;"> -0.2663847 </td>
   <td style="text-align:right;"> 0.5691887 </td>
   <td style="text-align:right;"> -0.0345530 </td>
   <td style="text-align:right;"> 0.5311746 </td>
   <td style="text-align:right;"> 0.4405867 </td>
   <td style="text-align:right;"> 0.5194436 </td>
   <td style="text-align:right;"> -0.0061813 </td>
   <td style="text-align:right;"> 0.3248869 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.1786405 </td>
   <td style="text-align:right;"> 0.1222643 </td>
   <td style="text-align:right;"> -0.1220613 </td>
   <td style="text-align:right;"> 0.1450651 </td>
   <td style="text-align:right;"> -0.2332098 </td>
   <td style="text-align:right;"> 0.3463839 </td>
   <td style="text-align:right;"> 0.4570297 </td>
   <td style="text-align:right;"> 0.5850148 </td>
   <td style="text-align:right;"> 0.5541040 </td>
   <td style="text-align:right;"> 0.1458865 </td>
   <td style="text-align:right;"> 0.3509271 </td>
   <td style="text-align:right;"> 0.5384227 </td>
   <td style="text-align:right;"> 0.2103815 </td>
   <td style="text-align:right;"> 0.0847033 </td>
   <td style="text-align:right;"> 0.0527266 </td>
   <td style="text-align:right;"> 0.6067756 </td>
   <td style="text-align:right;"> 0.1833791 </td>
   <td style="text-align:right;"> 0.2992517 </td>
   <td style="text-align:right;"> 0.3983270 </td>
   <td style="text-align:right;"> 0.4590767 </td>
   <td style="text-align:right;"> -0.0127906 </td>
   <td style="text-align:right;"> 0.4377405 </td>
   <td style="text-align:right;"> 0.4437190 </td>
   <td style="text-align:right;"> -0.0210200 </td>
   <td style="text-align:right;"> 0.6007881 </td>
   <td style="text-align:right;"> 0.6001283 </td>
   <td style="text-align:right;"> -0.0252173 </td>
   <td style="text-align:right;"> 0.5605500 </td>
   <td style="text-align:right;"> 0.2121094 </td>
   <td style="text-align:right;"> 0.6329001 </td>
   <td style="text-align:right;"> -0.0149174 </td>
   <td style="text-align:right;"> 0.1112550 </td>
   <td style="text-align:right;"> 0.2662364 </td>
   <td style="text-align:right;"> 0.3002933 </td>
   <td style="text-align:right;"> 0.3078119 </td>
   <td style="text-align:right;"> 0.5142097 </td>
   <td style="text-align:right;"> -0.0111774 </td>
   <td style="text-align:right;"> 0.7067664 </td>
   <td style="text-align:right;"> 0.1859939 </td>
   <td style="text-align:right;"> 0.2154662 </td>
   <td style="text-align:right;"> 0.4102892 </td>
   <td style="text-align:right;"> -0.7690763 </td>
   <td style="text-align:right;"> -0.3032435 </td>
   <td style="text-align:right;"> 0.5629398 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LDLC </td>
   <td style="text-align:right;"> 0.3326801 </td>
   <td style="text-align:right;"> -0.1602412 </td>
   <td style="text-align:right;"> -0.3812499 </td>
   <td style="text-align:right;"> 0.4485102 </td>
   <td style="text-align:right;"> 0.1923222 </td>
   <td style="text-align:right;"> 0.2902029 </td>
   <td style="text-align:right;"> 0.9593977 </td>
   <td style="text-align:right;"> 0.3393993 </td>
   <td style="text-align:right;"> -0.1786405 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.0891320 </td>
   <td style="text-align:right;"> -0.1355682 </td>
   <td style="text-align:right;"> -0.4106565 </td>
   <td style="text-align:right;"> -0.0544225 </td>
   <td style="text-align:right;"> 0.3260744 </td>
   <td style="text-align:right;"> 0.0350296 </td>
   <td style="text-align:right;"> -0.2475307 </td>
   <td style="text-align:right;"> -0.0976329 </td>
   <td style="text-align:right;"> -0.3567807 </td>
   <td style="text-align:right;"> -0.4746744 </td>
   <td style="text-align:right;"> -0.0871632 </td>
   <td style="text-align:right;"> -0.2921999 </td>
   <td style="text-align:right;"> -0.0337229 </td>
   <td style="text-align:right;"> 0.2875625 </td>
   <td style="text-align:right;"> -0.5154323 </td>
   <td style="text-align:right;"> -0.1501120 </td>
   <td style="text-align:right;"> -0.0980013 </td>
   <td style="text-align:right;"> 0.0404411 </td>
   <td style="text-align:right;"> -0.5812730 </td>
   <td style="text-align:right;"> -0.1625928 </td>
   <td style="text-align:right;"> -0.4352023 </td>
   <td style="text-align:right;"> -0.6632347 </td>
   <td style="text-align:right;"> 0.2746022 </td>
   <td style="text-align:right;"> -0.2174518 </td>
   <td style="text-align:right;"> -0.3518337 </td>
   <td style="text-align:right;"> -0.4947465 </td>
   <td style="text-align:right;"> -0.2437681 </td>
   <td style="text-align:right;"> -0.4824073 </td>
   <td style="text-align:right;"> 0.0707855 </td>
   <td style="text-align:right;"> -0.2134242 </td>
   <td style="text-align:right;"> -0.3033643 </td>
   <td style="text-align:right;"> -0.2346563 </td>
   <td style="text-align:right;"> -0.9008284 </td>
   <td style="text-align:right;"> -0.4513308 </td>
   <td style="text-align:right;"> -0.5790447 </td>
   <td style="text-align:right;"> -0.4927250 </td>
   <td style="text-align:right;"> -0.5088037 </td>
   <td style="text-align:right;"> -0.3817900 </td>
   <td style="text-align:right;"> -0.4135695 </td>
   <td style="text-align:right;"> -0.2810713 </td>
   <td style="text-align:right;"> 0.4334539 </td>
   <td style="text-align:right;"> 0.6838157 </td>
   <td style="text-align:right;"> -0.1798006 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Statin </td>
   <td style="text-align:right;"> 0.3563483 </td>
   <td style="text-align:right;"> 0.3056611 </td>
   <td style="text-align:right;"> 0.2222222 </td>
   <td style="text-align:right;"> 0.0430331 </td>
   <td style="text-align:right;"> 0.0430331 </td>
   <td style="text-align:right;"> 0.1490712 </td>
   <td style="text-align:right;"> -0.1938006 </td>
   <td style="text-align:right;"> 0.4891688 </td>
   <td style="text-align:right;"> 0.1222643 </td>
   <td style="text-align:right;"> -0.0891320 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.0463875 </td>
   <td style="text-align:right;"> -0.1207768 </td>
   <td style="text-align:right;"> -0.0524142 </td>
   <td style="text-align:right;"> 0.1761176 </td>
   <td style="text-align:right;"> -0.0870389 </td>
   <td style="text-align:right;"> -0.2115733 </td>
   <td style="text-align:right;"> -0.4010337 </td>
   <td style="text-align:right;"> -0.5003344 </td>
   <td style="text-align:right;"> 0.1089970 </td>
   <td style="text-align:right;"> -0.3874548 </td>
   <td style="text-align:right;"> -0.1211334 </td>
   <td style="text-align:right;"> -0.2047578 </td>
   <td style="text-align:right;"> -0.0638842 </td>
   <td style="text-align:right;"> -0.1482268 </td>
   <td style="text-align:right;"> -0.4124044 </td>
   <td style="text-align:right;"> 0.0342380 </td>
   <td style="text-align:right;"> -0.5003482 </td>
   <td style="text-align:right;"> -0.0153457 </td>
   <td style="text-align:right;"> 0.0926246 </td>
   <td style="text-align:right;"> -0.0502045 </td>
   <td style="text-align:right;"> -0.1376247 </td>
   <td style="text-align:right;"> -0.2987831 </td>
   <td style="text-align:right;"> 0.3016631 </td>
   <td style="text-align:right;"> -0.2313635 </td>
   <td style="text-align:right;"> -0.0117934 </td>
   <td style="text-align:right;"> -0.2877258 </td>
   <td style="text-align:right;"> 0.1802598 </td>
   <td style="text-align:right;"> -0.3708440 </td>
   <td style="text-align:right;"> -0.5915026 </td>
   <td style="text-align:right;"> 0.0647371 </td>
   <td style="text-align:right;"> -0.4376760 </td>
   <td style="text-align:right;"> 0.1975705 </td>
   <td style="text-align:right;"> -0.1371942 </td>
   <td style="text-align:right;"> 0.1242043 </td>
   <td style="text-align:right;"> 0.0776946 </td>
   <td style="text-align:right;"> -0.1519739 </td>
   <td style="text-align:right;"> -0.1673150 </td>
   <td style="text-align:right;"> -0.4853268 </td>
   <td style="text-align:right;"> -0.5084465 </td>
   <td style="text-align:right;"> -0.0242666 </td>
   <td style="text-align:right;"> 0.2045341 </td>
   <td style="text-align:right;"> -0.3159206 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NIHSS_baseline </td>
   <td style="text-align:right;"> -0.7066620 </td>
   <td style="text-align:right;"> -0.2519262 </td>
   <td style="text-align:right;"> 0.4638749 </td>
   <td style="text-align:right;"> -0.1796580 </td>
   <td style="text-align:right;"> 0.2155896 </td>
   <td style="text-align:right;"> 0.1037256 </td>
   <td style="text-align:right;"> -0.1526244 </td>
   <td style="text-align:right;"> -0.1133281 </td>
   <td style="text-align:right;"> -0.1220613 </td>
   <td style="text-align:right;"> -0.1355682 </td>
   <td style="text-align:right;"> -0.0463875 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.5502278 </td>
   <td style="text-align:right;"> 0.5105867 </td>
   <td style="text-align:right;"> -0.6396876 </td>
   <td style="text-align:right;"> 0.3480218 </td>
   <td style="text-align:right;"> 0.2875622 </td>
   <td style="text-align:right;"> 0.2952877 </td>
   <td style="text-align:right;"> -0.1371178 </td>
   <td style="text-align:right;"> 0.6009155 </td>
   <td style="text-align:right;"> 0.1169844 </td>
   <td style="text-align:right;"> 0.2043706 </td>
   <td style="text-align:right;"> 0.0886097 </td>
   <td style="text-align:right;"> 0.2871463 </td>
   <td style="text-align:right;"> 0.0643033 </td>
   <td style="text-align:right;"> 0.4202140 </td>
   <td style="text-align:right;"> 0.7153340 </td>
   <td style="text-align:right;"> 0.2547761 </td>
   <td style="text-align:right;"> 0.3144246 </td>
   <td style="text-align:right;"> 0.4148064 </td>
   <td style="text-align:right;"> 0.4232019 </td>
   <td style="text-align:right;"> 0.4115959 </td>
   <td style="text-align:right;"> 0.4122684 </td>
   <td style="text-align:right;"> 0.2847252 </td>
   <td style="text-align:right;"> 0.4049609 </td>
   <td style="text-align:right;"> 0.5805498 </td>
   <td style="text-align:right;"> 0.2601642 </td>
   <td style="text-align:right;"> 0.4318794 </td>
   <td style="text-align:right;"> 0.2683951 </td>
   <td style="text-align:right;"> 0.3571332 </td>
   <td style="text-align:right;"> 0.1083986 </td>
   <td style="text-align:right;"> 0.5148074 </td>
   <td style="text-align:right;"> 0.3108288 </td>
   <td style="text-align:right;"> 0.5611951 </td>
   <td style="text-align:right;"> 0.2051371 </td>
   <td style="text-align:right;"> 0.6076428 </td>
   <td style="text-align:right;"> 0.0874129 </td>
   <td style="text-align:right;"> 0.5512961 </td>
   <td style="text-align:right;"> 0.5589009 </td>
   <td style="text-align:right;"> 0.1721538 </td>
   <td style="text-align:right;"> -0.1939736 </td>
   <td style="text-align:right;"> 0.1952591 </td>
   <td style="text-align:right;"> 0.3317746 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NIHSS_3mo </td>
   <td style="text-align:right;"> -0.4790877 </td>
   <td style="text-align:right;"> -0.0086494 </td>
   <td style="text-align:right;"> 0.0275456 </td>
   <td style="text-align:right;"> -0.3873435 </td>
   <td style="text-align:right;"> 0.1903892 </td>
   <td style="text-align:right;"> -0.0596986 </td>
   <td style="text-align:right;"> -0.2710240 </td>
   <td style="text-align:right;"> -0.5286057 </td>
   <td style="text-align:right;"> 0.1450651 </td>
   <td style="text-align:right;"> -0.4106565 </td>
   <td style="text-align:right;"> -0.1207768 </td>
   <td style="text-align:right;"> 0.5502278 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.7869718 </td>
   <td style="text-align:right;"> -0.2829038 </td>
   <td style="text-align:right;"> 0.5513979 </td>
   <td style="text-align:right;"> 0.3449104 </td>
   <td style="text-align:right;"> 0.3977183 </td>
   <td style="text-align:right;"> -0.1596144 </td>
   <td style="text-align:right;"> 0.7236294 </td>
   <td style="text-align:right;"> 0.6603669 </td>
   <td style="text-align:right;"> 0.2436845 </td>
   <td style="text-align:right;"> 0.2180131 </td>
   <td style="text-align:right;"> -0.2535170 </td>
   <td style="text-align:right;"> 0.5369207 </td>
   <td style="text-align:right;"> 0.6087493 </td>
   <td style="text-align:right;"> 0.6580000 </td>
   <td style="text-align:right;"> 0.1151279 </td>
   <td style="text-align:right;"> 0.7060497 </td>
   <td style="text-align:right;"> 0.8603869 </td>
   <td style="text-align:right;"> 0.8797579 </td>
   <td style="text-align:right;"> 0.6949345 </td>
   <td style="text-align:right;"> 0.5122723 </td>
   <td style="text-align:right;"> 0.2356105 </td>
   <td style="text-align:right;"> 0.5663707 </td>
   <td style="text-align:right;"> 0.7355172 </td>
   <td style="text-align:right;"> 0.2596322 </td>
   <td style="text-align:right;"> 0.4218165 </td>
   <td style="text-align:right;"> 0.1570567 </td>
   <td style="text-align:right;"> 0.4552276 </td>
   <td style="text-align:right;"> 0.4204540 </td>
   <td style="text-align:right;"> 0.5986834 </td>
   <td style="text-align:right;"> 0.4577798 </td>
   <td style="text-align:right;"> 0.7515435 </td>
   <td style="text-align:right;"> 0.4459148 </td>
   <td style="text-align:right;"> 0.6385606 </td>
   <td style="text-align:right;"> 0.2123547 </td>
   <td style="text-align:right;"> 0.8763335 </td>
   <td style="text-align:right;"> 0.7575520 </td>
   <td style="text-align:right;"> 0.3730915 </td>
   <td style="text-align:right;"> -0.4496292 </td>
   <td style="text-align:right;"> -0.1352114 </td>
   <td style="text-align:right;"> 0.3191727 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mRS </td>
   <td style="text-align:right;"> -0.2241327 </td>
   <td style="text-align:right;"> -0.2877577 </td>
   <td style="text-align:right;"> -0.1397713 </td>
   <td style="text-align:right;"> -0.4736655 </td>
   <td style="text-align:right;"> -0.1759329 </td>
   <td style="text-align:right;"> -0.3516054 </td>
   <td style="text-align:right;"> -0.0279343 </td>
   <td style="text-align:right;"> -0.3632492 </td>
   <td style="text-align:right;"> -0.2332098 </td>
   <td style="text-align:right;"> -0.0544225 </td>
   <td style="text-align:right;"> -0.0524142 </td>
   <td style="text-align:right;"> 0.5105867 </td>
   <td style="text-align:right;"> 0.7869718 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.3751604 </td>
   <td style="text-align:right;"> 0.5716619 </td>
   <td style="text-align:right;"> -0.0206494 </td>
   <td style="text-align:right;"> 0.1220707 </td>
   <td style="text-align:right;"> -0.2212478 </td>
   <td style="text-align:right;"> 0.3494943 </td>
   <td style="text-align:right;"> 0.3820224 </td>
   <td style="text-align:right;"> 0.1506733 </td>
   <td style="text-align:right;"> 0.2105947 </td>
   <td style="text-align:right;"> 0.0243512 </td>
   <td style="text-align:right;"> 0.1898888 </td>
   <td style="text-align:right;"> 0.4577474 </td>
   <td style="text-align:right;"> 0.5559495 </td>
   <td style="text-align:right;"> 0.1217595 </td>
   <td style="text-align:right;"> 0.5036646 </td>
   <td style="text-align:right;"> 0.8048582 </td>
   <td style="text-align:right;"> 0.6416481 </td>
   <td style="text-align:right;"> 0.2619831 </td>
   <td style="text-align:right;"> 0.6044156 </td>
   <td style="text-align:right;"> 0.1680644 </td>
   <td style="text-align:right;"> 0.2225495 </td>
   <td style="text-align:right;"> 0.5123814 </td>
   <td style="text-align:right;"> -0.0626878 </td>
   <td style="text-align:right;"> 0.3595042 </td>
   <td style="text-align:right;"> -0.1498550 </td>
   <td style="text-align:right;"> 0.4045777 </td>
   <td style="text-align:right;"> 0.3919728 </td>
   <td style="text-align:right;"> 0.3860824 </td>
   <td style="text-align:right;"> 0.1661860 </td>
   <td style="text-align:right;"> 0.6558529 </td>
   <td style="text-align:right;"> 0.3172330 </td>
   <td style="text-align:right;"> 0.4125514 </td>
   <td style="text-align:right;"> -0.0410795 </td>
   <td style="text-align:right;"> 0.7180336 </td>
   <td style="text-align:right;"> 0.5349442 </td>
   <td style="text-align:right;"> 0.0616950 </td>
   <td style="text-align:right;"> -0.0046871 </td>
   <td style="text-align:right;"> 0.1623411 </td>
   <td style="text-align:right;"> 0.1074513 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CEC </td>
   <td style="text-align:right;"> 0.5278085 </td>
   <td style="text-align:right;"> 0.2426529 </td>
   <td style="text-align:right;"> -0.0891744 </td>
   <td style="text-align:right;"> 0.2678994 </td>
   <td style="text-align:right;"> 0.2719383 </td>
   <td style="text-align:right;"> 0.2159655 </td>
   <td style="text-align:right;"> 0.3942284 </td>
   <td style="text-align:right;"> 0.2252365 </td>
   <td style="text-align:right;"> 0.3463839 </td>
   <td style="text-align:right;"> 0.3260744 </td>
   <td style="text-align:right;"> 0.1761176 </td>
   <td style="text-align:right;"> -0.6396876 </td>
   <td style="text-align:right;"> -0.2829038 </td>
   <td style="text-align:right;"> -0.3751604 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.1256450 </td>
   <td style="text-align:right;"> -0.1189595 </td>
   <td style="text-align:right;"> 0.0703113 </td>
   <td style="text-align:right;"> -0.3453354 </td>
   <td style="text-align:right;"> -0.3169741 </td>
   <td style="text-align:right;"> -0.0147856 </td>
   <td style="text-align:right;"> -0.4012123 </td>
   <td style="text-align:right;"> -0.2443690 </td>
   <td style="text-align:right;"> -0.1345758 </td>
   <td style="text-align:right;"> 0.1504463 </td>
   <td style="text-align:right;"> -0.0340426 </td>
   <td style="text-align:right;"> -0.4871177 </td>
   <td style="text-align:right;"> -0.2978615 </td>
   <td style="text-align:right;"> -0.2985350 </td>
   <td style="text-align:right;"> -0.0067555 </td>
   <td style="text-align:right;"> -0.0511688 </td>
   <td style="text-align:right;"> -0.0195514 </td>
   <td style="text-align:right;"> -0.1747329 </td>
   <td style="text-align:right;"> -0.1931675 </td>
   <td style="text-align:right;"> -0.2337801 </td>
   <td style="text-align:right;"> -0.6300267 </td>
   <td style="text-align:right;"> -0.1204253 </td>
   <td style="text-align:right;"> -0.5352327 </td>
   <td style="text-align:right;"> -0.1742889 </td>
   <td style="text-align:right;"> -0.1006438 </td>
   <td style="text-align:right;"> 0.2449507 </td>
   <td style="text-align:right;"> -0.1928564 </td>
   <td style="text-align:right;"> -0.4105992 </td>
   <td style="text-align:right;"> -0.3347895 </td>
   <td style="text-align:right;"> -0.3669278 </td>
   <td style="text-align:right;"> -0.6577265 </td>
   <td style="text-align:right;"> -0.2544973 </td>
   <td style="text-align:right;"> -0.4589151 </td>
   <td style="text-align:right;"> -0.4692880 </td>
   <td style="text-align:right;"> -0.0466410 </td>
   <td style="text-align:right;"> -0.1409265 </td>
   <td style="text-align:right;"> 0.1967519 </td>
   <td style="text-align:right;"> -0.1741065 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> -0.4786708 </td>
   <td style="text-align:right;"> 0.0791748 </td>
   <td style="text-align:right;"> -0.0915318 </td>
   <td style="text-align:right;"> 0.1282102 </td>
   <td style="text-align:right;"> 0.0338402 </td>
   <td style="text-align:right;"> -0.0084306 </td>
   <td style="text-align:right;"> 0.1313135 </td>
   <td style="text-align:right;"> 0.1260281 </td>
   <td style="text-align:right;"> 0.4570297 </td>
   <td style="text-align:right;"> 0.0350296 </td>
   <td style="text-align:right;"> -0.0870389 </td>
   <td style="text-align:right;"> 0.3480218 </td>
   <td style="text-align:right;"> 0.5513979 </td>
   <td style="text-align:right;"> 0.5716619 </td>
   <td style="text-align:right;"> -0.1256450 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.6351921 </td>
   <td style="text-align:right;"> 0.7259423 </td>
   <td style="text-align:right;"> 0.1072889 </td>
   <td style="text-align:right;"> 0.5168249 </td>
   <td style="text-align:right;"> 0.5689523 </td>
   <td style="text-align:right;"> 0.4332171 </td>
   <td style="text-align:right;"> 0.5755308 </td>
   <td style="text-align:right;"> 0.4353681 </td>
   <td style="text-align:right;"> 0.3212961 </td>
   <td style="text-align:right;"> 0.6147224 </td>
   <td style="text-align:right;"> 0.7329522 </td>
   <td style="text-align:right;"> 0.4747425 </td>
   <td style="text-align:right;"> 0.7616066 </td>
   <td style="text-align:right;"> 0.3335474 </td>
   <td style="text-align:right;"> 0.7354998 </td>
   <td style="text-align:right;"> 0.3606876 </td>
   <td style="text-align:right;"> 0.6060614 </td>
   <td style="text-align:right;"> 0.7300562 </td>
   <td style="text-align:right;"> 0.7155161 </td>
   <td style="text-align:right;"> 0.2461515 </td>
   <td style="text-align:right;"> 0.6500264 </td>
   <td style="text-align:right;"> 0.6374733 </td>
   <td style="text-align:right;"> 0.4722285 </td>
   <td style="text-align:right;"> 0.2651067 </td>
   <td style="text-align:right;"> 0.1083762 </td>
   <td style="text-align:right;"> 0.7081633 </td>
   <td style="text-align:right;"> 0.1785471 </td>
   <td style="text-align:right;"> 0.6681524 </td>
   <td style="text-align:right;"> 0.6712421 </td>
   <td style="text-align:right;"> 0.1979961 </td>
   <td style="text-align:right;"> 0.4971349 </td>
   <td style="text-align:right;"> 0.7533810 </td>
   <td style="text-align:right;"> 0.5305691 </td>
   <td style="text-align:right;"> 0.5980351 </td>
   <td style="text-align:right;"> -0.3982757 </td>
   <td style="text-align:right;"> -0.1852456 </td>
   <td style="text-align:right;"> 0.7621499 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> -0.7524823 </td>
   <td style="text-align:right;"> 0.4523301 </td>
   <td style="text-align:right;"> 0.2390248 </td>
   <td style="text-align:right;"> 0.3694524 </td>
   <td style="text-align:right;"> 0.1818629 </td>
   <td style="text-align:right;"> 0.1611085 </td>
   <td style="text-align:right;"> -0.0916029 </td>
   <td style="text-align:right;"> -0.0708605 </td>
   <td style="text-align:right;"> 0.5850148 </td>
   <td style="text-align:right;"> -0.2475307 </td>
   <td style="text-align:right;"> -0.2115733 </td>
   <td style="text-align:right;"> 0.2875622 </td>
   <td style="text-align:right;"> 0.3449104 </td>
   <td style="text-align:right;"> -0.0206494 </td>
   <td style="text-align:right;"> -0.1189595 </td>
   <td style="text-align:right;"> 0.6351921 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.8758582 </td>
   <td style="text-align:right;"> 0.2442549 </td>
   <td style="text-align:right;"> 0.6147289 </td>
   <td style="text-align:right;"> 0.4215484 </td>
   <td style="text-align:right;"> 0.5597743 </td>
   <td style="text-align:right;"> 0.5964651 </td>
   <td style="text-align:right;"> 0.3206183 </td>
   <td style="text-align:right;"> 0.2006985 </td>
   <td style="text-align:right;"> 0.4727005 </td>
   <td style="text-align:right;"> 0.5145967 </td>
   <td style="text-align:right;"> 0.2612014 </td>
   <td style="text-align:right;"> 0.6075399 </td>
   <td style="text-align:right;"> -0.0862428 </td>
   <td style="text-align:right;"> 0.5497778 </td>
   <td style="text-align:right;"> 0.5277019 </td>
   <td style="text-align:right;"> 0.3816114 </td>
   <td style="text-align:right;"> 0.5539599 </td>
   <td style="text-align:right;"> 0.7519443 </td>
   <td style="text-align:right;"> 0.2467659 </td>
   <td style="text-align:right;"> 0.9535400 </td>
   <td style="text-align:right;"> 0.5624739 </td>
   <td style="text-align:right;"> 0.6679847 </td>
   <td style="text-align:right;"> 0.1186489 </td>
   <td style="text-align:right;"> -0.2070180 </td>
   <td style="text-align:right;"> 0.7495221 </td>
   <td style="text-align:right;"> 0.3825978 </td>
   <td style="text-align:right;"> 0.4436520 </td>
   <td style="text-align:right;"> 0.6138754 </td>
   <td style="text-align:right;"> 0.2370297 </td>
   <td style="text-align:right;"> 0.5995786 </td>
   <td style="text-align:right;"> 0.5839850 </td>
   <td style="text-align:right;"> 0.5034716 </td>
   <td style="text-align:right;"> 0.8332376 </td>
   <td style="text-align:right;"> -0.6684651 </td>
   <td style="text-align:right;"> -0.5166206 </td>
   <td style="text-align:right;"> 0.8825878 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> -0.6881281 </td>
   <td style="text-align:right;"> 0.2244215 </td>
   <td style="text-align:right;"> 0.2691527 </td>
   <td style="text-align:right;"> 0.2851676 </td>
   <td style="text-align:right;"> 0.1759026 </td>
   <td style="text-align:right;"> 0.0617825 </td>
   <td style="text-align:right;"> 0.0658435 </td>
   <td style="text-align:right;"> -0.0757404 </td>
   <td style="text-align:right;"> 0.5541040 </td>
   <td style="text-align:right;"> -0.0976329 </td>
   <td style="text-align:right;"> -0.4010337 </td>
   <td style="text-align:right;"> 0.2952877 </td>
   <td style="text-align:right;"> 0.3977183 </td>
   <td style="text-align:right;"> 0.1220707 </td>
   <td style="text-align:right;"> 0.0703113 </td>
   <td style="text-align:right;"> 0.7259423 </td>
   <td style="text-align:right;"> 0.8758582 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2565816 </td>
   <td style="text-align:right;"> 0.4750401 </td>
   <td style="text-align:right;"> 0.4940548 </td>
   <td style="text-align:right;"> 0.4552300 </td>
   <td style="text-align:right;"> 0.5300990 </td>
   <td style="text-align:right;"> 0.4451348 </td>
   <td style="text-align:right;"> 0.3846242 </td>
   <td style="text-align:right;"> 0.7464376 </td>
   <td style="text-align:right;"> 0.4647209 </td>
   <td style="text-align:right;"> 0.4088404 </td>
   <td style="text-align:right;"> 0.5927055 </td>
   <td style="text-align:right;"> 0.0797110 </td>
   <td style="text-align:right;"> 0.6610746 </td>
   <td style="text-align:right;"> 0.6385938 </td>
   <td style="text-align:right;"> 0.5155371 </td>
   <td style="text-align:right;"> 0.4978122 </td>
   <td style="text-align:right;"> 0.7278452 </td>
   <td style="text-align:right;"> 0.0402958 </td>
   <td style="text-align:right;"> 0.9018618 </td>
   <td style="text-align:right;"> 0.3775674 </td>
   <td style="text-align:right;"> 0.5985135 </td>
   <td style="text-align:right;"> 0.4878040 </td>
   <td style="text-align:right;"> 0.1241041 </td>
   <td style="text-align:right;"> 0.8989520 </td>
   <td style="text-align:right;"> 0.1997049 </td>
   <td style="text-align:right;"> 0.5586852 </td>
   <td style="text-align:right;"> 0.4717185 </td>
   <td style="text-align:right;"> -0.0039871 </td>
   <td style="text-align:right;"> 0.5452433 </td>
   <td style="text-align:right;"> 0.5587800 </td>
   <td style="text-align:right;"> 0.5975853 </td>
   <td style="text-align:right;"> 0.9054749 </td>
   <td style="text-align:right;"> -0.6607967 </td>
   <td style="text-align:right;"> -0.3760829 </td>
   <td style="text-align:right;"> 0.9194690 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.1918124 </td>
   <td style="text-align:right;"> -0.0907414 </td>
   <td style="text-align:right;"> 0.0100892 </td>
   <td style="text-align:right;"> 0.0466718 </td>
   <td style="text-align:right;"> -0.3222404 </td>
   <td style="text-align:right;"> -0.2385781 </td>
   <td style="text-align:right;"> -0.3586046 </td>
   <td style="text-align:right;"> 0.0260635 </td>
   <td style="text-align:right;"> 0.1458865 </td>
   <td style="text-align:right;"> -0.3567807 </td>
   <td style="text-align:right;"> -0.5003344 </td>
   <td style="text-align:right;"> -0.1371178 </td>
   <td style="text-align:right;"> -0.1596144 </td>
   <td style="text-align:right;"> -0.2212478 </td>
   <td style="text-align:right;"> -0.3453354 </td>
   <td style="text-align:right;"> 0.1072889 </td>
   <td style="text-align:right;"> 0.2442549 </td>
   <td style="text-align:right;"> 0.2565816 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.0590291 </td>
   <td style="text-align:right;"> 0.1824074 </td>
   <td style="text-align:right;"> 0.5265322 </td>
   <td style="text-align:right;"> 0.4313938 </td>
   <td style="text-align:right;"> 0.1732468 </td>
   <td style="text-align:right;"> 0.2834743 </td>
   <td style="text-align:right;"> 0.2640349 </td>
   <td style="text-align:right;"> 0.0520705 </td>
   <td style="text-align:right;"> 0.5316522 </td>
   <td style="text-align:right;"> 0.3682907 </td>
   <td style="text-align:right;"> -0.4789534 </td>
   <td style="text-align:right;"> -0.0693717 </td>
   <td style="text-align:right;"> 0.1192397 </td>
   <td style="text-align:right;"> -0.0731678 </td>
   <td style="text-align:right;"> 0.3288737 </td>
   <td style="text-align:right;"> 0.3435548 </td>
   <td style="text-align:right;"> -0.0739376 </td>
   <td style="text-align:right;"> 0.4436714 </td>
   <td style="text-align:right;"> 0.2567102 </td>
   <td style="text-align:right;"> 0.3597064 </td>
   <td style="text-align:right;"> 0.3292876 </td>
   <td style="text-align:right;"> -0.0419796 </td>
   <td style="text-align:right;"> 0.3068174 </td>
   <td style="text-align:right;"> 0.2012880 </td>
   <td style="text-align:right;"> 0.0750472 </td>
   <td style="text-align:right;"> 0.3547808 </td>
   <td style="text-align:right;"> -0.0514128 </td>
   <td style="text-align:right;"> 0.6863451 </td>
   <td style="text-align:right;"> 0.1387526 </td>
   <td style="text-align:right;"> 0.3824594 </td>
   <td style="text-align:right;"> 0.5171061 </td>
   <td style="text-align:right;"> -0.0576496 </td>
   <td style="text-align:right;"> -0.7536652 </td>
   <td style="text-align:right;"> 0.4812776 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> -0.5732583 </td>
   <td style="text-align:right;"> 0.0267910 </td>
   <td style="text-align:right;"> 0.2390585 </td>
   <td style="text-align:right;"> -0.2144546 </td>
   <td style="text-align:right;"> 0.4662228 </td>
   <td style="text-align:right;"> 0.3068168 </td>
   <td style="text-align:right;"> -0.4066107 </td>
   <td style="text-align:right;"> -0.0678586 </td>
   <td style="text-align:right;"> 0.3509271 </td>
   <td style="text-align:right;"> -0.4746744 </td>
   <td style="text-align:right;"> 0.1089970 </td>
   <td style="text-align:right;"> 0.6009155 </td>
   <td style="text-align:right;"> 0.7236294 </td>
   <td style="text-align:right;"> 0.3494943 </td>
   <td style="text-align:right;"> -0.3169741 </td>
   <td style="text-align:right;"> 0.5168249 </td>
   <td style="text-align:right;"> 0.6147289 </td>
   <td style="text-align:right;"> 0.4750401 </td>
   <td style="text-align:right;"> -0.0590291 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.4617281 </td>
   <td style="text-align:right;"> 0.2304898 </td>
   <td style="text-align:right;"> 0.2808565 </td>
   <td style="text-align:right;"> -0.1597475 </td>
   <td style="text-align:right;"> 0.4291118 </td>
   <td style="text-align:right;"> 0.4658743 </td>
   <td style="text-align:right;"> 0.7621979 </td>
   <td style="text-align:right;"> 0.1281356 </td>
   <td style="text-align:right;"> 0.6972039 </td>
   <td style="text-align:right;"> 0.4508817 </td>
   <td style="text-align:right;"> 0.6901806 </td>
   <td style="text-align:right;"> 0.6207686 </td>
   <td style="text-align:right;"> 0.1631298 </td>
   <td style="text-align:right;"> 0.5148754 </td>
   <td style="text-align:right;"> 0.8224194 </td>
   <td style="text-align:right;"> 0.7928389 </td>
   <td style="text-align:right;"> 0.5731809 </td>
   <td style="text-align:right;"> 0.7169611 </td>
   <td style="text-align:right;"> 0.4675432 </td>
   <td style="text-align:right;"> 0.0561774 </td>
   <td style="text-align:right;"> -0.0063348 </td>
   <td style="text-align:right;"> 0.5821295 </td>
   <td style="text-align:right;"> 0.5400491 </td>
   <td style="text-align:right;"> 0.5088258 </td>
   <td style="text-align:right;"> 0.5017794 </td>
   <td style="text-align:right;"> 0.8090469 </td>
   <td style="text-align:right;"> 0.3272334 </td>
   <td style="text-align:right;"> 0.7719014 </td>
   <td style="text-align:right;"> 0.5069060 </td>
   <td style="text-align:right;"> 0.4824979 </td>
   <td style="text-align:right;"> -0.4294212 </td>
   <td style="text-align:right;"> -0.2841066 </td>
   <td style="text-align:right;"> 0.4670299 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -0.2815303 </td>
   <td style="text-align:right;"> 0.1821326 </td>
   <td style="text-align:right;"> -0.4006884 </td>
   <td style="text-align:right;"> 0.1617256 </td>
   <td style="text-align:right;"> 0.4413547 </td>
   <td style="text-align:right;"> 0.3426146 </td>
   <td style="text-align:right;"> 0.1260855 </td>
   <td style="text-align:right;"> -0.2759877 </td>
   <td style="text-align:right;"> 0.5384227 </td>
   <td style="text-align:right;"> -0.0871632 </td>
   <td style="text-align:right;"> -0.3874548 </td>
   <td style="text-align:right;"> 0.1169844 </td>
   <td style="text-align:right;"> 0.6603669 </td>
   <td style="text-align:right;"> 0.3820224 </td>
   <td style="text-align:right;"> -0.0147856 </td>
   <td style="text-align:right;"> 0.5689523 </td>
   <td style="text-align:right;"> 0.4215484 </td>
   <td style="text-align:right;"> 0.4940548 </td>
   <td style="text-align:right;"> 0.1824074 </td>
   <td style="text-align:right;"> 0.4617281 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3357546 </td>
   <td style="text-align:right;"> 0.3404110 </td>
   <td style="text-align:right;"> -0.3142665 </td>
   <td style="text-align:right;"> 0.6471009 </td>
   <td style="text-align:right;"> 0.5472780 </td>
   <td style="text-align:right;"> 0.5901962 </td>
   <td style="text-align:right;"> 0.4561498 </td>
   <td style="text-align:right;"> 0.5902370 </td>
   <td style="text-align:right;"> 0.4819358 </td>
   <td style="text-align:right;"> 0.6006924 </td>
   <td style="text-align:right;"> 0.4390801 </td>
   <td style="text-align:right;"> 0.4792509 </td>
   <td style="text-align:right;"> 0.3256641 </td>
   <td style="text-align:right;"> 0.6610850 </td>
   <td style="text-align:right;"> 0.3994200 </td>
   <td style="text-align:right;"> 0.4025222 </td>
   <td style="text-align:right;"> 0.1363309 </td>
   <td style="text-align:right;"> 0.6304528 </td>
   <td style="text-align:right;"> 0.3405765 </td>
   <td style="text-align:right;"> 0.1959483 </td>
   <td style="text-align:right;"> 0.4929873 </td>
   <td style="text-align:right;"> 0.0675309 </td>
   <td style="text-align:right;"> 0.4151545 </td>
   <td style="text-align:right;"> 0.3334066 </td>
   <td style="text-align:right;"> 0.3157012 </td>
   <td style="text-align:right;"> 0.4755774 </td>
   <td style="text-align:right;"> 0.6556621 </td>
   <td style="text-align:right;"> 0.6926927 </td>
   <td style="text-align:right;"> 0.5027962 </td>
   <td style="text-align:right;"> -0.4927414 </td>
   <td style="text-align:right;"> -0.2652725 </td>
   <td style="text-align:right;"> 0.4725342 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> -0.3909040 </td>
   <td style="text-align:right;"> 0.4705084 </td>
   <td style="text-align:right;"> 0.2603396 </td>
   <td style="text-align:right;"> 0.3569149 </td>
   <td style="text-align:right;"> -0.2692017 </td>
   <td style="text-align:right;"> -0.2154853 </td>
   <td style="text-align:right;"> -0.2125722 </td>
   <td style="text-align:right;"> -0.1900140 </td>
   <td style="text-align:right;"> 0.2103815 </td>
   <td style="text-align:right;"> -0.2921999 </td>
   <td style="text-align:right;"> -0.1211334 </td>
   <td style="text-align:right;"> 0.2043706 </td>
   <td style="text-align:right;"> 0.2436845 </td>
   <td style="text-align:right;"> 0.1506733 </td>
   <td style="text-align:right;"> -0.4012123 </td>
   <td style="text-align:right;"> 0.4332171 </td>
   <td style="text-align:right;"> 0.5597743 </td>
   <td style="text-align:right;"> 0.4552300 </td>
   <td style="text-align:right;"> 0.5265322 </td>
   <td style="text-align:right;"> 0.2304898 </td>
   <td style="text-align:right;"> 0.3357546 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.8392239 </td>
   <td style="text-align:right;"> 0.2398648 </td>
   <td style="text-align:right;"> 0.0659796 </td>
   <td style="text-align:right;"> 0.4646061 </td>
   <td style="text-align:right;"> 0.5058931 </td>
   <td style="text-align:right;"> 0.0863236 </td>
   <td style="text-align:right;"> 0.6153463 </td>
   <td style="text-align:right;"> -0.0869148 </td>
   <td style="text-align:right;"> 0.2848650 </td>
   <td style="text-align:right;"> 0.3013061 </td>
   <td style="text-align:right;"> 0.5908215 </td>
   <td style="text-align:right;"> 0.5911818 </td>
   <td style="text-align:right;"> 0.3322368 </td>
   <td style="text-align:right;"> 0.0761906 </td>
   <td style="text-align:right;"> 0.5926821 </td>
   <td style="text-align:right;"> 0.3834425 </td>
   <td style="text-align:right;"> 0.3019536 </td>
   <td style="text-align:right;"> 0.2500276 </td>
   <td style="text-align:right;"> 0.0760422 </td>
   <td style="text-align:right;"> 0.5631244 </td>
   <td style="text-align:right;"> 0.2680727 </td>
   <td style="text-align:right;"> 0.2799497 </td>
   <td style="text-align:right;"> 0.5553014 </td>
   <td style="text-align:right;"> 0.0746563 </td>
   <td style="text-align:right;"> 0.6153451 </td>
   <td style="text-align:right;"> 0.5665743 </td>
   <td style="text-align:right;"> 0.6355951 </td>
   <td style="text-align:right;"> 0.5908747 </td>
   <td style="text-align:right;"> -0.4547195 </td>
   <td style="text-align:right;"> -0.6732367 </td>
   <td style="text-align:right;"> 0.6968867 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.2448976 </td>
   <td style="text-align:right;"> 0.2405681 </td>
   <td style="text-align:right;"> 0.0612116 </td>
   <td style="text-align:right;"> 0.2531168 </td>
   <td style="text-align:right;"> -0.2078159 </td>
   <td style="text-align:right;"> -0.2158945 </td>
   <td style="text-align:right;"> 0.0143034 </td>
   <td style="text-align:right;"> -0.1055000 </td>
   <td style="text-align:right;"> 0.0847033 </td>
   <td style="text-align:right;"> -0.0337229 </td>
   <td style="text-align:right;"> -0.2047578 </td>
   <td style="text-align:right;"> 0.0886097 </td>
   <td style="text-align:right;"> 0.2180131 </td>
   <td style="text-align:right;"> 0.2105947 </td>
   <td style="text-align:right;"> -0.2443690 </td>
   <td style="text-align:right;"> 0.5755308 </td>
   <td style="text-align:right;"> 0.5964651 </td>
   <td style="text-align:right;"> 0.5300990 </td>
   <td style="text-align:right;"> 0.4313938 </td>
   <td style="text-align:right;"> 0.2808565 </td>
   <td style="text-align:right;"> 0.3404110 </td>
   <td style="text-align:right;"> 0.8392239 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2419668 </td>
   <td style="text-align:right;"> -0.0872801 </td>
   <td style="text-align:right;"> 0.5474467 </td>
   <td style="text-align:right;"> 0.4994299 </td>
   <td style="text-align:right;"> 0.0165612 </td>
   <td style="text-align:right;"> 0.5747526 </td>
   <td style="text-align:right;"> -0.0994736 </td>
   <td style="text-align:right;"> 0.2601044 </td>
   <td style="text-align:right;"> 0.1108993 </td>
   <td style="text-align:right;"> 0.6333788 </td>
   <td style="text-align:right;"> 0.5338646 </td>
   <td style="text-align:right;"> 0.4044251 </td>
   <td style="text-align:right;"> 0.0584903 </td>
   <td style="text-align:right;"> 0.6628341 </td>
   <td style="text-align:right;"> 0.4499366 </td>
   <td style="text-align:right;"> 0.3101765 </td>
   <td style="text-align:right;"> 0.1156390 </td>
   <td style="text-align:right;"> -0.1467836 </td>
   <td style="text-align:right;"> 0.6298970 </td>
   <td style="text-align:right;"> -0.0227967 </td>
   <td style="text-align:right;"> 0.1120483 </td>
   <td style="text-align:right;"> 0.4110367 </td>
   <td style="text-align:right;"> 0.0417548 </td>
   <td style="text-align:right;"> 0.3618197 </td>
   <td style="text-align:right;"> 0.6119762 </td>
   <td style="text-align:right;"> 0.4566401 </td>
   <td style="text-align:right;"> 0.7000831 </td>
   <td style="text-align:right;"> -0.1683674 </td>
   <td style="text-align:right;"> -0.6052790 </td>
   <td style="text-align:right;"> 0.6842991 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -0.4082899 </td>
   <td style="text-align:right;"> -0.0403214 </td>
   <td style="text-align:right;"> 0.3484210 </td>
   <td style="text-align:right;"> 0.3045794 </td>
   <td style="text-align:right;"> -0.3685675 </td>
   <td style="text-align:right;"> -0.2519418 </td>
   <td style="text-align:right;"> 0.2282246 </td>
   <td style="text-align:right;"> 0.3859239 </td>
   <td style="text-align:right;"> 0.0527266 </td>
   <td style="text-align:right;"> 0.2875625 </td>
   <td style="text-align:right;"> -0.0638842 </td>
   <td style="text-align:right;"> 0.2871463 </td>
   <td style="text-align:right;"> -0.2535170 </td>
   <td style="text-align:right;"> 0.0243512 </td>
   <td style="text-align:right;"> -0.1345758 </td>
   <td style="text-align:right;"> 0.4353681 </td>
   <td style="text-align:right;"> 0.3206183 </td>
   <td style="text-align:right;"> 0.4451348 </td>
   <td style="text-align:right;"> 0.1732468 </td>
   <td style="text-align:right;"> -0.1597475 </td>
   <td style="text-align:right;"> -0.3142665 </td>
   <td style="text-align:right;"> 0.2398648 </td>
   <td style="text-align:right;"> 0.2419668 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.2725051 </td>
   <td style="text-align:right;"> 0.1818532 </td>
   <td style="text-align:right;"> 0.1264939 </td>
   <td style="text-align:right;"> 0.3653275 </td>
   <td style="text-align:right;"> 0.0634117 </td>
   <td style="text-align:right;"> -0.3487176 </td>
   <td style="text-align:right;"> 0.0428218 </td>
   <td style="text-align:right;"> -0.0372837 </td>
   <td style="text-align:right;"> 0.2873455 </td>
   <td style="text-align:right;"> 0.3908089 </td>
   <td style="text-align:right;"> 0.0704608 </td>
   <td style="text-align:right;"> -0.3967991 </td>
   <td style="text-align:right;"> 0.3893285 </td>
   <td style="text-align:right;"> 0.2651898 </td>
   <td style="text-align:right;"> 0.1246656 </td>
   <td style="text-align:right;"> 0.2136055 </td>
   <td style="text-align:right;"> -0.0508816 </td>
   <td style="text-align:right;"> 0.3511412 </td>
   <td style="text-align:right;"> -0.0259897 </td>
   <td style="text-align:right;"> 0.2975372 </td>
   <td style="text-align:right;"> 0.2685467 </td>
   <td style="text-align:right;"> -0.3643828 </td>
   <td style="text-align:right;"> 0.2372717 </td>
   <td style="text-align:right;"> -0.0181701 </td>
   <td style="text-align:right;"> 0.0373242 </td>
   <td style="text-align:right;"> 0.2464706 </td>
   <td style="text-align:right;"> -0.0762704 </td>
   <td style="text-align:right;"> 0.0943284 </td>
   <td style="text-align:right;"> 0.5219961 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> -0.1873223 </td>
   <td style="text-align:right;"> -0.0071323 </td>
   <td style="text-align:right;"> 0.0324905 </td>
   <td style="text-align:right;"> -0.1924477 </td>
   <td style="text-align:right;"> 0.3282038 </td>
   <td style="text-align:right;"> 0.1714115 </td>
   <td style="text-align:right;"> -0.3809930 </td>
   <td style="text-align:right;"> -0.0668582 </td>
   <td style="text-align:right;"> 0.6067756 </td>
   <td style="text-align:right;"> -0.5154323 </td>
   <td style="text-align:right;"> -0.1482268 </td>
   <td style="text-align:right;"> 0.0643033 </td>
   <td style="text-align:right;"> 0.5369207 </td>
   <td style="text-align:right;"> 0.1898888 </td>
   <td style="text-align:right;"> 0.1504463 </td>
   <td style="text-align:right;"> 0.3212961 </td>
   <td style="text-align:right;"> 0.2006985 </td>
   <td style="text-align:right;"> 0.3846242 </td>
   <td style="text-align:right;"> 0.2834743 </td>
   <td style="text-align:right;"> 0.4291118 </td>
   <td style="text-align:right;"> 0.6471009 </td>
   <td style="text-align:right;"> 0.0659796 </td>
   <td style="text-align:right;"> -0.0872801 </td>
   <td style="text-align:right;"> -0.2725051 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.5219083 </td>
   <td style="text-align:right;"> 0.2927442 </td>
   <td style="text-align:right;"> 0.4538686 </td>
   <td style="text-align:right;"> 0.6091247 </td>
   <td style="text-align:right;"> 0.4652722 </td>
   <td style="text-align:right;"> 0.6502975 </td>
   <td style="text-align:right;"> 0.7457794 </td>
   <td style="text-align:right;"> -0.0501904 </td>
   <td style="text-align:right;"> 0.3532298 </td>
   <td style="text-align:right;"> 0.6055645 </td>
   <td style="text-align:right;"> 0.2602985 </td>
   <td style="text-align:right;"> 0.2756478 </td>
   <td style="text-align:right;"> 0.1669929 </td>
   <td style="text-align:right;"> 0.3212792 </td>
   <td style="text-align:right;"> 0.5341332 </td>
   <td style="text-align:right;"> 0.6273665 </td>
   <td style="text-align:right;"> 0.3964797 </td>
   <td style="text-align:right;"> 0.4348529 </td>
   <td style="text-align:right;"> 0.5459679 </td>
   <td style="text-align:right;"> 0.3843371 </td>
   <td style="text-align:right;"> 0.2213689 </td>
   <td style="text-align:right;"> 0.5857368 </td>
   <td style="text-align:right;"> 0.3929039 </td>
   <td style="text-align:right;"> 0.5433404 </td>
   <td style="text-align:right;"> 0.3945620 </td>
   <td style="text-align:right;"> -0.5786438 </td>
   <td style="text-align:right;"> -0.3473329 </td>
   <td style="text-align:right;"> 0.3354050 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> -0.4092369 </td>
   <td style="text-align:right;"> -0.1489583 </td>
   <td style="text-align:right;"> 0.3275359 </td>
   <td style="text-align:right;"> -0.1009288 </td>
   <td style="text-align:right;"> 0.1444066 </td>
   <td style="text-align:right;"> -0.1120529 </td>
   <td style="text-align:right;"> -0.0534254 </td>
   <td style="text-align:right;"> -0.2296133 </td>
   <td style="text-align:right;"> 0.1833791 </td>
   <td style="text-align:right;"> -0.1501120 </td>
   <td style="text-align:right;"> -0.4124044 </td>
   <td style="text-align:right;"> 0.4202140 </td>
   <td style="text-align:right;"> 0.6087493 </td>
   <td style="text-align:right;"> 0.4577474 </td>
   <td style="text-align:right;"> -0.0340426 </td>
   <td style="text-align:right;"> 0.6147224 </td>
   <td style="text-align:right;"> 0.4727005 </td>
   <td style="text-align:right;"> 0.7464376 </td>
   <td style="text-align:right;"> 0.2640349 </td>
   <td style="text-align:right;"> 0.4658743 </td>
   <td style="text-align:right;"> 0.5472780 </td>
   <td style="text-align:right;"> 0.4646061 </td>
   <td style="text-align:right;"> 0.5474467 </td>
   <td style="text-align:right;"> 0.1818532 </td>
   <td style="text-align:right;"> 0.5219083 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.5578647 </td>
   <td style="text-align:right;"> 0.2767281 </td>
   <td style="text-align:right;"> 0.6401258 </td>
   <td style="text-align:right;"> 0.4466485 </td>
   <td style="text-align:right;"> 0.6833425 </td>
   <td style="text-align:right;"> 0.6812417 </td>
   <td style="text-align:right;"> 0.6332077 </td>
   <td style="text-align:right;"> 0.4040170 </td>
   <td style="text-align:right;"> 0.5949567 </td>
   <td style="text-align:right;"> 0.1421149 </td>
   <td style="text-align:right;"> 0.5903817 </td>
   <td style="text-align:right;"> 0.2562439 </td>
   <td style="text-align:right;"> 0.3067591 </td>
   <td style="text-align:right;"> 0.7691030 </td>
   <td style="text-align:right;"> 0.5093678 </td>
   <td style="text-align:right;"> 0.9145146 </td>
   <td style="text-align:right;"> 0.0654532 </td>
   <td style="text-align:right;"> 0.5283532 </td>
   <td style="text-align:right;"> 0.2372470 </td>
   <td style="text-align:right;"> 0.0862332 </td>
   <td style="text-align:right;"> 0.3226870 </td>
   <td style="text-align:right;"> 0.6861475 </td>
   <td style="text-align:right;"> 0.7550992 </td>
   <td style="text-align:right;"> 0.7804699 </td>
   <td style="text-align:right;"> -0.4128823 </td>
   <td style="text-align:right;"> -0.3199132 </td>
   <td style="text-align:right;"> 0.6821804 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> -0.5374950 </td>
   <td style="text-align:right;"> 0.0168897 </td>
   <td style="text-align:right;"> 0.0870459 </td>
   <td style="text-align:right;"> 0.0893930 </td>
   <td style="text-align:right;"> 0.3551515 </td>
   <td style="text-align:right;"> 0.3075865 </td>
   <td style="text-align:right;"> -0.0446908 </td>
   <td style="text-align:right;"> 0.0742155 </td>
   <td style="text-align:right;"> 0.2992517 </td>
   <td style="text-align:right;"> -0.0980013 </td>
   <td style="text-align:right;"> 0.0342380 </td>
   <td style="text-align:right;"> 0.7153340 </td>
   <td style="text-align:right;"> 0.6580000 </td>
   <td style="text-align:right;"> 0.5559495 </td>
   <td style="text-align:right;"> -0.4871177 </td>
   <td style="text-align:right;"> 0.7329522 </td>
   <td style="text-align:right;"> 0.5145967 </td>
   <td style="text-align:right;"> 0.4647209 </td>
   <td style="text-align:right;"> 0.0520705 </td>
   <td style="text-align:right;"> 0.7621979 </td>
   <td style="text-align:right;"> 0.5901962 </td>
   <td style="text-align:right;"> 0.5058931 </td>
   <td style="text-align:right;"> 0.4994299 </td>
   <td style="text-align:right;"> 0.1264939 </td>
   <td style="text-align:right;"> 0.2927442 </td>
   <td style="text-align:right;"> 0.5578647 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3754631 </td>
   <td style="text-align:right;"> 0.6928877 </td>
   <td style="text-align:right;"> 0.4354274 </td>
   <td style="text-align:right;"> 0.5873501 </td>
   <td style="text-align:right;"> 0.3749863 </td>
   <td style="text-align:right;"> 0.5607952 </td>
   <td style="text-align:right;"> 0.7130496 </td>
   <td style="text-align:right;"> 0.7285554 </td>
   <td style="text-align:right;"> 0.6068560 </td>
   <td style="text-align:right;"> 0.5160752 </td>
   <td style="text-align:right;"> 0.6250553 </td>
   <td style="text-align:right;"> 0.5836229 </td>
   <td style="text-align:right;"> 0.1730824 </td>
   <td style="text-align:right;"> 0.0376492 </td>
   <td style="text-align:right;"> 0.6109257 </td>
   <td style="text-align:right;"> 0.2375350 </td>
   <td style="text-align:right;"> 0.5254704 </td>
   <td style="text-align:right;"> 0.4954871 </td>
   <td style="text-align:right;"> 0.6254837 </td>
   <td style="text-align:right;"> 0.4127620 </td>
   <td style="text-align:right;"> 0.8226371 </td>
   <td style="text-align:right;"> 0.6558231 </td>
   <td style="text-align:right;"> 0.4120752 </td>
   <td style="text-align:right;"> -0.3455088 </td>
   <td style="text-align:right;"> -0.1209037 </td>
   <td style="text-align:right;"> 0.5947196 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -0.4572435 </td>
   <td style="text-align:right;"> -0.2822302 </td>
   <td style="text-align:right;"> -0.2494951 </td>
   <td style="text-align:right;"> 0.1204422 </td>
   <td style="text-align:right;"> 0.1252168 </td>
   <td style="text-align:right;"> 0.2033120 </td>
   <td style="text-align:right;"> 0.0914941 </td>
   <td style="text-align:right;"> 0.2375354 </td>
   <td style="text-align:right;"> 0.3983270 </td>
   <td style="text-align:right;"> 0.0404411 </td>
   <td style="text-align:right;"> -0.5003482 </td>
   <td style="text-align:right;"> 0.2547761 </td>
   <td style="text-align:right;"> 0.1151279 </td>
   <td style="text-align:right;"> 0.1217595 </td>
   <td style="text-align:right;"> -0.2978615 </td>
   <td style="text-align:right;"> 0.4747425 </td>
   <td style="text-align:right;"> 0.2612014 </td>
   <td style="text-align:right;"> 0.4088404 </td>
   <td style="text-align:right;"> 0.5316522 </td>
   <td style="text-align:right;"> 0.1281356 </td>
   <td style="text-align:right;"> 0.4561498 </td>
   <td style="text-align:right;"> 0.0863236 </td>
   <td style="text-align:right;"> 0.0165612 </td>
   <td style="text-align:right;"> 0.3653275 </td>
   <td style="text-align:right;"> 0.4538686 </td>
   <td style="text-align:right;"> 0.2767281 </td>
   <td style="text-align:right;"> 0.3754631 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3187656 </td>
   <td style="text-align:right;"> -0.0792166 </td>
   <td style="text-align:right;"> 0.2227632 </td>
   <td style="text-align:right;"> 0.1586953 </td>
   <td style="text-align:right;"> 0.0706148 </td>
   <td style="text-align:right;"> 0.3851757 </td>
   <td style="text-align:right;"> 0.5778641 </td>
   <td style="text-align:right;"> 0.1250655 </td>
   <td style="text-align:right;"> 0.3737274 </td>
   <td style="text-align:right;"> 0.2634514 </td>
   <td style="text-align:right;"> 0.6791056 </td>
   <td style="text-align:right;"> 0.3917704 </td>
   <td style="text-align:right;"> -0.0105170 </td>
   <td style="text-align:right;"> 0.3368321 </td>
   <td style="text-align:right;"> 0.1034316 </td>
   <td style="text-align:right;"> 0.4521939 </td>
   <td style="text-align:right;"> 0.3549600 </td>
   <td style="text-align:right;"> 0.1245253 </td>
   <td style="text-align:right;"> 0.6215820 </td>
   <td style="text-align:right;"> 0.2266872 </td>
   <td style="text-align:right;"> 0.4350536 </td>
   <td style="text-align:right;"> 0.3261663 </td>
   <td style="text-align:right;"> -0.1615376 </td>
   <td style="text-align:right;"> -0.0799177 </td>
   <td style="text-align:right;"> 0.5011157 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> -0.4776855 </td>
   <td style="text-align:right;"> 0.1650874 </td>
   <td style="text-align:right;"> 0.1219456 </td>
   <td style="text-align:right;"> -0.1546067 </td>
   <td style="text-align:right;"> -0.0607473 </td>
   <td style="text-align:right;"> -0.1603868 </td>
   <td style="text-align:right;"> -0.4812312 </td>
   <td style="text-align:right;"> -0.1155071 </td>
   <td style="text-align:right;"> 0.4590767 </td>
   <td style="text-align:right;"> -0.5812730 </td>
   <td style="text-align:right;"> -0.0153457 </td>
   <td style="text-align:right;"> 0.3144246 </td>
   <td style="text-align:right;"> 0.7060497 </td>
   <td style="text-align:right;"> 0.5036646 </td>
   <td style="text-align:right;"> -0.2985350 </td>
   <td style="text-align:right;"> 0.7616066 </td>
   <td style="text-align:right;"> 0.6075399 </td>
   <td style="text-align:right;"> 0.5927055 </td>
   <td style="text-align:right;"> 0.3682907 </td>
   <td style="text-align:right;"> 0.6972039 </td>
   <td style="text-align:right;"> 0.5902370 </td>
   <td style="text-align:right;"> 0.6153463 </td>
   <td style="text-align:right;"> 0.5747526 </td>
   <td style="text-align:right;"> 0.0634117 </td>
   <td style="text-align:right;"> 0.6091247 </td>
   <td style="text-align:right;"> 0.6401258 </td>
   <td style="text-align:right;"> 0.6928877 </td>
   <td style="text-align:right;"> 0.3187656 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3926584 </td>
   <td style="text-align:right;"> 0.8224318 </td>
   <td style="text-align:right;"> 0.6698242 </td>
   <td style="text-align:right;"> 0.3625210 </td>
   <td style="text-align:right;"> 0.7504284 </td>
   <td style="text-align:right;"> 0.7816907 </td>
   <td style="text-align:right;"> 0.4954554 </td>
   <td style="text-align:right;"> 0.6441189 </td>
   <td style="text-align:right;"> 0.7641897 </td>
   <td style="text-align:right;"> 0.3401732 </td>
   <td style="text-align:right;"> 0.3507908 </td>
   <td style="text-align:right;"> 0.3137763 </td>
   <td style="text-align:right;"> 0.7063244 </td>
   <td style="text-align:right;"> 0.6217246 </td>
   <td style="text-align:right;"> 0.7289069 </td>
   <td style="text-align:right;"> 0.8491263 </td>
   <td style="text-align:right;"> 0.4589792 </td>
   <td style="text-align:right;"> 0.7154843 </td>
   <td style="text-align:right;"> 0.8724481 </td>
   <td style="text-align:right;"> 0.7189703 </td>
   <td style="text-align:right;"> 0.6622139 </td>
   <td style="text-align:right;"> -0.5543335 </td>
   <td style="text-align:right;"> -0.6347501 </td>
   <td style="text-align:right;"> 0.7051819 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -0.0913140 </td>
   <td style="text-align:right;"> -0.1311561 </td>
   <td style="text-align:right;"> -0.0282349 </td>
   <td style="text-align:right;"> -0.4478025 </td>
   <td style="text-align:right;"> 0.2372384 </td>
   <td style="text-align:right;"> -0.0277883 </td>
   <td style="text-align:right;"> -0.0734986 </td>
   <td style="text-align:right;"> -0.4100490 </td>
   <td style="text-align:right;"> -0.0127906 </td>
   <td style="text-align:right;"> -0.1625928 </td>
   <td style="text-align:right;"> 0.0926246 </td>
   <td style="text-align:right;"> 0.4148064 </td>
   <td style="text-align:right;"> 0.8603869 </td>
   <td style="text-align:right;"> 0.8048582 </td>
   <td style="text-align:right;"> -0.0067555 </td>
   <td style="text-align:right;"> 0.3335474 </td>
   <td style="text-align:right;"> -0.0862428 </td>
   <td style="text-align:right;"> 0.0797110 </td>
   <td style="text-align:right;"> -0.4789534 </td>
   <td style="text-align:right;"> 0.4508817 </td>
   <td style="text-align:right;"> 0.4819358 </td>
   <td style="text-align:right;"> -0.0869148 </td>
   <td style="text-align:right;"> -0.0994736 </td>
   <td style="text-align:right;"> -0.3487176 </td>
   <td style="text-align:right;"> 0.4652722 </td>
   <td style="text-align:right;"> 0.4466485 </td>
   <td style="text-align:right;"> 0.4354274 </td>
   <td style="text-align:right;"> -0.0792166 </td>
   <td style="text-align:right;"> 0.3926584 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.7079528 </td>
   <td style="text-align:right;"> 0.4992616 </td>
   <td style="text-align:right;"> 0.4131319 </td>
   <td style="text-align:right;"> 0.0130131 </td>
   <td style="text-align:right;"> 0.2040231 </td>
   <td style="text-align:right;"> 0.5093248 </td>
   <td style="text-align:right;"> -0.1728693 </td>
   <td style="text-align:right;"> 0.0799592 </td>
   <td style="text-align:right;"> -0.1631011 </td>
   <td style="text-align:right;"> 0.3960579 </td>
   <td style="text-align:right;"> 0.6136143 </td>
   <td style="text-align:right;"> 0.2813408 </td>
   <td style="text-align:right;"> 0.1877935 </td>
   <td style="text-align:right;"> 0.5571466 </td>
   <td style="text-align:right;"> 0.1035302 </td>
   <td style="text-align:right;"> 0.4132668 </td>
   <td style="text-align:right;"> -0.1432812 </td>
   <td style="text-align:right;"> 0.5747898 </td>
   <td style="text-align:right;"> 0.4613204 </td>
   <td style="text-align:right;"> -0.0201596 </td>
   <td style="text-align:right;"> -0.2457316 </td>
   <td style="text-align:right;"> 0.2534905 </td>
   <td style="text-align:right;"> -0.0701946 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> -0.5470097 </td>
   <td style="text-align:right;"> 0.1385272 </td>
   <td style="text-align:right;"> 0.1630347 </td>
   <td style="text-align:right;"> -0.2299454 </td>
   <td style="text-align:right;"> 0.1097366 </td>
   <td style="text-align:right;"> -0.1066563 </td>
   <td style="text-align:right;"> -0.2771120 </td>
   <td style="text-align:right;"> -0.3023832 </td>
   <td style="text-align:right;"> 0.4377405 </td>
   <td style="text-align:right;"> -0.4352023 </td>
   <td style="text-align:right;"> -0.0502045 </td>
   <td style="text-align:right;"> 0.4232019 </td>
   <td style="text-align:right;"> 0.8797579 </td>
   <td style="text-align:right;"> 0.6416481 </td>
   <td style="text-align:right;"> -0.0511688 </td>
   <td style="text-align:right;"> 0.7354998 </td>
   <td style="text-align:right;"> 0.5497778 </td>
   <td style="text-align:right;"> 0.6610746 </td>
   <td style="text-align:right;"> -0.0693717 </td>
   <td style="text-align:right;"> 0.6901806 </td>
   <td style="text-align:right;"> 0.6006924 </td>
   <td style="text-align:right;"> 0.2848650 </td>
   <td style="text-align:right;"> 0.2601044 </td>
   <td style="text-align:right;"> 0.0428218 </td>
   <td style="text-align:right;"> 0.6502975 </td>
   <td style="text-align:right;"> 0.6833425 </td>
   <td style="text-align:right;"> 0.5873501 </td>
   <td style="text-align:right;"> 0.2227632 </td>
   <td style="text-align:right;"> 0.8224318 </td>
   <td style="text-align:right;"> 0.7079528 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.8242375 </td>
   <td style="text-align:right;"> 0.4587244 </td>
   <td style="text-align:right;"> 0.4602229 </td>
   <td style="text-align:right;"> 0.6752230 </td>
   <td style="text-align:right;"> 0.4803340 </td>
   <td style="text-align:right;"> 0.4951427 </td>
   <td style="text-align:right;"> 0.5205709 </td>
   <td style="text-align:right;"> 0.2125266 </td>
   <td style="text-align:right;"> 0.5038684 </td>
   <td style="text-align:right;"> 0.5088342 </td>
   <td style="text-align:right;"> 0.7258697 </td>
   <td style="text-align:right;"> 0.5517272 </td>
   <td style="text-align:right;"> 0.8747534 </td>
   <td style="text-align:right;"> 0.6478655 </td>
   <td style="text-align:right;"> 0.3910145 </td>
   <td style="text-align:right;"> 0.4358423 </td>
   <td style="text-align:right;"> 0.8131381 </td>
   <td style="text-align:right;"> 0.6973238 </td>
   <td style="text-align:right;"> 0.5581929 </td>
   <td style="text-align:right;"> -0.6671117 </td>
   <td style="text-align:right;"> -0.2484420 </td>
   <td style="text-align:right;"> 0.5668516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -0.5895328 </td>
   <td style="text-align:right;"> 0.1802071 </td>
   <td style="text-align:right;"> 0.5183584 </td>
   <td style="text-align:right;"> -0.1961637 </td>
   <td style="text-align:right;"> 0.1879125 </td>
   <td style="text-align:right;"> -0.0486015 </td>
   <td style="text-align:right;"> -0.5105520 </td>
   <td style="text-align:right;"> -0.3606594 </td>
   <td style="text-align:right;"> 0.4437190 </td>
   <td style="text-align:right;"> -0.6632347 </td>
   <td style="text-align:right;"> -0.1376247 </td>
   <td style="text-align:right;"> 0.4115959 </td>
   <td style="text-align:right;"> 0.6949345 </td>
   <td style="text-align:right;"> 0.2619831 </td>
   <td style="text-align:right;"> -0.0195514 </td>
   <td style="text-align:right;"> 0.3606876 </td>
   <td style="text-align:right;"> 0.5277019 </td>
   <td style="text-align:right;"> 0.6385938 </td>
   <td style="text-align:right;"> 0.1192397 </td>
   <td style="text-align:right;"> 0.6207686 </td>
   <td style="text-align:right;"> 0.4390801 </td>
   <td style="text-align:right;"> 0.3013061 </td>
   <td style="text-align:right;"> 0.1108993 </td>
   <td style="text-align:right;"> -0.0372837 </td>
   <td style="text-align:right;"> 0.7457794 </td>
   <td style="text-align:right;"> 0.6812417 </td>
   <td style="text-align:right;"> 0.3749863 </td>
   <td style="text-align:right;"> 0.1586953 </td>
   <td style="text-align:right;"> 0.6698242 </td>
   <td style="text-align:right;"> 0.4992616 </td>
   <td style="text-align:right;"> 0.8242375 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2237020 </td>
   <td style="text-align:right;"> 0.3221978 </td>
   <td style="text-align:right;"> 0.5943551 </td>
   <td style="text-align:right;"> 0.3544901 </td>
   <td style="text-align:right;"> 0.5182924 </td>
   <td style="text-align:right;"> 0.3166905 </td>
   <td style="text-align:right;"> 0.2117597 </td>
   <td style="text-align:right;"> 0.6448734 </td>
   <td style="text-align:right;"> 0.6107666 </td>
   <td style="text-align:right;"> 0.7102881 </td>
   <td style="text-align:right;"> 0.6493740 </td>
   <td style="text-align:right;"> 0.7312017 </td>
   <td style="text-align:right;"> 0.4751019 </td>
   <td style="text-align:right;"> 0.3065580 </td>
   <td style="text-align:right;"> 0.5019571 </td>
   <td style="text-align:right;"> 0.5800730 </td>
   <td style="text-align:right;"> 0.7117702 </td>
   <td style="text-align:right;"> 0.6092980 </td>
   <td style="text-align:right;"> -0.7911389 </td>
   <td style="text-align:right;"> -0.4073193 </td>
   <td style="text-align:right;"> 0.5412015 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> -0.3670954 </td>
   <td style="text-align:right;"> 0.1856233 </td>
   <td style="text-align:right;"> 0.0294007 </td>
   <td style="text-align:right;"> 0.2603343 </td>
   <td style="text-align:right;"> -0.0487667 </td>
   <td style="text-align:right;"> -0.1467302 </td>
   <td style="text-align:right;"> 0.3955915 </td>
   <td style="text-align:right;"> -0.3761779 </td>
   <td style="text-align:right;"> -0.0210200 </td>
   <td style="text-align:right;"> 0.2746022 </td>
   <td style="text-align:right;"> -0.2987831 </td>
   <td style="text-align:right;"> 0.4122684 </td>
   <td style="text-align:right;"> 0.5122723 </td>
   <td style="text-align:right;"> 0.6044156 </td>
   <td style="text-align:right;"> -0.1747329 </td>
   <td style="text-align:right;"> 0.6060614 </td>
   <td style="text-align:right;"> 0.3816114 </td>
   <td style="text-align:right;"> 0.5155371 </td>
   <td style="text-align:right;"> -0.0731678 </td>
   <td style="text-align:right;"> 0.1631298 </td>
   <td style="text-align:right;"> 0.4792509 </td>
   <td style="text-align:right;"> 0.5908215 </td>
   <td style="text-align:right;"> 0.6333788 </td>
   <td style="text-align:right;"> 0.2873455 </td>
   <td style="text-align:right;"> -0.0501904 </td>
   <td style="text-align:right;"> 0.6332077 </td>
   <td style="text-align:right;"> 0.5607952 </td>
   <td style="text-align:right;"> 0.0706148 </td>
   <td style="text-align:right;"> 0.3625210 </td>
   <td style="text-align:right;"> 0.4131319 </td>
   <td style="text-align:right;"> 0.4587244 </td>
   <td style="text-align:right;"> 0.2237020 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2358929 </td>
   <td style="text-align:right;"> 0.1922827 </td>
   <td style="text-align:right;"> 0.0680747 </td>
   <td style="text-align:right;"> 0.3345599 </td>
   <td style="text-align:right;"> 0.0663694 </td>
   <td style="text-align:right;"> 0.1943130 </td>
   <td style="text-align:right;"> 0.4483846 </td>
   <td style="text-align:right;"> 0.2114937 </td>
   <td style="text-align:right;"> 0.6332613 </td>
   <td style="text-align:right;"> -0.1875635 </td>
   <td style="text-align:right;"> 0.3956055 </td>
   <td style="text-align:right;"> 0.1668547 </td>
   <td style="text-align:right;"> -0.0197972 </td>
   <td style="text-align:right;"> 0.0784152 </td>
   <td style="text-align:right;"> 0.6163482 </td>
   <td style="text-align:right;"> 0.6310694 </td>
   <td style="text-align:right;"> 0.4317252 </td>
   <td style="text-align:right;"> -0.2952387 </td>
   <td style="text-align:right;"> 0.0171207 </td>
   <td style="text-align:right;"> 0.5112659 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> -0.3271235 </td>
   <td style="text-align:right;"> 0.2472663 </td>
   <td style="text-align:right;"> 0.2122196 </td>
   <td style="text-align:right;"> 0.2876959 </td>
   <td style="text-align:right;"> 0.0765954 </td>
   <td style="text-align:right;"> 0.1658216 </td>
   <td style="text-align:right;"> -0.2052591 </td>
   <td style="text-align:right;"> 0.4991277 </td>
   <td style="text-align:right;"> 0.6007881 </td>
   <td style="text-align:right;"> -0.2174518 </td>
   <td style="text-align:right;"> 0.3016631 </td>
   <td style="text-align:right;"> 0.2847252 </td>
   <td style="text-align:right;"> 0.2356105 </td>
   <td style="text-align:right;"> 0.1680644 </td>
   <td style="text-align:right;"> -0.1931675 </td>
   <td style="text-align:right;"> 0.7300562 </td>
   <td style="text-align:right;"> 0.5539599 </td>
   <td style="text-align:right;"> 0.4978122 </td>
   <td style="text-align:right;"> 0.3288737 </td>
   <td style="text-align:right;"> 0.5148754 </td>
   <td style="text-align:right;"> 0.3256641 </td>
   <td style="text-align:right;"> 0.5911818 </td>
   <td style="text-align:right;"> 0.5338646 </td>
   <td style="text-align:right;"> 0.3908089 </td>
   <td style="text-align:right;"> 0.3532298 </td>
   <td style="text-align:right;"> 0.4040170 </td>
   <td style="text-align:right;"> 0.7130496 </td>
   <td style="text-align:right;"> 0.3851757 </td>
   <td style="text-align:right;"> 0.7504284 </td>
   <td style="text-align:right;"> 0.0130131 </td>
   <td style="text-align:right;"> 0.4602229 </td>
   <td style="text-align:right;"> 0.3221978 </td>
   <td style="text-align:right;"> 0.2358929 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.6625044 </td>
   <td style="text-align:right;"> 0.1487245 </td>
   <td style="text-align:right;"> 0.6375605 </td>
   <td style="text-align:right;"> 0.7045650 </td>
   <td style="text-align:right;"> 0.4842155 </td>
   <td style="text-align:right;"> 0.0270460 </td>
   <td style="text-align:right;"> 0.0956780 </td>
   <td style="text-align:right;"> 0.4812212 </td>
   <td style="text-align:right;"> 0.3570380 </td>
   <td style="text-align:right;"> 0.4250268 </td>
   <td style="text-align:right;"> 0.7271529 </td>
   <td style="text-align:right;"> 0.2142938 </td>
   <td style="text-align:right;"> 0.7158626 </td>
   <td style="text-align:right;"> 0.5282233 </td>
   <td style="text-align:right;"> 0.3625363 </td>
   <td style="text-align:right;"> 0.4572592 </td>
   <td style="text-align:right;"> -0.4522680 </td>
   <td style="text-align:right;"> -0.4096391 </td>
   <td style="text-align:right;"> 0.7181758 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> -0.6289787 </td>
   <td style="text-align:right;"> -0.0042868 </td>
   <td style="text-align:right;"> 0.0347987 </td>
   <td style="text-align:right;"> 0.0000615 </td>
   <td style="text-align:right;"> 0.4021136 </td>
   <td style="text-align:right;"> 0.3135620 </td>
   <td style="text-align:right;"> -0.2419448 </td>
   <td style="text-align:right;"> 0.0797681 </td>
   <td style="text-align:right;"> 0.6001283 </td>
   <td style="text-align:right;"> -0.3518337 </td>
   <td style="text-align:right;"> -0.2313635 </td>
   <td style="text-align:right;"> 0.4049609 </td>
   <td style="text-align:right;"> 0.5663707 </td>
   <td style="text-align:right;"> 0.2225495 </td>
   <td style="text-align:right;"> -0.2337801 </td>
   <td style="text-align:right;"> 0.7155161 </td>
   <td style="text-align:right;"> 0.7519443 </td>
   <td style="text-align:right;"> 0.7278452 </td>
   <td style="text-align:right;"> 0.3435548 </td>
   <td style="text-align:right;"> 0.8224194 </td>
   <td style="text-align:right;"> 0.6610850 </td>
   <td style="text-align:right;"> 0.3322368 </td>
   <td style="text-align:right;"> 0.4044251 </td>
   <td style="text-align:right;"> 0.0704608 </td>
   <td style="text-align:right;"> 0.6055645 </td>
   <td style="text-align:right;"> 0.5949567 </td>
   <td style="text-align:right;"> 0.7285554 </td>
   <td style="text-align:right;"> 0.5778641 </td>
   <td style="text-align:right;"> 0.7816907 </td>
   <td style="text-align:right;"> 0.2040231 </td>
   <td style="text-align:right;"> 0.6752230 </td>
   <td style="text-align:right;"> 0.5943551 </td>
   <td style="text-align:right;"> 0.1922827 </td>
   <td style="text-align:right;"> 0.6625044 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.5465050 </td>
   <td style="text-align:right;"> 0.8016590 </td>
   <td style="text-align:right;"> 0.6941032 </td>
   <td style="text-align:right;"> 0.7650474 </td>
   <td style="text-align:right;"> 0.2414734 </td>
   <td style="text-align:right;"> -0.0289782 </td>
   <td style="text-align:right;"> 0.7203265 </td>
   <td style="text-align:right;"> 0.4241936 </td>
   <td style="text-align:right;"> 0.5633014 </td>
   <td style="text-align:right;"> 0.6186419 </td>
   <td style="text-align:right;"> 0.5485451 </td>
   <td style="text-align:right;"> 0.6595585 </td>
   <td style="text-align:right;"> 0.7293122 </td>
   <td style="text-align:right;"> 0.6098670 </td>
   <td style="text-align:right;"> 0.7369762 </td>
   <td style="text-align:right;"> -0.4924795 </td>
   <td style="text-align:right;"> -0.4530231 </td>
   <td style="text-align:right;"> 0.7560358 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> -0.4727843 </td>
   <td style="text-align:right;"> -0.1563835 </td>
   <td style="text-align:right;"> -0.0816276 </td>
   <td style="text-align:right;"> -0.4753756 </td>
   <td style="text-align:right;"> 0.2467022 </td>
   <td style="text-align:right;"> 0.1262534 </td>
   <td style="text-align:right;"> -0.4570208 </td>
   <td style="text-align:right;"> -0.3528549 </td>
   <td style="text-align:right;"> -0.0252173 </td>
   <td style="text-align:right;"> -0.4947465 </td>
   <td style="text-align:right;"> -0.0117934 </td>
   <td style="text-align:right;"> 0.5805498 </td>
   <td style="text-align:right;"> 0.7355172 </td>
   <td style="text-align:right;"> 0.5123814 </td>
   <td style="text-align:right;"> -0.6300267 </td>
   <td style="text-align:right;"> 0.2461515 </td>
   <td style="text-align:right;"> 0.2467659 </td>
   <td style="text-align:right;"> 0.0402958 </td>
   <td style="text-align:right;"> -0.0739376 </td>
   <td style="text-align:right;"> 0.7928389 </td>
   <td style="text-align:right;"> 0.3994200 </td>
   <td style="text-align:right;"> 0.0761906 </td>
   <td style="text-align:right;"> 0.0584903 </td>
   <td style="text-align:right;"> -0.3967991 </td>
   <td style="text-align:right;"> 0.2602985 </td>
   <td style="text-align:right;"> 0.1421149 </td>
   <td style="text-align:right;"> 0.6068560 </td>
   <td style="text-align:right;"> 0.1250655 </td>
   <td style="text-align:right;"> 0.4954554 </td>
   <td style="text-align:right;"> 0.5093248 </td>
   <td style="text-align:right;"> 0.4803340 </td>
   <td style="text-align:right;"> 0.3544901 </td>
   <td style="text-align:right;"> 0.0680747 </td>
   <td style="text-align:right;"> 0.1487245 </td>
   <td style="text-align:right;"> 0.5465050 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.1522551 </td>
   <td style="text-align:right;"> 0.5946744 </td>
   <td style="text-align:right;"> 0.2560886 </td>
   <td style="text-align:right;"> -0.0146505 </td>
   <td style="text-align:right;"> -0.1233475 </td>
   <td style="text-align:right;"> 0.2470535 </td>
   <td style="text-align:right;"> 0.5493776 </td>
   <td style="text-align:right;"> 0.4546046 </td>
   <td style="text-align:right;"> 0.3912704 </td>
   <td style="text-align:right;"> 0.9826091 </td>
   <td style="text-align:right;"> 0.1433290 </td>
   <td style="text-align:right;"> 0.6942590 </td>
   <td style="text-align:right;"> 0.4766386 </td>
   <td style="text-align:right;"> 0.1131029 </td>
   <td style="text-align:right;"> -0.1043352 </td>
   <td style="text-align:right;"> -0.1437649 </td>
   <td style="text-align:right;"> 0.0755796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -0.6846201 </td>
   <td style="text-align:right;"> 0.2829309 </td>
   <td style="text-align:right;"> 0.2850472 </td>
   <td style="text-align:right;"> 0.3225066 </td>
   <td style="text-align:right;"> 0.1504787 </td>
   <td style="text-align:right;"> 0.1216399 </td>
   <td style="text-align:right;"> -0.1275314 </td>
   <td style="text-align:right;"> 0.0548097 </td>
   <td style="text-align:right;"> 0.5605500 </td>
   <td style="text-align:right;"> -0.2437681 </td>
   <td style="text-align:right;"> -0.2877258 </td>
   <td style="text-align:right;"> 0.2601642 </td>
   <td style="text-align:right;"> 0.2596322 </td>
   <td style="text-align:right;"> -0.0626878 </td>
   <td style="text-align:right;"> -0.1204253 </td>
   <td style="text-align:right;"> 0.6500264 </td>
   <td style="text-align:right;"> 0.9535400 </td>
   <td style="text-align:right;"> 0.9018618 </td>
   <td style="text-align:right;"> 0.4436714 </td>
   <td style="text-align:right;"> 0.5731809 </td>
   <td style="text-align:right;"> 0.4025222 </td>
   <td style="text-align:right;"> 0.5926821 </td>
   <td style="text-align:right;"> 0.6628341 </td>
   <td style="text-align:right;"> 0.3893285 </td>
   <td style="text-align:right;"> 0.2756478 </td>
   <td style="text-align:right;"> 0.5903817 </td>
   <td style="text-align:right;"> 0.5160752 </td>
   <td style="text-align:right;"> 0.3737274 </td>
   <td style="text-align:right;"> 0.6441189 </td>
   <td style="text-align:right;"> -0.1728693 </td>
   <td style="text-align:right;"> 0.4951427 </td>
   <td style="text-align:right;"> 0.5182924 </td>
   <td style="text-align:right;"> 0.3345599 </td>
   <td style="text-align:right;"> 0.6375605 </td>
   <td style="text-align:right;"> 0.8016590 </td>
   <td style="text-align:right;"> 0.1522551 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.5821847 </td>
   <td style="text-align:right;"> 0.6953924 </td>
   <td style="text-align:right;"> 0.2254210 </td>
   <td style="text-align:right;"> -0.1491895 </td>
   <td style="text-align:right;"> 0.8126103 </td>
   <td style="text-align:right;"> 0.3134236 </td>
   <td style="text-align:right;"> 0.3979449 </td>
   <td style="text-align:right;"> 0.5820092 </td>
   <td style="text-align:right;"> 0.1616629 </td>
   <td style="text-align:right;"> 0.6620807 </td>
   <td style="text-align:right;"> 0.5591358 </td>
   <td style="text-align:right;"> 0.5167963 </td>
   <td style="text-align:right;"> 0.9169482 </td>
   <td style="text-align:right;"> -0.5767370 </td>
   <td style="text-align:right;"> -0.6042877 </td>
   <td style="text-align:right;"> 0.9429335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> -0.5248808 </td>
   <td style="text-align:right;"> -0.0389185 </td>
   <td style="text-align:right;"> 0.1547731 </td>
   <td style="text-align:right;"> -0.2537164 </td>
   <td style="text-align:right;"> -0.1482612 </td>
   <td style="text-align:right;"> -0.1454950 </td>
   <td style="text-align:right;"> -0.5135293 </td>
   <td style="text-align:right;"> 0.1588479 </td>
   <td style="text-align:right;"> 0.2121094 </td>
   <td style="text-align:right;"> -0.4824073 </td>
   <td style="text-align:right;"> 0.1802598 </td>
   <td style="text-align:right;"> 0.4318794 </td>
   <td style="text-align:right;"> 0.4218165 </td>
   <td style="text-align:right;"> 0.3595042 </td>
   <td style="text-align:right;"> -0.5352327 </td>
   <td style="text-align:right;"> 0.6374733 </td>
   <td style="text-align:right;"> 0.5624739 </td>
   <td style="text-align:right;"> 0.3775674 </td>
   <td style="text-align:right;"> 0.2567102 </td>
   <td style="text-align:right;"> 0.7169611 </td>
   <td style="text-align:right;"> 0.1363309 </td>
   <td style="text-align:right;"> 0.3834425 </td>
   <td style="text-align:right;"> 0.4499366 </td>
   <td style="text-align:right;"> 0.2651898 </td>
   <td style="text-align:right;"> 0.1669929 </td>
   <td style="text-align:right;"> 0.2562439 </td>
   <td style="text-align:right;"> 0.6250553 </td>
   <td style="text-align:right;"> 0.2634514 </td>
   <td style="text-align:right;"> 0.7641897 </td>
   <td style="text-align:right;"> 0.0799592 </td>
   <td style="text-align:right;"> 0.5205709 </td>
   <td style="text-align:right;"> 0.3166905 </td>
   <td style="text-align:right;"> 0.0663694 </td>
   <td style="text-align:right;"> 0.7045650 </td>
   <td style="text-align:right;"> 0.6941032 </td>
   <td style="text-align:right;"> 0.5946744 </td>
   <td style="text-align:right;"> 0.5821847 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2865680 </td>
   <td style="text-align:right;"> -0.0686263 </td>
   <td style="text-align:right;"> -0.1752869 </td>
   <td style="text-align:right;"> 0.4702436 </td>
   <td style="text-align:right;"> 0.6398505 </td>
   <td style="text-align:right;"> 0.5477372 </td>
   <td style="text-align:right;"> 0.8084788 </td>
   <td style="text-align:right;"> 0.6261250 </td>
   <td style="text-align:right;"> 0.5043378 </td>
   <td style="text-align:right;"> 0.6844641 </td>
   <td style="text-align:right;"> 0.3391044 </td>
   <td style="text-align:right;"> 0.4257947 </td>
   <td style="text-align:right;"> -0.1802631 </td>
   <td style="text-align:right;"> -0.4545788 </td>
   <td style="text-align:right;"> 0.5387702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.5440953 </td>
   <td style="text-align:right;"> 0.1520278 </td>
   <td style="text-align:right;"> -0.1526560 </td>
   <td style="text-align:right;"> 0.5055409 </td>
   <td style="text-align:right;"> 0.5784861 </td>
   <td style="text-align:right;"> 0.6538657 </td>
   <td style="text-align:right;"> 0.1963514 </td>
   <td style="text-align:right;"> 0.2121357 </td>
   <td style="text-align:right;"> 0.6329001 </td>
   <td style="text-align:right;"> 0.0707855 </td>
   <td style="text-align:right;"> -0.3708440 </td>
   <td style="text-align:right;"> 0.2683951 </td>
   <td style="text-align:right;"> 0.1570567 </td>
   <td style="text-align:right;"> -0.1498550 </td>
   <td style="text-align:right;"> -0.1742889 </td>
   <td style="text-align:right;"> 0.4722285 </td>
   <td style="text-align:right;"> 0.6679847 </td>
   <td style="text-align:right;"> 0.5985135 </td>
   <td style="text-align:right;"> 0.3597064 </td>
   <td style="text-align:right;"> 0.4675432 </td>
   <td style="text-align:right;"> 0.6304528 </td>
   <td style="text-align:right;"> 0.3019536 </td>
   <td style="text-align:right;"> 0.3101765 </td>
   <td style="text-align:right;"> 0.1246656 </td>
   <td style="text-align:right;"> 0.3212792 </td>
   <td style="text-align:right;"> 0.3067591 </td>
   <td style="text-align:right;"> 0.5836229 </td>
   <td style="text-align:right;"> 0.6791056 </td>
   <td style="text-align:right;"> 0.3401732 </td>
   <td style="text-align:right;"> -0.1631011 </td>
   <td style="text-align:right;"> 0.2125266 </td>
   <td style="text-align:right;"> 0.2117597 </td>
   <td style="text-align:right;"> 0.1943130 </td>
   <td style="text-align:right;"> 0.4842155 </td>
   <td style="text-align:right;"> 0.7650474 </td>
   <td style="text-align:right;"> 0.2560886 </td>
   <td style="text-align:right;"> 0.6953924 </td>
   <td style="text-align:right;"> 0.2865680 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.0483182 </td>
   <td style="text-align:right;"> -0.3429842 </td>
   <td style="text-align:right;"> 0.4508250 </td>
   <td style="text-align:right;"> 0.0267657 </td>
   <td style="text-align:right;"> 0.1787624 </td>
   <td style="text-align:right;"> 0.2911605 </td>
   <td style="text-align:right;"> 0.2908171 </td>
   <td style="text-align:right;"> 0.5864228 </td>
   <td style="text-align:right;"> 0.3518794 </td>
   <td style="text-align:right;"> 0.4255995 </td>
   <td style="text-align:right;"> 0.5509849 </td>
   <td style="text-align:right;"> -0.4017788 </td>
   <td style="text-align:right;"> -0.2316881 </td>
   <td style="text-align:right;"> 0.6655955 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.3904141 </td>
   <td style="text-align:right;"> -0.2682738 </td>
   <td style="text-align:right;"> 0.3214070 </td>
   <td style="text-align:right;"> -0.2421058 </td>
   <td style="text-align:right;"> -0.1251021 </td>
   <td style="text-align:right;"> -0.3514561 </td>
   <td style="text-align:right;"> -0.1255198 </td>
   <td style="text-align:right;"> -0.4267040 </td>
   <td style="text-align:right;"> -0.0149174 </td>
   <td style="text-align:right;"> -0.2134242 </td>
   <td style="text-align:right;"> -0.5915026 </td>
   <td style="text-align:right;"> 0.3571332 </td>
   <td style="text-align:right;"> 0.4552276 </td>
   <td style="text-align:right;"> 0.4045777 </td>
   <td style="text-align:right;"> -0.1006438 </td>
   <td style="text-align:right;"> 0.2651067 </td>
   <td style="text-align:right;"> 0.1186489 </td>
   <td style="text-align:right;"> 0.4878040 </td>
   <td style="text-align:right;"> 0.3292876 </td>
   <td style="text-align:right;"> 0.0561774 </td>
   <td style="text-align:right;"> 0.3405765 </td>
   <td style="text-align:right;"> 0.2500276 </td>
   <td style="text-align:right;"> 0.1156390 </td>
   <td style="text-align:right;"> 0.2136055 </td>
   <td style="text-align:right;"> 0.5341332 </td>
   <td style="text-align:right;"> 0.7691030 </td>
   <td style="text-align:right;"> 0.1730824 </td>
   <td style="text-align:right;"> 0.3917704 </td>
   <td style="text-align:right;"> 0.3507908 </td>
   <td style="text-align:right;"> 0.3960579 </td>
   <td style="text-align:right;"> 0.5038684 </td>
   <td style="text-align:right;"> 0.6448734 </td>
   <td style="text-align:right;"> 0.4483846 </td>
   <td style="text-align:right;"> 0.0270460 </td>
   <td style="text-align:right;"> 0.2414734 </td>
   <td style="text-align:right;"> -0.0146505 </td>
   <td style="text-align:right;"> 0.2254210 </td>
   <td style="text-align:right;"> -0.0686263 </td>
   <td style="text-align:right;"> 0.0483182 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.7075767 </td>
   <td style="text-align:right;"> 0.6578923 </td>
   <td style="text-align:right;"> 0.1371019 </td>
   <td style="text-align:right;"> 0.5938416 </td>
   <td style="text-align:right;"> 0.0964395 </td>
   <td style="text-align:right;"> -0.0990521 </td>
   <td style="text-align:right;"> 0.2808560 </td>
   <td style="text-align:right;"> 0.3627018 </td>
   <td style="text-align:right;"> 0.7272263 </td>
   <td style="text-align:right;"> 0.4843486 </td>
   <td style="text-align:right;"> -0.3420923 </td>
   <td style="text-align:right;"> -0.1552865 </td>
   <td style="text-align:right;"> 0.4060250 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 0.0708207 </td>
   <td style="text-align:right;"> -0.0293808 </td>
   <td style="text-align:right;"> 0.3526735 </td>
   <td style="text-align:right;"> -0.2736705 </td>
   <td style="text-align:right;"> -0.0853650 </td>
   <td style="text-align:right;"> -0.2874374 </td>
   <td style="text-align:right;"> -0.2468028 </td>
   <td style="text-align:right;"> -0.2104401 </td>
   <td style="text-align:right;"> 0.1112550 </td>
   <td style="text-align:right;"> -0.3033643 </td>
   <td style="text-align:right;"> 0.0647371 </td>
   <td style="text-align:right;"> 0.1083986 </td>
   <td style="text-align:right;"> 0.4204540 </td>
   <td style="text-align:right;"> 0.3919728 </td>
   <td style="text-align:right;"> 0.2449507 </td>
   <td style="text-align:right;"> 0.1083762 </td>
   <td style="text-align:right;"> -0.2070180 </td>
   <td style="text-align:right;"> 0.1241041 </td>
   <td style="text-align:right;"> -0.0419796 </td>
   <td style="text-align:right;"> -0.0063348 </td>
   <td style="text-align:right;"> 0.1959483 </td>
   <td style="text-align:right;"> 0.0760422 </td>
   <td style="text-align:right;"> -0.1467836 </td>
   <td style="text-align:right;"> -0.0508816 </td>
   <td style="text-align:right;"> 0.6273665 </td>
   <td style="text-align:right;"> 0.5093678 </td>
   <td style="text-align:right;"> 0.0376492 </td>
   <td style="text-align:right;"> -0.0105170 </td>
   <td style="text-align:right;"> 0.3137763 </td>
   <td style="text-align:right;"> 0.6136143 </td>
   <td style="text-align:right;"> 0.5088342 </td>
   <td style="text-align:right;"> 0.6107666 </td>
   <td style="text-align:right;"> 0.2114937 </td>
   <td style="text-align:right;"> 0.0956780 </td>
   <td style="text-align:right;"> -0.0289782 </td>
   <td style="text-align:right;"> -0.1233475 </td>
   <td style="text-align:right;"> -0.1491895 </td>
   <td style="text-align:right;"> -0.1752869 </td>
   <td style="text-align:right;"> -0.3429842 </td>
   <td style="text-align:right;"> 0.7075767 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2648826 </td>
   <td style="text-align:right;"> 0.2154943 </td>
   <td style="text-align:right;"> 0.4766351 </td>
   <td style="text-align:right;"> 0.0857621 </td>
   <td style="text-align:right;"> -0.1844332 </td>
   <td style="text-align:right;"> 0.1405545 </td>
   <td style="text-align:right;"> 0.1762827 </td>
   <td style="text-align:right;"> 0.3971086 </td>
   <td style="text-align:right;"> 0.0661701 </td>
   <td style="text-align:right;"> -0.4042761 </td>
   <td style="text-align:right;"> -0.0117582 </td>
   <td style="text-align:right;"> 0.0440284 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -0.6964293 </td>
   <td style="text-align:right;"> -0.0006991 </td>
   <td style="text-align:right;"> 0.3792686 </td>
   <td style="text-align:right;"> -0.0052876 </td>
   <td style="text-align:right;"> 0.0742457 </td>
   <td style="text-align:right;"> -0.1252244 </td>
   <td style="text-align:right;"> -0.1176711 </td>
   <td style="text-align:right;"> -0.2513375 </td>
   <td style="text-align:right;"> 0.2662364 </td>
   <td style="text-align:right;"> -0.2346563 </td>
   <td style="text-align:right;"> -0.4376760 </td>
   <td style="text-align:right;"> 0.5148074 </td>
   <td style="text-align:right;"> 0.5986834 </td>
   <td style="text-align:right;"> 0.3860824 </td>
   <td style="text-align:right;"> -0.1928564 </td>
   <td style="text-align:right;"> 0.7081633 </td>
   <td style="text-align:right;"> 0.7495221 </td>
   <td style="text-align:right;"> 0.8989520 </td>
   <td style="text-align:right;"> 0.3068174 </td>
   <td style="text-align:right;"> 0.5821295 </td>
   <td style="text-align:right;"> 0.4929873 </td>
   <td style="text-align:right;"> 0.5631244 </td>
   <td style="text-align:right;"> 0.6298970 </td>
   <td style="text-align:right;"> 0.3511412 </td>
   <td style="text-align:right;"> 0.3964797 </td>
   <td style="text-align:right;"> 0.9145146 </td>
   <td style="text-align:right;"> 0.6109257 </td>
   <td style="text-align:right;"> 0.3368321 </td>
   <td style="text-align:right;"> 0.7063244 </td>
   <td style="text-align:right;"> 0.2813408 </td>
   <td style="text-align:right;"> 0.7258697 </td>
   <td style="text-align:right;"> 0.7102881 </td>
   <td style="text-align:right;"> 0.6332613 </td>
   <td style="text-align:right;"> 0.4812212 </td>
   <td style="text-align:right;"> 0.7203265 </td>
   <td style="text-align:right;"> 0.2470535 </td>
   <td style="text-align:right;"> 0.8126103 </td>
   <td style="text-align:right;"> 0.4702436 </td>
   <td style="text-align:right;"> 0.4508250 </td>
   <td style="text-align:right;"> 0.6578923 </td>
   <td style="text-align:right;"> 0.2648826 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.2575785 </td>
   <td style="text-align:right;"> 0.6368137 </td>
   <td style="text-align:right;"> 0.4502532 </td>
   <td style="text-align:right;"> 0.1978331 </td>
   <td style="text-align:right;"> 0.4615546 </td>
   <td style="text-align:right;"> 0.7589218 </td>
   <td style="text-align:right;"> 0.7856971 </td>
   <td style="text-align:right;"> 0.9012553 </td>
   <td style="text-align:right;"> -0.5245674 </td>
   <td style="text-align:right;"> -0.4224012 </td>
   <td style="text-align:right;"> 0.8600432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> -0.5499636 </td>
   <td style="text-align:right;"> 0.2488211 </td>
   <td style="text-align:right;"> 0.3632190 </td>
   <td style="text-align:right;"> -0.3282950 </td>
   <td style="text-align:right;"> -0.2055219 </td>
   <td style="text-align:right;"> -0.2364365 </td>
   <td style="text-align:right;"> -0.8482983 </td>
   <td style="text-align:right;"> -0.2206988 </td>
   <td style="text-align:right;"> 0.3002933 </td>
   <td style="text-align:right;"> -0.9008284 </td>
   <td style="text-align:right;"> 0.1975705 </td>
   <td style="text-align:right;"> 0.3108288 </td>
   <td style="text-align:right;"> 0.4577798 </td>
   <td style="text-align:right;"> 0.1661860 </td>
   <td style="text-align:right;"> -0.4105992 </td>
   <td style="text-align:right;"> 0.1785471 </td>
   <td style="text-align:right;"> 0.3825978 </td>
   <td style="text-align:right;"> 0.1997049 </td>
   <td style="text-align:right;"> 0.2012880 </td>
   <td style="text-align:right;"> 0.5400491 </td>
   <td style="text-align:right;"> 0.0675309 </td>
   <td style="text-align:right;"> 0.2680727 </td>
   <td style="text-align:right;"> -0.0227967 </td>
   <td style="text-align:right;"> -0.0259897 </td>
   <td style="text-align:right;"> 0.4348529 </td>
   <td style="text-align:right;"> 0.0654532 </td>
   <td style="text-align:right;"> 0.2375350 </td>
   <td style="text-align:right;"> 0.1034316 </td>
   <td style="text-align:right;"> 0.6217246 </td>
   <td style="text-align:right;"> 0.1877935 </td>
   <td style="text-align:right;"> 0.5517272 </td>
   <td style="text-align:right;"> 0.6493740 </td>
   <td style="text-align:right;"> -0.1875635 </td>
   <td style="text-align:right;"> 0.3570380 </td>
   <td style="text-align:right;"> 0.4241936 </td>
   <td style="text-align:right;"> 0.5493776 </td>
   <td style="text-align:right;"> 0.3134236 </td>
   <td style="text-align:right;"> 0.6398505 </td>
   <td style="text-align:right;"> 0.0267657 </td>
   <td style="text-align:right;"> 0.1371019 </td>
   <td style="text-align:right;"> 0.2154943 </td>
   <td style="text-align:right;"> 0.2575785 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.6640750 </td>
   <td style="text-align:right;"> 0.7693593 </td>
   <td style="text-align:right;"> 0.5484572 </td>
   <td style="text-align:right;"> 0.5833896 </td>
   <td style="text-align:right;"> 0.4390556 </td>
   <td style="text-align:right;"> 0.4117307 </td>
   <td style="text-align:right;"> 0.2266379 </td>
   <td style="text-align:right;"> -0.5384446 </td>
   <td style="text-align:right;"> -0.5039197 </td>
   <td style="text-align:right;"> 0.2911579 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.7051383 </td>
   <td style="text-align:right;"> 0.0228269 </td>
   <td style="text-align:right;"> 0.1998870 </td>
   <td style="text-align:right;"> -0.2596329 </td>
   <td style="text-align:right;"> -0.1385398 </td>
   <td style="text-align:right;"> -0.2653618 </td>
   <td style="text-align:right;"> -0.3368757 </td>
   <td style="text-align:right;"> -0.2791451 </td>
   <td style="text-align:right;"> 0.3078119 </td>
   <td style="text-align:right;"> -0.4513308 </td>
   <td style="text-align:right;"> -0.1371942 </td>
   <td style="text-align:right;"> 0.5611951 </td>
   <td style="text-align:right;"> 0.7515435 </td>
   <td style="text-align:right;"> 0.6558529 </td>
   <td style="text-align:right;"> -0.3347895 </td>
   <td style="text-align:right;"> 0.6681524 </td>
   <td style="text-align:right;"> 0.4436520 </td>
   <td style="text-align:right;"> 0.5586852 </td>
   <td style="text-align:right;"> 0.0750472 </td>
   <td style="text-align:right;"> 0.5088258 </td>
   <td style="text-align:right;"> 0.4151545 </td>
   <td style="text-align:right;"> 0.2799497 </td>
   <td style="text-align:right;"> 0.1120483 </td>
   <td style="text-align:right;"> 0.2975372 </td>
   <td style="text-align:right;"> 0.5459679 </td>
   <td style="text-align:right;"> 0.5283532 </td>
   <td style="text-align:right;"> 0.5254704 </td>
   <td style="text-align:right;"> 0.4521939 </td>
   <td style="text-align:right;"> 0.7289069 </td>
   <td style="text-align:right;"> 0.5571466 </td>
   <td style="text-align:right;"> 0.8747534 </td>
   <td style="text-align:right;"> 0.7312017 </td>
   <td style="text-align:right;"> 0.3956055 </td>
   <td style="text-align:right;"> 0.4250268 </td>
   <td style="text-align:right;"> 0.5633014 </td>
   <td style="text-align:right;"> 0.4546046 </td>
   <td style="text-align:right;"> 0.3979449 </td>
   <td style="text-align:right;"> 0.5477372 </td>
   <td style="text-align:right;"> 0.1787624 </td>
   <td style="text-align:right;"> 0.5938416 </td>
   <td style="text-align:right;"> 0.4766351 </td>
   <td style="text-align:right;"> 0.6368137 </td>
   <td style="text-align:right;"> 0.6640750 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.7229519 </td>
   <td style="text-align:right;"> 0.3773718 </td>
   <td style="text-align:right;"> 0.5342736 </td>
   <td style="text-align:right;"> 0.7006260 </td>
   <td style="text-align:right;"> 0.7253313 </td>
   <td style="text-align:right;"> 0.4127562 </td>
   <td style="text-align:right;"> -0.5830683 </td>
   <td style="text-align:right;"> -0.1659296 </td>
   <td style="text-align:right;"> 0.5468311 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -0.5467648 </td>
   <td style="text-align:right;"> 0.3595014 </td>
   <td style="text-align:right;"> 0.0663445 </td>
   <td style="text-align:right;"> -0.0025225 </td>
   <td style="text-align:right;"> -0.2857829 </td>
   <td style="text-align:right;"> -0.2244579 </td>
   <td style="text-align:right;"> -0.4947810 </td>
   <td style="text-align:right;"> 0.0023138 </td>
   <td style="text-align:right;"> 0.5142097 </td>
   <td style="text-align:right;"> -0.5790447 </td>
   <td style="text-align:right;"> 0.1242043 </td>
   <td style="text-align:right;"> 0.2051371 </td>
   <td style="text-align:right;"> 0.4459148 </td>
   <td style="text-align:right;"> 0.3172330 </td>
   <td style="text-align:right;"> -0.3669278 </td>
   <td style="text-align:right;"> 0.6712421 </td>
   <td style="text-align:right;"> 0.6138754 </td>
   <td style="text-align:right;"> 0.4717185 </td>
   <td style="text-align:right;"> 0.3547808 </td>
   <td style="text-align:right;"> 0.5017794 </td>
   <td style="text-align:right;"> 0.3334066 </td>
   <td style="text-align:right;"> 0.5553014 </td>
   <td style="text-align:right;"> 0.4110367 </td>
   <td style="text-align:right;"> 0.2685467 </td>
   <td style="text-align:right;"> 0.3843371 </td>
   <td style="text-align:right;"> 0.2372470 </td>
   <td style="text-align:right;"> 0.4954871 </td>
   <td style="text-align:right;"> 0.3549600 </td>
   <td style="text-align:right;"> 0.8491263 </td>
   <td style="text-align:right;"> 0.1035302 </td>
   <td style="text-align:right;"> 0.6478655 </td>
   <td style="text-align:right;"> 0.4751019 </td>
   <td style="text-align:right;"> 0.1668547 </td>
   <td style="text-align:right;"> 0.7271529 </td>
   <td style="text-align:right;"> 0.6186419 </td>
   <td style="text-align:right;"> 0.3912704 </td>
   <td style="text-align:right;"> 0.5820092 </td>
   <td style="text-align:right;"> 0.8084788 </td>
   <td style="text-align:right;"> 0.2911605 </td>
   <td style="text-align:right;"> 0.0964395 </td>
   <td style="text-align:right;"> 0.0857621 </td>
   <td style="text-align:right;"> 0.4502532 </td>
   <td style="text-align:right;"> 0.7693593 </td>
   <td style="text-align:right;"> 0.7229519 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3720631 </td>
   <td style="text-align:right;"> 0.7992582 </td>
   <td style="text-align:right;"> 0.6415606 </td>
   <td style="text-align:right;"> 0.4952203 </td>
   <td style="text-align:right;"> 0.4592053 </td>
   <td style="text-align:right;"> -0.5565559 </td>
   <td style="text-align:right;"> -0.5824528 </td>
   <td style="text-align:right;"> 0.6438321 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> -0.4553864 </td>
   <td style="text-align:right;"> -0.1659305 </td>
   <td style="text-align:right;"> -0.0089270 </td>
   <td style="text-align:right;"> -0.4294680 </td>
   <td style="text-align:right;"> 0.2944537 </td>
   <td style="text-align:right;"> 0.2081899 </td>
   <td style="text-align:right;"> -0.4874364 </td>
   <td style="text-align:right;"> -0.2150772 </td>
   <td style="text-align:right;"> -0.0111774 </td>
   <td style="text-align:right;"> -0.4927250 </td>
   <td style="text-align:right;"> 0.0776946 </td>
   <td style="text-align:right;"> 0.6076428 </td>
   <td style="text-align:right;"> 0.6385606 </td>
   <td style="text-align:right;"> 0.4125514 </td>
   <td style="text-align:right;"> -0.6577265 </td>
   <td style="text-align:right;"> 0.1979961 </td>
   <td style="text-align:right;"> 0.2370297 </td>
   <td style="text-align:right;"> -0.0039871 </td>
   <td style="text-align:right;"> -0.0514128 </td>
   <td style="text-align:right;"> 0.8090469 </td>
   <td style="text-align:right;"> 0.3157012 </td>
   <td style="text-align:right;"> 0.0746563 </td>
   <td style="text-align:right;"> 0.0417548 </td>
   <td style="text-align:right;"> -0.3643828 </td>
   <td style="text-align:right;"> 0.2213689 </td>
   <td style="text-align:right;"> 0.0862332 </td>
   <td style="text-align:right;"> 0.6254837 </td>
   <td style="text-align:right;"> 0.1245253 </td>
   <td style="text-align:right;"> 0.4589792 </td>
   <td style="text-align:right;"> 0.4132668 </td>
   <td style="text-align:right;"> 0.3910145 </td>
   <td style="text-align:right;"> 0.3065580 </td>
   <td style="text-align:right;"> -0.0197972 </td>
   <td style="text-align:right;"> 0.2142938 </td>
   <td style="text-align:right;"> 0.5485451 </td>
   <td style="text-align:right;"> 0.9826091 </td>
   <td style="text-align:right;"> 0.1616629 </td>
   <td style="text-align:right;"> 0.6261250 </td>
   <td style="text-align:right;"> 0.2908171 </td>
   <td style="text-align:right;"> -0.0990521 </td>
   <td style="text-align:right;"> -0.1844332 </td>
   <td style="text-align:right;"> 0.1978331 </td>
   <td style="text-align:right;"> 0.5484572 </td>
   <td style="text-align:right;"> 0.3773718 </td>
   <td style="text-align:right;"> 0.3720631 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.1556419 </td>
   <td style="text-align:right;"> 0.6299097 </td>
   <td style="text-align:right;"> 0.4046092 </td>
   <td style="text-align:right;"> 0.0730052 </td>
   <td style="text-align:right;"> -0.0724221 </td>
   <td style="text-align:right;"> -0.1400568 </td>
   <td style="text-align:right;"> 0.0681694 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -0.5117469 </td>
   <td style="text-align:right;"> 0.3662092 </td>
   <td style="text-align:right;"> 0.0896694 </td>
   <td style="text-align:right;"> 0.2803789 </td>
   <td style="text-align:right;"> -0.0612862 </td>
   <td style="text-align:right;"> 0.0426049 </td>
   <td style="text-align:right;"> -0.3927382 </td>
   <td style="text-align:right;"> 0.1091440 </td>
   <td style="text-align:right;"> 0.7067664 </td>
   <td style="text-align:right;"> -0.5088037 </td>
   <td style="text-align:right;"> -0.1519739 </td>
   <td style="text-align:right;"> 0.0874129 </td>
   <td style="text-align:right;"> 0.2123547 </td>
   <td style="text-align:right;"> -0.0410795 </td>
   <td style="text-align:right;"> -0.2544973 </td>
   <td style="text-align:right;"> 0.4971349 </td>
   <td style="text-align:right;"> 0.5995786 </td>
   <td style="text-align:right;"> 0.5452433 </td>
   <td style="text-align:right;"> 0.6863451 </td>
   <td style="text-align:right;"> 0.3272334 </td>
   <td style="text-align:right;"> 0.4755774 </td>
   <td style="text-align:right;"> 0.6153451 </td>
   <td style="text-align:right;"> 0.3618197 </td>
   <td style="text-align:right;"> 0.2372717 </td>
   <td style="text-align:right;"> 0.5857368 </td>
   <td style="text-align:right;"> 0.3226870 </td>
   <td style="text-align:right;"> 0.4127620 </td>
   <td style="text-align:right;"> 0.6215820 </td>
   <td style="text-align:right;"> 0.7154843 </td>
   <td style="text-align:right;"> -0.1432812 </td>
   <td style="text-align:right;"> 0.4358423 </td>
   <td style="text-align:right;"> 0.5019571 </td>
   <td style="text-align:right;"> 0.0784152 </td>
   <td style="text-align:right;"> 0.7158626 </td>
   <td style="text-align:right;"> 0.6595585 </td>
   <td style="text-align:right;"> 0.1433290 </td>
   <td style="text-align:right;"> 0.6620807 </td>
   <td style="text-align:right;"> 0.5043378 </td>
   <td style="text-align:right;"> 0.5864228 </td>
   <td style="text-align:right;"> 0.2808560 </td>
   <td style="text-align:right;"> 0.1405545 </td>
   <td style="text-align:right;"> 0.4615546 </td>
   <td style="text-align:right;"> 0.5833896 </td>
   <td style="text-align:right;"> 0.5342736 </td>
   <td style="text-align:right;"> 0.7992582 </td>
   <td style="text-align:right;"> 0.1556419 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.4223671 </td>
   <td style="text-align:right;"> 0.5771987 </td>
   <td style="text-align:right;"> 0.5763982 </td>
   <td style="text-align:right;"> -0.6559887 </td>
   <td style="text-align:right;"> -0.6834933 </td>
   <td style="text-align:right;"> 0.7464062 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> -0.5739959 </td>
   <td style="text-align:right;"> 0.0393083 </td>
   <td style="text-align:right;"> 0.0395150 </td>
   <td style="text-align:right;"> -0.2168310 </td>
   <td style="text-align:right;"> 0.0517411 </td>
   <td style="text-align:right;"> -0.1131923 </td>
   <td style="text-align:right;"> -0.2725973 </td>
   <td style="text-align:right;"> -0.3603460 </td>
   <td style="text-align:right;"> 0.1859939 </td>
   <td style="text-align:right;"> -0.3817900 </td>
   <td style="text-align:right;"> -0.1673150 </td>
   <td style="text-align:right;"> 0.5512961 </td>
   <td style="text-align:right;"> 0.8763335 </td>
   <td style="text-align:right;"> 0.7180336 </td>
   <td style="text-align:right;"> -0.4589151 </td>
   <td style="text-align:right;"> 0.7533810 </td>
   <td style="text-align:right;"> 0.5839850 </td>
   <td style="text-align:right;"> 0.5587800 </td>
   <td style="text-align:right;"> 0.1387526 </td>
   <td style="text-align:right;"> 0.7719014 </td>
   <td style="text-align:right;"> 0.6556621 </td>
   <td style="text-align:right;"> 0.5665743 </td>
   <td style="text-align:right;"> 0.6119762 </td>
   <td style="text-align:right;"> -0.0181701 </td>
   <td style="text-align:right;"> 0.3929039 </td>
   <td style="text-align:right;"> 0.6861475 </td>
   <td style="text-align:right;"> 0.8226371 </td>
   <td style="text-align:right;"> 0.2266872 </td>
   <td style="text-align:right;"> 0.8724481 </td>
   <td style="text-align:right;"> 0.5747898 </td>
   <td style="text-align:right;"> 0.8131381 </td>
   <td style="text-align:right;"> 0.5800730 </td>
   <td style="text-align:right;"> 0.6163482 </td>
   <td style="text-align:right;"> 0.5282233 </td>
   <td style="text-align:right;"> 0.7293122 </td>
   <td style="text-align:right;"> 0.6942590 </td>
   <td style="text-align:right;"> 0.5591358 </td>
   <td style="text-align:right;"> 0.6844641 </td>
   <td style="text-align:right;"> 0.3518794 </td>
   <td style="text-align:right;"> 0.3627018 </td>
   <td style="text-align:right;"> 0.1762827 </td>
   <td style="text-align:right;"> 0.7589218 </td>
   <td style="text-align:right;"> 0.4390556 </td>
   <td style="text-align:right;"> 0.7006260 </td>
   <td style="text-align:right;"> 0.6415606 </td>
   <td style="text-align:right;"> 0.6299097 </td>
   <td style="text-align:right;"> 0.4223671 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.8068349 </td>
   <td style="text-align:right;"> 0.6105340 </td>
   <td style="text-align:right;"> -0.4036570 </td>
   <td style="text-align:right;"> -0.3993933 </td>
   <td style="text-align:right;"> 0.6102114 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -0.6761977 </td>
   <td style="text-align:right;"> 0.0647647 </td>
   <td style="text-align:right;"> 0.1729750 </td>
   <td style="text-align:right;"> -0.0689355 </td>
   <td style="text-align:right;"> 0.0357834 </td>
   <td style="text-align:right;"> -0.1173486 </td>
   <td style="text-align:right;"> -0.2588695 </td>
   <td style="text-align:right;"> -0.5131760 </td>
   <td style="text-align:right;"> 0.2154662 </td>
   <td style="text-align:right;"> -0.4135695 </td>
   <td style="text-align:right;"> -0.4853268 </td>
   <td style="text-align:right;"> 0.5589009 </td>
   <td style="text-align:right;"> 0.7575520 </td>
   <td style="text-align:right;"> 0.5349442 </td>
   <td style="text-align:right;"> -0.4692880 </td>
   <td style="text-align:right;"> 0.5305691 </td>
   <td style="text-align:right;"> 0.5034716 </td>
   <td style="text-align:right;"> 0.5975853 </td>
   <td style="text-align:right;"> 0.3824594 </td>
   <td style="text-align:right;"> 0.5069060 </td>
   <td style="text-align:right;"> 0.6926927 </td>
   <td style="text-align:right;"> 0.6355951 </td>
   <td style="text-align:right;"> 0.4566401 </td>
   <td style="text-align:right;"> 0.0373242 </td>
   <td style="text-align:right;"> 0.5433404 </td>
   <td style="text-align:right;"> 0.7550992 </td>
   <td style="text-align:right;"> 0.6558231 </td>
   <td style="text-align:right;"> 0.4350536 </td>
   <td style="text-align:right;"> 0.7189703 </td>
   <td style="text-align:right;"> 0.4613204 </td>
   <td style="text-align:right;"> 0.6973238 </td>
   <td style="text-align:right;"> 0.7117702 </td>
   <td style="text-align:right;"> 0.6310694 </td>
   <td style="text-align:right;"> 0.3625363 </td>
   <td style="text-align:right;"> 0.6098670 </td>
   <td style="text-align:right;"> 0.4766386 </td>
   <td style="text-align:right;"> 0.5167963 </td>
   <td style="text-align:right;"> 0.3391044 </td>
   <td style="text-align:right;"> 0.4255995 </td>
   <td style="text-align:right;"> 0.7272263 </td>
   <td style="text-align:right;"> 0.3971086 </td>
   <td style="text-align:right;"> 0.7856971 </td>
   <td style="text-align:right;"> 0.4117307 </td>
   <td style="text-align:right;"> 0.7253313 </td>
   <td style="text-align:right;"> 0.4952203 </td>
   <td style="text-align:right;"> 0.4046092 </td>
   <td style="text-align:right;"> 0.5771987 </td>
   <td style="text-align:right;"> 0.8068349 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.6441088 </td>
   <td style="text-align:right;"> -0.5656363 </td>
   <td style="text-align:right;"> -0.4231725 </td>
   <td style="text-align:right;"> 0.6533322 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> -0.5786772 </td>
   <td style="text-align:right;"> 0.1457133 </td>
   <td style="text-align:right;"> 0.2423324 </td>
   <td style="text-align:right;"> 0.1344105 </td>
   <td style="text-align:right;"> 0.0641090 </td>
   <td style="text-align:right;"> -0.0723769 </td>
   <td style="text-align:right;"> -0.1433470 </td>
   <td style="text-align:right;"> -0.1967517 </td>
   <td style="text-align:right;"> 0.4102892 </td>
   <td style="text-align:right;"> -0.2810713 </td>
   <td style="text-align:right;"> -0.5084465 </td>
   <td style="text-align:right;"> 0.1721538 </td>
   <td style="text-align:right;"> 0.3730915 </td>
   <td style="text-align:right;"> 0.0616950 </td>
   <td style="text-align:right;"> -0.0466410 </td>
   <td style="text-align:right;"> 0.5980351 </td>
   <td style="text-align:right;"> 0.8332376 </td>
   <td style="text-align:right;"> 0.9054749 </td>
   <td style="text-align:right;"> 0.5171061 </td>
   <td style="text-align:right;"> 0.4824979 </td>
   <td style="text-align:right;"> 0.5027962 </td>
   <td style="text-align:right;"> 0.5908747 </td>
   <td style="text-align:right;"> 0.7000831 </td>
   <td style="text-align:right;"> 0.2464706 </td>
   <td style="text-align:right;"> 0.3945620 </td>
   <td style="text-align:right;"> 0.7804699 </td>
   <td style="text-align:right;"> 0.4120752 </td>
   <td style="text-align:right;"> 0.3261663 </td>
   <td style="text-align:right;"> 0.6622139 </td>
   <td style="text-align:right;"> -0.0201596 </td>
   <td style="text-align:right;"> 0.5581929 </td>
   <td style="text-align:right;"> 0.6092980 </td>
   <td style="text-align:right;"> 0.4317252 </td>
   <td style="text-align:right;"> 0.4572592 </td>
   <td style="text-align:right;"> 0.7369762 </td>
   <td style="text-align:right;"> 0.1131029 </td>
   <td style="text-align:right;"> 0.9169482 </td>
   <td style="text-align:right;"> 0.4257947 </td>
   <td style="text-align:right;"> 0.5509849 </td>
   <td style="text-align:right;"> 0.4843486 </td>
   <td style="text-align:right;"> 0.0661701 </td>
   <td style="text-align:right;"> 0.9012553 </td>
   <td style="text-align:right;"> 0.2266379 </td>
   <td style="text-align:right;"> 0.4127562 </td>
   <td style="text-align:right;"> 0.4592053 </td>
   <td style="text-align:right;"> 0.0730052 </td>
   <td style="text-align:right;"> 0.5763982 </td>
   <td style="text-align:right;"> 0.6105340 </td>
   <td style="text-align:right;"> 0.6441088 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.5222243 </td>
   <td style="text-align:right;"> -0.6705987 </td>
   <td style="text-align:right;"> 0.8849769 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 0.5551287 </td>
   <td style="text-align:right;"> -0.6841765 </td>
   <td style="text-align:right;"> -0.3563901 </td>
   <td style="text-align:right;"> -0.3800747 </td>
   <td style="text-align:right;"> -0.2431921 </td>
   <td style="text-align:right;"> -0.1965567 </td>
   <td style="text-align:right;"> 0.2165828 </td>
   <td style="text-align:right;"> 0.1782288 </td>
   <td style="text-align:right;"> -0.7690763 </td>
   <td style="text-align:right;"> 0.4334539 </td>
   <td style="text-align:right;"> -0.0242666 </td>
   <td style="text-align:right;"> -0.1939736 </td>
   <td style="text-align:right;"> -0.4496292 </td>
   <td style="text-align:right;"> -0.0046871 </td>
   <td style="text-align:right;"> -0.1409265 </td>
   <td style="text-align:right;"> -0.3982757 </td>
   <td style="text-align:right;"> -0.6684651 </td>
   <td style="text-align:right;"> -0.6607967 </td>
   <td style="text-align:right;"> -0.0576496 </td>
   <td style="text-align:right;"> -0.4294212 </td>
   <td style="text-align:right;"> -0.4927414 </td>
   <td style="text-align:right;"> -0.4547195 </td>
   <td style="text-align:right;"> -0.1683674 </td>
   <td style="text-align:right;"> -0.0762704 </td>
   <td style="text-align:right;"> -0.5786438 </td>
   <td style="text-align:right;"> -0.4128823 </td>
   <td style="text-align:right;"> -0.3455088 </td>
   <td style="text-align:right;"> -0.1615376 </td>
   <td style="text-align:right;"> -0.5543335 </td>
   <td style="text-align:right;"> -0.2457316 </td>
   <td style="text-align:right;"> -0.6671117 </td>
   <td style="text-align:right;"> -0.7911389 </td>
   <td style="text-align:right;"> -0.2952387 </td>
   <td style="text-align:right;"> -0.4522680 </td>
   <td style="text-align:right;"> -0.4924795 </td>
   <td style="text-align:right;"> -0.1043352 </td>
   <td style="text-align:right;"> -0.5767370 </td>
   <td style="text-align:right;"> -0.1802631 </td>
   <td style="text-align:right;"> -0.4017788 </td>
   <td style="text-align:right;"> -0.3420923 </td>
   <td style="text-align:right;"> -0.4042761 </td>
   <td style="text-align:right;"> -0.5245674 </td>
   <td style="text-align:right;"> -0.5384446 </td>
   <td style="text-align:right;"> -0.5830683 </td>
   <td style="text-align:right;"> -0.5565559 </td>
   <td style="text-align:right;"> -0.0724221 </td>
   <td style="text-align:right;"> -0.6559887 </td>
   <td style="text-align:right;"> -0.4036570 </td>
   <td style="text-align:right;"> -0.5656363 </td>
   <td style="text-align:right;"> -0.5222243 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3678249 </td>
   <td style="text-align:right;"> -0.6307091 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.2306563 </td>
   <td style="text-align:right;"> -0.3090736 </td>
   <td style="text-align:right;"> -0.1558881 </td>
   <td style="text-align:right;"> 0.0122501 </td>
   <td style="text-align:right;"> 0.2751016 </td>
   <td style="text-align:right;"> 0.2813786 </td>
   <td style="text-align:right;"> 0.6138286 </td>
   <td style="text-align:right;"> 0.2026866 </td>
   <td style="text-align:right;"> -0.3032435 </td>
   <td style="text-align:right;"> 0.6838157 </td>
   <td style="text-align:right;"> 0.2045341 </td>
   <td style="text-align:right;"> 0.1952591 </td>
   <td style="text-align:right;"> -0.1352114 </td>
   <td style="text-align:right;"> 0.1623411 </td>
   <td style="text-align:right;"> 0.1967519 </td>
   <td style="text-align:right;"> -0.1852456 </td>
   <td style="text-align:right;"> -0.5166206 </td>
   <td style="text-align:right;"> -0.3760829 </td>
   <td style="text-align:right;"> -0.7536652 </td>
   <td style="text-align:right;"> -0.2841066 </td>
   <td style="text-align:right;"> -0.2652725 </td>
   <td style="text-align:right;"> -0.6732367 </td>
   <td style="text-align:right;"> -0.6052790 </td>
   <td style="text-align:right;"> 0.0943284 </td>
   <td style="text-align:right;"> -0.3473329 </td>
   <td style="text-align:right;"> -0.3199132 </td>
   <td style="text-align:right;"> -0.1209037 </td>
   <td style="text-align:right;"> -0.0799177 </td>
   <td style="text-align:right;"> -0.6347501 </td>
   <td style="text-align:right;"> 0.2534905 </td>
   <td style="text-align:right;"> -0.2484420 </td>
   <td style="text-align:right;"> -0.4073193 </td>
   <td style="text-align:right;"> 0.0171207 </td>
   <td style="text-align:right;"> -0.4096391 </td>
   <td style="text-align:right;"> -0.4530231 </td>
   <td style="text-align:right;"> -0.1437649 </td>
   <td style="text-align:right;"> -0.6042877 </td>
   <td style="text-align:right;"> -0.4545788 </td>
   <td style="text-align:right;"> -0.2316881 </td>
   <td style="text-align:right;"> -0.1552865 </td>
   <td style="text-align:right;"> -0.0117582 </td>
   <td style="text-align:right;"> -0.4224012 </td>
   <td style="text-align:right;"> -0.5039197 </td>
   <td style="text-align:right;"> -0.1659296 </td>
   <td style="text-align:right;"> -0.5824528 </td>
   <td style="text-align:right;"> -0.1400568 </td>
   <td style="text-align:right;"> -0.6834933 </td>
   <td style="text-align:right;"> -0.3993933 </td>
   <td style="text-align:right;"> -0.4231725 </td>
   <td style="text-align:right;"> -0.6705987 </td>
   <td style="text-align:right;"> 0.3678249 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> -0.5371928 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -0.7017066 </td>
   <td style="text-align:right;"> 0.2794976 </td>
   <td style="text-align:right;"> 0.2739832 </td>
   <td style="text-align:right;"> 0.3695822 </td>
   <td style="text-align:right;"> 0.0401972 </td>
   <td style="text-align:right;"> 0.0282879 </td>
   <td style="text-align:right;"> -0.0548690 </td>
   <td style="text-align:right;"> 0.0483113 </td>
   <td style="text-align:right;"> 0.5629398 </td>
   <td style="text-align:right;"> -0.1798006 </td>
   <td style="text-align:right;"> -0.3159206 </td>
   <td style="text-align:right;"> 0.3317746 </td>
   <td style="text-align:right;"> 0.3191727 </td>
   <td style="text-align:right;"> 0.1074513 </td>
   <td style="text-align:right;"> -0.1741065 </td>
   <td style="text-align:right;"> 0.7621499 </td>
   <td style="text-align:right;"> 0.8825878 </td>
   <td style="text-align:right;"> 0.9194690 </td>
   <td style="text-align:right;"> 0.4812776 </td>
   <td style="text-align:right;"> 0.4670299 </td>
   <td style="text-align:right;"> 0.4725342 </td>
   <td style="text-align:right;"> 0.6968867 </td>
   <td style="text-align:right;"> 0.6842991 </td>
   <td style="text-align:right;"> 0.5219961 </td>
   <td style="text-align:right;"> 0.3354050 </td>
   <td style="text-align:right;"> 0.6821804 </td>
   <td style="text-align:right;"> 0.5947196 </td>
   <td style="text-align:right;"> 0.5011157 </td>
   <td style="text-align:right;"> 0.7051819 </td>
   <td style="text-align:right;"> -0.0701946 </td>
   <td style="text-align:right;"> 0.5668516 </td>
   <td style="text-align:right;"> 0.5412015 </td>
   <td style="text-align:right;"> 0.5112659 </td>
   <td style="text-align:right;"> 0.7181758 </td>
   <td style="text-align:right;"> 0.7560358 </td>
   <td style="text-align:right;"> 0.0755796 </td>
   <td style="text-align:right;"> 0.9429335 </td>
   <td style="text-align:right;"> 0.5387702 </td>
   <td style="text-align:right;"> 0.6655955 </td>
   <td style="text-align:right;"> 0.4060250 </td>
   <td style="text-align:right;"> 0.0440284 </td>
   <td style="text-align:right;"> 0.8600432 </td>
   <td style="text-align:right;"> 0.2911579 </td>
   <td style="text-align:right;"> 0.5468311 </td>
   <td style="text-align:right;"> 0.6438321 </td>
   <td style="text-align:right;"> 0.0681694 </td>
   <td style="text-align:right;"> 0.7464062 </td>
   <td style="text-align:right;"> 0.6102114 </td>
   <td style="text-align:right;"> 0.6533322 </td>
   <td style="text-align:right;"> 0.8849769 </td>
   <td style="text-align:right;"> -0.6307091 </td>
   <td style="text-align:right;"> -0.5371928 </td>
   <td style="text-align:right;"> 1.0000000 </td>
  </tr>
</tbody>
</table></div>

<br>  
<br>  



## Supplemental Table 5. Linear regression statistics for HDL protein logFC and stroke recovery {.tabset}

All linear models assessing protein log2FC following stroke and their relation to stroke recovery are recorded here. Models are tab separated into 3 main categories: NIHSS models - these assess protein relation to NIHSS scores at 3 months; mRS models - these assess protein relation to mRS scores at 3 months; and Age and Sex models - these are additional parameters examined in the review process. 

### NIHSS models



__Model: NIHSS_3mo ~ NIHSS_baseline + tPA + Protein Log2FC__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + tPa + i, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_NIH3_baseline.tPa.Prot <- out_df
rownames(lm_NIH3_baseline.tPa.Prot) <- rn[,1]
lm_NIH3_baseline.tPa.Prot$p.adj <- p.adjust(lm_NIH3_baseline.tPa.Prot$p, 
                                            method = "BH")

kable(lm_NIH3_baseline.tPa.Prot) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> -1.2130203 </td>
   <td style="text-align:right;"> 5.916280 </td>
   <td style="text-align:right;"> 0.2096030 </td>
   <td style="text-align:right;"> 0.1017329 </td>
   <td style="text-align:right;"> 2.502726 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> -0.0080804 </td>
   <td style="text-align:right;"> 6.282966 </td>
   <td style="text-align:right;"> 0.1612593 </td>
   <td style="text-align:right;"> 0.1476270 </td>
   <td style="text-align:right;"> 2.089494 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.8035443 </td>
   <td style="text-align:right;"> 5.999027 </td>
   <td style="text-align:right;"> 0.2316193 </td>
   <td style="text-align:right;"> 0.0850359 </td>
   <td style="text-align:right;"> 2.708150 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -2.0631617 </td>
   <td style="text-align:right;"> 6.276763 </td>
   <td style="text-align:right;"> 0.1575452 </td>
   <td style="text-align:right;"> 0.1517343 </td>
   <td style="text-align:right;"> 2.059708 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 2.9933515 </td>
   <td style="text-align:right;"> 5.424179 </td>
   <td style="text-align:right;"> 0.3967288 </td>
   <td style="text-align:right;"> 0.0176284 </td>
   <td style="text-align:right;"> 4.726565 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1116468 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -6.4993059 </td>
   <td style="text-align:right;"> 4.213293 </td>
   <td style="text-align:right;"> 0.7469416 </td>
   <td style="text-align:right;"> 0.0013438 </td>
   <td style="text-align:right;"> 12.806627 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.0209971 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> -2.4246885 </td>
   <td style="text-align:right;"> 6.445506 </td>
   <td style="text-align:right;"> 0.1561935 </td>
   <td style="text-align:right;"> 0.1532512 </td>
   <td style="text-align:right;"> 2.048933 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -2.6128163 </td>
   <td style="text-align:right;"> 6.358492 </td>
   <td style="text-align:right;"> 0.1685174 </td>
   <td style="text-align:right;"> 0.1398531 </td>
   <td style="text-align:right;"> 2.148469 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -1.1692406 </td>
   <td style="text-align:right;"> 6.025861 </td>
   <td style="text-align:right;"> 0.1799668 </td>
   <td style="text-align:right;"> 0.1282499 </td>
   <td style="text-align:right;"> 2.243623 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> -0.3846900 </td>
   <td style="text-align:right;"> 6.054828 </td>
   <td style="text-align:right;"> 0.1844325 </td>
   <td style="text-align:right;"> 0.1239367 </td>
   <td style="text-align:right;"> 2.281461 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 3.0997628 </td>
   <td style="text-align:right;"> 6.127625 </td>
   <td style="text-align:right;"> 0.2844073 </td>
   <td style="text-align:right;"> 0.0538647 </td>
   <td style="text-align:right;"> 3.252177 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1462042 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 9.2829523 </td>
   <td style="text-align:right;"> 6.128948 </td>
   <td style="text-align:right;"> 0.4504761 </td>
   <td style="text-align:right;"> 0.0094853 </td>
   <td style="text-align:right;"> 5.645289 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0901104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -1.2519457 </td>
   <td style="text-align:right;"> 6.120398 </td>
   <td style="text-align:right;"> 0.1552174 </td>
   <td style="text-align:right;"> 0.1543538 </td>
   <td style="text-align:right;"> 2.041174 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 2.1355272 </td>
   <td style="text-align:right;"> 5.552717 </td>
   <td style="text-align:right;"> 0.3527780 </td>
   <td style="text-align:right;"> 0.0280118 </td>
   <td style="text-align:right;"> 4.088701 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1182722 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -7.4112534 </td>
   <td style="text-align:right;"> 3.693201 </td>
   <td style="text-align:right;"> 0.7221328 </td>
   <td style="text-align:right;"> 0.0000923 </td>
   <td style="text-align:right;"> 15.726771 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0035079 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 1.2980299 </td>
   <td style="text-align:right;"> 4.377695 </td>
   <td style="text-align:right;"> 0.5759572 </td>
   <td style="text-align:right;"> 0.0016577 </td>
   <td style="text-align:right;"> 8.696764 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0209971 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -0.6916747 </td>
   <td style="text-align:right;"> 5.309327 </td>
   <td style="text-align:right;"> 0.3640139 </td>
   <td style="text-align:right;"> 0.0249695 </td>
   <td style="text-align:right;"> 4.243381 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1182722 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> -1.1150148 </td>
   <td style="text-align:right;"> 6.228018 </td>
   <td style="text-align:right;"> 0.1309215 </td>
   <td style="text-align:right;"> 0.1838574 </td>
   <td style="text-align:right;"> 1.853649 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1838574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> -0.5616944 </td>
   <td style="text-align:right;"> 6.105875 </td>
   <td style="text-align:right;"> 0.1682405 </td>
   <td style="text-align:right;"> 0.1401436 </td>
   <td style="text-align:right;"> 2.146200 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 1.6567386 </td>
   <td style="text-align:right;"> 5.917846 </td>
   <td style="text-align:right;"> 0.2720873 </td>
   <td style="text-align:right;"> 0.0601334 </td>
   <td style="text-align:right;"> 3.118150 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1480206 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.9675912 </td>
   <td style="text-align:right;"> 5.596336 </td>
   <td style="text-align:right;"> 0.3177269 </td>
   <td style="text-align:right;"> 0.0395378 </td>
   <td style="text-align:right;"> 3.638902 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1365850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -0.0816117 </td>
   <td style="text-align:right;"> 6.447519 </td>
   <td style="text-align:right;"> 0.1473384 </td>
   <td style="text-align:right;"> 0.1634841 </td>
   <td style="text-align:right;"> 1.979189 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1725665 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.1889865 </td>
   <td style="text-align:right;"> 6.457399 </td>
   <td style="text-align:right;"> 0.1543615 </td>
   <td style="text-align:right;"> 0.1553259 </td>
   <td style="text-align:right;"> 2.034384 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 1.3604303 </td>
   <td style="text-align:right;"> 6.641231 </td>
   <td style="text-align:right;"> 0.1771497 </td>
   <td style="text-align:right;"> 0.1310315 </td>
   <td style="text-align:right;"> 2.219965 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 0.4438287 </td>
   <td style="text-align:right;"> 6.039290 </td>
   <td style="text-align:right;"> 0.2129540 </td>
   <td style="text-align:right;"> 0.0990344 </td>
   <td style="text-align:right;"> 2.533252 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 2.5121178 </td>
   <td style="text-align:right;"> 8.103455 </td>
   <td style="text-align:right;"> 0.2075250 </td>
   <td style="text-align:right;"> 0.1151493 </td>
   <td style="text-align:right;"> 2.396637 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -0.0543307 </td>
   <td style="text-align:right;"> 5.574148 </td>
   <td style="text-align:right;"> 0.3051929 </td>
   <td style="text-align:right;"> 0.0445047 </td>
   <td style="text-align:right;"> 3.489075 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1382499 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.0545733 </td>
   <td style="text-align:right;"> 5.963175 </td>
   <td style="text-align:right;"> 0.2164386 </td>
   <td style="text-align:right;"> 0.0962896 </td>
   <td style="text-align:right;"> 2.565271 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 1.6503142 </td>
   <td style="text-align:right;"> 5.191217 </td>
   <td style="text-align:right;"> 0.4162325 </td>
   <td style="text-align:right;"> 0.0141804 </td>
   <td style="text-align:right;"> 5.040394 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1077711 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -2.3859519 </td>
   <td style="text-align:right;"> 5.616571 </td>
   <td style="text-align:right;"> 0.2986463 </td>
   <td style="text-align:right;"> 0.0472960 </td>
   <td style="text-align:right;"> 3.412947 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1382499 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 1.1731638 </td>
   <td style="text-align:right;"> 6.057722 </td>
   <td style="text-align:right;"> 0.2333804 </td>
   <td style="text-align:right;"> 0.0838022 </td>
   <td style="text-align:right;"> 2.725091 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -1.5935663 </td>
   <td style="text-align:right;"> 6.144839 </td>
   <td style="text-align:right;"> 0.1589702 </td>
   <td style="text-align:right;"> 0.1501480 </td>
   <td style="text-align:right;"> 2.071105 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> -0.9096423 </td>
   <td style="text-align:right;"> 5.283493 </td>
   <td style="text-align:right;"> 0.3696189 </td>
   <td style="text-align:right;"> 0.0235576 </td>
   <td style="text-align:right;"> 4.322604 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1182722 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -1.2324058 </td>
   <td style="text-align:right;"> 5.410324 </td>
   <td style="text-align:right;"> 0.3389244 </td>
   <td style="text-align:right;"> 0.0321792 </td>
   <td style="text-align:right;"> 3.905223 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1222810 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> -0.5285131 </td>
   <td style="text-align:right;"> 6.061566 </td>
   <td style="text-align:right;"> 0.1790117 </td>
   <td style="text-align:right;"> 0.1291877 </td>
   <td style="text-align:right;"> 2.235583 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> -1.3330810 </td>
   <td style="text-align:right;"> 8.654650 </td>
   <td style="text-align:right;"> 0.3054611 </td>
   <td style="text-align:right;"> 0.0623245 </td>
   <td style="text-align:right;"> 3.199021 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.1480206 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> -0.9701808 </td>
   <td style="text-align:right;"> 6.010835 </td>
   <td style="text-align:right;"> 0.1843196 </td>
   <td style="text-align:right;"> 0.1240443 </td>
   <td style="text-align:right;"> 2.280499 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1686395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -0.6927444 </td>
   <td style="text-align:right;"> 6.224055 </td>
   <td style="text-align:right;"> 0.1432865 </td>
   <td style="text-align:right;"> 0.1683406 </td>
   <td style="text-align:right;"> 1.947758 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1728904 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_NIH3_baseline.tPa.Prot, file="lm_3mos_bs_Prot_tpa.csv")
```
<br>  

__Model: NIHSS_3mo ~ tPA + Protein LFC__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ tPa + i, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_NIH3_tPa.Prot <- out_df
rownames(lm_NIH3_tPa.Prot) <- rn[,1]
lm_NIH3_tPa.Prot$p.adj <- p.adjust(lm_NIH3_tPa.Prot$p, method = "BH")

kable(lm_NIH3_tPa.Prot) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 9.9572937 </td>
   <td style="text-align:right;"> 2.945714 </td>
   <td style="text-align:right;"> 0.0270856 </td>
   <td style="text-align:right;"> 0.3183299 </td>
   <td style="text-align:right;"> 1.2366366 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.5498425 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 11.5252795 </td>
   <td style="text-align:right;"> 2.865387 </td>
   <td style="text-align:right;"> -0.0105529 </td>
   <td style="text-align:right;"> 0.4231655 </td>
   <td style="text-align:right;"> 0.9112374 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6068641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 11.4742080 </td>
   <td style="text-align:right;"> 2.713229 </td>
   <td style="text-align:right;"> 0.0870669 </td>
   <td style="text-align:right;"> 0.1975193 </td>
   <td style="text-align:right;"> 1.8106498 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4415137 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 8.0562118 </td>
   <td style="text-align:right;"> 4.108854 </td>
   <td style="text-align:right;"> -0.0130881 </td>
   <td style="text-align:right;"> 0.4311929 </td>
   <td style="text-align:right;"> 0.8901884 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6068641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 10.3602621 </td>
   <td style="text-align:right;"> 2.304339 </td>
   <td style="text-align:right;"> 0.3478916 </td>
   <td style="text-align:right;"> 0.0158386 </td>
   <td style="text-align:right;"> 5.5346434 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1504664 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 9.8753335 </td>
   <td style="text-align:right;"> 3.217725 </td>
   <td style="text-align:right;"> 0.2897539 </td>
   <td style="text-align:right;"> 0.0726337 </td>
   <td style="text-align:right;"> 3.4477758 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.2300067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 11.2086544 </td>
   <td style="text-align:right;"> 2.958066 </td>
   <td style="text-align:right;"> -0.0883478 </td>
   <td style="text-align:right;"> 0.7380321 </td>
   <td style="text-align:right;"> 0.3100035 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7414067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 11.1480616 </td>
   <td style="text-align:right;"> 2.989055 </td>
   <td style="text-align:right;"> -0.0884068 </td>
   <td style="text-align:right;"> 0.7383322 </td>
   <td style="text-align:right;"> 0.3095802 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7414067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 10.7817540 </td>
   <td style="text-align:right;"> 2.910757 </td>
   <td style="text-align:right;"> -0.0294898 </td>
   <td style="text-align:right;"> 0.4863879 </td>
   <td style="text-align:right;"> 0.7565173 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6160914 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 10.2780322 </td>
   <td style="text-align:right;"> 2.871270 </td>
   <td style="text-align:right;"> 0.0308307 </td>
   <td style="text-align:right;"> 0.3092537 </td>
   <td style="text-align:right;"> 1.2703974 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.5498425 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 12.2151796 </td>
   <td style="text-align:right;"> 2.560776 </td>
   <td style="text-align:right;"> 0.2067984 </td>
   <td style="text-align:right;"> 0.0688168 </td>
   <td style="text-align:right;"> 3.2160646 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2300067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 10.8924970 </td>
   <td style="text-align:right;"> 2.037128 </td>
   <td style="text-align:right;"> 0.4842473 </td>
   <td style="text-align:right;"> 0.0027266 </td>
   <td style="text-align:right;"> 8.9807684 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0345375 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 11.1752066 </td>
   <td style="text-align:right;"> 2.935441 </td>
   <td style="text-align:right;"> -0.0720561 </td>
   <td style="text-align:right;"> 0.6590961 </td>
   <td style="text-align:right;"> 0.4286898 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7414067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 11.5623180 </td>
   <td style="text-align:right;"> 2.465028 </td>
   <td style="text-align:right;"> 0.2462362 </td>
   <td style="text-align:right;"> 0.0469442 </td>
   <td style="text-align:right;"> 3.7767426 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2229850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -0.6816439 </td>
   <td style="text-align:right;"> 2.780696 </td>
   <td style="text-align:right;"> 0.6329183 </td>
   <td style="text-align:right;"> 0.0002128 </td>
   <td style="text-align:right;"> 15.6556102 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0080876 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 8.6273733 </td>
   <td style="text-align:right;"> 2.082438 </td>
   <td style="text-align:right;"> 0.5058874 </td>
   <td style="text-align:right;"> 0.0019770 </td>
   <td style="text-align:right;"> 9.7025551 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0345375 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 10.6127123 </td>
   <td style="text-align:right;"> 2.594966 </td>
   <td style="text-align:right;"> 0.1710716 </td>
   <td style="text-align:right;"> 0.0957624 </td>
   <td style="text-align:right;"> 2.7542027 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2698245 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 11.0558086 </td>
   <td style="text-align:right;"> 3.121169 </td>
   <td style="text-align:right;"> -0.0879354 </td>
   <td style="text-align:right;"> 0.7359373 </td>
   <td style="text-align:right;"> 0.3129640 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7414067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 11.5239895 </td>
   <td style="text-align:right;"> 2.920159 </td>
   <td style="text-align:right;"> -0.0420750 </td>
   <td style="text-align:right;"> 0.5327946 </td>
   <td style="text-align:right;"> 0.6568026 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6531031 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 9.2124216 </td>
   <td style="text-align:right;"> 2.625197 </td>
   <td style="text-align:right;"> 0.2236120 </td>
   <td style="text-align:right;"> 0.0586015 </td>
   <td style="text-align:right;"> 3.4481340 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2300067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 6.4928953 </td>
   <td style="text-align:right;"> 2.866008 </td>
   <td style="text-align:right;"> 0.3035438 </td>
   <td style="text-align:right;"> 0.0259431 </td>
   <td style="text-align:right;"> 4.7046441 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1971672 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 11.7734330 </td>
   <td style="text-align:right;"> 2.937070 </td>
   <td style="text-align:right;"> -0.0279535 </td>
   <td style="text-align:right;"> 0.4809705 </td>
   <td style="text-align:right;"> 0.7688569 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6160914 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 11.5408719 </td>
   <td style="text-align:right;"> 2.846449 </td>
   <td style="text-align:right;"> 0.0021651 </td>
   <td style="text-align:right;"> 0.3848201 </td>
   <td style="text-align:right;"> 1.0184431 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6068641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 12.6775418 </td>
   <td style="text-align:right;"> 2.965716 </td>
   <td style="text-align:right;"> 0.0398746 </td>
   <td style="text-align:right;"> 0.2882554 </td>
   <td style="text-align:right;"> 1.3530098 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.5476852 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 12.4254524 </td>
   <td style="text-align:right;"> 2.995944 </td>
   <td style="text-align:right;"> 0.0087478 </td>
   <td style="text-align:right;"> 0.3661834 </td>
   <td style="text-align:right;"> 1.0750127 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6049987 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 12.6673691 </td>
   <td style="text-align:right;"> 2.851310 </td>
   <td style="text-align:right;"> 0.1634196 </td>
   <td style="text-align:right;"> 0.1126186 </td>
   <td style="text-align:right;"> 2.5627392 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2853004 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 9.9695280 </td>
   <td style="text-align:right;"> 2.649895 </td>
   <td style="text-align:right;"> 0.1669307 </td>
   <td style="text-align:right;"> 0.0994090 </td>
   <td style="text-align:right;"> 2.7032333 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2698245 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 10.7240310 </td>
   <td style="text-align:right;"> 2.757719 </td>
   <td style="text-align:right;"> 0.0649333 </td>
   <td style="text-align:right;"> 0.2363943 </td>
   <td style="text-align:right;"> 1.5902603 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4932164 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 11.8384106 </td>
   <td style="text-align:right;"> 2.427682 </td>
   <td style="text-align:right;"> 0.2733048 </td>
   <td style="text-align:right;"> 0.0356829 </td>
   <td style="text-align:right;"> 4.1967878 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2229850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 11.5197608 </td>
   <td style="text-align:right;"> 2.854867 </td>
   <td style="text-align:right;"> -0.0043021 </td>
   <td style="text-align:right;"> 0.4039247 </td>
   <td style="text-align:right;"> 0.9635885 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6068641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 7.8071392 </td>
   <td style="text-align:right;"> 2.906551 </td>
   <td style="text-align:right;"> 0.2056851 </td>
   <td style="text-align:right;"> 0.0695445 </td>
   <td style="text-align:right;"> 3.2010454 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2300067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 11.1002208 </td>
   <td style="text-align:right;"> 3.227587 </td>
   <td style="text-align:right;"> -0.0890100 </td>
   <td style="text-align:right;"> 0.7414067 </td>
   <td style="text-align:right;"> 0.3052544 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7414067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 8.2747488 </td>
   <td style="text-align:right;"> 2.702081 </td>
   <td style="text-align:right;"> 0.2489358 </td>
   <td style="text-align:right;"> 0.0456978 </td>
   <td style="text-align:right;"> 3.8172747 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2229850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 10.3035325 </td>
   <td style="text-align:right;"> 2.672512 </td>
   <td style="text-align:right;"> 0.1365997 </td>
   <td style="text-align:right;"> 0.1299897 </td>
   <td style="text-align:right;"> 2.3447968 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3087254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 11.4507506 </td>
   <td style="text-align:right;"> 2.885092 </td>
   <td style="text-align:right;"> -0.0272789 </td>
   <td style="text-align:right;"> 0.4786082 </td>
   <td style="text-align:right;"> 0.7742869 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6160914 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 16.6675601 </td>
   <td style="text-align:right;"> 4.487634 </td>
   <td style="text-align:right;"> 0.0697266 </td>
   <td style="text-align:right;"> 0.2466082 </td>
   <td style="text-align:right;"> 1.5621462 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.4932164 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 11.7282198 </td>
   <td style="text-align:right;"> 3.430009 </td>
   <td style="text-align:right;"> -0.0828692 </td>
   <td style="text-align:right;"> 0.7106204 </td>
   <td style="text-align:right;"> 0.3495165 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7414067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 11.4633205 </td>
   <td style="text-align:right;"> 2.962263 </td>
   <td style="text-align:right;"> -0.0659528 </td>
   <td style="text-align:right;"> 0.6314693 </td>
   <td style="text-align:right;"> 0.4740869 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7414067 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_NIH3_tPa.Prot, file="lm_3mos_Prot_tpa.csv")
```

<br>  

### mRS models

__Model: mRS ~ NIHSS_baseline + Protein LFC__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(mRS ~ NIHSS_baseline +i, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_mRS_baseline.Prot <- out_df
rownames(lm_mRS_baseline.Prot) <- rn[,1]
lm_mRS_baseline.Prot$p.adj <- p.adjust(lm_mRS_baseline.Prot$p, 
                                       method = "BH")

kable(lm_mRS_baseline.Prot) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 2.371108 </td>
   <td style="text-align:right;"> 0.8724370 </td>
   <td style="text-align:right;"> 0.1134909 </td>
   <td style="text-align:right;"> 0.1395460 </td>
   <td style="text-align:right;"> 2.2161901 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 2.563239 </td>
   <td style="text-align:right;"> 0.9229100 </td>
   <td style="text-align:right;"> -0.0072841 </td>
   <td style="text-align:right;"> 0.4132398 </td>
   <td style="text-align:right;"> 0.9313011 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 2.593084 </td>
   <td style="text-align:right;"> 0.9270040 </td>
   <td style="text-align:right;"> -0.0036709 </td>
   <td style="text-align:right;"> 0.4008083 </td>
   <td style="text-align:right;"> 0.9652537 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 2.437534 </td>
   <td style="text-align:right;"> 0.9546943 </td>
   <td style="text-align:right;"> 0.0048882 </td>
   <td style="text-align:right;"> 0.3726672 </td>
   <td style="text-align:right;"> 1.0466661 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 2.877345 </td>
   <td style="text-align:right;"> 0.8925587 </td>
   <td style="text-align:right;"> 0.1068357 </td>
   <td style="text-align:right;"> 0.1487054 </td>
   <td style="text-align:right;"> 2.1363413 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 3.017866 </td>
   <td style="text-align:right;"> 1.1414515 </td>
   <td style="text-align:right;"> -0.0903771 </td>
   <td style="text-align:right;"> 0.6664684 </td>
   <td style="text-align:right;"> 0.4197974 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.6664684 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 2.747893 </td>
   <td style="text-align:right;"> 0.9439024 </td>
   <td style="text-align:right;"> 0.0216918 </td>
   <td style="text-align:right;"> 0.3224434 </td>
   <td style="text-align:right;"> 1.2106415 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 2.580663 </td>
   <td style="text-align:right;"> 0.9594457 </td>
   <td style="text-align:right;"> -0.0070242 </td>
   <td style="text-align:right;"> 0.4123344 </td>
   <td style="text-align:right;"> 0.9337352 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 2.538652 </td>
   <td style="text-align:right;"> 0.8852775 </td>
   <td style="text-align:right;"> 0.0683624 </td>
   <td style="text-align:right;"> 0.2128188 </td>
   <td style="text-align:right;"> 1.6970977 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 2.587239 </td>
   <td style="text-align:right;"> 0.9593540 </td>
   <td style="text-align:right;"> -0.0067958 </td>
   <td style="text-align:right;"> 0.4115400 </td>
   <td style="text-align:right;"> 0.9358758 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 3.234854 </td>
   <td style="text-align:right;"> 0.9286697 </td>
   <td style="text-align:right;"> 0.1479559 </td>
   <td style="text-align:right;"> 0.0996180 </td>
   <td style="text-align:right;"> 2.6496579 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 3.398774 </td>
   <td style="text-align:right;"> 0.9991048 </td>
   <td style="text-align:right;"> 0.1280788 </td>
   <td style="text-align:right;"> 0.1211899 </td>
   <td style="text-align:right;"> 2.3954806 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 2.493922 </td>
   <td style="text-align:right;"> 0.9133365 </td>
   <td style="text-align:right;"> 0.0197196 </td>
   <td style="text-align:right;"> 0.3280105 </td>
   <td style="text-align:right;"> 1.1911051 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 2.927182 </td>
   <td style="text-align:right;"> 0.9163853 </td>
   <td style="text-align:right;"> 0.0891771 </td>
   <td style="text-align:right;"> 0.1756300 </td>
   <td style="text-align:right;"> 1.9301285 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 2.411458 </td>
   <td style="text-align:right;"> 0.9191981 </td>
   <td style="text-align:right;"> 0.0331369 </td>
   <td style="text-align:right;"> 0.2917511 </td>
   <td style="text-align:right;"> 1.3255897 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 2.997781 </td>
   <td style="text-align:right;"> 0.8127055 </td>
   <td style="text-align:right;"> 0.2526007 </td>
   <td style="text-align:right;"> 0.0327045 </td>
   <td style="text-align:right;"> 4.2107427 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 2.628085 </td>
   <td style="text-align:right;"> 0.9112088 </td>
   <td style="text-align:right;"> 0.0224371 </td>
   <td style="text-align:right;"> 0.3203615 </td>
   <td style="text-align:right;"> 1.2180444 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 2.562557 </td>
   <td style="text-align:right;"> 0.9203967 </td>
   <td style="text-align:right;"> -0.0072746 </td>
   <td style="text-align:right;"> 0.4132065 </td>
   <td style="text-align:right;"> 0.9313907 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 2.609774 </td>
   <td style="text-align:right;"> 0.9128739 </td>
   <td style="text-align:right;"> 0.0156375 </td>
   <td style="text-align:right;"> 0.3398037 </td>
   <td style="text-align:right;"> 1.1509166 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 2.626605 </td>
   <td style="text-align:right;"> 0.9398368 </td>
   <td style="text-align:right;"> -0.0013972 </td>
   <td style="text-align:right;"> 0.3931556 </td>
   <td style="text-align:right;"> 0.9867450 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 2.833198 </td>
   <td style="text-align:right;"> 0.9226066 </td>
   <td style="text-align:right;"> 0.0596947 </td>
   <td style="text-align:right;"> 0.2302480 </td>
   <td style="text-align:right;"> 1.6031014 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 2.509002 </td>
   <td style="text-align:right;"> 0.9311220 </td>
   <td style="text-align:right;"> -0.0007586 </td>
   <td style="text-align:right;"> 0.3910296 </td>
   <td style="text-align:right;"> 0.9927985 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 3.044055 </td>
   <td style="text-align:right;"> 0.9149563 </td>
   <td style="text-align:right;"> 0.1183176 </td>
   <td style="text-align:right;"> 0.1332183 </td>
   <td style="text-align:right;"> 2.2748548 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 2.550909 </td>
   <td style="text-align:right;"> 0.9493431 </td>
   <td style="text-align:right;"> -0.0071548 </td>
   <td style="text-align:right;"> 0.4127892 </td>
   <td style="text-align:right;"> 0.9325119 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 2.912482 </td>
   <td style="text-align:right;"> 0.9072688 </td>
   <td style="text-align:right;"> 0.0960492 </td>
   <td style="text-align:right;"> 0.1646800 </td>
   <td style="text-align:right;"> 2.0094220 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 2.586282 </td>
   <td style="text-align:right;"> 1.1776220 </td>
   <td style="text-align:right;"> 0.0388747 </td>
   <td style="text-align:right;"> 0.2838048 </td>
   <td style="text-align:right;"> 1.3640240 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 2.653959 </td>
   <td style="text-align:right;"> 0.8720581 </td>
   <td style="text-align:right;"> 0.1004226 </td>
   <td style="text-align:right;"> 0.1580294 </td>
   <td style="text-align:right;"> 2.0605143 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 2.612220 </td>
   <td style="text-align:right;"> 0.9129570 </td>
   <td style="text-align:right;"> 0.0160405 </td>
   <td style="text-align:right;"> 0.3386231 </td>
   <td style="text-align:right;"> 1.1548689 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 3.189082 </td>
   <td style="text-align:right;"> 0.8054971 </td>
   <td style="text-align:right;"> 0.2930221 </td>
   <td style="text-align:right;"> 0.0203873 </td>
   <td style="text-align:right;"> 4.9374788 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 2.382790 </td>
   <td style="text-align:right;"> 0.8980815 </td>
   <td style="text-align:right;"> 0.0684448 </td>
   <td style="text-align:right;"> 0.2126588 </td>
   <td style="text-align:right;"> 1.6979996 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 2.914178 </td>
   <td style="text-align:right;"> 0.9647777 </td>
   <td style="text-align:right;"> 0.0465031 </td>
   <td style="text-align:right;"> 0.2591935 </td>
   <td style="text-align:right;"> 1.4633254 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 2.507233 </td>
   <td style="text-align:right;"> 0.9124729 </td>
   <td style="text-align:right;"> 0.0179718 </td>
   <td style="text-align:right;"> 0.3330151 </td>
   <td style="text-align:right;"> 1.1738564 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 2.742828 </td>
   <td style="text-align:right;"> 0.8252821 </td>
   <td style="text-align:right;"> 0.1988061 </td>
   <td style="text-align:right;"> 0.0590446 </td>
   <td style="text-align:right;"> 3.3573039 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 2.610161 </td>
   <td style="text-align:right;"> 0.9057822 </td>
   <td style="text-align:right;"> 0.0285499 </td>
   <td style="text-align:right;"> 0.3037277 </td>
   <td style="text-align:right;"> 1.2791950 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 2.579671 </td>
   <td style="text-align:right;"> 0.9147338 </td>
   <td style="text-align:right;"> 0.0063604 </td>
   <td style="text-align:right;"> 0.3680068 </td>
   <td style="text-align:right;"> 1.0608105 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 1.726190 </td>
   <td style="text-align:right;"> 1.1874912 </td>
   <td style="text-align:right;"> 0.0637891 </td>
   <td style="text-align:right;"> 0.2385724 </td>
   <td style="text-align:right;"> 1.5791509 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 2.562017 </td>
   <td style="text-align:right;"> 0.9219178 </td>
   <td style="text-align:right;"> -0.0072941 </td>
   <td style="text-align:right;"> 0.4132747 </td>
   <td style="text-align:right;"> 0.9312076 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 2.557071 </td>
   <td style="text-align:right;"> 0.9194479 </td>
   <td style="text-align:right;"> -0.0047234 </td>
   <td style="text-align:right;"> 0.4043949 </td>
   <td style="text-align:right;"> 0.9553385 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.4244443 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_mRS_baseline.Prot, file="lm_mRS_baseline_Prot.csv")
```
<br>  

__Model: mRS ~ NIHSS_baseline + tPA + Protein LFC__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(mRS ~ NIHSS_baseline + tPa + i, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_mRS_baseline.tPa.Prot <- out_df
rownames(lm_mRS_baseline.tPa.Prot) <- rn[,1]
lm_mRS_baseline.tPa.Prot$p.adj <- p.adjust(lm_mRS_baseline.tPa.Prot$p, 
                                           method = "BH")

kable(lm_mRS_baseline.tPa.Prot) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 2.599040 </td>
   <td style="text-align:right;"> 0.9066326 </td>
   <td style="text-align:right;"> 0.1089221 </td>
   <td style="text-align:right;"> 0.1925623 </td>
   <td style="text-align:right;"> 1.7741638 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 2.820346 </td>
   <td style="text-align:right;"> 0.9610369 </td>
   <td style="text-align:right;"> -0.0101267 </td>
   <td style="text-align:right;"> 0.4460106 </td>
   <td style="text-align:right;"> 0.9365070 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 2.937106 </td>
   <td style="text-align:right;"> 0.9695143 </td>
   <td style="text-align:right;"> 0.0116167 </td>
   <td style="text-align:right;"> 0.3878463 </td>
   <td style="text-align:right;"> 1.0744371 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 2.723904 </td>
   <td style="text-align:right;"> 1.0157311 </td>
   <td style="text-align:right;"> -0.0091507 </td>
   <td style="text-align:right;"> 0.4432822 </td>
   <td style="text-align:right;"> 0.9425709 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 3.029771 </td>
   <td style="text-align:right;"> 0.9258550 </td>
   <td style="text-align:right;"> 0.0840670 </td>
   <td style="text-align:right;"> 0.2329336 </td>
   <td style="text-align:right;"> 1.5812921 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 3.211201 </td>
   <td style="text-align:right;"> 1.1194368 </td>
   <td style="text-align:right;"> -0.0303175 </td>
   <td style="text-align:right;"> 0.4891425 </td>
   <td style="text-align:right;"> 0.8626817 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.4891425 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 2.944288 </td>
   <td style="text-align:right;"> 0.9740824 </td>
   <td style="text-align:right;"> 0.0105317 </td>
   <td style="text-align:right;"> 0.3906183 </td>
   <td style="text-align:right;"> 1.0674106 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 2.801541 </td>
   <td style="text-align:right;"> 0.9892435 </td>
   <td style="text-align:right;"> -0.0121489 </td>
   <td style="text-align:right;"> 0.4516986 </td>
   <td style="text-align:right;"> 0.9239802 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 2.764370 </td>
   <td style="text-align:right;"> 0.9201917 </td>
   <td style="text-align:right;"> 0.0619529 </td>
   <td style="text-align:right;"> 0.2740962 </td>
   <td style="text-align:right;"> 1.4182822 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 2.770438 </td>
   <td style="text-align:right;"> 0.9799717 </td>
   <td style="text-align:right;"> -0.0110167 </td>
   <td style="text-align:right;"> 0.4485079 </td>
   <td style="text-align:right;"> 0.9309882 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 3.519085 </td>
   <td style="text-align:right;"> 0.9522339 </td>
   <td style="text-align:right;"> 0.1641795 </td>
   <td style="text-align:right;"> 0.1225462 </td>
   <td style="text-align:right;"> 2.2440511 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 3.662448 </td>
   <td style="text-align:right;"> 1.0233954 </td>
   <td style="text-align:right;"> 0.1369027 </td>
   <td style="text-align:right;"> 0.1539478 </td>
   <td style="text-align:right;"> 2.0045800 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 2.748489 </td>
   <td style="text-align:right;"> 0.9356212 </td>
   <td style="text-align:right;"> 0.0330208 </td>
   <td style="text-align:right;"> 0.3359559 </td>
   <td style="text-align:right;"> 1.2162735 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 3.095302 </td>
   <td style="text-align:right;"> 0.9467999 </td>
   <td style="text-align:right;"> 0.0722055 </td>
   <td style="text-align:right;"> 0.2543747 </td>
   <td style="text-align:right;"> 1.4928911 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 2.651620 </td>
   <td style="text-align:right;"> 0.9954240 </td>
   <td style="text-align:right;"> 0.0028053 </td>
   <td style="text-align:right;"> 0.4107545 </td>
   <td style="text-align:right;"> 1.0178168 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 3.175109 </td>
   <td style="text-align:right;"> 0.8414866 </td>
   <td style="text-align:right;"> 0.2434810 </td>
   <td style="text-align:right;"> 0.0595074 </td>
   <td style="text-align:right;"> 3.0383444 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 2.826001 </td>
   <td style="text-align:right;"> 0.9477675 </td>
   <td style="text-align:right;"> 0.0062303 </td>
   <td style="text-align:right;"> 0.4017425 </td>
   <td style="text-align:right;"> 1.0397059 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 2.805541 </td>
   <td style="text-align:right;"> 0.9557694 </td>
   <td style="text-align:right;"> -0.0111243 </td>
   <td style="text-align:right;"> 0.4488106 </td>
   <td style="text-align:right;"> 0.9303213 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 2.888692 </td>
   <td style="text-align:right;"> 0.9440561 </td>
   <td style="text-align:right;"> 0.0255667 </td>
   <td style="text-align:right;"> 0.3534259 </td>
   <td style="text-align:right;"> 1.1661710 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 2.865849 </td>
   <td style="text-align:right;"> 0.9740485 </td>
   <td style="text-align:right;"> -0.0056774 </td>
   <td style="text-align:right;"> 0.4336619 </td>
   <td style="text-align:right;"> 0.9642460 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 2.921937 </td>
   <td style="text-align:right;"> 0.9536312 </td>
   <td style="text-align:right;"> 0.0213904 </td>
   <td style="text-align:right;"> 0.3634939 </td>
   <td style="text-align:right;"> 1.1384335 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 2.775871 </td>
   <td style="text-align:right;"> 0.9813330 </td>
   <td style="text-align:right;"> -0.0114285 </td>
   <td style="text-align:right;"> 0.4496669 </td>
   <td style="text-align:right;"> 0.9284371 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 3.141334 </td>
   <td style="text-align:right;"> 0.9450901 </td>
   <td style="text-align:right;"> 0.0852010 </td>
   <td style="text-align:right;"> 0.2309593 </td>
   <td style="text-align:right;"> 1.5898635 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 2.996944 </td>
   <td style="text-align:right;"> 1.0310725 </td>
   <td style="text-align:right;"> 0.0025994 </td>
   <td style="text-align:right;"> 0.4113005 </td>
   <td style="text-align:right;"> 1.0165059 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 3.240607 </td>
   <td style="text-align:right;"> 0.9327769 </td>
   <td style="text-align:right;"> 0.1225892 </td>
   <td style="text-align:right;"> 0.1728416 </td>
   <td style="text-align:right;"> 1.8848743 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 3.025399 </td>
   <td style="text-align:right;"> 1.2714326 </td>
   <td style="text-align:right;"> 0.0315846 </td>
   <td style="text-align:right;"> 0.3450550 </td>
   <td style="text-align:right;"> 1.1956883 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 2.935963 </td>
   <td style="text-align:right;"> 0.8946232 </td>
   <td style="text-align:right;"> 0.1207738 </td>
   <td style="text-align:right;"> 0.1753644 </td>
   <td style="text-align:right;"> 1.8699708 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 2.812435 </td>
   <td style="text-align:right;"> 0.9500323 </td>
   <td style="text-align:right;"> -0.0001737 </td>
   <td style="text-align:right;"> 0.4187038 </td>
   <td style="text-align:right;"> 0.9989004 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 3.316955 </td>
   <td style="text-align:right;"> 0.8355452 </td>
   <td style="text-align:right;"> 0.2729161 </td>
   <td style="text-align:right;"> 0.0444242 </td>
   <td style="text-align:right;"> 3.3772612 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 2.637106 </td>
   <td style="text-align:right;"> 0.9155258 </td>
   <td style="text-align:right;"> 0.0872376 </td>
   <td style="text-align:right;"> 0.2274462 </td>
   <td style="text-align:right;"> 1.6053103 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 2.984045 </td>
   <td style="text-align:right;"> 0.9879683 </td>
   <td style="text-align:right;"> 0.0122803 </td>
   <td style="text-align:right;"> 0.3861577 </td>
   <td style="text-align:right;"> 1.0787421 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 2.778757 </td>
   <td style="text-align:right;"> 0.9299009 </td>
   <td style="text-align:right;"> 0.0414744 </td>
   <td style="text-align:right;"> 0.3169123 </td>
   <td style="text-align:right;"> 1.2740364 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 2.964653 </td>
   <td style="text-align:right;"> 0.8527924 </td>
   <td style="text-align:right;"> 0.2004404 </td>
   <td style="text-align:right;"> 0.0890728 </td>
   <td style="text-align:right;"> 2.5876936 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 2.875869 </td>
   <td style="text-align:right;"> 0.9359740 </td>
   <td style="text-align:right;"> 0.0362768 </td>
   <td style="text-align:right;"> 0.3285249 </td>
   <td style="text-align:right;"> 1.2384014 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 2.878099 </td>
   <td style="text-align:right;"> 0.9437903 </td>
   <td style="text-align:right;"> 0.0235880 </td>
   <td style="text-align:right;"> 0.3581710 </td>
   <td style="text-align:right;"> 1.1529994 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 2.148371 </td>
   <td style="text-align:right;"> 1.3985758 </td>
   <td style="text-align:right;"> 0.0226317 </td>
   <td style="text-align:right;"> 0.3702280 </td>
   <td style="text-align:right;"> 1.1312161 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 2.806360 </td>
   <td style="text-align:right;"> 0.9582220 </td>
   <td style="text-align:right;"> -0.0117308 </td>
   <td style="text-align:right;"> 0.4505186 </td>
   <td style="text-align:right;"> 0.9265662 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 2.799434 </td>
   <td style="text-align:right;"> 0.9585086 </td>
   <td style="text-align:right;"> -0.0121473 </td>
   <td style="text-align:right;"> 0.4516940 </td>
   <td style="text-align:right;"> 0.9239901 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.4639066 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_mRS_baseline.tPa.Prot, file="lm_mRS_baseline_Prot_tpa.csv")
```
<br>  

__Model: mRS ~ tPA + Protein LFC__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(mRS ~ tPa + i, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_mRS_tPa.Prot <- out_df
rownames(lm_mRS_tPa.Prot) <- rn[,1]
lm_mRS_tPa.Prot$p.adj <- p.adjust(lm_mRS_tPa.Prot$p, method = "BH")

kable(lm_mRS_tPa.Prot) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 3.735864 </td>
   <td style="text-align:right;"> 0.4283955 </td>
   <td style="text-align:right;"> 0.0568937 </td>
   <td style="text-align:right;"> 0.2361434 </td>
   <td style="text-align:right;"> 1.5730954 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 3.995009 </td>
   <td style="text-align:right;"> 0.4136660 </td>
   <td style="text-align:right;"> -0.0585448 </td>
   <td style="text-align:right;"> 0.6301424 </td>
   <td style="text-align:right;"> 0.4745849 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 4.011792 </td>
   <td style="text-align:right;"> 0.4053159 </td>
   <td style="text-align:right;"> -0.0163105 </td>
   <td style="text-align:right;"> 0.4457944 </td>
   <td style="text-align:right;"> 0.8475374 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 3.912173 </td>
   <td style="text-align:right;"> 0.5577727 </td>
   <td style="text-align:right;"> -0.0635910 </td>
   <td style="text-align:right;"> 0.6561374 </td>
   <td style="text-align:right;"> 0.4320048 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 3.865702 </td>
   <td style="text-align:right;"> 0.3927744 </td>
   <td style="text-align:right;"> 0.0843813 </td>
   <td style="text-align:right;"> 0.1836473 </td>
   <td style="text-align:right;"> 1.8754985 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 4.347656 </td>
   <td style="text-align:right;"> 0.4946422 </td>
   <td style="text-align:right;"> -0.0538252 </td>
   <td style="text-align:right;"> 0.5431637 </td>
   <td style="text-align:right;"> 0.6424678 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 3.970894 </td>
   <td style="text-align:right;"> 0.4047451 </td>
   <td style="text-align:right;"> -0.0090899 </td>
   <td style="text-align:right;"> 0.4195792 </td>
   <td style="text-align:right;"> 0.9144242 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 3.957631 </td>
   <td style="text-align:right;"> 0.4203016 </td>
   <td style="text-align:right;"> -0.0510361 </td>
   <td style="text-align:right;"> 0.5931437 </td>
   <td style="text-align:right;"> 0.5387004 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 4.018509 </td>
   <td style="text-align:right;"> 0.4038288 </td>
   <td style="text-align:right;"> -0.0077529 </td>
   <td style="text-align:right;"> 0.4148772 </td>
   <td style="text-align:right;"> 0.9269145 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 3.986013 </td>
   <td style="text-align:right;"> 0.4171233 </td>
   <td style="text-align:right;"> -0.0621199 </td>
   <td style="text-align:right;"> 0.6484633 </td>
   <td style="text-align:right;"> 0.4443761 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 4.133040 </td>
   <td style="text-align:right;"> 0.3663607 </td>
   <td style="text-align:right;"> 0.1892264 </td>
   <td style="text-align:right;"> 0.0653216 </td>
   <td style="text-align:right;"> 3.2172041 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.5256422 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 3.926937 </td>
   <td style="text-align:right;"> 0.3644606 </td>
   <td style="text-align:right;"> 0.1837568 </td>
   <td style="text-align:right;"> 0.0691635 </td>
   <td style="text-align:right;"> 3.1386880 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.5256422 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 3.982167 </td>
   <td style="text-align:right;"> 0.4085878 </td>
   <td style="text-align:right;"> -0.0305462 </td>
   <td style="text-align:right;"> 0.5017456 </td>
   <td style="text-align:right;"> 0.7184127 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 4.016405 </td>
   <td style="text-align:right;"> 0.3887660 </td>
   <td style="text-align:right;"> 0.0647531 </td>
   <td style="text-align:right;"> 0.2199295 </td>
   <td style="text-align:right;"> 1.6577456 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 3.733676 </td>
   <td style="text-align:right;"> 0.5831437 </td>
   <td style="text-align:right;"> -0.0417587 </td>
   <td style="text-align:right;"> 0.5500865 </td>
   <td style="text-align:right;"> 0.6191941 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 3.730182 </td>
   <td style="text-align:right;"> 0.3582387 </td>
   <td style="text-align:right;"> 0.2641845 </td>
   <td style="text-align:right;"> 0.0286382 </td>
   <td style="text-align:right;"> 4.4108460 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.5256422 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 3.967805 </td>
   <td style="text-align:right;"> 0.4123323 </td>
   <td style="text-align:right;"> -0.0389123 </td>
   <td style="text-align:right;"> 0.5374410 </td>
   <td style="text-align:right;"> 0.6441789 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 3.999278 </td>
   <td style="text-align:right;"> 0.4347152 </td>
   <td style="text-align:right;"> -0.0670520 </td>
   <td style="text-align:right;"> 0.6745090 </td>
   <td style="text-align:right;"> 0.4030336 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 4.046296 </td>
   <td style="text-align:right;"> 0.4097979 </td>
   <td style="text-align:right;"> -0.0221687 </td>
   <td style="text-align:right;"> 0.4681146 </td>
   <td style="text-align:right;"> 0.7939648 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 3.867302 </td>
   <td style="text-align:right;"> 0.4362202 </td>
   <td style="text-align:right;"> -0.0243846 </td>
   <td style="text-align:right;"> 0.4768107 </td>
   <td style="text-align:right;"> 0.7738608 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 3.627435 </td>
   <td style="text-align:right;"> 0.4812485 </td>
   <td style="text-align:right;"> 0.0364575 </td>
   <td style="text-align:right;"> 0.2833432 </td>
   <td style="text-align:right;"> 1.3594507 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 4.008773 </td>
   <td style="text-align:right;"> 0.4183304 </td>
   <td style="text-align:right;"> -0.0654205 </td>
   <td style="text-align:right;"> 0.6657930 </td>
   <td style="text-align:right;"> 0.4166670 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 4.011095 </td>
   <td style="text-align:right;"> 0.3845743 </td>
   <td style="text-align:right;"> 0.0844028 </td>
   <td style="text-align:right;"> 0.1836107 </td>
   <td style="text-align:right;"> 1.8757420 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 4.117225 </td>
   <td style="text-align:right;"> 0.4282277 </td>
   <td style="text-align:right;"> -0.0220118 </td>
   <td style="text-align:right;"> 0.4675040 </td>
   <td style="text-align:right;"> 0.7953921 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 4.226214 </td>
   <td style="text-align:right;"> 0.4006122 </td>
   <td style="text-align:right;"> 0.1038762 </td>
   <td style="text-align:right;"> 0.1529461 </td>
   <td style="text-align:right;"> 2.1012137 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 4.135377 </td>
   <td style="text-align:right;"> 0.4349347 </td>
   <td style="text-align:right;"> 0.0398075 </td>
   <td style="text-align:right;"> 0.2816088 </td>
   <td style="text-align:right;"> 1.3731207 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 3.814210 </td>
   <td style="text-align:right;"> 0.3920519 </td>
   <td style="text-align:right;"> 0.1109639 </td>
   <td style="text-align:right;"> 0.1429635 </td>
   <td style="text-align:right;"> 2.1857302 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 3.944679 </td>
   <td style="text-align:right;"> 0.4205394 </td>
   <td style="text-align:right;"> -0.0442134 </td>
   <td style="text-align:right;"> 0.5612018 </td>
   <td style="text-align:right;"> 0.5977570 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 4.096485 </td>
   <td style="text-align:right;"> 0.3448605 </td>
   <td style="text-align:right;"> 0.2708407 </td>
   <td style="text-align:right;"> 0.0265095 </td>
   <td style="text-align:right;"> 4.5287031 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.5256422 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 4.029561 </td>
   <td style="text-align:right;"> 0.4049689 </td>
   <td style="text-align:right;"> -0.0099539 </td>
   <td style="text-align:right;"> 0.4226427 </td>
   <td style="text-align:right;"> 0.9063700 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 3.713030 </td>
   <td style="text-align:right;"> 0.4527160 </td>
   <td style="text-align:right;"> 0.0301444 </td>
   <td style="text-align:right;"> 0.2995162 </td>
   <td style="text-align:right;"> 1.2952730 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 4.096971 </td>
   <td style="text-align:right;"> 0.4415212 </td>
   <td style="text-align:right;"> -0.0449862 </td>
   <td style="text-align:right;"> 0.5647418 </td>
   <td style="text-align:right;"> 0.5910291 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 3.650963 </td>
   <td style="text-align:right;"> 0.3850194 </td>
   <td style="text-align:right;"> 0.2091049 </td>
   <td style="text-align:right;"> 0.0528957 </td>
   <td style="text-align:right;"> 3.5117071 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.5256422 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 3.961666 </td>
   <td style="text-align:right;"> 0.4034427 </td>
   <td style="text-align:right;"> 0.0000349 </td>
   <td style="text-align:right;"> 0.3884020 </td>
   <td style="text-align:right;"> 1.0003313 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 4.010645 </td>
   <td style="text-align:right;"> 0.4059233 </td>
   <td style="text-align:right;"> -0.0194918 </td>
   <td style="text-align:right;"> 0.4577961 </td>
   <td style="text-align:right;"> 0.8183679 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 4.141181 </td>
   <td style="text-align:right;"> 0.6760056 </td>
   <td style="text-align:right;"> -0.0799705 </td>
   <td style="text-align:right;"> 0.6964768 </td>
   <td style="text-align:right;"> 0.3705856 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6964768 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 3.947568 </td>
   <td style="text-align:right;"> 0.4898421 </td>
   <td style="text-align:right;"> -0.0645240 </td>
   <td style="text-align:right;"> 0.6610460 </td>
   <td style="text-align:right;"> 0.4241763 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 4.002177 </td>
   <td style="text-align:right;"> 0.4168396 </td>
   <td style="text-align:right;"> -0.0668525 </td>
   <td style="text-align:right;"> 0.6734378 </td>
   <td style="text-align:right;"> 0.4046985 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.6927390 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_mRS_tPa.Prot, file="lm_mRS_tpa_Prot.csv")
```
<br>  

### Age and Sex models {.tabset}

The influence of other covariates, such as age and sex, and their impact on the association with protein change and stroke recovery scores at 3 months were of interest to reviewers. Here are various combinations of age and sex with our previously examined covariates - tPA status and NIHSS score at baseline.

#### Models with just protein, age, and sex

__Statistics for Model: NIHSS_3mo ~ Protein LFC + Sex__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ i + Sex, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.sex <- out_df
rownames(lm_Prot.sex) <- rn[,1]
lm_Prot.sex$p.adj <- p.adjust(lm_Prot.sex$p, method = "BH")

kable(lm_Prot.sex) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 11.889776 </td>
   <td style="text-align:right;"> 4.079504 </td>
   <td style="text-align:right;"> 0.0434214 </td>
   <td style="text-align:right;"> 0.2803642 </td>
   <td style="text-align:right;"> 1.385835 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 13.809569 </td>
   <td style="text-align:right;"> 4.102943 </td>
   <td style="text-align:right;"> 0.0122480 </td>
   <td style="text-align:right;"> 0.3565963 </td>
   <td style="text-align:right;"> 1.105399 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 12.318600 </td>
   <td style="text-align:right;"> 4.323884 </td>
   <td style="text-align:right;"> 0.0251588 </td>
   <td style="text-align:right;"> 0.3230885 </td>
   <td style="text-align:right;"> 1.219369 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 10.325211 </td>
   <td style="text-align:right;"> 4.081838 </td>
   <td style="text-align:right;"> 0.1088929 </td>
   <td style="text-align:right;"> 0.1647373 </td>
   <td style="text-align:right;"> 2.038696 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3294746 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 11.272666 </td>
   <td style="text-align:right;"> 2.829747 </td>
   <td style="text-align:right;"> 0.3590366 </td>
   <td style="text-align:right;"> 0.0139177 </td>
   <td style="text-align:right;"> 5.761287 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1322180 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 11.127434 </td>
   <td style="text-align:right;"> 4.458111 </td>
   <td style="text-align:right;"> 0.2455735 </td>
   <td style="text-align:right;"> 0.0982151 </td>
   <td style="text-align:right;"> 2.953061 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.2332609 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 14.414293 </td>
   <td style="text-align:right;"> 3.497576 </td>
   <td style="text-align:right;"> 0.0465999 </td>
   <td style="text-align:right;"> 0.2734521 </td>
   <td style="text-align:right;"> 1.415460 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 13.700603 </td>
   <td style="text-align:right;"> 3.349643 </td>
   <td style="text-align:right;"> 0.0282114 </td>
   <td style="text-align:right;"> 0.3155776 </td>
   <td style="text-align:right;"> 1.246758 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 15.374829 </td>
   <td style="text-align:right;"> 3.035557 </td>
   <td style="text-align:right;"> 0.2550524 </td>
   <td style="text-align:right;"> 0.0429794 </td>
   <td style="text-align:right;"> 3.910198 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1633218 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 13.788597 </td>
   <td style="text-align:right;"> 3.061161 </td>
   <td style="text-align:right;"> 0.1849520 </td>
   <td style="text-align:right;"> 0.0843707 </td>
   <td style="text-align:right;"> 2.928834 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2332609 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 13.563177 </td>
   <td style="text-align:right;"> 2.950891 </td>
   <td style="text-align:right;"> 0.2414192 </td>
   <td style="text-align:right;"> 0.0492415 </td>
   <td style="text-align:right;"> 3.705134 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1701071 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 8.846064 </td>
   <td style="text-align:right;"> 3.022773 </td>
   <td style="text-align:right;"> 0.4023450 </td>
   <td style="text-align:right;"> 0.0082356 </td>
   <td style="text-align:right;"> 6.722252 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1043179 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 13.950823 </td>
   <td style="text-align:right;"> 3.699377 </td>
   <td style="text-align:right;"> 0.0155301 </td>
   <td style="text-align:right;"> 0.3478051 </td>
   <td style="text-align:right;"> 1.134088 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 12.941490 </td>
   <td style="text-align:right;"> 2.886179 </td>
   <td style="text-align:right;"> 0.2804902 </td>
   <td style="text-align:right;"> 0.0331202 </td>
   <td style="text-align:right;"> 4.313598 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1488432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 3.138915 </td>
   <td style="text-align:right;"> 2.741424 </td>
   <td style="text-align:right;"> 0.6667697 </td>
   <td style="text-align:right;"> 0.0001030 </td>
   <td style="text-align:right;"> 18.007882 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0039146 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 10.266759 </td>
   <td style="text-align:right;"> 2.433574 </td>
   <td style="text-align:right;"> 0.5394185 </td>
   <td style="text-align:right;"> 0.0011671 </td>
   <td style="text-align:right;"> 10.954933 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0221753 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 11.050231 </td>
   <td style="text-align:right;"> 3.424055 </td>
   <td style="text-align:right;"> 0.1708348 </td>
   <td style="text-align:right;"> 0.0959678 </td>
   <td style="text-align:right;"> 2.751275 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2332609 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 14.109763 </td>
   <td style="text-align:right;"> 3.792406 </td>
   <td style="text-align:right;"> 0.0176783 </td>
   <td style="text-align:right;"> 0.3421531 </td>
   <td style="text-align:right;"> 1.152970 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 13.373137 </td>
   <td style="text-align:right;"> 3.456990 </td>
   <td style="text-align:right;"> 0.0155377 </td>
   <td style="text-align:right;"> 0.3477848 </td>
   <td style="text-align:right;"> 1.134155 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 9.279308 </td>
   <td style="text-align:right;"> 3.919805 </td>
   <td style="text-align:right;"> 0.1805006 </td>
   <td style="text-align:right;"> 0.0878886 </td>
   <td style="text-align:right;"> 2.872186 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2332609 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 7.904220 </td>
   <td style="text-align:right;"> 3.583077 </td>
   <td style="text-align:right;"> 0.3109268 </td>
   <td style="text-align:right;"> 0.0239501 </td>
   <td style="text-align:right;"> 4.835410 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1488432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 13.724230 </td>
   <td style="text-align:right;"> 3.718181 </td>
   <td style="text-align:right;"> 0.0121871 </td>
   <td style="text-align:right;"> 0.3567614 </td>
   <td style="text-align:right;"> 1.104868 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 13.217823 </td>
   <td style="text-align:right;"> 3.333999 </td>
   <td style="text-align:right;"> 0.0497396 </td>
   <td style="text-align:right;"> 0.2667701 </td>
   <td style="text-align:right;"> 1.444917 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 13.254403 </td>
   <td style="text-align:right;"> 3.523495 </td>
   <td style="text-align:right;"> 0.0173674 </td>
   <td style="text-align:right;"> 0.3429662 </td>
   <td style="text-align:right;"> 1.150232 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 14.348071 </td>
   <td style="text-align:right;"> 3.319181 </td>
   <td style="text-align:right;"> 0.0843326 </td>
   <td style="text-align:right;"> 0.2019997 </td>
   <td style="text-align:right;"> 1.782846 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 15.766860 </td>
   <td style="text-align:right;"> 3.175529 </td>
   <td style="text-align:right;"> 0.3035417 </td>
   <td style="text-align:right;"> 0.0312120 </td>
   <td style="text-align:right;"> 4.486689 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1488432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 10.052467 </td>
   <td style="text-align:right;"> 4.040366 </td>
   <td style="text-align:right;"> 0.1261811 </td>
   <td style="text-align:right;"> 0.1422257 </td>
   <td style="text-align:right;"> 2.227416 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3002543 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 12.049603 </td>
   <td style="text-align:right;"> 3.525259 </td>
   <td style="text-align:right;"> 0.0849925 </td>
   <td style="text-align:right;"> 0.2009104 </td>
   <td style="text-align:right;"> 1.789542 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 12.362871 </td>
   <td style="text-align:right;"> 2.923936 </td>
   <td style="text-align:right;"> 0.2781796 </td>
   <td style="text-align:right;"> 0.0339263 </td>
   <td style="text-align:right;"> 4.275782 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1488432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 13.056175 </td>
   <td style="text-align:right;"> 3.481932 </td>
   <td style="text-align:right;"> 0.0291188 </td>
   <td style="text-align:right;"> 0.3133743 </td>
   <td style="text-align:right;"> 1.254933 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 9.802748 </td>
   <td style="text-align:right;"> 3.482846 </td>
   <td style="text-align:right;"> 0.2314996 </td>
   <td style="text-align:right;"> 0.0542810 </td>
   <td style="text-align:right;"> 3.560502 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1718900 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 13.330285 </td>
   <td style="text-align:right;"> 3.405352 </td>
   <td style="text-align:right;"> 0.0218619 </td>
   <td style="text-align:right;"> 0.3313742 </td>
   <td style="text-align:right;"> 1.189980 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 9.746144 </td>
   <td style="text-align:right;"> 3.319675 </td>
   <td style="text-align:right;"> 0.2744800 </td>
   <td style="text-align:right;"> 0.0352523 </td>
   <td style="text-align:right;"> 4.215735 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1488432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 11.549761 </td>
   <td style="text-align:right;"> 3.368307 </td>
   <td style="text-align:right;"> 0.1533617 </td>
   <td style="text-align:right;"> 0.1122152 </td>
   <td style="text-align:right;"> 2.539707 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2508340 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 13.406537 </td>
   <td style="text-align:right;"> 3.695251 </td>
   <td style="text-align:right;"> 0.0123410 </td>
   <td style="text-align:right;"> 0.3563447 </td>
   <td style="text-align:right;"> 1.106209 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 15.062508 </td>
   <td style="text-align:right;"> 3.896719 </td>
   <td style="text-align:right;"> 0.0475449 </td>
   <td style="text-align:right;"> 0.2874273 </td>
   <td style="text-align:right;"> 1.374387 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 14.035763 </td>
   <td style="text-align:right;"> 3.724677 </td>
   <td style="text-align:right;"> 0.0170245 </td>
   <td style="text-align:right;"> 0.3438648 </td>
   <td style="text-align:right;"> 1.147214 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 14.901588 </td>
   <td style="text-align:right;"> 3.882924 </td>
   <td style="text-align:right;"> 0.0395291 </td>
   <td style="text-align:right;"> 0.2890342 </td>
   <td style="text-align:right;"> 1.349825 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.3567614 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.sex, file="lm_Prot_sex.csv")
```

__Statistics for Model: NIHSS_3mo ~ Protein LFC + Age__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ i + Age_blood, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.age <- out_df
rownames(lm_Prot.age) <- rn[,1]
lm_Prot.age$p.adj <- p.adjust(lm_Prot.age$p, method = "BH")

kable(lm_Prot.age) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 23.424502 </td>
   <td style="text-align:right;"> 13.736442 </td>
   <td style="text-align:right;"> 0.0488416 </td>
   <td style="text-align:right;"> 0.2686667 </td>
   <td style="text-align:right;"> 1.4364718 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4907894 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 27.883156 </td>
   <td style="text-align:right;"> 15.130305 </td>
   <td style="text-align:right;"> 0.0121839 </td>
   <td style="text-align:right;"> 0.3567699 </td>
   <td style="text-align:right;"> 1.1048406 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.5811336 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 23.334978 </td>
   <td style="text-align:right;"> 13.736815 </td>
   <td style="text-align:right;"> 0.0476386 </td>
   <td style="text-align:right;"> 0.2712257 </td>
   <td style="text-align:right;"> 1.4251833 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4907894 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 16.281749 </td>
   <td style="text-align:right;"> 14.439916 </td>
   <td style="text-align:right;"> -0.0201499 </td>
   <td style="text-align:right;"> 0.4542526 </td>
   <td style="text-align:right;"> 0.8321090 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6639076 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 10.097358 </td>
   <td style="text-align:right;"> 11.867237 </td>
   <td style="text-align:right;"> 0.3232394 </td>
   <td style="text-align:right;"> 0.0209209 </td>
   <td style="text-align:right;"> 5.0598321 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1397931 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 18.211988 </td>
   <td style="text-align:right;"> 15.915635 </td>
   <td style="text-align:right;"> 0.2016917 </td>
   <td style="text-align:right;"> 0.1303008 </td>
   <td style="text-align:right;"> 2.5158931 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.3300953 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 20.696002 </td>
   <td style="text-align:right;"> 15.790830 </td>
   <td style="text-align:right;"> -0.0949422 </td>
   <td style="text-align:right;"> 0.7722386 </td>
   <td style="text-align:right;"> 0.2629668 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7780796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 20.008798 </td>
   <td style="text-align:right;"> 15.050521 </td>
   <td style="text-align:right;"> -0.0960429 </td>
   <td style="text-align:right;"> 0.7780796 </td>
   <td style="text-align:right;"> 0.2551712 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7780796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 18.139256 </td>
   <td style="text-align:right;"> 14.386156 </td>
   <td style="text-align:right;"> -0.0491372 </td>
   <td style="text-align:right;"> 0.5604793 </td>
   <td style="text-align:right;"> 0.6018957 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6960109 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 19.465135 </td>
   <td style="text-align:right;"> 13.471881 </td>
   <td style="text-align:right;"> 0.0588178 </td>
   <td style="text-align:right;"> 0.2482391 </td>
   <td style="text-align:right;"> 1.5311954 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4907894 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 15.267424 </td>
   <td style="text-align:right;"> 12.584015 </td>
   <td style="text-align:right;"> 0.1987829 </td>
   <td style="text-align:right;"> 0.0742069 </td>
   <td style="text-align:right;"> 3.1088594 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2550054 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 13.554468 </td>
   <td style="text-align:right;"> 10.810165 </td>
   <td style="text-align:right;"> 0.4100596 </td>
   <td style="text-align:right;"> 0.0074710 </td>
   <td style="text-align:right;"> 6.9082363 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0946326 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 19.238521 </td>
   <td style="text-align:right;"> 15.296677 </td>
   <td style="text-align:right;"> -0.0945860 </td>
   <td style="text-align:right;"> 0.7703564 </td>
   <td style="text-align:right;"> 0.2654931 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7780796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 25.308248 </td>
   <td style="text-align:right;"> 11.724467 </td>
   <td style="text-align:right;"> 0.3042321 </td>
   <td style="text-align:right;"> 0.0257514 </td>
   <td style="text-align:right;"> 4.7167183 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1397931 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 3.208628 </td>
   <td style="text-align:right;"> 9.049416 </td>
   <td style="text-align:right;"> 0.6253190 </td>
   <td style="text-align:right;"> 0.0002482 </td>
   <td style="text-align:right;"> 15.1859665 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0049930 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 26.791434 </td>
   <td style="text-align:right;"> 8.626796 </td>
   <td style="text-align:right;"> 0.6224521 </td>
   <td style="text-align:right;"> 0.0002628 </td>
   <td style="text-align:right;"> 15.0137006 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.0049930 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 27.233829 </td>
   <td style="text-align:right;"> 12.223755 </td>
   <td style="text-align:right;"> 0.2627410 </td>
   <td style="text-align:right;"> 0.0397620 </td>
   <td style="text-align:right;"> 4.0291913 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1678840 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 20.963131 </td>
   <td style="text-align:right;"> 14.644331 </td>
   <td style="text-align:right;"> -0.0826214 </td>
   <td style="text-align:right;"> 0.7094015 </td>
   <td style="text-align:right;"> 0.3513135 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7780796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 22.766568 </td>
   <td style="text-align:right;"> 14.513102 </td>
   <td style="text-align:right;"> -0.0404950 </td>
   <td style="text-align:right;"> 0.5267657 </td>
   <td style="text-align:right;"> 0.6691889 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6902447 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 16.062013 </td>
   <td style="text-align:right;"> 12.607025 </td>
   <td style="text-align:right;"> 0.1900020 </td>
   <td style="text-align:right;"> 0.0805280 </td>
   <td style="text-align:right;"> 2.9938536 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2550054 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 12.990256 </td>
   <td style="text-align:right;"> 11.704334 </td>
   <td style="text-align:right;"> 0.3167745 </td>
   <td style="text-align:right;"> 0.0224672 </td>
   <td style="text-align:right;"> 4.9409875 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1397931 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 24.417230 </td>
   <td style="text-align:right;"> 14.864832 </td>
   <td style="text-align:right;"> -0.0326801 </td>
   <td style="text-align:right;"> 0.4978070 </td>
   <td style="text-align:right;"> 0.7310100 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6755952 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 21.585793 </td>
   <td style="text-align:right;"> 13.723926 </td>
   <td style="text-align:right;"> 0.0301482 </td>
   <td style="text-align:right;"> 0.3108908 </td>
   <td style="text-align:right;"> 1.2642256 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.5369932 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 19.539596 </td>
   <td style="text-align:right;"> 14.240563 </td>
   <td style="text-align:right;"> -0.0509536 </td>
   <td style="text-align:right;"> 0.5677983 </td>
   <td style="text-align:right;"> 0.5878927 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6960109 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 20.346354 </td>
   <td style="text-align:right;"> 13.828538 </td>
   <td style="text-align:right;"> 0.0084419 </td>
   <td style="text-align:right;"> 0.3670318 </td>
   <td style="text-align:right;"> 1.0723674 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.5811336 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 14.249313 </td>
   <td style="text-align:right;"> 14.331453 </td>
   <td style="text-align:right;"> 0.0868600 </td>
   <td style="text-align:right;"> 0.2078825 </td>
   <td style="text-align:right;"> 1.7609784 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4503615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 20.010261 </td>
   <td style="text-align:right;"> 12.727042 </td>
   <td style="text-align:right;"> 0.1596268 </td>
   <td style="text-align:right;"> 0.1061351 </td>
   <td style="text-align:right;"> 2.6145534 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2940162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 22.956586 </td>
   <td style="text-align:right;"> 13.279772 </td>
   <td style="text-align:right;"> 0.0996667 </td>
   <td style="text-align:right;"> 0.1779681 </td>
   <td style="text-align:right;"> 1.9409484 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.4226743 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 24.671454 </td>
   <td style="text-align:right;"> 11.458090 </td>
   <td style="text-align:right;"> 0.3310526 </td>
   <td style="text-align:right;"> 0.0191760 </td>
   <td style="text-align:right;"> 5.2065296 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1397931 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 21.226462 </td>
   <td style="text-align:right;"> 14.048072 </td>
   <td style="text-align:right;"> -0.0164112 </td>
   <td style="text-align:right;"> 0.4419146 </td>
   <td style="text-align:right;"> 0.8627572 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6639076 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 14.047523 </td>
   <td style="text-align:right;"> 12.498414 </td>
   <td style="text-align:right;"> 0.2191329 </td>
   <td style="text-align:right;"> 0.0611851 </td>
   <td style="text-align:right;"> 3.3853349 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2325036 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 20.403971 </td>
   <td style="text-align:right;"> 15.931740 </td>
   <td style="text-align:right;"> -0.0956828 </td>
   <td style="text-align:right;"> 0.7761648 </td>
   <td style="text-align:right;"> 0.2577195 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7780796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 17.980926 </td>
   <td style="text-align:right;"> 11.823147 </td>
   <td style="text-align:right;"> 0.2772785 </td>
   <td style="text-align:right;"> 0.0342452 </td>
   <td style="text-align:right;"> 4.2611006 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.1626647 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 21.805121 </td>
   <td style="text-align:right;"> 12.775023 </td>
   <td style="text-align:right;"> 0.1573386 </td>
   <td style="text-align:right;"> 0.1083218 </td>
   <td style="text-align:right;"> 2.5870887 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.2940162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 23.725653 </td>
   <td style="text-align:right;"> 14.575890 </td>
   <td style="text-align:right;"> -0.0274497 </td>
   <td style="text-align:right;"> 0.4792053 </td>
   <td style="text-align:right;"> 0.7729115 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.6744371 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 41.383586 </td>
   <td style="text-align:right;"> 19.115526 </td>
   <td style="text-align:right;"> 0.0902442 </td>
   <td style="text-align:right;"> 0.2133291 </td>
   <td style="text-align:right;"> 1.7439702 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.4503615 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 21.972609 </td>
   <td style="text-align:right;"> 15.176577 </td>
   <td style="text-align:right;"> -0.0821600 </td>
   <td style="text-align:right;"> 0.7071373 </td>
   <td style="text-align:right;"> 0.3546607 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7780796 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 22.903689 </td>
   <td style="text-align:right;"> 15.134971 </td>
   <td style="text-align:right;"> -0.0686826 </td>
   <td style="text-align:right;"> 0.6436991 </td>
   <td style="text-align:right;"> 0.4537183 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7643927 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.age, file="lm_Prot_age.csv")
```

__Statistics for model: NIHSS_3mo ~ Protein LFC + Age + Sex__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ i + Age_blood + Sex, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.sex.age <- out_df
rownames(lm_Prot.sex.age) <- rn[,1]
lm_Prot.sex.age$p.adj <- p.adjust(lm_Prot.sex.age$p, method = "BH")

kable(lm_Prot.sex.age) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 25.586256 </td>
   <td style="text-align:right;"> 13.922390 </td>
   <td style="text-align:right;"> 0.0471265 </td>
   <td style="text-align:right;"> 0.3195578 </td>
   <td style="text-align:right;"> 1.2802577 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 26.538030 </td>
   <td style="text-align:right;"> 15.340935 </td>
   <td style="text-align:right;"> -0.0050082 </td>
   <td style="text-align:right;"> 0.4337038 </td>
   <td style="text-align:right;"> 0.9717616 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 24.772558 </td>
   <td style="text-align:right;"> 14.104653 </td>
   <td style="text-align:right;"> 0.0160629 </td>
   <td style="text-align:right;"> 0.3847075 </td>
   <td style="text-align:right;"> 1.0925092 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 20.695160 </td>
   <td style="text-align:right;"> 13.935093 </td>
   <td style="text-align:right;"> 0.0849353 </td>
   <td style="text-align:right;"> 0.2513968 </td>
   <td style="text-align:right;"> 1.5259741 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 13.456701 </td>
   <td style="text-align:right;"> 12.506458 </td>
   <td style="text-align:right;"> 0.3148324 </td>
   <td style="text-align:right;"> 0.0406422 </td>
   <td style="text-align:right;"> 3.6038157 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2980429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 25.315413 </td>
   <td style="text-align:right;"> 16.819779 </td>
   <td style="text-align:right;"> 0.2275516 </td>
   <td style="text-align:right;"> 0.1603432 </td>
   <td style="text-align:right;"> 2.1783391 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.3584141 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 22.241249 </td>
   <td style="text-align:right;"> 15.132842 </td>
   <td style="text-align:right;"> -0.0012314 </td>
   <td style="text-align:right;"> 0.4246369 </td>
   <td style="text-align:right;"> 0.9930306 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 23.467366 </td>
   <td style="text-align:right;"> 14.604789 </td>
   <td style="text-align:right;"> -0.0071743 </td>
   <td style="text-align:right;"> 0.4389597 </td>
   <td style="text-align:right;"> 0.9596352 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 23.180098 </td>
   <td style="text-align:right;"> 12.528756 </td>
   <td style="text-align:right;"> 0.2247326 </td>
   <td style="text-align:right;"> 0.0900018 </td>
   <td style="text-align:right;"> 2.6426391 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.3337357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 24.216509 </td>
   <td style="text-align:right;"> 12.967800 </td>
   <td style="text-align:right;"> 0.1675065 </td>
   <td style="text-align:right;"> 0.1409161 </td>
   <td style="text-align:right;"> 2.1401936 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.3584141 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 18.725051 </td>
   <td style="text-align:right;"> 13.085111 </td>
   <td style="text-align:right;"> 0.1966725 </td>
   <td style="text-align:right;"> 0.1127027 </td>
   <td style="text-align:right;"> 2.3873269 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.3337357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 13.803219 </td>
   <td style="text-align:right;"> 11.876513 </td>
   <td style="text-align:right;"> 0.3680969 </td>
   <td style="text-align:right;"> 0.0239342 </td>
   <td style="text-align:right;"> 4.3009534 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2980429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 27.938548 </td>
   <td style="text-align:right;"> 15.574692 </td>
   <td style="text-align:right;"> 0.0059420 </td>
   <td style="text-align:right;"> 0.4077573 </td>
   <td style="text-align:right;"> 1.0338728 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 27.498239 </td>
   <td style="text-align:right;"> 11.883768 </td>
   <td style="text-align:right;"> 0.3077518 </td>
   <td style="text-align:right;"> 0.0434512 </td>
   <td style="text-align:right;"> 3.5192215 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2980429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 6.669984 </td>
   <td style="text-align:right;"> 9.129042 </td>
   <td style="text-align:right;"> 0.6471365 </td>
   <td style="text-align:right;"> 0.0004747 </td>
   <td style="text-align:right;"> 11.3924214 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0116834 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 28.697482 </td>
   <td style="text-align:right;"> 8.645601 </td>
   <td style="text-align:right;"> 0.6334972 </td>
   <td style="text-align:right;"> 0.0006149 </td>
   <td style="text-align:right;"> 10.7947887 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0116834 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 27.595542 </td>
   <td style="text-align:right;"> 12.694235 </td>
   <td style="text-align:right;"> 0.2140422 </td>
   <td style="text-align:right;"> 0.0981706 </td>
   <td style="text-align:right;"> 2.5432200 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.3337357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 24.302494 </td>
   <td style="text-align:right;"> 14.359793 </td>
   <td style="text-align:right;"> -0.0131981 </td>
   <td style="text-align:right;"> 0.4537886 </td>
   <td style="text-align:right;"> 0.9261848 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 25.401352 </td>
   <td style="text-align:right;"> 14.393306 </td>
   <td style="text-align:right;"> -0.0016974 </td>
   <td style="text-align:right;"> 0.4257489 </td>
   <td style="text-align:right;"> 0.9903976 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 18.513356 </td>
   <td style="text-align:right;"> 13.564179 </td>
   <td style="text-align:right;"> 0.1526652 </td>
   <td style="text-align:right;"> 0.1572665 </td>
   <td style="text-align:right;"> 2.0209696 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.3584141 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 14.873494 </td>
   <td style="text-align:right;"> 12.724179 </td>
   <td style="text-align:right;"> 0.2785660 </td>
   <td style="text-align:right;"> 0.0567667 </td>
   <td style="text-align:right;"> 3.1880600 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2980429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 25.159025 </td>
   <td style="text-align:right;"> 14.727790 </td>
   <td style="text-align:right;"> -0.0117639 </td>
   <td style="text-align:right;"> 0.4502297 </td>
   <td style="text-align:right;"> 0.9341129 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 24.295987 </td>
   <td style="text-align:right;"> 14.013744 </td>
   <td style="text-align:right;"> 0.0279217 </td>
   <td style="text-align:right;"> 0.3588422 </td>
   <td style="text-align:right;"> 1.1627679 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 24.010919 </td>
   <td style="text-align:right;"> 14.397801 </td>
   <td style="text-align:right;"> -0.0099166 </td>
   <td style="text-align:right;"> 0.4456717 </td>
   <td style="text-align:right;"> 0.9443579 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 24.240603 </td>
   <td style="text-align:right;"> 13.812014 </td>
   <td style="text-align:right;"> 0.0557187 </td>
   <td style="text-align:right;"> 0.3030131 </td>
   <td style="text-align:right;"> 1.3343697 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 17.544632 </td>
   <td style="text-align:right;"> 13.081106 </td>
   <td style="text-align:right;"> 0.2511040 </td>
   <td style="text-align:right;"> 0.0825050 </td>
   <td style="text-align:right;"> 2.7882606 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3337357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 21.754451 </td>
   <td style="text-align:right;"> 13.495138 </td>
   <td style="text-align:right;"> 0.1159791 </td>
   <td style="text-align:right;"> 0.2040570 </td>
   <td style="text-align:right;"> 1.7434382 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4307871 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 25.016047 </td>
   <td style="text-align:right;"> 13.614285 </td>
   <td style="text-align:right;"> 0.0833066 </td>
   <td style="text-align:right;"> 0.2540895 </td>
   <td style="text-align:right;"> 1.5149713 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 25.487053 </td>
   <td style="text-align:right;"> 11.971696 </td>
   <td style="text-align:right;"> 0.2912478 </td>
   <td style="text-align:right;"> 0.0506218 </td>
   <td style="text-align:right;"> 3.3286048 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2980429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 24.504377 </td>
   <td style="text-align:right;"> 14.146578 </td>
   <td style="text-align:right;"> 0.0091736 </td>
   <td style="text-align:right;"> 0.4003000 </td>
   <td style="text-align:right;"> 1.0524651 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 17.069269 </td>
   <td style="text-align:right;"> 13.327321 </td>
   <td style="text-align:right;"> 0.1950170 </td>
   <td style="text-align:right;"> 0.1141727 </td>
   <td style="text-align:right;"> 2.3728198 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.3337357 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 23.758291 </td>
   <td style="text-align:right;"> 15.486979 </td>
   <td style="text-align:right;"> -0.0134537 </td>
   <td style="text-align:right;"> 0.4544246 </td>
   <td style="text-align:right;"> 0.9247743 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 20.652383 </td>
   <td style="text-align:right;"> 12.276080 </td>
   <td style="text-align:right;"> 0.2672610 </td>
   <td style="text-align:right;"> 0.0627459 </td>
   <td style="text-align:right;"> 3.0668734 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2980429 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 24.410288 </td>
   <td style="text-align:right;"> 13.059940 </td>
   <td style="text-align:right;"> 0.1555350 </td>
   <td style="text-align:right;"> 0.1539943 </td>
   <td style="text-align:right;"> 2.0436966 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.3584141 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 25.402418 </td>
   <td style="text-align:right;"> 14.492706 </td>
   <td style="text-align:right;"> -0.0055117 </td>
   <td style="text-align:right;"> 0.4349219 </td>
   <td style="text-align:right;"> 0.9689382 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 37.390104 </td>
   <td style="text-align:right;"> 19.903650 </td>
   <td style="text-align:right;"> 0.0695432 </td>
   <td style="text-align:right;"> 0.2979425 </td>
   <td style="text-align:right;"> 1.3737044 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 26.498696 </td>
   <td style="text-align:right;"> 14.898588 </td>
   <td style="text-align:right;"> 0.0001787 </td>
   <td style="text-align:right;"> 0.4212835 </td>
   <td style="text-align:right;"> 1.0010126 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 23.157246 </td>
   <td style="text-align:right;"> 14.675101 </td>
   <td style="text-align:right;"> -0.0045730 </td>
   <td style="text-align:right;"> 0.4326526 </td>
   <td style="text-align:right;"> 0.9742046 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.4544246 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.sex.age, file="lm_Prot_sex_age.csv")
```

#### Models + NIHSS baseline scores

__Statistics for Model: NIHSS_3mo ~ NIHSS_baseline + Protein LFC + Sex__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + i + Sex, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.sex <- out_df
rownames(lm_Prot.sex) <- rn[,1]
lm_Prot.sex$p.adj <- p.adjust(lm_Prot.sex$p, method = "BH")

kable(lm_Prot.sex) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> -0.2088169 </td>
   <td style="text-align:right;"> 7.492393 </td>
   <td style="text-align:right;"> 0.1798943 </td>
   <td style="text-align:right;"> 0.1283209 </td>
   <td style="text-align:right;"> 2.243012 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 1.7200633 </td>
   <td style="text-align:right;"> 7.611666 </td>
   <td style="text-align:right;"> 0.1469100 </td>
   <td style="text-align:right;"> 0.1639922 </td>
   <td style="text-align:right;"> 1.975853 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.2304292 </td>
   <td style="text-align:right;"> 7.675928 </td>
   <td style="text-align:right;"> 0.1602443 </td>
   <td style="text-align:right;"> 0.1487407 </td>
   <td style="text-align:right;"> 2.081328 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 0.8632579 </td>
   <td style="text-align:right;"> 7.157509 </td>
   <td style="text-align:right;"> 0.1890321 </td>
   <td style="text-align:right;"> 0.1196153 </td>
   <td style="text-align:right;"> 2.320869 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 3.9202329 </td>
   <td style="text-align:right;"> 6.278070 </td>
   <td style="text-align:right;"> 0.3876863 </td>
   <td style="text-align:right;"> 0.0194491 </td>
   <td style="text-align:right;"> 4.587849 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1178289 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -9.9298757 </td>
   <td style="text-align:right;"> 8.168441 </td>
   <td style="text-align:right;"> 0.5575740 </td>
   <td style="text-align:right;"> 0.0154102 </td>
   <td style="text-align:right;"> 6.041061 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.1178289 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.6926167 </td>
   <td style="text-align:right;"> 8.878236 </td>
   <td style="text-align:right;"> 0.1474121 </td>
   <td style="text-align:right;"> 0.1633967 </td>
   <td style="text-align:right;"> 1.979764 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.2955347 </td>
   <td style="text-align:right;"> 8.313964 </td>
   <td style="text-align:right;"> 0.1573236 </td>
   <td style="text-align:right;"> 0.1519821 </td>
   <td style="text-align:right;"> 2.057939 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 5.9061402 </td>
   <td style="text-align:right;"> 7.002618 </td>
   <td style="text-align:right;"> 0.3107911 </td>
   <td style="text-align:right;"> 0.0422265 </td>
   <td style="text-align:right;"> 3.555320 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1337173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 4.0684245 </td>
   <td style="text-align:right;"> 7.119851 </td>
   <td style="text-align:right;"> 0.2475263 </td>
   <td style="text-align:right;"> 0.0744078 </td>
   <td style="text-align:right;"> 2.864051 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1602214 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 4.6090730 </td>
   <td style="text-align:right;"> 6.926569 </td>
   <td style="text-align:right;"> 0.2894267 </td>
   <td style="text-align:right;"> 0.0514690 </td>
   <td style="text-align:right;"> 3.308115 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1504479 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 6.1935215 </td>
   <td style="text-align:right;"> 6.635332 </td>
   <td style="text-align:right;"> 0.3688814 </td>
   <td style="text-align:right;"> 0.0237395 </td>
   <td style="text-align:right;"> 4.312100 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1178289 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 1.2023077 </td>
   <td style="text-align:right;"> 7.815265 </td>
   <td style="text-align:right;"> 0.1465319 </td>
   <td style="text-align:right;"> 0.1644419 </td>
   <td style="text-align:right;"> 1.972910 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 3.4097114 </td>
   <td style="text-align:right;"> 6.428693 </td>
   <td style="text-align:right;"> 0.3531485 </td>
   <td style="text-align:right;"> 0.0279068 </td>
   <td style="text-align:right;"> 4.093717 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1178289 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -4.1787404 </td>
   <td style="text-align:right;"> 4.235839 </td>
   <td style="text-align:right;"> 0.7300634 </td>
   <td style="text-align:right;"> 0.0000756 </td>
   <td style="text-align:right;"> 16.325919 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0028745 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 2.8485672 </td>
   <td style="text-align:right;"> 5.108437 </td>
   <td style="text-align:right;"> 0.5850775 </td>
   <td style="text-align:right;"> 0.0014304 </td>
   <td style="text-align:right;"> 8.990501 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0271773 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -2.6834687 </td>
   <td style="text-align:right;"> 6.609730 </td>
   <td style="text-align:right;"> 0.3603549 </td>
   <td style="text-align:right;"> 0.0259286 </td>
   <td style="text-align:right;"> 4.192413 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1178289 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 1.8987412 </td>
   <td style="text-align:right;"> 7.583684 </td>
   <td style="text-align:right;"> 0.1487640 </td>
   <td style="text-align:right;"> 0.1618015 </td>
   <td style="text-align:right;"> 1.990320 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 1.0180081 </td>
   <td style="text-align:right;"> 7.364705 </td>
   <td style="text-align:right;"> 0.1547749 </td>
   <td style="text-align:right;"> 0.1548558 </td>
   <td style="text-align:right;"> 2.037662 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 1.1828052 </td>
   <td style="text-align:right;"> 6.951526 </td>
   <td style="text-align:right;"> 0.2286309 </td>
   <td style="text-align:right;"> 0.0871629 </td>
   <td style="text-align:right;"> 2.679579 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 1.8555001 </td>
   <td style="text-align:right;"> 6.521507 </td>
   <td style="text-align:right;"> 0.3210847 </td>
   <td style="text-align:right;"> 0.0382876 </td>
   <td style="text-align:right;"> 3.679981 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1337173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 1.6930876 </td>
   <td style="text-align:right;"> 7.369200 </td>
   <td style="text-align:right;"> 0.1491284 </td>
   <td style="text-align:right;"> 0.1613737 </td>
   <td style="text-align:right;"> 1.993171 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 1.8586517 </td>
   <td style="text-align:right;"> 7.295734 </td>
   <td style="text-align:right;"> 0.1601900 </td>
   <td style="text-align:right;"> 0.1488005 </td>
   <td style="text-align:right;"> 2.080891 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 1.4579540 </td>
   <td style="text-align:right;"> 7.309820 </td>
   <td style="text-align:right;"> 0.1461831 </td>
   <td style="text-align:right;"> 0.1648574 </td>
   <td style="text-align:right;"> 1.970198 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 2.3429112 </td>
   <td style="text-align:right;"> 7.030618 </td>
   <td style="text-align:right;"> 0.2195798 </td>
   <td style="text-align:right;"> 0.0938680 </td>
   <td style="text-align:right;"> 2.594379 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 10.8393031 </td>
   <td style="text-align:right;"> 11.039468 </td>
   <td style="text-align:right;"> 0.2623588 </td>
   <td style="text-align:right;"> 0.0754072 </td>
   <td style="text-align:right;"> 2.896921 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.1602214 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -1.5293030 </td>
   <td style="text-align:right;"> 7.127669 </td>
   <td style="text-align:right;"> 0.2566818 </td>
   <td style="text-align:right;"> 0.0687964 </td>
   <td style="text-align:right;"> 2.956807 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1602214 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.7992054 </td>
   <td style="text-align:right;"> 7.097751 </td>
   <td style="text-align:right;"> 0.2020762 </td>
   <td style="text-align:right;"> 0.1080090 </td>
   <td style="text-align:right;"> 2.435098 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 0.8332060 </td>
   <td style="text-align:right;"> 6.064005 </td>
   <td style="text-align:right;"> 0.4133841 </td>
   <td style="text-align:right;"> 0.0146456 </td>
   <td style="text-align:right;"> 4.993260 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1178289 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -2.9397836 </td>
   <td style="text-align:right;"> 7.599848 </td>
   <td style="text-align:right;"> 0.2451950 </td>
   <td style="text-align:right;"> 0.0758943 </td>
   <td style="text-align:right;"> 2.840791 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1602214 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 2.9243051 </td>
   <td style="text-align:right;"> 6.949753 </td>
   <td style="text-align:right;"> 0.2465475 </td>
   <td style="text-align:right;"> 0.0750290 </td>
   <td style="text-align:right;"> 2.854268 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1602214 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 0.7143441 </td>
   <td style="text-align:right;"> 7.617995 </td>
   <td style="text-align:right;"> 0.1527188 </td>
   <td style="text-align:right;"> 0.1572050 </td>
   <td style="text-align:right;"> 2.021392 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> -0.0933342 </td>
   <td style="text-align:right;"> 6.337102 </td>
   <td style="text-align:right;"> 0.3662088 </td>
   <td style="text-align:right;"> 0.0244085 </td>
   <td style="text-align:right;"> 4.274238 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1178289 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -1.1016023 </td>
   <td style="text-align:right;"> 6.703505 </td>
   <td style="text-align:right;"> 0.3127894 </td>
   <td style="text-align:right;"> 0.0414369 </td>
   <td style="text-align:right;"> 3.579228 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1337173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.8297137 </td>
   <td style="text-align:right;"> 7.578155 </td>
   <td style="text-align:right;"> 0.1515598 </td>
   <td style="text-align:right;"> 0.1585413 </td>
   <td style="text-align:right;"> 2.012256 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> -5.1529841 </td>
   <td style="text-align:right;"> 10.117311 </td>
   <td style="text-align:right;"> 0.2504934 </td>
   <td style="text-align:right;"> 0.0947576 </td>
   <td style="text-align:right;"> 2.671056 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 1.1947170 </td>
   <td style="text-align:right;"> 7.129631 </td>
   <td style="text-align:right;"> 0.1892430 </td>
   <td style="text-align:right;"> 0.1194201 </td>
   <td style="text-align:right;"> 2.322686 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 2.7098151 </td>
   <td style="text-align:right;"> 7.832467 </td>
   <td style="text-align:right;"> 0.1569556 </td>
   <td style="text-align:right;"> 0.1523945 </td>
   <td style="text-align:right;"> 2.055004 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1648574 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.sex, file="lm_3mos_Prot_sex.csv")
```

__Statistics for Model: NIHSS_3mo ~ NIHSS_baseline + Protein LFC + Age__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + i + Age_blood, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.age <- out_df
rownames(lm_Prot.age) <- rn[,1]
lm_Prot.age$p.adj <- p.adjust(lm_Prot.age$p, method = "BH")

kable(lm_Prot.age) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 3.3803291 </td>
   <td style="text-align:right;"> 17.038993 </td>
   <td style="text-align:right;"> 0.1702848 </td>
   <td style="text-align:right;"> 0.1380099 </td>
   <td style="text-align:right;"> 2.162986 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 4.1865183 </td>
   <td style="text-align:right;"> 20.328334 </td>
   <td style="text-align:right;"> 0.1129516 </td>
   <td style="text-align:right;"> 0.2083490 </td>
   <td style="text-align:right;"> 1.721561 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 3.4935179 </td>
   <td style="text-align:right;"> 17.383667 </td>
   <td style="text-align:right;"> 0.1556355 </td>
   <td style="text-align:right;"> 0.1538807 </td>
   <td style="text-align:right;"> 2.044495 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -2.5729292 </td>
   <td style="text-align:right;"> 16.834634 </td>
   <td style="text-align:right;"> 0.1216777 </td>
   <td style="text-align:right;"> 0.1961621 </td>
   <td style="text-align:right;"> 1.785027 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> -2.6365226 </td>
   <td style="text-align:right;"> 14.173471 </td>
   <td style="text-align:right;"> 0.3764397 </td>
   <td style="text-align:right;"> 0.0219293 </td>
   <td style="text-align:right;"> 4.420934 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1171942 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -7.6593291 </td>
   <td style="text-align:right;"> 14.708705 </td>
   <td style="text-align:right;"> 0.5545841 </td>
   <td style="text-align:right;"> 0.0158642 </td>
   <td style="text-align:right;"> 5.980372 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.1171942 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> -0.0874362 </td>
   <td style="text-align:right;"> 17.145150 </td>
   <td style="text-align:right;"> 0.1185617 </td>
   <td style="text-align:right;"> 0.2004495 </td>
   <td style="text-align:right;"> 1.762220 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -2.1691207 </td>
   <td style="text-align:right;"> 16.739605 </td>
   <td style="text-align:right;"> 0.1300351 </td>
   <td style="text-align:right;"> 0.1850108 </td>
   <td style="text-align:right;"> 1.847006 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -3.6362296 </td>
   <td style="text-align:right;"> 16.742945 </td>
   <td style="text-align:right;"> 0.1394963 </td>
   <td style="text-align:right;"> 0.1729844 </td>
   <td style="text-align:right;"> 1.918624 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.7582079 </td>
   <td style="text-align:right;"> 16.481638 </td>
   <td style="text-align:right;"> 0.1748731 </td>
   <td style="text-align:right;"> 0.1333141 </td>
   <td style="text-align:right;"> 2.200963 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 0.2049807 </td>
   <td style="text-align:right;"> 15.389130 </td>
   <td style="text-align:right;"> 0.2695753 </td>
   <td style="text-align:right;"> 0.0614818 </td>
   <td style="text-align:right;"> 3.091377 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1801844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 10.2209216 </td>
   <td style="text-align:right;"> 15.044179 </td>
   <td style="text-align:right;"> 0.3727977 </td>
   <td style="text-align:right;"> 0.0227867 </td>
   <td style="text-align:right;"> 4.368165 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1171942 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -5.1045332 </td>
   <td style="text-align:right;"> 18.014313 </td>
   <td style="text-align:right;"> 0.1087595 </td>
   <td style="text-align:right;"> 0.2144055 </td>
   <td style="text-align:right;"> 1.691512 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 10.1459879 </td>
   <td style="text-align:right;"> 15.321200 </td>
   <td style="text-align:right;"> 0.3536810 </td>
   <td style="text-align:right;"> 0.0277565 </td>
   <td style="text-align:right;"> 4.100933 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1171942 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -10.4758684 </td>
   <td style="text-align:right;"> 9.640469 </td>
   <td style="text-align:right;"> 0.7184663 </td>
   <td style="text-align:right;"> 0.0001010 </td>
   <td style="text-align:right;"> 15.461176 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0038387 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 17.4666053 </td>
   <td style="text-align:right;"> 11.733918 </td>
   <td style="text-align:right;"> 0.6308098 </td>
   <td style="text-align:right;"> 0.0006463 </td>
   <td style="text-align:right;"> 10.682244 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0122802 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 8.5539801 </td>
   <td style="text-align:right;"> 14.755990 </td>
   <td style="text-align:right;"> 0.3787167 </td>
   <td style="text-align:right;"> 0.0214071 </td>
   <td style="text-align:right;"> 4.454239 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1171942 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> -1.4484791 </td>
   <td style="text-align:right;"> 17.507269 </td>
   <td style="text-align:right;"> 0.0951431 </td>
   <td style="text-align:right;"> 0.2350034 </td>
   <td style="text-align:right;"> 1.595834 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.9894503 </td>
   <td style="text-align:right;"> 17.371584 </td>
   <td style="text-align:right;"> 0.1237849 </td>
   <td style="text-align:right;"> 0.1933029 </td>
   <td style="text-align:right;"> 1.800543 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 2.7678602 </td>
   <td style="text-align:right;"> 16.135830 </td>
   <td style="text-align:right;"> 0.2230381 </td>
   <td style="text-align:right;"> 0.0912588 </td>
   <td style="text-align:right;"> 2.626698 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2039903 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 3.3306857 </td>
   <td style="text-align:right;"> 15.024160 </td>
   <td style="text-align:right;"> 0.3188879 </td>
   <td style="text-align:right;"> 0.0391018 </td>
   <td style="text-align:right;"> 3.653060 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1458151 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.5566866 </td>
   <td style="text-align:right;"> 19.266904 </td>
   <td style="text-align:right;"> 0.0989384 </td>
   <td style="text-align:right;"> 0.2291180 </td>
   <td style="text-align:right;"> 1.622212 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 2.0192169 </td>
   <td style="text-align:right;"> 17.408655 </td>
   <td style="text-align:right;"> 0.1356829 </td>
   <td style="text-align:right;"> 0.1777565 </td>
   <td style="text-align:right;"> 1.889569 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -1.2372784 </td>
   <td style="text-align:right;"> 17.173972 </td>
   <td style="text-align:right;"> 0.1017284 </td>
   <td style="text-align:right;"> 0.2248631 </td>
   <td style="text-align:right;"> 1.641744 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.9738580 </td>
   <td style="text-align:right;"> 16.205920 </td>
   <td style="text-align:right;"> 0.1870098 </td>
   <td style="text-align:right;"> 0.1215003 </td>
   <td style="text-align:right;"> 2.303487 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -0.6312117 </td>
   <td style="text-align:right;"> 17.003583 </td>
   <td style="text-align:right;"> 0.1596215 </td>
   <td style="text-align:right;"> 0.1619041 </td>
   <td style="text-align:right;"> 2.013013 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 2.3719248 </td>
   <td style="text-align:right;"> 15.632456 </td>
   <td style="text-align:right;"> 0.2608059 </td>
   <td style="text-align:right;"> 0.0663837 </td>
   <td style="text-align:right;"> 2.999340 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1801844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 4.3049786 </td>
   <td style="text-align:right;"> 16.789474 </td>
   <td style="text-align:right;"> 0.1959662 </td>
   <td style="text-align:right;"> 0.1133280 </td>
   <td style="text-align:right;"> 2.381130 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 7.7424146 </td>
   <td style="text-align:right;"> 14.073098 </td>
   <td style="text-align:right;"> 0.4227080 </td>
   <td style="text-align:right;"> 0.0131685 </td>
   <td style="text-align:right;"> 5.149279 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1171942 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -3.1868406 </td>
   <td style="text-align:right;"> 15.635926 </td>
   <td style="text-align:right;"> 0.2422924 </td>
   <td style="text-align:right;"> 0.0777787 </td>
   <td style="text-align:right;"> 2.812032 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1970393 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 2.7529892 </td>
   <td style="text-align:right;"> 16.029961 </td>
   <td style="text-align:right;"> 0.2312499 </td>
   <td style="text-align:right;"> 0.0852965 </td>
   <td style="text-align:right;"> 2.704606 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2025793 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -0.0806894 </td>
   <td style="text-align:right;"> 17.222666 </td>
   <td style="text-align:right;"> 0.1153023 </td>
   <td style="text-align:right;"> 0.2050107 </td>
   <td style="text-align:right;"> 1.738534 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 2.0674459 </td>
   <td style="text-align:right;"> 14.444470 </td>
   <td style="text-align:right;"> 0.3613681 </td>
   <td style="text-align:right;"> 0.0256600 </td>
   <td style="text-align:right;"> 4.206467 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1171942 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 1.8240285 </td>
   <td style="text-align:right;"> 15.014874 </td>
   <td style="text-align:right;"> 0.3108334 </td>
   <td style="text-align:right;"> 0.0422096 </td>
   <td style="text-align:right;"> 3.555825 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.1458151 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 1.9107780 </td>
   <td style="text-align:right;"> 17.616197 </td>
   <td style="text-align:right;"> 0.1274695 </td>
   <td style="text-align:right;"> 0.1883807 </td>
   <td style="text-align:right;"> 1.827853 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 12.8438233 </td>
   <td style="text-align:right;"> 21.082848 </td>
   <td style="text-align:right;"> 0.3025112 </td>
   <td style="text-align:right;"> 0.0638041 </td>
   <td style="text-align:right;"> 3.168573 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.1801844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> -0.1685610 </td>
   <td style="text-align:right;"> 16.698350 </td>
   <td style="text-align:right;"> 0.1472864 </td>
   <td style="text-align:right;"> 0.1635457 </td>
   <td style="text-align:right;"> 1.978785 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -0.3574103 </td>
   <td style="text-align:right;"> 18.364217 </td>
   <td style="text-align:right;"> 0.0975294 </td>
   <td style="text-align:right;"> 0.2312899 </td>
   <td style="text-align:right;"> 1.612393 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2350034 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.age, file="lm_3mos_Prot_age.csv")
```

__Statistics for model: NIHSS_3mo ~ NIHSS_baseline + Protein LFC + Age + Sex__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + i + Age_blood + Sex, 
            data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.sex.age <- out_df
rownames(lm_Prot.sex.age) <- rn[,1]
lm_Prot.sex.age$p.adj <- p.adjust(lm_Prot.sex.age$p, method = "BH")

kable(lm_Prot.sex.age) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 6.5097261 </td>
   <td style="text-align:right;"> 18.33139 </td>
   <td style="text-align:right;"> 0.1277691 </td>
   <td style="text-align:right;"> 0.2275744 </td>
   <td style="text-align:right;"> 1.622563 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 4.2550396 </td>
   <td style="text-align:right;"> 20.67419 </td>
   <td style="text-align:right;"> 0.0825299 </td>
   <td style="text-align:right;"> 0.2936666 </td>
   <td style="text-align:right;"> 1.382303 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 5.3628052 </td>
   <td style="text-align:right;"> 18.50299 </td>
   <td style="text-align:right;"> 0.1021598 </td>
   <td style="text-align:right;"> 0.2636221 </td>
   <td style="text-align:right;"> 1.483582 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 4.8658994 </td>
   <td style="text-align:right;"> 18.13700 </td>
   <td style="text-align:right;"> 0.1305602 </td>
   <td style="text-align:right;"> 0.2238595 </td>
   <td style="text-align:right;"> 1.638205 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 0.7343024 </td>
   <td style="text-align:right;"> 15.84489 </td>
   <td style="text-align:right;"> 0.3430422 </td>
   <td style="text-align:right;"> 0.0482517 </td>
   <td style="text-align:right;"> 3.219213 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2413173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -10.6084203 </td>
   <td style="text-align:right;"> 19.96062 </td>
   <td style="text-align:right;"> 0.5023594 </td>
   <td style="text-align:right;"> 0.0444953 </td>
   <td style="text-align:right;"> 4.028446 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.2413173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 4.4049891 </td>
   <td style="text-align:right;"> 18.59534 </td>
   <td style="text-align:right;"> 0.0855393 </td>
   <td style="text-align:right;"> 0.2889243 </td>
   <td style="text-align:right;"> 1.397548 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 3.1185659 </td>
   <td style="text-align:right;"> 18.75026 </td>
   <td style="text-align:right;"> 0.0954265 </td>
   <td style="text-align:right;"> 0.2736913 </td>
   <td style="text-align:right;"> 1.448347 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 8.2659142 </td>
   <td style="text-align:right;"> 16.87574 </td>
   <td style="text-align:right;"> 0.2591424 </td>
   <td style="text-align:right;"> 0.0950099 </td>
   <td style="text-align:right;"> 2.486595 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 9.0991416 </td>
   <td style="text-align:right;"> 17.77418 </td>
   <td style="text-align:right;"> 0.1956267 </td>
   <td style="text-align:right;"> 0.1487575 </td>
   <td style="text-align:right;"> 2.033616 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 4.0820901 </td>
   <td style="text-align:right;"> 17.01005 </td>
   <td style="text-align:right;"> 0.2348360 </td>
   <td style="text-align:right;"> 0.1134829 </td>
   <td style="text-align:right;"> 2.304365 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 10.4496127 </td>
   <td style="text-align:right;"> 16.21822 </td>
   <td style="text-align:right;"> 0.3246918 </td>
   <td style="text-align:right;"> 0.0564710 </td>
   <td style="text-align:right;"> 3.043423 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2413173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 4.2993692 </td>
   <td style="text-align:right;"> 21.94547 </td>
   <td style="text-align:right;"> 0.0825080 </td>
   <td style="text-align:right;"> 0.2937013 </td>
   <td style="text-align:right;"> 1.382193 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 13.7067465 </td>
   <td style="text-align:right;"> 16.50551 </td>
   <td style="text-align:right;"> 0.3273023 </td>
   <td style="text-align:right;"> 0.0552395 </td>
   <td style="text-align:right;"> 3.067845 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2413173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -7.0157353 </td>
   <td style="text-align:right;"> 10.66958 </td>
   <td style="text-align:right;"> 0.7111884 </td>
   <td style="text-align:right;"> 0.0003307 </td>
   <td style="text-align:right;"> 11.465474 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0125685 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 21.1508350 </td>
   <td style="text-align:right;"> 12.49968 </td>
   <td style="text-align:right;"> 0.6258139 </td>
   <td style="text-align:right;"> 0.0016560 </td>
   <td style="text-align:right;"> 8.107986 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0314646 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 7.2346682 </td>
   <td style="text-align:right;"> 15.90079 </td>
   <td style="text-align:right;"> 0.3353944 </td>
   <td style="text-align:right;"> 0.0515548 </td>
   <td style="text-align:right;"> 3.144770 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2413173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 4.2037110 </td>
   <td style="text-align:right;"> 18.65942 </td>
   <td style="text-align:right;"> 0.0845909 </td>
   <td style="text-align:right;"> 0.2904135 </td>
   <td style="text-align:right;"> 1.392733 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 5.4217854 </td>
   <td style="text-align:right;"> 18.65113 </td>
   <td style="text-align:right;"> 0.0944200 </td>
   <td style="text-align:right;"> 0.2752176 </td>
   <td style="text-align:right;"> 1.443125 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 5.1168844 </td>
   <td style="text-align:right;"> 17.69048 </td>
   <td style="text-align:right;"> 0.1730647 </td>
   <td style="text-align:right;"> 0.1723586 </td>
   <td style="text-align:right;"> 1.889459 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 5.0489969 </td>
   <td style="text-align:right;"> 16.60093 </td>
   <td style="text-align:right;"> 0.2713464 </td>
   <td style="text-align:right;"> 0.0866417 </td>
   <td style="text-align:right;"> 2.582675 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 3.4170158 </td>
   <td style="text-align:right;"> 19.69206 </td>
   <td style="text-align:right;"> 0.0843123 </td>
   <td style="text-align:right;"> 0.2908519 </td>
   <td style="text-align:right;"> 1.391321 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 5.6687399 </td>
   <td style="text-align:right;"> 18.62282 </td>
   <td style="text-align:right;"> 0.0990619 </td>
   <td style="text-align:right;"> 0.2682242 </td>
   <td style="text-align:right;"> 1.467305 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 4.4513993 </td>
   <td style="text-align:right;"> 18.62289 </td>
   <td style="text-align:right;"> 0.0826897 </td>
   <td style="text-align:right;"> 0.2934135 </td>
   <td style="text-align:right;"> 1.383110 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 4.0944683 </td>
   <td style="text-align:right;"> 17.82068 </td>
   <td style="text-align:right;"> 0.1602970 </td>
   <td style="text-align:right;"> 0.1868430 </td>
   <td style="text-align:right;"> 1.811313 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 11.7029722 </td>
   <td style="text-align:right;"> 19.12270 </td>
   <td style="text-align:right;"> 0.2011013 </td>
   <td style="text-align:right;"> 0.1575540 </td>
   <td style="text-align:right;"> 2.006892 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 2.9675365 </td>
   <td style="text-align:right;"> 17.37421 </td>
   <td style="text-align:right;"> 0.2045032 </td>
   <td style="text-align:right;"> 0.1401495 </td>
   <td style="text-align:right;"> 2.092573 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 6.9608237 </td>
   <td style="text-align:right;"> 18.09889 </td>
   <td style="text-align:right;"> 0.1497634 </td>
   <td style="text-align:right;"> 0.1994242 </td>
   <td style="text-align:right;"> 1.748609 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 7.3654291 </td>
   <td style="text-align:right;"> 15.37229 </td>
   <td style="text-align:right;"> 0.3785947 </td>
   <td style="text-align:right;"> 0.0350116 </td>
   <td style="text-align:right;"> 3.589337 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2413173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -1.2039974 </td>
   <td style="text-align:right;"> 18.05626 </td>
   <td style="text-align:right;"> 0.1878464 </td>
   <td style="text-align:right;"> 0.1566136 </td>
   <td style="text-align:right;"> 1.983000 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 6.1103559 </td>
   <td style="text-align:right;"> 17.53323 </td>
   <td style="text-align:right;"> 0.1910613 </td>
   <td style="text-align:right;"> 0.1533317 </td>
   <td style="text-align:right;"> 2.003797 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 5.1400974 </td>
   <td style="text-align:right;"> 18.61468 </td>
   <td style="text-align:right;"> 0.0923529 </td>
   <td style="text-align:right;"> 0.2783696 </td>
   <td style="text-align:right;"> 1.432437 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 4.7812418 </td>
   <td style="text-align:right;"> 15.99626 </td>
   <td style="text-align:right;"> 0.3232629 </td>
   <td style="text-align:right;"> 0.0571541 </td>
   <td style="text-align:right;"> 3.030134 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2413173 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 4.3060714 </td>
   <td style="text-align:right;"> 16.64600 </td>
   <td style="text-align:right;"> 0.2671208 </td>
   <td style="text-align:right;"> 0.0894730 </td>
   <td style="text-align:right;"> 2.549046 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 5.3906995 </td>
   <td style="text-align:right;"> 18.72198 </td>
   <td style="text-align:right;"> 0.0913264 </td>
   <td style="text-align:right;"> 0.2799436 </td>
   <td style="text-align:right;"> 1.427147 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 12.6908331 </td>
   <td style="text-align:right;"> 22.12635 </td>
   <td style="text-align:right;"> 0.2394339 </td>
   <td style="text-align:right;"> 0.1383589 </td>
   <td style="text-align:right;"> 2.180538 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 5.7521546 </td>
   <td style="text-align:right;"> 18.17935 </td>
   <td style="text-align:right;"> 0.1319012 </td>
   <td style="text-align:right;"> 0.2220894 </td>
   <td style="text-align:right;"> 1.645756 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 3.2589506 </td>
   <td style="text-align:right;"> 18.80154 </td>
   <td style="text-align:right;"> 0.0921793 </td>
   <td style="text-align:right;"> 0.2786354 </td>
   <td style="text-align:right;"> 1.431541 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2937013 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.sex.age, file="lm_3mos_Prot_sex_age.csv")
```


#### Models + NIHSS baseline scores + tPA status

__Statistics for Model: NIHSS_3mo ~ NIHSS_baseline + tPA + Protein LFC + Sex__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + tPa +  i + Sex, 
            data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_tpa.Prot.sex <- out_df
rownames(lm_tpa.Prot.sex) <- rn[,1]
lm_tpa.Prot.sex$p.adj <- p.adjust(lm_tpa.Prot.sex$p, method = "BH")

kable(lm_tpa.Prot.sex) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 1.2391775 </td>
   <td style="text-align:right;"> 7.737380 </td>
   <td style="text-align:right;"> 0.1656433 </td>
   <td style="text-align:right;"> 0.1806768 </td>
   <td style="text-align:right;"> 1.843744 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 2.6859362 </td>
   <td style="text-align:right;"> 7.825703 </td>
   <td style="text-align:right;"> 0.1214018 </td>
   <td style="text-align:right;"> 0.2362057 </td>
   <td style="text-align:right;"> 1.587251 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 1.2252788 </td>
   <td style="text-align:right;"> 7.670245 </td>
   <td style="text-align:right;"> 0.1730759 </td>
   <td style="text-align:right;"> 0.1723462 </td>
   <td style="text-align:right;"> 1.889529 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 2.3214416 </td>
   <td style="text-align:right;"> 7.516845 </td>
   <td style="text-align:right;"> 0.1638149 </td>
   <td style="text-align:right;"> 0.1827691 </td>
   <td style="text-align:right;"> 1.832607 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 5.1094467 </td>
   <td style="text-align:right;"> 6.572878 </td>
   <td style="text-align:right;"> 0.3679158 </td>
   <td style="text-align:right;"> 0.0386416 </td>
   <td style="text-align:right;"> 3.473788 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -8.3103055 </td>
   <td style="text-align:right;"> 6.527107 </td>
   <td style="text-align:right;"> 0.7202987 </td>
   <td style="text-align:right;"> 0.0051431 </td>
   <td style="text-align:right;"> 8.725727 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.0651458 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 1.8973687 </td>
   <td style="text-align:right;"> 9.139320 </td>
   <td style="text-align:right;"> 0.1225512 </td>
   <td style="text-align:right;"> 0.2346315 </td>
   <td style="text-align:right;"> 1.593588 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 1.1109129 </td>
   <td style="text-align:right;"> 8.621379 </td>
   <td style="text-align:right;"> 0.1332813 </td>
   <td style="text-align:right;"> 0.2202779 </td>
   <td style="text-align:right;"> 1.653552 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 7.5761012 </td>
   <td style="text-align:right;"> 7.273351 </td>
   <td style="text-align:right;"> 0.3029699 </td>
   <td style="text-align:right;"> 0.0675663 </td>
   <td style="text-align:right;"> 2.847297 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 4.5102218 </td>
   <td style="text-align:right;"> 7.459538 </td>
   <td style="text-align:right;"> 0.1973194 </td>
   <td style="text-align:right;"> 0.1470870 </td>
   <td style="text-align:right;"> 2.044758 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 5.4505957 </td>
   <td style="text-align:right;"> 7.258478 </td>
   <td style="text-align:right;"> 0.2529130 </td>
   <td style="text-align:right;"> 0.0995105 </td>
   <td style="text-align:right;"> 2.438762 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 9.0247918 </td>
   <td style="text-align:right;"> 6.735093 </td>
   <td style="text-align:right;"> 0.4088173 </td>
   <td style="text-align:right;"> 0.0261824 </td>
   <td style="text-align:right;"> 3.938979 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 2.3321852 </td>
   <td style="text-align:right;"> 8.044901 </td>
   <td style="text-align:right;"> 0.1237131 </td>
   <td style="text-align:right;"> 0.2330474 </td>
   <td style="text-align:right;"> 1.600010 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 4.3699985 </td>
   <td style="text-align:right;"> 6.772640 </td>
   <td style="text-align:right;"> 0.3221654 </td>
   <td style="text-align:right;"> 0.0576831 </td>
   <td style="text-align:right;"> 3.019966 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -5.0906670 </td>
   <td style="text-align:right;"> 4.616843 </td>
   <td style="text-align:right;"> 0.7166123 </td>
   <td style="text-align:right;"> 0.0002937 </td>
   <td style="text-align:right;"> 11.747124 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0111607 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 3.4035953 </td>
   <td style="text-align:right;"> 5.432908 </td>
   <td style="text-align:right;"> 0.5590249 </td>
   <td style="text-align:right;"> 0.0045326 </td>
   <td style="text-align:right;"> 6.387733 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0651458 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -1.8466503 </td>
   <td style="text-align:right;"> 7.187156 </td>
   <td style="text-align:right;"> 0.3183539 </td>
   <td style="text-align:right;"> 0.0595499 </td>
   <td style="text-align:right;"> 2.984907 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 3.7895023 </td>
   <td style="text-align:right;"> 8.022679 </td>
   <td style="text-align:right;"> 0.1275613 </td>
   <td style="text-align:right;"> 0.2278526 </td>
   <td style="text-align:right;"> 1.621403 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 2.4986723 </td>
   <td style="text-align:right;"> 7.658823 </td>
   <td style="text-align:right;"> 0.1355396 </td>
   <td style="text-align:right;"> 0.2173355 </td>
   <td style="text-align:right;"> 1.666362 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 2.8959695 </td>
   <td style="text-align:right;"> 7.212644 </td>
   <td style="text-align:right;"> 0.2223815 </td>
   <td style="text-align:right;"> 0.1239294 </td>
   <td style="text-align:right;"> 2.215405 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 2.0790219 </td>
   <td style="text-align:right;"> 7.009730 </td>
   <td style="text-align:right;"> 0.2696892 </td>
   <td style="text-align:right;"> 0.0877438 </td>
   <td style="text-align:right;"> 2.569440 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 2.9253860 </td>
   <td style="text-align:right;"> 7.683232 </td>
   <td style="text-align:right;"> 0.1197396 </td>
   <td style="text-align:right;"> 0.2384949 </td>
   <td style="text-align:right;"> 1.578117 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 3.0096175 </td>
   <td style="text-align:right;"> 7.659656 </td>
   <td style="text-align:right;"> 0.1239613 </td>
   <td style="text-align:right;"> 0.2327099 </td>
   <td style="text-align:right;"> 1.601384 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 3.5819963 </td>
   <td style="text-align:right;"> 7.676784 </td>
   <td style="text-align:right;"> 0.1393748 </td>
   <td style="text-align:right;"> 0.2124007 </td>
   <td style="text-align:right;"> 1.688270 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 3.5531438 </td>
   <td style="text-align:right;"> 7.393821 </td>
   <td style="text-align:right;"> 0.1878702 </td>
   <td style="text-align:right;"> 0.1565890 </td>
   <td style="text-align:right;"> 1.983154 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 17.8810839 </td>
   <td style="text-align:right;"> 12.027333 </td>
   <td style="text-align:right;"> 0.3003582 </td>
   <td style="text-align:right;"> 0.0804016 </td>
   <td style="text-align:right;"> 2.717211 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0192104 </td>
   <td style="text-align:right;"> 7.333552 </td>
   <td style="text-align:right;"> 0.2517616 </td>
   <td style="text-align:right;"> 0.1003597 </td>
   <td style="text-align:right;"> 2.430008 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 2.1643708 </td>
   <td style="text-align:right;"> 7.495096 </td>
   <td style="text-align:right;"> 0.1714382 </td>
   <td style="text-align:right;"> 0.1741578 </td>
   <td style="text-align:right;"> 1.879370 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 1.3101393 </td>
   <td style="text-align:right;"> 6.522342 </td>
   <td style="text-align:right;"> 0.3717404 </td>
   <td style="text-align:right;"> 0.0373089 </td>
   <td style="text-align:right;"> 3.514720 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -1.5361288 </td>
   <td style="text-align:right;"> 7.720165 </td>
   <td style="text-align:right;"> 0.2463251 </td>
   <td style="text-align:right;"> 0.1044433 </td>
   <td style="text-align:right;"> 2.389036 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 3.3725849 </td>
   <td style="text-align:right;"> 7.354469 </td>
   <td style="text-align:right;"> 0.1937698 </td>
   <td style="text-align:right;"> 0.1506059 </td>
   <td style="text-align:right;"> 2.021447 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 2.0540097 </td>
   <td style="text-align:right;"> 7.877860 </td>
   <td style="text-align:right;"> 0.1326683 </td>
   <td style="text-align:right;"> 0.2210813 </td>
   <td style="text-align:right;"> 1.650086 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.8581251 </td>
   <td style="text-align:right;"> 6.768259 </td>
   <td style="text-align:right;"> 0.3310623 </td>
   <td style="text-align:right;"> 0.0535027 </td>
   <td style="text-align:right;"> 3.103357 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.3372548 </td>
   <td style="text-align:right;"> 7.016046 </td>
   <td style="text-align:right;"> 0.2954783 </td>
   <td style="text-align:right;"> 0.0717599 </td>
   <td style="text-align:right;"> 2.782462 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 2.0952866 </td>
   <td style="text-align:right;"> 7.780833 </td>
   <td style="text-align:right;"> 0.1367274 </td>
   <td style="text-align:right;"> 0.2157988 </td>
   <td style="text-align:right;"> 1.673126 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> -0.5870736 </td>
   <td style="text-align:right;"> 11.267151 </td>
   <td style="text-align:right;"> 0.2431660 </td>
   <td style="text-align:right;"> 0.1352563 </td>
   <td style="text-align:right;"> 2.204851 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 2.6939662 </td>
   <td style="text-align:right;"> 7.470356 </td>
   <td style="text-align:right;"> 0.1670618 </td>
   <td style="text-align:right;"> 0.1790652 </td>
   <td style="text-align:right;"> 1.852419 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 3.3604148 </td>
   <td style="text-align:right;"> 8.056240 </td>
   <td style="text-align:right;"> 0.1214789 </td>
   <td style="text-align:right;"> 0.2360999 </td>
   <td style="text-align:right;"> 1.587675 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2384949 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_tpa.Prot.sex, file="lm_3mos_tpa_Prot_sex.csv")
```

__Statistics for Model: NIHSS_3mo ~ NIHSS_baseline +tPA +  Protein LFC + Age__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline +tPa +  i + Age_blood, 
            data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.tpa.age <- out_df
rownames(lm_Prot.tpa.age) <- rn[,1]
lm_Prot.tpa.age$p.adj <- p.adjust(lm_Prot.tpa.age$p, method = "BH")

kable(lm_Prot.tpa.age) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.2895624 </td>
   <td style="text-align:right;"> 17.67012 </td>
   <td style="text-align:right;"> 0.1493413 </td>
   <td style="text-align:right;"> 0.1999404 </td>
   <td style="text-align:right;"> 1.746128 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 1.9428131 </td>
   <td style="text-align:right;"> 20.66649 </td>
   <td style="text-align:right;"> 0.0974278 </td>
   <td style="text-align:right;"> 0.2706726 </td>
   <td style="text-align:right;"> 1.458765 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.7773535 </td>
   <td style="text-align:right;"> 17.37499 </td>
   <td style="text-align:right;"> 0.1725133 </td>
   <td style="text-align:right;"> 0.1729671 </td>
   <td style="text-align:right;"> 1.886034 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -5.5426051 </td>
   <td style="text-align:right;"> 17.50373 </td>
   <td style="text-align:right;"> 0.0959289 </td>
   <td style="text-align:right;"> 0.2729314 </td>
   <td style="text-align:right;"> 1.450958 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> -5.6136043 </td>
   <td style="text-align:right;"> 14.60526 </td>
   <td style="text-align:right;"> 0.3699856 </td>
   <td style="text-align:right;"> 0.0379158 </td>
   <td style="text-align:right;"> 3.495878 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2401332 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -29.7260912 </td>
   <td style="text-align:right;"> 10.75620 </td>
   <td style="text-align:right;"> 0.8275571 </td>
   <td style="text-align:right;"> 0.0007930 </td>
   <td style="text-align:right;"> 15.397058 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.0150675 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> -3.1964598 </td>
   <td style="text-align:right;"> 17.87763 </td>
   <td style="text-align:right;"> 0.0914367 </td>
   <td style="text-align:right;"> 0.2797742 </td>
   <td style="text-align:right;"> 1.427715 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -5.2088689 </td>
   <td style="text-align:right;"> 17.39361 </td>
   <td style="text-align:right;"> 0.1063454 </td>
   <td style="text-align:right;"> 0.2574873 </td>
   <td style="text-align:right;"> 1.505753 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -7.1275925 </td>
   <td style="text-align:right;"> 17.32821 </td>
   <td style="text-align:right;"> 0.1260123 </td>
   <td style="text-align:right;"> 0.2299340 </td>
   <td style="text-align:right;"> 1.612769 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> -1.3003930 </td>
   <td style="text-align:right;"> 17.78615 </td>
   <td style="text-align:right;"> 0.1219011 </td>
   <td style="text-align:right;"> 0.2355209 </td>
   <td style="text-align:right;"> 1.590002 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> -2.1623853 </td>
   <td style="text-align:right;"> 16.17407 </td>
   <td style="text-align:right;"> 0.2367002 </td>
   <td style="text-align:right;"> 0.1119778 </td>
   <td style="text-align:right;"> 2.317930 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2836771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 6.8428766 </td>
   <td style="text-align:right;"> 14.80164 </td>
   <td style="text-align:right;"> 0.4097175 </td>
   <td style="text-align:right;"> 0.0259498 </td>
   <td style="text-align:right;"> 3.949943 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2401332 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -11.4192068 </td>
   <td style="text-align:right;"> 19.00283 </td>
   <td style="text-align:right;"> 0.1121767 </td>
   <td style="text-align:right;"> 0.2490992 </td>
   <td style="text-align:right;"> 1.536989 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 8.3421729 </td>
   <td style="text-align:right;"> 16.52187 </td>
   <td style="text-align:right;"> 0.3114866 </td>
   <td style="text-align:right;"> 0.0630315 </td>
   <td style="text-align:right;"> 2.922719 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2661331 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -9.6184852 </td>
   <td style="text-align:right;"> 10.07870 </td>
   <td style="text-align:right;"> 0.7020426 </td>
   <td style="text-align:right;"> 0.0004020 </td>
   <td style="text-align:right;"> 11.013783 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0150675 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 18.4911230 </td>
   <td style="text-align:right;"> 12.87915 </td>
   <td style="text-align:right;"> 0.6041619 </td>
   <td style="text-align:right;"> 0.0023415 </td>
   <td style="text-align:right;"> 7.486713 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0296591 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 7.9326784 </td>
   <td style="text-align:right;"> 16.12005 </td>
   <td style="text-align:right;"> 0.3316996 </td>
   <td style="text-align:right;"> 0.0532126 </td>
   <td style="text-align:right;"> 3.109416 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2661331 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> -5.4815654 </td>
   <td style="text-align:right;"> 18.50813 </td>
   <td style="text-align:right;"> 0.0686079 </td>
   <td style="text-align:right;"> 0.3162444 </td>
   <td style="text-align:right;"> 1.313062 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3162444 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> -2.0611928 </td>
   <td style="text-align:right;"> 17.93156 </td>
   <td style="text-align:right;"> 0.1048092 </td>
   <td style="text-align:right;"> 0.2597278 </td>
   <td style="text-align:right;"> 1.497591 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> -0.4696029 </td>
   <td style="text-align:right;"> 16.55290 </td>
   <td style="text-align:right;"> 0.2172459 </td>
   <td style="text-align:right;"> 0.1284386 </td>
   <td style="text-align:right;"> 2.179547 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3050418 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 3.2595920 </td>
   <td style="text-align:right;"> 16.39024 </td>
   <td style="text-align:right;"> 0.2665058 </td>
   <td style="text-align:right;"> 0.0898909 </td>
   <td style="text-align:right;"> 2.544183 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2836771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -1.4264303 </td>
   <td style="text-align:right;"> 19.58134 </td>
   <td style="text-align:right;"> 0.0821260 </td>
   <td style="text-align:right;"> 0.2943068 </td>
   <td style="text-align:right;"> 1.380265 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3083894 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> -1.3199694 </td>
   <td style="text-align:right;"> 18.89416 </td>
   <td style="text-align:right;"> 0.0898232 </td>
   <td style="text-align:right;"> 0.2822589 </td>
   <td style="text-align:right;"> 1.419422 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -4.7589415 </td>
   <td style="text-align:right;"> 17.22786 </td>
   <td style="text-align:right;"> 0.1239538 </td>
   <td style="text-align:right;"> 0.2327202 </td>
   <td style="text-align:right;"> 1.601342 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -3.6854876 </td>
   <td style="text-align:right;"> 16.94365 </td>
   <td style="text-align:right;"> 0.1568710 </td>
   <td style="text-align:right;"> 0.1908719 </td>
   <td style="text-align:right;"> 1.790747 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -4.0759859 </td>
   <td style="text-align:right;"> 17.42770 </td>
   <td style="text-align:right;"> 0.1545716 </td>
   <td style="text-align:right;"> 0.2076717 </td>
   <td style="text-align:right;"> 1.731329 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -0.7907884 </td>
   <td style="text-align:right;"> 16.10374 </td>
   <td style="text-align:right;"> 0.2518844 </td>
   <td style="text-align:right;"> 0.1002688 </td>
   <td style="text-align:right;"> 2.430940 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2836771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 1.4977657 </td>
   <td style="text-align:right;"> 17.84272 </td>
   <td style="text-align:right;"> 0.1566471 </td>
   <td style="text-align:right;"> 0.1911372 </td>
   <td style="text-align:right;"> 1.789409 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 7.2994895 </td>
   <td style="text-align:right;"> 15.30789 </td>
   <td style="text-align:right;"> 0.3787432 </td>
   <td style="text-align:right;"> 0.0349631 </td>
   <td style="text-align:right;"> 3.590972 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2401332 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -7.0565606 </td>
   <td style="text-align:right;"> 15.96497 </td>
   <td style="text-align:right;"> 0.2503836 </td>
   <td style="text-align:right;"> 0.1013831 </td>
   <td style="text-align:right;"> 2.419567 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2836771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 1.6695105 </td>
   <td style="text-align:right;"> 17.53925 </td>
   <td style="text-align:right;"> 0.1744680 </td>
   <td style="text-align:right;"> 0.1708171 </td>
   <td style="text-align:right;"> 1.898195 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -3.1548056 </td>
   <td style="text-align:right;"> 17.81222 </td>
   <td style="text-align:right;"> 0.0948891 </td>
   <td style="text-align:right;"> 0.2745056 </td>
   <td style="text-align:right;"> 1.445557 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.4506881 </td>
   <td style="text-align:right;"> 15.36942 </td>
   <td style="text-align:right;"> 0.3215964 </td>
   <td style="text-align:right;"> 0.0579588 </td>
   <td style="text-align:right;"> 3.014708 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2661331 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -0.8403701 </td>
   <td style="text-align:right;"> 15.67524 </td>
   <td style="text-align:right;"> 0.2881117 </td>
   <td style="text-align:right;"> 0.0760764 </td>
   <td style="text-align:right;"> 2.720038 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2836771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> -0.9783384 </td>
   <td style="text-align:right;"> 18.01834 </td>
   <td style="text-align:right;"> 0.1159070 </td>
   <td style="text-align:right;"> 0.2438301 </td>
   <td style="text-align:right;"> 1.557187 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 13.5131578 </td>
   <td style="text-align:right;"> 21.43560 </td>
   <td style="text-align:right;"> 0.2800980 </td>
   <td style="text-align:right;"> 0.1071591 </td>
   <td style="text-align:right;"> 2.459042 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.2836771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> -3.1969214 </td>
   <td style="text-align:right;"> 17.37419 </td>
   <td style="text-align:right;"> 0.1228470 </td>
   <td style="text-align:right;"> 0.2342275 </td>
   <td style="text-align:right;"> 1.595221 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3064525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -2.7901684 </td>
   <td style="text-align:right;"> 18.78143 </td>
   <td style="text-align:right;"> 0.0783880 </td>
   <td style="text-align:right;"> 0.3002739 </td>
   <td style="text-align:right;"> 1.361485 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3083894 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.tpa.age, file="lm_3mos_Prot_tpa_age.csv")
```

__Statistics for model: NIHSS_3mo ~ NIHSS_baseline + tPA +  Protein LFC + Age + Sex__


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline +tPa +  i + Age_blood + Sex, 
            data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Prot.tpa.sex.age <- out_df
rownames(lm_Prot.tpa.sex.age) <- rn[,1]
lm_Prot.tpa.sex.age$p.adj <- p.adjust(lm_Prot.tpa.sex.age$p, method = "BH")

kable(lm_Prot.tpa.sex.age) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 3.2541492 </td>
   <td style="text-align:right;"> 19.15239 </td>
   <td style="text-align:right;"> 0.0971250 </td>
   <td style="text-align:right;"> 0.3036095 </td>
   <td style="text-align:right;"> 1.365748 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 2.2811602 </td>
   <td style="text-align:right;"> 21.23062 </td>
   <td style="text-align:right;"> 0.0482191 </td>
   <td style="text-align:right;"> 0.3776589 </td>
   <td style="text-align:right;"> 1.172251 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 1.2804776 </td>
   <td style="text-align:right;"> 18.91532 </td>
   <td style="text-align:right;"> 0.1041664 </td>
   <td style="text-align:right;"> 0.2936599 </td>
   <td style="text-align:right;"> 1.395347 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> 1.8040734 </td>
   <td style="text-align:right;"> 19.03453 </td>
   <td style="text-align:right;"> 0.0941999 </td>
   <td style="text-align:right;"> 0.3077972 </td>
   <td style="text-align:right;"> 1.353588 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> -2.5245036 </td>
   <td style="text-align:right;"> 16.44575 </td>
   <td style="text-align:right;"> 0.3297319 </td>
   <td style="text-align:right;"> 0.0756556 </td>
   <td style="text-align:right;"> 2.672597 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -47.0356757 </td>
   <td style="text-align:right;"> 12.08665 </td>
   <td style="text-align:right;"> 0.8805105 </td>
   <td style="text-align:right;"> 0.0006382 </td>
   <td style="text-align:right;"> 18.685446 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 0.0211151 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 1.2010483 </td>
   <td style="text-align:right;"> 19.48334 </td>
   <td style="text-align:right;"> 0.0495633 </td>
   <td style="text-align:right;"> 0.3755139 </td>
   <td style="text-align:right;"> 1.177303 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.1607546 </td>
   <td style="text-align:right;"> 19.62329 </td>
   <td style="text-align:right;"> 0.0614700 </td>
   <td style="text-align:right;"> 0.3567736 </td>
   <td style="text-align:right;"> 1.222686 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 4.7944239 </td>
   <td style="text-align:right;"> 17.45899 </td>
   <td style="text-align:right;"> 0.2468450 </td>
   <td style="text-align:right;"> 0.1336637 </td>
   <td style="text-align:right;"> 2.114343 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 7.6080213 </td>
   <td style="text-align:right;"> 19.51656 </td>
   <td style="text-align:right;"> 0.1325919 </td>
   <td style="text-align:right;"> 0.2554195 </td>
   <td style="text-align:right;"> 1.519723 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 1.6413023 </td>
   <td style="text-align:right;"> 17.94032 </td>
   <td style="text-align:right;"> 0.1943316 </td>
   <td style="text-align:right;"> 0.1832853 </td>
   <td style="text-align:right;"> 1.820098 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 6.1692020 </td>
   <td style="text-align:right;"> 16.09679 </td>
   <td style="text-align:right;"> 0.3616166 </td>
   <td style="text-align:right;"> 0.0592129 </td>
   <td style="text-align:right;"> 2.925953 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -4.1562644 </td>
   <td style="text-align:right;"> 24.61008 </td>
   <td style="text-align:right;"> 0.0568624 </td>
   <td style="text-align:right;"> 0.3639695 </td>
   <td style="text-align:right;"> 1.204988 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 11.9427542 </td>
   <td style="text-align:right;"> 17.83345 </td>
   <td style="text-align:right;"> 0.2784872 </td>
   <td style="text-align:right;"> 0.1087330 </td>
   <td style="text-align:right;"> 2.312321 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -6.0740429 </td>
   <td style="text-align:right;"> 11.16329 </td>
   <td style="text-align:right;"> 0.6932401 </td>
   <td style="text-align:right;"> 0.0011113 </td>
   <td style="text-align:right;"> 8.683588 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.0211151 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 22.3535138 </td>
   <td style="text-align:right;"> 13.71833 </td>
   <td style="text-align:right;"> 0.5970666 </td>
   <td style="text-align:right;"> 0.0051089 </td>
   <td style="text-align:right;"> 6.038119 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.0647132 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 6.8579829 </td>
   <td style="text-align:right;"> 17.19214 </td>
   <td style="text-align:right;"> 0.2804006 </td>
   <td style="text-align:right;"> 0.1073398 </td>
   <td style="text-align:right;"> 2.324851 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.1748363 </td>
   <td style="text-align:right;"> 19.59127 </td>
   <td style="text-align:right;"> 0.0581207 </td>
   <td style="text-align:right;"> 0.3619972 </td>
   <td style="text-align:right;"> 1.209804 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 2.1681669 </td>
   <td style="text-align:right;"> 19.44904 </td>
   <td style="text-align:right;"> 0.0635283 </td>
   <td style="text-align:right;"> 0.3535823 </td>
   <td style="text-align:right;"> 1.230649 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 1.4316975 </td>
   <td style="text-align:right;"> 18.33755 </td>
   <td style="text-align:right;"> 0.1581173 </td>
   <td style="text-align:right;"> 0.2237687 </td>
   <td style="text-align:right;"> 1.638567 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 4.8599648 </td>
   <td style="text-align:right;"> 17.90391 </td>
   <td style="text-align:right;"> 0.2107319 </td>
   <td style="text-align:right;"> 0.1666494 </td>
   <td style="text-align:right;"> 1.907789 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 1.3171987 </td>
   <td style="text-align:right;"> 20.31195 </td>
   <td style="text-align:right;"> 0.0469737 </td>
   <td style="text-align:right;"> 0.3796516 </td>
   <td style="text-align:right;"> 1.167583 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 2.2510024 </td>
   <td style="text-align:right;"> 19.98364 </td>
   <td style="text-align:right;"> 0.0510936 </td>
   <td style="text-align:right;"> 0.3730790 </td>
   <td style="text-align:right;"> 1.183072 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.5433447 </td>
   <td style="text-align:right;"> 19.50273 </td>
   <td style="text-align:right;"> 0.0718112 </td>
   <td style="text-align:right;"> 0.3408879 </td>
   <td style="text-align:right;"> 1.263048 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 1.2792526 </td>
   <td style="text-align:right;"> 18.73149 </td>
   <td style="text-align:right;"> 0.1214903 </td>
   <td style="text-align:right;"> 0.2699829 </td>
   <td style="text-align:right;"> 1.470191 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 10.0604557 </td>
   <td style="text-align:right;"> 18.45775 </td>
   <td style="text-align:right;"> 0.2587444 </td>
   <td style="text-align:right;"> 0.1393289 </td>
   <td style="text-align:right;"> 2.116999 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -0.7244520 </td>
   <td style="text-align:right;"> 18.04080 </td>
   <td style="text-align:right;"> 0.1895481 </td>
   <td style="text-align:right;"> 0.1883358 </td>
   <td style="text-align:right;"> 1.795190 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 4.1542830 </td>
   <td style="text-align:right;"> 19.22326 </td>
   <td style="text-align:right;"> 0.1033497 </td>
   <td style="text-align:right;"> 0.2948042 </td>
   <td style="text-align:right;"> 1.391891 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 6.9818456 </td>
   <td style="text-align:right;"> 16.59335 </td>
   <td style="text-align:right;"> 0.3272340 </td>
   <td style="text-align:right;"> 0.0770717 </td>
   <td style="text-align:right;"> 2.653763 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -6.1794083 </td>
   <td style="text-align:right;"> 18.71291 </td>
   <td style="text-align:right;"> 0.1886159 </td>
   <td style="text-align:right;"> 0.1893306 </td>
   <td style="text-align:right;"> 1.790370 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 4.9205155 </td>
   <td style="text-align:right;"> 19.00086 </td>
   <td style="text-align:right;"> 0.1271601 </td>
   <td style="text-align:right;"> 0.2624851 </td>
   <td style="text-align:right;"> 1.495331 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 1.8768738 </td>
   <td style="text-align:right;"> 19.43689 </td>
   <td style="text-align:right;"> 0.0603985 </td>
   <td style="text-align:right;"> 0.3584405 </td>
   <td style="text-align:right;"> 1.218555 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 3.1522511 </td>
   <td style="text-align:right;"> 17.02596 </td>
   <td style="text-align:right;"> 0.2766375 </td>
   <td style="text-align:right;"> 0.1100921 </td>
   <td style="text-align:right;"> 2.300271 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 1.4863948 </td>
   <td style="text-align:right;"> 17.45610 </td>
   <td style="text-align:right;"> 0.2371021 </td>
   <td style="text-align:right;"> 0.1420769 </td>
   <td style="text-align:right;"> 2.056691 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 2.1545184 </td>
   <td style="text-align:right;"> 19.42533 </td>
   <td style="text-align:right;"> 0.0647889 </td>
   <td style="text-align:right;"> 0.3516350 </td>
   <td style="text-align:right;"> 1.235543 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 13.3802704 </td>
   <td style="text-align:right;"> 22.59206 </td>
   <td style="text-align:right;"> 0.2083797 </td>
   <td style="text-align:right;"> 0.2026591 </td>
   <td style="text-align:right;"> 1.789696 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 2.5953189 </td>
   <td style="text-align:right;"> 19.05779 </td>
   <td style="text-align:right;"> 0.0976527 </td>
   <td style="text-align:right;"> 0.3028574 </td>
   <td style="text-align:right;"> 1.367951 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.7802630 </td>
   <td style="text-align:right;"> 19.60864 </td>
   <td style="text-align:right;"> 0.0499457 </td>
   <td style="text-align:right;"> 0.3749048 </td>
   <td style="text-align:right;"> 1.178743 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3796516 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Prot.tpa.sex.age, file="lm_3mos_Prot_tpa_sex_age.csv")
```

<br>  
<br>  

## Supplemental Figure 4. Plots of HDL protein logFC correlation with stroke recovery

Visualizations of the relationship between proteins with significant relation to stroke recovery at 3 months,


```r
ggplotRegression <- function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], 
                             y = names(fit$model)[1])) + 
     geom_point(size=2, alpha=0.8) +
     theme_bw() +
     stat_smooth(method = "lm", col = "Dodgerblue2", alpha=0.3)

}

r1 <- ggplotRegression(lm(APOF ~ NIHSS_3mo, data=lfc_df))
r2 <- ggplotRegression(lm(APOL1 ~ NIHSS_3mo, data=lfc_df))
r3 <- ggplotRegression(lm(APMAP ~ NIHSS_3mo, data=lfc_df))
r4 <- ggplotRegression(lm(LPA ~ NIHSS_3mo, data=lfc_df))
r5 <- ggplotRegression(lm(APOC4 ~ NIHSS_3mo, data=lfc_df))
r6 <- ggplotRegression(lm(APOM ~ NIHSS_3mo, data=lfc_df))
r7 <- ggplotRegression(lm(PCYOX1 ~ NIHSS_3mo, data=lfc_df))
r8 <- ggplotRegression(lm(PON1 ~ NIHSS_3mo, data=lfc_df))
r9 <- ggplotRegression(lm(APOE ~ NIHSS_3mo, data=lfc_df))
r10 <- ggplotRegression(lm(PPBP ~ NIHSS_3mo, data=lfc_df))

grid.arrange(r1, r2, r3, r4, r5, r6, r7, r8, r9,
             top = "Significantly associated protein changes with recovery", 
             ncol = 5)
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-42-1.png)<!-- -->

<br>  
<br>

<br>


## Supplemental Figure 5. 

Visualization of the linear regression of plasma lipid measures and stroke severity and recovery scores.


```r
ch1 <- ggplotRegression(lm(Cholesterol ~ NIHSS_baseline, data=lfc_df))
ch2 <- ggplotRegression(lm(Cholesterol ~ NIHSS_3mo, data=lfc_df))
ch3 <- ggplotRegression(lm(Cholesterol ~ mRS, data=lfc_df))
h1 <- ggplotRegression(lm(HDLC ~ NIHSS_baseline, data=lfc_df))
h2 <- ggplotRegression(lm(HDLC ~ NIHSS_3mo, data=lfc_df))
h3 <- ggplotRegression(lm(HDLC ~ mRS, data=lfc_df))
l1 <- ggplotRegression(lm(LDLC ~ NIHSS_baseline, data=lfc_df))
l2 <- ggplotRegression(lm(LDLC ~ NIHSS_3mo, data=lfc_df))
l3 <- ggplotRegression(lm(LDLC ~ mRS, data=lfc_df))
t1 <- ggplotRegression(lm(Triglycerides ~ NIHSS_baseline, data=lfc_df))
t2 <- ggplotRegression(lm(Triglycerides ~ NIHSS_3mo, data=lfc_df))
t3 <- ggplotRegression(lm(Triglycerides ~ mRS, data=lfc_df))
ce1 <- ggplotRegression(lm(CEC ~ NIHSS_baseline, data=lfc_df))
ce2 <- ggplotRegression(lm(CEC ~ NIHSS_3mo, data=lfc_df))
ce3 <- ggplotRegression(lm(CEC ~ mRS, data=lfc_df))

grid.arrange(ch1,ch2,ch3,h1,h2,h3,l1,l2,l3,t1,t2,t3,ce1,ce2,ce3, ncol = 3)
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-43-1.png)<!-- -->

<br>  
<br>  

## Supplemental Table 6.

Table with statistics of linear regressions of plasma lipid measures and stroke severity and recovery scores.


```r
cec_b <- lm_out(lm(CEC ~ NIHSS_baseline, data = lfc_df))
cec_3 <- lm_out(lm(CEC ~ NIHSS_3mo, data = lfc_df))
cec_m <- lm_out(lm(CEC ~ mRS, data = lfc_df))

hdlc_b <- lm_out(lm(HDLC ~ NIHSS_baseline, data = lfc_df))
hdlc_3 <- lm_out(lm(HDLC ~ NIHSS_3mo, data = lfc_df))
hdlc_m <- lm_out(lm(HDLC ~ mRS, data = lfc_df))

ldlc_b <- lm_out(lm(LDLC ~ NIHSS_baseline, data = lfc_df))
ldlc_3 <- lm_out(lm(LDLC ~ NIHSS_3mo, data = lfc_df))
ldlc_m <- lm_out(lm(LDLC ~ mRS, data = lfc_df))

chol_b <- lm_out(lm(Cholesterol ~ NIHSS_baseline, data = lfc_df))
chol_3 <- lm_out(lm(Cholesterol  ~ NIHSS_3mo, data = lfc_df))
chol_m <- lm_out(lm(Cholesterol  ~ mRS, data = lfc_df))

trig_b <- lm_out(lm(Triglycerides ~ NIHSS_baseline, data = lfc_df))
trig_3 <- lm_out(lm(Triglycerides ~ NIHSS_3mo, data = lfc_df))
trig_m <- lm_out(lm(Triglycerides ~ mRS, data = lfc_df))


lm_lipids <- as.data.frame(rbind(cec_b, cec_3, cec_m, 
                                 hdlc_b, hdlc_3, hdlc_m, 
                                 ldlc_b, ldlc_3, ldlc_m, 
                                 chol_b, chol_3, chol_m, 
                                 trig_b, trig_3, trig_m))

row.names(lm_lipids) <- c("cec_b", "cec_3", "cec_m", 
                          "hdlc_b", "hdlc_3", "hdlc_m", 
                          "ldlc_b", "ldlc_3", "ldlc_m", 
                          "chol_b", "chol_3", "chol_m", 
                          "trig_b", "trig_3", "trig_m")
#write.csv(lm_lipids, file="lm_lipids.csv")

kable((lm_lipids)) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> cec_b </td>
   <td style="text-align:right;"> 13.03154 </td>
   <td style="text-align:right;"> 1.6418388 </td>
   <td style="text-align:right;"> 0.1312226 </td>
   <td style="text-align:right;"> 0.0647695 </td>
   <td style="text-align:right;"> 3.8698151 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cec_3 </td>
   <td style="text-align:right;"> 10.27024 </td>
   <td style="text-align:right;"> 0.9166501 </td>
   <td style="text-align:right;"> -0.0359908 </td>
   <td style="text-align:right;"> 0.5313244 </td>
   <td style="text-align:right;"> 0.4094116 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cec_m </td>
   <td style="text-align:right;"> 10.04931 </td>
   <td style="text-align:right;"> 1.7968614 </td>
   <td style="text-align:right;"> -0.0554495 </td>
   <td style="text-align:right;"> 0.9665424 </td>
   <td style="text-align:right;"> 0.0018090 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hdlc_b </td>
   <td style="text-align:right;"> 37.58491 </td>
   <td style="text-align:right;"> 7.5632716 </td>
   <td style="text-align:right;"> -0.0452372 </td>
   <td style="text-align:right;"> 0.6783559 </td>
   <td style="text-align:right;"> 0.1776914 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hdlc_3 </td>
   <td style="text-align:right;"> 38.37018 </td>
   <td style="text-align:right;"> 3.2842553 </td>
   <td style="text-align:right;"> -0.0580929 </td>
   <td style="text-align:right;"> 0.7995816 </td>
   <td style="text-align:right;"> 0.0666419 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> hdlc_m </td>
   <td style="text-align:right;"> 35.37778 </td>
   <td style="text-align:right;"> 7.4335169 </td>
   <td style="text-align:right;"> -0.0241061 </td>
   <td style="text-align:right;"> 0.4667848 </td>
   <td style="text-align:right;"> 0.5527644 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldlc_b </td>
   <td style="text-align:right;"> 175.71924 </td>
   <td style="text-align:right;"> 38.2922736 </td>
   <td style="text-align:right;"> 0.0187504 </td>
   <td style="text-align:right;"> 0.2582403 </td>
   <td style="text-align:right;"> 1.3630654 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldlc_3 </td>
   <td style="text-align:right;"> 142.23138 </td>
   <td style="text-align:right;"> 20.9505274 </td>
   <td style="text-align:right;"> -0.0362448 </td>
   <td style="text-align:right;"> 0.5333231 </td>
   <td style="text-align:right;"> 0.4053896 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ldlc_m </td>
   <td style="text-align:right;"> 134.70667 </td>
   <td style="text-align:right;"> 39.4337072 </td>
   <td style="text-align:right;"> -0.0554878 </td>
   <td style="text-align:right;"> 0.9732560 </td>
   <td style="text-align:right;"> 0.0011556 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> chol_b </td>
   <td style="text-align:right;"> 178.18662 </td>
   <td style="text-align:right;"> 24.2827204 </td>
   <td style="text-align:right;"> 0.0308232 </td>
   <td style="text-align:right;"> 0.2214450 </td>
   <td style="text-align:right;"> 1.6042662 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> chol_3 </td>
   <td style="text-align:right;"> 155.24378 </td>
   <td style="text-align:right;"> 13.1459087 </td>
   <td style="text-align:right;"> -0.0148308 </td>
   <td style="text-align:right;"> 0.3987975 </td>
   <td style="text-align:right;"> 0.7515610 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> chol_m </td>
   <td style="text-align:right;"> 145.15556 </td>
   <td style="text-align:right;"> 25.1433825 </td>
   <td style="text-align:right;"> -0.0539416 </td>
   <td style="text-align:right;"> 0.8699886 </td>
   <td style="text-align:right;"> 0.0275640 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trig_b </td>
   <td style="text-align:right;"> 175.58765 </td>
   <td style="text-align:right;"> 70.1964478 </td>
   <td style="text-align:right;"> -0.0221870 </td>
   <td style="text-align:right;"> 0.4532869 </td>
   <td style="text-align:right;"> 0.5875964 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trig_3 </td>
   <td style="text-align:right;"> 126.78892 </td>
   <td style="text-align:right;"> 38.2253928 </td>
   <td style="text-align:right;"> -0.0611397 </td>
   <td style="text-align:right;"> 0.8879095 </td>
   <td style="text-align:right;"> 0.0205102 </td>
   <td style="text-align:right;"> 16 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> trig_m </td>
   <td style="text-align:right;"> 124.64444 </td>
   <td style="text-align:right;"> 70.8288325 </td>
   <td style="text-align:right;"> -0.0555555 </td>
   <td style="text-align:right;"> 0.9993468 </td>
   <td style="text-align:right;"> 0.0000007 </td>
   <td style="text-align:right;"> 18 </td>
  </tr>
</tbody>
</table></div>

<br>  
<br>  


<br>


## Supplemental Table 7. Linear regression statistics for plasma lipid levels and HDL protein Log2FC {.tabset}

### Protein ~ CEC


```r
out_df <- data.frame()
 
for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(i ~ CEC, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_CEC.Prot <- out_df
rownames(lm_CEC.Prot) <- rn[,1]
lm_CEC.Prot$p.adj <- p.adjust(lm_CEC.Prot$p, method = "BH")

kable(lm_CEC.Prot, caption="Model: protein FC ~ CEC") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Model: protein FC ~ CEC</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.0779021 </td>
   <td style="text-align:right;"> 0.0391957 </td>
   <td style="text-align:right;"> 0.0494306 </td>
   <td style="text-align:right;"> 0.1755938 </td>
   <td style="text-align:right;"> 1.9880197 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8039673 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.0699900 </td>
   <td style="text-align:right;"> 0.0538070 </td>
   <td style="text-align:right;"> 0.0135398 </td>
   <td style="text-align:right;"> 0.2762563 </td>
   <td style="text-align:right;"> 1.2607876 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8510934 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.0165648 </td>
   <td style="text-align:right;"> 0.0323453 </td>
   <td style="text-align:right;"> -0.0520875 </td>
   <td style="text-align:right;"> 0.8103052 </td>
   <td style="text-align:right;"> 0.0593347 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.0496082 </td>
   <td style="text-align:right;"> 0.0356342 </td>
   <td style="text-align:right;"> -0.0451033 </td>
   <td style="text-align:right;"> 0.6763825 </td>
   <td style="text-align:right;"> 0.1800216 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 0.0038084 </td>
   <td style="text-align:right;"> 0.0413596 </td>
   <td style="text-align:right;"> -0.0554473 </td>
   <td style="text-align:right;"> 0.9661985 </td>
   <td style="text-align:right;"> 0.0018464 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.0394094 </td>
   <td style="text-align:right;"> 0.0801729 </td>
   <td style="text-align:right;"> -0.0756027 </td>
   <td style="text-align:right;"> 0.9014063 </td>
   <td style="text-align:right;"> 0.0159584 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.0092009 </td>
   <td style="text-align:right;"> 0.0161757 </td>
   <td style="text-align:right;"> -0.0348565 </td>
   <td style="text-align:right;"> 0.5559669 </td>
   <td style="text-align:right;"> 0.3600336 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.0096234 </td>
   <td style="text-align:right;"> 0.0337279 </td>
   <td style="text-align:right;"> -0.0553887 </td>
   <td style="text-align:right;"> 0.9580387 </td>
   <td style="text-align:right;"> 0.0028465 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 0.0241574 </td>
   <td style="text-align:right;"> 0.0386544 </td>
   <td style="text-align:right;"> -0.0227962 </td>
   <td style="text-align:right;"> 0.4575068 </td>
   <td style="text-align:right;"> 0.5765251 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9559702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> -0.0262709 </td>
   <td style="text-align:right;"> 0.0454338 </td>
   <td style="text-align:right;"> -0.0375938 </td>
   <td style="text-align:right;"> 0.5835790 </td>
   <td style="text-align:right;"> 0.3115973 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> -0.0241311 </td>
   <td style="text-align:right;"> 0.0264288 </td>
   <td style="text-align:right;"> -0.0218936 </td>
   <td style="text-align:right;"> 0.4512749 </td>
   <td style="text-align:right;"> 0.5929345 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9559702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.0773456 </td>
   <td style="text-align:right;"> 0.0326758 </td>
   <td style="text-align:right;"> 0.1771433 </td>
   <td style="text-align:right;"> 0.0367370 </td>
   <td style="text-align:right;"> 5.0902913 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.4653358 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.0029549 </td>
   <td style="text-align:right;"> 0.0150714 </td>
   <td style="text-align:right;"> -0.0554791 </td>
   <td style="text-align:right;"> 0.9715988 </td>
   <td style="text-align:right;"> 0.0013033 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.0125890 </td>
   <td style="text-align:right;"> 0.0187469 </td>
   <td style="text-align:right;"> -0.0146956 </td>
   <td style="text-align:right;"> 0.4057493 </td>
   <td style="text-align:right;"> 0.7248275 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9559702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.0427291 </td>
   <td style="text-align:right;"> 0.0239833 </td>
   <td style="text-align:right;"> -0.0141270 </td>
   <td style="text-align:right;"> 0.4024356 </td>
   <td style="text-align:right;"> 0.7353267 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9559702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> -0.0014935 </td>
   <td style="text-align:right;"> 0.0221411 </td>
   <td style="text-align:right;"> -0.0462312 </td>
   <td style="text-align:right;"> 0.6934822 </td>
   <td style="text-align:right;"> 0.1604221 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.0014112 </td>
   <td style="text-align:right;"> 0.0220518 </td>
   <td style="text-align:right;"> -0.0552125 </td>
   <td style="text-align:right;"> 0.9398662 </td>
   <td style="text-align:right;"> 0.0058521 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.0418872 </td>
   <td style="text-align:right;"> 0.0257649 </td>
   <td style="text-align:right;"> 0.0379505 </td>
   <td style="text-align:right;"> 0.2024987 </td>
   <td style="text-align:right;"> 1.7495045 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8039673 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.0076901 </td>
   <td style="text-align:right;"> 0.0316581 </td>
   <td style="text-align:right;"> -0.0505627 </td>
   <td style="text-align:right;"> 0.7732624 </td>
   <td style="text-align:right;"> 0.0855460 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0534084 </td>
   <td style="text-align:right;"> 0.0291147 </td>
   <td style="text-align:right;"> 0.0344501 </td>
   <td style="text-align:right;"> 0.2115703 </td>
   <td style="text-align:right;"> 1.6779057 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8039673 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.1638576 </td>
   <td style="text-align:right;"> 0.0431242 </td>
   <td style="text-align:right;"> 0.3242193 </td>
   <td style="text-align:right;"> 0.0051797 </td>
   <td style="text-align:right;"> 10.1156297 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.1968268 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.0819470 </td>
   <td style="text-align:right;"> 0.0749594 </td>
   <td style="text-align:right;"> 0.0085807 </td>
   <td style="text-align:right;"> 0.2947942 </td>
   <td style="text-align:right;"> 1.1644439 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8510934 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.0966391 </td>
   <td style="text-align:right;"> 0.0622115 </td>
   <td style="text-align:right;"> 0.1193409 </td>
   <td style="text-align:right;"> 0.0748733 </td>
   <td style="text-align:right;"> 3.5747504 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.7112967 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 0.0103599 </td>
   <td style="text-align:right;"> 0.0488062 </td>
   <td style="text-align:right;"> -0.0536973 </td>
   <td style="text-align:right;"> 0.8605806 </td>
   <td style="text-align:right;"> 0.0317440 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.0292633 </td>
   <td style="text-align:right;"> 0.0218049 </td>
   <td style="text-align:right;"> 0.0039287 </td>
   <td style="text-align:right;"> 0.3135607 </td>
   <td style="text-align:right;"> 1.0749401 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8510934 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -0.0259670 </td>
   <td style="text-align:right;"> 0.0713428 </td>
   <td style="text-align:right;"> -0.0567283 </td>
   <td style="text-align:right;"> 0.8565044 </td>
   <td style="text-align:right;"> 0.0337066 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0266037 </td>
   <td style="text-align:right;"> 0.0385518 </td>
   <td style="text-align:right;"> -0.0524089 </td>
   <td style="text-align:right;"> 0.8191635 </td>
   <td style="text-align:right;"> 0.0538189 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.0567286 </td>
   <td style="text-align:right;"> 0.0367219 </td>
   <td style="text-align:right;"> 0.0578860 </td>
   <td style="text-align:right;"> 0.1582331 </td>
   <td style="text-align:right;"> 2.1674104 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8039673 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.0038018 </td>
   <td style="text-align:right;"> 0.0204652 </td>
   <td style="text-align:right;"> -0.0550794 </td>
   <td style="text-align:right;"> 0.9291828 </td>
   <td style="text-align:right;"> 0.0081226 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -0.0155150 </td>
   <td style="text-align:right;"> 0.0381950 </td>
   <td style="text-align:right;"> -0.0451839 </td>
   <td style="text-align:right;"> 0.6775677 </td>
   <td style="text-align:right;"> 0.1786199 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.1536059 </td>
   <td style="text-align:right;"> 0.0495389 </td>
   <td style="text-align:right;"> 0.2456271 </td>
   <td style="text-align:right;"> 0.0152621 </td>
   <td style="text-align:right;"> 7.1864823 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.2899802 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -0.0348875 </td>
   <td style="text-align:right;"> 0.0220580 </td>
   <td style="text-align:right;"> 0.0349203 </td>
   <td style="text-align:right;"> 0.2103263 </td>
   <td style="text-align:right;"> 1.6874923 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8039673 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.0096549 </td>
   <td style="text-align:right;"> 0.0144388 </td>
   <td style="text-align:right;"> -0.0544844 </td>
   <td style="text-align:right;"> 0.8939372 </td>
   <td style="text-align:right;"> 0.0182849 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -0.0008887 </td>
   <td style="text-align:right;"> 0.0200831 </td>
   <td style="text-align:right;"> -0.0512740 </td>
   <td style="text-align:right;"> 0.7896597 </td>
   <td style="text-align:right;"> 0.0733083 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.0229026 </td>
   <td style="text-align:right;"> 0.0415556 </td>
   <td style="text-align:right;"> -0.0480433 </td>
   <td style="text-align:right;"> 0.7236273 </td>
   <td style="text-align:right;"> 0.1290221 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9715988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 0.2625887 </td>
   <td style="text-align:right;"> 0.1340848 </td>
   <td style="text-align:right;"> 0.0260088 </td>
   <td style="text-align:right;"> 0.2454236 </td>
   <td style="text-align:right;"> 1.4539571 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.8478270 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.1767261 </td>
   <td style="text-align:right;"> 0.1171787 </td>
   <td style="text-align:right;"> -0.0256350 </td>
   <td style="text-align:right;"> 0.4779851 </td>
   <td style="text-align:right;"> 0.5251089 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9559702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.0590932 </td>
   <td style="text-align:right;"> 0.0407263 </td>
   <td style="text-align:right;"> 0.0507315 </td>
   <td style="text-align:right;"> 0.1727963 </td>
   <td style="text-align:right;"> 2.0154122 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8039673 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_CEC.Prot, file="lm_CEC_prot.csv")
```

### Protein ~ HDLC


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(i ~ HDLC, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_HDLC.Prot <- out_df
rownames(lm_HDLC.Prot) <- rn[,1]
lm_HDLC.Prot$p.adj <- p.adjust(lm_HDLC.Prot$p, method = "BH")

kable(lm_HDLC.Prot, caption="Model: protein FC ~ HDLC") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Model: protein FC ~ HDLC</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> -0.0201788 </td>
   <td style="text-align:right;"> 0.0386085 </td>
   <td style="text-align:right;"> 0.0213724 </td>
   <td style="text-align:right;"> 0.2496919 </td>
   <td style="text-align:right;"> 1.4149449 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3953455 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> -0.0832703 </td>
   <td style="text-align:right;"> 0.0488755 </td>
   <td style="text-align:right;"> 0.1363653 </td>
   <td style="text-align:right;"> 0.0608202 </td>
   <td style="text-align:right;"> 4.0000430 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> -0.0475301 </td>
   <td style="text-align:right;"> 0.0282915 </td>
   <td style="text-align:right;"> 0.1459415 </td>
   <td style="text-align:right;"> 0.0540785 </td>
   <td style="text-align:right;"> 4.2467204 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.0871811 </td>
   <td style="text-align:right;"> 0.0323568 </td>
   <td style="text-align:right;"> 0.0856739 </td>
   <td style="text-align:right;"> 0.1127320 </td>
   <td style="text-align:right;"> 2.7803326 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 0.0334864 </td>
   <td style="text-align:right;"> 0.0395729 </td>
   <td style="text-align:right;"> -0.0252407 </td>
   <td style="text-align:right;"> 0.4750565 </td>
   <td style="text-align:right;"> 0.5322340 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6943134 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.0220041 </td>
   <td style="text-align:right;"> 0.0985482 </td>
   <td style="text-align:right;"> -0.0764023 </td>
   <td style="text-align:right;"> 0.9379984 </td>
   <td style="text-align:right;"> 0.0062892 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.9379984 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> -0.0032976 </td>
   <td style="text-align:right;"> 0.0158418 </td>
   <td style="text-align:right;"> -0.0531924 </td>
   <td style="text-align:right;"> 0.8429757 </td>
   <td style="text-align:right;"> 0.0403888 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9152308 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.0051558 </td>
   <td style="text-align:right;"> 0.0327387 </td>
   <td style="text-align:right;"> -0.0551234 </td>
   <td style="text-align:right;"> 0.9325205 </td>
   <td style="text-align:right;"> 0.0073730 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9379984 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -0.0514782 </td>
   <td style="text-align:right;"> 0.0363519 </td>
   <td style="text-align:right;"> 0.0401805 </td>
   <td style="text-align:right;"> 0.1969421 </td>
   <td style="text-align:right;"> 1.7953895 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3563715 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> -0.0695093 </td>
   <td style="text-align:right;"> 0.0413081 </td>
   <td style="text-align:right;"> 0.0899085 </td>
   <td style="text-align:right;"> 0.1070780 </td>
   <td style="text-align:right;"> 2.8770225 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> -0.0552428 </td>
   <td style="text-align:right;"> 0.0229503 </td>
   <td style="text-align:right;"> 0.1823442 </td>
   <td style="text-align:right;"> 0.0344204 </td>
   <td style="text-align:right;"> 5.2371626 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> -0.0110141 </td>
   <td style="text-align:right;"> 0.0356933 </td>
   <td style="text-align:right;"> -0.0418125 </td>
   <td style="text-align:right;"> 0.6319377 </td>
   <td style="text-align:right;"> 0.2374469 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8004544 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.0013342 </td>
   <td style="text-align:right;"> 0.0146224 </td>
   <td style="text-align:right;"> -0.0542070 </td>
   <td style="text-align:right;"> 0.8810760 </td>
   <td style="text-align:right;"> 0.0230265 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9300247 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> -0.0306107 </td>
   <td style="text-align:right;"> 0.0172926 </td>
   <td style="text-align:right;"> 0.0838992 </td>
   <td style="text-align:right;"> 0.1151904 </td>
   <td style="text-align:right;"> 2.7400757 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.0581788 </td>
   <td style="text-align:right;"> 0.0221353 </td>
   <td style="text-align:right;"> 0.0833746 </td>
   <td style="text-align:right;"> 0.1159275 </td>
   <td style="text-align:right;"> 2.7282053 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> -0.0355506 </td>
   <td style="text-align:right;"> 0.0189285 </td>
   <td style="text-align:right;"> 0.1886539 </td>
   <td style="text-align:right;"> 0.0317955 </td>
   <td style="text-align:right;"> 5.4178724 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -0.0255403 </td>
   <td style="text-align:right;"> 0.0205053 </td>
   <td style="text-align:right;"> 0.0318811 </td>
   <td style="text-align:right;"> 0.2185129 </td>
   <td style="text-align:right;"> 1.6256888 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3774314 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> -0.0271609 </td>
   <td style="text-align:right;"> 0.0246917 </td>
   <td style="text-align:right;"> 0.0624597 </td>
   <td style="text-align:right;"> 0.1496045 </td>
   <td style="text-align:right;"> 2.2657966 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3511419 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> -0.0525070 </td>
   <td style="text-align:right;"> 0.0281711 </td>
   <td style="text-align:right;"> 0.1173166 </td>
   <td style="text-align:right;"> 0.0767420 </td>
   <td style="text-align:right;"> 3.5252716 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> -0.0207994 </td>
   <td style="text-align:right;"> 0.0280975 </td>
   <td style="text-align:right;"> 0.0458149 </td>
   <td style="text-align:right;"> 0.1836263 </td>
   <td style="text-align:right;"> 1.9122791 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3511419 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.0196052 </td>
   <td style="text-align:right;"> 0.0522512 </td>
   <td style="text-align:right;"> -0.0526922 </td>
   <td style="text-align:right;"> 0.8273727 </td>
   <td style="text-align:right;"> 0.0489605 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9152308 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -0.1314570 </td>
   <td style="text-align:right;"> 0.0675217 </td>
   <td style="text-align:right;"> 0.1464348 </td>
   <td style="text-align:right;"> 0.0537515 </td>
   <td style="text-align:right;"> 4.2595770 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> -0.1048661 </td>
   <td style="text-align:right;"> 0.0626041 </td>
   <td style="text-align:right;"> 0.0537263 </td>
   <td style="text-align:right;"> 0.1665348 </td>
   <td style="text-align:right;"> 2.0787582 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3511419 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.0411157 </td>
   <td style="text-align:right;"> 0.0462507 </td>
   <td style="text-align:right;"> -0.0040332 </td>
   <td style="text-align:right;"> 0.3492451 </td>
   <td style="text-align:right;"> 0.9236770 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.5308525 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.0388759 </td>
   <td style="text-align:right;"> 0.0203932 </td>
   <td style="text-align:right;"> 0.0755167 </td>
   <td style="text-align:right;"> 0.1275594 </td>
   <td style="text-align:right;"> 2.5520207 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 0.0052287 </td>
   <td style="text-align:right;"> 0.0680526 </td>
   <td style="text-align:right;"> -0.0539213 </td>
   <td style="text-align:right;"> 0.7819492 </td>
   <td style="text-align:right;"> 0.0790747 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.9004263 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> -0.0390296 </td>
   <td style="text-align:right;"> 0.0348231 </td>
   <td style="text-align:right;"> 0.0888809 </td>
   <td style="text-align:right;"> 0.1084233 </td>
   <td style="text-align:right;"> 2.8534747 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.0187182 </td>
   <td style="text-align:right;"> 0.0375707 </td>
   <td style="text-align:right;"> -0.0464006 </td>
   <td style="text-align:right;"> 0.6961509 </td>
   <td style="text-align:right;"> 0.1574822 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8266792 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.0308242 </td>
   <td style="text-align:right;"> 0.0188988 </td>
   <td style="text-align:right;"> 0.0452958 </td>
   <td style="text-align:right;"> 0.1848115 </td>
   <td style="text-align:right;"> 1.9014516 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3511419 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 0.0161327 </td>
   <td style="text-align:right;"> 0.0370593 </td>
   <td style="text-align:right;"> -0.0440413 </td>
   <td style="text-align:right;"> 0.6612358 </td>
   <td style="text-align:right;"> 0.1985134 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.8105471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> -0.0043681 </td>
   <td style="text-align:right;"> 0.0564440 </td>
   <td style="text-align:right;"> -0.0391425 </td>
   <td style="text-align:right;"> 0.6004146 </td>
   <td style="text-align:right;"> 0.2843066 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.7965301 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -0.0432881 </td>
   <td style="text-align:right;"> 0.0205887 </td>
   <td style="text-align:right;"> 0.1078619 </td>
   <td style="text-align:right;"> 0.0860956 </td>
   <td style="text-align:right;"> 3.2971517 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> -0.0081322 </td>
   <td style="text-align:right;"> 0.0134800 </td>
   <td style="text-align:right;"> 0.0247751 </td>
   <td style="text-align:right;"> 0.2390774 </td>
   <td style="text-align:right;"> 1.4826846 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3949974 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -0.0054182 </td>
   <td style="text-align:right;"> 0.0193898 </td>
   <td style="text-align:right;"> -0.0398012 </td>
   <td style="text-align:right;"> 0.6078783 </td>
   <td style="text-align:right;"> 0.2727238 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.7965301 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> -0.0746266 </td>
   <td style="text-align:right;"> 0.0350776 </td>
   <td style="text-align:right;"> 0.2076371 </td>
   <td style="text-align:right;"> 0.0249909 </td>
   <td style="text-align:right;"> 5.9789112 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 0.2659861 </td>
   <td style="text-align:right;"> 0.1165678 </td>
   <td style="text-align:right;"> 0.0576632 </td>
   <td style="text-align:right;"> 0.1724136 </td>
   <td style="text-align:right;"> 2.0402600 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.3511419 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.1575315 </td>
   <td style="text-align:right;"> 0.1143744 </td>
   <td style="text-align:right;"> -0.0368077 </td>
   <td style="text-align:right;"> 0.5753851 </td>
   <td style="text-align:right;"> 0.3254819 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.7965301 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> -0.0569407 </td>
   <td style="text-align:right;"> 0.0390545 </td>
   <td style="text-align:right;"> 0.0737535 </td>
   <td style="text-align:right;"> 0.1303292 </td>
   <td style="text-align:right;"> 2.5128980 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3301674 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_HDLC.Prot, file="lm_HDLC_prot.csv")
```

### Protein ~ LDLC


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(i ~ LDLC, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_LDLC.Prot <- out_df
rownames(lm_LDLC.Prot) <- rn[,1]
lm_LDLC.Prot$p.adj <- p.adjust(lm_LDLC.Prot$p, method = "BH")

kable(lm_LDLC.Prot, caption="Model: protein FC ~ LDLC") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Model: protein FC ~ LDLC</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.0356560 </td>
   <td style="text-align:right;"> 0.0262373 </td>
   <td style="text-align:right;"> -0.0428495 </td>
   <td style="text-align:right;"> 0.6451859 </td>
   <td style="text-align:right;"> 0.2193116 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.0356063 </td>
   <td style="text-align:right;"> 0.0350317 </td>
   <td style="text-align:right;"> -0.0237694 </td>
   <td style="text-align:right;"> 0.4643724 </td>
   <td style="text-align:right;"> 0.5588667 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.0169871 </td>
   <td style="text-align:right;"> 0.0206034 </td>
   <td style="text-align:right;"> -0.0451582 </td>
   <td style="text-align:right;"> 0.6771894 </td>
   <td style="text-align:right;"> 0.1790666 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.0200022 </td>
   <td style="text-align:right;"> 0.0225684 </td>
   <td style="text-align:right;"> -0.0263643 </td>
   <td style="text-align:right;"> 0.4834766 </td>
   <td style="text-align:right;"> 0.5119458 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> -0.0209638 </td>
   <td style="text-align:right;"> 0.0255573 </td>
   <td style="text-align:right;"> 0.0132877 </td>
   <td style="text-align:right;"> 0.2771647 </td>
   <td style="text-align:right;"> 1.2558665 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.0481044 </td>
   <td style="text-align:right;"> 0.0785289 </td>
   <td style="text-align:right;"> -0.0720100 </td>
   <td style="text-align:right;"> 0.8109708 </td>
   <td style="text-align:right;"> 0.0595801 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.0059977 </td>
   <td style="text-align:right;"> 0.0103200 </td>
   <td style="text-align:right;"> -0.0312943 </td>
   <td style="text-align:right;"> 0.5234407 </td>
   <td style="text-align:right;"> 0.4234507 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.0126770 </td>
   <td style="text-align:right;"> 0.0215221 </td>
   <td style="text-align:right;"> -0.0521593 </td>
   <td style="text-align:right;"> 0.8122461 </td>
   <td style="text-align:right;"> 0.0581017 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -0.0105100 </td>
   <td style="text-align:right;"> 0.0250465 </td>
   <td style="text-align:right;"> -0.0513811 </td>
   <td style="text-align:right;"> 0.7922495 </td>
   <td style="text-align:right;"> 0.0714686 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.0416785 </td>
   <td style="text-align:right;"> 0.0271245 </td>
   <td style="text-align:right;"> 0.0945414 </td>
   <td style="text-align:right;"> 0.1012185 </td>
   <td style="text-align:right;"> 2.9838427 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6718986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> -0.0014635 </td>
   <td style="text-align:right;"> 0.0171502 </td>
   <td style="text-align:right;"> -0.0535663 </td>
   <td style="text-align:right;"> 0.8557976 </td>
   <td style="text-align:right;"> 0.0339860 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.0087796 </td>
   <td style="text-align:right;"> 0.0236400 </td>
   <td style="text-align:right;"> -0.0544875 </td>
   <td style="text-align:right;"> 0.8940917 </td>
   <td style="text-align:right;"> 0.0182314 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -0.0106196 </td>
   <td style="text-align:right;"> 0.0089371 </td>
   <td style="text-align:right;"> 0.0913140 </td>
   <td style="text-align:right;"> 0.1052652 </td>
   <td style="text-align:right;"> 2.9093140 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6718986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.0067822 </td>
   <td style="text-align:right;"> 0.0119685 </td>
   <td style="text-align:right;"> -0.0125879 </td>
   <td style="text-align:right;"> 0.3936519 </td>
   <td style="text-align:right;"> 0.7638037 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.0221217 </td>
   <td style="text-align:right;"> 0.0156365 </td>
   <td style="text-align:right;"> -0.0554349 </td>
   <td style="text-align:right;"> 0.9643253 </td>
   <td style="text-align:right;"> 0.0020569 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 0.0220569 </td>
   <td style="text-align:right;"> 0.0136934 </td>
   <td style="text-align:right;"> 0.0202252 </td>
   <td style="text-align:right;"> 0.2533916 </td>
   <td style="text-align:right;"> 1.3922109 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.0066617 </td>
   <td style="text-align:right;"> 0.0139857 </td>
   <td style="text-align:right;"> -0.0391979 </td>
   <td style="text-align:right;"> 0.6010351 </td>
   <td style="text-align:right;"> 0.2833319 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.0283302 </td>
   <td style="text-align:right;"> 0.0165164 </td>
   <td style="text-align:right;"> 0.0320580 </td>
   <td style="text-align:right;"> 0.2180269 </td>
   <td style="text-align:right;"> 1.6292750 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.0003206 </td>
   <td style="text-align:right;"> 0.0202762 </td>
   <td style="text-align:right;"> -0.0551220 </td>
   <td style="text-align:right;"> 0.9324118 </td>
   <td style="text-align:right;"> 0.0073969 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0344130 </td>
   <td style="text-align:right;"> 0.0189315 </td>
   <td style="text-align:right;"> 0.0004631 </td>
   <td style="text-align:right;"> 0.3285017 </td>
   <td style="text-align:right;"> 1.0088032 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.0897848 </td>
   <td style="text-align:right;"> 0.0309932 </td>
   <td style="text-align:right;"> 0.1453785 </td>
   <td style="text-align:right;"> 0.0544540 </td>
   <td style="text-align:right;"> 4.2320649 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6718986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.0517620 </td>
   <td style="text-align:right;"> 0.0478757 </td>
   <td style="text-align:right;"> 0.0098256 </td>
   <td style="text-align:right;"> 0.2900046 </td>
   <td style="text-align:right;"> 1.1885380 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.0106171 </td>
   <td style="text-align:right;"> 0.0429359 </td>
   <td style="text-align:right;"> -0.0270283 </td>
   <td style="text-align:right;"> 0.4885641 </td>
   <td style="text-align:right;"> 0.4999773 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.0436825 </td>
   <td style="text-align:right;"> 0.0289760 </td>
   <td style="text-align:right;"> 0.0906722 </td>
   <td style="text-align:right;"> 0.1060892 </td>
   <td style="text-align:right;"> 2.8945542 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6718986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.0041576 </td>
   <td style="text-align:right;"> 0.0143225 </td>
   <td style="text-align:right;"> -0.0522010 </td>
   <td style="text-align:right;"> 0.8133821 </td>
   <td style="text-align:right;"> 0.0573865 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -0.0065644 </td>
   <td style="text-align:right;"> 0.0453032 </td>
   <td style="text-align:right;"> -0.0572237 </td>
   <td style="text-align:right;"> 0.8744659 </td>
   <td style="text-align:right;"> 0.0257243 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0218790 </td>
   <td style="text-align:right;"> 0.0246543 </td>
   <td style="text-align:right;"> -0.0537982 </td>
   <td style="text-align:right;"> 0.8643839 </td>
   <td style="text-align:right;"> 0.0300175 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.0096709 </td>
   <td style="text-align:right;"> 0.0248034 </td>
   <td style="text-align:right;"> -0.0523308 </td>
   <td style="text-align:right;"> 0.8169688 </td>
   <td style="text-align:right;"> 0.0551590 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.0033631 </td>
   <td style="text-align:right;"> 0.0130697 </td>
   <td style="text-align:right;"> -0.0535681 </td>
   <td style="text-align:right;"> 0.8558633 </td>
   <td style="text-align:right;"> 0.0339547 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -0.0462463 </td>
   <td style="text-align:right;"> 0.0214956 </td>
   <td style="text-align:right;"> 0.1894973 </td>
   <td style="text-align:right;"> 0.0314594 </td>
   <td style="text-align:right;"> 5.4422410 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6718986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.0924692 </td>
   <td style="text-align:right;"> 0.0332422 </td>
   <td style="text-align:right;"> 0.1683331 </td>
   <td style="text-align:right;"> 0.0410029 </td>
   <td style="text-align:right;"> 4.8456857 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6718986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -0.0137679 </td>
   <td style="text-align:right;"> 0.0146445 </td>
   <td style="text-align:right;"> -0.0414964 </td>
   <td style="text-align:right;"> 0.6280196 </td>
   <td style="text-align:right;"> 0.2429826 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.0033247 </td>
   <td style="text-align:right;"> 0.0091629 </td>
   <td style="text-align:right;"> -0.0397360 </td>
   <td style="text-align:right;"> 0.6071315 </td>
   <td style="text-align:right;"> 0.2738687 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -0.0004125 </td>
   <td style="text-align:right;"> 0.0128027 </td>
   <td style="text-align:right;"> -0.0460157 </td>
   <td style="text-align:right;"> 0.6901271 </td>
   <td style="text-align:right;"> 0.1641632 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.0265429 </td>
   <td style="text-align:right;"> 0.0262499 </td>
   <td style="text-align:right;"> -0.0238831 </td>
   <td style="text-align:right;"> 0.4651847 </td>
   <td style="text-align:right;"> 0.5568060 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 0.1450677 </td>
   <td style="text-align:right;"> 0.0812280 </td>
   <td style="text-align:right;"> -0.0439621 </td>
   <td style="text-align:right;"> 0.6013416 </td>
   <td style="text-align:right;"> 0.2841166 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.0971948 </td>
   <td style="text-align:right;"> 0.0759688 </td>
   <td style="text-align:right;"> -0.0554595 </td>
   <td style="text-align:right;"> 0.9681671 </td>
   <td style="text-align:right;"> 0.0016375 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.0232381 </td>
   <td style="text-align:right;"> 0.0269577 </td>
   <td style="text-align:right;"> -0.0183148 </td>
   <td style="text-align:right;"> 0.4277679 </td>
   <td style="text-align:right;"> 0.6582766 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9681671 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_LDLC.Prot, file="lm_LDLC_prot.csv")
```

### Protein ~ Cholesterol


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(i ~ Cholesterol, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Chol.Prot <- out_df
rownames(lm_Chol.Prot) <- rn[,1]
lm_Chol.Prot$p.adj <- p.adjust(lm_Chol.Prot$p, method = "BH")

kable(lm_Chol.Prot, caption="Model: protein FC ~ Cholesterol") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Model: protein FC ~ Cholesterol</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.0309864 </td>
   <td style="text-align:right;"> 0.0438915 </td>
   <td style="text-align:right;"> -0.0541171 </td>
   <td style="text-align:right;"> 0.8772047 </td>
   <td style="text-align:right;"> 0.0245633 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.0339686 </td>
   <td style="text-align:right;"> 0.0589340 </td>
   <td style="text-align:right;"> -0.0465438 </td>
   <td style="text-align:right;"> 0.6984293 </td>
   <td style="text-align:right;"> 0.1549974 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.0121426 </td>
   <td style="text-align:right;"> 0.0344432 </td>
   <td style="text-align:right;"> -0.0550147 </td>
   <td style="text-align:right;"> 0.9245350 </td>
   <td style="text-align:right;"> 0.0092274 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.0195548 </td>
   <td style="text-align:right;"> 0.0378991 </td>
   <td style="text-align:right;"> -0.0454572 </td>
   <td style="text-align:right;"> 0.6816283 </td>
   <td style="text-align:right;"> 0.1738670 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> -0.0161707 </td>
   <td style="text-align:right;"> 0.0436683 </td>
   <td style="text-align:right;"> -0.0404928 </td>
   <td style="text-align:right;"> 0.6159262 </td>
   <td style="text-align:right;"> 0.2605784 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.0266968 </td>
   <td style="text-align:right;"> 0.1064127 </td>
   <td style="text-align:right;"> -0.0768572 </td>
   <td style="text-align:right;"> 0.9779298 </td>
   <td style="text-align:right;"> 0.0007953 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.9779298 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.0113707 </td>
   <td style="text-align:right;"> 0.0171442 </td>
   <td style="text-align:right;"> -0.0280408 </td>
   <td style="text-align:right;"> 0.4964888 </td>
   <td style="text-align:right;"> 0.4817567 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.0176353 </td>
   <td style="text-align:right;"> 0.0357905 </td>
   <td style="text-align:right;"> -0.0509773 </td>
   <td style="text-align:right;"> 0.7826549 </td>
   <td style="text-align:right;"> 0.0784111 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -0.0324890 </td>
   <td style="text-align:right;"> 0.0411967 </td>
   <td style="text-align:right;"> -0.0273949 </td>
   <td style="text-align:right;"> 0.4914097 </td>
   <td style="text-align:right;"> 0.4933759 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.0580168 </td>
   <td style="text-align:right;"> 0.0465384 </td>
   <td style="text-align:right;"> 0.0372473 </td>
   <td style="text-align:right;"> 0.2042865 </td>
   <td style="text-align:right;"> 1.7350779 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> -0.0157060 </td>
   <td style="text-align:right;"> 0.0284313 </td>
   <td style="text-align:right;"> -0.0458432 </td>
   <td style="text-align:right;"> 0.6874720 </td>
   <td style="text-align:right;"> 0.1671596 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 0.0206517 </td>
   <td style="text-align:right;"> 0.0391905 </td>
   <td style="text-align:right;"> -0.0467767 </td>
   <td style="text-align:right;"> 0.7021790 </td>
   <td style="text-align:right;"> 0.1509578 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -0.0192057 </td>
   <td style="text-align:right;"> 0.0150567 </td>
   <td style="text-align:right;"> 0.0684163 </td>
   <td style="text-align:right;"> 0.1390963 </td>
   <td style="text-align:right;"> 2.3953771 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.0122269 </td>
   <td style="text-align:right;"> 0.0200002 </td>
   <td style="text-align:right;"> -0.0213392 </td>
   <td style="text-align:right;"> 0.4475104 </td>
   <td style="text-align:right;"> 0.6030270 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.0338146 </td>
   <td style="text-align:right;"> 0.0258815 </td>
   <td style="text-align:right;"> -0.0444245 </td>
   <td style="text-align:right;"> 0.6666035 </td>
   <td style="text-align:right;"> 0.1918362 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 0.0214057 </td>
   <td style="text-align:right;"> 0.0233945 </td>
   <td style="text-align:right;"> -0.0329491 </td>
   <td style="text-align:right;"> 0.5381191 </td>
   <td style="text-align:right;"> 0.3939366 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.0110104 </td>
   <td style="text-align:right;"> 0.0232950 </td>
   <td style="text-align:right;"> -0.0413575 </td>
   <td style="text-align:right;"> 0.6263154 </td>
   <td style="text-align:right;"> 0.2454154 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.0237802 </td>
   <td style="text-align:right;"> 0.0284692 </td>
   <td style="text-align:right;"> -0.0387599 </td>
   <td style="text-align:right;"> 0.5961656 </td>
   <td style="text-align:right;"> 0.2910406 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.0020880 </td>
   <td style="text-align:right;"> 0.0337346 </td>
   <td style="text-align:right;"> -0.0549346 </td>
   <td style="text-align:right;"> 0.9191572 </td>
   <td style="text-align:right;"> 0.0105944 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0516383 </td>
   <td style="text-align:right;"> 0.0312543 </td>
   <td style="text-align:right;"> 0.0160072 </td>
   <td style="text-align:right;"> 0.2675481 </td>
   <td style="text-align:right;"> 1.3090852 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.1543161 </td>
   <td style="text-align:right;"> 0.0488934 </td>
   <td style="text-align:right;"> 0.2317773 </td>
   <td style="text-align:right;"> 0.0183025 </td>
   <td style="text-align:right;"> 6.7324102 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3477483 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.0678416 </td>
   <td style="text-align:right;"> 0.0807587 </td>
   <td style="text-align:right;"> -0.0176691 </td>
   <td style="text-align:right;"> 0.4237170 </td>
   <td style="text-align:right;"> 0.6701154 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.0347716 </td>
   <td style="text-align:right;"> 0.0713135 </td>
   <td style="text-align:right;"> -0.0233703 </td>
   <td style="text-align:right;"> 0.4615379 </td>
   <td style="text-align:right;"> 0.5661044 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.0857394 </td>
   <td style="text-align:right;"> 0.0474123 </td>
   <td style="text-align:right;"> 0.1206325 </td>
   <td style="text-align:right;"> 0.0737044 </td>
   <td style="text-align:right;"> 3.6064398 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9335889 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.0153616 </td>
   <td style="text-align:right;"> 0.0237899 </td>
   <td style="text-align:right;"> -0.0485544 </td>
   <td style="text-align:right;"> 0.7328530 </td>
   <td style="text-align:right;"> 0.1201857 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 0.0254078 </td>
   <td style="text-align:right;"> 0.0743213 </td>
   <td style="text-align:right;"> -0.0412809 </td>
   <td style="text-align:right;"> 0.5994668 </td>
   <td style="text-align:right;"> 0.2864021 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0104952 </td>
   <td style="text-align:right;"> 0.0410171 </td>
   <td style="text-align:right;"> -0.0535283 </td>
   <td style="text-align:right;"> 0.8544414 </td>
   <td style="text-align:right;"> 0.0346362 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 0.0452826 </td>
   <td style="text-align:right;"> 0.0401240 </td>
   <td style="text-align:right;"> 0.0053155 </td>
   <td style="text-align:right;"> 0.3078178 </td>
   <td style="text-align:right;"> 1.1015341 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.0076837 </td>
   <td style="text-align:right;"> 0.0217613 </td>
   <td style="text-align:right;"> -0.0549838 </td>
   <td style="text-align:right;"> 0.9224162 </td>
   <td style="text-align:right;"> 0.0097546 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -0.0613655 </td>
   <td style="text-align:right;"> 0.0380036 </td>
   <td style="text-align:right;"> 0.0849363 </td>
   <td style="text-align:right;"> 0.1137472 </td>
   <td style="text-align:right;"> 2.7635823 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.1613506 </td>
   <td style="text-align:right;"> 0.0528181 </td>
   <td style="text-align:right;"> 0.2416294 </td>
   <td style="text-align:right;"> 0.0160876 </td>
   <td style="text-align:right;"> 7.0537140 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.3477483 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -0.0287797 </td>
   <td style="text-align:right;"> 0.0239624 </td>
   <td style="text-align:right;"> -0.0071916 </td>
   <td style="text-align:right;"> 0.3648297 </td>
   <td style="text-align:right;"> 0.8643349 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> -0.0059185 </td>
   <td style="text-align:right;"> 0.0150005 </td>
   <td style="text-align:right;"> -0.0065007 </td>
   <td style="text-align:right;"> 0.3613444 </td>
   <td style="text-align:right;"> 0.8772841 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> -0.0102774 </td>
   <td style="text-align:right;"> 0.0211029 </td>
   <td style="text-align:right;"> -0.0265072 </td>
   <td style="text-align:right;"> 0.4845641 </td>
   <td style="text-align:right;"> 0.5093693 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.0191487 </td>
   <td style="text-align:right;"> 0.0442715 </td>
   <td style="text-align:right;"> -0.0519363 </td>
   <td style="text-align:right;"> 0.8062871 </td>
   <td style="text-align:right;"> 0.0619309 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 0.2108644 </td>
   <td style="text-align:right;"> 0.1328598 </td>
   <td style="text-align:right;"> -0.0198253 </td>
   <td style="text-align:right;"> 0.4252376 </td>
   <td style="text-align:right;"> 0.6695211 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.0827745 </td>
   <td style="text-align:right;"> 0.1263789 </td>
   <td style="text-align:right;"> -0.0550367 </td>
   <td style="text-align:right;"> 0.9260830 </td>
   <td style="text-align:right;"> 0.0088515 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.0303569 </td>
   <td style="text-align:right;"> 0.0451851 </td>
   <td style="text-align:right;"> -0.0333622 </td>
   <td style="text-align:right;"> 0.5418967 </td>
   <td style="text-align:right;"> 0.3865829 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9511123 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Chol.Prot, file="lm_Chol_prot.csv")
```

### Protein ~ Triglycerides


```r
out_df <- data.frame()

for (i in lfc_df[,25:ncol(lfc_df)]){
  fit <- lm(i ~ Triglycerides, data = lfc_df)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

lm_Tri.Prot <- out_df
rownames(lm_Tri.Prot) <- rn[,1]
lm_Tri.Prot$p.adj <- p.adjust(lm_Tri.Prot$p, method = "BH")

kable(lm_Tri.Prot, caption="Model: protein FC ~ Triglycerides") %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
<caption>Model: protein FC ~ Triglycerides</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.0294655 </td>
   <td style="text-align:right;"> 0.0161079 </td>
   <td style="text-align:right;"> -0.0458212 </td>
   <td style="text-align:right;"> 0.6871355 </td>
   <td style="text-align:right;"> 0.1675416 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.0180462 </td>
   <td style="text-align:right;"> 0.0217153 </td>
   <td style="text-align:right;"> -0.0466624 </td>
   <td style="text-align:right;"> 0.7003322 </td>
   <td style="text-align:right;"> 0.1529397 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.0090616 </td>
   <td style="text-align:right;"> 0.0126936 </td>
   <td style="text-align:right;"> -0.0555439 </td>
   <td style="text-align:right;"> 0.9889268 </td>
   <td style="text-align:right;"> 0.0001980 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.0305318 </td>
   <td style="text-align:right;"> 0.0139685 </td>
   <td style="text-align:right;"> -0.0461625 </td>
   <td style="text-align:right;"> 0.6924082 </td>
   <td style="text-align:right;"> 0.1616141 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> -0.0165170 </td>
   <td style="text-align:right;"> 0.0147919 </td>
   <td style="text-align:right;"> 0.1205552 </td>
   <td style="text-align:right;"> 0.0737739 </td>
   <td style="text-align:right;"> 3.6045393 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6133293 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> 0.0780910 </td>
   <td style="text-align:right;"> 0.0524074 </td>
   <td style="text-align:right;"> 0.0004454 </td>
   <td style="text-align:right;"> 0.3341126 </td>
   <td style="text-align:right;"> 1.0062379 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 0.0021315 </td>
   <td style="text-align:right;"> 0.0063618 </td>
   <td style="text-align:right;"> -0.0427717 </td>
   <td style="text-align:right;"> 0.6441700 </td>
   <td style="text-align:right;"> 0.2206708 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> -0.0092046 </td>
   <td style="text-align:right;"> 0.0132096 </td>
   <td style="text-align:right;"> -0.0545965 </td>
   <td style="text-align:right;"> 0.8996121 </td>
   <td style="text-align:right;"> 0.0163696 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> -0.0066410 </td>
   <td style="text-align:right;"> 0.0153695 </td>
   <td style="text-align:right;"> -0.0533822 </td>
   <td style="text-align:right;"> 0.8493417 </td>
   <td style="text-align:right;"> 0.0371381 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.0177875 </td>
   <td style="text-align:right;"> 0.0169778 </td>
   <td style="text-align:right;"> 0.0561378 </td>
   <td style="text-align:right;"> 0.1616676 </td>
   <td style="text-align:right;"> 2.1300566 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> -0.0030451 </td>
   <td style="text-align:right;"> 0.0105162 </td>
   <td style="text-align:right;"> -0.0540056 </td>
   <td style="text-align:right;"> 0.8725695 </td>
   <td style="text-align:right;"> 0.0264704 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> -0.0006567 </td>
   <td style="text-align:right;"> 0.0143688 </td>
   <td style="text-align:right;"> -0.0365370 </td>
   <td style="text-align:right;"> 0.5726152 </td>
   <td style="text-align:right;"> 0.3302673 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> -0.0053435 </td>
   <td style="text-align:right;"> 0.0052773 </td>
   <td style="text-align:right;"> 0.1569617 </td>
   <td style="text-align:right;"> 0.0472102 </td>
   <td style="text-align:right;"> 4.5375282 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6133293 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> -0.0016916 </td>
   <td style="text-align:right;"> 0.0074827 </td>
   <td style="text-align:right;"> -0.0531026 </td>
   <td style="text-align:right;"> 0.8400577 </td>
   <td style="text-align:right;"> 0.0419264 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> 0.0221746 </td>
   <td style="text-align:right;"> 0.0095849 </td>
   <td style="text-align:right;"> -0.0551790 </td>
   <td style="text-align:right;"> 0.9370047 </td>
   <td style="text-align:right;"> 0.0064235 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 0.0137767 </td>
   <td style="text-align:right;"> 0.0084805 </td>
   <td style="text-align:right;"> 0.0001372 </td>
   <td style="text-align:right;"> 0.3299519 </td>
   <td style="text-align:right;"> 1.0026076 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> -0.0005936 </td>
   <td style="text-align:right;"> 0.0086406 </td>
   <td style="text-align:right;"> -0.0553828 </td>
   <td style="text-align:right;"> 0.9573037 </td>
   <td style="text-align:right;"> 0.0029472 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.0229094 </td>
   <td style="text-align:right;"> 0.0096919 </td>
   <td style="text-align:right;"> 0.1131831 </td>
   <td style="text-align:right;"> 0.0807012 </td>
   <td style="text-align:right;"> 3.4249413 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6133293 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> -0.0096693 </td>
   <td style="text-align:right;"> 0.0121761 </td>
   <td style="text-align:right;"> -0.0123744 </td>
   <td style="text-align:right;"> 0.3924547 </td>
   <td style="text-align:right;"> 0.7677595 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.0173666 </td>
   <td style="text-align:right;"> 0.0119258 </td>
   <td style="text-align:right;"> -0.0553573 </td>
   <td style="text-align:right;"> 0.9542696 </td>
   <td style="text-align:right;"> 0.0033815 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 0.0490938 </td>
   <td style="text-align:right;"> 0.0203898 </td>
   <td style="text-align:right;"> 0.0158416 </td>
   <td style="text-align:right;"> 0.2681224 </td>
   <td style="text-align:right;"> 1.3058358 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> 0.0159183 </td>
   <td style="text-align:right;"> 0.0300730 </td>
   <td style="text-align:right;"> -0.0395193 </td>
   <td style="text-align:right;"> 0.6046608 </td>
   <td style="text-align:right;"> 0.2776789 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> -0.0195974 </td>
   <td style="text-align:right;"> 0.0266775 </td>
   <td style="text-align:right;"> -0.0549504 </td>
   <td style="text-align:right;"> 0.9201883 </td>
   <td style="text-align:right;"> 0.0103249 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> -0.0253704 </td>
   <td style="text-align:right;"> 0.0172963 </td>
   <td style="text-align:right;"> 0.1379170 </td>
   <td style="text-align:right;"> 0.0596752 </td>
   <td style="text-align:right;"> 4.0396411 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.6133293 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> -0.0039960 </td>
   <td style="text-align:right;"> 0.0087373 </td>
   <td style="text-align:right;"> -0.0418695 </td>
   <td style="text-align:right;"> 0.6326498 </td>
   <td style="text-align:right;"> 0.2364493 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> -0.0209425 </td>
   <td style="text-align:right;"> 0.0280072 </td>
   <td style="text-align:right;"> -0.0512458 </td>
   <td style="text-align:right;"> 0.7305929 </td>
   <td style="text-align:right;"> 0.1225423 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.0176433 </td>
   <td style="text-align:right;"> 0.0151269 </td>
   <td style="text-align:right;"> -0.0555220 </td>
   <td style="text-align:right;"> 0.9811797 </td>
   <td style="text-align:right;"> 0.0005722 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> -0.0070348 </td>
   <td style="text-align:right;"> 0.0148454 </td>
   <td style="text-align:right;"> -0.0030272 </td>
   <td style="text-align:right;"> 0.3444626 </td>
   <td style="text-align:right;"> 0.9426560 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> -0.0056590 </td>
   <td style="text-align:right;"> 0.0080200 </td>
   <td style="text-align:right;"> -0.0555487 </td>
   <td style="text-align:right;"> 0.9915227 </td>
   <td style="text-align:right;"> 0.0001161 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -0.0317976 </td>
   <td style="text-align:right;"> 0.0115850 </td>
   <td style="text-align:right;"> 0.3736102 </td>
   <td style="text-align:right;"> 0.0024908 </td>
   <td style="text-align:right;"> 12.3325484 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.0946486 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 0.0456621 </td>
   <td style="text-align:right;"> 0.0220862 </td>
   <td style="text-align:right;"> 0.0231932 </td>
   <td style="text-align:right;"> 0.2439466 </td>
   <td style="text-align:right;"> 1.4511349 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> -0.0116262 </td>
   <td style="text-align:right;"> 0.0089351 </td>
   <td style="text-align:right;"> -0.0315909 </td>
   <td style="text-align:right;"> 0.5260204 </td>
   <td style="text-align:right;"> 0.4181532 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.0057842 </td>
   <td style="text-align:right;"> 0.0056288 </td>
   <td style="text-align:right;"> -0.0439537 </td>
   <td style="text-align:right;"> 0.6600238 </td>
   <td style="text-align:right;"> 0.2000402 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.0034655 </td>
   <td style="text-align:right;"> 0.0078797 </td>
   <td style="text-align:right;"> -0.0542577 </td>
   <td style="text-align:right;"> 0.8833205 </td>
   <td style="text-align:right;"> 0.0221588 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 0.0145279 </td>
   <td style="text-align:right;"> 0.0162366 </td>
   <td style="text-align:right;"> -0.0422719 </td>
   <td style="text-align:right;"> 0.6377302 </td>
   <td style="text-align:right;"> 0.2294078 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> 0.1358450 </td>
   <td style="text-align:right;"> 0.0504026 </td>
   <td style="text-align:right;"> -0.0233998 </td>
   <td style="text-align:right;"> 0.4457159 </td>
   <td style="text-align:right;"> 0.6112988 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 0.1155952 </td>
   <td style="text-align:right;"> 0.0461366 </td>
   <td style="text-align:right;"> -0.0357637 </td>
   <td style="text-align:right;"> 0.5648403 </td>
   <td style="text-align:right;"> 0.3439516 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.0074240 </td>
   <td style="text-align:right;"> 0.0167737 </td>
   <td style="text-align:right;"> -0.0489912 </td>
   <td style="text-align:right;"> 0.7410375 </td>
   <td style="text-align:right;"> 0.1126409 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.9915227 </td>
  </tr>
</tbody>
</table></div>

```r
#write.csv(lm_Tri.Prot, file="lm_Tri_prot.csv")
```

<br>  
<br>  

## Supplemental Figure 6. CEC distributions of stroke tPA status

__Distribution of cholesterol efflux capacity at 96 h post stroke with and without tPA treatment__


```r
cols <- c("No"="#202f52", "Yes" = "goldenrod1")

ggplot(PS96, aes(x = factor(tPa), y = PS96$CEC, fill = tPa)) +
    geom_violin(alpha=0.15, width=0.8) +
    geom_boxplot(width=0.4, alpha = 0.8) +
    geom_jitter(fill = "black", size = 3, shape = 21, 
                position = position_jitter(0.075), alpha=0.6) +
    scale_fill_manual(values=cols) +
    scale_x_discrete(limits=c("No", "Yes")) +
    labs(x="tPa", y="cholesterol efflux capacity (%)") + 
    theme_bw() +
    expand_limits(y=c(5,20))
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-50-1.png)<!-- -->

<br>  
__Distribution of cholesterol efflux capacity at 24 h post stroke with and without tPA treatment__

*(Figure not included in manuscript)*

```r
cols <- c("No"="#202f52", "Yes" = "goldenrod1")

ggplot(S24, aes(x = factor(tPa), y = S24$CEC, fill = tPa)) +
    geom_violin(alpha=0.15, width=0.8) +
    geom_boxplot(width=0.4, alpha = 0.8) +
    geom_jitter(fill = "black", size = 3, shape = 21, 
                position = position_jitter(0.075), alpha=0.6) +
    scale_fill_manual(values=cols) +
    scale_x_discrete(limits=c("No", "Yes")) +
    labs(x="tPa", y="cholesterol efflux capacity (%)") + 
    theme_bw() +
    expand_limits(y=c(5,20))
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-51-1.png)<!-- -->


## Supplemental Table 8. T-test statistics for plasma lipid levels by stroke tPA status

tPA significance: Just 24 h, just 96 h, and log2FC 24 - 96 h

__T-test for 24 h post stroke cholesterol efflux capacity with and without tPA__

```r
#I am too lazy to format the output:
t.test(S24$CEC ~ S24$tPa)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  S24$CEC by S24$tPa
## t = -0.50385, df = 32.952, p-value = 0.6177
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -2.521311  1.520420
## sample estimates:
##  mean in group No mean in group Yes 
##          10.31823          10.81868
```

__T-test for 96 h post stroke cholesterol efflux capacity with and without tPA__

```r
t.test(PS96$CEC ~ PS96$tPa)
```

```
## 
## 	Welch Two Sample t-test
## 
## data:  PS96$CEC by PS96$tPa
## t = -2.3625, df = 15.604, p-value = 0.03151
## alternative hypothesis: true difference in means is not equal to 0
## 95 percent confidence interval:
##  -5.9112956 -0.3138944
## sample estimates:
##  mean in group No mean in group Yes 
##          9.374977         12.487572
```

<br>
<br>


# 8. Additional info for reviewers

## APOA1 differences between groups

APOA1 is the major protein component of HDL particles. Therefore it was of interest to reviewers to know whether APOA1 showed a trend between our sample groups. Although not significant, APOA1 trends slightly higher in both stroke time points, inverse of the plasma HDL-C levels. 


```r
allViolinplot(logAll, logAll$APOA1) +
  ggtitle("APOA1") 
```

![](Stroke_manuscript_Ranalysis_files/figure-html/unnamed-chunk-54-1.png)<!-- -->



## Influence of outlier NIHSS on protein & recovery correlations

Since proteins are found to be associated with NIHSS recovery scores, we re-analyzed their relationship after removal of the individual with a very high NIHSS score (=37). The protein APOF still is significantly correlated after removal of the outlier, and after multiple testing correction.

__Statistics without outlier for model: NIHSS_3mo ~ NIHSS_baseline + Protein LFC__


```r
#dataset without the outlier individual
no_out <- data.frame(subset(lfc_df, NIHSS_3mo!=37))

out_df <- data.frame()

for (i in no_out[,25:ncol(no_out)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + i, data = no_out)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

olm_Prot.base <- out_df
rn <- data.frame(colnames(lfc_df[,25:ncol(no_out)]))
rownames(olm_Prot.base) <- rn[,1]
olm_Prot.base$p.adj <- p.adjust(olm_Prot.base$p, method = "BH")

kable(olm_Prot.base) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 0.8144147 </td>
   <td style="text-align:right;"> 4.163289 </td>
   <td style="text-align:right;"> 0.0809889 </td>
   <td style="text-align:right;"> 0.2174211 </td>
   <td style="text-align:right;"> 1.705009 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> 0.5539283 </td>
   <td style="text-align:right;"> 3.836684 </td>
   <td style="text-align:right;"> 0.2110991 </td>
   <td style="text-align:right;"> 0.0746816 </td>
   <td style="text-align:right;"> 3.140691 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.8680136 </td>
   <td style="text-align:right;"> 4.120343 </td>
   <td style="text-align:right;"> 0.0878042 </td>
   <td style="text-align:right;"> 0.2063824 </td>
   <td style="text-align:right;"> 1.770047 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.4904576 </td>
   <td style="text-align:right;"> 3.779554 </td>
   <td style="text-align:right;"> 0.2623291 </td>
   <td style="text-align:right;"> 0.0466761 </td>
   <td style="text-align:right;"> 3.844945 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 2.3133828 </td>
   <td style="text-align:right;"> 4.082918 </td>
   <td style="text-align:right;"> 0.1691456 </td>
   <td style="text-align:right;"> 0.1073324 </td>
   <td style="text-align:right;"> 2.628643 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -2.7107679 </td>
   <td style="text-align:right;"> 5.016232 </td>
   <td style="text-align:right;"> 0.2451942 </td>
   <td style="text-align:right;"> 0.1143096 </td>
   <td style="text-align:right;"> 2.786642 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 1.2727777 </td>
   <td style="text-align:right;"> 4.506484 </td>
   <td style="text-align:right;"> 0.0797137 </td>
   <td style="text-align:right;"> 0.2195417 </td>
   <td style="text-align:right;"> 1.692947 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 1.1188579 </td>
   <td style="text-align:right;"> 4.477085 </td>
   <td style="text-align:right;"> 0.0781307 </td>
   <td style="text-align:right;"> 0.2221989 </td>
   <td style="text-align:right;"> 1.678019 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 0.9239649 </td>
   <td style="text-align:right;"> 4.158255 </td>
   <td style="text-align:right;"> 0.0773938 </td>
   <td style="text-align:right;"> 0.2234451 </td>
   <td style="text-align:right;"> 1.671089 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 0.8423450 </td>
   <td style="text-align:right;"> 4.189520 </td>
   <td style="text-align:right;"> 0.0786209 </td>
   <td style="text-align:right;"> 0.2213731 </td>
   <td style="text-align:right;"> 1.682637 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 1.4946210 </td>
   <td style="text-align:right;"> 4.524750 </td>
   <td style="text-align:right;"> 0.0833480 </td>
   <td style="text-align:right;"> 0.2135443 </td>
   <td style="text-align:right;"> 1.727412 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 4.5943333 </td>
   <td style="text-align:right;"> 4.675745 </td>
   <td style="text-align:right;"> 0.1908274 </td>
   <td style="text-align:right;"> 0.0891959 </td>
   <td style="text-align:right;"> 2.886642 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 0.6394007 </td>
   <td style="text-align:right;"> 4.076051 </td>
   <td style="text-align:right;"> 0.1139508 </td>
   <td style="text-align:right;"> 0.1683686 </td>
   <td style="text-align:right;"> 2.028844 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 0.9421442 </td>
   <td style="text-align:right;"> 4.267120 </td>
   <td style="text-align:right;"> 0.0773584 </td>
   <td style="text-align:right;"> 0.2235052 </td>
   <td style="text-align:right;"> 1.670756 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -3.3333485 </td>
   <td style="text-align:right;"> 2.521519 </td>
   <td style="text-align:right;"> 0.6927753 </td>
   <td style="text-align:right;"> 0.0001014 </td>
   <td style="text-align:right;"> 19.039571 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.0038551 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 1.4677398 </td>
   <td style="text-align:right;"> 3.749809 </td>
   <td style="text-align:right;"> 0.2481330 </td>
   <td style="text-align:right;"> 0.0533388 </td>
   <td style="text-align:right;"> 3.640180 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 0.8796850 </td>
   <td style="text-align:right;"> 4.140807 </td>
   <td style="text-align:right;"> 0.0805077 </td>
   <td style="text-align:right;"> 0.2182193 </td>
   <td style="text-align:right;"> 1.700453 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 0.9913720 </td>
   <td style="text-align:right;"> 4.140878 </td>
   <td style="text-align:right;"> 0.0807411 </td>
   <td style="text-align:right;"> 0.2178318 </td>
   <td style="text-align:right;"> 1.702663 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 0.9233404 </td>
   <td style="text-align:right;"> 4.133208 </td>
   <td style="text-align:right;"> 0.0807599 </td>
   <td style="text-align:right;"> 0.2178005 </td>
   <td style="text-align:right;"> 1.702841 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.3717694 </td>
   <td style="text-align:right;"> 4.118828 </td>
   <td style="text-align:right;"> 0.1153067 </td>
   <td style="text-align:right;"> 0.1665734 </td>
   <td style="text-align:right;"> 2.042682 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 1.8781474 </td>
   <td style="text-align:right;"> 4.110813 </td>
   <td style="text-align:right;"> 0.1387306 </td>
   <td style="text-align:right;"> 0.1380478 </td>
   <td style="text-align:right;"> 2.288615 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -0.0105591 </td>
   <td style="text-align:right;"> 3.668521 </td>
   <td style="text-align:right;"> 0.2873763 </td>
   <td style="text-align:right;"> 0.0366504 </td>
   <td style="text-align:right;"> 4.226120 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.1593853 </td>
   <td style="text-align:right;"> 4.382672 </td>
   <td style="text-align:right;"> 0.0937358 </td>
   <td style="text-align:right;"> 0.1971696 </td>
   <td style="text-align:right;"> 1.827448 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 0.6078341 </td>
   <td style="text-align:right;"> 4.205270 </td>
   <td style="text-align:right;"> 0.0869776 </td>
   <td style="text-align:right;"> 0.2076951 </td>
   <td style="text-align:right;"> 1.762107 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 1.5312190 </td>
   <td style="text-align:right;"> 4.156768 </td>
   <td style="text-align:right;"> 0.1089198 </td>
   <td style="text-align:right;"> 0.1751757 </td>
   <td style="text-align:right;"> 1.977868 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 2.7590971 </td>
   <td style="text-align:right;"> 5.049204 </td>
   <td style="text-align:right;"> 0.1826453 </td>
   <td style="text-align:right;"> 0.1063426 </td>
   <td style="text-align:right;"> 2.675943 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 0.9306476 </td>
   <td style="text-align:right;"> 4.145314 </td>
   <td style="text-align:right;"> 0.0773667 </td>
   <td style="text-align:right;"> 0.2234910 </td>
   <td style="text-align:right;"> 1.670834 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 1.0736540 </td>
   <td style="text-align:right;"> 4.160416 </td>
   <td style="text-align:right;"> 0.0823612 </td>
   <td style="text-align:right;"> 0.2151586 </td>
   <td style="text-align:right;"> 1.718027 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 2.6143496 </td>
   <td style="text-align:right;"> 3.699334 </td>
   <td style="text-align:right;"> 0.2980753 </td>
   <td style="text-align:right;"> 0.0329678 </td>
   <td style="text-align:right;"> 4.397234 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> -0.2142649 </td>
   <td style="text-align:right;"> 4.024293 </td>
   <td style="text-align:right;"> 0.1725021 </td>
   <td style="text-align:right;"> 0.1043337 </td>
   <td style="text-align:right;"> 2.667698 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 1.7603363 </td>
   <td style="text-align:right;"> 4.326166 </td>
   <td style="text-align:right;"> 0.0995311 </td>
   <td style="text-align:right;"> 0.1885113 </td>
   <td style="text-align:right;"> 1.884260 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 1.3861492 </td>
   <td style="text-align:right;"> 4.228634 </td>
   <td style="text-align:right;"> 0.0908367 </td>
   <td style="text-align:right;"> 0.2016274 </td>
   <td style="text-align:right;"> 1.799299 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.7386214 </td>
   <td style="text-align:right;"> 4.011571 </td>
   <td style="text-align:right;"> 0.1361138 </td>
   <td style="text-align:right;"> 0.1410107 </td>
   <td style="text-align:right;"> 2.260479 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.5411645 </td>
   <td style="text-align:right;"> 4.084948 </td>
   <td style="text-align:right;"> 0.1156547 </td>
   <td style="text-align:right;"> 0.1661153 </td>
   <td style="text-align:right;"> 2.046240 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 1.5519329 </td>
   <td style="text-align:right;"> 3.810622 </td>
   <td style="text-align:right;"> 0.2262034 </td>
   <td style="text-align:right;"> 0.0652295 </td>
   <td style="text-align:right;"> 3.338634 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> -2.9595088 </td>
   <td style="text-align:right;"> 5.369377 </td>
   <td style="text-align:right;"> 0.1582030 </td>
   <td style="text-align:right;"> 0.1411118 </td>
   <td style="text-align:right;"> 2.315544 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 1.4434337 </td>
   <td style="text-align:right;"> 3.991002 </td>
   <td style="text-align:right;"> 0.1538817 </td>
   <td style="text-align:right;"> 0.1219198 </td>
   <td style="text-align:right;"> 2.454942 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.9657528 </td>
   <td style="text-align:right;"> 3.991743 </td>
   <td style="text-align:right;"> 0.1425617 </td>
   <td style="text-align:right;"> 0.1338063 </td>
   <td style="text-align:right;"> 2.330117 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.2235052 </td>
  </tr>
</tbody>
</table></div>

<br>  
<br>  

__Statistics without outlier for model: NIHSS_3mo ~ NIHSS_baseline + tPA + Protein LFC__


```r
out_df <- data.frame()

for (i in no_out[,25:ncol(no_out)]){
  fit <- lm(NIHSS_3mo ~ NIHSS_baseline + i + tPa, data = no_out)
  out <- data.frame(lm_out(fit))
  out_df <- data.frame(rbind(out_df, out))
}

olm_Prot.base <- out_df
rn <- data.frame(colnames(lfc_df[,25:ncol(no_out)]))
rownames(olm_Prot.base) <- rn[,1]
olm_Prot.base$p.adj <- p.adjust(olm_Prot.base$p, method = "BH")

kable(olm_Prot.base) %>%
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:300px; overflow-x: scroll; width:100%; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> i </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> se </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> r2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> f </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> df </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> p.adj </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> SERPINA1 </td>
   <td style="text-align:right;"> 1.1611163 </td>
   <td style="text-align:right;"> 4.464745 </td>
   <td style="text-align:right;"> 0.0168314 </td>
   <td style="text-align:right;"> 0.3875221 </td>
   <td style="text-align:right;"> 1.091304 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ALB </td>
   <td style="text-align:right;"> -0.0615344 </td>
   <td style="text-align:right;"> 4.198236 </td>
   <td style="text-align:right;"> 0.1626123 </td>
   <td style="text-align:right;"> 0.1586143 </td>
   <td style="text-align:right;"> 2.035680 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AMBP </td>
   <td style="text-align:right;"> 0.9809250 </td>
   <td style="text-align:right;"> 4.535403 </td>
   <td style="text-align:right;"> 0.0180550 </td>
   <td style="text-align:right;"> 0.3849668 </td>
   <td style="text-align:right;"> 1.098064 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ANTXR2 </td>
   <td style="text-align:right;"> -0.2929500 </td>
   <td style="text-align:right;"> 4.090548 </td>
   <td style="text-align:right;"> 0.2073063 </td>
   <td style="text-align:right;"> 0.1153354 </td>
   <td style="text-align:right;"> 2.394781 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APMAP </td>
   <td style="text-align:right;"> 2.7137104 </td>
   <td style="text-align:right;"> 4.401286 </td>
   <td style="text-align:right;"> 0.1122511 </td>
   <td style="text-align:right;"> 0.2212869 </td>
   <td style="text-align:right;"> 1.674372 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LPA </td>
   <td style="text-align:right;"> -3.0196796 </td>
   <td style="text-align:right;"> 4.326625 </td>
   <td style="text-align:right;"> 0.4391587 </td>
   <td style="text-align:right;"> 0.0558434 </td>
   <td style="text-align:right;"> 3.871130 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA1 </td>
   <td style="text-align:right;"> 1.5085328 </td>
   <td style="text-align:right;"> 4.781126 </td>
   <td style="text-align:right;"> 0.0128478 </td>
   <td style="text-align:right;"> 0.3959191 </td>
   <td style="text-align:right;"> 1.069413 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA2 </td>
   <td style="text-align:right;"> 1.3804323 </td>
   <td style="text-align:right;"> 4.763326 </td>
   <td style="text-align:right;"> 0.0115617 </td>
   <td style="text-align:right;"> 0.3986553 </td>
   <td style="text-align:right;"> 1.062384 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOA4 </td>
   <td style="text-align:right;"> 1.2317363 </td>
   <td style="text-align:right;"> 4.474806 </td>
   <td style="text-align:right;"> 0.0112491 </td>
   <td style="text-align:right;"> 0.3993221 </td>
   <td style="text-align:right;"> 1.060678 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOB </td>
   <td style="text-align:right;"> 1.1738418 </td>
   <td style="text-align:right;"> 4.471922 </td>
   <td style="text-align:right;"> 0.0145438 </td>
   <td style="text-align:right;"> 0.3923296 </td>
   <td style="text-align:right;"> 1.078712 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC3 </td>
   <td style="text-align:right;"> 1.7575707 </td>
   <td style="text-align:right;"> 4.822147 </td>
   <td style="text-align:right;"> 0.0168739 </td>
   <td style="text-align:right;"> 0.3874332 </td>
   <td style="text-align:right;"> 1.091539 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOC4 </td>
   <td style="text-align:right;"> 6.1070616 </td>
   <td style="text-align:right;"> 5.189360 </td>
   <td style="text-align:right;"> 0.1627832 </td>
   <td style="text-align:right;"> 0.1584279 </td>
   <td style="text-align:right;"> 2.036980 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOD </td>
   <td style="text-align:right;"> 1.0964344 </td>
   <td style="text-align:right;"> 4.367067 </td>
   <td style="text-align:right;"> 0.0567413 </td>
   <td style="text-align:right;"> 0.3099432 </td>
   <td style="text-align:right;"> 1.320824 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOE </td>
   <td style="text-align:right;"> 1.2612422 </td>
   <td style="text-align:right;"> 4.601103 </td>
   <td style="text-align:right;"> 0.0110701 </td>
   <td style="text-align:right;"> 0.3997044 </td>
   <td style="text-align:right;"> 1.059702 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOF </td>
   <td style="text-align:right;"> -4.4405590 </td>
   <td style="text-align:right;"> 2.673525 </td>
   <td style="text-align:right;"> 0.6995848 </td>
   <td style="text-align:right;"> 0.0002819 </td>
   <td style="text-align:right;"> 13.419873 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.0107127 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOL1 </td>
   <td style="text-align:right;"> 1.7004439 </td>
   <td style="text-align:right;"> 4.044498 </td>
   <td style="text-align:right;"> 0.1929546 </td>
   <td style="text-align:right;"> 0.1280765 </td>
   <td style="text-align:right;"> 2.275135 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> APOM </td>
   <td style="text-align:right;"> 1.1760332 </td>
   <td style="text-align:right;"> 4.475752 </td>
   <td style="text-align:right;"> 0.0138330 </td>
   <td style="text-align:right;"> 0.3938314 </td>
   <td style="text-align:right;"> 1.074811 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CAMP </td>
   <td style="text-align:right;"> 1.3601375 </td>
   <td style="text-align:right;"> 4.478451 </td>
   <td style="text-align:right;"> 0.0160713 </td>
   <td style="text-align:right;"> 0.3891152 </td>
   <td style="text-align:right;"> 1.087114 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLU </td>
   <td style="text-align:right;"> 1.1918348 </td>
   <td style="text-align:right;"> 4.472946 </td>
   <td style="text-align:right;"> 0.0133243 </td>
   <td style="text-align:right;"> 0.3949084 </td>
   <td style="text-align:right;"> 1.072023 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> C3 </td>
   <td style="text-align:right;"> 0.4089940 </td>
   <td style="text-align:right;"> 4.543170 </td>
   <td style="text-align:right;"> 0.0472962 </td>
   <td style="text-align:right;"> 0.3272323 </td>
   <td style="text-align:right;"> 1.264769 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PPBP </td>
   <td style="text-align:right;"> 1.8297837 </td>
   <td style="text-align:right;"> 4.371577 </td>
   <td style="text-align:right;"> 0.0726613 </td>
   <td style="text-align:right;"> 0.2822813 </td>
   <td style="text-align:right;"> 1.417892 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AHSG </td>
   <td style="text-align:right;"> -1.1046025 </td>
   <td style="text-align:right;"> 4.012901 </td>
   <td style="text-align:right;"> 0.2634698 </td>
   <td style="text-align:right;"> 0.0747339 </td>
   <td style="text-align:right;"> 2.907827 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FGA </td>
   <td style="text-align:right;"> 0.4954864 </td>
   <td style="text-align:right;"> 4.615849 </td>
   <td style="text-align:right;"> 0.0340870 </td>
   <td style="text-align:right;"> 0.3525216 </td>
   <td style="text-align:right;"> 1.188213 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> HPR </td>
   <td style="text-align:right;"> 0.7059121 </td>
   <td style="text-align:right;"> 4.857755 </td>
   <td style="text-align:right;"> 0.0169048 </td>
   <td style="text-align:right;"> 0.3873686 </td>
   <td style="text-align:right;"> 1.091709 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> IHH </td>
   <td style="text-align:right;"> 1.7870760 </td>
   <td style="text-align:right;"> 4.467451 </td>
   <td style="text-align:right;"> 0.0437753 </td>
   <td style="text-align:right;"> 0.3338462 </td>
   <td style="text-align:right;"> 1.244156 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ITIH4 </td>
   <td style="text-align:right;"> 3.6680923 </td>
   <td style="text-align:right;"> 5.721650 </td>
   <td style="text-align:right;"> 0.1255557 </td>
   <td style="text-align:right;"> 0.2163271 </td>
   <td style="text-align:right;"> 1.717917 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LCAT </td>
   <td style="text-align:right;"> 1.2520097 </td>
   <td style="text-align:right;"> 4.467496 </td>
   <td style="text-align:right;"> 0.0114693 </td>
   <td style="text-align:right;"> 0.3988522 </td>
   <td style="text-align:right;"> 1.061880 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SELL </td>
   <td style="text-align:right;"> 1.3791919 </td>
   <td style="text-align:right;"> 4.484490 </td>
   <td style="text-align:right;"> 0.0163140 </td>
   <td style="text-align:right;"> 0.3886060 </td>
   <td style="text-align:right;"> 1.088451 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PCYOX1 </td>
   <td style="text-align:right;"> 2.5459895 </td>
   <td style="text-align:right;"> 3.958682 </td>
   <td style="text-align:right;"> 0.2443703 </td>
   <td style="text-align:right;"> 0.0869995 </td>
   <td style="text-align:right;"> 2.724798 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GPLD1 </td>
   <td style="text-align:right;"> 0.2931842 </td>
   <td style="text-align:right;"> 4.265244 </td>
   <td style="text-align:right;"> 0.1250618 </td>
   <td style="text-align:right;"> 0.2038282 </td>
   <td style="text-align:right;"> 1.762335 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PF4 </td>
   <td style="text-align:right;"> 1.8104755 </td>
   <td style="text-align:right;"> 4.559238 </td>
   <td style="text-align:right;"> 0.0305584 </td>
   <td style="text-align:right;"> 0.3594975 </td>
   <td style="text-align:right;"> 1.168116 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLTP </td>
   <td style="text-align:right;"> 1.5535706 </td>
   <td style="text-align:right;"> 4.508899 </td>
   <td style="text-align:right;"> 0.0227981 </td>
   <td style="text-align:right;"> 0.3751673 </td>
   <td style="text-align:right;"> 1.124427 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON1 </td>
   <td style="text-align:right;"> 0.9806576 </td>
   <td style="text-align:right;"> 4.336050 </td>
   <td style="text-align:right;"> 0.0724076 </td>
   <td style="text-align:right;"> 0.2827077 </td>
   <td style="text-align:right;"> 1.416319 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PON3 </td>
   <td style="text-align:right;"> 0.9278551 </td>
   <td style="text-align:right;"> 4.385091 </td>
   <td style="text-align:right;"> 0.0553638 </td>
   <td style="text-align:right;"> 0.3124237 </td>
   <td style="text-align:right;"> 1.312579 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RBP4 </td>
   <td style="text-align:right;"> 1.1008838 </td>
   <td style="text-align:right;"> 4.074769 </td>
   <td style="text-align:right;"> 0.1776615 </td>
   <td style="text-align:right;"> 0.1428356 </td>
   <td style="text-align:right;"> 2.152236 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA1 </td>
   <td style="text-align:right;"> -3.8299588 </td>
   <td style="text-align:right;"> 6.725735 </td>
   <td style="text-align:right;"> 0.0861916 </td>
   <td style="text-align:right;"> 0.2839599 </td>
   <td style="text-align:right;"> 1.440166 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAA2 </td>
   <td style="text-align:right;"> 1.5482505 </td>
   <td style="text-align:right;"> 4.296471 </td>
   <td style="text-align:right;"> 0.0893798 </td>
   <td style="text-align:right;"> 0.2551932 </td>
   <td style="text-align:right;"> 1.523481 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTR </td>
   <td style="text-align:right;"> 0.8328318 </td>
   <td style="text-align:right;"> 4.336121 </td>
   <td style="text-align:right;"> 0.0773619 </td>
   <td style="text-align:right;"> 0.2744642 </td>
   <td style="text-align:right;"> 1.447193 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.3997044 </td>
  </tr>
</tbody>
</table></div>

<br>  



<br>
<br>  
<br>  

... I hope you found this document helpful! Sorry for any terrible syntax or poor commenting - feel free to reach out with comments or questions!

-Deanna



```r
sessionInfo()
```

```
## R version 3.6.3 (2020-02-29)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18363)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] xlsx_0.6.3         statmod_1.4.34     kableExtra_1.1.0   knitr_1.28        
##  [5] stargazer_5.2.2    scales_1.1.1       limma_3.42.2       magrittr_1.5      
##  [9] Gmisc_1.10.0       htmlTable_1.13.3   Rcpp_1.0.4         plyr_1.8.6        
## [13] tidyr_1.0.3        purrr_0.3.3        mclust_5.4.6       corrplot_0.84     
## [17] Hmisc_4.4-0        Formula_1.2-3      survival_3.1-8     lattice_0.20-38   
## [21] psych_1.9.12.31    pheatmap_1.0.12    gridExtra_2.3      dplyr_0.8.5       
## [25] RColorBrewer_1.1-2 ggplot2_3.3.0     
## 
## loaded via a namespace (and not attached):
##  [1] httr_1.4.1          viridisLite_0.3.0   splines_3.6.3      
##  [4] assertthat_0.2.1    highr_0.8           latticeExtra_0.6-29
##  [7] xlsxjars_0.6.1      yaml_2.2.1          pillar_1.4.4       
## [10] backports_1.1.5     glue_1.3.2          digest_0.6.25      
## [13] checkmate_2.0.0     rvest_0.3.5         colorspace_1.4-1   
## [16] htmltools_0.4.0     Matrix_1.2-18       XML_3.99-0.3       
## [19] pkgconfig_2.0.3     webshot_0.5.2       jpeg_0.1-8.1       
## [22] tibble_3.0.1        mgcv_1.8-31         farver_2.0.3       
## [25] generics_0.0.2      ellipsis_0.3.0      withr_2.2.0        
## [28] nnet_7.3-12         mnormt_1.5-7        crayon_1.3.4       
## [31] evaluate_0.14       nlme_3.1-144        xml2_1.3.2         
## [34] foreign_0.8-75      tools_3.6.3         data.table_1.12.8  
## [37] hms_0.5.3           lifecycle_0.2.0     stringr_1.4.0      
## [40] munsell_0.5.0       cluster_2.1.0       compiler_3.6.3     
## [43] rlang_0.4.6         rstudioapi_0.11     htmlwidgets_1.5.1  
## [46] labeling_0.3        base64enc_0.1-3     rmarkdown_2.1      
## [49] gtable_0.3.0        abind_1.4-5         forestplot_1.9     
## [52] R6_2.4.1            lubridate_1.7.8     readr_1.3.1        
## [55] rJava_0.9-12        stringi_1.4.6       parallel_3.6.3     
## [58] vctrs_0.2.4         rpart_4.1-15        acepack_1.4.1      
## [61] png_0.1-7           tidyselect_1.0.0    xfun_0.13
```


