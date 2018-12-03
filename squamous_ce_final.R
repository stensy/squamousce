#####
'''
squamous cell comparative effectiveness for manuscript

    specifically, this is using the NCDB 2014 dataset

multiple parts to analysis.

PART 1: ALL TREATMENTS
I. importing data, cleaning
    includes coding of histology
II. excluding cases
    note that some exclusions (eg no treatment) come after
    coding for type of treatment received
III. multinomial propensity scoring via twang
IV. cox proportional hazards modeling via survey design

PART 2: RC VS RC + ADJUVANT CHEMO
?exclude patients who died within 30 days perioperatively?
I. import and clean data
II. exclude cases
III. propensity scoring
IV. cox proportional hazards modeling

'''
##### importing the info
rm(list=ls())
setwd("~/Desktop/ncdb_2014/bladder")

library(readr)
library(twang)
library(MatchIt)
library(plyr)
library(dplyr)

ncdb <- read_csv("~/Desktop/ncdb_2014/bladder/ncdb_bladder_delimited.txt")

##### variable cleaning #####
ncdb$TNM_CLIN_M <- factor(ncdb$TNM_CLIN_M)
ncdb$TNM_CLIN_T <- factor(ncdb$TNM_CLIN_T)
ncdb$TNM_CLIN_N <- factor(ncdb$TNM_CLIN_N)
ncdb$SEX <- factor(ncdb$SEX)

##### coding histology #####
ncdb$HISTOLOGY <- mapvalues(ncdb$HISTOLOGY, 
                         from = c('8000','8001','8002','8003','8004','8005','8010','8011','8012','8013','8014','8020',
                                    '8021','8022','8030','8031','8032','8033','8035','8041','8042','8044','8045','8046',
                                    '8050','8051','8052','8070','8071','8072','8073','8074','8076','8082','8083','8084',
                                    '8103','8120','8122','8123','8124','8130','8131','8140','8141','8142','8144','8145',
                                    '8160','8170','8210','8230','8240','8243','8244','8245','8246','8255','8260','8261',
                                    '8262','8263','8270','8310','8312','8315','8318','8320','8323','8340','8341','8380',
                                    '8440','8460','8470','8480','8481','8490','8500','8503','8507','8523','8542','8550',
                                    '8560','8572','8574','8575','8680','8693','8700','8720','8730','8772','8800','8801',
                                    '8802','8803','8804','8805','8806','8810','8811','8815','8825','8830','8840','8850',
                                    '8854','8857','8890','8891','8895','8896','8900','8901','8910','8912','8920','8935',
                                    '8936','8940','8950','8951','8963','8980','8981','8982','8990','9061','9071','9100',
                                    '9110','9120','9130','9150','9170','9180','9220','9231','9260','9364','9473','9540',
                                    '9580'),
                         to = c('Neoplasm, malignant','Tumor cells, malignant','Malignant tumor, small cell type',
                                    'Malignant tumor, giant cell type','Malignant tumor, spindle cell type',
                                    'Malignant tumor, clear cell type','Carcinoma, NOS','Epithelioma, malignant',
                                    'Large cell carcinoma, NOS','Large cell neuroendocrine carcinoma',
                                    'Large cell carcinoma with rhabdoid phenotype','Carcinoma, undifferentiated type, NOS',
                                    'Carcinoma, anaplastic type, NOS','Pleomorphic carcinoma',
                                    'Giant cell and spindle cell carcinoma','Giant cell carcinoma','Spindle cell carcinoma',
                                    'Pseudosarcomatous carcinoma','Carcinoma with osteoclast-like giant cells',
                                    'Small cell carcinoma, NOS','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','Papillary carcinoma, NOS','Verrucous carcinoma, NOS',
                                    'Papillary squamous cell carcinoma','Squamous cell carcinoma, NOS',
                                    'Sq. cell carcinoma, keratinizing, NOS','Sq. cell carcinoma, lg. cell, non-ker.',
                                    'Sq. cell carcinoma, sm. cell, non-ker.','Sq. cell carcinoma, spindle cell',
                                    'Sq. cell carcinoma, micro-invasive','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','Transitional cell carcinoma, NOS','Trans. cell carcinoma, spindle cell',
                                    'Basaloid carcinoma','Cloacogenic carcinoma','Papillary trans. cell carcinoma',
                                    'Transitional cell carcinoma, micropapillary','Adenocarcinoma, NOS',
                                    'Scirrhous adenocarcinoma','unknown_path','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','Solid carcinoma, NOS','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','unknown_path','Adenocarcinoma with mixed subtypes',
                                    'Papillary adenocarcinoma, NOS','Adenocarcinoma in villous adenoma','unknown_path',
                                    'unknown_path','unknown_path','Clear cell adenocarcinoma, NOS','unknown_path',
                                    'unknown_path','unknown_path','Granular cell carcinoma','Mixed cell adenocarcinoma',
                                    'unknown_path','unknown_path','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','Mucinous adenocarcinoma','Mucin-producing adenocarcinoma',
                                    'Signet ring cell carcinoma','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','Paraganglioma, malignant','unknown_path',
                                    'unknown_path','unknown_path','unknown_path','unknown_path','Sarcoma, NOS',
                                    'Spindle cell sarcoma','Giant cell sarcoma','Small cell sarcoma',
                                    'Epithelioid sarcoma','Undifferentiated sarcoma','Desmoplastic small round cell tumor',
                                    'Fibrosarcoma, NOS','Fibromyxosarcoma','Solitary fibrous tumor, malignant',
                                    'unknown_path','Fibrous histiocytoma, malignant','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','Leiomyosarcoma, NOS','Epithelioid leiomyosarcoma',
                                    'Myosarcoma','Myxoid leiomyosarcoma','Rhabdomyosarcoma, NOS',
                                    'Pleomorphic rhabdomyosarcoma, adult type','Embryonal rhabdomyosarcoma',
                                    'Spindle cell rhabdomyosarcoma','Alveolar rhabdomyosarcoma','unknown_path',
                                    'unknown_path','unknown_path','Mullerian mixed tumor','Mesodermal mixed tumor',
                                    'unknown_path','Carcinosarcoma, NOS','Carcinosarcoma, embryonal type',
                                    'Malignant myoepithelioma','Mesenchymoma, malignant','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','unknown_path','unknown_path','unknown_path',
                                    'unknown_path','unknown_path','unknown_path','unknown_path'))

ncdb$HISTOLOGY <- factor(ncdb$HISTOLOGY)
histlevs <- levels(ncdb$HISTOLOGY)

ncdb$hist_rc <- mapvalues(ncdb$HISTOLOGY, from=c('Large cell carcinoma, NOS','Mesenchymoma, malignant','Giant cell sarcoma','Papillary trans. cell carcinoma','Trans. cell carcinoma, spindle cell','Papillary carcinoma, NOS','Malignant tumor, small cell type','Giant cell carcinoma','Pseudosarcomatous carcinoma','Epithelioid sarcoma','Sarcoma, NOS','Adenocarcinoma, NOS','Fibromyxosarcoma','Epithelioid leiomyosarcoma','Large cell neuroendocrine carcinoma','Granular cell carcinoma','Desmoplastic small round cell tumor','Mucinous adenocarcinoma','Mixed cell adenocarcinoma','Spindle cell carcinoma','Alveolar rhabdomyosarcoma','Malignant tumor, spindle cell type','Spindle cell sarcoma','Mucin-producing adenocarcinoma','Mullerian mixed tumor','Solitary fibrous tumor, malignant','Tumor cells, malignant','Malignant tumor, giant cell type','Signet ring cell carcinoma','Adenocarcinoma with mixed subtypes','Pleomorphic rhabdomyosarcoma, adult type','Small cell sarcoma','Leiomyosarcoma, NOS','Myxoid leiomyosarcoma','Adenocarcinoma in villous adenoma','Transitional cell carcinoma, NOS','Papillary squamous cell carcinoma','Malignant tumor, clear cell type','Carcinoma, NOS','Pleomorphic carcinoma','Mesodermal mixed tumor','Fibrosarcoma, NOS','Papillary adenocarcinoma, NOS','Clear cell adenocarcinoma, NOS','Paraganglioma, malignant','Carcinoma with osteoclast-like giant cells','Small cell carcinoma, NOS','Cloacogenic carcinoma','Rhabdomyosarcoma, NOS','Neoplasm, malignant','Squamous cell carcinoma, NOS','Scirrhous adenocarcinoma','Large cell carcinoma with rhabdoid phenotype','Sq. cell carcinoma, spindle cell','Malignant myoepithelioma','Sq. cell carcinoma, sm. cell, non-ker.','Fibrous histiocytoma, malignant','Sq. cell carcinoma, keratinizing, NOS','Carcinoma, anaplastic type, NOS','Myosarcoma','Sq. cell carcinoma, micro-invasive','Verrucous carcinoma, NOS','Giant cell and spindle cell carcinoma','Transitional cell carcinoma, micropapillary','Sq. cell carcinoma, lg. cell, non-ker.','Carcinoma, undifferentiated type, NOS','Solid carcinoma, NOS','Carcinosarcoma, NOS','Carcinosarcoma, embryonal type','Basaloid carcinoma','Spindle cell rhabdomyosarcoma','unknown_path','Undifferentiated sarcoma','Embryonal rhabdomyosarcoma','Epithelioma, malignant'),
                       to=c('large_cell','other','giant_cell','transitional_cell','spindle_cell','transitional_cell','small_cell','giant_cell','other','other','nonspecific','adenocarcinoma','other','myosarcoma','large_cell','other','other','adenocarcinoma','other','spindle_cell','myosarcoma','spindle_cell','spindle_cell','other','other','other','nonspecific','giant_cell','signet_ring','adenocarcinoma','myosarcoma','small_cell','myosarcoma','myosarcoma','adenocarcinoma','transitional_cell','squamous_cell','clear_cell','nonspecific','other','other','other','adenocarcinoma','clear_cell','other','other','small_cell','other','myosarcoma','nonspecific','squamous_cell','adenocarcinoma','large_cell','squamous_cell','other','squamous_cell','other','squamous_cell','other','myosarcoma','squamous_cell','other','other','micropapillary','squamous_cell','nonspecific','nonspecific','other','other','other','myosarcoma','nonspecific','other','myosarcoma','other'))

ncdb$hist_rc <- factor(ncdb$hist_rc)

### recode race
ncdb$RACE <- factor(ncdb$RACE)
ncdb$race_rc <- ncdb$RACE
ncdb$race_rc <- recode(ncdb$race_rc, "1: White"="white", 
                       "10: Vietnamese"="other", 
                       "15: Asian Indian or Pakistani"="other",
                       "16: Asian Indian"="other",
                       "2: Black"="black",
                       "20: Micronesian, NOS"="other",
                       "3: American Indian, Aleutian, or Eskimo"="other",
                       "4: Chinese"="other",
                       "6: Filipino"="other",
                       "7: Hawaiian"="other",
                       "8: Korean"="other",
                       "11: Laotian"="other",
                       "12: Hmong"="other",
                       "13: Kampuchean"="other",
                       "14: Thai"="other",
                       "17: Pakistani"="other",
                       "21: Chamorran"="other",
                       "22: Guamanian, NOS"="other",
                       "25: Polynesian, NOS"="other",
                       "27: Samoan"="other",
                       "28: Tongan"="other",
                       "30: Melanesian, NOS"="other",
                       "31: Fiji Islander"="other",
                       "32: New Guinean"="other",
                       "5: Japanese"="other",
                       "97: Pacific Islander, NOS"="other",
                       "96: Other Asian, including Asian, NOS and Oriental, NOS"="other",
                       "98: Other"="other",
                       "99: Unknown"="other")

ncdb$race_rc <- factor(ncdb$race_rc)

### recode insurance status
ncdb$INSURANCE_STATUS <- factor(ncdb$INSURANCE_STATUS)
ncdb$insurance_rc <- ncdb$INSURANCE_STATUS
ncdb$insurance_rc <- recode(ncdb$insurance_rc, "0: Not Insured" = "uninsured",
                            "1: Private Insurance"="private", "2: Medicaid" = "medicaid",
                            "3: Medicare" = "medicare", "4: Other Government"="other",
                            "9: Insurance Status Unknown"="other")

ncdb$insurance_rc <- factor(ncdb$insurance_rc)

### factor income
ncdb$MED_INC_QUAR_12 <- factor(ncdb$MED_INC_QUAR_12)

### factor education
ncdb$NO_HSD_QUAR_12 <- factor(ncdb$NO_HSD_QUAR_12)

### recode urban rural
ncdb$UR_CD_13 <- factor(ncdb$UR_CD_13)
ncdb$rural_rc <- ncdb$UR_CD_13
ncdb$rural_rc <- recode(ncdb$rural_rc, 
                "1: Counties in metro areas of 1 million population or more" = "metro >1 million",
                "2: Counties in metro areas of 250,000 to 1 million population" = "metro 250k-1million",
                "3: Counties in metro areas of fewer than 250,000 population" = "metro < 250k",
                "4: Urban population of 20,000 or more adjacent to a metro area" = "pop<20k",
                "5: Urban population of 20,000 or more not adjacent to a metro area" = "pop<20k",
                "6: Urban population of 2,500 to 19,999, adjacent to a metro area" = "pop<20k",
                "7: Urban population of 2,500 to 19,999, not adjacent to a metro area" = "pop<20k",
                "8: Completely rural or less than 2,500 urban population, adjacent to a metro area" = "pop<20k",
                "9: Completely rural or less than 2,500 urban population, not adjacent to a metro area" = "pop<20k")

ncdb$rural_rc <- factor(ncdb$rural_rc)

### factor charlson
ncdb$CDCC_TOTAL <- factor(ncdb$CDCC_TOTAL)

### factor year
ncdb$YEAR_OF_DIAGNOSIS <- factor(ncdb$YEAR_OF_DIAGNOSIS)

### factor facility type
ncdb$FACILITY_TYPE_CD <- factor(ncdb$FACILITY_TYPE_CD)

### factor facility location
ncdb$FACILITY_LOCATION_CD <- factor(ncdb$FACILITY_LOCATION_CD)

##### cutting to just squamous cell; descriptives #####

## first, include only squamous cell

ncdb <- subset(ncdb, ncdb$hist_rc == "squamous_cell")

ncdb$hist_rc <- factor(ncdb$hist_rc)

##### descriptive statistics for all squamous cell carcinoma in NCDB #####

## clinical staging
table(ncdb$TNM_CLIN_T)
prop.table(table(ncdb$TNM_CLIN_T))
table(ncdb$TNM_CLIN_N)
prop.table(table(ncdb$TNM_CLIN_N))
table(ncdb$TNM_CLIN_M)
prop.table(table(ncdb$TNM_CLIN_M))

## pathologic staging
table(ncdb$TNM_PATH_T)
prop.table(table(ncdb$TNM_PATH_T))
table(ncdb$TNM_PATH_N)
prop.table(table(ncdb$TNM_PATH_N))
table(ncdb$TNM_PATH_M)
prop.table(table(ncdb$TNM_PATH_M))

##### excluding cases #####

## missing data
# missing followup time
ncdb <- subset(ncdb, ncdb$DX_LASTCONTACT_DEATH_MONTHS != "NA")

# exclude for metastatic disease

ncdb <- subset(ncdb, ncdb$TNM_CLIN_M != "1" & ncdb$TNM_CLIN_M != "88" & ncdb$TNM_CLIN_M != "NA")

ncdb <- subset(ncdb, ncdb$TNM_CLIN_T == "2" | ncdb$TNM_CLIN_T == "2A" | 
                   ncdb$TNM_CLIN_T == "2B" | ncdb$TNM_CLIN_T == "3" |
                   ncdb$TNM_CLIN_T == "3A" | ncdb$TNM_CLIN_T == "3B")

ncdb <- subset(ncdb, ncdb$TNM_CLIN_N == "0")
ncdb$TNM_CLIN_T <- factor(ncdb$TNM_CLIN_T)


##### collapse the T stages into just 2 v 3
ncdb$TNM_CLIN_T <- factor(ncdb$TNM_CLIN_T)
ncdb$TNM_CLIN_T[ncdb$TNM_CLIN_T == "2A"] <- "2"
ncdb$TNM_CLIN_T[ncdb$TNM_CLIN_T == "2B"] <- "2"
ncdb$TNM_CLIN_T[ncdb$TNM_CLIN_T == "3A"] <- "3"
ncdb$TNM_CLIN_T[ncdb$TNM_CLIN_T == "3B"] <- "3"

##### code for treatment types #####
## code for radical cystectomy
ncdb$rad_cyst <- 0
ncdb$rad_cyst[is.na(ncdb$rad_cyst)] <- 0
ncdb$rad_cyst[ncdb$RX_SUMM_SURG_PRIM_SITE > 49 & ncdb$RX_SUMM_SURG_PRIM_SITE < 90] <- 1
ncdb$rad_cyst <- factor(ncdb$rad_cyst)

## code for neoadjuvant chemo
ncdb$nac <- 0
ncdb$chemsurgtime <- ncdb$DX_DEFSURG_STARTED_DAYS - ncdb$DX_CHEMO_STARTED_DAYS
ncdb$nac[ncdb$rad_cyst == 1 & (ncdb$RX_SUMM_CHEMO == "1: Chemotherapy administered, type and number of agents not documented" | ncdb$RX_SUMM_CHEMO == "2: Single-agent chemotherapy" | ncdb$RX_SUMM_CHEMO == "3: Multiagent chemotherapy") & (ncdb$chemsurgtime > -181) & (ncdb$chemsurgtime < 0)] <- 1
ncdb$nac <- factor(ncdb$nac)

ncdb$RX_SUMM_SYSTEMIC_SUR_SEQ <- factor(ncdb$RX_SUMM_SYSTEMIC_SUR_SEQ)
ncdb$RX_SUMM_RADIATION <- factor(ncdb$RX_SUMM_RADIATION)

## code for adjuvant chemo
ncdb$ac <- 0
ncdb$ac[ncdb$rad_cyst == 1 & (ncdb$RX_SUMM_CHEMO == "1: Chemotherapy administered, type and number of agents not documented" | ncdb$RX_SUMM_CHEMO == "2: Single-agent chemotherapy" | ncdb$RX_SUMM_CHEMO == "3: Multiagent chemotherapy") & (ncdb$chemsurgtime > 0) & (ncdb$chemsurgtime < 180)] <- 1
ncdb$ac <- factor(ncdb$ac)

## code for chemo alone
ncdb$only_chemo <- 0
ncdb$only_chemo[(ncdb$RX_SUMM_CHEMO == "1: Chemotherapy administered, type and number of agents not documented" | ncdb$RX_SUMM_CHEMO == "2: Single-agent chemotherapy" | ncdb$RX_SUMM_CHEMO == "3: Multiagent chemotherapy") & ncdb$rad_cyst == 0 & ncdb$RX_SUMM_RADIATION == "0: None"] <- 1
ncdb$only_chemo <- factor(ncdb$only_chemo)

## code for radiation alone
ncdb$only_rads <- 0
ncdb$only_rads[ncdb$RX_SUMM_RADIATION == "1: Beam radiation" & ncdb$RX_SUMM_CHEMO == "0: None" & ncdb$rad_cyst == 0] <- 1
ncdb$only_rads <- factor(ncdb$only_rads)

## code for cystectomy alone
ncdb$only_cyst <- 0
ncdb$only_cyst[ncdb$rad_cyst == 1 & ncdb$RX_SUMM_RADIATION == "0: None" & ncdb$RX_SUMM_CHEMO == "0: None"] <- 1
ncdb$only_cyst <- factor(ncdb$only_cyst)

## code for chemorads
ncdb$chemorads <- 0
ncdb$chemorads[ncdb$rad_cyst == 0 & ncdb$RX_SUMM_RADIATION == "1: Beam radiation" & (ncdb$RX_SUMM_CHEMO == "1: Chemotherapy administered, type and number of agents not documented" | ncdb$RX_SUMM_CHEMO == "2: Single-agent chemotherapy" | ncdb$RX_SUMM_CHEMO == "3: Multiagent chemotherapy")] <- 1
ncdb$chemorads <- factor(ncdb$chemorads)

## code for radiation and surgery
ncdb$surg_rads <- 0
ncdb$surg_rads[ncdb$rad_cyst == 1 & ncdb$RX_SUMM_RADIATION == "1: Beam radiation"] <- 1
ncdb$surg_rads <- factor(ncdb$surg_rads)

## code for chemo, radiation, cystectomy
ncdb$trimodal <- 0
ncdb$trimodal[ncdb$rad_cyst == 1 & ncdb$RX_SUMM_RADIATION == "1: Beam radiation" & (ncdb$RX_SUMM_CHEMO == "1: Chemotherapy administered, type and number of agents not documented" | ncdb$RX_SUMM_CHEMO == "2: Single-agent chemotherapy" | ncdb$RX_SUMM_CHEMO == "3: Multiagent chemotherapy")] <- 1
ncdb$trimodal <- factor(ncdb$trimodal)

### now combine all those treatment options into one variable
ncdb$treatment <- 0
ncdb$treatment[ncdb$only_cyst == 1] <- 1
ncdb$treatment[ncdb$nac == 1] <- 2
ncdb$treatment[ncdb$ac == 1] <- 3
ncdb$treatment[ncdb$only_chemo == 1] <- 4
ncdb$treatment[ncdb$only_rads == 1] <- 5
ncdb$treatment[ncdb$chemorads == 1] <- 6
ncdb$treatment[ncdb$surg_rads == 1] <- 7
ncdb$treatment[ncdb$trimodal == 1] <- 8
ncdb$treatment <- factor(ncdb$treatment, levels = c(0:8), labels=c("None", "Cystectomy Alone", "Neoadjuvant Chemo", "Adjuvant Chemo", "Chemo Alone", "Radiation Alone", "Chemorads", "Cystectomy and radiation", "Trimodal"))

ncdb <- subset(ncdb, ncdb$treatment != "Cystectomy and radiation" & ncdb$treatment != "Trimodal")

ncdb$treatment <- factor(ncdb$treatment)

## drop those not receiving treatment
ncdb <- subset(ncdb, ncdb$treatment != "None")
ncdb$treatment <- factor(ncdb$treatment)

table(ncdb$treatment)
prop.table(table(ncdb$treatment))

##### multinomial propensity score via twang #####

ncdb$TNM_CLIN_T <- factor(ncdb$TNM_CLIN_T)
ncdb$YEAR_OF_DIAGNOSIS <- factor(ncdb$YEAR_OF_DIAGNOSIS)
ncdb <- as.data.frame(ncdb)

# make multinomial ps model
mnps.scc <- mnps(treatment ~ AGE + TNM_CLIN_T +
                     SEX + CDCC_TOTAL + race_rc + insurance_rc + 
                     MED_INC_QUAR_12 + NO_HSD_QUAR_12 + rural_rc + 
                     YEAR_OF_DIAGNOSIS + FACILITY_TYPE_CD + 
                     FACILITY_LOCATION_CD,
                 data=ncdb,
                 n.trees = 3000,
                 interaction.depth = 2,
                 shrinkage = 0.01,
                 perm.test.iters = 0,
                 stop.method = c("es.mean", "ks.mean"),
                 estimand = "ATE",
                 verbose = FALSE)

# create balance plots
plot(mnps.scc, plots=4)
plot(mnps.scc, plots=3)
plot(mnps.scc, plots=2, subset="ks.mean")

# create balance table
bal.table(mnps.scc, stop.method = "ks.mean")
bal.table(mnps.scc, stop.method = "ks.mean", collapse.to='covariate')

# give estimates
summary(mnps.scc)

# create survival variable
ncdb$surv <- Surv(ncdb$DX_LASTCONTACT_DEATH_MONTHS, ncdb$PUF_VITAL_STATUS=="0: Dead")

# extract weights from multinomial ps model
ncdb$w <- get.weights(mnps.scc, stop.method= "ks.mean")


###### kaplan meier for overall survival strat by treatment #####
survfit(ncdb$surv ~ ncdb$treatment)

##### cox proportional hazards model via survey design #####
# use design package
design.mnps <- svydesign(ids=~1, weights=~w, data=ncdb)

# use cox model from design package
coxmodel <- svycoxph(surv ~ treatment + AGE + TNM_CLIN_T + SEX + CDCC_TOTAL +
                    race_rc + insurance_rc + MED_INC_QUAR_12 + NO_HSD_QUAR_12 +
                        rural_rc + YEAR_OF_DIAGNOSIS + FACILITY_TYPE_CD +
                        FACILITY_LOCATION_CD, 
         design=design.mnps)

summary(coxmodel)

##### subset for conditional RC vs RC + AC analysis ####
ncdb_ac <- subset(ncdb, ncdb$treatment == "Cystectomy Alone" | ncdb$treatment == "Adjuvant Chemo")
ncdb_ac$treatment <- factor(ncdb_ac$treatment)

### recode TNM_PATH_T, TNM_PATH_N, TNM_PATH_M
ncdb_ac$TNM_PATH_T <- factor(ncdb_ac$TNM_PATH_T)
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="IS"] <- "IS"
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="X"] <- NA
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="2A"] <- "2"
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="2B"] <- "2"
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="3A"] <- "3"
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="3B"] <- "3"
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="4A"] <- "4"
ncdb_ac$TNM_PATH_T[ncdb_ac$TNM_PATH_T=="4B"] <- "4"
ncdb_ac$TNM_PATH_T <- factor(ncdb_ac$TNM_PATH_T)
ncdb_ac$TNM_PATH_N <- factor(ncdb_ac$TNM_PATH_N)
ncdb_ac$TNM_PATH_N[ncdb_ac$TNM_PATH_N == "X"] <- NA
ncdb_ac$TNM_PATH_N <- factor(ncdb_ac$TNM_PATH_N)

## path stage unavailable for some of these folks; drop them
ncdb_ac <- subset(ncdb_ac, is.na(ncdb_ac$TNM_PATH_T) == FALSE)
ncdb_ac <- subset(ncdb_ac, is.na(ncdb_ac$TNM_PATH_N) == FALSE)

## do the propensity score MATCH this time
# first only keep the necessary variables
keepvars <- c("AGE", "SEX", "CDCC_TOTAL", "TNM_PATH_T",
              "TNM_PATH_N", "PUF_30_DAY_MORT_CD", "treatment",
              "PUF_VITAL_STATUS", "DX_LASTCONTACT_DEATH_MONTHS")
ncdb_ac <- ncdb_ac[,keepvars]
ncdb_ac$adjuvant <- 0
ncdb_ac$adjuvant[ncdb_ac$treatment == "Adjuvant Chemo"] <- 1
scc_match <- matchit(adjuvant ~ AGE + TNM_PATH_T + TNM_PATH_N +
                         SEX + CDCC_TOTAL, data = ncdb_ac, 
                     method = "nearest", ratio = 3)
summary(scc_match)
plot(scc_match, type = "hist")
plot(scc_match, type = "jitter")

ac_matched <- match.data(scc_match)

ac_matched$surv <- Surv(ac_matched$DX_LASTCONTACT_DEATH_MONTHS, ac_matched$PUF_VITAL_STATUS == "0: Dead")

ac_matched$adjuvant <- factor(ac_matched$adjuvant)
ac_matched$TNM_PATH_T <- relevel(ac_matched$TNM_PATH_T, ref = "2")
ac_matched$TNM_PATH_N <- factor(ac_matched$TNM_PATH_N)
ac_matched$CDCC_TOTAL <- factor(ac_matched$CDCC_TOTAL)

ac_cox <- coxph(surv ~ adjuvant + AGE + TNM_PATH_T + TNM_PATH_N +
                    SEX + CDCC_TOTAL, data = ac_matched)
summary(ac_cox)

##### secondary analysis lumping all RC treatments into one #####

##### multinomial propensity score via twang #####

ncdb$TNM_CLIN_T <- factor(ncdb$TNM_CLIN_T)
ncdb <- as.data.frame(ncdb)

ncdb$treatment2 <- 0
ncdb$treatment2[ncdb$treatment == "Cystectomy Alone"] <- "cystectomy"
ncdb$treatment2[ncdb$treatment == "Neoadjuvant Chemo"] <- "cystectomy"
ncdb$treatment2[ncdb$treatment == "Adjuvant Chemo"] <- "cystectomy"
ncdb$treatment2[ncdb$treatment == "Chemo Alone"] <- "chemo alone"
ncdb$treatment2[ncdb$treatment == "Radiation Alone"] <- "rads alone"
ncdb$treatment2[ncdb$treatment == "Chemorads"] <- "chemorads"

ncdb$treatment2 <- factor(ncdb$treatment2)


mnps2.scc <- mnps(treatment2 ~ AGE + TNM_CLIN_T +
                     SEX + CDCC_TOTAL,
                 data=ncdb,
                 n.trees = 3000,
                 interaction.depth = 2,
                 shrinkage = 0.01,
                 perm.test.iters = 0,
                 stop.method = c("es.mean", "ks.mean"),
                 estimand = "ATE",
                 verbose = FALSE)

#plot(mnps2.scc, plots=4)
#bal.table(mnps.scc, stop.method = "ks.mean")
#bal.table(mnps.scc, collapse.to='covariate')

summary(mnps2.scc)

ncdb$w2 <- get.weights(mnps2.scc, stop.method= "ks.mean")


##### cox proportional hazards model via survey design #####
ncdb$treatment2 <- relevel(ncdb$treatment2, ref = "cystectomy")
design2.mnps <- svydesign(ids=~1, weights=~w2, data=ncdb)

coxmodel2 <- svycoxph(surv ~ treatment2 + AGE + TNM_CLIN_T + SEX + CDCC_TOTAL, 
                     design=design2.mnps)
summary(coxmodel2)


##### exploratory: how many folks got adjuvant radiation?
table(ncdb$rad_cyst, ncdb$RX_SUMM_RADIATION)
ncdb$cyst_rad <- 0
ncdb$cyst_rad[ncdb$rad_cyst == "1" & ncdb$RX_SUMM_RADIATION == "1: Beam radiation"] <- 1
ncdb$cyst_rad <- factor(ncdb$cyst_rad)

summary(ncdb$cyst_rad)


##### individual pairwise propensity score models #####

### neoadjuvant
require(MatchIt)

neo_match <- subset(ncdb, ncdb$treatment == "Neoadjuvant Chemo" | ncdb$treatment == "Cystectomy Alone")
neo_match <- neo_match[,c("AGE", "SEX", "CDCC_TOTAL", "TNM_CLIN_T", "PUF_VITAL_STATUS", "treatment", "DX_LASTCONTACT_DEATH_MONTHS")]

neo_match$treatment <- factor(neo_match$treatment)

neo_match <- na.omit(neo_match)
neo_match$neo <- 0
neo_match$neo[neo_match$treatment == "Neoadjuvant Chemo"] <- 1
neo_matched <- matchit(neo ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T, 
                       data=neo_match, method = "nearest", ratio=2)
summary(neo_matched)
plot(neo_matched, type = "hist")

postmatch_neo <- match.data(neo_matched)
postmatch_neo$surv <- Surv(postmatch_neo$DX_LASTCONTACT_DEATH_MONTHS, postmatch_neo$PUF_VITAL_STATUS == "0: Dead")

neo_model <- coxph(surv ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T + treatment, data=postmatch_neo)
summary(neo_model)

### adjuvant
adj_match <- subset(ncdb, ncdb$treatment == "Adjuvant Chemo" | ncdb$treatment == "Cystectomy Alone")
adj_match <- adj_match[,c("AGE", "SEX", "CDCC_TOTAL", "TNM_CLIN_T", "PUF_VITAL_STATUS", "treatment", "DX_LASTCONTACT_DEATH_MONTHS")]

adj_match$treatment <- factor(adj_match$treatment)

adj_match <- na.omit(adj_match)
adj_match$adj <- 0
adj_match$adj[adj_match$treatment == "Adjuvant Chemo"] <- 1
adj_matched <- matchit(adj ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T, 
                       data=adj_match, method = "nearest", ratio=2)
summary(adj_matched)
plot(adj_matched, type = "hist")

postmatch_adj <- match.data(adj_matched)
postmatch_adj$surv <- Surv(postmatch_adj$DX_LASTCONTACT_DEATH_MONTHS, postmatch_adj$PUF_VITAL_STATUS == "0: Dead")

adj_model <- coxph(surv ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T + treatment, data=postmatch_adj)
summary(adj_model)

### chemotherapy
chemo_match <- subset(ncdb, ncdb$treatment == "Chemo Alone" | ncdb$treatment == "Cystectomy Alone")
chemo_match <- chemo_match[,c("AGE", "SEX", "CDCC_TOTAL", "TNM_CLIN_T", "PUF_VITAL_STATUS", "treatment", "DX_LASTCONTACT_DEATH_MONTHS")]

chemo_match$treatment <- factor(chemo_match$treatment)

chemo_match <- na.omit(chemo_match)
chemo_match$chemo <- 0
chemo_match$chemo[chemo_match$treatment == "Chemo Alone"] <- 1
chemo_matched <- matchit(chemo ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T, 
                       data=chemo_match, method = "nearest", ratio=2)
summary(chemo_matched)
plot(chemo_matched, type = "hist")

postmatch_chemo <- match.data(chemo_matched)
postmatch_chemo$surv <- Surv(postmatch_chemo$DX_LASTCONTACT_DEATH_MONTHS, postmatch_chemo$PUF_VITAL_STATUS == "0: Dead")

chemo_model <- coxph(surv ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T + treatment, data=postmatch_chemo)
summary(chemo_model)


### radiation alone
rads_match <- subset(ncdb, ncdb$treatment == "Radiation Alone" | ncdb$treatment == "Cystectomy Alone")
rads_match <- rads_match[,c("AGE", "SEX", "CDCC_TOTAL", "TNM_CLIN_T", "PUF_VITAL_STATUS", "treatment", "DX_LASTCONTACT_DEATH_MONTHS")]

rads_match$treatment <- factor(rads_match$treatment)

rads_match <- na.omit(rads_match)
rads_match$rads <- 0
rads_match$rads[rads_match$treatment == "Radiation Alone"] <- 1
rads_matched <- matchit(rads ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T, 
                         data=rads_match, method = "nearest", ratio=2)
summary(rads_matched)
plot(rads_matched, type = "hist")

postmatch_rads <- match.data(rads_matched)
postmatch_rads$surv <- Surv(postmatch_rads$DX_LASTCONTACT_DEATH_MONTHS, postmatch_rads$PUF_VITAL_STATUS == "0: Dead")

rads_model <- coxph(surv ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T + treatment, data=postmatch_rads)
summary(rads_model)

### chemorads
chemorads_match <- subset(ncdb, ncdb$treatment == "Chemorads" | ncdb$treatment == "Cystectomy Alone")
chemorads_match <- chemorads_match[,c("AGE", "SEX", "CDCC_TOTAL", "TNM_CLIN_T", "PUF_VITAL_STATUS", "treatment", "DX_LASTCONTACT_DEATH_MONTHS")]

chemorads_match$treatment <- factor(chemorads_match$treatment)

chemorads_match <- na.omit(chemorads_match)
chemorads_match$chemorads <- 0
chemorads_match$chemorads[chemorads_match$treatment == "Chemorads"] <- 1
chemorads_matched <- matchit(chemorads ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T, 
                        data=chemorads_match, method = "nearest", ratio=2)
summary(chemorads_matched)
plot(chemorads_matched, type = "hist")

postmatch_chemorads <- match.data(chemorads_matched)
postmatch_chemorads$surv <- Surv(postmatch_chemorads$DX_LASTCONTACT_DEATH_MONTHS, postmatch_chemorads$PUF_VITAL_STATUS == "0: Dead")

chemorads_model <- coxph(surv ~ AGE + SEX + CDCC_TOTAL + TNM_CLIN_T + treatment, data=postmatch_chemorads)
summary(chemorads_model)


#### table 1 re-do #####
tapply(ncdb$AGE, ncdb$treatment, mean)

table(ncdb$SEX, ncdb$treatment)
prop.table(table(ncdb$SEX, ncdb$treatment),2)

table(ncdb$FACILITY_TYPE_CD, ncdb$treatment)
prop.table(table(ncdb$FACILITY_TYPE_CD, ncdb$treatment),2)

ncdb$race_rc <- 0
ncdb$race_rc[ncdb$RACE == "1: White"] <- 1
ncdb$race_rc[ncdb$RACE == "2: Black"] <- 2
ncdb$race_rc <- factor(ncdb$race_rc, levels=c(0:2), labels=c("Other", "White", "Black"))

table(ncdb$race_rc, ncdb$treatment)
prop.table(table(ncdb$race_rc, ncdb$treatment),2)

table(ncdb$INSURANCE_STATUS, ncdb$treatment)
prop.table(table(ncdb$INSURANCE_STATUS, ncdb$treatment),2)

table(ncdb$MED_INC_QUAR_12, ncdb$treatment)
prop.table(table(ncdb$MED_INC_QUAR_12, ncdb$treatment),2)


table(ncdb$CDCC_TOTAL, ncdb$treatment)
prop.table(table(ncdb$CDCC_TOTAL, ncdb$treatment),2)

table(ncdb$YEAR_OF_DIAGNOSIS, ncdb$treatment)
prop.table(table(ncdb$YEAR_OF_DIAGNOSIS, ncdb$treatment),2)

table(ncdb$TNM_CLIN_T, ncdb$treatment)
prop.table(table(ncdb$TNM_CLIN_T, ncdb$treatment),2)


## how many events?
table(ncdb$PUF_VITAL_STATUS, ncdb$treatment)






### surgical approach
ncdb$RX_HOSP_SURG_APPR_2010 <- factor(ncdb$RX_HOSP_SURG_APPR_2010)
table(ncdb$RX_HOSP_SURG_APPR_2010, ncdb$YEAR_OF_DIAGNOSIS)
prop.table(table(ncdb$RX_HOSP_SURG_APPR_2010[ncdb$rad_cyst == 1]))
