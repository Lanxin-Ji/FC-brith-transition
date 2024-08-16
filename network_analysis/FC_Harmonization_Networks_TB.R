##########################################################
# FC harmonization for network data
# Lanxin Ji for Birth transition paper R1 
# Last Edited: 13 Aug 2024 by Tanya Bhatia
###########################################################

#################################
# install longCombat package
#################################
install.packages('devtools')
devtools::install_github("jcbeer/longCombat")

install.packages("Matrix", dependencies = TRUE)

#################################
# load longCombat package
#################################
library(longCombat)
# check documentation
?longCombat

#################################
# install and load invgamma & lmer package
#################################
# install.packages('invgamma')
library(invgamma)
library(lme4)
library(ggplot2)

# data_age = read.csv("/Users/jil02/Documents/Projects_Working/Ji_Birth_transition/R1_analysis/Combat_analysis/Data_for_R_n184_final.csv")
# data_age_for_combat = read.csv("/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/Data_for_R_n184_final_batch.csv")
# filtered_out_data <- subset(data_age_for_combat, abs(mean_pry)>0.5 | abs(mean_xyz)>0.5 | abs(mean_max_pry)>1 | abs(mean_max_xyz)>1)
# filtered_data <- subset(data_age_for_combat, abs(mean_pry)<=0.5 & abs(mean_xyz)<=0.5 & abs(mean_max_pry)<=1 & abs(mean_max_xyz)<=1)
# filtered_data$ID <- as.factor(filtered_data$ID)
# filtered_data <- filtered_data[!is.na(filtered_data$GA_scan), ]

# Load data containing connectivity values for each network pair for each subject (rows:subjects, cols: subject ID, GA at scan, sequence.ID, and network pairs)
filtered_data = read.csv("/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/network_avg_df_input_for_combat.csv")
filtered_data$ID <- as.factor(filtered_data$ID)
filtered_data <- filtered_data[!is.na(filtered_data$GA_scan), ]
#################################
# batchTimeViz() -- visualize change in batch over time
#################################
batchTimeViz(batchvar='sequence.ID',
             timevar='GA_scan',
             data=filtered_data)


#################################
# batchBoxplot() -- to visualize residuals across batches
# can do for each feature you are interested in
#################################

# make batch boxplot for V1 (network?), do not adjust for batch 
batchBoxplot(idvar='ID', 
             batchvar='sequence.ID', 
             feature='V1', 
             formula='GA_scan',
             ranef='(1|ID)',
             data=filtered_data,
             colors=1:5)

# make batch boxplot for V1 (network?), DO adjust for batch 
batchBoxplot(idvar='ID', 
             batchvar='sequence.ID', 
             feature='V1', 
             formula='GA_scan',
             ranef='(1|ID)',
             data=filtered_data,
             adjustBatch=TRUE,
             orderby='var',
             colors=1:5)

#################################
# addTest() -- test for additive scanner effects
#################################
featurenames <- paste0('V',1:36)

addTestTable <- addTest(idvar='ID', 
                        batchvar='sequence.ID', 
                        features=featurenames, 
                        formula='GA_scan',
                        ranef='(1|ID)',
                        data=filtered_data)

#################################
# multTest() -- test for multiplicative scanner effects
#################################
multTestTable <- multTest(idvar='ID', 
                          batchvar='sequence.ID', 
                          features=featurenames, 
                          formula='GA_scan',
                          ranef='(1|ID)',
                          data=filtered_data)


#################################
# longCombat() -- apply longitudinal ComBat
#################################
df <- filtered_data
#featurenames <- paste0('V',1:36)

df <- df[, !(names(df) %in% c("V_ID"))] 
#df <- df[, !(names(df) %in% c("GA.at.Birth", "Preterm", "fetal_pos", "fetal_neg", "infant_pos", "infant_neg", "Fetal_var", "Infant_var"))] 
df_combat <- longCombat(idvar='ID', 
                        timevar='GA_scan',
                        batchvar='sequence.ID', 
                        features=featurenames, 
                        formula='GA_scan',
                        ranef='(1|ID)',
                        data=df)

#################################
# get the harmonized data
filtered_data_harmonized <- df_combat$data_combat
# write harmonized data to csv
write.csv(filtered_data_harmonized,"/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/networks_data_combat_output.csv")
# save combat feature names
featurenames.combat <- names(filtered_data_harmonized)[4:ncol(filtered_data_harmonized)]
# merge with original dataframe
filtered_data <- merge(filtered_data, filtered_data_harmonized, by=c('ID', 'GA_scan','sequence.ID'))

#################################
# test for additive scanner effects in combatted data
#################################
addTestTableCombat <- addTest(idvar='ID', 
                              batchvar='sequence.ID', 
                              features=featurenames.combat, 
                              formula='GA_scan',
                              ranef='(1|ID)',
                              data=filtered_data)

boxplot(-log(as.numeric(addTestTable$`KR p-value`), base=10),
        -log(as.numeric(addTestTableCombat$`KR p-value`), base=10),
        ylim=c(0, 18),
        las=1,
        ylab='additive batch effect -log10(p-value)',
        names=c('before ComBat', 'after ComBat'))


###########################################


