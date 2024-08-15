##########################################################
# FC harmonization and plot
# Lanxin Ji for Birth transition paper R1 
# 13 Aug 2024
###########################################################

#################################
# load basic library
packages <- c("utils", "tidyverse", "downloadthis",                             # packages for data management
              "foreign", "MplusAutomation",                                     # packages for writing data         
              "sjPlot", "broom", "kableExtra",                                  # packages for generating tables
              "nlme", "lme4", "lmerTest", "stats","lmer","lmerTest",            # packages for MLMs
              "mgcv", "gamm4", "itsadug", "splines","npreg",                    # packages for GAMMs
              "lavaan",                                                         # packages for SEMs
              "ggplot2", "semPlot", "ggeffects", "ggstatsplot", "RColorBrewer", # packages for visualization   
              "interactions")                                                   # packages for probing interactions
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, 
                           rownames(installed.packages())), 
                   repos = "http://cran.us.r-project.org")
}
invisible(lapply(packages, 
                 library, 
                 character.only = TRUE))



#################################
# install longCombat package
#################################
# install.packages('devtools')
# devtools::install_github("jcbeer/longCombat")

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

data_age = read.csv("/Users/jil02/Documents/Projects_Working/Ji_Birth_transition/R1_analysis/Combat_analysis/Data_for_R_n184_final.csv")
filtered_out_data <- subset(data_age, abs(mean_pry)>0.5 | abs(mean_xyz)>0.5 | abs(mean_max_pry)>1 | abs(mean_max_xyz)>1)
filtered_data <- subset(data_age, abs(mean_pry)<=0.5 & abs(mean_xyz)<=0.5 & abs(mean_max_pry)<=1 & abs(mean_max_xyz)<=1)
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
# make batch boxplot for all_pos_fc, do not adjust for batch 
batchBoxplot(idvar='ID', 
             batchvar='sequence.ID', 
             feature='all_pos_fc', 
             formula='GA_scan',
             ranef='(1|ID)',
             data=filtered_data,
             colors=1:5)

batchBoxplot(idvar='ID', 
             batchvar='sequence.ID', 
             feature='all_neg_fc', 
             formula='GA_scan',
             ranef='(1|ID)',
             data=filtered_data,
             colors=1:5)

# make batch boxplot for all_pos_fc, DO adjust for batch 
# order by increasing batch variance
# (centers boxplot means on the zero line)
batchBoxplot(idvar='ID', 
             batchvar='sequence.ID', 
             feature='all_neg_fc', 
             formula='GA_scan',
             ranef='(1|ID)',
             data=filtered_data,
             adjustBatch=TRUE,
             orderby='var',
             colors=1:5)

#################################
# trajPlot() -- visualize trajectories
#################################
# for everyone
trajPlot(idvar='ID', 
         timevar='GA_scan',
         feature='all_pos_fc', 
         batchvar='sequence.ID',  
         data=filtered_data,
         point.col=filtered_data$sequence.ID)

trajPlot(idvar='ID', 
         timevar='GA_scan',
         feature='all_neg_fc', 
         batchvar='sequence.ID',  
         data=filtered_data,
         point.col=filtered_data$sequence.ID)

#################################
# addTest() -- test for additive scanner effects
#################################
featurenames <- c('all_pos_fc','all_neg_fc','GE_pos', 'GE_neg', 'LE_pos')
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

df <- df[, !(names(df) %in% c("GA.at.Birth", "Preterm", "fetal_pos", "fetal_neg", "infant_pos", "infant_neg", "Fetal_var", "Infant_var"))] 
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
# save combat feature names
featurenames.combat <- names(filtered_data_harmonized)[4:8]
# merge with original dataframe
filtered_data <- merge(filtered_data, filtered_data_harmonized[,c(1,2,4:8)], by=c('ID', 'GA_scan'))

#################################
# test for additive scanner effects in combatted data
#################################
addTestTableCombat <- addTest(idvar='ID', 
                              batchvar='sequence.ID', 
                              features=featurenames.combat, 
                              formula='GA_scan',
                              ranef='(1|ID)',
                              data=filtered_data)

# there are still some significant additive batch effects (p<0.05)
# but p-values tend to be larger (-log10(p-values are smaller)) in overall distribution
boxplot(-log(as.numeric(addTestTable$`KR p-value`), base=10),
        -log(as.numeric(addTestTableCombat$`KR p-value`), base=10),
        ylim=c(0, 18),
        las=1,
        ylab='additive batch effect -log10(p-value)',
        names=c('before ComBat', 'after ComBat'))

####################################
# plot features before and after combat
#####!!!!!!!!! should we adjust batch in box plots??????
#######################################
par(mfrow=c(1,2))
batchBoxplot(idvar='ID', 
             batchvar='sequence.ID', 
             feature='all_pos_fc', 
             formula='GA_scan',
             ranef='(1|ID)',
             data=filtered_data,
             colors=1:5,
             adjustBatch=TRUE,
             orderby='var',
             title='pos_fc before ComBat')

# check feature3 boxplot after combat
batchBoxplot(idvar='ID', 
             batchvar='sequence.ID', 
             feature='all_pos_fc.combat', 
             formula='GA_scan',
             ranef='(1|ID)',
             data=filtered_data,
             colors=1:5,
             adjustBatch=TRUE,
             orderby='var',
             title='pos_fc after ComBat')

###########################################

#################################
# plot trajectories before and after combat
#################################
par(mfrow=c(1,2))
ggplot(filtered_data,aes(x=GA_scan,y=all_pos_fc)) +
  geom_point(aes(shape=as.factor(Fetal)), alpha=.8)+
  geom_line(aes(group=ID),alpha=.2)+
  geom_smooth(method='gam',color='firebrick')+
  theme_bw()+ylim(c(-0.1,0.4))

ggplot(filtered_data,aes(x=GA_scan,y=all_pos_fc.combat)) +
  geom_point(aes(shape=as.factor(Fetal)), alpha=.8)+
  geom_line(aes(group=ID),alpha=.2)+
  geom_smooth(method='gam',color='firebrick')+
  theme_bw()+ylim(c(-0.1,0.4))

ggplot(filtered_data,aes(x=GA_scan,y=all_neg_fc)) +
  geom_point(aes(shape=as.factor(Fetal)), alpha=.8)+
  geom_line(aes(group=ID),alpha=.2)+
  geom_smooth(method='gam',color='deepskyblue2')+
  theme_bw()+ylim(c(-0.1,0.4))

ggplot(filtered_data,aes(x=GA_scan,y=all_neg_fc.combat)) +
  geom_point(aes(shape=as.factor(Fetal)),alpha=.8)+
  geom_line(aes(group=ID),alpha=.2)+
  geom_smooth(method='gam',color='deepskyblue2')+
  theme_bw()+ylim(c(-0.1,0.4))

#################################
# plot trajectories for graph measure after combat
#################################

ggplot(filtered_data,aes(x=GA_scan,y=GE_pos.combat)) +
  geom_point(aes(shape=as.factor(Fetal)), alpha=.8)+
  geom_line(aes(group=ID),alpha=.2)+
  geom_smooth(method='gam',color='firebrick')+
  theme_bw()+ylim(c(0.4,0.55))

ggplot(filtered_data,aes(x=GA_scan,y=LE_pos.combat)) +
  geom_point(aes(shape=as.factor(Fetal)), alpha=.8)+
  geom_line(aes(group=ID),alpha=.2)+
  geom_smooth(method='gam',color='firebrick')+
  theme_bw()+ylim(c(0.5,0.9))

ggplot(filtered_data,aes(x=GA_scan,y=GE_neg.combat)) +
  geom_point(aes(shape=as.factor(Fetal)), alpha=.8)+
  geom_line(aes(group=ID),alpha=.2)+
  geom_smooth(method='gam',color='navy')+
  theme_bw()+ylim(c(0.4,0.7))


##################################
# plot heatmap for change rate
#################################

#gamm fit on all_pos_fc.combat
gamm <- gamm4(all_pos_fc.combat ~ 1 + s(GA_scan),
              random = ~ (1 | ID),
              data = filtered_data)
summary(gamm$gam)

plot.gam(gamm$gam, se = TRUE, rug = TRUE, shade = TRUE,
         xlab = "Age", ylab = "Fitted FC Values")

# infant: seq(42, 55, length.out = 65)
# cross-birth seq(25, 55, length.out = 150)
# fetal: seq(24, 40, length.out = 80)
x0 <- seq(25, 55, length.out = 150)
y_pred <- predict(gamm$gam, newdata=data.frame(GA_scan = x0))
dy_dx <- diff(y_pred) / diff(x0)

# Convert the matrix to a data frame
df <- data.frame(values = dy_dx, position = 1:149)
df <- drop_na(df)
ggplot(df , aes(x = position, y = 1)) +
  geom_raster(aes(fill = values), interpolate=TRUE) +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(-0.05, 0.05)) +
  theme_classic()

#gamm fit on all_neg_fc.combat
gamm <- gamm4(all_neg_fc.combat ~ 1 + s(GA_scan),
              random = ~ (1 | ID),
              data = filtered_data)
summary(gamm$gam)
x0 <- seq(25, 55, length.out = 150)
y_pred <- predict(gamm$gam, newdata=data.frame(GA_scan = x0))
dy_dx <- diff(y_pred) / diff(x0)
df <- data.frame(values = dy_dx, position = 1:149)
df <- drop_na(df)
ggplot(df , aes(x = position, y = 1)) +
  geom_raster(aes(fill = values), interpolate=TRUE) +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(-0.05, 0.05)) +
  theme_classic()

#gamm fit on GE_pos.combat
gamm <- gamm4(GE_pos.combat ~ 1 + s(GA_scan),
              random = ~ (1 | ID),
              data = filtered_data)
summary(gamm$gam)
x0 <- seq(25, 55, length.out = 150)
y_pred <- predict(gamm$gam, newdata=data.frame(GA_scan = x0))
dy_dx <- diff(y_pred) / diff(x0)
df <- data.frame(values = dy_dx, position = 1:149)
df <- drop_na(df)
ggplot(df , aes(x = position, y = 1)) +
  geom_raster(aes(fill = values), interpolate=TRUE) +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(-0.02, 0.02)) +
  theme_classic()

#gamm fit on GE_neg.combat
gamm <- gamm4(GE_neg.combat ~ 1 + s(GA_scan),
              random = ~ (1 | ID),
              data = filtered_data)
summary(gamm$gam)
x0 <- seq(25, 55, length.out = 150)
y_pred <- predict(gamm$gam, newdata=data.frame(GA_scan = x0))
dy_dx <- diff(y_pred) / diff(x0)
df <- data.frame(values = dy_dx, position = 1:149)
df <- drop_na(df)
ggplot(df , aes(x = position, y = 1)) +
  geom_raster(aes(fill = values), interpolate=TRUE) +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(-0.02, 0.02)) +
  theme_classic()

#gamm fit on LE_pos.combat
gamm <- gamm4(LE_pos.combat ~ 1 + s(GA_scan),
              random = ~ (1 | ID),
              data = filtered_data)
summary(gamm$gam)
x0 <- seq(25, 55, length.out = 150)
y_pred <- predict(gamm$gam, newdata=data.frame(GA_scan = x0))
dy_dx <- diff(y_pred) / diff(x0)
df <- data.frame(values = dy_dx, position = 1:149)
df <- drop_na(df)
ggplot(df , aes(x = position, y = 1)) +
  geom_raster(aes(fill = values), interpolate=TRUE) +
  scale_fill_gradient2(low="navy", mid="white", high="red", 
                       midpoint=0, limits=c(-0.02, 0.02)) +
  theme_classic()
