library(R.matlab)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(cowplot)

#matfile = readMat("/Users/jil02/Dropbox (NYU Langone Health)/TANYA_LANXIN/birth_transition_project_0923/data_and_notes/network_FC_vec_new.mat")
matfile = readMat("/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/network_FC_vec_new.mat")

df = matfile$network.avg
df_new= as.data.frame(df)
df_t = t(df_new) # x axis is network pairs, y axis is subject
network_avg = df_t

#age_data = read.csv('/Users/jil02/Dropbox (NYU Langone Health)/TANYA_LANXIN/birth_transition_project_0923/data_and_notes/Data_for_R_n184_final.csv')
age_data = read.csv('/Users/tb2322/NYU Langone Health Dropbox/Tanya Bhatia/TANYA_LANXIN/birth_transtion_project_0824/data_and_notes/Data_for_R_n184_final.csv')
conn_age_df = cbind.data.frame(age_data, network_avg)

filtered_data <- subset(conn_age_df, abs(mean_pry)<=0.5 & abs(mean_xyz)<=0.5 & abs(mean_max_pry)<=1 & abs(mean_max_xyz)<=1)
col_name <- as.character(i)

### gamm fit
filtered_data$ID <- as.factor(filtered_data$ID)

colnames(filtered_data)[colnames(filtered_data) == "36"] <- "var36"
gamm <- gamm4(var36 ~ 1 + s(GA_scan),
              random = ~ (1 | ID),
              data = filtered_data)
summary(gamm$gam)

plot.gam(gamm$gam, se = TRUE, rug = TRUE, shade = TRUE,
         xlab = "Age", ylab = "Fitted Modularity Values")


x0 <- seq(25, 55, length.out = 150)
y_pred <- predict(gamm$gam, newdata=data.frame(GA_scan = x0))
dy_dx <- diff(y_pred) / diff(x0)


# Convert the matrix to a data frame
df <- data.frame(values = dy_dx, position = 1:149)
df <- drop_na(df)
df$values[df$values > 0.004] <- 0.004
df$values[df$values < -0.004] <- -0.004
# Create a ggplot for the current subplot
ggplot(df , aes(x = position, y = 1)) +
  geom_raster(aes(fill = values), interpolate=TRUE) +
  scale_fill_gradient2(low="navy", mid="white", high="red", midpoint=0, limits=c(-0.004,0.004)) +
  theme_classic()


