#This script was used to construct the SEM models and subsequently predict the effects of climate change scenarios. 
#piecewiseSEM 

library(piecewiseSEM)
library(dplyr)
library(nlme)

#first we imported the table with the predictors and response variables, which were the  percentage of plant-growth promoting taxa (in fraction from 0 to 1), and the Shannon diversity. These two response varibles were in the first two columns of the table (we import two seperate tables for Bacteria and Fungi).
dataset1 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")

#We proceeded to log the table to normalize the data  (we specify a range of columns to log-transform because the last column contains non-continuous data)
dataset2 = log(dataset1[,1:20] + 1)

#We then construct the first model, using latitude and longitude for spatial autocorrection, and country for random effect (these data are also in the imported table)

#Calculate all the direct correlations in the model
direct = lme(Shannon ~  C_Percent + Na + PGPT + EVI2 + pH, random = ~ 1 | Country, correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes),  data = dataset2)
direct2 = lme(PGPT ~  pH + EVI2 + N_Percent , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct3 = lme(pH ~  EVI2 + N_Percent + Mean_Ann_Temp, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct4 = lme(Na ~ Mean_Ann_Temp , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct5 = lme(C_Percent ~  Mean_Ann_Temp = Mean_Ann_Precip  , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct6 = lme(EVI2 ~  Mean_Ann_Temp + Mean_Ann_Precip , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct7 = lme(N_Percent ~  Mean_Ann_Temp + Mean_Ann_Precip, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)

#After calculating all the direct correlations in the model, we calculate the combined model. In this case, we excluded the correlations between C_Percent, N_Percent, and pH (represented by %~~%), because we do not see them as true correlations in the model.  
modelfinalfung = psem(direct,direct2,direct3, direct4,direct5, direct6,  direct7, N_Percent %~~% C_Percent,pH %~~% C_Percent)
model_final_data_fung = summary(modelfinalfung, standardize = "scale")

#Once the model is calculated, piecewiseSEM suggests additional correlations that will increase the predictability of the model (this is represented by increasing the p-value of the model and decreasing the AIC). 
#After these additional correlations are accounted form, we consider the model to be optimal (i.e. with the lowest AIC and non-significant p-value)

#Optimal model for Bacteria 
direct = lme(Shannon ~   C_Percent + pH + Mean_Ann_Temp  , random = ~ 1 | Country, correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes),  data = dataset2)
direct2 = lme(PGPR ~  pH + Shannon + Na + Mean_Ann_Precip + Mean_Ann_Temp + N_Percent, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct3 = lme(pH ~   Mean_Ann_Precip + Na + N_Percent + Mean_Ann_Temp + EVI2, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct4 = lme(Na ~ Mean_Ann_Temp , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct5 = lme(C_Percent ~  Mean_Ann_Temp + Mean_Ann_Precip + Na  , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct6 = lme(EVI2 ~ Mean_Ann_Precip, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct7 = lme(N_Percent ~  Mean_Ann_Temp + Mean_Ann_Precip + Na, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)

summary(direct)
model_final_bact = psem(direct,direct2,direct3, direct4, direct5, direct6,  direct7, N_Percent %~~% C_Percent)
model_final_data = summary(model_final_bact, standardize = "scale")
model_final_data$coefficients

#Calculate total standardized effect of Mean Annual Temperature (MAT) and Mean Annual Precipitation (MAP) on Shannon diversity and abundance of PGP taxa. 
#total effect of MAP on PGPB
MAP_PGPR_total = model_final_data$coefficients[10,8] + (model_final_data$coefficients[25,8]* model_final_data$coefficients[15,8]* model_final_data$coefficients[7,8]) + (model_final_data$coefficients[13,8]* model_final_data$coefficients[7,8]) + (model_final_data$coefficients[23,8]* model_final_data$coefficients[12,8]*model_final_data$coefficients[7,8])

#total effect of MAP on Shannon
MAP_Shannon_total = (model_final_data$coefficients[19,8]* model_final_data$coefficients[1,8]) + (model_final_data$coefficients[13,8]* model_final_data$coefficients[2,8]) + (model_final_data$coefficients[25,8]* model_final_data$coefficients[15,8]* model_final_data$coefficients[2,8]) + (model_final_data$coefficients[23,8]* model_final_data$coefficients[12,8]* model_final_data$coefficients[2,8])

#total effect of MAT on PGPB
MAT_PGPR_total = (model_final_data$coefficients[11,8]) + (model_final_data$coefficients[16,8]* model_final_data$coefficients[14,8]* model_final_data$coefficients[7,8]) + (model_final_data$coefficients[26,8]* model_final_data$coefficients[15,8]* model_final_data$coefficients[7,8]) 

#total effect of MAT on Shannon
MAT_Shannon_total = (model_final_data$coefficients[26,8]* model_final_data$coefficients[15,8]* model_final_data$coefficients[2,8]) + (model_final_data$coefficients[16,8]* model_final_data$coefficients[14,8]* model_final_data$coefficients[2,8]) + (model_final_data$coefficients[1,8]*model_final_data$coefficients[20,8])



#Optimal model for Fungi 
direct = lme(Shannon_Fungi ~  Mean_Ann_Precip + Na + PGPF + Mean_Annual_Temp, random = ~ 1 | Country, correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes),  data = dataset2)
direct2 = lme(PGPF ~  pH + EVI2 + Na + Mean_Ann_Precip , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct3 = lme(pH ~  EVI2 + Mean_Ann_Precip + Na + N_Percent + Mean_Ann_Temp, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct4 = lme(Na ~ Mean_Ann_Temp , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct5 = lme(C_Percent ~  Mean_Ann_Precip + Mean_Ann_Temp + Na  , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct6 = lme(EVI2 ~  Mean_Ann_Precip , random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)
direct7 = lme(N_Percent ~  Mean_Ann_Precip + Mean_Ann_Temp + Na, random = ~ 1 | Country,correlation = corSpatial(form = ~ Corrected_Longitude + Corrected_Latitudes), data = dataset2)


modelfinalfung = psem(direct,direct2,direct3, direct4,direct5, direct6,  direct7, N_Percent %~~% C_Percent,pH %~~% C_Percent)
model_final_data_fung = summary(modelfinalfung, standardize = "scale")
model_final_data_fung$coefficients

#total effect of MAP on PGPF
MAP_PGPF_total_fungi = model_final_data_fung$coefficients[8,8] + (model_final_data_fung$coefficients[20,8]*model_final_data_fung$coefficients[6,8]) + (model_final_data_fung$coefficients[23,8]*model_final_data_fung$coefficients[13,8]*model_final_data_fung$coefficients[5,8])
#total effect of MAP on Shannon
MAP_Shannon_total_fungi = model_final_data_fung$coefficients[2,8] + (model_final_data_fung$coefficients[24,8]*model_final_data_fung$coefficients[13,8]*model_final_data_fung$coefficients[5,8]) + (model_final_data_fung$coefficients[23,8]*model_final_data_fung$coefficients[13,8]*model_final_data_fung$coefficients[5,8]*model_final_data_fung$coefficients[4,8])


#total effect of MAT on PGPF
MAP_PGPF_total_fungi = model_final_data_fung$coefficients[9,8] + (model_final_data_fung$coefficients[20,8]*model_final_data_fung$coefficients[6,8]*model_final_data_fung$coefficients[4,8]) + (model_final_data_fung$coefficients[23,8]*model_final_data_fung$coefficients[13,8]*model_final_data_fung$coefficients[5,8]*model_final_data_fung$coefficients[4,8])
#total effect of MAT on Shannon
MAT_Shannon_total = (model_final_data$coefficients[26,8]* model_final_data$coefficients[15,8]* model_final_data$coefficients[2,8]) + (model_final_data$coefficients[16,8]* model_final_data$coefficients[14,8]* model_final_data$coefficients[2,8]) + (model_final_data$coefficients[1,8]*model_final_data$coefficients[20,8])


######################################################################################################################################
#After the models are constructed, we used them to predect the efect of different future MAT and MAP values on the response variables. The MAT and MAP values were estimated with the MIROC6 model, for two time periods (2040-2060, 2080-2100) and two scenarios (SSP126, SSP585)


#Predictions using the Bacterial SEM 
#timne period 2080-2100, SSP585
#Import the table with estimated MAP and MAT values 
dataset3 = read.table(file.choose(), header = TRUE, row.names = 1, sep = ",")

#create a sub-table with just the values for a specific time period and scenario
predictdata_100_585 = data.frame(Mean_Ann_Temp= dataset3$X2100_585_b1, Mean_Ann_Precip = dataset3$X2100_585_b12)

#Predictions
Na_100_585_2 = as.data.frame( predict(direct4,newdata = predictdata_100_585, level = 0))
predict_100_585_Na= cbind(predictdata_100_585, "Na" = Na_100_585_2[,1])
EVI2_100_585_2 = as.data.frame( predict(direct6,newdata = predictdata_100_585, level = 0))
predict_100_585_Na_EVI2 = cbind(predict_100_585_Na, "EVI2" = EVI2_100_585_2[,1])
N_Percent_100_585 = as.data.frame(predict(direct7, newdata = predict_100_585_Na_EVI2, level = 0))
predict_100_585_Na_EVI2_N = cbind(predict_100_585_Na_EVI2, "N_Percent" = N_Percent_100_585[,1])
C_Percent_100_585 = as.data.frame(predict(direct5, newdata = predict_100_585_Na_EVI2_N, level = 0))
predict_100_585_Na_EVI2_N_C = cbind(predict_100_585_Na_EVI2_N, "C_Percent" = C_Percent_100_585[,1])
pH_100_585 = as.data.frame(predict(direct3, newdata = predict_100_585_Na_EVI2_N_C, level = 0))
predict_100_585_Na_EVI2_N_C_pH = cbind(predict_100_585_Na_EVI2_N_C , "pH" = pH_100_585[,1])
Shannon_100_585 = as.data.frame(predict(direct, newdata = predict_100_585_Na_EVI2_N_C_pH, level = 0))
predict_100_585_Na_EVI2_N_C_pH_Shannon = cbind(predict_100_585_Na_EVI2_N_C_pH , "Shannon" = Shannon_100_585[,1])
PGPB_100_585 = as.data.frame(predict(direct2, newdata = predict_100_585_Na_EVI2_N_C_pH_Shannon, level = 0))
predict_100_585_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_100_585_Na_EVI2_N_C_pH_Shannon , "PGPB" = PGPB_100_585[,1])
write.csv(predict_100_585_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_bact_100years_585.csv")

#timne period 2040-2060, SSP585
predictdata_50_585 = data.frame(Mean_Ann_Temp= dataset3$X2050_585_b1, Mean_Ann_Precip = dataset3$X2050_585_b12)

Na_50_585_2 = as.data.frame( predict(direct4,newdata = predictdata_50_585, level = 0))
predict_50_585_Na= cbind(predictdata_50_585, "Na" = Na_50_585_2[,1])
EVI2_50_585_2 = as.data.frame( predict(direct6,newdata = predictdata_50_585, level = 0))
predict_50_585_Na_EVI2 = cbind(predict_50_585_Na, "EVI2" = EVI2_50_585_2[,1])
N_Percent_50_585 = as.data.frame(predict(direct7, newdata = predict_50_585_Na_EVI2, level = 0))
predict_50_585_Na_EVI2_N = cbind(predict_50_585_Na_EVI2, "N_Percent" = N_Percent_50_585[,1])
C_Percent_50_585 = as.data.frame(predict(direct5, newdata = predict_50_585_Na_EVI2_N, level = 0))
predict_50_585_Na_EVI2_N_C = cbind(predict_50_585_Na_EVI2_N, "C_Percent" = C_Percent_50_585[,1])
pH_50_585 = as.data.frame(predict(direct3, newdata = predict_50_585_Na_EVI2_N_C, level = 0))
predict_50_585_Na_EVI2_N_C_pH = cbind(predict_50_585_Na_EVI2_N_C , "pH" = pH_50_585[,1])
Shannon_50_585 = as.data.frame(predict(direct, newdata = predict_50_585_Na_EVI2_N_C_pH, level = 0))
predict_50_585_Na_EVI2_N_C_pH_Shannon = cbind(predict_50_585_Na_EVI2_N_C_pH , "Shannon" = Shannon_50_585[,1])
PGPB_50_585 = as.data.frame(predict(direct2, newdata = predict_50_585_Na_EVI2_N_C_pH_Shannon, level = 0))
predict_50_585_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_50_585_Na_EVI2_N_C_pH_Shannon , "PGPB" = PGPB_50_585[,1])
write.csv(predict_50_585_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_bact_50years_585.csv")

#timne period 2080-2100, SSP126
predictdata_100_126 = data.frame(Mean_Ann_Temp= dataset3$X2100_126_b1, Mean_Ann_Precip = dataset3$X2100_126_b12)

Na_100_126_2 = as.data.frame( predict(direct4,newdata = predictdata_100_126, level = 0))
predict_100_126_Na= cbind(predictdata_100_126, "Na" = Na_100_126_2[,1])
EVI2_100_126_2 = as.data.frame( predict(direct6,newdata = predictdata_100_126, level = 0))
predict_100_126_Na_EVI2 = cbind(predict_100_126_Na, "EVI2" = EVI2_100_126_2[,1])
N_Percent_100_126 = as.data.frame(predict(direct7, newdata = predict_100_126_Na_EVI2, level = 0))
predict_100_126_Na_EVI2_N = cbind(predict_100_126_Na_EVI2, "N_Percent" = N_Percent_100_126[,1])
C_Percent_100_126 = as.data.frame(predict(direct5, newdata = predict_100_126_Na_EVI2_N, level = 0))
predict_100_126_Na_EVI2_N_C = cbind(predict_100_126_Na_EVI2_N, "C_Percent" = C_Percent_100_126[,1])
pH_100_126 = as.data.frame(predict(direct3, newdata = predict_100_126_Na_EVI2_N_C, level = 0))
predict_100_126_Na_EVI2_N_C_pH = cbind(predict_100_126_Na_EVI2_N_C , "pH" = pH_100_126[,1])
Shannon_100_126 = as.data.frame(predict(direct, newdata = predict_100_126_Na_EVI2_N_C_pH, level = 0))
predict_100_126_Na_EVI2_N_C_pH_Shannon = cbind(predict_100_126_Na_EVI2_N_C_pH , "Shannon" = Shannon_100_126[,1])
PGPB_100_126 = as.data.frame(predict(direct2, newdata = predict_100_126_Na_EVI2_N_C_pH_Shannon, level = 0))
predict_100_126_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_100_126_Na_EVI2_N_C_pH_Shannon , "PGPB" = PGPB_100_126[,1])
write.csv(predict_100_126_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_bact_100years_126.csv")


#timne period 2040-2060, SSP126
predictdata_50_126 = data.frame(Mean_Ann_Temp= dataset3$X2050_126_b1, Mean_Ann_Precip = dataset3$X2050_126_b12)

Na_50_126_2 = as.data.frame( predict(direct4,newdata = predictdata_50_126, level = 0))
predict_50_126_Na= cbind(predictdata_50_126, "Na" = Na_50_126_2[,1])
EVI2_50_126_2 = as.data.frame( predict(direct6,newdata = predictdata_50_126, level = 0))
predict_50_126_Na_EVI2 = cbind(predict_50_126_Na, "EVI2" = EVI2_50_126_2[,1])
N_Percent_50_126 = as.data.frame(predict(direct7, newdata = predict_50_126_Na_EVI2, level = 0))
predict_50_126_Na_EVI2_N = cbind(predict_50_126_Na_EVI2, "N_Percent" = N_Percent_50_126[,1])
C_Percent_50_126 = as.data.frame(predict(direct5, newdata = predict_50_126_Na_EVI2_N, level = 0))
predict_50_126_Na_EVI2_N_C = cbind(predict_50_126_Na_EVI2_N, "C_Percent" = C_Percent_50_126[,1])
pH_50_126 = as.data.frame(predict(direct3, newdata = predict_50_126_Na_EVI2_N_C, level = 0))
predict_50_126_Na_EVI2_N_C_pH = cbind(predict_50_126_Na_EVI2_N_C , "pH" = pH_50_126[,1])
Shannon_50_126 = as.data.frame(predict(direct, newdata = predict_50_126_Na_EVI2_N_C_pH, level = 0))
predict_50_126_Na_EVI2_N_C_pH_Shannon = cbind(predict_50_126_Na_EVI2_N_C_pH , "Shannon" = Shannon_50_126[,1])
PGPB_50_126 = as.data.frame(predict(direct2, newdata = predict_50_126_Na_EVI2_N_C_pH_Shannon, level = 0))
predict_50_126_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_50_126_Na_EVI2_N_C_pH_Shannon , "PGPB" = PGPB_50_126[,1])
write.csv(predict_50_126_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_bact_50years_126.csv")



#Fungal Predictions(starting point is the same table used for the bacterial predictions, but the SEM model is the one for Fungi)
#timne period 2080-2100, SSP585
predictdata_100_585 = data.frame(Mean_Ann_Temp= dataset3$X2100_585_b1, Mean_Ann_Precip = dataset3$X2100_585_b12)

Na_100_585_2 = as.data.frame( predict(direct4,newdata = predictdata_100_585, level = 0))
predict_100_585_Na= cbind(predictdata_100_585, "Na" = Na_100_585_2[,1])
EVI2_100_585_2 = as.data.frame( predict(direct6,newdata = predictdata_100_585, level = 0))
predict_100_585_Na_EVI2 = cbind(predict_100_585_Na, "EVI2" = EVI2_100_585_2[,1])
N_Percent_100_585 = as.data.frame(predict(direct7, newdata = predict_100_585_Na_EVI2, level = 0))
predict_100_585_Na_EVI2_N = cbind(predict_100_585_Na_EVI2, "N_Percent" = N_Percent_100_585[,1])
C_Percent_100_585 = as.data.frame(predict(direct5, newdata = predict_100_585_Na_EVI2_N, level = 0))
predict_100_585_Na_EVI2_N_C = cbind(predict_100_585_Na_EVI2_N, "C_Percent" = C_Percent_100_585[,1])
pH_100_585 = as.data.frame(predict(direct3, newdata = predict_100_585_Na_EVI2_N_C, level = 0))
predict_100_585_Na_EVI2_N_C_pH = cbind(predict_100_585_Na_EVI2_N_C , "pH" = pH_100_585[,1])
PGPF_100_585 = as.data.frame(predict(direct2, newdata = predict_100_585_Na_EVI2_N_C_pH, level = 0))
predict_100_585_Na_EVI2_N_C_pH_PGPF = cbind(predict_100_585_Na_EVI2_N_C_pH , "PGBF" = PGPF_100_585[,1])
Shannon_100_585 = as.data.frame(predict(direct, newdata = predict_100_585_Na_EVI2_N_C_pH_PGPF, level = 0))
predict_100_585_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_100_585_Na_EVI2_N_C_pH_PGPF, "Shannon" = Shannon_100_585[,1])
write.csv(predict_100_585_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_fungi_100years_585_v2.csv")

#timne period 2040-2060, SSP585
predictdata_50_585 = data.frame(Mean_Ann_Temp= dataset3$X2050_585_b1, Mean_Ann_Precip = dataset3$X2050_585_b12)

Na_50_585_2 = as.data.frame( predict(direct4,newdata = predictdata_50_585, level = 0))
predict_50_585_Na= cbind(predictdata_50_585, "Na" = Na_50_585_2[,1])
EVI2_50_585_2 = as.data.frame( predict(direct6,newdata = predictdata_50_585, level = 0))
predict_50_585_Na_EVI2 = cbind(predict_50_585_Na, "EVI2" = EVI2_50_585_2[,1])
N_Percent_50_585 = as.data.frame(predict(direct7, newdata = predict_50_585_Na_EVI2, level = 0))
predict_50_585_Na_EVI2_N = cbind(predict_50_585_Na_EVI2, "N_Percent" = N_Percent_50_585[,1])
C_Percent_50_585 = as.data.frame(predict(direct5, newdata = predict_50_585_Na_EVI2_N, level = 0))
predict_50_585_Na_EVI2_N_C = cbind(predict_50_585_Na_EVI2_N, "C_Percent" = C_Percent_50_585[,1])
pH_50_585 = as.data.frame(predict(direct3, newdata = predict_50_585_Na_EVI2_N_C, level = 0))
predict_50_585_Na_EVI2_N_C_pH = cbind(predict_50_585_Na_EVI2_N_C , "pH" = pH_50_585[,1])
PGPF_50_585 = as.data.frame(predict(direct2, newdata = predict_50_585_Na_EVI2_N_C_pH, level = 0))
predict_50_585_Na_EVI2_N_C_pH_PGPF = cbind(predict_50_585_Na_EVI2_N_C_pH , "PGBF" = PGPF_50_585[,1])
Shannon_50_585 = as.data.frame(predict(direct, newdata = predict_50_585_Na_EVI2_N_C_pH_PGPF, level = 0))
predict_50_585_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_50_585_Na_EVI2_N_C_pH_PGPF, "Shannon" = Shannon_50_585[,1])
write.csv(predict_50_585_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_fungi_50years_585_v2.csv")

#timne period 2080-2100, SSP126
predictdata_100_126 = data.frame(Mean_Ann_Temp= dataset3$X2100_126_b1, Mean_Ann_Precip = dataset3$X2100_126_b12)

Na_100_126_2 = as.data.frame( predict(direct4,newdata = predictdata_100_126, level = 0))
predict_100_126_Na= cbind(predictdata_100_126, "Na" = Na_100_126_2[,1])
EVI2_100_126_2 = as.data.frame( predict(direct6,newdata = predictdata_100_126, level = 0))
predict_100_126_Na_EVI2 = cbind(predict_100_126_Na, "EVI2" = EVI2_100_126_2[,1])
N_Percent_100_126 = as.data.frame(predict(direct7, newdata = predict_100_126_Na_EVI2, level = 0))
predict_100_126_Na_EVI2_N = cbind(predict_100_126_Na_EVI2, "N_Percent" = N_Percent_100_126[,1])
C_Percent_100_126 = as.data.frame(predict(direct5, newdata = predict_100_126_Na_EVI2_N, level = 0))
predict_100_126_Na_EVI2_N_C = cbind(predict_100_126_Na_EVI2_N, "C_Percent" = C_Percent_100_126[,1])
pH_100_126 = as.data.frame(predict(direct3, newdata = predict_100_126_Na_EVI2_N_C, level = 0))
predict_100_126_Na_EVI2_N_C_pH = cbind(predict_100_126_Na_EVI2_N_C , "pH" = pH_100_126[,1])
PGPF_100_126 = as.data.frame(predict(direct2, newdata = predict_100_126_Na_EVI2_N_C_pH, level = 0))
predict_100_126_Na_EVI2_N_C_pH_PGPF = cbind(predict_100_126_Na_EVI2_N_C_pH , "PGBF" = PGPF_100_126[,1])
Shannon_100_126 = as.data.frame(predict(direct, newdata = predict_100_126_Na_EVI2_N_C_pH_PGPF, level = 0))
predict_100_126_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_100_126_Na_EVI2_N_C_pH_PGPF, "Shannon" = Shannon_100_126[,1])
write.csv(predict_100_126_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_fungi_100years_126_v2.csv")

#timne period 2040-2060, SSP126
predictdata_50_126 = data.frame(Mean_Ann_Temp= dataset3$X2050_126_b1, Mean_Ann_Precip = dataset3$X2050_126_b12)

Na_50_126_2 = as.data.frame( predict(direct4,newdata = predictdata_50_126, level = 0))
predict_50_126_Na= cbind(predictdata_50_126, "Na" = Na_50_126_2[,1])
EVI2_50_126_2 = as.data.frame( predict(direct6,newdata = predictdata_50_126, level = 0))
predict_50_126_Na_EVI2 = cbind(predict_50_126_Na, "EVI2" = EVI2_50_126_2[,1])
N_Percent_50_126 = as.data.frame(predict(direct7, newdata = predict_50_126_Na_EVI2, level = 0))
predict_50_126_Na_EVI2_N = cbind(predict_50_126_Na_EVI2, "N_Percent" = N_Percent_50_126[,1])
C_Percent_50_126 = as.data.frame(predict(direct5, newdata = predict_50_126_Na_EVI2_N, level = 0))
predict_50_126_Na_EVI2_N_C = cbind(predict_50_126_Na_EVI2_N, "C_Percent" = C_Percent_50_126[,1])
pH_50_126 = as.data.frame(predict(direct3, newdata = predict_50_126_Na_EVI2_N_C, level = 0))
predict_50_126_Na_EVI2_N_C_pH = cbind(predict_50_126_Na_EVI2_N_C , "pH" = pH_50_126[,1])
PGPF_50_126 = as.data.frame(predict(direct2, newdata = predict_50_126_Na_EVI2_N_C_pH, level = 0))
predict_50_126_Na_EVI2_N_C_pH_PGPF = cbind(predict_50_126_Na_EVI2_N_C_pH , "PGBF" = PGPF_50_126[,1])
Shannon_50_126 = as.data.frame(predict(direct, newdata = predict_50_126_Na_EVI2_N_C_pH_PGPF, level = 0))
predict_50_126_Na_EVI2_N_C_pH_Shannon_PGPB = cbind(predict_50_126_Na_EVI2_N_C_pH_PGPF, "Shannon" = Shannon_50_126[,1])
write.csv(predict_50_126_Na_EVI2_N_C_pH_Shannon_PGPB, "table_predicted_fungi_50years_126_v2.csv")