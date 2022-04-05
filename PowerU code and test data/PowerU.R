# PowerU.R is prepared to

# 1. Select features which are suitable for signal calibration
# 2. Perform cross-validation for each selected metabolic feature to choose the best regression model
# 3. Construct calibration curve using serial QC injections under the best regression model
# 4. Convert MS signal intensities in real samples to the corresponding QC relative concentration using the built model
# 5. Perform Shapiro-Wilk normality test on converted data
# 6. Select running t-test or U test according to data normality (designed for two groups)
# 7. Add notations for each metabolic feature

# Created: 2021-03-18
# Created by: Huaxu Yu
# Version: 1.0.0
# Copy right 2020 Huan lab @ University of British Columbia

# Input: QC data file; sample data file in .csv format
# Output: One sample data file in .csv format

##################################################
#Parameter setting

#Import Library for polynomial regression
library(polynom)

# File input
# the data path for the metabolite intensity table
Calibration_datapath = ""
# the file name of peak intensity table of real samples
sample_FileName = "sample_input.csv"
# the file name of peak intensity table of serial QC samples
QC_FileName = "QC_input.csv"

# QC information. Relative loading amount of serial QC samples
QC_conc = c(10,20,30,40,50,60,70,80,90,100)

# Sample sizes in two groups
sample_size = c(6,6)

# Filters for feature selection based on overall linearity
R2_threshold = 0.8              # the feature with R2 value lager than this threshold will be qualified for signal calibration
k_threshold = 0                 # the feature with k value lager than this threshold will be qualified for signal calibration
int_threshold = 0               # an intensity larger than this value will be known as valid intensity
LR_QC_points = 6                # the number of selected QC samples less than this value will only be calibrated using linear model
QR_QC_points = 9                # the number of selected QC samples less than this value but no less than LR_QC_points will be calibrated using linear or quadratic model

# the number of selected QC samples no less than QR_QC_points will be calibrated using linear quadratic, or cubic model

##################################################
#Load files

setwd(Calibration_datapath)
#Load sample file
sample_table = read.csv(sample_FileName, stringsAsFactors = FALSE)
#Load QC file
QC_table = read.csv(QC_FileName, stringsAsFactors = FALSE)

##################################################
#Create function for cross-validation
cross_validation = function(intensity_seq,conc_seq,order_number) {
  comparison_result = c()
  for (p in 2:(length(intensity_seq)-3)) {
    for (q in (p+2):(length(intensity_seq)-1)) {
      real_FC = q/p

      valid_data1 = as.numeric(intensity_seq[p])
      valid_data2 = as.numeric(intensity_seq[q])

      training_intensity_seq = intensity_seq[-c(p,q)]
      training_conc_seq = conc_seq[-c(p,q)]

      #Uncalibrated Ratio
      Uncali_FC = as.numeric(valid_data2/valid_data1)

      #Linear regression
      Linear_valid_data = calibrate_intensity(training_intensity_seq,training_conc_seq,1,c(valid_data1,valid_data2))
      Linear_calibrated_FC = Linear_valid_data[[1]][2]/Linear_valid_data[[1]][1]

      #Quadratic regression
      if(order_number >= 2){
        Quadratic_valid_data = calibrate_intensity(training_intensity_seq,training_conc_seq,2,c(valid_data1,valid_data2))
        Quadratic_Calibrated_FC = Quadratic_valid_data[[1]][2]/Quadratic_valid_data[[1]][1]
      } else{Quadratic_Calibrated_FC = 10000}


      #Cubic regression if order number is 3
      if(order_number >= 3){
        Cubic_valid_data = calibrate_intensity(training_intensity_seq,training_conc_seq,3,c(valid_data1,valid_data2))
        Cubic_Calibrated_FC = Cubic_valid_data[[1]][2]/Cubic_valid_data[[1]][1]
      } else{Cubic_Calibrated_FC = 10000}

      FC_diff = abs(c(Uncali_FC,Linear_calibrated_FC,Quadratic_Calibrated_FC,Cubic_Calibrated_FC)-real_FC)
      comparison_result = c(comparison_result, match(min(FC_diff),FC_diff))
    }
  }
  #Use lower order of regression if two models show same performance in cross-cvalidation
  return(as.numeric(names(sort(table(comparison_result),decreasing=TRUE)[1]))-1)
}

#####################################################
# Create function for MS signal calibration
# Need to input selected QC intensities, selected QC concentrations, regression model, and intensities for calibration
calibrate_intensity = function(s_QC_int, s_QC_conc, model_level, real_intensities){

  if(model_level != 0){
    calibrated_intensity = c()

    Re_model = lm(as.numeric(s_QC_int) ~ poly(s_QC_conc, model_level, raw = T))
    Re_coeff = Re_model$coefficients

    for (int in 1:length(real_intensities)) {
      cali_int = 0
      Re_equation = polynomial(c(Re_coeff[1]-real_intensities[int],Re_coeff[2:length(Re_coeff)]))
      All_solutions = solve(Re_equation)
      All_solutions = All_solutions[which(Im(All_solutions) == 0)]

      pre_cali_int = (Re(All_solutions)[Re(All_solutions) < tail(s_QC_conc,1) & Re(All_solutions) > 0])

      if(length(pre_cali_int) == 0){pre_cali_int = (Re(All_solutions)[Re(All_solutions) < 1.5*tail(s_QC_conc,1)])}
      if(length(pre_cali_int) == 0){pre_cali_int = (Re(All_solutions)[Re(All_solutions) < 2*tail(s_QC_conc,1)])}
      if(length(pre_cali_int) != 0){cali_int = max(pre_cali_int)}

      if(model_level == 1){cali_int = All_solutions}

      if(cali_int<0){cali_int = real_intensities[int]/s_QC_int[1]*10}
      calibrated_intensity = c(calibrated_intensity,cali_int)
    }
    calibrated_intensity = list(calibrated_intensity)} else {calibrated_intensity = list(real_intensities)}

  return(calibrated_intensity)
}

######################################################
#Acquire the number of serial diluted QC samples
QC_points = length(QC_conc)
#Prepare calibrated sample table
calibrated_table = sample_table

calibrated_table$notation = NA
calibrated_table$model = NA
calibrated_table$QC_number = NA
model_pool = c("Uncali.","Linear", "Quadratic", "Cubic")

for (i in 1:nrow(sample_table)) {
  QC_int = QC_table[i,4:(3+QC_points)]
  valid_int = which(QC_int > int_threshold)
  QC_int = QC_int[valid_int]
  selected_QC_conc = QC_conc[valid_int]
  selected_int_points = length(selected_QC_conc)

  #perform linear regression to select the features good for signal calibration
  if (selected_int_points <= 3) {
    calibrated_table$notation[i] = "Insufficient_QC_points"
  } else {
    filter_regression = lm(as.numeric(QC_int) ~ selected_QC_conc)
    slope = as.numeric(filter_regression$coefficients[2])
    cor_coeff = as.numeric(summary(filter_regression)[8])
    if(slope > k_threshold & cor_coeff > R2_threshold & selected_int_points > 3){
      calibrated_table$notation[i] = "Good for calibration"
      if(selected_int_points < LR_QC_points){
        best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,1)
      } else if(selected_int_points>=LR_QC_points & selected_int_points<QR_QC_points){
        best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,2)
      } else {
        best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,3)
      }
      calibrated_int = calibrate_intensity(QC_int,selected_QC_conc,best_model,sample_table[i,4:ncol(sample_table)])[[1]]

      for (j in 4:(ncol(sample_table))) {
        if(sample_table[i,j] != 0 & as.numeric(calibrated_int[j-3]) == 0){
          best_model = cross_validation(as.numeric(QC_int),selected_QC_conc,1)
          calibrated_int = calibrate_intensity(QC_int,selected_QC_conc,best_model,sample_table[i,4:ncol(sample_table)])[[1]]
          break
        }
      }
      calibrated_table[i,4:ncol(sample_table)] = calibrated_int
      calibrated_table$model[i] = model_pool[best_model+1]
      calibrated_table$QC_number[i] = selected_int_points

    } else {
      calibrated_table$notation[i] = "Not suitable for calibration"
    }
  }
}

calibrated_table$p.value = 0
calibrated_table$test = 0

for (i in 1:nrow(calibrated_table)) {
  g1 = as.numeric(calibrated_table[i,4:(3+sample_size[1])])
  g2 = as.numeric(calibrated_table[i,(4+sample_size[1]):(3+sample_size[1]+sample_size[2])])
  if (var(g1) == 0 | var(g2) == 0) {
    calibrated_table$test[i] = "t-test, identical values observed within group."
    calibrated_table$p.value[i] = t.test(g1, g2, conf.level = 0.95, var.equal = F)$p.value
    next
  }
  n1 = shapiro.test(g1)$p.value
  n2 = shapiro.test(g2)$p.value
  f = var.test(g1,g2)$p.value
  if(n1 > 0.1 & n2 > 0.1){
    calibrated_table$test[i] = "t-test"
    if(f > 0.05){
      calibrated_table$p.value[i] = t.test(g1, g2, conf.level = 0.95, var.equal = T)$p.value
    } else{
      calibrated_table$p.value[i] = t.test(g1, g2, conf.level = 0.95, var.equal = F)$p.value
      }
  } else{
    calibrated_table$test[i] = "U-test"
    calibrated_table$p.value[i] = wilcox.test(g1, g2, exact = F)$p.value}
}

write.csv(calibrated_table,"calibrated_sample_table.csv", row.names = F)




