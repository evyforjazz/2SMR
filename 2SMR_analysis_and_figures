#install package
install.packages("devtools")
devtools::install_github("MRCIEU/TwoSampleMR")

#library package
library(TwoSampleMR)

#exposure data
extdata <- read.csv("/Users/wangyiwei/PhD_thesis/MR_gur_skin/UK-twins-after-clumping.csv", header = TRUE)
mwas_file <- format_data(extdata, type="exposure")
mwas_file_cluster <- clump_data(mwas_file, clump_r2 = 0.1, clump_kb = 500, clump_p1 = 0.05)

#outcome data
ao <- available_outcomes()
allergy_ID = c('ukb-b-5911', 'ukb-b-9841', 'ukb-b-4601', 'ukb-b-16702', 'ukb-b-10351', 'ukb-b-18787', 'ukb-b-9039',
              'ukb-b-20296', 'ukb-b-16207', 'ukb-b-17241', 'ukb-a-93' ,'ukb-a-446', 'ukb-a-444', 'ukb-a-447', 
              'finn-a-CHILDHOOD_ALLERGY','finn-b-POLLENALLERGY', 'ieu-a-996') 
allergy_outcome <- extract_outcome_data(snps = mwas_file_cluster$SNP, outcomes = allergy_ID, proxies = FALSE, maf_threshold = 0.01)

#harmonise data
dat_llergy <- harmonise_data(mwas_file_cluster, allergy_outcome, action = 2)

#perform 2smr 
res_allergy <- mr(dat_allergy)

#sensitive analysis
mr_heter_allergy <- mr_heterogeneity(dat_allergy, method_list=c("mr_egger_regression", "mr_ivw"))
mr_ple_allergy <- mr_pleiotropy_test(dat_allergy)
res_single_allergy <- mr_singlesnp(dat_allergy)
res_loo_allergy <- mr_leaveoneout(dat_allergy)

#calculate odds ratio
res_or_allergy <- generate_odds_ratios(res_allergy)
res_single_allergy_or <- generate_odds_ratios(res_single_allergy)

#save results
write.csv(allergy_outcome, file = "/Users/wangyiwei/Desktop/allergy_outcome.csv")
write.csv(res_allergy, file = "/Users/wangyiwei/Desktop/res_allergy.csv")
write.csv(mr_heter_allergy, file = "/Users/wangyiwei/Desktop/me_heter_allergy.csv")
write.csv(mr_ple_allergy, file = "/Users/wangyiwei/Desktop/mr_ple_allergy.csv")
write.csv(res_single_allergy, file = "/Users/wangyiwei/Desktop/res_single_allergy.csv")
write.csv(res_loo_allergy, file = "/Users/wangyiwei/Desktop/res_loo_allergy.csv")
write.csv(res_or_allergy, file = "/Users/wangyiwei/Desktop/res_or_allergy.csv")
write.csv(res_single_allergy_or, file = "/Users/wangyiwei/Desktop/res_single_allergy_or.csv")

#prepare forest plot 

library(forestplot)

#figure1
AD <- read.csv("/Users/wangyiwei/Desktop/AD_single.csv", header = F)
forestplot(as.matrix(AD[,c(1:5)]), mean = AD$V6, lower = AD$V7, upper = AD$V8, 
           zero = 1, 
           xlog = F, 
           fn.ci_norm = fpDrawCircleCI, 
           boxsize = 0.15, 
           col=fpColors(line = "#CC79A7",
                        box="#D55E00"), 
           lty.ci = 7,  
           lwd.ci = 3, 
           ci.vertices.height = 0.05,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), 
           xlab = gpar(cex = 0.2), 
           cex = 1.2), 
           lineheight = "auto", 
           graph.pos = 5, 
           xticks = c(0.98,0.99,1,1.01, 1.02, 1.03)）

#figure3a
single_446 <- read.csv("/Users/wangyiwei/Desktop/single_446.csv", header = F)
forestplot(as.matrix(single_446[,c(1:5)]), mean = single_446$V6, lower = single_446$V7, upper = single_446$V8, 
           zero = 1, 
           xlog = F, 
           fn.ci_norm = fpDrawCircleCI, 
           boxsize = 0.15, 
           col=fpColors(line = "#CC79A7",
                        box="#D55E00"), 
           lty.ci = 7,  
           lwd.ci = 3, 
           ci.vertices.height = 0.05,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), 
           xlab = gpar(cex = 0.2), 
           cex = 1.2), 
           lineheight = "auto", 
           graph.pos = 5, 
           xticks = c(0.98,0.99,1,1.01, 1.02, 1.03)）
           
#figure3b
single_20296 <- read.csv("/Users/wangyiwei/Desktop/single_202962.csv", header = F)
forestplot(as.matrix(single_20296[,c(1:5)]), mean = single_20296$V6, lower = single_20296$V7, upper = single_20296$V8, 
           zero = 1, 
           xlog = F, 
           fn.ci_norm = fpDrawCircleCI, 
           boxsize = 0.15, 
           col=fpColors(line = "#CC79A7",
                        box="#D55E00"), 
           lty.ci = 7,  
           lwd.ci = 3, 
           ci.vertices.height = 0.05,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), 
           xlab = gpar(cex = 0.2), 
           cex = 1.2), 
           lineheight = "auto", 
           graph.pos = 5, 
           xticks = c(0.98,0.99,1,1.01, 1.02, 1.03)）

#figure2
single_allergy <- read.csv("/Users/wangyiwei/Desktop/single_allergy.csv", header = F)
forestplot(as.matrix(single_allergy[,c(1:6)]), mean = single_allergy$V7, lower = single_allergy$V8, upper = single_allergy$V9, 
           zero = 1, 
           xlog = F, 
           fn.ci_norm = fpDrawCircleCI, 
           boxsize = 0.15, 
           col=fpColors(line = "#CC79A7",
                        box="#D55E00"), 
           lty.ci = 7,  
           lwd.ci = 3, 
           ci.vertices.height = 0.05,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.7), 
           xlab = gpar(cex = 0.2), 
           cex = 1.2), 
           lineheight = "auto", graph.pos = 5, 
           xticks = c(0.98,0.99,1,1.01, 1.02, 1.03)）
           
          
