########### Notes and Project Background ####
# The purpose of this modelling is to investigate what is affecting breeding success of the eastern quoll (Dasyurus viverrinus),
# in an effort to understand drivers of reproductive sucess and population-level patterns of genetic diversity. 
# Prior investigation using single-nucleotide polymorphisms (genetic data) has provided us with a pedigree for all known
# quolls in Mulligan's Flat Woodland Sanctuary, from which we have been able to determine some reproductive patterns. 
# What remains uncertain is what makes a 'good' quoll: that is, what traits lead to higher reproductive success/output?

# Generally, we're asking the question "What affects breeding in the population" in a hurdle-model approach.
#   1. Of the morphological and genetic factors possessed by an individual (colour, heterozygosity, etc), which affect breeding success?
#   2. For individuals that breed successfully, do any traits influence the number of offspring produced?

# Note that in this analysis "successful" production of offspring is equivalent to recruiting offspring into the
# next generation, rather than conception/birth being the endpoint. This is because individuals can only be captured and 
# microchipped/sampled if they survive from birth to dispersal from the natal den. 

#### SNP filtering set-up ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  +   install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("qvalue")
library(plotly)
library(SNPRelate)
library(adegenet)
library(hierfstat)
library(dartR)
gl.install.vanilla.dartR() # This installs all missing and required packages of dartR (my version, dartR 2.9.7)
# set working directory # 
dartdata <- "SNP_2_bettongANDduplicates_removed.csv"
ind_metadata <- "idpop_quoll_allinfo.csv"
gl.data <- gl.read.dart(dartdata, ind.metafile=ind_metadata)
names(gl.data$other$loc.metrics) # check data object column names
#### Add SNP quality metrics #####
#We can add a metric of mean depth across both alleles from DArTSeq data (default provided is depth per allele)
depth_bothalleles <- gl.data$other$loc.metrics$AvgCountRef + gl.data$other$loc.metrics$AvgCountSnp
#Also add a ratio of seq depth of ref to SNP allele that I like to use to check for PCR bias towards one allele (may be indicative of null allele propensity at that SNP)
alleleDepthRatio <- gl.data$other$loc.metrics$AvgCountRef/gl.data$other$loc.metrics$AvgCountSnp

#### Plot basic SNP quality metrics (save your plot)####
glPlot(gl.data)
par(mfrow = c(2,2))
hist(gl.data$other$loc.metrics$CallRate, xlab = "Call Rate", main = "Call Rate")
hist(depth_bothalleles, xlab = "Depth", main = "Depth")
hist(alleleDepthRatio, xlab = "Allele Seq Depth Ratio", main = "Allele Seq Depth Ratio")
hist(gl.data$other$loc.metrics$RepAvg, xlab = "Repeatability", main = "Repeatability")

#### Calculate basic pop gen metrics ####
#Overall allele frequencies, assuming of bi-allelic SNPs)
gl.data$other$loc.metrics$alf <- gl.alf(gl.data)[,2]

# Calculate allele frequencies and heterozygosities within sub-pops
allelefreq1 <- aggregate(gl.data, list(gl.data$pop), function(x) {(mean(x, na.rm=TRUE)/2)})
hetcount <- aggregate(gl.data, list(gl.data$pop), function(x) {sum(x[x=="1"], na.rm=TRUE)})
samplsize <- aggregate(gl.data, list(gl.data$pop), function (x) {length(x[!is.na(x)])})
Hs <- 1-(allelefreq1[,-1]^2) - ((1-allelefreq1[,-1])^2)
Ho <- hetcount[,-1]/samplsize[,-1]
Fis <- (Hs-Ho)/Hs

#Append means across sub-pops to genlight object
gl.data$other$loc.metrics$Ht <- 2*gl.data$other$loc.metrics$alf*(1-gl.data$other$loc.metrics$alf)
gl.data$other$loc.metrics$meanHs <- colMeans(Hs, na.rm=TRUE)
gl.data$other$loc.metrics$meanHo <- colMeans(Ho, na.rm=TRUE)
gl.data$other$loc.metrics$meanFis <- colMeans(Fis, na.rm=TRUE)
gl.data$other$loc.metrics$mean_p <- colMeans(allelefreq1[-1], na.rm=TRUE)
gl.data$other$loc.metrics$Ho <- colSums(hetcount[,-1],na.rm=TRUE)/colSums(samplsize[,-1])
gl.data$other$loc.metrics$Fit <-(gl.data$other$loc.metrics$Ht- gl.data$other$loc.metrics$Ho)/gl.data$other$loc.metrics$Ht

#Estimate FST as (Ht-Hs)/Ht. 
gl.data$other$loc.metrics$Fst <-(gl.data$other$loc.metrics$Ht-gl.data$other$loc.metrics$meanHs)/gl.data$other$loc.metrics$Ht

#### Plot pop gen stats ####
par(mfrow = c(2,4))
hist(gl.data$other$loc.metrics$mean_p, xlab = "Allele Freq (mean within pops)", main = NULL, ylab = NULL)
hist(gl.data$other$loc.metrics$meanHo, xlab = "Ho (mean within pops)", main = NULL, ylab = NULL)
hist(gl.data$other$loc.metrics$meanHs, xlab = "Mean Hs", main = NULL, ylab = NULL)
hist(gl.data$other$loc.metrics$Ht, xlab = "Ht", main = NULL, ylab = NULL)
plot(gl.data$other$loc.metrics$meanHs, gl.data$other$loc.metrics$meanHo, cex = 0.5, col = "blue", xlab = "Mean Hs", ylab = "Ho")
hist(gl.data$other$loc.metrics$meanFis, xlab = "Fis", main = NULL, ylab = NULL)
hist(gl.data$other$loc.metrics$Fit, xlab = "Fit", main = NULL, ylab = NULL)
hist(gl.data$other$loc.metrics$Fst, xlab = "Fst", main = NULL, ylab = NULL)

#### Plot all these pop gen metrics against SNP quality####
par(mfrow = c(2,3))
plot(gl.data$other$loc.metrics$CallRate,gl.data$other$loc.metrics$meanFis, cex = 0.4, col = "green", xlab = "Call Rate", ylab = "Fis")
plot(gl.data$other$loc.metrics$RepAvg,gl.data$other$loc.metrics$meanFis, cex = 0.4, col = "pink", xlab = "Repeatability", ylab = "Fis")
plot(depth_bothalleles,gl.data$other$loc.metrics$meanFis, cex = 0.4, col = "blue", xlab = "Depth", ylab = "Fis")
plot(depth_bothalleles,gl.data$other$loc.metrics$meanHo, cex = 0.4, col = "orange", xlab = "Depth", ylab = "Ho")
plot(gl.data$other$loc.metrics$CallRate,gl.data$other$loc.metrics$meanHo, cex = 0.4, col = "red", xlab = "Call Rate", ylab = "Ho")
plot(depth_bothalleles,gl.data$other$loc.metrics$CallRate, cex = 0.4, col = "purple", xlab = "Depth", ylab = "Call Rate")

#### PCoA before filtering ####

pcoa <- gl.pcoa(gl.data)
gl.pcoa.scree(pcoa)
gl.pcoa.plot(pcoa, gl.data)

#3d plot is nice too (shape can be cube or sphere, tetrahaedron doesn't work)
gl.pcoa.plot.3d(pcoa, gl.data, title = "PCoA", xaxis = 1, yaxis = 2, zaxis = 3, shape = "sphere", radius = 2, legend = "topright")

#This probably won't change much with SNP filtering as the broad pattrns of structure are pretty stable, but some if the individual 
#and pop-level diversity metrics are sensitive to genotype quality then it may

#### filtering based on loci ####
#Call Rate threshold (proportion of samples typed per SNP) - base on frequency rug plot from earlier
gl.data.filt <- gl.filter.callrate(gl.data, method = "loc", threshold = 0.95, mono.rm= TRUE, plot.out = TRUE)
#Repeatability threshold (again have a look at earlier plots)
gl.data.filt <- gl.filter.reproducibility(gl.data.filt, threshold = 0.95, v = 3)
# Test affect of multiple MAFs on population structure, as it can have noticeable affect (see Linck and Battey, 2019' Brockett et al 2022)
gl.data.filt <- gl.filter.maf(gl.data.filt, threshold = 0.05, v = 3)
gl.filter.monomorphs(gl.data.filt) # no monomorphs to remove 

#### filtering based on individuals ####
ind_callrate <- gl.report.callrate(gl.data.filt, method = "ind", v = 3)
gl.data.filt <- gl.filter.callrate(gl.data.filt, method = "ind", threshold = 0.90, mono.rm = TRUE, recalc = FALSE, plot.out = TRUE, v = 3)
gl.filter.monomorphs(gl.data.filt)
#### save filtered DartR object ####
save(gl.data.filt, file = "gl_maf05_allinfo.R")

#### Check rug plot ####
# plot rug plot before and after, together
par(mfrow = c(2,1))
glPlot(gl.data)
glPlot(gl.data.filt)
#### Plot pop gen of filtered dataset####

#First recalculate the pop gen metrics if you have removed any individuals
gl.data.filt$other$loc.metrics$alf <- gl.alf(gl.data.filt)[,2]
allelefreq1 <- aggregate(gl.data.filt, list(gl.data.filt$pop), function(x) {(mean(x, na.rm=TRUE)/2)})
hetcount <- aggregate(gl.data.filt, list(gl.data.filt$pop), function(x) {sum(x[x=="1"], na.rm=TRUE)})
samplsize <- aggregate(gl.data.filt, list(gl.data.filt$pop), function (x) {length(x[!is.na(x)])})
Hs <- 1-(allelefreq1[,-1]^2) - ((1-allelefreq1[,-1])^2)
Ho <- hetcount[,-1]/samplsize[,-1]
Fis <- (Hs-Ho)/Hs
gl.data.filt$other$loc.metrics$Ht <- 2*gl.data.filt$other$loc.metrics$alf*(1-gl.data.filt$other$loc.metrics$alf)
gl.data.filt$other$loc.metrics$meanHs <- colMeans(Hs, na.rm=TRUE)
gl.data.filt$other$loc.metrics$meanHo <- colMeans(Ho, na.rm=TRUE)
gl.data.filt$other$loc.metrics$meanFis <- colMeans(Fis, na.rm=TRUE)
gl.data.filt$other$loc.metrics$Ho <- colSums(hetcount[,-1],na.rm=TRUE)/colSums(samplsize[,-1])
gl.data.filt$other$loc.metrics$Fit <-(gl.data.filt$other$loc.metrics$Ht-gl.data.filt$other$loc.metrics$Ho)/gl.data.filt$other$loc.metrics$Ht
gl.data.filt$other$loc.metrics$mean_p <- colMeans(allelefreq1[-1], na.rm=TRUE)
gl.data.filt$other$loc.metrics$Fst <-(gl.data.filt$other$loc.metrics$Ht - gl.data.filt$other$loc.metrics$meanHs)/gl.data.filt$other$loc.metrics$Ht

depth_bothalleles <- gl.data.filt$other$loc.metrics$AvgCountRef + gl.data.filt$other$loc.metrics$AvgCountSnp
alleleDepthRatio <- gl.data.filt$other$loc.metrics$AvgCountRef/gl.data.filt$other$loc.metrics$AvgCountSnp

#Second, re plot for post-filt dataset
par(mfrow = c(2,4))
hist(gl.data.filt$other$loc.metrics$mean_p, xlab = "Allele Freq - mean within pops", main = NULL, ylab = NULL)
hist(gl.data.filt$other$loc.metrics$meanHo, xlab = "Mean Ho", main = NULL, ylab = NULL)
hist(gl.data.filt$other$loc.metrics$meanHs, xlab = "Mean Hs", main = NULL, ylab = NULL)
hist(gl.data.filt$other$loc.metrics$Ht, xlab = "Ht", main = NULL, ylab = NULL)
plot(gl.data.filt$other$loc.metrics$meanHs, gl.data.filt$other$loc.metrics$Ho, cex = 0.5, col = "blue", xlab = "Mean Hs", ylab = "Ho")
hist(gl.data.filt$other$loc.metrics$meanFis, xlab = "Fis", main = NULL, ylab = NULL)
hist(gl.data.filt$other$loc.metrics$Fit, xlab = "Fit", main = NULL, ylab = NULL)
hist(gl.data.filt$other$loc.metrics$Fst, xlab = "Fst", main = NULL, ylab = NULL)

#Plot basic SNP metrics
par(mfrow = c(2,2))
hist(gl.data.filt$other$loc.metrics$CallRate, xlab = "Call Rate", main = NULL)
hist(depth_bothalleles, xlab = "Depth", main = NULL)
hist(alleleDepthRatio, xlab = "Allele Seq Depth Ratio", main = NULL)
hist(gl.data.filt$other$loc.metrics$RepAvg, xlab = "Repeatability", main = NULL)

#### Summaries of basic genetic diversity before and after filtering #####
postfilt_Fst <- mean(gl.data.filt$other$loc.metrics$Fst, na.rm=TRUE)
postfilt_Ho <- mean(gl.data.filt$other$loc.metrics$meanHo, na.rm=TRUE)
postfilt_Hs <- mean(gl.data.filt$other$loc.metrics$meanHs, na.rm=TRUE)
postfilt_Ht <- mean(gl.data.filt$other$loc.metrics$Ht, na.rm=TRUE)
postfilt_Fis <- mean(gl.data.filt$other$loc.metrics$meanFis, na.rm=TRUE)
popgen_postfilt <- c(postfilt_Hs,postfilt_Ho,postfilt_Ht,postfilt_Fis,postfilt_Fst) 

#What changed compared to pre-filtering?
prefilt_Fst <- mean(gl.data$other$loc.metrics$Fst, na.rm=TRUE)
prefilt_Ho <- mean(gl.data$other$loc.metrics$meanHo, na.rm=TRUE)
prefilt_Hs <- mean(gl.data$other$loc.metrics$meanHs, na.rm=TRUE)
prefilt_Ht <- mean(gl.data$other$loc.metrics$Ht, na.rm=TRUE)
prefilt_Fis <- mean(gl.data$other$loc.metrics$meanFis, na.rm=TRUE)
popgen_prefilt <- c(prefilt_Hs,prefilt_Ho,prefilt_Ht,prefilt_Fis,prefilt_Fst) 

popgen_summary <- rbind(popgen_prefilt, popgen_postfilt)
colnames(popgen_summary) <- c("Hs", "Ho", "Ht", "Fis", "Fst")
popgen_summary
write.csv(popgen_summary, file="popgen_summary_prepost_maf05.csv")


########################################################################
#### Tests of reproductive change over time ####
library(trend)      # For Mann-Kendall test
library(lmtest)     # For regression tests
library(ggplot2)    # For plotting
library(tidyr)
library(dplyr)      # For data manipulation
library(RColorBrewer)

repro_data <- read.csv("Parentage_over_time.csv")

# Mann-Kendall tests for multivariate change over time ####
mk_offspring <- mk.test(repro_data$offspring)
mk_pct_females <- mk.test(repro_data$pct_females_mothers)
mk_offspring_per_mother <- mk.test(repro_data$avg_offspring_per_mother)
mk_pct_males <- mk.test(repro_data$pct_males_fathers)
mk_offspring_per_father <- mk.test(repro_data$avg_offspring_per_father)

cat("=== MANN-KENDALL TREND TESTS ===\n")
mk_offspring
mk_pct_females
mk_offspring_per_mother
mk_pct_males
mk_offspring_per_father

# linear regression for change over time ####

run_regression <- function(y_var, metric_name) {
  model <- lm(y_var ~ year, data = repro_data)
  model_summary <- summary(model)
  
  cat("\nLinear Regression:", metric_name, "vs Year\n")
  cat("========================\n")
  cat("Slope:", model_summary$coefficients[2,1], "per year\n")
  cat("p-value:", model_summary$coefficients[2,4], "\n")
  cat("R-squared:", model_summary$r.squared, "\n")
  cat("Interpretation:", 
      ifelse(model_summary$coefficients[2,4] < 0.05, 
             "Significant trend", 
             ifelse(model_summary$coefficients[2,4] < 0.1, 
                    "Marginally significant trend", 
                    "No significant trend")),
      "\n\n")
  
  return(model)
}
cat("\n=== LINEAR REGRESSION ANALYSIS ===\n")
reg_offspring <- run_regression(repro_data$offspring, "Total Offspring")
reg_pct_females <- run_regression(repro_data$pct_females_mothers, "% Females that are Mothers")
reg_offspring_per_mother <- run_regression(repro_data$avg_offspring_per_mother, "Avg Offspring per Mother")
reg_pct_males <- run_regression(repro_data$pct_males_fathers, "% Males that are Fathers")
reg_offspring_per_father <- run_regression(repro_data$avg_offspring_per_father, "Avg Offspring per Father")

cat("\n=== SPEARMAN RANK CORRELATION ===\n")
sp_offspring <- run_spearman(repro_data$offspring, "Total Offspring")
sp_pct_females <- run_spearman(repro_data$pct_females_mothers, "% Females that are Mothers")
sp_offspring_per_mother <- run_spearman(repro_data$avg_offspring_per_mother, "Avg Offspring per Mother")
sp_pct_males <- run_spearman(repro_data$pct_males_fathers, "% Males that are Fathers")
sp_offspring_per_father <- run_spearman(repro_data$avg_offspring_per_father, "Avg Offspring per Father")

# Chi-square test for consecutive year-year comparisons ####
run_chi_square <- function(success1, total1, success2, total2, year1, year2, metric_name) {
  # Create contingency table
  contingency <- matrix(c(success1, total1 - success1, success2, total2 - success2), 
                        nrow = 2, byrow = TRUE)
  
  # Run chi-square test
  chi_result <- chisq.test(contingency, correct = TRUE)
  
  # Calculate proportions
  prop1 <- success1 / total1
  prop2 <- success2 / total2
  
  cat("\nChi-square Test:", metric_name, "from", year1, "to", year2, "\n")
  cat("========================\n")
  cat(year1, ":", round(prop1 * 100, 1), "% (", success1, "/", total1, ")\n", sep = "")
  cat(year2, ":", round(prop2 * 100, 1), "% (", success2, "/", total2, ")\n", sep = "")
  cat("Change:", round((prop2 - prop1) * 100, 1), "percentage points\n", sep = "")
  cat("Chi-square:", round(chi_result$statistic, 4), "\n")
  cat("p-value:", round(chi_result$p.value, 4), "\n")
  cat("Interpretation:", 
      ifelse(chi_result$p.value < 0.05, 
             "Significant difference", 
             "No significant difference"),
      "\n\n")
  
  return(chi_result)
}

# Female reproductive participation
cat("\n--- Percentage of Females that are Mothers ---\n")
female_2016_2017 <- run_chi_square(repro_data$mothers[1], repro_data$adult_females[1], 
                                   repro_data$mothers[2], repro_data$adult_females[2],
                                   2016, 2017, "% Females as Mothers")
female_2017_2018 <- run_chi_square(repro_data$mothers[2], repro_data$adult_females[2], 
                                   repro_data$mothers[3], repro_data$adult_females[3],
                                   2017, 2018, "% Females as Mothers")
female_2018_2019 <- run_chi_square(repro_data$mothers[3], repro_data$adult_females[3], 
                                   repro_data$mothers[4], repro_data$adult_females[4],
                                   2018, 2019, "% Females as Mothers")
female_2019_2020 <- run_chi_square(repro_data$mothers[4], repro_data$adult_females[4], 
                                   repro_data$mothers[5], repro_data$adult_females[5],
                                   2019, 2020, "% Females as Mothers")

# Male reproductive participation
cat("\n--- Percentage of Males that are Fathers ---\n")
male_2016_2017 <- run_chi_square(repro_data$fathers[1], repro_data$adult_males[1], 
                                 repro_data$fathers[2], repro_data$adult_males[2],
                                 2016, 2017, "% Males as Fathers")
male_2017_2018 <- run_chi_square(repro_data$fathers[2], repro_data$adult_males[2], 
                                 repro_data$fathers[3], repro_data$adult_males[3],
                                 2017, 2018, "% Males as Fathers")
male_2018_2019 <- run_chi_square(repro_data$fathers[3], repro_data$adult_males[3], 
                                 repro_data$fathers[4], repro_data$adult_males[4],
                                 2018, 2019, "% Males as Fathers")
male_2019_2020 <- run_chi_square(repro_data$fathers[4], repro_data$adult_males[4], 
                                 repro_data$fathers[5], repro_data$adult_males[5],
                                 2019, 2020, "% Males as Fathers")





# Test plots ####
ggplot(repro_data, aes(x = year)) +
  geom_line(aes(y = pct_females_mothers * 100, color = "Females"), size = 1.2) +
  geom_point(aes(y = pct_females_mothers * 100, color = "Females"), size = 2) +
  geom_line(aes(y = pct_males_fathers * 100, color = "Males"), size = 1.2) +
  geom_point(aes(y = pct_males_fathers * 100, color = "Males"), size = 2) +
  labs(title = "Reproductive Participation by Sex",
       x = "Year",
       y = "Percentage of adult population of given sex",
       color = "Sex") +
  scale_color_manual(values = c("Females" = "plum2", "Males" = "skyblue3")) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Per Capita Reproductive Success
ggplot(repro_data, aes(x = year)) +
  geom_line(aes(y = avg_offspring_per_mother, color = "Per Mother"), size = 1.2) +
  geom_point(aes(y = avg_offspring_per_mother, color = "Per Mother"), size = 3) +
  geom_line(aes(y = avg_offspring_per_father, color = "Per Father"), size = 1.2) +
  geom_point(aes(y = avg_offspring_per_father, color = "Per Father"), size = 3) +
  labs(title = "Per Capita Reproductive Success",
       subtitle = "Average Offspring per Reproductive Individual",
       x = "Year",
       y = "Average Offspring",
       color = "Metric") +
  scale_color_manual(values = c("Per Mother" = "deeppink3", "Per Father" = "darkblue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

repro_data$pct_females_mothers_pct <- repro_data$pct_females_mothers * 100
repro_data$pct_males_fathers_pct <- repro_data$pct_males_fathers * 100

par(mar = c(5, 5, 4, 5) + 0.1)  # Adjust margins for labels

# Total Offspring Production
ggplot(repro_data, aes(x = factor(year), y = offspring)) +
  geom_bar(stat = "identity", fill = "darkolivegreen4", alpha = 0.7) +
  geom_text(aes(label = offspring), vjust = -0.5) +
  labs(title = "Total Offspring Production by Year",
       x = "Year",
       y = "Number of Offspring") +
  theme_minimal()

# Dual-axis plot ####

# Create long format data for breeding success
breeding_data <- repro_data %>%
  select(year, pct_females_mothers, pct_males_fathers) %>%
  pivot_longer(cols = c(pct_females_mothers, pct_males_fathers),
               names_to = "sex_metric",
               values_to = "breeding_percentage") %>%
  mutate(sex = ifelse(grepl("females", sex_metric), "Female", "Male"),
         metric_type = "Breeding Success")

# Create long format data for population sizes
population_data <- repro_data %>%
  select(year, adult_females, adult_males) %>%
  pivot_longer(cols = c(adult_females, adult_males),
               names_to = "sex_pop",
               values_to = "population_size") %>%
  mutate(sex = ifelse(grepl("females", sex_pop), "Female", "Male"),
         metric_type = "Population Size")

# Corrected ggplot code for percentage data as whole numbers
trends <- ggplot(repro_data) +
  geom_line(aes(x=year, y=pct_females_mothers, 
                color="% females adults"), size=1.2) +
  geom_point(aes(x=year, y=pct_females_mothers, 
                 color="% females adults"), size=2) +
  geom_line(aes(x=year, y=pct_males_fathers, 
                color="% males adults"), size=1.2) +
  geom_point(aes(x=year, y=pct_males_fathers, 
                 color="% males adults"), size=2) +
  geom_line(aes(x=year, y=adult_females*2, 
                color="Number of adult females"), linetype="dashed", size=1.2) +
  geom_point(aes(x=year, y=adult_females*2, 
                 color="Number of adult females"),  shape=15, size=3) +
  geom_line(aes(x=year, y=adult_males*2, 
                color="Number of adult males"), linetype="dashed", size=1.2) +
  geom_point(aes(x=year, y=adult_males*2, 
                 color="Number of adult males"), shape=18, size=3) +
  # Scales and labels
  scale_y_continuous(
    name="Reproductive Success Rate (%)",
    limits=c(0, 120),
    sec.axis=sec_axis(~./2, name="Adult Population Size")) +
  scale_color_manual(
    name="Metric",
    values=c("% females adults"="purple4", 
             "% males adults"="deepskyblue1",
             "Number of adult females"="maroon", 
             "Number of adult males"="darkslategray3")) +
  scale_x_continuous(breaks = repro_data$year) +
  labs(title="Reproductive Success and Population Size by Year", 
       x="Year") +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.25),
    axis.title.y = element_text(color = "black", size = 12),
    axis.title.y.right = element_text(color = "black", size = 12),
    axis.title.x = element_text(size = 12),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.title = element_blank())


trends

ggsave("Reproductive success trends per year.png",trends)

########################################################################
#### Modelling set-up ####
library(tidyverse)    # For data manipulation and visualization
library(glmmTMB)      # For GLMM modelling
library(DHARMa)       # For model diagnostics
library(MuMIn)        # For model selection
library(performance)  # for running model performance checks
library(bbmle)        # AICtab 
library(AICcmodavg)   # AIC with averaged model
library(dplyr)
library(car)          # For VIF calculation
library(data.table)   # tabulating model outputs
library(effects)
library(webshot)
library(ggplot2)      # Plotting
library(ggeffects)    
library(sjPlot)       
library(patchwork)
library(grid)
library(gridExtra)

set.seed(1234)
dat <-read.csv(file="Parentage_data_formattedforglmm.csv")
head(dat)
summary(dat)
dat[dat == "NA"] <- NA
# Set categorical variables as factors
dat$sex <- as.factor(dat$sex)
dat$morph <- as.factor(dat$morph)
dat$cohort <- as.factor(dat$cohort)
dat$cluster <- as.factor(dat$cluster)
dat$year <- as.factor(dat$year)
dat$id <- as.factor(dat$id)
dat$breedingsuccess <- as.numeric(dat$breedingsuccess)
dat$offspring <- as.numeric(dat$offspring)
dat$age <- as.numeric(dat$age)
dat$weight <- as.numeric(dat$weight)
dat$indhet <- as.numeric(dat$indhet)

# Scale continuous variables so they're comparable to each other 
# Note that for plotting you will need to make sure that actual numbers are used not the scaled versions
dat$age_scaled<-scale(dat$age) 
dat$weight_scaled<-scale(dat$weight)
dat$indhet_scaled<-scale(dat$indhet)
dat$parentgensim_scaled <- scale(dat$parentgensim)

#look at adjusted dataset
head(dat)

# Create a dataset with most NAs removed, as not all information is always collected every year for every individual 
nona <- subset(dat, weight!="NA")
nona <- subset(nona, morph!="NA")
nona <- subset(nona, age!="NA")
nona <- subset(nona, year!="NA")
nona <- subset(nona, sex!="NA")
nona <- subset(nona, parentgensim!="NA")
nona <- subset(nona, indhet!="NA")
nona <- subset(nona, cohort!="NA")
nona <- subset(nona, cluster!="NA")

# after some preliminary investigation, seems as though individual 88C7E may be an outlier, with reproductive success at 200g.
# this could be human error, inaccurate parentage assignment, or incorrect field measurements being recorded. 
pruned <- subset(nona, id!="88C7E") # potential outlier and erroneous point.

# Take a subset of the no-NA dataset containing either males or females, for within-sex tests
males <- subset(pruned, sex!="F")
females <- subset(pruned, sex!="M")

# Subset the dataset into individuals that successfully produced offspring, 
# to reduce the impact of unsuccessful individuals for tests on number of offspring
success <- subset(pruned, breedingsuccess==1)

#### Basic summary statistics ####
breeding_summary <- pruned %>%
  group_by(sex) %>%
  summarize(
    n = n(),
    success_rate = mean(breedingsuccess, na.rm = TRUE),
    avg_offspring = mean(offspring, na.rm = TRUE)
  )
print(breeding_summary)

sexbreed_plot <- ggplot(breeding_summary, aes(x = sex, y = success_rate, fill = sex)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(x = "Sex", y = "Breeding Success Rate", title = "Breeding Success by Sex") +
  ylim(0, 0.7) +
  theme_minimal() +
  theme(legend.position = "none")

sexbreed_plot

breeding_by_age <- nona %>%
  group_by(age) %>%
  summarize(
    n = n(),
    success_rate = mean(breedingsuccess, na.rm = TRUE)
  )
print(breeding_by_age)

bredage_plot <- ggplot(breeding_by_age, aes(x = age, y = success_rate)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(x = "Age", y = "Breeding Success Rate", title = "Breeding Success by Age") +
  ylim(0, 0.7) +
  theme_minimal()

bredage_plot

breeding_by_age_sex <- nona %>%
  group_by(age, sex) %>%
  summarize(
    n = n(),
    success_rate = mean(breedingsuccess, na.rm = TRUE)
  )
print(breeding_by_age_sex)

agesex_plot <- ggplot(breeding_by_age_sex, aes(x = age, y = success_rate, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(x = "Age", y = "Breeding Success Rate", title = "Breeding Success by Age and Sex") +
  ylim(0, 0.7) +
  theme_minimal()
agesex_plot

breeding_by_morph <- nona %>%
  filter(!is.na(morph)) %>%
  group_by(morph) %>%
  summarize(
    n = n(),
    success_rate = mean(breedingsuccess, na.rm = TRUE)
  )
print(breeding_by_morph)

morph_plot <- ggplot(breeding_by_morph, aes(x = morph, y = success_rate, fill = morph)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(x = "Color Morph", y = "Breeding Success Rate", title = "Breeding Success by Morph") +
  theme_minimal() +
  theme(legend.position = "none")
morph_plot
#### Investigate distribution of the data ########
# weight ####

# violin plot to see whether success are different to failures, in comparison to weight
success_violin <- ggplot(pruned, aes(x = factor(breedingsuccess), y = weight)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1)+
  xlab("Breeding Success")+
  ylab("Weight (kg)")+
  theme_bw()
success_violin


success_by_weight<- ggplot(data=pruned)+
  aes(x=weight,y=breedingsuccess)+
  geom_point(alpha = 0.5)+
  geom_smooth(method = glm, 
              method.args = list(family=binomial),
              se=TRUE, 
              colour="black", 
              fill="forestgreen")+
  xlab("Weight (kg)")+
  ylab("")
success_by_weight

offspring_by_weight<- ggplot(data=pruned)+
  aes(x=weight,y=offspring)+
  geom_smooth(method = glm, 
              method.args = list(family=poisson),
              se=TRUE, 
              colour="black", 
              fill="forestgreen")+
  xlab("Weight (kg)")+
  ylab("")
offspring_by_weight

# sex ####
offspring_by_sex<- ggplot(nona,aes(x=sex,y=offspring))+
  geom_boxplot(outlier.fill="black",
               outlier.shape=1)+
  stat_summary(fun=mean)+
  xlab("Sex")+
  ylab("Number of Offspring")
offspring_by_sex


# morph ####
offspring_by_morph<- ggplot(nona,aes(x=morph,y=offspring))+
  geom_boxplot(outlier.fill="black",
               outlier.shape=1)+
  stat_summary(fun=mean)+
  xlab("Morph")+
  ylab("")
offspring_by_morph

success_by_morph<- ggplot(nona,aes(x=morph,y=breedingsuccess))+
  geom_violin(alpha = 0.7) +
  geom_boxplot(width=0.05)+
  xlab("Morph")+
  ylab("")
success_by_morph

empirical_rates <- pruned %>%
  group_by(morph) %>%
  summarize(
    success_rate = mean(breedingsuccess, na.rm = TRUE),
    sample_size = n(),
    se = sqrt(success_rate * (1 - success_rate) / sample_size),
    lower = pmax(0, success_rate - 1.96 * se),
    upper = pmin(1, success_rate + 1.96 * se))

# Create empirical plot
empirical_plot <- ggplot(empirical_rates, aes(x = morph, y = success_rate, fill = morph)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_text(aes(label = paste("n =", sample_size)), vjust = -0.8) +
  scale_fill_manual(values = c("darkgreen", "tan4")) +
  labs(x = "Morph", 
       y = "Observed Breeding Success Rate", 
       title = "Observed Breeding Success Rate by Morph") +
  ylim(0, max(empirical_rates$upper) * 1.1) +  # Add space for sample size text
  theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size = 14, colour = "black"))

# age ####
success_by_age_violin <- ggplot(nona, aes(x = factor(breedingsuccess), y = age)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1)+
  xlab("Breeding Success")+
  ylab("Age (years)")+
  theme_bw()
success_by_age_violin


success_by_age<- ggplot(data=nona)+
  aes(x=age,y=breedingsuccess)+
  geom_smooth(method = glm,
              method.args =list(family="binomial"),
              se=TRUE,
              colour="black", 
              fill="darkgoldenrod1")+
  xlab("Age (years)")+
  ylab("Likelihood of Breeding Success")
success_by_age

offspring_by_age<- ggplot(data=nona)+
  aes(x=age,y=offspring)+
  geom_smooth(method = glm,
              method.args =list(family="poisson"),
              se=TRUE,
              colour="black", 
              fill="darkgoldenrod1")+
  xlab("Age (years)")+
  ylab("Predicted number of offspring")
offspring_by_age

success_by_age_pruned<- ggplot(data=pruned)+
  aes(x=age,y=breedingsuccess)+
  geom_smooth(method = glm,
              method.args =list(family="binomial"),
              se=TRUE,
              colour="black", 
              fill="darkgoldenrod1")+
  xlab("Age (years)")+
  ylab("Likelihood of Breeding Success")
success_by_age_pruned

offspring_by_age_pruned<- ggplot(data=pruned)+
  aes(x=age,y=offspring)+
  geom_smooth(method = glm,
              method.args =list(family="poisson"),
              se=TRUE,
              colour="black", 
              fill="darkgoldenrod1")+
  xlab("Age (years)")+
  ylab("")
offspring_by_age_pruned


# individual heterozygosity ####
success_by_het_violin <- ggplot(nona, aes(x = factor(breedingsuccess), y = indhet)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1)+
  xlab("Breeding Success")+
  ylab("Individual Heterozygosity")+
  theme_bw()
success_by_het_violin

success_by_het<- ggplot(data=nona)+
  aes(x=indhet,y=breedingsuccess)+
  geom_smooth(method = glm,
              method.args=list("binomial"),
              se=TRUE,
              colour="black", 
              fill="firebrick3")+
  xlab("Individual Heterozygosity")+
  ylab("")
success_by_het

offspring_by_het<- ggplot(data=nona)+
  aes(x=indhet,y=offspring)+
  geom_smooth(method = glm,
              method.args=list("poisson"),
              se=TRUE,
              colour="black", 
              fill="firebrick4")+
  xlab("Individual Heterozygosity")+
  ylab("")
offspring_by_het

#### plotting trends ####

trends <- (success_by_weight|offspring_by_weight)/
  (success_by_age|offspring_by_age)/
  (success_by_morph|offspring_by_morph)/
  (success_by_het|offspring_by_het) & theme_classic()
trends

#### Modelling predictors of breeding success ####
binmod_pruned<-glmmTMB(breedingsuccess ~ (scale(weight)+I((scale(weight))^2))*scale(age)+scale(indhet)+sex+morph+
                         (1|id)+
                         (1|year),
                       family="binomial", 
                       data=pruned, 
                       na.action=na.fail) 


# check model fit in GLOBAL MODEL (as per AICmodav package instructions by Mazerolle, 2020)
check_overdispersion(binmod_pruned)
check_collinearity(binmod_pruned)
ch_bin_nona <- check_model(binmod_pruned)
ch_bin_nona

summary(binmod_pruned)
plot_model(binmod_pruned, type="pred",
           terms="weight[all]",
           ci.lvl = 0.95)

# test whether a subset of factors make a better model
modsel_pruned<-dredge(binmod_pruned, trace=2)
modsel_sub_pruned<-subset(modsel_pruned, delta<5)
modsel_sub_pruned

# top model has lowest degrees of freedom, and later models do not improve 
# fit much. Therefore, select top model rather than averaging. 
topmod_pruned<-get.models(modsel_sub_pruned, subset = 1)[[1]]
topmod_pruned
summary(topmod_pruned)

check_overdispersion(topmod_pruned)
check_collinearity(topmod_pruned)
ch_bin_topmod_pruned <- check_model(topmod_pruned)

# summary shows that weight has the greatest effect on breeding success (increasing 
# weight linked to increased success), and that age and fawn morph are both negatively 
# associated with breeding success. Slight negative impact of weight quadratic likely 
# reflective of age overriding positive effects of weight 
# (noting that individuals get heavier as they get older)
plot_model(topmod_pruned)

# model effects table ####
modeltabletest <- tab_model(topmod_pruned)
model_table <- tab_model(topmod_pruned,
                         transform = "exp",
                         show.est = TRUE,           
                         show.stat = TRUE,          
                         show.p = TRUE,             
                         show.ci = 0.95,            
                         show.re.var = TRUE,        
                         show.aic = FALSE,           
                         show.r2 = TRUE,            
                         show.icc = FALSE,
                         collapse.ci = TRUE,
                         title = "Factors Affecting Breeding Success",
                         dv.labels = "Predicted impact", 
                         pred.labels = c(           
                           "(Intercept)" = "Intercept",
                           "scale(weight)" = "Weight (kg)",
                           "I((scale(weight))^2)" = "Weight (kg) quadratic",
                           "scale(age)" = "Age (years)",
                           "morphFawn" = "Morph (fawn)"),
                         string.est = "Odds Ratios (95% CI)",
                         string.stat = "z-statistic",
                         string.p = "p-value",
                         file = "breeding_success_model.html")

model_table
webshot("breeding_success_model.html", 
        "breeding_success_model.png",
        vwidth = 800, vheight = 600)

#### Plotting predictors: breeding success ####
summary(topmod_pruned)
topmod_pruned<-glmmTMB(breedingsuccess ~ scale(weight)+I((scale(weight))^2)+scale(age)+morph+(1|id),
                       family="binomial", data=pruned, na.action=na.fail) 
# weight ####
pred_data <- data.frame(
  weight = seq(min(pruned$weight, na.rm=TRUE), max(pruned$weight, na.rm=TRUE), length.out=100),
  age = mean(pruned$age, na.rm=TRUE),
  morph = levels(pruned$morph)[1],
  id = NA)

preds <- predict(topmod_pruned, newdata=pred_data, se.fit=TRUE, re.form=NA)

# Create a data frame with predictions and confidence intervals
plotdat <- data.frame(
  weight = pred_data$weight,  # Original weight scale
  fit_link = preds$fit,       # Predictions on the link (logit) scale
  se_link = preds$se.fit)      # Standard errors on the link scale
# Calculate confidence intervals on the link scale (95% CI)
plotdat <- plotdat %>%
  mutate(
    lower_link = fit_link - 1.96 * se_link,
    upper_link = fit_link + 1.96 * se_link)
# Convert from link scale (logit) to response scale (probability)
plotdat <- plotdat %>%
  mutate(
    fit_response = plogis(fit_link),         # Predicted probability
    lower_response = plogis(lower_link),     # Lower CI bound
    upper_response = plogis(upper_link))      # Upper CI bound
# Find the weight with maximum breeding success probability
plotdat <- plotdat %>%
  mutate(bestweight = weight[which.max(fit_response)])

plotweightbreed <- {ggplot(data = plotdat, 
                           aes(x = weight,
                               y = fit_response, 
                               ymin = lower_response,
                               ymax = upper_response, 
                               xmin=min(pruned$weight), 
                               xmax=max(pruned$weight))) +
    # Add confidence interval ribbon
    geom_ribbon(alpha = 0.2,
                colour = "darkgreen", 
                fill = "darkgreen", 
                linetype = "dotted") +
    # Add the main prediction line
    geom_line(linetype = "solid", 
              colour = "green4",
              linewidth = 1) +
    # Add axis labels
    xlab("Weight (kg)") +
    ylab("Likelihood of Breeding Success") +
    # Set plot theme
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_text(angle = 0, vjust = 0.5, color = "black"),
          axis.text.x = element_text(color = "black"),
          strip.text.x = element_text(colour = "black"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          text = element_text(size = 14, colour = "black", face = "plain"))
}
plotweightbreed
ggsave("breeding_success_by_weight.png", plotweightbreed, width = 8, height = 6, dpi = 300)


# age ####

# Calculate age mean and sd
age_mean <- mean(pruned$age, na.rm=TRUE)
age_sd <- sd(pruned$age, na.rm=TRUE)
# Create sequence of age values
age_values <- seq(min(pruned$age, na.rm=TRUE), 
                  max(pruned$age, na.rm=TRUE), 
                  length.out=100)
# Get model coefficients
age_coef <- simple_coefs["scale(age)", "Estimate"]
age_se <- simple_coefs["scale(age)", "Std. Error"]  # Get standard error for CI
# Get intercept
intercept <- simple_coefs["(Intercept)", "Estimate"]
# Get mean values for other predictors to hold constant
weight_mean_scaled <- mean(scale(pruned$weight), na.rm=TRUE)
weight_sq_mean <- mean((scale(pruned$weight))^2, na.rm=TRUE)
# Get coefficients for other terms
weight_coef <- simple_coefs["scale(weight)", "Estimate"]
weight_sq_coef <- simple_coefs["I((scale(weight))^2)", "Estimate"]

# Calculate linear predictor for each age value
manual_pred <- data.frame(
  age = age_values,
  age_scaled = (age_values - age_mean) / age_sd)
# Calculate linear predictor using coefficient
manual_pred$linear_pred <- intercept + 
  weight_coef * weight_mean_scaled +
  weight_sq_coef * weight_sq_mean +
  age_coef * manual_pred$age_scaled
# Convert to probability
manual_pred$probability <- plogis(manual_pred$linear_pred)
# Calculate 95% confidence intervals on the link scale
manual_pred$lower_linear <- manual_pred$linear_pred - 1.96 * age_se * abs(manual_pred$age_scaled)
manual_pred$upper_linear <- manual_pred$linear_pred + 1.96 * age_se * abs(manual_pred$age_scaled)
# Convert to probability scale
manual_pred$lower_prob <- plogis(manual_pred$lower_linear)
manual_pred$upper_prob <- plogis(manual_pred$upper_linear)

# Find optimal age based on coefficient direction
best_age <- manual_pred$age[which.max(manual_pred$probability)]

# Create age effect plot
age_plot <- {ggplot(manual_pred, aes(x=age, 
                                     y=probability)) +
    geom_ribbon(aes(ymin=lower_prob, 
                    ymax=upper_prob), 
                alpha=0.2, 
                fill="darkblue") +
    geom_line(color="blue4", 
              size=1) +
    #geom_vline(xintercept=best_age, 
            #   linetype="dashed", 
            #   color="grey41") +
   # annotate("text", 
           #  x=best_age + 0.1,
           #  y=1,
          #   label=paste0(round(best_age, 1), " years"),
           #  color="gray41") +
    labs(x="Age (years)", 
         y="Likelihood of Breeding Success",
         title="Effect of Age on Breeding Success (manual calculation)") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_text(angle = 0, vjust = 0.5, color = "black"),
          axis.text.x = element_text(color = "black"),
          strip.text.x = element_text(colour = "black"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          text = element_text(size = 14, colour = "black", face = "plain"))
}
age_plot

ggsave("age_effect_on_breeding_success.png", plot=age_plot, width = 8, height = 6, dpi = 300 )

# morph ####
morph_violin <- ggplot(pruned, 
                       aes(x = morph, 
                           y = breedingsuccess, 
                           fill = morph)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.07, 
               alpha = 0.7) +
  geom_jitter(width = 0.15, 
              height = 0.02, 
              alpha = 0.5,
              colour = "black") +
  scale_fill_manual(values = c("grey30", "tan4")) + # Adjust colors to match morph colors
  labs(x = "Morph", 
       y = "Breeding Success", 
       title = "Raw Breeding Success by Morph") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size = 14, colour = "black"))
morph_violin

summary(topmod_pruned)$coefficients

morph_coef <- -2.4320712
morph_se <- 0.8101449
morph_p <- -3.002020
intercept <- 1.1873517

# Extract mean values for other predictors to hold constant
weight_mean <- mean(pruned$weight, na.rm=TRUE)
age_mean <- mean(pruned$age, na.rm=TRUE)
weight_scaled <- 0  # Mean of scaled variable is 0
age_scaled <- 0     # Mean of scaled variable is 0

# Create basis for linear predictor (everything except morph effect)
base_linear_pred <- intercept
# Add effects of other variables
base_linear_pred <- base_linear_pred + (2.5050649  * weight_scaled)
base_linear_pred <- base_linear_pred + (-0.6925422 * weight_scaled^2)
base_linear_pred <- base_linear_pred + (-1.1900311 * age_scaled)
# Create prediction data
morph_levels <- levels(pruned$morph)
pred_data <- data.frame(
  morph = factor(morph_levels, levels = morph_levels))
# Calculate linear predictors for each morph
pred_data$linear_pred <- base_linear_pred
# For the Fawn morph, add the morph coefficient
fawn_idx <- which(pred_data$morph == "Fawn")
pred_data$linear_pred[fawn_idx] <- base_linear_pred + morph_coef
# Convert to probability scale
pred_data$predicted_success <- plogis(pred_data$linear_pred)
# Calculate standard errors and confidence intervals
pred_data$se_link <- 0  # Initialize
pred_data$se_link[pred_data$morph == levels(pruned$morph)[1]] <- sqrt(0.5862252^2)
pred_data$se_link[pred_data$morph == "Fawn"] <- sqrt(0.5862252^2 + 0.8101449^2)
# Calculate 95% confidence intervals on link scale
pred_data$lower_link <- pred_data$linear_pred - 1.96 * pred_data$se_link
pred_data$upper_link <- pred_data$linear_pred + 1.96 * pred_data$se_link
# Convert to probability scale
pred_data$lower <- plogis(pred_data$lower_link)
pred_data$upper <- plogis(pred_data$upper_link)

morph_plot <- ggplot(pred_data, aes(x = morph, 
                                    y = predicted_success, 
                                    fill = morph)) +
  geom_bar(stat = "identity", 
           width = 0.6, 
           alpha = 0.7) +
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper), 
                width = 0.2) +
  scale_fill_manual(values = c("black", "tan3")) +
  labs(x = "Morph", 
       y = "Predicted Breeding Success", 
       title = "Predicted Breeding Success by Morph") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(size = 14, colour = "black"))
morph_plot

ggsave("predicted_success_by_morph_manual.png", morph_plot, width = 8, height = 6, dpi = 300)


#### Create combined figure ####
morph_plot2 <- {ggplot(pred_data, aes(x = morph, 
                                      y = predicted_success, 
                                      fill = morph)) +
    geom_bar(stat = "identity", 
             width = 0.6, 
             alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, 
                      ymax = upper), 
                  width = 0.2, 
                  color = "black") +
    scale_fill_manual(values = c("black", "tan3")) +
    labs(x = "Morph",
         y = NULL, 
         title = "C)") +
    ylim(0, 1) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "gray90"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(color = "black", size = 9))}
age_plot2<- {ggplot(manual_pred, aes(x = age, 
                                     y = probability)) +
    geom_ribbon(aes(ymin = lower_prob, 
                    ymax = upper_prob), 
                fill = "blue", alpha = 0.2) +
    geom_line(color = "blue4", size = 1) +
    labs(x = "Age (years)",
         y = NULL, 
         title = "A)") +
    ylim(0, 1) +  # Set consistent y-axis limits
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "gray90"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(color = "black", size = 9))}
weight_plot2 <- {ggplot(plotdat, aes(x = weight, 
                                     y = fit_response)) +
    geom_ribbon(aes(ymin = lower_response, 
                    ymax = upper_response), 
                fill = "darkgreen", alpha = 0.2) +
    geom_line(color = "green4", size = 1) +
    labs(x = "Weight (kilograms)",
         y = NULL, 
         title = "B)") +
    ylim(0, 1) +  # Set consistent y-axis limits
    theme_minimal() +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "gray90"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(color = "black", size = 9))}

combined_figure <- ((age_plot2/weight_plot2) | morph_plot2 )+ 
  plot_layout(ncol=2, widths = c(1.5, 0.7))+
  plot_annotation(title = "Factors Affecting Breeding Success",
                  theme = theme(plot.title = element_text(size = 16, 
                                                          face = "bold", 
                                                          hjust = 0.5)))
combined_figure

combined_grob <- patchworkGrob(combined_figure)
y_label <- textGrob("Predicted likelihood of breeding success", 
                    rot = 90, 
                    gp = gpar(fontsize = 14),
                    vjust = 1)
final_figure <- grid.arrange(
  y_label, combined_grob,
  ncol = 2, 
  widths = c(0.05, 0.95))

ggsave("Factors_affecting_breeding_success.png", final_figure, width = 8, height = 6, dpi = 300)

#### Modelling predictors for number of offspring ####

summary(success)

# 31 observations across 4 years
# 24 unique individuals - 12 males, 12 females


# Running a pseudo hurdle model. First stage was modelling drivers of reproductive
# success/failure. Second stage is modelling, for the successful individuals, 
# drivers of higher fecundity. 

psn_success<-glmmTMB(offspring~ 
                       I((scale(weight))^2)*scale(age)+
                       scale(weight)+
                       scale(age)+
                       scale(indhet)+
                       morph+
                       sex+
                       (1|id)+
                       (1|year),
                     data=success,
                     family="poisson", 
                     na.action=na.fail)
check_overdispersion(psn_success)
check_collinearity((psn_success))
check_model(psn_success) 
summary(psn_success)

# try a nbinom2 model to see how it affects fit
nb_success <- glmmTMB(offspring ~ I((scale(weight))^2)*scale(age)+
                        scale(weight) + 
                        scale(age) + 
                        scale(indhet) + 
                        morph + 
                        sex + 
                        (1|id) + 
                        (1|year),
                      data = success,
                      family = "nbinom2")
check_overdispersion(nb_success)
check_overdispersion(nb_success)
check_model(nb_success)

# Compare
AIC(psn_success, nb_success)

# poisson model has lower AIC and one df less, most parsimonious
# will run dredge though no significant predictors in full model, just to check
psnsel<-dredge(psn_success)
psnsel_sub<-subset(psnsel, delta<5)
psnsel_sub

# Convert to data frame and export
write.csv(psnsel_sub, "model_selection_table.csv", row.names = FALSE)

# select top model because it has the lowest df
topmod_psn<-get.models(psnsel_sub, subset = 1)[[1]]
topmod_psn
summary(topmod_psn)

check_overdispersion(topmod_psn)
check_collinearity(topmod_psn)
check_model(topmod_psn) 

# Model table ####
model_table_offspring <- tab_model(topmod_psn,
                                   transform = "exp",
                                   show.est = TRUE,           
                                   show.stat = TRUE,         
                                   show.p = TRUE,            
                                   show.ci = 0.95,          
                                   show.re.var = TRUE,      
                                   show.aic = FALSE,        
                                   show.r2 = TRUE,       
                                   show.icc = FALSE, 
                                   collapse.ci = TRUE,
                                   pred.labels = c(           # Custom names for predictors
                                     "(Intercept)" = "Intercept",
                                     "scale(age)" = "Age (years)",
                                     "morphFawn" = "Morph (fawn)"),
                                   dv.labels = "Predicted impact",
                                   string.est = "Odds Ratios (95% CI)",
                                   string.stat = "z-statistic",
                                   string.p = "p-value",
                                   file = "offspring_model.html") # Save to file
model_table_offspring

webshot("offspring_model.html", 
        "offspring_model.png",
        vwidth = 800, vheight = 600)

#### Modelling females  ####
# binomial ####
summary(females)
binmod_f<-glmmTMB(breedingsuccess ~ scale(weight)+
                    scale(age)+
                    scale(indhet)+
                    morph+
                    (1|id)+(1|year),
                  family="binomial",
                  data=females,
                  na.action=na.fail) 
check_model(binmod_f)
summary(binmod_f)

plot_model(binmod_f, type="pred",
           ci.lvl = 0.95)

modsel_f<-dredge(binmod_f)
modsel_f<-subset(modsel_f, delta < 5)
modsel_f

binmod_f_avg <- model.avg(subset(modsel_f, delta<2), fit=TRUE)
summary(binmod_f_avg)

check_overdispersion(binmod_f) # can't run these checks on averaged models, but since components the same, is comparable
check_collinearity(binmod_f)
check_model(binmod_f)
# plot binomial for females only ####

weight_mean <- mean(females$weight, na.rm = TRUE) 
weight_sd <- sd(females$weight, na.rm = TRUE)
weight_kg_min <- min(females$weight, na.rm = TRUE)
weight_kg_max <- max(females$weight, na.rm = TRUE)
weight_kg_seq <- seq(from = weight_kg_min, to = weight_kg_max, length.out = 100)

weight_std_seq <- (weight_kg_seq - weight_mean) / weight_sd

intercept <- -0.8555      
weight_coef <- 1.4213      
intercept_se <- 0.5955     
weight_se <- 0.5121        

pred_link <- intercept + weight_coef * weight_std_seq
z_68 <- qnorm(0.84)

pred_se <- sqrt(intercept_se^2 + (weight_se * weight_std_seq)^2)
lower_link <- pred_link - z_68 * pred_se
upper_link <- pred_link + z_68 * pred_se
pred_response <- plogis(pred_link)
lower_response <- plogis(lower_link)
upper_response <- plogis(upper_link)

plotdat_f <- data.frame(
  weight_kg = weight_kg_seq,
  weight_std = weight_std_seq,
  fit = pred_response,
  lower = lower_response,
  upper = upper_response)


f_weight_plot <- ggplot(plotdat_f, 
                        aes(x = weight_kg, 
                            y = fit)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper), 
              alpha = 0.2, 
              fill = "tomato3") +
  labs(title = "Effect of weight on female breeding success",
       x = "Weight (kilograms)",
       y = "Probability of Breeding Success") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10))
f_weight_plot


# model tables ####
model_table_f <- tab_model(binmod_f_avg,
                           transform = "exp",
                           show.est = TRUE,           
                           show.stat = TRUE,          
                           show.p = TRUE,             
                           show.ci = 0.95,            
                           show.re.var = TRUE,        
                           show.aic = FALSE,           
                           show.r2 = TRUE,            
                           show.icc = FALSE,
                           collapse.ci = TRUE,
                           title = "Factors Affecting Female Breeding Success",
                           dv.labels = "Predicted impact", 
                           pred.labels = c(           
                             "cond((Int))" = "Intercept",
                             "cond(morphFawn)" = "Morph (fawn)",
                             "cond(scale(age))" = "Age (years)",
                             "cond(scale(weight))" = "Weight (kg)",
                             "cond(scale(indhet))" = "Heterozygosity"),
                           string.est = "Odds Ratios (95% CI)",
                           string.stat = "z-statistic",
                           string.p = "p-value",
                  file = "fems_model.html"
                  )

model_table_f

webshot("fems_model.html", 
        "fems_model.png",
        vwidth = 800, vheight = 600)

# Get the top model
topmod_f <- get.models(modsel_sub_f, subset = 1)[[1]]
# Or extract them specifically
VarCorr(topmod_f)

# psn ####
fems <- subset(females, females$breedingsuccess==1)
summary(fems) # 13 success - this is a tiny dataset, will model off all females and test success only just to see 


mod_fems<-glmmTMB(offspring~ 
                    scale(weight)+
                    scale(age)+
                    scale(indhet)+
                    morph+
                    (1|id)+
                    (1|year),
                  family="poisson",
                  data=fems, 
                  na.action=na.fail)
summary(mod_fems)
check_overdispersion(mod_fems)
check_collinearity(mod_fems)
ch_f <- check_model(mod_fems)
ch_f # possibly too few observations to really be able to tell anything, 
# observations not lining up with predictions

mod_f<-glmmTMB(offspring~ 
                 scale(weight)+
                 scale(age)+
                 morph+
                 (1|id)+
                 (1|year),
               family="poisson",
               data=females,
               ziformula = ~.,
               na.action=na.fail)
summary(mod_f)
check_model(mod_f)# much better fit of residuals 

# # # # # # # #
modself<-dredge(mod_f)
modsel_sub_f<-subset(modself, delta < 5)
modsel_sub_f

# there are multiple models with low delta, so will run model averaging to 
# account for all significant options. No suddden changes in df or AIC. 
topmod_f<-get.models(modsel_sub_f, subset = 1)[[1]]
summary (topmod_f)

# no effects on the number of offspring that a female has, outside of 
# effects on whether she successfully reproduces or not (see binomial)

#### Modelling male ####
summary(males)
mals <- subset(males, males$breedingsuccess==1)
summary(mals) # 18 success

# binomial ####
binmod_m<-glmmTMB(breedingsuccess ~ scale(weight)+
                    scale(age)+
                    morph+
                    (1|id)+(1|year),
                  family="binomial",
                  data=males,
                  na.action=na.fail) 
check_model(binmod_m)
summary(binmod_m)

plot_model(binmod_m, type="pred",
           ci.lvl = 0.95)

modsel_m<-dredge(binmod_m)
modsel_m<-subset(modsel_m, delta < 5)
modsel_m
topmod_m <- get.models(modsel_m, subset = 1)[[1]]
summary(topmod_m)


check_overdispersion(topmod_m)
check_collinearity(topmod_m)
# no significant results, although weight has strong support for affecting 
# breeding success

# plot binomial for males only ####

pred_data_m <- data.frame(
  weight = seq(min(males$weight, na.rm=TRUE), max(males$weight, na.rm=TRUE), length.out=100),
  age = mean(males$age, na.rm=TRUE),
  morph = levels(males$morph)[1],
  id = NA)

preds_m <- predict(topmod_m, newdata=pred_data_m, se.fit=TRUE, re.form=NA)

# Create a data frame with predictions and confidence intervals
plotdat_m <- data.frame(
  weight_m = pred_data_m$weight,  # Original weight scale
  fit_link_m = preds_m$fit,       # Predictions on the link (logit) scale
  se_link_m = preds_m$se.fit)      # Standard errors on the link scale

# Calculate confidence intervals on the link scale (95% CI)
plotdat_m <- plotdat_m %>%
  mutate(
    lower_link_m = fit_link_m - 1.96 * se_link_m,
    upper_link_m = fit_link_m + 1.96 * se_link_m)
# Convert from link scale (logit) to response scale (probability)
plotdat_m <- plotdat_m %>%
  mutate(
    fit_response_m = plogis(fit_link_m),         # Predicted probability
    lower_response_m = plogis(lower_link_m),     # Lower CI bound
    upper_response_m = plogis(upper_link_m))      # Upper CI bound
# Find the weight with maximum breeding success probability
plotdat_m <- plotdat_m %>%
  mutate(bestweight_m = weight_m[which.max(fit_response_m)])

plotbreed_m <- {ggplot(data = plotdat_m, 
                       aes(x = weight_m,
                           y = fit_response_m, 
                           ymin = lower_response_m,
                           ymax = upper_response_m, 
                           xmin=min(males$weight), 
                           xmax=max(males$weight))) +
    # Add confidence interval ribbon
    geom_ribbon(alpha = 0.2,
                colour = "darkgreen", 
                fill = "darkgoldenrod", 
                linetype = "dotted") +
    # Add the main prediction line
    geom_line(linetype = "solid", 
              colour = "black",
              linewidth = 1) +
    # Add axis labels
    xlab("Weight (kg)") +
    ylab("Likelihood of Breeding Success") +
    # Set plot theme
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_text(angle = 0, vjust = 0.5, color = "black"),
          axis.text.x = element_text(color = "black"),
          strip.text.x = element_text(colour = "black"),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          text = element_text(size = 14, colour = "black", face = "plain"))
}
plotbreed_m

ggsave("breeding_success_males.png", plotbreed_m, width = 8, height = 6, dpi = 300)




# model tables ####
summary(topmod_m)

model_table_m <- tab_model(topmod_m, 
                           transform = "exp",
                           show.est = TRUE,           
                           show.stat = TRUE,          
                           show.p = TRUE,             
                           show.ci = 0.95,            
                           show.re.var = TRUE,        
                           show.aic = FALSE,           
                           show.r2 = TRUE,            
                           show.icc = FALSE,
                           collapse.ci = TRUE,
                           title = "Factors Affecting Male Breeding Success",
                           dv.labels = "Predicted impact", 
                           pred.labels = c(           
                             "Inercept" = "Intercept",
                             "morphFawn" = "Morph (fawn)",
                             "scale(age)" = "Age (years)",
                             "scale(weight)" = "Weight (kg)"),
                           string.est = "Odds Ratios (95% CI)",
                           string.stat = "z-statistic",
                           string.p = "p-value",
                           file = "males_model.html")

model_table_m

webshot("males_model.html", 
        "males_model.png",
        vwidth = 800, vheight = 600)




# psn ####
mod_mals<-glmmTMB(offspring~ 
                    scale(weight)+
                    scale(age)+
                    scale(indhet)+
                    morph+
                    (1|id)+
                    (1|year),
                  family="poisson",
                  data=mals, 
                  na.action=na.fail)
summary(mod_mals)
ch_m <- check_model(mod_mals)
ch_m # possibly too few observations to really be able to tell anything, 
# observations not lining up with predictions

mod_m<-glmmTMB(offspring~ scale(age)+morph+(1|id)+(1|year),
               family="poisson",
               data=males,
               ziformula = ~.,
               na.action=na.fail)
summary(mod_m)
check_model(mod_m)# much better fit of residuals 

# no effects on the number of offspring that a male has, possibly 
# a result of the small dataset. Does potentially indicate that, 
#population level effects of morph, weight, age affect females more that males

#### combined outputs for separated sexes ####
# plot ####
x_scale <- scale_x_continuous(
  limits = c(0.25, 1.75),
  breaks = seq(0, 1.75, by = 0.25),
  expand = c(0, 0))

f_weight_plot <- ggplot(
  plotdat_f,
  aes(x = weight_kg, y = fit)
) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.6,
    fill = "maroon"
  ) +
  x_scale +
  labs(
    title = "A)",
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    plot.title = element_text(face = "bold", size = 12),
    text = element_text(size = 12)
  )
f_weight_plot

plotbreed_m <- ggplot(
  plotdat_m,
  aes(
    x = weight_m,
    y = fit_response_m,
    ymin = lower_response_m,
    ymax = upper_response_m
  )
) +
  geom_ribbon(
    alpha = 0.3,
    colour = "darkgreen",
    fill = "cadetblue1",
    linetype = "dotted"
  ) +
  geom_line(
    colour = "black",
    linewidth = 1
  ) +
  x_scale +
  labs(
    title = "B)",
    x = "Weight (kilograms)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    text = element_text(size = 12)
  )
plotbreed_m


combined_figure_sexes <- (f_weight_plot/plotbreed_m)+ 
  plot_annotation(title = "Factors Affecting Breeding Success",
                  theme = theme(plot.title = element_text(size = 16, 
                                                          face = "bold", 
                                                          hjust = 0.5)))
combined_figure_sexes

combined_grob <- patchworkGrob(combined_figure_sexes)
y_label <- textGrob("Predicted likelihood of breeding success", 
                    rot = 90, 
                    gp = gpar(fontsize = 14),
                    vjust = 1)

final_figure <- grid.arrange(
  y_label, combined_grob,
  ncol = 2, 
  widths = c(0.05, 0.95))

ggsave("Factors__breeding_success+sexes.png", final_figure, width = 8, height = 6, dpi = 300)

#### model disgnostic plots ####
library(gtable)
library(gridExtra)
library(grid)

# Function to create a formatted table grob
create_table_grob <- function(data, title) {
  # Create table
  table_grob <- tableGrob(data, rows = NULL,
                          theme = ttheme_default(
                            base_size = 8,
                            core = list(fg_params = list(hjust = 0, x = 0.05),
                                        padding = unit(c(1.5, 1.5), "mm")),
                            colhead = list(fg_params = list(hjust = 0, x = 0.05),
                                           padding = unit(c(1.5, 1.5), "mm"))
                          ))
  
  # Add title
  title_grob <- textGrob(title, gp = gpar(fontsize = 9, fontface = "bold"),
                         just = "left", x = 0.02)
  
  # Combine title and table
  combined <- gtable_add_rows(table_grob, heights = unit(0.5, "cm"), pos = 0)
  combined <- gtable_add_grob(combined, title_grob, 
                              t = 1, l = 1, r = ncol(combined))
  
  return(combined)
}

# Extract collinearity results for each model
binom_collin <- as.data.frame(check_collinearity(topmod_pruned))
binom_collin[sapply(binom_collin, is.numeric)] <- round(binom_collin[sapply(binom_collin, is.numeric)], 3)

poisson_collin <- as.data.frame(check_collinearity(psn_success))
poisson_collin[sapply(poisson_collin, is.numeric)] <- round(poisson_collin[sapply(poisson_collin, is.numeric)], 3)

fem_collin <- as.data.frame(check_collinearity(binmod_f))
fem_collin[sapply(fem_collin, is.numeric)] <- round(fem_collin[sapply(fem_collin, is.numeric)], 3)

male_collin <- as.data.frame(check_collinearity(topmod_m))
male_collin[sapply(male_collin, is.numeric)] <- round(male_collin[sapply(male_collin, is.numeric)], 3)

# Create collinearity table grobs
col_table1 <- create_table_grob(binom_collin, "Binomial model of breeding success - Collinearity")
col_table2 <- create_table_grob(poisson_collin, "Poisson model of individual fecundity - Collinearity")
col_table3 <- create_table_grob(fem_collin, "Binomial model of female breeding success - Collinearity")
col_table4 <- create_table_grob(male_collin, "Binomial model of male breeding success - Collinearity")

# Create overdispersion data frames and tables
# Model 1: Binomial breeding success
binom_overdisp1 <- check_overdispersion(topmod_pruned)
overdisp_df1 <- data.frame(
  Metric = c("Dispersion ratio", "p-value"),
  Value = c(
    round(binom_overdisp1$dispersion_ratio, 3),
    round(binom_overdisp1$p_value, 3)
  )
)
od_table1 <- create_table_grob(overdisp_df1, "Overdispersion")

# Model 2: Poisson individual fecundity
psn_overdisp <- check_overdispersion(psn_success)
overdisp_df2 <- data.frame(
  Metric = c("Dispersion ratio", "Pearson's Chi-Squared", "p-value"),
  Value = c(
    round(psn_overdisp$dispersion_ratio, 3),
    round(psn_overdisp$chisq, 3),
    round(psn_overdisp$p_value, 3)
  )
)
od_table2 <- create_table_grob(overdisp_df2, "Overdispersion")

# Model 3: Female breeding success
fem_overdisp <- check_overdispersion(binmod_f)
overdisp_df3 <- data.frame(
  Metric = c("Dispersion ratio", "p-value"),
  Value = c(
    round(fem_overdisp$dispersion_ratio, 3),
    round(fem_overdisp$p_value, 3)
  )
)
od_table3 <- create_table_grob(overdisp_df3, "Overdispersion")

# Model 4: Male breeding success
male_overdisp <- check_overdispersion(topmod_m)
overdisp_df4 <- data.frame(
  Metric = c("Dispersion ratio", "p-value"),
  Value = c(
    round(male_overdisp$dispersion_ratio, 3),
    round(male_overdisp$p_value, 3)
  )
)
od_table4 <- create_table_grob(overdisp_df4, "Overdispersion")

# Combine each model's collinearity and overdispersion side-by-side
# Then stack all models vertically
combined_plot <- grid.arrange(
  col_table1, od_table1,
  col_table2, od_table2,
  col_table3, od_table3,
  col_table4, od_table4,
  ncol = 2, nrow = 4,
  widths = c(0.80, 0.20),  # Collinearity gets 75% width, overdispersion 25%
  heights = c(1, 1, 1, 1)
)

# Save
ggsave("model_diagnostics_combined.png", combined_plot, 
       width = 210, height = 120, units = "mm", dpi = 300, bg = "white")
