library(tidyverse)
library(cvequality)
library(car)
library(readxl)
###Germination data

germdata <- read_xlsx("Price et al-AJB 2020-AppendixS3.xlsx", sheet = "Germination")


#Create Inbred, a boolean value indicating if a seed is inbred or not
germdata$Inbred <- 0
for(i in 1: length(rownames(germdata))){
  if(germdata[i,7] == 0){
    germdata[i, 9] <- "N"
  }else{germdata[i, 9] <- "Y"}
}
#Format
germdata$Time <- as.numeric(germdata$Time)
germdata$Length_mm <- as.numeric(germdata$Length_mm)
#calculate estimated time to germination
germdata$Est_time <- germdata$Time - germdata$Length_mm * 2


#Calculate group means and Mann-whitney U test for significant change in trait value due to inbreeding
ExpFam_germ_Mean <- germdata %>% group_by(Family, Inbred) %>% summarise(mean = mean(Est_time, na.rm = T))
germdata %>% group_by(Family, Parent, Inbred) %>% summarise(mean = mean(Est_time, na.rm = T))




G_p <- germdata %>% group_by(Family) %>% summarise(p = wilcox.test(Est_time~Inbred)$p.value)
germdata %>% group_by(Family, Parent) %>% summarise(p = wilcox.test(Est_time~Inbred)$p.value)


#Calculate inbreeding
IG <- ExpFam_germ_Mean %>% filter(Inbred == "Y")
OG <- ExpFam_germ_Mean %>% filter(Inbred == "N")

IBG <- 1-(IG[3]/OG[3])
IBG$Family <- c("FS1", "HS4", "HS7")

#Count observations
germdata %>% group_by(Family, Inbred) %>% summarise(count= length(which(!is.na(Est_time))))

#Logisitc regression for germination
summary(glm(Germinate~Inbred, family = "binomial", data = filter(germdata, Family == "FS1")))
summary(glm(Germinate~Inbred, family = "binomial", data = filter(germdata, Family == "HS4")))
summary(glm(Germinate~Inbred, family = "binomial", data = filter(germdata, Family == "HS7")))




#######Adult plant data
dataA <- read_xlsx("Price et al-AJB 2020-AppendixS3.xlsx", sheet = "Adult Plant")

dataA$TD1 <- as.numeric(dataA$TD1)    
dataA$TD2 <- as.numeric(dataA$TD2) 
dataA$H1 <- as.numeric(dataA$H1) 
dataA$H2 <- as.numeric(dataA$H2) 
dataA$BD1 <- as.numeric(dataA$BD1)
dataA$BD2 <- as.numeric(dataA$BD2) 
dataA$Stalks <- as.numeric(dataA$Stalks)
dataA$Flowers <- as.numeric(dataA$Flowers)
dataA$Total_seed_weight <- as.numeric(dataA$Total_seed_weight)
dataA$weight_per_seed <- as.numeric(dataA$weight_per_seed)
dataA$Seed_fill <- as.numeric(dataA$Seed_fill)
dataA$FT19 <- as.numeric(dataA$FT19)
dataA$FT219 <- as.numeric(dataA$FT219)
#Calculate means and derived values
dataA$TDmean <- rowMeans(dataA[,c(11,12)], na.rm = T)/1000
dataA$Hmean <- rowMeans(dataA[,c(13,14)], na.rm = T)
dataA$BDmean <- rowMeans(dataA[,c(16,17)], na.rm = T)/1000
dataA$Hcm <- dataA$Hmean *2.54
dataA$FPS <- dataA$Flowers/dataA$Stalks
dataA$est_y <- dataA$Flowers*dataA$Total_seed_weight

#Change estimated yield to 0 for plants with 0 flowers (rather than "NA")
for(i in 1: length(dataA$ID)){
  if(dataA[i,7] == 0 & is.na(dataA[i,7])==F){
    dataA[i, 30] <- 0
  } 
}

#Recode Inbred sign to better represent signified. Now a Boolean

dataA$Inbred <- as.character(dataA$Inbred)

for(i in 1: length(dataA$ID)){
  if(dataA[i,15] == 1){
    dataA[i, 15] <- "N"
  }else{dataA[i, 15] <- "Y"}
}
dataA$Inbred <- as.factor(dataA$Inbred)

#Logistic regression for survival


summary(glm(Survival~Inbred, family = "binomial", data = filter(dataA, Family == "FS1")))
summary(glm(Survival~Inbred, family = "binomial", data = filter(dataA, Family == "HS4")))
summary(glm(Survival~Inbred, family = "binomial", data = filter(dataA, Family == "HS7")))

#Percent survival 
mean(dataA$Survival)


#Give a unique code for each maternal family
dataA$Mother <- paste(dataA$Family, "-", dataA$Parent)

#Trait means
ExpFam_mean <- dataA %>% group_by(Family, Inbred) %>% summarise(FT = mean(FT19, na.rm = T), H = mean(Hcm, na.rm = T), BD= mean(BDmean, na.rm = T), TD = mean(TDmean, na.rm = T), S = mean(Stalks, na.rm = T), Fll = mean(Flowers, na.rm = T), FPSq= mean(FPS, na.rm = T), est= mean(est_y, na.rm = T), wps= mean(weight_per_seed, na.rm = T), ss= mean(Seed_fill, na.rm = T) )
MatFam_means <- dataA %>% group_by(Family, Parent, Inbred) %>% summarise(FT = mean(FT19, na.rm = T), H = mean(Hcm, na.rm = T), BD= mean(BDmean, na.rm = T), TD = mean(TDmean, na.rm = T), S = mean(Stalks, na.rm = T), Fll = mean(Flowers, na.rm = T), FPSq= mean(FPS, na.rm = T), est= mean(est_y, na.rm = T), wps= mean(weight_per_seed, na.rm = T), ss= mean(Seed_fill, na.rm = T) )

#Calculate Inbreeding 
IS <- ExpFam_mean %>% filter(Inbred == "Y")
FR <- ExpFam_mean %>% filter(Inbred == "N")

inbreeding <- 1-(IS[c(3:12)]/FR[c(3:12)])
inbreeding$Family <- c("FS1", "HS4", "HS7")
#Count observations
dataA %>% group_by(Family, Inbred) %>% summarise(FT = length(which(!is.na(FT19))), H = length(which(!is.na(Hcm))), BD= length(which(!is.na(BDmean))), TD = length(which(!is.na(TDmean))), S = length(which(!is.na(Stalks))), Fll = length(which(!is.na(Flowers))), FPSq= length(which(!is.na(FPS))), est= length(which(!is.na(est_y))), wps= length(which(!is.na(weight_per_seed))), ss= length(which(!is.na(Seed_fill))))

#Mann-whitney U test for significant change in trait value due to inbreeding
#Experimental Family level
ExpFam_p <- dataA %>% group_by(Family) %>% summarise(FT = wilcox.test(FT19~Inbred)$p.value, H = wilcox.test(Hcm~Inbred)$p.value, BD= wilcox.test(BDmean~Inbred)$p.value, TD = wilcox.test(TDmean~Inbred)$p.value, S = wilcox.test(Stalks~Inbred)$p.value, Fll = wilcox.test(Flowers~Inbred)$p.value, FPSq= wilcox.test(FPS~Inbred)$p.value, est= wilcox.test(est_y~Inbred)$p.value, wps= wilcox.test(weight_per_seed~Inbred)$p.value, ss= wilcox.test(Seed_fill~Inbred)$p.value )


##Difference in trait value due to type of treatment

Type_p <- dataA %>% filter(Inbred == "N") %>% group_by(Family) %>% summarise(FT = wilcox.test(FT19~Type)$p.value, H = wilcox.test(Hcm~Type)$p.value, BD= wilcox.test(BDmean~Type)$p.value, TD = wilcox.test(TDmean~Type)$p.value, S = wilcox.test(Stalks~Type)$p.value, Fll = wilcox.test(Flowers~Type)$p.value, FPSq= wilcox.test(FPS~Type)$p.value, est= wilcox.test(est_y~Type)$p.value, wps= wilcox.test(weight_per_seed~Type)$p.value, ss= wilcox.test(Seed_fill~Type)$p.value )
dataA %>% filter(Inbred == "N") %>% group_by(Family, Type) %>% summarise(Stalk_mean = mean(Stalks, na.rm = T), Flower_mean = mean(Flowers, na.rm = T), est_y_mean = mean(est_y, na.rm = T))



#Maternal Family level significance of change due to inbreeding
MaternalFam_p <- dataA %>% group_by(Family, Parent) %>% summarise(FT = wilcox.test(FT19~Inbred)$p.value, H = wilcox.test(Hcm~Inbred)$p.value, BD= wilcox.test(BDmean~Inbred)$p.value, TD = wilcox.test(TDmean~Inbred)$p.value, S = wilcox.test(Stalks~Inbred)$p.value, Fll = wilcox.test(Flowers~Inbred)$p.value, FPSq= wilcox.test(FPS~Inbred)$p.value, est= wilcox.test(est_y~Inbred)$p.value, wps= wilcox.test(weight_per_seed~Inbred)$p.value, ss= wilcox.test(Seed_fill~Inbred)$p.value )

#Calculate coefficient of variation and test for significant change

ExpFam_cv <- dataA %>% group_by(Family, Inbred) %>% summarise(FT = (sd(FT19, na.rm = T)/ mean(FT19, na.rm = T)), H = (sd(Hcm, na.rm = T)/mean(Hcm, na.rm = T)), BD= (sd(BDmean, na.rm = T)/mean(BDmean, na.rm = T)), TD = (sd(TDmean, na.rm = T)/ mean(TDmean, na.rm = T)), S = (sd(Stalks, na.rm = T)/mean(Stalks, na.rm = T)), Fll = (sd(Flowers, na.rm = T)/ mean(Flowers, na.rm = T)), FPSq = (sd(FPS, na.rm = T)/ mean(FPS, na.rm = T)), est = (sd(est_y, na.rm = T)/mean(est_y, na.rm = T)) , wps = (sd(weight_per_seed, na.rm = T)/ mean(weight_per_seed, na.rm = T)) , ss = (sd(Seed_fill, na.rm = T)/mean(Seed_fill, na.rm = T)))

dataA %>% group_by(Family) %>% filter(is.na(Hcm)==F) %>% summarise(asymptotic_test(Hcm, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(BDmean)==F) %>% summarise(asymptotic_test(BDmean, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(TDmean)==F) %>% summarise(asymptotic_test(TDmean, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(FT19)==F) %>% summarise(asymptotic_test(FT19, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(Stalks)==F) %>% summarise(asymptotic_test(Stalks, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(Flowers)==F) %>% summarise(asymptotic_test(Flowers, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(FPS)==F) %>% summarise(asymptotic_test(FPS, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(Seed_fill)==F) %>% summarise(asymptotic_test(Seed_fill, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(est_y)==F) %>% summarise(asymptotic_test(est_y, Inbred)$p_value)
dataA %>% group_by(Family) %>% filter(is.na(weight_per_seed)==F) %>% summarise(asymptotic_test(weight_per_seed, Inbred)$p_value)

#Shaprio-Wilkes test for deviation from normality for seed set
dataA %>% group_by(Family, Inbred) %>% filter(is.na(Seed_fill)==F) %>% summarise(shapiro.test(Seed_fill)$p.value)


#correlations between trait values for outbreds

FScortable_out <-  dataA %>% filter(Family == "FS1", Inbred == "N") %>% select(Hmean, TDmean, BDmean, Stalks, Flowers, est_y) %>%
  cor(use="pairwise.complete.obs")

HScortable_out <-  dataA %>% filter(Family == "HS7" | Family == "HS4", Inbred == "N") %>% select(Hmean, TDmean, BDmean, Stalks, Flowers, est_y) %>%
  cor(use="pairwise.complete.obs")

#Changes in correlation coefficient due to inbreeding depression
FScortable <- dataA %>% filter(Family == "FS1", Inbred == "Y") %>% select(Hmean, TDmean, BDmean, Stalks, Flowers, est_y) %>%
  cor(use="pairwise.complete.obs") - dataA %>% filter(Family == "FS1", Inbred == "N") %>% select(Hmean, TDmean, BDmean, Stalks, Flowers, est_y) %>%
  cor(use="pairwise.complete.obs")

HScortable <- dataA %>% filter(Family == "HS7" | Family == "HS4", Inbred == "Y") %>% select(Hmean, TDmean, BDmean, Stalks, Flowers, est_y) %>%
  cor(use="pairwise.complete.obs") - dataA %>% filter(Family == "HS7" | Family == "HS4", Inbred == "N") %>% select(Hmean, TDmean, BDmean, Stalks, Flowers, est_y) %>%
  cor(use="pairwise.complete.obs")

##Permutation test for significance of change in correlation coefficient 

corfun  <- function(data, trait1, trait2, to_perm = TRUE){
  data <- data  %>%
    filter(Family  == "FS1" & !is.na(get(trait1)) & !is.na(get(trait2)))
  if(to_perm){data <-   mutate(data,Inbred = sample(Inbred))}  
  #
  data %>%
    group_by(Inbred)                        %>%
    summarise(corAB = cor(get(trait1), get(trait2)))   %>%
    summarise(diff = diff(corAB)) %>%
    pull()
}


perm.dist <- replicate(n = 500, corfun(data = dataA, "TDmean", "BDmean"))
obs.perm  <- corfun(data = dataA, "TDmean", "BDmean", to_perm = FALSE)
p.value   <- mean(abs(perm.dist) >= abs(obs.perm))

traits_to_perm <- c("Hcm", "BDmean", "TDmean", "Stalks", "weight_per_seed", "Flowers")
cor_p_values <- as.data.frame(matrix(nrow = 6, ncol = 6))

for( i in 1:6){
  for( j in 1:6){
    perm.dist <- replicate(n = 500, corfun(data = dataA, traits_to_perm[i], traits_to_perm[j]))
    obs.perm  <- corfun(data = dataA, traits_to_perm[i], traits_to_perm[j], to_perm = FALSE)
    p.value   <- mean(abs(perm.dist) >= abs(obs.perm))
    cor_p_values[i,j] <- p.value
  }
}

colnames(cor_p_values) <- traits_to_perm
row.names(cor_p_values) <- traits_to_perm


###HS

corfunhs  <- function(data, trait1, trait2, to_perm = TRUE){
  data <- data  %>%
    filter(Family  != "FS1" & !is.na(get(trait1)) & !is.na(get(trait2)))
  if(to_perm){data <-   mutate(data,Inbred = sample(Inbred))}  
  #
  data %>%
    group_by(Inbred)                        %>%
    summarise(corAB = cor(get(trait1), get(trait2)))   %>%
    summarise(diff = diff(corAB)) %>%
    pull()
}


traits_to_perm <- c("Hcm", "BDmean", "TDmean", "Stalks", "weight_per_seed", "Flowers")
cor_p_values_HS <- as.data.frame(matrix(nrow = 6, ncol = 6))

colnames(cor_p_values_HS) <- traits_to_perm
row.names(cor_p_values_HS) <- traits_to_perm

for( i in 1:6){
  for( j in 1:6){
    perm.dist <- replicate(n = 500, corfunhs(data = dataA, traits_to_perm[i], traits_to_perm[j]))
    obs.perm  <- corfunhs(data = dataA, traits_to_perm[i], traits_to_perm[j], to_perm = FALSE)
    p.value   <- mean(abs(perm.dist) >= abs(obs.perm))
    cor_p_values_HS[i,j] <- p.value
  }
}
