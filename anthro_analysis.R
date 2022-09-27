library(lme4)
library(tidyverse)
library(magrittr)
library(gridExtra)
load("~/Library/CloudStorage/OneDrive-HarvardUniversity/School/Harvard/HaSET/Code/anthro_data.RData")

##### Exploratory analysis and data wrangling #####
# Modify dataset
dat <- anthro %>% subset(age >= 0 & age < 1000 & !is.na(length)) # Remove invalid age and missing lengths - 13418 observations

# Categorizing age - 10787 observations out of 13418 - 80%
dat %<>% mutate(visit = case_when(age < 4 ~ "Birth", 
                                  age > 5 & age < 12 ~ "Day 6",
                                  age > 20 & age < 40 ~ "Day 28",
                                  age > 40 & age < 60 ~ "Day 42",
                                  age > 165 & age < 200 ~ "6 months",
                                  age > 350 & age < 380 ~ "12 months",
                                  age > 715  ~ "24 months")) %>% 
  mutate(visit = factor(visit,levels = c("Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months")))

# Convert age to months
dat %<>% mutate(age=age/30.437)
# Create standardized age variable
dat %<>% mutate(age_s = scale(age))
# Remove unnecessary variables
dat %<>% select(-c(cdob, intdt, dcid)) 
# Number of observations per child + relabel ID
dat %<>% group_by(cuuid) %>% mutate(n_obs = n(), cuuid=factor(cur_group_id())) %>% ungroup()
# Subset WHO standard data
who_line <- who %>% subset(age<25) 
who_sd <- who %>% filter(who=="-2 SD") # 2 standard deviations from WHO median

# Add stunted column
isStunted <- Vectorize(function(age_in, length_in){
  who_standard <- who_sd[which.min(abs(who_sd$age-age_in)) ,"standard"] # Get corresponding standard
  return(length_in < who_standard)
})
dat %<>% mutate(stunt=isStunted(age,length))




##### Modeling #####
# LME with random intercept
mod1 <- lmer(length ~ age_s + (1|cuuid), data=dat) 

# LME with random intercept+slope
mod2 <- lmer(length ~ age_s + (1+age_s|cuuid), data=dat) 

# LME with random intercept and quadratic age
mod3 <- lmer(length ~ age_s + I(age_s^2) + (1|cuuid), data=dat)

# LME with random intercept+slope and quadratic age
mod4 <- lmer(length ~ age_s + I(age_s^2) + (1 + age_s|cuuid), data=dat)

# Create dataset with both observed and predicted values
pred.v.obs <- dat %>% select(age,age_s,length,visit,cuuid,n_obs) %>% 
  mutate(p_length1=predict(mod1), p_length2=predict(mod2), p_length3=predict(mod3), p_length4=predict(mod4))




##### Removing unreliable points #####
# Obtain dataset values with more than a single observation: n=12587
dat.long <- pred.v.obs %>% filter(n_obs >= 2)

# Monotonically increasing:
# Function to remove non-monotonically increasing points
removeMI <- function(dat.in){
  id.rm <- which(dat.in$length - lag(dat.in$length) < 0) # Function to remove non-increasing rows
  
  while (length(id.rm) > 0){
    dat.in <- slice(dat.in,-id.rm)
    id.rm <- which(dat.in$length - lag(dat.in$length) < 0) # Function to remove non-increasing rows
  }
  return(dat.in)
}
dat.mi <- dat %>% group_by(cuuid) %>% group_modify(~ removeMI(.)) # n=11657 after removal
dat.mi.nb <- dat %>% filter(visit != "Birth" | is.na(visit)) %>% group_by(cuuid) %>% group_modify(~ removeMI(.)) # No birth version

# Function that removes points not within 1.5sd, given the model
removesd <- function(mod.in, dat.in){
  sd_val <- summary(mod.in)$sigma
  dat.sd <- dat.in %>% select(cuuid,age,age_s,length,visit,n_obs) %>% 
    mutate(p_length = predict(mod.in)) %>%
    filter(length >= (p_length-1.5*sd_val) & length <= (p_length+1.5*sd_val))
  return(dat.sd)
}

# Function that calculates prevalence for the given data
prevCalc <- function(dat.in, modelname=""){
  dat.in <- dat.in %>% group_by(cuuid,visit) %>% group_modify(~ averageMulti(.)) %>% ungroup()
  pred.stunt <- dat.in %>% select(age,length,visit) %>% mutate(stunt=isStunted(age,length))
  prop <- table(pred.stunt$visit, pred.stunt$stunt) %>% prop.table(1) %>% .[,2]
  out <- tibble(c("Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months"), prop, modelname)
  colnames(out) <- c("Category", "Proportion", "Model")
  #cat("Total observations:", sum(out$Frequency), " Missing observations:",sum(is.na(pred.stunt$visit)),"\n")
  return(out)
}

# Function that averages out values if more than one observation per individual in a given time window
averageMulti <- function(dat.in){
  # If more than one observation in a window, average them out
  if (nrow(dat.in) > 1){
    avg <- dat.in %>% select(age,length) %>% apply(2,mean)
    dat.out <- data.frame(t(avg))
    return(dat.out)
  }else{
    return(dat.in)
  }
}

# Re-run 6 models but dataset excludes birth category
dat.nobirth <- dat %>% filter(visit != "Birth" | is.na(visit))

mod1.nb <- lmer(length ~ age_s + (1|cuuid), data=dat.nobirth) 
mod2.nb <- lmer(length ~ age_s + (1+age_s|cuuid), data=dat.nobirth) 
mod3.nb <- lmer(length ~ age_s + I(age_s^2) + (1|cuuid), data=dat.nobirth)
mod4.nb <- lmer(length ~ age_s + I(age_s^2) + (1 + age_s|cuuid), data=dat.nobirth)

allPrev.2sd <- rbind(prevCalc(removesd(mod1,dat),"Mod1"),
                     prevCalc(removesd(mod2,dat),"Mod2"),
                     prevCalc(removesd(mod3,dat),"Mod3"),
                     prevCalc(removesd(mod4,dat),"Mod4"),
                     prevCalc(removesd(mod1.nb,dat.nobirth),"Mod1 NB"),
                     prevCalc(removesd(mod2.nb,dat.nobirth),"Mod2 NB"),
                     prevCalc(removesd(mod3.nb,dat.nobirth),"Mod3 NB"),
                     prevCalc(removesd(mod4.nb,dat.nobirth),"Mod4 NB"),
                     prevCalc(dat.mi,"Monotonically increasing"),
                     prevCalc(dat.mi.nb,"Monotonically increasing NB"),
                     prevCalc(dat,"Observed"))
allPrev.2sd %<>% pivot_wider(names_from = Category, values_from = Proportion)

# 10787 out of 13418 individuals who fall in a visit category
table(dat$visit)

# Numbers from each model removed
removedPoints <- rbind(table(removesd(mod1,dat)$visit), 
                       table(removesd(mod2,dat)$visit),
                       table(removesd(mod3,dat)$visit),
                       table(removesd(mod4,dat)$visit),
                       table(removesd(mod1.nb,dat.nobirth)$visit),
                       table(removesd(mod2.nb,dat.nobirth)$visit),
                       table(removesd(mod3.nb,dat.nobirth)$visit),
                       table(removesd(mod4.nb,dat.nobirth)$visit),
                       table(dat.mi$visit),
                       table(dat.mi.nb$visit),
                       table(dat$visit))

removedPoints <- removedPoints %>% data.frame()
removedPoints$total <- apply(removedPoints, 1, sum)
colnames(removedPoints) <- c("Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months","Total")
rownames(removedPoints) <- allPrev.2sd$Model
removedPoints

# Incidence + reversal condition only on the previous time point - who was not stunted, look to see how many became stunted
incCalc <- function(dat.in, modelname=""){
  # Average out multiple instances per window
  dat.in <- dat.in %>% group_by(cuuid,visit) %>% group_modify(~ averageMulti(.)) %>% ungroup()
  pred.stunt <- dat.in %>% select(cuuid,age,length,visit) %>% mutate(stunt=isStunted(age,length))

  # Calculate incidence by category
  selectID <- pred.stunt %>% pull(cuuid)
  categories <- c("Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months")
  out <- data.frame(matrix(ncol = length(categories), nrow = 0))
  colnames(out) <- categories
  for (i in categories){
    stuntCount <- pred.stunt %>% filter(visit==i & stunt==TRUE & cuuid %in% selectID) %>% nrow
    if (stuntCount == 0 | i == "Birth"){
      out[1,i] <- NA
      out[2,i] <- NA
    }
    else{
      prevN <- pred.stunt %>% filter(visit==i & cuuid %in% selectID) %>% nrow()
      out[1,i] <- stuntCount/prevN
      out[2,i] <- prevN
    }
    
    # Only keep individuals in previous cohort
    selectID <- pred.stunt %>% filter(visit==i & stunt==FALSE) %>% pull(cuuid)
  }
  
  return(tibble(Model=modelname,round(out,3)))
}

allInc.2sd <- rbind(incCalc(removesd(mod1,dat),"Mod1"),
                     incCalc(removesd(mod2,dat),"Mod2"),
                     incCalc(removesd(mod3,dat),"Mod3"),
                     incCalc(removesd(mod4,dat),"Mod4"),
                     incCalc(removesd(mod1.nb,dat.nobirth),"Mod1 NB"),
                     incCalc(removesd(mod2.nb,dat.nobirth),"Mod2 NB"),
                     incCalc(removesd(mod3.nb,dat.nobirth),"Mod3 NB"),
                     incCalc(removesd(mod4.nb,dat.nobirth),"Mod4 NB"),
                     incCalc(dat.mi,"Monotonically increasing"),
                     incCalc(dat.mi.nb,"Monotonically increasing NB"),
                     incCalc(dat,"Observed"))

allInc.2sd[rep(c(T,F),11),] # Get incidence
allInc.2sd[rep(c(F,T),11),] # Get incidence counts




# Remove facility visits after birth
# Reversal numbers
# Use weight at all? More accurate measurements





##### Plotting #####
# nrow(anthro) # 19078 observations
# unique(anthro$cuuid) %>% length # 4273 children
# anthro %>% count(cuuid) %>% pull(n) %>% table() # Number of observations per child

# # Scatter plot of age*length vs. WHO median
al_scatter <- dat %>% mutate(age=age*30.437) %>% ggplot(aes(x=age, y=length)) + geom_point(shape=1, color="orange") +
  scale_x_continuous(name="Age (days)", breaks=c(0, 28, 180, 365, 730))
al_scatter +  theme(legend.position = "none") + ylab("Length") + geom_vline(xintercept = 165, linetype="dotted", color="Blue") +
  geom_vline(xintercept = 200, linetype="dotted", color="Blue") + geom_vline(xintercept = 350, linetype="dotted", color="Red") +
  geom_vline(xintercept = 380, linetype="dotted", color="Red") + geom_vline(xintercept = 715, linetype="dotted", color="green") +
  geom_vline(xintercept = 745, linetype="dotted", color="green") + geom_vline(xintercept = 20, linetype="dotted", color="black") +
  geom_vline(xintercept = 40, linetype="dotted", color="black") + geom_vline(xintercept = 0, linetype="dotted", color="cyan") +
  geom_vline(xintercept = 4, linetype="dotted", color="cyan") + ggtitle("Scatter plot of Age*Length with age windows ")

# # Scatter plot function of predicted vs. observed length
# plotOP <- function(dat_in, modelText){
#   dat_in %>% ggplot(aes(x=age, y=value, color=name)) + geom_point() +
#     ggtitle(paste("Predicted vs. observed length",modelText)) +
#     #geom_line(data=who_line, aes(x=age, y=standard, group=who), inherit.aes = FALSE, size=1) +
#     scale_color_manual(labels = c("Observed", "Predicted"), values = c("orange", "darkturquoise")) +
#     labs(y="Length", x="Age",color="")
# }
# 
# # Plot observed vs predicted length for all models
# plot.mod1 <- pred.v.obs %>% pivot_longer(c("length", "p_length1")) %>% plotOP("Model1")
# plot.mod2 <- pred.v.obs %>% pivot_longer(c("length", "p_length2")) %>% plotOP("Model2")
# plot.mod3 <- pred.v.obs %>% pivot_longer(c("length", "p_length3")) %>% plotOP("Model3")
# plot.mod4 <- pred.v.obs %>% pivot_longer(c("length", "p_length4")) %>% plotOP("Model4")
# plot.mod5 <- pred.v.obs %>% pivot_longer(c("length", "p_length5")) %>% plotOP("Model5")
# plot.mod6 <- pred.v.obs %>% pivot_longer(c("length", "p_length6")) %>% plotOP("Model6")
# grid.arrange(plot.mod1,plot.mod2,plot.mod3,plot.mod4,plot.mod5,plot.mod6)

# # Sampled longitudinal trajectories and their predicted lines
# plotL <- function(dat.in,plot.id){
#   #sampled_ids <- sample(unique(dat.in$cuuid),7)
#   sampled_ids <- plot.id
#   plot.long <- dat.in %>% filter(cuuid %in% sampled_ids) %>%
#     select(cuuid, age, length, p_length4) %>%
#     pivot_longer(c("length", "p_length4"))
#   plot.long %>% ggplot(aes(x=age, y=value, color=cuuid, linetype=name)) + geom_line() +
#     xlab("Age") + ylab("Length")
# }
# plot.id <- c(3313,2528,2622,333,3841)
# plot.f <- plotL(dat.long,plot.id) + ggtitle("Full data")
# plot.2sd <- plotL(dat.2sd,plot.id) + ggtitle("Data with 2 SD removed")
# plot.2sdmi <- plotL(dat.mi,plot.id) + ggtitle("Data with 2 SD removed and MI")
# library(ggpubr)
# ggarrange(plot.f,plot.2sd,plot.2sdmi, ncol=3, common.legend = TRUE)

# # Only take one observation per individual in each window
# checkMulti <- function(dat.in){
#   empty <- data.frame()
#   if (nrow(dat.in) > 1)
#     return(dat.in)
#   else
#     return(empty)
# }
# 
# dat.12 <- dat %>% filter(visit == "12 months") %>% group_by(cuuid) %>% group_modify(~ checkMulti(.))
# test <- dat.12 %>% group_by(cuuid,visit) %>% group_modify(~ averageMulti(.)) %>% ungroup()
# ggplot(aes(x=age, y=length, color=cuuid), data = test) + geom_point() +
#   xlab("Age (months)") + ylab("Length") + ggtitle("Individuals with >1 measurement around 12 months (n=135)") +
#   theme(legend.position = "none") + xlim(c(0,24)) + ylim(c(38,100)) + geom_vline(xintercept = 10.22, linetype="dotted") + geom_vline(xintercept = 13.44, linetype="dotted")
# View(dat.2sd %>% group_by(cuuid,visit) %>% filter(!is.na(visit)) %>% group_modify(~ checkMulti(.))) #5781 entries!





##### OLD CODE #####
# # Categorizing age into birth, 28 days, 6 months, 12 months, 24 months
# dat %<>% mutate(visit = case_when(age < 10 ~ "Birth", 
#                                   age >= 15 & age < 50 ~ "1 month",
#                                   age > 130 & age < 230 ~ "6 months",
#                                   age > 310 & age < 410 ~ "12 months",
#                                   age > 670 & age < 800 ~ "24 months")) %>% 
#   mutate(visit = factor(visit,levels = c("Birth","1 month","6 months","12 months","24 months")))
# 
# # Append splines to dataset
# dat_splines <- dat %>% mutate(spline1 = case_when(age > 0 & age <= 1 ~ age, age==0 ~ 0, TRUE ~ 1),
#                               spline2 = case_when(age > 1 & age <= 6 ~ age-1, age <= 1 ~ 0, TRUE ~ 5),
#                               spline3 = case_when(age > 6 & age <= 12 ~ age-6, age <= 6 ~ 0, TRUE ~ 6),
#                               spline4 = case_when(age > 12 & age <= 29 ~ age-12, age <= 12 ~ 0, TRUE ~ 12))

# # LME with random intercept and piecewise-linear age
# mod5 <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1|cuuid), data=dat_splines) 
# 
# # LME with random intercept+slope and piecewise-linear age
# mod6 <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1 + age_s|cuuid), data=dat_splines) 

# dat_splines.nobirth <- dat_splines %>% filter(visit != "Birth" | is.na(visit))
# mod5.nb <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1|cuuid), data=dat_splines.nobirth) 
# mod6.nb <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1 + age_s|cuuid), data=dat_splines.nobirth) 
