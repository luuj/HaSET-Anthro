library(lme4)
library(tidyverse)
library(magrittr)
library(gridExtra)
load("~/Library/CloudStorage/OneDrive-HarvardUniversity/School/Harvard/HaSET/Code/anthro_data.RData")

##### Exploratory analysis and data wrangling #####
# nrow(anthro) # 19078 observations
# unique(anthro$cuuid) %>% length # 4273 children
# anthro %>% count(cuuid) %>% pull(n) %>% table() # Number of observations per child

# Modify dataset
dat <- anthro %>% subset(age >= 0 & age < 1000 & !is.na(length)) # Remove invalid age and missing lengths - 13418 observations
# Categorizing age into birth, 28 days, 6 months, 12 months, 24 months
dat %<>% mutate(visit = case_when(age < 10 ~ "Birth", 
                                 age >= 15 & age < 50 ~ "1 month",
                                 age > 130 & age < 230 ~ "6 months",
                                 age > 310 & age < 410 ~ "12 months",
                                 age > 670 & age < 800 ~ "24 months")) %>% 
         mutate(visit = factor(visit,levels = c("Birth","1 month","6 months","12 months","24 months")))
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
# Append splines to dataset
dat_splines <- dat %>% mutate(spline1 = case_when(age > 0 & age <= 1 ~ age, age==0 ~ 0, TRUE ~ 1),
                              spline2 = case_when(age > 1 & age <= 6 ~ age-1, age <= 1 ~ 0, TRUE ~ 5),
                              spline3 = case_when(age > 6 & age <= 12 ~ age-6, age <= 6 ~ 0, TRUE ~ 6),
                              spline4 = case_when(age > 12 & age <= 29 ~ age-12, age <= 12 ~ 0, TRUE ~ 12))
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

# LME with random intercept and piecewise-linear age
mod5 <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1|cuuid), data=dat_splines) 

# LME with random intercept+slope and piecewise-linear age
mod6 <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1 + age_s|cuuid), data=dat_splines) 

# Create dataset with both observed and predicted values
pred.v.obs <- dat %>% select(age,age_s,length,visit,cuuid,n_obs) %>% 
  mutate(p_length1=predict(mod1), p_length2=predict(mod2), p_length3=predict(mod3), 
         p_length4=predict(mod4), p_length5=predict(mod5), p_length6=predict(mod6))




##### Plotting #####
# # Scatter plot of age*length vs. WHO median
# al_scatter <- dat %>% ggplot(aes(x=age, y=length)) + geom_point(shape=3) +
#   scale_x_continuous(name="Age (months)", breaks=c(0, 1, 6, 12, 24)) +
#   geom_line(data=who_line, aes(x=age, y=standard, color=who, linetype=who))
# al_scatter

# # Scatter plot function of predicted vs. observed length 
# plotOP <- function(dat_in, modelText){
#   dat_in %>% ggplot(aes(x=age, y=value, color=name)) + geom_point() +
#     ggtitle(paste("Predicted vs. observed length",modelText)) + 
#     geom_line(data=who_line, aes(x=age, y=standard, group=who), inherit.aes = FALSE, size=1) +
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
# 



##### Removing unreliable points #####
# Obtain dataset values with more than a single observation: n=12587
dat.long <- pred.v.obs %>% filter(n_obs >= 2)

# Removing points outside of 2 SD (model4): n=12089
# 262 values higher than upper boundary, 256 values lower than lower boundary
sd_val <- summary(mod4)$sigma
dat.2sd <- dat.long %>% select(cuuid,age,age_s,length,visit,n_obs,p_length4) %>% 
  filter(length >= (p_length4-2*sd_val) & length <= (p_length4+2*sd_val))

# Monotonically increasing after 2SD: n=7139
checkMI <- function(dat.in){
  all(dat.in == cummax(dat.in))
}
length.id <- dat.2sd %>% group_by(cuuid) %>% group_map(~ pull(.,length))
length.id <- tibble("cuuid" = unique(dat.2sd$cuuid), "isMI" = sapply(length.id,checkMI))
dat.mi <- left_join(dat.2sd, length.id, by="cuuid") %>% filter(isMI == T)

# Not monotonically increasing data: n=4950
removeMI <- function(dat.in){
  id.rm <-which(dat.in$length - lag(dat.in$length) < 0) # Function to remove non-increasing rows
  slice(dat.in,-id.rm)
}
dat.nmi <- left_join(dat.2sd, length.id, by="cuuid") %>% filter(isMI == F)
dat.nmi <- dat.2sd %>% group_by(cuuid) %>% group_modify(~ removeMI(.)) # n=3744 after removal
dat.mi.combined <- full_join(dat.mi, dat.nmi) %>% select(-isMI) # n=10,883 in total

# # Sampled longitudinal trajectories and their predicted lines
# plotL <- function(dat.in){
#   sampled_ids <- sample(unique(dat.in$cuuid),5)
#   plot.long <- dat.in %>% filter(cuuid %in% sampled_ids) %>% 
#     select(cuuid, age, length, p_length4) %>%
#     pivot_longer(c("length", "p_length4"))
#   plot.long %>% ggplot(aes(x=age, y=value, color=cuuid, linetype=name)) + geom_line()
# }
# plotL(dat.long)
# plotL(dat.nobirth)
# plotL(dat.2sd)
# plotL(dat.mi)
# plotL(dat.mi.combined)


# Re-run 6 models but dataset excludes birth category
dat.nobirth <- dat %>% filter(visit != "Birth" | is.na(visit))
dat_splines.nobirth <- dat_splines %>% filter(visit != "Birth" | is.na(visit))

mod1.nb <- lmer(length ~ age_s + (1|cuuid), data=dat.nobirth) 
mod2.nb <- lmer(length ~ age_s + (1+age_s|cuuid), data=dat.nobirth) 
mod3.nb <- lmer(length ~ age_s + I(age_s^2) + (1|cuuid), data=dat.nobirth)
mod4.nb <- lmer(length ~ age_s + I(age_s^2) + (1 + age_s|cuuid), data=dat.nobirth)
mod5.nb <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1|cuuid), data=dat_splines.nobirth) 
mod6.nb <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1 + age_s|cuuid), data=dat_splines.nobirth) 

# 12 sets of prevalences - present them graphically for 2SD algorithm with CI
# Function that removes points not within 2sd, given the model
remove2sd <- function(mod.in, dat.in){
  sd_val <- summary(mod.in)$sigma
  dat.2sd <- dat.in %>% select(cuuid,age,age_s,length,visit,n_obs) %>% 
    mutate(p_length = predict(mod.in)) %>%
    filter(length >= (p_length-2*sd_val) & length <= (p_length+2*sd_val))
  return(dat.2sd)
}

# Function that calculates prevalence for the given data
prevCalc <- function(dat.in, modelname=""){
  pred.stunt <- dat.in %>% select(age,length,visit) %>% mutate(stunt=isStunted(age,length))
  prop <- table(pred.stunt$visit, pred.stunt$stunt) %>% prop.table(1) %>% .[,2]
  out <- tibble(c("Birth","1 month","6 months","12 months","24 months"), prop, modelname)
  colnames(out) <- c("Category", "Proportion", "Model")
  #cat("Total observations:", sum(out$Frequency), " Missing observations:",sum(is.na(pred.stunt$visit)),"\n")
  return(out)
}

allPrev.2sd <- rbind(prevCalc(remove2sd(mod1,dat),"Mod1"),
      prevCalc(remove2sd(mod2,dat),"Mod2"),
      prevCalc(remove2sd(mod3,dat),"Mod3"),
      prevCalc(remove2sd(mod4,dat),"Mod4"),
      prevCalc(remove2sd(mod5,dat),"Mod5"),
      prevCalc(remove2sd(mod6,dat),"Mod6"),
      prevCalc(remove2sd(mod1.nb,dat.nobirth),"Mod1 NB"),
      prevCalc(remove2sd(mod2.nb,dat.nobirth),"Mod2 NB"),
      prevCalc(remove2sd(mod3.nb,dat.nobirth),"Mod3 NB"),
      prevCalc(remove2sd(mod4.nb,dat.nobirth),"Mod4 NB"),
      prevCalc(remove2sd(mod5.nb,dat.nobirth),"Mod5 NB"),
      prevCalc(remove2sd(mod6.nb,dat.nobirth),"Mod6 NB"),
      prevCalc(dat,"Full"))
allPrev.2sd %>% pivot_wider(names_from = Category, values_from = Proportion)

