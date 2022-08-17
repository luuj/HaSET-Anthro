library(lme4)
library(tidyverse)
library(magrittr)
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
dat %<>% group_by(cuuid) %>% mutate(n_obs = n(), cuuid=cur_group_id()) %>% ungroup()




##### Modeling #####
# LME with random intercept
mod1 <- lmer(length ~ age_s + (1|cuuid), data=dat) 

# LME with random intercept and quadratic age
mod2 <- lmer(length ~ age_s + I(age_s^2) + (1|cuuid), data=dat)

# LME with random intercept+slope and quadratic age
mod3 <- lmer(length ~ age_s + I(age_s^2) + (1 + age_s|cuuid), data=dat)

# Create dataset with both observed and predicted values
pred.v.obs <- dat %>% select(age,length,visit,cuuid,n_obs) %>% 
  mutate(p_length1=predict(mod1), p_length2=predict(mod2), p_length3=predict(mod3))



##### Plotting #####
# Scatter plot of age*length vs. WHO median
who_line <- who %>% subset(age<25) 
al_scatter <- dat %>% ggplot(aes(x=age, y=length)) + geom_point(shape=3) + 
  scale_x_continuous(name="Age (months)", breaks=c(0, 1, 6, 12, 24)) + 
  geom_line(data=who_line, aes(x=age, y=standard, color=who, linetype=who))
al_scatter

# Scatter plot of predicted vs. observed length from model1
# Plot observed vs predicted length (model1)
pred.v.obs %>% pivot_longer(c("length", "p_length1")) %>% ggplot(aes(x=age, y=value, color=name)) + geom_point() +
  ggtitle("Predicted vs. observed length (model1)")
# Plot observed vs predicted length (model2)
pred.v.obs %>% pivot_longer(c("length", "p_length2")) %>% ggplot(aes(x=age, y=value, color=name)) + geom_point() +
  ggtitle("Predicted vs. observed length (model2)")
pred.v.obs %>% pivot_longer(c("length", "p_length3")) %>% ggplot(aes(x=age, y=value, color=name)) + geom_point() +
  ggtitle("Predicted vs. observed length (model3)") + 
  geom_line(data=who_line, aes(x=age, y=standard, color=who, linetype=who))

# Histogram of age_s - very non-normal
hist(dat$age_s)



##### Prevalence of stunting #####  
who_sd <- who %>% filter(who=="-2 SD") # 2 standard deviations from WHO median
isStunted <- Vectorize(function(age_in, length_in){
  who_standard <- who_sd[which.min(abs(who_sd$age-age_in)) ,"standard"] # Get corresponding standard
  return(length_in < who_standard)
})

# Add stunted column
dat %<>% mutate(stunt=isStunted(age,length))
# Proportion stunted by group
table(dat$visit, dat$stunt) %>% prop.table(1)

# Use model 3 to predict prevalence of stunting
pred.stunt <- pred.v.obs %>% select(age,p_length3,visit) %>% mutate(stunt=isStunted(age,p_length3))
table(pred.stunt$visit, pred.stunt$stunt) %>% prop.table(1)



##### Removing unreliable points #####
# Sampled longitudinal trajectories
long_dat <- pred.v.obs %>% filter(n_obs > 2)
sampled_ids <- sample(unique(long_dat$cuuid),5)
long_dat %<>% filter(cuuid %in% sampled_ids) %>% 
  select(cuuid, age, length, p_length3) %>%
  pivot_longer(c("length", "p_length3"))

long_dat %>% ggplot(aes(x=age, y=value, col=cuuid, linetype=name)) + 
  geom_line()









