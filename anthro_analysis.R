library(lme4)
library(tidyverse)
library(magrittr)
library(gridExtra)
load("~/Library/CloudStorage/OneDrive-HarvardUniversity/School/Harvard/HaSET/Code/anthro_data.RData")

##### Exploratory analysis and data wrangling #####
dat <- anthro %>% subset(age >= 0 & age < 1000 & !is.na(length)) # Remove invalid age and missing lengths - 13418 observations

# Categorizing age
dat %<>% mutate(visit = case_when(age <= 3 ~ "Birth", 
                                  age > 3 & age < 10 ~ "Day 6",
                                  age > 20 & age < 35 ~ "Day 28",
                                  age > 35 & age < 50 ~ "Day 42",
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
who_sd <- who %>% filter(who=="-2 SD")

# Add stunted column
isStunted <- Vectorize(function(age_in, length_in){
  who_standard <- who_sd[which.min(abs(who_sd$age-age_in)) ,"standard"] # Get corresponding standard
  return(length_in < who_standard)
})
dat %<>% mutate(stunt=isStunted(age,length))

# Remove facility visits after birth
dat %<>% filter(location=="community" | visit=="Birth")

# Append splines to dataset
dat %<>% mutate(spline1 = case_when(age > 0 & age <= 1 ~ age, age==0 ~ 0, TRUE ~ 1),
                spline2 = case_when(age > 1 & age <= 6 ~ age-1, age <= 1 ~ 0, TRUE ~ 5),
                spline3 = case_when(age > 6 & age <= 12 ~ age-6, age <= 6 ~ 0, TRUE ~ 6),
                spline4 = case_when(age > 12 & age <= 29 ~ age-12, age <= 12 ~ 0, TRUE ~ 12))




##### Modeling #####
# LME with random intercept
mod1 <- lmer(length ~ age_s + (1|cuuid), data=dat, REML=F) 

# LME with random intercept+slope
mod2 <- lmer(length ~ age_s + (1+age_s|cuuid), data=dat, REML=F) 

# LME with random intercept and quadratic age
mod3 <- lmer(length ~ age_s + I(age_s^2) + (1|cuuid), data=dat, REML=F)

# LME with random intercept+slope and quadratic age
mod4 <- lmer(length ~ age_s + I(age_s^2) + (1 + age_s|cuuid), data=dat, REML=F)

# LME with random intercept and piecewise-linear age
mod5 <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1|cuuid), data=dat, REML=F)

# LME with random intercept+slope and piecewise-linear age
mod6 <- lmer(length ~ spline1 + spline2 + spline3 + spline4 + (1 + age_s|cuuid), data=dat, REML=F)





##### Removing unreliable points #####
# Function to remove non-monotonically increasing points
removeMI <- function(dat.in){
  id.rm <- which(dat.in$length - lag(dat.in$length) < 0) # Function to remove non-increasing rows
  
  while (length(id.rm) > 0){
    dat.in <- slice(dat.in,-id.rm)
    id.rm <- which(dat.in$length - lag(dat.in$length) < 0) 
  }
  return(dat.in)
}
dat.mi <- dat %>% group_by(cuuid) %>% group_modify(~ removeMI(.)) # n=11657 after removal

# Function that removes points not within 1.5sd, given the model
removesd <- function(mod.in, dat.in){
  sd_val <- summary(mod.in)$sigma
  dat.sd <- dat.in %>% select(cuuid,age,age_s,length,visit,n_obs) %>% 
    mutate(p_length = predict(mod.in)) %>%
    filter(length >= (p_length-1.5*sd_val) & length <= (p_length+1.5*sd_val))
  return(dat.sd)
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

# Function that calculates prevalence for the given data
prevCalc <- function(dat.in, modelname=""){
  dat.in <- dat.in %>% group_by(cuuid,visit) %>% group_modify(~ averageMulti(.)) %>% ungroup()
  pred.stunt <- dat.in %>% select(age,length,visit) %>% mutate(stunt=isStunted(age,length))
  prop <- table(pred.stunt$visit, pred.stunt$stunt) %>% prop.table(1) %>% .[,2]
  out <- tibble(c("Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months"), prop, modelname)
  colnames(out) <- c("Category", "Proportion", "Model")
  return(out)
}

# Function to calculate incidence and reversal
# Baseline determines if we are calculating incidence/reversal based on shifting time points (6 days vs. 6 months, 6 months vs 12 months) or a baseline time point (6 days vs 6 months, 6 days vs 12 months)
# To select a baseline, pass in one of the categories: "Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months"
inc_and_rev_Calc <- function(dat.in, modelname="", baseline=NULL){
  # Average out multiple instances per window
  dat.in <- dat.in %>% group_by(cuuid,visit) %>% group_modify(~ averageMulti(.)) %>% ungroup()
  pred.stunt <- dat.in %>% select(cuuid,age,length,visit) %>% mutate(stunt=isStunted(age,length))
  
  # Setup storage output
  selectID.inc <- selectID.rev <- pred.stunt %>% pull(cuuid)
  # Choose the categories we are interested in
  categories <- c("Day 6","6 months","12 months","24 months")
  #categories <- c("Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months")
  
  out.inc <- data.frame(matrix(ncol = length(categories), nrow = 0))
  out.rev <- data.frame(matrix(ncol = length(categories), nrow = 0))
  colnames(out.inc) <- categories
  colnames(out.rev) <- categories
  first <- T
  
  # Calculate incidence and reversal by category
  for (i in categories){
    # Obtain number of individuals stunted/not stunted
    stuntCount.inc <- pred.stunt %>% filter(visit==i & stunt==TRUE & cuuid %in% selectID.inc) %>% nrow
    stuntCount.rev <- pred.stunt %>% filter(visit==i & stunt==FALSE & cuuid %in% selectID.rev) %>% nrow
    
    # Ignore first column as we need two time points to calculate incidence/reversal
    if (first==T){
      out.inc[c(1,2),i] <- NA
      out.rev[c(1,2),i] <- NA
      first <- F
    }
    else{ # Obtain denominator and calculate incidence/reversal
      prevN.inc <- pred.stunt %>% filter(visit==i & cuuid %in% selectID.inc) %>% nrow()
      out.inc[1,i] <- stuntCount.inc/prevN.inc
      out.inc[2,i] <- prevN.inc
      
      prevN.rev <- pred.stunt %>% filter(visit==i & cuuid %in% selectID.rev) %>% nrow()
      out.rev[1,i] <- stuntCount.rev/prevN.rev
      out.rev[2,i] <- prevN.rev
    }
    
    # If baseline is null, only keep individuals in previous cohort 
    if (is.null(baseline)){
      selectID.inc <- pred.stunt %>% filter(visit==i & stunt==FALSE) %>% pull(cuuid)
      selectID.rev <- pred.stunt %>% filter(visit==i & stunt==TRUE) %>% pull(cuuid)
    }else{ # If baseline is specified, keep individuals from baseline
      selectID.inc <- pred.stunt %>% filter(visit==baseline & stunt==FALSE) %>% pull(cuuid)
      selectID.rev <- pred.stunt %>% filter(visit==baseline & stunt==TRUE) %>% pull(cuuid)
    }
  }
  return(cbind(tibble(Model=modelname,round(out.inc,3)), tibble(Model=modelname,round(out.rev,3))))
}




##### Generate tables #####
# Prevalence table
allPrev.2sd <- rbind(prevCalc(removesd(mod6,dat),"Mod6"),
                     prevCalc(dat,"Observed"))
allPrev.2sd %<>% pivot_wider(names_from = Category, values_from = Proportion)
allPrev.2sd # Get prevalence

# How many points from each model removed
removedPoints <- rbind(table(removesd(mod6,dat)$visit),
                       table(dat$visit))
removedPoints <- removedPoints %>% data.frame()
removedPoints$total <- apply(removedPoints, 1, sum)
colnames(removedPoints) <- c("Birth","Day 6","Day 28","Day 42","6 months","12 months","24 months","Total")
rownames(removedPoints) <- allPrev.2sd$Model
removedPoints # Get window counts

# Incidence and reversal tables
allIR.2sd <- rbind(inc_and_rev_Calc(removesd(mod6,dat),"Mod6"),
                   inc_and_rev_Calc(dat,"Observed"))

IR_n <- nrow(allIR.2sd)/2
allIR.2sd[rep(c(T,F),IR_n),1:5] # Get incidence
allIR.2sd[rep(c(F,T),IR_n),1:5] # Get incidence counts
allIR.2sd[rep(c(T,F),IR_n),6:10] # Get reversal
allIR.2sd[rep(c(F,T),IR_n),6:10] # Get reversal counts

# Incidence and reversal tables w/ common baseline of Day 6
allIR.2sd.baseline <- rbind(inc_and_rev_Calc(removesd(mod6,dat),"Mod6", "Day 6"),
                   inc_and_rev_Calc(dat,"Observed", "Day 6"))

IR_n <- nrow(allIR.2sd.baseline)/2
allIR.2sd.baseline[rep(c(T,F),IR_n),1:5] # Get incidence w/ common baseline
allIR.2sd.baseline[rep(c(F,T),IR_n),1:5] # Get incidence counts w/ common baseline
allIR.2sd.baseline[rep(c(T,F),IR_n),6:10] # Get reversal w/ common baseline
allIR.2sd.baseline[rep(c(F,T),IR_n),6:10] # Get reversal counts w/ common baseline




##### Model diagnostics #####
AIC(mod1,mod2,mod3,mod4,mod5,mod6)
BIC(mod1,mod2,mod3,mod4,mod5,mod6)
sapply(c(mod1,mod2,mod3,mod4,mod5,mod6),deviance)




##### Plotting #####
# Scatter plot of age*length with age windows
al_scatter <- dat %>% mutate(age=age*30.437) %>% ggplot(aes(x=age, y=length)) + geom_point(shape=1, color="orange") +
  scale_x_continuous(name="Age (days)", breaks=c(0, 28, 180, 365, 730))
al_scatter +  theme(legend.position = "none") + ylab("Length") + geom_vline(xintercept = 165, linetype="dotted", color="Blue") +
  geom_vline(xintercept = 200, linetype="dotted", color="Blue") + geom_vline(xintercept = 350, linetype="dotted", color="Red") +
  geom_vline(xintercept = 380, linetype="dotted", color="Red") + geom_vline(xintercept = 715, linetype="dotted", color="green") +
  geom_vline(xintercept = 745, linetype="dotted", color="green") + geom_vline(xintercept = 40, linetype="dotted", color="yellow") +
  geom_vline(xintercept = 60, linetype="dotted", color="yellow") +geom_vline(xintercept = 20, linetype="dotted", color="black") +
  geom_vline(xintercept = 40, linetype="dotted", color="black") + geom_vline(xintercept = 12, linetype="dotted", color="cyan") +
  geom_vline(xintercept = 5, linetype="dotted", color="cyan") + geom_vline(xintercept = 0, linetype="dotted", color="pink") +
  geom_vline(xintercept = 4, linetype="dotted", color="pink") + ggtitle("Scatter plot of Age*Length with age windows ")

# # Scatter plot function of predicted vs. observed length
plotOP <- function(dat_in, modelText){
  dat_in %>% ggplot(aes(x=age, y=value, color=name)) + geom_point() +
    ggtitle(paste("Predicted vs. observed length",modelText)) +
    scale_color_manual(labels = c("Observed", "Predicted"), values = c("orange", "darkturquoise")) +
    labs(y="Length", x="Age",color="")
}

# Plot observed vs predicted length for model 6
# Create dataset with both observed and predicted values
pred.v.obs <- dat %>% select(age,age_s,length,visit,cuuid,n_obs) %>% mutate(p_length6=predict(mod6))
pred.v.obs %>% pivot_longer(c("length", "p_length6")) %>% plotOP("Model 6")

# Plot removed points for model 6
removed_pts <- anti_join(dat, removesd(mod6,dat)) %>% mutate(removed=T) %>% select(age,length,removed)
plot_pts <- removesd(mod6,dat) %>% mutate(removed=F) %>% select(age,length,removed) %>% rbind(removed_pts)
ggplot(data=plot_pts, aes(x=age,y=length, color=removed)) + geom_point() + 
  xlab("Age") + ylab("Length") + ggtitle("Points removed from model") 

