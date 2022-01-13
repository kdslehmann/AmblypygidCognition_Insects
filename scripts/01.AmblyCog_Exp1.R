#######################################################################################
##                          Amblypygid Cognition Analysis                            ##        
##                            1. Experiment 1                                        ##
##                             by - Colton Watts and                                 ##
##                                  Kenna D. S. Lehmann                              ##
##                                                                                   ##
##                            last updated: 13 Jan 2022                              ##
#######################################################################################
## DESCRIPTION: This script loads, prepares, and analyzes the data from Fiona's 
## experiment with Phrynus marginemaculatus. Experiment 1 in Insects Manuscript
## Lehmann KDS, Shogren FG, Fallick M, Watts JC, Schoenberg D, Wiegmann DD, Bingman VP, Hebets EA. 
## Exploring Higher-Order Conceptual Learning in an Arthropod with a Large Multisensory Processing Center. 
## Insects. 2022; 13(1):81. https://doi.org/10.3390/insects13010081
#######################################################################################

#### ENVIRONMENT PREP ###
rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(4386)

# load libraries
library(dplyr)
library(survival)
library(gplots)
library(viridis)
library(ggplot2)
library(tidyr)
library(here)
library(lubridate)
library(betareg) 
library(plyr)
library(lme4)
library(fitdistrplus)
library(glmmTMB)
library(performance)
library(ggpubr)
library(cowplot)
library(ggfortify)

# Setup -------------------------------------------------------------------
## load data sets 
NonNovelDat = read.csv(here::here('data/Amblypygid Data  (Non-Novel Test Trials)  - Sheet1.csv')) 
NovelDat = read.csv(here::here("data/Amblypygid Data (Novel Test Odors)  - Sheet1.csv")) 

# checking what variables are present in each, usually the names get scrambled a lot during import.
names(NonNovelDat)
names(NovelDat)

## checking to see how much variability there was in the "bookkeeping" variables, e.g., pokes and whatnot
## ignore issues with variable coding at the moment (e.g., there is "y" and "y  ", which R reads differently)
table(NonNovelDat$Probe.Shelters..Y.N.)
table(NovelDat$Probe.Shelters..Y.N.)  ## they almost always probe shelters, so there's not enough variation here to include in subsequent analyses

table(NonNovelDat$Leave.Shelter..Voluntarily..Y.N.)
table(NovelDat$Leave.Shelter..Voluntarily..Y.N.) ## and almost never leave voluntarily, so again not much variation to explore/control for here.

### Clean and standardize data
## Get Individual IDs correctly ordered.
NonNovelDat$X <- factor(NonNovelDat$X, levels= c('Ind 1', 'Ind 2', 'Ind 3', 'Ind 4', 'Ind 5',
                                                 'Ind 6', 'Ind 7', 'Ind 8', 'Ind 9', 'Ind 10',
                                                 'Ind 11', 'Ind 12', 'Ind 13', 'Ind 14'))
NovelDat$X <- factor(NovelDat$X, levels= c('Ind 1', 'Ind 2', 'Ind 3', 'Ind 4', 'Ind 5',
                                                 'Ind 6', 'Ind 7', 'Ind 8', 'Ind 9', 'Ind 10',
                                                 'Ind 11', 'Ind 12', 'Ind 13', 'Ind 14'))

## There's a problem though, in that the times are in minutes:seconds, which isn't conducive to analysis
## need to get it in seconds
## This is why we loaded lubridate up above.
## create a new variable to store the converted times
NonNovelDat$TimeInPlus.Sec = ms(NonNovelDat$Time.in...1) ## step one, converts to a weird format
NonNovelDat$TimeInPlus.Sec = as.numeric(NonNovelDat$TimeInPlus.Sec) ## step two, apply as.numeric() function
NonNovelDat$TimeInPlus.Sec  ## double check. Now it's in seconds! Woo!

## go ahead and do this for other time measures too
NonNovelDat$TimeInMinus.Sec = ms(NonNovelDat$Time.in..) ## step one, converts to a weird format
NonNovelDat$TimeInMinus.Sec = as.numeric(NonNovelDat$TimeInMinus.Sec) ## step two, apply as.numeric() function
NonNovelDat$TimeInMinus.Sec  ## double check. Now it's in seconds! Woo!

NonNovelDat$TimeInNeutral.Sec = ms(NonNovelDat$Time.in.neutral) ## step one, converts to a weird format
NonNovelDat$TimeInNeutral.Sec = as.numeric(NonNovelDat$TimeInNeutral.Sec) ## step two, apply as.numeric()
NonNovelDat$TimeInNeutral.Sec  ## double check. Now it's in seconds! Woo!

## repeat for the other data set
NovelDat$TimeInPlus.Sec = ms(NovelDat$Time.in...1) # convert to weird format
NovelDat$TimeInPlus.Sec = as.numeric(NovelDat$TimeInPlus.Sec) # convert from weird format to seconds
NovelDat$TimeInPlus.Sec # double-check as well.

NovelDat$TimeInMinus.Sec = ms(NovelDat$Time.in..) # convert to weird format
NovelDat$TimeInMinus.Sec = as.numeric(NovelDat$TimeInMinus.Sec) # convert from weird format to seconds
NovelDat$TimeInMinus.Sec # double-check as well.

NovelDat$TimeInNeutral.Sec = ms(NovelDat$Time.in.neutral) # convert to weird format
NovelDat$TimeInNeutral.Sec = as.numeric(NovelDat$TimeInNeutral.Sec) # convert from weird format to seconds
NovelDat$TimeInNeutral.Sec # double-check as well.


# Chi-Squared Test of Test Trials -----------------------------------------
### Okay, so the simplest way to look at these data is to ask whether the first choice was more frequently the + than the -.
## Since each individual is only measured once, can just ask whether the proportion of +s across all trials is significantly different from random (0.5) using a chi-square test.
### Check the same thing for the NonNovelDat
NonNovel.table = table(NonNovelDat$X1st.choice)
NonNovel.table #double-check. Looks right, but definitely looks like no difference.
chisq.test(NonNovel.table) ## and chi-square test confirms there's no real difference.

### Chi-squared with the final, Novel trials
Novel.table = table(NovelDat$X1st.choice)
Novel.table ## double-check, looks good.

chisq.test(Novel.table)  ## all right, get a p-value of 0.03251! They are significantly more likely to choose the positive side.
## can either make a plot of this or just show the table, up to you.

choice_table <- data.frame(Frequency=c(Novel.table, NonNovel.table),
                           trial_type= c('novel','novel','trained','trained'),
                           choice=c('incorrect','correct', 'incorrect','correct'))
choice_table$trial_type <- factor(choice_table$trial_type, levels=c('trained','novel'))
choice_table$choice <- factor(choice_table$choice, levels=c('incorrect', 'correct'))

#### PLOT barchart results ####
choice_ratios_olf <- ggplot(choice_table, aes(fill=choice, y=Frequency, x=trial_type))+
  geom_bar(position="dodge", stat="identity", color="black") + 
  xlab('') +
  ylab('Number of trials') + 
  ylim(c(0,12)) +
  theme_classic()  + 
  scale_fill_manual(values=c("incorrect"='gray', 'correct' ='grey25')) +
  coord_cartesian(ylim=c(0,12)) +
  scale_y_continuous(breaks=seq(0,12,2)) +
  theme(text=element_text(size=14), axis.text=element_text(size=13)) 
choice_ratios_olf

# Proportion of time spent during test trials -----------------------------
### All right, but you probably want to know whether they spent significantly more time on one side than the 
### other. This turns into a bit of a nightmare because you can't directly compare the time spent in one zone 
## to the time spent in another -- an increase in time spent on one side necessarily decreases the time spent 
## on the other side or the neutral zone. I think the "right" way to analyze this
### type of data is with a beta regression, which accounts for the fact that you're modeling a proportion.

## calculate our response variable --  proportion of time in + by dividing by total time spent in either + or -. 
## Compare to null expectation of 0.5.
#For nonnovel dataset
NonNovelDat$PropTimePlus = NonNovelDat$TimeInPlus.Sec/(NonNovelDat$TimeInPlus.Sec + NonNovelDat$TimeInMinus.Sec)
NonNovelDat$PropTimeMinus = NonNovelDat$TimeInMinus.Sec/(NonNovelDat$TimeInPlus.Sec + NonNovelDat$TimeInMinus.Sec)

NonNovelDat$PropTimePlus.all <- NonNovelDat$TimeInPlus.Sec/(NonNovelDat$TimeInPlus.Sec + NonNovelDat$TimeInMinus.Sec + NonNovelDat$TimeInNeutral.Sec)
NonNovelDat$PropTimeMinus.all <- NonNovelDat$TimeInMinus.Sec/(NonNovelDat$TimeInPlus.Sec + NonNovelDat$TimeInMinus.Sec + NonNovelDat$TimeInNeutral.Sec)
NonNovelDat$PropTimeNeutral.all <- NonNovelDat$TimeInNeutral.Sec/(NonNovelDat$TimeInPlus.Sec + NonNovelDat$TimeInMinus.Sec + NonNovelDat$TimeInNeutral.Sec)

#For Novel stim dataset
NovelDat$PropTimePlus <- NovelDat$TimeInPlus.Sec/(NovelDat$TimeInPlus.Sec + NovelDat$TimeInMinus.Sec)
NovelDat$PropTimeMinus <- NovelDat$TimeInMinus.Sec/(NovelDat$TimeInPlus.Sec + NovelDat$TimeInMinus.Sec)

NovelDat$PropTimePlus.all <- NovelDat$TimeInPlus.Sec/(NovelDat$TimeInPlus.Sec + NovelDat$TimeInMinus.Sec + NovelDat$TimeInNeutral.Sec)
NovelDat$PropTimeMinus.all <- NovelDat$TimeInMinus.Sec/(NovelDat$TimeInPlus.Sec + NovelDat$TimeInMinus.Sec + NovelDat$TimeInNeutral.Sec)
NovelDat$PropTimeNeutral.all <- NovelDat$TimeInNeutral.Sec/(NovelDat$TimeInPlus.Sec + NovelDat$TimeInMinus.Sec + NovelDat$TimeInNeutral.Sec)

## double-check
NovelDat$PropTimePlus
NonNovelDat$PropTimePlus

## create ggplotable dataframe to look at spread of proportions
dat1 <- data.frame(treatment='Nov', PropTimePlus=NovelDat$PropTimePlus)
dat2 <- data.frame(treatment='NoNov', PropTimePlus=NonNovelDat$PropTimePlus)
dat <- rbind(dat1,dat2)
ggplot(dat, aes(x=treatment, y=PropTimePlus)) +
  geom_violin(trim=TRUE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)

se <- function(x) sqrt(var(x) / length(x))
se(NonNovelDat$PropTimeMinus.all)

### Create a ggplottable dataframe to show the mean proportions in barplot form
prop_barplot <- data.frame(proportion=c(mean(NonNovelDat$PropTimeMinus.all), 
                                           mean(NonNovelDat$PropTimeNeutral.all),
                                           mean(NonNovelDat$PropTimePlus.all),
                                           mean(NovelDat$PropTimeMinus.all),
                                           mean(NovelDat$PropTimeNeutral.all),
                                           mean(NovelDat$PropTimePlus.all)),
                           trial_type=c('trained','trained','trained', 'novel','novel','novel'),
                           area=c('incorrect','neutral','correct'),
                           sterr=c(se(NonNovelDat$PropTimeMinus.all), 
                                   se(NonNovelDat$PropTimeNeutral.all),
                                   se(NonNovelDat$PropTimePlus.all),
                                   se(NovelDat$PropTimeMinus.all),
                                   se(NovelDat$PropTimeNeutral.all),
                                   se(NovelDat$PropTimePlus.all))
                           )
prop_barplot$trial_type <- factor(prop_barplot$trial_type, levels=c('trained','novel'))
prop_barplot$area <- factor(prop_barplot$area, levels=c('incorrect', 'neutral', 'correct'))

### Save this
time_prop_olf <- ggplot(prop_barplot, aes(fill=area, y=proportion, x=trial_type))+
  geom_bar(position="dodge", stat="identity", color="black") + 
  xlab('') +
  ylab('Proportion of trial time in area') +
  theme_classic()  + 
  scale_fill_manual(name='arena area', values=c("incorrect"='gray', 'neutral'='white', 'correct' ='grey25')) +
  theme(text=element_text(size=14), axis.text=element_text(size=13)) +
  geom_errorbar(aes(ymin=proportion-sterr, ymax=proportion+sterr), width=.2, position=position_dodge(.9)) +
  scale_y_continuous(breaks=seq(0,.6,0.2), limits=c(0,0.6))


# beta regressions now
# so we use Smithson and Verkuilen 2006 to transform the data
NonNovelDat$PropTimePlus.t <- (NonNovelDat$PropTimePlus*(nrow(NonNovelDat)-1) + 0.5)/nrow(NonNovelDat)
NovelDat$PropTimePlus.t <- (NovelDat$PropTimePlus*(nrow(NovelDat)-1) + 0.5)/nrow(NovelDat)

## model proportion of time spent on either side, ignoring Neutral area, for non-novel stim test trials.
nonnovel.timemod <- betareg(PropTimePlus.t ~ 1, data = NonNovelDat)

## Back transform to proportion
NNMean <- 1/(1 + exp(-(nonnovel.timemod$coefficients$mean)))  ## mean is 0.464
NNMean <- (NNMean*nrow(NonNovelDat)-0.5)/(nrow(NonNovelDat)-1)

# Gethe conf intervals
NNCIs <- 1/(1 + exp(-(confint(nonnovel.timemod)[c(1,3)]))) ### This overlaps 0.5, so they spend as much time in + as -, when they're not in neutral
NNCIs <- (NNCIs*nrow(NonNovelDat)-0.5)/(nrow(NonNovelDat)-1)

# now for novel data
novel.timemod <- betareg(PropTimePlus.t ~ 1, data = NovelDat) 
NovMean <- 1/(1 + exp(-novel.timemod$coefficients$mean))  ## mean is 0.59859
NovMean <- (NovMean*nrow(NovelDat)-0.5)/(nrow(NovelDat)-1)
NovCIs <- 1/(1 + exp(-(confint(novel.timemod)[c(1,3)]))) ### This also overlaps 0.5.
NovCIs <- (NovCIs*nrow(NovelDat)-0.5)/(nrow(NovelDat)-1)

# use values to create plottable data frame
toy.dat.3 = data.frame("PropTimePlus" = c(NNMean, NovMean),
                       "StimulusType" = c("trained", "novel"),
                       "UpperCI" = c(NNCIs[2], NovCIs[2]),
                       "LowerCI" = c(NNCIs[1], NovCIs[1]))
toy.dat.3$StimulusType <- factor(toy.dat.3$StimulusType, levels=c('trained', 'novel'))

###### Save figure panel for publication
modeled_values_olf <- ggplot(data = toy.dat.3, aes(x = StimulusType, y = PropTimePlus)) + 
  ylab("Proportion of time on correct side") + xlab("Test Stimuli") +
  geom_point() +
  geom_errorbar(aes(x = StimulusType, ymin = LowerCI, ymax = UpperCI, width = 0.25)) +
  ylim(0,1) + 
  theme_classic() +
  theme(text=element_text(size=14), axis.text=element_text(size=13)) +
  geom_hline(yintercept = 0.50, lty = 2)

png(here::here('output/04.olfactorytests.png'), width=300, height=800)
cowplot::plot_grid(choice_ratios_olf, time_prop_olf, modeled_values_olf, 
          ncol=1, nrow=3, labels=c("A", "B", "C"), vjust=1, 
          align='v', axis='lr')
dev.off()


# Analysis of training trials ---------------------------------------------
# two trials per day for each individual, 7 days of training (5 initially, 2 more after 1st test trials)
# I've combined the trial days/trial no. on each day (1 or 2) into a single "TrialNo" variable 
# Each individual scored for time in the "+" side of arena, and also for whether or not they went into
# the correct shelter ("achieved success").
# Those that achieved success also have a value for time to success.

tt.dat = read.csv(here::here("data/Amblypygid_TrainTrials_FINAL_20211120.csv"))
names(tt.dat)
str(tt.dat)
unique(tt.dat$ID)

## lets first look at whether the likelihood of achieving success increased over trials
## need to run a mixed effects glm where we control for individual ID
## First, have to convert Achieved.Success.Y.N. to 0s and 1s
str(tt.dat$Achieved.Success.Y.N.)
tt.dat$Achieved.Success.Y.N.=revalue(tt.dat$Achieved.Success.Y.N., c("Y"=1,"N"=0))
str(tt.dat$Achieved.Success.Y.N.)
# now make sure it's numeric...
tt.dat$Achieved.Success.Y.N. = as.integer(tt.dat$Achieved.Success.Y.N.)
str(tt.dat$Achieved.Success.Y.N. )
tt.dat$ID <- factor(tt.dat$ID, levels= c('Ind 1', 'Ind 2', 'Ind 3', 'Ind 4', 'Ind 5',
                                                   'Ind 6', 'Ind 7', 'Ind 8', 'Ind 9', 'Ind 10',
                                                   'Ind 11', 'Ind 12', 'Ind 13', 'Ind 14'))

# Create and compare models
success.mod = glmer(Achieved.Success.Y.N. ~ TrialNo + (1|ID), data = tt.dat, family = "binomial")
anova(success.mod,
      update(success.mod, .~.-TrialNo), test = "LRT")
summary(success.mod)

# Print model info to nice table 
sjPlot::tab_model(success.mod, show.se = T, show.ci = F, show.re.var = F, show.intercept = T, 
                  show.icc = FALSE, show.ngroups = F, show.r2 = F,
                  dv.labels = c("Success (Y/N)"), 
                  string.se = "SE", transform = NULL, file = here::here("output/2.table_succ1.doc"))

# plot the model
# need to create a smoothing function for binomial response
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

# by individual, with mean included
B_mod.succYN <- ggplot(tt.dat, aes(x = TrialNo, y = Achieved.Success.Y.N., color = ID)) + 
  binomial_smooth(se=F) + binomial_smooth(color = "black", se=T)+
  geom_point() + geom_jitter(width = 0.25, height = 0.05)+
  ylab("Achieved Success Y/N") + xlab("") + 
  scale_x_continuous(breaks=seq(1,14,1)) +
  theme_classic() + 
  theme(text=element_text(size=14), axis.text=element_text(size=13))


# Model time to enter correct door in training trials ----------------------------------------
## let's see if they get quicker at doing so! ##
str(tt.dat$Time.to.Success.)
# this is unfortunately in mm:ss format, so need to use the lubridate stuff from above to get it into seconds
tt.dat$Time.to.Success.Sec <- ms(tt.dat$Time.to.Success.) ## step one, converts to a weird format
tt.dat$Time.to.Success.Sec <- as.numeric(tt.dat$Time.to.Success.Sec) ## step two, apply as.numeric() function
tt.dat$Time.to.Success.Sec

tt.dat.succ$Time.to.Success.Sec = ms(tt.dat.succ$Time.to.Success.) ## step one, converts to a weird format
tt.dat.succ$Time.to.Success.Sec = as.numeric(tt.dat.succ$Time.to.Success.Sec) ## step two, apply as.numeric() function
tt.dat.succ$Time.to.Success.Sec  ## double check. Now it's in seconds! Woo!

#### Instead of a linear model, we are going to use a survival analysis.
# create time to success and censoring data
tt.dat$Time.to.Success.Sec[is.na(tt.dat$Time.to.Success.Se)] <- 1201
tt.dat$Event <- 1  ## Indicates success was achieved before 20 minutes
tt.dat$Event[tt.dat$Time.to.Success.Sec>1200] <- 0

# Estimate survival curves with Kaplan-Meier method
tsucc <- survfit(Surv(Time.to.Success.Sec, Event) ~ TrialNo, data=tt.dat)
names(tsucc)
plot(tsucc)

# Cox regression
tsuccC <- coxph(Surv(Time.to.Success.Sec, Event) ~ TrialNo, data=tt.dat)

# Make a plot that is understandable 
tt.datpretty <- dplyr::mutate(tt.dat, Trial=ifelse((TrialNo<8), '1 to 7', '8 to 14'),
                       Trial=factor(Trial))
ttsuccP <- survfit(Surv(Time.to.Success.Sec, Event) ~Trial, data=tt.datpretty)

timetosuccess <- autoplot(ttsuccP, censor=FALSE)
D_modtimetosucc_olf <- timetosuccess + 
  theme_classic() + 
  theme(text=element_text(size=14), axis.text=element_text(size=13))  +
  ylab('% not entered shelter') +
  xlab('time (sec)') + 
  labs(fill='trials', color='trials') +
  coord_cartesian(xlim=c(0,1150))

#plots combined
not_used <- ggplot(tt.dat.succ, aes(x = TrialNo, y = log(Time.to.Success.Sec), color = ID)) + 
  geom_smooth(method = "lm", se = F) + geom_smooth(method = "lm", color = "black")+
  geom_point() +
  ylab("log(Time to Success)") + xlab("Training Trial Number") + 
  theme_classic() + 
  theme(text=element_text(size=14), axis.text=element_text(size=13))



# Training Trials,  proportion of time spent on the plus vs minus  --------
#### Now consider the amount of time spent on the + side of the arena, regardless of whether they eventually went into the + refuge ####
names(tt.dat)
str(tt.dat$Time.in...1)
# convert this time to seconds using the function in lubridate package
tt.dat$Time.in...1.Sec = lubridate::ms(tt.dat$Time.in...1) ## step one, converts to a weird format
tt.dat$Time.in...1.Sec = as.numeric(tt.dat$Time.in...1.Sec) ## step two, apply as.numeric() function
tt.dat$Time.in...1.Sec  ## double check. Now it's in seconds! Woo!

tt.dat$Time.in..Sec = lubridate::ms(tt.dat$Time.in..) ## step one, converts to a weird format
tt.dat$Time.in..Sec = as.numeric(tt.dat$Time.in..Sec) ## step two, apply as.numeric() function
tt.dat$Time.in..Sec  ## double check. Now it's in seconds! Woo!

tt.dat$Prop.in.Plus = tt.dat$Time.in...1.Sec/(tt.dat$Time.in...1.Sec + tt.dat$Time.in..Sec)
tt.dat$Prop.in.Plus

# use a beta regression to model this proportion as a function of trial number
# Need to transform first.
tt.dat$Prop.in.Plus.t <- (tt.dat$Prop.in.Plus*(nrow(tt.dat)-1) + 0.5)/nrow(tt.dat)
tt.dat$ID <- factor(tt.dat$ID, levels= c('Ind 1', 'Ind 2', 'Ind 3', 'Ind 4', 'Ind 5',
                                         'Ind 6', 'Ind 7', 'Ind 8', 'Ind 9', 'Ind 10',
                                         'Ind 11', 'Ind 12', 'Ind 13', 'Ind 14'))

prop.mod <- glmmTMB(Prop.in.Plus.t ~ TrialNo + (1|ID), tt.dat, family=beta_family)
prop.mod2 <- glmmTMB(Prop.in.Plus.t ~  (1|ID), tt.dat, family=beta_family)
anova(prop.mod, prop.mod2, test = "LRT")
summary(prop.mod)

sjPlot::tab_model(prop.mod, show.se = T, show.ci = F, show.re.var = F, show.intercept = T, 
                  show.icc = FALSE, show.ngroups = F, show.r2 = F,
                  dv.labels = c("Proportion of time spent on correct side"), 
                  string.se = "SE", transform = NULL, file = here::here("output/3.table_proptime.doc"))

C_proptimeplus <- ggplot(tt.dat, aes(x = TrialNo, y = Prop.in.Plus, color = ID)) + binomial_smooth(se = F) + 
  geom_point() + geom_jitter(width = 0.2, height = 0) + binomial_smooth(color = "black")+
  ylab("Proportion of Time in +") + xlab("Training Trial Number") + 
  theme_classic() + 
  scale_x_continuous(breaks=seq(1,14,1)) +
  theme(text=element_text(size=14), axis.text=element_text(size=13))


# Does training trial affect what side of the arena they go to 1st? --------
# get data in correct formats
tt.dat$X1st.choice
str(tt.dat$X1st.choice)
tt.dat$X1st.choice = as.numeric(as.factor(tt.dat$X1st.choice))-1
str(tt.dat$X1st.choice)

# make model
first.mod = glmer(X1st.choice ~ TrialNo + (1|ID), data =  tt.dat, family = "binomial")
summary(first.mod)
anova(first.mod,
      update(first.mod, .~. - TrialNo), test = "LRT")

# export and save model info to table
sjPlot::tab_model(first.mod, show.se = T, show.ci = F, show.re.var = F, show.intercept = T, 
                  show.icc = FALSE, show.ngroups = F, show.r2 = F,
                  dv.labels = c("Correct first choice (Y/N)"), 
                  string.se = "SE", transform = NULL, file = here::here("output/1.table_firstchoice.doc"))


#plot 
A_modfirstchoice_olf <- ggplot(tt.dat, aes(x = TrialNo, y = X1st.choice, color = ID)) + 
  geom_point() + geom_jitter(width = 0.2, height = 0)+
  binomial_smooth(se = F) + binomial_smooth(color = "black")+
  ylab("Correct first choice Y/N") + xlab("") + 
  scale_x_continuous(breaks=seq(1,14,1)) +
  theme_classic() + 
  theme(text=element_text(size=14), axis.text=element_text(size=13))

# save the legend separately
mod_leg <- get_legend(A_modfirstchoice_olf)

# compine plots for publication
modelplots_olf <- cowplot::plot_grid(A_modfirstchoice_olf + theme(legend.position = 'none'),
                                     B_mod.succYN + theme(legend.position = 'none'),
                                     C_proptimeplus + theme(legend.position = 'none'),
                                     axis = 'lrbt',
                                     align = 'vh',
                                     labels=c('AUTO'),
                                     nrow=3, ncol=1)


D_modtimetosucc_olf

# Print plots for publication
png(here::here('output/02.modeledbytrial_olf.png'), height=900, width=400)
cowplot::plot_grid(modelplots_olf, mod_leg,
          ncol=2, rel_widths=c(1,.2))
dev.off()

png(here::here('output/03.modeledsurival.png'), height=300, width=400)
D_modtimetosucc_olf
dev.off()


