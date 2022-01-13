#######################################################################################
##                          Amblypygid Cognition Analysis                            ##        
##                            2. Experiment 2                                        ##
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
library(scales)
library(grid)
library(glue)
library(rlang)
library(sf)
library(png)
library(ggpattern)
library(sjPlot)
library(ggfortify)


# Setup -------------------------------------------------------------------

## load data sets 
tests = read.csv(here::here('data/MariahDataTests.csv')) 
trials = read.csv(here::here("data/MariahDataTrials.csv"), blank.lines.skip = TRUE) 
trials <- trials[!is.na(trials$ID),]

#### Let's look at the test trials and separate out the tests into days for analysis and figures
day1 <- as.data.frame(table(tests$FirstChoice[tests$Day==1], tests$SameORdiff[tests$Day==1]))
day2 <- as.data.frame(table(tests$FirstChoice[tests$Day==2], tests$SameORdiff[tests$Day==2]))
day3 <- as.data.frame(table(tests$FirstChoice[tests$Day==3], tests$SameORdiff[tests$Day==3]))
day3 <- day3[day3$Var1!='none',]
alldays <- as.data.frame(table(tests$FirstChoice, tests$SameORdiff))
alldays <- alldays[alldays$Var1!='none',]

## Chi squared tests of the days, excluding those tests where the individual didn't choose.
chisq.test(table(tests$FirstChoice[tests$Day==1]))
chisq.test(table(tests$FirstChoice[tests$Day==2]))
chisq.test(table(tests$FirstChoice[tests$Day==3 & tests$FirstChoice!='none']))
chisq.test(table(tests$FirstChoice[tests$FirstChoice!='none']))

# Start creating plotable data
day1$day <- 'test 1'
day2$day <- 'test 2'
day3$day <- 'novel test'
alldays$day <- 'combined'

# combine plotable data
testsby_day <- rbind(day1,day2,day3,alldays)
testsby_day <- dplyr::rename(testsby_day, choice=Var1, treatment=Var2)
testsby_day$day <- factor(testsby_day$day, levels=c('test 1', 'test 2', 'novel test', 'combined'))
testsby_day$choice <- factor(testsby_day$choice, levels=c('incorrect', 'correct'))
testsby_day$choice[is.na(testsby_day$choice)] <- 'incorrect'
testsby_day$treatment <- factor(testsby_day$treatment, levels=c("same", 'dif'))

# make plot
bars_freq_days <- ggplot(testsby_day, aes(fill=choice, y=Freq, x=choice, pattern=treatment)) +
  geom_bar_pattern(position='stack', stat="identity", color='black',
                   pattern_fill='white',
                   pattern_color='white',
                   pattern_alpha=1,
                   pattern_size=.4,
                   pattern_spacing=0.15,
                   pattern_key_scale_factor=.1) + 
  facet_wrap(~day, nrow=1, strip.position = 'top') +
  theme(panel.spacing = unit(0, "cm")) +
  ylab('Number of trials') +
  xlab('') +
  theme_classic()  + 
  scale_fill_manual(values=c('incorrect'='gray', 'correct' ='grey25'), 
                    guide=guide_legend(override.aes = list(pattern='none'))) +
  scale_pattern_manual(values=c('dif'='stripe', 'same'='none'), labels=c('same','different')) +
  scale_y_continuous(expand=c(0,0), limits= c(0,15), breaks=seq(0,14,2)) +
  theme(text=element_text(size=14), axis.text=element_text(size=13), 
        axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

# create data and plots where tests are separated by treatment 
# (i.e. whether they were trained to  learn"same" or "different" concepts)
testsSame <- subset(tests, SameORdiff=='same')
testsDiff <- subset(tests, SameORdiff=='dif')

same <- as.data.frame(table(testsSame$FirstChoice))
same <- same[same$Var1!='none',]
same$treatment <- 'same'
diff <- as.data.frame(table(testsDiff$FirstChoice))
diff <- diff[diff$Var1!='none',]
diff$treatment <- 'diff'

# chisquared tests for these, first combining the first two days of testing
chisq.test(table(testsSame$FirstChoice[testsSame$Day!=3]))
chisq.test(table(testsDiff$FirstChoice[testsDiff$Day!=3]))
# checking the final day of testing for a difference
chisq.test(table(testsSame$FirstChoice[testsSame$Day==3]))
chisq.test(table(testsDiff$FirstChoice[testsDiff$Day==3]))
# checking the different treatments for a difference between correct and incorrect choices
chisq.test(same$Freq)
chisq.test(diff$Freq)

# make something plottable
samediff <- rbind(same,diff)
samediff <- dplyr::rename(samediff, choice=Var1, Frequency=Freq)
samediff$treatment <- factor(samediff$treatment, levels=c('same','different'))
samediff$treatment[is.na(samediff$treatment)] <- 'different'
samediff$choice <- factor(samediff$choice, levels=c('incorrect', 'correct'))
samediff$choice[is.na(samediff$choice)] <- 'incorrect'

# plot
bars_freq_SD <- ggplot(samediff, aes(fill=choice, y=Frequency, x=treatment, pattern=treatment)) +
  geom_bar_pattern(position="dodge", stat="identity", color='grey25',
                   pattern_fill='white',
                   pattern_color='white',
                   pattern_alpha=1,
                   pattern_size=.4,
                   pattern_spacing=0.08,
                   pattern_key_scale_factor=.1) + 
  xlab('') +
  ylab('') +
  theme_classic()  + 
  scale_pattern_manual(values=c('different'='stripe', 'same'='none')) +
  scale_fill_manual(values=c('incorrect'='gray', 'correct' ='grey25')) +
  scale_y_continuous(expand=c(0,0), limits= c(0,15)) +
  theme(text=element_text(size=14), axis.text=element_text(size=13), 
        axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
bars_freq_SD

# save the legend separately 
freq_leg <- get_legend(bars_freq_days)

##### Plot the these babies
freq_bars <- cowplot::plot_grid(bars_freq_days + theme(legend.position = 'none'),
                       bars_freq_SD + theme(legend.position = 'none'),
                       freq_leg,
                       align='h', axis='bt',
                       nrow=1, labels=c('A','B'), rel_widths = c(4,2,1))


##### Let's look at the time spent on either side
# calculate proportions
tests$PropPlusTime <- tests$PlusTime/(tests$PlusTime + tests$MinusTime)
tests$PropMinusTime <- tests$MinusTime/(tests$PlusTime + tests$MinusTime)
tests$PropPlusTimeAll <- tests$PlusTime/900
tests$PropMinusTimeAll <- tests$MinusTime/900
tests$PropCenterTimeAll <- tests$CenterTime/900

df1 <- data.frame(time='plus', prop=tests$PropPlusTimeAll)
df2 <- data.frame(time='center', prop=tests$PropCenterTimeAll)
df3 <- data.frame(time='minus', prop=tests$PropMinusTimeAll)
proptimes <- rbind(df1,df2,df3)


## reorganize data to the get the means for all these proportions
se <- function(x) sqrt(var(x) / length(x))
tests$Day <- as.factor(tests$Day)
proptime_byday <- tests %>%
  dplyr::group_by(Day) %>%
  dplyr::summarise(plusprop=mean(PropPlusTimeAll), cenprop= mean(PropCenterTimeAll), minprop=mean(PropMinusTimeAll))

# also get the standard errors so we can add error bars
setime_byday <- tests %>%
  dplyr::group_by(Day) %>%
  dplyr::summarise(plusse=se(PropPlusTimeAll), cense= se(PropCenterTimeAll), minse=se(PropMinusTimeAll))


proptime_byday$Day <- factor(proptime_byday$Day, levels=c('1','2','3','combined'))
proptime_byday <- rbind(proptime_byday, c('combined', mean(tests$PropPlusTimeAll), mean(tests$PropCenterTimeAll),
                                          mean(tests$PropMinusTimeAll)))

setime_byday$Day <- factor(setime_byday$Day, levels=c('1','2','3','combined'))
setime_byday <- rbind(setime_byday, c('combined', se(tests$PropPlusTimeAll), se(tests$PropCenterTimeAll),
                                          se(tests$PropMinusTimeAll)))

proptime_byday <- pivot_longer(proptime_byday, cols=c('plusprop','cenprop','minprop'), names_to = 'proportiontype',
                               values_to='proportion')
setime_byday <- pivot_longer(setime_byday, cols=c('plusse','cense','minse'), names_to='setype',
                             values_to='stderr')
proptime_byday$sterr <- setime_byday$stderr

proptime_byexp <- tests %>%
  dplyr::group_by(SameORdiff) %>%
  dplyr::summarise(plusprop=mean(PropPlusTimeAll), cenprop= mean(PropCenterTimeAll), minprop=mean(PropMinusTimeAll))

setime_byexp <- tests %>%
  dplyr::group_by(SameORdiff) %>%
  dplyr::summarise(plusse=se(PropPlusTimeAll), cense= se(PropCenterTimeAll), minse=se(PropMinusTimeAll))

proptime_byexp <- pivot_longer(proptime_byexp, cols=c('plusprop','cenprop','minprop'), names_to = 'proportiontype',
                               values_to='proportion')

setime_byexp <- pivot_longer(setime_byexp, cols=c('plusse','cense','minse'), names_to = 'proportiontype',
                               values_to='sterr')

proptime_byexp$sterr <- setime_byexp$sterr

# cleanup for plotting 
proptime_byday$proportiontype <- factor(proptime_byday$proportiontype, levels=c('minprop', 'cenprop','plusprop'))
proptime_byday$proportion <- as.numeric(proptime_byday$proportion)
proptime_byexp$proportiontype <- factor(proptime_byexp$proportiontype, levels=c('minprop', 'cenprop','plusprop'))
proptime_byexp$SameORdiff <- factor(proptime_byexp$SameORdiff, levels=c('same','different'))
proptime_byexp$SameORdiff[is.na(proptime_byexp$SameORdiff)] <- 'different'

proptime_byday$proportion <- as.numeric(proptime_byday$proportion)
proptime_byday$sterr <- as.numeric(proptime_byday$sterr)

# plot
bars_time_days <- ggplot(proptime_byday, aes(fill=proportiontype, y=proportion, x=Day)) +
  geom_bar(position='dodge', stat='identity', color='black') + 
  geom_errorbar(aes(ymin=proportion-sterr, ymax=proportion+sterr), 
                width=.2, position=position_dodge(.9)) +
  ylab('Proportion of time in area') +
  xlab('') +
  theme_classic()  + 
  scale_fill_manual("area", values=c("minprop"='gray', "cenprop"='white', "plusprop"='grey25')) +
  scale_y_continuous(expand=c(0,0), limits= c(0,0.7)) +
  theme(text=element_text(size=14), axis.text=element_text(size=13), 
        axis.text.x=element_blank(), axis.ticks.x=element_blank()) + theme(legend.position = 'none')

bars_time_SD <- ggplot(proptime_byexp, aes(fill=proportiontype, y=proportion, x=SameORdiff)) +
  geom_bar(position='dodge', stat='identity', color="black") + 
  geom_errorbar(aes(ymin=proportion-sterr, ymax=proportion+sterr), 
                width=.2, position=position_dodge(.9)) +
  ylab('') +
  xlab('') +
  theme_classic()  + 
  scale_fill_manual("arena area", values=c("minprop"='gray', "cenprop"='white', "plusprop"='grey25'),
                    labels=c('incorrect', 'neutral', 'correct')) +
  scale_y_continuous(expand=c(0,0), limits= c(0,0.7)) +
  theme(text=element_text(size=14), axis.text=element_text(size=13), 
        axis.text.x=element_blank(), axis.ticks.x=element_blank())

#save the legend separately
time_leg <- get_legend(bars_time_SD)

#save the plot 
time_bars <- cowplot::plot_grid(bars_time_days + theme(legend.position = 'none'),
                       bars_time_SD + theme(legend.position = 'none'),
                       time_leg,
                       align='h', axis='bt',
                       nrow=1, labels=c('c','d'), rel_widths = c(4,2,1))
  

# Model proportion of time in front of correct door -----------------------
# beta regressions now and get info from them to plot the calculated coefficients
# and SEs to see if they're significantly different from 0.5
# so we use Smithson and Verkuilen 2006 to transform the data
tests$PropTimePlus.t <- (tests$PropPlusTime*(nrow(tests)-1) + 0.5)/nrow(tests)
tests$backttest <- (tests$PropTimePlus.t*nrow(tests)-0.5)/(nrow(tests)-1)
test.timepropmod.all <- betareg(PropTimePlus.t ~ 1, data = tests)
test.timepropmod.1 <- betareg(PropTimePlus.t ~ 1, data = tests[tests$Day==1,])
test.timepropmod.2 <- betareg(PropTimePlus.t ~ 1, data = tests[tests$Day==2,])
test.timepropmod.3 <- betareg(PropTimePlus.t ~ 1, data = tests[tests$Day==3,])
test.timepropmod.dif <- betareg(PropTimePlus.t ~ 1, data = tests[tests$SameORdiff=='dif',])
test.timepropmod.same <- betareg(PropTimePlus.t ~ 1, data = tests[tests$SameORdiff=='same',])

## Back transform to proportion then undo other transform
betamean.all <- 1/(1 + exp(-(test.timepropmod.all$coefficients$mean)))  
betamean.all <- (betamean.all*nrow(tests)-0.5)/(nrow(tests)-1)
# Get the conf intervals
betaCIs.all <- 1/(1 + exp(-(confint(test.timepropmod.all)[c(1,3)]))) 
betaCIs.all <- (betaCIs.all*nrow(tests)-0.5)/(nrow(tests)-1)

## Back transform to proportion then undo other transform
betamean.1 <- 1/(1 + exp(-(test.timepropmod.1$coefficients$mean)))  
betamean.1 <- (betamean.1*nrow(tests)-0.5)/(nrow(tests)-1)
# Get the conf intervals
betaCIs.1 <- 1/(1 + exp(-(confint(test.timepropmod.1)[c(1,3)]))) 
betaCIs.1 <- (betaCIs.1*nrow(tests)-0.5)/(nrow(tests)-1)

## Back transform to proportion then undo other transform
betamean.2 <- 1/(1 + exp(-(test.timepropmod.2$coefficients$mean)))  
betamean.2 <- (betamean.2*nrow(tests)-0.5)/(nrow(tests)-1)
# Get the conf intervals
betaCIs.2 <- 1/(1 + exp(-(confint(test.timepropmod.2)[c(1,3)]))) 
betaCIs.2 <- (betaCIs.2*nrow(tests)-0.5)/(nrow(tests)-1)

## Back transform to proportion then undo other transform
betamean.3 <- 1/(1 + exp(-(test.timepropmod.3$coefficients$mean)))  
betamean.3 <- (betamean.3*nrow(tests)-0.5)/(nrow(tests)-1)
# Get the conf intervals
betaCIs.3 <- 1/(1 + exp(-(confint(test.timepropmod.3)[c(1,3)]))) 
betaCIs.3 <- (betaCIs.3*nrow(tests)-0.5)/(nrow(tests)-1)

## Back transform to proportion then undo other transform
betamean.dif <- 1/(1 + exp(-(test.timepropmod.dif$coefficients$mean)))  
betamean.dif <- (betamean.dif*nrow(tests)-0.5)/(nrow(tests)-1)
# Get the conf intervals
betaCIs.dif <- 1/(1 + exp(-(confint(test.timepropmod.dif)[c(1,3)]))) 
betaCIs.dif <- (betaCIs.dif*nrow(tests)-0.5)/(nrow(tests)-1)

## Back transform to proportion then undo other transform
betamean.same <- 1/(1 + exp(-(test.timepropmod.same$coefficients$mean)))  
betamean.same <- (betamean.same*nrow(tests)-0.5)/(nrow(tests)-1)
# Get the conf intervals
betaCIs.same <- 1/(1 + exp(-(confint(test.timepropmod.same)[c(1,3)]))) 
betaCIs.same <- (betaCIs.same*nrow(tests)-0.5)/(nrow(tests)-1)

# Combine in a dataframe for plotting
modtimeprop = data.frame('grouping' = c('combined', 'test 1', 'test 2', 'novel test', 'same', 'different'),
                        "propmean" = c(betamean.all, betamean.1, betamean.2, betamean.3, betamean.same, betamean.dif),
                       "UpperCI" = c(betaCIs.all[2], betaCIs.1[2],betaCIs.2[2],betaCIs.3[2],betaCIs.same[2],betaCIs.dif[2]),
                       "LowerCI" = c(betaCIs.all[1], betaCIs.1[1],betaCIs.2[1],betaCIs.3[1],betaCIs.same[1],betaCIs.dif[1]))

modtimeprop$grouping <- factor(modtimeprop$grouping, levels=c('test 1', 'test 2', 'novel test', 'combined', 'same', 'different'))

modtimeprop_days <- modtimeprop[1:4,]
modtimeprop_exp <- modtimeprop[5:6,]

# plot results by day
CImodprop_days <- ggplot(data = modtimeprop_days, aes(x = grouping, y = propmean)) + 
  ylab("Proportion of time on correct side  ") + xlab("") +
  geom_point() +
  geom_errorbar(aes(x = grouping, ymin = LowerCI, ymax = UpperCI, width = 0.25)) +
  ylim(0,1) + 
  theme_classic() +
  theme(text=element_text(size=14), axis.text=element_text(size=13)) +
  geom_hline(yintercept = 0.50, lty = 2)

# plot results by treatment (same/different)
CImodprop_exp <- ggplot(data = modtimeprop_exp, aes(x = grouping, y = propmean)) + 
  ylab('') + xlab("") +
  geom_point() +
  geom_errorbar(aes(x = grouping, ymin = LowerCI, ymax = UpperCI, width = 0.25)) +
  ylim(0,1) + 
  theme_classic() +
  theme(text=element_text(size=14), axis.text=element_text(size=13)) +
  geom_hline(yintercept = 0.50, lty = 2)

CIsplots <- cowplot::plot_grid(CImodprop_days, CImodprop_exp,
          align='h', axis='bt',
          nrow=1, labels=c('e','f'), rel_widths = c(4,2,1))

box <- rectGrob(x=.77, y=.99, width=.2, height=.01, gp=gpar(fill='white', alpha=0))

# Print the plot of everything for publication
png(here::here('output/06.texturetests.png'), width=610, height=800)
cowplot::plot_grid(bars_freq_days + theme(legend.position = 'none'), bars_freq_SD + theme(legend.position = 'none'), freq_leg,
          bars_time_days, bars_time_SD + theme(legend.position = 'none'), time_leg,
          CImodprop_days, CImodprop_exp,
          nrow=3, ncol=3, rel_widths = c(4,2.1,1),
          rel_heights = c(1,1,1),
          axis = 'lrbt',
          align = 'vh',
          labels=c('A', 'B', '', 'C','D','', 'E', 'F',''), hjust=-2)

dev.off()



# TRAINING TRIALS ---------------------------------------------------------

#### Looking at the training trials
trains <- trials[trials$Trial != 'TEST',]
trains <- trains[trains$Trial != 'TEST NOVEL',]

### create a variable that is "success"
trains$TimeToEnter <- as.numeric(trains$TimeToEnter)
trains$success <- NA
trains$success[trains$TimeToEnter>900] <- 0
trains$success[is.na(trains$success)] <- 1
trains$TrialNo <- as.integer(trains$TrialNo)
trains$ID <- as.factor(trains$ID)

success.mod = glmer(success ~ TrialNo + (1|ID), data = trains, family = "binomial")
anova(success.mod,
      update(success.mod, .~.-TrialNo), test = "LRT")
summary(success.mod)

# export and save info about the model of success liklihood
sjPlot::tab_model(success.mod, show.se = T, show.ci = F, show.re.var = F, show.intercept = T, 
                  show.icc = FALSE, show.ngroups = F, show.r2 = F,
                  dv.labels = c("log(Time to success in seconds)"), 
                  string.se = "SE", transform = NULL, file = here::here("output/5.table_succ2.docx"))

# plot it
# need to create a smoothing function for binomial response
binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

# now by individual, with mean included
mod.succes <- ggplot(trains, aes(x = TrialNo, y = success, color = ID)) + 
  binomial_smooth(se=F) + binomial_smooth(color = "black", se=T)+
  geom_point() + geom_jitter(width = 0.25, height = 0.05)+
  ylab("Achieved Success Y/N") + xlab("Training trial number") + 
  theme_classic() + 
  xlim(1,12) + 
  scale_x_continuous(breaks=seq(1,12,1)) +
  theme(text=element_text(size=14), axis.text=element_text(size=13),
        plot.margin = unit(c(.5,.5,0,.5), 'cm'))


### See if time to success decreases over time
# now fit a lme model of the time to success
trains$TimeToEnter <- as.numeric(trains$TimeToEnter)
trains$TrialNo <- as.integer(trains$TrialNo)
trains$ID <- as.factor(trains$ID)
trains.succ <- subset(trains, success==1)
hist(trains$TimeToEnter) # may need a different distribution...
hist(log(trains$TimeToEnter))

#### Instead of linear model, we are going to use a survival analysis.
# create time to success and censoring data
trains$Event <- 1  ## Indicates success was achieved before 20 minutes
trains$Event[trains$TimeToEnter>900] <- 0

# Estimate survival curves with Kaplan-Meier method
tsucc <- survfit(Surv(TimeToEnter, Event) ~ TrialNo, data=trains)
tsucc1 <- survfit(Surv(TimeToEnter, Event) ~ 1, data=trains)
names(tsucc)
plot(tsucc)

# Cox regression
tsuccC <- coxph(Surv(TimeToEnter, Event) ~ TrialNo, data=trains)
tsuccC1 <- coxph(Surv(TimeToEnter, Event) ~ 1, data=trains)
anova(tsuccC, tsuccC1, test='LRT')

# plot
autoplot(tsucc)
# pretty plot
trainspretty <- dplyr::mutate(trains, Trial=ifelse((TrialNo<8), '1 to 6', '7 to 12'),
                              Trial=factor(Trial))
ttsuccP <- survfit(Surv(TimeToEnter) ~Trial, data=trainspretty)

timetosuccess <- autoplot(ttsuccP)

# save the plot
survivalplotexp2 <- timetosuccess + 
  theme_classic() + 
  theme(text=element_text(size=14), axis.text=element_text(size=13))  +
  ylab('% not entered shelter') +
  xlab('time (sec)') + 
  labs(fill='trials', color='trials') +
  coord_cartesian(xlim=c(0,900)) +
  scale_x_continuous(breaks=seq(0,900,150)) +
  geom_vline(aes(xintercept=900), linetype='dotted')

# organize the plots
mod.plots <- cowplot::plot_grid(mod.succes +  theme(legend.position = 'none'),
          mod.time.to.success + theme(legend.position = 'none'), nrow=2, labels='AUTO', align = 'vh')

mod.panel <- cowplot::plot_grid(mod.succes, survivalplotexp2, nrow=2, labels='AUTO')          

png(here::here('output/05.textmods.png'), width=400, height=650)
mod.panel
dev.off()
  
  
  
  
  
  
  
  