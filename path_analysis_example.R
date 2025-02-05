### Kormos et al. (2016) Trends and sensitivities of low streamflow extremes to discharge timing 
### and magnitude in Pacific Northwest mountain streams. Water Resources Research, 52, doi:10.1002/2015WR018125
### https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2015WR018125
################################################################################################
### Example below to reproduce path analysis results presented in Figure 2 in the paper referenced above
################################################################################################

library(EGRET) # assist with USGS flow data retrieval and 7-day moving average flow metric calculation
library(dataRetrieval) # to retrieve flow data
library(tidyverse)
library(lubridate)
library(knitr)

library(fasstr) # to calculate flow center of timing (calculates sooo many more metrics)

# path analysis
# https://rpubs.com/Agrele/SEM
# https://stats.stackexchange.com/questions/297347/estimation-of-direct-and-total-effects-with-regressions-and-sem-lavaan

library(lavaan)
# install.packages("semPlot")
# install.packages("OpenMx")
library(semPlot)
# library(lavaanPlot)
# library(OpenMx)

options(scipen=999)

################################################################################################
################################################################################################
### Get Boise River daily flow data from USGS
################################################################################################

staid <- c('13185000')
staname <- c('Boise nr Twin Springs')

boise <- NULL
basin_info <- NULL

for (i in 1:length(staid)){
  
  cat("\nProcessing ", staname[i], "\n")
  Daily <- readNWISDaily(staid[i], "00060")
  INFO<- readNWISInfo(staid[i],'00060',interactive = FALSE)
  # eList <- as.egret(INFO,Daily,NA,NA)
  Daily$site_no <- INFO$site_no[1]
  Daily$station_nm <- INFO$station_nm[1]
  Daily$staname <- staname[i]
  Daily$drainSqKm <- INFO$drainSqKm[1]
  
  boise <- rbind(boise,Daily)
  basin_info <- bind_rows(INFO,basin_info)
}

# Note there are some missing data (16 NAs) in 2025
filter(boise,is.na(Q))
boise <- filter(boise,waterYear!=2025)

ggplot(boise,aes(Date,Q)) + geom_line() + geom_smooth()

# Note that the first waterYear (2011) is incomplete
glimpse(filter(boise,waterYear==min(waterYear)))
boise <- filter(boise,waterYear>1911)

################################################################################################
################################################################################################
### Calculate final annual flow metrics used in path analysis
################################################################################################

### calculate minimum 7-day moving average flow (min7q) metric for low flow dry period (June 1 - Nov 15)
boise.min7q <- boise %>% mutate(Year = year(Date), pDate = as.Date(paste(2000,Month,day(Date),sep="-"))) %>% 
  filter(between(pDate,as.Date("2000-06-01"),as.Date("2000-11-15"))) %>% group_by(Year,site_no) %>% summarize(min7q = min(Q7)) %>% ungroup()
### calculate mean annual flow (MAF) metric for waterYear (Oct 1 - Sep 30)
boise.wy <- boise %>% mutate(Year = waterYear) %>% group_by(Year,site_no) %>% summarize(MAF = mean(Q)) %>% ungroup()
### calculate waterYear flow center of timing (COT) metric for waterYear (Oct 1 - Sep 30) - DoY_50pct_TotalQ = waterYear Julian date 50% of flow passed gauge
boise.cot <- calc_annual_flow_timing(data=boise, dates = Date, values = Q, groups = site_no, 
                                     percent_total = c(25, 33.3, 50, 75), water_year_start = 10, 
                                     start_year=1912, end_year=2024, months = 1:12, transpose = FALSE)
boise.cot <- transmute(boise.cot, site_no = site_no, Year = Year, COT = DoY_50pct_TotalQ)

### pull the metrics together for modeling
boise.summary <- left_join(boise.wy,boise.min7q, by = c("Year","site_no"))
boise.summary <- left_join(boise.summary,boise.cot, by = c("Year","site_no"))

### Notice apparent step change in COT ~corresponding to 1970s PDO shift
ggplot(boise.summary %>% gather(-Year,-site_no,key=metric,value=Value), aes(Year,Value)) +
  theme_bw() +
  geom_line() +
  geom_smooth() +
  facet_wrap(~metric,ncol=1,scales='free_y')

################################################################################################
################################################################################################
### Perform path analysis using sem function in lavaan package
################################################################################################

### set up path model - based on the YouTube video tutorial this is the way to specify the indirect and total effects...
### https://www.youtube.com/watch?v=ab3Y-xnXigA R tutorial: Structural equation modelling part 2 (indirect effects)
### however the results don't seem to square with the example in the paper or they might but it is still unclear what to add to 
### the direct effect of 0.78 between MAF ~ min7q to get the net effect...should be equal to 0.07
mod<-'
COT ~ a*MAF
min7q ~ b*COT + c*MAF
# indirect and total effects
indirect:= a*b
total:= c + (a*b)
'

### the path model execution
mod.res <- sem(mod,data=boise.summary %>% filter(between(Year,1948,2013)))
# mod.res <- sem(mod,data=boise.summary %>% filter(between(Year,1948,2013)),se='bootstrap',bootstrap = 10000)

################################################################################################
### evaluate model and extract relevant statistics - rho (correlation coefficient), beta (direct effect), and total (net) effect
summary(mod.res, fit.measures=T,standardized=T,rsq=T,ci=T)

parameterestimates(mod.res) %>% 
  kable()

# parameterestimates(mod.res, boot.ci.type = "bca.simple", standardized = TRUE) %>% 
#  kable()

################################################################################################
### Still not sure how to calculate 'net effect'
# fitMeasures(mod.res)
### get rho and beta - Kormos' correlation coefficient and standardized regression coefficient (direct effect), respectively...I think
### the results correspond almost exactly to results in Figure 2
R <- lavInspect(mod.res, what = 'cor.ov') # rho
R <- data.frame(Stat = "ρ", COT_min7q = round(R["COT","min7q"],2), MAF_COT = round(R["MAF","COT"],2), MAF_min7q = round(R["MAF","min7q"],2))
B = lavInspect(mod.res, what = 'std')$beta
B <- data.frame(Stat = "ß", COT_min7q = round(B["min7q","COT"],2), MAF_COT = round(B["COT","MAF"],2), MAF_min7q = round(B["min7q","MAF"],2))
NE <- NULL # ? Total (Net) Effect

################################################################################################
### simple path model plot....illustrates that min7q affected much more strongly by mean annual flow (MAF) vs flow center of timing (COT)
### would like to try and see how close to a publication ready graph this could be made into
semPaths(mod.res, what='std', 
         edge.label.cex=1.25, curvePivot = TRUE, 
         fade=FALSE)
