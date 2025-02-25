# Path analysis of PSP Vital Sign low flow indicator

library(EGRET)
library(dataRetrieval)
library(tidyverse)
library(lubridate)
library(knitr)

library(mapview)
library(sf)

library(ggh4x)
library(plotly)
library(ggrepel)

library(fasstr) # to calculate flow center of timing (calculates sooo many more metrics)

# path analysis
# https://rpubs.com/Agrele/SEM
# https://stats.stackexchange.com/questions/297347/estimation-of-direct-and-total-effects-with-regressions-and-sem-lavaan

library(lavaan)
# install.packages("semPlot")
# install.packages("OpenMx")
library(semPlot)
library(lavaanPlot)
# library(OpenMx)

options(scipen=8,digits=7)

# https://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-r
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}


######################################
### get daily flow data
######################################

staid <- c('12134500','12167000','12149000','12048000')
# staname <- c('Skykomish','NF Stillaguamish','Snoqualmie nr Carnation','Dungeness nr Sequim')
staname <- c('Skykomish','NF Stillaguamish','Snoqualmie nr Carnation','Dungeness')

staname <- c('Snoqualmie nr Carnation','NF Snoqualmie','Sauk abv White Chuck','NF Stillaguamish',
             'Puyallup at Puyallup','Sauk nr Sauk','Naselle','NF Skokomish',
             'Satsop','Skykomish','Dungeness')
staid <- c('12149000','12142000','12186000','12167000',
           '12101500','12189500','12010000','12056500',
           '12035000','12134500','12048000')
df <- NULL
PS_drainSqKm <- 0
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
  Daily$lat <- INFO$dec_lat_va
  Daily$lon <- INFO$dec_long_va
  
  df <- rbind(df,Daily)
  basin_info <- bind_rows(INFO,basin_info)
}
###

# Note there are some substantial data gaps

###############################################################
### patch in Naselle data
###############################################################
tmp <- readNWISdv('12010000', "00060")
tmp <- renameNWISColumns(tmp)
tmp <- tmp %>% filter(Date>=as.Date("2003-09-30")) %>% transmute(site_no=site_no, 
                                                                 lat = attr(.,"siteInfo")$dec_lat_va, 
                                                                 lon = attr(.,"siteInfo")$dec_lon_va, 
                                                                 Date=Date,Q=Flow/35.3147,Month=month(Date),waterYear=ifelse(month(Date)>9,year(Date)+1,year(Date)),
                                                                 Qualifier=Flow_cd,station_nm='NASELLE RIVER NEAR NASELLE, WA',
                                                                 staname = "Naselle", Day = as.numeric(strftime(Date, 
                                                                                                                format = "%j")))
# do something different so I can extend Q7 and Q30 stats for Naselle
### fill in a few metrics that readNWISdv does not deliver
tmp2 <- filter(df,site_no == '12010000'&!Date>=as.Date("2003-09-30"))
tmp3 <- bind_rows(tmp2,tmp)

ma30 <- function(x,n=30){stats::filter(x,rep(1/n,n), sides=1)}
ma07 <- function(x,n=7){stats::filter(x,rep(1/n,n), sides=1)}
tmp3 <- mutate(tmp3, LogQ = ifelse(is.na(LogQ),log(Q),LogQ), Q7 = as.numeric(ma07(Q)), Q30 = as.numeric(ma30(Q)), drainSqKm = 141.9313)

df <- filter(df,(!site_no == '12010000'))
df <- bind_rows(df,tmp3)

###############################################################

df <- mutate(df,Year=year(Date),Datem=as.Date(paste(2000,1,1,sep='-'))+Day-1)

ggplot(df,aes(Date,Q)) +
  geom_line() +
  facet_wrap(~staname,scales="free_y")

# Test effect of changing baseline....  
# summaryQ <- df %>% filter(between(Year,1948,1978)) %>%
summaryQ <- df %>% filter(between(Year,1948,1998)) %>%
  # summaryQ <- df %>% filter(between(Year,1925,2023)) %>% 
  group_by(site_no,Day) %>%
  summarize(p75 = quantile(Q, probs = .75, na.rm = TRUE),
            p25 = quantile(Q, probs = .25, na.rm = TRUE),
            p10 = quantile(Q, probs = 0.1, na.rm = TRUE),
            p05 = quantile(Q, probs = 0.05, na.rm = TRUE),
            p00 = quantile(Q, probs = 0, na.rm = TRUE)) 

df <- left_join(df,summaryQ,by=c("site_no","Day"))

df <- mutate(df,score = ifelse(Q<p25,1,0))

summary.score <- df %>% filter(between(Datem,as.Date("2000-07-15"),as.Date("2000-09-15"))) %>% 
  group_by(site_no,station_nm,staname,lat,lon,Year) %>% 
  summarize(ndays = length(Date), score = sum(score), pdays = score/ndays*100, sQ = sum(Q,na.rm = T), mQ = mean(Q,na.rm = T)) %>% ungroup()
summary.Q <- df %>% group_by(site_no,waterYear) %>% 
  summarize(tQ = sum(Q,na.rm = T), Q = mean(Q,na.rm = T)) %>% rename(Year = waterYear) %>% ungroup()

summary.score <- left_join(summary.score,summary.Q,by=c("site_no","Year")) %>% mutate(summary.score, pQ = sQ/tQ*100)

### calculate winter mean flow metric for waterYear (Nov 1 - Mar 31)
summary.Q.winter <- df %>% filter(Month %in% c(11,12,1,2,3)) %>% group_by(site_no,waterYear) %>% 
  summarize(wQ = mean(Q,na.rm = T)) %>% rename(Year = waterYear) %>% ungroup()
summary.score <- left_join(summary.score,summary.Q.winter,by=c("site_no","Year")) 

### calculate mean annual flow (MAF) metric for waterYear (Oct 1 - Sep 30)
summary.Q.MAF <- df %>% mutate(Year = waterYear) %>% group_by(Year,site_no) %>% summarize(MAF = mean(Q)) %>% ungroup()
summary.score <- left_join(summary.score,summary.Q.MAF,by=c("site_no","Year")) 

# remove scores in years with less than complete data
summary.score <- filter(summary.score,ndays>=63)
# bin results                        
summary.score$pcat <- cut(summary.score$pdays, breaks = c(0,25,50,75,100),right=F,include.lowest = T, labels = c("<25%","25-<50%","50-<75%",">75%"))
head(summary.score,20)

summary.score <- mutate(summary.score, staname = factor(staname, levels = rev(c('Naselle','NF Skokomish','NF Stillaguamish','Satsop',
                                                                                'Puyallup at Puyallup','Skykomish','Snoqualmie nr Carnation','NF Snoqualmie',
                                                                                'Dungeness','Sauk abv White Chuck','Sauk nr Sauk'
))))

# riverPalette <- c(rep("#4059AD",4),rep("#6B9AC4",3),rep("#97D8C4",4))
riverPalette <- rev(c(rep("#4059AD",4),rep("#6B9AC4",3),rep("#363636",4)))
# The palette
# cbPalette <- c("#80CDC1","#F5F5F5","#DFC27D","#A6611A")
# cbPalette <- c("#80CDC1","#EBEBEB","#DFC27D","#A6611A")
cbPalette <- c("#80CDC1","#E6E6E6","#DFC27D","#A6611A")

numvec <- seq(1920,2025,1)
every_nth(numvec,10,inverse = T)

ggplot(summary.score %>% filter(between(Year,1920,2024)), aes(Year,staname,fill=pcat)) +
  theme_bw() +
  geom_tile(color='gray') +
  scale_fill_manual(values = cbPalette) +
  # labs(x="",y="USGS gage",title="% of days below the baseline 25th percentile") +
  # labs(x="",y="",title="% of days July 15-September 15 flow below the baseline (1948-1998) 25th percentile") +
  labs(x="",y="") +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text = element_text(size = 16), plot.title = element_text(size = 25), legend.text = element_text(size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(color = riverPalette)) +
  # panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(1920,2025,1), 
                     labels = every_nth(numvec,10,inverse = T),
                     limits = c(1920,2025), expand = c(0,0))

# https://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set
strip <- strip_themed(background_x = elem_list_rect(fill = rev(riverPalette)))

ggplot(summary.score %>% filter(between(Year,1920,2024)),aes(Year,pdays)) +
  theme_bw() +
  geom_point(color='gray') +
  geom_smooth() +
  # facet_wrap(~staname) +
  facet_wrap2(~staname, strip = strip) +
  scale_y_continuous(limits = c(0,100)) +
  labs(x="",y="%Days below the Jul 15-Sep 15 baseline 25th Percentile",title="% of days below the baseline 25th percentile") +
  theme(strip.text = element_text(color = "white"))

ggplot(summary.score %>% filter(between(Year,1920,2024)),aes(Year,mQ)) +
  theme_bw() +
  geom_point(color='gray') +
  geom_smooth() +
  # facet_wrap(~staname, scales = "free_y") +
  facet_wrap2(~staname, strip = strip, scales = 'free_y') +
  labs(x="",y="Jul 15-Sep 15 mean flow",title="Jul 15-Sep 15 mean flow") +
  theme(strip.text = element_text(color = "white"))

ggplot(summary.score %>% filter(between(Year,1920,2024)),aes(Year,Q)) +
  theme_bw() +
  geom_point(color='gray') +
  geom_smooth() +
  # facet_wrap(~staname, scales = "free_y") +
  facet_wrap2(~staname, strip = strip, scales = "free_y") +
  labs(x="",y="Water Year (Oct 1-Sep 30) mean flow") +
  theme(strip.text = element_text(color = "white")) 

summary.min7q <- df %>% mutate(pDate = as.Date(paste(2000,Month,day(Date),sep="-"))) %>% 
  filter(between(pDate,as.Date("2000-06-01"),as.Date("2000-11-15"))) %>% group_by(Year,staname,site_no) %>% summarize(min7q = min(Q7),
                                                                                                                      minDate = last(Date[which(Q7 == min(Q7))])) %>% ungroup() %>% 
  mutate(staname = factor(staname, levels = rev(c('Naselle','NF Skokomish','NF Stillaguamish','Satsop',
                                                  'Puyallup at Puyallup','Skykomish','Snoqualmie nr Carnation','NF Snoqualmie',
                                                  'Dungeness','Sauk abv White Chuck','Sauk nr Sauk'
  )))) %>% arrange(staname,Year) %>% select(-staname)
summary.score <- left_join(summary.score,summary.min7q,by=c("site_no","Year")) 

cot <- calc_annual_flow_timing(data=df, dates = Date, values = Q, groups = site_no, 
                               percent_total = c(25, 33.3, 50, 75), water_year_start = 10, 
                               start_year=1920, end_year=2024, months = 1:12, transpose = FALSE)

cot <- left_join(cot,summary.score,by=c("site_no","Year")) %>% filter(!is.na(staname)) %>% rename(COT = DoY_50pct_TotalQ)

ggplot(cot,aes(Year,COT)) +
  theme_bw() +
  geom_point(color='gray') +
  geom_smooth() +
  facet_wrap(~staname) +
  facet_wrap2(~staname, strip = strip) +
  labs(x="",y="Julian Day (Oct 1 = 1)",title="Flow Center of Timing") +
  theme(strip.text = element_text(color = "white")) 

# ggplot(cot %>% mutate(Date_50pct_TotalQ = ifelse(!is.na(Date_50pct_TotalQ),as.Date(paste(2000,month(Date_50pct_TotalQ),day(Date_50pct_TotalQ),sep="-")),NA)),aes(`Date_50pct_TotalQ`,Year)) +
ggplot(cot %>% mutate(COT = as.Date("2000-10-01")+COT-1),aes(Year,COT)) +
  theme_bw() +
  geom_point(color='gray') +
  geom_smooth() +
  facet_wrap(~staname) +
  labs(x="",y="Calendar Day",title="Flow Center of Timing") 

#######################################################################################################################################
### data now ready for path analysis
#######################################################################################################################################
### 
mod1<-'
COT ~ a*wQ
pdays ~ b*COT + c*wQ
# indirect and total effects
a1 :=a
b1 :=b
c1 :=c
indirect:= a*b
total:= c + (a*b)
'

mod2<-'
COT ~ a*MAF
mQ ~ b*COT + c*MAF
# indirect and total effects
a1 :=a
b1 :=b
c1 :=c
indirect:= a*b
total:= c + (a*b)
'

################################################################################################
################################################################################################
### Perform path analysis using sem function in lavaan package
################################################################################################
yearstart = 1977
yearend = 2024

path.res <- NULL
path.res2 <- NULL

for(sta in staid){
# tmp <- filter(cot, site_no == sta&between(Year,1930,2024))
tmp <- filter(cot, site_no == sta&between(Year,yearstart,yearend))

### the path model execution
mod.res <- sem(mod2,data=tmp)

################################################################################################
### Results for rho, beta, and net effects in Figure 2
# fitMeasures(mod.res)
### get rho and beta - Kormos' correlation coefficient and standardized regression coefficient (direct effect), respectively...I think
### the results correspond almost exactly to results in Figure 2
R <- lavInspect(mod.res, what = 'cor.ov') # rho
R <- data.frame(Stat = "ρ", COT_mQ = round(R["COT","mQ"],2), MAF_COT = round(R["MAF","COT"],2), MAF_mQ = round(R["MAF","mQ"],2))
B = lavInspect(mod.res, what = 'std')$beta
B <- data.frame(Stat = "ß", COT_mQ = round(B["mQ","COT"],2), MAF_COT = round(B["COT","MAF"],2), MAF_mQ = round(B["mQ","MAF"],2))
# 'total' effect in summary(mod.res)$pe$std.all[11] 
NE = lavInspect(mod.res, what = 'std')$beta
NE <- data.frame(Stat = "NE", COT_mQ = round(NE["mQ","COT"],2), MAF_COT = round(NE["COT","MAF"],2), MAF_mQ = round(summary(mod.res, fit.measures=T,standardized=T,rsq=T,ci=T)$pe$std.all[11],2))

################################################################################################
# results table for Figure 2
stats <- bind_rows(R,B,NE) %>% mutate(site_no = tmp$site_no[1], staname = tmp$staname[1], lat = tmp$lat[1], lon = tmp$lon[1])
# stats %>% kable()
path.res <- bind_rows(path.res,stats)
################################################################################################
# long format results
R <- gather(R,-Stat, key = Path, value = Value)
B <- gather(B,-Stat, key = Path, value = Value)
NE <- gather(NE,-Stat, key = Path, value = Value)
stats2 <- bind_rows(R,B,NE) %>% mutate(site_no = tmp$site_no[1], staname = tmp$staname[1], lat = tmp$lat[1], lon = tmp$lon[1])
path.res2 <- bind_rows(path.res2,stats2)
################################################################################################

semPaths(mod.res, what='std', 
         edge.label.cex=1.25, curvePivot = TRUE, 
         fade=FALSE)
title(paste(tmp$staname[1],sta,sep=" "))

}

saveRDS(path.res,paste0('path.',yearstart,'.',yearend,'.rds'))
path.1930_1976 <- readRDS('path.1930.1976.rds')
path.1977_2024 <- readRDS('path.1977.2024.rds')

# # having trouble figuring out how to normalize symbols for each map layer
# mapview(path.res %>% filter(Stat == 'NE') %>% mutate(`COT_mQ` = scale(`COT_mQ`)), cex = 'COT_mQ', col.regions = "#56B4E9", xcol = 'lon', ycol = 'lat', crs = 4269, layer.name = 'NE COT') +
#   # mapview(path.res %>% filter(Stat == 'NE') %>% st_as_sf(coords = c("lon", "lat")) %>% st_jitter(factor = 0.01),
#   mapview(path.res %>% filter(Stat == 'NE') %>% mutate(`MAF_mQ` = scale(`MAF_mQ`), lat = lat-0.05, lon = lon-0.05), 
#           cex = 'MAF_mQ', col.regions = "#E69F00", xcol = 'lon', ycol = 'lat', crs = 4269, layer.name = 'NE MAF')   
# 
# #### try ggplot2
# wash_map <-  st_read("//kc.kingcounty.lcl/dnrp/WLRD/Users/UsersSTS1/degaspec/R/GLM/gis/counties.shp")
# st_crs(wash_map) = 2926
# wash_map <- st_transform(wash_map,crs=4326) %>% 
#   filter(!COUNTYLABE %in% c("Okanogan","Ferry","Pend Oreille","Spokane","Whitman","Garfield","Columbia","Asotin","Walla Walla","Franklin","Benton","Adams","Lincoln","Douglas","Grant","Chelan","Kittitas","Yakima","Klickitat","Stevens"))

path.resx <- st_as_sf(path.res2, coords = c("lon", "lat"), crs = 4326, remove = F)

tmp <- filter(path.resx, Stat == 'NE'&!Path == 'MAF_COT') %>% mutate(Path = factor(Path, levels = c("MAF_mQ","COT_mQ")))
ggplot(wash_map) + geom_sf(fill="transparent") + 
  theme_bw() +
  geom_sf(data = tmp %>% st_jitter(factor = 0.01), aes(size=Value, color = Path), alpha = 0.5) +
  ggtitle(paste0('Mean Summer Flow (Jul 15 - Sep 15) ',yearstart,' to ',yearend)) +
  lims(size = c(0,1)) +
  scale_color_manual(labels = c("MAF","COT"), values = c("#56B4E9","#E69F00")) +
  guides(size = guide_legend(title = "Net Effect"), color = guide_legend(override.aes = list(size=5)))

# p <- ggplot(path.res %>% filter(Stat == 'NE'), aes(MAF_mQ,COT_mQ,text=staname)) +
  ggplot(path.res %>% filter(Stat == 'NE'), aes(MAF_mQ,COT_mQ,text=staname,label=staname)) +
  theme_bw() +
  geom_point() +
  geom_text_repel(max.overlaps = Inf) +
  ggtitle(paste0('Mean Summer Flow (Jul 15 - Sep 15) ',yearstart,' to ',yearend)) +  
  expand_limits(x = c(0,1), y = c(0,1)) +
  geom_abline(slope = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
# ggplotly(p)



xlm <- range(c(path.1930_1976$MAF_mQ[path.1930_1976$Stat=='NE'],path.1977_2024$MAF_mQ[path.1977_2024$Stat=='NE']))
ylm <- range(c(path.1930_1976$COT_mQ[path.1930_1976$Stat=='NE'],path.1977_2024$COT_mQ[path.1977_2024$Stat=='NE']))

plot(path.1930_1976$COT_mQ[path.1930_1976$Stat=='NE']~path.1930_1976$MAF_mQ[path.1930_1976$Stat=='NE'],col=4,cex=1.1,xlim=xlm,ylim=ylm)
points(path.1977_2024$COT_mQ[path.1977_2024$Stat=='NE']~path.1977_2024$MAF_mQ[path.1977_2024$Stat=='NE'],pch=16,cex=0.9)
arrows(path.1930_1976$MAF_mQ[path.1930_1976$Stat=='NE'],path.1930_1976$COT_mQ[path.1930_1976$Stat=='NE'],path.1977_2024$MAF_mQ[path.1977_2024$Stat=='NE'],path.1977_2024$COT_mQ[path.1977_2024$Stat=='NE'],length=0.12)

path.1930_1976.x <- path.1930_1976 %>% filter(Stat=='NE') %>% select(-`MAF_COT`,-Stat) %>% rename(COT_mQ.1930_1976 = COT_mQ, MAF_mQ.1930_1976 = MAF_mQ)
path.1977_2024.x <- path.1977_2024 %>% filter(Stat=='NE') %>% select(-`MAF_COT`,-Stat) %>% rename(COT_mQ.1977_2024 = COT_mQ, MAF_mQ.1977_2024 = MAF_mQ)

path.mQ.w <- left_join(path.1930_1976.x,path.1977_2024.x,by=c('site_no','staname','lat','lon'))

path.1930_1976.y <- path.1930_1976 %>% filter(Stat=='NE') %>% select(-`MAF_COT`,-Stat) %>% mutate(Period = '1930-1976', staname = NA)
path.1977_2024.y <- path.1977_2024 %>% filter(Stat=='NE') %>% select(-`MAF_COT`,-Stat) %>% mutate(Period = '1977-2024')

path.mQ.l <-bind_rows(path.1930_1976.y,path.1977_2024.y)

ggplot( ) +
  theme_bw() +
  geom_point(data=path.mQ.l, aes(MAF_mQ,COT_mQ,color=Period),size=4) +
  scale_color_manual(values = c("#56B4E9","#E69F00")) +
  #geom_point(aes(MAF_mQ.1930_1976,COT_mQ.1930_1976), color = "#56B4E9", size = 4) +
  #geom_point(aes(MAF_mQ.1977_2024,COT_mQ.1977_2024), color = "#E69F00", size = 4) +
  geom_segment(data = path.mQ.w, aes(x=MAF_mQ.1930_1976, y=COT_mQ.1930_1976, xend=MAF_mQ.1977_2024, yend=COT_mQ.1977_2024), 
               arrow = arrow(), inherit.aes = F) +
  geom_text_repel(data = path.mQ.l,aes(MAF_mQ,COT_mQ,label=staname),max.overlaps = Inf) +
  ggtitle(paste0('Net Effects on Mean Summer Flow (Jul 15 - Sep 15) ')) +  
  expand_limits(x = c(0,1), y = c(0,1)) +
  geom_abline(slope = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
  



################################################################################################
### evaluate model and extract relevant statistics - rho (correlation coefficient), beta (direct effect), and total (net) effect
summary(mod.res, fit.measures=T,standardized=T,rsq=T,ci=T)

parameterestimates(mod.res) %>% 
  kable()

################################################################################################
### Results for rho, beta, and net effects in Figure 2
# fitMeasures(mod.res)
### get rho and beta - Kormos' correlation coefficient and standardized regression coefficient (direct effect), respectively...I think
### the results correspond almost exactly to results in Figure 2
R <- lavInspect(mod.res, what = 'cor.ov') # rho
R <- data.frame(Stat = "ρ", COT_mQ = round(R["COT","mQ"],2), MAF_COT = round(R["MAF","COT"],2), MAF_mQ = round(R["MAF","mQ"],2))
B = lavInspect(mod.res, what = 'std')$beta
B <- data.frame(Stat = "ß", COT_mQ = round(B["mQ","COT"],2), MAF_COT = round(B["COT","MAF"],2), MAF_mQ = round(B["mQ","MAF"],2))
# 'total' effect in summary(mod.res)$pe$std.all[11] 
NE = lavInspect(mod.res, what = 'std')$beta
NE <- data.frame(Stat = "NE", COT_mQ = round(NE["mQ","COT"],2), MAF_COT = round(NE["COT","MAF"],2), MAF_mQ = round(summary(mod.res, fit.measures=T,standardized=T,rsq=T,ci=T)$pe$std.all[11],2))

# results table for Figure 2
stats <- bind_rows(R,B,NE)
stats %>% kable()

################################################################################################
### simple path model plot....
### would like to try and see how close to a publication ready graph this could be made into
semPaths(mod.res, what='std', 
         edge.label.cex=1.25, curvePivot = TRUE, 
         fade=FALSE)

lavaanPlot(name = 'test', mod.res, coefs = T)

