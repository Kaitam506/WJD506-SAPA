#Software for the Automated Processing of Aerosol SAPA V.0.3.4
# New feautres: 
#Removes short, fat peaks which may be noise
#Processes data, with quantification
#Plot time series
#Can check where bottlenecks are using profvis command
#Cluster to speed up work present but not intergrated
#Source apportionment for CHOS & CHON components

#Install needed packages and loaded needed premade assets before progressing
install.packages("devtools")
devtools::install_github('davidcarslaw/openair')
devtools::install_github("kassambara/factoextra")
install.packages("ggplot2")
install.packages("ape")
install.packages("dendextend")
install.packages("doParallel")
install.packages("lubridate")
install.packages("openair")
install.packages("tidyr")
install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
install.packages.2('dendextend')
install.packages.2('colorspace')
install.packages("plotly")
install.packages("plotrix")


library("ape")
library(plotly)
library(ggplot2)
library(tidyr)
library(dplyr)
library(grid)
library(stringr)
library(plyr)
library(plotrix)
library(RgoogleMaps)
library(reshape2)
library(threadr)
library(openair)
library(lubridate)
library("devtools")
library("factoextra")
library(dendextend)
library(colorspace)
library(data.table)
library(plotrix)

#check for bottlenecks, put code in profvis command
#profvis({})
#Build multiplot function to put many plots together on one page, multiplot by Stuart Grange
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))}}}
#create clusters to speed up code
#cl<-makeCluster(4)
#clusterExport(cl, c(("")))
#clusterEvalQ(cl, {source()})
#stopCluster(cl)

#Input Tracefinder data, manipulate to work with code
CHINA.2016 <- read.csv("C:/Users/William/Documents/WACL stuff/06112017Chinaquant4.csv", na.strings="N/F")
CHINA.2016$Area <- as.numeric(as.character(CHINA.2016$Area))
CHINA.2016$Height <- as.numeric(as.character(CHINA.2016$Height))
#area to height ratio (removes short, wide peaks)
CHINA.2016<-CHINA.2016[which((CHINA.2016$Area/CHINA.2016$Height)<10 ),]
colnames(CHINA.2016)[1] <- "Compound"
colnames(CHINA.2016)[2] <- "RT"
allcompounds2 <- CHINA.2016
allcompounds2 <- allcompounds2 %>% 
  select(Sample.ID,
         Compound, 
         Area,
         RT)

allcompounds2 <- allcompounds2[!is.na(allcompounds2$Area),]
#input extra data, massloading for filter info, compound list for elemental composition, PM25 data from AMS, MetGas is from Leeds meteorlogical 
massloading <- read.csv("C:/Users/William/Documents/WACL stuff/massloading.csv", na.strings="N/F")
Compound_list <- read.csv("C:/Users/William/Documents/WACL stuff/compounds.csv", na.strings=0)
PM25 <- read.csv("C:/Users/William/Documents/WACL stuff/PM2.5mass.csv", na.strings=0)
summerMETGAS <- read.csv("C:/Users/William/Documents/WACL stuff/summerGasPhaseMET.csv", na.strings=0)
colnames(summerMETGAS)[1] <- "datetime"
summerMETGAS$datetime <- dmy_hms(summerMETGAS$datetime, tz = "UTC")

#Manipulate and marge extra files with orbitrap data
massloading <- massloading %>% 
  select(Sample.ID,
         datetime,
         Volume)

allcompounds2 <- allcompounds2 %>% 
  left_join(massloading, by = "Sample.ID")

allcompounds2 <- allcompounds2 %>% 
  left_join(Compound_list, by = "Compound")

#Calculate aerosol metrics & group by heteroatoms
allcompounds2$Area <- as.numeric(allcompounds2$Area)
allcompounds2$Volume <- as.numeric(allcompounds2$Volume)
allcompounds2$Sulphur[is.na(allcompounds2$Sulphur)] = 0
allcompounds2$Nitrogen[is.na(allcompounds2$Nitrogen)] = 0
allcompounds2$Sulphur <- as.numeric(allcompounds2$Sulphur)
allcompounds2$Nitrogen <- as.numeric(allcompounds2$Nitrogen)
allcompounds2$datetime <- dmy_hm(allcompounds2$datetime, tz = "UTC")
allcompounds2$Quant <- 0
allcompounds2$Group = "CHO"
allcompounds2$Group[allcompounds2$Nitrogen >= 1] = paste(allcompounds2$Group[allcompounds2$Nitrogen == 1],"N",sep = "")
allcompounds2$Group[allcompounds2$Sulphur == 1] = paste(allcompounds2$Group[allcompounds2$Sulphur == 1],"S",sep = "")
for(x in 1:length(allcompounds2$Compound)){
  allcompounds2$DBE[x] <- (allcompounds2$Carbon[x]+1-((allcompounds2$Hydrogen[x]/2)+(allcompounds2$Nitrogen[x]/2)))}
for(x in 1:length(allcompounds2$Compound)){
  allcompounds2$AI[x] <- ((1+allcompounds2$Carbon[x]-allcompounds2$Oxygen[x]-allcompounds2$Sulphur[x]-(allcompounds2$Hydrogen[x]/2))/(allcompounds2$Carbon[x]-allcompounds2$Oxygen[x]-allcompounds2$Sulphur[x]-allcompounds2$Nitrogen[x]))}
for(x in 1:length(allcompounds2$Compound))
{allcompounds2$HCratio[x] <- (allcompounds2$Hydrogen[x]/allcompounds2$Carbon[x])}
for(x in 1:length(allcompounds2$Compound))
{allcompounds2$OCratio[x] <- (allcompounds2$Oxygen[x]/allcompounds2$Carbon[x])}
{for(x in 1:length(allcompounds2$Area))
  allcompounds2$Quant[x] <- ((allcompounds2$Area[x]*8)/allcompounds2$Volume[x])}   #area times by 8, as 1/8 of filter taken
allcompounds2[allcompounds2<0] <- 0
allcompounds2 <- allcompounds2[!is.na(allcompounds2$Area),]
for(x in 1:length(allcompounds2$Compound))
{allcompounds2$OSc[x] <- (((allcompounds2$Oxygen[x]/allcompounds2$Carbon[x])*2)-(allcompounds2$Hydrogen[x]/allcompounds2$Carbon[x]))}

#pick summer or winter samples
allcompounds<-allcompounds2
allcompounds$datetime <- ymd_hms(allcompounds$datetime, tz = "UTC")
allcompoundsSummer <- allcompounds[which(allcompounds$datetime>"2016-12-10 23:59:00" ),]
allcompoundsWinter <- allcompounds[which(allcompounds$datetime<"2016-12-10 23:59:00" ),]
allcompounds2 <- allcompoundsWinter

Samples <- allcompounds2[!grepl("ppb", allcompounds2$Sample.ID),]
Samples <- Samples[!grepl("Blank", Samples$Sample.ID),]
Standards <- CHINA.2016[grepl("ppb", CHINA.2016$Sample.ID),]
Standards$Conc <- 0
#Input values for standard concentrations
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHON100ppb') Standards$Conc[x] <- 100}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHON50ppb') Standards$Conc[x] <- 50}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHON25ppb') Standards$Conc[x] <- 25}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHON12ppb') Standards$Conc[x] <- 12.5}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHON6ppb') Standards$Conc[x] <- 6.25}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHON3ppb') Standards$Conc[x] <- 3.125}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHON1ppb') Standards$Conc[x] <- 1.5625}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHOS100ppb') Standards$Conc[x] <- 100}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHOS50ppb') Standards$Conc[x] <- 50}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHOS25ppb') Standards$Conc[x] <- 25}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHOS12ppb') Standards$Conc[x] <- 12.5}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHOS6ppb') Standards$Conc[x] <- 6.25}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHOS3ppb') Standards$Conc[x] <- 3.125}
{for(x in 1:length(Standards$Sample.ID))
  if (Standards$Sample.ID[x] == 'CHOS1ppb') Standards$Conc[x] <- 1.5625}
pinonicSTD <- Standards[grepl("cis", Standards$Compound),]
pinonicSTD <- pinonicSTD[grepl("CHO", pinonicSTD$Sample.ID),]
STDcali24np <- Standards[ which(Standards$Compound=='2,4-Dinitrophenol'),]
STDcali24np <- STDcali24np[ which(STDcali24np$Area > 10000),]
ggplot(STDcali24np, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("2,4-Dinitrophenol concentration (ppb)")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("2,4-Dinitrophenol concentration (ppb)") + geom_smooth(method=lm)
STDcalinaph <- Standards[ which(Standards$Compound=='2-Nitro-1-Naphthol'),]
STDcalinaph <- STDcalinaph[ which(STDcalinaph$Area > 10000),]
ggplot(STDcalinaph, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("2-Nitro-1-Naphthol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcali3n2p <- Standards[ which(Standards$Compound=='3-Methyl-2-nitrophenol'),]
STDcali3n2p <- STDcali3n2p[ which(STDcali3n2p$Area > 10000),]
ggplot(STDcali3n2p, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("3-Methyl-2-nitrophenol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcali26dim4np <- Standards[ which(Standards$Compound=='2,6-dimethyl-4-nitrophenol'),]
STDcali26dim4np <- STDcali26dim4np[ which(STDcali26dim4np$Area > 10000),]
ggplot(STDcali26dim4np, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("2,6-dimethyl-4-nitrophenol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcali2n3p <- Standards[ which(Standards$Compound=='2-Methyl-3-nitrophenol'),]
STDcali2n3p <- STDcali2n3p[ which(STDcali2n3p$Area > 100000),]
ggplot(STDcali2n3p, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("2-Methyl-3-nitrophenol concentration (ppb)")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("2-Methyl-3-nitrophenol concentration (ppb)") + geom_smooth(method=lm)
STDcali4m2np <- Standards[ which(Standards$Compound=='4-Methyl-2-nitrophenol'),]
STDcali4m2np <- STDcali4m2np[ which(STDcali4m2np$Area > 1),]
ggplot(STDcali4m2np, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("4-Methyl-2-nitrophenol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcali2np <- Standards[ which(Standards$Compound=='2-Nitrophenol'),]
STDcali2np <- STDcali2np[ which(STDcali2np$Area > 1),]
ggplot(STDcali2np, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("2-Nitrophenol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcali3np <- Standards[ which(Standards$Compound=='3-Nitrophenol'),]
STDcali3np <- STDcali3np[ which(STDcali3np$Area > 1),]
ggplot(STDcali3np, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("3-Nitrophenol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcali4np <- Standards[ which(Standards$Compound=='4-Nitrophenol'),]
STDcali4np <- STDcali4np[ which(STDcali4np$Area > 1),]
ggplot(STDcali4np, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("4-Nitrophenol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcali2n5p <- Standards[ which(Standards$Compound=='2-Methyl-5-nitrophenol'),]
STDcali2n5p <- STDcali2n5p[ which(STDcali2n5p$Area > 100000),]
ggplot(STDcali2n5p, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("2-Methyl-5-nitrophenol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcalidodecyl <- Standards[ which(Standards$Compound=='Dodecyl Sulfate'),]
STDcalidodecyl <- STDcalidodecyl[ which(STDcalidodecyl$Area > 6800000),]
ggplot(STDcalidodecyl, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("Dodecyl Sulfate")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcalioctyl <- Standards[ which(Standards$Compound=='Octyl Sulfate'),]
STDcalioctyl <- STDcalioctyl[ which(STDcalioctyl$Area > 1000000),]
ggplot(STDcalioctyl, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("Octyl Sulfate")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcalinaphyl <- Standards[ which(Standards$Compound=='2-Napthyl Sulfate'),]
STDcalinaphyl <- STDcalinaphyl[ which(STDcalinaphyl$Area > 100000),]
ggplot(STDcalinaphyl, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("2-Naphthyl Sulfate")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcalipinon <- Standards[ which(Standards$Compound=='cis-pinonic acid (pinene SOA)'),]
STDcalipinon <- STDcalipinon[ which(STDcalipinon$Area > 10000),]
ggplot(STDcalipinon, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("cis-pinonic acid (pinene SOA)")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcalicres <- Standards[ which(Standards$Compound=='C7H7NO3 - RT6.86'),]
STDcalicres <- STDcalicres[ which(STDcalicres$Area > 100000),]
ggplot(STDcalicres, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("4-nitr-m-cresol")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)
STDcalimann <- Standards[ which(Standards$Compound=='C6H10O5SO4 - RT0.71'),]
STDcalimann <- STDcalimann[ which(STDcalimann$Area > 100000),]
ggplot(STDcalimann, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) +ggtitle("6-mannose-d-sulfate")  + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("Concentration ppb") + geom_smooth(method=lm)

STDcalinitrosulf <- CHINA.2016[ which(CHINA.2016$Compound=='C6H5NO6S - RT2.08'),]
STDcalinitrosulf <- STDcalinitrosulf[ which(STDcalinitrosulf$Area < 30000),]



rsq <- function (x, y) cor(x, y) ^ 2
rsq(STDcalipinon$Conc, STDcalipinon$Area)
lm.pinon<-lm(formula= STDcalipinon$Area~STDcalipinon$Conc)
summary(lm.pinon)

rsq(STDcali2n3p$Conc, STDcali2n3p$Area)
lm.2n3p<-lm(formula= STDcali2n3p$Area~STDcali2n3p$Conc)
summary(lm.2n3p)

rsq(STDcali2n5p$Conc, STDcali2n5p$Area)
lm.2n5p<-lm(formula= STDcali2n5p$Area~STDcali2n5p$Conc)
summary(lm.2n5p)

rsq(STDcali4m2np$Conc, STDcali4m2np$Area)
lm.4m2np<-lm(formula= STDcali4m2np$Area~STDcali4m2np$Conc)
summary(lm.4m2np)

rsq(STDcali24np$Conc, STDcali24np$Area)
lm.24np<-lm(formula= STDcali24np$Area~STDcali24np$Conc)
summary(lm.24np)

rsq(STDcali26dim4np$Conc, STDcali26dim4np$Area)
lm.26dim4np<-lm(formula= STDcali26dim4np$Area~STDcali26dim4np$Conc)
summary(lm.26dim4np)

rsq(STDcali3n2p$Conc, STDcali3n2p$Area)
lm.3n2p<-lm(formula= STDcali3n2p$Area~STDcali3n2p$Conc)
summary(lm.3n2p)

rsq(STDcali4np$Conc, STDcali4np$Area)
lm.4np<-lm(formula= STDcali4np$Area~STDcali4np$Conc)
summary(lm.4np)

rsq(STDcali3np$Conc, STDcali3np$Area)
lm.3np<-lm(formula= STDcali3np$Area~STDcali3np$Conc)
summary(lm.3np)

rsq(STDcali2np$Conc, STDcali2np$Area)
lm.2np<-lm(formula= STDcali2np$Area~STDcali2np$Conc)
summary(lm.2np)

rsq(STDcalinaph$Conc, STDcalinaph$Area)
lm.naph<-lm(formula= STDcalinaph$Area~STDcalinaph$Conc)
summary(lm.naph)

rsq(STDcalioctyl$Conc, STDcalioctyl$Area)
lm.octyl<-lm(formula= STDcalioctyl$Area~STDcalioctyl$Conc)
summary(lm.octyl)

rsq(STDcalinaphyl$Conc, STDcalinaphyl$Area)
lm.naphyl<-lm(formula= STDcalinaphyl$Area~STDcalinaphyl$Conc)
summary(lm.naphyl)

rsq(STDcalidodecyl$Conc, STDcalidodecyl$Area)
lm.dodecyl<-lm(formula= STDcalidodecyl$Area~STDcalidodecyl$Conc)
summary(lm.dodecyl)

rsq(STDcalicres$Conc, STDcalicres$Area)
lm.cres<-lm(formula= STDcalicres$Area~STDcalicres$Conc)
summary(lm.cres)

rsq(STDcalimann$Conc, STDcalimann$Area)
lm.mann<-lm(formula= STDcalimann$Area~STDcalimann$Conc)
summary(lm.mann)


SummedOC <- aggregate(Quant ~ Sample.ID, allcompounds2, sum)
Trends <- Samples[ which(Samples$Compound=='C7H7NO3 - RT6.86'),]
Trends2 <- Samples[ which(Samples$Compound=='2,6-dimethyl-4-nitrophenol'),]
Trends <- rbind(Trends, Trends2)
Trends3 <- Samples[ which(Samples$Compound=='C5H12SO4 - RT4.32'),]
Trends4 <- Samples[ which(Samples$Compound=='Octyl Sulfate'),]
Trends3 <- rbind(Trends3, Trends4)
Trends5 <- Samples[ which(Samples$Compound=='C5H9COOH - RT1.79'),]
Trends6 <- Samples[ which(Samples$Compound=='C10H16O2 - RT6.83 (a-pinene SOA)'),]
Trends5 <- rbind(Trends5, Trends6)
Quant <- Samples[ which(Samples$Compound=='Octyl Sulfate'),]
lm.oct <- lm(octSTD$Area ~ octSTD$Conc)
oct.slope <- unname(lm.oct$coefficients[2])
oct.intercept <- unname(lm.oct$coefficients[1])
Quant2 <- Samples[ which(Samples$Compound=='cis-pinonic acid (pinene SOA)'),]
Quant <- rbind(Quant, Quant2)
Quant2 <- Samples[ which(Samples$Compound=='2,6-dimethyl-4-nitrophenol'),]
Quant <- rbind(Quant, Quant2)
Quant$Conc <- 0
lm.cis <- lm(pinonicSTD$Area ~ pinonicSTD$Conc)
cis.slope <- unname(lm.cis$coefficients[2])
cis.intercept <- unname(lm.cis$coefficients[1])
lm.dimeth <- lm(dimethSTD$Area ~ dimethSTD$Conc)
dimeth.slope <- unname(lm.dimeth$coefficients[2])
dimeth.intercept <- unname(lm.dimeth$coefficients[1])
for(x in 1:length(Quant$Compound))
  if (Quant$Compound[x] == 'Octyl Sulfate') {Quant$Conc[x] <- (((Quant$Area[x]-780972)/234125)/Quant$Volume[x])}
for(x in 1:length(Quant$Compound))
  if (Quant$Compound[x] == 'cis-pinonic acid (pinene SOA)') {Quant$Conc[x] <- (((Quant$Area[x]-25915)/5151.9)/Quant$Volume[x])}
for(x in 1:length(Quant$Compound))
  if (Quant$Compound[x] == '2,6-dimethyl-4-nitrophenol') {Quant$Conc[x] <- (((Quant$Area[x]-1448624)/273600)/Quant$Volume[x])}
Avg <- aggregate(Quant ~ Compound, Samples, mean)
Avg <- Avg %>% 
  right_join(Samples, by = "Compound") 

Avg <- Avg %>% 
  select(Compound,
         Quant.x,
         Carbon,
         Hydrogen,
         Nitrogen,
         Sulphur,
         Oxygen,
         DBE,
         OCratio,
         HCratio,
         Group,
         OSc,
         RT)
colnames(Avg)[2] <- "Quant"
ggplot(Avg, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                panel.grid.major = element_blank(),
                                                                                                panel.grid.minor = element_blank(),
                                                                                                panel.background = element_blank()) +ggtitle("Summed compounds") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + xlab("Carbon number")
ggplot(Trends, aes(x=datetime, y=Quant, group=Compound))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", values = c("red", "sienna3", "skyblue2", "dodgerblue3", "darkred", "steelblue4")) + theme(legend.position="right") + scale_y_log10() + labs(x = "Date", y = "Log of peak area per volume sampled (a.u. m-3)")
ggplot(Quant, aes(x=datetime, y=Conc, group=Compound))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ppb)")
Day18 <- Samples[grepl("11-18", Samples$datetime),]
Day22 <- Samples[grepl("11-22", Samples$datetime),]
Day27 <- Samples[grepl("11-27", Samples$datetime),]
Day04 <- Samples[grepl("12-04", Samples$datetime),]
HetGroup <- Samples[!grepl("Sulphate", Samples$Compound),]
HetGroupCHO <- HetGroup[ which(HetGroup$Group=='CHO'),]
HetGroupCHO <- aggregate(Quant ~ datetime, HetGroupCHO, mean)
HetGroupCHO$Group <- "CHO"
HetGroupCHON <- HetGroup[ which(HetGroup$Group=='CHON'),]
HetGroupCHON <- aggregate(Quant ~ datetime, HetGroupCHON, mean)
HetGroupCHON$Group <- "CHON"
HetGroupCHOS <- HetGroup[ which(HetGroup$Group=='CHOS'),]
HetGroupCHOS <- aggregate(Quant ~ datetime, HetGroupCHOS, mean)
HetGroupCHOS$Group <- "CHOS"
HetGroupCHONS <- HetGroup[ which(HetGroup$Group=='CHONS'),]
HetGroupCHONS <- aggregate(Quant ~ datetime, HetGroupCHONS, mean)
HetGroupCHONS$Group <- "CHONS"
HetGroup<-rbind(HetGroupCHO, HetGroupCHON, HetGroupCHOS)
plot18<-ggplot(Day18, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                          panel.grid.major = element_blank(),
                                                                                                          panel.grid.minor = element_blank(),
                                                                                                          panel.background = element_blank()) +ggtitle("18/11/2016") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + xlab("Carbon number")
plot22<-ggplot(Day22, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                          panel.grid.major = element_blank(),
                                                                                                          panel.grid.minor = element_blank(),
                                                                                                          panel.background = element_blank()) +ggtitle("22/11/2016") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + xlab("Carbon number")
plot27<-ggplot(Day27, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                          panel.grid.major = element_blank(),
                                                                                                          panel.grid.minor = element_blank(),
                                                                                                          panel.background = element_blank()) +ggtitle("27/11/2016") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + xlab("Carbon number")
plot04<-ggplot(Day04, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                          panel.grid.major = element_blank(),
                                                                                                          panel.grid.minor = element_blank(),
                                                                                                          panel.background = element_blank()) +ggtitle("04/12/2016") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + xlab("Carbon number")
multiplot(plot18, plot27, plot22, plot04, cols=2)
ggplot(HetGroup, aes(x=datetime, y=Quant, group=Group))  + geom_line(aes(colour = Group), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="right") +  labs(x = "Date", y = "Peak area per volume sampled (a.u. m-3)")
ggplot(HetGroup, aes(x=datetime, y=Quant, group=Group))  + geom_line(aes(colour = Group), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position="right") +  labs(x = "Date", y = "Peak area per volume sampled (a.u. m-3)") + scale_y_log10()

plotCHO<-ggplot(HetGroupCHO, aes(x=datetime, y=Quant))  + geom_line(color="red", size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual() + theme(legend.position="right") +  labs(x = "Date", y = "Peak area (a.u. m-3)") +ggtitle("Total CHO")
plotCHOS<-ggplot(HetGroupCHOS, aes(x=datetime, y=Quant))  + geom_line(color="blue", size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual() + theme(legend.position="right") +  labs(x = "Date", y = "Peak area (a.u. m-3)") +ggtitle("Total CHOS")
plotCHON<-ggplot(HetGroupCHON, aes(x=datetime, y=Quant))  + geom_line(color="green", size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual() + theme(legend.position="right") +  labs(x = "Date", y = "Peak area (a.u. m-3)") +ggtitle("Total CHON")
plottrend3<-ggplot(Trends, aes(x=datetime, y=Quant, group=Compound))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", values = c("red", "darkred", "steelblue4")) +   labs(x = "Date", y = "Peak area (a.u. m-3)")
plottrend2<-ggplot(Trends3, aes(x=datetime, y=Quant, group=Compound))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", values = c("skyblue2", "dodgerblue3", "darkred", "steelblue4")) +   labs(x = "Date", y = "Peak area (a.u. m-3)")
plottrend1<-ggplot(Trends5, aes(x=datetime, y=Quant, group=Compound))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", values = c("green", "darkgreen", "darkred", "steelblue4")) +  labs(x = "Date", y = "Peak area (a.u. m-3)")

multiplot(plotCHO, plotCHOS, plotCHON, plottrend1, plottrend2, plottrend3, cols=2)
WinAMS <- read.csv("C:/Users/William/Documents/WACL stuff/AMSedit.csv", na.strings=0)
SumAMS <- read.csv("C:/Users/William/Documents/WACL stuff/SummerAMSEdited.csv", na.strings=0)
AMS <- WinAMS
AMS$Sample.ID <- as.character(AMS$Sample.ID)
AMS$datetime <- dmy_hm(AMS$datetime, tz = "UTC")

sulfate <- HetGroupCHOS %>% 
  left_join(AMS, by = "datetime")

nitrate <- HetGroupCHON %>% 
  left_join(AMS, by = "datetime")

oxy <- HetGroupCHO %>% 
  left_join(AMS, by = "datetime")

oxy2 <- SummedOC

oxy2 <- oxy2 %>% 
  left_join(massloading, by = "Sample.ID") %>%
  select(-Volume)

oxy2$datetime <- dmy_hm(oxy2$datetime, tz = "UTC")

oxy2 <- oxy2 %>% 
  left_join(AMS, by = "datetime")

sulfate$SO4 <- as.numeric(sulfate$SO4)
plotsulfate <- ggplot(sulfate, aes(x = datetime))
plotsulfate <- plotsulfate + geom_line(aes(y = Quant, colour = "Total CHOS")) 
plotsulfate <- plotsulfate + geom_line(aes(y = SO4*5000, colour = "AMS SO4"))
plotsulfate <- plotsulfate + scale_y_continuous(sec.axis = sec_axis(~.*0.0002, name = "SO4 ppb"))
plotsulfate <- plotsulfate + scale_colour_manual(values = c("blue", "red"))
plotsulfate <- plotsulfate + labs(y = "Peak area (a.u. m-3)",
                                  x = "Date",
                                  colour = "Parameter")
plotsulfate <- plotsulfate + theme(axis.text=element_text(size=8), axis.title=element_text(size=8), legend.position = c(0.05, 0.8), legend.title=element_text(size=8), legend.text=element_text(size=7))
nitrate$NO3 <- as.numeric(nitrate$NO3)
plotnitrate <- ggplot(nitrate, aes(x = datetime))
plotnitrate <- plotnitrate + geom_line(aes(y = Quant, colour = "Total CHON"))
plotnitrate <- plotnitrate + geom_line(aes(y = NO3*25000, colour = "AMS NO3"))
plotnitrate <- plotnitrate + scale_y_continuous(sec.axis = sec_axis(~.*0.00004, name = "NO3 ppb"))
plotnitrate <- plotnitrate + scale_colour_manual(values = c("blue", "red"))
plotnitrate <- plotnitrate + labs(y = "Peak area (a.u. m-3)",
                                  x = "Date",
                                  colour = "Parameter")
plotnitrate <- plotnitrate + theme(axis.text=element_text(size=8), axis.title=element_text(size=8), legend.position = c(0.05, 0.8), legend.title=element_text(size=8), legend.text=element_text(size=7))
oxy$Org <- as.numeric(oxy$Org)
plotoxy <- ggplot(oxy, aes(x = datetime))
plotoxy <- plotoxy + geom_line(aes(y = Quant, colour = "Total CHO"))
plotoxy <- plotoxy + geom_line(aes(y = Org*2500, colour = "Organic AMS"))
plotoxy <- plotoxy + scale_y_continuous(sec.axis = sec_axis(~.*0.0004, name = "Organic ppb"))
plotoxy <- plotoxy + scale_colour_manual(values = c("blue", "red"))
plotoxy <- plotoxy + labs(y = "Peak area (a.u. m-3)",
                          x = "Date",
                          colour = "Parameter")
plotoxy <- plotoxy + theme(axis.text=element_text(size=8), axis.title=element_text(size=8), legend.position = c(0.05, 0.8), legend.title=element_text(size=8), legend.text=element_text(size=7))
oxy2$Org <- as.numeric(oxy2$Org)
plotoxy2 <- ggplot(oxy2, aes(x = datetime))
plotoxy2 <- plotoxy2 + geom_line(aes(y = Quant, colour = "Total OC"))
plotoxy2 <- plotoxy2 + geom_line(aes(y = Org*2500000, colour = "Organic AMS"))
plotoxy2 <- plotoxy2 + scale_y_continuous(sec.axis = sec_axis(~.*0.0000004, name = "Organic ppb"))
plotoxy2 <- plotoxy2 + scale_colour_manual(values = c("blue", "red"))
plotoxy2 <- plotoxy2 + labs(y = "Peak area (a.u. m-3)",
                            x = "Date",
                            colour = "Parameter")
plotoxy2 <- plotoxy2 + theme(axis.text=element_text(size=8), axis.title=element_text(size=8), legend.position = c(0.05, 0.8), legend.title=element_text(size=8), legend.text=element_text(size=7))
multiplot(plotsulfate, plotnitrate, plotoxy, plotoxy2, cols=1)
# Do some interactive plotting
nitrate %>% 
  rename(date=datetime) %>% 
  mutate(NO3 = NO3 * 50000) %>% 
  time_dygraph(variable = c("Quant", "NO3"))
# Create timeseries object
nitrate_ts <- nitrate 
names(nitrate_ts)[1] <- "date"
nitrate_ts <- nitrate_ts%>% 
  select(date,
         NO3,
         Quant) %>% 
  data_frame_to_timeseries()
# Clean up
nitrate_ts <- cbind(nitrate_ts$NO3, nitrate_ts$Quant)
names(nitrate_ts) <- c("NO3", "Quant")
# Plot two-axes
dygraph(nitrate_ts) %>%
  dySeries() %>% 
  dyAxis("y", label = "NO3") %>%
  dyAxis("y2", label = "Peak area per volume (a.u. m-3)", axisLabelFontSize = 10, axisLabelWidth = 70, independentTicks = TRUE) %>%
  dySeries("Quant", axis = 'y2')%>% 
  dyRangeSelector()

R2 <- AMS %>% 
  left_join(HetGroupCHO, by = "datetime") %>%
  select(-Group,
         -Sample.ID) 
names(R2)[13] <- "CHO"
R2 <- R2 %>% 
  left_join(HetGroupCHOS, by = "datetime") %>%
  select(-Group) 
names(R2)[14] <- "CHOS"
R2 <- R2 %>% 
  left_join(HetGroupCHON, by = "datetime") %>%
  select(-Group) 
names(R2)[15] <- "CHON"
corPlot(R2, pollutants = c("CHO","CHOS","CHON","Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS")
names(PM25)[1] <- "datetime"
PM25$datetime <- dmy_hm(PM25$datetime, tz = "UTC")
R2low <- PM25[which(PM25$PM2.5 < 50 ),]
R2high <- PM25[which(PM25$PM2.5 > 100 ),]
R2low <- R2low %>% 
  left_join(R2, by = "datetime") 
R2low <- subset(R2low, select = -c(X,Sample.ID) )
R2high <- R2high %>% 
  left_join(R2, by = "datetime") 
R2high <- subset(R2high, select = -c(X,Sample.ID) )
corPlot(R2low, pollutants = c("CHO","CHOS","CHON","Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS, clean days", order="alphabet", cluster=FALSE)
corPlot(R2high, pollutants = c("CHO","CHOS","CHON","Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS, polluted days", order="alphabet", cluster=FALSE)
AMSnoedit <- read.csv("C:/Users/William/Documents/WACL stuff/AMSnoedit.csv", na.strings=0)
AMSnoeditlow <- AMSnoedit[which(AMSnoedit$Org < 50 ),]
AMSnoedithigh <- AMSnoedit[which(AMSnoedit$Org > 50 ),]
corPlot(AMSnoeditlow, pollutants = c("Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with full AMS, clean days", order="alphabet", cluster=FALSE)
corPlot(AMSnoedithigh, pollutants = c("Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with full AMS, polluted days", order="alphabet", cluster=FALSE)
orgsulf <- Samples[which(Samples$Group == "CHOS" ),]
orgsulf <- orgsulf %>%
  select(Compound,
         Sample.ID,
         Quant,
         datetime) 
library(reshape)
md<-melt(orgsulf, id=c("Sample.ID", "Compound"))
orgsulf2<-cast(md, id+Compound~value)
class(orgsulf$Compound)
orgsulf <- orgsulf %>%
  select(Compound,
         Quant) 
orgsulf2 <- unstack(orgsulf)
md2<-melt(orgsulf, id=c("Compound"))
orgsulf3<-t(md2)
md2$num <- 1:(unique(md2[1]))
md2 <- md2[c(1,2,4,3)]
library("reshape2", lib.loc="~/R/win-library/3.2")
md2 <- md2 %>% 
  select(-variable) 
md3 <- dcast(md2, Compound ~ num, value.var="value", fill=0)
md3<-reshape(md2, idvar=c("Compound","num"), direction="wide")
md3<-data.frame(unique(md2[1]), 
                as.data.frame.matrix(xtabs(value ~ do.call(paste, md2[1]) + num, md2)))
#correlated just CHOS
orgsulf <- Samples[which(Samples$Group == "CHOS" ),]
orgsulf <- orgsulf %>%
  select(Compound,
         Sample.ID,
         Quant,
         datetime) 

n = unique(orgsulf[,1])

timeseries = seq(1478735940,1481241540,60)
timeseries = data.frame(unixtime = timeseries)
orgsulf$unixtime = as.numeric(orgsulf$datetime)

temp_compound = orgsulf[orgsulf$Compound == n[1],]
temp_compound = temp_compound[,c("Quant","unixtime")]
temp_compound = data.frame(temp_compound)
names(temp_compound) = c(n[1],"unixtime")
timeseries = left_join(timeseries,temp_compound,"unixtime")

for (i in 2:length(n)){
  temp_compound = orgsulf[orgsulf$Compound == n[i],]
  temp_compound = temp_compound[,c("Quant","unixtime")]
  temp_compound = data.frame(temp_compound)
  names(temp_compound) = c(n[i],"unixtime")
  
  timeseries = left_join(timeseries,temp_compound,"unixtime")
}
#possibly kill computer with 191 corplot of CHOS
filepath = paste(getwd(),"/CHOS_corplot.pdf",sep = "")
pdf(width=60, height=60, file=filepath)
corPlot(timeseries,
        pollutants = names(timeseries)[c(2:191)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors with AMS, CHOS",
        order="alphabet",
        cluster=TRUE)
dev.off()
DBE22<-ggplot(Day22, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                         panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(),
                                                                                                         panel.background = element_blank()) +ggtitle("22/11/2016, cleaner day") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + xlab("Carbon number")

DBE27<-ggplot(Day27, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                         panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(),
                                                                                                         panel.background = element_blank()) +ggtitle("27/11/2016, polluted day") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + xlab("Carbon number")

VK22<-ggplot(Day22, aes(OCratio, HCratio, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                             panel.grid.major = element_blank(),
                                                                                                             panel.grid.minor = element_blank(),
                                                                                                             panel.background = element_blank()) +ggtitle("22/11/2016, cleaner day") + xlim(0, 2.5) + ylim(0, 2.5) + guides(size=FALSE) + xlab("O : C ratio") + ylab("H : C ratio")

VK27<-ggplot(Day27, aes(OCratio, HCratio, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                                                             panel.grid.major = element_blank(),
                                                                                                             panel.grid.minor = element_blank(),
                                                                                                             panel.background = element_blank()) +ggtitle("27/11/2016, polluted day") + xlim(0, 2.5) + ylim(0, 2.5) + guides(size=FALSE) + xlab("O : C ratio") + ylab("H : C ratio")
multiplot(VK22, DBE22, VK27, DBE27, cols=2)
SummedOC <- SummedOC %>% 
  right_join(orgsulf, by = "Sample.ID")
SummedOC <- SummedOC %>%
  select(datetime,
         Quant.x)
names(SummedOC)[2] <- "Quant"
SummedOC$Group <- "Total OC"
#remove blocked filters
SummedOC1 <- SummedOC[which(SummedOC$datetime<"2016-11-17 23:59:00" ),]
SummedOC2 <- SummedOC[which(SummedOC$datetime>="2016-11-17 23:59:00" ),]
SummedOC1$num <- 1
SummedOC2$num <- 2
SummedOC<-rbind(SummedOC1,SummedOC2)
HetGroup <- rbind(HetGroup, HetGroupCHONS)
HetGroup <- HetGroup[which(HetGroup$datetime!="2016-12-03 23:59:00" ),]
HetGroup <- HetGroup[which(HetGroup$datetime!="2016-11-25 23:59:00" ),]
HetGroupCHO1 <- HetGroupCHO[which(HetGroupCHO$datetime<"2016-11-17 23:59:00" ),]
HetGroupCHO1$num <- 1
HetGroupCHO2 <- HetGroupCHO[which(HetGroupCHO$datetime>="2016-11-17 23:59:00" ),]
HetGroupCHO2$num <- 2
HetGroupCHO<-rbind(HetGroupCHO1,HetGroupCHO2)
HetGroupCHON1 <- HetGroupCHON[which(HetGroupCHON$datetime<"2016-11-17 23:59:00" ),]
HetGroupCHON1$num <- 1
HetGroupCHON2 <- HetGroupCHON[which(HetGroupCHON$datetime>="2016-11-17 23:59:00" ),]
HetGroupCHON2$num <- 2
HetGroupCHON<-rbind(HetGroupCHON1,HetGroupCHON2)
HetGroupCHOS1 <- HetGroupCHOS[which(HetGroupCHOS$datetime<"2016-11-17 23:59:00" ),]
HetGroupCHOS1$num <- 1
HetGroupCHOS2 <- HetGroupCHOS[which(HetGroupCHOS$datetime>="2016-11-17 23:59:00" ),]
HetGroupCHOS2$num <- 2
HetGroupCHOS<-rbind(HetGroupCHOS1,HetGroupCHOS2)
HetGroupCHONS1 <- HetGroupCHONS[which(HetGroupCHONS$datetime<"2016-11-17 23:59:00" ),]
HetGroupCHONS1$num <- 1
HetGroupCHONS2 <- HetGroupCHONS[which(HetGroupCHONS$datetime>="2016-11-17 23:59:00" ),]
HetGroupCHONS2$num <- 2
HetGroupCHONS<-rbind(HetGroupCHONS1,HetGroupCHONS2)
SummedOC1 <- SummedOC[which(SummedOC$datetime<"2016-11-17 23:59:00" ),]
SummedOC2 <- SummedOC[which(SummedOC$datetime>="2016-11-17 23:59:00" ),]
HetGroup2 <- rbind(HetGroup, SummedOC)
trend1<-ggplot(HetGroupCHO, aes(x=datetime, y=Quant, group=num)) + geom_line(color='green4', size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHO Peak area (a.u. m-3)") +theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                                                                                           axis.title=element_text(size=8))
trend2<-ggplot(HetGroupCHON, aes(x=datetime, y=Quant, group=num)) + geom_line(color='navyblue', size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHON Peak area (a.u. m-3)") +theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                                                                                               axis.title=element_text(size=8))
trend3<-ggplot(HetGroupCHOS, aes(x=datetime, y=Quant, group=num)) + geom_line(color='red2', size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHOS Peak area (a.u. m-3)")+theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                                                                                          axis.title=element_text(size=8))
trend4<-ggplot(HetGroupCHONS, aes(x=datetime, y=Quant, group=num)) + geom_line(color='magenta4', size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHONS Peak area (a.u. m-3)")+theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                                                                                                axis.title=element_text(size=8))
trend5<-ggplot(SummedOC, aes(x=datetime, y=Quant, group=num)) + geom_line(size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "Total OC Peak area (a.u. m-3)") +theme(axis.text=element_text(size=8),
                                                                                                                                                                                                                                                                                                                             axis.title=element_text(size=8))
multiplot(trend1, trend2,trend3, trend4, trend5, cols=1)
R22 <- R2 %>% 
  right_join(SummedOC, by = "datetime")
names(R22)[16] <- "TotalOC"
corPlot(R22, pollutants = c("Org","SO4","TotalOC","CHOS"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS", order="alphabet", cluster=FALSE)
trend6<-ggplot(R2, aes(x=CHOS, y=SO4)) + geom_point(size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "LC-MS CHOS (a.u. m-3)", y = "AMS Sulfate ug m-3") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+ geom_smooth(method = "lm", se=FALSE)
trend7<-ggplot(R22, aes(x=TotalOC, y=Org)) + geom_point(size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "LC-MS total OC (a.u. m-3)", y = "AMS Organic fraction ug m-3") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+ geom_smooth(method = "lm", se=FALSE)
multiplot(trend7, trend6, cols=2)

Trendsbio3 <- Samples[ which(Samples$Compound=='2-(2-carboxyethyl)-3,3-dimethylcyclobutanecarboxy acid (B-Caryophyllene SOA C10H16O4)'),]
Trendsbio31 <- Trendsbio3[which(Trendsbio3$datetime<"2016-11-17 23:59:00" ),]
Trendsbio31$num <- 1
Trendsbio32 <- Trendsbio3[which(Trendsbio3$datetime>="2016-11-17 23:59:00" ),]
Trendsbio32$num <- 2
Trendsbio3<-rbind(Trendsbio31,Trendsbio32)
Trendsbio4 <- Samples[ which(Samples$Compound=='C10H16O3 - RT3.92 (Pinene SOA)'),]
Trendsbio41 <- Trendsbio4[which(Trendsbio4$datetime<"2016-11-17 23:59:00" ),]
Trendsbio41$num <- 1
Trendsbio42 <- Trendsbio4[which(Trendsbio4$datetime>="2016-11-17 23:59:00" ),]
Trendsbio42$num <- 2
Trendsbio4<-rbind(Trendsbio41,Trendsbio42)
trendbiob<-ggplot(Trendsbio3, aes(x=datetime, y=Quant, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("C10H16O4 - RT4.75 (BSOA tracer)", "C10H16O3 - RT3.92 (BSOA tracer)"), values = c("darkblue")) + theme(legend.position="right") + labs(x = "Date", y = "Peak area per volume sampled (a.u. m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
trendbiob2<-ggplot(Trendsbio4, aes(x=datetime, y=Quant, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("C10H16O3 - RT3.92 (BSOA tracer)"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Peak area per volume sampled (a.u. m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio <- Samples[ which(Samples$Compound=='Octyl Sulfate'),]
Trendsbio1 <- Trendsbio[which(Trendsbio$datetime<"2016-11-17 23:59:00" ),]
Trendsbio1$num <- 1
Trendsbio2 <- Trendsbio[which(Trendsbio$datetime>="2016-11-17 23:59:00" ),]
Trendsbio2$num <- 2
Trendsbio<-rbind(Trendsbio1,Trendsbio2)
Trendsbio2 <- Samples[ which(Samples$Compound=='C10H7NO3 - RT8.25'),]
Trendsbio21 <- Trendsbio2[which(Trendsbio2$datetime<"2016-11-17 23:59:00" ),]
Trendsbio21$num <- 1
Trendsbio22 <- Trendsbio2[which(Trendsbio2$datetime>="2016-11-17 23:59:00" ),]
Trendsbio22$num <- 2
Trendsbio2<-rbind(Trendsbio21,Trendsbio22)
Trendsbio2$Conc <- 0
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio$Compound[x] == 'Octyl Sulfate') {Trendsbio$Conc[x] <- ((((Trendsbio$Area[x]-780972)/234125)*(8/3))/Trendsbio$Volume[x])}
for(x in 1:length(Trendsbio2$Compound))
  if (Trendsbio2$Compound[x] == 'C10H7NO3 - RT8.25') {Trendsbio2$Conc[x] <- ((((Trendsbio2$Area[x]-5291)/1321.2)*(8/3))/Trendsbio2$Volume[x])}
trendbio<-ggplot(Trendsbio, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("Octyl Sulfate"), values = c("darkblue")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
trendbio2<-ggplot(Trendsbio2, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-Nitro-1-Naphthol", "Octyl Sulfate"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
multiplot(trendbio2, trendbiob, trendbiob2, cols=1)
names(R22)[1] <- "date"
timeVariation(R22, pollutant = c("CHO","CHOS","CHON"), name.pol = c("Total CHO", "Total CHOS","Total CHON") , ylab = "Heterogroup totals (a.u./m3)")

#correlated just CHOS
Samples2 <- Samples
Samples2 <- Samples2 %>%
  select(Compound,
         Sample.ID,
         Quant,
         datetime) 

n = unique(Samples2[,2])

timeseries2 = seq(1478735940,1481241540,60)
timeseries3 = seq(1478735940,1481241540,60)
timeseries2 = data.frame(unixtime = timeseries2)
timeseries3 = data.frame(unixtime = timeseries3)
Samples2$unixtime = as.numeric(Samples2$datetime)
Samples2$Quant = as.numeric(Samples2$Quant)

for (i in 2:length(n)){
  temp_compound2 = Samples2[Samples2$Compound == n[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(n[i],"unixtime")
  
  timeseries2 = right_join(temp_compound2,timeseries2,"unixtime")
}
timeseries3 <- time %>% 
  left_join(timeseries2, by = "unixtime")

timeseries3[is.na(timeseries3)] <- 0
write.csv(timeseries3, file="winterfilterstimeseries.csv")

time<-HetGroupCHO %>% 
  select(datetime)
names(time)[1] <- "unixtime"

time$unixtime <- as.numeric(time$unixtime)
time2 <- merge(time,timeseries2,by="unixtime")

ROC <- SummedOC %>% 
  right_join(PM25, by = "datetime")
ROC <- ROC[which(ROC$datetime<="2016-12-04 23:59:00" ),]
ROC <- ROC[which(ROC$datetime>="2016-11-14 23:59:00" ),]
plotOC <- ggplot(ROC, aes(x = datetime))
plotOC <- plotOC + geom_line(aes(y = Quant, colour = "Total OC"))
plotOC <- plotOC + geom_line(aes(y = PM2.5*1000000, colour = "PM2.5"))
plotOC <- plotOC + scale_y_continuous(sec.axis = sec_axis(~.*0.000001, name = "PM2.5 (ug m-3)"))
plotOC <- plotOC + scale_colour_manual(values = c("blue", "red"))
plotOC <- plotOC + labs(y = "Peak area (a.u. m-3)",
                        x = "Date",
                        colour = "Parameter")
plotOC
timeseries2<-timeseries2[,-4]
timeseries2<-timeseries2[,-5]
time2 <- merge(time,timeseries2,by="unixtime")
time3 <- sulfate %>% 
  select(datetime)
time3$unixtime <- 0
time3$unixtime = as.numeric(time3$datetime)
time3 <- merge(time3,time2,by="unixtime")
time3 <- time3 %>% 
  select(-unixtime)
time3[(time3<2)] <- 0
write.csv(time3, file="winterfilters.csv")

ROC2 <- SummedOC %>% 
  right_join(PM25, by = "datetime")
ROC2 <- ROC2[which(ROC2$datetime<="2016-12-04 23:58:00" ),]
plotOC <- ggplot(ROC2, aes(x = datetime))
plotOC <- plotOC + geom_line(aes(y = Quant, colour = "Total OC"))
plotOC <- plotOC + geom_line(aes(y = PM2.5*1000000, colour = "PM2.5"))
plotOC <- plotOC + scale_y_continuous(sec.axis = sec_axis(~.*0.000001, name = "PM2.5 (ug m-3)"))
plotOC <- plotOC + scale_colour_manual(values = c("blue", "red"))
plotOC <- plotOC + labs(y = "Peak area (a.u. m-3)",
                        x = "Date",
                        colour = "Parameter")
plotOC
corPlot(ROC2, pollutants = c("PM2.5","Quant"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with ACSA", order="alphabet", cluster=FALSE) #Aerosol composition speciation analyser for PM2.5

R22 <- R22[which(R22$date<="2016-12-04 23:58:00" ),]
corPlot(R22, pollutants = c("Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS", order="alphabet", cluster=FALSE)
R2 <- R2[which(R2$datetime<="2016-12-04 23:58:00" ),]
corPlot(R2, pollutants = c("CHO","CHOS","CHON","Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS", order="alphabet", cluster=FALSE)
R2high <- R2high[which(R2high$date<="2016-12-04 23:58:00" ),]
corPlot(R2high, pollutants = c("Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS on polluted days", order="alphabet", cluster=FALSE)
R2low <- R2low[which(R2low$date<="2016-12-04 23:58:00" ),]
corPlot(R2low, pollutants = c("Org","SO4","NO3","BBOA","CCOA","LOOOA","MOOOA", "OPOA", "COA"), cols = "jet", dendrogram = FALSE, main = "Correlations of library factors with AMS on clean days", order="alphabet", cluster=FALSE)
# log transform 
R2high[is.na(R2high)] <- 1
log.r2h <- log(R2high[, 2:16])
r2.timeh <- R2high[, 1]
# apply PCA - scale. = TRUE is highly  GROUPS high PM
# advisable, but default is FALSE. 
r2.pcah <- prcomp(log.r2h,
                  center = TRUE,
                  scale. = TRUE) 
plot(r2.pcah, type = "l", main="high pollution days")
biplot(r2.pcah)
summary(r2.pcah)
# log transform 
R2low[is.na(R2low)] <- 1
log.r2l <- log(R2low[, 2:16])
r2.timel <- R2low[, 1]
# apply PCA - scale. = TRUE is highly  GROUPS low PM
# advisable, but default is FALSE. 
r2.pcal <- prcomp(log.r2l,
                  center = TRUE,
                  scale. = TRUE) 
plot(r2.pcal, type = "l", main="low pollution days")
biplot(r2.pcal)
summary(r2.pcal)
# log transform 
time3[is.na(time3)] <- 1
log.r22 <- log(time3[, 2:611])
r2.time2 <- time3[, 1]
# apply PCA - scale. = TRUE is highly   COMPOUNDS
# advisable, but default is FALSE. 
r2.pca2 <- prcomp(log.r22,
                  center = TRUE,
                  scale. = TRUE) 
plot(r2.pca2, type = "l")
biplot(r2.pca2)
summary(r2.pca2)
#correlated just CHON
orgnit <- Samples[which(Samples$Group == "CHON" ),]
orgnit <- orgnit %>%
  select(Compound,
         Sample.ID,
         Quant,
         datetime) 

n = unique(orgsulf[,1])

timeseries = seq(1478735940,1481241540,60)
timeseries = data.frame(unixtime = timeseries)
orgnit$unixtime = as.numeric(orgnit$datetime)

temp_compound = orgsulf[orgsulf$Compound == n[1],]
temp_compound = temp_compound[,c("Quant","unixtime")]
temp_compound = data.frame(temp_compound)
names(temp_compound) = c(n[1],"unixtime")
timeseries = left_join(timeseries,temp_compound,"unixtime")

for (i in 2:length(n)){
  temp_compound = orgsulf[orgsulf$Compound == n[i],]
  temp_compound = temp_compound[,c("Quant","unixtime")]
  temp_compound = data.frame(temp_compound)
  names(temp_compound) = c(n[i],"unixtime")
  
  timeseries = left_join(timeseries,temp_compound,"unixtime")
}
#possibly kill computer with 191 corplot of CHOS
filepath = paste(getwd(),"/CHOS_corplot.pdf",sep = "")
pdf(width=60, height=60, file=filepath)
corPlot(timeseries,
        pollutants = names(timeseries)[c(2:188)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors with AMS, CHOS",
        order="alphabet",
        cluster=TRUE)
dev.off()
#correlated just CHON
orgCHON <- Samples[which(Samples$Group == "CHON" ),]
orgCHON <- orgCHON %>%
  select(Compound,
         Sample.ID,
         Quant,
         datetime) 

n = unique(orgCHON[,1])

timeseries = seq(1478735940,1481241540,60)
timeseries = data.frame(unixtime = timeseries)
orgCHON$unixtime = as.numeric(orgCHON$datetime)

temp_compound = orgCHON[orgCHON$Compound == n[1],]
temp_compound = temp_compound[,c("Quant","unixtime")]
temp_compound = data.frame(temp_compound)
names(temp_compound) = c(n[1],"unixtime")
timeseries = left_join(timeseries,temp_compound,"unixtime")

for (i in 2:length(n)){
  temp_compound = orgCHON[orgCHON$Compound == n[i],]
  temp_compound = temp_compound[,c("Quant","unixtime")]
  temp_compound = data.frame(temp_compound)
  names(temp_compound) = c(n[i],"unixtime")
  
  timeseries = left_join(timeseries,temp_compound,"unixtime")
}

#possibly kill computer with 191 corplot of CHONS
filepath = paste(getwd(),"/CHON_corplot.pdf",sep = "")
pdf(width=60, height=60, file=filepath)
corPlot(timeseries,
        pollutants = names(timeseries)[c(2:41)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors with AMS, CHON",
        order="alphabet",
        cluster=TRUE)
dev.off()

timeseries[is.na(timeseries)] <- 0
colSums(timeseries != 0)
timeseries2 <- timeseries[, !colSums(timeseries != 0) < 20]
timeseries[timeseries < 0.1] <- NA
timeseries2[timeseries2 < 0.1] <- NA
timeseries2$'2-Methyl-3-nitrophenol' <- NULL
timeseries2$'C7H7NO3 - RT6.86' <- NULL
outCHON <- corPlot(timeseries2,
                   pollutants = names(timeseries2)[c(2:28)],
                   cols = "jet",
                   dendrogram = TRUE,
                   main = "Correlations of library factors with AMS, CHON",
                   order="alphabet",
                   cluster=TRUE)
mydendCHON <- as.dendrogram(outCHON$clust)
dend <- plot(as.phylo(outCHON$clust), cex=0.7)
dend <- color_labels(dend, k = 3)

timeseries_CHON_new<-timeseries2

ggplot(Avg, aes(Carbon, OSc, colour = RT, size = Quant)) + geom_point(aes(color=RT))+ scale_color_gradient(low="green3", high="red") + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE) + xlab("Carbon number") + ylab("OSc") + xlim(2, 20) + ylim(-2, 4)
ggplot(Avg, aes(Carbon, OCratio, colour = RT, size = Quant)) + geom_point(aes(color=RT))+ scale_color_gradient(low="green3", high="red") + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE) + xlab("Carbon number") + ylab("Oxygen : Carbon ratio") + xlim(2, 20) + ylim(0, 3)

Sulfcomps <- Avg[which(Avg$Group == "CHOS" ),]
SulfOSc <- ggplot(Sulfcomps, aes(Carbon, OSc, colour = RT, size = Quant)) + geom_point(aes(color=RT))+ scale_color_gradient(low="green3", high="red") + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE) +ggtitle("CHOS") + xlab("Carbon number") + ylab("OSc") + xlim(3, 20) + ylim(-2, 2.5)
Nitcomps <- Avg[which(Avg$Group == "CHON" ),]
NitOSc <- ggplot(Nitcomps, aes(Carbon, OSc, colour = RT, size = Quant)) + geom_point(aes(color=RT))+ scale_color_gradient(low="green3", high="red") + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE) +ggtitle("CHON") + xlab("Carbon number") + ylab("OSc") + xlim(3, 20) + ylim(-2, 2.5)
Orgcomps <- Avg[which(Avg$Group == "CHO" ),]
OrgOSc <- ggplot(Orgcomps, aes(Carbon, OSc, colour = RT, size = Quant)) + geom_point(aes(color=RT))+ scale_color_gradient(low="green3", high="red") + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE) +ggtitle("CHO") + xlab("Carbon number") + ylab("OSc") + xlim(3, 20) + ylim(-2, 2.5)
multiplot(SulfOSc, NitOSc, OrgOSc, cols=2)
ggplot(Avg, aes(OCratio, HCratio, colour = Group, size = Quant)) + geom_point()+ theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE)  + xlab("O : C ratio") + ylab("H : C ratio") + xlim(0, 3) + ylim(0, 3) + geom_hline(yintercept = 2, color="red") + geom_abline(intercept = 2, slope = -1, color="red") + geom_abline(intercept = 2, slope = -2, color="red") + geom_abline(intercept = 2, slope = 2, color="blue")


R2high2 <- R2high
R2high2$unixtime = as.numeric(R2high2$datetime)
R2high2 <- R2high2 %>%
  select(unixtime) 

n = unique(orgsulf[,1])
temp_compound = orgsulf[orgsulf$Compound == n[1],]
temp_compound = temp_compound[,c("Quant","unixtime")]
temp_compound = data.frame(temp_compound)
names(temp_compound) = c(n[1],"unixtime")
R2high2 = left_join(R2high2,temp_compound,"unixtime")

for (i in 2:length(n)){
  temp_compound = orgsulf[orgsulf$Compound == n[i],]
  temp_compound = temp_compound[,c("Quant","unixtime")]
  temp_compound = data.frame(temp_compound)
  names(temp_compound) = c(n[i],"unixtime")
  
  R2high2 = left_join(R2high2,temp_compound,"unixtime")
}
#possibly kill computer with 191 corplot of CHOS
filepath = paste(getwd(),"/CHOS_corplothigh.pdf",sep = "")
pdf(width=60, height=60, file=filepath)
corPlot(R2high2,
        pollutants = names(R2high2)[c(2:191)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors with high pollution AMS, CHOS",
        order="alphabet",
        cluster=TRUE)
dev.off()

R2low2 <- R2low
R2low2$unixtime = as.numeric(R2low2$datetime)
R2low2 <- R2low2 %>%
  select(unixtime) 

n = unique(orgsulf[,1])
temp_compound = orgsulf[orgsulf$Compound == n[1],]
temp_compound = temp_compound[,c("Quant","unixtime")]
temp_compound = data.frame(temp_compound)
names(temp_compound) = c(n[1],"unixtime")
R2low2 = left_join(R2low2,temp_compound,"unixtime")

for (i in 2:length(n)){
  temp_compound = orgsulf[orgsulf$Compound == n[i],]
  temp_compound = temp_compound[,c("Quant","unixtime")]
  temp_compound = data.frame(temp_compound)
  names(temp_compound) = c(n[i],"unixtime")
  
  R2low2 = left_join(R2low2,temp_compound,"unixtime")
}
#possibly kill computer with 191 corplot of CHOS
filepath = paste(getwd(),"/CHOS_corplotlow.pdf",sep = "")
pdf(width=60, height=60, file=filepath)
corPlot(R2low2,
        pollutants = names(R2low2)[c(2:191)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors with low pollution AMS, CHOS",
        order="alphabet",
        cluster=TRUE,
        numbers=FALSE)
dev.off()

Trendsbio3 <- Samples[ which(Samples$Compound=='2,4-Dinitrophenol'),]
Trendsbio31 <- Trendsbio3[which(Trendsbio3$datetime<"2016-11-17 23:59:00" ),]
Trendsbio31$num <- 1
Trendsbio32 <- Trendsbio3[which(Trendsbio3$datetime>="2016-11-17 23:59:00" ),]
Trendsbio32$num <- 2
Trendsbio3<-rbind(Trendsbio31,Trendsbio32)
Trendsbio4 <- Samples[ which(Samples$Compound=='2-Methyl-3-nitrophenol'),]
Trendsbio41 <- Trendsbio4[which(Trendsbio4$datetime<"2016-11-17 23:59:00" ),]
Trendsbio41$num <- 1
Trendsbio42 <- Trendsbio4[which(Trendsbio4$datetime>="2016-11-17 23:59:00" ),]
Trendsbio42$num <- 2
Trendsbio4<-rbind(Trendsbio41,Trendsbio42)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio4$Compound[x] == '2-Methyl-3-nitrophenol') {Trendsbio4$Conc[x] <- ((((Trendsbio4$Area[x]-298958)/172182)*(8/3))/Trendsbio4$Volume[x])}
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio3$Compound[x] == '2,4-Dinitrophenol') {Trendsbio3$Conc[x] <- ((((Trendsbio3$Area[x]-97500)/85439)*(8/3))/Trendsbio3$Volume[x])}
trendbiob<-ggplot(Trendsbio3, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2,4-Dinitrophenol", "2-Methyl-3-nitrophenol"), values = c("green3")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
trendbiob2<-ggplot(Trendsbio4, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-Methyl-3-nitrophenol"), values = c("purple3")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio <- Samples[ which(Samples$Compound=='Octyl Sulfate'),]
Trendsbio1 <- Trendsbio[which(Trendsbio$datetime<"2016-11-17 23:59:00" ),]
Trendsbio1$num <- 1
Trendsbio2 <- Trendsbio[which(Trendsbio$datetime>="2016-11-17 23:59:00" ),]
Trendsbio2$num <- 2
Trendsbio<-rbind(Trendsbio1,Trendsbio2)
Trendsbio2 <- Samples[ which(Samples$Compound=='C10H7NO3 - RT8.25'),]
Trendsbio21 <- Trendsbio2[which(Trendsbio2$datetime<"2016-11-17 23:59:00" ),]
Trendsbio21$num <- 1
Trendsbio22 <- Trendsbio2[which(Trendsbio2$datetime>="2016-11-17 23:59:00" ),]
Trendsbio22$num <- 2
Trendsbio2<-rbind(Trendsbio21,Trendsbio22)
Trendsbio2$Conc <- 0
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio$Compound[x] == 'Octyl Sulfate') {Trendsbio$Conc[x] <- ((((Trendsbio$Area[x]-780972)/234125)*(8/3))/Trendsbio$Volume[x])}
for(x in 1:length(Trendsbio2$Compound))
  if (Trendsbio2$Compound[x] == 'C10H7NO3 - RT8.25') {Trendsbio2$Conc[x] <- ((((Trendsbio2$Area[x]-4669)/1790.6)*(8/3))/Trendsbio2$Volume[x])}
trendbio<-ggplot(Trendsbio, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("Octyl Sulfate"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
trendbio2<-ggplot(Trendsbio2, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-Nitro-1-Naphthol", "Octyl Sulfate"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio5 <- Samples[ which(Samples$Compound=='2,6-dimethyl-4-nitrophenol'),]
Trendsbio51 <- Trendsbio5[which(Trendsbio5$datetime<"2016-11-17 23:59:00" ),]
Trendsbio51$num <- 1
Trendsbio52 <- Trendsbio5[which(Trendsbio5$datetime>="2016-11-17 23:59:00" ),]
Trendsbio52$num <- 2
Trendsbio5<-rbind(Trendsbio51,Trendsbio52)
Trendsbio6 <- Samples[ which(Samples$Compound=='4-Nitrophenol'),]
Trendsbio61 <- Trendsbio6[which(Trendsbio6$datetime<"2016-11-17 23:59:00" ),]
Trendsbio61$num <- 1
Trendsbio62 <- Trendsbio6[which(Trendsbio6$datetime>="2016-11-17 23:59:00" ),]
Trendsbio62$num <- 2
Trendsbio6<-rbind(Trendsbio61,Trendsbio62)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio6$Compound[x] == '4-Nitrophenol') {Trendsbio6$Conc[x] <- ((((Trendsbio6$Area[x]-538611)/112702)*(8/3))/Trendsbio6$Volume[x])}
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio5$Compound[x] == '2,6-dimethyl-4-nitrophenol') {Trendsbio5$Conc[x] <- ((((Trendsbio5$Area[x]-690764)/230716)*(8/3))/Trendsbio5$Volume[x])}
trendbioc<-ggplot(Trendsbio5, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2,6-dimethyl-4-nitrophenol", "4-Nitrophenol"), values = c("green")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
trendbioc2<-ggplot(Trendsbio6, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("4-Nitrophenol"), values = c("purple")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio7 <- Samples[ which(Samples$Compound=='2-Nitrophenol'),]
Trendsbio71 <- Trendsbio7[which(Trendsbio7$datetime<"2016-11-17 23:59:00" ),]
Trendsbio71$num <- 1
Trendsbio72 <- Trendsbio7[which(Trendsbio7$datetime>="2016-11-17 23:59:00" ),]
Trendsbio72$num <- 2
Trendsbio7<-rbind(Trendsbio71,Trendsbio72)
Trendsbio8 <- Samples[ which(Samples$Compound=='3-Methyl-2-nitrophenol'),]
Trendsbio81 <- Trendsbio8[which(Trendsbio8$datetime<"2016-11-17 23:59:00" ),]
Trendsbio81$num <- 1
Trendsbio82 <- Trendsbio8[which(Trendsbio8$datetime>="2016-11-17 23:59:00" ),]
Trendsbio82$num <- 2
Trendsbio8<-rbind(Trendsbio81,Trendsbio82)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio8$Compound[x] == '3-Methyl-2-nitrophenol') {Trendsbio8$Conc[x] <- ((((Trendsbio8$Area[x]-298958)/172182)*(8/3))/Trendsbio8$Volume[x])}
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio7$Compound[x] == '2-Nitrophenol') {Trendsbio7$Conc[x] <- ((((Trendsbio7$Area[x]-20120)/1881)*(8/3))/Trendsbio7$Volume[x])}
trendbiod<-ggplot(Trendsbio7, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-Nitrophenol", "3-Methyl-2-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
trendbiod2<-ggplot(Trendsbio8, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("3-Methyl-2-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio9 <- Samples[ which(Samples$Compound=='3-Nitrophenol'),]
Trendsbio91 <- Trendsbio9[which(Trendsbio9$datetime<"2016-11-17 23:59:00" ),]
Trendsbio91$num <- 1
Trendsbio92 <- Trendsbio9[which(Trendsbio9$datetime>="2016-11-17 23:59:00" ),]
Trendsbio92$num <- 2
Trendsbio9<-rbind(Trendsbio91,Trendsbio92)
Trendsbio10 <- Samples[ which(Samples$Compound=='2-Methyl-5-nitrophenol'),]
Trendsbio101 <- Trendsbio10[which(Trendsbio10$datetime<"2016-11-17 23:59:00" ),]
Trendsbio101$num <- 1
Trendsbio102 <- Trendsbio10[which(Trendsbio10$datetime>="2016-11-17 23:59:00" ),]
Trendsbio102$num <- 2
Trendsbio10<-rbind(Trendsbio101,Trendsbio102)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio10$Compound[x] == '2-Methyl-5-nitrophenol') {Trendsbio10$Conc[x] <- ((((Trendsbio10$Area[x]-286389)/126443)*(8/3))/Trendsbio10$Volume[x])}
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio9$Compound[x] == '3-Nitrophenol') {Trendsbio9$Conc[x] <- ((((Trendsbio9$Area[x]-190625)/92675)*(8/3))/Trendsbio9$Volume[x])}
trendbioe<-ggplot(Trendsbio9, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("3-Nitrophenol", "2-Methyl-5-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
trendbioe2<-ggplot(Trendsbio10, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-Methyl-5-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio12 <- Samples[ which(Samples$Compound=='4-Methyl-3-nitrophenol'),]
Trendsbio121 <- Trendsbio12[which(Trendsbio12$datetime<"2016-11-17 23:59:00" ),]
Trendsbio121$num <- 1
Trendsbio122 <- Trendsbio12[which(Trendsbio12$datetime>="2016-11-17 23:59:00" ),]
Trendsbio122$num <- 2
Trendsbio12<-rbind(Trendsbio121,Trendsbio122)
Trendsbio11 <- Samples[ which(Samples$Compound=='cis-pinonic acid (pinene SOA)'),]
Trendsbio111 <- Trendsbio11[which(Trendsbio11$datetime<"2016-11-17 23:59:00" ),]
Trendsbio111$num <- 1
Trendsbio112 <- Trendsbio11[which(Trendsbio11$datetime>="2016-11-17 23:59:00" ),]
Trendsbio112$num <- 2
Trendsbio11<-rbind(Trendsbio111,Trendsbio112)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio12$Compound[x] == '4-Methyl-3-nitrophenol') {Trendsbio12$Conc[x] <- ((((Trendsbio12$Area[x]-298958)/172182)*(8/3))/Trendsbio12$Volume[x])}
trendbiof2<-ggplot(Trendsbio12, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("4-Methyl-3-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio11$Compound[x] == 'cis-pinonic acid (pinene SOA)') {Trendsbio11$Conc[x] <- ((((Trendsbio11$Area[x]-97500)/85439)*(8/3))/Trendsbio11$Volume[x])}
trendbiof<-ggplot(Trendsbio11, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("cis-pinonic acid (pinene SOA)"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio13 <- Samples[ which(Samples$Compound=='2-Napthyl Sulfate'),]
Trendsbio131 <- Trendsbio13[which(Trendsbio13$datetime<"2016-11-17 23:59:00" ),]
Trendsbio131$num <- 1
Trendsbio132 <- Trendsbio13[which(Trendsbio13$datetime>="2016-11-17 23:59:00" ),]
Trendsbio132$num <- 2
Trendsbio13<-rbind(Trendsbio131,Trendsbio132)
Trendsbio14 <- Samples[ which(Samples$Compound=='2-(2-carboxyethyl)-3,3-dimethylcyclobutanecarboxy acid (B-Caryophyllene SOA C10H16O4)'),]
Trendsbio141 <- Trendsbio14[which(Trendsbio14$datetime<"2016-11-17 23:59:00" ),]
Trendsbio141$num <- 1
Trendsbio142 <- Trendsbio14[which(Trendsbio14$datetime>="2016-11-17 23:59:00" ),]
Trendsbio142$num <- 2
Trendsbio14<-rbind(Trendsbio141,Trendsbio142)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio13$Compound[x] == '2-Napthyl Sulfate') {Trendsbio13$Conc[x] <- ((((Trendsbio13$Area[x]-10822)/2683)*(8/3))/Trendsbio13$Volume[x])}
trendbiog2<-ggplot(Trendsbio13, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-Napthyl Sulfate"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio14$Compound[x] == '2-(2-carboxyethyl)-3,3-dimethylcyclobutanecarboxy acid (B-Caryophyllene SOA C10H16O4)') {Trendsbio14$Conc[x] <- ((((Trendsbio14$Area[x]-90833)/250508)*(8/3))/Trendsbio14$Volume[x])}
trendbiog<-ggplot(Trendsbio14, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("B-Caryophyllene SOA C10H16O4"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

Trendsbio15 <- Samples[ which(Samples$Compound=='Dodecyl Sulfate'),]
Trendsbio151 <- Trendsbio15[which(Trendsbio15$datetime<"2016-11-17 23:59:00" ),]
Trendsbio151$num <- 1
Trendsbio152 <- Trendsbio15[which(Trendsbio15$datetime>="2016-11-17 23:59:00" ),]
Trendsbio152$num <- 2
Trendsbio15<-rbind(Trendsbio151,Trendsbio152)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio15$Compound[x] == 'Dodecyl Sulfate') {Trendsbio15$Conc[x] <- ((((Trendsbio15$Area[x]-872083)/189985)*(8/3))/Trendsbio15$Volume[x])}
trendbioh2<-ggplot(Trendsbio15, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("Dodecyl Sulfate"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
Trendsbio16 <- Samples[ which(Samples$Compound=='4-Methoxy-2-nitrophenol'),]
Trendsbio161 <- Trendsbio16[which(Trendsbio16$datetime<"2016-11-17 23:59:00" ),]
Trendsbio161$num <- 1
Trendsbio162 <- Trendsbio16[which(Trendsbio16$datetime>="2016-11-17 23:59:00" ),]
Trendsbio162$num <- 2
Trendsbio16<-rbind(Trendsbio161,Trendsbio162)
for(x in 1:length(Trendsbio$Compound))
  if (Trendsbio16$Compound[x] == '4-Methoxy-2-nitrophenol') {Trendsbio16$Conc[x] <- ((((Trendsbio16$Area[x]-359)/9564)*(8/3))/Trendsbio16$Volume[x])}
trendbioh<-ggplot(Trendsbio16, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("4-Methoxy-2-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))

multiplot(trendbio, trendbio2, trendbiob, trendbiob2, trendbioc, trendbioc2, trendbiod, trendbiod2, trendbioe2, trendbiof2, trendbiof, trendbiog2, trendbiog, trendbioh2, trendbioh, cols=4)
multiplot(trendbio2, trendbiob, trendbiob2, trendbioc, trendbioc2, trendbiod, trendbiod2, trendbioe2, trendbiof2, trendbiog2, trendbioh, cols=4)
multiplot(trendbio2, trendbiob, trendbiob2, trendbioc, trendbioc2, trendbiod2, trendbiof2, trendbioh, cols=3)
multiplot(trendbio, trendbiog2, cols=1)

STDcorr <- Trendsbio2 %>%
  select(datetime,
         Conc) 
colnames(STDcorr)[2] <- "2-Nitro-1-naphthol"
STDcorr <- STDcorr %>% 
  right_join(Trendsbio, by = "datetime")%>%
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[3] <- "OctylSulfate" 
STDcorr <- STDcorr %>% 
  right_join(Trendsbio3, by = "datetime")%>%
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[4] <- "2,4-Dinitrophenol"
STDcorr <- STDcorr%>% 
  right_join(Trendsbio5, by = "datetime")%>%
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[5] <- "2,6-dimethyl-4-nitrophenol"
STDcorr <- STDcorr%>% 
  right_join(Trendsbio6, by = "datetime")%>%
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[6] <- "4-Nitrophenol"
STDcorr <- STDcorr%>%
  right_join(Trendsbio8, by = "datetime")%>% 
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[7] <- "3-Methyl-2-nitrophenol"
STDcorr <- STDcorr%>%
  right_join(Trendsbio12, by = "datetime")%>% 
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[8] <- "4-Methyl-3-nitrophenol"
STDcorr <- STDcorr%>%
  right_join(Trendsbio11, by = "datetime")%>% 
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[9] <- "cis-pinonicacid"
STDcorr <- STDcorr%>%
  right_join(Trendsbio13, by = "datetime")%>% 
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr)[10] <- "2-NapthylSulfate"
STDcorr2 <- STDcorr%>%
  right_join(Trendsbio15, by = "datetime")%>% 
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr2)[11] <- "DodecylSulfate"
STDcorr2 <- STDcorr2%>%
  right_join(Trendsbio16, by = "datetime")%>% 
  select(-num,
         -OSc,
         -OCratio,
         -HCratio,
         -AI,
         -DBE,
         -Group,
         -Quant,
         -Sulphur,
         -Carbon,
         -Nitrogen,
         -Oxygen,
         -Hydrogen,
         -Area,
         -Sample.ID,
         -Compound,
         -Expected.RT,
         -RT,
         -Volume)
colnames(STDcorr2)[12] <- "4-Methoxy-2-nitrophenol"

corPlot(STDcorr,
        pollutants = names(STDcorr)[c(2:13)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of nitrophenol compounds",
        order="alphabet",
        cluster=TRUE)

STDcorr3<-STDcorr[c(1,5)]
STDcorr4<-STDcorr2[c(1,12)]
STDcorr3 <- STDcorr3%>%
  left_join(STDcorr4, by = "datetime")
colnames(STDcorr3)[2] <- "dimethyl"
colnames(STDcorr3)[3] <- "Methoxy"

plotTWO2 <- ggplot(STDcorr3, aes(datetime, dimethyl)) + 
  geom_line(colour = "#68382C", size = 1.5) + 
  ggtitle("2,6-dimethyl-4-nitrophenol") +
  labs(x=NULL,y=NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "gray50", size = 0.5),
    panel.grid.major.x = element_blank(),
    text = element_text(family="ITCOfficinaSans LT Book"),
    axis.text.y = element_text(colour="#68382C", size = 14),
    axis.text.x = element_text(size = 14),
    axis.ticks = element_line(colour = 'gray50'),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = -0.135, vjust=2.12, colour="#68382C", size = 14, family = "OfficinaSanITCMedium")) 
plotTWO <- ggplot(STDcorr3, aes(datetime, Methoxy)) + 
  geom_line(colour = "#00A4E6", size = 1.5) +  
  ggtitle("4-Methoxy-2-nitrophenol") +
  labs(x=NULL,y=NULL) +
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "gray50", size = 0.5),
    panel.grid.major.x = element_blank(),
    text = element_text(family="ITCOfficinaSans LT Book"),
    axis.text.y = element_text(colour="#00A4E6", size=14),
    axis.text.x = element_text(size = 14),
    axis.ticks = element_line(colour = 'gray50'),
    axis.ticks.length = unit(.25, "cm"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.85, vjust=2.12, colour = "#00a4e6", size = 14, family = "OfficinaSanITCMedium"))

plotTWO3 <- ggplot(STDcorr3, aes(x=datetime))
plotTWO3 <- plotTWO3 + geom_line(aes(y = Methoxy/5, colour = "4-Methoxy-2-nitrophenol"))
plotTWO3 <- plotTWO3 + geom_line(aes(y = Naphth, colour = "2,6-dimethyl-4-nitrophenol"))
plotTWO3 <- plotTWO3 + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "4-Methoxy-2-nitrophenol Concentration ng/m3"))
plotTWO3 <- plotTWO3 + scale_colour_manual(values = c("blue", "red"))
plotTWO3 <- plotTWO3 + labs(y = "2,6-dimethyl-4-nitrophenol Concentration ng/m3",
                            x = "Date")
plotTWO3

STDcorr5 <- STDcorr3%>%
  left_join(AMS, by = "datetime")
corPlot(STDcorr5,
        pollutants = names(STDcorr5)[c(2:14)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of nitrophenol compounds",
        order="alphabet",
        cluster=TRUE)

plotTWO4 <- ggplot(STDcorr5, aes(x=datetime))
plotTWO4 <- plotTWO4 + geom_line(aes(y = Methoxy/3, colour = "4-Methoxy-2-nitrophenol"))
plotTWO4 <- plotTWO4 + geom_line(aes(y = BBOA, colour = "Biomass burning"))
plotTWO4 <- plotTWO4 + scale_y_continuous(sec.axis = sec_axis(~.*3, name = "4-Methoxy-2-nitrophenol Concentration ng/m3"))
plotTWO4 <- plotTWO4 + scale_colour_manual(values = c("blue", "red"))
plotTWO4 <- plotTWO4 + labs(y = "Biomass burning ug/m3",
                            x = "Date")
plotTWO4

plotTWO5 <- ggplot(STDcorr5, aes(x=datetime))
plotTWO5 <- plotTWO5 + geom_line(aes(y = dimethyl/2, colour = "2,6-dimethyl-4-nitrophenol"))
plotTWO5 <- plotTWO5 + geom_line(aes(y = OPOA, colour = "Oxygenated primary emissions"))
plotTWO5 <- plotTWO5 + scale_y_continuous(sec.axis = sec_axis(~.*2, name = "2,6-dimethyl-4-nitrophenol Concentration ng/m3"))
plotTWO5 <- plotTWO5 + scale_colour_manual(values = c("blue", "red"))
plotTWO5 <- plotTWO5 + labs(y = "OPOA ug/m3",
                            x = "Date")
plotTWO5

Box1 <- Trendsbio %>%
  select(Compound,
         Conc) 
Box2 <- Trendsbio2 %>%
  select(Compound,
         Conc) 
Box3 <- Trendsbio3 %>%
  select(Compound,
         Conc) 
Box4 <- Trendsbio4 %>%
  select(Compound,
         Conc) 
Box5 <- Trendsbio5 %>%
  select(Compound,
         Conc) 
Box6 <- Trendsbio6 %>%
  select(Compound,
         Conc) 
Box7 <- Trendsbio7 %>%
  select(Compound,
         Conc)
Box8 <- Trendsbio8 %>%
  select(Compound,
         Conc) 
Box9 <- Trendsbio9 %>%
  select(Compound,
         Conc) 
Box10 <- Trendsbio10 %>%
  select(Compound,
         Conc) 
Box11 <- Trendsbio11 %>%
  select(Compound,
         Conc) 
Box12 <- Trendsbio12 %>%
  select(Compound,
         Conc) 
Box13 <- Trendsbio13 %>%
  select(Compound,
         Conc)
Box14 <- Trendsbio14 %>%
  select(Compound,
         Conc)
Box15 <- Trendsbio15 %>%
  select(Compound,
         Conc) 
Box16 <- Trendsbio16 %>%
  select(Compound,
         Conc) 
Box <- rbind(Box1, Box3, Box4, Box5, Box6, Box7, Box8, Box9, Box10, Box11, Box12, Box13, Box15, Box16)
plotBox <- ggplot(Box, aes(x = Compound, y = Conc)) +
  geom_boxplot(alpha = 0.7,
               outlier.colour = "#1F3552", outlier.shape = 20) +
  scale_y_continuous(name = "Concentration",
                     breaks = seq(0, 250, 20),
                     limits=c(0, 160)) +
  scale_x_discrete(name = "Compound") +
  ggtitle("Boxplot of mean compound concentration")
plotBox

Daynitro <- Day27[grepl("itro", Day27$Compound),]
Daysulf <- Day27[grepl("CHOS", Day27$Group),]
DayCHONS <- Day27[grepl("CHONS", Day27$Group),]
sum(Daynitro$Quant)
sum(Daysulf$Quant)
sum(DayCHONS$Quant)

time4 <- Trendsbio2
hourly<-ggplot(time4, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-nitro-1-naphthol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
time5 <- Trendsbio3
hourly2<-ggplot(time5, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2,4-dinitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
time6 <- Trendsbio4
hourly3<-ggplot(time6, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-methyl-3-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) 
time7 <- Trendsbio16
hourly4<-ggplot(time7, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("4-Methoxy-2-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) 
multiplot(hourly, hourly2, hourly3, hourly4, cols=2)

slices <- c(80.42, 12.56, 1.12, 5.88) 
lbls <- c("CHON 80.42%", "CHOS 12.56%", "CHONS 1.12%", "CHO 5.88%")
pie3D(slices,labels=lbls,explode=0.1,
      main="Proportion of groups in averaged aerosol sample")

HetGroup2 <- read.csv("C:/Users/William/Documents/WACL stuff/HetGroup.csv", na.strings="N/F")
HetGroup2$datetime <- HetGroup$datetime
ggplot(HetGroup2, aes(fill=Group, y=Quant, x=datetime)) + 
  geom_area(color='black', size=.2, alpha=.4) + labs(x = "Date", y = "% Composition")

Averages <- data.frame(datetime=as.Date(character()),
                       OCratio=numeric(), 
                       HCratio=numeric(), 
                       DBE=numeric(), 
                       OSc=numeric(), 
                       RT=numeric(), 
                       stringsAsFactors=FALSE) 
Averages2 <- data.frame(datetime="0",
                        OCratio="0", 
                        HCratio="0", 
                        DBE="0", 
                        OSc="0", 
                        RT="0", 
                        stringsAsFactors=FALSE) 
for(i in 1:length(unique(Samples$datetime))){
  this.a <- unique(Samples$datetime)[i]
  hold <- Samples[ which(Samples$datetime == this.a),]
  hold[is.na(hold)] <- 0
  avgOC <- mean(hold$OCratio)
  avgHC <- mean(hold$HCratio)
  avgOSc <- mean(hold$OSc)
  avgDBE <- mean(hold$DBE)
  avgRT <- mean(hold$RT)
  
  Averages2$OCratio <- rbind(weighted.mean(hold$OCratio, hold$Area, na.rm=TRUE))
  Averages2$HCratio <- rbind(weighted.mean(hold$HCratio, hold$Area, na.rm=TRUE))
  Averages2$OSc <- rbind(weighted.mean(hold$OSc, hold$Area, na.rm=TRUE))
  Averages2$DBE <- rbind(weighted.mean(hold$DBE, hold$Area, na.rm=TRUE))
  Averages2$RT <- rbind(weighted.mean(hold$RT, hold$Area, na.rm=TRUE))
  Averages <- rbind(Averages, Averages2)
}
Averages$datetime<-unique(Samples$datetime)

avgOC<-ggplot(Averages,aes(x=datetime,y=OCratio))+geom_line(color="red", size=1.1) +theme_bw()+labs(y="O : C ratio", x="Date")+coord_cartesian(ylim=c(0.1,0.8))
avgHC<-ggplot(Averages,aes(x=datetime,y=HCratio))+geom_line(color="red", size=1.1) +theme_bw()+labs(y="H : C ratio", x="Date")+coord_cartesian(ylim=c(0.2,1.4))
avgOSc<-ggplot(Averages,aes(x=datetime,y=OSc))+geom_line(color="darkgreen", size=1.1) +theme_bw()+labs(y="OSc", x="Date")+coord_cartesian(ylim=c(-0.7,-0.1))
avgDBE<-ggplot(Averages,aes(x=datetime,y=DBE))+geom_line(color="blue", size=1.1) +theme_bw()+labs(y="DBE", x="Date")+coord_cartesian(ylim=c(0,4))
avgRT<-ggplot(Averages,aes(x=datetime,y=RT))+geom_line(color="purple", size=1.1) +theme_bw()+labs(y="Retention time", x="Date")+coord_cartesian(ylim=c(1,6))
multiplot(avgOC, avgHC, avgOSc, avgDBE, avgRT, cols=1)

winteravgOC<-Averages

colnames(winteravgOC)[1] <- "datetime"

winteravgOC$datetime <- ymd_hms(winteravgOC$datetime, tz = "UTC")
HetGroupCHO$datetime <- ymd_hms(HetGroupCHO$datetime, tz = "UTC")

winteravgOC2 <- winteravgOC %>%  
  select(datetime,
         OCratio) 

winteravgOC2$OCratio<-as.numeric(winteravgOC2$OCratio)

winteravgOC2 <- winteravgOC2 %>%
  right_join(HetGroupCHO, by = "datetime")

winterAVG_OC<-weighted.mean(winteravgOC2$OCratio, winteravgOC2$Quant, na.rm=TRUE)
summerAVG_OC<-weighted.mean(summeravgOC2$OCratio, summeravgOC2$Quant, na.rm=TRUE)

scatter<-AMS
scatter <- scatter %>%
  right_join(HetGroupCHO, by = "datetime")%>%
  right_join(HetGroupCHON, by = "datetime")%>%
  right_join(HetGroupCHOS, by = "datetime")%>%
  right_join(SummedOC, by = "datetime")%>%
  right_join(Trendsbio16, by = "datetime")
lm_fit <- lm(Quant ~ Org, data=scatter)
summary(lm_fit)
predicted_df <- data.frame(mpg_pred = predict(lm_fit, df))
OCCHOS<-ggplot(scatter,aes(x=Quant,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon", x="CHOS")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
OCCHO<-ggplot(scatter,aes(x=Quant.x,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon", x="CHO")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
OCCHON<-ggplot(scatter,aes(x=Quant.y,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon", x="CHON")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
BBOACHOS<-ggplot(scatter,aes(x=Quant,y=BBOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Biomass burning organic aerosol", x="CHOS")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
BBOACHO<-ggplot(scatter,aes(x=Quant.x,y=BBOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Biomass burning organic aerosol", x="CHO")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
BBOACHON<-ggplot(scatter,aes(x=Quant.y,y=BBOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Biomass burning organic aerosol", x="CHON")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
OPOACHOS<-ggplot(scatter,aes(x=Quant,y=OPOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Oxygenated primary organic aerosol", x="CHOS")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
OPOACHO<-ggplot(scatter,aes(x=Quant.x,y=OPOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Oxygenated primary organic aerosol", x="CHO")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
OPOACHON<-ggplot(scatter,aes(x=Quant.y,y=OPOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Oxygenated primary organic aerosol", x="CHON")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
multiplot(OCCHO,OCCHON,OCCHOS,BBOACHO,BBOACHON,BBOACHOS,OPOACHO,OPOACHON,OPOACHOS, cols=3)

totalorgAMS<-ggplot(scatter,aes(x=Quant.y.y,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon (AMS) ug/m3", x="Total organic carbon peak area (Library)")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 

lm_eqn <- function(df){
  m <- lm(Org ~ Quant.y.y, scatter);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
ggplot(scatter,aes(x=Quant.y.y,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon (AMS) ug/m3", x="Total organic carbon peak area (Library)")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_text(x = 50000000, y = 100, label = lm_eqn(scatter), parse = TRUE)

scatter3 <- scatter[which(scatter$'Quant.y.y'>20000000 ),]
scatter3 <- scatter3[which(scatter3$Org<100 ),]
lm_eqn <- function(df){
  m <- lm(Org ~ Quant.y.y, scatter3);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq));                 
}
ggplot(scatter3,aes(x=Quant.y.y,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon (AMS) ug/m3", x="Total organic carbon peak area (Library)")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_text(x = 50000000, y = 100, label = lm_eqn(scatter3), parse = TRUE)

Samples3 <- Samples2

OStrends <- Samples3[which(Samples3$Compound == "C11H12O7S - RT1.41" ),]
OStrends2 <- Samples3[which(Samples3$Compound == "C12H20O6SO4 - RT9.82" ),]
OStrends3 <- Samples3[which(Samples3$Compound == "Octyl Sulfate" ),]
OStrends4 <- Samples3[which(Samples3$Compound == "C10H17NO3SO4 - RT6.98" ),]
OS<-ggplot(OStrends,aes(x=datetime,y=Quant))+geom_line(color="red") +theme_bw()+labs(y="C11H12O7S - RT1.41 Peak area", x="Date")
OS2<-ggplot(OStrends2,aes(x=datetime,y=Quant))+geom_line(color="red") +theme_bw()+labs(y="C12H20O6SO4 - RT9.82 Peak area", x="Date")
OS3<-ggplot(OStrends3,aes(x=datetime,y=Quant))+geom_line(color="red") +theme_bw()+labs(y="Octyl Sulfate Peak area", x="Date")
OS4<-ggplot(OStrends4,aes(x=datetime,y=Quant))+geom_line(color="red") +theme_bw()+labs(y="C10H17NO3SO4 - RT6.98 Peak area", x="Date")
multiplot(OS3, OS2, OS, OS4, cols=2)
ggplot(scatter,aes(x=Quant.x,y=MOOOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="MOOOA", x="CHO")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggplot(scatter,aes(x=Conc,y=BBOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Biomass Burning organic aerosol  ug/m3", x="4-Methoxy-2-nitrophenol ng/m3")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

scatter2 <- scatter[ which(scatter$datetime != "2016-11-29 13:04:00"),]
ggplot(scatter2,aes(x=Conc,y=BBOA))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Biomass Burning organic aerosol  ug/m3", x="4-Methoxy-2-nitrophenol ng/m3")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

STDcorr3 <- STDcorr2 %>%
  right_join(AMS, by = "datetime")
STDcorr4 <- STDcorr3%>%
  right_join(Trendsbio7, by = "datetime")
STDcorr4 <- STDcorr4%>%  
  select(datetime,
         Conc)
STDcorr3 <- STDcorr3%>%
  left_join(STDcorr4, by = "datetime")
STDcorr5 <- STDcorr2%>%
  right_join(Trendsbio4, by = "datetime")
STDcorr5 <- STDcorr5%>%  
  select(datetime,
         Conc)
STDcorr3 <- STDcorr3%>%
  left_join(STDcorr5, by = "datetime")
STDcorr6 <- STDcorr2%>%
  right_join(Trendsbio10, by = "datetime")
STDcorr6 <- STDcorr6%>%  
  select(datetime,
         Conc)
STDcorr3 <- STDcorr3%>%
  left_join(STDcorr6, by = "datetime")

for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc[x] <- ((((STDcorr3$`2,4-Dinitrophenol`[x]/1000))/STDcorr3$Org[x])*100)}
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc2[x] <- ((((STDcorr3$`2-Nitro-1-naphthol`[x]/1000))/STDcorr3$Org[x])*100)}
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc3[x] <- ((((STDcorr3$`2,6-dimethyl-4-nitrophenol`[x]/1000))/STDcorr3$Org[x])*100)}
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc4[x] <- ((((STDcorr3$`4-Nitrophenol`[x]/1000))/STDcorr3$Org[x])*100)}
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc5[x] <- ((((STDcorr3$`3-Methyl-2-nitrophenol`[x]/1000))/STDcorr3$Org[x])*100)}
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc6[x] <- ((((STDcorr3$`4-Methyl-3-nitrophenol`[x]/1000))/STDcorr3$Org[x])*100)}
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc7[x] <- ((((STDcorr3$`4-Methoxy-2-nitrophenol`[x]/1000))/STDcorr3$Org[x])*100)}
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc8[x] <- ((((STDcorr3$`Conc.x`[x]/1000))/STDcorr3$Org[x])*100)}
STDcorr3$perc8<-as.numeric(STDcorr3$perc8)
STDcorr3$perc8[is.na(STDcorr3$perc8)] <- 0
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc9[x] <- ((((STDcorr3$`Conc.y`[x]/1000))/STDcorr3$Org[x])*100)}
STDcorr3$perc9<-as.numeric(STDcorr3$perc9)
STDcorr3$perc9[is.na(STDcorr3$perc9)] <- 0
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perc10[x] <- ((((STDcorr3$`Conc`[x]/1000))/STDcorr3$Org[x])*100)}
STDcorr3$perc10<-as.numeric(STDcorr3$perc10)
STDcorr3$perc10[is.na(STDcorr3$perc10)] <- 0
for(x in 1:length(STDcorr3$datetime)) {STDcorr3$perctotal[x] <- (STDcorr3$perc[x]+STDcorr3$perc2[x]+STDcorr3$perc3[x]+STDcorr3$perc4[x]+STDcorr3$perc5[x]+STDcorr3$perc6[x]+STDcorr3$perc7[x]+STDcorr3$perc8[x]+STDcorr3$perc9[x]+STDcorr3$perc10[x])}

STDcorr4 <- STDcorr3%>%  
  select(-perctotal,
         -perc,
         -perc2,
         -perc3,
         -perc4,
         -perc5,
         -perc6,
         -perc7,
         -perc8,
         -perc9,
         -perc10,
         -Conc,
         -Conc.x,
         -Conc.y,
         -Sample.ID,
         -Chl,
         -NH4,
         -`2-NapthylSulfate`,
         -`OctylSulfate`,
         -`2,6-dimethyl-4-nitrophenol`,
         -`4-Nitrophenol`,
         -`2,6-dimethyl-4-nitrophenol`,
         -`4-Methyl-3-nitrophenol`,
         -`cis-pinonicacid`,
         -`DodecylSulfate`,
         -CCOA,
         -COA)

STDcorr5 <- Samples[which(Samples$Compound == "C10H17NO3SO4 - RT6.98" ),]
STDcorr5 <- STDcorr5%>%  
  select(datetime,
         Quant)
colnames(STDcorr5)[2] <- "C10H17NO3SO4 - RT6.98"
STDcorr4 <- STDcorr4%>%
  left_join(STDcorr5, by = "datetime")
corPlot(STDcorr4,
        pollutants = names(STDcorr4)[c(2:13)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of nitroaromatics with AMS",
        order="alphabet",
        cluster=TRUE)

massloading2 <- read.csv("C:/Users/William/Documents/WACL stuff/massloading.csv", na.strings="N/F")
colnames(massloading2)[15] <- "Potassium ugm-3"
massloading2 <- massloading2%>%  
  select(datetime,
         `Potassium ugm-3`)
massloading2$datetime <- dmy_hm(massloading2$datetime, tz = "UTC")
STDcorr6 <- STDcorr4%>%
  left_join(massloading2, by = "datetime")
STDcorr6 <- STDcorr6[ which(STDcorr6$datetime != "2016-11-29 13:04:00"),]
STDcorr6 <- STDcorr6[ which(STDcorr6$datetime != "2016-12-04 12:00:00"),]
ggplot(STDcorr6,aes(x=`4-Methoxy-2-nitrophenol`,y=`Potassium ugm-3`))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Potassium  ug/m3", x="4-Methoxy-2-nitrophenol ng/m3")+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+coord_cartesian(ylim=c(0,4))

slices <- c(0.175, 1.865, 0.104, 0.627, 0.075, 0.188, 0.126, 0.018, 0.204, 0.318, 93.5) 
lbls <- c("2,4-Dinitrophenol 0.175", "2-nitro-1-naphthol 1.865", "2,6-dimethyl-4-nitrophenol 0.104", "4-nitrophenol 0.627", "3-methyl-2-nitrophenol 0.075", "4-methyl-3-nitrophenol 0.188", "4-methoxy-2-nitrophenol 0.126", "2-Methyl-5-nitrophenol 0.018", "2-Methyl-3-nitrophenol 0.204", "2-nitrophenol 0.318", "Others")
pie3D(slices,labels=lbls,explode=0.1,
      main="Proportion of nitro-aromatics")
colnames(STDcorr3)[1] <- "date"
timePlot(selectByDate(STDcorr3), pollutant=c("2-Nitro-1-naphthol","2,4-Dinitrophenol","4-Nitrophenol","2,6-dimethyl-4-nitrophenol"))

OStrends5 <- Samples3[which(Samples3$Compound == "C6H5NO6S - RT3.27" ),]
OS5<-ggplot(OStrends5,aes(x=datetime,y=Quant))+geom_line(color="red") +theme_bw()+labs(y="C6H5NO6S - RT3.27 Peak area", x="Date")

time4 <- Trendsbio2
hourly<-ggplot(time4, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-nitro-1-naphthol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
time5 <- Trendsbio3
hourly2<-ggplot(time5, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2,4-dinitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8))
time6 <- Trendsbio4
hourly3<-ggplot(time6, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("2-methyl-3-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) 
time7 <- Trendsbio16
hourly4<-ggplot(time7, aes(x=datetime, y=Conc, group=num))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("4-Methoxy-2-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) 
multiplot(hourly, hourly2, hourly3, hourly4, cols=2)


Boxsummer1<-rbind(time4, time5, time6, time7)
Boxsummer1 <- Boxsummer1[which(Boxsummer1$Conc<80 ),]
boxplot(Conc~Compound, data=Boxsummer1, main="Winter", xlab="Compound", ylab="Concentration ug/m3", ylim(0,80))


METGAS128 <- summerMETGAS[ which(summerMETGAS$datetime > "2017-05-18 12:59:00" & summerMETGAS$datetime < "2017-05-18 17:31:00"),]

SummerTime1<-ggplot(time3, aes(x=datetime, y='C10H10O7S - RT1.15'))  + geom_line(aes(colour = Compound), size=1.1)  + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_colour_manual("Compound", labels = c("4-Methoxy-2-nitrophenol"), values = c("darkred")) + theme(legend.position="right") + labs(x = "Date", y = "Concentration (ng m-3)") +theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) 

Samples2<-Samples
Samples21 <- Samples2[which(Samples2$datetime<"2016-11-17 23:59:00" ),]
Samples21$num <- 1
Samples22 <- Samples2[which(Samples2$datetime>="2016-11-17 23:59:00" ),]
Samples22$num <- 2
Samples2<-rbind(Samples21,Samples22)

OStrends <- Samples2[which(Samples2$Compound == "C7H6OHNO3 - RT6.55" ),]
OStrends2 <- Samples2[which(Samples2$Compound == "4-Methoxy-2-nitrophenol" ),]
OStrends3 <- Samples2[which(Samples2$Compound == "C10H19OH(NO3)2 - RT10.05" ),]
OStrends4 <- Samples2[which(Samples2$Compound == "C6H5NO3 - RT5.14" ),]
OStrends5 <- Samples2[which(Samples2$Compound == "C12H11NO3 - RT4.80" ),]
OStrends6 <- Samples2[which(Samples2$Compound == "C10H12OHNO3 - RT7.07" ),]
OStrends7 <- Samples2[which(Samples2$Compound == "C7H7NO3 - RT11.74" ),]
OStrends8 <- Samples2[which(Samples2$Compound == "C7H7NO4 - RT11.58" ),]
OStrends9 <- Samples2[which(Samples2$Compound == "C6H4OHNO3 - RT11.44" ),]
OStrends10 <- Samples2[which(Samples2$Compound == "C6H4OHNO3 - RT10.82" ),]
OStrends11 <- Samples2[which(Samples2$Compound == "C6H4OHNO3 - RT4.11" ),]
OStrends12 <- Samples2[which(Samples2$Compound == "C7H7NO4 - RT5.62" ),]
OStrends13 <- Samples2[which(Samples2$Compound == "2,4-Dinitrophenol" ),]
OStrends14 <- Samples2[which(Samples2$Compound == "C11H9NO3 - RT8.64" ),]
OStrends15 <- Samples2[which(Samples2$Compound == "C10H7NO3 - RT8.25" ),]
OStrends16 <- Samples2[which(Samples2$Compound == "C18H27O3NO3 - RT7.02" ),]
OStrends17 <- Samples2[which(Samples2$Compound == "C18H25O3NO3 - RT7.17" ),]
OStrends18 <- Samples2[which(Samples2$Compound == "C14H11NO3 - RT4.59" ),]
OStrends19 <- Samples2[which(Samples2$Compound == "C7H7NO4 - RT5.07" ),]
OStrends20 <- Samples2[which(Samples2$Compound == "C7H7NO4 - RT3.97" ),]
OStrends21 <- Samples2[which(Samples2$Compound == "C7H7NO3 - RT6.86" ),]
OStrends22 <- Samples2[which(Samples2$Compound == "2-Methyl-3-nitrophenol" ),]
OStrends23 <- Samples2[which(Samples2$Compound == "4-Methyl-3-nitrophenol" ),]
OStrends24 <- Samples2[which(Samples2$Compound == "2,6-dimethyl-4-nitrophenol" ),]
OStrends25 <- Samples2[which(Samples2$Compound == "3-Methyl-2-nitrophenol" ),]
OStrends26 <- Samples2[which(Samples2$Compound == "4-Nitrophenol" ),]
OStrends27 <- Samples2[which(Samples2$Compound == "C6H5COOHNO3 - RT6.43" ),]
OStrends28 <- Samples2[which(Samples2$Compound == "C7H5ON" ),]
OS<-ggplot(OStrends,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H6OHNO3 - RT6.55 Peak area", x="Date") 
OS2<-ggplot(OStrends2,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="4-Methoxy-2-nitrophenol Peak area", x="Date")
OS3<-ggplot(OStrends3,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C10H19OH(NO3)2 - RT10.05 Peak area", x="Date")
OS4<-ggplot(OStrends4,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C6H5NO3 - RT5.14 Peak area", x="Date")
OS5<-ggplot(OStrends5,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C12H11NO3 - RT4.80", x="Date")
OS6<-ggplot(OStrends6,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C10H12OHNO3 - RT7.07", x="Date")
OS7<-ggplot(OStrends7,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H7NO3 - RT11.74", x="Date")
OS8<-ggplot(OStrends8,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H7NO4 - RT11.58", x="Date")
OS9<-ggplot(OStrends9,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C6H4OHNO3 - RT11.44", x="Date")
OS10<-ggplot(OStrends10,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C6H4OHNO3 - RT10.82", x="Date")
OS11<-ggplot(OStrends11,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C6H4OHNO3 - RT4.11 Peak area", x="Date")
OS12<-ggplot(OStrends12,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H7NO4 - RT5.62 Peak area", x="Date")
OS13<-ggplot(OStrends13,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="2,4-Dinitrophenol Peak area", x="Date")
OS14<-ggplot(OStrends14,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C11H9NO3 - RT8.64 Peak area", x="Date")
OS15<-ggplot(OStrends15,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C10H7NO3 - RT8.25 Peak area", x="Date")
OS16<-ggplot(OStrends16,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C18H27O3NO3 - RT7.02 Peak area", x="Date")
OS17<-ggplot(OStrends17,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C18H25O3NO3 - RT7.17 Peak area", x="Date")
OS18<-ggplot(OStrends18,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C14H11NO3 - RT4.59 Peak area", x="Date")
OS19<-ggplot(OStrends19,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H7NO4 - RT5.07 Peak area", x="Date")
OS20<-ggplot(OStrends20,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H7NO4 - RT3.97 Peak area", x="Date")
OS21<-ggplot(OStrends21,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H7NO3 - RT6.86 Peak area", x="Date")
OS22<-ggplot(OStrends22,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="2-Methyl-3-nitrophenol Peak area", x="Date")
OS23<-ggplot(OStrends23,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="4-Methyl-3-nitrophenol Peak area", x="Date")
OS24<-ggplot(OStrends24,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="2,6-dimethyl-4-nitrophenol Peak area", x="Date")
OS25<-ggplot(OStrends25,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="3-Methyl-2-nitrophenol Peak area", x="Date")
OS26<-ggplot(OStrends26,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="4-Nitrophenol Peak area", x="Date")
OS27<-ggplot(OStrends27,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C6H5COOHNO3 - RT6.43 Peak area", x="Date")
OS28<-ggplot(OStrends28,aes(x=datetime,y=Quant, group=num))+geom_line(color="red") +theme_bw()+labs(y="C7H5ON Peak area", x="Date")
multiplot(OS3, OS2, OS, OS4, OS5, OS6, OS7, OS8, OS9, OS10, OS11, OS12, OS13, OS14, OS15, OS16, OS17, OS18, OS19, OS20, OS23, OS24, OS25, OS26, OS27, OS28, cols=5)

nitrate$Quant <- NULL
colnames(timeseries2)[1] <- "datetime"
timeseries2$datetime <- as.POSIXct(timeseries2$datetime, origin="1970-01-01")
AMSnit <- nitrate%>%
  left_join(timeseries2, by = "datetime")

AMSnit$'4-Nitrophenol' <- as.numeric(AMSnit$'4-Nitrophenol')
AMSnit$datetime <- as.POSIXct(AMSnit$datetime, origin="1970-01-01")
ggplot(STDcorr,aes(x=datetime,y=OctylSulfate))+geom_line(color="red") +theme_bw()+labs(y='Octyl Sulfate', x="Date")

n = unique(orgsulf[,1])
timeseries = seq(1478735940,1481241540,60)
timeseries = data.frame(unixtime = timeseries)
orgnit$unixtime = as.numeric(orgnit$datetime)
temp_compound = orgsulf[orgsulf$Compound == n[1],]
temp_compound = temp_compound[,c("Quant","unixtime")]
temp_compound = data.frame(temp_compound)
names(temp_compound) = c(n[1],"unixtime")
timeseries = left_join(timeseries,temp_compound,"unixtime")
for (i in 2:length(n)){
  temp_compound = orgsulf[orgsulf$Compound == n[i],]
  temp_compound = temp_compound[,c("Quant","unixtime")]
  temp_compound = data.frame(temp_compound)
  names(temp_compound) = c(n[i],"unixtime")
  
  timeseries = left_join(timeseries,temp_compound,"unixtime")
}
timeseries[is.na(timeseries)] <- 0
colSums(timeseries != 0)
timeseries2 <- timeseries[, !colSums(timeseries != 0) < 20]
timeseries[timeseries < 0.1] <- NA
timeseries2[timeseries2 < 0.1] <- NA
#possibly kill computer with 191 corplot of CHOS
filepath = paste(getwd(),"/CHOS_corplot.pdf",sep = "")
pdf(width=60, height=60, file=filepath)
corPlot(timeseries,
        pollutants = names(timeseries)[c(2:188)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors with AMS, CHOS",
        order="alphabet",
        cluster=TRUE)
dev.off()
#possibly kill computer with 191 corplot of CHOS
filepath = paste(getwd(),"/CHOS_corplot2.pdf",sep = "")
pdf(width=60, height=60, file=filepath)
corPlot(timeseries2,
        pollutants = names(timeseries2)[c(2:129)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors with AMS, CHOS",
        order="alphabet",
        cluster=TRUE)
dev.off()

outCHON <-corPlot(STDcorr4,
                  pollutants = names(STDcorr4)[c(2:13)],
                  cols = "jet",
                  dendrogram = TRUE,
                  main = "Correlations of nitroaromatics with AMS",
                  order="alphabet",
                  cluster=TRUE)
mydendCHON <- as.dendrogram(outCHON$clust)
dend <- plot(as.phylo(outCHON$clust), cex=0.7)

Standardpoint <- data.frame('Sample.ID'='Standard','Compound'='Standard', 'Area'=10000000, 'RT'=0, 'datetime'='2018-01-04 09:00:00', 'Volume'=0,'Expected.RT'=0,'Carbon'=20,'Hydrogen'=0,'Oxygen'=0,'Nitrogen'=0,'Sulphur'=0,'Quant'=10000000,'Group'='Standard','DBE'=0,'AI'=0,'HCratio'=0,'OCratio'=3,'OSc'=-2,'num'=0)

Standardpoint$'Compound' <- as.character((Standardpoint$'Compound'))
Standardpoint$'Group' <- as.character((Standardpoint$'Group'))
Standardpoint$'Sample.ID' <- as.character((Standardpoint$'Sample.ID'))
Standardpoint$Area <- as.numeric((Standardpoint$Area))
Standardpoint$Quant <- as.numeric((Standardpoint$Quant))
Samples3 <- rbind(Samples2, Standardpoint)

ggplot(Samples3, aes(OCratio, HCratio, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                          panel.grid.major = element_blank(),
                                                                                                          panel.grid.minor = element_blank(),
                                                                                                          panel.background = element_blank()) +ggtitle("Average winter") + xlim(0, 3) + ylim(0, 3) + guides(size=FALSE) + ylab("H:C ratio") + xlab("O:C ratio") + scale_color_manual(breaks = c("CHO", "CHON", "CHONS", "CHOS", "Standard"), values = c("green3", "blue", "purple", "red", "grey22"))

ggplot(Samples3, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                     panel.grid.major = element_blank(),
                                                                                                     panel.grid.minor = element_blank(),
                                                                                                     panel.background = element_blank()) +ggtitle("Average winter") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + ylab("DBE") + xlab("Carbon number") + scale_color_manual(breaks = c("CHO", "CHON", "CHONS", "CHOS", "Standard"), values = c("green3", "blue", "purple", "red", "grey22"))




HetGroup2$datetime <- as.POSIXct(HetGroup2$datetime, origin="1970-01-01")
HetGroup3 <- HetGroup2[ which(HetGroup2$datetime < "2016-12-05 00:00:00"),]
HetGroupCHO$datetime <- as.POSIXct(HetGroupCHO$datetime, origin="1970-01-01")
HetGroupCHO <- HetGroupCHO[ which(HetGroupCHO$datetime < "2016-12-05 00:00:00"),]
HetGroupCHON$datetime <- as.POSIXct(HetGroupCHON$datetime, origin="1970-01-01")
HetGroupCHON <- HetGroupCHON[ which(HetGroupCHON$datetime < "2016-12-05 00:00:00"),]
HetGroupCHOS$datetime <- as.POSIXct(HetGroupCHOS$datetime, origin="1970-01-01")
HetGroupCHOS <- HetGroupCHOS[ which(HetGroupCHOS$datetime < "2016-12-05 00:00:00"),]
HetGroupCHONS$datetime <- as.POSIXct(HetGroupCHONS$datetime, origin="1970-01-01")
HetGroupCHONS <- HetGroupCHONS[ which(HetGroupCHONS$datetime < "2016-12-05 00:00:00"),]
SummedOC$datetime <- as.POSIXct(SummedOC$datetime, origin="1970-01-01")
SummedOC <- SummedOC[ which(SummedOC$datetime < "2016-12-05 00:00:00"),]

WinCHO<-ggplot(HetGroupCHO, aes(x=datetime, y=Quant, group=num))+ylim(0,400000) + geom_line(color='green4', size=1.1)  + theme_bw() + theme(text = element_text(size=20),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHO Peak area (a.u. m-3)") +theme(axis.text=element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                                        axis.title=element_text(size=20))
WinCHON<-ggplot(HetGroupCHON, aes(x=datetime, y=Quant, group=num))+ylim(0,4500000) + geom_line(color='navyblue', size=1.1)  + theme_bw() + theme(text = element_text(size=20),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHON Peak area (a.u. m-3)") +theme(axis.text=element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                                             axis.title=element_text(size=20))
WinCHOS<-ggplot(HetGroupCHOS, aes(x=datetime, y=Quant, group=num))+ylim(0,650000) + geom_line(color='red2', size=1.1)  + theme_bw() + theme(text = element_text(size=20),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHOS Peak area (a.u. m-3)")+theme(axis.text=element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                                        axis.title=element_text(size=20))
WinCHONS<-ggplot(HetGroupCHONS, aes(x=datetime, y=Quant, group=num))+ylim(0,200000) + geom_line(color='magenta4', size=1.1)  + theme_bw() + theme(text = element_text(size=20),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "CHONS Peak area (a.u. m-3)")+theme(axis.text=element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                                               axis.title=element_text(size=20))
WinOC<-ggplot(SummedOC, aes(x=datetime, y=Quant, group=num))+ylim(0,100000000) + geom_line(size=1.1)  + theme_bw() + theme(text = element_text(size=20),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +  labs(x = "Date", y = "Total OC Peak area (a.u. m-3)") +theme(axis.text=element_text(size=20),
                                                                                                                                                                                                                                                                                                                                                         axis.title=element_text(size=20))
multiplot(WinCHO, WinCHON, WinCHOS, WinCHONS, WinOC, cols=1)

AMSComp<-AMS
AMSComp <- AMSComp%>%
  right_join(SummedOC, by = "datetime")

AMSComp <- AMSComp[ which(AMSComp$datetime > "2016-11-15 00:00:00"),]

lm_eqn <- function(AMSComp){
  m <- lm(Org ~ Quant, AMSComp);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
WinAMSScat<-ggplot(AMSComp,aes(x=Quant,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon (AMS) ug/m3", x="Total organic carbon peak area (Library)")+ xlim(0, 120000000) + ylim(0, 120) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_text(x = 25000000, y = 115, label = lm_eqn(AMSComp), parse = TRUE) +theme(axis.text=element_text(size=20), axis.title=element_text(size=20))
HetGroupCHO$Conc<-0
for(x in 1:length(HetGroupCHO$datetime)) {HetGroupCHO$Conc[x] <- ((HetGroupCHO$Quant[x]-97500)/85439)}
HetGroupCHO[HetGroupCHO < 0] <- 0
HetGroupCHON$Conc<-0
for(x in 1:length(HetGroupCHON$datetime)) {HetGroupCHON$Conc[x] <- ((HetGroupCHON$Quant[x]-308980)/104986)}
HetGroupCHON[HetGroupCHON < 0] <- 0
HetGroupCHOS$Conc<-0
for(x in 1:length(HetGroupCHOS$datetime)) {HetGroupCHOS$Conc[x] <- ((HetGroupCHOS$Quant[x]-324579)/147725)}
HetGroupCHOS[HetGroupCHOS < 0] <- 0
SummedOCNormal <- rbind(HetGroupCHO, HetGroupCHON, HetGroupCHOS)
SummedOCNormal2 <- aggregate(Conc ~ datetime, SummedOCNormal, sum)
AMSComp <- AMSComp%>%
  right_join(SummedOCNormal2, by = "datetime")

AMSComp <- AMSComp[ which(AMSComp$datetime > "2016-11-15 00:00:00"),]

lm_eqn <- function(AMSComp){
  m <- lm(Org ~ Conc, AMSComp);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}
WinAMSNormScat<-ggplot(AMSComp,aes(x=Conc,y=Org))+geom_point(color="blue", size=1.1)+ geom_smooth(method="lm", se=FALSE, colour="red")+labs(y="Total organic carbon (AMS) ug/m3", x="Ionisation corrected total organic carbon (Library) ug/m3")+ xlim(0, 30) + ylim(0, 120) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_text(x = 5, y = 115, label = lm_eqn(AMSComp), parse = TRUE) +theme(axis.text=element_text(size=20), axis.title=element_text(size=20))
multiplot(WinAMSScat, WinAMSNormScat, cols=1)
multiplot(WinAMSScat, WinAMSNormScat, SumAMSScat, SumAMSNormScat, cols=2)


WinTimeAMS <- ggplot(AMSComp, aes(x = datetime))
WinTimeAMS <- WinTimeAMS + geom_line(aes(y = Org, colour = "TOC (AMS) ug/m3"))
WinTimeAMS <- WinTimeAMS + geom_line(aes(y = Conc/0.2, colour = "Ionisation corrected TOC ug/m3"))
WinTimeAMS <- WinTimeAMS + scale_y_continuous(sec.axis = sec_axis(~.*0.2, name = "Ionisation corrected TOC ug/m3"))
WinTimeAMS <- WinTimeAMS + scale_colour_manual(values = c("blue", "red"))
WinTimeAMS <- WinTimeAMS + labs(y = "TOC ug/m3",
                                x = "Date",
                                colour = "Winter")
WinTimeAMS <- WinTimeAMS + theme(legend.position = c(0.4, 0.9))
WinTimeAMS

WinTimeAMS2 <- ggplot(AMSComp, aes(x = datetime))
WinTimeAMS2 <- WinTimeAMS2 + geom_line(aes(y = Org, colour = "Total organic carbon (AMS) ug/m3"))
WinTimeAMS2 <- WinTimeAMS2 + geom_line(aes(y = Quant/2500000, colour = "Total organic carbon (library) a.u./m3"))
WinTimeAMS2 <- WinTimeAMS2 + scale_y_continuous(sec.axis = sec_axis(~.*2500000, name = "Peak area a.u./m3"))
WinTimeAMS2 <- WinTimeAMS2 + scale_colour_manual(values = c("blue", "red"))
WinTimeAMS2 <- WinTimeAMS2 + labs(y = "TOC ug/m3",
                                  x = "Date",
                                  colour = "Winter")
WinTimeAMS2 <- WinTimeAMS2 + theme(legend.position = c(0.4, 0.9))
WinTimeAMS2

multiplot(WinTimeAMS2, WinTimeAMS, cols=1)
multiplot(WinTimeAMS2, WinTimeAMS, SumTimeAMS2, SumTimeAMS, cols=2)

mean(HetGroupCHON$Quant, na.rm=TRUE)

timeseries$unixtime<-as.POSIXct(timeseries$unixtime, origin="1970-01-01 00:00:00")

timeTemp <- timeseries[!is.na(timeseries[i]),]
timeTemp[is.na(timeTemp)] <- 0
filepath = paste(getwd(),"test.pdf",sep = "")
pdf(width=10, height=10, file=filepath)
for (i in 2:ncol(timeTemp)){
  
  plot(timeTemp$unixtime,timeTemp[,i],type="l",xlab = "Date")
  
}
dev.off() 

ncol(timeTemp)

CIMS <- read.csv("C:/Users/William/Documents/WACL stuff/nitrophenolaerosoldata.csv", na.strings=0)
timeseries = seq(1478735940,1481241540,60)
timeseries = data.frame(unixtime = timeseries)
timeseries$unixtime<-as.POSIXct(timeseries$unixtime, origin="1970-01-01 00:00:00")
colnames(CIMS)[1] <- "datetime"
colnames(timeseries)[1] <- "datetime"
CIMS$datetime <- dmy_hm(CIMS$datetime, tz = "UTC")
timeseries2 <- CIMS%>%
  right_join(timeseries, by = "datetime")

WinTimeCIMS <- ggplot(AMSComp, aes(x = datetime))
WinTimeCIMS <- WinTimeCIMS + geom_line(aes(y = Org, colour = "Total organic carbon (AMS) ug/m3"))
WinTimeCIMS <- WinTimeCIMS + geom_line(aes(y = Quant/2500000, colour = "Total organic carbon (library) a.u./m3"))
WinTimeCIMS <- WinTimeCIMS + scale_y_continuous(sec.axis = sec_axis(~.*2500000, name = "Peak area a.u./m3"))
WinTimeCIMS <- WinTimeCIMS + scale_colour_manual(values = c("blue", "red"))
WinTimeCIMS <- WinTimeCIMS + labs(y = "TOC ug/m3",
                                  x = "Date",
                                  colour = "Winter")
WinTimeCIMS <- WinTimeCIMS + theme(legend.position = c(0.4, 0.9))
WinTimeCIMS

mean(Trendsbio7$Conc, na.rm=TRUE)

ggplot(Samples3, aes(Carbon, DBE, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                     panel.grid.major = element_blank(),
                                                                                                     panel.grid.minor = element_blank(),
                                                                                                     panel.background = element_blank()) +ggtitle("Average winter") + xlim(0, 20) + ylim(0, 8) + guides(size=FALSE) + ylab("DBE") + xlab("Carbon number")



ggplot(Samples3, aes(Carbon, OSc, colour = Group, size = Quant)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                     panel.grid.major = element_blank(),
                                                                                                     panel.grid.minor = element_blank(),
                                                                                                     panel.background = element_blank()) +ggtitle("Average winter") + xlim(0, 20) + ylim(-2, 5) + guides(size=FALSE) + ylab("OSc") + xlab("Carbon number")



STDcali24np <- Standards[ which(Standards$Compound=='2,4-Dinitrophenol'),]
STDcali24np <- STDcali24np[ which(STDcali24np$Area > 10000),]
ggplot(STDcali24np, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("2,4-Nitrophenol Concentration ppb") + geom_smooth(method=lm)
calib<-lm(Area ~ Conc, data = STDcali24np)
summary.lm(calib)

STDcalinaph <- Standards[ which(Standards$Compound=='2-Nitro-1-Naphthol'),]
STDcalinaph <- STDcalinaph[ which(STDcalinaph$Area > 10000),]
ggplot(STDcalinaph, aes(Conc, Area)) + geom_point() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank()) + guides(size=FALSE) + ylab("Peak area a.u.") + xlab("2-Nitro-1-Naphthol Concentration ppb") + geom_smooth(method=lm)
calib2<-lm(Area ~ Conc, data = STDcalinaph)
summary.lm(calib2)

avgmetwinter <- read.csv("C:/Users/William/Documents/WACL stuff/aphh_winter_filter_aggregate_merge.csv", na.strings=0)
wintertimeseries <- read.csv("C:/Users/William/Documents/WACL stuff/wintertimeres.csv", na="NA", check.names = F)
colnames(avgmetwinter)[4] <- "Sample.ID"
colnames(avgmetwinter)[2] <- "datetime"
avgmetwinter$Sample.ID<-as.character(avgmetwinter$Sample.ID)
avgmetwinter_8<-avgmetwinter
colnames(avgmetwinter)[93] <- "ws"
colnames(avgmetwinter)[92] <- "wd"
colnames(avgmetwinter_8)[91] <- "ws"
colnames(avgmetwinter_8)[90] <- "wd"
avgmetwinter$wd<-as.numeric(avgmetwinter$wd)
avgmetwinter$ws<-as.numeric(avgmetwinter$ws)
avgmetwinter_8$wd<-as.numeric(avgmetwinter_8$wd)
avgmetwinter_8$ws<-as.numeric(avgmetwinter_8$ws)

avgmetwinter2 <- avgmetwinter%>%
  select(wd,ws,datetime) 

avgmetwinter2$wd<-as.numeric(avgmetwinter2$wd)
avgmetwinter2$ws<-as.numeric(avgmetwinter2$ws)

avgmetwinter2$datetime  <- ymd_hms(avgmetwinter2$datetime, tz = "UTC")
wintertimeseries$datetime  <- dmy_hm(wintertimeseries$datetime, tz = "UTC")
avgmetwinter3<-avgmetwinter2 %>% 
  right_join(wintertimeseries, by = "datetime")

colnames(avgmetwinter3)[3] <- "date"
polar_group2<-polarPlot(avgmetwinter3, pollutant = c("2-Nitro-1-Naphthol"), k = 50, key.header = "Peak area (a.u. / m3)")
polar_group3<-polarPlot(avgmetwinter3, pollutant = c("2-Methyl-5-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")
polar_group1<-polarPlot(avgmetwinter3, pollutant = c("4-Methyl-3-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")
polar_group4<-polarPlot(avgmetwinter3, pollutant = c("4-Methoxy-2-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")

polarCluster(avgmetwinter3, pollutant = c("2-Nitro-1-Naphthol"), k = 50, statistic = "cpf", percentile = 90, n.clusters = 2)


#all r corr plot, this may take a while
filepath = paste(getwd(),"/allRcorr.pdf",sep = "")
pdf(width=80, height=80, file=filepath)
corPlot(avgmetwinter3,
        pollutants = names(avgmetwinter3)[c(4:621)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

nits = unique(Nitcomps[,1])

CHONs<- avgmetwinter3[,grep(c("itro"), names(avgmetwinter3))]

avgmetwinter4 <- avgmetwinter3%>%
  select(wd,ws,date,'2,6-dimethyl-4-nitrophenol')

CHONs<-CHONs %>% 
  right_join(avgmetwinter4, by = "2,6-dimethyl-4-nitrophenol")

#remove 3np due to low data content
CHONs<-CHONs[,-c(2,6)] 

filepath = paste(getwd(),"/nitrophenRcorr.pdf",sep = "")
pdf(width=10, height=10, file=filepath)
corPlot(CHONs,
        pollutants = names(CHONs)[c(1:10)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified nitro-aromatics",
        order="alphabet",
        cluster=TRUE)
dev.off()

CHONs_r_plot <- corPlot(CHONs,
                      pollutants = names(CHONs)[c(1:10)],
                      cols = "jet",
                      dendrogram = TRUE,
                      main = "Correlations of identified nitro-aromatics",
                      order="alphabet",
                      cluster=TRUE)
CHONs_r_dend <- as.dendrogram(CHONs_r_plot$clust)

CHOS_samples <- Samples[which(Samples$Group=="CHOS" ),]
CHON_samples <- Samples[which(Samples$Group=="CHON" ),]
CHONS_samples <- Samples[which(Samples$Group=="CHONS" ),]
CHO_samples <- Samples[which(Samples$Group=="CHO" ),]

CHOS_compounds = unique(CHOS_samples[,2])
CHONS_compounds = unique(CHONS_samples[,2])
CHON_compounds = unique(CHON_samples[,2])
CHO_compounds = unique(CHO_samples[,2])

timeseries_CHON2 = seq(1478735940,1481241540,60)
timeseries_CHON3 = seq(1478735940,1481241540,60)
timeseries_CHON2 = data.frame(unixtime = timeseries_CHON2)
timeseries_CHON3 = data.frame(unixtime = timeseries_CHON3)
CHON_samples$unixtime = as.numeric(CHON_samples$datetime)
CHON_samples$Quant = as.numeric(CHON_samples$Quant)

for (i in 2:length(CHON_compounds)){
  temp_compound2 = CHON_samples[CHON_samples$Compound == CHON_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHON_compounds[i],"unixtime")
  
  timeseries_CHON2 = right_join(temp_compound2,timeseries_CHON2,"unixtime")
}
timeseries_CHON3 <- time %>% 
  left_join(timeseries_CHON2, by = "unixtime")

timeseries_CHONS2 = seq(1478735940,1481241540,60)
timeseries_CHONS3 = seq(1478735940,1481241540,60)
timeseries_CHONS2 = data.frame(unixtime = timeseries_CHONS2)
timeseries_CHONS3 = data.frame(unixtime = timeseries_CHONS3)
CHONS_samples$unixtime = as.numeric(CHONS_samples$datetime)
CHONS_samples$Quant = as.numeric(CHONS_samples$Quant)

for (i in 2:length(CHONS_compounds)){
  temp_compound2 = CHONS_samples[CHONS_samples$Compound == CHONS_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHONS_compounds[i],"unixtime")
  
  timeseries_CHONS2 = right_join(temp_compound2,timeseries_CHONS2,"unixtime")
}
timeseries_CHONS3 <- time %>% 
  left_join(timeseries_CHONS2, by = "unixtime")

timeseries_CHO2 = seq(1478735940,1481241540,60)
timeseries_CHO3 = seq(1478735940,1481241540,60)
timeseries_CHO2 = data.frame(unixtime = timeseries_CHO2)
timeseries_CHO3 = data.frame(unixtime = timeseries_CHO3)
CHO_samples$unixtime = as.numeric(CHO_samples$datetime)
CHO_samples$Quant = as.numeric(CHO_samples$Quant)

for (i in 2:length(CHO_compounds)){
  temp_compound2 = CHO_samples[CHO_samples$Compound == CHO_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHO_compounds[i],"unixtime")
  
  timeseries_CHO2 = right_join(temp_compound2,timeseries_CHO2,"unixtime")
}
timeseries_CHO3 <- time %>% 
  left_join(timeseries_CHO2, by = "unixtime")

timeseries_CHOS2 = seq(1478735940,1481241540,60)
timeseries_CHOS3 = seq(1478735940,1481241540,60)
timeseries_CHOS2 = data.frame(unixtime = timeseries_CHOS2)
timeseries_CHOS3 = data.frame(unixtime = timeseries_CHOS3)
CHOS_samples$unixtime = as.numeric(CHOS_samples$datetime)
CHOS_samples$Quant = as.numeric(CHOS_samples$Quant)

for (i in 2:length(CHOS_compounds)){
  temp_compound2 = CHOS_samples[CHOS_samples$Compound == CHOS_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHOS_compounds[i],"unixtime")
  
  timeseries_CHOS2 = right_join(temp_compound2,timeseries_CHOS2,"unixtime")
}
timeseries_CHOS3 <- time %>% 
  left_join(timeseries_CHOS2, by = "unixtime")

filepath = paste(getwd(),"/CHOSRcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHOS3,
        pollutants = names(timeseries_CHOS3)[c(2:188)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHOS",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

timeseries_CHON4<-timeseries_CHON3 %>% 
  select(-contains('C10H19OH(NO3)2 - RT10'),
         -contains('C19H29O3NO3 - RT7.4'),
         -contains('C18H25O3NO3 - RT8.5'),
         -contains('18H25O3NO3 - RT6.2'),
         -contains('C11H19NO9 - RT5.3'),
         -contains('2-Methyl-3-nitropheno'))

filepath = paste(getwd(),"/CHONRcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHON4,
        pollutants = names(timeseries_CHON4)[c(2:40)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHON",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

filepath = paste(getwd(),"/CHORcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHO3,
        pollutants = names(timeseries_CHO3)[c(2:231)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHO",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

filepath = paste(getwd(),"/CHONSRcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHONS3,
        pollutants = names(timeseries_CHONS3)[c(2:105)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHONS",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

CHONSs_r_plot <-corPlot(timeseries_CHONS3,
                       pollutants = names(timeseries_CHONS3)[c(2:45)],
                       cols = "jet",
                       dendrogram = TRUE,
                       main = "Correlations of identified CHON",
                       order="alphabet",
                       cluster=TRUE,
                       text.col = NA)
dev.off()
CHONSs_r_dend <- as.dendrogram(CHONSs_r_plot$clust)

timeseries_CHOS4<-timeseries_CHOS3
timeseries_CHOS4[is.na(timeseries_CHOS4)] <- 1
log_CHOS <- log(timeseries_CHOS4[, 2:188])

timeseries_CHON4<-timeseries_CHON3
timeseries_CHON4[is.na(timeseries_CHON4)] <- 1
log_CHON <- log(timeseries_CHON4[, 2:45])

timeseries_CHONS4<-timeseries_CHONS3
timeseries_CHONS4[is.na(timeseries_CHONS4)] <- 1
log_CHONS <- log(timeseries_CHONS4[, 2:105])

timeseries_CHO4<-timeseries_CHO3
timeseries_CHO4[is.na(timeseries_CHO4)] <- 1
log_CHO <- log(timeseries_CHO4[, 2:231])

CHO_pca <- prcomp(log_CHO,
                  center = TRUE,
                  scale. = TRUE) 
plot(CHO_pca, type = "l")

CHON_pca <- prcomp(log_CHON,
                  center = TRUE,
                  scale. = TRUE) 
plot(CHON_pca, type = "l")

CHOS_pca <- prcomp(log_CHOS,
                  center = TRUE,
                  scale. = TRUE) 
plot(CHOS_pca, type = "l")

CHONS_pca <- prcomp(log_CHONS,
                  center = TRUE,
                  scale. = TRUE) 
plot(CHONS_pca, type = "l")

fviz_pca_var(CHON_pca, col.var="contrib")+
  scale_color_gradient2(low="blue", mid="green",
                        high="red", midpoint=3) +
  theme_minimal()

# log transform 
avgmetwinter4<-avgmetwinter3
avgmetwinter4[is.na(avgmetwinter4)] <- 1
log_all <- log(avgmetwinter4[, 4:621])
#ir.species <- iris[, 5]

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
all_pca <- prcomp(log_all,
                 center = TRUE,
                 scale. = TRUE) 
princomp(log_all, cor = FALSE, scores = TRUE)

colnames(avgmetwinter2)[3] <- "unixtime"

CHON_PCA_group1_source<-polarPlot(avgmetwinter3, pollutant = c("2,6-dimethyl-4-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")
CHON_PCA_group2_source<-polarPlot(avgmetwinter3, pollutant = c("C7H7NO4 - RT5.07"), k = 50, key.header = "Peak area (a.u. / m3)")
CHON_PCA_group3_source<-polarPlot(avgmetwinter3, pollutant = c("2-Nitro-1-Naphthol"), k = 50, key.header = "Peak area (a.u. / m3)")

CHOS_PCA_group1_source<-polarPlot(avgmetwinter3, pollutant = c("C9H12O5S - RT3.86"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group2_source<-polarPlot(avgmetwinter3, pollutant = c("C5H8SO5 - RT0.94"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group3_source<-polarPlot(avgmetwinter3, pollutant = c("C16H28O4SO4 - RT6.26"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group5_source<-polarPlot(avgmetwinter3, pollutant = c("C5H6O4SO4 - RT0.71"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group5_source<-polarCluster(avgmetwinter3, pollutant = c("C5H6O4SO4 - RT0.71"), k = 50, key.header = "Peak area (a.u. / m3)", n.clusters = 3)
CHOS_PCA_group6_source<-polarPlot(avgmetwinter3, pollutant = c("C10H20OSO4 - RT6.71"), k = 50, key.header = "Peak area (a.u. / m3)")

CHO_PCA_group1_source<-polarPlot(avgmetwinter3, pollutant = c("C9H8O4 - RT5.25"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group2_source<-polarPlot(avgmetwinter3, pollutant = c("C7H6O2 - RT3.64"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group3_source<-polarPlot(avgmetwinter3, pollutant = c("C8H6O4 - RT4.06"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group4_source<-polarPlot(avgmetwinter3, pollutant = c("C6H10O4 - RT1.18 (Pinene SOA)"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group5_source<-polarPlot(avgmetwinter3, pollutant = c("nor-caryophyllenic acid (B-Caryophyllene SOA C8H12O4)"), k = 50, key.header = "Peak area (a.u. / m3)")


wintertimeseries_transposed <- read.csv("C:/Users/William/Documents/WACL stuff/wintertimeres_transposed.csv", na="NA", check.names = F)

wintertimeseries_transposed_1<-wintertimeseries_transposed
wintertimeseries_transposed_1[is.na(wintertimeseries_transposed_1)] <- 1
log_CHOS_t <- log(wintertimeseries_transposed_1[, 2:124])
CHOS_t_pca <- prcomp(log_CHOS_t,
                   center = TRUE,
                   scale. = TRUE) 

CHON_dend <- corPlot(timeseries_CHON3,
        pollutants = names(timeseries_CHON3)[c(2:45)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHON",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

CHON_dend_table<-Map(as.data.frame(CHON_dend))

mydendCHON <- as.dendrogram(CHON_dend$clust)
plot(mydendCHON)
k <- 6
mydendCHON <- color_branches(mydendCHON, k = k)
plot(mydendCHON)
labels_dend <- labels(mydendCHON)
groups <- cutree(mydendCHON, k=6, order_clusters_as_data = FALSE)
dends <- list()
for(i in 1:k) {
  labels_to_keep <- labels_dend[i != groups]
  dends[[i]] <- prune(dend, labels_to_keep)}

CHON_dend_table <-as.data.frame(CHON_dend$data)
CHON_dend_table <- CHON_dend_table[ which(CHON_dend_table$cor > -1),]
CHON_dend_table <- CHON_dend_table[ which(CHON_dend_table$cor < 1),]

n = unique(CHON_dend_table[,1])
CHON_dend_table2 <- CHON_dend_table %>% 
  select(x)

for (i in 2:length(n)){
  temp_compound2 = CHON_dend_table[CHON_dend_table$y == n[i],]
  temp_compound2 = temp_compound2[,c("cor","x")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(n[i],"x")
  
  timeseries3 = right_join(temp_compound2,CHON_dend_table2,"x")
}

timeseries3 <- time %>% 
  left_join(timeseries2, by = "unixtime")


m_dist<-dist(log_CHOS_t,diag = FALSE )
m_hclust<-hclust(m_dist, method= "complete")
plot(m_hclust)

# I'll do this to just 4 clusters for illustrative purposes
k <- 6
cols <- rainbow_hcl(k)
dend <- as.dendrogram(m_hclust)
dend <- color_branches(dend, k = k)
plot(dend)
labels_dend <- labels(dend)
groups <- cutree(dend, k=6, order_clusters_as_data = FALSE)
dends <- list()
for(i in 1:k) {
  labels_to_keep <- labels_dend[i != groups]
  dends[[i]] <- prune(dend, labels_to_keep)
}

length(which(timeseries_CHON3$`C8H8OHNO3 - RT8.52` != "NA")) 

timeseries_CHON5<-timeseries_CHON3
timeseries_CHON5[is.na(timeseries_CHON5)] <- 0
colSums(timeseries_CHON5 != 0)
timeseries_CHON5 <- timeseries_CHON5[, !colSums(timeseries_CHON5 != 0) < 5]
timeseries_CHON5[timeseries_CHON5 < 0.1] <- NA

timeseries_CHON6<-timeseries_CHON5 %>% 
  select(-contains('C10H19OH(NO3)2 - RT10'),
         -contains('C19H29O3NO3 - RT7.4'),
         -contains('C18H25O3NO3 - RT8.5'),
         -contains('18H25O3NO3 - RT6.2'),
         -contains('C11H19NO9 - RT5.3'),
         -contains('2-Methyl-3-nitropheno'))

filepath = paste(getwd(),"/CHONRcorr2.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHON6,
        pollutants = names(timeseries_CHON6)[c(2:27)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHON",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

dend <- plot(as.phylo(outCHON$clust), cex=0.7)

wintertimeseries_transposed <- read.csv("C:/Users/William/Documents/WACL stuff/Wintertimeres.csv", na="NA", check.names = F)

wintertimeseries_transposed2 <- wintertimeseries_transposed %>% 
  select(datetime,
         '4-Methoxy-2-nitrophenol',
         '4-Methyl-2-nitrophenol',
         '2-Methyl-3-nitrophenol',
         '4-Nitrocatechol',
         '2-Nitrophenol')

wintertimeseries_transposed2$datetime <- dmy_hm(wintertimeseries_transposed2$datetime, tz = "UTC")

AMS_STD_Comp <- AMSnit %>% 
  select(datetime,
         CCOA,
         COA,
         BBOA,
         OPOA,
         LOOOA,
         MOOOA,
         Chl,
         NH4,
         NO3,
         SO4,
         Org,
         '2,4-Dinitrophenol',
         '2,6-dimethyl-4-nitrophenol',
         '3-Methyl-2-nitrophenol',
         '4-Methyl-3-nitrophenol',
         '4-Nitrophenol',
         'C10H7NO3 - RT8.25')

wintertimeseries_transposed3 <- wintertimeseries_transposed2 %>% 
  right_join(AMS_STD_Comp, by = "datetime")

colnames(wintertimeseries_transposed3)[23] <- "2-Nitro-1-naphthol"


filepath = paste(getwd(),"/CHONRcorrAMS.pdf",sep = "")
pdf(width=20, height=20, file=filepath)
corPlot(wintertimeseries_transposed3,
        pollutants = names(wintertimeseries_transposed3)[c(2:23)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHON with AMS",
        order="alphabet",
        cluster=TRUE)
dev.off()
  
wintertimeseries_transposed <- read.csv("C:/Users/William/Documents/WACL stuff/aphh_winter_filter_aggregate_merge.csv", na="NA", check.names = F)

timeseries_CHON5_group1 <- timeseries_CHON5 %>% 
  select(unixtime,
         '2-Nitrophenol',
         'C10H12OHNO3 - RT7.07')

timeseries_CHON5_group4 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C7H7NO4 - RT11.58',
         'C8H8OHNO3 - RT8.52',
         'C8H8OHNO3 - RT7.87',
         '4-Methoxy-2-nitrophenol')

timeseries_CHON5_group2 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C6H4OHNO3 - RT11.44',
         'C7H6OHNO3 - RT6.55',
         'C16H10O4N2 - RT6.98',
         'C18H25O3NO3 - RT6.35')

timeseries_CHON5_group3 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C6H4OHNO3 - RT10.82',
         'C8H8OHNO3 - RT5.00',
         'C7H7NO4 - RT5.07',
         'C8H8OHNO3 - RT6.50')

timeseries_CHON5_group6 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C10H7NO3 - RT8.25',
         'C11H9NO3 - RT8.64',
         'C9H5O4N - RT3.09',
         'C9H5O4N - RT3.95')

timeseries_CHON5_group5 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C7H7NO3 - RT6.86',
         'C8H5NO3 - RT3.65',
         'C7H7NO4 - RT3.97',
         '4-Nitrophenol',
         'C6H5NO3 - RT5.14',
         'C7H7NO4 - RT5.62',
         '4-Methyl-3-nitrophenol',
         '2,6-dimethyl-4-nitrophenol')

#timeseries_CHON5_group1[is.na(timeseries_CHON5_group1)] <- 0
timeseries_CHON5_group1$sum <- rowSums( timeseries_CHON5_group1[2:3] , na.rm=TRUE)
timeseries_CHON5_group1$unixtime<-as.POSIXct(timeseries_CHON5_group1$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group1, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group1") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
timeseries_CHON5_group1 <- timeseries_CHON5_group1 %>% 
  left_join(avgmetwinter2, by = "unixtime")

CHON_PCA_group1_source<-polarPlot(timeseries_CHON5_group1, pollutant = c("sum"), k = 10, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group2[is.na(timeseries_CHON5_group2)] <- 0
timeseries_CHON5_group2$sum <- rowSums( timeseries_CHON5_group2[2:5] , na.rm=TRUE)
timeseries_CHON5_group2$unixtime<-as.POSIXct(timeseries_CHON5_group2$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group2, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group2") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
timeseries_CHON5_group2 <- timeseries_CHON5_group2 %>% 
  left_join(avgmetwinter2, by = "unixtime")

CHON_PCA_group2_source<-polarPlot(timeseries_CHON5_group2, pollutant = c("sum"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group3[is.na(timeseries_CHON5_group3)] <- 0
timeseries_CHON5_group3$sum <- rowSums( timeseries_CHON5_group3[2:4], na.rm=TRUE)
timeseries_CHON5_group3$unixtime<-as.POSIXct(timeseries_CHON5_group3$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group3, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group3") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
timeseries_CHON5_group3 <- timeseries_CHON5_group3 %>% 
  left_join(avgmetwinter2, by = "unixtime")

CHON_PCA_group3_source<-polarPlot(timeseries_CHON5_group3, pollutant = c("sum"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group4[is.na(timeseries_CHON5_group4)] <- 0
timeseries_CHON5_group4$sum <- rowSums( timeseries_CHON5_group4[2:5], na.rm=TRUE )
timeseries_CHON5_group4$unixtime<-as.POSIXct(timeseries_CHON5_group4$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group4, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group4") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
timeseries_CHON5_group4 <- timeseries_CHON5_group4 %>% 
  left_join(avgmetwinter2, by = "unixtime")

CHON_PCA_group4_source<-polarPlot(timeseries_CHON5_group4, pollutant = c("sum"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group5[is.na(timeseries_CHON5_group5)] <- 0
timeseries_CHON5_group5$sum <- rowSums( timeseries_CHON5_group5[2:9] , na.rm=TRUE)
timeseries_CHON5_group5$unixtime<-as.POSIXct(timeseries_CHON5_group5$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group5, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                     panel.grid.major = element_blank(),
                                                                                                     panel.grid.minor = element_blank(),
                                                                                                     panel.background = element_blank()) +ggtitle("Average group5") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
timeseries_CHON5_group5 <- timeseries_CHON5_group5 %>% 
  left_join(avgmetwinter2, by = "unixtime")

CHON_PCA_group5_source<-polarPlot(timeseries_CHON5_group5, pollutant = c("sum"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group6[is.na(timeseries_CHON5_group6)] <- 0
timeseries_CHON5_group6$sum <- rowSums( timeseries_CHON5_group6[2:5], na.rm=TRUE )
timeseries_CHON5_group6$unixtime<-as.POSIXct(timeseries_CHON5_group6$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group6, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group6") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
ggplot(timeseries_CHON5_group6, aes(unixtime, "C11H9NO3 - RT8.64")) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("C11H9NO3 - RT8.64 timeseries") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
timeseries_CHON5_group6 <- timeseries_CHON5_group6 %>% 
  left_join(avgmetwinter2, by = "unixtime")

CHON_PCA_group6_source<-polarPlot(timeseries_CHON5_group6, pollutant = c("sum"), k = 50, key.header = "Peak area (a.u. / m3)")


wintertimeseries_transposed_BTNratio <- wintertimeseries_transposed %>% 
  select('date_mid',
         'Benzene_ptrtof',
         'Toluene_ptrtof',
         'MLH_m_agl')

colnames(wintertimeseries_transposed_BTNratio)[4] <- "BLheight"

wintertimeseries_transposed_BTNratio$BTratio <- 0
{for(x in 1:length(wintertimeseries_transposed_BTNratio$BTratio)) wintertimeseries_transposed_BTNratio$BTratio[x] <- ((wintertimeseries_transposed_BTNratio$Benzene_ptrtof[x])/(wintertimeseries_transposed_BTNratio$Toluene_ptrtof[x]))}
wintertimeseries_transposed_BTNratio$date_mid <- dmy_hm(wintertimeseries_transposed_BTNratio$date_mid, tz = "UTC")
wintertimeseries_transposed_BTNratio <- wintertimeseries_transposed_BTNratio[ which(wintertimeseries_transposed_BTNratio$date_mid > "2016-11-23 00:00:00"),]
ggplot(wintertimeseries_transposed_BTNratio, aes(date_mid, BTratio)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ylim(0,2) +ggtitle("Benzene : Toluene Ratio") + guides(size=FALSE) + ylab("B:T ratio") + xlab("Date")

ggplot(wintertimeseries_transposed_BTNratio, aes(date_mid, BLheight)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                        panel.grid.major = element_blank(),
                                                                                                        panel.grid.minor = element_blank(),
                                                                                                        panel.background = element_blank()) +ggtitle("Boundary layer height (m)") + guides(size=FALSE) + ylab("BL height (m)") + xlab("Date")


timeseries_CHON_proportion1 <- timeseries_CHON5_group1 %>% 
  select(unixtime,
         sum) 

colnames(timeseries_CHON_proportion1)[2] <- "Group1"

timeseries_CHON_proportion2 <- timeseries_CHON5_group2 %>% 
  select(unixtime,
         sum) 

colnames(timeseries_CHON_proportion2)[2] <- "Group2"

timeseries_CHON_proportion3 <- timeseries_CHON5_group3 %>% 
  select(unixtime,
         sum) 

colnames(timeseries_CHON_proportion3)[2] <- "Group3"

timeseries_CHON_proportion4 <- timeseries_CHON5_group4 %>% 
  select(unixtime,
         sum) 

colnames(timeseries_CHON_proportion4)[2] <- "Group4"

timeseries_CHON_proportion5 <- timeseries_CHON5_group5 %>% 
  select(unixtime,
         sum) 

colnames(timeseries_CHON_proportion5)[2] <- "Group5"

timeseries_CHON_proportion6 <- timeseries_CHON5_group6 %>% 
  select(unixtime,
         sum) 

colnames(timeseries_CHON_proportion6)[2] <- "Group6"

timeseries_CHON_proportion<-timeseries_CHON_proportion6 %>% 
  right_join(timeseries_CHON_proportion5, by = "unixtime") %>% 
  right_join(timeseries_CHON_proportion4, by = "unixtime") %>% 
  right_join(timeseries_CHON_proportion3, by = "unixtime") %>% 
  right_join(timeseries_CHON_proportion2, by = "unixtime") %>% 
  right_join(timeseries_CHON_proportion1, by = "unixtime")

timeseries_CHON_proportion$sum <- rowSums( timeseries_CHON_proportion[2:7] )

timeseries_CHON_proportion$group1 <- (timeseries_CHON_proportion$Group1/timeseries_CHON_proportion$sum)*100
timeseries_CHON_proportion$group2 <- (timeseries_CHON_proportion$Group2/timeseries_CHON_proportion$sum)*100
timeseries_CHON_proportion$group3 <- (timeseries_CHON_proportion$Group3/timeseries_CHON_proportion$sum)*100
timeseries_CHON_proportion$group4 <- (timeseries_CHON_proportion$Group4/timeseries_CHON_proportion$sum)*100
timeseries_CHON_proportion$group5 <- (timeseries_CHON_proportion$Group5/timeseries_CHON_proportion$sum)*100
timeseries_CHON_proportion$group6 <- (timeseries_CHON_proportion$Group6/timeseries_CHON_proportion$sum)*100

timeseries_CHON_proportion$sum2 <- rowSums( timeseries_CHON_proportion[9:12] )

write.csv(timeseries_CHON_proportion, file="CHONclusterproportion.csv")
CHONclusterproportion_2 <- read.csv("C:/Users/William/Documents/WACL stuff/CHONclusterproportion.csv", na="NA", check.names = F)
CHONclusterproportion_2$unixtime <- dmy_hm(CHONclusterproportion_2$unixtime, tz = "UTC")

gg <- ggplot(CHONclusterproportion_2, aes(unixtime,Prop), y=Prop)
gg <- gg + geom_area(aes(colour=Group, fill=Group))
gg

slices <- c(0.4216, 95.4896, 4.08879) 
lbls <- c("Groups 1, 2, 3 & 4", "Group 5", "Group 6")
pie(slices,labels=lbls,explode=0.1,
      main="Proportion of CHON by cluster")

timeseries_CHON5_group1$prop1 <- (timeseries_CHON5_group1$`2-Nitrophenol`/timeseries_CHON5_group1$sum)*100
timeseries_CHON5_group1$prop2 <- (timeseries_CHON5_group1$`C10H12OHNO3 - RT7.07`/timeseries_CHON5_group1$sum)*100

mean(timeseries_CHON5_group1$prop1, na.rm=TRUE)

slices <- c(8.633243, 8.1667) 
lbls <- c("2-Nitrophenol", "C10H12OHNO3 - RT7.07")
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 1")

timeseries_CHON5_group2$prop1 <- (timeseries_CHON5_group2$`C6H4OHNO3 - RT11.44`/timeseries_CHON5_group2$sum)*100
timeseries_CHON5_group2$prop2 <- (timeseries_CHON5_group2$`C7H6OHNO3 - RT6.55`/timeseries_CHON5_group2$sum)*100
timeseries_CHON5_group2$prop3 <- (timeseries_CHON5_group2$`C16H10O4N2 - RT6.98`/timeseries_CHON5_group2$sum)*100
timeseries_CHON5_group2$prop4 <- (timeseries_CHON5_group2$`C18H25O3NO3 - RT6.35`/timeseries_CHON5_group2$sum)*100

mean(timeseries_CHON5_group2$prop1, na.rm=TRUE)
mean(timeseries_CHON5_group2$prop2, na.rm=TRUE)
mean(timeseries_CHON5_group2$prop3, na.rm=TRUE)
mean(timeseries_CHON5_group2$prop4, na.rm=TRUE)

slices <- c(66.31, 23.025, 7.151, 3.513) 
lbls <- c("C6H4OHNO3 - RT11.44", "C7H6OHNO3 - RT6.55", "C16H10O4N2 - RT6.98", "C18H25O3NO3 - RT6.35")
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 2")

timeseries_CHON5_group3$prop1 <- (timeseries_CHON5_group3$`C6H4OHNO3 - RT10.82`/timeseries_CHON5_group3$sum)*100
timeseries_CHON5_group3$prop2 <- (timeseries_CHON5_group3$`C8H8OHNO3 - RT5.00`/timeseries_CHON5_group3$sum)*100
timeseries_CHON5_group3$prop3 <- (timeseries_CHON5_group3$`C7H7NO4 - RT5.07`/timeseries_CHON5_group3$sum)*100
timeseries_CHON5_group3$prop4 <- (timeseries_CHON5_group3$`C8H8OHNO3 - RT6.50`/timeseries_CHON5_group3$sum)*100


timeseries_CHON5_group3 <- timeseries_CHON5_group3[which(timeseries_CHON5_group3$sum>0 ),]

mean(timeseries_CHON5_group3$prop1, na.rm=TRUE)
mean(timeseries_CHON5_group3$prop2, na.rm=TRUE)
mean(timeseries_CHON5_group3$prop3, na.rm=TRUE)
mean(timeseries_CHON5_group3$prop4, na.rm=TRUE)

slices <- c(37.765, 59.47, 40.27, 30.953) 
lbls <- c("C6H4OHNO3 - RT10.82", "C8H8OHNO3 - RT5.00", "C7H7NO4 - RT5.07", "C8H8OHNO3 - RT6.50")
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 3")

timeseries_CHON5_group4$prop1 <- (timeseries_CHON5_group4$`C7H7NO4 - RT11.58`/timeseries_CHON5_group4$sum)*100
timeseries_CHON5_group4$prop2 <- (timeseries_CHON5_group4$`C8H8OHNO3 - RT8.52`/timeseries_CHON5_group4$sum)*100
timeseries_CHON5_group4$prop3 <- (timeseries_CHON5_group4$`C8H8OHNO3 - RT7.87`/timeseries_CHON5_group4$sum)*100
timeseries_CHON5_group4$prop4 <- (timeseries_CHON5_group4$`4-Methoxy-2-nitrophenol`/timeseries_CHON5_group4$sum)*100

timeseries_CHON5_group4 <- timeseries_CHON5_group4[which(timeseries_CHON5_group4$sum>0 ),]

mean(timeseries_CHON5_group4$prop1, na.rm=TRUE)
mean(timeseries_CHON5_group4$prop2, na.rm=TRUE)
mean(timeseries_CHON5_group4$prop3, na.rm=TRUE)
mean(timeseries_CHON5_group4$prop4, na.rm=TRUE)

slices <- c(59.137, 26.231, 26.3, 64.9) 
lbls <- c("C7H7NO4 - RT11.58", "C8H8OHNO3 - RT8.52", "C8H8OHNO3 - RT7.87", "4-Methoxy-2-nitrophenol")
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 4")


timeseries_CHON5_group5$prop1 <- (timeseries_CHON5_group5$`C7H7NO3 - RT6.86`/timeseries_CHON5_group5$sum)*100
timeseries_CHON5_group5$prop2 <- (timeseries_CHON5_group5$`C8H5NO3 - RT3.65`/timeseries_CHON5_group5$sum)*100
timeseries_CHON5_group5$prop3 <- (timeseries_CHON5_group5$`C7H7NO4 - RT3.97`/timeseries_CHON5_group5$sum)*100
timeseries_CHON5_group5$prop4 <- (timeseries_CHON5_group5$`4-Nitrophenol`/timeseries_CHON5_group5$sum)*100
timeseries_CHON5_group5$prop5 <- (timeseries_CHON5_group5$`C6H5NO3 - RT5.14`/timeseries_CHON5_group5$sum)*100
timeseries_CHON5_group5$prop6 <- (timeseries_CHON5_group5$`C7H7NO4 - RT5.62`/timeseries_CHON5_group5$sum)*100
timeseries_CHON5_group5$prop7 <- (timeseries_CHON5_group5$`4-Methyl-3-nitrophenol`/timeseries_CHON5_group5$sum)*100
timeseries_CHON5_group5$prop8 <- (timeseries_CHON5_group5$`2,6-dimethyl-4-nitrophenol`/timeseries_CHON5_group5$sum)*100

timeseries_CHON5_group5 <- timeseries_CHON5_group5[which(timeseries_CHON5_group5$sum>0 ),]

mean(timeseries_CHON5_group5$prop1, na.rm=TRUE)
mean(timeseries_CHON5_group5$prop2, na.rm=TRUE)
mean(timeseries_CHON5_group5$prop3, na.rm=TRUE)
mean(timeseries_CHON5_group5$prop4, na.rm=TRUE)
mean(timeseries_CHON5_group5$prop5, na.rm=TRUE)
mean(timeseries_CHON5_group5$prop6, na.rm=TRUE)
mean(timeseries_CHON5_group5$prop7, na.rm=TRUE)
mean(timeseries_CHON5_group5$prop8, na.rm=TRUE)

slices <- c(14.099, 0.519, 0.601, 32.857, 32.282, 0.828, 14.09, 7.42) 
lbls <- c("C7H7NO3 - RT6.86", "C8H5NO3 - RT3.65", "C7H7NO4 - RT3.97", "4-Nitrophenol", "C6H5NO3 - RT5.14", "C7H7NO4 - RT5.62", "4-Methyl-3-nitrophenol", "2,6-dimethyl-4-nitrophenol")
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 5", radius=1.2, cex=0.4)

timeseries_CHON5_group6$prop1 <- (timeseries_CHON5_group6$`C10H7NO3 - RT8.25`/timeseries_CHON5_group6$sum)*100
timeseries_CHON5_group6$prop2 <- (timeseries_CHON5_group6$`C11H9NO3 - RT8.64`/timeseries_CHON5_group6$sum)*100
timeseries_CHON5_group6$prop3 <- (timeseries_CHON5_group6$`C9H5O4N - RT3.09`/timeseries_CHON5_group6$sum)*100
timeseries_CHON5_group6$prop4 <- (timeseries_CHON5_group6$`C9H5O4N - RT3.95`/timeseries_CHON5_group6$sum)*100

timeseries_CHON5_group6 <- timeseries_CHON5_group6[which(timeseries_CHON5_group6$sum>0 ),]

mean(timeseries_CHON5_group6$prop1, na.rm=TRUE)
mean(timeseries_CHON5_group6$prop2, na.rm=TRUE)
mean(timeseries_CHON5_group6$prop3, na.rm=TRUE)
mean(timeseries_CHON5_group6$prop4, na.rm=TRUE)

slices <- c(65.699, 22.459, 9.508, 4.15) 
lbls <- c("2-Nitro1-naphthol", "C11H9NO3 - RT8.64", "C9H5O4N - RT3.09", "C9H5O4N - RT3.95")
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 6")

avgmetwinter_82 <- avgmetwinter_8%>%
  select(wd,ws,datetime) 

avgmetwinter_82$wd<-as.numeric(avgmetwinter_82$wd)
avgmetwinter_82$ws<-as.numeric(avgmetwinter_82$ws)

avgmetwinter_82$datetime  <- ymd_hms(avgmetwinter_82$datetime, tz = "UTC")
avgmetwinter_83<-avgmetwinter_82 %>% 
  right_join(wintertimeseries, by = "datetime")

colnames(avgmetwinter_83)[3] <- "date"
polar_group2<-polarPlot(avgmetwinter_83, pollutant = c("2-Nitro-1-Naphthol"), k = 50, key.header = "Peak area (a.u. / m3)")
polar_group3<-polarPlot(avgmetwinter_83, pollutant = c("2-Methyl-5-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")
polar_group1<-polarPlot(avgmetwinter_83, pollutant = c("4-Methyl-3-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")
polar_group4<-polarPlot(avgmetwinter_83, pollutant = c("4-Methoxy-2-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")

polarCluster(avgmetwinter_83, pollutant = c("2-Nitro-1-Naphthol"), k = 50, statistic = "cpf", percentile = 90, n.clusters = 2)


#all r corr plot, this may take a while
filepath = paste(getwd(),"/allRcorr.pdf",sep = "")
pdf(width=80, height=80, file=filepath)
corPlot(avgmetwinter_83,
        pollutants = names(avgmetwinter_83)[c(4:621)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of library factors",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

nits = unique(Nitcomps[,1])

CHONs<- avgmetwinter_83[,grep(c("itro"), names(avgmetwinter_83))]

avgmetwinter_84 <- avgmetwinter_83%>%
  select(wd,ws,date,'2,6-dimethyl-4-nitrophenol')

CHONs<-CHONs %>% 
  right_join(avgmetwinter_84, by = "2,6-dimethyl-4-nitrophenol")

#remove 3np due to low data content
CHONs<-CHONs[,-c(2,6)] 

filepath = paste(getwd(),"/nitrophenRcorr.pdf",sep = "")
pdf(width=10, height=10, file=filepath)
corPlot(CHONs,
        pollutants = names(CHONs)[c(1:10)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified nitro-aromatics",
        order="alphabet",
        cluster=TRUE)
dev.off()

CHONs_r_plot <- corPlot(CHONs,
                        pollutants = names(CHONs)[c(1:10)],
                        cols = "jet",
                        dendrogram = TRUE,
                        main = "Correlations of identified nitro-aromatics",
                        order="alphabet",
                        cluster=TRUE)
CHONs_r_dend <- as.dendrogram(CHONs_r_plot$clust)

CHOS_samples <- Samples[which(Samples$Group=="CHOS" ),]
CHON_samples <- Samples[which(Samples$Group=="CHON" ),]
CHONS_samples <- Samples[which(Samples$Group=="CHONS" ),]
CHO_samples <- Samples[which(Samples$Group=="CHO" ),]

CHOS_compounds = unique(CHOS_samples[,2])
CHONS_compounds = unique(CHONS_samples[,2])
CHON_compounds = unique(CHON_samples[,2])
CHO_compounds = unique(CHO_samples[,2])

timeseries_CHON2 = seq(1478735940,1481241540,60)
timeseries_CHON3 = seq(1478735940,1481241540,60)
timeseries_CHON2 = data.frame(unixtime = timeseries_CHON2)
timeseries_CHON3 = data.frame(unixtime = timeseries_CHON3)
CHON_samples$unixtime = as.numeric(CHON_samples$datetime)
CHON_samples$Quant = as.numeric(CHON_samples$Quant)

for (i in 2:length(CHON_compounds)){
  temp_compound2 = CHON_samples[CHON_samples$Compound == CHON_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHON_compounds[i],"unixtime")
  
  timeseries_CHON2 = right_join(temp_compound2,timeseries_CHON2,"unixtime")
}
timeseries_CHON3 <- time %>% 
  left_join(timeseries_CHON2, by = "unixtime")

timeseries_CHONS2 = seq(1478735940,1481241540,60)
timeseries_CHONS3 = seq(1478735940,1481241540,60)
timeseries_CHONS2 = data.frame(unixtime = timeseries_CHONS2)
timeseries_CHONS3 = data.frame(unixtime = timeseries_CHONS3)
CHONS_samples$unixtime = as.numeric(CHONS_samples$datetime)
CHONS_samples$Quant = as.numeric(CHONS_samples$Quant)

for (i in 2:length(CHONS_compounds)){
  temp_compound2 = CHONS_samples[CHONS_samples$Compound == CHONS_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHONS_compounds[i],"unixtime")
  
  timeseries_CHONS2 = right_join(temp_compound2,timeseries_CHONS2,"unixtime")
}
timeseries_CHONS3 <- time %>% 
  left_join(timeseries_CHONS2, by = "unixtime")

timeseries_CHO2 = seq(1478735940,1481241540,60)
timeseries_CHO3 = seq(1478735940,1481241540,60)
timeseries_CHO2 = data.frame(unixtime = timeseries_CHO2)
timeseries_CHO3 = data.frame(unixtime = timeseries_CHO3)
CHO_samples$unixtime = as.numeric(CHO_samples$datetime)
CHO_samples$Quant = as.numeric(CHO_samples$Quant)

for (i in 2:length(CHO_compounds)){
  temp_compound2 = CHO_samples[CHO_samples$Compound == CHO_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHO_compounds[i],"unixtime")
  
  timeseries_CHO2 = right_join(temp_compound2,timeseries_CHO2,"unixtime")
}
timeseries_CHO3 <- time %>% 
  left_join(timeseries_CHO2, by = "unixtime")

timeseries_CHOS2 = seq(1478735940,1481241540,60)
timeseries_CHOS3 = seq(1478735940,1481241540,60)
timeseries_CHOS2 = data.frame(unixtime = timeseries_CHOS2)
timeseries_CHOS3 = data.frame(unixtime = timeseries_CHOS3)
CHOS_samples$unixtime = as.numeric(CHOS_samples$datetime)
CHOS_samples$Quant = as.numeric(CHOS_samples$Quant)

for (i in 2:length(CHOS_compounds)){
  temp_compound2 = CHOS_samples[CHOS_samples$Compound == CHOS_compounds[i],]
  temp_compound2 = temp_compound2[,c("Quant","unixtime")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(CHOS_compounds[i],"unixtime")
  
  timeseries_CHOS2 = right_join(temp_compound2,timeseries_CHOS2,"unixtime")
}
timeseries_CHOS3 <- time %>% 
  left_join(timeseries_CHOS2, by = "unixtime")

filepath = paste(getwd(),"/CHOSRcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHOS3,
        pollutants = names(timeseries_CHOS3)[c(2:188)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHOS",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

timeseries_CHON4<-timeseries_CHON3 %>% 
  select(-contains('C10H19OH(NO3)2 - RT10'),
         -contains('C19H29O3NO3 - RT7.4'),
         -contains('C18H25O3NO3 - RT8.5'),
         -contains('18H25O3NO3 - RT6.2'),
         -contains('C11H19NO9 - RT5.3'),
         -contains('2-Methyl-3-nitropheno'))

filepath = paste(getwd(),"/CHONRcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHON4,
        pollutants = names(timeseries_CHON4)[c(2:40)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHON",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

filepath = paste(getwd(),"/CHORcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHO3,
        pollutants = names(timeseries_CHO3)[c(2:231)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHO",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

filepath = paste(getwd(),"/CHONSRcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHONS3,
        pollutants = names(timeseries_CHONS3)[c(2:105)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHONS",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

CHONSs_r_plot <-corPlot(timeseries_CHONS3,
                        pollutants = names(timeseries_CHONS3)[c(2:45)],
                        cols = "jet",
                        dendrogram = TRUE,
                        main = "Correlations of identified CHON",
                        order="alphabet",
                        cluster=TRUE,
                        text.col = NA)
dev.off()
CHONSs_r_dend <- as.dendrogram(CHONSs_r_plot$clust)

timeseries_CHOS4<-timeseries_CHOS3
timeseries_CHOS4[is.na(timeseries_CHOS4)] <- 1
log_CHOS <- log(timeseries_CHOS4[, 2:188])

timeseries_CHON4<-timeseries_CHON3
timeseries_CHON4[is.na(timeseries_CHON4)] <- 1
log_CHON <- log(timeseries_CHON4[, 2:45])

timeseries_CHONS4<-timeseries_CHONS3
timeseries_CHONS4[is.na(timeseries_CHONS4)] <- 1
log_CHONS <- log(timeseries_CHONS4[, 2:105])

timeseries_CHO4<-timeseries_CHO3
timeseries_CHO4[is.na(timeseries_CHO4)] <- 1
log_CHO <- log(timeseries_CHO4[, 2:231])

CHO_pca <- prcomp(log_CHO,
                  center = TRUE,
                  scale. = TRUE) 
plot(CHO_pca, type = "l")

CHON_pca <- prcomp(log_CHON,
                   center = TRUE,
                   scale. = TRUE) 
plot(CHON_pca, type = "l")

CHOS_pca <- prcomp(log_CHOS,
                   center = TRUE,
                   scale. = TRUE) 
plot(CHOS_pca, type = "l")

CHONS_pca <- prcomp(log_CHONS,
                    center = TRUE,
                    scale. = TRUE) 
plot(CHONS_pca, type = "l")

fviz_pca_var(CHON_pca, col.var="contrib")+
  scale_color_gradient2(low="blue", mid="green",
                        high="red", midpoint=3) +
  theme_minimal()

# log transform 
avgmetwinter_84<-avgmetwinter_83
avgmetwinter_84[is.na(avgmetwinter_84)] <- 1
log_all <- log(avgmetwinter_84[, 4:621])
#ir.species <- iris[, 5]

# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
all_pca <- prcomp(log_all,
                  center = TRUE,
                  scale. = TRUE) 
princomp(log_all, cor = FALSE, scores = TRUE)

CHON_PCA_group1_source<-polarPlot(avgmetwinter_83, pollutant = c("2,6-dimethyl-4-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")
CHON_PCA_group2_source<-polarPlot(avgmetwinter_83, pollutant = c("C7H7NO4 - RT5.07"), k = 50, key.header = "Peak area (a.u. / m3)")
CHON_PCA_group3_source<-polarPlot(avgmetwinter_83, pollutant = c("2-Nitro-1-Naphthol"), k = 50, key.header = "Peak area (a.u. / m3)")

CHOS_PCA_group1_source<-polarPlot(avgmetwinter_83, pollutant = c("C9H12O5S - RT3.86"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group2_source<-polarPlot(avgmetwinter_83, pollutant = c("C5H8SO5 - RT0.94"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group3_source<-polarPlot(avgmetwinter_83, pollutant = c("C16H28O4SO4 - RT6.26"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group5_source<-polarPlot(avgmetwinter_83, pollutant = c("C5H6O4SO4 - RT0.71"), k = 50, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group5_source<-polarCluster(avgmetwinter_83, pollutant = c("C5H6O4SO4 - RT0.71"), k = 50, key.header = "Peak area (a.u. / m3)", n.clusters = 3)
CHOS_PCA_group6_source<-polarPlot(avgmetwinter_83, pollutant = c("C10H20OSO4 - RT6.71"), k = 50, key.header = "Peak area (a.u. / m3)")

CHO_PCA_group1_source<-polarPlot(avgmetwinter_83, pollutant = c("C9H8O4 - RT5.25"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group2_source<-polarPlot(avgmetwinter_83, pollutant = c("C7H6O2 - RT3.64"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group3_source<-polarPlot(avgmetwinter_83, pollutant = c("C8H6O4 - RT4.06"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group4_source<-polarPlot(avgmetwinter_83, pollutant = c("C6H10O4 - RT1.18 (Pinene SOA)"), k = 50, key.header = "Peak area (a.u. / m3)")
CHO_PCA_group5_source<-polarPlot(avgmetwinter_83, pollutant = c("nor-caryophyllenic acid (B-Caryophyllene SOA C8H12O4)"), k = 50, key.header = "Peak area (a.u. / m3)")


wintertimeseries_transposed <- read.csv("C:/Users/William/Documents/WACL stuff/wintertimeres_transposed.csv", na="NA", check.names = F)

wintertimeseries_transposed_1<-wintertimeseries_transposed
wintertimeseries_transposed_1[is.na(wintertimeseries_transposed_1)] <- 1
log_CHOS_t <- log(wintertimeseries_transposed_1[, 2:124])
CHOS_t_pca <- prcomp(log_CHOS_t,
                     center = TRUE,
                     scale. = TRUE) 

CHON_dend <- corPlot(timeseries_CHON3,
                     pollutants = names(timeseries_CHON3)[c(2:45)],
                     cols = "jet",
                     dendrogram = TRUE,
                     main = "Correlations of identified CHON",
                     order="alphabet",
                     cluster=TRUE,
                     text.col = NA)
dev.off()

CHON_dend_table<-Map(as.data.frame(CHON_dend))

mydendCHON <- as.dendrogram(CHON_dend$clust)
plot(mydendCHON)
k <- 6
mydendCHON <- color_branches(mydendCHON, k = k)
plot(mydendCHON)
labels_dend <- labels(mydendCHON)
groups <- cutree(mydendCHON, k=6, order_clusters_as_data = FALSE)
dends <- list()
for(i in 1:k) {
  labels_to_keep <- labels_dend[i != groups]
  dends[[i]] <- prune(dend, labels_to_keep)}

CHON_dend_table <-as.data.frame(CHON_dend$data)
CHON_dend_table <- CHON_dend_table[ which(CHON_dend_table$cor > -1),]
CHON_dend_table <- CHON_dend_table[ which(CHON_dend_table$cor < 1),]

n = unique(CHON_dend_table[,1])
CHON_dend_table2 <- CHON_dend_table %>% 
  select(x)

for (i in 2:length(n)){
  temp_compound2 = CHON_dend_table[CHON_dend_table$y == n[i],]
  temp_compound2 = temp_compound2[,c("cor","x")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(n[i],"x")
  
  timeseries3 = right_join(temp_compound2,CHON_dend_table2,"x")
}

timeseries3 <- time %>% 
  left_join(timeseries2, by = "unixtime")


m_dist<-dist(log_CHOS_t,diag = FALSE )
m_hclust<-hclust(m_dist, method= "complete")
plot(m_hclust)

# Source apportionment (saucy portion)
k <- 6
cols <- rainbow_hcl(k)
dend <- as.dendrogram(m_hclust)
dend <- color_branches(dend, k = k)
plot(dend)
labels_dend <- labels(dend)
groups <- cutree(dend, k=6, order_clusters_as_data = FALSE)
dends <- list()
for(i in 1:k) {
  labels_to_keep <- labels_dend[i != groups]
  dends[[i]] <- prune(dend, labels_to_keep)
}

length(which(timeseries_CHON3$`C8H8OHNO3 - RT8.52` != "NA")) 

timeseries_CHON5<-timeseries_CHON3
timeseries_CHON5[is.na(timeseries_CHON5)] <- 0
colSums(timeseries_CHON5 != 0)
timeseries_CHON5 <- timeseries_CHON5[, !colSums(timeseries_CHON5 != 0) < 5]
timeseries_CHON5[timeseries_CHON5 < 0.1] <- NA

timeseries_CHON6<-timeseries_CHON5 %>% 
  select(-contains('C10H19OH(NO3)2 - RT10'),
         -contains('C19H29O3NO3 - RT7.4'),
         -contains('C18H25O3NO3 - RT8.5'),
         -contains('18H25O3NO3 - RT6.2'),
         -contains('C11H19NO9 - RT5.3'),
         -contains('2-Methyl-3-nitropheno'))

filepath = paste(getwd(),"/CHONRcorr2.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHON6,
        pollutants = names(timeseries_CHON6)[c(2:27)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHON",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

dend <- plot(as.phylo(outCHON$clust), cex=0.7)

wintertimeseries_transposed <- read.csv("C:/Users/William/Documents/WACL stuff/Wintertimeres.csv", na="NA", check.names = F)

wintertimeseries_transposed2 <- wintertimeseries_transposed %>% 
  select(datetime,
         '4-Methoxy-2-nitrophenol',
         '4-Methyl-2-nitrophenol',
         '2-Methyl-3-nitrophenol', 
         '4-Nitrocatechol',
         '2-Nitrophenol')

wintertimeseries_transposed2$datetime <- dmy_hm(wintertimeseries_transposed2$datetime, tz = "UTC")

AMS_STD_Comp <- AMSnit %>% 
  select(datetime,
         CCOA,
         COA,
         BBOA,
         OPOA,
         LOOOA,
         MOOOA,
         Chl,
         NH4,
         NO3,
         SO4,
         Org,
         '2,4-Dinitrophenol',
         '2,6-dimethyl-4-nitrophenol',
         '3-Methyl-2-nitrophenol',
         '4-Methyl-3-nitrophenol',
         '4-Nitrophenol',
         'C10H7NO3 - RT8.25')

wintertimeseries_transposed3 <- wintertimeseries_transposed2 %>% 
  right_join(AMS_STD_Comp, by = "datetime")

colnames(wintertimeseries_transposed3)[23] <- "2-Nitro-1-naphthol"


filepath = paste(getwd(),"/CHONRcorrAMS.pdf",sep = "")
pdf(width=20, height=20, file=filepath)
corPlot(wintertimeseries_transposed3,
        pollutants = names(wintertimeseries_transposed3)[c(2:23)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHON with AMS",
        order="alphabet",
        cluster=TRUE)
dev.off()

wintertimeseries_transposed <- read.csv("C:/Users/William/Documents/WACL stuff/aphh_winter_filter_aggregate_merge.csv", na="NA", check.names = F)

timeseries_CHON5_group1 <- timeseries_CHON5 %>% 
  select(unixtime,
         '2-Nitrophenol',
         'C10H12OHNO3 - RT7.07')

timeseries_CHON5_group4 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C7H7NO4 - RT11.58',
         'C8H8OHNO3 - RT8.52',
         'C8H8OHNO3 - RT7.87',
         '4-Methoxy-2-nitrophenol')

timeseries_CHON5_group2 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C6H4OHNO3 - RT11.44',
         'C7H6OHNO3 - RT6.55',
         'C16H10O4N2 - RT6.98',
         'C18H25O3NO3 - RT6.35')

timeseries_CHON5_group3 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C6H4OHNO3 - RT10.82',
         'C8H8OHNO3 - RT5.00',
         'C7H7NO4 - RT5.07',
         'C8H8OHNO3 - RT6.50')

timeseries_CHON5_group6 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C10H7NO3 - RT8.25',
         'C11H9NO3 - RT8.64',
         'C9H5O4N - RT3.09',
         'C9H5O4N - RT3.95')

timeseries_CHON5_group5 <- timeseries_CHON5 %>% 
  select(unixtime,
         'C7H7NO3 - RT6.86',
         'C8H5NO3 - RT3.65',
         'C7H7NO4 - RT3.97',
         '4-Nitrophenol',
         'C6H5NO3 - RT5.14',
         'C7H7NO4 - RT5.62',
         '4-Methyl-3-nitrophenol',
         '2,6-dimethyl-4-nitrophenol')

#timeseries_CHON5_group1[is.na(timeseries_CHON5_group1)] <- 0
timeseries_CHON5_group1$sum <- rowSums( timeseries_CHON5_group1[2:3] , na.rm=TRUE)
timeseries_CHON5_group1$unixtime<-as.POSIXct(timeseries_CHON5_group1$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group1, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group1") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
CHON_PCA_group1_source<-polarPlot(avgmetwinter_83, pollutant = c("2-Nitrophenol"), k = 10, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group2[is.na(timeseries_CHON5_group2)] <- 0
timeseries_CHON5_group2$sum <- rowSums( timeseries_CHON5_group2[2:5] , na.rm=TRUE)
timeseries_CHON5_group2$unixtime<-as.POSIXct(timeseries_CHON5_group2$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group2, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group2") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
CHON_PCA_group2_source<-polarPlot(avgmetwinter_83, pollutant = c("C6H4OHNO3 - RT11.44"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group3[is.na(timeseries_CHON5_group3)] <- 0
timeseries_CHON5_group3$sum <- rowSums( timeseries_CHON5_group3[2:4], na.rm=TRUE)
timeseries_CHON5_group3$unixtime<-as.POSIXct(timeseries_CHON5_group3$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group3, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group3") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
CHON_PCA_group3_source<-polarPlot(avgmetwinter_83, pollutant = c("C8H8OHNO3 - RT5.00"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group4[is.na(timeseries_CHON5_group4)] <- 0
timeseries_CHON5_group4$sum <- rowSums( timeseries_CHON5_group4[2:5], na.rm=TRUE )
timeseries_CHON5_group4$unixtime<-as.POSIXct(timeseries_CHON5_group4$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group4, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group4") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
CHON_PCA_group4_source<-polarPlot(avgmetwinter_83, pollutant = c("C8H8OHNO3 - RT8.52"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group5[is.na(timeseries_CHON5_group5)] <- 0
timeseries_CHON5_group5$sum <- rowSums( timeseries_CHON5_group5[2:9] , na.rm=TRUE)
timeseries_CHON5_group5$unixtime<-as.POSIXct(timeseries_CHON5_group5$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group5, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group5") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
CHON_PCA_group5_source<-polarPlot(avgmetwinter_83, pollutant = c("4-Methyl-3-nitrophenol"), k = 50, key.header = "Peak area (a.u. / m3)")

#timeseries_CHON5_group6[is.na(timeseries_CHON5_group6)] <- 0
timeseries_CHON5_group6$sum <- rowSums( timeseries_CHON5_group6[2:5], na.rm=TRUE )
timeseries_CHON5_group6$unixtime<-as.POSIXct(timeseries_CHON5_group6$unixtime, origin="1970-01-01 00:00:00")
ggplot(timeseries_CHON5_group6, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Average group6") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
ggplot(timeseries_CHON5_group6, aes(unixtime, "C11H9NO3 - RT8.64")) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                       panel.grid.major = element_blank(),
                                                                                                       panel.grid.minor = element_blank(),
                                                                                                       panel.background = element_blank()) +ggtitle("C11H9NO3 - RT8.64 timeseries") + guides(size=FALSE) + ylab("Sum peak area") + xlab("Date")
CHON_PCA_group6_source<-polarPlot(avgmetwinter_83, pollutant = c("C11H9NO3 - RT8.64"), k = 50, key.header = "Peak area (a.u. / m3)")

GC_aggregate <- read.csv("C:/Users/William/Documents/WACL stuff/aphh_winter_gc_aggregate_merge.csv", na="NA", check.names = F)
colnames(GC_aggregate)[66] <- "BLheight"
GC_aggregate$date_mid <- ymd_hms(GC_aggregate$date_mid, tz = "UTC")
ggplot(GC_aggregate, aes(date_mid, BLheight)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                                         panel.grid.major = element_blank(),
                                                                                                         panel.grid.minor = element_blank(),
                                                                                                         panel.background = element_blank()) +ggtitle("Mixing layer height") + guides(size=FALSE) + ylab("ML height (m)") + xlab("Date")

log_CHOS_t <- log(wintertimeseries_transposed_1[, 2:124])
CHOS_t_pca <- prcomp(log_CHOS_t,
                     center = TRUE,
                     scale. = TRUE) 

timeseries_CHOS3_QA<-timeseries_CHOS3
timeseries_CHOS3_QA[is.na(timeseries_CHOS3_QA)] <- 0
timeseries_CHOS3_QA<-timeseries_CHOS3_QA[, !colSums(timeseries_CHOS3_QA != 0) < 40]
timeseries_CHOS3_QA[timeseries_CHOS3_QA < 0.1] <- NA

timeseries_CHOS3_QA <- timeseries_CHOS3_QA %>% 
  select(-'C5H12SO7 - RT0.75',
         -'C3H6SO6 - RT0.78')

CHOS_dend <- corPlot(timeseries_CHOS3_QA,
                     pollutants = names(timeseries_CHOS3_QA)[c(2:82)],
                     cols = "jet",
                     dendrogram = TRUE,
                     main = "Correlations of identified CHOS",
                     order="alphabet",
                     cluster=TRUE,
                     text.col = NA)
dev.off()



filepath = paste(getwd(),"/CHOSRcorr.pdf",sep = "")
pdf(width=20, height=20, file=filepath)
corPlot(timeseries_CHOS3_QA,
        pollutants = names(timeseries_CHOS3_QA)[c(2:82)],
        cols = "jet",
        dendrogram = TRUE,
        main = "Correlations of identified CHOS",
        order="alphabet",
        cluster=TRUE,
        text.col = NA)
dev.off()

CHOS_dend_table<-Map(as.data.frame(CHOS_dend))

mydendCHOS <- as.dendrogram(CHOS_dend$clust)
plot(mydendCHOS)
k <- 12
mydendCHOS <- color_branches(mydendCHOS, k = k)
plot(mydendCHOS)
labels_dend <- labels(mydendCHOS)
groups <- cutree(mydendCHOS, k=6, order_clusters_as_data = FALSE)
dends <- list()
for(i in 1:k) {
  labels_to_keep <- labels_dend[i != groups]
  dends[[i]] <- prune(dend, labels_to_keep)}

colors = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000", "FF6666", "3333FF")
clus12 = cutree(mydendCHOS, 11)
plot(as.phylo(mydendCHOS), type = "fan", tip.color = colors[clus12],
     label.offset = 0.3, cex = 0.7)

CHOS_dend_table <-as.data.frame(CHOS_dend$data)
CHOS_dend_table <- CHOS_dend_table[ which(CHOS_dend_table$cor > -1),]
CHOS_dend_table <- CHOS_dend_table[ which(CHOS_dend_table$cor < 1),]

n = unique(CHOS_dend_table[,1])
CHOS_dend_table2 <- CHOS_dend_table %>% 
  select(x)

for (i in 2:length(n)){
  temp_compound2 = CHOS_dend_table[CHOS_dend_table$y == n[i],]
  temp_compound2 = temp_compound2[,c("cor","x")]
  temp_compound2 = data.frame(temp_compound2)
  names(temp_compound2) = c(n[i],"x")
  
  timeseries3 = right_join(temp_compound2,CHOS_dend_table2,"x")
}

timeseries3 <- time %>% 
  left_join(timeseries2, by = "unixtime")


m_dist<-dist(log_CHOS_t,diag = FALSE )
m_hclust<-hclust(m_dist, method= "complete")
plot(m_hclust)

CHOS_PCA_group1_source<-polarPlot(avgmetwinter3, pollutant = c("C13H20O8SO4 - RT1.38"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group1_source
CHOS_PCA_group2_source<-polarPlot(avgmetwinter3, pollutant = c("C9H12O5S - RT2.01"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group2_source
CHOS_PCA_group3_source<-polarPlot(avgmetwinter3, pollutant = c("C6H12OSO4 - RT1.36"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group3_source
CHOS_PCA_group4_source<-polarPlot(avgmetwinter3, pollutant = c("C11H12O7S - RT1.41"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group4_source
CHOS_PCA_group5_source<-polarPlot(avgmetwinter3, pollutant = c("C4H10SO4 - RT1.70"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group5_source
CHOS_PCA_group6_source<-polarPlot(avgmetwinter3, pollutant = c("C5H11O3SO4 - RT5.19"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group6_source
CHOS_PCA_group7_source<-polarPlot(avgmetwinter3, pollutant = c("C8H18SO4 - RT7.43"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group7_source

timeseries_CHOS4group1 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C11H12O7S - RT3.09',
         'C8H12O6SO4 - RT0.77',
         'C10H16O4SO4 @ 1.28',
         'C10H10O7S - RT1.15',
         'C8H14SO10 - RT0.79',
         'C5H8SO8 - RT0.77',
         'C6H10O5SO4 - RT0.71',
         'C5H6O3SO4 - RT0.71',
         'C10H10O6S - RT1.73',
         'C10H10O6S - RT2.06',
         'C6H8SO6 - RT0.79',
         'C3H6SO6 - RT0.78',
         'Lactic acid sulfate',
         'C6H10O3SO4 - RT0.87',
         'C10H18SO5 - RT5.69',
         'C8H14O3SO4 - RT0.89',
         'C10H16O6SO4 @ 0.77',
         'C9H16O5SO4 - RT0.73',
         'C8H12O5SO4 - RT0.77',
         'C6H10O2SO4 - RT1.16')

timeseries_CHOS4group2 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C9H10O5S - RT3.24',
         'C7H14SO5 - RT2.81',
         'C9H16O4SO4 - RT0.83',
         'C3H6OSO4 - RT0.71',
         'C5H12SO7 - RT0.75',
         'C5H12O3SO4 - RT0.70',
         'C5H10SO6 - RT0.80',
         'C6H10O4SO4 - RT0.77',
         'C5H10SO7 - RT0.78',
         'C4H6SO6 - RT0.76',
         'C5H6OSO4 - RT0.83',
         'C4H8SO7 - RT0.75',
         'C5H6O2SO4 - RT0.73',
         'C5H12OSO4 - RT0.98',
         'C5H10SO5 - RT1.43',
         'C10H10O6S - RT0.84',
         'C9H12O5S - RT3.86',
         'C7H8O4S - RT1.07',
         'C10H12O7S - RT0.98',
         'C9H8O5S - RT2.27',
         'C9H12O5S - RT6.29',
         'C11H20OSO4 - RT5.21',
         'C7H6O4S - RT1.85',
         'C6H14SO4 - RT6.10',
         'C5H12SO4 - RT4.32',
         'C7H8O4S - RT3.22',
         'C14H24OSO4 - RT8.41',
         'C14H27OHSO4 - RT10.09',
         'C13H26SO6 - RT8.83')

timeseries_CHOS4group3 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C12H26SO4 - RT11.50',
         'C9H20OSO4 - RT6.39',
         'C10H22OSO4 - RT6.64',
         'C8H14O2SO4 - RT1.24',
         'C9H20SO4 - RT9.05',
         'C10H10O7S - RT3.18')

timeseries_CHOS4group4 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C9H16SO6 - RT3.61',
         'C9H16O3SO4 - RT1.14',
         'C8H14SO10 - RT3.07',
         'C8H14O2SO4 - RT2.12')

timeseries_CHOS4group5 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C7H16SO4 - RT7.48',
         'C13H28OSO4 - RT9.59',
         'C8H18OSO4 - RT5.05',
         'C16H32O5SO4 - RT10.51')

timeseries_CHOS4group6 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C9H18SO6 - RT5.88',
         'C7H14SO5 - RT3.29',
         'C10H20OSO4 - RT6.71')

timeseries_CHOS4group7 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C4H10SO4 - RT1.70',
         'C9H8O5S - RT4.82',
         'C7H16SO4 - RT6.89',
         'C6H14SO4 - RT5.62')

timeseries_CHOS4group8 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C8H14O3SO4 - RT1.49',
         'C10H18SO6 - RT4.22',
         'C13H26OSO4 - RT9.34',
         'C10H17NSO9 - RT4.24',
         'Octyl Sulfate',
         'C5H6O4SO4 - RT0.71')

timeseries_CHOS4group9 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C5H10SO5 - RT0.94')

timeseries_CHOS4group10 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'Ethyl sulfate',
         'C8H8O4S - RT1.37',
         'Dodecyl Sulfate',
         'C12H20O6SO4 - RT9.82')

timeseries_CHOS4group11 <- timeseries_CHOS4 %>% 
  select(unixtime,
         'C12H26OSO4 - RT8.81',
         'C11H22OSO4 - RT7.37')

timeseries_CHOS4group1$sum <- rowSums( timeseries_CHOS4group1[2:21] , na.rm=TRUE)
timeseries_CHOS4group1$unixtime<-as.POSIXct(timeseries_CHOS4group1$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group1 <- timeseries_CHOS4group1 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group1, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Summed group 1") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")


timeseries_CHOS4group2$sum <- rowSums( timeseries_CHOS4group2[2:30] , na.rm=TRUE)
timeseries_CHOS4group2$unixtime<-as.POSIXct(timeseries_CHOS4group2$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group2 <- timeseries_CHOS4group2 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group2, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Summed group 2") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group3$sum <- rowSums( timeseries_CHOS4group3[2:7] , na.rm=TRUE)
timeseries_CHOS4group3$unixtime<-as.POSIXct(timeseries_CHOS4group3$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group3 <- timeseries_CHOS4group3 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group3, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Summed group 3") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group4$sum <- rowSums( timeseries_CHOS4group4[2:5] , na.rm=TRUE)
timeseries_CHOS4group4$unixtime<-as.POSIXct(timeseries_CHOS4group4$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group4 <- timeseries_CHOS4group4 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group4, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Summed group 4") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group5$sum <- rowSums( timeseries_CHOS4group5[2:5] , na.rm=TRUE)
timeseries_CHOS4group5$unixtime<-as.POSIXct(timeseries_CHOS4group5$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group5 <- timeseries_CHOS4group5 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group5, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Summed group 5") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group6$sum <- rowSums( timeseries_CHOS4group6[2:4] , na.rm=TRUE)
timeseries_CHOS4group6$unixtime<-as.POSIXct(timeseries_CHOS4group6$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group6 <- timeseries_CHOS4group6 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group6, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Summed group 6") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group7$sum <- rowSums( timeseries_CHOS4group7[2:5] , na.rm=TRUE)
timeseries_CHOS4group7$unixtime<-as.POSIXct(timeseries_CHOS4group7$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group7 <- timeseries_CHOS4group7 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group7, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank(),
                                                                                       panel.background = element_blank()) +ggtitle("Summed group 7") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group8$sum <- rowSums( timeseries_CHOS4group8[2:7] , na.rm=TRUE)
timeseries_CHOS4group8$unixtime<-as.POSIXct(timeseries_CHOS4group8$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group8 <- timeseries_CHOS4group8 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group8, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                      panel.grid.major = element_blank(),
                                                                                      panel.grid.minor = element_blank(),
                                                                                      panel.background = element_blank()) +ggtitle("Summed group 8") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group9$sum <- rowSums( timeseries_CHOS4group9[2:2] , na.rm=TRUE)
timeseries_CHOS4group9$unixtime<-as.POSIXct(timeseries_CHOS4group9$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group9 <- timeseries_CHOS4group9 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group9, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                      panel.grid.major = element_blank(),
                                                                                      panel.grid.minor = element_blank(),
                                                                                      panel.background = element_blank()) +ggtitle("Summed group 9") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group10$sum <- rowSums( timeseries_CHOS4group10[2:5] , na.rm=TRUE)
timeseries_CHOS4group10$unixtime<-as.POSIXct(timeseries_CHOS4group10$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group10 <- timeseries_CHOS4group10 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group10, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                      panel.grid.major = element_blank(),
                                                                                      panel.grid.minor = element_blank(),
                                                                                      panel.background = element_blank()) +ggtitle("Summed group 10") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

timeseries_CHOS4group11$sum <- rowSums( timeseries_CHOS4group11[2:3] , na.rm=TRUE)
timeseries_CHOS4group11$unixtime<-as.POSIXct(timeseries_CHOS4group11$unixtime, origin="1970-01-01 00:00:00")
timeseries_CHOS4group11 <- timeseries_CHOS4group11 %>% 
  left_join(avgmetwinter2, by = "unixtime")
ggplot(timeseries_CHOS4group11, aes(unixtime, sum)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                      panel.grid.major = element_blank(),
                                                                                      panel.grid.minor = element_blank(),
                                                                                      panel.background = element_blank()) +ggtitle("Summed group 11") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")


CHOS_PCA_group1_source<-polarPlot(timeseries_CHOS4group1, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group1_source
CHOS_PCA_group2_source<-polarPlot(timeseries_CHOS4group2, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group2_source
CHOS_PCA_group3_source<-polarPlot(timeseries_CHOS4group3, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group3_source
CHOS_PCA_group4_source<-polarPlot(timeseries_CHOS4group4, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group4_source
CHOS_PCA_group5_source<-polarPlot(timeseries_CHOS4group5, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group5_source
CHOS_PCA_group6_source<-polarPlot(timeseries_CHOS4group6, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group6_source
CHOS_PCA_group7_source<-polarPlot(timeseries_CHOS4group7, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group7_source
CHOS_PCA_group8_source<-polarPlot(timeseries_CHOS4group8, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group8_source
CHOS_PCA_group9_source<-polarPlot(timeseries_CHOS4group9, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group9_source
CHOS_PCA_group10_source<-polarPlot(timeseries_CHOS4group10, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group10_source
CHOS_PCA_group11_source<-polarPlot(timeseries_CHOS4group11, pollutant = c("sum"), k = 25, key.header = "Peak area (a.u. / m3)")
CHOS_PCA_group11_source


filepath = paste(getwd(),"/CHOS_cluster5_Rcorr.pdf",sep = "")
pdf(width=40, height=40, file=filepath)
corPlot(timeseries_CHOS4group5,
                     pollutants = names(timeseries_CHOS4group5)[c(2:141)],
                     cols = "jet",
                     dendrogram = TRUE,
                     main = "Correlations of identified CHOS, cluster 5",
                     order="alphabet",
                     cluster=TRUE,
                     text.col = NA)
dev.off()

timeseries_CHOS5_7 <- timeseries_CHOS4group7 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_7)[2] <- "G7"

timeseries_CHOS5_9 <- timeseries_CHOS4group9 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_9)[2] <- "G9"

timeseries_CHOS5_10 <- timeseries_CHOS4group10 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_10)[2] <- "G10"

timeseries_CHOS5_11 <- timeseries_CHOS4group11 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_11)[2] <- "G11"

timeseries_CHOS5_1 <- timeseries_CHOS4group1 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_1)[2] <- "G1"

timeseries_CHOS5_2 <- timeseries_CHOS4group2 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_2)[2] <- "G2"

timeseries_CHOS5_3 <- timeseries_CHOS4group3 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_3)[2] <- "G3"

timeseries_CHOS5_4 <- timeseries_CHOS4group4 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_4)[2] <- "G4"

timeseries_CHOS5_6 <- timeseries_CHOS4group6 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_6)[2] <- "G6"

timeseries_CHOS5_5 <- timeseries_CHOS4group5 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_5)[2] <- "G5"

timeseries_CHOS5_8 <- timeseries_CHOS4group8 %>% 
  select(unixtime,
         sum)

colnames(timeseries_CHOS5_8)[2] <- "G8"

timeseries_CHOS5 <- timeseries_CHOS5_7

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_9, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_10, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_11, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_1, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_2, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_3, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_4, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_6, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_5, by = "unixtime")

timeseries_CHOS5 <- timeseries_CHOS5 %>% 
  left_join(timeseries_CHOS5_8, by = "unixtime")

timeseries_CHOS5$sum <- rowSums( timeseries_CHOS5[2:11] , na.rm=TRUE)

timeseries_CHOS5$prop1 <- (timeseries_CHOS5$G1/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop2 <- (timeseries_CHOS5$G2/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop3 <- (timeseries_CHOS5$G3/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop4 <- (timeseries_CHOS5$G4/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop5 <- (timeseries_CHOS5$G5/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop6 <- (timeseries_CHOS5$G6/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop7 <- (timeseries_CHOS5$G7/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop8 <- (timeseries_CHOS5$G8/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop9 <- (timeseries_CHOS5$G9/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop10 <- (timeseries_CHOS5$G10/timeseries_CHOS5$sum)*100
timeseries_CHOS5$prop11 <- (timeseries_CHOS5$G11/timeseries_CHOS5$sum)*100

write.csv(timeseries_CHOS5, file = "CHOS3.csv")

timeseries_CHOS5_trend <- read.csv("C:/Users/William/Documents/WACL stuff/CHOS4.csv", na.strings=0)
colnames(timeseries_CHOS5_trend)[4] <- "Cluster"
timeseries_CHOS5_trend$unixtime <- dmy_hm(timeseries_CHOS5_trend$unixtime, tz = "UTC")
ggplot(timeseries_CHOS5_trend, aes(x=unixtime,y=G7,group=Cluster,fill=Cluster)) + theme(panel.background = element_blank(), axis.line = element_line(colour = "Black"))+ geom_area() +ggtitle("Contribution to total CHOS from cluster groups") + ylab("Peak area a.u.") + xlab("Date")+scale_fill_discrete(breaks = c("Group 1", "Group 2","Group 3", "Group 4","Group 5", "Group 6","Group 7", "Group 8","Group 9", "Group 10", "Group 11") )


meansvec <- colMeans(timeseries_CHOS4group1[ , c(2:21)] ) 
slices <- c(meansvec) 
lbls <- c('C11H12O7S - RT3.09',
          'C8H12O6SO4 - RT0.77',
          'C10H16O4SO4 @ 1.28',
          'C10H10O7S - RT1.15',
          'C8H14SO10 - RT0.79',
          'C5H8SO8 - RT0.77',
          'C6H10O5SO4 - RT0.71',
          'C5H6O3SO4 - RT0.71',
          'C10H10O6S - RT1.73',
          'C10H10O6S - RT2.06',
          'C6H8SO6 - RT0.79',
          'C3H6SO6 - RT0.78',
          'Lactic acid sulfate',
          'C6H10O3SO4 - RT0.87',
          'C10H18SO5 - RT5.69',
          'C8H14O3SO4 - RT0.89',
          'C10H16O6SO4 @ 0.77',
          'C9H16O5SO4 - RT0.73',
          'C8H12O5SO4 - RT0.77',
          'C6H10O2SO4 - RT1.16')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 1")
timeseries_CHOS4group1_1<-timeseries_CHOS4group1
timeseries_CHOS4group1_1$'C11H12O7S - RT3.09prop' <- (timeseries_CHOS4group1_1$'C11H12O7S - RT3.09'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C8H12O6SO4 - RT0.77prop' <- (timeseries_CHOS4group1_1$'C8H12O6SO4 - RT0.77'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C10H16O4SO4 @ 1.28prop' <- (timeseries_CHOS4group1_1$'C10H16O4SO4 @ 1.28'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C10H10O7S - RT1.15prop' <- (timeseries_CHOS4group1_1$'C10H10O7S - RT1.15'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C8H14SO10 - RT0.79prop' <- (timeseries_CHOS4group1_1$'C8H14SO10 - RT0.79'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C5H8SO8 - RT0.77prop' <- (timeseries_CHOS4group1_1$'C5H8SO8 - RT0.77'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C6H10O5SO4 - RT0.71prop' <- (timeseries_CHOS4group1_1$'C6H10O5SO4 - RT0.71'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C5H6O3SO4 - RT0.71' <- (timeseries_CHOS4group1_1$'C5H6O3SO4 - RT0.71'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C10H10O6S - RT1.73prop' <- (timeseries_CHOS4group1_1$'C10H10O6S - RT1.73'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C10H10O6S - RT2.06prop' <- (timeseries_CHOS4group1_1$'C10H10O6S - RT2.06'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C6H8SO6 - RT0.79prop' <- (timeseries_CHOS4group1_1$'C6H8SO6 - RT0.79'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C3H6SO6 - RT0.78prop' <- (timeseries_CHOS4group1_1$'C3H6SO6 - RT0.78'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'Lactic acid sulfateprop' <- (timeseries_CHOS4group1_1$'Lactic acid sulfate'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C6H10O3SO4 - RT0.87prop' <- (timeseries_CHOS4group1_1$'C6H10O3SO4 - RT0.87'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C10H18SO5 - RT5.69prop' <- (timeseries_CHOS4group1_1$'C10H18SO5 - RT5.69'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C8H14O3SO4 - RT0.89prop' <- (timeseries_CHOS4group1_1$'C8H14O3SO4 - RT0.89'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C10H16O6SO4 @ 0.77prop' <- (timeseries_CHOS4group1_1$'C10H16O6SO4 @ 0.77'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C9H16O5SO4 - RT0.73prop' <- (timeseries_CHOS4group1_1$'C9H16O5SO4 - RT0.73'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C8H12O5SO4 - RT0.77prop' <- (timeseries_CHOS4group1_1$'C8H12O5SO4 - RT0.77'/timeseries_CHOS4group1_1$sum)*100
timeseries_CHOS4group1_1$'C6H10O2SO4 - RT1.16prop' <- (timeseries_CHOS4group1_1$'C6H10O2SO4 - RT1.16'/timeseries_CHOS4group1_1$sum)*100


meansvec <- colMeans(timeseries_CHOS4group2[ , c(2:30)] ) 
slices <- c(meansvec) 
lbls <- c('C9H10O5S - RT3.24',
          'C7H14SO5 - RT2.81',
          'C9H16O4SO4 - RT0.83',
          'C3H6OSO4 - RT0.71',
          'C5H12SO7 - RT0.75',
          'C5H12O3SO4 - RT0.70',
          'C5H10SO6 - RT0.80',
          'C6H10O4SO4 - RT0.77',
          'C5H10SO7 - RT0.78',
          'C4H6SO6 - RT0.76',
          'C5H6OSO4 - RT0.83',
          'C4H8SO7 - RT0.75',
          'C5H6O2SO4 - RT0.73',
          'C5H12OSO4 - RT0.98',
          'C5H10SO5 - RT1.43',
          'C10H10O6S - RT0.84',
          'C9H12O5S - RT3.86',
          'C7H8O4S - RT1.07',
          'C10H12O7S - RT0.98',
          'C9H8O5S - RT2.27',
          'C9H12O5S - RT6.29',
          'C11H20OSO4 - RT5.21',
          'C7H6O4S - RT1.85',
          'C6H14SO4 - RT6.10',
          'C5H12SO4 - RT4.32',
          'C7H8O4S - RT3.22',
          'C14H24OSO4 - RT8.41',
          'C14H27OHSO4 - RT10.09',
          'C13H26SO6 - RT8.83')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 2")
timeseries_CHOS4group2_1<-timeseries_CHOS4group2
timeseries_CHOS4group2_1$'C9H10O5S - RT3.24prop' <- (timeseries_CHOS4group2_1$'C9H10O5S - RT3.24'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C7H14SO5 - RT2.81prop' <- (timeseries_CHOS4group2_1$'C7H14SO5 - RT2.81'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C9H16O4SO4 - RT0.83prop' <- (timeseries_CHOS4group2_1$'C9H16O4SO4 - RT0.83'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C3H6OSO4 - RT0.71prop' <- (timeseries_CHOS4group2_1$'C3H6OSO4 - RT0.71'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H12SO7 - RT0.75prop' <- (timeseries_CHOS4group2_1$'C5H12SO7 - RT0.75'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H12O3SO4 - RT0.70prop' <- (timeseries_CHOS4group2_1$'C5H12O3SO4 - RT0.70'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H10SO6 - RT0.80prop' <- (timeseries_CHOS4group2_1$'C5H10SO6 - RT0.80'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C6H10O4SO4 - RT0.77prop' <- (timeseries_CHOS4group2_1$'C6H10O4SO4 - RT0.77'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H10SO7 - RT0.78prop' <- (timeseries_CHOS4group2_1$'C5H10SO7 - RT0.78'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C4H6SO6 - RT0.76prop' <- (timeseries_CHOS4group2_1$'C4H6SO6 - RT0.76'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H6OSO4 - RT0.83prop' <- (timeseries_CHOS4group2_1$'C5H6OSO4 - RT0.83'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C4H8SO7 - RT0.75prop' <- (timeseries_CHOS4group2_1$'C4H8SO7 - RT0.75'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H6O2SO4 - RT0.73prop' <- (timeseries_CHOS4group2_1$'C5H6O2SO4 - RT0.73'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H12OSO4 - RT0.98prop' <- (timeseries_CHOS4group2_1$'C5H12OSO4 - RT0.98'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H10SO5 - RT1.43prop' <- (timeseries_CHOS4group2_1$'C5H10SO5 - RT1.43'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C10H10O6S - RT0.84prop' <- (timeseries_CHOS4group2_1$'C10H10O6S - RT0.84'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C9H12O5S - RT3.86prop' <- (timeseries_CHOS4group2_1$'C9H12O5S - RT3.86'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C7H8O4S - RT1.07prop' <- (timeseries_CHOS4group2_1$'C7H8O4S - RT1.07'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C10H12O7S - RT0.98prop' <- (timeseries_CHOS4group2_1$'C10H12O7S - RT0.98'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C9H8O5S - RT2.27prop' <- (timeseries_CHOS4group2_1$'C9H8O5S - RT2.27'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C9H12O5S - RT6.29prop' <- (timeseries_CHOS4group2_1$'C9H12O5S - RT6.29'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C11H20OSO4 - RT5.21prop' <- (timeseries_CHOS4group2_1$'C11H20OSO4 - RT5.21'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C7H6O4S - RT1.85prop' <- (timeseries_CHOS4group2_1$'C7H6O4S - RT1.85'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C6H14SO4 - RT6.10prop' <- (timeseries_CHOS4group2_1$'C6H14SO4 - RT6.10'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C5H12SO4 - RT4.32prop' <- (timeseries_CHOS4group2_1$'C5H12SO4 - RT4.32'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C7H8O4S - RT3.22prop' <- (timeseries_CHOS4group2_1$'C7H8O4S - RT3.22'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C14H24OSO4 - RT8.41prop' <- (timeseries_CHOS4group2_1$'C14H24OSO4 - RT8.41'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C14H27OHSO4 - RT10.09prop' <- (timeseries_CHOS4group2_1$'C14H27OHSO4 - RT10.09'/timeseries_CHOS4group2_1$sum)*100
timeseries_CHOS4group2_1$'C13H26SO6 - RT8.83prop' <- (timeseries_CHOS4group2_1$'C13H26SO6 - RT8.83'/timeseries_CHOS4group2_1$sum)*100

meansvec <- colMeans(timeseries_CHOS4group3[ , c(2:7)] ) 
slices <- c(meansvec) 
lbls <- c('C12H26SO4 - RT11.50',
          'C9H20OSO4 - RT6.39',
          'C10H22OSO4 - RT6.64',
          'C8H14O2SO4 - RT1.24',
          'C9H20SO4 - RT9.05',
          'C10H10O7S - RT3.18')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 3")
timeseries_CHOS4group3_1<-timeseries_CHOS4group3
timeseries_CHOS4group3_1$'C12H26SO4 - RT11.50prop' <- (timeseries_CHOS4group3_1$'C12H26SO4 - RT11.50'/timeseries_CHOS4group3_1$sum)*100
timeseries_CHOS4group3_1$'C9H20OSO4 - RT6.39prop' <- (timeseries_CHOS4group3_1$'C9H20OSO4 - RT6.39'/timeseries_CHOS4group3_1$sum)*100
timeseries_CHOS4group3_1$'C10H22OSO4 - RT6.64prop' <- (timeseries_CHOS4group3_1$'C10H22OSO4 - RT6.64'/timeseries_CHOS4group3_1$sum)*100
timeseries_CHOS4group3_1$'C8H14O2SO4 - RT1.24prop' <- (timeseries_CHOS4group3_1$'C8H14O2SO4 - RT1.24'/timeseries_CHOS4group3_1$sum)*100
timeseries_CHOS4group3_1$'C9H20SO4 - RT9.05prop' <- (timeseries_CHOS4group3_1$'C9H20SO4 - RT9.05'/timeseries_CHOS4group3_1$sum)*100
timeseries_CHOS4group3_1$'C10H10O7S - RT3.18prop' <- (timeseries_CHOS4group3_1$'C10H10O7S - RT3.18'/timeseries_CHOS4group3_1$sum)*100

meansvec <- colMeans(timeseries_CHOS4group4[ , c(2:5)] ) 
slices <- c(meansvec) 
lbls <- c('C9H16SO6 - RT3.61',
          'C9H16O3SO4 - RT1.14',
          'C8H14SO10 - RT3.07',
          'C8H14O2SO4 - RT2.12')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 4")
timeseries_CHOS4group4_1<-timeseries_CHOS4group4
timeseries_CHOS4group4_1$'C9H16SO6 - RT3.61prop' <- (timeseries_CHOS4group4_1$'C9H16SO6 - RT3.61'/timeseries_CHOS4group4_1$sum)*100
timeseries_CHOS4group4_1$'C9H16O3SO4 - RT1.14prop' <- (timeseries_CHOS4group4_1$'C9H16O3SO4 - RT1.14'/timeseries_CHOS4group4_1$sum)*100
timeseries_CHOS4group4_1$'C8H14SO10 - RT3.07prop' <- (timeseries_CHOS4group4_1$'C8H14SO10 - RT3.07'/timeseries_CHOS4group4_1$sum)*100
timeseries_CHOS4group4_1$'C8H14O2SO4 - RT2.12prop' <- (timeseries_CHOS4group4_1$'C8H14O2SO4 - RT2.12'/timeseries_CHOS4group4_1$sum)*100

meansvec <- colMeans(timeseries_CHOS4group5[ , c(2:5)] ) 
slices <- c(meansvec) 
lbls <- c('C7H16SO4 - RT7.48',
          'C13H28OSO4 - RT9.59',
          'C8H18OSO4 - RT5.05',
          'C16H32O5SO4 - RT10.51')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 5")
timeseries_CHOS4group5_1<-timeseries_CHOS4group5
timeseries_CHOS4group5_1$'C7H16SO4 - RT7.48prop' <- (timeseries_CHOS4group5_1$'C7H16SO4 - RT7.48'/timeseries_CHOS4group5_1$sum)*100
timeseries_CHOS4group5_1$'C13H28OSO4 - RT9.59prop' <- (timeseries_CHOS4group5_1$'C13H28OSO4 - RT9.59'/timeseries_CHOS4group5_1$sum)*100
timeseries_CHOS4group5_1$'C8H18OSO4 - RT5.05prop' <- (timeseries_CHOS4group5_1$'C8H18OSO4 - RT5.05'/timeseries_CHOS4group5_1$sum)*100
timeseries_CHOS4group5_1$'C16H32O5SO4 - RT10.51prop' <- (timeseries_CHOS4group5_1$'C16H32O5SO4 - RT10.51'/timeseries_CHOS4group5_1$sum)*100

meansvec <- colMeans(timeseries_CHOS4group6[ , c(2:4)] ) 
slices <- c(meansvec) 
lbls <- c('C9H18SO6 - RT5.88',
          'C7H14SO5 - RT3.29',
          'C10H20OSO4 - RT6.71')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 6")
timeseries_CHOS4group6_1<-timeseries_CHOS4group6
timeseries_CHOS4group6_1$'C9H18SO6 - RT5.88prop' <- (timeseries_CHOS4group6_1$'C9H18SO6 - RT5.88'/timeseries_CHOS4group6_1$sum)*100
timeseries_CHOS4group6_1$'C7H14SO5 - RT3.29prop' <- (timeseries_CHOS4group6_1$'C7H14SO5 - RT3.29'/timeseries_CHOS4group6_1$sum)*100
timeseries_CHOS4group6_1$'C10H20OSO4 - RT6.71prop' <- (timeseries_CHOS4group6_1$'C10H20OSO4 - RT6.71'/timeseries_CHOS4group6_1$sum)*100

meansvec <- colMeans(timeseries_CHOS4group7[ , c(2:5)] ) 
slices <- c(meansvec) 
lbls <- c('C4H10SO4 - RT1.70',
          'C9H8O5S - RT4.82',
          'C7H16SO4 - RT6.89',
          'C6H14SO4 - RT5.62')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 7")
timeseries_CHOS4group7_1<-timeseries_CHOS4group7
timeseries_CHOS4group7_1$'C4H10SO4 - RT1.70prop' <- (timeseries_CHOS4group7_1$'C4H10SO4 - RT1.70'/timeseries_CHOS4group7_1$sum)*100
timeseries_CHOS4group7_1$'C9H8O5S - RT4.82prop' <- (timeseries_CHOS4group7_1$'C9H8O5S - RT4.82'/timeseries_CHOS4group7_1$sum)*100
timeseries_CHOS4group7_1$'C7H16SO4 - RT6.89prop' <- (timeseries_CHOS4group7_1$'C7H16SO4 - RT6.89'/timeseries_CHOS4group7_1$sum)*100
timeseries_CHOS4group7_1$'C6H14SO4 - RT5.62prop' <- (timeseries_CHOS4group7_1$'C6H14SO4 - RT5.62'/timeseries_CHOS4group7_1$sum)*100

meansvec <- colMeans(timeseries_CHOS4group8[ , c(2:7)] ) 
slices <- c(meansvec) 
lbls <- c('C8H14O3SO4 - RT1.49',
          'C10H18SO6 - RT4.22',
          'C13H26OSO4 - RT9.34',
          'C10H17NSO9 - RT4.24',
          'Octyl Sulfate',
          'C5H6O4SO4 - RT0.71')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 8")
timeseries_CHOS4group8_1<-timeseries_CHOS4group8
timeseries_CHOS4group8_1$'C8H14O3SO4 - RT1.49prop' <- (timeseries_CHOS4group8_1$'C8H14O3SO4 - RT1.49'/timeseries_CHOS4group8_1$sum)*100
timeseries_CHOS4group8_1$'C10H18SO6 - RT4.22prop' <- (timeseries_CHOS4group8_1$'C10H18SO6 - RT4.22'/timeseries_CHOS4group8_1$sum)*100
timeseries_CHOS4group8_1$'C13H26OSO4 - RT9.34prop' <- (timeseries_CHOS4group8_1$'C13H26OSO4 - RT9.34'/timeseries_CHOS4group8_1$sum)*100
timeseries_CHOS4group8_1$'C10H17NSO9 - RT4.24prop' <- (timeseries_CHOS4group8_1$'C10H17NSO9 - RT4.24'/timeseries_CHOS4group8_1$sum)*100
timeseries_CHOS4group8_1$'Octyl Sulfateprop' <- (timeseries_CHOS4group8_1$'Octyl Sulfate'/timeseries_CHOS4group8_1$sum)*100
timeseries_CHOS4group8_1$'C5H6O4SO4 - RT0.71prop' <- (timeseries_CHOS4group8_1$'C5H6O4SO4 - RT0.71'/timeseries_CHOS4group8_1$sum)*100

meansvec <- 100 
slices <- c(meansvec) 
lbls <- c('C5H10SO5 - RT0.94')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 9")

meansvec <- colMeans(timeseries_CHOS4group10[ , c(2:5)] ) 
slices <- c(meansvec) 
lbls <- c('Ethyl sulfate',
          'C8H8O4S - RT1.37',
          'Dodecyl Sulfate',
          'C12H20O6SO4 - RT9.82')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 10")
timeseries_CHOS4group10_1<-timeseries_CHOS4group10
timeseries_CHOS4group10_1$'Ethyl sulfateprop' <- (timeseries_CHOS4group10_1$'Ethyl sulfate'/timeseries_CHOS4group10_1$sum)*100
timeseries_CHOS4group10_1$'C8H8O4S - RT1.37prop' <- (timeseries_CHOS4group10_1$'C8H8O4S - RT1.37'/timeseries_CHOS4group10_1$sum)*100
timeseries_CHOS4group10_1$'Dodecyl Sulfateprop' <- (timeseries_CHOS4group10_1$'Dodecyl Sulfate'/timeseries_CHOS4group10_1$sum)*100
timeseries_CHOS4group10_1$'C12H20O6SO4 - RT9.82prop' <- (timeseries_CHOS4group10_1$'C12H20O6SO4 - RT9.82'/timeseries_CHOS4group10_1$sum)*100

meansvec <- colMeans(timeseries_CHOS4group11[ , c(2:3)] ) 
slices <- c(meansvec) 
lbls <- c('C12H26OSO4 - RT8.81',
          'C11H22OSO4 - RT7.37')
pie(slices,labels=lbls,explode=0.1,
    main="Proportion of cluster 11")
timeseries_CHOS4group11_1<-timeseries_CHOS4group11
timeseries_CHOS4group11_1$'C12H26OSO4 - RT8.81prop' <- (timeseries_CHOS4group11_1$'C12H26OSO4 - RT8.81'/timeseries_CHOS4group11_1$sum)*100
timeseries_CHOS4group11_1$'C11H22OSO4 - RT7.37prop' <- (timeseries_CHOS4group11_1$'C11H22OSO4 - RT7.37'/timeseries_CHOS4group11_1$sum)*100

write.csv(timeseries_CHOS4group1_1, file = "CHOSgroup1.csv")
write.csv(timeseries_CHOS4group2_1, file = "CHOSgroup2.csv")
write.csv(timeseries_CHOS4group3_1, file = "CHOSgroup3.csv")
write.csv(timeseries_CHOS4group4_1, file = "CHOSgroup4.csv")
write.csv(timeseries_CHOS4group5_1, file = "CHOSgroup5.csv")
write.csv(timeseries_CHOS4group6_1, file = "CHOSgroup6.csv")
write.csv(timeseries_CHOS4group7_1, file = "CHOSgroup7.csv")
write.csv(timeseries_CHOS4group8_1, file = "CHOSgroup8.csv")
write.csv(timeseries_CHOS4group9_1, file = "CHOSgroup9.csv")
write.csv(timeseries_CHOS4group10_1, file = "CHOSgroup10.csv")
write.csv(timeseries_CHOS4group11_1, file = "CHOSgroup11.csv")

HetGroup2<-HetGroup2 %>% 
  left_join(AMS, by = "datetime")

HetGroup2<-HetGroup2 %>% 
  left_join(PM25, by = "datetime")

HetGroup2_CHO <- HetGroup2[which(HetGroup2$Group=="CHO" ),]
HetGroup2_CHOS <- HetGroup2[which(HetGroup2$Group=="CHOS" ),]

HetGroup2_CHO <- HetGroup2_CHO[ which(HetGroup2_CHO$datetime > "2016-11-19 00:00:00"),]
HetGroup2_CHO <- HetGroup2_CHO[ which(HetGroup2_CHO$datetime < "2016-11-28 00:00:00"),]
HetGroup2_CHO <- HetGroup2_CHO[ which(HetGroup2_CHO$'PM2.5' > 20),]

HetGroup2_CHOS <- HetGroup2_CHOS[ which(HetGroup2_CHOS$datetime > "2016-11-19 00:00:00"),]
HetGroup2_CHOS <- HetGroup2_CHOS[ which(HetGroup2_CHOS$datetime < "2016-11-28 00:00:00"),]
HetGroup2_CHOS <- HetGroup2_CHOS[ which(HetGroup2_CHOS$'SO4' > 3),]

plot(HetGroup2_CHO$Quant, HetGroup2_CHO$'PM2.5')
plot(HetGroup2_CHOS$Quant, HetGroup2_CHOS$'SO4')

lm.CHO<-lm(formula= HetGroup2_CHO$Quant~HetGroup2_CHO$'PM2.5')
summary(lm.CHO)

lm.CHOS<-lm(formula= HetGroup2_CHOS$Quant~HetGroup2_CHOS$'SO4')
summary(lm.CHOS)

WinterHetGroup <- HetGroup
WinterHetGroup$Season <- 'Winter'

ggplot(HetGroup, aes(x=Group, y=Quant, color=Group)) +
  geom_boxplot()

boxplot<-rbind(WinterHetGroup, SummerHetGroup)

ggplot(boxplot, aes(x=Group, y=Quant, color=Season)) +
  geom_boxplot() + ylim(0, 4000000)

boxplot_CHO <- boxplot[ which(boxplot$Group == 'CHO'),]
boxplot_CHOS <- boxplot[ which(boxplot$Group == 'CHOS'),]
boxplot_CHON <- boxplot[ which(boxplot$Group == 'CHON'),]
boxplot_CHONS <- boxplot[ which(boxplot$Group == 'CHONS'),]

ggplot(boxplot_CHO, aes(x=Season, y=Quant, color=Season)) +
  geom_boxplot() + ylim(0, 300000) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Peak area a.u.") +stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE) + ggtitle("CHO")
ggplot(boxplot_CHON, aes(x=Season, y=Quant, color=Season)) +
  geom_boxplot() + ylim(0, 4000000) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                            panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Peak area a.u.") +stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE) + ggtitle("CHON")
ggplot(boxplot_CHOS, aes(x=Season, y=Quant, color=Season)) +
  geom_boxplot() + ylim(0, 200000) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                          panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Peak area a.u.") +stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE) + ggtitle("CHOS")
ggplot(boxplot_CHONS, aes(x=Season, y=Quant, color=Season)) +
  geom_boxplot() + ylim(0, 100000) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Peak area a.u.") +stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE) + ggtitle("CHONS")




ggplot(HetGroupCHONS, aes(datetime, Quant)) + geom_line() + theme_bw() + theme(text = element_text(size=20), axis.line = element_line(colour = "black"),
                                                                                      panel.grid.major = element_blank(),
                                                                                      panel.grid.minor = element_blank(),
                                                                                      panel.background = element_blank()) +ggtitle("CHONS") + guides(size=FALSE) + ylab("Summed peak area") + xlab("Date")

colnames(AMS)[1] <- "unixtime"
colnames(avgmetwinter)[2] <- "unixtime"

avgmetwinter$unixtime <- ymd_hms(avgmetwinter$unixtime, tz = "UTC")
AMS$unixtime <- ymd_hms(AMS$unixtime, tz = "UTC")

timeseries_CHOS4group1_2 <- timeseries_CHOS4group1%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group2_2 <- timeseries_CHOS4group2%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group3_2 <- timeseries_CHOS4group3%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group4_2 <- timeseries_CHOS4group4%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group5_2 <- timeseries_CHOS4group5%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group6_2 <- timeseries_CHOS4group6%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group7_2 <- timeseries_CHOS4group7%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group8_2 <- timeseries_CHOS4group8%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group9_2 <- timeseries_CHOS4group9%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group10_2 <- timeseries_CHOS4group10%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")
timeseries_CHOS4group11_2 <- timeseries_CHOS4group11%>%
  left_join(avgmetwinter, by = "unixtime") %>% 
  left_join(AMS, by = "unixtime")

group1O3<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("O3")
group1RH<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("Relative humidity %")
group1Temp<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("Temperature oC")
group1SO2<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("SO2")
group1NOx<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("NOx")
group1_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("AMS CCOA")
group1_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("AMS COA")
group1_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("AMS BBOA (ug m-3)")
group1_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("AMS LOOOA")
group1_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("AMS MOOOA (ug m-3)")
multiplot(group1O3, group1SO2, group1NOx, group1RH, group1Temp, group1_AMS_BBOA, group1_AMS_CCOA, group1_AMS_COA, group1_AMS_MOOOA, group1_AMS_LOOOA, cols = 5)
lm.G1BBOA<-lm(formula= timeseries_CHOS4group1_2$BBOA~timeseries_CHOS4group1_2$sum)
summary(lm.G1BBOA)
lm.G1MOOOA<-lm(formula= timeseries_CHOS4group1_2$MOOOA~timeseries_CHOS4group1_2$sum)
summary(lm.G1MOOOA)
ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G1SO4<-lm(formula= timeseries_CHOS4group1_2$SO4~timeseries_CHOS4group1_2$sum)
summary(lm.G1SO4)
ggplot(data=subset(timeseries_CHOS4group1_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 1 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G1NO3<-lm(formula= timeseries_CHOS4group1_2$NO3~timeseries_CHOS4group1_2$sum)
summary(lm.G1NO3)

group2O3<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("O3")
group2RH<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("Relative humidity %")
group2Temp<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("Temperature oC")
group2SO2<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("SO2")
group2NOx<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("NOx")
group2_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("AMS CCOA")
group2_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("AMS COA")
group2_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("AMS BBOA")
group2_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("AMS LOOOA")
group2_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("AMS MOOOA")
multiplot(group2O3, group2SO2, group2NOx, group2RH, group2Temp, group2_AMS_BBOA, group2_AMS_CCOA, group2_AMS_COA, group2_AMS_MOOOA, group2_AMS_LOOOA, cols = 5)
lm.G2SO2<-lm(formula= timeseries_CHOS4group2_2$so2_ppbv~timeseries_CHOS4group2_2$sum)
summary(lm.G2SO2)
lm.G2LOOOA<-lm(formula= timeseries_CHOS4group2_2$LOOOA~timeseries_CHOS4group2_2$sum)
summary(lm.G2LOOOA)
lm.G2BBOA<-lm(formula= timeseries_CHOS4group2_2$BBOA~timeseries_CHOS4group2_2$sum)
summary(lm.G2BBOA)
ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G2SO4<-lm(formula= timeseries_CHOS4group2_2$SO4~timeseries_CHOS4group2_2$sum)
summary(lm.G2SO4)
ggplot(data=subset(timeseries_CHOS4group2_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 2 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G2NO3<-lm(formula= timeseries_CHOS4group2_2$NO3~timeseries_CHOS4group2_2$sum)
summary(lm.G2NO3)


group3O3<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("O3")
group3RH<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("Relative humidity %")
group3Temp<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("Temperature oC")
group3SO2<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("SO2")
group3NOx<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("NOx")
group3_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("AMS CCOA")
group3_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("AMS COA")
group3_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("AMS BBOA")
group3_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("AMS LOOOA")
group3_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("AMS MOOOA")
multiplot(group3O3, group3SO2, group3NOx, group3RH, group3Temp, group3_AMS_BBOA, group3_AMS_CCOA, group3_AMS_COA, group3_AMS_MOOOA, group3_AMS_LOOOA, cols = 5)
lm.G3CCOA<-lm(formula= timeseries_CHOS4group3_2$CCOA~timeseries_CHOS4group3_2$sum)
summary(lm.G3CCOA)
ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G3SO4<-lm(formula= timeseries_CHOS4group3_2$SO4~timeseries_CHOS4group3_2$sum)
summary(lm.G3SO4)
ggplot(data=subset(timeseries_CHOS4group3_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 3 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G3NO3<-lm(formula= timeseries_CHOS4group3_2$NO3~timeseries_CHOS4group3_2$sum)
summary(lm.G3NO3)


group4O3<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("O3")
group4RH<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("Relative humidity %")
group4Temp<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("Temperature oC")
group4SO2<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("SO2")
group4NOx<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("NOx")
group4_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("AMS CCOA")
group4_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("AMS COA")
group4_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("AMS BBOA")
group4_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("AMS LOOOA")
group4_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("AMS MOOOA")
multiplot(group4O3, group4SO2, group4NOx, group4RH, group4Temp, group4_AMS_BBOA, group4_AMS_CCOA, group4_AMS_COA, group4_AMS_MOOOA, group4_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G4SO4<-lm(formula= timeseries_CHOS4group4_2$SO4~timeseries_CHOS4group4_2$sum)
summary(lm.G4SO4)
ggplot(data=subset(timeseries_CHOS4group4_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 4 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G4NO3<-lm(formula= timeseries_CHOS4group4_2$NO3~timeseries_CHOS4group4_2$sum)
summary(lm.G4NO3)


group5O3<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("O3")
group5RH<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("Relative humidity %")
group5Temp<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("Temperature oC")
group5SO2<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("SO2")
group5NOx<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("NOx")
group5_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("AMS CCOA")
group5_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("AMS COA")
group5_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("AMS BBOA")
group5_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("AMS LOOOA")
group5_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("AMS MOOOA")
multiplot(group5O3, group5SO2, group5NOx, group5RH, group5Temp, group5_AMS_BBOA, group5_AMS_CCOA, group5_AMS_COA, group5_AMS_MOOOA, group5_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G5SO4<-lm(formula= timeseries_CHOS4group5_2$SO4~timeseries_CHOS4group5_2$sum)
summary(lm.G5SO4)
ggplot(data=subset(timeseries_CHOS4group5_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 5 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G5NO3<-lm(formula= timeseries_CHOS4group5_2$NO3~timeseries_CHOS4group5_2$sum)
summary(lm.G5NO3)


group6O3<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("O3")
group6RH<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("Relative humidity %")
group6Temp<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("Temperature oC")
group6SO2<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("SO2")
group6NOx<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("NOx")
group6_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("AMS CCOA")
group6_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("AMS COA")
group6_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("AMS BBOA")
group6_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("AMS LOOOA")
group6_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("AMS MOOOA")
multiplot(group6O3, group6SO2, group6NOx, group6RH, group6Temp, group6_AMS_BBOA, group6_AMS_CCOA, group6_AMS_COA, group6_AMS_MOOOA, group6_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G6SO4<-lm(formula= timeseries_CHOS4group6_2$SO4~timeseries_CHOS4group6_2$sum)
summary(lm.G6SO4)
ggplot(data=subset(timeseries_CHOS4group6_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 6 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G6NO3<-lm(formula= timeseries_CHOS4group6_2$NO3~timeseries_CHOS4group6_2$sum)
summary(lm.G6NO3)


group7O3<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("O3")
group7RH<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("Relative humidity %")
group7Temp<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("Temperature oC")
group7SO2<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("SO2")
group7NOx<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("NOx")
group7_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("AMS CCOA")
group7_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("AMS COA")
group7_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("AMS BBOA")
group7_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("AMS LOOOA")
group7_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("AMS MOOOA")
multiplot(group7O3, group7SO2, group7NOx, group7RH, group7Temp, group7_AMS_BBOA, group7_AMS_CCOA, group7_AMS_COA, group7_AMS_MOOOA, group7_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G7SO4<-lm(formula= timeseries_CHOS4group7_2$SO4~timeseries_CHOS4group7_2$sum)
summary(lm.G7SO4)
ggplot(data=subset(timeseries_CHOS4group7_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 7 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G7NO3<-lm(formula= timeseries_CHOS4group7_2$NO3~timeseries_CHOS4group7_2$sum)
summary(lm.G7NO3)


group8O3<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("O3")
group8RH<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("Relative humidity %")
group8Temp<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("Temperature oC")
group8SO2<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("SO2")
group8NOx<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("NOx")
group8_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("AMS CCOA")
group8_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("AMS COA")
group8_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("AMS BBOA")
group8_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("AMS LOOOA")
group8_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("AMS MOOOA")
multiplot(group8O3, group8SO2, group8NOx, group8RH, group8Temp, group8_AMS_BBOA, group8_AMS_CCOA, group8_AMS_COA, group8_AMS_MOOOA, group8_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G8SO4<-lm(formula= timeseries_CHOS4group8_2$SO4~timeseries_CHOS4group8_2$sum)
summary(lm.G8SO4)
ggplot(data=subset(timeseries_CHOS4group8_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 8 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G8NO3<-lm(formula= timeseries_CHOS4group8_2$NO3~timeseries_CHOS4group8_2$sum)
summary(lm.G8NO3)


group9O3<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("O3")
group9RH<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("Relative humidity %")
group9Temp<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("Temperature oC")
group9SO2<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("SO2")
group9NOx<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("NOx")
group9_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("AMS CCOA")
group9_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("AMS COA")
group9_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("AMS BBOA")
group9_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("AMS LOOOA")
group9_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("AMS MOOOA")
multiplot(group9O3, group9SO2, group9NOx, group9RH, group9Temp, group9_AMS_BBOA, group9_AMS_CCOA, group9_AMS_COA, group9_AMS_MOOOA, group9_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G9SO4<-lm(formula= timeseries_CHOS4group9_2$SO4~timeseries_CHOS4group9_2$sum)
summary(lm.G9SO4)
ggplot(data=subset(timeseries_CHOS4group9_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 9 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G9NO3<-lm(formula= timeseries_CHOS4group9_2$NO3~timeseries_CHOS4group9_2$sum)
summary(lm.G9NO3)


group10O3<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("O3")
group10RH<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("Relative humidity %")
group10Temp<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("Temperature oC")
group10SO2<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("SO2")
group10NOx<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("NOx")
group10_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("AMS CCOA")
group10_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("AMS COA")
group10_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("AMS BBOA")
group10_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("AMS LOOOA")
group10_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("AMS MOOOA")
multiplot(group10O3, group10SO2, group10NOx, group10RH, group10Temp, group10_AMS_BBOA, group10_AMS_CCOA, group10_AMS_COA, group10_AMS_MOOOA, group10_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G10SO4<-lm(formula= timeseries_CHOS4group10_2$SO4~timeseries_CHOS4group10_2$sum)
summary(lm.G10SO4)
ggplot(data=subset(timeseries_CHOS4group10_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 10 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G10NO3<-lm(formula= timeseries_CHOS4group10_2$NO3~timeseries_CHOS4group10_2$sum)
summary(lm.G10NO3)


group11O3<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(o3_ppbv)), aes(x=sum, y=o3_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("O3")
group11RH<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(rh_8m.y)), aes(x=sum, y=rh_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("Relative humidity %")
group11Temp<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(temp_8m.y)), aes(x=sum, y=temp_8m.y)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("Temperature oC")
group11SO2<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(so2_ppbv)), aes(x=sum, y=so2_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("SO2")
group11NOx<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(nox_ppbv)), aes(x=sum, y=nox_ppbv)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("NOx")
group11_AMS_CCOA<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(CCOA)), aes(x=sum, y=CCOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("AMS CCOA")
group11_AMS_COA<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(COA)), aes(x=sum, y=COA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("AMS COA")
group11_AMS_BBOA<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(BBOA)), aes(x=sum, y=BBOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("AMS BBOA")
group11_AMS_LOOOA<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(LOOOA)), aes(x=sum, y=LOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("AMS LOOOA")
group11_AMS_MOOOA<-ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(MOOOA)), aes(x=sum, y=MOOOA)) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("AMS MOOOA")
multiplot(group11O3, group11SO2, group11NOx, group11RH, group11Temp, group11_AMS_BBOA, group11_AMS_CCOA, group11_AMS_COA, group11_AMS_MOOOA, group11_AMS_LOOOA, cols = 5)
ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(SO4)), aes(x=sum, y=SO4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("AMS SO4 (ug m-3)")
lm.G11SO4<-lm(formula= timeseries_CHOS4group11_2$SO4~timeseries_CHOS4group11_2$sum)
summary(lm.G11SO4)
ggplot(data=subset(timeseries_CHOS4group11_2, !is.na(NO3)), aes(x=sum, y=NO3)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_point() + geom_smooth(method=lm) + xlab("Group 11 Peak area") + ylab("AMS NO3 (ug m-3)")
lm.G11NO3<-lm(formula= timeseries_CHOS4group11_2$NO3~timeseries_CHOS4group11_2$sum)
summary(lm.G11NO3)


box_CHO_w <- boxplot_CHO[ which(boxplot_CHO$Season=='Winter'),]
mean(box_CHO_w$Quant, na.rm=TRUE)
median.default(box_CHO_w$Quant, na.rm = TRUE)
box_CHO_s <- boxplot_CHO[ which(boxplot_CHO$Season=='Summer'),]
meadian(box_CHO_s$Quant, na.rm=TRUE)
median.default(box_CHO_s$Quant, na.rm = TRUE)
85774/97062

box_CHOS_w <- boxplot_CHOS[ which(boxplot_CHOS$Season=='Winter'),]
mean(box_CHOS_w$Quant, na.rm=TRUE)
median.default(box_CHOS_w$Quant, na.rm = TRUE)
box_CHOS_s <- boxplot_CHOS[ which(boxplot_CHOS$Season=='Summer'),]
mean(box_CHOS_s$Quant, na.rm=TRUE)
median.default(box_CHOS_s$Quant, na.rm = TRUE)
157537/176518

box_CHON_w <- boxplot_CHON[ which(boxplot_CHON$Season=='Winter'),]
mean(box_CHON_w$Quant, na.rm=TRUE)
median.default(box_CHON_w$Quant, na.rm = TRUE)
box_CHON_s <- boxplot_CHON[ which(boxplot_CHON$Season=='Summer'),]
mean(box_CHON_s$Quant, na.rm=TRUE)
median.default(box_CHON_s$Quant, na.rm = TRUE)
1432207/685064

box_CHONS_w <- boxplot_CHONS[ which(boxplot_CHONS$Season=='Winter'),]
mean(box_CHONS_w$Quant, na.rm=TRUE)
median.default(box_CHONS_w$Quant, na.rm = TRUE)
box_CHONS_s <- boxplot_CHONS[ which(boxplot_CHONS$Season=='Summer'),]
mean(box_CHONS_s$Quant, na.rm=TRUE)
median.default(box_CHONS_s$Quant, na.rm = TRUE)
8024/17994


write.csv(Averages, file = "Avgmetricswinter.csv")

massloading2 <- massloading %>% 
  select(datetime,
         Volume)

massloading2$datetime <- dmy_hm(massloading2$datetime, tz = "UTC")

massloading3 <- massloading2%>%
  left_join(Averages, by = "datetime")

for(x in 1:length(massloading3$OCratio))
  if (massloading3$Volume[x] != 'NA') {massloading3$minutes[x] <- (massloading3$Volume[x]/1.3)} else{massloading3$minutes[x] <- NA}

weighted.mean(massloading3$OCratio, massloading3$minutes, na.rm = TRUE)



