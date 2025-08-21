library(tidyverse)
mydata=read.csv("pfpr_intvn.csv")
newData=filter(mydata, year==2000)
newData[,4]
newData[,c(4)]
newData[,c(3,4)]#selecting/viewing multiple columns
sum(newData[,4])
names(mydata)
sum(mydata[,4])
min(mydata$itn_coverage)
max(mydata$itn_coverage)
mean(mydata$itn_coverage)
tidyverse_packages()
names(mydata)
newDataa=mydata %>% select(county,residence,year,pfpr)
newDataa=mydata %>% select(-county,-residence,year)
head(newDataa)
newDataa=mydata |> filter(county %in% c("Kisumu","Baringo","Machakos"))
newDataa
newDatao=mydata |> select(-county,-year)
newDatao
newDataa=mydata |> filter(county %in% c("Kisumu","Baringo","Machakos"),residence=="urban")
newDataa
#mutate function
newDataa=mutate(mydata,pfprRate=pfpr*100)
#group_by() and summarise()
newDataa=mydata %>% group_by(county)%>% summarise(pfpr_Avg=mean(pfpr))
newDataa
#grouping using more than one column
newDataa=mydata %>% group_by(county,year)%>% summarise(pfpr_Avg=mean(pfpr))
newDataa
newDataa=mydata %>% group_by(county,year)%>% summarise(pfpr_Avg=round(mean(pfpr),2))# get the mean per county in each year in 2dp
newDataa
