library(tidyverse)
library(data.table)
library(devtools)
load_all("/dcs04/scharpf/data/annapragada/useful.stuff.aa")
library(readxl)

meta<-fread("../data/All_LCr_Clinical_Metadata_withscores.csv",header=T) %>% select(-V1)
meta2<-fread("../data/All_Cross_Reactivity_Metadata_withscores.csv",header=T) %>% select(-V1) %>% filter(existing_val != "Yes")


#DELFI Feature Matrix output from delfi pipeline
delfi<-fread("allfeatures.5mb.hg19.csv")  %>% filter(id %in% meta$id|id %in% meta2$id)
delfi$arm<-paste0("zscore_",delfi$arm)
delfi<-delfi %>% spread(key=arm,value=zscore)


#Repeat element feature matrix output from artemis pipeline -- normalized counts
artemis<-fread("artemis.csv")

artemis<-artemis %>% select(-V1)
artemis<-tibble(artemis)
artemis<-artemis %>% filter(id %in% meta$id|id %in% meta2$id)

#Epigenetic bins feature matrix output from artemis pipeline
epi<-fread("epi_bins.csv")
epi<-epi %>% select(-V1)
epi$id<-gsub(".hg19","",epi$id)
epi<-epi %>% filter(id %in% meta$id|id %in% meta2$id)

e<-fread("Expected.csv") #This reference is provided in ARTEMIS pipeline, use to center and scale the features
e<-e %>% filter(total_kmers>1000)

e$fam<-sapply(str_split(e$feature,"#"),"[",2)
e$fam<-sapply(str_split(e$fam,"_"),"[",1)
e<-e %>% mutate(fam=if_else(is.na(fam),"Satellite",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("rRNA","snRNA","scRNA","tRNA","srpRNA"),"RNA_DNA Elements",fam))
e<-e %>% mutate(fam=if_else(fam %in% c("DNA","DNA?","RC","Retroposon"),"RNA_DNA Elements",fam))


artemis<-artemis %>% select(id,e$feature)
test<-artemis %>% gather(key=feature,value=count,-id)
test<-inner_join(test,e %>% select(feature,fam),by="feature")
test<-test %>% group_by(id,fam) %>% summarize(c=scale(count)[,1],f=feature)

#add "class_fam"
test$f<-paste0("class_",test$fam,"_",test$f)

test<-test %>% ungroup() %>% select(-fam)
test<-test %>% spread(key=f,value=c)

artemis<-test
delfi<-delfi %>% select(-starts_with("cov"))

data<-inner_join(delfi,artemis,by="id")


#epi feature classes - H3K27me3, H3K36me3, H3K9me3, H4K20me1, states_1_5, states_10_13, states_7_9


data<-inner_join(data,epi,by="id")

#######

data2<-data
meta1<-meta %>% filter(training=="Discovery")
meta1<-meta1 %>% mutate(type=if_else(type=="No Known Liver Disease","healthy","cancer"))
data2<-inner_join(data2,meta1 %>% select(type,id),by="id")

write.csv(data2,"Cirrhosis_Train.csv")

data2<-data
meta1<-meta %>% filter(training=="Validation")
meta1<-meta1 %>% mutate(type=if_else(type=="No Known Liver Disease","healthy","cancer"))
meta1<-meta1 %>% select(id,type)
meta2<-meta2 %>% mutate(type="healthy") %>% select(id,type)
meta2<-rbind(meta1,meta2)



data2<-inner_join(data2,meta2 %>% select(type,id),by="id")

write.csv(data2,"Cirrhosis_Test.csv")



