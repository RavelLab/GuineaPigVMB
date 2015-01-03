"Barplots of the relative abundance of the top 30 phylotypes found in the guinea pig vaginal microbiota"
------------------------
  
  
#set working directory
setwd("../../data")
# set a CRAN mirror
local({r <- getOption("repos")
       r["CRAN"] <- "http://lib.stat.cmu.edu/R/CRAN/"
       options(repos=r)})
#install plyr and ggplots
install.packages(c("ggplot2","plyr"),dep=T)

library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("knitr", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("lattice", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("labeling", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("memoise", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("plyr", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")
library("scales", lib.loc="/Library/Frameworks/R.framework/Versions/3.1/Resources/library")

#Upload Data (Healthy)

Data <- read.csv("n_in_barplot.csv", header=TRUE)
head(Data)

#Select Taxon (Healthy)
Data$taxa <- factor(Data$taxa,levels=c("Corynebacterium","Anaerococcus","Clostridiaceae(Family)","Peptoniphilus","Aerococcus","Facklamia","Acinetobacter","Actinomyces","Actinomycetales.Order.","Aerococcaceae(Family)","Allobaculum","Arcanobacterium","Bacteroidales_S24-7(Family)","Clostridiaceae.otu1","Enterococcus","Erysipelotrichacea_PSB-M-3","Escherichia","Helcococcus","Porphyromonadaceae.Family.","Porphyromonas","Proteus","Peptococcus","Ruminococcaceae.Family.","Rummeliibacillus","Solibacillus","Staphylococcus","Streptococcaceae.Family.","Weissella","Chlamydia","Other"))

#Plot mean (Healthy)

p <-ggplot(Data, aes(x = Sample_ID)) + geom_bar(aes(weight=RA, fill = taxa), 
                                                position = 'fill') + scale_fill_manual(values=c("#F69799","#FACBCC","#FCC199","#FCDDCA","#FCE897","#FFF3CC","#E7EC99","#F1F4CE","#C5DD91","#E0EDC8","#A6D38A","#CFE7C3","#A6D6B2","#CFE9D7","#A8DBD6","#D1EBE9","#A3DDF2","#CDEBF5","#9DBDE2","#CCDBEF","#9494C8","#CBC8E3","#B098C7","#DACAE3","#CC9DC8","#E5CBE2","#E39CC4","#F0CCE1","#EB2627","#FFFFFF")) + 
  facet_grid(.~ animal + cycle, scale="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p

#note colors for scale_fill_manual
#scale_fill_manual(values=c("#F69799","#FACBCC","#FCC199","#FCDDCA","#FCE897","#FFF3CC"
#,"#E7EC99","#F1F4CE","#C5DD91","#E0EDC8","#A6D38A","#CFE7C3","#A6D6B2","#CFE9D7","#A8DBD6"
#,"#D1EBE9","#A3DDF2","#CDEBF5","#9DBDE2","#CCDBEF","#9494C8","#CBC8E3","#B098C7","#DACAE3"
#,"#CC9DC8","#E5CBE2","#E39CC4","#F0CCE1","#EB2627","#FFFFFF"))

#note script for taxa/phlotypes
#Data$taxa <- factor(Data$taxa,levels=c("Corynebacterium","Anaerococcus","Clostridiaceae
# (Family)","Peptoniphilus","Aerococcus","Facklamia","Acinetobacter","Actinomyces",
# "Actinomycetales.Order.","Aerococcaceae(Family)","Allobaculum","Arcanobacterium",
# "Bacteroidales_S24-7(Family)","Clostridiaceae.otu1","Enterococcus",
#"Erysipelotrichacea_PSB-M-3","Escherichia","Helcococcus","Porphyromonadaceae.Family.",
#"Porphyromonas","Proteus","Peptococcus","Ruminococcaceae.Family.","Rummeliibacillus",
#"Solibacillus","Staphylococcus","Streptococcaceae.Family.","Weissella","Chlamydia","Other"))


#Upload Data (Mock-infected)
Data <- read.csv("mk_barplot.csv", header=TRUE)
#Select Taxon (Mock-infected)
Data$taxa <- factor(Data$taxa,levels=c("Corynebacterium","Anaerococcus","Clostridiaceae(Family)","Peptoniphilus","Aerococcus","Facklamia","Acinetobacter","Actinomyces","Actinomycetales.Order.","Aerococcaceae(Family)","Allobaculum","Arcanobacterium","Bacteroidales_S24-7(Family)","Clostridiaceae.otu1","Enterococcus","Erysipelotrichacea_PSB-M-3","Escherichia","Helcococcus","Porphyromonadaceae.Family.","Porphyromonas","Proteus","Peptococcus","Ruminococcaceae.Family.","Rummeliibacillus","Solibacillus","Staphylococcus","Streptococcaceae.Family.","Weissella","Chlamydia","Other"))

#Plot mean (Mock infected)

p <-ggplot(Data, aes(x = Sample_ID)) + geom_bar(aes(weight=RA, fill = taxa), position = 'fill')+ 
  scale_fill_manual(values=c("#F69799","#FACBCC","#FCC199","#FCDDCA","#FCE897","#FFF3CC","#E7EC99","#F1F4CE","#C5DD91","#E0EDC8","#A6D38A","#CFE7C3","#A6D6B2","#CFE9D7","#A8DBD6","#D1EBE9","#A3DDF2","#CDEBF5","#9DBDE2","#CCDBEF","#9494C8","#CBC8E3","#B098C7","#DACAE3","#CC9DC8","#E5CBE2","#E39CC4","#F0CCE1","#EB2627","#FFFFFF")) +
  facet_grid(.~ animal + cycle, scale="free") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
#note colors for scale_fill_manual
#scale_fill_manual(values=c("#F69799","#FACBCC","#FCC199","#FCDDCA","#FCE897","#FFF3CC"
#,"#E7EC99","#F1F4CE","#C5DD91","#E0EDC8","#A6D38A","#CFE7C3","#A6D6B2","#CFE9D7","#A8DBD6"
#,"#D1EBE9","#A3DDF2","#CDEBF5","#9DBDE2","#CCDBEF","#9494C8","#CBC8E3","#B098C7","#DACAE3"
#,"#CC9DC8","#E5CBE2","#E39CC4","#F0CCE1","#EB2627","#FFFFFF"))

#note script for taxa/phlotypes
#Data$taxa <- factor(Data$taxa,levels=c("Corynebacterium","Anaerococcus","Clostridiaceae
# (Family)","Peptoniphilus","Aerococcus","Facklamia","Acinetobacter","Actinomyces",
# "Actinomycetales.Order.","Aerococcaceae(Family)","Allobaculum","Arcanobacterium",
# "Bacteroidales_S24-7(Family)","Clostridiaceae.otu1","Enterococcus",
#"Erysipelotrichacea_PSB-M-3","Escherichia","Helcococcus","Porphyromonadaceae.Family.",
#"Porphyromonas","Proteus","Peptococcus","Ruminococcaceae.Family.","Rummeliibacillus",
#"Solibacillus","Staphylococcus","Streptococcaceae.Family.","Weissella","Chlamydia","Other"))




#Upload Data (Infected)
Data <- read.csv("in_barplot.csv", header=TRUE)

Data$taxa <- factor(Data$taxa,levels=c("Corynebacterium","Anaerococcus","Clostridiaceae(Family)","Peptoniphilus","Aerococcus","Facklamia","Acinetobacter","Actinomyces","Actinomycetales.Order.","Aerococcaceae(Family)","Allobaculum","Arcanobacterium","Bacteroidales_S24-7(Family)","Clostridiaceae.otu1","Enterococcus","Erysipelotrichacea_PSB-M-3","Escherichia","Helcococcus","Porphyromonadaceae.Family.","Porphyromonas","Proteus","Peptococcus","Ruminococcaceae.Family.","Rummeliibacillus","Solibacillus","Staphylococcus","Streptococcaceae.Family.","Weissella","Chlamydia","Other"))

#Plot mean (Infected)

p <-ggplot(Data, aes(x = Sample_ID)) + geom_bar(aes(weight=RA, fill = taxa), position = 'fill') + scale_fill_manual(values=c("#F69799","#FACBCC","#FCC199","#FCDDCA","#FCE897","#FFF3CC","#E7EC99","#F1F4CE","#C5DD91","#E0EDC8","#A6D38A","#CFE7C3","#A6D6B2","#CFE9D7","#A8DBD6","#D1EBE9","#A3DDF2","#CDEBF5","#9DBDE2","#CCDBEF","#9494C8","#CBC8E3","#B098C7","#DACAE3","#CC9DC8","#E5CBE2","#E39CC4","#F0CCE1","#EB2627","#FFFFFF")) + facet_grid(.~ animal + cycle, scale="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
#note colors for scale_fill_manual
#scale_fill_manual(values=c("#F69799","#FACBCC","#FCC199","#FCDDCA","#FCE897","#FFF3CC"
#,"#E7EC99","#F1F4CE","#C5DD91","#E0EDC8","#A6D38A","#CFE7C3","#A6D6B2","#CFE9D7","#A8DBD6"
#,"#D1EBE9","#A3DDF2","#CDEBF5","#9DBDE2","#CCDBEF","#9494C8","#CBC8E3","#B098C7","#DACAE3"
#,"#CC9DC8","#E5CBE2","#E39CC4","#F0CCE1","#EB2627","#FFFFFF"))

#note script for taxa/phlotypes
#Data$taxa <- factor(Data$taxa,levels=c("Corynebacterium","Anaerococcus","Clostridiaceae
# (Family)","Peptoniphilus","Aerococcus","Facklamia","Acinetobacter","Actinomyces",
# "Actinomycetales.Order.","Aerococcaceae(Family)","Allobaculum","Arcanobacterium",
# "Bacteroidales_S24-7(Family)","Clostridiaceae.otu1","Enterococcus",
#"Erysipelotrichacea_PSB-M-3","Escherichia","Helcococcus","Porphyromonadaceae.Family.",
#"Porphyromonas","Proteus","Peptococcus","Ruminococcaceae.Family.","Rummeliibacillus",
#"Solibacillus","Staphylococcus","Streptococcaceae.Family.","Weissella","Chlamydia","Other"))
