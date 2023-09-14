#NMDS 16S 
library(ggplot2)
library(vegan)
library(readxl)
library(mixOmics)

table_all_16 <- read_delim("16S_all.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
meta_all_16 <-read_delim("meta_all_16S.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)


#I will formate the table....
tab_all16 <-table_all_16[,-1]

sums<-rowSums(tab_all16)
min(sums)
#[1] 3574

tabN_all_16 <-as.data.frame.array(rrarefy(tab_all16,sample=3574))

library(ape)
pco_all16<-pcoa(vegdist(tabN_all_16,"bray"))
biplot.pcoa(pco_all16)


#PCA

pcoa_all16 <- data.frame(pco_all16$vectors)

dat_pcoall16 <- pcoa_all16[, c("Axis.1", "Axis.2")]

dat_pcoall16
meta_all_16
dat_pcoall16 <- cbind(dat_pcoall16, meta_all_16$Sample)


pcoa_all16 <- data.frame(pco_all16$vectors)


ggplot(dat_pcoall16,aes(x = Axis.1, y = Axis.2,color = meta_all_16$Sample )) +
  geom_jitter() +
  scale_shape_manual(values = c(15,16,17)) + theme_minimal()+
  stat_ellipse(data=dat_pcoall16, aes(x=Axis.1, y = Axis.2, color = meta_all_16$Sample ))+
  scale_color_manual(values=c("#99CCFF","#66CC99","#CC3300")) +
  labs(title = " 16S", color ="Samples")


adonis2(tabN_all_16 ~ Sample, data = meta_all_16, permutations = 999, method= "bray")
#adonis2(formula = tabN_all_16 ~ Sample, data = meta_all_16, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#Sample     2    4.756 0.04155 5.8314  0.001 ***
#Residual 269  109.697 0.95845                  
#Total    271  114.453 1.00000                  


#NMDS itS 
library(ggplot2)
library(vegan)
library(readxl)

all_moy_ITS <- read_delim("ITS_all.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
metaall_its <- read_delim("meta_all_ITS.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)

library(mixOmics)
library(vegan)


#I will formate the table....

View(tab_allits)
tab2_allITS <-all_moy_ITS[,-1]

tab_allITS2 <- as.data.frame(tab2_allITS)
tab2

sums1<-rowSums(tab_allITS2)
min(sums1)
#[1] 3475 

tabN_all_its <-as.data.frame.array(rrarefy(tab_allITS2,sample=3475))

library(ape)
pco_allits<-pcoa(vegdist(tabN_all_its,"bray"))
biplot.pcoa(pco_allits)

#PCA
pcoa_allits <- data.frame(pco_allits$vectors)

dat_pcoallits <- pcoa_allits[, c("Axis.1", "Axis.2")]
dat_pcoallits <- cbind(dat_pcoallits, metaall_its$Sample, metaall_its$PC1)

ggplot(dat_pcoallits,aes(x = Axis.1, y = Axis.2,color = metaall_its$Sample )) +
  geom_jitter() +
  scale_shape_manual(values = c(15,16,17)) + theme_minimal()+
  stat_ellipse(data=dat_pcoallits, aes(x=Axis.1, y = Axis.2, color = metaall_its$Sample ))+
  scale_color_manual(values=c("#99CCFF","#66CC99","#CC3300")) +
  labs(title = " ITS", color ="Samples")


adonis2(tabN_all_its ~ Sample, data = metaall_its, permutations = 999, method= "bray")
#adonis2(formula = tabN_all_its ~ Sample, data = metaall_its, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#Sample     2    2.301 0.02639 2.8322  0.001 ***
#  Residual 209   84.886 0.97361                  
#Total    211   87.186 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#NMDS 18S 

library(ggplot2)
library(vegan)
library(readxl)

table_all_18 <- read_delim("18S_all.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)

meta_all_18 <-read_delim("meta_all_ITS.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)

library(mixOmics)
library(vegan)

#I will formate the table....

tab_all18 <-table_all_18[,-1]

sums<-rowSums(tab_all18)
min(sums)
#[1] 6702

tabN_all_18 <-as.data.frame.array(rrarefy(tab_all18,sample=6702))

library(ape)
pco_all18<-pcoa(vegdist(tabN_all_18,"bray"))
biplot.pcoa(pco_all18)
  
#PCA
pcoa_all18 <- data.frame(pco_all18$vectors)

dat_pcoall18 <- pcoa_all18[, c("Axis.1", "Axis.2")]
dat_pcoall18 <- cbind(dat_pcoall18, meta_all_18$condition_trie, meta_all_18$X.PC1)



ggplot(dat_pcoall18,aes(x = Axis.1, y = Axis.2,color = meta_all_18$condition_trie )) +
  geom_jitter() +
  scale_shape_manual(values = c(15,17)) + theme_minimal()+
  stat_ellipse(data=dat_pcoall18, aes(x=Axis.1, y = Axis.2, color = meta_all_18$condition_trie ))+
  scale_color_manual(values=c("#99CCFF","#CC3300")) +
  labs(title = "18S", color ="Samples")

ggplot(dat_pcoall18) +
  aes(
    x = Axis.1,
    y = Axis.2,
    colour = `meta_all_18$condition_trie`) +
  geom_point(shape = "circle", size = 1.5) +
  scale_color_hue(direction = 1) +
  theme_minimal()

adonis2(tabN_all_18 ~ condition_trie, data = meta_all_18, permutations = 999, method= "bray")
#adonis2(formula = tabN_all_18 ~ condition_trie, data = meta_all_18, permutations = 999, method = "bray")
#Df SumOfSqs      R2     F Pr(>F)
#condition_trie  1   0.5396 0.04557 1.146   0.12
#Residual       24  11.2992 0.95443             
#Total          25  11.8388 1.00000




