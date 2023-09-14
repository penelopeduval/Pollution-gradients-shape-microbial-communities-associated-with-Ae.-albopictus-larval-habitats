library(readxl)
library(mixOmics)
library(ade4)


#physicochemical analysis combined with pollution variables

table <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/toute_var.xlsx", 
                    col_types = c("text", "text", "text", 
                                  "text", "text", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric", 
                                  "numeric", "numeric", "numeric"))

NC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/toute_var.xlsx", 
                 sheet = "NC", col_types = c("text", "text", 
                                             "text", "text", "text", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric"))

EC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/toute_var.xlsx", 
                 sheet = "EC", col_types = c("text", "text", 
                                             "text", "text", "text", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric"))

pca <- pca(table[,c(7:31)],scale=TRUE, center=TRUE)


plot(pca)
Colonisation<-table$Colonisation
biplot(pca, group=Colonisation, ind.names = FALSE, legend = FALSE, labels = FALSE)

#PCA contruction for colonized and uncolonized water samples 
pca_nc <- pca(NC[,c(7:31)],scale=TRUE, center=TRUE)
biplot(pca_nc, ind.names = FALSE, legend = FALSE)
pca_ec <- pca(EC[,c(7:31)],scale=TRUE, center=TRUE)
biplot(pca_ec, ind.names = FALSE)

#procurste analysis 
pro<-procuste(pca_nc$variates$X,pca_ec$variates$X)
plot(pro)
X<-pca_ec$variates$X

randtest(pro)

biplot(pca_nc)

biplot(pca_ec)

pls1 <- plsda(table[,c(7:31)], table$Colonisation)
plotIndiv(pls1,legend=T, ellipse = T, star=T)
plotVar(pls1)
biplot(pls1)
plotLoadings(pls1, contrib='max', method='mean', legend = FALSE, title =  )

lm1<-lm(1/(1+table$N2O_Vol)~table$Colonisation)
library(car)
library(MASS)
plot(lm1)
qqPlot(residuals(lm1))
Anova(lm1)

wilcox.test(table$N2O_Vol~table$Colonisation)

library(ggplot2)

ggplot(table, aes(x=Colonisation, y=N2O_Vol))+geom_boxplot()


library(car)
library(MASS)
wilcox.test(table$CH4b_Vol~table$Colonisation)

pca_nc$variates->PC_nc
pca_ec$variates->PC_ec


#impact on 16S communities

table_EC16 <- read_delim("table_EC_moy16S.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
meta_EC_16 <- read_delim("meta_EC_moy16S.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)


table_NC16 <- read_delim("table_NC_moy16S.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
meta_NC_16 <- read_delim("table_NC_moy16S.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)

tab2_EC<-table_EC16[,-1]
tab2_NC<-table_NC16[,-1]

library(vegan)
tab2_NC <-as.data.frame.array(rrarefy(tab2_NC,sample= 3574))
tab2_EC <-as.data.frame.array(rrarefy(tab2_EC,sample= 3574))


library(ape)
pco_nc<-pcoa(vegdist(tab2_NC,"bray"))
biplot.pcoa(pco_nc)

pco_ec<-pcoa(vegdist(tab2_EC,"bray"))
biplot.pcoa(pco_ec)


library(ade4)

pro_M <- procuste(pco_nc$vectors[,c(1,2)], pco_ec$vectors[,c(1,2)])
plot(pro_M)
randtest(pro_M)

PC_ec1<-as.data.frame.array(PC_ec$X)
db_ec <- capscale(tab2_EC~PC_ec1$PC1+PC_ec1$PC2, dist="bray")

adonis2(tab2_EC~PC_ec1$PC1+PC_ec1$PC2, dist="bray")

ordistep(db_ec)->oec
oec

#colonized water samples
plot(oec)
coorec<-as.data.frame.array(oec$CCA$v)
coorec[which(coorec$CAP1<=-0.3),]

df_ec<-cbind(PC_ec1$PC1,tab2_EC$Cluster16S_13,tab2_EC$Cluster16S_14, PC_ec1$PC2)
df_ec<-as.data.frame.array(df_ec)
library(ggplot2)

ggplot(df_ec, aes(x =V1, y=log10(1+V2)))+geom_point()
ggplot(df_ec, aes(x =V1, y=log10(1+V3)))+geom_point()

###V1= PC1, V2=Cluster16S_13, V3=Cluster16S_14

cor.test(df_ec$V1,df_ec$V2, method="spearman")


cor.test(df_ec$V4,df_ec$V2, method="spearman") #otu 13 corrélée avec PC2 

cor.test(df_ec$V1,df_ec$V3, method="spearman")

coorec<-as.data.frame.array(oec$CCA$v)
coorec[which(coorec$CAP2<=-0.3),]

df_ecPC2<-cbind(PC_ec1$PC2,tab2_EC$Cluster16S_13,tab2_EC$Cluster16S_14)
df_ecPC2<-as.data.frame.array(df_ecPC2)

cor.test(df_ecPC2$V1,df_ecPC2$V2, method="spearman")

###V1= PC1, V2=Cluster16S_13, V3=Cluster16S_14

cor.test(df_ecPC2$V1,df_ecPC2$V3, method="spearman")
cor.test(df_ec$V1,df_ec$V3, method="spearman")

# uncolonized water samples 

PC_nc1<-as.data.frame.array(PC_nc$X)
db_nc <- capscale(tab2_NC~PC_nc1$PC1+PC_nc1$PC2, dist="bray")
plot(db_nc)

adonis2(tab2_NC~PC_nc1$PC1+PC_nc1$PC2, dist="bray")

ordistep(db_nc)->onc

onc <- capscale(tab2_NC~PC_nc1$PC1+PC_nc1$PC2, dist="bray")

onc <- capscale(tab2_NC~PC_nc1$PC1, dist="bray")
plot(onc)
coornc<-as.data.frame.array(onc$CCA$v)
rownames(coornc[which(coornc$CAP1<=-0.1),])

df_nc<-cbind(PC_nc1$PC1,tab2_NC$Cluster16S_50,tab2_NC$Cluster16S_142)
df_nc<-as.data.frame.array(df_nc)
library(ggplot2)

ggplot(df_nc, aes(x =V1, y=log10(1+V2)))+geom_point()
ggplot(df_nc, aes(x =V1, y=log10(1+V3)))+geom_point()


cor.test(df_nc$V1,df_nc$V2, method="spearman")

cor.test(df_nc$V1,df_nc$V3, method="spearman")


#impact on ITS communities

table_Ec_ITS <- read_delim("table_EC_moyITS.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
meta_Ec_ITS <- read_delim("metadata_EC_ITS.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)


table_Nc_ITS <- read_delim("table_NC_moyITS.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
meta_Nc_ITS <- read_delim("metadata_NC_ITS.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)


library(mixOmics)
library(vegan)


#I will formate the table....
tab6_EC<-table_Ec_ITS[,-1]
tab6_NC<-table_Nc_ITS[,-1]

sums<-rowSums(tab6_EC)
min(sums)
#[1] 5341
sums<-rowSums(tab6_NC)
min(sums)
#[1] 3475

norm_tab_EC<-as.data.frame.array(rrarefy(tab6_EC,sample=3475))
tab7_EC<-as.data.frame.array(norm_tab_EC)

norm_tab_NC<-as.data.frame.array(rrarefy(tab6_NC,sample=3475))
tab7_NC<-as.data.frame.array(norm_tab_NC)


library(ape)
pco_nc<-pcoa(vegdist(tab7_NC,"bray"))
biplot.pcoa(pco_nc)


pco_ec<-pcoa(vegdist(tab7_EC,"bray"))
biplot.pcoa(pco_ec)


library(ade4)


PC_poll <- write.table(PC_ec1)

pro_M <- procuste(pco_nc$vectors[,c(1,2)], pco_ec$vectors[,c(1,2)])
plot(pro_M)
randtest(pro_M)


# colonized water samples 

Pollution_EC_ITS <- read_excel("Pollution_E.xlsx", sheet = "EC_ITS", col_types = c("text","numeric", "numeric"))
View(Pollution_EC_ITS)

adonis2(tab7_EC~PC_ec1$X.PC1+PC_ec1$X.PC2)

PC_poll1 <- Pollution_EC_ITS[,-1]
PC_ec1<-as.data.frame(PC_poll1)
db_ec <- capscale(tab7_EC~PC_ec1$X.PC1+PC_ec1$X.PC2, dist="bray")
plot(db_ec)

db_ec

ordistep(db_ec)->oec
plot(oec)


ordistep(db_ec)-> oec
oec


db_ec1 <- capscale(tab7_EC~PC_ec1$X.PC1, dist="bray")
plot(db_ec1)


### correlation with PC1 

coorec<-as.data.frame.array(db_ec$CCA$v)

rownames(coorec[which(db_ec$CCA$v<=-0.3),])
#"ClusterITS_11"

rownames(coorec[which(coorec$CAP1>0.3),])
# "ClusterITS_92" "ClusterITS_6" 

df_ec<-cbind(PC_ec1$X.PC1,tab7_EC$ClusterITS_11,tab7_EC$ClusterITS_92, tab7_EC$ClusterITS_6, PC_ec1$X.PC2)
df_ec<-as.data.frame.array(df_ec)
library(ggplot2)

ggplot(df_ec, aes(x =V1, y=log10(1+V2)))+geom_point()
ggplot(df_ec, aes(x =V1, y=log10(1+V3)))+geom_point()
ggplot(df_ec, aes(x =V1, y=log(1+V4)))+geom_point()

###V1= PC1, V2=ClusterITS_11, V3=ClusterITS_92, V4=ClusterITS_6

cor.test(df_ec$V1,df_ec$V2, method="spearman")

cor.test(df_ec$V5,df_ec$V2, method="spearman")

cor.test(df_ec$V1,df_ec$V3, method="spearman")

cor.test(df_ec$V1,df_ec$V4, method="spearman")


# uncolonized water samples 
Pollution_NC_ITS <- read_excel("Pollution_E.xlsx", 
                               sheet = "NC_ITS", col_types = c("text", 
                                                               "numeric", "numeric"))

PC_poll_NC <- Pollution_NC_ITS[,-1]

PC_NC4<-as.data.frame(PC_poll_NC)
db_nc_its <- capscale(tab7_NC~PC_NC4$X.PC1+PC_NC4$X.PC2, dist="bray")
plot(db_nc_its)

db_nc_its


ordistep(db_nc_its)->onc
plot(onc)


adonis2(tab7_NC~PC_NC1$X.PC1+PC_NC1$X.PC2)
         

### correlation with PC1 
db_nc1 <- capscale(tab7_NC~PC_NC4$X.PC1, dist="bray")
plot(db_nc1)
coornc<-as.data.frame.array(db_nc1$CCA$v)

coornc$CAP1

rownames(coornc[which(coornc$CAP1>0.3),])
#[1] "ClusterITS_30"

rownames(coornc[which(coornc$CAP1<=-0.3),])
#[1] "ClusterITS_60" "ClusterITS_24"

coord_PC1_NC <- data.frame(coornc$CAP1)
 
write.table(coord_PC1_NC,"coord_PC1_NC.txt",row.names=F,col.names=F,quote=F,sep='\t')

df_nc<-cbind(PC_NC4$X.PC1,tab7_NC$ClusterITS_30,tab7_NC$ClusterITS_60, tab7_NC$ClusterITS_24)

df_nc<-as.data.frame.array(df_nc)
library(ggplot2)

ggplot(df_nc, aes(x =V1, y=log10(1+V2)))+geom_point()
ggplot(df_nc, aes(x =V1, y=log10(1+V3)))+geom_point()
ggplot(df_nc, aes(x =V1, y=log10(1+V4)))+geom_point()


cor.test(df_nc$V1,df_nc$V2, method="spearman")

cor.test(df_nc$V1,df_nc$V3, method="spearman")

cor.test(df_nc$V1,df_nc$V4, method="spearman")


library(ggplot2)

ggplot(dat_pcoNC1) +
  aes(x = Axis.1, y = Axis.2, colour = `as.numeric(meta_Nc_ITS$X.PC1)`) +
  geom_point(shape = "circle", 
             size = 4L) +
  scale_color_gradient(low = "#F3D43B", high = "#132B43") +
  labs(title = "Non colonized water ITS ", 
       color = "Pollution (PC1) ") +
  theme_minimal() +
  xlim(-0.4, 0.4)

ggplot(dat_pcoEC1) +
  aes(x = Axis.1, y = Axis.2, colour = dat_pcoEC1$`Pollution_EC_ITS$X.PC1`) +
  geom_point(shape = "circle", 
             size = 4L) +
  scale_color_gradient(low = "#F3D43B", high = "#132B43") +
  labs(title = "Colonized water ITS ", 
       color = "Pollution (PC1) ") +
  theme_minimal() +
  xlim(-0.4, 0.4)

#impact on 18S communities
library(readxl)

table_EC_18S  <- read_delim("table_EC_moy18S.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
table_NC_18S  <- read_delim("table_NC_moy18S.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)


poll_NC_18S <- read_excel("Pollution_E.xlsx", sheet = "NC_18S", col_types = c("text", "numeric", "numeric"))
poll_EC_18S <- read_excel("Pollution_E.xlsx", sheet = "EC_18S", col_types = c("text", "numeric", "numeric"))


library(mixOmics)
library(vegan)



#I will formate the table....
T18S_EC<-table_EC_18S[,-1]
T18S_NC<-table_NC_18S[,-1]

sums<-rowSums(T18S_EC)
min(sums)
#[1] 6702
sums<-rowSums(T18S_NC)
min(sums)
#[1] 7195

norm_18S_EC<-as.data.frame.array(rrarefy(T18S_EC,sample=6702))
tab8_EC<-as.data.frame.array(norm_18S_EC)

norm_18S_NC<-as.data.frame.array(rrarefy(T18S_NC,sample=6702))
tab8_NC<-as.data.frame.array(norm_18S_NC)


library(ape)
pco_18Snc<-pcoa(vegdist(tab8_NC,"bray"))
biplot.pcoa(pco_18Snc)


pco_18Sec<-pcoa(vegdist(tab8_EC,"bray"))
biplot.pcoa(pco_18Sec)


#Procruste analysis
library(ade4)
pro_18S <- procuste(pco_18Snc$vectors[,c(1,2)], pco_18Sec$vectors[,c(1,2)])
plot(pro_18S)
randtest(pro_18S)


# colonized water samples 

poll_EC18s <- poll_EC_18S[,-1]
P_EC18s<-as.data.frame(poll_EC18s)

db_ec_18S <- capscale(tab8_EC~P_EC18s$X.PC1+P_EC18s$X.PC2, dist="bray")
plot(db_ec_18S)

db_ec_18S

adonis2(tab8_EC~P_EC18s$X.PC1+P_EC18s$X.PC2, dist="bray")


ordistep(db_ec_18S)->oec_18S
plot(oec_18S)


df_ec18<-cbind(P_EC18s$X.PC1, tab8_EC$Cluster18S_10, tab8_EC$Cluster18S_5, tab8_EC$Cluster18S_17, tab8_EC$Cluster18S_3)
df_ec18s<-as.data.frame.array(df_ec18)
library(ggplot2)



###V2 = Cluster18S_10, V3 = Cluster18S_5, V4 = Cluster18S_17, V5 = Cluster18S_3

cor.test(df_ec18s$V1,df_ec18s$V2, method="spearman")

cor.test(df_ec18s$V1,df_ec18s$V3, method="spearman")

cor.test(df_ec18s$V1,df_ec18s$V4, method="spearman")

cor.test(df_ec18s$V1,df_ec18s$V5, method="spearman")
 

########NC

poll_NC18s <- poll_NC_18S[,-1]
P_NC18s<-as.data.frame(poll_NC18s)

db_nc_18S <- capscale(tab8_NC~P_NC18s$X.PC1+P_NC18s$X.PC2, dist="bray")
plot(db_nc_18S)

adonis2(tab8_NC~P_NC18s$X.PC1+P_NC18s$X.PC2, dist="bray")

View(db_nc_18S)

#jolie DBRDA 
bp_18S_v <- data.frame(db_nc_18S$CCA$v) #coordonnées pour les OTU sur CAP1 et CAP2 
bp_18S_u <- data.frame(db_nc_18S$CCA$u) #coordonnées pour les échantillons 
View(bp_18S_u)
bp_18S_w <- data.frame(db_nc_18S$CCA$w) #coordonnées pour les fleches PC1 
View(bp_18S_w)

library(vegan)
library(dplyr)
library(ggplot2)
library(grid)


ordistep(db_nc_18S)->onc_18s
plot(onc_18s)


# correlation with PC1  
db_nc18S1 <- capscale(tab8_NC~P_NC18s$X.PC1, dist="bray")
plot(db_nc18S1)
coornc_18S<-as.data.frame.array(db_nc18S1$CCA$v)

coornc_18S1 <- as.data.frame(coornc_18S)



coord_PC1_NC18 <- data.frame(coornc_18S$CAP1)
write.table(coornc_18S,"coord_PC1_NC18.txt",row.names=T,col.names=F,quote=F,sep='\t')


df_nc18<-cbind(P_NC18s$X.PC1,tab8_NC$Cluster18S_1,tab8_NC$Cluster18S_4, tab8_NC$Cluster18S_11)


df_nc18s<-as.data.frame.array(df_nc18)
library(ggplot2)

ggplot(df_nc18s, aes(x =V1, y=log10(1+V2)))+geom_point()
ggplot(df_nc18s, aes(x =V1, y=log10(1+V3)))+geom_point()
ggplot(df_nc18s, aes(x =V1, y=log10(1+V4)))+geom_point()


cor.test(df_nc18s$V1,df_nc18s$V2, method="spearman")

cor.test(df_nc18s$V1,df_nc18s$V3, method="spearman")

cor.test(df_nc18s$V1,df_nc18s$V4, method="spearman")




#PCoA 

library(esquisse)
esquisser()
library(ggplot2)

#NC
library(readxl)
Pollution_NC_18S <- read_excel("Pollution_E.xlsx",sheet = "C_18S", col_types = c("text","numeric", "numeric"))
pcoa_nc_18S <- data.frame(pco_18Snc$vectors)

dat_pcoNC18 <- pcoa_nc_18S[, c("Axis.1", "Axis.2")]
dat_pcoNC18 <- cbind(dat_pcoNC18, Pollution_NC_18S$X.PC1)

ggplot(dat_pcoNC18) +
  aes(x = Axis.1, y = Axis.2, colour = `Pollution_NC_18S$X.PC1`) +
  geom_point(shape = "circle", 
             size = 4L) +
  scale_color_gradient(low = "#F3D43B", high = "#132B43") +
  labs(title = "Non colonized water 18S", 
       color = "Pollution (PC1)") +
  theme_minimal() +
  xlim(-0.4, 0.6) +
  ylim(-0.4, 0.5)

#EC
Pollution_EC_18S <- read_excel("Pollution_E.xlsx",sheet = "EC_18S", col_types = c("text","numeric", "numeric"))
pcoa_ec_18S <- data.frame(pco_18Sec$vectors)

dat_pcoEC18 <- pcoa_ec_18S[, c("Axis.1", "Axis.2")]
dat_pcoEC18 <- cbind(dat_pcoEC18, Pollution_EC_18S$X.PC1)

ggplot(dat_pcoEC18) +
  aes(x = Axis.1, y = Axis.2, colour = `Pollution_EC_18S$X.PC1`) +
  geom_point(shape = "circle", 
             size = 4L) +
  scale_color_gradient(low = "#F3D43B", high = "#132B43") +
  labs(title = "Colonized water 18S", 
       color = "Pollution (PC1)") +
  theme_minimal() +
  xlim(-0.4, 0.6) +
  ylim(-0.4, 0.5)




