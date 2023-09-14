library(readxl)
library(mixOmics)

MICROPOL_all <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/MICROPOLLUANTS_moy.xlsx")
MICROPOL_EC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/MICROPOLLUANTS_moy.xlsx", sheet = "EC")
MICROPOL_NC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/MICROPOLLUANTS_moy.xlsx", sheet = "NC")

pca <- pca(MICROPOL_all[,c(3:149)],scale=TRUE, center=TRUE)

plot(pca) 
Colonisation<- MICROPOL_all$COLONISATION
biplot(pca, group=Colonisation,  ind.names = FALSE, legend = FALSE, labels = FALSE)



pcaM_nc <- pca(MICROPOL_NC[,c(2:148)],scale=TRUE, center=TRUE)

biplot(pcaM_nc, ind.names = FALSE, legend = FALSE, labels = FALSE)

pcaM_ec <- pca(MICROPOL_EC[,c(2:148)],scale=TRUE, center=TRUE)
biplot(pcaM_ec, ind.names = FALSE, legend = FALSE, labels = FALSE)

#procruste analysis
library(ade4)
proMP<-procuste(pcaM_nc$variates$X,pcaM_ec$variates$X)
plot(proMP)
randtest(proMP)

plsMP <- plsda(MICROPOL_all[,c(3:149)], MICROPOL_all$COLONISATION)

plotIndiv(plsMP,legend=T, ellipse = T, star=T)
plotVar(plsMP)
biplot(plsMP)
plotLoadings(plsMP, contrib='max', method='mean', legend = FALSE)

#test the most discriminate parameters
wilcox.test(MICROPOL_all$M61T49~MICROPOL_all$COLONISATION)

#coordinate extraction for further analysis 
pcaM_nc$variates->cbind(PCmp_nc+MICROPOL_NC$name)
pcaM_nc$variates->cbind(+MICROPOL_NC$name)
PCmp_nc <- cbind(pcaM_nc+MICROPOL_NC$name)
pcaM_ec$variates->PCmp_ec

write.table(PCmp_ec,"PCmp_ec.xls",row.names=T,col.names=F,quote=F,sep='\t')
write.table(PCmp_nc,"PCmp_nc.xls",row.names=T,col.names=F,quote=F,sep='\t')



#impact on 16S communities

table_EC16 <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/16S/DATA_BRUT/EC_NC_moy/table_E_moy16S.xlsx", 
                         sheet = "EC")
meta_EC_16 <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/16S/DATA_BRUT/EC_NC_moy/meta_E_moy16S.xlsx", 
                         sheet = "EC")

tab2_EC<-table_EC16[,-1]

table_NC16 <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/16S/DATA_BRUT/EC_NC_moy/table_E_moy16S.xlsx", 
                         sheet = "NC")
meta_NC_16 <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/16S/DATA_BRUT/EC_NC_moy/meta_E_moy16S.xlsx", sheet = "NC", col_types = c("text", "text", "text", "text", "text", "text", "text",  "text", "numeric", "numeric", "numeric"))

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

#procruste analysis
pro_M <- procuste(pco_nc$vectors[,c(1,2)], pco_ec$vectors[,c(1,2)])
plot(pro_M)
randtest(pro_M)

#test on colonized water samples 
PCmp_ec1<-as.data.frame.array(PCmp_ec$X)
db_ecMP <- capscale(tab2_EC~PCmp_ec1$PC1+PCmp_ec1$PC2, dist="bray")

adonis2(tab2_EC~PCmp_ec1$PC1+PCmp_ec1$PC2, dist="bray")


ordistep(db_ecMP)->oecMP
oecMP

plot(db_ecMP)
coorec_16mp<-as.data.frame.array(db_ecMP$CCA$v)
rownames(coorec[which(coorec_16mp$CAP1<0),])

write.table(coorec_16mp)
write.table(coorec_16mp,"coordEC_16sVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

df_ecPC2mp<-cbind(PCmp_ec1$PC2,tab2_EC$Cluster16S_39,tab2_EC$Cluster16S_21, tab2_EC$Cluster16S_142, PCmp_ec1$PC1)
df_ecPC2mp<-as.data.frame.array(df_ecPC2mp)


cor.test(df_ecPC2mp$V1,df_ecPC2mp$V2, method="spearman")

cor.test(df_ecPC2mp$V5,df_ecPC2mp$V2, method="spearman")

cor.test(df_ecPC2mp$V1,df_ecPC2mp$V3, method="spearman")

cor.test(df_ecPC2mp$V1, df_ecPC2mp$V4, method = "spearman")


#test on uncolonized waters samples

adonis2(tab2_NC~meta_NC_16$PC1_mp +meta_NC_16$PC2_mp, dist="bray")

db_ncMP <- capscale(tab2_NC~meta_NC_16$PC1_mp+meta_NC_16$PC2_mp, dist="bray")
db_ncMP
ordistep(db_ncMP)->oncMP

plot(oncMP)

coornc_16mp<-as.data.frame.array(db_ncMP$CCA$v)
write.table(coornc_16mp,"coord_16sVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')


#impact on ITS communities

table_Ec_ITS <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/ITS/data_brut/EC_NC_moyITS/table_OTU_E_moyITS.xlsx", sheet = "EC_T")
meta_Ec_ITS <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/ITS/data_brut/EC_NC_moyITS/metadata_E_ITS.xlsx", sheet = "EC")
table_Nc_ITS <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/ITS/data_brut/EC_NC_moyITS/table_OTU_E_moyITS.xlsx", sheet = "NC_T")
meta_Nc_ITS <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/ITS/data_brut/EC_NC_moyITS/metadata_E_ITS.xlsx", sheet = "NC")

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
pco_ncits<-pcoa(vegdist(tab7_NC,"bray"))
biplot.pcoa(pco_ncits)

pco_ecits<-pcoa(vegdist(tab7_EC,"bray"))
biplot.pcoa(pco_ecits)


#procruste analysis
pro_mpITS <- procuste(pco_nc$vectors[,c(1,2)], pco_ec$vectors[,c(1,2)])
plot(pro_mpITS)
randtest(pro_mpITS)

#test on colonized water samples 
adonis2(tab7_EC~meta_Ec_ITS$PC1_mp + meta_Ec_ITS$PC2_mp, dist="bray")

db_ecMPITS <- capscale(tab7_EC~meta_Ec_ITS$PC1_mp + meta_Ec_ITS$PC2_mp, dist="bray")
ordistep(db_ecMPITS)->oecMPITS

oecMPITS
capscale(formula = tab7_EC ~ meta_Ec_ITS$PC2_mp, distance = "bray")


plot(db_ecMPITS)

coorec_itsmp<-as.data.frame.array(db_ecMPITS$CCA$v)

write.table(coorec_itsmp,"coordEC_itssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

df_ecPC2mpITS<-cbind(meta_Ec_ITS$PC2_mp,tab7_EC$ClusterITS_8,tab7_EC$ClusterITS_6, tab7_EC$ClusterITS_3, tab7_EC$ClusterITS_44, tab7_EC$ClusterITS_14, meta_Ec_ITS$PC1_mp)
df_ecPC2mpITS<-as.data.frame.array(df_ecPC2mpITS)


cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V2, method="spearman")


cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V3, method="spearman")


cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V4, method="spearman")

cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V5, method="spearman")

cor.test(df_ecPC2mpITS$V7,df_ecPC2mpITS$V5, method="spearman")

cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V6, method="spearman")



#test on uncolonized water samples 

adonis2(tab7_NC~meta_Nc_ITS$PC1_mp + meta_Nc_ITS$PC2_mp, dist="bray")
db_ncMPITS <- capscale(tab7_NC~meta_Nc_ITS$PC1_mp + meta_Nc_ITS$PC2_mp, dist="bray")
ordistep(db_ncMPITS)->oncMPITS
oncMPITS

capscale(formula = tab7_NC ~ meta_Nc_ITS$PC2_mp, distance = "bray")

plot(db_ncMPITS)

coornc_itsmp<-as.data.frame.array(db_ncMPITS$CCA$v)

write.table(coornc_itsmp,"coordNC_itssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

df_ncPC1mpITS<-cbind(meta_Nc_ITS$PC1_mp,tab7_NC$ClusterITS_3,tab7_NC$ClusterITS_30, tab7_NC$ClusterITS_61, tab7_NC$ClusterITS_163, tab7_NC$ClusterITS_107, meta_Nc_ITS$PC2_mp)
df_ncPC1mpITS<-as.data.frame.array(df_ncPC1mpITS)

#test PC1
cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V2, method="spearman")

cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V3, method="spearman")

#test PC2

cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V4, method="spearman")

cor.test(df_ncPC1mpITS$V7,df_ncPC1mpITS$V5, method="spearman")

cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V6, method="spearman")


#impact on ITS communities

library(readxl)

table_EC_18S <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/18S/DATA_BRUT/table_OTU_moy18S.xlsx", sheet = "EC")
table_NC_18S <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/18S/DATA_BRUT/table_OTU_moy18S.xlsx", sheet = "NC")

met_18SEC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/18S/DATA_BRUT/metadata_moy18S.xlsx", sheet = "ec")
met_18SNC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/MICROBIO/meta/NGS/VF_analyses/18S/DATA_BRUT/metadata_moy18S.xlsx", sheet = "nc")

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


#test on colonized water samples 

adonis2(tab8_EC~met_18SEC$PC1_mp + met_18SEC$PC2_mp, dist="bray")

db_ecMP18S <- capscale(tab8_EC~met_18SEC$PC1_mp + met_18SEC$PC2_mp, dist="bray")

ordistep(db_ecMP18S)->oecMP18S

oecMP18S
capscale(formula = tab8_EC ~ met_18SEC$PC1_mp, distance = "bray")

plot(db_ecMP18S)

coorec_18smp<-as.data.frame.array(db_ecMP18S$CCA$v)

write.table(coorec_18smp,"coordEC_18ssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

df_ecPC2mp18S<-cbind(met_18SEC$PC2_mp,tab8_EC$Cluster18S_970,tab8_EC$Cluster18S_16, tab8_EC$Cluster18S_10, tab8_EC$Cluster18S_17)
df_ecPC2mp18S<-as.data.frame.array(df_ecPC2mp18S)


cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V2, method="spearman")


cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V3, method="spearman")

cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V4, method="spearman")


cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V5, method="spearman")


#test on uncolonized water samples 

adonis2(tab8_NC~met_18SNC$PC1_mp + met_18SNC$PC2_mp, dist="bray")

db_ncMP18S <- capscale(tab8_NC~met_18SNC$PC1_mp + met_18SNC$PC2_mp, dist="bray")

ordistep(db_ncMP18S)->oncMP18S

oncMP18S

plot(db_ncMP18S)

coornc_18smp<-as.data.frame.array(db_ncMP18S$CCA$v)

write.table(coornc_18smp,"coordNC_18ssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

df_ncPC1mp18S<-cbind(met_18SNC$PC1_mp,tab8_NC$Cluster18S_24,tab8_NC$Cluster18S_11, tab8_NC$Cluster18S_1)
df_ncPC1mp18S<-as.data.frame.array(df_ncPC1mp18S)

cor.test(df_ncPC1mp18S$V1,df_ncPC1mp18S$V2, method="spearman")

cor.test(df_ncPC1mp18S$V1,df_ncPC1mp18S$V3, method="spearman")

cor.test(df_ncPC1mp18S$V1,df_ncPC1mp18S$V4, method="spearman")
