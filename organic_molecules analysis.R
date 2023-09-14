library(readxl)
MICROPOL_all <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/MICROPOLLUANTS_moy.xlsx")
MICROPOL_EC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/MICROPOLLUANTS_moy.xlsx", sheet = "EC")
MICROPOL_NC <- read_excel("C:/Users/Penelope/Desktop/THESE/RESULTATS/Bioinfo/MICROPOLLUANTS_moy.xlsx", sheet = "NC")


library(mixOmics)

pca <- pca(MICROPOL_all[,c(3:149)],scale=TRUE, center=TRUE)

plot(pca) 
Colonisation<- MICROPOL_all$COLONISATION
biplot(pca, group=Colonisation,  ind.names = FALSE, legend = FALSE, labels = FALSE)



pcaM_nc <- pca(MICROPOL_NC[,c(2:148)],scale=TRUE, center=TRUE)

biplot(pcaM_nc, ind.names = FALSE, legend = FALSE, labels = FALSE)

pcaM_ec <- pca(MICROPOL_EC[,c(2:148)],scale=TRUE, center=TRUE)
biplot(pcaM_ec, ind.names = FALSE, legend = FALSE, labels = FALSE)

library(ade4)

proMP<-procuste(pcaM_nc$variates$X,pcaM_ec$variates$X)

plot(proMP)

#Monte-Carlo test
#Call: procuste.randtest(df1 = df1, df2 = df2, nrepet = nrepet)
#Observation: 0.0606304 
#Based on 999 replicates
#Simulated p-value: 0.985 
#Alternative hypothesis: greater 
#Std.Obs Expectation    Variance 
#-1.36885035  0.21914880  0.01341058 


randtest(proMP)

plsMP <- plsda(MICROPOL_all[,c(3:149)], MICROPOL_all$COLONISATION)

plotIndiv(plsMP,legend=T, ellipse = T, star=T)
plotVar(plsMP)
biplot(plsMP)
plotLoadings(plsMP, contrib='max', method='mean', legend = FALSE)

wilcox.test(MICROPOL_all$M61T49~MICROPOL_all$COLONISATION)
#Wilcoxon rank sum exact test
#data:  MICROPOL_all$M61T49 by MICROPOL_all$COLONISATION
#W = 211, p-value = 0.2469
#alternative hypothesis: true location shift is not equal to 0

#Prendre les coordonées des points pour EC et NC

pcaM_nc$variates->cbind(PCmp_nc+MICROPOL_NC$name)
pcaM_nc$variates->cbind(+MICROPOL_NC$name)
PCmp_nc <- cbind(pcaM_nc+MICROPOL_NC$name)
pcaM_ec$variates->PCmp_ec

write.table(PCmp_ec,"PCmp_ec.xls",row.names=T,col.names=F,quote=F,sep='\t')
write.table(PCmp_nc,"PCmp_nc.xls",row.names=T,col.names=F,quote=F,sep='\t')



#tester l'impact des micropolluants sur le 16S 
#début

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

pro_M <- procuste(pco_nc$vectors[,c(1,2)], pco_ec$vectors[,c(1,2)])
plot(pro_M)
randtest(pro_M)

#Monte-Carlo test
#Call: procuste.randtest(df1 = df1, df2 = df2, nrepet = nrepet)

#Observation: 0.3056967 

#Based on 999 replicates
#Simulated p-value: 0.229 
#Alternative hypothesis: greater 

#Std.Obs Expectation    Variance 
#0.731283661 0.241624980 0.007676461 

PCmp_ec1<-as.data.frame.array(PCmp_ec$X)
db_ecMP <- capscale(tab2_EC~PCmp_ec1$PC1+PCmp_ec1$PC2, dist="bray")

adonis2(tab2_EC~PCmp_ec1$PC1+PCmp_ec1$PC2, dist="bray")
#adonis2(formula = tab2_EC ~ PCmp_ec1$PC1 + PCmp_ec1$PC2, dist = "bray")
#Df SumOfSqs      R2      F Pr(>F)   
#PCmp_ec1$PC1  1   0.4376 0.05124 1.1580  0.167   
#PCmp_ec1$PC2  1   0.5446 0.06377 1.4412  0.007 **
#  Residual     20   7.5582 0.88499                 
#Total        22   8.5405 1.00000 

ordistep(db_ecMP)->oecMP
#Start: tab2_EC ~ PCmp_ec1$PC1 + PCmp_ec1$PC2 
#Df    AIC      F Pr(>F)   
# PCmp_ec1$PC1  1 50.797 1.1616   0.18   
# PCmp_ec1$PC2  1 51.098 1.4400   0.01 **
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Step: tab2_EC ~ PCmp_ec1$PC2 
#Df    AIC     F Pr(>F)  
#PCmp_ec1$PC2  1 50.311 1.429  0.015 *
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


oecMP

#Call: capscale(formula = tab2_EC ~ PCmp_ec1$PC2, distance = "bray")
#Inertia Proportion Rank
#Total         8.54162    1.00000     
#Constrained   0.54421    0.06371    1
#Unconstrained 7.99740    0.93629   21
#Inertia is squared Bray distance 
#Species scores projected from ‘tab2_EC’ 

#Eigenvalues for constrained axes:
#  CAP1 
#0.5442 

#Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8 
#0.8907 0.6605 0.6033 0.5658 0.5570 0.4632 0.4506 0.4232 
#(Showing 8 of 21 unconstrained eigenvalues)

plot(db_ecMP)

coorec_16mp<-as.data.frame.array(db_ecMP$CCA$v)
rownames(coorec[which(coorec_16mp$CAP1<0),])

write.table(coorec_16mp)
write.table(coorec_16mp,"coordEC_16sVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

#tester si les plus différentes si significative
df_ecPC2mp<-cbind(PCmp_ec1$PC2,tab2_EC$Cluster16S_39,tab2_EC$Cluster16S_21, tab2_EC$Cluster16S_142, PCmp_ec1$PC1)
df_ecPC2mp<-as.data.frame.array(df_ecPC2mp)


cor.test(df_ecPC2mp$V1,df_ecPC2mp$V2, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mp$V1 and df_ecPC2mp$V2
#S = 930.23, p-value = 0.007765
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 0.5404006 

#résultats sur pc1 
cor.test(df_ecPC2mp$V5,df_ecPC2mp$V2, method="spearman")


cor.test(df_ecPC2mp$V1,df_ecPC2mp$V3, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mp$V1 and df_ecPC2mp$V3
#S = 2757.9, p-value = 0.08905
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho -0.3626045 

cor.test(df_ecPC2mp$V1, df_ecPC2mp$V4, method = "spearman")
#	Spearman's rank correlation rho
#data:  df_ecPC2mp$V1 and df_ecPC2mp$V4
#S = 1752.5, p-value = 0.5417
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.1341588 

#NC

adonis2(tab2_NC~meta_NC_16$PC1_mp +meta_NC_16$PC2_mp, dist="bray")
#adonis2(formula = tab2_NC ~ meta_NC_16$PC1_mp + meta_NC_16$PC2_mp, dist = "bray")
#Df SumOfSqs      R2      F Pr(>F)  
#meta_NC_16$PC1_mp  1   0.4778 0.05463 1.2260  0.095 .
#meta_NC_16$PC2_mp  1   0.4742 0.05421 1.2167  0.117  
#Residual          20   7.7950 0.89116                
#Total             22   8.7470 1.00000 

db_ncMP <- capscale(tab2_NC~meta_NC_16$PC1_mp+meta_NC_16$PC2_mp, dist="bray")
db_ncMP
ordistep(db_ncMP)->oncMP
#Start: tab2_NC ~ PCmp_nc1$PC1 + PCmp_nc1$PC2 
#Df    AIC      F Pr(>F)  
# PCmp_nc1$PC2  1 51.552 1.2101  0.110  
# PCmp_nc1$PC1  1 51.569 1.2263  0.095 .
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Step: tab2_NC ~ PCmp_nc1$PC1 
#Df    AIC      F Pr(>F)
# PCmp_nc1$PC1  1 50.845 1.2141  0.125
#Step: tab2_NC ~ 1 


plot(oncMP)

coornc_16mp<-as.data.frame.array(db_ncMP$CCA$v)
write.table(coornc_16mp,"coord_16sVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')


#test de l'effet des micropolluants sur its 

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


library(ade4)
#Std.Obs  Expectation     Variance 
#-0.041201095  0.282547071  0.009808011 


pro_mpITS <- procuste(pco_nc$vectors[,c(1,2)], pco_ec$vectors[,c(1,2)])
plot(pro_mpITS)
randtest(pro_mpITS)

adonis2(tab7_EC~meta_Ec_ITS$PC1_mp + meta_Ec_ITS$PC2_mp, dist="bray")
#adonis2(formula = tab7_EC ~ meta_Ec_ITS$PC1_mp + meta_Ec_ITS$PC2_mp, dist = "bray")
#Df SumOfSqs      R2      F Pr(>F)  
#meta_Ec_ITS$PC1_mp  1   0.3946 0.05534 0.9483  0.681  
#meta_Ec_ITS$PC2_mp  1   0.4943 0.06932 1.1879  0.068 .
#Residual           15   6.2418 0.87534                
#Total              17   7.1307 1.00000   

db_ecMPITS <- capscale(tab7_EC~meta_Ec_ITS$PC1_mp + meta_Ec_ITS$PC2_mp, dist="bray")

ordistep(db_ecMPITS)->oecMPITS
#Start: tab7_EC ~ meta_Ec_ITS$PC1_mp + meta_Ec_ITS$PC2_mp 
#Df    AIC      F Pr(>F)  
#meta_Ec_ITS$PC1_mp  1 37.200 1.1020   0.24  
#meta_Ec_ITS$PC2_mp  1 37.303 1.1945   0.05 *
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Step: tab7_EC ~ meta_Ec_ITS$PC2_mp 
#Df    AIC      F Pr(>F)
#meta_Ec_ITS$PC2_mp  1 36.332 1.0391  0.335


oecMPITS
capscale(formula = tab7_EC ~ meta_Ec_ITS$PC2_mp, distance = "bray")

#Call: capscale(formula = tab7_EC ~ meta_Ec_ITS$PC2_mp, distance = "bray")
#Inertia Proportion Rank
#Total         7.13145    1.00000     
#Constrained   0.43490    0.06098    1
#Unconstrained 6.69655    0.93902   16
#Inertia is squared Bray distance 
#Species scores projected from ‘tab7_EC’ 

#Eigenvalues for constrained axes:
#  CAP1 
#0.4349 
#Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8   MDS9  MDS10  MDS11  MDS12  MDS13  MDS14  MDS15  MDS16 
#0.6604 0.6271 0.5725 0.5555 0.4883 0.4731 0.4533 0.4154 0.3889 0.3787 0.3454 0.3209 0.2842 0.2655 0.2546 0.2128 

plot(db_ecMPITS)

coorec_itsmp<-as.data.frame.array(db_ecMPITS$CCA$v)

write.table(coorec_itsmp,"coordEC_itssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

#tester si les plus différentes si significative
df_ecPC2mpITS<-cbind(meta_Ec_ITS$PC2_mp,tab7_EC$ClusterITS_8,tab7_EC$ClusterITS_6, tab7_EC$ClusterITS_3, tab7_EC$ClusterITS_44, tab7_EC$ClusterITS_14, meta_Ec_ITS$PC1_mp)
df_ecPC2mpITS<-as.data.frame.array(df_ecPC2mpITS)


cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V2, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mpITS$V1 and df_ecPC2mpITS$V2
#S = 988.29, p-value = 0.9375
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho -0.01990396 

cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V3, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mpITS$V1 and df_ecPC2mpITS$V3
#S = 1152.1, p-value = 0.4527
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho -0.1889757 


cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V4, method="spearman")
#	Spearman's rank correlation rho
#data:  df_ecPC2mpITS$V1 and df_ecPC2mpITS$V4
#S = 1048, p-value = 0.7482
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho -0.08152735 

cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V5, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mpITS$V1 and df_ecPC2mpITS$V5
#S = 499.83, p-value = 0.04174
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 0.4841805 

cor.test(df_ecPC2mpITS$V7,df_ecPC2mpITS$V5, method="spearman")


cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V6, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mpITS$V1 and df_ecPC2mpITS$V6
#S = 583.96, p-value = 0.1025
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 0.3973597

#NC

adonis2(tab7_NC~meta_Nc_ITS$PC1_mp + meta_Nc_ITS$PC2_mp, dist="bray")
#adonis2(formula = tab7_NC ~ meta_Nc_ITS$PC1_mp + meta_Nc_ITS$PC2_mp, dist = "bray")
#Df SumOfSqs      R2      F Pr(>F)  
#meta_Nc_ITS$PC1_mp  1   0.5549 0.08017 1.3918  0.025 *
#meta_Nc_ITS$PC2_mp  1   0.3861 0.05578 0.9684  0.577  
#Residual           15   5.9802 0.86405                
#Total              17   6.9211 1.00000 

db_ncMPITS <- capscale(tab7_NC~meta_Nc_ITS$PC1_mp + meta_Nc_ITS$PC2_mp, dist="bray")

ordistep(db_ncMPITS)->oncMPITS
#Start: tab7_NC ~ meta_Nc_ITS$PC1_mp + meta_Nc_ITS$PC2_mp 
#Df    AIC      F Pr(>F)  
# meta_Nc_ITS$PC2_mp  1 36.228 0.9701  0.615  
# meta_Nc_ITS$PC1_mp  1 36.669 1.3662  0.040 *
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Step: tab7_NC ~ meta_Nc_ITS$PC1_mp 
#Df    AIC      F Pr(>F)  
#meta_Nc_ITS$PC1_mp  1 35.711 1.3745  0.035 *
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


oncMPITS
capscale(formula = tab7_EC ~ meta_Ec_ITS$PC2_mp, distance = "bray")

#Call: capscale(formula = tab7_NC ~ meta_Nc_ITS$PC1_mp, distance = "bray")
#nertia Proportion Rank
#Total         6.88955    1.00000     
#Constrained   0.54505    0.07911    1
#Unconstrained 6.34450    0.92089   16
#Inertia is squared Bray distance 
#Species scores projected from ‘tab7_NC’ 

#Eigenvalues for constrained axes:
#  CAP1 
#0.545 

#Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8   MDS9  MDS10  MDS11  MDS12  MDS13  MDS14  MDS15  MDS16 
#0.8264 0.6930 0.5641 0.4884 0.4798 0.4431 0.4416 0.4013 0.3865 0.3417 0.3009 0.2615 0.2246 0.2122 0.1546 0.1249 

plot(db_ncMPITS)

coornc_itsmp<-as.data.frame.array(db_ncMPITS$CCA$v)

write.table(coornc_itsmp,"coordNC_itssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

#tester si les plus différentes si significative
df_ncPC1mpITS<-cbind(meta_Nc_ITS$PC1_mp,tab7_NC$ClusterITS_3,tab7_NC$ClusterITS_30, tab7_NC$ClusterITS_61, tab7_NC$ClusterITS_163, tab7_NC$ClusterITS_107, meta_Nc_ITS$PC2_mp)
df_ncPC1mpITS<-as.data.frame.array(df_ncPC1mpITS)


cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V2, method="spearman")
#Spearman's rank correlation rho
#data:  df_ncPC1mpITS$V1 and df_ncPC1mpITS$V2
#S = 584, p-value = 0.1036
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.3973168 

cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V3, method="spearman")
#Spearman's rank correlation rho
#data:  df_ncPC1mpITS$V1 and df_ncPC1mpITS$V3
#S = 478.26, p-value = 0.03198
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.5064353 

#test PC2
cor.test(df_ncPC1mpITS$V7,df_ncPC1mpITS$V3, method="spearman")
#data:  df_ncPC1mpITS$V7 and df_ncPC1mpITS$V3
#S = 1015.2, p-value = 0.8509
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho -0.04771218 

cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V4, method="spearman")
#	Spearman's rank correlation rho
#data:  df_ncPC1mpITS$V1 and df_ncPC1mpITS$V4
#S = 602.08, p-value = 0.1213
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.3786552

cor.test(df_ecPC2mpITS$V1,df_ecPC2mpITS$V5, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mpITS$V1 and df_ecPC2mpITS$V5
#S = 499.83, p-value = 0.04174
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 0.4841805 

#test PC2
cor.test(df_ncPC1mpITS$V7,df_ncPC1mpITS$V5, method="spearman")
#data:  df_ncPC1mpITS$V7 and df_ncPC1mpITS$V5
#S = 549.82, p-value = 0.07297
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.4325905 


cor.test(df_ncPC1mpITS$V1,df_ncPC1mpITS$V6, method="spearman")
#Spearman's rank correlation rho
#data:  df_ncPC1mpITS$V1 and df_ncPC1mpITS$V6
#S = 1362.7, p-value = 0.09432
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho -0.4062821 

#18S
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


#EC

adonis2(tab8_EC~met_18SEC$PC1_mp + met_18SEC$PC2_mp, dist="bray")
#Number of permutations: 999
#adonis2(formula = tab8_EC ~ met_18SEC$PC1_mp + met_18SEC$PC2_mp, dist = "bray")
#Df SumOfSqs      R2      F Pr(>F)
#met_18SEC$PC1_mp  1   0.4722 0.08325 1.0011  0.489
#met_18SEC$PC2_mp  1   0.4835 0.08523 1.0250  0.343
#Residual         10   4.7169 0.83153              
#Total            12   5.6726 1.00000

db_ecMP18S <- capscale(tab8_EC~met_18SEC$PC1_mp + met_18SEC$PC2_mp, dist="bray")

ordistep(db_ecMP18S)->oecMP18S
#                   Df   AIC      F Pr(>F)
#met_18SEC$PC1_mp  1 24.36 1.0004  0.435
#met_18SEC$PC2_mp  1 24.39 1.0260  0.365

#Step: tab8_EC ~ met_18SEC$PC2_mp 
#Df    AIC      F Pr(>F)
#met_18SEC$PC2_mp  1 23.519 1.0257    0.4
#Step: tab8_EC ~ 1 


oecMP18S
capscale(formula = tab8_EC ~ met_18SEC$PC1_mp, distance = "bray")

#Inertia Proportion Rank
#Total         5.67095    1.00000     
#Constrained   0.47163    0.08317    1
#Unconstrained 5.19932    0.91683   11
#Inertia is squared Bray distance 
#Species scores projected from ‘tab8_EC’ 
#Eigenvalues for constrained axes:
#  CAP1 
#0.4716 

#Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8   MDS9  MDS10  MDS11 
#0.6497 0.5502 0.5429 0.5322 0.4860 0.4610 0.4472 0.4380 0.3891 0.3627 0.3402 

plot(db_ecMP18S)

coorec_18smp<-as.data.frame.array(db_ecMP18S$CCA$v)

write.table(coorec_18smp,"coordEC_18ssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

#au niveau du 18S, les eaux colonisées ne sont pas corrélées avec PC1 ni PC2  
df_ecPC2mp18S<-cbind(met_18SEC$PC2_mp,tab8_EC$Cluster18S_970,tab8_EC$Cluster18S_16, tab8_EC$Cluster18S_10, tab8_EC$Cluster18S_17)
df_ecPC2mp18S<-as.data.frame.array(df_ecPC2mp18S)


cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V2, method="spearman")
#	Spearman's rank correlation rho
#data:  df_ecPC2mp18S$V1 and df_ecPC2mp18S$V2
#S = 348.1, p-value = 0.8873
#alternative hypothesis: true rho is not equal to 0
#sample estimates:  rho 0.04367853 

cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V3, method="spearman")
#	Spearman's rank correlation rho
#data:  df_ecPC2mp18S$V1 and df_ecPC2mp18S$V3
#S = 364, p-value = 1
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0

cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V4, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mp18S$V1 and df_ecPC2mp18S$V4
#S = 460.91, p-value = 0.3793
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#       rho -0.2662498 

cor.test(df_ecPC2mp18S$V1,df_ecPC2mp18S$V5, method="spearman")
#Spearman's rank correlation rho
#data:  df_ecPC2mp18S$V1 and df_ecPC2mp18S$V5
#S = 195.5, p-value = 0.1112
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.46291 

#NC

adonis2(tab8_NC~met_18SNC$PC1_mp + met_18SNC$PC2_mp, dist="bray")
#adonis2(formula = tab8_NC ~ met_18SNC$PC1_mp + met_18SNC$PC2_mp, dist = "bray")
#Df SumOfSqs      R2      F Pr(>F)  
#met_18SNC$PC1_mp  1   0.5955 0.10581 1.3184  0.065 .
#met_18SNC$PC2_mp  1   0.5156 0.09161 1.1415  0.167  
#Residual         10   4.5171 0.80258                
#Total            12   5.6282 1.00000 

db_ncMP18S <- capscale(tab8_NC~met_18SNC$PC1_mp + met_18SNC$PC2_mp, dist="bray")

ordistep(db_ncMP18S)->oncMP18S
#Start: tab8_NC ~ met_18SNC$PC1_mp + met_18SNC$PC2_mp 
#Df    AIC      F Pr(>F)  
#met_18SNC$PC2_mp  1 23.962 1.1406  0.170  
#met_18SNC$PC1_mp  1 24.194 1.3408  0.075 .
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Step: tab8_NC ~ met_18SNC$PC1_mp 
#Df    AIC     F Pr(>F)  
#met_18SNC$PC1_mp  1 23.417 1.302   0.07 .
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

oncMP18S

#Inertia Proportion Rank
#Total          5.6265     1.0000     
#Constrained    0.5955     0.1058    1
#Unconstrained  5.0310     0.8942   11
#Inertia is squared Bray distance 
#Species scores projected from ‘tab8_NC’ 

#Eigenvalues for constrained axes:
#  CAP1 0.5955 

#Eigenvalues for unconstrained axes:
#  MDS1   MDS2   MDS3   MDS4   MDS5   MDS6   MDS7   MDS8   MDS9  MDS10  MDS11 
#0.8434 0.6351 0.6309 0.5540 0.4978 0.4724 0.4400 0.4004 0.2558 0.1887 0.1124

plot(db_ncMP18S)

coornc_18smp<-as.data.frame.array(db_ncMP18S$CCA$v)

write.table(coornc_18smp,"coordNC_18ssVSmp.xls",row.names=T,col.names=F,quote=F,sep='\t')

#au niveau du 18S, les eaux colonisées ne sont pas corrélées avec PC1 ni PC2  
df_ncPC1mp18S<-cbind(met_18SNC$PC1_mp,tab8_NC$Cluster18S_24,tab8_NC$Cluster18S_11, tab8_NC$Cluster18S_1)
df_ncPC1mp18S<-as.data.frame.array(df_ncPC1mp18S)

cor.test(df_ncPC1mp18S$V1,df_ncPC1mp18S$V2, method="spearman")
#	Spearman's rank correlation rho
#data:  df_ncPC1mp18S$V1 and df_ncPC1mp18S$V2
#S = 194.94, p-value = 0.1098
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.4644434

cor.test(df_ncPC1mp18S$V1,df_ncPC1mp18S$V3, method="spearman")
#	Spearman's rank correlation rho
#data:  df_ncPC1mp18S$V1 and df_ncPC1mp18S$V3
#S = 448.25, p-value = 0.4467
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho -0.231455 

cor.test(df_ncPC1mp18S$V1,df_ncPC1mp18S$V4, method="spearman")
#	Spearman's rank correlation rho
#data:  df_ncPC1mp18S$V1 and df_ncPC1mp18S$V4
#S = 166.13, p-value = 0.05484
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 0.5435874 