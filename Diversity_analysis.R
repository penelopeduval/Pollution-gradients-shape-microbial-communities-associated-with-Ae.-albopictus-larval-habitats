library(ggplot2)
library(car)
library(nlme)

#ITS
shan <- read_delim("ITS_shannon_all.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)

ggplot(data = shan, )

ggplot(shan, aes(x=condition, y=shan, color= condition, shape = condition)) +
  geom_violin(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(), size=2)+
  scale_shape_manual(values=c(17, 16, 15))+
  theme_classic() + facet_wrap(~facet) +
  scale_color_manual(values=c("#3399FF", "#00CC33", "#FF3333"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5) +
  ylim(0, 2.5)

lm_shan <- lm((shan$shan)^2~shan$condition)

qqPlot(lm_shan)
plot(lm_shan) #fit parfaitement à la normalité et hétéroscédasticité
leveneTest(lm_shan)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value Pr(>F)
#group  2  0.8301 0.4418

Anova(lm_shan)
#Response: (shan$shan)^2
#Sum Sq Df F value  Pr(>F)  
#shan$condition 10.245  2  3.3632 0.04246 *
#  Residuals      77.680 51   

library(emmeans)

em_shan <- emmeans(lm_shan, pairwise~condition)
em_shan
#$contrasts
#contrast estimate    SE df t.ratio p.value
#EC - L     -0.787 0.411 51  -1.914  0.1452
#EC - NC     0.230 0.411 51   0.559  0.8422
#L - NC      1.017 0.411 51   2.473  0.0436

#16S

shan2 <- read_delim("16S_shannon_all.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
shan2$shan <- shan2$`indice de Shannon`

ggplot(shan2, aes(x=Condition, y=shan, color= Condition, shape = Condition)) +
  geom_violin(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(), size=2)+
  scale_shape_manual(values=c(17, 16, 15))+
  theme_classic() + facet_wrap(~facet) +
  scale_color_manual(values=c("#3399FF", "#00CC33", "#FF3333"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5)  +
  ylim(0, 5)

lm_shan2 <- lm(shan2$shan~shan2$Condition)

qqPlot(lm_shan2)
plot(lm_shan2) #fit pas et problème d'hétéroscédasticité
leveneTest(lm_shan2)

kruskal.test(shan2$shan~shan2$Condition)
#Kruskal-Wallis chi-squared = 23.343, df = 2, p-value = 8.532e-06

pairwise.wilcox.test(shan2$shan, shan2$Condition,
                     p.adjust.method = "BH")

#Pairwise comparisons using Wilcoxon rank sum test with continuity correction 

#data:  shan2$shan and shan2$Condition 

#EC     L     
#L  0.0012 -     
#  NC 0.5276 0.0005

#18S 

shan3 <- read_delim("18S_shannon_all.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)

ggplot(shan3, aes(x=condition, y=shan, color= condition, shape = condition)) +
  geom_violin(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(), size=2)+
  scale_shape_manual(values=c(17, 15))+
  theme_classic() + facet_wrap(~facet) +
  scale_color_manual(values=c("#3399FF", "#FF3333"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5) +
  ylim(0, 4)

lm_shan3 <- lm((shan3$shan)^2~shan3$condition)

qqPlot(lm_shan3)
leveneTest(lm_shan3)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value Pr(>F)
#group  1  0.0123 0.9128


Anova(lm_shan3)
#Response: (shan3$shan)^2
#Sum Sq Df F value Pr(>F)
#shan3$condition   1.606  1  0.2807 0.6011
#Residuals       137.296 24
