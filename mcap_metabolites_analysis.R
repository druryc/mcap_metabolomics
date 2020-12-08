setwd("~/CRD_GBS/metabolites/data")
library(tidyverse);library(readxl);library(cowplot);library(RColorBrewer);library(janitor);library(viridis);library(vegan);library(car);library(magicfor);library(grid);library(factoextra);library(egg)
########################################## COLOR SCHEME ################################################################### #####
#list<-read_excel("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(molecular_family)%>%
#     distinct()%>%arrange(molecular_family)%>%mutate(molecular_family=str_to_title(molecular_family))
#gg_color_hue <- function(n) {hues = seq(15, 375, length = n + 1)
#     hcl(h = hues, l = 65, c = 100)[1:n]}
#n = nrow(list);cols = gg_color_hue(n);as.list(list)
cols <- c("Amino Acid"="#F8766D",
          "Betaine"="#ED813E",
          "Bile Acid"="#DE8C00",
          "Carnitine"="#CD9600",
          "Chlorophyll"="#B79F00",
          "Diterpenoid"="#9DA700",
          "Eicosanoid"="#7CAE00",
          "Endocannabinoid"="#49B500",
          "Fatty Acid"="#00BA38",
          "Indole"="#00BE67",
          "Lactone"="#00C08B",
          "Microbial Natural Product"="#00C1A9",
          "Monoacylglyceride"="#00BFC4",
          "Nucleotide"="#00BBDC",
          "Peptide"="#00B4F0",
          "Phosphatidic Acids"="#00A9FF",
          "Phosphocholine"="#FF6C91",
          "Phosphoethanolamine"="#9F8CFF",
          "Phosphoserine"="#C77CFF",
          "Prostaglandin"="#E36EF6",
          "Steroid"="#F564E3",
          "Triterpenoid"="#FF61CC",
          "Unknown"="lightgray",
          "Xanthin"="#619CFF",
          "Other"="lightgray")



############################## DIVERSITY (SF) ############################################################################# #####
raw<-read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%filter(time_point=="T1")
div<-as.matrix(raw[,17:3575])
list<-read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%filter(time_point=="T1")%>%select(colony,bleaching_history_phenotype)

shannon<-as.data.frame(diversity(div, index = "shannon", MARGIN = 1, base = exp(1)));names(shannon)[1]<-"shannon"
richness<-as.data.frame(div)%>%mutate(richness=rowSums(.!=0))%>%select(richness)
metrics<-bind_cols(list,shannon,richness)%>%mutate(evenness=shannon/log(richness))

leveneTest(shannon~bleaching_history_phenotype,data=metrics);t.test(shannon~bleaching_history_phenotype,data=metrics)
leveneTest(richness~bleaching_history_phenotype,data=metrics);t.test(richness~bleaching_history_phenotype,data=metrics)
leveneTest(evenness~bleaching_history_phenotype,data=metrics);t.test(evenness~bleaching_history_phenotype,data=metrics)

metrics$bleaching_history_phenotype <- factor(metrics$bleaching_history_phenotype,labels=c("B","NB"))


a<-ggplot(metrics)+
     geom_boxplot(aes(bleaching_history_phenotype,shannon,fill=bleaching_history_phenotype))+
     geom_jitter(aes(bleaching_history_phenotype,shannon),width=0.02,alpha=0.2)+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     theme_classic(base_size=8)+
     ylab("Shannon's Diversity Index")+
     xlab("")+
     theme(legend.position="none")+
     annotate("text",label="p=0.010",x=2,y=5.05,size=2,hjust=1,fontface = 'italic')

b<-ggplot(metrics)+geom_boxplot(aes(bleaching_history_phenotype,richness,fill=bleaching_history_phenotype))+
     geom_jitter(aes(bleaching_history_phenotype,richness),width=0.02,alpha=0.2)+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     theme_classic(base_size=8)+
     ylab("Richness")+
     xlab("")+
     theme(legend.position="none")+
     annotate("text",label="p=0.021",x=2,y=3240,size=2,hjust=1,fontface = 'italic')

c<-ggplot(metrics)+geom_boxplot(aes(bleaching_history_phenotype,evenness,fill=bleaching_history_phenotype))+
     geom_jitter(aes(bleaching_history_phenotype,evenness),width=0.02,alpha=0.2)+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     theme_classic(base_size=8)+
     ylab("Evenness")+
     xlab("")+
     annotate("text",label="p=0.013",x=2,y=.62,size=2,hjust=1,fontface = 'italic')+
     theme(legend.position="none",legend.key.size = unit(0.5,"cm"))

metrics$bleaching_history_phenotype <- factor(metrics$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
lg<-get_legend(ggplot(metrics)+geom_boxplot(aes(bleaching_history_phenotype,evenness,fill=bleaching_history_phenotype))+
                    scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
                    theme_classic(base_size=8))

quartz(width=4,height=2.5)
plot_grid(a,b,c,lg,nrow=1,align="v",rel_widths=c(1,1,1,1),labels=c('A','B','C'),label_size=8,label_x=c(0.36,0.42,0.4))






############################## RANDOM FOREST VARIABLE IMPORTANCES ######################################################### #####
data<-read_xlsx("bucket_table.xlsx",sheet="T1_T2")
coral_RF <- randomForest(as.factor(data$ATTRIBUTE_Phenotype) ~., data[,11:2307], importance = TRUE, ntree = 1000)
VI<-varImp(coral_RF)
#write.table(VI, "variable_importance.xlsx",quote=FALSE,row.names=FALSE,sep="\t")

############################## PCOAs (F1) ################################################################################# #####
#using loadings from JMP, see github.com/druryc/mcap_metabolomics/JMP
all<-read_xlsx("pcoa_coordinates.xlsx",sheet="Lab_All")%>%clean_names()
all$bleaching_history_phenotype <- factor(all$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
all$colony <- factor(all$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
all$colony <- as.factor(all$colony)
data.env<-all%>%select(sample_number,bleaching_history_phenotype,colony,prin1,prin2)
coord<-all%>%select(prin1,prin2)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ bleaching_history_phenotype*colony, data = data.env) #PERMANOVA 

a<-ggplot(all)+geom_point(aes(prin1,prin2,fill=bleaching_history_phenotype,color=colony),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-80,80),breaks=seq(-80,80,20))+
     scale_x_continuous(limits=c(-80,80),breaks=seq(-80,80,20))+
     ylab("Component 2 (12.4%)")+
     xlab("Component 1 (15.6%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="left")+
     annotate("text",x=-80,y=-80,hjust=0,vjust=1,label="~Phenotype p<0.001, ~Genotype p<0.001",size=2,fontface = 'italic');a

val<-read_xlsx("pcoa_coordinates.xlsx",sheet="Val_All")%>%clean_names()%>%mutate(prin1=prin1*(-1))
val$bleaching_history_phenotype <- factor(val$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
val$colony <- factor(val$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
val$colony <- as.factor(val$colony)
data.env<-val%>%select(sample_number,bleaching_history_phenotype,prin1,prin2)
coord<-val%>%select(prin1,prin2)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ bleaching_history_phenotype, data = data.env) #PERMANOVA 

b<-ggplot(val)+geom_point(aes(prin1,prin2,fill=bleaching_history_phenotype),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-80,80),breaks=seq(-80,80,20))+
     scale_x_continuous(limits=c(-80,80),breaks=seq(-80,80,20))+
     ylab("Component 2 (15.6%)")+
     xlab("Component 1 (26.6%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.position="none")+
     annotate("text",x=-80,y=-80,hjust=0,vjust=1,label="~Phenotype p=0.006",size=2,fontface = 'italic');b

VIP<-left_join(transform(read_xlsx("variable_importance.xlsx",sheet="All")%>%clean_names()%>%dplyr::rename(molecule=1), molecule=reorder(molecule, -mean_decrease_accuracy))%>%
     head(500)%>%rownames_to_column("n")%>%mutate(n=as.numeric(n))%>%
     separate(molecule,into=c("cluster_id","other"),sep="_",remove=FALSE)%>%
     separate(cluster_id,into=c("trash","cluster_id"),sep=1)%>%select(-other,-trash)%>%
     mutate(cluster_id=as.numeric(cluster_id)),read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(cluster_id,molecular_family)%>%mutate(cluster_id=as.numeric(cluster_id)),by="cluster_id")%>%
     mutate(molecular_family=str_to_title(molecular_family))

c<-ggplot(VIP)+geom_bar(aes(n,mean_decrease_accuracy,fill=molecular_family),stat="identity",width=1)+
     scale_fill_manual(values=cols,name="")+
     theme_classic(base_size=8)+
     ylab("Variable Importance")+
     scale_x_continuous(limits=c(0,500),breaks=seq(0,500,100))+
     theme(axis.title.x=element_blank(),
           legend.position="left",
          legend.key.size = unit(0.15,"cm"),
          legend.text=element_text(size=6))+
     annotate("segment",x=1,xend=30,y=9.08,yend=9.08)+annotate("text",x=32,y=9.08,label="Betaine Lipid",hjust=0,size=2)+
     annotate("segment",x=20,xend=30,y=7.2,yend=8)+annotate("text",x=32,y=8,label="Glutamine-Phenylalanine",hjust=0,size=2)+
     annotate("segment",x=124,xend=124,y=3.45,yend=8)+annotate("text",x=126,y=8,label="1-Arachidoyl-2-hydroxy-sn-glycero-3-phosphocholine",hjust=0,size=2)+
     annotate("segment",x=132,xend=132,y=3.30,yend=7.3)+annotate("text",x=134,y=7.3,label="Octadecatrienoic acid, 3-(hexopyranosyloxy)-2-hydroxypropyl ester",hjust=0,size=2)+
     annotate("segment",x=149,xend=149,y=3.03,yend=6.6)+annotate("text",x=151,y=6.6,label="Octadecatrienoic acid, 3-(hexopyranosyloxy)-2-hydroxypropyl ester",hjust=0,size=2)+
     annotate("segment",x=161,xend=161,y=2.72,yend=5.9)+annotate("text",x=163,y=5.9,label="PAF",hjust=0,size=2)+
     annotate("segment",x=171,xend=171,y=2.58,yend=5.2)+annotate("text",x=173,y=5.2,label="Lyso PC",hjust=0,size=2)+
     annotate("segment",x=198,xend=198,y=2.32,yend=4.5)+annotate("text",x=200,y=4.5,label="1-(1Z-Hexadecenyl)-sn-glycero-3-phosphocholine",hjust=0,size=2)+
     annotate("segment",x=228,xend=228,y=2.07,yend=3.7)+annotate("text",x=230,y=3.7,label="Lyso-PAF C-18",hjust=0,size=2)+
     annotate("segment",x=230,xend=230,y=2.07,yend=3.0)+annotate("text",x=232,y=3.0,label="Tryptamine",hjust=0,size=2)+
     annotate("segment",x=242,xend=242,y=2,yend=2.3)+annotate("text",x=244,y=2.3,label="Monolinolein",hjust=0,size=2)+
     annotate("segment",x=276,xend=300,y=1.80,yend=2.3)+annotate("text",x=302,y=2.3,label="Pheophorbide A",hjust=0,size=2);c
     
d<-VIP%>%
     ggplot(.)+geom_boxplot(aes(reorder(molecular_family,-n,FUN=median),n,fill=molecular_family),color=NA,outlier.shape=NA)+
     geom_point(aes(reorder(molecular_family,-n,FUN=median),n),size=0.5,alpha=0.3)+
     theme_classic(base_size=8)+
     coord_flip()+
     xlab("")+
     scale_fill_manual(values=cols)+
     ylab("Molecule (Sorted by VI)")+
     theme(axis.title.y=element_blank(),
           legend.position="none",
           legend.key.size = unit(0.4,"cm"),
           legend.text=element_text(size=6),
           panel.background = element_rect(fill = "transparent",colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.grid.major.y=element_line(color="lightgray",linetype="dotted"),
           axis.text.y=element_text(size=6))+
     guides(fill=guide_legend(nrow=2));d

plots <- align_plots(a,c,d, align = 'v', axis = 'l')
z<-plot_grid(plots[[1]],b,align="h",axis="lr",nrow=1,rel_widths=c(2,1.5),labels=c("A","B"),label_size=8,label_x=c(0.4,0.12))
quartz(width=(183/25.4),height=5.5)
x<-plot_grid(c,d,nrow=2,align="v",axis="l",labels=c("C","D"))
plot_grid(z,plots[[2]],plots[[3]],align="v",axis="r",nrow=3,ncol=1,rel_heights=c(1.8,1,1),labels=c("","C","D"),label_size=8,label_x=c(0.23,0.23))



############################## W/IN NETWORK INITIAL DIFFERENCES (F2) ) #################################################### #####
setwd("~/CRD_GBS/metabolites/data")
all<-read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%filter(time_point=="T1")%>%select(-filename,-feature,-mac,-well_id,-sample_id,-experiment,-sum,-sample_number,-plate_number,-nubbin_number,-tank)%>%gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))
val<-read_xlsx("bucket_table.xlsx",sheet="Val")%>%clean_names%>%select(-filename,-feature,-mac,-well_id,-sample_id,-experiment,-sum,-sample_number,-plate_number,-nubbin_number,-tank)%>%gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))

a1<-val%>%filter(stringr::str_detect(molecule, "x384_"))
a1$bleaching_history_phenotype <- factor(a1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=a1)
wilcox.test(aoc~bleaching_history_phenotype,data=a1)
my_y_title <- expression(paste( italic("m/z"), "490.368"))
a<-ggplot(a1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.title.x=element_blank())+
     scale_y_continuous(limits=c(16,23),breaks=seq(16,23,1))+
     annotate("text",x=2,y=16,label="p<0.001",size=1.5,fontface = 'italic',vjust=0);a

b1<-val%>%filter(stringr::str_detect(molecule, "x384_"))
b1$bleaching_history_phenotype <- factor(b1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=b1)
t.test(aoc~bleaching_history_phenotype,data=b1)
b<-ggplot(b1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(16,23),breaks=seq(16,23,1))+
     annotate("text",x=2,y=16,label="p<0.001",size=1.5,fontface = 'italic',vjust=0);b

c1<-val%>%filter(stringr::str_detect(molecule, "x2475_"))
c1$bleaching_history_phenotype <- factor(c1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=c1)
wilcox.test(aoc~bleaching_history_phenotype,data=c1)
my_y_title <- expression(paste( italic("m/z"), "484.327"))

c<-ggplot(c1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.title.x=element_blank())+
     scale_y_continuous(limits=c(-1,16),breaks=seq(0,15,5))+
     annotate("text",x=2,y=-1,label="p<0.001",size=1.5,fontface = 'italic',vjust=0);c

quartz(w=4,h=2)
plot_grid(a,c,labels=c("A","B"),label_size=8,label_x=c(0.12,0.12))

d1<-val%>%filter(stringr::str_detect(molecule, "x2475_"))
d1$bleaching_history_phenotype <- factor(d1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=d1)
wilcox.test(aoc~bleaching_history_phenotype,data=d1)
d<-ggplot(d1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(-1,16),breaks=seq(0,15,5))+
     annotate("text",x=2,y=-1,label="p=0.002",size=1.5,fontface = 'italic',vjust=0);d


e1<-all%>%filter(stringr::str_detect(molecule, "x3078_"))
e1$bleaching_history_phenotype <- factor(e1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=e1)
t.test(aoc~bleaching_history_phenotype,data=e1)
my_y_title <- expression(paste( italic("m/z"), "518.404"))
e<-ggplot(e1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(13,17),breaks=seq(13,17,1))+
     annotate("text",x=2,y=13,label="p=0.004",size=1.5,fontface = 'italic',vjust=0);e

f1<-val%>%filter(stringr::str_detect(molecule, "x3078_"))
f1$bleaching_history_phenotype <- factor(f1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=f1)
t.test(aoc~bleaching_history_phenotype,data=f1)
my_y_title <- expression(paste( italic("m/z"), "518.404"))
f<-ggplot(f1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(9,16),breaks=seq(9,16,2))+
     annotate("text",x=2,y=9,label="p<0.001",size=1.5,fontface = 'italic',vjust=0);f

g1<-all%>%filter(stringr::str_detect(molecule, "x775_"))
g1$bleaching_history_phenotype <- factor(e1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=g1)
wilcox.test(aoc~bleaching_history_phenotype,data=g1)
my_y_title <- expression(paste( italic("m/z"), "516.384"))
g<-ggplot(g1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(12,16),breaks=seq(12,16,1))+
     annotate("text",x=2,y=12,label="p=0.848",size=1.5,fontface = 'italic',vjust=0);g

h1<-val%>%filter(stringr::str_detect(molecule, "x775_"))
h1$bleaching_history_phenotype <- factor(h1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=h1)
t.test(aoc~bleaching_history_phenotype,data=h1)
my_y_title <- expression(paste( italic("m/z"), "516.384"))
h<-ggplot(h1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(12,16),breaks=seq(12,16,1))+
     annotate("text",x=2,y=12,label="p=0.008",size=1.5,fontface = 'italic',vjust=0);h

m1<-all%>%filter(stringr::str_detect(molecule, "x382_"))
m1$bleaching_history_phenotype <- factor(m1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=m1)
t.test(aoc~bleaching_history_phenotype,data=m1)
my_y_title <- expression(paste( italic("m/z"), "562.368"))
m<-ggplot(m1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(15.9,22),breaks=seq(16,22,1))+
     annotate("text",x=2,y=16,label="p<0.001",size=1.5,fontface = 'italic',vjust=0);m

n1<-val%>%filter(stringr::str_detect(molecule, "x382_"))
n1$bleaching_history_phenotype <- factor(n1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=n1)
t.test(aoc~bleaching_history_phenotype,data=n1)
my_y_title <- expression(paste( italic("m/z"), "562.368"))
n<-ggplot(n1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(18,22),breaks=seq(18,22,1))+
     annotate("text",x=2,y=18,label="p=0.002",size=1.5,fontface = 'italic',vjust=0);n

o1<-all%>%filter(stringr::str_detect(molecule, "x8356_"))
o1$bleaching_history_phenotype <- factor(o1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=o1)
t.test(aoc~bleaching_history_phenotype,data=o1)
my_y_title <- expression(paste( italic("m/z"), "822.584"))
o<-ggplot(m1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
    scale_y_continuous(limits=c(15.9,22),breaks=seq(16,22,1))+
     annotate("text",x=2,y=16,label="p=0.002",size=1.5,fontface = 'italic',vjust=0);o

p1<-val%>%filter(stringr::str_detect(molecule, "x8356_"))
p1$bleaching_history_phenotype <- factor(p1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(aoc~bleaching_history_phenotype,data=p1)
t.test(aoc~bleaching_history_phenotype,data=p1)
my_y_title <- expression(paste( italic("m/z"), "822.584"))
p<-ggplot(p1)+geom_boxplot(aes(bleaching_history_phenotype,aoc,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,aoc),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(my_y_title)+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+
     scale_y_continuous(limits=c(0,14.1),breaks=seq(0,14,1))+
     annotate("text",x=2,y=0,label="p=0.116",size=1.5,fontface = 'italic',vjust=0);p

#EDGE FIGURES
all_edge<-read_excel("edge.xlsx",sheet="T1")%>%clean_names()
val_edge<-read_excel("edge.xlsx",sheet="Validation")%>%clean_names()

i1<-all_edge%>%select(bleaching_history_phenotype,c4h8_gain)
i1$bleaching_history_phenotype <- factor(i1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(c4h8_gain~bleaching_history_phenotype,data=i1)
t.test(c4h8_gain~bleaching_history_phenotype,data=i1)
i<-ggplot(i1)+geom_boxplot(aes(bleaching_history_phenotype,c4h8_gain,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,c4h8_gain),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(expression("C"[4]*"H"[8]*" Addition"))+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.title.x=element_blank())+
     scale_y_continuous(limits=c(8.5,9.5),breaks=seq(8.5,9.5,0.2))+
     annotate("text",x=2,y=8.5,label="p<0.001",size=1.5,fontface = 'italic',vjust=0);i

j1<-val_edge%>%select(bleaching_history_phenotype,c4h8_gain)%>%mutate(c4h8_gain=log10(c4h8_gain))
j1$bleaching_history_phenotype <- factor(j1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(c4h8_gain~bleaching_history_phenotype,data=j1)
t.test(c4h8_gain~bleaching_history_phenotype,data=j1)
j<-ggplot(j1)+geom_boxplot(aes(bleaching_history_phenotype,c4h8_gain,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,c4h8_gain),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab(expression("C"[4]*"H"[8]*" Addition"))+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.title.x=element_blank())+
     scale_y_continuous(limits=c(6,7.2),breaks=seq(6,7.2,0.2))+
     annotate("text",x=2,y=6,label="p=0.015",size=1.5,fontface = 'italic',vjust=0);j

k1<-all_edge%>%select(bleaching_history_phenotype,h2_gain)
k1$bleaching_history_phenotype <- factor(k1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(h2_gain~bleaching_history_phenotype,data=k1)
t.test(h2_gain~bleaching_history_phenotype,data=k1)
k<-ggplot(k1)+geom_boxplot(aes(bleaching_history_phenotype,h2_gain,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,h2_gain),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab("Saturation")+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.title.x=element_blank())+
     scale_y_continuous(limits=c(9,9.6),breaks=seq(9,9.6,0.1))+
     annotate("text",x=2,y=9,label="p<0.001",size=1.5,fontface = 'italic',vjust=0);k

l1<-val_edge%>%select(bleaching_history_phenotype,h2_gain)%>%mutate(h2_gain=log10(h2_gain))
l1$bleaching_history_phenotype <- factor(l1$bleaching_history_phenotype,labels=c("B","NB"))
leveneTest(h2_gain~bleaching_history_phenotype,data=l1)
t.test(h2_gain~bleaching_history_phenotype,data=l1)
l<-ggplot(l1)+geom_boxplot(aes(bleaching_history_phenotype,h2_gain,fill=bleaching_history_phenotype),outlier.alpha = 0.1)+geom_jitter(aes(bleaching_history_phenotype,h2_gain),width=0.02,alpha=0.1)+
     theme_classic(base_size=6)+ylab("Saturation")+scale_fill_manual(values=c("white","darkgray"))+
     theme(legend.position="none",axis.title.x=element_blank())+
     scale_y_continuous(limits=c(8,9.6),breaks=seq(8,9.6,0.5))+
     annotate("text",x=2,y=8,label="p=0.029",size=1.5,fontface = 'italic',vjust=0);l

quartz(width=1.2,height=5)
plot_grid(c,a,m,o,e,g,k,ncol=1,labels=c("B","C","D","E","F","G","H"),align="v",axis="l",label_size=8,label_x=c(0.23,0.23,0.23,0.23,0.23,0.23,0.23,0.26))

############################## BETAINE INITIAL DIFFERENCES (F2) ) ######################################################### ##### 
quartz(width=1.2,height=3)
plot_grid(d,b,ncol=1,align="v",axis="l",label_size=8,label_x=c(0.23))

############################## DGTS INITIAL DIFFERENCES ################################################################### #####
library(tidyverse);library(magicfor);library(readxl);library(egg)
all<-read_xlsx("bucket_table.xlsx",sheet="DGTS")%>%select(-contains('...2'))%>%
     select(BleachingHistoryPhenotype,contains('m\\z'))%>%
     rename('m\\z522.372'=3)%>%
     select(-'m\\z474.368',-'m\\z496.356',-'m\\z498.373',-'m\\z524.389')%>%
     gather(molecule,auc,-BleachingHistoryPhenotype)%>%
     mutate(auc=log(auc))%>%mutate(auc=case_when(auc=="-Inf"~0,TRUE ~ as.numeric(auc)))%>%
     mutate(molecule=as.factor(molecule))

all$BleachingHistoryPhenotype <- factor(all$BleachingHistoryPhenotype,labels=c("B","NB"))
list<-unique(all$molecule)
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[1]),exact=FALSE)$p.value
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[2]),exact=FALSE)$p.value
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[3]),exact=FALSE)$p.value
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[4]),exact=FALSE)$p.value
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[5]),exact=FALSE)$p.value
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[6]),exact=FALSE)$p.value
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[7]),exact=FALSE)$p.value
wilcox.test(auc~BleachingHistoryPhenotype,data=all%>%filter(molecule==list[8]),exact=FALSE)$p.value


quartz(w=5,h=2)
plot<-ggplot(all)+geom_boxplot(aes(BleachingHistoryPhenotype,auc,fill=BleachingHistoryPhenotype),outlier.alpha = 0.1)+
     geom_jitter(aes(BleachingHistoryPhenotype,auc),width=0.02,alpha=0.1)+
     facet_wrap(~molecule,nrow=1)+
     theme_classic(base_size=8)+scale_fill_manual(values=c("white","darkgray"))+
     guides(fill=FALSE)+
     xlab("Historical Bleaching Phenotype")+
     scale_y_continuous(limits=c(0,20));plot


tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
     
     gb <- ggplot_build(p)
     lay <- gb$layout$layout
     tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
     p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                   vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

my_tag <- c("p=0.095","p=0.323","p=0.017","p=0.107","p=0.571","p=0.291","p=0.976","p=0.044")
tag_facet(plot, 
          x = 1, y = 19, 
          vjust = 0, hjust = -0.1,
          open = "", close = "",
          fontface = 'italic',
          size = 2,
          tag_pool = my_tag)




############################## PCOAs COMPARTMENTS (F3) #################################################################### #####
#shared
data<-read_xlsx("bucket_table.xlsx",sheet="Shared")%>%clean_names()%>%select(-filename,-attribute_coral_species,-attribute_phenotype,-attribute_replicate_id,-experiment,-attribute_sample_type)%>%
     gather(molecule,aoc,-bleaching_history_phenotype,-colony,-temperature_treatment)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0, TRUE ~ as.numeric(aoc)))%>%
     spread(molecule,aoc)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
data$colony <- as.factor(data$colony)

pca<-prcomp(data[,4:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:3],axes$data)
coord<-plotdata%>%select(x,y)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ bleaching_history_phenotype*colony, data = plotdata) #PERMANOVA 

axes$labels$x
axes$labels$y

shared<-ggplot(plotdata)+geom_point(aes(x,y,fill=bleaching_history_phenotype,color=colony),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-140,140,25),breaks=c(-140,-100,-50,0,50,100,140))+
     scale_x_continuous(limits=c(-140,140,25),breaks=c(-140,-100,-50,0,50,100,140))+     
     xlab("Component 1 (11.5%)")+
     ylab("Component 2 (7.9%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="left",
           legend.spacing.y = unit(0.0, "cm"),
           legend.margin= margin(0.2,0,0,0, unit="cm"),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.background = element_rect(fill = "transparent",colour = NA))+
     annotate("text",x=-140,y=-140,hjust=0,vjust=1,label="~Phenotype p<0.001,~Genotype p<0.001",size=2,fontface = 'italic');shared

#coral
data<-read_xlsx("bucket_table.xlsx",sheet="Coral_Only")%>%clean_names()%>%select(-filename,-attribute_coral_species,-attribute_phenotype,-attribute_replicate_id,-experiment,-attribute_sample_type)%>%
     gather(molecule,aoc,-bleaching_history_phenotype,-colony,-temperature_treatment)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0, TRUE ~ as.numeric(aoc)))%>%
     spread(molecule,aoc)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
data$colony <- as.factor(data$colony)

pca<-prcomp(data[,4:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:3],axes$data)
coord<-plotdata%>%select(x,y)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ bleaching_history_phenotype*colony, data = plotdata) #PERMANOVA 

axes$labels$x
axes$labels$y

coral<-ggplot(plotdata)+geom_point(aes(x,y,fill=bleaching_history_phenotype,color=colony),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-140,140,25),breaks=c(-140,-100,-50,0,50,100,140))+
     scale_x_continuous(limits=c(-140,140,25),breaks=c(-140,-100,-50,0,50,100,140))+
     xlab("Component 1 (13.2%)")+
     ylab("Component 2 (8.5%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none")+
     annotate("text",x=-140,y=-140,hjust=0,vjust=1,label="~Phenotype p<0.001,~Genotype p<0.001",size=2,fontface = 'italic');coral

#symbiont
data<-read_xlsx("bucket_table.xlsx",sheet="Symbiont_Only")%>%clean_names()%>%select(-filename,-attribute_coral_species,-attribute_phenotype,-attribute_replicate_id,-experiment,-attribute_sample_type)%>%
     gather(molecule,aoc,-bleaching_history_phenotype,-colony,-temperature_treatment)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0, TRUE ~ as.numeric(aoc)))%>%
     spread(molecule,aoc)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
data$colony <- as.factor(data$colony)

pca<-prcomp(data[,4:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:3],axes$data)
coord<-plotdata%>%select(x,y)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ bleaching_history_phenotype*colony, data = plotdata) #PERMANOVA 

axes$labels$x
axes$labels$y
symb<-ggplot(plotdata)+geom_point(aes(x,y,fill=bleaching_history_phenotype,color=colony),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-140,140,25),breaks=c(-140,-100,-50,0,50,100,140))+
     scale_x_continuous(limits=c(-140,140,25),breaks=c(-140,-100,-50,0,50,100,140))+
     xlab("Component 1 (13.9%)")+
     ylab("Component 2 (9.8%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none")+
     annotate("text",x=-140,y=-140,hjust=0,vjust=1,label="~Phenotype p<0.001,~Genotype p<0.001",size=2,fontface = 'italic');symb

#shared 1445
#coral 660
#symbiont 1783

library(VennDiagram)
a<-ggdraw(draw.pairwise.venn((660+1447), (1783+1447), 1447,
                             c("H", "S"),
                             fill = c("white", "white"),
                             alpha = c(0.8, 0.8),
                             fontfamily='sans',
                             cat.fontfamil='sans',
                             cex=0.4,0.4,
                             cat.cex=0.5,lwd=1,scaled=TRUE,
                             cat.prompts=TRUE,
                             cat.pos=c(230,-240),
                             cat.dist=c(0.04,0.04)))+
     theme_classic(base_size=7)+
     theme(axis.line.x=element_blank(),
           axis.line.y=element_blank());a

inset<-ggdraw()+draw_plot(a,0.65,0.74,0.3,0.3)+draw_plot(shared)


quartz(w=7.2,h=2)
plot_grid(inset,coral,symb,nrow=1,rel_widths=c(1.4,1,1),labels=c("A","B","C"),label_size=8,label_x=c(0.41,0.2,0.2))
detach(package:VennDiagram)
############################## PCOAs KNOWN (SUPPLEMENTAL) ################################################################# #####
library(tidyverse);library(readxl);library(cowplot);library(RColorBrewer);library(janitor);library(viridis);library(vegan)
known<-read_xlsx("pcoa_coordinates.xlsx",sheet="Lab_Known")%>%clean_names()%>%mutate(prin1=prin1*(-1))
val<-read_xlsx("pcoa_coordinates.xlsx",sheet="Validation_Known")%>%clean_names()%>%mutate(prin1=prin1*(-1))

known$bleaching_history_phenotype <- factor(known$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
known$colony <- factor(known$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
known$colony <- as.factor(known$colony)
data.env<-known%>%select(sample_number,bleaching_history_phenotype,colony,prin1,prin2)
coord<-known%>%select(prin1,prin2)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ bleaching_history_phenotype+colony, data = data.env) #PERMANOVA 

a<-ggplot(known)+geom_point(aes(prin1,prin2,fill=bleaching_history_phenotype,color=colony),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-15,15),breaks=seq(-15,15,3))+
     scale_x_continuous(limits=c(-15,15),breaks=seq(-15,15,3))+
     ylab("Component 2 (15.4%)")+
     xlab("Component 1 (19.2%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="left")+
     annotate("text",x=-15,y=-15,hjust=0,vjust=1,label="~Phenotype p<0.001, ~Genotype p<0.001",size=2,fontface = 'italic');a

val$bleaching_history_phenotype <- factor(val$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
val$colony <- factor(val$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
val$colony <- as.factor(val$colony)
data.env<-val%>%select(sample_number,bleaching_history_phenotype,prin1,prin2)
coord<-val%>%select(prin1,prin2)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ bleaching_history_phenotype, data = data.env) #PERMANOVA 

b<-ggplot(val)+geom_point(aes(prin1,prin2,fill=bleaching_history_phenotype),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(-15,15),breaks=seq(-15,15,3))+
     scale_x_continuous(limits=c(-15,15),breaks=seq(-15,15,3))+
     ylab("Component 2 (15.9%)")+
     xlab("Component 1 (42.7%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.position="none")+
     annotate("text",x=-15,y=-15,hjust=0,vjust=1,label="~Phenotype p=0.004",size=2,fontface = 'italic');b

VIP<-transform(read_xlsx("variable_importance.xlsx",sheet="Known")%>%clean_names()%>%rename(molecule=1), molecule=reorder(molecule, -mean_decrease_accuracy))%>%
     head(500)%>%rownames_to_column("n")%>%mutate(n=as.numeric(n))%>%mutate(molecular_family=str_to_title(molecular_family))%>%
     mutate(molecular_family=case_when(molecular_family=="Contaminant"~"Other",
                                       molecular_family=="Unknown"~"Other",
                                       TRUE~as.character(molecular_family)))


VIP$molecular_family<-factor(VIP$molecular_family, levels=c("Amino Acid","Carnitine","Chlorophyll","Eicosanoid","Fatty Acid","Monoacylglyceride","Nucleotide","Peptide","Phosphocholine","Prostaglandin","Steroid","Xanthin","Other"))

c<-ggplot(VIP)+geom_bar(aes(n,mean_decrease_accuracy,fill=molecular_family),stat="identity",width=1)+
     scale_fill_manual(values=cols,name="")+
     theme_classic(base_size=8)+
     ylab("Variable Importance")+
     scale_x_continuous(limits=c(0,139),breaks=seq(0,139,25))+
     theme(axis.title.x=element_blank(),
           legend.position="left",
           legend.key.size = unit(0.15,"cm"),
           legend.text=element_text(size=6));c

d<-VIP%>%
     ggplot(.)+geom_boxplot(aes(reorder(molecular_family,-n,FUN=median),n,fill=molecular_family),color=NA,outlier.shape=NA)+
     scale_fill_manual(values=cols,name="")+
     geom_point(aes(reorder(molecular_family,-n,FUN=median),n),size=0.5,alpha=0.3)+
     theme_classic(base_size=8)+
     scale_y_continuous(limits=c(0,139),breaks=seq(0,139,25))+
     coord_flip()+
     xlab("")+
     ylab("Molecule (Sorted by VI)")+
     theme(axis.title.y=element_blank(),
           legend.position="none",
           legend.key.size = unit(0.4,"cm"),
           legend.text=element_text(size=6),
           panel.background = element_rect(fill = "transparent",colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.grid.major.y=element_line(color="lightgray",linetype="dotted"),
           axis.text.y=element_text(size=5))+
     guides(fill=guide_legend(nrow=2));d

quartz(width=(160/25.4),height=4.5)

plots <- align_plots(a,c,d, align = 'v', axis = 'l')
x<-plot_grid(plots[[1]],b,nrow=1,rel_widths=c(2,1.5),align="h",axis="lr",labels=c("A","B"),label_size=8,label_x=c(0.39,0.14))
y<-plot_grid(x,plots[[2]],plots[[3]],align="v",axis="r",nrow=3,rel_heights=c(2,1,1),labels=c("","C","D"),label_size=8,label_x=c(0.22,0.22));y




############################## OVERALL RESPONSE TO TEMP STRESS (F5) ####################################################### #####
mf<-read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(cluster_id,mean_decrease_accuracy,molecular_family)
raw<-read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names()%>%select(-filename,-feature,-mac,-well_id,-sample_id,-experiment,-sum,-sample_number,-plate_number,-nubbin_number,-tank)%>%
     gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
     separate(molecule,into=c("cluster_id","trash"),remove=FALSE)%>%select(-trash)%>%separate(cluster_id,into=c("x","cluster_id"),sep=1)%>%select(-x)%>%
     mutate(cluster_id=as.numeric(cluster_id))%>%mutate(aoc=(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0, TRUE ~ as.numeric(aoc)))
data<-left_join(raw,mf,by="cluster_id")%>%group_by(colony,time_point,cluster_id)%>%mutate(mean_aoc=mean(aoc))%>%mutate(count=length(which(aoc==0)))%>%
     mutate(aoc=case_when(aoc==0&count==1 ~ mean_aoc,
                          TRUE ~ as.numeric(aoc)))%>%ungroup()

correction<-data%>%filter(temperature_treatment=="Control")%>%filter(time_point=="T1")%>%select(cluster_id,colony,aoc)%>%rename(initial=aoc)
value<-left_join(data%>%filter(temperature_treatment=="Constant"&time_point=="T1"),correction,by=c("colony","cluster_id"))%>%mutate(diff=aoc-initial)%>%select(cluster_id,colony,diff)
inta<-left_join(data%>%filter((temperature_treatment=="Constant")&time_point=="T2"),value,by=c("colony","cluster_id"))%>%mutate(aoc_stand=aoc-diff)%>%select(-diff)%>%arrange(colony,cluster_id)%>%rename(constant_aoc_stand=aoc_stand)
intb<-left_join(data%>%filter((temperature_treatment=="Control")&time_point=="T2"),value,by=c("colony","cluster_id"))%>%mutate(aoc_stand=aoc)%>%select(-diff)%>%arrange(colony,cluster_id)%>%rename(control_aoc_stand=aoc_stand)%>%select(control_aoc_stand)

out<-bind_cols(inta,intb)%>%
     mutate(rawchange=constant_aoc_stand-control_aoc_stand)%>%select(bleaching_history_phenotype,cluster_id,molecule,molecular_family,rawchange)%>%
     group_by(cluster_id,bleaching_history_phenotype)%>%mutate(change=mean(rawchange))%>%distinct(bleaching_history_phenotype,cluster_id,change,.keep_all=TRUE)%>%select(-rawchange)%>%
     spread(bleaching_history_phenotype,change)%>%rename(change_b=Bleached,change_nb='Not Bleached')


corr<-left_join(out,mf,by="cluster_id")%>%select(-molecular_family.y)%>%rename(molecular_family=molecular_family.x)%>%arrange(desc(mean_decrease_accuracy))%>%rownames_to_column()%>%
     mutate(rowname=as.numeric(rowname))%>%mutate(top=case_when(rowname<=100 ~"Top 100 VI",
                                                                rowname>100 ~ "Other"))%>%select(-rowname)%>%arrange(mean_decrease_accuracy)
corr$top<-factor(corr$top,levels=c("Top 100 VI","Other"))
lb<-quantile(corr$change_nb,0.25)-(1.5*IQR(corr$change_b)) #set bounds for classification
rb<-quantile(corr$change_nb,0.75)+(1.5*IQR(corr$change_b)) #set bounds for classification
tb<-quantile(corr$change_nb,0.75)+(1.5*IQR(corr$change_nb)) #set bounds for classification
bb<-quantile(corr$change_nb,0.25)-(1.5*IQR(corr$change_nb)) #set bounds for classification


########################################## FULL ANALYSIS ################################################################## #####
#slope is significantly different from 1
model<-lm(change_nb~change_b,data=corr);summary(model)
model0<-lm(change_nb~1+offset(change_b),data=corr)
anova(model,model0)
comp_data<-corr%>%select(cluster_id,change_b,change_nb)%>%gather(pheno,change,-cluster_id)
car::leveneTest(change~pheno,data=comp_data) #significant difference in variance
a<-var((comp_data%>%filter(pheno=="change_b"))$change)
b<-var((comp_data%>%filter(pheno=="change_nb"))$change)
(a-b)/a #32.5% larger variance in bleached corals
wilcox.test(change~pheno,data=comp_data) #no significant difference in mean change

########################################## PLOTS ########################################################################## #####
anno1=expression("R"^2~"= 0.2407; p<0.001")
a<-ggplot(corr)+
     geom_abline(slope=1,color="lightgray",linetype="dotted")+
     annotate("rect", xmin=-200000000, xmax=lb, ymin=bb, ymax=tb, alpha=0.2, fill="gray")+
     annotate("rect", xmin=rb, xmax=200000000, ymin=bb, ymax=tb, alpha=0.2, fill="gray")+
     annotate("rect", xmin=lb, xmax=rb, ymin=tb, ymax=200000000, alpha=0.2, fill="gray")+
     annotate("rect", xmin=lb, xmax=rb, ymin=-200000000, ymax=bb, alpha=0.2, fill="gray")+
     annotate("rect", xmin=lb, xmax=rb, ymin=bb, ymax=tb, fill="darkgray")+
     theme_classic(base_size=8)+
     scale_alpha_manual(values=c(0.5,0.2),name="")+
     scale_color_manual(values=c("orange","black"),name="")+
     scale_size_manual(values=c(1,0.5),name="")+
     geom_vline(xintercept=lb,linetype="dotted")+     #geom_vline(xintercept=lib,linetype="dotted")+
     geom_vline(xintercept=rb,linetype="dotted")+     #geom_vline(xintercept=rib,linetype="dotted")+
     geom_hline(yintercept=tb,linetype="dotted")+     #geom_hline(yintercept=tib,linetype="dotted")+
     geom_hline(yintercept=bb,linetype="dotted")+     #geom_hline(yintercept=bib,linetype="dotted")+
     ylab("Non-bleached NHR")+
     xlab("Bleached NHR")+
     scale_y_continuous(limits=c(-200000000,200000000))+
     scale_x_continuous(limits=c(-200000000,200000000))+
     geom_point(aes(change_b,change_nb,color=top,alpha=top,size=top))+
     theme(legend.position=c(0.8,0.2),legend.key.size = unit(0.5,"cm"))+
     geom_smooth(aes(change_b,change_nb),method="lm")+
     annotate("text",x=100000000,y=70000000,label=anno1,size=2,fontface = 'italic');a

x<-ggplot(corr%>%filter(molecular_family!="unknown"&molecular_family!="contaminant"&molecular_family!="ontaminant"))+
     theme_classic(base_size=8)+
     geom_vline(xintercept=0,linetype="dotted")+     geom_hline(yintercept=0,linetype="dotted")+
     coord_cartesian(xlim = c(-5,5),ylim=c(-5,5))+
     ylab("Non-bleached NHR")+
     xlab("Bleached NHR")+
     #geom_point(aes(change_b,change_nb,color=molecular_family),size=1)+
     stat_ellipse(aes(change_b,change_nb,color=molecular_family))+
     theme(legend.key.size = unit(0.5,"cm"));x

########################################## SUMMARY VALUES ################################################################# #####
#summary values
noresponse<-nrow(corr%>%filter(change_b>lb,change_b<rb,change_nb>bb,change_nb<tb,top=="Top 100 VI"));noresponse
bleached<-nrow(corr%>%filter((change_b<lb|change_b>rb),change_nb>bb,change_nb<tb,top=="Top 100 VI"));bleached
nonbleached<-nrow(corr%>%filter((change_nb<bb|change_nb>tb),change_b>lb,change_b<rb,top=="Top 100 VI"));nonbleached
cooperative1<-nrow(corr%>%filter(change_b>rb,change_nb>tb,top=="Top 100 VI"));cooperative1
cooperative2<-nrow(corr%>%filter(change_b<lb,change_nb<bb,top=="Top 100 VI"));cooperative2
antagonistic1<-nrow(corr%>%filter(change_b>rb,change_nb<bb,top=="Top 100 VI"));antagonistic1
antagonistic2<-nrow(corr%>%filter(change_b<lb,change_nb>tb,top=="Top 100 VI"));antagonistic2
ant<-sum(antagonistic1,antagonistic2);ant
coop<-sum(cooperative1,cooperative2);coop

#summary of all metabolites
noresponse_a<-nrow(corr%>%filter(change_b>lb,change_b<rb,change_nb>bb,change_nb<tb));noresponse_a
bleached_a<-nrow(corr%>%filter((change_b<lb|change_b>rb),change_nb>bb,change_nb<tb));bleached_a
nonbleached_a<-nrow(corr%>%filter((change_nb<bb|change_nb>tb),change_b>lb,change_b<rb));nonbleached_a
cooperative1_a<-nrow(corr%>%filter(change_b>rb,change_nb>tb));cooperative1_a
cooperative2_a<-nrow(corr%>%filter(change_b<lb,change_nb<bb));cooperative2_a
antagonistic1_a<-nrow(corr%>%filter(change_b>rb,change_nb<bb));antagonistic1_a
antagonistic2_a<-nrow(corr%>%filter(change_b<lb,change_nb>tb));antagonistic2_a
ant_a<-sum(antagonistic1_a,antagonistic2_a)
coop_a<-sum(cooperative1_a,cooperative2_a)
frame<-c("No Response","Bleached Specific","Non-bleached Specific", "Cooperative","Antagonistic","No Response","Bleached Specific","Non-bleached Specific", "Cooperative","Antagonistic")
n<-c(noresponse,bleached,nonbleached,ant,coop,noresponse_a,bleached_a,nonbleached_a,ant_a,coop_a)
dataset<-c("Top VI","Top VI","Top VI","Top VI","Top VI","All","All","All","All","All")
figdata<-data.frame(frame,n,dataset)%>%mutate(prop=case_when(dataset=="Top VI"~n/100,
                                                             dataset=="All"~n/3560))

figdata$frame<-factor(figdata$frame,levels=c("Bleached Specific","Non-bleached Specific","Cooperative","Antagonistic","No Response"),labels=c("B-Specific","NB-Specific","Cooperative","Antagonistic","No Response"))
b<-ggplot(figdata)+geom_col(aes(dataset,prop,fill=frame))+
     scale_fill_viridis(discrete=TRUE,name="",guide = guide_legend(nrow=2,reverse = TRUE,byrow=TRUE))+
     theme_classic(base_size=8)+ylab("Proportion")+xlab("")+
     theme(legend.key.size = unit(0.3,"cm"),
           legend.position="bottom")+coord_flip();b

fisher.test(rbind(c(bleached+nonbleached,bleached_a+nonbleached_a),c(100,3560)),alternative="greater")$p.value
figdata%>%filter(dataset!="All")%>%summarise(sum=sum(prop))
figdata<-figdata%>%mutate(prop=case_when(dataset=="Top VI" ~prop*(100/99.0099),
                                dataset=="All"~prop))%>%filter(dataset!="All")

########################################## FAMILY ANALYSIS ################################################################ #####
data<-corr%>%filter(molecular_family!="ontaminant"&molecular_family!="contaminant"&molecular_family!="unknown")
magic_for(print,silent=TRUE)

#test if change in each family is significantly different from 0 for each phenotype
for (i in unique(data$molecular_family)){
     print(t.test((data%>%filter(molecular_family==i))$change_b, mu = 0, alternative = "two.sided")$p.value)
     print(t.test((data%>%filter(molecular_family==i))$change_nb, mu = 0, alternative = "two.sided")$p.value)
}
out<-magic_result_as_dataframe()%>%select(1,2,4)%>%rename(molecular_family=1,b_p_value=2,nb_p_value=3)
sig<-format(bind_cols(out,as.data.frame(p.adjust(out$b_p_value,method="BH")),as.data.frame(p.adjust(out$nb_p_value,method="BH")))%>%
                 rename(b_adjust=4,nb_adjust=5)%>%mutate(b_adjust=round(b_adjust,3))%>%mutate(nb_adjust=round(nb_adjust,3))%>%
                 mutate(sig=case_when(b_adjust<=.05&nb_adjust<=.05~"Both",
                                      b_adjust<=.05&nb_adjust>.05~"Bleached",
                                      b_adjust>.05&nb_adjust<=.05~"Non-bleached",
                                      TRUE ~ "NS")),scientific=FALSE)%>%
     select(-b_p_value,-nb_p_value)
#write.table((sig%>%mutate(molecular_family=str_to_title(molecular_family))%>%rename('Molecular Family'=1,'Adjusted p-value (Bleached)'=2,'Adjusted p-value (Non-bleached)'=3,'Significance'=4)),"../figs/STx_family_change.txt",sep="\t",quote=FALSE,row.names=FALSE)
fam_viz<-left_join(corr,sig,by="molecular_family")%>%filter(sig!="NS")%>%filter(molecular_family!="unknown")%>%ungroup()%>%mutate(molecular_family=str_to_title(molecular_family))
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
     
     gb <- ggplot_build(p)
     lay <- gb$layout$layout
     tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
     p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                   vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}
my_tag <- c("Both p<0.003","Bleached p<0.001","Bleached p=0.010","Bleached p=0.002","Both p<0.002","Bleached p=0.002","Bleached p=0.027")
c<-ggplot(fam_viz)+
     theme_classic(base_size=8)+
     scale_color_manual(values=cols,name="")+
     geom_vline(xintercept=0,linetype="dotted")+
     geom_hline(yintercept=0,linetype="dotted")+
     scale_y_continuous(limits=c(-200000000,200000000))+
     scale_x_continuous(limits=c(-200000000,200000000))+
     ylab("Non-bleached NHR")+
     xlab("Bleached NHR")+
     geom_point(aes(change_b,change_nb,color=molecular_family),size=0.5,alpha=.5)+
     facet_wrap(~molecular_family,ncol=2)+
     theme(legend.position="none",strip.text=element_text(size=5),
           axis.text.x=element_blank(),
           axis.text.y=element_blank());c

c<-tag_facet2(c, 
              x = -2e+08, y = -2e+08, 
              vjust = 0, hjust = 0,
              open = "", close = "",
              fontface = 'italic',
              size = 2,
              tag_pool = my_tag)

quartz(w=5,h=4)
plots <- align_plots(b,c, align = 'h', axis = 'b')
p1<-plot_grid(a,plots[[1]],ncol=1,align="v",axis="l",rel_heights=c(2.5,1),labels=c("A","B"),label_size=8,label_x=c(0.16,0.16),label_y=c(1,1.1))
plot_grid(p1,plots[[2]],ncol=2,rel_widths=c(1.7,1),labels=c("","C"),label_size=8,label_x=c(1,0))

#take each category, reclassify as category + other, run MWU test on ranks
distvals<-corr%>%mutate(dist=sqrt((change_b^2)+(change_nb^2)))%>%filter(molecular_family!="ontaminant"&molecular_family!="contaminant"&molecular_family!="unknown")
magic_for(print,silent=TRUE)
for(i in unique(distvals$molecular_family)){     
     data<-distvals%>%mutate(split=case_when(molecular_family==i ~"a",molecular_family!=i ~"z"));
     print(wilcox.test(dist~split,data=data,alternative="greater")$p.value)
}
mwu<-magic_result_as_dataframe()
mwu_sig<-bind_cols(mwu,as.data.frame(p.adjust(mwu[,2],method="BH")))%>%rename(molecular_family=1,p_value=2,p_adjust=3)%>%
     mutate(sig=case_when(p_adjust<=.05~"Significant", TRUE ~ "NS"))%>%mutate(molecular_family=str_to_title(molecular_family))%>%select(-p_value)%>%mutate(p_adjust=round(p_adjust,3))%>%
     rename('Molecular Family'=1,'Adjusted p-value'=2,'Significance'=3)

#write.table(mwu_sig,"../figs/STx_family_magnitude.txt",sep="\t",quote=FALSE,row.names=FALSE)
########################################## OUTLIER ANALYIS ################################################################ #####
nrow(corr%>%filter(change_b>rb|change_b<lb|change_nb>tb|change_nb<bb)) #1014 outliers
nrow(corr) #3560 molecules, expected proportion 0.284
exp<-nrow(corr%>%filter(change_b>rb|change_b<lb|change_nb>tb|change_nb<bb))/nrow(corr)
1014/3560

#by top 5% distance from origin
nr1<-nrow(corr%>%ungroup()%>%mutate(dist=sqrt((change_b^2)+(change_nb^2)))%>%top_frac(0.05));nr2<-nrow(corr%>%ungroup()%>%mutate(dist=sqrt((change_b^2)+(change_nb^2)))%>%top_frac(-0.95))
outdata<-left_join(corr%>%ungroup()%>%mutate(dist=sqrt((change_b^2)+(change_nb^2)))%>%top_frac(0.05)%>%group_by(molecular_family)%>%tally()%>%rename(actual=n),
                   corr%>%ungroup()%>%mutate(dist=sqrt((change_b^2)+(change_nb^2)))%>%top_frac(-0.95)%>%group_by(molecular_family)%>%tally()%>%mutate(expected=n)%>%select(-n),by="molecular_family")%>%filter(molecular_family!="ontaminant"&molecular_family!="contaminant"&molecular_family!="unknown")%>% replace(is.na(.), 0)

magic_for(print,silent=TRUE)
for (i in 1:nrow(outdata)){
     print(fisher.test(rbind(c(outdata[[i,2]],outdata[[i,3]]),c(nr1,nr2)), alternative="greater")$p.value)
}

out<-magic_result_as_dataframe()%>%rename(p_value=3)%>%select(p_value)
sig<-bind_cols(outdata,out,as.data.frame(p.adjust(out$p_value,method="BH")))%>%rename(p_adjust=5)%>%
     mutate(sig=case_when(p_adjust<=.05~"Significant",TRUE ~ "NS"))%>%mutate(molecular_family=str_to_title(molecular_family))%>%select(-p_value)%>%mutate(p_adjust=round(p_adjust,3))%>%
rename('Molecular Family'=1,'Actual Count'=2,'Expected Count'=3,'Adjusted p-value'=4,'Significance'=5)

#write.table(sig,"../figs/STx_outlier_enrichment.txt",sep="\t",quote=FALSE,row.names=FALSE)
############################## TEMP STRESS PCOA SETUP (F4) ################################################################ #####
metadata<-read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%select(colony,bleaching_history_phenotype)%>%distinct()
#averaged initial for all fragments by genotype plus constant fragments at the end
raw<-left_join(bind_rows(read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%
                              select(-filename,-feature,-mac,-well_id,-sample_id,-experiment,-sum,-sample_number,-plate_number,-nubbin_number,-tank)%>%
                              filter(time_point=="T1")%>%group_by(colony)%>%summarise_if(is.numeric, mean, na.rm = TRUE)%>%mutate(time_point="T1"),
                         read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%
                              filter(temperature_treatment=="Constant")%>%filter(time_point=="T2")%>%
                              select(-filename,-feature,-mac,-well_id,-sample_id,-experiment,-sum,-sample_number,-plate_number,-nubbin_number,-tank)%>%
                              mutate(time_point="T2"))%>%
                    select(colony,time_point,everything(),-bleaching_history_phenotype),metadata,by="colony")%>%
     select(colony,time_point,bleaching_history_phenotype,temperature_treatment,everything())
list<-read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%group_by(molecular_family)%>%tally()

########################################## OVERALL (F4) ################################################################### #####
data<-left_join(raw%>%gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
                     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))%>%
                     separate(molecule,into=c("cluster_id","other"),sep="_",remove=FALSE)%>%
                     separate(cluster_id,into=c("trash","cluster_id"),sep=1)%>%select(-other,-trash)%>%
                     mutate(cluster_id=as.numeric(cluster_id)),read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(cluster_id,molecular_family)%>%mutate(cluster_id=as.numeric(cluster_id)),by="cluster_id")%>%
     select(-cluster_id,-molecular_family)%>%
     spread(molecule,aoc)%>%filter(colony!=203&colony!=214)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$time_point <- factor(data$time_point,labels=c("Initial","Heat"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"));data$colony <- as.factor(data$colony)

pca<-prcomp(data[,5:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:4],axes$data)
data.env<-plotdata%>%select(time_point,bleaching_history_phenotype,colony,x,y)
#quartz()
legend<-get_legend(ggplot(plotdata)+geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+geom_path(aes(x,y,group=colony),arrow=arrow(length = unit(0.2, "cm")),color="gray")+theme_classic(base_size=8)+
                        scale_fill_manual(values=c("white","darkgray"),name="")+scale_color_manual(values=c("black","red"),name="")+theme(legend.key.size = unit(0.25,"cm"),legend.position="bottom"))
set.seed(3839)
dist<-as.matrix(dist(plotdata%>%select(x,y)));adonis(dist ~ bleaching_history_phenotype*time_point, data = data.env)
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Non-bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Non-bleached"))) 
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Bleached"))) 
a<-ggplot(plotdata)+
     geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+
     geom_path(aes(x,y,group=colony),color="gray")+
     stat_ellipse(aes(x,y,group=bleaching_history_phenotype),color="gray",linetype="dotted")+
     theme_classic(base_size=8)+
     xlab(axes$labels$x)+
     ylab(axes$labels$y)+
     scale_fill_manual(values=c("white","darkgray"),name="")+
     scale_color_manual(values=c("black","red"),name="")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           plot.title = element_text(size = 8))+
     ggtitle("All Compounds")+
     annotate("text",x=-90,y=-145,hjust=0,vjust=0,label="~Phenotype*Time p=0.018",size=2,fontface = 'italic');a

#t_data<-left_join(bind_cols(plotdata%>%select(colony,time_point,x)%>%spread(time_point,x),
#                            plotdata%>%select(colony,time_point,y)%>%spread(time_point,y))%>%
#                       select(-colony1)%>%mutate(dist=sqrt((((Initial-Heat)^2)+((Initial1-Heat1)^2)))),metadata%>%mutate(colony=as.factor(colony)),by="colony")
#leveneTest(dist~bleaching_history_phenotype,data=t_data)     
#t.test(dist~bleaching_history_phenotype,data=t_data)     

########################################## SYMBIONT (F4) ################################################################## #####
metadata<-read_xlsx("bucket_table.xlsx",sheet="T1_T2_symbiont")%>%clean_names()%>%select(colony,bleach_status)%>%distinct()
#averaged initial for all fragments by genotype plus constant fragments at the end
raw<-bind_rows(read_xlsx("bucket_table.xlsx",sheet="T1_T2_symbiont")%>%clean_names()%>%
                    select(colony,bleach_status,everything(),-filename,-coral,-replicate,-experiment,-bnb)%>%
                    filter(sample_type=="T1Coral")%>%group_by(colony)%>%summarise_if(is.numeric, mean, na.rm = TRUE)%>%mutate(sample_type="T1")%>%ungroup(),
               read_xlsx("bucket_table.xlsx",sheet="T1_T2_symbiont")%>%clean_names()%>%
                    select(-filename,-coral,-replicate,-experiment,-bnb)%>%
                    filter(temp=="Constant")%>%filter(sample_type=="T2Coral")%>%
                    mutate(sample_type="T2"))%>%
     select(colony,sample_type,bleach_status,everything())%>%
     select(colony,sample_type,everything())

data<-raw%>%select(-bleach_status)%>%inner_join(.,metadata,by="colony")%>%select(bleach_status,colony,everything())
data$bleaching_history_phenotype <- factor(data$bleach_status,labels=c("Bleached","Non-bleached"))
data$time_point <- factor(data$sample_type,labels=c("Initial","Heat"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"));data$colony <- as.factor(data$colony)
data<-data%>%select(time_point,temp,colony,bleaching_history_phenotype,bleach_status,sample_type,everything(),-temp,-bleach_status,-sample_type)%>%
     gather(molecule,aoc,-colony,-bleaching_history_phenotype,-time_point)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))%>%
     spread(molecule,aoc)

pca<-prcomp(data[,4:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:4],axes$data)
set.seed(3839)
dist<-as.matrix(dist(plotdata%>%select(x,y)));adonis(dist ~ bleaching_history_phenotype*time_point, data = plotdata)
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Non-bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (plotdata%>%filter(bleaching_history_phenotype=="Non-bleached"))) 
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (plotdata%>%filter(bleaching_history_phenotype=="Bleached"))) 

b<-ggplot(plotdata)+
     geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+
     geom_path(aes(x,y,group=colony),color="gray")+
     stat_ellipse(aes(x,y,group=bleaching_history_phenotype),color="gray",linetype="dotted")+
     theme_classic(base_size=8)+
     xlab(axes$labels$x)+
     ylab(axes$labels$y)+
     scale_fill_manual(values=c("white","darkgray"),name="")+
     scale_color_manual(values=c("black","red"),name="")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           plot.title = element_text(size = 8))+
     ggtitle("Symbiont Compounds")+
     annotate("text",x=-270,y=-250,hjust=0,vjust=0,label="~Phenotype p<0.001\n~Phenotype|Time p<0.012",size=2,fontface = 'italic')

########################################## HOST (F4) ###################################################################### #####
metadata<-read_xlsx("bucket_table.xlsx",sheet="T1_T2_coral")%>%clean_names()%>%select(colony,bleach_status)%>%distinct()
#averaged initial for all fragments by genotype plus constant fragments at the end
raw<-bind_rows(read_xlsx("bucket_table.xlsx",sheet="T1_T2_coral")%>%clean_names()%>%
                    select(colony,bleach_status,everything(),-filename,-coral,-replicate,-experiment,-bnb)%>%
                    filter(sample_type=="T1Coral")%>%group_by(colony)%>%summarise_if(is.numeric, mean, na.rm = TRUE)%>%mutate(sample_type="T1")%>%ungroup(),
               read_xlsx("bucket_table.xlsx",sheet="T1_T2_coral")%>%clean_names()%>%
                    select(-filename,-coral,-replicate,-experiment,-bnb)%>%
                    filter(temp=="Constant")%>%filter(sample_type=="T2Coral")%>%
                    mutate(sample_type="T2"))%>%
     select(colony,sample_type,bleach_status,everything())%>%
     select(colony,sample_type,everything())

data<-raw%>%select(-bleach_status)%>%inner_join(.,metadata,by="colony")%>%select(bleach_status,colony,everything())
data$bleaching_history_phenotype <- factor(data$bleach_status,labels=c("Bleached","Non-bleached"))
data$time_point <- factor(data$sample_type,labels=c("Initial","Heat"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"));data$colony <- as.factor(data$colony)
data<-data%>%select(time_point,temp,colony,bleaching_history_phenotype,bleach_status,sample_type,everything(),-temp,-bleach_status,-sample_type)%>%
     gather(molecule,aoc,-colony,-bleaching_history_phenotype,-time_point)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))%>%
     spread(molecule,aoc)

pca<-prcomp(data[,4:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:4],axes$data)
set.seed(3839)
dist<-as.matrix(dist(plotdata%>%select(x,y)));adonis(dist ~ bleaching_history_phenotype*time_point, data = plotdata)
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Non-bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (plotdata%>%filter(bleaching_history_phenotype=="Non-bleached"))) 
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (plotdata%>%filter(bleaching_history_phenotype=="Bleached"))) 

c<-ggplot(plotdata)+
     geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+
     geom_path(aes(x,y,group=colony),color="gray")+
     stat_ellipse(aes(x,y,group=bleaching_history_phenotype),color="gray",linetype="dotted")+
     theme_classic(base_size=8)+
     xlab(axes$labels$x)+
     ylab(axes$labels$y)+
     scale_fill_manual(values=c("white","darkgray"),name="")+
     scale_color_manual(values=c("black","red"),name="")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           plot.title = element_text(size = 8))+
     ggtitle("Host Compounds")+
     annotate("text",x=-135,y=-145,hjust=0,vjust=0,label="~Phenotype p<0.001\n~Phenotype|Time p<0.027",size=2,fontface = 'italic')

########################################## PLOTS ########################################################################## #####
#plot for overall, host, symb
quartz(w=5.5,h=2) 
plot_grid(plot_grid(a,b,c,align="h",axis="l",rel_widths=c(1,1,1),nrow=1,labels=c("A","B","C"),label_size=8,label_x=c(0.23,0.23,0.23),label_y=c(0.9,0.9,0.9)),NULL,legend,nrow=3,rel_heights=c(5,-0.2,1))


########################################## MEMBRANE LIPIDS (SF) ########################################################### #####
metadata<-read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%select(colony,bleaching_history_phenotype)%>%distinct()
#averaged initial for all fragments by genotype plus constant fragments at the end
raw<-left_join(bind_rows(read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%
                              select(-filename,-feature,-mac,-well_id,-sample_id,-experiment,-sum,-sample_number,-plate_number,-nubbin_number,-tank)%>%
                              filter(time_point=="T1")%>%group_by(colony)%>%summarise_if(is.numeric, mean, na.rm = TRUE)%>%mutate(time_point="T1"),
                         read_xlsx("bucket_table.xlsx",sheet="T1_T2")%>%clean_names%>%
                              filter(temperature_treatment=="Constant")%>%filter(time_point=="T2")%>%
                              select(-filename,-feature,-mac,-well_id,-sample_id,-experiment,-sum,-sample_number,-plate_number,-nubbin_number,-tank)%>%
                              mutate(time_point="T2"))%>%
                    select(colony,time_point,everything(),-bleaching_history_phenotype),metadata,by="colony")%>%
     select(colony,time_point,bleaching_history_phenotype,temperature_treatment,everything())

data<-left_join(raw%>%gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
                     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))%>%
                     separate(molecule,into=c("cluster_id","other"),sep="_",remove=FALSE)%>%
                     separate(cluster_id,into=c("trash","cluster_id"),sep=1)%>%select(-other,-trash)%>%
                     mutate(cluster_id=as.numeric(cluster_id)),read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(cluster_id,molecular_family)%>%mutate(cluster_id=as.numeric(cluster_id)),by="cluster_id")%>%
     filter(molecular_family=="betaine"|molecular_family=="phosphocholine")%>%
     select(-cluster_id,-molecular_family)%>%
     spread(molecule,aoc)%>%filter(colony!=203&colony!=214)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$time_point <- factor(data$time_point,labels=c("Initial","Heat"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"));data$colony <- as.factor(data$colony)

pca<-prcomp(data[,5:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:4],axes$data)
data.env<-plotdata%>%select(time_point,bleaching_history_phenotype,colony,x,y)

legend<-get_legend(ggplot(plotdata)+geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+geom_path(aes(x,y,group=colony),arrow=arrow(length = unit(0.2, "cm")),color="gray")+theme_classic(base_size=8)+
                        scale_fill_manual(values=c("white","darkgray"),name="")+scale_color_manual(values=c("black","red"),name="")+theme(legend.key.size = unit(0.25,"cm"),legend.position="bottom"))
set.seed(3839)
dist<-as.matrix(dist(plotdata%>%select(x,y)));adonis(dist ~ bleaching_history_phenotype*time_point, data = data.env)
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Non-bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Non-bleached"))) 
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Bleached"))) 
memb<-ggplot(plotdata)+
     geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+
     geom_path(aes(x,y,group=colony),color="gray")+
     stat_ellipse(aes(x,y,group=bleaching_history_phenotype),color="gray",linetype="dotted")+
     theme_classic(base_size=8)+
     xlab(axes$labels$x)+
     ylab(axes$labels$y)+
     scale_fill_manual(values=c("white","darkgray"),name="")+
     scale_color_manual(values=c("black","red"),name="")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           plot.title = element_text(size = 8))+
     ggtitle("Membrane Lipids")+
     annotate("text",x=-33,y=-39,hjust=0,vjust=0,label="Phenotype p=0.002",size=2,fontface = 'italic');memb

########################################## MICROBIAL NATURAL PRODUCTS (SF) ################################################ #####
data<-left_join(raw%>%gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
                     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))%>%
                     separate(molecule,into=c("cluster_id","other"),sep="_",remove=FALSE)%>%
                     separate(cluster_id,into=c("trash","cluster_id"),sep=1)%>%select(-other,-trash)%>%
                     mutate(cluster_id=as.numeric(cluster_id)),read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(cluster_id,molecular_family)%>%mutate(cluster_id=as.numeric(cluster_id)),by="cluster_id")%>%
     filter(molecular_family=="indole"|molecular_family=="microbial natural product")%>%
     select(-cluster_id,-molecular_family)%>%
     spread(molecule,aoc)%>%filter(colony!=203&colony!=214)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$time_point <- factor(data$time_point,labels=c("Initial","Heat"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"));data$colony <- as.factor(data$colony)

pca<-prcomp(data[,5:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:4],axes$data)
data.env<-plotdata%>%select(time_point,bleaching_history_phenotype,colony,x,y)

set.seed(3839)
dist<-as.matrix(dist(plotdata%>%select(x,y)));adonis(dist ~ bleaching_history_phenotype*time_point, data = data.env)
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Non-bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Non-bleached"))) 
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Bleached"))) 
micro<-ggplot(plotdata)+
     geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+
     geom_path(aes(x,y,group=colony),color="gray")+
     stat_ellipse(aes(x,y,group=bleaching_history_phenotype),color="gray",linetype="dotted")+
     theme_classic(base_size=8)+
     xlab(axes$labels$x)+
     ylab(axes$labels$y)+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_manual(values=c("black","red"),name="Time Point")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           plot.title = element_text(size = 8))+
     ggtitle("Putatitve Microbial Products")+
     annotate("text",x=-7,y=-5,hjust=0,vjust=0,label="Phenotype p=0.005",size=2,fontface = 'italic')

########################################## BIOACTIVE MOLECULES (SF) ####################################################### #####
data<-left_join(raw%>%gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
                     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))%>%
                     separate(molecule,into=c("cluster_id","other"),sep="_",remove=FALSE)%>%
                     separate(cluster_id,into=c("trash","cluster_id"),sep=1)%>%select(-other,-trash)%>%
                     mutate(cluster_id=as.numeric(cluster_id)),read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(cluster_id,molecular_family)%>%mutate(cluster_id=as.numeric(cluster_id)),by="cluster_id")%>%
     filter(molecular_family=="prostaglandin"|molecular_family=="steroid"|molecular_family=="eicosanoid"|molecular_family=="phosphatidic acids"|grepl('paf_c|x7771_', molecule))%>%
     select(-cluster_id,-molecular_family)%>%
     spread(molecule,aoc)%>%filter(colony!=203&colony!=214)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$time_point <- factor(data$time_point,labels=c("Initial","Heat"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"));data$colony <- as.factor(data$colony)

pca<-prcomp(data[,5:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:4],axes$data)
data.env<-plotdata%>%select(time_point,bleaching_history_phenotype,colony,x,y)

set.seed(3839)
dist<-as.matrix(dist(plotdata%>%select(x,y)));adonis(dist ~ bleaching_history_phenotype*time_point, data = data.env)
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Non-bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Non-bleached"))) 
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Bleached"))) 
bio<-ggplot(plotdata)+
     geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+
     geom_path(aes(x,y,group=colony),color="gray")+
     stat_ellipse(aes(x,y,group=bleaching_history_phenotype),color="gray",linetype="dotted")+
     theme_classic(base_size=8)+
     xlab(axes$labels$x)+
     ylab(axes$labels$y)+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_manual(values=c("black","red"),name="Time Point")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           plot.title = element_text(size = 8))+
     ggtitle("Bioactive Molecules")+
     annotate("text",x=-33,y=-33,hjust=0,vjust=0,label="Phenotype p<0.001",size=2,fontface = 'italic');bio

########################################## PHOTOSYNTHETIC PRODUCTS (SF) ################################################### #####
data<-left_join(raw%>%gather(molecule,aoc,-colony,-bleaching_history_phenotype,-temperature_treatment,-time_point)%>%
                     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0,TRUE ~ as.numeric(aoc)))%>%
                     separate(molecule,into=c("cluster_id","other"),sep="_",remove=FALSE)%>%
                     separate(cluster_id,into=c("trash","cluster_id"),sep=1)%>%select(-other,-trash)%>%
                     mutate(cluster_id=as.numeric(cluster_id)),read_xlsx("component_ID.xlsx",sheet="Molecular Families")%>%clean_names()%>%select(cluster_id,molecular_family)%>%mutate(cluster_id=as.numeric(cluster_id)),by="cluster_id")%>%
     filter(molecular_family=="chlorophyll"|molecular_family=="xanthin")%>%
     select(-cluster_id,-molecular_family)%>%
     spread(molecule,aoc)%>%filter(colony!=203&colony!=214)

data$bleaching_history_phenotype <- factor(data$bleaching_history_phenotype,labels=c("Bleached","Non-bleached"))
data$time_point <- factor(data$time_point,labels=c("Initial","Heat"))
data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"));data$colony <- as.factor(data$colony)

pca<-prcomp(data[,5:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:4],axes$data)
data.env<-plotdata%>%select(time_point,bleaching_history_phenotype,colony,x,y)

set.seed(3839)
dist<-as.matrix(dist(plotdata%>%select(x,y)));adonis(dist ~ bleaching_history_phenotype*time_point, data = data.env)
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Non-bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Non-bleached"))) 
dist<-as.matrix(dist(plotdata%>%filter(bleaching_history_phenotype=="Bleached")%>%select(x,y)));adonis(dist ~ time_point, data = (data.env%>%filter(bleaching_history_phenotype=="Bleached"))) 
quartz()
d<-ggplot(plotdata)+
     geom_point(aes(x,y,fill=bleaching_history_phenotype,color=time_point),pch=21,size=2,stroke = 1)+
     geom_path(aes(x,y,group=colony),color="gray")+
     stat_ellipse(aes(x,y,group=bleaching_history_phenotype),color="gray",linetype="dotted")+
     theme_classic(base_size=8)+
     xlab(axes$labels$x)+
     ylab(axes$labels$y)+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_manual(values=c("black","red"),name="Time Point")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           plot.title = element_text(size = 8))+
     ggtitle("Photosynthetic Products")+
     annotate("text",x=-5,y=-5,hjust=0,vjust=0,label="Phenotype p=0.078",size=2,fontface = 'italic');d


########################################## PLOTS ########################################################################## #####
#plot for bioactive compounds, membrane lipids, microbial products
quartz(h=2,w=5.5)
plot_grid(plot_grid(memb,bio,micro,align="h",axis="l",labels=c("A","B","C"),nrow=1,label_size=8,label_x=c(0.19,0.19,0.21),label_y=c(0.9,0.9,0.9)),NULL,legend,nrow=3,rel_heights=c(5,-0.2,1))



############################## SYMBIONT TEMP RESPONSE (SF) ################################################################ #####
library(tidyverse);library(readxl);library(cowplot)
data<-read_excel("symbionts.xlsx",sheet="fvfm")%>%select(Treatment,Fragment,Colony,Phenotype,Tank,T0_fvfm,T5_fvfm)%>%filter(T5_fvfm>0.425)%>%filter(T0_fvfm>0.5)%>%filter(T0_fvfm<0.85)%>%
     #mutate(T0=1,T5=T5_fvfm/T0_fvfm)%>%select(-T5_fvfm,-T0_fvfm)%>%
     gather(time,fvfm,-Treatment,-Fragment,-Phenotype,-Tank,-Colony)%>%
     filter(Treatment=="Constant High"|Treatment=="Control")

data$Phenotype <- factor(data$Phenotype,labels=c("Bleached","Non-bleached"))
data$time <- factor(data$time,labels=c("Initial","Final"))

a<-ggplot(data)+
     #geom_point(aes(time,fvfm))+
     facet_wrap(~Treatment)+
     geom_line(aes(time,fvfm,group=Fragment,color=Phenotype),alpha=0.5)+
     theme_classic(base_size=8)+
     ylab("Fv/fm")+
     theme(legend.key.size=unit(0.2,"cm"),axis.title.x=element_blank(),legend.position=c(0.85,0.15))+
     scale_color_manual(values=c("orange","gray"),name="")+guides(linetype=FALSE);a

tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
     
     gb <- ggplot_build(p)
     lay <- gb$layout$layout
     tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
     p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                   vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

my_tag <- c("~Treament p<0.001", "")
b<-tag_facet2(a, 
             x = -Inf, y = -Inf, 
             vjust = -1, hjust = -0.1,
             open = "", close = "",
             fontface = 'italic',
             size = 2,
             tag_pool = my_tag);b

x<-data%>%spread(time,fvfm)%>%mutate(diff=Final-Initial)
summary(aov(diff~Treatment*Phenotype,data=x))
x%>%group_by(Treatment)%>%summarise(mean=mean(diff))
0.0447/0.00738

quartz(w=3,h=2);a

data<-read_excel("symbionts.xlsx",sheet="summary")%>%mutate(year=as.factor(Year))%>%filter(year==2018)%>%group_by(colony,phenotype)%>%
     summarise(prop_d=mean(prop_d),prop_c=mean(prop_c))%>%
     gather(symb,prop,-phenotype,-colony)

data$colony <- factor(data$colony, levels=c("222","214","202","20","12","2","210","212","220","204","4","203","201","19","11","211","1","221","209","213","219","3"))
data$symb<- factor(data$symb, labels=c("Cladocopium","Durusdinium"))
b<-ggplot(data)+geom_bar(aes(colony,prop,fill=symb),stat="identity",position="stack")+
     theme_classic(base_size=8)+scale_y_continuous(breaks=seq(0,1.1,0.25))+ylab("2018\nProportion")+
     scale_fill_manual(values=c("goldenrod1","navy"),labels=c(expression(italic("Cladocopium")),expression(italic("Durusdinium"))))+
     theme(legend.key.size=unit(0.2,"cm"),axis.title.x=element_blank(),legend.title=element_blank(),
     legend.position="bottom",axis.text.x=element_text(angle=90));b

data<-read_excel("symbionts.xlsx",sheet="summary")%>%mutate(year=as.factor(Year))%>%group_by(colony,year,phenotype)%>%
     summarise(prop_d=mean(prop_d),prop_c=mean(prop_c))%>%
     gather(symb,prop,-phenotype,-colony,-year)%>%filter(year==2019)

data$colony <- factor(data$colony, levels=c("222","214","202","20","12","2","210","212","220","204","4","203","201","19","11","211","1","221","209","213","219","3"))
data$symb<- factor(data$symb, labels=c("Cladocopium","Durusdinium"))
c<-ggplot(data)+geom_bar(aes(colony,prop,fill=symb),stat="identity",position="stack")+
     theme_classic(base_size=8)+scale_y_continuous(breaks=seq(0,1,0.25))+ylab("2019\nProportion")+
     scale_fill_manual(values=c("goldenrod1","navy"),labels=c(expression(italic("Cladocopium")),expression(italic("Durusdinium"))))+
     theme(legend.key.size=unit(0.2,"cm"),axis.title.x=element_blank(),legend.title=element_blank(),
           legend.position="none",axis.text.x=element_text(angle=90))+
     annotate("segment", x = 1, xend = 11, y = 1.03, yend = 1.03,colour = "black",size=0.5)+
     annotate("segment", x = 12, xend = 22, y = 1.03, yend = 1.03,colour = "black",size=0.5)+
     annotate("text", x = 6, y = 1.11,colour = "black",label="Non-bleached",size=2)+
     annotate("text", x = 17.5, y = 1.11,colour = "black",label="Bleached",size=2)


quartz(w=3,h=4)
plot_grid(c,b,ncol=1,rel_heights=c(1,1.1),labels=c("A","B"),label_size=8)




############################## PCOAs SYMBIONT PELLET ###################################################################### #####
data<-read_xlsx("bucket_table.xlsx",sheet="Symbiont_pellet")%>%clean_names()%>%select(-attribute_coral_species,-attribute_replicate_id,-experiment,-attribute_sample_type,-attribute_time_point,-attribute_colony,-temperature_treatment)%>%
     gather(molecule,aoc,-attribute_phenotype,-filename,-wash)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0, TRUE ~ as.numeric(aoc)))%>%distinct()%>%
     spread(molecule,aoc)%>%select(-filename)%>%filter(wash==3)%>%select(-wash)

data$attribute_phenotype <- factor(data$attribute_phenotype,levels=c("B","NB"),labels=c("Bleached","Non-bleached"))
#data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
#data$colony <- as.factor(data$colony)

pca<-prcomp(data[,2:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1],axes$data)
coord<-plotdata%>%select(x,y)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ attribute_phenotype, data = plotdata) #PERMANOVA 

axes$labels$x
axes$labels$y
pellet<-ggplot(plotdata)+geom_point(aes(x,y,fill=attribute_phenotype),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     #scale_y_continuous(limits=c(-300,300),breaks=seq(-300,300,100))+
     #scale_x_continuous(limits=c(-300,300),breaks=seq(-300,300,100))+     
     xlab("Component 1 (55.7%)")+
     ylab("Component 2 (26.6%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="left",
           legend.spacing.y = unit(0.0, "cm"),
           legend.margin= margin(0.2,0,0,0, unit="cm"),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.background = element_rect(fill = "transparent",colour = NA))+
     annotate("text",x=-200,y=-130,hjust=0,vjust=1,label="PERMANOVA ~Phenotype p=0.1",size=1.5,fontface = 'italic');pellet

############################## PCOAs BLEACHED FRAG ######################################################################## #####
data<-read_xlsx("bucket_table.xlsx",sheet="Bleached_Frag")%>%clean_names()%>%filter(attribute_time_point=="Bleached")%>%
     select(-attribute_coral_species,-attribute_replicate_id,-experiment,-attribute_sample_type,-attribute_time_point,-attribute_colony,-temperature_treatment)%>%
     gather(molecule,aoc,-attribute_phenotype,-filename,-location)%>%
     mutate(aoc=log(aoc))%>%mutate(aoc=case_when(aoc=="-Inf"~0, TRUE ~ as.numeric(aoc)))%>%distinct()%>%
     spread(molecule,aoc)%>%select(-filename)%>%filter(location=="field")%>%select(-location)

data$attribute_phenotype <- factor(data$attribute_phenotype,levels=c("B","NB"),labels=c("Bleached","Non-bleached"))
#data$colony <- factor(data$colony, levels=c("222","214","202","20","12","211","203","201","19","11"))
#data$colony <- as.factor(data$colony)

pca<-prcomp(data[,2:ncol(data)])
axes<-fviz_pca_ind(pca,axes = c(1,2))
plotdata<-bind_cols(data[,1:3],axes$data)
coord<-plotdata%>%select(x,y)
data<-as.matrix(dist(coord))
set.seed(3839);adonis2(data ~ attribute_phenotype, data = plotdata,nperm=100) #PERMANOVA 

axes$labels$x
axes$labels$y
frag<-ggplot(plotdata)+geom_point(aes(x,y,fill=attribute_phenotype),pch=21,size=2,stroke = 1)+
     theme_classic(base_size=8)+
     #scale_y_continuous(limits=c(-300,300),breaks=seq(-300,300,50))+
     #scale_x_continuous(limits=c(-300,300),breaks=seq(-300,300,50))+     
     xlab("Component 1 (40.9%)")+
     ylab("Component 2 (19.4%)")+
     scale_fill_manual(values=c("white","darkgray"),name="Historical\nBleaching\nPhenotype")+
     scale_color_viridis(discrete=TRUE,name="Colony")+
     theme(legend.key.size = unit(0.25,"cm"),
           legend.position="none",
           legend.spacing.y = unit(0.0, "cm"),
           legend.margin= margin(0.2,0,0,0, unit="cm"),
           plot.background = element_rect(fill = "transparent",colour = NA),
           panel.background = element_rect(fill = "transparent",colour = NA))+
     annotate("text",x=-110,y=-125,hjust=0,vjust=1,label="PERMANOVA ~Phenotype p=0.1",size=1.5,fontface = 'italic');frag
quartz(w=5.5,h=2)
plot_grid(pellet,frag,nrow=1,rel_widths = c(1.3,1),labels=c("A","B"),label_size=8,label_x=c(0.39,0.17))



