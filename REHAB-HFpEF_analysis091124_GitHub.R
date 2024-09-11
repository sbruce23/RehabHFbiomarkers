rm(list=ls())

library(tidyverse)
library(XLConnect)
library(janitor)
library(gt)
library(gtsummary)
library(reshape2)

############################################
## 0.1 Data dictionary
#############################################
#age (years)
#sex (1 for female, 0 for male)
#race___4 (1 for white, 0 for non-white)
#hf_cat (1 for hfpef (EF >=45%), 0 for hfref (EF < 45%))
#creatinine_value (mg/dL)
#bl_sppb, fu_sppb (baseline (bl) and follow up (fu) SPPB score)
#bl_smw, fu_smw (baseline (bl) and follow up (fu) six minute walk test distance (m))
#rehosp (rehospitalization)
#death (death)
#event (rehospitalization or death)
#lfu (lost to follow up)
#eGFR (derived from age, sex, and creatinine level 
#(https://www.kidney.org/content/ckd-epi-creatinine-equation-2021))

############################################
## 0.2 Set theme for ggplot2 plots ##
############################################
hw <- theme_gray()+ theme(
  plot.title=element_text(hjust=0.5,size=18,face="bold"),
  plot.subtitle=element_text(hjust=0.5,size=12),
  plot.caption=element_text(hjust=-.5,size=10),
  strip.background=element_rect(fill=rgb(.9,.95,1),
                                colour=gray(.5), linewidth=.2),
  panel.border=element_rect(fill=FALSE,colour=gray(.70)),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.spacing.x = unit(0.2,"cm"),
  panel.spacing.y = unit(0.2,"cm"),
  axis.text=element_text(colour="black",size=10),
  axis.text.y=element_text(margin=ggplot2::margin(0,3,0,3)),
  axis.text.x=element_text(margin=ggplot2::margin(-1,0,3,0)),
  axis.title=element_text(size=16,face="bold"),
  legend.text=element_text(size=14),
  legend.title = element_blank(),
);

############################################
## 1. Read in and initial data pre-processing ##
############################################

## Read in data
## OMITTED
## two data frames created: 
#1) rhfds contains participant clinical and demographic data
#2) rbiods contains participant biomarker measurements

## Clean up column names and data types
rbiods <- clean_names(rbiods);
colnames(rbiods)[3] <- "sample_id";
colnames(rbiods)[4] <- "sample_type";
colnames(rbiods)[5] <- "sample_volume_ml";
colnames(rbiods)[9] <- "hf_category_1_hfpef_0_hfref";
colnames(rbiods)[12] <- "creatinine_mg_dl";
colnames(rbiods)[13] <- "troponin_i_pg_ml";

rhfds$sex <- factor(rhfds$sex,levels=c(0,1));
rhfds$race___4 <- factor(rhfds$race___4,levels=c(0,1));
rhfds$hf_cat <- factor(rhfds$hf_cat,levels=c(0,1));
rhfds$creatinine_value <- as.numeric(rhfds$creatinine_value)
rhfds$bl_sppb <- as.numeric(rhfds$bl_sppb)
rhfds$fu_sppb <- as.numeric(rhfds$fu_sppb)
rhfds$bl_smw <- as.numeric(rhfds$bl_smw)
rhfds$fu_smw <- as.numeric(rhfds$fu_smw)


## Create eGFR
kpa = ifelse(rhfds$sex==1,0.7,0.9);
alpha = ifelse(rhfds$sex==1,-0.241,-0.302)
const = ifelse(rhfds$sex==1,1.012,1)
rhfds$egfr = 142*
  (pmin(as.numeric(rhfds$creatinine_value)/kpa,1)^alpha)*
  (pmax(as.numeric(rhfds$creatinine_value)/kpa,1)^-1.2)*
  (0.9938^rhfds$age)*const

## Clean up entries in rbiods

## Find and remove all (D) indicators
rbiods <- apply(rbiods,2,function(x) gsub("(D)","",x,fixed=TRUE));

## Find and remove all empty spaces in biomarker columns 12-16
rbiods[,12:16] <- apply(rbiods[,12:16],2,function(x) gsub(" ","",x,fixed=TRUE));

## Change QNS to NA and add QNS indicator for all biomarkers (columns 12-16)
qns_ind <- apply(rbiods[,12:16],2,function(x) as.numeric(grepl("QNS",x,fixed=TRUE)));
colnames(qns_ind) <- paste(colnames(qns_ind),"_qns_ind",sep="");
rbiods <- as.data.frame(cbind(rbiods,qns_ind));
rbiods[,12:16] <- apply(rbiods[,12:16],2,function(x) gsub("QNS",NA,x,fixed=TRUE));

## Add censoring indicator columns and change values to NA
cens_ind <- matrix(0,nrow=nrow(rbiods),ncol=5)
cens_ind[apply(rbiods[,12:16],2,function(x) grepl("[><]",x))] <- 
  unlist(apply(rbiods[,12:16],2,function(x) x[grep("[><]",x)]));
colnames(cens_ind) <- paste(colnames(rbiods)[12:16],"_cens_ind",sep="");
rbiods <- as.data.frame(cbind(rbiods,cens_ind));
rbiods[,12:16] <- apply(rbiods[,12:16],2,function(x) gsub("[><]",NA,x));

## Impute censored observations: 
## Halfway between 0 and lower limit of detection for left censored observations
## Upper limit of detection for right censored observations
rbiods$creatinine_mg_dl[rbiods$creatinine_mg_dl_cens_ind != 0] <- 0.15/2;
rbiods$hs_crp_mg_l[rbiods$hs_crp_mg_l_cens_ind != 0] <- 90;
rbiods$nt_pro_bnp[rbiods$nt_pro_bnp_cens_ind != 0] <- 35000;
rbiods$troponin_t_ng_l[rbiods$troponin_t_ng_l_cens_ind != 0] <- 6/2;

## Convert biomarker fields and censoring and QNS indicators to numeric
rbiods[,12:21] <- apply(rbiods[,12:21],2,function(x) as.numeric(x));

############################################
## 2. Supplemental Figure 1 - Consort diagram numbers ##
############################################

#total number of participants 
length(unique(rhfds$study_id)) #349 

#number of participants with some baseline biomarkers measured 
length(rbiods$subject_id[rbiods$timepoint=="Baseline"]) #242

#missingness by biomarker at baseline
apply(rbiods[rbiods$timepoint=="Baseline",12:16],2,function(x) sum(is.na(x))) 
# creatinine_mg_dl troponin_i_pg_ml      hs_crp_mg_l       nt_pro_bnp  troponin_t_ng_l 
# 9                8                9                2                0 

#out of the 242 participants with some baseline biomarkers measured,
#how many were lost to follow up and are missing both follow-up SPPB and 6MWD measurements?
sum(rhfds$study_id %in% rbiods$subject_id[rbiods$timepoint=="Baseline"]& rhfds$lfu==1
    & is.na(rhfds$fu_sppb) & is.na(rhfds$fu_smw)) #16

#out of the 242 participants with some baseline biomarkers measured,
#how many died and are missing both follow-up SPPB and 6MWD measurements?
sum(rhfds$study_id %in% rbiods$subject_id[rbiods$timepoint=="Baseline"] & rhfds$death==1
    & is.na(rhfds$fu_sppb) & is.na(rhfds$fu_smw)) #14

#out of the 242 participants with some baseline biomarkers measured,
#how many have a follow-up SPPB or 6MWD measurement (or both)?
sum(rhfds$study_id %in% rbiods$subject_id[rbiods$timepoint=="Baseline"] 
    & (!is.na(rhfds$fu_smw) | !is.na(rhfds$fu_sppb))) #212

#out of the 242 participants with some baseline biomarkers measured,
#how many have a follow-up SPPB measurement?
sum(rhfds$study_id %in% rbiods$subject_id[rbiods$timepoint=="Baseline"] 
    & !is.na(rhfds$fu_sppb)) #212

#out of the 242 participants with some baseline biomarkers measured,
#how many have both follow-up SPPB and 6MWD measurements?
sum(rhfds$study_id %in% rbiods$subject_id[rbiods$timepoint=="Baseline"] 
    & (!is.na(rhfds$fu_smw) & !is.na(rhfds$fu_sppb))) #185

#out of the 242 participants with some baseline biomarkers measured,
#how many have both follow-up SPPB and 6MWD measurements and all 5 baseline biomarkers?

sum(rhfds$study_id %in% rbiods$subject_id[rbiods$timepoint=="Baseline"] 
    & !is.na(rhfds$fu_smw) & !is.na(rhfds$fu_sppb)
    & apply(merge(rhfds[,'study_id',drop=FALSE],
                  rbiods[rbiods$timepoint=="Baseline",1:16],
                  by.x='study_id',by.y='subject_id',all.x=TRUE)[,12:16],1,
            function(x) !anyNA(x))) #178 

#out of the 242 participants with some baseline biomarkers measured,
#how many have some follow-up biomarkers measured?
sum(rbiods$subject_id[rbiods$timepoint=="Baseline"] %in% 
      rbiods$subject_id[rbiods$timepoint=="Follow Up"]) #157

#out of the 157 with some baseline and some follow up biomarker measurements, 
#how many at follow-up were missing for each biomarker?
colSums(is.na(rbiods[rbiods$subject_id %in%
                       (rbiods$subject_id[rbiods$timepoint=="Baseline"]
                        [rbiods$subject_id[rbiods$timepoint=="Baseline"] %in% 
                            rbiods$subject_id[rbiods$timepoint=="Follow Up"]]) & 
                       rbiods$timepoint=="Follow Up",12:16]))
# creatinine_mg_dl troponin_i_pg_ml      hs_crp_mg_l       nt_pro_bnp  troponin_t_ng_l 
# 7                6                7                2                2 

############################################
## 3. Baseline biomarker outlier removal and final pre-processing ##
############################################

## Identify and remove outliers in baseline biomarker levels

## How many outliers are there for each biomarker?
## More than 3 standard deviations beyond the mean after log (base 2) transformation
apply(rbiods[rbiods$timepoint=='Baseline',12:16],2,
      function(x) sum((log2(x+0.001)>mean(log2(x+0.001),na.rm=TRUE)+3*sd(log2(x+0.001),na.rm=TRUE))|
                        (log2(x+0.001)<mean(log2(x+0.001),na.rm=TRUE)-3*sd(log2(x+0.001),na.rm=TRUE)),na.rm=TRUE))

## Which rows correspond to outliers for each biomarker?
bmoutliers <- apply(rbiods[rbiods$timepoint=='Baseline',12:16],2,
                    function(x) rbiods[rbiods$timepoint=='Baseline',c("sample_id","timepoint")]
                    [which((log2(x+0.001)>mean(log2(x+0.001),na.rm=TRUE)+
                              3*sd(log2(x+0.001),na.rm=TRUE))|
                             (log2(x+0.001)<mean(log2(x+0.001),na.rm=TRUE)-
                                3*sd(log2(x+0.001),na.rm=TRUE))),])

## Replace outlier values with NA
for (i in 1:length(bmoutliers)){
  bmtmp <- names(bmoutliers)[i]
  if (nrow(bmoutliers[[i]]>0)){
    for (j in 1:nrow(bmoutliers[[i]])){
      print(bmoutliers[[i]][j,])
      idtmp <- bmoutliers[[i]][j,1]
      rbiods[rbiods$timepoint=='Baseline' & rbiods$sample_id == idtmp,bmtmp] <- NA
    }
  }
}

## Merge rhfds and rbiods
rcombds <- merge(x=rbiods,y=rhfds,
                 by.x='subject_id',by.y='study_id',all.x=TRUE);

############################################
## 4. Table 1 (Baseline Patient Characteristics) ##
############################################
varlist = c("intervention_1_control_0","age","sex","race___4",
            "hf_cat","egfr",
            "troponin_t_ng_l","troponin_i_pg_ml","nt_pro_bnp",
            "hs_crp_mg_l","creatinine_mg_dl",
            "bl_sppb","bl_smw")

tbl_summary(rcombds[rcombds$timepoint=='Baseline',varlist],
            by=intervention_1_control_0,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)"
            ),label = list(age ~ "Age",
                           sex ~ "Sex (Female=1)",
                           race___4 ~ "Race(White=1)",
                           hf_cat ~ "HF Category (HFpEF=1) (Study)",
                           egfr ~ "EGFR (Study)",
                           troponin_t_ng_l ~ "Hs-cTnT (Inova)",
                           troponin_i_pg_ml ~ "Hs-cTnI (Inova)",
                           nt_pro_bnp ~ "NT-proBNP (Inova)",
                           hs_crp_mg_l ~ "Hs-CRP (Inova)",
                           creatinine_mg_dl ~ "Creatinine (Inova)"
            ),
            missing_text = "(Missing)"
) %>% add_p(list(all_continuous() ~ "wilcox.test",all_categorical() ~ "fisher.test")) %>% add_n(
) %>% modify_header(label ~ "**Variable**") %>% 
  modify_caption("**Table 1. Baseline Patient Characteristics**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Control (0) vs. Intervention (1)**") %>%
  as_gt() %>% 
  gtsave(filename = "table1_baselinecharacteristics.html")

############################################
## 5. Spearman correlation of baseline biomarkers ##
############################################
scmat <- round(cor(rbiods[rbiods$timepoint=="Baseline",12:16],
    method="spearman",use="pairwise.complete.obs"),2)
scmat[lower.tri(scmat)] <- NA
rownames(scmat) <- colnames(scmat) <- c("Creatinine (mg/dL)",
                                        "Hs-cTnI (ng/L)","Hs-CRP (mg/dL)",
                                        "NT-proBNP (ng/L)","Hs-cTnT (ng/L)")
scmat_melt <- melt(scmat, na.rm = TRUE)

#compute p-values
for (i in 1:nrow(scmat_melt)){
    b1 <- which(rownames(scmat)==as.character(scmat_melt$Var1[i]))
    b2 <- which(rownames(scmat)==as.character(scmat_melt$Var2[i]))
    scmat_melt$pval[i] <- cor.test(rbiods[rbiods$timepoint=="Baseline",11+b1],
                    rbiods[rbiods$timepoint=="Baseline",11+b2],
                    method="spearman",use="pairwise.complete.obs")$p.value
}
scmat_melt$sig = scmat_melt$pval<0.01
scmat_melt$text = ifelse(scmat_melt$sig==TRUE & scmat_melt$value != 1,paste0(scmat_melt$value,"*"),scmat_melt$value)
ggheatmap <- ggplot(scmat_melt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = text), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.5, 0.7),
    legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
# Print the heatmap
pdf("SpearmanCorrelation.pdf",width=5,height=5)
print(ggheatmap)
dev.off()

ggsave(ggheatmap,file='SpearmanCorrelation.eps',width=5,height=5,device="eps")

############################################
## 6. Table 2 (Linear Regression Results) ##
############################################
bm.name=c("troponin_t_ng_l",
          "troponin_i_pg_ml",
          "nt_pro_bnp",
          "hs_crp_mg_l",
          "creatinine_mg_dl")

for (i in 1:length(bm.name)){
  
  fname=paste("SPPBvs",bm.name[i],"_regressionanalysis.html",sep="")
  ttl=paste("Follow up SPPB ~ Baseline SPPB + log2(",bm.name[i],") + Intervention + More",sep="")
  print(tab_model(summary(lm(data=rcombds[rcombds$timepoint=="Baseline",],
                             formula=paste0("fu_sppb ~ bl_sppb + 
                               log2(",bm.name[i],"+0.001)+
                               intervention_1_control_0"))),
                  summary(lm(data=rcombds[rcombds$timepoint=="Baseline",],
                             formula=paste0("fu_sppb ~ bl_sppb + 
                               log2(",bm.name[i],"+0.001)*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"))),
                  title=ttl,file=fname,
                  digits=2,digits.p=4,show.obs=TRUE))
  
  fname=paste("6MWDvs",bm.name[i],"_regressionanalysis.html",sep="")
  ttl=paste("Follow up 6MWD ~ Baseline 6MWD + log2(",bm.name[i],") + Intervention + More",sep="")
  print(tab_model(summary(lm(data=rcombds[rcombds$timepoint=="Baseline",],
                             formula=paste0("fu_smw ~ bl_smw + 
                               log2(",bm.name[i],"+0.001)+
                               intervention_1_control_0"))),
                  summary(lm(data=rcombds[rcombds$timepoint=="Baseline",],
                             formula=paste0("fu_smw ~ bl_smw + 
                               log2(",bm.name[i],"+0.001)*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"))),
                  title=ttl,file=fname,
                  digits=2,digits.p=4,show.obs=TRUE))
  
}  

############################################
## 7. Figures 1 and 2 (Linear Regression Results) ##
############################################

#create long form dataset
rcombds_long_bl=pivot_longer(data=rcombds[rcombds$timepoint=='Baseline',
                                       c("subject_id","timepoint","intervention_1_control_0",
                                          bm.name,"bl_smw","fu_smw","bl_sppb","fu_sppb")],
                             cols=c('creatinine_mg_dl':'troponin_t_ng_l'),
                          names_to='Biomarker',values_to='ExpressionLevel')
rcombds_long_bl$intervention_1_control_0 <- factor(rcombds_long_bl$intervention_1_control_0,levels=c(0,1))
rcombds_long_bl$sppb_chg = rcombds_long_bl$fu_sppb - rcombds_long_bl$bl_sppb
rcombds_long_bl$smw_chg = rcombds_long_bl$fu_smw - rcombds_long_bl$bl_smw
colnames(rcombds_long_bl)[3] <- "Intervention"
rcombds_long_bl$Biomarker = factor(rcombds_long_bl$Biomarker,levels=c("creatinine_mg_dl","hs_crp_mg_l",
                                            "troponin_t_ng_l","nt_pro_bnp","troponin_i_pg_ml"))
bm.labs <- c("Creatinine (mg/dL)", "Hs-CRP (mg/dL)","Hs-cTnT (ng/L)", "NT-proBNP (ng/L)",
             "Hs-cTnI (ng/L)")
names(bm.labs) <- levels(rcombds_long_bl$Biomarker)

#SPPB
fig1 <- ggplot(data=rcombds_long_bl,mapping=aes(x=ExpressionLevel,y=sppb_chg,color=Intervention))+
  geom_smooth(method='lm')+
  geom_hline(yintercept=0,color='black',linetype='dashed')+
  facet_wrap('Biomarker',scales='free',ncol=3,nrow=2,
             labeller = labeller(Biomarker = bm.labs))+
  labs(x="Expression Level",y="Change in SPPB Score")+
  scale_x_continuous(trans='log2',labels = function(x) sprintf("%.0f", x))+
  scale_color_discrete(breaks=c(1,0),
                       labels=c("Rehabilitation\nIntervention","Attention\nControl"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85,0.25),
        legend.box.background = element_rect(color = "black"))

#6MWD
fig2 <- ggplot(data=rcombds_long_bl,mapping=aes(x=ExpressionLevel,y=smw_chg,color=Intervention))+
  geom_smooth(method='lm')+
  geom_hline(yintercept=0,color='black',linetype='dashed')+
  facet_wrap('Biomarker',scales='free',ncol=3,nrow=2,
             labeller = labeller(Biomarker = bm.labs))+
  labs(x="Expression Level",y="Change in 6MWD (meters)")+
  scale_x_continuous(trans='log2',labels = function(x) sprintf("%.0f", x))+
  scale_color_discrete(breaks=c(1,0),
                       labels=c("Rehabilitation\nIntervention","Attention\nControl"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.85,0.25),
        legend.box.background = element_rect(color = "black"))


# Print figures
pdf("Figure1_ChangeinSPPB.pdf",width=6,height=3)
print(fig1)
dev.off()
ggsave(fig1,file='Figure1_ChangeinSPPB.eps',width=6,height=3,device="eps")

pdf("Figure2_Changein6MWD.pdf",width=6,height=3)
print(fig2)
dev.off()
ggsave(fig2,file='Figure2_Changein6MWD.eps',width=6,height=3,device="eps")

############################################
## 8. Supplemental Table 1 ##
############################################
#supplemental table 1: comparing those with baseline biomarkers measured vs. no baseline biomarkers
rhfds$bm_measured = (rhfds$study_id %in% 
                       rcombds$subject_id[rcombds$timepoint=='Baseline'])
rhfds$bm_measured = factor(rhfds$bm_measured,levels=c("TRUE","FALSE"))
tbl_summary(rhfds[,c("bm_measured","age","sex","race___4","hf_cat","egfr","bl_sppb","bl_smw")],
            by=bm_measured,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)"
            ),label = list(age ~ "Age",
                           sex ~ "Sex (Female=1)",
                           race___4 ~ "Race(White=1)",
                           hf_cat ~ "HF Category (HFpEF=1) (Study)",
                           egfr ~ "EGFR (Study)",
                           bl_sppb ~ "SPPB Score (Baseline)",
                           bl_smw ~ "6MWD (meters) (Baseline)"
            ),
            missing_text = "(Missing)"
) %>% add_p(list(all_continuous() ~ "wilcox.test",all_categorical() ~ "fisher.test")) %>% add_n(
) %>% modify_header(label ~ "**Variable**") %>% 
  modify_caption("**Supplemental Table 1. Baseline Patient Characteristics**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Measured vs. Unmeasured**") %>%
  as_gt() %>% 
  gtsave(filename = "SuppTable1_measuredvsunmeasured.html")

############################################
## 9. Supplemental Table 3 ##
############################################

tbl_summary(rcombds[rcombds$timepoint=='Baseline',
                    c("intervention_1_control_0","event","rehosp","death","lfu")],
            by=intervention_1_control_0,
            type = list(all_continuous() ~ "continuous2"
            ),
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)"
            ),label = list(rehosp ~ "Rehospitalization",
                           death ~ "Death",
                           event ~ "Rehospitalization or Death",
                           lfu ~ "Lost to follow up"
            ),
            missing_text = "(Missing)"
) %>% add_p(list(all_continuous() ~ "wilcox.test",all_categorical() ~ "fisher.test")) %>% add_n(
) %>% modify_header(label ~ "**Variable**") %>% 
  modify_caption("**Supplemental Table 2**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Control (0) vs. Intervention (1)**") %>%
  as_gt() %>% 
  gtsave(filename = "SuppTable3_Secondaryoutcomes.html")

############################################
## 10. Supplemental Table 4 (Regression Results - Secondary Outcomes) ##
############################################
bm.name=c("troponin_t_ng_l",
          "troponin_i_pg_ml",
          "nt_pro_bnp",
          "hs_crp_mg_l",
          "creatinine_mg_dl")

for (i in 1:length(bm.name)){
  
  fname=paste("DeathorRehospvs",bm.name[i],"_regressionanalysis.html",sep="")
  ttl=paste("logit(Death or Rehosp) ~ log2(",bm.name[i],") + Intervention + More",sep="")
  print(tab_model(glm(data=rcombds[rcombds$timepoint=="Baseline",],
                      formula=paste0("event ~  
                               log2(",bm.name[i],"+0.001)+
                               intervention_1_control_0"), family=binomial(link="logit")),
                  glm(data=rcombds[rcombds$timepoint=="Baseline",],
                      formula=paste0("event ~  
                               log2(",bm.name[i],"+0.001)*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"), family=binomial(link="logit")),
                  title=ttl,
                  file=fname,
                  digits=4,digits.p=4,show.obs=TRUE))
  
  fname=paste("Deathvs",bm.name[i],"_regressionanalysis.html",sep="")
  ttl=paste("logit(Death) ~ log2(",bm.name[i],") + Intervention + More",sep="")
  print(tab_model(glm(data=rcombds[rcombds$timepoint=="Baseline",],
                      formula=paste0("death ~  
                               log2(",bm.name[i],"+0.001)+
                               intervention_1_control_0"), family=binomial(link="logit")),
                  glm(data=rcombds[rcombds$timepoint=="Baseline",],
                      formula=paste0("death ~  
                               log2(",bm.name[i],"+0.001)*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"), family=binomial(link="logit")),
                  title=ttl,
                  file=fname,
                  digits=4,digits.p=4,show.obs=TRUE))
  
  fname=paste("Rehospvs",bm.name[i],"_regressionanalysis.html",sep="")
  ttl=paste("logit(Rehosp) ~ log2(",bm.name[i],") + Intervention + More",sep="")
  print(tab_model(glm(data=rcombds[rcombds$timepoint=="Baseline",],
                      formula=paste0("rehosp ~  
                               log2(",bm.name[i],"+0.001)+
                               intervention_1_control_0"), family=binomial(link="logit")),
                  glm(data=rcombds[rcombds$timepoint=="Baseline",],
                      formula=paste0("rehosp ~  
                               log2(",bm.name[i],"+0.001)*intervention_1_control_0+
                                            age+sex+race___4+hf_cat*intervention_1_control_0"), family=binomial(link="logit")),
                  title=ttl,
                  file=fname,
                  digits=4,digits.p=4,show.obs=TRUE))
  
}  

############################################
## 11. Propensity matching ##
############################################

ps.df=rcombds[rcombds$timepoint=="Baseline",]
ps.df$trt=as.numeric(ps.df$intervention_1_control_0)
ps.df=ps.df[complete.cases(ps.df[,c("creatinine_mg_dl",'troponin_i_pg_ml',
                                    'hs_crp_mg_l','nt_pro_bnp','troponin_t_ng_l','bl_smw','fu_smw')]),]
#table(is.na(ps.df$fu_sppb),is.na(ps.df$fu_smw))


#propensity score for treatment assignment
ppty=glm(trt~age+sex+race___4+hf_cat
         ,family=binomial(link="logit"),data=ps.df)
prop.score=predict.glm(ppty,type="response",na.action=na.exclude)
ps.df=as.data.frame(cbind(ps.df,prop.score))
ppty.distance=match_on(ppty)

set.trt=ps.df[ps.df[, "trt"] == 1,]
set.cont=ps.df[ps.df[,"trt"]==0,]
match.set=c()
calp.sd=1
for (i in 1:dim(set.trt)[1]){
  pair.set1=c()
  var.select1=ifelse(abs(log2(set.cont$creatinine_mg_dl)-
                           log2(set.trt$creatinine_mg_dl[i]))<
                       calp.sd*sd(log2(set.trt$creatinine_mg_dl)),1,0)
  var.select2=ifelse(abs(log2(set.cont$troponin_i_pg_ml)-
                           log2(set.trt$troponin_i_pg_ml[i]))<
                       calp.sd*sd(log2(set.trt$troponin_i_pg_ml)),1,0)
  var.select3=ifelse(abs(log2(set.cont$troponin_t_ng_l)-
                           log2(set.trt$troponin_t_ng_l[i]))<
                       calp.sd*sd(log2(set.trt$troponin_t_ng_l)),1,0)
  var.select4=ifelse(abs(log2(set.cont$nt_pro_bnp)-
                           log2(set.trt$nt_pro_bnp[i]))<
                       calp.sd*sd(log2(set.trt$nt_pro_bnp)),1,0)
  var.select5=ifelse(abs(log2(set.cont$hs_crp_mg_l)-
                           log2(set.trt$hs_crp_mg_l[i]))<
                       calp.sd*sd(log2(set.trt$hs_crp_mg_l)),1,0)
  var.selectf=ifelse(var.select1+var.select2+var.select3+
                       var.select4+var.select5==5,1,0)
  id=c(1:length(var.select1))
  lenset.cont2=as.data.frame(cbind(
    id,set.cont,var.selectf))
  lenset.cont3=lenset.cont2[lenset.cont2[,"var.selectf"]==1,]
  if (dim(lenset.cont3)[1]>0){
    cont.index=which.min(abs(lenset.cont3$prop.score-set.trt$prop.score[i]))
    pair.set=rbind(set.trt[i,],
                   lenset.cont3[cont.index,c(-1,-dim(lenset.cont3)[2])])
    pairin=rep(i,2)
    pair.set1=cbind(pair.set,pairin)
    set.cont = set.cont[!set.cont$prop.score == 
                          lenset.cont3$prop.score[cont.index],]
  } else {
    pair.set1=c()
    set.cont = set.cont
  }
  match.set=rbind(match.set,pair.set1)
}


######################
## 12. Matching Tree Functions
######################

#tree pruning and selection functions from Matching Tree paper (Zhang et al., 2021)
#https://doi.org/10.1016/j.csda.2021.107188

#MT (m)------
#prune tree function
parent <- function(x) {
  if (x[1] != 1)
    c(Recall(if (x %% 2 == 0L) x / 2 else (x - 1) / 2), x) else x
}

jesse.tree.select=function(input,set1.caliper.tree){
  fit.overlap=TRUE
  ptree1=input
  while (fit.overlap==TRUE){
    if (dim(ptree1$frame)[1]==1) break
    rownames=as.numeric(rownames(ptree1$frame))
    size=dim(ptree1$frame)[1]
    rownamesmax=rownames[which.max(rownames)]
    rownamessecond=rownames[which(rownames == 
                                    sort(unique(rownames),partial=size-1)[size-1])]
    
    subset1=set1.caliper.tree[ptree1$where==which(rownames==rownamesmax),]
    subset2=set1.caliper.tree[ptree1$where==which(rownames==rownamessecond),]
    dummy=c(rep("1",dim(subset1)[1]),rep("2",dim(subset2)[1]))
    subsetf=cbind(rbind(subset1,subset2),dummy)
    
    subsetf$dummy=as.factor(subsetf$dummy)
    subsetf$trt=as.factor(subsetf$trt)
    ls.lm=lm(outcome~dummy*trt,data=subsetf)
    ls.list=lsmeans(ls.lm,list(pairwise~dummy|trt,pairwise~trt|dummy))
    out1=as.matrix(summary(ls.list[[4]]))
    
    q.value=1.96 #could change q.value for different methods
    CI1=c(as.numeric(out1[1,3])-q.value*as.numeric(out1[1,4]),
          as.numeric(out1[1,3])+q.value*as.numeric(out1[1,4]))
    CI2=c(as.numeric(out1[2,3])-q.value*as.numeric(out1[2,4]),
          as.numeric(out1[2,3])+q.value*as.numeric(out1[2,4]))  
    fit.overlap=max(CI1[1],CI2[1])<min(CI1[2],CI2[2])
    fit.overlap[is.na(fit.overlap)] <- FALSE
    if (fit.overlap==FALSE) break
    ptree1=snip.rpart(ptree1,
                      toss=tail(head(parent(rownamesmax), -1),n=1))}
  if (dim(ptree1$frame)[1]>1){
    return(ptree1)
  } else {
    return (input)
  }
}


######################
## 13. SPPB Matching Tree
######################

#calculate the difference between treated and control
treat.set=match.set[match.set[,"trt"]==1,]
control.set=match.set[match.set[,"trt"]==0,]
outcome.diff=(treat.set$fu_sppb-treat.set$bl_sppb)-
  (control.set$fu_sppb-control.set$bl_sppb)
set.total=as.data.frame(rbind(cbind(treat.set,outcome.diff),cbind(control.set,outcome.diff)))

#build matching tree using all data
tree.rpart=rpart(outcome.diff~creatinine_mg_dl+troponin_i_pg_ml+hs_crp_mg_l+
                   nt_pro_bnp+troponin_t_ng_l ,method="anova",data=set.total,
                 control=rpart.control(minbucket = 20))
#this is the result after prune
set.total$outcome = set.total$fu_smw - set.total$bl_smw
tree.all=jesse.tree.select(tree.rpart,set.total)

pdf("Matchingtree_SPPB.pdf")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()
setEPS()
postscript("Matchingtree_SPPB.eps")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()

#get confidence intervals
set.total$sppb_node=rpart.predict.leaves(tree.all,set.total,type="where")
modscore.sppb=lm(data=set.total,formula= outcome.diff~0+as.factor(sppb_node))
confint(modscore.sppb)
# 2.5 %   97.5 %
# as.factor(sppb_node)2 0.1271625 1.551409
# as.factor(sppb_node)4 0.5409918 2.763356
# as.factor(sppb_node)5 2.0441062 3.750766

summary(modscore.sppb)
# Call:
#   lm(formula = outcome.diff ~ 0 + as.factor(sppb_node), data = set.total)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.6522 -1.8393 -0.8393  2.1026  7.1607 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# as.factor(sppb_node)2   0.8393     0.3595   2.335  0.02130 *  
#   as.factor(sppb_node)4   1.6522     0.5610   2.945  0.00391 ** 
#   as.factor(sppb_node)5   2.8974     0.4308   6.726 7.13e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.69 on 115 degrees of freedom
# Multiple R-squared:  0.3404,	Adjusted R-squared:  0.3232 
# F-statistic: 19.79 on 3 and 115 DF,  p-value: 2.066e-10

#Leave one out stability analysis
set.seed(67)
folds <- cut(seq(1,nrow(treat.set)),breaks=nrow(treat.set),labels=FALSE)
#shuffle folds
folds<- folds[sample(nrow(treat.set))]

#store trees and subpops
tree.2=list()
tree.plot=list()
tree.leaves=list()
tree.mod.train=list()
tree.mod.test=list()
for (i in 1:max(folds)){
  
  #identify pairs in the train
  ids <- treat.set$subject_id[folds != i]
  pairs.train <- set.total$pairin[set.total$subject_id %in% ids]
  
  #this is the result from CART (update with different rpart.control params)
  tmp=rpart(outcome.diff~creatinine_mg_dl+troponin_i_pg_ml+hs_crp_mg_l+
                     nt_pro_bnp+troponin_t_ng_l,
                   method="anova",data=set.total[set.total$pairin %in% pairs.train,],
                   control=rpart.control(minbucket = 20))

  #this is the result after prune
  set.total$outcome = set.total$fu_smw - set.total$bl_smw
  tree.2[[i]]=jesse.tree.select(tmp,set.total)
  tree.plot[[i]]=rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
  
  #numbers are average outcome diff
  #get leaves for all observations
  tree.leaves[[i]]=rpart.predict.leaves(tree.2[[i]],set.total,type="where")
  
  df=cbind(set.total,tree.leaves[[i]])
  colnames(df)
  tree.mod.train[[i]]=summary(lm(outcome.diff ~ 0 + as.factor(`tree.leaves[[i]]`), data=df[df$pairin %in% pairs.train,]))

}

pdf("Matchingtree_SPPB_LOO.pdf", onefile = TRUE)
for (i in 1:length(tree.2)) {
  rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
}
dev.off()

######################
## 14. 6MWD Matching Tree
######################

#calculate the difference between treated and control
treat.set=match.set[match.set[,"trt"]==1,]
control.set=match.set[match.set[,"trt"]==0,]
outcome.diff=(treat.set$fu_smw-treat.set$bl_smw)-
  (control.set$fu_smw-control.set$bl_smw)
set.total=as.data.frame(rbind(cbind(treat.set,outcome.diff),cbind(control.set,outcome.diff)))

#build matching tree using all data
tree.rpart=rpart(outcome.diff~creatinine_mg_dl+troponin_i_pg_ml+hs_crp_mg_l+
                   nt_pro_bnp+troponin_t_ng_l ,method="anova",data=set.total,
                 control=rpart.control(minbucket = 20))
#this is the result after prune
set.total$outcome = set.total$fu_smw - set.total$bl_smw
tree.all=jesse.tree.select(tree.rpart,set.total)

pdf("Matchingtree_6MWD.pdf")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()
setEPS()
postscript("Matchingtree_6MWD.eps")
rpart.plot(tree.all,box.palette = 'Reds',extra=1)
dev.off()

#get confidence intervals
set.total$smw_node=rpart.predict.leaves(tree.all,set.total,type="where")
modscore.6mwd=lm(data=set.total,formula= outcome.diff~0+as.factor(smw_node))
confint(modscore.6mwd)
# 2.5 %    97.5 %
# as.factor(smw_node)2 -0.4550698  60.19986
# as.factor(smw_node)4 20.4393863  90.47766
# as.factor(smw_node)5 78.0548485 179.54995

summary(modscore.6mwd)
# Call:
#   lm(formula = outcome.diff ~ 0 + as.factor(smw_node), data = set.total)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -192.01  -74.02   11.63   54.96  499.82 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# as.factor(smw_node)2    29.87      15.31   1.951  0.05348 .  
# as.factor(smw_node)4    55.46      17.68   3.137  0.00217 ** 
#   as.factor(smw_node)5   128.80      25.62   5.027 1.84e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 114.6 on 115 degrees of freedom
# Multiple R-squared:  0.2529,	Adjusted R-squared:  0.2334 
# F-statistic: 12.97 on 3 and 115 DF,  p-value: 2.33e-07

#Leave one out stability analysis
set.seed(34)
folds <- cut(seq(1,nrow(treat.set)),breaks=nrow(treat.set),labels=FALSE)
#shuffle folds
folds<- folds[sample(nrow(treat.set))]

#store trees and subpops
tree.2=list()
tree.plot=list()
tree.leaves=list()
tree.mod.train=list()
tree.mod.test=list()
for (i in 1:max(folds)){
  
  #identify pairs in the train
  ids <- treat.set$subject_id[folds != i]
  pairs.train <- set.total$pairin[set.total$subject_id %in% ids]
  
  #this is the result from CART (update with different rpart.control params)
  tmp=rpart(outcome.diff~creatinine_mg_dl+troponin_i_pg_ml+hs_crp_mg_l+
              nt_pro_bnp+troponin_t_ng_l,
            method="anova",data=set.total[set.total$pairin %in% pairs.train,],
            control=rpart.control(minbucket = 20))
  
  #this is the result after prune
  set.total$outcome = set.total$fu_smw - set.total$bl_smw
  tree.2[[i]]=jesse.tree.select(tmp,set.total)
  tree.plot[[i]]=rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
  
  #numbers are average outcome diff
  #get leaves for all observations
  tree.leaves[[i]]=rpart.predict.leaves(tree.2[[i]],set.total,type="where")
  
  df=cbind(set.total,tree.leaves[[i]])
  colnames(df)
  tree.mod.train[[i]]=summary(lm(outcome.diff ~ 0 + as.factor(`tree.leaves[[i]]`), data=df[df$pairin %in% pairs.train,]))
  
}

pdf("Matchingtree_6MWD_LOO.pdf", onefile = TRUE)
for (i in 1:length(tree.2)) {
  rpart.plot(tree.2[[i]],box.palette = 'Reds',extra=1)
}
dev.off()

######################
## 15. Supplemental table 2 matched vs. unmatched
######################
varlist=c("intervention_1_control_0","age","sex","race___4",
          "hf_cat","egfr","troponin_t_ng_l","troponin_i_pg_ml","nt_pro_bnp","hs_crp_mg_l",
          "creatinine_mg_dl")
tmpdf = rbind(data.frame(set.total[,varlist],dataset="matched"),
              data.frame(setdiff(ps.df[,varlist],set.total[,varlist]),dataset="unmatched"))
tmpdf$dataset=factor(tmpdf$dataset,levels=c("unmatched","matched"))
tmpdf$grp=str_c(tmpdf$intervention_1_control_0,tmpdf$dataset)
tmpdf$grp=as.factor(tmpdf$grp)
tmpdf$grp=fct_collapse(tmpdf$grp, Unmatched = c("0unmatched","1unmatched"), MatchedControl = "0matched",MatchedIntervention = "1matched")
tmpdf$grp = factor(tmpdf$grp,levels = c("MatchedControl","MatchedIntervention","Unmatched"))
tmpdf=tmpdf[,-which(names(tmpdf) %in% c("intervention_1_control_0","dataset"))]

tbl_summary(tmpdf,
            by=grp,
            type = all_continuous() ~ "continuous2",
            statistic = list(all_continuous() ~ c("{mean} ({sd})",
                                                  "{median} ({p25}, {p75})", 
                                                  "{min}, {max}"),
                             all_categorical() ~ "{n} / {N} ({p}%)"
            ),label = list(age ~ "Age",
                           sex ~ "Sex (Female=1)",
                           race___4 ~ "Race(White=1)",
                           hf_cat ~ "HF Category (HFpEF=1) (Study)",
                           egfr ~ "EGFR (Study)",
                           troponin_t_ng_l ~ "Hs-cTnT (Inova)",
                           troponin_i_pg_ml ~ "Hs-cTnI (Inova)",
                           nt_pro_bnp ~ "NT-proBNP (Inova)",
                           hs_crp_mg_l ~ "Hs-CRP (Inova)",
                           creatinine_mg_dl ~ "Creatinine (Inova)"
            ),
            missing_text = "(Missing)"
)  %>% modify_header(label ~ "**Variable**") %>% 
  modify_caption("**Supplemental Table 2. Baseline Measurements (Unmatched vs. Matched)**") %>%
  as_gt() %>% 
  gtsave(filename = "Supptable2_matchedvsunmatched.html")

######################
## 16. Supplemental figure 2 
######################
#create long form dataset
rcombds_long_fu=pivot_longer(data=rcombds[rcombds$timepoint=='Follow Up',
                                          c("subject_id","timepoint","intervention_1_control_0",
                                            bm.name,"bl_smw","fu_smw","bl_sppb","fu_sppb")],
                             cols=c('creatinine_mg_dl':'troponin_t_ng_l'),
                             names_to='Biomarker',values_to='ExpressionLevel')
rcombds_long_fu$intervention_1_control_0 <- factor(rcombds_long_fu$intervention_1_control_0,levels=c(0,1))
rcombds_long_fu$sppb_chg = rcombds_long_fu$fu_sppb - rcombds_long_fu$bl_sppb
rcombds_long_fu$smw_chg = rcombds_long_fu$fu_smw - rcombds_long_fu$bl_smw
colnames(rcombds_long_fu)[3] <- "Intervention"
rcombds_long_fu$Biomarker = factor(rcombds_long_fu$Biomarker,levels=c("creatinine_mg_dl","hs_crp_mg_l",
                                                                      "troponin_t_ng_l","nt_pro_bnp","troponin_i_pg_ml"))
bm.labs <- c("Creatinine (mg/dL)", "Hs-CRP (mg/dL)","Hs-cTnT (ng/L)", "NT-proBNP (ng/L)",
             "Hs-cTnI (ng/L)")
names(bm.labs) <- levels(rcombds_long_fu$Biomarker)

rcombds_wide=merge(rcombds_long_bl,
                   rcombds_long_fu,by=c('subject_id','Biomarker','Intervention',"bl_smw","fu_smw","bl_sppb","fu_sppb","sppb_chg","smw_chg"),
                   all.x=TRUE)
rcombds_wide$Biomarker_change = log2(rcombds_wide$ExpressionLevel.y+0.001) - log2(rcombds_wide$ExpressionLevel.x+0.001)

supfig2 <- ggplot(data=rcombds_wide,
             mapping=aes(x=Intervention, y=ExpressionLevel.y/ExpressionLevel.x, fill=Intervention))+
        geom_violin()+facet_wrap_paginate('Biomarker',scales='free',ncol=2,nrow=3,
                                          page=1,labeller = labeller(Biomarker = bm.labs))+hw+
        labs(y="Ratio of change from baseline to 12-weeks")+
        geom_hline(yintercept=1,linetype=2)+
        scale_y_continuous(trans='log2',labels = function(x) sprintf("%.3f", x))+
        scale_x_discrete(breaks=c(1,0),
                         labels=c("RI","AC"))+
        scale_fill_discrete(breaks=c(0,1),
                            labels=c("Attention Control (AC)","Rehabilitation Intervention (RI)"))+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_text(face="bold"),
              legend.title = element_blank(),
              legend.position = c(0.75,0.15),
              legend.box.background = element_rect(color = "black"))


pdf("SuppFigure2_biomarkerchange.pdf",height=9,width=7,onefile = TRUE)
supfig2
dev.off()
ggsave(supfig2,file='SuppFigure2_biomarkerchange.eps',height=9,width=7,device="eps")


