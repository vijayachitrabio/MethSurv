#load the required libararies
library(survMisc)
library(data.table)
library(shinythemes)
library(shiny)
library(shinyBS)
library(ggplot2)
library(base)
library(grid)
library(pheatmap)
path = "/srv/shiny-server/"

###load the necessary files####
load(paste0(path, "LUAD_meth.RData"))
LUAD_meth=as.data.table(LUAD_meth)
load(paste0(path, "LUAD_clinic_new.RData"))
load(paste0(path,"LUAD_annot.RData"))

load(paste0(path,"UCS_meth.RData"))
UCS_meth=as.data.table(UCS_meth)
load(paste0(path,"UCS_clinic_new.RData"))
load(paste0(path,"UCS_annot.RData"))

load(paste0(path,"UCEC_meth.RData"))
UCEC_meth=as.data.table(UCEC_meth)
load(paste0(path,"UCEC_clinic_new.RData"))
load(paste0(path,"UCEC_annot.RData"))

load(paste0(path, "MESO_meth.RData"))
MESO_meth=as.data.table(MESO_meth)
load(paste0(path, "MESO_clinic_new.RData"))
load(paste0(path, "MESO_annot.RData"))

load(paste0(path, "BRCA_meth.RData"))
BRCA_meth=as.data.table(BRCA_meth)
load(paste0(path, "BRCA_clinic_new.RData"))
load(paste0(path, "BRCA_annot.RData"))

load(paste0(path, "BLCA_meth.RData"))
BLCA_meth=as.data.table(BLCA_meth)
load(paste0(path, "BLCA_clinic_new.RData"))
load(paste0(path, "BLCA_annot.RData"))

load(paste0(path, "ACC_meth.RData"))
ACC_meth=as.data.table(ACC_meth)
load(paste0(path, "ACC_clinic_new.RData"))
load(paste0(path, "ACC_annot.RData"))


load(paste0(path, "CESC_meth.RData"))
CESC_meth=as.data.table(CESC_meth)
load(paste0(path, "CESC_clinic_new.RData"))
load(paste0(path, "CESC_annot.RData"))

load(paste0(path, "COAD_meth.RData"))
COAD_meth=as.data.table(COAD_meth)
load(paste0(path, "COAD_clinic_new.RData"))
load(paste0(path, "COAD_annot.RData"))

load(paste0(path, "ESCA_meth.RData"))
ESCA_meth=as.data.table(ESCA_meth)
load(paste0(path, "ESCA_clinic_new.RData"))
load(paste0(path, "ESCA_annot.RData"))


load(paste0(path, "LGG_meth.RData"))
LGG_meth=as.data.table(LGG_meth)
load(paste0(path, "LGG_clinic_new.RData"))
load(paste0(path,"LGG_annot.RData"))


load(paste0(path, "PAAD_meth.RData"))
PAAD_meth=as.data.table(PAAD_meth)
load(paste0(path, "PAAD_clinic_new.RData"))
load(paste0(path,"PAAD_annot.RData"))


#load(paste0(path, "PCPG_meth.RData"))
#load(paste0(path, "PCPG_clinic_new.RData"))
#load(paste0(path,"PCPG_annot.RData"))

#load(paste0(path, "PRAD_meth.RData"))
#load(paste0(path, "PRAD_clinic_new.RData"))
#load(paste0(path,"PRAD_annot.RData"))

load(paste0(path, "GBM_meth.RData"))
GBM_meth=as.data.table(GBM_meth)
load(paste0(path, "GBM_clinic_new.RData"))
load(paste0(path,"GBM_annot.RData"))

load(paste0(path, "READ_meth.RData"))
READ_meth=as.data.table(READ_meth)
load(paste0(path, "READ_clinic_new.RData"))
load(paste0(path,"READ_annot.RData"))

load(paste0(path, "HNSC_meth.RData"))
HNSC_meth=as.data.table(HNSC_meth)
load(paste0(path, "HNSC_clinic_new.RData"))
load(paste0(path,"HNSC_annot.RData"))

#load(paste0(path, "TGCT_meth.RData"))
#load(paste0(path, "TGCT_clinic_new.RData"))
#load(paste0(path,"TGCT_annot.RData"))

load(paste0(path, "SKCM_meth.RData"))
SKCM_meth=as.data.table(SKCM_meth)
load(paste0(path, "SKCM_clinic_new.RData"))
load(paste0(path,"SKCM_annot.RData"))


load(paste0(path, "SARC_meth.RData"))
SARC_meth=as.data.table(SARC_meth)
load(paste0(path, "SARC_clinic_new.RData"))
load(paste0(path,"SARC_annot.RData"))





load(paste0(path, "LIHC_meth.RData"))
LIHC_meth=as.data.table(LIHC_meth)
load(paste0(path, "LIHC_clinic_new.RData"))
load(paste0(path,"LIHC_annot.RData"))


load(paste0(path, "LUSC_meth.RData"))
LUSC_meth=as.data.table(LUSC_meth)
load(paste0(path, "LUSC_clinic_new.RData"))
load(paste0(path,"LUSC_annot.RData"))

load(paste0(path, "STAD_meth.RData"))
STAD_meth=as.data.table(STAD_meth)
load(paste0(path, "STAD_clinic_new.RData"))
load(paste0(path,"STAD_annot.RData"))


#load(paste0(path, "THCA_meth.RData"))
#load(paste0(path, "THCA_clinic_new.RData"))
#load(paste0(path,"THCA_annot.RData"))

#load(paste0(path, "THYM_meth.RData"))
#load(paste0(path, "THYM_clinic_new.RData"))
#load(paste0(path,"THYM_annot.RData"))

load(paste0(path, "KICH_meth.RData"))
KICH_meth=as.data.table(KICH_meth)
load(paste0(path, "KICH_clinic_new.RData"))
load(paste0(path,"KICH_annot.RData"))

load(paste0(path, "KIRP_meth.RData"))
KIRP_meth=as.data.table(KIRP_meth)
load(paste0(path, "KIRP_clinic_new.RData"))
load(paste0(path,"KIRP_annot.RData"))


load(paste0(path, "KIRC_meth.RData"))
KIRC_meth=as.data.table(KIRC_meth)
load(paste0(path, "KIRC_clinic_new.RData"))
load(paste0(path,"KIRC_annot.RData"))

load(paste0(path, "UVM_meth.RData"))
UVM_meth=as.data.table(UVM_meth)
load(paste0(path, "UVM_clinic_new.RData"))
load(paste0(path,"UVM_annot.RData"))


load(paste0(path, "LAML_meth.RData"))
LAML_meth=as.data.table(LAML_meth)
load(paste0(path, "LAML_clinic_new.RData"))
load(paste0(path,"LAML_annot.RData"))
UCS_meth=na.omit(UCS_meth)
UCS_meth=subset(UCS_meth,UCSC_RefGene_Name!="NA")
#row.names(UCS_meth)=UCS_meth$Name 

LUAD_meth=na.omit(LUAD_meth)
LUAD_meth=subset(LUAD_meth,UCSC_RefGene_Name!="NA")
#row.names(LUAD_meth)=as.character(LUAD_meth$Name) 

UCEC_meth=na.omit(UCEC_meth)
UCEC_meth=subset(UCEC_meth,UCSC_RefGene_Name!="NA")
#row.names(UCEC_meth)=as.character(UCEC_meth$Name)

MESO_meth=na.omit(MESO_meth)
MESO_meth=subset(MESO_meth,UCSC_RefGene_Name!="NA")
#row.names(MESO_meth)=MESO_meth$Name 


BRCA_meth=na.omit(BRCA_meth)
BRCA_meth=subset(BRCA_meth,UCSC_RefGene_Name!="NA")
#row.names(BRCA_meth)=BRCA_meth$Name 


BLCA_meth=na.omit(BLCA_meth)
BLCA_meth=subset(BLCA_meth,UCSC_RefGene_Name!="NA")
#row.names(BLCA_meth)=BLCA_meth$Name

ACC_meth=na.omit(ACC_meth)
ACC_meth=subset(ACC_meth,UCSC_RefGene_Name!="NA")
#row.names(ACC_meth)=ACC_meth$Name

CESC_meth=na.omit(CESC_meth)
CESC_meth=subset(CESC_meth,UCSC_RefGene_Name!="NA")
#row.names(CESC_meth)=CESC_meth$Name

COAD_meth=na.omit(COAD_meth)
COAD_meth=subset(COAD_meth,UCSC_RefGene_Name!="NA")
#row.names(COAD_meth)=COAD_meth$Name


ESCA_meth=na.omit(ESCA_meth)
ESCA_meth=subset(ESCA_meth,UCSC_RefGene_Name!="NA")
#row.names(ESCA_meth)=ESCA_meth$Name


LGG_meth=na.omit(LGG_meth)
LGG_meth=subset(LGG_meth,UCSC_RefGene_Name!="NA")
#row.names(LGG_meth)=as.character(LGG_meth$Name) 

PAAD_meth=na.omit(PAAD_meth)
PAAD_meth=subset(PAAD_meth,UCSC_RefGene_Name!="NA")
#row.names(PAAD_meth)=as.character(PAAD_meth$Name) 



#PCPG_meth=na.omit(PCPG_meth)
#PCPG_meth=subset(PCPG_meth,UCSC_RefGene_Name!="NA")
##row.names(PCPG_meth)=as.character(PCPG_meth$Name) 


#PRAD_meth=na.omit(PRAD_meth)
#PRAD_meth=subset(PRAD_meth,UCSC_RefGene_Name!="NA")
##row.names(PRAD_meth)=as.character(PRAD_meth$Name) 


GBM_meth=na.omit(GBM_meth)
GBM_meth=subset(GBM_meth,UCSC_RefGene_Name!="NA")
#row.names(GBM_meth)=as.character(GBM_meth$Name) 


READ_meth=na.omit(READ_meth)
READ_meth=subset(READ_meth,UCSC_RefGene_Name!="NA")
#row.names(READ_meth)=as.character(READ_meth$Name) 

HNSC_meth=na.omit(HNSC_meth)
HNSC_meth=subset(HNSC_meth,UCSC_RefGene_Name!="NA")
#row.names(HNSC_meth)=as.character(HNSC_meth$Name) 

#TGCT_meth=na.omit(TGCT_meth)
#TGCT_meth=subset(TGCT_meth,UCSC_RefGene_Name!="NA")
##row.names(TGCT_meth)=as.character(TGCT_meth$Name) 


SKCM_meth=na.omit(SKCM_meth)
SKCM_meth=subset(SKCM_meth,UCSC_RefGene_Name!="NA")
#row.names(SKCM_meth)=as.character(SKCM_meth$Name) 


SARC_meth=na.omit(SARC_meth)
SARC_meth=subset(SARC_meth,UCSC_RefGene_Name!="NA")
#row.names(SARC_meth)=as.character(SARC_meth$Name) 


LIHC_meth=na.omit(LIHC_meth)
LIHC_meth=subset(LIHC_meth,UCSC_RefGene_Name!="NA")
#row.names(LIHC_meth)=as.character(LIHC_meth$Name) 

LUSC_meth=na.omit(LUSC_meth)
LUSC_meth=subset(LUSC_meth,UCSC_RefGene_Name!="NA")
#row.names(LUSC_meth)=as.character(LUSC_meth$Name) 

STAD_meth=na.omit(STAD_meth)
STAD_meth=subset(STAD_meth,UCSC_RefGene_Name!="NA")
#row.names(STAD_meth)=as.character(STAD_meth$Name) 

#THCA_meth=na.omit(THCA_meth)
#THCA_meth=subset(THCA_meth,UCSC_RefGene_Name!="NA")
##row.names(THCA_meth)=as.character(THCA_meth$Name) 

#THYM_meth=na.omit(THYM_meth)
#THYM_meth=subset(THYM_meth,UCSC_RefGene_Name!="NA")
##row.names(THYM_meth)=as.character(THYM_meth$Name) 

KICH_meth=na.omit(KICH_meth)
KICH_meth=subset(KICH_meth,UCSC_RefGene_Name!="NA")
#row.names(KICH_meth)=as.character(KICH_meth$Name) 


KIRP_meth=na.omit(KIRP_meth)
KIRP_meth=subset(KIRP_meth,UCSC_RefGene_Name!="NA")
#row.names(KIRP_meth)=as.character(KIRP_meth$Name) 


KIRC_meth=na.omit(KIRC_meth)
KIRC_meth=subset(KIRC_meth,UCSC_RefGene_Name!="NA")
#row.names(KIRC_meth)=as.character(KIRC_meth$Name) 


UVM_meth=na.omit(UVM_meth)
UVM_meth=subset(UVM_meth,UCSC_RefGene_Name!="NA")
#row.names(UVM_meth)=as.character(UVM_meth$Name) 


LAML_meth=na.omit(LAML_meth)
LAML_meth=subset(LAML_meth,UCSC_RefGene_Name!="NA")
#row.names(LAML_meth)=as.character(LAML_meth$Name)

#load top biomarkers files

load("ACC_top_2017.RData")
load("BLCA_top_2017.RData")
load("BRCA_top_2017.RData")
load("CESC_top_2017.RData")
load("COAD_top_2017.RData")
load("ESCA_top_2017.RData")
load("GBM_top_2017.RData")
load("HNSC_top_2017.RData")
load("KICH_top_2017.RData")
load("KIRC_top_2017.RData")
load("KIRP_top_2017.RData")
load("LAML_top_2017.RData")
load("LGG_top_2017.RData")
load("LIHC_top_2017.RData")
load("LUAD_top_2017.RData")
load("LUSC_top_2017.RData")
load("MESO_top_2017.RData")
load("PAAD_top_2017.RData")
load("READ_top_2017.RData")
load("SARC_top_2017.RData")
load("SKCM_top_2017.RData")
load("STAD_top_2017.RData")
load("UCEC_top_2017.RData")
load("UCS_top_2017.RData")
load("UVM_top_2017.RData")

# assign the available covariates for each cancer type

UCS_covariate= c("BMI","age","stage")
LUAD_covariate=c("age","sex","stage")
UCEC_covariate= c("BMI","age","stage")
MESO_covariate=c("age","sex","stage")
BRCA_covariate=c("age","her2_status_by_ihc","er_status_by_ihc","pr_status_by_ihc","stage")
BLCA_covariate= c("BMI","age","sex","grade")
CESC_covariate= c("BMI","age")
COAD_covariate= c("BMI","age","sex")
ESCA_covariate= c("BMI","age", "sex","stage","grade")
ACC_covariate=c("age", "sex","stage")
GBM_covariate=c("age", "sex")
HNSC_covariate=c("age","sex","stage")
LGG_covariate=c("age", "sex","grade")
PAAD_covariate=c("age","sex","stage")
SARC_covariate=c("age", "sex")
SKCM_covariate=c("age", "sex")
READ_covariate=c("BMI","age","sex")
LIHC_covariate= c("BMI","age","sex","stage","grade")
LUSC_covariate=c("age","sex","stage")
STAD_covariate=c("age","sex","stage")
KICH_covariate=c("age", "sex")
KIRP_covariate= c("BMI","age","sex","stage")
KIRC_covariate=c("age","sex","stage")
UVM_covariate= c("BMI","age","sex","pathologic_N")
LAML_covariate=c("age", "sex")

#how the tool tip should be positioned

tooltipplace="top"

###survival plot function
cox_model_best_split_plot<-function(cancer_data1,cancer_data2,gene,cov_vector=NULL,threshold=0.25,region,group,probe,split){
  
  #set user selected gene, probe and regions as current gene data and extract it
  current_gene_data=cancer_data1[which( cancer_data1$UCSC_RefGene_Name %in% gene & cancer_data1$UCSC_RefGene_Group %in% region & cancer_data1$Relation_to_UCSC_CpG_Island %in% group & cancer_data1$Name %in% probe) , ]
 
  current_gene_data=current_gene_data[,!colnames(current_gene_data) %in% c("Name", "CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"),with=FALSE]
  methylation<-as.numeric(current_gene_data)
  time_to_event<-as.numeric(cancer_data2$time_to_event)#time up to death from data 
  event<-as.numeric(cancer_data2$event)
  age<-as.numeric(cancer_data2$age)
  BMI<-as.numeric(cancer_data2$BMI)
sex=cancer_data2$sex
grade=cancer_data2$grade
stage=cancer_data2$clinical_stage
her2_status_by_ihc=cancer_data2$her2_status_by_ihc
er_status_by_ihc=cancer_data2$er_status_by_ihc
pr_status_by_ihc=cancer_data2$pr_status_by_ihc
pathologic_N=cancer_data2$pathologic_N
  
  # optimal splitpoint by log-rank test
  
  covariates<-c("methylation",cov_vector)
  
  t_curr=data.frame(t(current_gene_data))   
 colnames(t_curr)[1]="meth"
 t_curr$event=cancer_data2$event
    t_curr$time_to_event=cancer_data2$time_to_event
    t_curr=data.frame(t_curr)
    
    res.cut <- surv_cutpoint(t_curr, time = "time_to_event", event = "event",
                         variables = c("meth" ))
                         
 maxstat_split<- res.cut$cutpoint[,1]
                         
  
  ### splitting by mean, median, upper and lower quantiles
  
  mean<-factor(ifelse(methylation<=mean(methylation),"hypo","hyper")) 
  mean<-relevel(mean, ref = "hypo")
  median<-factor(ifelse(methylation<=median(methylation),"hypo","hyper")) 
  median<-relevel(median, ref = "hypo")
  q25<-factor(ifelse(methylation<=quantile(methylation,0.25),"hypo","hyper")) 
  q25<-relevel(q25, ref = "hypo")
  q75<-factor(ifelse(methylation<=quantile(methylation,0.75),"hypo","hyper")) 
  q75<-relevel(q75, ref = "hypo")
  maxstat<-factor(ifelse(methylation<=maxstat_split,"hypo","hyper"))
  maxstat<-relevel(maxstat, ref = "hypo")
  
  #split_df<-data.frame(mean,median,q25,q75)
  split_list<-list(mean,median,q25,q75,maxstat)
  names(split_list)<-c("mean","median","q25","q75","maxstat")
    
    
    
  split_threshold<-function(x, threshold){
    t1<-table(x)#how many hyper and hypo in 1 type of variable
    does_suit<-"yes"
    if(min(t1)/sum(t1)<threshold){does_suit<-"no"}
    return(does_suit)
  }
  
  #does_suit_threshold<-apply(split_df,2,split_threshold,threshold=threshold)
  does_suit_threshold<-sapply(split_list,split_threshold,threshold=threshold)
  split_list_yes<-split_list[which(does_suit_threshold %in% "yes")]
  
  #meth_v<-split_list_yes[[1]]
  
  HR_fun0<-function(meth_v,te,e,covs){ 
    surv_fun<-Surv(time = as.numeric(te), event = as.numeric(e)) # survival function part of cox model
    cov_vector1<-c("meth_v",covs)#add the covariates
    if(is.null(covs)){ 
      KM_formula_1<-as.formula(surv_fun ~ meth_v)#if 0 then regress only on methylation
    }
    else{
      KM_formula_1<-as.formula(paste("surv_fun ~ ", paste(cov_vector1, collapse= "+")))#other collapse with other covariates
    }
    
    KM_for_HR<-coxph(KM_formula_1)#formula with the conditions, execute the cox model
    result<-summary(KM_for_HR)$coefficients
    res<-summary(KM_for_HR)
    meth_ind<-which(rownames(result) %in% "meth_vhyper")#estimate for methylation releveling hyper 
    HR<-round(result[meth_ind,2],3)#second coloumn  from res1
    LR_test_pvalue<-as.numeric(signif(res$logtest[3],2))
    return(c(HR,LR_test_pvalue))
    #return(HR)
  }
  
  results<-sapply(split_list_yes,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
  HRs<-results[1,]
  #LR_tests<-results[-1,]
  #HRs<-sapply(split_list_yes,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
  for_best_split<-ifelse(HRs<1,1/HRs,HRs)
  best_split<-names(for_best_split)[which.max(for_best_split)]
  
  
  if (split=="best"){split<-best_split}
  
  
  
 if (does_suit_threshold[names(does_suit_threshold) %in% split] %in% "no"){

    
    plot(c(),xlim=c(1,10),ylim=c(1,10),axes=0,ylab="",xlab="")
    leg<-legend("left",
           text.col = "red",bty = "n",text.font=4,cex=1.5,
           legend=c("There are not enough observations
in one out of two groups. 
Try another threshold or/and split."),plot=FALSE)
    leftx <- leg$rect$left
    rightx <- (leg$rect$left + leg$rect$w) * 1.7
    topy <- leg$rect$top
    bottomy <- (leg$rect$top - leg$rect$h)
    legend(x = c(leftx, rightx), y = c(topy, bottomy),
           text.col = "red",bty = "n",text.font=4,cex=1.5,
           legend=c("There are not enough observations
in one out of two groups. 
Try another split."))
#return(NULL)
  }
  
  else{ 
    meth_status_l<-split_list[names(split_list) %in% split]
    meth_status<-meth_status_l[[1]]
    results_for_legend<-sapply(meth_status_l,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
    HR<-results_for_legend[-2,]
    LR_test_pvalue<-results_for_legend[-1,] 
    #LR_test_pvalue<-gsub("\\b0\\b","< 2.2e-16",LR_test_pvalue)
    hyper_n<-length(meth_status[meth_status %in% "hyper"])#number of obervations of hyper meth
    hypo_n<-length(meth_status[meth_status %in% "hypo"])#number of observations of hypo meth
    
  ### PLOT
  
    surv_fun<-Surv(time = as.numeric(time_to_event), event = as.numeric(event)) # 
    if(is.null(cov_vector)){ 
      KM_formula_2<-as.formula(surv_fun ~ strata(meth_status))#strata to get 2 different curves for plot (2 grps) instead of estimates
    }
    else{
      KM_formula_2<-as.formula(paste("surv_fun ~ ", paste(c(cov_vector,"strata(meth_status)"), collapse= "+")))
    }
    
    
    KM_for_plot<-coxph(KM_formula_2)
    # summary(survfit(KM_for_plot))
    #surv fit is built in command for producing survival plots
    plot(survfit(KM_for_plot),conf.int=FALSE,mark.time=FALSE,col=c("blue","red"),xlab="Survival time (days)",
         ylab="Survival Probability",main=paste0(gene," - ",region,"-",group,"-",probe)) 
    
    legend("bottomright",inset=c(0.01,0.01),cex=1.0,#lty=5:1,
           box.lwd = 0,box.col = "white",bg = "white",lwd=c(2,2),col=c("blue","red"),
           legend=c(paste0("Lower (n=",hypo_n,")"),paste0("Higher (n=",hyper_n,")")) )
    legend("topright",inset=c(0.01,0.01),cex=1.0,#lty=5:1,                         ### HR<-result[meth_ind,2]  ## (after getting the proper location put it into    #the window of plot)
           box.lwd = 0,box.col = "white",bg = "white",
           legend=c(paste0("LR test p-value=",LR_test_pvalue,"\n","HR=",HR))) 
    
  }
  
}

### density plot

distr_plot<-function(cancer_data1,cancer_data2,gene,cov_vector=NULL,threshold=0.25,region,group,probe,split){
  #head(rownames(cancer_data1)) 
  #inds<-grep(paste0("^",cpg,"$"),cancer_data1$Name,ignore.case=TRUE)[1] #,value=TRUE
  #current_gene_data<-cancer_data1[inds,-(1:4)]
  #region<-data.frame(cancer_data1[inds,(2:4)])
  #region[] <- lapply(region, as.character)
  #region_text<-paste0(region[1]," - ",region[2]," - ",region[3],".")
  
  current_gene_data=cancer_data1[which( cancer_data1$UCSC_RefGene_Name %in% gene & cancer_data1$UCSC_RefGene_Group %in% region & cancer_data1$Relation_to_UCSC_CpG_Island %in% group & cancer_data1$Name %in% probe) , ]
  # #current_gene_data<-cancer_data1[inds,]
  #rownames(current_gene_data)=current_gene_data$Name
  # #head(current_gene_data)
  current_gene_data=current_gene_data[,!colnames(current_gene_data) %in% c("Name", "CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"),with=FALSE]
  methylation<-as.numeric(current_gene_data)
  time_to_event<-as.numeric(cancer_data2$time_to_event)#time up to death from data 
  event<-as.numeric(cancer_data2$event)
  age<-as.numeric(cancer_data2$age)
  BMI<-as.numeric(cancer_data2$BMI)
sex=cancer_data2$sex
grade=cancer_data2$grade
stage=cancer_data2$clinical_stage
her2_status_by_ihc=cancer_data2$her2_status_by_ihc
er_status_by_ihc=cancer_data2$er_status_by_ihc
pr_status_by_ihc=cancer_data2$pr_status_by_ihc
pathologic_N=cancer_data2$pathologic_N
  
  # optimal splitpoint by log-rank test
  
 
  covariates<-c("methylation",cov_vector)
  #surv_fun<-Surv(time = time_to_event, event = event)
  #cph <- coxph(as.formula(paste("surv_fun ~ ", paste(covariates, collapse= "+"))))
  #op <- cutp(cph)$methylation
  #maxstat_split<-op$methylation[which.min(op$p)]
  #surv_cutpoint
  
  #current_gene=subset(cancer_data1,Name==cpg)
  t_curr=data.frame(t(current_gene_data))   
 colnames(t_curr)[1]="meth"
 t_curr$event=cancer_data2$event
    t_curr$time_to_event=cancer_data2$time_to_event
    t_curr=data.frame(t_curr)
    
    res.cut <- surv_cutpoint(t_curr, time = "time_to_event", event = "event",
                         variables = c("meth" ))
                         
 maxstat_split<- res.cut$cutpoint[,1]
  
  
  mean<-factor(ifelse(methylation<=mean(methylation),"hypo","hyper")) 
  mean<-relevel(mean, ref = "hypo")
  median<-factor(ifelse(methylation<=median(methylation),"hypo","hyper")) 
  median<-relevel(median, ref = "hypo")
  q25<-factor(ifelse(methylation<=quantile(methylation,0.25),"hypo","hyper")) 
  q25<-relevel(q25, ref = "hypo")
  q75<-factor(ifelse(methylation<=quantile(methylation,0.75),"hypo","hyper")) 
  q75<-relevel(q75, ref = "hypo")
  maxstat<-factor(ifelse(methylation<=maxstat_split,"hypo","hyper"))
  maxstat<-relevel(maxstat, ref = "hypo")
  
  #split_df<-data.frame(mean,median,q25,q75)
  split_list<-list(mean,median,q25,q75,maxstat)
  names(split_list)<-c("mean","median","q25","q75","maxstat")
  
  
  
  split_threshold<-function(x, threshold){
    t1<-table(x)#how many hyper and hypo in 1 type of variable
    does_suit<-"yes"
    if(min(t1)/sum(t1)<threshold){does_suit<-"no"}
    return(does_suit)
  }
  
  #does_suit_threshold<-apply(split_df,2,split_threshold,threshold=threshold)
  does_suit_threshold<-sapply(split_list,split_threshold,threshold=threshold)
  split_list_yes<-split_list[which(does_suit_threshold %in% "yes")]
  
  
  #meth_v<-split_list_yes[[2]]
  
  HR_fun0<-function(meth_v,te,e,covs){ 
    surv_fun<-Surv(time = as.numeric(te), event = as.numeric(e)) # survival function part of cox model
    cov_vector1<-c("meth_v",covs)#add the covariates
    if(is.null(covs)){ 
      KM_formula_1<-as.formula(surv_fun ~ meth_v)#if 0 then regress only on methylation
    }
    else{
      KM_formula_1<-as.formula(paste("surv_fun ~ ", paste(cov_vector1, collapse= "+")))#other collapse with other covariates
    }
    
    KM_for_HR<-coxph(KM_formula_1)#formula with the conditions, execute the cox model
    result<-summary(KM_for_HR)$coefficients
    res<-summary(KM_for_HR)
    meth_ind<-which(rownames(result) %in% "meth_vhyper")#estimate for methylation releveling hyper 
    HR<-round(result[meth_ind,2],3)#second coloumn  from res1
    LR_test_pvalue<-as.numeric(round(res$logtest[3],3))
    return(c(HR,LR_test_pvalue))
    #return(HR)
  }
  
  results<-sapply(split_list_yes,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
  HRs<-results[-2,]
  #LR_tests<-results[-1,]
  #HRs<-sapply(split_list_yes,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
  for_best_split<-ifelse(HRs<1,1/HRs,HRs)
  best_split<-names(for_best_split)[which.max(for_best_split)]
  
  if (split=="best"){split<-best_split}
  
  
  
  if (does_suit_threshold[names(does_suit_threshold) %in% split] %in% "no"){
#     print(paste0("There are not enough observations in one out of 
#           two groups. Try another threshold or/and split options."))
    return(NULL)
    
  }
  
  else{  
    
    
    #split_points
    meth<-as.numeric(current_gene_data)
    split_points<-round(data.frame("mean"=mean(meth,na.rm=TRUE),"median"=median(meth,na.rm=TRUE),
                             "q25"=quantile(meth,0.25),"q75"=quantile(meth,0.75),"maxstat"=maxstat_split),3)
    
    
    meth_status_l<-split_list[names(split_list) %in% split]
    meth_status<-meth_status_l[[1]]
    
    if(length(age)==0){age=rep(0,length(meth_status))}
    if(length(BMI)==0){BMI=rep(0,length(meth_status))}
    sp<-split_points[colnames(split_points) %in% split]
    q1<-split_points[colnames(split_points) %in% "q25"]
    q2<-split_points[colnames(split_points) %in% "median"]
    q3<-split_points[colnames(split_points) %in% "q75"]
    m<-split_points[colnames(split_points) %in% "mean"]
    opt<-split_points[colnames(split_points) %in% "maxstat"]
    
    #n_pars<-length(cov_vector)+1
    vars_df<-data.frame(meth_status,"meth"=meth,age,BMI)
    
    #if (n_pars>2){layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))}
    #if (n_pars==2){par(mfrow=c(2,1))}
    
    
    #vars_df_current<-data.frame(vars_df[,colnames(vars_df) %in% cov_vector])
    #names(vars_df_current)<-cov_vector
    
    #max(density(vars_df$meth)$y)/2
    #par(cex.axis=1.0)
    par(font.axis=1,font.lab=1,cex.lab=1.3,cex.axis=1.2)
    plot(density(vars_df$meth),xlab="Beta",main=paste0("Splitting of methylation by ",split))
    abline(v=q1,col="#FF1493",lty=3) 
    
    text(q1+(0.07*(max(vars_df$meth,na.rm=TRUE)-min(vars_df$meth,na.rm=TRUE))), 0+0.8*(max(density(vars_df$meth)$y,na.rm=TRUE)), "q25",cex = .8,col="#FF1493")
    #median
    abline(v=q2,col="blue",lty=3) 
    text(q2+(0.07*(max(vars_df$meth,na.rm=TRUE)-min(vars_df$meth,na.rm=TRUE))), 0+0.8*(max(density(vars_df$meth)$y,na.rm=TRUE)), "median",cex = .8,col="blue")
    #q75
    abline(v=q3,col="orange",lty=3) 
    text(q3+(0.07*(max(vars_df$meth,na.rm=TRUE)-min(vars_df$meth,na.rm=TRUE))), 0+0.8*(max(density(vars_df$meth)$y,na.rm=TRUE)), "q75",cex = .8,col="orange")
    #Mean
    abline(v=m,col="#006400",lty=3) 
    text(m+(0.07*(max(vars_df$meth,na.rm=TRUE)-min(vars_df$meth,na.rm=TRUE))), 0+0.9*(max(density(vars_df$meth)$y,na.rm=TRUE)), "mean",cex = .8,col="#006400")
    #opt_p
    abline(v=opt,col="#700000",lty=3) 
    text(opt+(0.07*(max(vars_df$meth,na.rm=TRUE)-min(vars_df$meth,na.rm=TRUE))), 0+0.7*(max(density(vars_df$meth)$y,na.rm=TRUE)), "maxstat",cex = .8,col="#700000")
    #current split
    abline(v=sp,col="red") 
    text(sp+(0.07*(max(vars_df$meth,na.rm=TRUE)-min(vars_df$meth,na.rm=TRUE))), 0+0.07*(max(density(vars_df$meth)$y,na.rm=TRUE)), paste0(sp),cex = .8,col="red")
    
    
    
    
     
}
}




### distribution plot end


###table start

cox_model_best_split_table<-function(cancer_data1,cancer_data2,gene,cov_vector=NULL,threshold=0.25,region,group,probe,split){
  #head(rownames(cancer_data1)) 
  
  #inds<-grep(paste0("^",cpg,"$"),cancer_data1$Name,ignore.case=TRUE)[1] #,value=TRUE
  #current_gene_data<-cancer_data1[inds,-(1:4)]
  #region<-data.frame(cancer_data1[inds,(2:4)])
  #region[] <- lapply(region, as.character)
  #region_text<-paste0(region[1]," - ",region[2]," - ",region[3],".")
  
  current_gene_data=cancer_data1[which( cancer_data1$UCSC_RefGene_Name %in% gene & cancer_data1$UCSC_RefGene_Group %in% region & cancer_data1$Relation_to_UCSC_CpG_Island %in% group & cancer_data1$Name %in% probe) , ]
  # #current_gene_data<-cancer_data1[inds,]
  #rownames(current_gene_data)=current_gene_data$Name
  # #head(current_gene_data)
   current_gene_data=current_gene_data[,!colnames(current_gene_data) %in% c("Name", "CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"),with=FALSE]
   methylation<-as.numeric(current_gene_data)
  time_to_event<-as.numeric(cancer_data2$time_to_event)#time up to death from data 
  event<-as.numeric(cancer_data2$event)
  age<-as.numeric(cancer_data2$age)
  BMI<-as.numeric(cancer_data2$BMI)
sex=cancer_data2$sex
grade=cancer_data2$grade
stage=cancer_data2$clinical_stage
her2_status_by_ihc=cancer_data2$her2_status_by_ihc
er_status_by_ihc=cancer_data2$er_status_by_ihc
pr_status_by_ihc=cancer_data2$pr_status_by_ihc
pathologic_N=cancer_data2$pathologic_N
  
  # optimal splitpoint by log-rank test
  
 # optimal splitpoint by log-rank test
  
  covariates<-c("methylation",cov_vector)
  surv_fun<-Surv(time = time_to_event, event = event)
  #cph <- coxph(as.formula(paste("surv_fun ~ ", paste(covariates, collapse= "+"))))
  #op <- cutp(cph)$methylation
  #maxstat_split<-op$methylation[which.min(op$p)]
  #surv_fun<-Surv(time = time_to_event, event = event)
  #cph <- coxph(as.formula(paste("surv_fun ~ ", paste(covariates, collapse= "+"))))
  #op <- cutp(cph)$methylation
  #maxstat_split<-op$methylation[which.min(op$p)]
  #surv_cutpoint
  
  #current_gene=subset(cancer_data1,Name==cpg)
  t_curr=data.frame(t(current_gene_data))   
 colnames(t_curr)[1]="meth"
 t_curr$event=cancer_data2$event
    t_curr$time_to_event=cancer_data2$time_to_event
    t_curr=data.frame(t_curr)
    
    res.cut <- surv_cutpoint(t_curr, time = "time_to_event", event = "event",
                         variables = c("meth" ))
                         
 maxstat_split<- res.cut$cutpoint[,1]
  
  
  mean<-factor(ifelse(methylation<=mean(methylation),"hypo","hyper")) 
  mean<-relevel(mean, ref = "hypo")
  median<-factor(ifelse(methylation<=median(methylation),"hypo","hyper")) 
  median<-relevel(median, ref = "hypo")
  q25<-factor(ifelse(methylation<=quantile(methylation,0.25),"hypo","hyper")) 
  q25<-relevel(q25, ref = "hypo")
  q75<-factor(ifelse(methylation<=quantile(methylation,0.75),"hypo","hyper")) 
  q75<-relevel(q75, ref = "hypo")
  maxstat<-factor(ifelse(methylation<=maxstat_split,"hypo","hyper"))
  maxstat<-relevel(maxstat, ref = "hypo")
  
  #split_df<-data.frame(mean,median,q25,q75)
  split_list<-list(mean,median,q25,q75,maxstat)
  names(split_list)<-c("mean","median","q25","q75","maxstat")
  
  
  
  split_threshold<-function(x, threshold){
    t1<-table(x)#how many hyper and hypo in 1 type of variable
    does_suit<-"yes"
    if(min(t1)/sum(t1)<threshold){does_suit<-"no"}
    return(does_suit)
  }
  
  #does_suit_threshold<-apply(split_df,2,split_threshold,threshold=threshold)
  does_suit_threshold<-sapply(split_list,split_threshold,threshold=threshold)
  split_list_yes<-split_list[which(does_suit_threshold %in% "yes")]
  
  
  #meth_v<-split_list_yes$mean
  
  HR_fun0<-function(meth_v,te,e,covs){ 
    
    surv_fun<-Surv(time = as.numeric(te), event = as.numeric(e)) # survival function part of cox model
    cov_vector1<-c("meth_v",covs)#add the covariates
    if(is.null(covs)){ 
      KM_formula_1<-as.formula(surv_fun ~ meth_v)#if 0 then regress only on methylation
    }
    else{
      KM_formula_1<-as.formula(paste("surv_fun ~ ", paste(cov_vector1, collapse= "+")))#other collapse with other covariates
    }
    
    KM_for_HR<-coxph(KM_formula_1)#formula with the conditions, execute the cox model
    phtest<-cox.zph(KM_for_HR)$table 
    phtest_p_local<-round(phtest[!rownames(phtest) %in% "GLOBAL",3],3)
    #phtest_p_global<-round(phtest[dim(phtest)[1],3],3)
    res<-summary(KM_for_HR)
    result<-res$coefficients
    #meth_ind<-which(rownames(result) %in% "meth_vhyper")
    #HR_meth<-round(result[meth_ind,2],3)
    HR<-round(result[,2],3)
    #HR_p<-round(result[,5],3)
    HR_p<-signif(result[,5],2)
    CI_df<-round(data.frame(res$conf.int),3)
    CI<-paste0("(",CI_df[,3],";",CI_df[,4],")")
    names(CI)<-rownames(CI_df)
    #LR_test_pvalue<-as.numeric(round(res$logtest[3],3))
    return(c(HR,CI,HR_p,phtest_p_local))#res,phtest))
    #return(HR)
  }
  
  #cov_vector<-c("age")
  #cov_vector<-NULL
  results<-sapply(split_list_yes,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
  HRs<-results[1,]
  for_best_split<-ifelse(as.numeric(HRs)<1,1/as.numeric(HRs),as.numeric(HRs))
  best_split<-colnames(results)[which.max(for_best_split)]
#   if(length(HRs)>1){best_split<-names(HRs)[which.max(for_best_split)]}
#   else{best_split<-names(split_list_yes)}
  
  if (split=="best"){split<-best_split}
  
 
  
  if (does_suit_threshold[names(does_suit_threshold) %in% split] %in% "no"){
#     print(paste0("There are not enough observations in one out of 
#           two groups according to threshold=",threshold," and split=",split,". 
#                  Try another threshold or/and split"))
    #plot(c(),xlim=c(1,10),ylim=c(1,10),axes=0,ylab="",xlab="")
return(NULL)
  }
  
  else{  
    
    
    
    meth_status_l<-split_list[names(split_list) %in% split]
    meth_status<-meth_status_l[[1]]
    results_for_table<-sapply(meth_status_l,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
    res_names<-unique(rownames(results_for_table))
    res_names<-res_names[!res_names %in% ""]
    n<-length(res_names)
    res_tab<-data.frame(matrix(results_for_table,nrow=n))
    rownames(res_tab)<-res_names
    rownames(res_tab)[1]<-"methylation(hyper)"
    colnames(res_tab)<-c("HR","CI","Wald_Pvalue","PH_test_Pvalue")
    
    ######
    # updated - add the additional colums with mean median and range:
    mean_meth<-round(mean(as.numeric(current_gene_data)),3)
    median_meth<-round(median(as.numeric(current_gene_data)),3)
    range_meth<-paste0("(",round(range(as.numeric(current_gene_data)),3)[1],";",round(range(as.numeric(current_gene_data)),3)[2],")")
    
    mean_vector<-c()
    median_vector<-c()
    range_vector<-c()
    
    mean_vector[1]<-mean_meth
    median_vector[1]<-median_meth
    range_vector[1]<-range_meth
    
    if ("age" %in% rownames(res_tab)){
      ind_age<-which(rownames(res_tab) %in% "age")
      mean_vector[ind_age]<-floor(mean(age,na.rm=TRUE))
      median_vector[ind_age]<-floor(median(age,na.rm=TRUE))
      range_vector[ind_age]<-paste0("(",floor(range(age,na.rm=TRUE))[1],";",floor(range(age,na.rm=TRUE))[2],")")
      #mean_age<-floor(mean(age,na.rm=TRUE));median_age<-floor(median(age,na.rm=TRUE));range_age<-paste0("(",floor(range(age,na.rm=TRUE))[1],";",floor(range(age,na.rm=TRUE))[2],")")}
    }
    
    if ("BMI" %in% rownames(res_tab)){
      ind_BMI<-which(rownames(res_tab) %in% "BMI")
      mean_vector[ind_BMI]<-floor(mean(BMI,na.rm=TRUE))
      median_vector[ind_BMI]<-floor(median(BMI,na.rm=TRUE))
      range_vector[ind_BMI]<-paste0("(",floor(range(BMI,na.rm=TRUE))[1],";",floor(range(BMI,na.rm=TRUE))[2],")")
    }   
    
   ######
    
     res_tab<-data.frame(res_tab,"Current_split"=split,"Best_split"=best_split,
                       "Mean"=mean_vector,"Median"=median_vector,"Range"=range_vector)
  res_tab=cbind(Name = rownames(res_tab), res_tab)
     res_tab=subset(res_tab,row.names(res_tab)=="methylation(hyper)")
   return(res_tab)
  }
  #return(res_tab)
}




##table function ends

###table end







###heatmap
drawhmap <- function(mat, gene,annot){
  
  mat=mat[which(mat$UCSC_RefGene_Name %in% gene) , ]
  mat=as.data.frame(mat)
  #annotation_col=subset(annot,select=-c(bcr_patient_barcode ,PID ,time_to_event) )
  #annotation_col=subset(annotation_data,select=c("race","age") )
  annotation_row=subset(mat,select=c("UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"))
  #mat=as.matrix(mat)
  rownames(mat)=mat$Name
  mat=mat[,!colnames(mat) %in% c("Name", "CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island")]
  #rownames(annotation_col) = colnames(mat)
  rownames(annotation_row) = rownames(mat)
  rowNamesSize = min(16, max(1, floor(650 / nrow(mat))))
  colNamesSize = min(16, max(1, floor(350 / ncol(mat))))
  pheatmap(mat, annotation_col = annot, annotation_row = annotation_row,fontsize_row=rowNamesSize,fontsize_col=colNamesSize)
  
} 





#### new lines###

load("ACC_vio.RData")   
load("BRCA_vio.RData")  
load("COAD_vio.RData")  
load("GBM_vio.RData")   
load("KICH_vio.RData")  
load("KIRP_vio.RData")  
load("LGG_vio.RData")   
load("LUAD_vio.RData")  
load("MESO_vio.RData")  
load("READ_vio.RData")  
load("SKCM_vio.RData")  
load("UCEC_vio.RData")  
load("UVM_vio.RData")
load("BLCA_vio.RData")  
load("CESC_vio.RData")  
load("ESCA_vio.RData")  
load("HNSC_vio.RData")  
load("KIRC_vio.RData")  
load("LAML_vio.RData")  
load("LIHC_vio.RData")  
load("LUSC_vio.RData")  
load("PAAD_vio.RData")  
load("SARC_vio.RData")  
load("STAD_vio.RData")  
load("UCS_vio.RData")

 drawvp <- function(meth_mat, var, gene, cpg, annot){

 meth_mat=meth_mat[which( meth_mat$UCSC_RefGene_Name %in% gene & meth_mat$Name %in% cpg) , ]
  #meth_mat=subset(meth_mat,UCSC_RefGene_Name==gene & Name==cpg)
  meth_mat=meth_mat[,!colnames(meth_mat)%in% c("Name", "CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"),with=FALSE]
  meth_mat=t(meth_mat)
  colnames(meth_mat)="Beta"
  annot=merge(annot,meth_mat,by="row.names")
  var_tab=subset(annot,select=var)
  xlabs=paste("N=",table(var_tab),sep="")
  ggplot(na.omit(annot), aes_string(x = var,y="Beta",fill=var)) +geom_violin(alpha=0.25)+geom_boxplot(width=0.1, fill="white",outlier.shape=NA)+ theme_bw()+scale_x_discrete(labels=xlabs)+
  theme(
 legend.title=element_text(size=10) , 
 legend.text=element_text(size=10),
 axis.text = element_text(size=10),axis.text.x=element_text(size=10),axis.text.y=element_text(size=10)

 )
}

###for buttons

shinyInput <- function(FUN, len, id, ...) {inputs <- character(len)
for (i in seq_len(len)) {
  inputs[i] <- as.character(FUN(paste0(id, i), ...))}
inputs
}

###for table based cpgs

  
 ###survival plot for region based analysis
cox_model_best_split_plot_v2<-function(cancer_data1,cancer_data2,gene,cov_vector=NULL,threshold=0.25,probe,split){
  #head(rownames(cancer_data1)) 
  
  #inds<-grep(paste0("^",cpg,"$"),cancer_data1$Name,ignore.case=TRUE)[1] #,value=TRUE
  #current_gene_data<-cancer_data1[inds,-(1:4)]
  #region<-data.frame(cancer_data1[inds,(2:4)])
  #region[] <- lapply(region, as.character)
  #region_text<-paste0(region[1]," - ",region[2]," - ",region[3],".")
  
  current_gene_data=cancer_data1[which(cancer_data1$Name %in% probe) , ]
gene=as.character(cancer_data1$UCSC_RefGene_Name)
region=as.character(cancer_data1$UCSC_RefGene_Group)
group=as.character(cancer_data1$Relation_to_UCSC_CpG_Island)

  # #current_gene_data<-cancer_data1[inds,]
  #rownames(current_gene_data)=current_gene_data$Name
  # #head(current_gene_data)
  current_gene_data=current_gene_data[,!colnames(current_gene_data) %in% c("Name", "CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_UCSC_CpG_Island"),with=FALSE]
  methylation<-as.numeric(current_gene_data)
  time_to_event<-as.numeric(cancer_data2$time_to_event)#time up to death from data 
  event<-as.numeric(cancer_data2$event)
  age<-as.numeric(cancer_data2$age)
  BMI<-as.numeric(cancer_data2$BMI)
sex=cancer_data2$sex
grade=cancer_data2$grade
stage=cancer_data2$clinical_stage
her2_status_by_ihc=cancer_data2$her2_status_by_ihc
er_status_by_ihc=cancer_data2$er_status_by_ihc
pr_status_by_ihc=cancer_data2$pr_status_by_ihc
pathologic_N=cancer_data2$pathologic_N
  
  # optimal splitpoint by log-rank test
  
  covariates<-c("methylation",cov_vector)
  #surv_fun<-Surv(time = time_to_event, event = event)
  #cph <- coxph(as.formula(paste("surv_fun ~ ", paste(covariates, collapse= "+"))))
  #op <- cutp(cph)$methylation
  #maxstat_split<-op$methylation[which.min(op$p)]
  #surv_cutpoint
  
  #current_gene=subset(cancer_data1,Name==cpg)
  t_curr=data.frame(t(current_gene_data))   
 colnames(t_curr)[1]="meth"
 t_curr$event=cancer_data2$event
    t_curr$time_to_event=cancer_data2$time_to_event
    t_curr=data.frame(t_curr)
    
    res.cut <- surv_cutpoint(t_curr, time = "time_to_event", event = "event",
                         variables = c("meth" ))
                         
 maxstat_split<- res.cut$cutpoint[,1]
                         
  
  ###
  
  mean<-factor(ifelse(methylation<=mean(methylation),"hypo","hyper")) 
  mean<-relevel(mean, ref = "hypo")
  median<-factor(ifelse(methylation<=median(methylation),"hypo","hyper")) 
  median<-relevel(median, ref = "hypo")
  q25<-factor(ifelse(methylation<=quantile(methylation,0.25),"hypo","hyper")) 
  q25<-relevel(q25, ref = "hypo")
  q75<-factor(ifelse(methylation<=quantile(methylation,0.75),"hypo","hyper")) 
  q75<-relevel(q75, ref = "hypo")
  maxstat<-factor(ifelse(methylation<=maxstat_split,"hypo","hyper"))
  maxstat<-relevel(maxstat, ref = "hypo")
  
  #split_df<-data.frame(mean,median,q25,q75)
  split_list<-list(mean,median,q25,q75)
  names(split_list)<-c("mean","median","q25","q75")
    
    
    
  split_threshold<-function(x, threshold){
    t1<-table(x)#how many hyper and hypo in 1 type of variable
    does_suit<-"yes"
    if(min(t1)/sum(t1)<threshold){does_suit<-"no"}
    return(does_suit)
  }
  
  #does_suit_threshold<-apply(split_df,2,split_threshold,threshold=threshold)
  does_suit_threshold<-sapply(split_list,split_threshold,threshold=threshold)
  split_list_yes<-split_list[which(does_suit_threshold %in% "yes")]
  
  #meth_v<-split_list_yes[[1]]
  
  HR_fun0<-function(meth_v,te,e,covs){ 
    surv_fun<-Surv(time = as.numeric(te), event = as.numeric(e)) # survival function part of cox model
    cov_vector1<-c("meth_v",covs)#add the covariates
    if(is.null(covs)){ 
      KM_formula_1<-as.formula(surv_fun ~ meth_v)#if 0 then regress only on methylation
    }
    else{
      KM_formula_1<-as.formula(paste("surv_fun ~ ", paste(cov_vector1, collapse= "+")))#other collapse with other covariates
    }
    
    KM_for_HR<-coxph(KM_formula_1)#formula with the conditions, execute the cox model
    result<-summary(KM_for_HR)$coefficients
    res<-summary(KM_for_HR)
    meth_ind<-which(rownames(result) %in% "meth_vhyper")#estimate for methylation releveling hyper 
    HR<-round(result[meth_ind,2],3)#second coloumn  from res1
    LR_test_pvalue<-as.numeric(signif(res$logtest[3],2))
    return(c(HR,LR_test_pvalue))
    #return(HR)
  }
  
  results<-sapply(split_list_yes,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
  HRs<-results[1,]
  #LR_tests<-results[-1,]
  #HRs<-sapply(split_list_yes,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
  for_best_split<-ifelse(HRs<1,1/HRs,HRs)
  best_split<-names(for_best_split)[which.max(for_best_split)]
  
  
  if (split=="best"){split<-best_split}
  
  
  
 if (does_suit_threshold[names(does_suit_threshold) %in% split] %in% "no"){

    
    plot(c(),xlim=c(1,10),ylim=c(1,10),axes=0,ylab="",xlab="")
    leg<-legend("left",
           text.col = "red",bty = "n",text.font=4,cex=1.5,
           legend=c("There are not enough observations
in one out of two groups. 
Try another threshold or/and split."),plot=FALSE)
    leftx <- leg$rect$left
    rightx <- (leg$rect$left + leg$rect$w) * 1.7
    topy <- leg$rect$top
    bottomy <- (leg$rect$top - leg$rect$h)
    legend(x = c(leftx, rightx), y = c(topy, bottomy),
           text.col = "red",bty = "n",text.font=4,cex=1.5,
           legend=c("There are not enough observations
in one out of two groups. 
Try another split."))
#return(NULL)
  }
  
  else{ 
    meth_status_l<-split_list[names(split_list) %in% split]
    meth_status<-meth_status_l[[1]]
    results_for_legend<-sapply(meth_status_l,HR_fun0,te=time_to_event,e=event,covs=cov_vector)
    HR<-results_for_legend[-2,]
    LR_test_pvalue<-results_for_legend[-1,] 
    #LR_test_pvalue<-gsub("\\b0\\b","< 2.2e-16",LR_test_pvalue)
    hyper_n<-length(meth_status[meth_status %in% "hyper"])#number of obervations of hyper meth
    hypo_n<-length(meth_status[meth_status %in% "hypo"])#number of observations of hypo meth
    
  ### PLOT
  
    surv_fun<-Surv(time = as.numeric(time_to_event), event = as.numeric(event)) # 
    if(is.null(cov_vector)){ 
      KM_formula_2<-as.formula(surv_fun ~ strata(meth_status))#strata to get 2 different curves for plot (2 grps) instead of estimates
    }
    else{
      KM_formula_2<-as.formula(paste("surv_fun ~ ", paste(c(cov_vector,"strata(meth_status)"), collapse= "+")))
    }
    
    
     KM_for_plot<-coxph(KM_formula_2)
    # summary(survfit(KM_for_plot))
    #surv fit is built in command for producing survival plots
    plot(survfit(KM_for_plot),conf.int=FALSE,mark.time=FALSE,col=c("blue","red"),xlab="Survival time (days)",
         ylab="Survival Probability",main=paste0(probe)) 
    
    legend("bottomright",inset=c(0.01,0.01),cex=1.0,#lty=5:1,
           box.lwd = 0,box.col = "white",bg = "white",lwd=c(2,2),col=c("blue","red"),
           legend=c(paste0("Lower (n=",hypo_n,")"),paste0("Higher (n=",hyper_n,")")) )
    legend("topright",inset=c(0.01,0.01),cex=1.0,#lty=5:1,                         ### HR<-result[meth_ind,2]  ## (after getting the proper location put it into    #the window of plot)
           box.lwd = 0,box.col = "white",bg = "white",
           legend=c(paste0("LR test p-value=",LR_test_pvalue,"\n","HR=",HR))) 
    
  }
  
}



####for region based analysis tab

load("ACC_null_2017.RData")   
load("BRCA_null_2017.RData")  
load("COAD_null_2017.RData")  
load("GBM_null_2017.RData")   
load("KICH_null_2017.RData")  
load("KIRP_null_2017.RData")  
load("LGG_null_2017.RData")   
load("LUAD_null_2017.RData")  
load("MESO_null_2017.RData")  
load("READ_null_2017.RData")  
load("SKCM_null_2017.RData")  
load("UCEC_null_2017.RData")  
load("UVM_null_2017.RData")
load("BLCA_null_2017.RData")  
load("CESC_null_2017.RData")  
load("ESCA_null_2017.RData")  
load("HNSC_null_2017.RData")  
load("KIRC_null_2017.RData")  
load("LAML_null_2017.RData")  
load("LIHC_null_2017.RData")  
load("LUSC_null_2017.RData")  
load("PAAD_null_2017.RData")  
load("SARC_null_2017.RData")  
load("STAD_null_2017.RData")  
load("UCS_null_2017.RData")
