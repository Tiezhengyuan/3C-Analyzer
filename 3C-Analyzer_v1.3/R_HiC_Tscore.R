#version 1.1

################
#R packages
if (!require("MASS")) {
  install.packages("MASS", dependencies = TRUE)
}
if (!require("nlmrt")) {
  install.packages("nlmrt", dependencies = TRUE)
}
if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
}
library(plyr)
library(MASS)
library(nlmrt)
###################

####################################
#options
####################################

R_ligation_frequency_dir="/home/yuan/3C_analyzer/result/ligation_frequency"
R_RC_csv_name="SizeBin_trans_ligation_GSE18199_SRR027958_GM22_RC.csv"
R_size_bin_csv="/home/yuan/3C_analyzer/result/genome_size_bin.csv"
R_enzyme_bin_csv="/home/yuan/3C_analyzer/result/genome_enzyme_bin.csv"
R_p_level="0.01"
R_RC_noise="2"
R_permutation_num="500"
R_sampling_ratio="0.4"

ligations_dir=R_ligation_frequency_dir
RC_csv_name=R_RC_csv_name
in_file=paste(ligations_dir, RC_csv_name, sep="/")
size_bin_csv=R_size_bin_csv
enzyme_bin_csv=R_enzyme_bin_csv
p.level=as.numeric(R_p_level)
RC.noise=as.numeric(R_RC_noise)
permutation.num=as.numeric(R_permutation_num)
sampling.ratio=as.numeric(R_sampling_ratio)


##################
#subroutines
##################
######
cis_regression<-function(RC_df){
  nls_class<-list()
  max(RC_df$RC)

  #get maximum RC grouped by chromosomal distance
  RC_group<- ddply(RC_df, ~dist, summarise, max=max(RC), mean=mean(RC), median=median(RC), 
                   min=min(RC), sum=sum(RC), sd=sd(RC), num=length(RC))
  head(RC_group)
  dim(RC_group)
  #scaling sum of RCs
  RC.fold<-round(log(max(RC_group[,'max'])/1e4, base=10))
  if(RC.fold<0){ RC.fold=0 }
  RC_group[,'rRC']<-RC_group[,'max']/(10**RC.fold) 
  
  #estimate parameters of double exponential functions
  cis_lm<-lm(rRC~dist, data=RC_group)
  cis_nls<-nlxb(rRC~RC0*exp(-k0*dist)+RC1*exp(-k1*dist), data=RC_group, 
                start=list(RC0=abs(cis_lm$coeff[1]), k0=abs(cis_lm$coeff[2]), RC1=abs(cis_lm$coeff[1]), 
                           k1=abs(cis_lm$coeff[2])), control=list(maxiter=1000) )
  cis_nls
  RC0<-abs(cis_nls$coeff[1])
  k0<-abs(cis_nls$coeff[2])
  RC1<-abs(cis_nls$coeff[3])
  k1<-abs(cis_nls$coeff[4])
  
  #calculate distRC
  mRC<-round( (RC0*exp(-k0*RC_df[,'dist'])+RC1*exp(-k1*RC_df[,'dist']))*(10**RC.fold) )
    head(sort(mRC,decreasing=T))
  distRC<-round(RC_df[,'RC']-mRC)
  distRC[distRC<0]<-0
    head(sort(distRC,decreasing=T))
  
  #
  nls_class[['double_exponential']]<-cis_nls
  nls_class[['mRC']]<-mRC
  nls_class[['distRC']]<-distRC
  return(nls_class)
}

##########
sig_analysis<-function(RC, RC.noise, permutation.num, p.level){
  
  #initiate
  exp_fit_rate=0
  p.values<-rep(0, length=length(RC))
  names(p.values)<-names(RC)
  p.adj<-p.values
  Tscore<-p.values
  sigRC.limit<-0
  
  #estimate lamda by permutation
  limits_num<-round(sampling.ratio*length(RC))
  samples_num<-ifelse(limits_num<99999, limits_num, 99999) 
  exp_fit_rates<-rep(0,length=permutation.num)
  for(i in 1:permutation.num){
    random_RC<-sample(RC, samples_num, replace=F)
    #plot(density(random_RC), xlim=c(0,10))
    #qqnorm(random_RC)
    #fit parameter of exponential function for significance analysis
    exp_fit_rates[i]<-fitdistr(random_RC, densfun="exponential")$estimate
  }
  #lamda of exponential distribution density funciton
  exp_fit_rate<-mean(exp_fit_rates)
  hist(exp_fit_rates, main=exp_fit_rate)
  abline(v=exp_fit_rate, col='red')
  
  #calculate cumulative probability:  probabilities are P[X > x]
  p.values<-pexp(RC, rate=exp_fit_rate, lower.tail=F)
  p.values[p.values==0|p.values<1e-320]<-1e-320
  #BH-adjusted p values
  p.adj<-p.adjust(p.values, "BH")
  #convert p.values
  Tscore<--10*log10(p.values)
  #RC limit based on exp_fit_df
  sigRC.limit<-qexp(p.level, rate=exp_fit_rate, lower.tail=F)

  
  sig.class<-list()
  sig.class[['exp_fit_rates']]<-exp_fit_rates
  sig.class[['exp.rate']]<-exp_fit_rate
  sig.class[['p.values']]<-p.values
  sig.class[['p.adj']]<-p.adj
  sig.class[['Tscore']]<-Tscore
  sig.class[['sigRC.limit']]<-sigRC.limit
  return(sig.class)
}




###############################
#main programming
################################
#options for Tscore cal
bin_type<-ifelse(grepl('Size', RC_csv_name), 'size_bin_name', 'enzyme_bin_name')
interactions_type<-ifelse(grepl('cis_ligation', RC_csv_name), 'cis', 'trans')

#Tscore calculation
print(paste('Calculate Tscores of ', interactions_type, ' interactions of ', RC_csv_name, sep="") )
#read table into RC_df!!!!!!!!!!!!!!!!!!1
RC_df<-read.table(file=in_file, header=T, sep=',')
head(RC_df)

###############
if(interactions_type== 'cis'){
  #
  RC_df[,'dist']<-abs(RC_df[,'bin1']-RC_df[,'bin2'])
  
  #cis non-linear regression
  nls_class<-cis_regression(RC_df)
  RC_df[,'mRC']<-nls_class[['mRC']]
  RC_df[,'distRC']<-nls_class[['distRC']]
  
  #exponential signficance analysis
  #fit parameter of exponential function for significance analysis of cis-interactions
  sig_class<-sig_analysis(RC=RC_df[,'distRC'], RC.noise, permutation.num, p.level)


}else {
  #exponential signficance analysis
  #fit parameter of exponential function for significance analysis of cis-interactions
  sig_class<-sig_analysis(RC=RC_df[,'RC'], RC.noise, permutation.num, p.level)
  
  
}
RC_df[,'p.values']<-sig_class[['p.values']]  
RC_df[,'p.adj']<-sig_class[['p.adj']]  
RC_df[,'Tscore']<-sig_class[['Tscore']]  

#export R data
print("Save R object Tscore_class to the file")
RData_file_name<-sub(pattern="_RC.csv", replacement="_Tscore.RData", RC_csv_name)
RData_file<-paste(ligations_dir, RData_file_name, sep='/')
save(RC_df, file = RData_file )

#export data frome into csv file
print("Save the data frame Tscore into the file")
Tscore_csv_name<-sub(pattern="_RC.csv", replacement="_Tscore.csv", RC_csv_name)
Tscore_csv<-paste(ligations_dir, Tscore_csv_name, sep='/')
write.csv(RC_df, file=Tscore_csv, row.names=F, quote=F)


####################################
#end
####################################
