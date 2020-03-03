#version 1.1

################
#R packages
library(plyr)
library(MASS)
library(nlmrt)
###################

####################################
#options
####################################

R_ligation_frequency_dir="/home/yuan/data_1/HiC/result/ligation_frequency"
R_RC_csv_name="Size_bin_H1_hESC_trans_ligation_RC.csv"
R_size_bin_csv="/home/yuan/data_1/HiC/result/genome_size_bin.csv"
R_enzyme_bin_csv="/home/yuan/data_1/HiC/result/genome_enzyme_bin.csv"
R_p_level="0.01"
R_RC_noise="1"
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
#Tscore of cis interactions
cis_Tscore_calculation<-function(in_file, RC.noise=2, permutation.num, sampling.ratio){
  Tscore_class<-vector(mode='list')
  
  #read RC file
  RC_1<-read.table(file=in_file, header=T, sep=',')
  #max(RC_1$RC)
  #get mean of RCs group by chromosomal distance
  RC_2<-RC_1[RC_1$RC>=RC.noise,]
  RC_group<- ddply(RC_2, ~dist, summarise, max=max(RC), mean=mean(RC), median=median(RC), 
                   min=min(RC), sum=sum(RC), sd=sd(RC), num=length(RC))
  RC.fold<-round(log(max(RC_group$max)/1e3, base=10)) #scaling sum of RCs
  if(RC.fold<0){RC.fold=0}
  RC_group$rRC<-RC_group$max/(10**RC.fold) 
  head(RC_group, n=20)
  dim(RC_group)
  Tscore_class[['RC_group']]<-RC_group
  
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
  Tscore_class[['double_exponential']]<-cis_nls

  #calculate distRC
  theoRC<-round( (RC0*exp(-k0*RC_1$dist)+RC1*exp(-k1*RC_1$dist))*(10**RC.fold) )
  head(sort(theoRC,decreasing=T))
  RC_1[,'mRC']<-theoRC
  distRC<-round(RC_1$RC-theoRC)
  distRC[distRC<0]<-0
  RC_1[,'distRC']<-round(distRC)
  head(sort(distRC,decreasing=T))
  #get max of distRCs group by chromosomal distance
  RC_3<-RC_1[RC_1$distRC>=RC.noise,]
  distRC_group<- ddply(RC_3, ~dist, summarise, mean=mean(distRC), median=median(distRC),
                        max=max(distRC), sum=sum(distRC), sd=sd(distRC), num=length(distRC))
  head(distRC_group)
  dim(distRC_group)
  Tscore_class[['distRC_group']]<-distRC_group
  ######
  #fit parameter of exponential function for significance analysis of cis-interactions
  #estimate lamda by permutation
  limits_num<-round(sampling.ratio*length(RC_2$RC))
  samples_num<-ifelse(limits_num<99999, limits_num, 99999) 
  exp_fit_rates<-rep(0,length=permutation.num)
  for(i in 1:permutation.num){
    random_RC<-sample(RC_3$distRC, samples_num, replace=F)
    #plot(density(random_RC), xlim=c(0,10))
    #qqnorm(random_RC)
    #fit parameter of exponential function for significance analysis
    exp_fit_rates[i]<-fitdistr(random_RC, densfun="exponential")$estimate
  }
  Tscore_class[['exp_fit_rates']]<-exp_fit_rates
  #lamda of exponential distribution density funciton
  exp_fit_rate<-mean(exp_fit_rates)
  hist(exp_fit_rates, main=exp_fit_rate)
  abline(v=exp_fit_rate, col='red')
  
  #calculate probability:  probabilities are P[X > x]
  P<-pexp(RC_1$distRC, rate=exp_fit_rate, lower.tail=F) 
  P[P<1e-320|P==0]<-1e-320
  RC_1[,'pvalue']<-P
  #Tscore calculation
  Tscore<-round(-10*log(P,base=10), digits=1)
  head(sort(Tscore,decreasing=T))
  RC_1[,'Tscore']<-Tscore
  #BH-adjusted p values
  BH.P<-p.adjust(P, "BH")
  RC_1[,'BH.pvalue']<-BH.P
  #Tscore calculation
  BH.Tscore<-round(-10*log(BH.P,base=10), digits=1)
  head(sort(BH.Tscore,decreasing=T))
  RC_1[,'BH.Tscore']<-BH.Tscore
  
  #sig_RC<-subset(RC_1, Tscore>=20&dist<5)
  #sig_RC<-sig_RC[order(sig_RC$Tscore,decreasing=T),]
  #head(sig_RC, n=10)
  #dim(sig_RC) 
 
  Tscore_class[['Tscore_df']]<-RC_1
  return(Tscore_class)
}


#####################
#Tscore of trans-interactions
trans_Tscore_calculation<-function(in_file, RC.noise=2, permutation.num, sampling.ratio){
  Tscore_class<-vector(mode='list')
  
  #read RC file
  RC_1<-read.table(file=in_file, header=T, sep=',')
  RC_2<-RC_1[RC_1$RC>=RC.noise,]
  RC_group<- ddply(RC_2, ~RC, summarise, max=max(RC), mean=mean(RC), median=median(RC), 
                   min=min(RC), sum=sum(RC), sd=sd(RC), num=length(RC))
  head(RC_group, n=10)
  dim(RC_group)
  Tscore_class[['RC_group']]<-RC_group
  #plot(density(RC_2$RC), xlim=c(0,100))
  #qqnorm(RC_2$RC)
  
  #estimate lamda by permutation
  limits_num<-round(sampling.ratio*length(RC_2$RC))
  samples_num<-ifelse(limits_num<99999, limits_num, 99999) 
  exp_fit_rates<-rep(0,length=permutation.num)
  for(i in 1:permutation.num){
    random_RC<-sample(RC_2$RC, samples_num, replace=F)
    #plot(density(random_RC), xlim=c(0,50))
    #qqnorm(random_RC)
    #fit parameter of exponential function for significance analysis
    exp_fit_rates[i]<-fitdistr(random_RC, densfun="exponential")$estimate
  }
  Tscore_class[['exp_fit_rates']]<-exp_fit_rates
  exp_fit_rate<-mean(exp_fit_rates)
  #hist(exp_fit_rates, main=exp_fit_rate)
  #abline(v=exp_fit_rate, col='red')
  
  #calculate probability:  probabilities are P[X > x]
  P<-pexp(RC_1$RC, rate=exp_fit_rate, lower.tail=F) 
  P[P<1e-320|P==0]<-1e-320
  RC_1[,'pvalue']<-P
  #Tscore calculation
  Tscore<-round(-10*log(P,base=10), digits=1)
  head(sort(Tscore,decreasing=T))
  RC_1[,'Tscore']<-Tscore
  #BH-adjusted p values
  BH.P<-p.adjust(P, "BH")
  RC_1[,'BH.pvalue']<-BH.P
  #Tscore calculation
  BH.Tscore<-round(-10*log(BH.P,base=10), digits=1)
  head(sort(BH.Tscore,decreasing=T))
  RC_1[,'BH.Tscore']<-BH.Tscore
  
  #sig_RC<-subset(RC_1, Tscore>=20)
  #sig_RC<-sig_RC[order(sig_RC$Tscore,decreasing=T),]
  #dim(sig_RC) 
  #head(sig_RC, n=10)
  
  
  Tscore_class[['Tscore_df']]<-RC_1
  return(Tscore_class)
}



###############################
#main programming
################################
#options for Tscore cal
size_bin_info<-read.csv(size_bin_csv, header=T, row.names=1)
bin_type<-ifelse(grepl('Size', RC_csv_name), 'size_bin_name', 'enzyme_bin_name')
interactions_type<-ifelse(grepl('cis_ligation', RC_csv_name), 'cis', 'trans')

#Tscore calculation
print(paste('Calculate Tscores of ', interactions_type, ' interactions of ', RC_csv_name, sep="") )
if(interactions_type== 'cis'){
  Tscore_class<-cis_Tscore_calculation(in_file=in_file, RC.noise=RC.noise,
                     permutation.num=permutation.num, sampling.ratio=sampling.ratio)
}else {
  Tscore_class<-trans_Tscore_calculation(in_file=in_file, RC.noise=RC.noise, 
                    permutation.num=permutation.num, sampling.ratio=sampling.ratio)
}

#export R data
print("Save R object Tscore_class to the file")
RData_file_name<-sub(pattern="_RC.csv", replacement="_Tscore.RData", RC_csv_name)
RData_file<-paste(ligations_dir, RData_file_name, sep='/')
save(Tscore_class, file = RData_file )

#export data frome into csv file
print("Save the data frame Tscore into the file")
Tscore_csv_name<-sub(pattern="_RC.csv", replacement="_Tscore.csv", RC_csv_name)
Tscore_csv<-paste(ligations_dir, Tscore_csv_name, sep='/')
write.csv(Tscore_class$Tscore_df, file=Tscore_csv, row.names=F, quote=F)


####################################
#end
####################################
