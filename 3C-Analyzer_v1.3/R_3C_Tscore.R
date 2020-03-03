
######3C-analyzer#######

####################################
#version 1.1, 2015.04.15
#Developer: Tiezheng Yuan Ph.D. tiezhengyuan@hotmail.com
####################################

####################################
#library
####################################
if (!require("MASS")) {
  install.packages("MASS", dependencies = TRUE)
}
if (!require("nlmrt")) {
  install.packages("nlmrt", dependencies = TRUE)
}
library(MASS)
library(nlmrt)
#library(sqldf)
#library(qvalue)
#library(plotrix)

####################################
#functions
####################################
read_RC<-function(in_file){
	#read genomic-interactions
	bin<-read.table(in_file, header=T, sep=",")
	#bin<-read.csv.sql(in_file, sql="select * from file", header=T, row.names=F, sep=',')

	#order based chromosome
	#index <- with(bin, order(as.numeric(bin[,1]), as.numeric(bin[,2])))
	bin<-bin[order(as.numeric(bin[,1]), as.numeric(bin[,2])), ]
	names(bin)[1]<-'chr'
	names(bin)[2]<-'bin'
	
	return(bin)
}
#################################
cis_regression<-function(raw_cis, view_bin, RC.noise){
  cis.reg.class<-list()
  
  #
  ordered_rawRC<-sort(raw_cis[,'RC'], decreasing=T)
  min_rows<-ifelse(length(ordered_rawRC)>1000, 1000, length(ordered_rawRC))
  RC.level<-ifelse(ordered_rawRC[min_rows]>RC.noise, RC.noise, ordered_rawRC[min_rows]) #at least 1000 rows were reserved
  #
  cis=subset(raw_cis, RC>=RC.level)
  #plot(x=cis[,2], y=cis[,3], xlim=c(view_bin-100, view_bin+100), pch=20)
  
  #filter large numbers
  RC.fold<-round(log(max(cis[,'RC'])/1e4, base=10))
  if(RC.fold<0){RC.fold=0}
  cis[,'rRC']<-cis[,'RC']/(10**RC.fold)
  cis[,'dist']<-abs(cis[,'bin']-view_bin)
  #calculate RC0 and k of fit_RC by power distribution
  single_cis<-aggregate(rRC~dist, data=cis, max)
  #remove the viewpoint bin if lower signal
  if(single_cis[1,'rRC']<single_cis[2,'rRC']){		single_cis<-single_cis[-1,]		} 
  #non-linear regression
  cis_lm<-lm(rRC~dist, data=single_cis)
  cis_nls<-nlxb(rRC~RC0*exp(-k0*dist)+RC1*exp(-k1*dist), data=single_cis, start=list(RC0=abs(cis_lm$coeff[1]), k0=abs(cis_lm$coeff[2]), RC1=abs(cis_lm$coeff[1]), k1=abs(cis_lm$coeff[2])), control=list(maxiter=1000) )
  cis_RC0<-cis_nls$coeff[1]
  cis_k0<-cis_nls$coeff[2]
  cis_RC1<-cis_nls$coeff[3]
  cis_k1<-cis_nls$coeff[4]

  #calculate distRC normalized by dist
  theoRC<-round( (cis_RC0*exp(-cis_k0*cis$dist)+cis_RC1*exp(-cis_k1*cis$dist))*(10**RC.fold) )
  names(theoRC)<-rownames(cis)
  cis[,'theoRC']<-theoRC
  #distRC for regression
  distRC<-cis[,'RC']-theoRC
  names(distRC)<-rownames(cis)
  cis[,'distRC']<-distRC
  
  #export theoRC and distRC based on raw_cis
  raw_theoRC<-round( (cis_RC0*exp(-cis_k0*raw_cis$dist)+cis_RC1*exp(-cis_k1*raw_cis$dist))*(10**RC.fold) )
  names(raw_theoRC)<-rownames(raw_cis)
  raw_distRC<-raw_cis[,'RC']-raw_theoRC
  names(raw_distRC)<-rownames(raw_cis)
  
  #subset(cis, dist<25)
  #cis.reg.class[['cis']]<-cis
  cis.reg.class[['cis_est']]<-c(cis_RC0, cis_k0, cis_RC1, cis_k1,round(mean(cis[,'RC'])))
  cis.reg.class[['theoRC']]<-theoRC
  cis.reg.class[['distRC']]<-distRC
  cis.reg.class[['raw_theoRC']]<-raw_theoRC
  cis.reg.class[['raw_distRC']]<-raw_distRC
  return(cis.reg.class)
}

##########
sig_analysis<-function(RC, RC.noise){
  
  #initiate
  exp_fit_rate=0
  p.values<-rep(0, length=length(RC))
  names(p.values)<-names(RC)
  p.adj<-p.values
  sigRC.limit<-0
  
  #significance analysis
  reg.RC=RC[RC>RC.noise]
  if(length(reg.RC)>=10){#2
    #fit log2 of RC using exponential distribution
    exp_fit_rate<-fitdistr(reg.RC, densfun="exponential")$estimate
    
    #calculate cumulative probability:  probabilities are P[X > x]
    p.values<-pexp(RC, rate=exp_fit_rate, lower.tail=F)
    p.values[p.values==0|p.values<1e-320]<-1e-320
    #BH-adjusted p values
    p.adj<-p.adjust(p.values, "BH")
    # q values
    #qval<-qvalue(P, fdr.level=FDR.level, pi0.method="bootsrap")$qvalues			
    #RC limit based on exp_fit_df
    sigRC.limit<-qexp(p.level, rate=exp_fit_rate, lower.tail=F)
  }#2

  sig.class<-list()
  sig.class[['exp.rate']]<-exp_fit_rate
  sig.class[['p.values']]<-p.values
  sig.class[['p.adj']]<-p.adj
  sig.class[['sigRC.limit']]<-sigRC.limit
  return(sig.class)
}
###############
#Tscore calculation
convert_Tscore<-function(pval_df){

  Tscore_df<-apply(pval_df[,-(1:4)],  MARGIN=2, FUN=function(x){
    x[x==0]=1
    Tscore<-round(-10*log(x,base=10), digits=1)
    return(Tscore)
  })
  Tscore_df<-cbind(pval_df[,(1:4)], Tscore_df)
  
  return(Tscore_df)
}
##########################################
###Tscore calculation
Tscore_calculation<-function(in_file, viewpoints, bin_type, p.level=0.01, FDR.level=0.05, RC.noise=5){
	if(RC.noise<1){	RC.noise<-1	}
	print('read RC data in .RData')
	bin_df<-read_RC(in_file)
	bin_chr<-bin_df[,1]
	bin_no<-bin_df[,2]
	#
	estimate<-as.data.frame(matrix(NA, nrow=10, ncol=ncol(bin_df[,-(1:4)])))
	colnames(estimate)<-colnames(bin_df[,-(1:4)])

  #
	print('setup data frame')
	theoRC_df<-as.data.frame(matrix(0, nrow=nrow(bin_df), ncol=ncol(bin_df)))
	colnames(theoRC_df)<-colnames(bin_df)
	rownames(theoRC_df)<-rownames(bin_df)
	theoRC_df[,1:4]<-bin_df[,1:4]
	distRC_df<-theoRC_df
	pval_df<-as.data.frame(matrix(1, nrow=nrow(bin_df), ncol=ncol(bin_df)))
	colnames(pval_df)<-colnames(bin_df)
	rownames(pval_df)<-rownames(bin_df)
	pval_df[,1:4]<-bin_df[,1:4]
	padj_df<-pval_df
	#qval_df<-pval_df

	print('Enter loop')
	for(i in 5:ncol(bin_df) ){#2
		view<-names(bin_df[i])
		view_chr<-viewpoints[view,'chr']
		view_bin<-viewpoints[view, bin_type]
		#print(paste(i-2,': ', view, sep=""))
		#cis-interaction
		raw_cis<-subset(bin_df[,c('chr','bin',view)], chr==view_chr)
		colnames(raw_cis)<-c('chr', 'bin','RC')
		raw_cis[,'dist']<-abs(raw_cis[,2]-view_bin)
		min_raw_cis<-subset(raw_cis, RC>=RC.noise)
		min_judging<-var(min_raw_cis[,'RC'])>mean(min_raw_cis[,'RC'])
		if(nrow(min_raw_cis)>20 & min_judging){#3

      #function of cis regression
      cis.reg.class<-cis_regression(raw_cis, view_bin, RC.noise)
			#refresh data frames
      estimate[c('cis_RC0','cis_k0','cis_RC1','cis_k1','cis_mean'),i-4]<-cis.reg.class[['cis_est']]
			theoRC_df[rownames(raw_cis),i]<-cis.reg.class[['raw_theoRC']]
			distRC_df[rownames(raw_cis),i]<-cis.reg.class[['raw_distRC']]
      
      #function of significance analysis 
      cis.sig.class<-sig_analysis(cis.reg.class[['raw_distRC']], RC.noise)
      #refresh data frames
			estimate['cis_rate',i-4]<-cis.sig.class[['exp.rate']]
			estimate['cisRC_limit',i-4]<-cis.sig.class[['sigRC.limit']]
      pval_df[rownames(raw_cis),i]<-cis.sig.class[['p.values']]
			padj_df[rownames(raw_cis),i]<-cis.sig.class[['p.adj']]
			
		}#3
		
		#trans-interaction
		#function of significance analysis 
		trans<-subset(bin_df[,c(1,2,i)], chr!=view_chr)
    trans.sig.class<-sig_analysis(trans[,3], RC.noise)
		#refresh data frames
		estimate['trans_mean',i-4]<-round(mean(trans[,3]))
    estimate['trans_rate',i-4]<-trans.sig.class[['exp.rate']]
		estimate['transRC_limit',i-4]<-trans.sig.class[['sigRC.limit']]
		pval_df[rownames(trans),i]<-trans.sig.class[['p.values']]
		padj_df[rownames(trans),i]<-trans.sig.class[['p.adj']]
	}#2
	
  #Tscore conversion
	Tscore_df<-convert_Tscore(pval_df)

  #Tscore_class
	Tscore_class=vector(mode="list")
	Tscore_class$RC=bin_df #mode is data frame
	Tscore_class$theoRC=theoRC_df
	Tscore_class$distRC=distRC_df
	Tscore_class$pval=pval_df
	Tscore_class$padj=padj_df
	#Tscore_class$qval=qval_df
	Tscore_class$Tscore=Tscore_df
	Tscore_class$estimate=estimate
	
	return(Tscore_class)
}


####################################
#options
####################################

R_ligation_frequency_dir="/home/yuan/3C_analyzer/result/ligation_frequency"
R_RC_csv_name="EnzymeBin_Viewpoints_Total_RC.csv"
R_viewpoints_csv="/home/yuan/3C_analyzer/result/viewpoints.csv"
R_size_bin_csv="/home/yuan/3C_analyzer/result/genome_size_bin.csv"
R_enzyme_bin_csv="/home/yuan/3C_analyzer/result/genome_enzyme_bin.csv"
R_p_level="0.01"
R_RC_noise="1"

ligations_dir=R_ligation_frequency_dir
RC_csv_name=R_RC_csv_name
in_file=paste(ligations_dir, RC_csv_name, sep="/")
viewpoints_csv=R_viewpoints_csv
size_bin_csv=R_size_bin_csv
enzyme_bin_csv=R_enzyme_bin_csv
p.level=as.numeric(R_p_level)
RC.noise=as.numeric(R_RC_noise)
####################################
#main program
####################################
#options for Tscore cal
viewpoints<-read.csv(viewpoints_csv, header=T, row.names=1)
bin_type<-ifelse(grepl('Size', RC_csv_name), 'size_bin_name', 'enzyme_bin_name')

#
print('Tscore calculation')

Tscore_class<-Tscore_calculation(in_file=in_file, viewpoints=viewpoints, 
                    bin_type=bin_type, p.level=p.level, RC.noise=RC.noise)

#export
print("Save R object Tscore_class to the file")
RData_file_name<-sub(pattern="_RC.csv", replacement="_Tscore.RData", RC_csv_name)
RData_file<-paste(ligations_dir, RData_file_name, sep='/')
save(Tscore_class, file = RData_file )

print("Save the data frame Tscore into the file")
Tscore_csv_name<-sub(pattern="_RC.csv", replacement="_Tscore.csv", RC_csv_name)
Tscore_csv<-paste(ligations_dir, Tscore_csv_name, sep='/')
write.csv(Tscore_class$Tscore, file=Tscore_csv, row.names=F, quote=F)


####################################
#end
####################################





