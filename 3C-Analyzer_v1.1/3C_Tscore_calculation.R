#
####################################
#version 1.0
####################################

####################################
#library
####################################
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

##########################################
###Tscore calculation
Tscore_calculation<-function(in_file, viewpoints, bin_type, p.level=0.01, FDR.level=0.05, RC.noise=5){
	if(RC.noise<1){	RC.noise<-1	}
	print('read RC data in .RData')
	bin_df<-read_RC(in_file)
	bin_chr<-bin_df[,1]
	bin_no<-bin_df[,2]
	#
	estimate<-as.data.frame(matrix(NA, nrow=10, ncol=ncol(bin_df[,-(1:2)])))
	colnames(estimate)<-colnames(bin_df[,-(1:2)])
	rownames(estimate)<-c('cis_mean', 'cis_RC0', 'cis_k0', 'cis_RC1', 'cis_k1', 'cis_rate', 'cis_theoRC',
												 'trans_mean', 'trans_rate', 'trans_theoRC')
	#
	print('setup data frame')
	theoRC_df<-as.data.frame(matrix(0, nrow=nrow(bin_df), ncol=ncol(bin_df)))
	colnames(theoRC_df)<-colnames(bin_df)
	rownames(theoRC_df)<-rownames(bin_df)
	theoRC_df[,1]<-bin_df[,1]
	theoRC_df[,2]<-bin_df[,2]
	distRC_df<-theoRC_df
	pval_df<-as.data.frame(matrix(1, nrow=nrow(bin_df), ncol=ncol(bin_df)))
	colnames(pval_df)<-colnames(bin_df)
	rownames(pval_df)<-rownames(bin_df)
	pval_df[,1]<-bin_df[,1]
	pval_df[,2]<-bin_df[,2]
	padj_df<-pval_df
	#qval_df<-pval_df

	print('Enter loop')
	for(i in 3:ncol(bin_df) ){#2
		view<-names(bin_df[i])
		view_chr<-viewpoints[view,'chr']
		view_bin<-viewpoints[view, bin_type]
		#print(paste(i-2,': ', view, sep=""))
		#cis-interaction
		raw_cis<-subset(bin_df[,c('chr','bin',view)], chr==view_chr)
		colnames(raw_cis)<-c('chr', 'bin','rawRC')
		min_raw_cis<-subset(raw_cis, rawRC>=RC.noise)
		min_judging<-var(min_raw_cis$rawRC)>mean(min_raw_cis$rawRC)
		if(nrow(min_raw_cis)>20 & min_judging){#3
			ordered_rawRC<-sort(raw_cis[,3], decreasing=T)
 			min_rows<-ifelse(length(ordered_rawRC)>1000, 1000, length(ordered_rawRC))
			RC.level<-ifelse(ordered_rawRC[min_rows]>RC.noise, RC.noise, ordered_rawRC[min_rows]) #at least 1000 rows were reserved
			cis=subset(raw_cis, rawRC>=RC.level)
			#plot(x=cis[,2], y=cis[,3], xlim=c(view_bin-100, view_bin+100), pch=20)
		
			#filter large numbers
			RC.fold<-round(log(max(cis[,3])/1e3, base=10))
			if(RC.fold<0){RC.fold=0}
			cis[,4]<-cis[,3]/(10**RC.fold)
			cis[,5]<-abs(cis[,2]-view_bin)
			colnames(cis)<-c('chr', 'bin','rawRC', 'RC', 'dist')
			estimate['cis_mean',i-2]<-round(mean(cis$rawRC))
			#calculate RC0 and k of fit_RC by power distribution
			single_cis<-aggregate(RC~dist, data=cis, max)
			if(single_cis[1,'RC']<single_cis[2,'RC']){		single_cis<-single_cis[-1,]		} 
			cis_lm<-lm(RC~dist, data=single_cis)
			cis_nls<-nlxb(RC~RC0*exp(-k0*dist)+RC1*exp(-k1*dist), data=single_cis, start=list(RC0=abs(cis_lm$coeff[1]), k0=abs(cis_lm$coeff[2]), RC1=abs(cis_lm$coeff[1]), k1=abs(cis_lm$coeff[2])), control=list(maxiter=1000) )
			cis_RC0<-cis_nls$coeff[1]
			cis_k0<-cis_nls$coeff[2]
			cis_RC1<-cis_nls$coeff[3]
			cis_k1<-cis_nls$coeff[4]
			estimate['cis_RC0',i-2]<-cis_RC0
			estimate['cis_k0',i-2]<-cis_k0
			estimate['cis_RC1',i-2]<-cis_RC1
			estimate['cis_k1',i-2]<-cis_k1
			#export theoRC_df
			raw_cis[,4]<-abs(raw_cis[,2]-view_bin)
			colnames(raw_cis)<-c('chr', 'bin','rawRC', 'dist')
			raw_theoRC<-round( (cis_RC0*exp(-cis_k0*raw_cis$dist)+cis_RC1*exp(-cis_k1*raw_cis$dist))*(10**RC.fold) )
			names(raw_theoRC)<-rownames(raw_cis)
			theoRC_df[names(raw_theoRC),i]<-raw_theoRC
			#calculate distRC normalized by dist
			theoRC<-round( (cis_RC0*exp(-cis_k0*cis$dist)+cis_RC1*exp(-cis_k1*cis$dist))*(10**RC.fold) )
			names(theoRC)<-rownames(cis)
			distRC<-cis$rawRC-theoRC
			names(distRC)<-rownames(cis)
			distRC_df[names(distRC),i]<-distRC
			cis<-cbind(cis,theoRC, distRC)
			#fitted by exponential distribution
			distRC=distRC[distRC>RC.noise]
			if(length(distRC)>=10){#4
				#fit log2 of RC using exponential distribution
				exp_fit_rate<-fitdistr(distRC, densfun="exponential")$estimate
				estimate['cis_rate',i-2]<-exp_fit_rate
				#calculate cumulative probability:  probabilities are P[X > x]
				P<-pexp(distRC, rate=exp_fit_rate, lower.tail=F)
				names(P)<-names(distRC)
				P[P==0]<-1e-320
				P[P<1e-320]<-1e-320
				pval_df[names(distRC),i]<-P
				#BH-adjusted p values
				Padj<-p.adjust(P, "BH")
				padj_df[names(distRC),i]<-Padj
				# q values
				#qval<-qvalue(P, fdr.level=FDR.level, pi0.method="bootsrap")$qvalues			
				#RC limit based on exp_fit_df
				cis_theoRC<-qexp(p.level, rate=exp_fit_rate, lower.tail=F)
				estimate['cis_theoRC',i-2]<-round(cis_theoRC)
			}#4
			subset(cis, dist<25)
		}#3
		
		#trans-interaction
		trans<-subset(bin_df[,c(1,2,i)], chr!=view_chr)
		trans=trans[trans[,3]>RC.noise, ]
		estimate['trans_mean',i-2]<-round(mean(trans[,3]))
		if(nrow(trans)>=10){#3
			#fit log2 of RC using exponential distribution
			exp_fit_rate<-fitdistr(trans[,3], densfun="exponential")$estimate
			estimate['trans_rate',i-2]<-exp_fit_rate
			#calculate cumulative probability:  probabilities are P[X > x]
			P<-pexp(trans[,3], rate=exp_fit_rate, lower.tail=F)
			P[P==0]<-1e-320
			P[P<1e-320]<-1e-320
			pval_df[rownames(trans),i]<-P
			#BH-adjusted p values
			padj_df[rownames(trans),i]<-p.adjust(P, "BH")
			# q values
			#qval_df[rownames(trans),i]<-qvalue(P, fdr.level=FDR.level, pi0.method="bootsrap")$qvalues			
			#RC limit based on exp_fit_df
			theo_RC<-qexp(p.level, rate=exp_fit_rate, lower.tail=F)
			estimate['trans_theoRC',i-2]<-round(theo_RC)
		}#3
	}#2
	#Tscore calculation
	Tscore_df<-apply(pval_df[,-(1:2)],  MARGIN=2, FUN=function(x){
					x[x==0]=1
					Tscore<-round(-10*log(x,base=10), digits=1)
					return(Tscore)
			})
	Tscore_df<-cbind(pval_df[,(1:2)], Tscore_df)
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
size_bin_info<-read.csv(size_bin_csv, header=T, row.names=1)
bin_type<-ifelse(grepl('Size', RC_csv_name), 'size_bin_name', 'enzyme_bin_name')

#
print('Tscore calculation')

Tscore_class<-Tscore_calculation(in_file=in_file, viewpoints=viewpoints, bin_type=bin_type, p.level=p.level, RC.noise=RC.noise)

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





