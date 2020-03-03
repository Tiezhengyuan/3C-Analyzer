#! /usr/bin/perl -w
use strict;
use warnings;
use List::Util;
use List::MoreUtils;
use threads; 
use Cwd;


#the file is the leading script

################################################
#subroutines
################################################
sub main_running{
	my ($variables_pointer, $raw_file)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $raw_to_sample_pointer=$variables{raw_to_sample_pointer};
	my %raw_to_sample=%$raw_to_sample_pointer;
	my @beginning_time=localtime(time);
	
	#refresh sample_dir
	my $raw_name=$raw_to_sample{$raw_file}->{'raw_file_name_head'};
	my $sample_name=$raw_to_sample{$raw_file}->{'sample_name'}; #sample name
	$variables{sample_name}=$sample_name;
	$variables{sample_dir}=$sample_info{$sample_name}->{'sample_dir'};
	$variables{sample_log}=$variables{sample_dir}.'/'.$raw_name.'.log';
	print "###Note: The result directory is temporarily changed to $variables{sample_dir} for $sample_name\n\n";
	
	#refresh sample_log
	sub_3C::refresh_log($variables{sample_log}, "status" , 'on');
	sub_3C::refresh_log($variables{sample_log}, "sample_name", $sample_name);
	sub_3C::refresh_log($variables{sample_log}, "beginning_time" , join(",", @beginning_time));
	
	
	print "\n###Analysis of $raw_file of $sample_name are triggered.\n";

	if($variables{'3C_trimming'} eq 'yes' and $variables{'raw_file_format'} eq 'FASTQ'){
		sub_3C::refresh_log($variables{sample_log}, "fastq_files", 2);
		my $R1_fastq_file=$raw_file;
		my $trimmed_R1_fastq_file=$variables{sample_dir}.'/'.$raw_name.'.fastq_trimmed';
		print"\n###Trim adapter and poor quality sequences of $R1_fastq_file into $trimmed_R1_fastq_file.\n";
		my $min_adapter_3=substr($variables{adapter_3}, 0, $variables{adapter_len});
		sub_3C::S1_trim_fastq(\%variables, $R1_fastq_file, $trimmed_R1_fastq_file, $min_adapter_3);
		
		my $R2_fastq_file=$raw_file;
		$R2_fastq_file=~s/_R1_/_R2_/;
		my $trimmed_R2_fastq_file=$trimmed_R1_fastq_file;
		$trimmed_R2_fastq_file=~s/_R1_/_R2_/;
		my $min_adapter_5=substr($variables{adapter_5}, 0, $variables{adapter_len});
		print"###Trim adapter and poor quality sequences of $R2_fastq_file into $trimmed_R2_fastq_file.\n";
		sub_3C::S1_trim_fastq(\%variables, $R2_fastq_file, $trimmed_R2_fastq_file, $min_adapter_5);
	}
	
	print "\n###Co-localization detetction begin:\n";
	if($variables{'3C_alignment'} eq 'yes' and $variables{'raw_file_format'} eq 'FASTQ'){
		print" \n###Pair-end fastq files alignment  of $raw_name.\n";
		my $trimmed_fastq_file=$variables{sample_dir}.'/'.$raw_name.'.fastq_trimmed';
		if(-f $trimmed_fastq_file){
			sub_3C::S2_Colocalization_detection(\%variables, $trimmed_fastq_file, $raw_name);
		}
		else{
			sub_3C::S2_Colocalization_detection(\%variables, $raw_file, $raw_name);
		}
	}
	
	if ($variables{'raw_file_format'} eq 'BAM'){
		print "\n###transfer BAM format of $raw_name into SAM format:\n";
		my $samtools_script=$variables{alignment_dir}.'/samtools';
		my $BAM_file=$variables{sample_dir}.'/'.$raw_name.".bam";
		my $SAM_file=$variables{sample_dir}.'/'.$raw_name.".sam";
		system("samtools_script view $BAM_file > $SAM_file");
	}
	
	if($variables{'3C_multiple_detection'} eq 'yes'){
		print" Extract information from multiple alignment from $raw_name:\n";
		sub_3C::S3_site_detection_multiple_alignment(\%variables, $raw_name);
	}
	if($variables{'3C_one_detection'} eq 'yes'){
		print" Extract information from one_alignment from $raw_name:\n";
		sub_3C::S3_site_detection_one_alignment(\%variables, $raw_name);
	}
	if($variables{'3C_no_detection'} eq 'yes'){
		print" Extract information from unalignment from $raw_name:\n";
		sub_3C::S3_site_detection_no_alignment(\%variables, $raw_name);
	}
	
	#record ending time of a sample analysis
	my @ending_time=localtime(time);
	my $running_time=sub_3C::get_time(\@beginning_time, \@ending_time);
	sub_3C::refresh_log($variables{sample_log}, "ending_time" , join(",", @ending_time));
	sub_3C::refresh_log($variables{sample_log}, "running_time", $running_time);
	sub_3C::refresh_log($variables{sample_log}, "status" , 'off');
	
	
	print "\n###Analysis of $raw_name is done.\n\n";
}
############


#####################################################
#main programe
#####################################################

#get the directory of perl scripts involved in Pscore
my $perl_dir=Cwd::getcwd();

#get subroutines 
require $perl_dir."/3C_subroutines.pm";

#read vairables from 3C_info.txt
my $var_file=$ARGV[0];
my $variables_pointer=sub_3C::process_info($var_file);
my %variables=%$variables_pointer;
$variables{ligations_dir}=$variables{result_dir}."/genomic_ligations";
$variables{ligation_frequency_dir}=$variables{result_dir}."/ligation_frequency";
$variables{log_dir}=$variables{result_dir}.'/sample_log';
mkdir($variables{log_dir}, 0755) unless -d $variables{log_dir};
$variables{genome_info_csv}=$variables{result_dir}.'/genome_info.csv'; 
$variables{probe_seq_csv}=$variables{result_dir}.'/probe_seq.csv'; #from $variables{site_info_file}
$variables{enzyme_bin_csv}=$variables{result_dir}.'/genome_enzyme_bin.csv'; 
$variables{enzyme_sites_file}=$variables{result_dir}.'/genome_enzyme_sites.txt'; 
$variables{size_bin_csv}=$variables{result_dir}.'/genome_size_bin.csv';
$variables{viewpoints_csv}=$variables{result_dir}.'/viewpoints.csv';
$variables{statistics_file}=$variables{result_dir}.'/statistics.txt';
$variables{Tscore_R_file}=$perl_dir.'/3C_Tscore_calculation.R';

#starting time
my @total_beginning_time=localtime(time);

#
print "\n\n########## ---Pre-procesing--- ##########\n";
#
print "Set connection between sample names and raw data from the file $variables{sample_info_file}.\n\n";
$variables_pointer=sub_3C::Pre_read_sample_info(\%variables);	#read sample_information
%variables=%$variables_pointer;
my @sample_names=split(',', $variables{sample_names});

print "initiate and clear the monitor.log\n\n";
sub_3C::Pre_initiate_monitor_log(\%variables);

#generate the information of reference chromosomes
print "Read the sequences of the reference genome!\n\n";
my $chr_info_pointer=sub_3C::Pre_size_bin(\%variables);
my %chr_info=%$chr_info_pointer;
$variables{chr_info_pointer}=\%chr_info;

#generate primary enzyme sites along reference chromosomes;
print "Position of primary enzyme sites of $variables{enzyme_site}\n\n";
my $enzyme_sites_pointer=sub_3C::Pre_enzyme_bin(\%variables);
my %enzyme_sites=%$enzyme_sites_pointer;
$variables{enzyme_sites_pointer}=\%enzyme_sites;

#generate the file of target site information and the file of chromosome-enzyme site information
print "Get sites information from $variables{site_info_file} and generate $variables{probe_seq_csv}\n\n";
my $site_info_pointer=sub_3C::Pre_site_info(\%variables);
my %site_info=%$site_info_pointer;
$variables{site_info_pointer}=\%site_info;

#check bowtie index
sub_3C::Pre_check_bowtie_index(\%variables, 'bowtie2');

print "############################\n";
########################################################################

print "\n\n\n########## ---Main running--- #########\n\n";
my $R1_raw_files_pointer=$variables{R1_raw_files_pointer};
my @raw_files_input=@$R1_raw_files_pointer;

if ($variables{'3C_alignment'} eq 'yes' or $variables{'3C_multiple_detection'} eq 'yes' or $variables{'3C_no_detection'} eq 'yes' or $variables{'3C_one_detection'} eq 'yes'){#1
	while(1){#2
		$variables{threads_num}=sub_3C::dynamic_threads(\%variables);
		if(threads->list()<$variables{threads_num} and @raw_files_input>0){
			my $raw_file=shift @raw_files_input;
			threads->create(\&main_running, \%variables, $raw_file);  #subroutine
			sub_3C::refresh_sample_log(\%variables);
		}
		#recover all threads
		foreach my $sub_thread( threads->list() ){#2
			$sub_thread->join() if $sub_thread->is_joinable();
		}#2
		last if threads->list()==0 and @raw_files_input==0;
		sleep 10;
	}#2
	sub_3C::refresh_sample_log(\%variables);
}#1

print "\n\n\n########## ---Export results--- ##########\n\n";
if ($variables{'3C_statistics'} eq 'yes'){
	print "Statistics begin!\n\n";
	sub_3C::S4_result_statistics(\%variables);
}


if($variables{'3C_RC_counting'} eq 'yes'){#1
	print "###export counting table by viewpoints index.\n";
	#0=displayid, 1=coordinate, 2=size_bin_no, 3=enzyme_bin_no, 4=co_size_bin_no, 5=right_co_size_bin_no, 
	#6=co_enzyme_bin_no, 7=right_co_enzyme_bin_no,
	#8=co_ref, 9=col_offset, 10=col_seq, 11=tag, 12=sample_name
	#12:sample_name is not in the file but will be added in function S5_Co_localization_counting 
	#col_1 always sample names
	my @export_ref=( 
			{row_1=>8, row_2=>4, row_2r=>5, row_attr=>'SizeBin', 
				col_1=>12, col_2=>0, col_attr=>'Viewpoints', thread=>1,	},
			{row_1=>8, row_2=>6, row_2r=>7, row_attr=>'EnzymeBin', 
				col_1=>12, col_2=>0, col_attr=>'Viewpoints', thread=>2,	},
		);
	
	#run
	while(1){#2
		if(threads->list()<$variables{threads_num} and @export_ref>0){
			my $export_pointer=shift @export_ref;
			my %export=%$export_pointer;
			threads->create(\&sub_3C::S5_Colocalization_counting, \%variables, \%export);  #subroutine
		}
		#recover all threads
		foreach my $sub_thread( threads->list() ){
			$sub_thread->join() if $sub_thread->is_joinable();
		}
		last if threads->list()==0 and @export_ref==0;
	}#2
	
}#1
#Tscore calculation
if($variables{'3C_Tscore_calculation'} eq 'yes'){#1
	sub_3C::S6_Tscore_calculation( \%variables);  #subroutine
}
#
my @total_ending_time=localtime(time);
my $total_running_time=sub_3C::get_time(\@total_beginning_time, \@total_ending_time);
sub_3C::refresh_log($variables{time_monitor_file}, "Total:ending_time" , join(",", @total_ending_time) );
sub_3C::refresh_log($variables{time_monitor_file}, "Total:running_time" , $total_running_time);
sub_3C::refresh_log($variables{time_monitor_file}, "Total:status" , 'off');

print "\n\n########## ---3C analysis is done!--- ##########\n\n";






