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
	print "\n###Analysis of $raw_file are triggered.\n";
		
	#initiate starting time
	$variables_pointer=sub_basic::initiate_starting_time(\%variables, $raw_file);
	%variables=%$variables_pointer;
	
	if($variables{'3C_trimming'} eq 'yes'){
		my $R1_fastq_file=$raw_file;
		my $trimmed_R1_fastq_file=$variables{sample_out_file}.'.fastq_trimmed';
		my $min_adapter_3=substr($variables{adapter_3}, 0, $variables{adapter_len});
		print"###Trim adapter and poor quality sequences of $R1_fastq_file into $trimmed_R1_fastq_file.\n";
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
	if($variables{'3C_colocalization'} eq 'yes'){
		print" \n###Pair-end fastq files alignment  of $variables{sample_out_file} begin:\n";
		sub_3C::S2_Colocalization_3CMTS(\%variables);
	}
	
	if($variables{'3C_other_detection'} eq 'yes'){
		#read one alignment file
		my $other_alignment_file=$variables{sample_out_file}.'.other_alignment';
		print" Extract co-localizations from $other_alignment_file:\n";
		my $pairs_pointer=sub_3C::read_pairs_output(\%variables, $other_alignment_file);
		my %pairs=%$pairs_pointer;
		#
		sub_3C::S3_detect_other_alignment_3CMTS(\%variables, \%pairs);
	}
	
	#record ending time of a sample analysis
	sub_basic::initiate_ending_time($variables{sample_log}, $variables{beginning_time});
		
	print "\n###Analysis of $raw_file is done.\n\n";
}
############


#####################################################
#main programe
#####################################################

#get the directory of perl scripts
my $perl_dir=Cwd::getcwd();

#get subroutines
require $perl_dir."/functions_basic.pm"; #sub_basic::
require $perl_dir."/functions_3C.pm";  #sub_3C::


#read vairables from 3C_info.txt
my $var_file=$ARGV[0];
my $variables_pointer=sub_3C::process_info($var_file);
my %variables=%$variables_pointer;
$variables{Tscore_R_file}=$perl_dir.'/R_3C_Tscore.R';
sub_basic::print_hash(\%variables);

#
print "\n\n########## ---Pre-procesing--- ##########\n";

#
print "Set connection between sample names and raw data from the file $variables{sample_info_file}.\n\n";
my $sample_info_pointer=sub_3C::Pre_sample_info(\%variables);	#read sample_information
%variables=%$sample_info_pointer;

#
print "initiate and clear the monitor.log\n\n";
my $total_beginning_time=sub_basic::initiate_log_files(\%variables);

#generate the information of reference chromosomes
print "Read the sequences of the reference genome!\n\n";
my $size_bin_pointer=sub_3C::Pre_size_bin(\%variables);
my %size_bins=%$size_bin_pointer;
$variables{size_bins_pointer}=\%size_bins;

#generate primary enzyme sites along reference chromosomes;
print "Position of primary enzyme sites of $variables{enzyme_site}\n\n";
my $enzyme_bins_pointer=sub_3C::Pre_enzyme_bin(\%variables);
my %enzyme_bins=%$enzyme_bins_pointer;
$variables{enzyme_bins_pointer}=\%enzyme_bins;

my $enzyme_bin_regions_pointer=sub_3C::Pre_enzyme_bin_regions(\%enzyme_bins);
my %enzyme_bin_regions=%$enzyme_bin_regions_pointer;
$variables{enzyme_bin_regions_pointer}=\%enzyme_bin_regions;

#generate the file of target site information and the file of chromosome-enzyme site information
print "Get sites information from $variables{site_info_file} and generate $variables{probe_seq_csv}\n\n";
my $viewpoint_info_pointer=sub_3C::Pre_viewpoint_info_3CMTS(\%variables);
my %viewpoint_info=%$viewpoint_info_pointer;
$variables{viewpoint_info_pointer}=\%viewpoint_info;

#check bowtie index
sub_basic::check_bowtie_index(\%variables, 'bowtie2');

print "############################\n";
#############
########################################################################

if ($variables{'3C_trimming'} eq 'yes' or $variables{'3C_colocalization'} eq 'yes' or 
				$variables{'3C_other_detection'} eq 'yes'){#1

	print "\n\n\n########## ---Main running--- #########\n\n";
	my $raw_files_input_pointer=$variables{raw_files_input_pointer};
	my @raw_files_input=@$raw_files_input_pointer;
	while(1){#2
		if(threads->list()<$variables{threads_num} and @raw_files_input>0){
			my $raw_file=shift @raw_files_input;
			threads->create(\&main_running, \%variables, $raw_file);  #subroutine
		}
		#recover all threads
		foreach my $sub_thread( threads->list() ){#2
			$sub_thread->join() if $sub_thread->is_joinable();
		}#2
		last if threads->list()==0 and @raw_files_input==0;
		#sleep 10;
	}#2
}#1
	
print "\n\n\n########## ---Export results--- ##########\n\n";
if ($variables{'3C_statistics'} eq 'yes'){
	print "Statistics begin!\n\n";
	sub_3C::S4_result_statistics(\%variables);
}


if($variables{'3C_RC_counting'} eq 'yes'){#1
	print "###export counting table by viewpoints index.\n";
	#0=coordinate, 1=vp_displayid, 2=vp_chr, 3=vp_pos, 4=vp_sbn, 5=vp_ebn, 
	#6=ref, 7=start, 8=end, 9=start_sbn, 10=end_sbn, 11=start_ebn, 12=end_ebn, 13=seq, 14=valid,
	my @export_ref=( 
									{ row=>'6,9', r_row=>'6,10', row_attr=>'SizeBin', col=>'1', 
										col_attr=>'Viewpoints', thread=>1,	},
									{ row=>'6,11', r_row=>'6,12', row_attr=>'EnzymeBin', col=>'1', 
										col_attr=>'Viewpoints', thread=>2,	},
		);
	
	#run
	while(1){#2
		if(threads->list()<$variables{threads_num} and @export_ref>0){
			my $export_pointer=shift @export_ref;
			my %export=%$export_pointer;
			threads->create(\&sub_3C::S5_RC_counting_3CMTS, \%variables, \%export);  #subroutine
		}
		#recover all threads
		foreach my $sub_thread( threads->list() ){
			$sub_thread->join() if $sub_thread->is_joinable();
		}
		last if threads->list()==0 and @export_ref==0;
	}#2
	
}#1
#Tscore calculation
if ($variables{'3C_Tscore'} eq 'yes'){
	sub_3C::S6_Tscore_calculation( \%variables);  #subroutine
}

#ending time
sub_basic::initiate_ending_time($variables{total_log_file}, $total_beginning_time);

print "\n\n########## ---3C analysis is done!--- ##########\n\n";






