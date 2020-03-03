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
	my $variables_time_pointer=sub_basic::initiate_starting_time(\%variables, $raw_file);
	%variables=%$variables_time_pointer;

	#########
	if($variables{'3C_colocalization'} eq 'yes'){
		sub_3C::S2_Colocalization_HiC(\%variables, $raw_file);
	}
	if($variables{'3C_other_detection'} eq 'yes'){
		#read other alignment file
		my $other_alignment_file=$variables{sample_out_file}.'.other_alignment';
		print" Extract co-localizations from $other_alignment_file:\n";
		my $pairs_pointer=sub_3C::read_pairs_output(\%variables, $other_alignment_file);
		my %pairs=%$pairs_pointer;
		#
		sub_3C::S3_detect_other_alignment_HiC(\%variables, \%pairs);
	}
	#########

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
$variables{Tscore_R_file}=$perl_dir.'/R_HiC_Tscore.R';
#print %variables
#sub_basic::print_hash(\%variables);
#
print "\n\n########## ---Pre-procesing--- ##########\n";
#
print "Set connection between sample names and raw data from the file $variables{sample_info_file}.\n\n";
my $variables_sample_pointer=sub_3C::Pre_sample_info(\%variables);	#read sample_information
%variables=%$variables_sample_pointer;

#
print "initiate and clear the monitor.log\n\n";
my $total_beginning_time=sub_basic::initiate_log_files(\%variables);

#generate information of reference chromosomes
print "Read the sequences of the reference genome!\n\n";
my $chr_info_pointer=sub_3C::Pre_chr_info(\%variables);
my %chr_info=%$chr_info_pointer;
$variables{chr_info_pointer}=\%chr_info;

#generate the size bins 
print "Length of sized bin is $variables{size_bin_len}\n";
my $size_bins_pointer=sub_3C::Pre_size_bin(\%variables);
my %size_bins=%$size_bins_pointer;
$variables{size_bins_pointer}=\%size_bins;

#generate primary enzyme sites along reference chromosomes;
print "Position of primary enzyme sites of $variables{enzyme_site}\n\n";
my $enzyme_bins_pointer=sub_3C::Pre_enzyme_bin(\%variables);
my %enzyme_bins=%$enzyme_bins_pointer;
$variables{enzyme_bins_pointer}=\%enzyme_bins;

my $enzyme_bin_regions_pointer=sub_3C::Pre_enzyme_bin_regions(\%enzyme_bins);
my %enzyme_bin_regions=%$enzyme_bin_regions_pointer;
$variables{enzyme_bin_regions_pointer}=\%enzyme_bin_regions;

#check bowtie index
sub_basic::check_bowtie_index(\%variables, 'bowtie2');

print "############################\n";
#############
########################################################################

print "\n\n\n########## ---Main running--- #########\n\n";
 my $raw_files_input_pointer=$variables{raw_files_input_pointer};
 my @raw_files_input=@$raw_files_input_pointer;
if ($variables{'3C_colocalization'} eq 'yes' or $variables{'3C_other_detection'} eq 'yes'){#1
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
	print "###export counting table.\n";
	my @sample_names=split(',', $variables{sample_names});
	my @export_ref=map { join(',', $_, 'cis_ligation', 'SizeBin') } @sample_names;
	push(@export_ref, map { join(',', $_, 'trans_ligation', 'SizeBin') } @sample_names);
	push(@export_ref, map { join(',', $_, 'cis_ligation', 'EnzymeBin') } @sample_names);
	push(@export_ref, map { join(',', $_, 'trans_ligation', 'EnzymeBin') } @sample_names);
	#run
	while(1){#2
		if(threads->list()<$variables{threads_num} and @export_ref>0){
			my $format=shift @export_ref;
			threads->create(\&sub_3C::S5_RC_counting_HiC, \%variables, $format );  #subroutine
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






