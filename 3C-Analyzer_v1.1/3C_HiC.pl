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
	$variables{sample_dir}=$sample_info{$sample_name}->{sample_dir};
	$variables{sample_log}=$variables{result_dir}.'/'.$sample_name.'/'.$raw_name.'.log';
	print "###Note: The result directory is temporarily changed to $variables{sample_dir} for $sample_name\n\n";
	#refresh sample_log
	sub_3C::refresh_log($variables{sample_log}, "status" , 'on');
	sub_3C::refresh_log($variables{sample_log}, "sample_name", $sample_name);
	sub_3C::refresh_log($variables{sample_log}, "beginning_time" , join(",", @beginning_time));
	
	
	print "\n###Analysis of $raw_file of $sample_name are triggered.\n";
	#########
	if($variables{'HiC_colocalization'} eq 'yes'){
		sub_3C::HiC_S1_colocalization_detection(\%variables, $raw_file, $raw_name);
	}
	#########
	
	#record ending time of a sample analysis
	my @ending_time=localtime(time);
	my $running_time=sub_3C::get_time(\@beginning_time, \@ending_time);
	sub_3C::refresh_log($variables{sample_log}, "ending_time" , join(",", @ending_time));
	sub_3C::refresh_log($variables{sample_log}, "running_time", $running_time);
	sub_3C::refresh_log($variables{sample_log}, "status" , 'off');
	
	print "\n###Analysis of $raw_file is done.\n\n";
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
$variables{enzyme_bin_csv}=$variables{result_dir}.'/genome_enzyme_bin.csv'; 
$variables{enzyme_sites_file}=$variables{result_dir}.'/genome_enzyme_sites.txt'; 
$variables{size_bin_csv}=$variables{result_dir}.'/genome_size_bin.csv';
$variables{statistics_file}=$variables{result_dir}.'/statistics.txt';
$variables{Tscore_R_file}=$perl_dir.'/HiC_Tscore_calculation.R';

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

#check bowtie index
sub_3C::Pre_check_bowtie_index(\%variables, 'bowtie2');

print "############################\n";
########################################################################

print "\n\n\n########## ---Main running--- #########\n\n";
my $R1_raw_files_pointer=$variables{R1_raw_files_pointer};
my @raw_files_input=@$R1_raw_files_pointer;

if ($variables{'HiC_colocalization'} eq 'yes'){#1
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
if ($variables{'HiC_statistics'} eq 'yes'){
	print "Statistics begin!\n\n";
	sub_3C::HiC_S2_statistics(\%variables);
}


if($variables{'HiC_RC_counting'} eq 'yes'){#1
	mkdir($variables{ligation_frequency_dir}, 0755) unless -d $variables{ligation_frequency_dir};
	foreach my $sample_name(@sample_names){
		#cis
		my $cis_file=$variables{ligations_dir}.'/'.$sample_name.'.cis_ligation';
		sub_3C::HiC_S3_cisRC_counting(\%variables, $cis_file);  #subroutine
		#trans
		my $trans_file=$variables{ligations_dir}.'/'.$sample_name.'.trans_ligation';
		sub_3C::HiC_S3_transRC_counting(\%variables, $trans_file);  #subroutine
	}
}#1
#Tscore calculation
if($variables{'HiC_Tscore_calculation'} eq 'yes'){#1
	sub_3C::HiC_S4_Tscore_calculation( \%variables);  #subroutine
}
#
my @total_ending_time=localtime(time);
my $total_running_time=sub_3C::get_time(\@total_beginning_time, \@total_ending_time);
sub_3C::refresh_log($variables{time_monitor_file}, "Total:ending_time" , join(",", @total_ending_time) );
sub_3C::refresh_log($variables{time_monitor_file}, "Total:running_time" , $total_running_time);
sub_3C::refresh_log($variables{time_monitor_file}, "Total:status" , 'off');

print "\n\n########## ---3C analysis is done!--- ##########\n\n";






