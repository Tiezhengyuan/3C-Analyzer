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
	my ($variables_pointer, $sample_name)=@_;
	my %variables=%$variables_pointer;
	
	#refresh sample dir
	$variables{result_dir}=$variables{result_dir}.'/'.$sample_name;
	mkdir($variables{result_dir}, 0755) unless -d $variables{result_dir};
	$variables{sample_log}=$variables{result_dir}.'/'.$sample_name.'.log';
	system("touch $variables{sample_log}") unless -f $variables{sample_log}; #create sample log file
	sub_3C::refresh_log($variables{sample_log}, "sample_name", $sample_name);
	print "###Note: The result directory is temporarily changed to $variables{result_dir} for $sample_name\n\n";
	
	#record beginning time of a sample analysis
	print "\n###Analysis of $sample_name are triggered.\n";
	my $status=sub_3C::read_log($variables{time_monitor_file}, "$sample_name:status");
	if ( $status eq 'no'){
		$status='on';
		sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:status", $status);
		my @beginning_time=localtime(time);
		sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:beginning_time" , join(",", @beginning_time));
	}
	#
	if($variables{'4C_colocalization'} eq 'yes'){
		print"\n###Co-localization of $sample_name.\n";
		sub_3C::C4_colocalization(\%variables, $sample_name);
	}
	
	#record ending time of a sample analysis
	if ($status eq 'on'){
		my $beginning_time_str=sub_3C::read_log($variables{time_monitor_file}, "$sample_name:beginning_time", $variables{lock_file});
		my @beginning_time= split(",", $beginning_time_str);
		my @ending_time=localtime(time);
		my $running_time=sub_3C::get_time(\@beginning_time, \@ending_time);
		sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:ending_time" , join(",", @ending_time));
		sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:running_time", $running_time);
		$status='yes';
		sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:status", $status);
	}
	
	print "\n###Analysis of $sample_name is done.\n\n";
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
$variables{lock_file}=$variables{result_dir}."/lock_file.log";
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
print "raw data directory: $variables{raw_data_dir}\n";  
print "The file list of raw data:\t";
my $raw_files_pointer=sub_3C::raw_files_list($variables{raw_data_dir}, $variables{raw_file_format});
my @raw_files=@$raw_files_pointer;
my $raw_files_num=@raw_files;
my @raw_files_input=@raw_files;
@raw_files_input=grep(/_R1_/, @raw_files_input) if $variables{raw_file_format} eq 'FASTQ' ;
print "A total of $raw_files_num will be analyzed.\n\n";

#
print "Set connection between sample names and raw data from the file $variables{sample_info_file}.\n\n";
my $out_pointer=sub_3C::Pre_read_sample_info($variables{sample_info_file}, \@raw_files, $variables{raw_file_format});	#read sample_information
my %out=%$out_pointer;
my $raw_to_sample_pointer=$out{raw_to_sample_pointer}; 
my %raw_to_sample=%$raw_to_sample_pointer;
$variables{raw_to_sample_pointer}=\%raw_to_sample;
my $sample_info_pointer=$out{sample_info_pointer}; 
my %sample_info=%$sample_info_pointer;
$variables{sample_info_pointer}=\%sample_info;
#sample names
my @sample_names=sort( keys %sample_info);
$variables{sample_names}=join(",", @sample_names);
sub_3C::refresh_log($variables{var_file}, 'sample_names', $variables{sample_names});

print "initiate and clear the monitor.log\n\n";
sub_3C::Pre_initiate_monitor_log($variables{time_monitor_file}, $variables{system_monitor_file}, \%sample_info);

#generate the information of reference chromosomes
print "Read the sequences of the reference genome!\n\n";
print "Length of sized bin is $variables{size_bin_len}\n";
my $chr_info_pointer=sub_3C::Pre_size_bin(\%variables);
my %chr_info=%$chr_info_pointer;
$variables{chr_info_pointer}=\%chr_info;

#generate primary enzyme sites along reference chromosomes;
print "Position of primary enzyme sites of $variables{enzyme_site}\n\n";
my $enzyme_sites_pointer=sub_3C::Pre_enzyme_bin(\%variables);
my %enzyme_sites=%$enzyme_sites_pointer;
$variables{enzyme_sites_pointer}=\%enzyme_sites;

#generate the file of target site information and the file of chromosome-enzyme site information
#######one 4C library and one fastq and one viewpoint for 4C-seq!!!!!!!!!!!!!
print "Get sites information from $variables{site_info_file} and generate $variables{viewpoints_csv}\n\n";
my $site_info_pointer=sub_3C::Pre_4C_site_info(\%variables);
my %site_info=%$site_info_pointer;
$variables{site_info_pointer}=\%site_info;

print "############################\n";

#
#print "Warning! Calculate the free space in $variables{result_dir}: ";
#my $space_info=split(",", sub_3C::Pre_free_space($variables{raw_data_dir}, $variables{result_dir}) );
#print "$space_info\n\n";

#check bowtie index
#build index
unless (-f $variables{genome_index}.'.1.bt2'){
	print "Build bowtie index of the genome reference---$variables{genome_fasta_file}.\n";
	my $bowtie_build_script=$variables{alignment_dir}.'/bowtie2-build';
	system("$bowtie_build_script $variables{genome_fasta_file} $variables{genome_index}");
}

#############
########################################################################

print "\n\n\n########## ---Main running--- #########\n\n";
while(1){
	if(threads->list()<$variables{threads_num} and @sample_names>0){
		my $sample_name=shift @sample_names;
		threads->create(\&main_running, \%variables, $sample_name);  #subroutine
	}
	#recover all threads
	foreach my $sub_thread( threads->list() ){#2
		$sub_thread->join() if $sub_thread->is_joinable();
	}#2
	last if threads->list()==0 and @sample_names==0;
	#sleep 10;
}#1

print "\n\n\n########## ---Export results--- ##########\n\n";
if ($variables{'4C_statistics'} eq 'yes'){
	print "Statistics begin!\n\n";
	sub_3C::C4_result_statistics(\%variables);
}


if($variables{'4C_RC_counting'} eq 'yes'){#1
	print "###export counting table by viewpoints index.\n";
	#0=displayid, 1=coordinate, 2=size_bin_no, 3=enzyme_bin_no
	#4=co_size_bin_no, 5=right_co_size_bin_no, 
	#6=co_enzyme_bin_no, 7=right_co_enzyme_bin_no,
	#8=co_ref, 9=co_offset, 10=co_seq, 11=tag,
	#displayid in 0 column is equal to sample_name, or viewpoints
	#col_1 always sample names
	my @export_ref=( 
		{ row_1=>8, row_2=>4, row_2r=>5, row_attr=>'SizeBin', col=>0, col_attr=>'Viewpoints', thread=>1,	},
		{ row_1=>8, row_2=>6, row_2r=>7, row_attr=>'EnzymeBin', col=>0, col_attr=>'Viewpoints', thread=>2,	},
		);
	
	#run
	while(1){#2
		if(threads->list()<$variables{threads_num} and @export_ref>0){
			my $export_pointer=shift @export_ref;
			my %export=%$export_pointer;
			threads->create(\&sub_3C::C4_Co_localization_counting, \%variables, \%export);  #subroutine
		}
		#recover all threads
		foreach my $sub_thread( threads->list() ){
			$sub_thread->join() if $sub_thread->is_joinable();
		}
		last if threads->list()==0 and @export_ref==0;
	}#2
	
}#1
#Tscore calculation
sub_3C::S6_Tscore_calculation( \%variables);  #subroutine

#
my @total_ending_time=localtime(time);
sub_3C::refresh_log($variables{time_monitor_file}, "Total:ending_time" , join(",", @total_ending_time) );
my $total_running_time=sub_3C::get_time(\@total_beginning_time, \@total_ending_time);
sub_3C::refresh_log($variables{time_monitor_file}, "Total:running_time" , $total_running_time);

print "\n\n########## ---3C analysis is done!--- ##########\n\n";






