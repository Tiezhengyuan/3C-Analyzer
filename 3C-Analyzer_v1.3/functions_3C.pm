#! /usr/bin/perl -w
use strict;
use warnings;
#use threads;
use List::MoreUtils;
use List::Util;
use File::Find;


#the file constains all subroutines required for running the pipeline of 3C_analyzer

#get basic subroutines
my $perl_dir=Cwd::getcwd();
require $perl_dir."/functions_basic.pm"; #sub_basic::

#
package sub_3C;

##############################################
#open parameter info file named "variables.txt" needed for processing. 
sub process_info{#1
	my $var_file=$_[0];
	
	my $variables_pointer=sub_basic::file_to_hash($var_file, '=');
	my %variables=%$variables_pointer;  # save parameters for processing
	#
	$variables{ligations_dir}=$variables{result_dir}."/genomic_ligations";
	mkdir($variables{ligations_dir}, 0755) unless -d $variables{ligations_dir};
	$variables{ligation_frequency_dir}=$variables{result_dir}."/ligation_frequency";
	mkdir($variables{ligation_frequency_dir}, 0755) unless -d $variables{ligation_frequency_dir};
	$variables{sample_log_dir}=$variables{result_dir}."/sample_log";
	mkdir($variables{sample_log_dir}, 0755) unless -d $variables{sample_log_dir};
	#parameters of 3C
	#$variables{var_file}=$variables{result_dir}.'/variables.txt';
	#sample info
	$variables{sample_info_file}=$variables{result_dir}.'/sample_info.csv';
	#viewpoints
	$variables{viewpoints_csv}=$variables{result_dir}.'/viewpoints.csv';
	$variables{site_info_file}=$variables{result_dir}.'/site_info.csv';
	#genome 
	$variables{genome_info_csv}=$variables{result_dir}.'/genome_info.csv'; 
	$variables{probe_seq_csv}=$variables{result_dir}.'/probe_seq.csv'; #from $variables{site_info_file}
	$variables{enzyme_bin_csv}=$variables{result_dir}.'/genome_enzyme_bin.csv'; 
	$variables{enzyme_sites_file}=$variables{result_dir}.'/genome_enzyme_sites.txt'; 
	$variables{size_bin_csv}=$variables{result_dir}.'/genome_size_bin.csv';
	#monitor log
	$variables{log_file}=$variables{result_dir}.'/3C_analyzer.log'; 
	$variables{total_log_file}=$variables{sample_log_dir}."/Total.log";
	$variables{time_monitor_file}=$variables{result_dir}.'/time_monitor.log';
	$variables{system_monitor_file}=$variables{result_dir}.'/system_monitor.log';
	$variables{statistics_file}=$variables{result_dir}.'/statistics.txt';
    #
	return(\%variables);
}#1

###################################################################
sub Pre_sample_info{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	   
	printf( "raw data directory: %s\n", $variables{raw_data_dir}); 
	#incrusive files and directories
	my $files_pointer=sub_basic::files_list($variables{raw_data_dir}, 'incrusive_file');
	my @incrusive_files=@$files_pointer;
	my @raw_files=grep(/\.fastq$|\.fq$/i, @incrusive_files);
	my $raw_files_num=@raw_files;
	print "A total of $raw_files_num will be analyzed.\n\n";
	
	#The default is paired-end sequencing
	my @raw_files_input=@raw_files;
	if ($raw_files[0]=~/R1_|_R1/){
		@raw_files_input=grep(/_R1|R1_/, @raw_files) ;
	}
	$variables{sequencing_end}=1 unless $raw_files[0]=~/R2_|_R2/;
	#print "@raw_files_input\n";
	
	print "The file list of raw data:\n";
	#read sample_info file
	my $samples_pointer=sub_basic::file_to_hash($variables{sample_info_file}, ',');
	my %samples=%$samples_pointer;
	#
	my (%sample_info, %raw_to_sample);
	foreach my $sample_name(keys %samples){#2
		my @raw_names=split(",", $samples{$sample_name});
		
		my $files_size=0;
		my @sample_raw_files;
		foreach my $raw_name( @raw_names ) {#3
			foreach my $raw_file(@raw_files){
				my $raw_file_name=sub_basic::file_operation($raw_file, 'file_name');
				if ($raw_file_name=~/^$raw_name/i){
					push(@sample_raw_files, $raw_file);
					$files_size += sub_basic::file_operation($raw_file, 'file_size');
					$raw_to_sample{$raw_file}=$sample_name;
				} 
			}
		}#3
		$sample_info{$sample_name}->{'raw_names'}=$samples{$sample_name};
		$sample_info{$sample_name}->{'raw_files'}=join(',', @sample_raw_files);
		$sample_info{$sample_name}->{'parallel_parts'}=@sample_raw_files/$variables{sequencing_end};
		$sample_info{$sample_name}->{'files_size'}=$files_size;
		$sample_info{$sample_name}->{'supposed_time'}=int($files_size/228383); #unit is second
		print "$sample_name (parallel parts=$sample_info{$sample_name}->{parallel_parts}):\n";
		print "\t$sample_info{$sample_name}->{raw_files}\n";
	}#2
	
	$variables{raw_files_input_pointer}=\@raw_files_input;
	$variables{sample_info_pointer}=\%sample_info;
	$variables{raw_to_sample_pointer}=\%raw_to_sample;
	my @sample_names=sort( keys %sample_info);
	$variables{sample_names}=join(",", @sample_names);
	
	return(\%variables);
}


#######################################
#genome DNA sequences
sub Pre_chr_info{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	
	#generate genome_info_csv file
	unless(-f $variables{genome_info_csv}){#2
		print "generate $variables{genome_info_csv}.\n";
		open my($OUT), ">", $variables{genome_info_csv} or die;
		print $OUT join(',', 'chr', 'chr_len', 'GC_perc', 'enzyme_num' ), "\n";
		open my($IN), "<",  $variables{genome_fasta_file} or die;
		my($chr,$seq);
		while( defined ($chr=<$IN>) && defined ($seq=<$IN>) ){#3
			chomp($chr, $seq);
			$chr=~s/^>//;
			my $seq_len=length($seq);
			my $GC_num=map {/G|C/gi} $seq;
			my $GC_cont=int($GC_num*100/$seq_len+0.5);
			my $enzyme_num=map {/$variables{enzyme_site}/gi} $seq;
			print $OUT join(',', $chr, $seq_len, $GC_cont, $enzyme_num), "\n";
			printf("\t...%s...\n", $chr);
		}#
		close($IN);
		close($OUT);
	}#2
	
	#get chromosome informations 
	print "Get genome information.\n";
	my $chr_info_pointer=sub_basic::file_to_hash($variables{genome_info_csv}, ',');
	my %chr_info=%$chr_info_pointer;
	#
	return(\%chr_info);
}

#############
#size bin list of genome DNA sequences
sub Pre_size_bin{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	
	#generate the size_bin file
	unless(-f $variables{size_bin_csv}){#2
		print "generate the file $variables{size_bin_csv}.\n";
		open my($OUT), ">", $variables{size_bin_csv} or die;
		print $OUT join(',', 'bin_name', 'chr', 'bin_no', 'bin_start', 'bin_end', 'enzyme_num', 'GC_cont'), "\n";
		open my($IN), "<", $variables{genome_fasta_file} or die;
		my($chr, $chr_seq);
		while( defined ($chr=<$IN>) && defined ($chr_seq=<$IN>) ){#3
			chomp($chr, $chr_seq);
			$chr=~s/^>//;
			my $chr_len=length($chr_seq);
			
			#export size bin along a chromosome
			my $bin_start=1;
			my $bin_no=1;
			while($bin_start < $chr_len){#4
				my $sub_seq=substr($chr_seq, $bin_start-1, $variables{size_bin_len});
				my $chr_sites_str=sub_3C::enzyme_sites_counting($sub_seq, $variables{enzyme_site}, $bin_start);
				my $enzyme_num=split(',', $chr_sites_str);
				my $bin_end=$bin_start+$variables{size_bin_len}-1;
				$bin_end=$chr_len if $bin_end>$chr_len;
				my $GC_num=map {/G|C/gi} $sub_seq;
				my $GC_cont=int($GC_num*100/length($sub_seq)+0.5);
				print $OUT join(',', $chr.'_'.$bin_no, $chr, $bin_no, $bin_start, $bin_end, $enzyme_num, $GC_cont), "\n";
				#print join(',', $chr.'_'.$bin_no, $chr, $bin_no, $bin_start, $bin_end, $enzyme_num, $GC_cont), "\n";
				
				#update for incrusive
				$bin_start += $variables{size_bin_len};
				$bin_no++;
			}#4
			printf("\t...%s...\n", $chr);
		}#3
		close($IN);
		close($OUT);
	}#2
	
	#get size bin informations
	print "Get size bin information.\n";
	my %size_bins;
	open my($IN), "<", $variables{size_bin_csv} or die;
	my $header_line=<$IN>;
	while(<$IN>){
		chomp($_);
		my @items=split(",", $_);
		my $chr=$items[1];
		my $bin_no=$items[2];
		my $bin_start=$items[3];
		my $bin_end=$items[4];
		$size_bins{$chr}->{$bin_no}->{bin_start}=$bin_start;
		$size_bins{$chr}->{$bin_no}->{bin_end}=$bin_end;
	}
	close($IN);
	
	#
	return(\%size_bins);
}

###################
#
sub get_enzyme_sites_file{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
		
	unless (-f $variables{enzyme_sites_file}	){#2
		print "get sequences of genome DNA and search enzyme sites:\n";
		open my($ENZ), ">", $variables{enzyme_sites_file} or die;
		open my($IN), "<",  $variables{genome_fasta_file} or die;
		my($chr,$seq);
		while( defined ($chr=<$IN>) && defined ($seq=<$IN>) ){#3
			chomp($chr, $seq);
			$chr=~s/^>//;
			my $seq_len=length($seq);
			#export enzyme sites along a chromosome
			my @sites=split(",", sub_3C::enzyme_sites_counting($seq, $variables{enzyme_site}, 1) );
			my $sites_num=@sites;
			my $chr_sites_str=join(',', @sites);
			print $ENZ join("\t", $chr, $seq_len, $sites_num, $chr_sites_str), "\n";
			print "\tChr$chr: $seq_len bp......\n";
		}#4
		close($ENZ);
	}#2

}

#######################################
#enzyme list of genome DNA sequences
sub Pre_enzyme_bin{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	
	## enzyme sites counting along chromosomes 
	sub_3C::get_enzyme_sites_file(\%variables);
	
	#count enzyme sites along chromosomes
	print "export $variables{enzyme_bin_csv}\n";
	my %enzyme_bins;
	open my($EN_bin), ">", $variables{enzyme_bin_csv} or die;
	#column names in enzyme_bin_file
	print $EN_bin join(',', 'bin_name', 'chr', 'bin_no', 'bin_start', 'bin_end', 'enzyme_pos', 'bin_len', 'size_bin_name'), "\n";
	open my($EN), "<", $variables{enzyme_sites_file} or die;
	while(<$EN>){#2
		chomp($_);
		my($chr, $chr_len, $sites_num, $sites_str)=split("\t", $_);
		my @sites=split(',', $sites_str);
		unshift (@sites, 1) if $sites[0]>1; #add the 5 end
		push (@sites, $chr_len) if $sites[-1]<$chr_len; #add the 3 end
		
		#calcaulate start/end position of enzyme fragment 
		my(@sites_1, @sites_2, $right_site);
		my $first_site=shift @sites;
		push(@sites_2,  $first_site);
		foreach my $site(@sites){#3
			push(@sites_2, $site); #combine neigbouring enzyme site
			my $bin_len=$sites_2[-1] - $sites_2[0];
			if( (@sites_2>=$variables{enzyme_bin_enzymes} and $bin_len>=$variables{enzyme_bin_lower} ) 
				or $bin_len>=$variables{enzyme_bin_upper} or $site==$sites[-1]){
					push(@sites_1, join('_', @sites_2) );
					@sites_2=$site;#clear @sites_2 and give the first element
			}
		}#3
		
		#
		my $bin_no=0;
		foreach my $start_end(@sites_1){#3
			$bin_no++;
			my @site_start_end=split('_', $start_end);
			my $start_site=$site_start_end[0];
			my $end_site=$site_start_end[-1];
			my $site_name=$chr.'_'.$bin_no;
			my $bin_len=$end_site-$start_site;
			my $size_bin_name=$chr.'_'.(int($start_site/$variables{size_bin_len})+1);
			#pop @site_start_end; #remove the last
			print $EN_bin join(',', $site_name, $chr, $bin_no, $start_site, $end_site, $start_end, $bin_len, $size_bin_name), "\n";
			#export to %enzyme_bin;
			$enzyme_bins{$chr}->{$bin_no}->{bin_start}=$start_site;
			$enzyme_bins{$chr}->{$bin_no}->{bin_end}=$end_site;
		}#3
	}#2
	close($EN);
	close($EN_bin);
	
	return(\%enzyme_bins);
}
################################
#for localizing enzyme bin
sub Pre_enzyme_bin_regions{
	my ($enzyme_bins_pointer)=@_;
	my %enzyme_bins=%$enzyme_bins_pointer;
	print "chromosome regions of enzyme bins";
	
	my %enzyme_bin_regions;
	foreach my $chr(keys %enzyme_bins){#2
		my $chr_pointer=$enzyme_bins{$chr};
		my %chr_bins=%$chr_pointer;
		foreach my $bin_no(keys %chr_bins){#3
			my $region_start=int($chr_bins{$bin_no}->{bin_start}/1e6)-1;
			my $region_end=int($chr_bins{$bin_no}->{bin_end}/1e6)+1;
			for (my $region=$region_start; $region<=$region_end; $region++){
				$enzyme_bin_regions{$chr}->{$region}->{$bin_no}->{bin_start}=$chr_bins{$bin_no}->{bin_start};
				$enzyme_bin_regions{$chr}->{$region}->{$bin_no}->{bin_end}=$chr_bins{$bin_no}->{bin_end};
			}
		}#3
	}#2
	
	return(\%enzyme_bin_regions);
}
##############################################################
#
sub Pre_viewpoint_info_3CMTS{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $enzyme_bin_regions_pointer=$variables{enzyme_bin_regions_pointer};
	my %enzyme_bin_regions=%$enzyme_bin_regions_pointer;
	my $enzyme_len=length($variables{enzyme_site});

	unless (-f $variables{probe_seq_csv}){#2
		print "read sequences of genome DNA:\n";
		my $genome_pointer=sub_basic::genome_sequences($variables{genome_fasta_file});
		my %genome=%$genome_pointer;

		print "specific known seq for detection.\n";
		my $num=1;
		open my($IN), "<", $variables{site_info_file} or die;
		open my($OUT), ">", $variables{probe_seq_csv} or die;
		while (<$IN>){#3
			chomp($_);
			my ($site_name, $chr, $site_position, $boundary_start, $boundary_end, $upstream_len, $downstream_len)=split(',', $_);
			if ($upstream_len>0){#4
				my $display_id=$site_name.'u';
				my $region_start=$site_position-$upstream_len;
				my $region_end=$site_position-1+$enzyme_len;
				my $fragment=substr($genome{$chr}, $site_position-$upstream_len-1, $upstream_len+$enzyme_len ); #EcorI is on the end
				my $shortest_len=sub_3C::shortest_unique_seq(\%genome, $enzyme_len, $chr, $site_position, 'up');
				print $OUT join(',', $display_id, $chr, $site_position, $region_start, $region_end, 
								$boundary_start, $boundary_end, $shortest_len, $fragment), "\n";
				#print "$display_id=", substr($genome{$chr}, $site_position-1, 20), "\n\n";
			}#4

			if ($downstream_len>0){#4	
				my $display_id=$site_name.'d';
				my $region_start=$site_position;
				my $region_end=$site_position+$downstream_len-1;
				my $revcom_fragment=substr($genome{$chr}, $site_position-1, $downstream_len+$enzyme_len ); #EcorI is on the head
				my $fragment=sub_basic::reverse_complement($revcom_fragment);#EcorI is on the end
				my $shortest_len=sub_3C::shortest_unique_seq(\%genome, $enzyme_len, $chr, $site_position, 'down');
				print $OUT join(',', $display_id, $chr, $site_position, $region_start, $region_end,
								$boundary_start, $boundary_end, $shortest_len, $fragment), "\n";
				#print "$display_id=", substr($genome{$chr}, $site_position-1, 20), "\n\n";
			}#4
			print "$site_name\n";
			$num++;
		}#3
		close($IN);
		close($OUT);
	}#2
	
	print "get site information. \n\n";
	my %viewpoint_info;
	open my($IN), "<", $variables{probe_seq_csv} or die;
	while(<$IN>){#2
		chomp($_);
		my @array=split(',', $_);
		my $display_id=$array[0];
		my $chr=$array[1];
		my $site_pos=$array[2];
		$viewpoint_info{$display_id}->{chr}=$chr;#chr
		$viewpoint_info{$display_id}->{site_position}=$site_pos;
		$viewpoint_info{$display_id}->{region_start}=$array[3];
		$viewpoint_info{$display_id}->{region_end}=$array[4];
		$viewpoint_info{$display_id}->{boundary_start}=$array[5];
		$viewpoint_info{$display_id}->{boundary_end}=$array[6];
		$viewpoint_info{$display_id}->{shortest_len}=$array[7]; 
		$viewpoint_info{$display_id}->{fragment}=$array[8];
		$viewpoint_info{$display_id}->{shortest_seq}=substr($viewpoint_info{$display_id}->{fragment}, -($viewpoint_info{$display_id}->{shortest_len}+$enzyme_len) );
		$viewpoint_info{$display_id}->{revcom_fragment}=sub_basic::reverse_complement( $viewpoint_info{$display_id}->{fragment} );
		#site region localization
		$viewpoint_info{$display_id}->{size_bin_name}=int($site_pos/$variables{size_bin_len})+1;
		my $corr_site_pos=($display_id=~/u$/) ? ($site_pos-1): ($site_pos+1);
		$viewpoint_info{$display_id}->{enzyme_bin_name}=sub_3C::enzyme_bin_localization($variables{enzyme_bin_regions_pointer}, $chr, $corr_site_pos);
	}#2
	close($IN);
	
	unless (-f $variables{viewpoints_csv}){
		print "export viewpoints information\n";
		open my($OUT), ">", $variables{viewpoints_csv} or die;
		print $OUT join(',', 'viewpoints', 'sample', 'display_id', 'chr', 'site_pos', 'region_start', 'region_end', 'size_bin_name', 'enzyme_bin_name'), "\n";
		my @sample_names=split(",", $variables{sample_names});
		foreach my $sample(@sample_names){
			foreach my $display_id(sort (keys %viewpoint_info) ){
				my $viewpoint=$sample.'_'.$display_id;
				print $OUT join(',', $viewpoint, $sample, $display_id, 
							$viewpoint_info{$display_id}->{chr}, $viewpoint_info{$display_id}->{site_position}, 
							$viewpoint_info{$display_id}->{region_start}, $viewpoint_info{$display_id}->{region_end}, 
							$viewpoint_info{$display_id}->{size_bin_name}, $viewpoint_info{$display_id}->{enzyme_bin_name}, ), "\n";
			}
		}
		close($OUT);
	}
	
	return(\%viewpoint_info);
}
#######################################################################
#
sub Pre_viewpoint_info_4C{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $enzyme_bin_regions_pointer=$variables{enzyme_bin_regions_pointer};
	my %enzyme_bin_regions=%$enzyme_bin_regions_pointer;
	my $enzyme_len=length($variables{enzyme_site});
	
	#read genome information
	unless(-f $variables{probe_seq_csv}){
		print "read sequences of genome DNA:\n";
		my $genome_pointer=sub_basic::genome_sequences($variables{genome_fasta_file});
		my %genome=%$genome_pointer;
		
		print "specific known seq for detection.\n";
		my $num=1;
		open my($IN), "<", $variables{site_info_file} or die;
		open my($OUT), ">", $variables{probe_seq_csv} or die;
		while (<$IN>){#3
			chomp($_);
			my($sample_name, $chr, $site_pos, $up_down, $view_len)=split(',', $_);
			
			my $size_bin_name=int($site_pos/$variables{size_bin_len})+1;
			my $enzyme_bin_name=sub_3C::enzyme_bin_localization($variables{enzyme_bin_regions_pointer}, $chr, $site_pos);
			#get sequence around viewpoints
			my $chr_seq=$genome{$chr};
			my $start=$site_pos+1-$view_len+length($variables{enzyme_site});
			my $end=$site_pos-1+length($variables{enzyme_site});
			my $view_seq=substr($chr_seq, $site_pos-1-$view_len+length($variables{enzyme_site}), $view_len);
			my $region_seq=substr($chr_seq, $site_pos-1-$view_len, $view_len*2);
			if($up_down eq 'Down'){
				$start=$site_pos;
				$end=$site_pos-1+$view_len;
				$view_seq=sub_basic::reverse_complement(substr($chr_seq, $site_pos-1, $view_len));
				$region_seq=sub_basic::reverse_complement($region_seq);
			}
			print $OUT join(',', $sample_name, $chr, $site_pos, $size_bin_name, $enzyme_bin_name,
										$start, $end, $view_seq, $region_seq), "\n";
		}
		close($IN);
		close($OUT);
	}
	
	print "get site information. \n\n";
	my %viewpoint_info;
	open my($IN), "<", $variables{probe_seq_csv} or die;
	while(<$IN>){#2
		chomp($_);
		my ($sample_name, $chr, $site_pos, $size_bin_name, $enzyme_bin_name, $start, $end, $view_seq, $region_seq)=split(',', $_);
		$viewpoint_info{$sample_name}->{site_name}=$sample_name;#chr
		$viewpoint_info{$sample_name}->{chr}=$chr;#chr
		$viewpoint_info{$sample_name}->{site_position}=$site_pos;
		$viewpoint_info{$sample_name}->{view_seq}=$view_seq;
		$viewpoint_info{$sample_name}->{region_seq}=$region_seq;
		$viewpoint_info{$sample_name}->{region_start}=$start;
		$viewpoint_info{$sample_name}->{region_end}=$end;
		$viewpoint_info{$sample_name}->{size_bin_name}=$size_bin_name;
		$viewpoint_info{$sample_name}->{enzyme_bin_name}=$enzyme_bin_name;
	}#2
	close($IN);
	
	unless (-f $variables{viewpoints_csv}){
		print "export viewpoints information\n";
		open my($OUT), ">", $variables{viewpoints_csv} or die;
		print $OUT join(',', 'viewpoints', 'chr', 'site_pos', 'region_start', 'region_end', 'size_bin_name', 'enzyme_bin_name'), "\n";
		my @sample_names=split(",", $variables{sample_names});
		foreach my $display_id(sort (keys %viewpoint_info) ){
			#$sample_name equal to $display_id;
			print $OUT join(',', $display_id, $viewpoint_info{$display_id}->{chr}, $viewpoint_info{$display_id}->{site_position}, 
							$viewpoint_info{$display_id}->{region_start}, $viewpoint_info{$display_id}->{region_end}, 
							$viewpoint_info{$display_id}->{size_bin_name}, $viewpoint_info{$display_id}->{enzyme_bin_name}, ), "\n";

		}
		close($OUT);
	}
	
	return(\%viewpoint_info);
}


############################################################
#calculate enzyme sites within a DNA sequence
sub enzyme_sites_counting{
	my($seq, $enzyme_site, $offset)=@_;

	my @chr_sites;
	my $pos=-1;
	do{
		$pos=index($seq, $enzyme_site, $pos+1);
		push(@chr_sites, $offset+$pos) unless $pos==-1;
	} until ($pos==-1);
	my $num=@chr_sites;

	my $out=join(',', @chr_sites);
	return($out);
}

#########################
#find the nearest enzyme site to offset of mappable read sequence against  chromosome
sub enzyme_offset{
	my ($chr, $sites_pointer, $offset)=@_;
	my %hash=%$sites_pointer;
	my $chr_len=$hash{chr_len};
	my $sites_num=$hash{sites_num};
	my @enzyme_sites=split(',', $hash{sites});
	
	my (%out, $left_offset, $right_offset, $fragment_len, $digestion_len);
	my @low_sites=grep {$_ < $offset} @enzyme_sites;
	if (@low_sites==0){
		$left_offset=1;
		$fragment_len=$enzyme_sites[0];
		$digestion_len=$offset;
	}
	elsif(@low_sites==@enzyme_sites){
		$left_offset=$enzyme_sites[-1];
		$fragment_len=$chr_len-$left_offset;
		$digestion_len=$offset-$left_offset+1;
	}
	else{
		$left_offset=$low_sites[-1];
		my $low_sites_num=@low_sites;
		$right_offset=$enzyme_sites[$low_sites_num];
		$fragment_len=$right_offset-$left_offset;
		$digestion_len=$offset-$left_offset+1;
	}
	$out{left_digestion_len}=$digestion_len;
	$out{left_viewpoint_info}=$chr.'_'.$left_offset.'_'.$fragment_len;
	
	return(\%out);
}
##########################
sub shortest_unique_seq{
	my($genome_seq_pointer, $enzyme_len, $chr, $site_pos, $direction)=@_;
	my %genome_seq=%$genome_seq_pointer;
	my $site_chr_seq=$genome_seq{$chr};

	my $len=15;
	#the identical chromosome with site_chr
	while($len<100){#3
		my $start=($direction eq 'up') ? ($site_pos-1-$len) : ($site_pos-1+$enzyme_len);
		my $seq=substr($site_chr_seq, $start, $len);
		my $index_1=index(substr($site_chr_seq, 0, $start+$enzyme_len), $seq);
		my $index_2=index(substr($site_chr_seq, $start+$enzyme_len), $seq);
		if ( $index_1>0 or $index_2>0){
			$len++;
		}
		else{
			last;
		}
	}#3
	delete $genome_seq{$chr};
	
	#the other chromosomes
	while( my($genome_chr, $parent_seq)=each (%genome_seq) ){#2
		while($len<100){#3
			my $start=($direction eq 'up') ? ($site_pos-1-$len) : ($site_pos-1+$enzyme_len);
			my $seq=substr($site_chr_seq, $start, $len);
			if ($parent_seq=~/$seq/){
				$len++;
			}
			else{
				last;
			}
		}#3
	}#2
	
	return($len);
}

################################################
#
sub hybrid_seq_assembly{
	my($enzyme_site, $seq_1, $seq_2, $share_len)=@_;
	my $enzyme_len=length($enzyme_site);
	my $len_1=length($seq_1);
	my $len_2=length($seq_2);
	my $max_len=($len_1<$len_2) ? $len_1 : $len_2;
	my $revcom_seq_2=sub_basic::reverse_complement($seq_2);
	#print "$seq_1\n$seq_2\n$revcom_seq_2\n";
	
	my %out=(ass_seq=>'NA', share_seq=>'NA', share_len=>0, seq_1r=>'NA', seq_2r=>'NA', hybrid=>'no');
	if($max_len>$share_len){#2
		foreach my $len($share_len..$max_len){#3
			my $seq_1_tail=substr($seq_1, -$len);
			if($revcom_seq_2=~/^$seq_1_tail/){
				$out{share_len}=$len;
				$out{share_seq}=$seq_1_tail;
			}
		}#3
		if($out{share_len}>0){
			$out{ass_seq}=$seq_1.substr($revcom_seq_2, -($len_2-$out{share_len}) );
		}
		my $enzyme_index=index($out{share_seq}, $enzyme_site);
		if($enzyme_index>=0){
			$out{seq_1r}=substr($seq_1, 0, $len_1-$out{share_len}+$enzyme_index+$enzyme_len);
			$out{seq_2r}=substr($revcom_seq_2, -($len_2-$enzyme_index) );
			$out{seq_2r}=sub_basic::reverse_complement($out{seq_2r});
			$out{hybrid}='yes';
		}
	}#2
	
	return(\%out);
}


##################################################
#trim adatper sequences and poor quality sequences from read seq in FASTQ format
sub S1_trim_fastq{
	my ($variables_pointer, $original_fastq_file, $trimmed_fastq_file, $min_adapter_seq)=@_;
	my %variables=%$variables_pointer;
			 
	open my($IN), "<", $original_fastq_file or die;
	open my($OUT), ">", $trimmed_fastq_file or die;
	my ($name_line, $seq_line, $third_line, $qual_line);
	while (defined($name_line = <$IN>) && defined($seq_line = <$IN>) && defined($third_line = <$IN>) && defined($qual_line = <$IN>) ){#2
		chomp($name_line, $seq_line, $third_line, $qual_line );
		#trim adapter sequences
		my $trim_len=$variables{min_trim_len};
		my $adapter_pos=index($seq_line, $min_adapter_seq);
		if($adapter_pos==-1){
			$trim_len=length($seq_line);
		}
		elsif ($adapter_pos>$variables{min_trim_len}){	#detect adapter and beyond min trim length
			$trim_len=$adapter_pos;
		}
		#export
		my $trimmed_seq=substr($seq_line, 0, $trim_len);
		my $trimmed_quality=substr($qual_line, 0, $trim_len);
		print $OUT join("\n", $name_line, $trimmed_seq, $third_line, $trimmed_quality), "\n";
	}#2
	close($IN);
	close($OUT);
	
}

####################################################
#
sub S2_Colocalization_3CMTS{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $R1_fastq=$variables{sample_out_file}.'.fastq_trimmed';
	my $R2_fastq=$R1_fastq;
	$R2_fastq=~s/_R1_/_R2_/;
	
	#sequences alignment
	my $sam_file=$variables{sample_out_file}.".sam";
	unless(-e $sam_file){
		print "$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $sam_file \n";
		system("$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $sam_file ");
	}
	
	print "\tread $sam_file\n";
	my %num;
	open my($ORI), ">", $variables{sample_out_file}.'.original_regions' or die;
	open my($cis_SITE), ">", $variables{sample_out_file}.'.cis_ligation' or die;
	open my($trans_SITE), ">", $variables{sample_out_file}.'.trans_ligation' or die;
	open my($OTH), ">", $variables{sample_out_file}.'.other_alignment' or die;
	open my($MULTI), ">", $variables{sample_out_file}.'.multiple_alignment' or die;
	open my($UN), ">", $variables{sample_out_file}.'.unknown_mapped_pairs' or die;
	open my($IN), "<", $sam_file or die;
	my ($line_1, $line_2);
	while( defined ($line_1 = <$IN>) && defined ($line_2 = <$IN>) ){ #2 read SAM
		chomp($line_1, $line_2);
		#read two lines per time
		my $single1_pointer=sub_3C::read_single_end(\%variables, $line_1);
		$single1_pointer=sub_3C::localize_viewpoint(\%variables, $single1_pointer);
		my %single1=%$single1_pointer;
		my $single2_pointer=sub_3C::read_single_end(\%variables, $line_2);
		$single2_pointer=sub_3C::localize_viewpoint(\%variables, $single2_pointer);
		my %single2=%$single2_pointer;
		#print "$single1{vp_displayid}:$single2{vp_displayid}\n";
		
		#get %pair
		my $pair_pointer=sub_3C::read_paired_end(\%variables, \%single1, \%single2);#subroutine
		my %pair=%$pair_pointer;
		#print "#$pair{valid_pair}\n";

		#judging
		if(exists $pair{localize_vp}){#3 #for viewpoint localization
			if($pair{valid_pair} eq 'cis_pair'){#
				print $ORI $pair{output}, "\n";
				$num{original_regions_num}++;
			}#
			elsif($pair{valid_pair} eq 'cis_unpair'){#
				print $cis_SITE $pair{output}, "\n";
				$num{cis_ligation_num}++;
				$num{cis_ligation_both_num}++;
			}#
			elsif($pair{valid_pair} eq 'trans'){#
				print $trans_SITE $pair{output}, "\n";
				$num{trans_ligation_num}++;
				$num{trans_ligation_both_num}++;
			}#
			$num{both_alignment_num}++;
		}#3
		else{#3 
			if ($pair{valid_pair} eq 'one_unique'){#4 one seq aligned
				print $OTH $pair{output}, "\n";
				$num{one_alignment_num}++;
			}#4
			elsif ($pair{valid_pair} eq 'unaligned'){#4 #unaligned
				print $OTH $pair{output}, "\n";
				$num{no_alignment_num}++;
			}#4
			elsif ($pair{valid_pair} eq 'multiple'){	#4
				print $MULTI $pair{output}, "\n";
				$num{multiple_alignment_num}++;
			}#4
			else{#4
				print $UN $pair{output}, "\n";
				$num{unknown_alignment_num}++;
				$num{both_alignment_num}++;
			}#4
		}#3
		$num{raw_reads_num}++;
	} #2 read SAM
	close($IN);
	close($ORI);
	close($cis_SITE);
	close($trans_SITE);
	close($OTH);
	close($UN);
	close($MULTI);
	
	#refresh sample log
	while(my($key, $value)=each(%num) ){
		sub_basic::refresh_log($variables{sample_log}, $key, $value);
	}

}

###############################################
#
sub S3_detect_other_alignment_3CMTS{
	my ($variables_pointer, $pairs_pointer, $alignment_type)=@_;
	my %variables=%$variables_pointer;
	my %pairs=%$pairs_pointer;
	
	#one seq alignment
	my $R1_pointer=sub_3C::one_seq_alignment(\%variables, $pairs{R1_pointer});
	my %R1=%$R1_pointer;
	my $R2_pointer=sub_3C::one_seq_alignment(\%variables, $pairs{R2_pointer});
	my %R2=%$R2_pointer;
	
	#export
	my %num;
	open my($ORI), ">>", $variables{sample_out_file}.'.original_regions' or die;
	open my($cis_SITE), ">>", $variables{sample_out_file}.'.cis_ligation' or die; # $found=2
	open my($trans_SITE), ">>", $variables{sample_out_file}.'.trans_ligation' or die; #$found=4
	foreach my $coordinate(keys %R1){#2
		my $single1_pointer=$R1{$coordinate};
		$single1_pointer=sub_3C::localize_viewpoint(\%variables, $single1_pointer);
		my $single2_pointer=$R2{$coordinate};
		$single2_pointer=sub_3C::localize_viewpoint(\%variables, $single2_pointer);
		#paired end combination
		my $pair_pointer=sub_3C::read_paired_end(\%variables, $single1_pointer, $single2_pointer);
		my %pair=%$pair_pointer;
		
		if(exists $pair{localize_vp}){#3
			if ($pair{valid_pair} eq 'cis_unpair'){	
				print $cis_SITE "$pair{output}\n";
				$num{cis_ligation_num}++;
			}
			elsif ($pair{valid_pair} eq 'cis_pair'){	
				print $ORI "$pair{output}\n";
				$num{original_regions_num}++;
			}
			elsif ($pair{valid_pair} eq 'trans'){	
				print $trans_SITE "$pair{output}\n";	
				$num{trans_ligation_num}++;
			}
		}#3
	}#2
	close($ORI);
	close($cis_SITE);
	close($trans_SITE);
	
	#refresh sample log
	while(my($key, $value)=each(%num) ){
		sub_basic::refresh_log($variables{sample_log}, $key, $value, 'add');
	}

}


###################################################
#
sub read_single_end{
	my ($variables_pointer, $sam_line)=@_;
	my %variables=%$variables_pointer;
	my @items=split("\t", $sam_line);
	
	my %single;
	#
	my @query_info=split(":", $items[0]);
	$single{coordinate}=(exists $query_info[-2]) ? $query_info[-2].":".$query_info[-1] : $query_info[-1];
	$single{len}=length($items[9]);
	$single{MapQ}=$items[4]; 
	$single{mate_ref}=$items[6]; 
	$single{mate_offset}=$items[7];
	#flag
	$single{flag}=sprintf("%08b", $items[1]); #convert to binary number
	if(substr($single{flag}, 5, 1)==1){#no reported alignment (+4)
		$single{ref}='*';
		$single{start}=0; 
		$single{start_sbn}=0; 
		$single{start_ebn}=0; 
		$single{end}=0;
		$single{end_sbn}=0;
		$single{end_ebn}=0;
	}
	else{
		$single{ref}=$items[2];
		$single{start}=$items[3]; 
		$single{start_sbn}=int($single{start}/$variables{size_bin_len})+1;
		$single{start_ebn}=enzyme_bin_localization($variables{enzyme_bin_regions_pointer}, $single{ref}, $single{start});
		$single{end}=$single{start}+$single{len};
		$single{end_sbn}=int($single{end}/$variables{size_bin_len})+1;
		$single{end_ebn}=enzyme_bin_localization($variables{enzyme_bin_regions_pointer}, $single{ref}, $single{end});
	}
	if(substr($single{flag}, 3, 1)==1){	#reversed reference strand	
		$single{seq}=sub_basic::reverse_complement($items[9]);
		$single{alignment}='revcom';
	}
	else{
		$single{seq}=$items[9];
		$single{alignment}='itself';
	}
	#option fields
	for (my $i=11; $i<@items; $i++){
		my($a, $b, $c)=split(':', $items[$i]);
		$single{$a}=$c;
	}
	#unique, multiple, no alignment
	if (exists $single{AS}){
		if (exists $single{XS} and $single{AS}<=$single{XS}){
			$single{valid}= 4;# multiple alignment
		}
		else{
			$single{valid}=2;#unique alignment
		}
	}
	else{	
		$single{valid}=1;#no alignment
	} 
	
	return(\%single);
}
####################################################
#all sequences from raw data without transfering revcom chain
#
#Note seq_1 and seq_2 is equal to R1_seq and R2_seq,
# however might be reversed-complementary in order to matching reference seq.
sub read_paired_end{
	my ($variables_pointer, $single1_pointer, $single2_pointer)=@_;
	my %variables=%$variables_pointer;
	my %s1=%$single1_pointer;
	my %s2=%$single2_pointer;
	
	#R1_seq is alway the upstream, and R2_seq is always the downstream
	if($s1{start} and $s2{start} and ($s1{start} > $s2{start}) ){
		my %tmp=%s1;
		%s1=%s2;
		%s2=%tmp;
	}

	my %pair;
	$pair{coordinate}=$s1{coordinate};
	#default output line
	$pair{output}=join("\t", $pair{coordinate}, 
			$s1{ref}, $s1{start}, $s1{end}, $s1{start_sbn}, $s1{end_sbn}, 
			$s1{start_ebn}, $s1{end_ebn}, $s1{seq}, $s1{valid},
			$s2{ref}, $s2{start}, $s2{end}, $s2{start_sbn}, $s2{end_sbn}, 
			$s2{start_ebn}, $s2{end_ebn}, $s2{seq}, $s2{valid},	);
	#valid_pair, coordinate, output
	#unaligned=1, unique=2, multiple=4
	#all possible: 2,3,4,5,6,8
	my $judging=$s1{valid}+$s2{valid};
	if ($judging==4){#2
		if ($s1{ref} eq $s2{ref}){
			$pair{valid_pair}='cis_unpair';
			$pair{valid_pair}='cis_pair' if ($s2{end}-$s1{start}) < $variables{colocalization_range};
		}
		else{
			$pair{valid_pair}= 'trans';
		}
		# the ouput line with viewpoint localization
		if(exists $s1{vp_displayid}){
			$pair{output}=join("\t", $pair{coordinate}, 
					$s1{vp_displayid}, $s1{vp_chr}, $s1{vp_pos}, $s1{vp_sbn}, $s1{vp_ebn}, 
					$s2{ref}, $s2{start}, $s2{end},$s2{start_sbn}, $s2{end_sbn},$s2{start_ebn}, 
					$s2{end_ebn}, $s2{seq}, $s2{valid}, 	);
			$pair{localize_vp}='R1';
		}
		elsif(exists $s2{vp_displayid}){
			$pair{output}=join("\t", $pair{coordinate}, 
					$s2{vp_displayid}, $s2{vp_chr}, $s2{vp_pos}, $s2{vp_sbn}, $s2{vp_ebn}, 
					$s1{ref}, $s1{start}, $s1{end},$s1{start_sbn}, $s1{end_sbn},$s1{start_ebn}, 
					$s1{end_ebn}, $s1{seq}, $s1{valid},	);
			$pair{localize_vp}='R2';
		}
	}#2
	elsif ($judging==3){#2 one sequence aligned
		$pair{valid_pair}='one_unique';
	}#2
	elsif ($judging==5 or $judging==6 or $judging==8){#2 multiple
		$pair{valid_pair}='multiple';
	}#2
	else{#2  both no alignment
		$pair{valid_pair}='unaligned';
	}#2
	
	return(\%pair);
}

######################3
#
sub read_pairs_output{
	my ($variables_pointer,$alignment_file)=@_;
	my %variables=%$variables_pointer;
	
	my (%pairs, %R1, %R2);
	open my($IN), "<", $alignment_file or die;
	while(<$IN>){
		chomp($_);
		my @items=split("\t", $_);
		my $coordinate=$items[0];
		my $len=@items;
		#
		print "$_\n" unless $len==19;
		$R1{$coordinate}->{coordinate}=$coordinate;
		$R1{$coordinate}->{ref}=$items[1];
		$R1{$coordinate}->{start}=$items[2];
		$R1{$coordinate}->{end}=$items[3];
		$R1{$coordinate}->{start_sbn}=$items[4];
		$R1{$coordinate}->{end_sbn}=$items[5];
		$R1{$coordinate}->{start_ebn} =$items[6];
		$R1{$coordinate}->{end_ebn} =$items[7];
		$R1{$coordinate}->{seq} = $items[8];
		$R1{$coordinate}->{trunc_seq}=sub_basic::split_seq($R1{$coordinate}->{seq}, $variables{enzyme_site});
		$R1{$coordinate}->{valid} =$items[9];
		$R2{$coordinate}->{coordinate}=$coordinate;
		$R2{$coordinate}->{ref}=$items[10];
		$R2{$coordinate}->{start}=$items[11];
		$R2{$coordinate}->{end}=$items[12];
		$R2{$coordinate}->{start_sbn}=$items[13];
		$R2{$coordinate}->{end_sbn}=$items[14];
		$R2{$coordinate}->{start_ebn} =$items[15];
		$R2{$coordinate}->{end_ebn} =$items[16];
		$R2{$coordinate}->{seq} =$items[17];
		$R2{$coordinate}->{trunc_seq}=sub_basic::split_seq($R2{$coordinate}->{seq}, $variables{enzyme_site});
		$R2{$coordinate}->{valid} =$items[18];
		
	}
	close($IN);
	
	#
	$pairs{R1_pointer}=\%R1;
	$pairs{R2_pointer}=\%R2;
	return(\%pairs);
}

######################3
#
sub read_unknown_output{
	my ($alignment_file, $sep)=@_;
	
	my (%pairs, %R1, %R2);
	open my($IN), "<", $alignment_file or die;
	while(my $line=<$IN>){#2 *.no_alignment circling
		chomp($line);
		my($coordinate, $seq_1, $seq_2)=split("\t", $line);
		#
		#print "$sep, $seq_1, $seq_2\n\n";
		my $out_pointer=hybrid_seq_assembly($sep, $seq_1, $seq_2, 15);
		my %out=%$out_pointer;
		if($out{hybrid} eq 'yes'){#3
			$R1{$coordinate}->{coordinate}=$coordinate;
			$R1{$coordinate}->{ref}='*';
			$R1{$coordinate}->{seq}=$seq_1;
			$R1{$coordinate}->{trunc_seq}=$out{seq_1r};
			$R2{$coordinate}->{coordinate}=$coordinate;
			$R2{$coordinate}->{ref}='*';
			$R2{$coordinate}->{seq}=$seq_2;
			$R2{$coordinate}->{trunc_seq}=$out{seq_2r};
		}#3
	}#2  *.no_alignment circling
	close($IN);
	#
	$pairs{R1_pointer}=\%R1;
	$pairs{R2_pointer}=\%R2;
	return(\%pairs);
}

###############################################
#
sub S3_detect_other_alignment_HiC{
	my ($variables_pointer, $pairs_pointer, $alignment_type)=@_;
	my %variables=%$variables_pointer;
	my %pairs=%$pairs_pointer;
	
	#one seq alignment
	my $R1_pointer=sub_3C::one_seq_alignment(\%variables, $pairs{R1_pointer});
	my %R1=%$R1_pointer;
	my $R2_pointer=sub_3C::one_seq_alignment(\%variables, $pairs{R2_pointer});
	my %R2=%$R2_pointer;
	
	#export
	my %num;
	open my($ORI), ">>", $variables{sample_out_file}.'.original_regions' or die;
	open my($cis_SITE), ">>", $variables{sample_out_file}.'.cis_ligation' or die; # $found=2
	open my($trans_SITE), ">>", $variables{sample_out_file}.'.trans_ligation' or die; #$found=4
	foreach my $coordinate(keys %R1){#2
		my $single1_pointer=$R1{$coordinate};
		my $single2_pointer=$R2{$coordinate};
		my $pair_pointer=sub_3C::read_paired_end(\%variables, $single1_pointer, $single2_pointer);
		my %pair=%$pair_pointer;
		
		my $output=join("\t", $pair{output}, $alignment_type );
		if ($pair{valid_pair} eq 'cis_unpair')		{	
			print $cis_SITE "$output\n";
			$num{cis_ligation_num}++;
		}
		if ($pair{valid_pair} eq 'cis_pair')		{	
			print $ORI "$output\n";
			$num{original_regions_num}++;
		}
		elsif ($pair{valid_pair} eq 'trans')	{	
			print $trans_SITE "$output\n";	
			$num{trans_ligation_num}++;
		}
	}#2
	close($ORI);
	close($cis_SITE);
	close($trans_SITE);
	
	#refresh sample log
	while(my($key, $value)=each(%num) ){
		sub_basic::refresh_log($variables{sample_log}, $key, $value, 'add');
	}

}

#######################################
sub one_seq_alignment{
	my($variables_pointer, $one_end_pointer)=@_;
	my %variables=%$variables_pointer;
	my %one_end=%$one_end_pointer;
	
	#export query fasta file
	my $fasta_file=$variables{sample_out_file}.'_one_alignment.fa';
	open my($OUT), ">", $fasta_file or die;
	foreach my $coordinate(keys %one_end){#2
		if ($one_end{$coordinate}->{ref} eq '*'){#3
			unless($one_end{$coordinate}->{seq} eq $one_end{$coordinate}->{trunc_seq}){
				print $OUT ">$coordinate\n", "$one_end{$coordinate}->{trunc_seq}\n";
			}
		}#3
	}#2
	close($OUT);
	
	#make alignment
	my $sam_file=$variables{sample_out_file}.'_one_alignment.sam';
	system("$variables{bowtie_options} -f $fasta_file -S $sam_file");
	
	#read sam file
	open my($IN), "<", $sam_file or die;
	while(my $line=<$IN>){#2
		chomp($line);
		my $single_pointer=read_single_end(\%variables, $line);
		my %single=%$single_pointer;
		
		#unmapped alignment: 1, unqiue alignment is 2, multiple alignmet is 4
		my $coordinate=$single{coordinate};
		$one_end{$coordinate}->{ref}=$single{ref};
		$one_end{$coordinate}->{start}=$single{start};
		$one_end{$coordinate}->{end}=$single{end};
		$one_end{$coordinate}->{start_sbn}=$single{start_sbn};
		$one_end{$coordinate}->{end_sbn}=$single{end_sbn};
		$one_end{$coordinate}->{start_ebn}=$single{start_ebn};
		$one_end{$coordinate}->{end_ebn}=$single{end_ebn};
		$one_end{$coordinate}->{valid}=$single{valid};
	}#2
	
	return(\%one_end);
}

#################################
#
sub localize_viewpoint{
	my($variables_pointer, $single_pointer)=@_;
	my %variables=%$variables_pointer;
	my %single=%$single_pointer;
	my $viewpoint_info_pointer=$variables{viewpoint_info_pointer};
	my %viewpoint_info=%$viewpoint_info_pointer;
	
	#print "#$single{ref}#\n";
	unless($single{ref} eq '*'){#2
		foreach my $display_id(keys %viewpoint_info){#3
			my $chr=$viewpoint_info{$display_id}->{chr};#
			my $r_start=$viewpoint_info{$display_id}->{region_start};#
			my $r_end=$viewpoint_info{$display_id}->{region_end};#
			if($chr eq $single{ref} and $single{start} > $r_start and $single{start} < $r_end){
				$single{vp_displayid}=$display_id;
				$single{vp_chr}=$chr;
				$single{vp_pos}=$viewpoint_info{$display_id}->{site_position};#
				$single{vp_sbn}=$viewpoint_info{$display_id}->{size_bin_name};#
				$single{vp_ebn}=$viewpoint_info{$display_id}->{enzyme_bin_name};#
				last;
			}
		}#3
	}#2
	
	return(\%single);
}
##########################################
sub enzyme_bin_localization{
	my($enzyme_bin_regions_pointer, $chr, $offset)=@_;
	my %enzyme_bin_regions=%$enzyme_bin_regions_pointer;
	
	my $enzyme_bin_no=0;
	my $region=int($offset/1e6);
	if(exists $enzyme_bin_regions{$chr}->{$region}){#2
		my $chr_sites_pointer=$enzyme_bin_regions{$chr}->{$region};
		my %chr_sites=%$chr_sites_pointer;
		foreach my $bin_no(keys %chr_sites){#3
			if ($offset>=$chr_sites{$bin_no}->{bin_start} and $offset<=$chr_sites{$bin_no}->{bin_end}){
				$enzyme_bin_no=$bin_no;
				last;
			}
		}#3
	}#2
	#print "fail to localizing enzyme bin ----$chr:$region:$offset\n" if $enzyme_bin_no==0;

	return($enzyme_bin_no);
}


##############################################################################
#export: statistics.txt 
sub S4_result_statistics{#1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});
	
	print "\texport ligations into ligation_dir\n";
	my @combined_file_tails=('unmapped_ligation', 'cis_ligation', 'trans_ligation', 'multiple_ligation', 'original_regions',);
	foreach my $sample_name(@sample_names){#2
		my $sample_dir=$variables{result_dir}.'/'.$sample_name;
		my $out_file_head=$variables{ligations_dir}.'/'.$sample_name;
		foreach my $file_tail(@combined_file_tails){#3
			sub_basic::combine_files($sample_dir, $file_tail, $out_file_head);
		}
	}#2
	
	print "\tcombine sample log files\n";
	sub_basic::combine_log_files(\%variables);
	
	#read sample log
	my %counting_hash;
	foreach my $sample_name(@sample_names){#2
		my $sample_log=$variables{sample_log_dir}.'/'.$sample_name.'.log';
		open my($IN), "<", $sample_log or die "Error: canot open $sample_log!\n";
		while(<$IN>){#3
			chomp($_);
			my($key, $value)=split("=", $_);
			$counting_hash{$key}->{$sample_name}=$value;
		}#3
		close($IN);
		$counting_hash{raw_files}->{$sample_name}=split(',', $sample_info{$sample_name}->{'raw_files'});
	}#2
	
	#generate statistics result 
	sub_basic::hash2_to_file(\%counting_hash, $variables{statistics_file});
	#
}


################################################
sub S6_Tscore_calculation{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
		print "###Tscore calculation.###\n";
		sub_basic::refresh_R_script($variables{Tscore_R_file}, "R_ligation_frequency_dir" , $variables{ligation_frequency_dir});
		sub_basic::refresh_R_script($variables{Tscore_R_file}, "R_viewpoints_csv" , $variables{viewpoints_csv});
		sub_basic::refresh_R_script($variables{Tscore_R_file}, "R_size_bin_csv" , $variables{size_bin_csv});
		sub_basic::refresh_R_script($variables{Tscore_R_file}, "R_enzyme_bin_csv" , $variables{enzyme_bin_csv});
		sub_basic::refresh_R_script($variables{Tscore_R_file}, "R_p_level" , $variables{RC_cumulative_p});
		sub_basic::refresh_R_script($variables{Tscore_R_file}, "R_RC_noise" , $variables{RC_noise});
		
		my $file_names_pointer=sub_basic::files_list($variables{ligation_frequency_dir}, 'file_name');
		my @file_names=@$file_names_pointer;
		my @RC_csv_names=grep(/_RC.csv$/, @file_names);
		#run
		foreach my $RC_csv_name(@RC_csv_names){#3
			sub_basic::refresh_R_script($variables{Tscore_R_file}, "R_RC_csv_name" , $RC_csv_name);
			#run R script
			print "\tRead $RC_csv_name, and run $variables{Tscore_R_file}\n";
			system("Rscript $variables{Tscore_R_file}");
		}#3
		#
}

#######################################################
#4C-seq
sub S2_Colocalization_4C{
	my ($variables_pointer, $raw_file)=@_;
	my %variables=%$variables_pointer;
	my $viewpoint_info_pointer=$variables{viewpoint_info_pointer};
	my %viewpoint_info=%$viewpoint_info_pointer;
	#viewpoints
	my $sample_name=$variables{sample_name};
	my $view_seq=$viewpoint_info{$sample_name}->{view_seq};
	my $region_seq=$viewpoint_info{$sample_name}->{region_seq};
	my $vp_chr=$viewpoint_info{$sample_name}->{chr};
	my $vp_pos=$viewpoint_info{$sample_name}->{site_position},
	my $vp_sbn=$viewpoint_info{$sample_name}->{size_bin_name};
	my $vp_ebn=$viewpoint_info{$sample_name}->{enzyme_bin_name},
	
	my (%end, %num);
	$num{raw_reads_num}=1;
	printf("\t###Read %s\n", $raw_file);
	open my($IN), "<", $raw_file or die;
	my ($name_line, $seq_line, $third_line, $qual_line);
	while (defined($name_line = <$IN>) && defined($seq_line = <$IN>) && defined($third_line = <$IN>) && defined($qual_line = <$IN>) ){#3
		chomp($name_line, $seq_line, $third_line, $qual_line );
		
		#find view seq
		#printf("####@@#%s\n", $region_seq);
		my $tail_seq=sub_basic::split_hybrid_seq($seq_line, $region_seq, $variables{enzyme_site}, 'tail');
		if ($tail_seq eq 'NA'){
			$num{original_regions_num}++ ;
		}
		elsif($tail_seq cmp $seq_line){#4
			$end{$num{raw_reads_num}}->{ref}='*';
			$end{$num{raw_reads_num}}->{trunc_seq}=$tail_seq;
			$end{$num{raw_reads_num}}->{seq}=$seq_line;
			$end{$num{raw_reads_num}}->{coordinate}=$num{raw_reads_num};
		}#4
		$num{raw_reads_num}++;
	}#3
	close($IN);
	
	printf("\t###genome mapping and co-localization of %s\n", $sample_name);
	my $end_pointer=sub_3C::one_seq_alignment(\%variables, \%end);
	%end=%$end_pointer;
	
	printf("\t###export co-localization of %s\n", $sample_name);
	open my($cis_SITE), ">", $variables{sample_out_file}.'.cis_ligation' or die; # $found=1
	open my($trans_SITE), ">", $variables{sample_out_file}.'.trans_ligation' or die; #$found=1
	open my($un_SITE), ">", $variables{sample_out_file}.'.unmapped_ligation' or die; # $found=2
	open my($multiple_SITE), ">", $variables{sample_out_file}.'.multiple_ligation' or die; # $found=4
	foreach my $coordinate(keys %end){#2
		my $output=join("\t", $coordinate, $sample_name, $vp_chr, $vp_pos, $vp_sbn, $vp_ebn, 
			$end{$coordinate}->{ref}, $end{$coordinate}->{start}, $end{$coordinate}->{end},
			$end{$coordinate}->{start_sbn}, $end{$coordinate}->{end_sbn},$end{$coordinate}->{start_ebn}, 
			$end{$coordinate}->{end_ebn}, $end{$coordinate}->{seq}, $end{$coordinate}->{valid}, );
		#
		my $found=$end{$coordinate}->{valid};
		if ($found==2){
			if($vp_chr eq $end{$coordinate}->{ref}){
				print $cis_SITE "$output\n";
				$num{cis_ligation_num}++;
			}
			else{	
				print $trans_SITE "$output\n";	
				$num{trans_ligation_num}++;
			}
		}
		elsif ($found==1){	
			print $un_SITE "$output\n";
			$num{unmapped_ligation_num}++;
		}
		elsif ($found==4){	
			print $multiple_SITE "$output\n";
			$num{multiple_ligation_num}++;
		}
	}#2
	close($cis_SITE);
	close($trans_SITE);
	close($un_SITE);
	close($multiple_SITE);

	#refresh sample log
	while(my($key, $value)=each(%num) ){
		sub_basic::refresh_log($variables{sample_log}, $key, $value, 'add');
	}
	#
}
####################
sub ligation_file_to_hash2{
	my ($export_format_pointer, $sample_file_head)=@_;
	
	#index of rows and columns for export
	my %export_format=%$export_format_pointer;
	my @row_index=split(',', $export_format{row});
	my @r_row_index=split(',', $export_format{r_row});
	my @col_index=split(',', $export_format{col});
	
	#
	my (%Total, %Unique, %seq_hash);
	foreach my $file_tail('cis_ligation', 'trans_ligation'){#2
		open my($IN), "<", $sample_file_head.'.'.$file_tail or die;
		while (my $line=<$IN>){#3
			chomp($line);
			my @items=split("\t", $line);
			my $row=join('_', map {$items[$_]} @row_index);
			my $r_row=join('_', map {$items[$_]} @r_row_index);
			my $col=join('_', map {$items[$_]} @col_index);
					
			#numeric sequence
			#my $seq=sub_basic::numeric_seq($items[10]);#numeric sequence
			my $seq=$items[10];
			
			#counting
			if($row eq $r_row){
				#update $row_name
				$Total{$row}->{$col}++;
				unless(exists $seq_hash{$seq}){
					$Unique{$row}->{$col}++;
					$seq_hash{$seq}=1;
				}
			}
			#row_2 cmp row_2r
			else{
				#update $row_name
				$Total{$row}->{$col} += 0.5;
				$Total{$r_row}->{$col} += 0.5;
				unless(exists $seq_hash{$seq}){
					$Unique{$row}->{$col} += 0.5;
					$Unique{$r_row}->{$col} += 0.5;
					$seq_hash{$seq}=1;
				}
			}
		}#3 
		close($IN);
	}#2
	my %out;
	$out{Total_pointer}=\%Total;
	$out{Unique_pointer}=\%Unique;
	return(\%out);
}

###############
sub counting_file{
	my ($variables_pointer, $export_pointer, $counting_pointer, $counting_file)=@_;
	my %variables=%$variables_pointer;
	my $size_bins_pointer=$variables{size_bins_pointer};
	my %size_bins=%$size_bins_pointer;
	my $enzyme_bins_pointer=$variables{enzyme_bins_pointer};
	my %enzyme_bins=%$enzyme_bins_pointer;
	my %export=%$export_pointer;
	my %counting=%$counting_pointer;
	#%size_bins or %enzyme_bins
	my %bins=($export{row_attr} eq 'SizeBin') ? %size_bins : %enzyme_bins;
	#
	my @row_names=keys %counting;
	@row_names= sort @row_names;
	my @col_names=map { keys $counting{$_} } @row_names;
	@col_names=List::MoreUtils::uniq @col_names;
	@col_names= sort @col_names;
	
	open my($OUT), ">", $counting_file;
		print $OUT join(',', 'chr', 'bin_no', 'bin_start', 'bin_end', @col_names), "\n";
	foreach my $row(@row_names){#2 one class, one file
		#interacted bins info
		my($chr, $bin_no)=split('_', $row);
		my($bin_start, $bin_end)= ($bins{$chr}->{$bin_no}) ? ($bins{$chr}->{$bin_no}->{bin_start}, $bins{$chr}->{$bin_no}->{bin_end}) : (0, 0);
		
		#count the number of reads
		my @RCs;
		foreach my $col(@col_names){#3
			#RC
			my $RC;
			if(exists $counting{$row}->{$col}) {
				$RC=$counting{$row}->{$col};
			}
			else{
				#rand number
				srand;
				while(1){
						#generate random number 0-1 default
						$RC= int(rand($variables{RC_noise})*100)/100;
						last if $RC>0;
				}
			}
			push(@RCs, $RC);
		}#3
		#export
		#print join(',', $chr, $bin_no, @RCs), "\n" unless $bin_start;
		print $OUT join(',', $chr, $bin_no, $bin_start, $bin_end, @RCs), "\n";
	}#2
	close($OUT);
}
#####################################
sub S5_RC_counting_4C{
	my ($variables_pointer, $export_pointer)=@_;
	my %variables=%$variables_pointer;
	my %export=%$export_pointer;
		
	#read statistics file
	my $statistics_pointer=sub_basic::file_to_hash2($variables{statistics_file});
	my %statistics=%$statistics_pointer; 
	
	#Read genomic ligations
	print "\tRead genomic ligations:\n";
	my(%allTotal, %allUnique, %allNormTotal);
	my @sample_names=split(",", $variables{sample_names});
	foreach my $sample_name(@sample_names){#2
		#print "\t\t$export_format{thread}: $sample_name\n";
		#read ligation file into hash2
		my $out_pointer=sub_3C::ligation_file_to_hash2($export_pointer, $variables{ligations_dir}.'/'.$sample_name);
		my %out=%$out_pointer;
		my $Total_pointer=$out{Total_pointer};
		my %Total=%$Total_pointer;
		my $Unique_pointer=$out{Unique_pointer};
		my %Unique=%$Unique_pointer;
		#normalization
		my $raw_reads_num=$statistics{raw_reads_num}->{$sample_name};
		my $NormTotal_pointer=sub_basic::scaling_normalization(\%Total, $raw_reads_num);
		my %NormTotal=%$NormTotal_pointer;
		
		#combine hash
		my $allTotal_pointer=sub_basic::combine_hash2(\%allTotal, \%Total);
		%allTotal=%$allTotal_pointer;
		my $allUnique_pointer=sub_basic::combine_hash2(\%allUnique, \%Total);
		%allUnique=%$allUnique_pointer;
		my $allNormTotal_pointer=sub_basic::combine_hash2(\%allNormTotal, \%Total);
		%allNormTotal=%$allNormTotal_pointer;
	}#2

	#
	my $file_head=$variables{ligation_frequency_dir}.'/'.$export{row_attr}.'_'.$export{col_attr};
	print "\t and then export genomic interactions by read counts into $file_head.\n";
	#total counting
	sub_3C::counting_file(\%variables, $export_pointer, \%allTotal, $file_head.'_Total_RC.csv');
	#normtotal counting
	sub_3C::counting_file(\%variables, $export_pointer, \%allNormTotal, $file_head.'_NormTotal_RC.csv');
	#unique counting
	sub_3C::counting_file(\%variables, $export_pointer, \%allUnique, $file_head.'_Unique_RC.csv');


}

#####################################
sub S5_RC_counting_3CMTS{
	my ($variables_pointer, $export_pointer)=@_;
	my %variables=%$variables_pointer;
	my %export=%$export_pointer;
		
	#read statistics file
	my $statistics_pointer=sub_basic::file_to_hash2($variables{statistics_file});
	my %statistics=%$statistics_pointer; 
	
	#Read genomic ligations
	print "\tRead genomic ligations:\n";
	my @sample_names=split(",", $variables{sample_names});
	foreach my $sample_name(@sample_names){#2
		#print "\t\t$export_format{thread}: $sample_name\n";
		#read ligation file into hash2
		my $out_pointer=sub_3C::ligation_file_to_hash2($export_pointer, $variables{ligations_dir}.'/'.$sample_name);
		my %out=%$out_pointer;
		my $Total_pointer=$out{Total_pointer};
		my %Total=%$Total_pointer;
		my $Unique_pointer=$out{Unique_pointer};
		my %Unique=%$Unique_pointer;
		#normalization
		my $raw_reads_num=$statistics{raw_reads_num}->{$sample_name};
		my $NormTotal_pointer=sub_basic::scaling_normalization(\%Total, $raw_reads_num);
		my %NormTotal=%$NormTotal_pointer;
		
		#
		my $file_head=$variables{ligation_frequency_dir}.'/'.$sample_name.'_'.$export{row_attr}.'_'.$export{col_attr};
		print "\t and then export genomic interactions by read counts into $file_head.\n";
		#total counting
		sub_3C::counting_file(\%variables, $export_pointer, \%Total, $file_head.'_Total_RC.csv');
		#normtotal counting
		sub_3C::counting_file(\%variables, $export_pointer, \%NormTotal, $file_head.'_NormTotal_RC.csv');
		#unique counting
		sub_3C::counting_file(\%variables, $export_pointer, \%Unique, $file_head.'_Unique_RC.csv');
	}#2



}
##################################################
sub S2_Colocalization_HiC{
	my ($variables_pointer, $R1_fastq)=@_;
	my %variables=%$variables_pointer;
	my $R2_fastq=$R1_fastq;
	$R2_fastq=~s/_R1_/_R2_/;

	#sequences alignment
	my $SAM_file=$variables{sample_out_file}.".sam";
	unless(-e $SAM_file){
		print" \n###Pair-end fastq files alignment  of $R1_fastq  begin:\n";
		print "$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $SAM_file \n";
		system("$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $SAM_file ");
		print "\nSequence alignment output: $SAM_file!\n";
	}
	
	print "read_alignment result from $SAM_file\n";
	my %num;
	open my($ORI), ">", $variables{sample_out_file}.'.original_regions' or die;
	open my($cis_SITE), ">", $variables{sample_out_file}.'.cis_ligation' or die;
	open my($trans_SITE), ">", $variables{sample_out_file}.'.trans_ligation' or die;
	open my($ONE), ">", $variables{sample_out_file}.'.one_alignment' or die;
	open my($NO), ">", $variables{sample_out_file}.'.no_alignment' or die;
	open my($MULTI), ">", $variables{sample_out_file}.'.multiple_alignment' or die;
	open my($IN), "<", $SAM_file or die;
	my ($line_1, $line_2);
	while( defined ($line_1 = <$IN>) && defined ($line_2 = <$IN>) ){ #2 read SAM
		chomp($line_1, $line_2);
		#read two lines per time
		my $single1_pointer=read_single_end(\%variables, $line_1);
		my %single1=%$single1_pointer;
		my $single2_pointer=read_single_end(\%variables, $line_2);
		my %single2=%$single2_pointer;
		#get %pair
		my $pair_pointer=sub_3C::read_paired_end(\%variables, \%single1, \%single2);#subroutine
		my %pair=%$pair_pointer;
		#print "$pair{valid_pair}\n";
		
		#judging
		if ($pair{valid_pair} eq 'one_unique'){#3
			print $ONE $pair{output}, "\n";
			$num{one_alignment_num}++;
		}#3
		elsif ($pair{valid_pair} eq 'multiple'){	#3
			print $MULTI $pair{output}, "\n";
			$num{multiple_alignment_num}++;
		}#3
		elsif ($pair{valid_pair} eq 'unaligned'){#3 #unaligned
			print $NO $pair{output}, "\n";
			$num{no_alignment_num}++;
		}#3
		else{#3#cis_ligation or trans_ligation
			my $output=join("\t", $pair{output}, 'both_unique', );
			if($pair{valid_pair} eq 'cis_pair'){#4
				print $ORI $output, "\n";
				$num{original_regions_num}++;
			}#4
			elsif($pair{valid_pair} eq 'cis_unpair'){#4
				print $cis_SITE $output, "\n";
				#print $output, "\n";
				$num{cis_ligation_num}++;
			}#4
			else{#4   $pair{valid_pair} eq 'trans'
				print $trans_SITE $output, "\n";
				$num{trans_ligation_num}++;
			}#4
			$num{both_alignment_num}++;
		}#3#cis_ligation or trans_ligation
		$num{raw_pairs_num}++;
	} #2 read SAM
	close($IN);
	close($ORI);
	close($cis_SITE);
	close($trans_SITE);
	close($ONE);
	close($NO);
	close($MULTI);
	#refresh sample log
	while(my($key, $value)=each(%num) ){
		sub_basic::refresh_log($variables{sample_log}, $key, $value);
	}
	
}

############
#0=coordinate, 1=chr, 2=start, 3=end, 4=size_bin_no, 5=right_size_bin_no, 
#6=enzyme_bin_no, 7=right_enzyme_bin_no, 8=seq, 9=valid, 
#10=co_chr, 11=co_start, 12=co_end, 13=co_size_bin_no, 14=co_right_size_bin_no, 
#15=co_enzyme_bin_no, 16=co_right_enzyme_bin_no, 17=co_seq, 18=co_valid,
#19=type
sub S5_RC_counting_HiC{
        my ($variables_pointer, $format)=@_;
        my %variables=%$variables_pointer;
        my ($sample_name, $file_tail, $type)=split(',', $format);
        my (%bin_info, $bin1_index, $bin2_index);
	if($type eq 'SizeBin'){
		my $size_bins_pointer=$variables{size_bins_pointer};
		%bin_info=%$size_bins_pointer;
		$bin1_index=4; #start_sbn
		$bin2_index=13; #start_sbn
	}
	else{
		my $enzyme_bins_pointer=$variables{enzyme_bins_pointer};
		%bin_info=%$enzyme_bins_pointer;
		$bin1_index=6;
		$bin2_index=15;
	}
        
        my $ligation_file=$variables{ligations_dir}.'/'.$sample_name.'.'.$file_tail;
	print "\tRead $ligation_file\n";
	my %RC_counting;
	open my($IN), "<", $ligation_file or die;
	while(my $line=<$IN>){#2
		chomp($line);
		my @items=split("\t", $line);
		#chromosome
		my $ref1=$items[1];
		my $bin1=$items[$bin1_index];
		my $ref2=$items[10];
		my $bin2=$items[$bin2_index];
                #counting
		#alway small values of chr or bin_no are on right side
		if($ref1<$ref2 or ($ref1==$ref2 and $bin1<=$bin2) ){
			$RC_counting{$ref1}->{$bin1}->{$ref2}->{$bin2}++;
		}
		else{
			$RC_counting{$ref2}->{$bin2}->{$ref1}->{$bin1}++;
		}
        }#2
        close($IN);

	#export size_bin
	my $out_file=$variables{ligation_frequency_dir}.'/'.$type.'_'.$file_tail.'_'.$sample_name.'_RC.csv';
	print "\t\texport RCs into $out_file\n";
	open my($OUT), ">", $out_file or die;
	print $OUT join(',', 'chr1', 'bin1', 'start1', 'end1', 'chr2', 'bin2', 'start2', 'end2', 'RC'), "\n";
	foreach my $ref1(sort {$a<=>$b} (keys %RC_counting)){#2
		my $pointer1=$RC_counting{$ref1};
		my %hash1=%$pointer1;
		foreach my $bin1(sort {$a<=>$b} (keys %hash1)){
			my $pointer2=$hash1{$bin1};
			my %hash2=%$pointer2;
			foreach my $ref2(sort {$a<=>$b} (keys %hash2)){
				my $pointer3=$hash2{$ref2};
				my %hash3=%$pointer3;
				foreach my $bin2(sort {$a<=>$b} (keys %hash3)){
					my $RC=$RC_counting{$ref1}->{$bin1}->{$ref2}->{$bin2};
					if($RC>=$variables{RC_noise}){
						my $start1=($bin_info{$ref1}->{$bin1}->{bin_start}) ? $bin_info{$ref1}->{$bin1}->{bin_start} : 0;
						my $end1=($bin_info{$ref1}->{$bin1}->{bin_end}) ? $bin_info{$ref1}->{$bin1}->{bin_end} : 0;
						my $start2=($bin_info{$ref2}->{$bin2}->{bin_start}) ? $bin_info{$ref2}->{$bin2}->{bin_start} : 0;
						my $end2=($bin_info{$ref2}->{$bin2}->{bin_end}) ? $bin_info{$ref2}->{$bin2}->{bin_end} : 0;
						print $OUT join(',', $ref1, $bin1, $start1, $end1, $ref2, $bin2, $start2, $end2, $RC), "\n";
						#print join(',', $ref1, $bin1, $start1, $end1, $ref2, $bin2, $start2, $end2, $RC), "\n";
					}
				}
			}
		}
	}#2
	close($OUT);
        #
}


#########################################################
1; 
# make sure the file returns true or require will not succeed!#



