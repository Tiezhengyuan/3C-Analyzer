#! /usr/bin/perl -w
use strict;
use warnings;
use threads;
use Bio::SeqIO;
use Bio::Seq::Quality;
#use Bio::Perl;
use List::MoreUtils;
use List::Util;
use File::Find;
use feature qw{switch};


#the file constains all subroutines required for running the pipeline of 3C_analyzer

#optimization note
#Save memory
#raw data reading were saved in files. no export %read_counts
#adapter removal were saved in files. not %insert_counts


package sub_3C;


################################################################
#get date and time
sub get_time{
	my ($begin_pointer, $end_pointer)=@_;
	my @begin=@$begin_pointer;
	my @end=@$end_pointer;
  
	my ($sec1,$min1,$hour1,$monthday1,$month1,$year1,$weekday1,$yearday1,$isdaylight1)=@begin;
	my ($sec2,$min2,$hour2,$monthday2,$month2,$year2,$weekday2,$yearday2,$isdaylight2)=@end;
	$year1 += 1900;
	$year2 += 1900;
	my $begin_time=$year1.'/'.$month1.'/'.$monthday1.', '.$hour1.':'.$min1.':'.$sec1;
	my $end_time=$year2.'/'.$month2.'/'.$monthday2.', '.$hour2.':'.$min2.':'.$sec2;
	my $duration_m=($year2-$year1)*365*24*60 + ($yearday2-$yearday1)*24*60 + ($hour2-$hour1)*60 + ($min2-$min1);
	my $my_time='Running time: '.$begin_time.'-----'.$end_time.'. Duration: '.$duration_m.'min';
	$my_time='NA' if $duration_m<0;
	return($my_time);
}
####################################
#refresh number of multi-threads
sub dynamic_threads{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	my $threads_num=sub_3C::read_log($variables{var_file}, 'threads_num');
	unless($threads_num=~/[^0-9]/){
		$threads_num=$variables{threads_num} if length($threads_num)==0 or $threads_num==0;
	}
	return($threads_num);
}
############################################
#open parameter info file named "variables.txt" needed for processing. 
sub process_info{#1

	my %variables;  # save parameters for processing
	open my ($INFO), "<", $_[0] or die;
	while (<$INFO>) {#2
		chomp($_);
		my ($name, $value) = split("=", $_); #split on the tabs
		$variables{$name} =$value  if $_=~/=/;  
	}#2
	close($INFO);

	return(\%variables);
}#1

###################################################################
sub Pre_read_sample_info{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	print "raw data directory: $variables{raw_data_dir}\n";  
	my $all_raw_files_pointer=sub_3C::raw_files_list($variables{raw_data_dir}, $variables{raw_file_format});
	my @all_raw_files=@$all_raw_files_pointer;
	
	#generate sample_info_file
	unless(-f $variables{sample_info_file}){#2
		my @sample_names;
		foreach my $file(@all_raw_files){#3
			my @array1=split("/", $file);
			my $file_name=$array1[-1];
			$file_name=~s/\.fastq$|\.fq$//;
			my @array2=split("_", $file_name);
			my $sample_name=join('_', @array2[0..$variables{name_slice}]);
			push(@sample_names, $sample_name) unless List::Util::first {$_ eq $sample_name} @sample_names;
		}#3
		@sample_names= sort @sample_names;
		
		print "generate $variables{sample_info_file} from raw files\n";
		open my($OUT), ">", $variables{sample_info_file} or die;
		foreach (@sample_names){
			print $OUT "$_;$_\n";
		}
		close($OUT);
	}#2
	
	#read sample_info file
	print "The file list of raw data:\n";
	my (%sample_info, %raw_to_sample, @matched_raw_files);
	open my ($IN), "<", $variables{sample_info_file} or die;
	while(<$IN>){
		chomp($_);
		my ($sample_name, $raw_names_str)=split(";", $_);
		my @raw_names=split(",", $raw_names_str);
		
		my $files_size=0;
		my (@raw_files, @raw_file_names, @raw_file_names_head);
		my (@fastq_R1_files, @R1_file_names, @R1_file_names_head, @fastq_R2_files);
		foreach my $raw_name( @raw_names ) {#3
			#print "##$sample_name: @raw_files\n";
			foreach my $raw_file(@all_raw_files){#4
				my @array=split("/", $raw_file);
				my $raw_file_name=$array[-1];
				my($raw_file_name_head, $raw_file_name_tail)=split(/\./, $raw_file_name);
				if ($raw_file_name=~/^$raw_name/i){#5
					push(@raw_files, $raw_file);
					push(@raw_file_names, $raw_file_name);
					push(@raw_file_names_head, $raw_file_name_head);
					if ($raw_file_name=~/_R1_/){
						push(@fastq_R1_files, $raw_file);
						push(@R1_file_names, $raw_file_name);
						push(@R1_file_names_head, $raw_file_name_head);
					}
					elsif ($raw_file_name=~/_R2_/){
						push(@fastq_R2_files, $raw_file);
					}
					my @stats=stat($raw_file);
					$files_size += $stats[7];
					$raw_to_sample{$raw_file}->{'sample_name'}=$sample_name;
					$raw_to_sample{$raw_file}->{'raw_file_name_head'}=$raw_file_name_head;
					push(@matched_raw_files, $raw_file) unless List::Util::first {$_ eq $raw_file} @matched_raw_files;
				}#5 
			}#4
		}#3
		$sample_info{$sample_name}->{'raw_names'}=$raw_names_str;
		$sample_info{$sample_name}->{'raw_files'}=join(',', @raw_files);
		$sample_info{$sample_name}->{'raw_file_names'}=join(',', @raw_file_names);
		$sample_info{$sample_name}->{'raw_file_names_head'}=join(',', @raw_file_names_head);
		$sample_info{$sample_name}->{'files_num'}=@raw_files;
		$sample_info{$sample_name}->{'files_size'}=$files_size;
		$sample_info{$sample_name}->{'fastq_R1_files'}=join(',', @fastq_R1_files);
		$sample_info{$sample_name}->{'R1_file_names'}=join(',', @R1_file_names);
		$sample_info{$sample_name}->{'R1_file_names_head'}=join(',', @R1_file_names_head);
		$sample_info{$sample_name}->{'fastq_R2_files'}=join(',', @fastq_R2_files);
		$sample_info{$sample_name}->{'sample_dir'}=$variables{result_dir}.'/'.$sample_name;
		mkdir($sample_info{$sample_name}->{'sample_dir'}, 0755) unless -d $sample_info{$sample_name}->{'sample_dir'};
		$sample_info{$sample_name}->{'sample_log'}=$variables{log_dir}.'/'.$sample_name.'.log';
		print "$sample_name (files=$sample_info{$sample_name}->{files_num}):\n";
		print "\t$sample_info{$sample_name}->{raw_file_names_head}\n";
	}
	close($IN);
	
	my @R1_raw_files=grep(/_R1_/, @matched_raw_files);
	my $all_files_num=@all_raw_files;
	my $matched_files_num=@matched_raw_files;
	print "#########\n";
	print "\tOf $all_files_num files, a total of $matched_files_num files will be analyzed.\n";
	if (@all_raw_files>@matched_raw_files){
		print "These raw files will not be analyzed:\n";
		foreach my $file(@all_raw_files){
			print "$file\n" unless List::Util::first {$_ eq $file} @matched_raw_files;
		}
	}
	
	my @sample_names=keys %sample_info;
	$variables{sample_names}=join(',', @sample_names);
	$variables{sample_info_pointer}=\%sample_info;
	$variables{raw_to_sample_pointer}=\%raw_to_sample;
	$variables{all_raw_files_pointer}=\@all_raw_files;
	$variables{matched_raw_files_pointer}=\@matched_raw_files;
	$variables{R1_raw_files_pointer}=\@R1_raw_files;
	return(\%variables);
}
###########################################
sub Pre_initiate_monitor_log{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(',', $variables{sample_names});
	my $raw_to_sample_pointer=$variables{raw_to_sample_pointer};
	my %raw_to_sample=%$raw_to_sample_pointer;
	
	my $total_supposed_time=0;
	my $total_files_num=0;
	open my($MON), ">", $variables{time_monitor_file} or die; 
	foreach my $sample_name(@sample_names) {
		my $files_num=$sample_info{$sample_name}->{'files_num'};
		$total_files_num += $files_num;
		my $supposed_time=int($sample_info{$sample_name}->{'files_size'}/228383); #unit is second
		$total_supposed_time += $supposed_time;
		print $MON $sample_name, ":supposed_time=$supposed_time\n";
		print $MON $sample_name, ":beginning_time=NA\n";
		print $MON $sample_name, ":ending_time=NA\n";
		print $MON $sample_name, ":running_time=NA\n";
		print $MON $sample_name, ":file_status=0\n";
		print $MON $sample_name, ":raw_files_num=$files_num\n";
		print $MON $sample_name, ":status=no\n";
	}
	print $MON "Total:supposed_time=$total_supposed_time\n";
	print $MON "Total:beginning_time=", join(",", localtime(time) ), "\n";
	print $MON "Total:ending_time=NA\n";
	print $MON "Total:running_time=NA\n";
	print $MON "Total:file_status=0\n";
	print $MON "Total:raw_files_num=$total_files_num\n";
	print $MON "Total:status=no\n";
	close($MON);
		
	#initiate and clear the monitor.log
	open my($SYS), ">", $variables{system_monitor_file} or die;
	print $SYS join("\t", 'Time', 'Duration(s)', 'CPU_usage(%)', 'Memory_usage(%)'), "\n";
	close($SYS);
	
	#initiate log files
	my $R1_raw_files_pointer=$variables{R1_raw_files_pointer};
	my @raw_files_input=@$R1_raw_files_pointer;
	foreach my $raw_file(@raw_files_input) {
		my $sample_name=$raw_to_sample{$raw_file}->{'sample_name'}; #sample name
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		my $sample_log=$sample_dir.'/'.$raw_to_sample{$raw_file}->{'raw_file_name_head'}.'.log';
		sub_3C::refresh_log($sample_log, "status" , 'no') if -f $sample_log;
	}

}

########################################
sub refresh_sample_log{
	my($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	
	#refresh sample_log
	foreach my $sample_name(keys %sample_info){#2
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		if(-d $sample_dir){#3
			#get all log files
			my $files_pointer=sub_3C::files_list($sample_dir, 'file');
			my @files=@$files_pointer;
			my @log_files=grep(/\.log$/, @files);
			my $log_files_num=@log_files;
			my $R1_files_num=split(',', $sample_info{$sample_name}->{'R1_file_names'} );
			my $all_log_files_num=($R1_files_num==0) ? $sample_info{$sample_name}->{'files_num'} : $R1_files_num;
			my %hash;
			 my %status=('no'=>0, 'on'=>0, 'off'=>0, 'all'=>$all_log_files_num);
			foreach my $log_file(@log_files){#4
				open my($IN), "<", $log_file or die;
				while(my $line=<$IN>){
					chomp($line);
					my($name, $value)=split('=', $line);
					if(exists $hash{$name}){
						if($name=~/_num$/){
							$hash{$name} += $value;
						}
						elsif ($name eq 'beginning_time'){
							my @old_time=split(',', $hash{$name});
							my @new_time=split(',', $value);
							$hash{$name}=$value if sub_3C::get_time(\@old_time, \@new_time) eq 'NA';
						}
						elsif ($name eq 'ending_time'){
							my @old_time=split(',', $hash{$name});
							my @new_time=split(',', $value);
							$hash{$name}=$value unless sub_3C::get_time(\@old_time, \@new_time) eq 'NA';
						}
					}
					else{
						$hash{$name}=$value;
					}
					if ($name eq 'status'){
						$status{$value}++;
					}
				}
				close($IN);
			}#4
			#judge running status: no, on, off
			#print "$status{no}=$status{on}=$status{off}=$status{all}\n";
			$hash{status}='no';
			$hash{running_time}='NA';
			if($status{on}>0){
				$hash{status}='on';
			}
			elsif($status{off}==$status{all}){
				my @beginning_time=split(',', $hash{'beginning_time'});
				my @ending_time=split(',', $hash{'ending_time'});
				$hash{running_time}=sub_3C::get_time(\@beginning_time, \@ending_time);
				$hash{status}='off';
			}
			
			#refresh sample_log and time_monitor_log
			unless ($hash{status} eq 'no'){
				#export to sample_log
				my $sample_log=$variables{log_dir}.'/'.$sample_name.'.log';
				#print "$sample_log\n";
				open my($OUT), ">", $sample_log or die;
				foreach(sort (keys %hash)){
					print $OUT "$_=$hash{$_}\n";
				}
				close($OUT);
				#refresh time_monitor_log
				sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:status", $hash{status});
				sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:file_status", $log_files_num);
				sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:beginning_time" , $hash{beginning_time}) if $hash{beginning_time};
				sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:ending_time" , $hash{ending_time}) if $hash{ending_time};
				sub_3C::refresh_log($variables{time_monitor_file}, "$sample_name:running_time" , $hash{running_time}) if $hash{running_time};
			}
		}#3
	}#2
}

#######################################
#size bin list of genome DNA sequences
sub Pre_size_bin{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	
	
	my %chr_info;
	if(-f $variables{genome_info_csv}){
		my $n=1;
		open my($IN), "<", $variables{genome_info_csv} or die;
		while(<$IN>){
			chomp($_);
			my @items=split(',', $_);
			unless($n==1){
				my $chr=$items[0];
				my $chr_len=$items[1];
				$chr_info{$chr}=$chr_len;
			}
			$n++;
		}
		close($IN);
	}
	else{
		print "Get genome information.\n";
		open my($OUT), ">", $variables{genome_info_csv} or die;
		print $OUT join(',', 'chr', 'chr_len', 'GC_cont', 'enzyme_num', ), "\n";
		my $in_obj = Bio::SeqIO->new(-file => $variables{genome_fasta_file}, -format => 'fasta');
		while (my $seq_obj = $in_obj->next_seq() ) {#3
			my $chr=$seq_obj->display_id;
			my $chr_seq=$seq_obj->seq;
			my $chr_len=$seq_obj->length;
			my $GC_num=map {/G|C/gi} $chr_seq;
			my $GC_cont=int($GC_num*100/$chr_len+0.5);
			my $enzyme_num=map {/$variables{enzyme_site}/gi} $chr_seq;
			# %chr_info
			$chr_info{$chr}=$chr_len;
			print $OUT join(',', $chr, $chr_len, $GC_cont, $enzyme_num, ), "\n";
		}
		close($OUT);
	}
	
	unless(-f $variables{size_bin_csv}){#2
		print "Get size bin information.\n";
		open my($OUT), ">", $variables{size_bin_csv} or die;
		print $OUT join(',', 'bin_name', 'chr', 'bin_no', 'bin_start', 'bin_end', 'enzyme_num', 'GC_cont'), "\n";
		my $in_obj = Bio::SeqIO->new(-file => $variables{genome_fasta_file}, -format => 'fasta');
		while (my $seq_obj = $in_obj->next_seq() ) {#3
			my $chr=$seq_obj->display_id;
			my $chr_seq=$seq_obj->seq;
			my $chr_len=$seq_obj->length;
			
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
		}#3
		close($OUT);
	}#2
	
	return(\%chr_info);
}

#######################################
#enzyme list of genome DNA sequences
sub Pre_enzyme_bin{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	
	## enzyme sites counting along chromosomes 
	unless (-f $variables{enzyme_sites_file}){#2
		print "get sequences of human genome DNA and search enzyme sites:\n";
		open my($ENZ), ">", $variables{enzyme_sites_file} or die;
		my $in_obj = Bio::SeqIO->new(-file => $variables{genome_fasta_file}, -format => 'fasta');
		while (my $seq_obj = $in_obj->next_seq() ) {#3
			my $chr=$seq_obj->display_id;
			my $seq=$seq_obj->seq;
			my $seq_len=$seq_obj->length;
			#export enzyme sites along a chromosome
			my @sites=split(",", sub_3C::enzyme_sites_counting($seq, $variables{enzyme_site}, 1) );
			my $sites_num=@sites;
			my $chr_sites_str=join(',', @sites);
			print $ENZ join("\t", $chr, $seq_len, $sites_num, $chr_sites_str), "\n";
			print "\t$chr......\n";
		}#4
		close($ENZ);
	}#2
	
	#count enzyme sites along chromosomes
	print "export $variables{enzyme_bin_csv}\n";
	my %enzyme_sites;
	open my($EN_bin), ">", $variables{enzyme_bin_csv} or die;
	#column names in enzyme_bin_file
	print $EN_bin join(',', 'bin_name', 'chr', 'chr_region', 'bin_no', 'bin_start', 'bin_end', 'enzyme_pos', 'bin_len', 'size_bin_name'), "\n";
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
			if( ($site-$sites_2[0])>=100 ){
				push(@sites_2, $site);
				push(@sites_1, join('_', @sites_2) );
				@sites_2=$site;#clear @sites_2 and give the first element
			}
			else{
				push(@sites_2, $site); #combine neigbouring enzyme site <100bp
			}
		}#3
		
		#
		my $i=0;
		foreach my $start_end(@sites_1){#3
			$i++;
			my @site_start_end=split('_', $start_end);
			my $start_site=$site_start_end[0];
			my $end_site=$site_start_end[-1];
			my $chr_region=int($start_site/1e5)+1;
			my $site_name=$chr.'_'.$i;
			my $bin_len=$end_site-$start_site;
			my $size_bin_name=$chr.'_'.(int($start_site/$variables{size_bin_len})+1);
			#pop @site_start_end; #remove the last
			print $EN_bin join(',', $site_name, $chr, $chr_region, $i, $start_site, $end_site, $start_end, $bin_len, $size_bin_name), "\n";
			#export %enzyme_sites
			$enzyme_sites{$chr}->{$chr_region}->{$i}->{left_site_pos}=$start_site;
			$enzyme_sites{$chr}->{$chr_region}->{$i}->{right_site_pos}=$end_site;
			my $r_chr_region=int($end_site/1e5)+1;
			unless($r_chr_region==$chr_region){
				$enzyme_sites{$chr}->{$r_chr_region}->{$i}->{left_site_pos}=$start_site;
				$enzyme_sites{$chr}->{$r_chr_region}->{$i}->{right_site_pos}=$end_site;
			}
		}#3
	}#2
	close($EN);
	close($EN_bin);
	
	return(\%enzyme_sites);
}
############################################################
sub Pre_check_bowtie_index{
	my ($variables_pointer, $version)=@_;
	my %variables=%$variables_pointer;
	
	#check index files
	my $found=1;
	if($version eq 'bowtie2'){#2
		my @tails=('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2');
		foreach (@tails){
			$found=0 unless -f $variables{genome_index}.$_;
		}
		#build index 
		if($found==0){
			my $bowtie_build_script=$variables{alignment_dir}.'/bowtie2-build';
			system("$bowtie_build_script $variables{genome_fasta_file} $variables{genome_index}");
		}
	}#2

}
#################################
#calculate free space of $raw_data_dir and $result_dir
sub Pre_free_space {
	my ($raw_data_dir, $result_dir)=@_;
	
	my (@incrusive_raw_data_files);
	my $raw_data_size=0;
	#incrusive files
	File::Find::find(sub {
			my $name=$File::Find::name;
			push(@incrusive_raw_data_files, $name) if -f $name;
		}, $raw_data_dir);
	#raw data size
	foreach my $raw_file(@incrusive_raw_data_files){
		my @stats=stat($raw_file);
		$raw_data_size+= $stats[7];
	}
	my $supposed_result_space=int($raw_data_size*3/(1024**3)+0.5);
	
	#free space of result dir
	open my ($SPA), "df $result_dir | " or die;
	my $result_free_space=0;
	my($l1, $l2);
	while( defined($l1=<$SPA>) && defined($l2=<$SPA>) ){
		chomp($l1, $l2);
		my ($disk_id, $total_space, $used_space, $free_space, $space_ratio, $mount_dir)=split(" ", $l2);
		$result_free_space=int($free_space/(1024**2)+0.5);
	}
	close($SPA);
	
	my $space_info="Free space of $result_dir is $result_free_space GB, and $supposed_result_space GB might be used.\n";
	return($space_info);
}

##############################################################
#
sub Pre_site_info{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	my $enzyme_len=length($variables{enzyme_site});

	unless (-f $variables{probe_seq_csv}){#2
		print "read sequences of genome DNA:\n";
		my %genome;
		my $in_obj = Bio::SeqIO->new(-file => $variables{genome_fasta_file}, -format => 'fasta');
		while (my $seq_obj = $in_obj->next_seq() ) {
			my $chr=$seq_obj->display_id;
			my $seq=$seq_obj->seq;
			my $seq_len=$seq_obj->length;
			$genome{$chr}=$seq;
			print "\t$chr......\n";
		}

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
				my $fragment=sub_3C::reverse_complement($revcom_fragment);#EcorI is on the end
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
	my %site_info;
	open my($IN), "<", $variables{probe_seq_csv} or die;
	while(<$IN>){#2
		chomp($_);
		my @array=split(',', $_);
		my $display_id=$array[0];
		my $chr=$array[1];
		my $site_pos=$array[2];
		$site_info{$display_id}->{chr}=$chr;#chr
		$site_info{$display_id}->{site_position}=$site_pos;
		$site_info{$display_id}->{region_start}=$array[3];
		$site_info{$display_id}->{region_end}=$array[4];
		$site_info{$display_id}->{boundary_start}=$array[5];
		$site_info{$display_id}->{boundary_end}=$array[6];
		$site_info{$display_id}->{shortest_len}=$array[7]; 
		$site_info{$display_id}->{fragment}=$array[8];
		$site_info{$display_id}->{shortest_seq}=substr($site_info{$display_id}->{fragment}, -($site_info{$display_id}->{shortest_len}+$enzyme_len) );
		$site_info{$display_id}->{revcom_fragment}=sub_3C::reverse_complement( $site_info{$display_id}->{fragment} );
		#site region localization
		$site_info{$display_id}->{size_bin_name}=int($site_pos/$variables{size_bin_len})+1;
		my $corr_site_pos=($display_id=~/u$/) ? ($site_pos-1): ($site_pos+1);
		$site_info{$display_id}->{enzyme_bin_name}=sub_3C::enzyme_bin_localization(\%variables, $chr, $corr_site_pos);
	}#2
	close($IN);
	
	unless (-f $variables{viewpoints_csv}){
		print "export viewpoints information\n";
		open my($OUT), ">", $variables{viewpoints_csv} or die;
		print $OUT join(',', 'viewpoints', 'sample', 'display_id', 'chr', 'site_pos', 'region_start', 'region_end', 'size_bin_name', 'enzyme_bin_name'), "\n";
		my @sample_names=split(",", $variables{sample_names});
		foreach my $sample(@sample_names){
			foreach my $display_id(sort (keys %site_info) ){
				my $viewpoint=$sample.'_'.$display_id;
				print $OUT join(',', $viewpoint, $sample, $display_id, 
							$site_info{$display_id}->{chr}, $site_info{$display_id}->{site_position}, 
							$site_info{$display_id}->{region_start}, $site_info{$display_id}->{region_end}, 
							$site_info{$display_id}->{size_bin_name}, $site_info{$display_id}->{enzyme_bin_name}, ), "\n";
			}
		}
		close($OUT);
	}
	
	return(\%site_info);
}
#######################################################################
#
sub Pre_4C_site_info{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	my $enzyme_len=length($variables{enzyme_site});
	
	#read genome information
	unless(-f $variables{probe_seq_csv}){
		print "read sequences of genome DNA:\n";
		my %genome;
		my $in_obj = Bio::SeqIO->new(-file => $variables{genome_fasta_file}, -format => 'fasta');
		while (my $seq_obj = $in_obj->next_seq() ) {
			my $chr=$seq_obj->display_id;
			my $seq=$seq_obj->seq;
			my $seq_len=$seq_obj->length;
			$genome{$chr}=$seq;
			print "\t$chr......\n";
		}
		print "specific known seq for detection.\n";
		my $num=1;
		open my($IN), "<", $variables{site_info_file} or die;
		open my($OUT), ">", $variables{probe_seq_csv} or die;
		while (<$IN>){#3
			chomp($_);
			my($sample_name, $chr, $site_pos, $up_down, $view_len)=split(',', $_);
			
			my $size_bin_name=int($site_pos/$variables{size_bin_len})+1;
			my $enzyme_bin_name=sub_3C::enzyme_bin_localization(\%variables, $chr, $site_pos);
			my $chr_region=int($site_pos/1e5)+1;
			#get sequence around viewpoints
			my $chr_seq=$genome{$chr};
			my $start=$site_pos+1-$view_len+length($variables{enzyme_site});
			my $end=$site_pos-1+length($variables{enzyme_site});
			my $view_seq=substr($chr_seq, $site_pos-1-$view_len+length($variables{enzyme_site}), $view_len);
			my $region_seq=substr($chr_seq, $site_pos-1-$view_len, $view_len*2);
			if($up_down eq 'Down'){
				$start=$site_pos;
				$end=$site_pos-1+$view_len;
				$view_seq=sub_3C::reverse_complement(substr($chr_seq, $site_pos-1, $view_len));
				$region_seq=sub_3C::reverse_complement($region_seq);
			}
			print $OUT join(',', $sample_name, $chr, $site_pos, $size_bin_name, $enzyme_bin_name, $start, $end, $view_seq, $region_seq), "\n";
		}
		close($IN);
		close($OUT);
	}
	
	
	print "get site information. \n\n";
	my %site_info;
	open my($IN), "<", $variables{probe_seq_csv} or die;
	while(<$IN>){#2
		chomp($_);
		my ($sample_name, $chr, $site_pos, $size_bin_name, $enzyme_bin_name, $start, $end, $view_seq, $region_seq)=split(',', $_);
		$site_info{$sample_name}->{site_name}=$sample_name;#chr
		$site_info{$sample_name}->{chr}=$chr;#chr
		$site_info{$sample_name}->{site_position}=$site_pos;
		$site_info{$sample_name}->{view_seq}=$view_seq;
		$site_info{$sample_name}->{region_seq}=$region_seq;
		$site_info{$sample_name}->{region_start}=$start;
		$site_info{$sample_name}->{region_end}=$end;
		$site_info{$sample_name}->{size_bin_name}=$size_bin_name;
		$site_info{$sample_name}->{enzyme_bin_name}=$enzyme_bin_name;
	}#2
	close($IN);
	
	unless (-f $variables{viewpoints_csv}){
		print "export viewpoints information\n";
		open my($OUT), ">", $variables{viewpoints_csv} or die;
		print $OUT join(',', 'viewpoints', 'chr', 'site_pos', 'region_start', 'region_end', 'size_bin_name', 'enzyme_bin_name'), "\n";
		my @sample_names=split(",", $variables{sample_names});
		foreach my $display_id(sort (keys %site_info) ){
			#$sample_name equal to $display_id;
			print $OUT join(',', $display_id, $site_info{$display_id}->{chr}, $site_info{$display_id}->{site_position}, 
							$site_info{$display_id}->{region_start}, $site_info{$display_id}->{region_end}, 
							$site_info{$display_id}->{size_bin_name}, $site_info{$display_id}->{enzyme_bin_name}, ), "\n";

		}
		close($OUT);
	}
	
	return(\%site_info);
}

########################################################################
#function: the specific files list in the specific direcotry
#list_type: files, incrusive_file, incrusive_file_name, file_name, dir, incrusive_dir 
sub files_list{#1
	my ($dir, $list_type)=@_;
	$dir .= '/' unless $dir=~/\/$/;
	my (@file_names, @dir_names, @incrusive_dir, @incrusive_files);
		
	# relative file name and directory name
	opendir(DIR, $dir);
	my @all=readdir(DIR);
	closedir(DIR);
	foreach (@all){
		unless ($_ =~ /^\.|~/){
			if (-d $dir.$_) {	push(@dir_names, $_) ;		}
			else {	push(@file_names, $_) ;		}
		}
	}
	#incrusive files and directories
	File::Find::find(sub {
			my $name=$File::Find::name;
			push(@incrusive_files, $name) if -f $name;
			push(@incrusive_dir, $name) if -d $name;
	}, $dir);

	my @files;
	if ($list_type eq 'dir_name'){
		@files=@dir_names;
	}
	elsif ($list_type eq 'dir'){
		@files=map{$dir.$_} @dir_names;
	}
	elsif ($list_type eq 'incrusive_dir'){
		@files=@incrusive_dir;
	}
	elsif ($list_type eq 'file_name'){
		@files=@file_names;
	}
	elsif ($list_type eq 'file'){
		@files=map {$dir.$_} @file_names;
	}
	elsif ($list_type eq 'incrusive_file'){
		@files=@incrusive_files;
	}
	elsif ($list_type eq 'incrusive_file_name'){
		foreach (@incrusive_files){
			my @array=split("/", $_);
			my $name=$array[-1];
			push(@files, $name);
		}
		@files=List::MoreUtils::uniq @files;
	}

	return(\@files);
}#1

####################################33
#function: raw data files list
#list_type: incrusive_files
sub raw_files_list{#1
	my ($raw_dir, $file_format)=@_;
	$raw_dir .= '/' unless $raw_dir=~/\/$/;
	
	#incrusive files and directories
	my @incrusive_files;
	File::Find::find(sub {
			my $name=$File::Find::name;
			push(@incrusive_files, $name) if -f $name;
	}, $raw_dir);

	my @raw_files;
	if ($file_format eq 'SAM'){
		@raw_files=grep(/\.sam$/i, @incrusive_files);
	}
	elsif ($file_format eq 'BAM'){
		@raw_files=grep(/\.bam$/i, @incrusive_files);
	}
	else{
		@raw_files=grep(/\.fastq$|\.fq$/i, @incrusive_files);
	}

	return(\@raw_files);
}#1
##########################
sub auto_sample_names{#1
	my ($dir, $file_format)=@_;
	$dir .= '/' unless $dir=~/\/$/;

	my (%sample_size, @sample_names);
	if ($file_format eq 'FASTQ'){
		File::Find::find(sub {
				my $file=$File::Find::name;
				if (-f $file and $file=~/\.fastq$|\.fq$/){#4
					my @array1=split("/", $file);
					my $file_name=$array1[-1];
					$file_name=~s/\.fastq$|\.fq$//;
					my @array2=split("_", $file_name);
					my $sample_name;
					if (@array2==1){	$sample_name=$file_name;	}
					elsif (@array2==2){	$sample_name=$array2[0]."_".$array2[1];	}
					else{	$sample_name=$array2[0]."_".$array2[1]."_".$array2[2];	}
					my @args=stat($file);
					$sample_size{$sample_name} += $args[7];
				}#4
		}, $dir);
	}
	elsif ($file_format eq 'SAM'){
		File::Find::find(sub {
				my $file=$File::Find::name;
				if (-f $file and $file=~/\.sam$/i){#4
					my @array1=split("/", $file);
					my $file_name=$array1[-1];
					$file_name=~s/\.sam$//i;
					my @array2=split("_", $file_name);
					my $sample_name;
					if (@array2==1){	$sample_name=$file_name;	}
					elsif (@array2==2){	$sample_name=$array2[0]."_".$array2[1];	}
					else{	$sample_name=$array2[0]."_".$array2[1]."_".$array2[2];	}
					my @args=stat($file);
					$sample_size{$sample_name} += $args[7];
				}#4
		}, $dir);
	}
	elsif ($file_format eq 'BAM'){
		File::Find::find(sub {
				my $file=$File::Find::name;
				if (-f $file and $file=~/\.bam$/i){#4
					my @array1=split("/", $file);
					my $file_name=$array1[-1];
					$file_name=~s/\.bam$//i;
					my @array2=split("_", $file_name);
					my $sample_name;
					if (@array2==1){	$sample_name=$file_name;	}
					elsif (@array2==2){	$sample_name=$array2[0]."_".$array2[1];	}
					else{	$sample_name=$array2[0]."_".$array2[1]."_".$array2[2];	}
					my @args=stat($file);
					$sample_size{$sample_name} += $args[7];
				}#4
		}, $dir);
	}
	#order output by file size
	foreach my $sample_name(sort { $sample_size{$b}<=>$sample_size{$a} } (keys %sample_size) ){
		push(@sample_names, $sample_name);
	}
	return(\@sample_names);
}#1


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
	$out{left_site_info}=$chr.'_'.$left_offset.'_'.$fragment_len;
	
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
############################################
sub refresh_R_script{
	my ($log_file, $r_name, $r_value)=@_;
	
	my @lines;
	#read old data
	open my($IN), "<", $log_file or die;
	while (<$IN>){
		chomp($_);
		if($_=~/^R\_/){
			my($name, $value)=split("=", $_);
			if($name eq $r_name){
				my $r_line=$name."=\"".$r_value."\"";
				push(@lines, $r_line);
			}
			else{	push(@lines, $_);	}
		}
		else{
			push(@lines, $_);
		}
	}
	close($IN);
	
	#export refreshed data
	open my($OUT), ">", $log_file or die;
	foreach (@lines){
		print $OUT "$_\n";
	}
	close($OUT);
}

#############################################
#type= 'replace' or 'add'
sub refresh_log{
	my ($log_file, $r_name, $r_value, $type)=@_;
	
	my %variables;
	#read old data
	if(-f $log_file){
		open my($IN), "<", $log_file or die;
		while (<$IN>){
			chomp($_);
			my($name, $value)=split("=", $_);
			$variables{$name}=$value if $name and $value and length($name)>0 and length($value)>0;
		}
		close($IN);
	}
	
	#refresh new data
	if (exists $variables{$r_name} and $type and $type eq 'add'){
		$variables{$r_name}=$variables{$r_name}+$r_value;
	}
	else{
		$variables{$r_name}=$r_value;
	}
	#export refreshed data
	open my($OUT), ">", $log_file or die;
	foreach my $key( sort {$a cmp $b} (keys %variables ) ){
		print $OUT "$key=$variables{$key}\n";
	}
	close($OUT);

}

#############################################
sub read_log{
	my ($log_file, $variable_name)=@_;
	
	my $variable_value='NA';
	open my($IN), "<", $log_file or die;
	flock $IN, 2;#lock file
	while (<$IN>){
		if($_=~/=/){
			chomp($_);
			my($name, $value)=split("=", $_);
			$variable_value=$value if $name eq $variable_name;
		}
	}
	flock $IN, 8; #unlock file
	close($IN);
	
	return($variable_value);
}
################################################
#
sub hybrid_seq_assembly{
	my($enzyme_site, $seq_1, $seq_2, $share_len)=@_;
	my $enzyme_len=length($enzyme_site);
	my $len_1=length($seq_1);
	my $len_2=length($seq_2);
	my $max_len=($len_1<$len_2) ? $len_1 : $len_2;
	my $revcom_seq_2=sub_3C::reverse_complement($seq_2);
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
			$out{seq_2r}=sub_3C::reverse_complement($out{seq_2r});
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
sub S2_Colocalization_detection{
	my ($variables_pointer, $R1_fastq, $name)=@_;
	my %variables=%$variables_pointer;
	my $site_info_pointer=$variables{site_info_pointer};
	my %site_info=%$site_info_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	my $R2_fastq=$R1_fastq;
	$R2_fastq=~s/_R1_/_R2_/;

	my $SAM_file=$variables{sample_dir}.'/'.$name.".sam";
	#sequences alignment
	unless(-e $SAM_file){
		print "$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $SAM_file \n";
		system("$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $SAM_file ");
	}
	print "\nSequence alignment output: $SAM_file!\n";
	
	print "read_alignment result from $SAM_file\n";
	my $original_regions=0;
	my $cis_ligation=0;
	my $trans_ligation=0;
	my $cis_ligation_both=0;
	my $trans_ligation_both=0;
	my $one_alignment=0;
	my $no_alignment=0;
	my $multiple_alignment=0;
	my $unknown_mapped_pairs=0;
	my $both_alignment=0;
	my $raw_pairs=0;
	open my($ORI), ">", $variables{sample_dir}.'/'.$name.'.original_regions' or die;
	open my($cis_SITE), ">", $variables{sample_dir}.'/'.$name.'.cis_ligation' or die;
	open my($trans_SITE), ">", $variables{sample_dir}.'/'.$name.'.trans_ligation' or die;
	open my($ONE), ">", $variables{sample_dir}.'/'.$name.'.one_alignment' or die;
	open my($NO), ">", $variables{sample_dir}.'/'.$name.'.no_alignment' or die;
	open my($MULTI), ">", $variables{sample_dir}.'/'.$name.'.multiple_alignment' or die;
	open my($UN), ">", $variables{sample_dir}.'/'.$name.'.unknown_mapped_pairs' or die;
	open my($IN), "<", $SAM_file or die;
	my ($line_1, $line_2);
	while( defined ($line_1 = <$IN>) && defined ($line_2 = <$IN>) ){ #2 read SAM
		chomp($line_1, $line_2);
		#read two lines per time
		my $pair_pointer=sub_3C::read_paired_end(\%variables, $line_1, $line_2);#subroutine
		my %pair=%$pair_pointer;
		my $output=join("\t", $pair{coordinate}, 
				$pair{ref_1}, $pair{offset_1}, $pair{seq_1}, $pair{alignment_1},
				$pair{ref_2}, $pair{offset_2}, $pair{seq_2}, $pair{alignment_2},
			);
		#judging
		if ($pair{valid_pair} eq 'one_unique'){#3
			print $ONE $output, "\n";
			$one_alignment++;
		}#3
		elsif ($pair{valid_pair}=~/multiple/){	#3
			print $MULTI "$output\t", "$pair{valid_pair}\n";
			$multiple_alignment++;
		}#3
		elsif ($pair{valid_pair} eq 'unaligned'){#3 #unaligned
			print $NO "$pair{coordinate}\t", "$pair{seq_1}\t", "$pair{seq_2}\n";
			$no_alignment++;
		}#3
		else{#3#cis_ligation or trans_ligation
			my $found=0;
			foreach my $display_id(keys %site_info){#4
				my $chr=$site_info{$display_id}->{chr};#
				my $r_start=$site_info{$display_id}->{region_start};#
				my $r_end=$site_info{$display_id}->{region_end};#
				if($chr eq $pair{ref_1} and $pair{offset_1} > $r_start and $pair{offset_1} < $r_end){
					my $co_out=sub_3C::co_seq_localization(\%variables, $display_id, $pair{coordinate}, $pair{ref_2}, $pair{offset_2}, $pair{seq_2});
					$output=join("\t", $co_out, 'both_unique', );
					$found=1;
				}
				elsif($chr eq $pair{ref_2} and $pair{offset_2} > $r_start and $pair{offset_2} < $r_end){
					my $co_out=sub_3C::co_seq_localization(\%variables, $display_id, $pair{coordinate}, $pair{ref_1}, $pair{offset_1}, $pair{seq_1});
					$output=join("\t", $co_out, 'both_unique', );
					$found=1;
				}
				last if $found==1;
			}#4
			#export
			
			if($found==0){#4
				print $UN $output, "\n";
				$unknown_mapped_pairs++;
			}#4
			elsif($found==1 and $pair{valid_pair} eq 'cis_pair'){#4
				print $ORI $output, "\n";
				$original_regions++;
			}#4
			elsif($found==1 and $pair{valid_pair} eq 'cis_unpair'){#4
				print $cis_SITE $output, "\n";
				#print $output, "\n";
				$cis_ligation++;
				$cis_ligation_both++;
			}#4
			elsif($found==1 and $pair{valid_pair} eq 'trans_unpair'){#4
				print $trans_SITE $output, "\n";
				$trans_ligation++;
				$trans_ligation_both++;
			}#4
			$both_alignment++;
		}#3#cis_ligation or trans_ligation
		$raw_pairs++;
	} #2 read SAM
	close($IN);
	close($ORI);
	close($cis_SITE);
	close($trans_SITE);
	close($ONE);
	close($UN);
	close($NO);
	close($MULTI);
	
	sub_3C::refresh_log($variables{sample_log}, "original_regions", $original_regions, 'add');
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation", $cis_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation", $trans_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation_both", $cis_ligation_both, 'add');
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation_both", $trans_ligation_both, 'add');
	sub_3C::refresh_log($variables{sample_log}, "one_alignment", $one_alignment, 'add');
	sub_3C::refresh_log($variables{sample_log}, "no_alignment", $no_alignment, 'add');
	sub_3C::refresh_log($variables{sample_log}, "multiple_alignment", $multiple_alignment, 'add');
	sub_3C::refresh_log($variables{sample_log}, "unknown_mapped_pairs", $unknown_mapped_pairs, 'add');
	sub_3C::refresh_log($variables{sample_log}, "both_alignment", $both_alignment, 'add');
	sub_3C::refresh_log($variables{sample_log}, "raw_pairs", $raw_pairs, 'add');
}

####################################################
#all sequences from raw data without transfering revcom chain
#Note seq_1 and seq_2 is equal to R1_seq and R2_seq,
# however might be reversed-complementary in order to matching reference seq.
sub read_paired_end{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my @items_1=split("\t", $_[1]);
	my @items_2=split("\t", $_[2]);
	
	if($items_1[1]>$items_2[1]){
		my @larger=@items_1;
		@items_1=@items_2;
		@items_2=@larger; #@item_1 is always R1_seq, and @item_2 is always R2_seq
	}
	
	my %pair;
	#R1
	my @query_info=split(":", $items_1[0]);
	$pair{coordinate}=$query_info[-2].":".$query_info[-1];
	$pair{flag_1}=sprintf("%08b", $items_1[1]); #convert to binary number
	if(substr($pair{flag_1}, 5, 1)==1){#no reported alignment (+4)
		$pair{ref_1}='*';
		$pair{offset_1}=0; 
	}
	else{
		$pair{ref_1}=$items_1[2];
		$pair{offset_1}=$items_1[3]; 
	}
	$pair{MapQ_1}=$items_1[4]; 
	$pair{mate_ref_1}=$items_1[6]; 
	$pair{mate_offset_1}=$items_1[7];
	if(substr($pair{flag_1}, 3, 1)==1){	#reversed reference strand	
		$pair{seq_1}=sub_3C::reverse_complement($items_1[9]);
		$pair{alignment_1}='revcom';
	}
	else{
		$pair{seq_1}=$items_1[9];
		$pair{alignment_1}='itself';
	}
	$pair{len_1}=length($items_1[9]);
	foreach (@items_1){
		my($a, $b, $c)=split(':', $_);
		$pair{$a.'_1'}=$c;
	}
	if (exists $pair{AS_1}){
		if (exists $pair{XS_1}){
			$pair{valid_1}=($pair{AS_1}>$pair{XS_1}) ? 2: 4;# unique : multiple alignment
		}
		else{
			$pair{valid_1}=2;#unique alignment
		}
	}
	else{	$pair{valid_1}=1;	} #no alignment
	
	#R2
	$pair{flag_2}=sprintf("%08b", $items_2[1]);
	if(substr($pair{flag_2}, 5, 1)==1){
		$pair{ref_2}='*';
		$pair{offset_2}=0; 
	}
	else{
		$pair{ref_2}=$items_2[2];
		$pair{offset_2}=$items_2[3]; 
	}
	$pair{MapQ_2}=$items_2[4]; 
	$pair{mate_ref_2}=$items_2[6]; 
	$pair{mate_offset_2}=$items_2[7];
	if(substr($pair{flag_2}, 3, 1)==1){		
		$pair{seq_2}=sub_3C::reverse_complement($items_2[9]);		
		$pair{alignment_2}='revcom';
	}
	else{
		$pair{seq_2}=$items_2[9];
		$pair{alignment_2}='itself';
	}
	$pair{len_2}=length($items_2[9]);
	
	foreach (@items_2){
		my($a, $b, $c)=split(':', $_);
		$pair{$a.'_2'}=$c;
	}
	if (exists $pair{AS_2}){
		if (exists $pair{XS_2}){
			$pair{valid_2}=($pair{AS_2}>$pair{XS_2}) ? 2: 4;# unique : multiple alignment
		}
		else{
			$pair{valid_2}=2;#unique alignment
		}
	}
	else{	$pair{valid_2}=1;	} #no alignment
	
	#unaligned=1, unique=2, multiple=4
	#all possible: 2,3,4,5,6,8
	my $judging=$pair{valid_1}+$pair{valid_2};
	if ($judging==3){	
		$pair{valid_pair}='one_unique';	
	}
	elsif ($judging==4){	
		if ($pair{ref_1} eq $pair{ref_2}){
			$pair{valid_pair}='cis_unpair';
			$pair{valid_pair}='cis_pair' if $pair{YT_2} eq 'CP';
		}
		else{
			$pair{valid_pair}= 'trans_unpair';
		}
	}
	elsif ($judging==5 or $judging==6 or $judging==8){	
		$pair{valid_pair}='multiple:'.$pair{valid_1}.':'.$pair{valid_2};
	}
	else{	
		$pair{valid_pair}='unaligned';
	}
	
	return(\%pair);
}


###############################################
#
sub S3_site_detection_one_alignment{
	my ($variables_pointer, $name)=@_;
	my %variables=%$variables_pointer;
	my $site_info_pointer=$variables{site_info_pointer};
	my %site_info=%$site_info_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	
	# $found=1 is unique alignment
	my %pair;
	my $cis_ligation=0;
	my $trans_ligation=0;
	my $cis_ligation_one=0;
	my $trans_ligation_one=0;
	my $multiple_ligation=0;
	my $unmapped_ligation=0;
	open my($cis_SITE), ">>", $variables{sample_dir}.'/'.$name.'.cis_ligation' or die; # $found=2
	open my($trans_SITE), ">>", $variables{sample_dir}.'/'.$name.'.trans_ligation' or die; #$found=4
	open my($un_SITE), ">>", $variables{sample_dir}.'/'.$name.'.unmapped_ligation' or die; # $found=8
	open my($multiple_SITE), ">>", $variables{sample_dir}.'/'.$name.'.multiple_ligation' or die; # $found=16
	#read one alignment file
	my $fasta_file=$variables{sample_dir}.'/'.$name.'_one_alignment.tmp_fasta';
	my $sam_file=$variables{sample_dir}.'/'.$name.'_one_alignment.tmp_sam';
	open my($OUT), ">", $fasta_file or die;
	open my($IN), "<", $variables{sample_dir}.'/'.$name.'.one_alignment' or die;
	while(my $line=<$IN>){#1 *.one_alignment circling
		chomp($line);
		my($coordinate, $ref_1, $offset_1, $seq_1, $alignment_1, $ref_2, $offset_2, $seq_2, $alignment_2)=split("\t", $line);
		
		my $found=0;
		my $out_line;
		#R1 is unmapped and R2 is mapped
		if($ref_1 eq '*'){#2 R1 is unmapped and R2 is mapped
			foreach my $display_id(keys %site_info){#3  #site_info
				my $chr=$site_info{$display_id}->{chr};
				my $region_start=$site_info{$display_id}->{region_start};
				my $region_end=$site_info{$display_id}->{region_end};
				#R1 is hybrid seq of viewpoint seq and co-seq, and R2 is a viewpoint seq
				if( $chr eq $ref_2 and $offset_2 > $region_start and $offset_2 < $region_end){#4
					$found=1;
					my $seq_trunc=end3_seq_truncation($site_info{$display_id}->{revcom_fragment}, $seq_1);
					if ($seq_1 cmp $seq_trunc and length($seq_trunc) >= $variables{query_length}){#5
						$pair{$coordinate}->{view_displayid}=$display_id;
						$pair{$coordinate}->{view_ref}=$ref_2;
						$pair{$coordinate}->{view_offset}=$offset_2;
						$pair{$coordinate}->{view_seq}=$seq_2;
						$pair{$coordinate}->{co_seq}=$seq_trunc;
						print $OUT ">$coordinate\n", "$seq_trunc\n";
					}#5
					else{
						$out_line=join("\t", , $coordinate, $seq_trunc, 'one_unique');
						$found=8; #default is un_SITE
					}
				}#4
				#R1 is a hybrid seq of viewpoint seq and co-seq, and R2 is a co-seq
				elsif($seq_1=~/$site_info{$display_id}->{shortest_seq}/){#4
					my $co_out=sub_3C::co_seq_localization(\%variables, $display_id, $coordinate, $ref_2, $offset_2, $seq_2);
					$out_line=join("\t", $co_out, 'one_unique', );
					$found=($ref_2 eq $chr) ? 2 : 4; #cis_SITE is 2 and trans_SITE is 4
				}#4
				last if $found>0;
			}#3 #site_info
		}#2 R1 is unmapped and R2 is mapped
		elsif($ref_2 eq '*'){#2 R2 is unmapped and R1 is mapped
			foreach my $display_id(keys %site_info){#3  #site_info
				my $chr=$site_info{$display_id}->{chr};
				my $region_start=$site_info{$display_id}->{region_start};
				my $region_end=$site_info{$display_id}->{region_end};
				#R1 is viewpoint seq, and R2 is hybrid seq of viewpoint seq and co-seq
				if( $chr eq $ref_1 and $offset_1 > $region_start and $offset_1 < $region_end){#4
					$found=1;
					my $seq_trunc=end3_seq_truncation($site_info{$display_id}->{revcom_fragment}, $seq_2);
					if ($seq_2 cmp $seq_trunc and length($seq_trunc) >= $variables{query_length}){#5
						$pair{$coordinate}->{view_displayid}=$display_id;
						$pair{$coordinate}->{view_ref}=$ref_1;
						$pair{$coordinate}->{view_offset}=$offset_1;
						$pair{$coordinate}->{view_seq}=$seq_1;
						$pair{$coordinate}->{co_seq}=$seq_trunc;
						print $OUT ">$coordinate\n", "$seq_trunc\n";
					}#5
					else{
						$out_line=join("\t", , $coordinate, $seq_trunc, 'one_unique');
						$found=8; #default is un_SITE
					}
				}#4
				#R2 is hybrid seq of unknown seq and known seq, and R1 is a co-regulated seq
				elsif($seq_2=~/$site_info{$display_id}->{shortest_seq}/){#4
					my $co_out=sub_3C::co_seq_localization(\%variables, $display_id, $coordinate, $ref_1, $offset_1, $seq_1);
					$out_line=join("\t", $co_out, 'one_unique', );
					$found= ($ref_1 eq $chr) ? 2 : 4 ##cis is 2 and trans is 4
				}#4
				last if $found>0;
			}#3 site_info
		}#2 R2 is unmapped and R1 is mapped
		
		#export
		if ($found==2)		{	
			print $cis_SITE "$out_line\n";
			$cis_ligation++;
			$cis_ligation_one++;
		}
		elsif ($found==4)	{	
			print $trans_SITE "$out_line\n";
			$trans_ligation++;
			$trans_ligation_one++;
		}
		elsif ($found==8)	{	
			print $un_SITE "$out_line\n";
			$unmapped_ligation++;
		}
	}#1  *.one_alignment circling
	close($IN);
	close($OUT);
	
	#export
	my $co_pair_pointer=sub_3C::one_seq_alignment(\%variables, \%pair, $name, $fasta_file, $sam_file, 'one_unique');
	my %co_pair=%$co_pair_pointer;
	foreach my $coordinate(keys %co_pair){#2
		my $found=$co_pair{$coordinate}->{valid_seq};
		if ($found==2)		{	
			print $cis_SITE "$co_pair{$coordinate}->{out_line}\n";
			$cis_ligation++;
			$cis_ligation_one++;
		}
		elsif ($found==4)	{	
			print $trans_SITE "$co_pair{$coordinate}->{out_line}\n";	
			$trans_ligation++;
			$trans_ligation_one++;
		}
		elsif ($found==8)	{	
			print $un_SITE "$co_pair{$coordinate}->{out_line}\n";
			$unmapped_ligation++;
		}
		elsif ($found==16)	{	
			print $multiple_SITE "$co_pair{$coordinate}->{out_line}\n";
			$multiple_ligation++;
		}
	}#2
	close($cis_SITE);
	close($trans_SITE);
	close($un_SITE);
	close($multiple_SITE);
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation", $cis_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation", $trans_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation_one", $cis_ligation_one, 'add');
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation_one", $trans_ligation_one, 'add');
	sub_3C::refresh_log($variables{sample_log}, "unmapped_ligation", $unmapped_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "multiple_ligation", $multiple_ligation, 'add');
	unlink($fasta_file, $sam_file);
}
#######################################
sub one_seq_alignment{
	my($variables_pointer, $pair_pointer, $sample_name, $fasta_file, $sam_file, $alignment_type)=@_;
	my %variables=%$variables_pointer;
	my %pair=%$pair_pointer;
	my $site_info_pointer=$variables{site_info_pointer};
	my %site_info=%$site_info_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	
	#make alignment
	print "$variables{bowtie_options} -f $fasta_file -S $sam_file";
	system("$variables{bowtie_options} -f $fasta_file -S $sam_file");
	
	#read sam file
	open my($IN), "<", $sam_file or die;
	while(my $line=<$IN>){#2
		chomp($line);
		my @items=split("\t", $line);
		my $coordinate=$items[0];
		
		my %out;
		#flag number
		#$out{flag}=sprintf("%08b", $items[1]); #convert to binary number
		#ref and offset
		#if(substr($out{flag}, 5, 1)==1){ #no reported alignment
		#	$out{ref}='*';
		#	$out{offset}=0; 
		#}
		#alignment direction: itself or revcom
		#$out{alignment}=(substr($out{flag}, 3, 1)==1)? 'revcom' : 'itself';
		#mapping quality
		#$out{MapQ}=$items[4]; 
		#if($out{MapQ} >= $variables{mapping_quality}){}
		
		#Optional fields
		foreach (@items){
			my($a, $b, $c)=split(':', $_);
			$out{$a}=$c;
		}
		#unmapped alignment: 1, unqiue alignment is 2, multiple alignmet is 16
		my $found=1;
		if (exists $out{AS}){#AS:i: is alignment score means there are alignments
			if (exists $out{XS}){#XS:i: is second alignment score means there might be multiple alignments
				$found=($out{AS}==$out{XS}) ? 16: 2;
			}
			else{
				$found=2;
			}
		}
				
		#export of valid_seq: cis alignment is 2, trans alignment 4, no alignment is 8, multiple alignmet is 16
		if($found==2){#unique alignment
			my $chr=$items[2];
			my $offset=$items[3]; 
			my $end_offset=$offset+length($pair{$coordinate}->{co_seq}); 
			$pair{$coordinate}->{co_ref}=$chr;
			$pair{$coordinate}->{co_offset}=$offset;
			#$pair{$coordinate}->{co_size_bin_name}=int($offset/$variables{size_bin_len})+1;
			#$pair{$coordinate}->{right_co_size_bin_name}=int($end_offset/$variables{size_bin_len})+1;
			#$pair{$coordinate}->{co_enzyme_bin_name}=sub_3C::enzyme_bin_localization(\%variables, $chr, $offset);
			#$pair{$coordinate}->{right_co_enzyme_bin_name}=sub_3C::enzyme_bin_localization(\%variables, $chr, $end_offset);
			#export
			$pair{$coordinate}->{valid_seq}= ($pair{$coordinate}->{view_ref} eq $pair{$coordinate}->{co_ref}) ? 2 : 4;
			my $co_out=sub_3C::co_seq_localization(\%variables, $pair{$coordinate}->{view_displayid}, $coordinate, 
				$pair{$coordinate}->{co_ref}, $pair{$coordinate}->{co_offset}, $pair{$coordinate}->{co_seq}); 
			$pair{$coordinate}->{out_line}=join("\t", $co_out, $alignment_type, );
		}
		elsif ($found==16){ #multiple alignment
			$pair{$coordinate}->{valid_seq}= 16;
			$pair{$coordinate}->{out_line}=join("\t", $pair{$coordinate}->{view_displayid}, 
							$coordinate, $pair{$coordinate}->{co_seq}, $alignment_type);
		}
		else{	#no alignment
			$pair{$coordinate}->{valid_seq}=8;
			$pair{$coordinate}->{out_line}=join("\t", $pair{$coordinate}->{view_displayid}, 
							$coordinate, $pair{$coordinate}->{co_seq}, $alignment_type);
		}
		#print "$coordinate:$pair{$coordinate}->{out_line}\n";
	}#2
	
	return(\%pair);
}
##########################################
sub reverse_complement {
	my $dna_seq = $_[0];

	# reverse the DNA sequence
	my $revcom_seq = reverse($dna_seq);
	# complement the reversed DNA sequence
	$revcom_seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

	return ($revcom_seq);
}
###############################################################
#
sub co_seq_localization{
	my($variables_pointer, $display_id, $coordinate, $ref, $offset, $seq)=@_;
	my %variables=%$variables_pointer;
	my $site_info_pointer=$variables{site_info_pointer};
	my %site_info=%$site_info_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	my $end_offset=$offset+length($seq);
	
	my $size_bin_name=$site_info{$display_id}->{size_bin_name};#
	my $enzyme_bin_name=$site_info{$display_id}->{enzyme_bin_name};#
	my $co_size_bin_name=int($offset/$variables{size_bin_len})+1;
	my $right_co_size_bin_name=int($end_offset/$variables{size_bin_len})+1;
	my $co_enzyme_bin_name=sub_3C::enzyme_bin_localization(\%variables, $ref, $offset);
	my $right_co_enzyme_bin_name=sub_3C::enzyme_bin_localization(\%variables, $ref, $end_offset);

	my $out=join("\t", $display_id, $coordinate, $size_bin_name, $enzyme_bin_name,
					$co_size_bin_name, $right_co_size_bin_name, 
					$co_enzyme_bin_name, $right_co_enzyme_bin_name, 
					$ref, $offset, $seq,);
	return($out);
}

##########################################
sub enzyme_bin_localization{
	my($variables_pointer, $chr, $offset)=@_;
	my %variables=%$variables_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	
	my $enzyme_bin_name=0;
	my $chr_region=int($offset/1e5)+1;
	until(exists $enzyme_sites{$chr}->{$chr_region}){
			$chr_region--;
	}
	my $chr_sites_pointer=$enzyme_sites{$chr}->{$chr_region};
	my %chr_sites=%$chr_sites_pointer;
	foreach my $site_no(keys %chr_sites){#3
		if ($offset>=$chr_sites{$site_no}->{left_site_pos} and $offset<$chr_sites{$site_no}->{right_site_pos}){
			$enzyme_bin_name=$site_no;
			#print "$chr:$chr_region=$site_no=";
			#print "$chr_sites{$site_no}->{left_site_pos}:$offset:$chr_sites{$site_no}->{right_site_pos}\n";
			last;
		}
	}#3
	return($enzyme_bin_name);
}
###############################################
#
sub S3_site_detection_no_alignment{
	my ($variables_pointer, $name)=@_;
	my %variables=%$variables_pointer;
	my $site_info_pointer=$variables{site_info_pointer};
	my %site_info=%$site_info_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	
	#read one alignment file
	my %pair;
	open my($IN), "<", $variables{sample_dir}.'/'.$name.'.no_alignment' or die;
	my $fasta_file=$variables{sample_dir}.'/'.$name.'_no_alignment.tmp_fasta';
	open my($OUT), ">", $fasta_file or die;
	while(my $line=<$IN>){#2 *.no_alignment circling
		chomp($line);
		my($coordinate, $seq_1, $seq_2)=split("\t", $line);
		
		my $out_pointer=hybrid_seq_assembly($variables{enzyme_site}, $seq_1, $seq_2, 15);
		my %out=%$out_pointer;
		if($out{hybrid} eq 'yes'){#3
			foreach my $display_id(keys %site_info){#4
				my $chr=$site_info{$display_id}->{chr};
				my $site_position=$site_info{$display_id}->{site_position};
				my $fragment =$site_info{$display_id}->{fragment};
				
				if($fragment=~/$out{seq_1r}/){
					$pair{$coordinate}->{view_displayid}=$display_id;
					$pair{$coordinate}->{view_ref}=$chr;
					$pair{$coordinate}->{view_offset}=$site_position;
					$pair{$coordinate}->{view_seq}=$out{seq_1r};
					$pair{$coordinate}->{co_seq}=$out{seq_2r};
					#print  "$seq_1\n$seq_2\n", "$out{seq_1r}\n$out{seq_2r}\n\n";
					print $OUT ">$coordinate\n", "$out{seq_2r}\n";
					last;
				}
				elsif($fragment=~/$out{seq_2r}/){
					$pair{$coordinate}->{view_displayid}=$display_id;
					$pair{$coordinate}->{view_ref}=$chr;
					$pair{$coordinate}->{view_offset}=$site_position;
					$pair{$coordinate}->{view_seq}=$out{seq_2r};
					$pair{$coordinate}->{co_seq}=$out{seq_1r};
					
					print $OUT ">$coordinate\n", "$out{seq_1r}\n";
					last;
				}
			}#4
		}#3
	}#2  *.no_alignment circling
	close($IN);
	close($OUT);
	
	#export
	my $cis_ligation=0;
	my $trans_ligation=0;
	my $cis_ligation_no=0;
	my $trans_ligation_no=0;
	my $multiple_ligation=0;
	my $unmapped_ligation=0;
	open my($cis_SITE), ">>", $variables{sample_dir}.'/'.$name.'.cis_ligation' or die;
	open my($trans_SITE), ">>", $variables{sample_dir}.'/'.$name.'.trans_ligation' or die;
	open my($un_SITE), ">>", $variables{sample_dir}.'/'.$name.'.unmapped_ligation' or die;
	open my($multiple_SITE), ">>", $variables{sample_dir}.'/'.$name.'.multiple_ligation' or die;
	my $sam_file=$variables{sample_dir}.'/'.$name.'_no_alignment.tmp_sam';
	my $co_pair_pointer=sub_3C::one_seq_alignment(\%variables, \%pair, $name, $fasta_file, $sam_file, 'no_alignment');
	my %co_pair=%$co_pair_pointer;
	foreach my $coordinate(keys %co_pair){#2
		my $found=$co_pair{$coordinate}->{valid_seq};
		if ($found==2)		{
			print $cis_SITE "$co_pair{$coordinate}->{out_line}\n";
			$cis_ligation++;
			$cis_ligation_no++;
		}
		elsif ($found==4)	{
			print $trans_SITE "$co_pair{$coordinate}->{out_line}\n";
			$trans_ligation++;
			$trans_ligation_no++;
		}
		elsif ($found==8)	{
			print $un_SITE "$co_pair{$coordinate}->{out_line}\n";
			$multiple_ligation++;
		}
		elsif ($found==16)	{	
			print $multiple_SITE "$co_pair{$coordinate}->{out_line}\n";
			$unmapped_ligation++;
		}
	}#2
	close($cis_SITE);
	close($trans_SITE);
	close($un_SITE);
	close($multiple_SITE);
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation", $cis_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation", $trans_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation_no", $cis_ligation_no, 'add');
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation_no", $trans_ligation_no, 'add');
	sub_3C::refresh_log($variables{sample_log}, "unmapped_ligation", $unmapped_ligation, 'add');
	sub_3C::refresh_log($variables{sample_log}, "multiple_ligation", $multiple_ligation, 'add');
	unlink($fasta_file, $sam_file);
}

###############################################
#
sub S3_site_detection_multiple_alignment{
	my ($variables_pointer, $name)=@_;
	my %variables=%$variables_pointer;
	my $site_info_pointer=$variables{site_info_pointer};
	my %site_info=%$site_info_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	
	my $multiple_ligation=0;
	open my($multiple_SITE), ">", $variables{sample_dir}.'/'.$name.'.multiple_ligation' or die;
	open my($IN), "<", $variables{sample_dir}.'/'.$name.'.multiple_alignment' or die;
	while(my $line=<$IN>){#2
		chomp($line);
		my($coordinate, $ref_1, $offset_1, $seq_1, $alignment_1, $ref_2, $offset_2, $seq_2, $alignment_2, $multiple_tag)=split("\t", $line);
		my ($multiple, $R1_judging, $R2_judging)=split(":", $multiple_tag);
		
		my $out_line;
		if($R1_judging==2 and $R2_judging==4){ #3 R1 unique alignment
			foreach my $display_id(keys %site_info){#4 site_info circling
				my $chr=$site_info{$display_id}->{chr};
				my $r_start=$site_info{$display_id}->{region_start};
				my $r_end=$site_info{$display_id}->{region_end};
				#viewpoint localization
				if($ref_1 eq $chr and $offset_1>=$r_start and $offset_1<=$r_end){
					$out_line=join("\t", $display_id, $coordinate, $seq_2, 'multiple_alignment');
					print $multiple_SITE "$out_line\n";
					$multiple_ligation++;
					last;
				}
			}#4
		}#3
		elsif($R1_judging==4 and $R2_judging==2){ #3 R2 unique alignment
			foreach my $display_id(keys %site_info){#4 site_info circling
				my $chr=$site_info{$display_id}->{chr};
				my $r_start=$site_info{$display_id}->{region_start};
				my $r_end=$site_info{$display_id}->{region_end};
				#viewpoint localization
				if($ref_2 eq $chr and $offset_2>=$r_start and $offset_2<=$r_end){
					$out_line=join("\t", $display_id, $coordinate, $seq_1, 'multiple_alignment');
					print $multiple_SITE "$out_line\n";
					$multiple_ligation++;
					last;
				}
			}#4
		}#3
		elsif($R1_judging==1 and $R2_judging==4){ #3 R1 no alignment
			my $right_enzyme_pos=rindex($seq_1, $variables{enzyme_site});
			if($right_enzyme_pos>0){#4
				my $sub_seq=substr($seq_1, 0, $right_enzyme_pos);
				my $sub_seq_len=length($sub_seq);
				foreach my $display_id(keys %site_info){#5 site_info circling
					my $shortest_len=$site_info{$display_id}->{shortest_len};
					my $fragment=$site_info{$display_id}->{fragment};
					my $revcom_fragment=$site_info{$display_id}->{revcom_fragment};
					#viewpoint localization
					if($sub_seq_len >=$shortest_len and ($fragment=~/$sub_seq/ or $revcom_fragment=~/$sub_seq/) ){
						$out_line=join("\t", $display_id, $coordinate, $seq_2, 'multiple_alignment');
						print $multiple_SITE "$out_line\n";
						$multiple_ligation++;
						last;
					}
				}#5
			}#4
		}#3
		elsif($R1_judging==4 and $R2_judging==1){ #3 R2 no alignment
			my $right_enzyme_pos=rindex($seq_2, $variables{enzyme_site});
			if($right_enzyme_pos>0){#4
				my $sub_seq=substr($seq_2, 0, $right_enzyme_pos);
				my $sub_seq_len=length($sub_seq);
				foreach my $display_id(keys %site_info){#5 site_info circling
					my $shortest_len=$site_info{$display_id}->{shortest_len};
					my $fragment=$site_info{$display_id}->{fragment};
					my $revcom_fragment=$site_info{$display_id}->{revcom_fragment};
					#viewpoint localization
					if($sub_seq_len >=$shortest_len and ($fragment=~/$sub_seq/ or $revcom_fragment=~/$sub_seq/) ){
						$out_line=join("\t", $display_id, $coordinate, $seq_1, 'multiple_alignment');
						print $multiple_SITE "$out_line\n";
						$multiple_ligation++;
						last;
					}
				}#5
			}#4
		}#3
	}#2
	close($IN);
	close($multiple_SITE);
	sub_3C::refresh_log($variables{sample_log}, "multiple_ligation", $multiple_ligation, 'add');
}


##################################################################
#truncate sequence of adapter from read
sub end3_seq_truncation{#1
	my ($adapter, $seq)=@_;
	my $adapter_len=length($adapter);
	my $seq_len=length($seq);
	
	my $seq_trunc=$seq;
	unless ($seq_len<=10 or $seq_trunc=~s/$adapter.*//){#2  all adapter sequence
		my $end_cut=($seq_len>$adapter_len) ? $adapter_len: $seq_len;
		for (my$i=$end_cut; $i>4; $i--){
			my $adapter_index=substr($adapter, 0, $i);
			last if $seq_trunc=~s/$adapter_index$//;
		}
	}#2

	return($seq_trunc);
}#1
###################################
#
sub numeric_seq{
	my $seq=$_[0];
	
	$seq=uc($seq);
	$seq=~s/A/1/g;
	$seq=~s/T/2/g;
	$seq=~s/G/3/g;
	$seq=~s/C/4/g;
	$seq=~s/[A-Z]/0/g;
	my @seq_arr=split("", $seq);
	my @power=(0..(@seq_arr-1));
	@power=reverse @power;
	@power=map {5**$_} @power;
	my @numeric_seq_arr=List::MoreUtils::pairwise {our($a, $b); $a*$b} @seq_arr, @power;
	my $numeric_seq=List::Util::sum(@numeric_seq_arr);
	
	return($numeric_seq);
}
#######################################################################################
#distance of two strings for comparison
#algorithm is Needleman-Wunsch
sub NW_dist{#2
  my $str1 = $_[0];
  my $str2 = $_[1];
  my $len1 = length($str1);
  my $len2 = length($str2);

  # We have to be a little careful with the index of the table we are using
  # Since we do not want to reconstruct the alignment, only two rows of the
  # matrix are sufficient. Obviously we will use pointers to these rows, so
  # we can exchange them freely without copying elements.

  # Initialize old_row as reference
  my $old_row;
  for (my $i = 0; $i <= $len1; $i++){#2
    $old_row->[$i] = $i;
  }#2

  my $new_row;
  for (my $i = 1; $i <= $len2; $i++) {#2
    $new_row = [];
    $new_row->[0] = $i;
    for (my $j = 1; $j <= $len1; $j++) {
      $new_row->[$j]=$old_row->[$j]+1;
      my $b=$new_row->[$j-1]+1;
      my $c=$old_row->[$j-1] + ((substr($str1, $j-1, 1) eq substr($str2, $i-1, 1)) ? 0 : 1);
      $new_row->[$j]=$b if $new_row->[$j] > $b;
      $new_row->[$j]=$c if $new_row->[$j] > $c;
    }
    $old_row = $new_row;
  }#2
  return($new_row->[$len1]);
}#2


##############################################################################
#export: statistics.txt 
sub S4_result_statistics{#1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});
	
	mkdir($variables{ligations_dir}, 0755) unless -d $variables{ligations_dir};
	#generate statistics result
	my %counting_hash;
	my @combined_file_tails=('unmapped_ligation', 'cis_ligation', 'trans_ligation', 'multiple_ligation', 'original_regions');
	foreach my $sample_name(@sample_names){#2
		my $sample_dir=$sample_info{$sample_name}->{sample_dir};
		#read sample log
		open my($IN), "<", $sample_info{$sample_name}->{sample_log} or die;
		while(<$IN>){#3
			chomp($_);
			my($key, $value)=split("=", $_);
			$counting_hash{$key}->{$sample_name}=$value;
		}#3
		$counting_hash{raw_files}->{$sample_name}=split(',', $sample_info{$sample_name}->{'raw_files'});
		
		#combine and export files from %combined_file_tails
		my @file_heads=split(',', $sample_info{$sample_name}->{'R1_file_names_head'});
		foreach my $file_tail(@combined_file_tails){#3
			my @sub_files=map {$sample_dir.'/'.$_.'.'.$file_tail} @file_heads;
			my $in_sub_files_str=join (' ', @sub_files);
			my $out_file=$variables{ligations_dir}.'/'.$sample_name.'.'.$file_tail;
			system("cat $in_sub_files_str > $out_file"); #combine
		}#3
	}#2
	
	#statistics
	open my($OUT), ">", $variables{statistics_file} or die;
	print $OUT join("\t", 'item_names', @sample_names), "\n";
	foreach my $item_name(sort (keys %counting_hash) ){
		my @counting_num=map { $counting_hash{$item_name}->{$_} } @sample_names;
		print $OUT join("\t", $item_name, @counting_num), "\n";
	}
	close($OUT);
}

##########################
sub read_statistics{
	my($statistics_file)=@_;
	
	#raw reads
	my (%col_names, %statistics);
	my $n=1;
	open my($STA), "<", $statistics_file or die;
	while (<$STA>){
		chomp($_);
		my @items=split("\t", $_);
		if($n==1){
			for(my $i=0;$i<@items; $i++){
				$col_names{$i}=$items[$i];
			}
		}
		else{
			for(my $i=1;$i<@items; $i++){
				my $col_name=$col_names{$i};
				my $row_name=$items[0];
				$statistics{$col_name}->{$row_name}=$items[$i];
				#print "$col_name:$row_name\n";
			}
		}
		$n++;
	}
	close($STA);
	
	return(\%statistics);
}
#####################################
sub S5_Colocalization_counting{
	my ($variables_pointer, $export_format_pointer)=@_;
	my %variables=%$variables_pointer;
	my %export_format=%$export_format_pointer;
	
	#read statistics file
	my $statistics_pointer=read_statistics($variables{statistics_file});
	my %statistics=%$statistics_pointer; 
	
	#Read genomic ligations
	print "\tRead genomic ligations:\n";
	my (%counting, %seq_hash, %view_points);
	my @sample_names=split(",", $variables{sample_names});
	foreach my $sample_name(@sample_names){#2
		print "\t\t$export_format{thread}: $sample_name\n";
		foreach my $file_tail('cis_ligation', 'trans_ligation'){#3
			open my($IN), "<", $variables{ligations_dir}.'/'.$sample_name.'.'.$file_tail or die;
			while (my $line=<$IN>){#4
				chomp($line);
				my @items=split("\t", $line);
				push(@items, $sample_name);
				my $row=$items[$export_format{row_1}].'_'.$items[$export_format{row_2}];
				my $r_row=$items[$export_format{row_1}].'_'.$items[$export_format{row_2r}];
				my $col=$items[$export_format{col_1}].'_'.$items[$export_format{col_2}];
				print "$sample_name:$file_tail:$col:@items\n" unless @items==13;
				my $row_col=$row.'_'.$col;
				my $r_row_col=$r_row.'_'.$col;
				#numeric sequence
				my $seq=sub_3C::numeric_seq($items[9]);#numeric sequence
				#my $seq=$items[9];
				if($row eq $r_row){
					#update $row_name
					$counting{$row}->{$col}->{total}++;
					unless(exists $seq_hash{$seq}){
						$counting{$row}->{$col}->{unique}++;
						$seq_hash{$seq}=1;
					}
				}
				#row_2 cmp row_2r
				else{
					#update $row_name
					$counting{$row}->{$col}->{total} += 0.5;
					$counting{$r_row}->{$col}->{total} += 0.5;
					unless(exists $seq_hash{$seq}){
						$counting{$row}->{$col}->{unique} += 0.5;
						$counting{$r_row}->{$col}->{unique} += 0.5;
						$seq_hash{$seq}=1;
					}
				}
				unless (exists $view_points{$col}){
					$view_points{$col}->{view_point_name}=$items[0];#display_id
					$view_points{$col}->{size_bin_name}=$items[1];
					$view_points{$col}->{enzyme_bin_name}=$items[2];
					$view_points{$col}->{sample_name}=$items[12];
				}
			}#4 
			close($IN);
		}#3
	}#2
	%seq_hash=(); 
	my @row_names=keys %counting;
	@row_names=sort @row_names;
	my @col_names=map { keys $counting{$_} } @row_names;
	@col_names=List::MoreUtils::uniq @col_names;
	@col_names= sort @col_names;
	
	#group @col_names into parts for minimizing the file size
	my (@col_groups, @tmp);
	my $n=1;
	foreach my $col_name(@col_names){
		push(@tmp, $col_name);
		if(int($n/100)==$n/100 or $col_name eq $col_names[-1]){
			push(@col_groups, [@tmp]);
			undef @tmp;
		}
		$n++;
	}
	
	
	print "\t and then export genomic interactions by read counts into files.\n\n";
	mkdir($variables{ligation_frequency_dir}, 0755) unless -d $variables{ligation_frequency_dir};
	for(my $i=0; $i<@col_groups; $i++){#2
		my $pointer=$col_groups[$i];
		my @sub_col_names=@$pointer;
		#exports
		my $file_head=$variables{ligation_frequency_dir}.'/'.$export_format{row_attr}.'_'.$export_format{col_attr}.'_'.$i;
		open my($TOT), ">", $file_head.'_Total_RC.csv';
		open my($UNI), ">", $file_head.'_Unique_RC.csv';
		open my($norm_TOT), ">", $file_head.'_NormTotal_RC.csv';
		open my($norm_UNI), ">", $file_head.'_NormUnique_RC.csv';
		print $TOT join(',', 'chr', 'bin', @sub_col_names), "\n";
		print $UNI join(',', 'chr', 'bin', @sub_col_names), "\n";
		print $norm_TOT join(',', 'chr', 'bin', @sub_col_names), "\n";
		print $norm_UNI join(',', 'chr', 'bin', @sub_col_names), "\n";
	
		foreach my $row(@row_names){#3 one class, one file
			my($chr, $bin)=split('_', $row);
			#count the number of reads
			my (@total_counting, @unique_counting, @norm_total_counting, @norm_unique_counting);
			foreach my $col(@sub_col_names){#4
				#RC
				my $total_counts=(exists $counting{$row}->{$col}->{total}) ? $counting{$row}->{$col}->{total} : 0;
				my $unique_counts=(exists $counting{$row}->{$col}->{unique}) ? $counting{$row}->{$col}->{unique} : 0; 
				#normalization by library size
				my $sample_name=$view_points{$col}->{sample_name};
				my $raw_pairs=$statistics{$sample_name}->{raw_pairs};
				my $norm_total_counts=int(($total_counts*1e9)/$raw_pairs + 0.5);
				my $norm_unique_counts=int(($unique_counts*1e9)/$raw_pairs + 0.5);
				#rand number
				srand;
				my $rand_RC;
				while(1){
					#generate random number 0-1 default
					$rand_RC= int(rand($variables{RC_noise})*100)/100;
					last if $rand_RC>0;
				}
				$total_counts=$rand_RC if $total_counts==0;
				$unique_counts=$rand_RC if $unique_counts==0;
				$norm_total_counts=$rand_RC if $norm_total_counts==0;
				$norm_unique_counts=$rand_RC if $norm_unique_counts==0;
				push(@total_counting, $total_counts);
				push(@unique_counting, $unique_counts);
				push(@norm_total_counting, $norm_total_counts);
				push(@norm_unique_counting, $norm_unique_counts);
			}#4
		
			print $TOT join(',', $chr, $bin, @total_counting), "\n";
			print $UNI join(',', $chr, $bin, @unique_counting), "\n";
			print $norm_TOT join(',', $chr, $bin, @norm_total_counting), "\n";
			print $norm_UNI join(',', $chr, $bin, @norm_unique_counting), "\n";
		}#3
		close($TOT);
		close($UNI);
		close($norm_TOT);
		close($norm_UNI);
	}#2
}
################################################
sub S6_Tscore_calculation{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	if($variables{'3C_Tscore_calculation'} eq 'yes'){#2
		print "###Tscore calculation.###\n";
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_ligation_frequency_dir" , $variables{ligation_frequency_dir});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_viewpoints_csv" , $variables{viewpoints_csv});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_size_bin_csv" , $variables{size_bin_csv});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_enzyme_bin_csv" , $variables{enzyme_bin_csv});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_p_level" , $variables{RC_cumulative_p});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_RC_noise" , $variables{RC_noise});
		
		my $file_names_pointer=sub_3C::files_list($variables{ligation_frequency_dir}, 'file_name');
		my @file_names=@$file_names_pointer;
		my @RC_csv_names=grep(/_RC.csv$/, @file_names);
		#run
		foreach my $RC_csv_name(@RC_csv_names){#3
			sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_RC_csv_name" , $RC_csv_name);
			#run R script
			print "\tRead $RC_csv_name, and run $variables{Tscore_R_file}\n";
			system("Rscript $variables{Tscore_R_file}");
		}#3
	
	}#2

}
#######################################################
#4C-seq
sub C4_S1_colocalization_detection{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_name=$variables{sample_name};
	my $enzyme_cutter=length($variables{enzyme_site});
	#
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @fastq_files=split(',', $sample_info{$sample_name}->{raw_files});
	#
	my $site_info_pointer=$variables{site_info_pointer};
	my %site_info=%$site_info_pointer;
	my $view_seq=$site_info{$sample_name}->{view_seq};
	my $region_seq=$site_info{$sample_name}->{region_seq};
	
	#
	my (%pair, %seq_index);
	my $raw_reads=1;
	my $query_reads=0;
	my $original_regions=0;
	my $fasta_file=$variables{sample_dir}.'/'.$sample_name.'.fa';
	open my ($OUT), ">", $fasta_file or die;
	foreach my $fastq_file(@fastq_files){#2
		print "\tRead $fastq_file\n";
		open my($IN), "<", $fastq_file or die;
		my ($name_line, $seq_line, $third_line, $qual_line);
		while (defined($name_line = <$IN>) && defined($seq_line = <$IN>) && defined($third_line = <$IN>) && defined($qual_line = <$IN>) ){#3
			chomp($name_line, $seq_line, $third_line, $qual_line );
			my $last_pos=length($seq_line);
			my $head_seq='NA';
			my $tail_seq='NA';
			#find view seq
			if($region_seq=~/$seq_line/){#4
				$original_regions++;
			}
			else{
				while($last_pos>0){#5
					my $tmp_seq=substr($seq_line, 0, $last_pos);
					$last_pos=rindex($tmp_seq, $variables{enzyme_site});
					if($last_pos>0){#6
						$head_seq=substr($tmp_seq, 0, $last_pos+$enzyme_cutter); #with enzyme site
						if( $view_seq=~/$head_seq$/ and ($last_pos+2)>=$variables{min_trim_len}){
							$tail_seq=substr($tmp_seq, $last_pos);#without enzyme site
							#print "$last_pos", substr($view_seq, -20), ":", $head_seq, "\n";
							$last_pos=-1;
						}
					}#6
				}#5
			}#4
			#export
			if(length($tail_seq)>=$variables{query_length}){#4
				if(exists $seq_index{$tail_seq}){
					$seq_index{$tail_seq} .= ','.$raw_reads;
				}
				else{
					$seq_index{$tail_seq}=$raw_reads;
					print $OUT ">$raw_reads\n", "$tail_seq\n";
					$query_reads++;
					$pair{$raw_reads}->{view_displayid}=$sample_name; #for 4C
					$pair{$raw_reads}->{view_ref}=$site_info{$sample_name}->{chr};
					$pair{$raw_reads}->{view_offset}=$site_info{$sample_name}->{site_position};
					$pair{$raw_reads}->{view_seq}='*';
					$pair{$raw_reads}->{co_seq}=$tail_seq;
				}
			}#4
			$raw_reads++;
		}#3
		close($IN);
	}#2
	close($OUT);
	sub_3C::refresh_log($variables{sample_log}, "raw_reads", $raw_reads);
	sub_3C::refresh_log($variables{sample_log}, "query_reads", $query_reads);
	sub_3C::refresh_log($variables{sample_log}, "original_regions", $original_regions);
	
	print "###genome mapping and co-localization\n";
	print "###one sequence alignment\n";
	my $sam_file=$variables{sample_dir}.'/'.$sample_name.".sam";
	my $co_pair_pointer=sub_3C::one_seq_alignment(\%variables, \%pair, $sample_name, $fasta_file, $sam_file, '4C');
	my %co_pair=%$co_pair_pointer;
	print "###export co-localization\n";
	my $cis_ligation=0;
	my $trans_ligation=0;
	my $unmapped_ligation=0;
	my $multiple_ligation=0;
	open my($cis_SITE), ">", $variables{sample_dir}.'/'.$sample_name.'.cis_ligation' or die; # $found=2
	open my($trans_SITE), ">", $variables{sample_dir}.'/'.$sample_name.'.trans_ligation' or die; #$found=4
	open my($un_SITE), ">", $variables{sample_dir}.'/'.$sample_name.'.unmapped_ligation' or die; # $found=8
	open my($multiple_SITE), ">", $variables{sample_dir}.'/'.$sample_name.'.multiple_ligation' or die; # $found=16
	foreach my $query_seq(keys %seq_index){#2
		my @raw_reads_arr=split(',', $seq_index{$query_seq});
		my $coordinate=$raw_reads_arr[0];
		
		foreach my $equal_raw_reads(@raw_reads_arr){#3
			#replace first_raw_reads with real raw_reads for unmapped/multiple_ligation
			my @out_arr=split("\t", $co_pair{$coordinate}->{out_line});
			#print "$coordinate:$out_arr[1]:$equal_raw_reads\n";
			$out_arr[1]=$equal_raw_reads;
			#print "$out_arr[1]\n";
			my $out_line=join("\t", @out_arr);
			#export co-seq
			my $found=$co_pair{$coordinate}->{valid_seq};
			if ($found==2)		{	
				print $cis_SITE "$out_line\n";
				$cis_ligation++;
			}
			elsif ($found==4)	{	
				print $trans_SITE "$out_line\n";	
				$trans_ligation++;
			}
			elsif ($found==8)	{	
				print $un_SITE "$out_line\n";
				$unmapped_ligation++;
			}
			elsif ($found==16)	{	
				print $multiple_SITE "$out_line\n";
				$multiple_ligation++;
			}
		}#3
	}#2
	close($cis_SITE);
	close($trans_SITE);
	close($un_SITE);
	close($multiple_SITE);
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation", $cis_ligation);
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation", $trans_ligation);
	sub_3C::refresh_log($variables{sample_log}, "unmapped_ligation", $unmapped_ligation);
	sub_3C::refresh_log($variables{sample_log}, "multiple_ligation", $multiple_ligation);
}

##############################################################################
#export: statistics.txt 
sub C4_S2_result_statistics{#1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});

	mkdir($variables{ligations_dir}, 0755) unless -d $variables{ligations_dir};
	#generate statistics result
	my %counting_hash;
	my @combined_file_tails=('cis_ligation', 'trans_ligation');
	foreach my $sample_name(@sample_names){#2
		#read sample log
		my $sample_dir=$variables{result_dir}.'/'.$sample_name;
		my $sample_log=$sample_dir.'/'.$sample_name.'.log';
		open my($IN), "<", $sample_log or die;
		while(<$IN>){#3
			chomp($_);
			my($key, $value)=split("=", $_);
			$counting_hash{$key}->{$sample_name}=$value if length($key)>=1 and length($value)>=1;
		}#3
		#copy *.cis_ligation and *.trans_ligation
		foreach my $file_tail(@combined_file_tails){
			my $in_file=$sample_dir.'/'.$sample_name.'.'.$file_tail;
			my $out_file=$variables{ligations_dir}.'/'.$sample_name.'.'.$file_tail;
			system("cp $in_file $out_file");
		}
	}#2
	
	#statistics
	open my($OUT), ">", $variables{statistics_file} or die;
	print $OUT join("\t", 'item_names', @sample_names), "\n";
	foreach my $item_name(sort (keys %counting_hash) ){
		my @counting_num=map { $counting_hash{$item_name}->{$_} } @sample_names;
		print $OUT join("\t", $item_name, @counting_num), "\n";
	}
	close($OUT);
}
#####################################
sub C4_S3_colocalization_counting{
	my ($variables_pointer, $export_format_pointer)=@_;
	my %variables=%$variables_pointer;
	my %export_format=%$export_format_pointer;
	
	#read statistics file
	my $statistics_pointer=read_statistics($variables{statistics_file});
	my %statistics=%$statistics_pointer; 
	
	#Read genomic ligations
	print "\tRead genomic ligations:\n";
	my (%counting, %seq_hash, %view_points);
	my @sample_names=split(",", $variables{sample_names});
	foreach my $sample_name(@sample_names){#2
		print "\t\t$export_format{thread}: $sample_name\n";
		foreach my $file_tail('cis_ligation', 'trans_ligation'){#3
			open my($IN), "<", $variables{ligations_dir}.'/'.$sample_name.'.'.$file_tail or die;
			while (my $line=<$IN>){#3
				chomp($line);
				my @items=split("\t", $line);
				push(@items, $sample_name);
				my $row=$items[$export_format{row_1}].'_'.$items[$export_format{row_2}];
				my $r_row=$items[$export_format{row_1}].'_'.$items[$export_format{row_2r}];
				#print "$sample_name:$file_tail:$r_row\n";
				my $col=$items[$export_format{col}];
				my $row_col=$row.'_'.$col;
				my $r_row_col=$r_row.'_'.$col;
			
				#numeric sequence
				my $seq=sub_3C::numeric_seq($items[10]);#numeric sequence
				#my $seq=$items[10];
				if($row eq $r_row){
					#update $row_name
					$counting{$row}->{$col}->{total}++;
					unless(exists $seq_hash{$seq}){
						$counting{$row}->{$col}->{unique}++;
						$seq_hash{$seq}=1;
					}
				}
				#row_2 cmp row_2r
				else{
					#update $row_name
					$counting{$row}->{$col}->{total} += 0.5;
					$counting{$r_row}->{$col}->{total} += 0.5;
					unless(exists $seq_hash{$seq}){
						$counting{$row}->{$col}->{unique} += 0.5;
						$counting{$r_row}->{$col}->{unique} += 0.5;
						$seq_hash{$seq}=1;
					}
				}
				unless (exists $view_points{$col}){
					$view_points{$col}->{view_point_name}=$items[0];
					$view_points{$col}->{size_bin_name}=$items[2];
					$view_points{$col}->{enzyme_bin_name}=$items[3];
					$view_points{$col}->{sample_name}=$items[0];
				}
			}#4 
			close($IN);
		}#3
	}#2
	%seq_hash=(); 
	my @row_names=keys %counting;
	@row_names= sort @row_names;
	my @col_names=map { keys $counting{$_} } @row_names;
	@col_names=List::MoreUtils::uniq @col_names;
	@col_names= sort @col_names;
	
	print "\t and then export genomic interactions by read counts into files.\n\n";
	mkdir($variables{ligation_frequency_dir}, 0755) unless -d $variables{ligation_frequency_dir};
	my $file_head=$variables{ligation_frequency_dir}.'/'.$export_format{row_attr}.'_'.$export_format{col_attr};
	open my($TOT), ">", $file_head.'_Total_RC.csv';
	open my($UNI), ">", $file_head.'_Unique_RC.csv';
	open my($norm_TOT), ">", $file_head.'_NormTotal_RC.csv';
	open my($norm_UNI), ">", $file_head.'_NormUnique_RC.csv';
	print $TOT join(',', 'chr', 'bin', @col_names), "\n";
	print $UNI join(',', 'chr', 'bin', @col_names), "\n";
	print $norm_TOT join(',', 'chr', 'bin', @col_names), "\n";
	print $norm_UNI join(',', 'chr', 'bin', @col_names), "\n";
	
	foreach my $row(@row_names){#2 one class, one file
		my($chr, $bin)=split('_', $row);
		#count the number of reads
		my (@total_counting, @unique_counting, @norm_total_counting, @norm_unique_counting);
		foreach my $col(@col_names){#3
			#RC
			my $total_counts=(exists $counting{$row}->{$col}->{total}) ? $counting{$row}->{$col}->{total} : 0;
			my $unique_counts=(exists $counting{$row}->{$col}->{unique}) ? $counting{$row}->{$col}->{unique} : 0; 
			#normalization by library size
			my $sample_name=$view_points{$col}->{sample_name};
			my $raw_reads=$statistics{$sample_name}->{raw_reads};
			my $norm_total_counts=int(($total_counts*1e9)/$raw_reads + 0.5);
			my $norm_unique_counts=int(($unique_counts*1e9)/$raw_reads + 0.5);
			#rand number
			srand;
			my $rand_RC;
			while(1){
				#generate random number 0-1 default
				$rand_RC= int(rand($variables{RC_noise})*100)/100;
				last if $rand_RC>0;
			}
			$total_counts=$rand_RC if $total_counts==0;
			$unique_counts=$rand_RC if $unique_counts==0;
			$norm_total_counts=$rand_RC if $norm_total_counts==0;
			$norm_unique_counts=$rand_RC if $norm_unique_counts==0;
			push(@total_counting, $total_counts);
			push(@unique_counting, $unique_counts);
			push(@norm_total_counting, $norm_total_counts);
			push(@norm_unique_counting, $norm_unique_counts);
		}#3
		
		#export
		#print "$row:$a=$b=$c=$d=$e\n";
		print $TOT join(',', $chr, $bin, @total_counting), "\n";
		print $UNI join(',', $chr, $bin, @unique_counting), "\n";
		print $norm_TOT join(',', $chr, $bin, @norm_total_counting), "\n";
		print $norm_UNI join(',', $chr, $bin, @norm_unique_counting), "\n";
	}#2
	close($TOT);
	close($UNI);
	close($norm_TOT);
	close($norm_UNI);
}
################################################
sub C4_S4_Tscore_calculation{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	if($variables{'4C_Tscore_calculation'} eq 'yes'){#2
		print "###Tscore calculation.###\n";
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_ligation_frequency_dir" , $variables{ligation_frequency_dir});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_viewpoints_csv" , $variables{viewpoints_csv});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_size_bin_csv" , $variables{size_bin_csv});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_enzyme_bin_csv" , $variables{enzyme_bin_csv});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_p_level" , $variables{RC_cumulative_p});
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_RC_noise" , $variables{RC_noise});
		
		my $file_names_pointer=sub_3C::files_list($variables{ligation_frequency_dir}, 'file_name');
		my @file_names=@$file_names_pointer;
		my @RC_csv_names=grep(/_RC.csv$/, @file_names);
		#run
		foreach my $RC_csv_name(@RC_csv_names){#3
			sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_RC_csv_name" , $RC_csv_name);
			#run R script
			print "\tRead $RC_csv_name, and run $variables{Tscore_R_file}\n";
			system("Rscript $variables{Tscore_R_file}");
		}#3
	
	}#2

}

##################################################
sub HiC_S1_colocalization_detection{
	my ($variables_pointer, $R1_fastq, $name)=@_;
	my %variables=%$variables_pointer;
	my $enzyme_sites_pointer=$variables{enzyme_sites_pointer};
	my %enzyme_sites=%$enzyme_sites_pointer;
	my $R2_fastq=$R1_fastq;
	$R2_fastq=~s/_R1_/_R2_/;

	my $SAM_file=$variables{sample_dir}.'/'.$name.".sam";
	#sequences alignment
	unless(-e $SAM_file){
		print" \n###Pair-end fastq files alignment  of $name begin:\n";
		print "$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $SAM_file \n";
		system("$variables{bowtie_options} -1 $R1_fastq -2 $R2_fastq -S $SAM_file ");
		print "\nSequence alignment output: $SAM_file!\n";
	}
	
	print "read_alignment result from $SAM_file\n";
	my $original_regions=0;
	my $cis_ligation=0;
	my $trans_ligation=0;
	my $one_alignment=0;
	my $no_alignment=0;
	my $multiple_alignment=0;
	my $unknown_mapped_pairs=0;
	my $both_alignment=0;
	my $raw_pairs=0;
	open my($ORI), ">", $variables{sample_dir}.'/'.$name.'.original_regions' or die;
	open my($cis_SITE), ">", $variables{sample_dir}.'/'.$name.'.cis_ligation' or die;
	open my($trans_SITE), ">", $variables{sample_dir}.'/'.$name.'.trans_ligation' or die;
	open my($ONE), ">", $variables{sample_dir}.'/'.$name.'.one_alignment' or die;
	open my($NO), ">", $variables{sample_dir}.'/'.$name.'.no_alignment' or die;
	open my($MULTI), ">", $variables{sample_dir}.'/'.$name.'.multiple_alignment' or die;
	open my($IN), "<", $SAM_file or die;
	my ($line_1, $line_2);
	while( defined ($line_1 = <$IN>) && defined ($line_2 = <$IN>) ){ #2 read SAM
		chomp($line_1, $line_2);
		#read two lines per time
		my $pair_pointer=sub_3C::read_paired_end(\%variables, $line_1, $line_2);#subroutine
		my %pair=%$pair_pointer;
		my $output=join("\t", $pair{coordinate}, 
				$pair{ref_1}, $pair{offset_1}, $pair{seq_1}, $pair{alignment_1},
				$pair{ref_2}, $pair{offset_2}, $pair{seq_2}, $pair{alignment_2},
			);
		#print "$pair{valid_pair}\n";
		#judging
		if ($pair{valid_pair} eq 'one_unique'){#3
			print $ONE $output, "\n";
			$one_alignment++;
		}#3
		elsif ($pair{valid_pair}=~/multiple/){	#3
			print $MULTI "$output\t", "$pair{valid_pair}\n";
			$multiple_alignment++;
		}#3
		elsif ($pair{valid_pair} eq 'unaligned'){#3 #unaligned
			print $NO "$pair{coordinate}\t", "$pair{seq_1}\t", "$pair{seq_2}\n";
			$no_alignment++;
		}#3
		else{#3#cis_ligation or trans_ligation
			my $ref=$pair{ref_1};
			my $offset=$pair{offset_1};
			my $end_offset=$offset+length($pair{seq_1});
			my $ref1_size_bin_name=int($offset/$variables{size_bin_len})+1;
			my $ref1_right_size_bin_name=int($end_offset/$variables{size_bin_len})+1;
			my $ref1_enzyme_bin_name=sub_3C::enzyme_bin_localization(\%variables, $ref, $offset);
			my $ref1_right_enzyme_bin_name=sub_3C::enzyme_bin_localization(\%variables, $ref, $end_offset);
			$ref=$pair{ref_2};
			$offset=$pair{offset_2};
			$end_offset=$offset+length($pair{seq_2});
			my $ref2_size_bin_name=int($offset/$variables{size_bin_len})+1;
			my $ref2_right_size_bin_name=int($end_offset/$variables{size_bin_len})+1;
			my $ref2_enzyme_bin_name=sub_3C::enzyme_bin_localization(\%variables, $ref, $offset);
			my $ref2_right_enzyme_bin_name=sub_3C::enzyme_bin_localization(\%variables, $ref, $end_offset);
			
			$output=join("\t", $pair{coordinate}, $pair{ref_1}, $ref1_size_bin_name, $ref1_right_size_bin_name, 
						$ref1_enzyme_bin_name, $ref1_right_enzyme_bin_name, $pair{seq_1},
						$pair{ref_2}, $ref2_size_bin_name, $ref2_right_size_bin_name, $ref2_enzyme_bin_name, 
						$ref2_right_enzyme_bin_name, $pair{seq_2}, 'both_unique');
					
			if($pair{valid_pair} eq 'cis_pair'){#4
				print $ORI $output, "\n";
				$original_regions++;
			}#4
			elsif($pair{valid_pair} eq 'cis_unpair'){#4
				print $cis_SITE $output, "\n";
				#print $output, "\n";
				$cis_ligation++;
			}#4
			else{#4   $pair{valid_pair} eq 'trans_unpair'
				print $trans_SITE $output, "\n";
				$trans_ligation++;
			}#4
			$both_alignment++;
		}#3#cis_ligation or trans_ligation
		$raw_pairs++;
	} #2 read SAM
	close($IN);
	close($ORI);
	close($cis_SITE);
	close($trans_SITE);
	close($ONE);
	close($NO);
	close($MULTI);
	
	sub_3C::refresh_log($variables{sample_log}, "original_regions_num", $original_regions);
	sub_3C::refresh_log($variables{sample_log}, "cis_ligation_num", $cis_ligation);
	sub_3C::refresh_log($variables{sample_log}, "trans_ligation_num", $trans_ligation);
	sub_3C::refresh_log($variables{sample_log}, "one_alignment_num", $one_alignment);
	sub_3C::refresh_log($variables{sample_log}, "no_alignment_num", $no_alignment);
	sub_3C::refresh_log($variables{sample_log}, "multiple_alignment_num", $multiple_alignment);
	sub_3C::refresh_log($variables{sample_log}, "both_alignment_num", $both_alignment);
	sub_3C::refresh_log($variables{sample_log}, "raw_pairs_num", $raw_pairs);
}
###################################################
sub HiC_S2_statistics{#1
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=keys %sample_info;
	@sample_names=sort @sample_names;
	
	mkdir($variables{ligations_dir}, 0755) unless -d $variables{ligations_dir};
	#generate statistics result
	my %counting_hash;
	my @combined_file_tails=('cis_ligation', 'trans_ligation', 'original_regions');
	foreach my $sample_name(@sample_names){#2
		#read sample log
		my $sample_dir=$variables{result_dir}.'/'.$sample_name;
		my $sample_log=$variables{log_dir}.'/'.$sample_name.'.log';
		open my($IN), "<", $sample_log or die;
		while(<$IN>){#3
			chomp($_);
			my($key, $value)=split("=", $_);
			$counting_hash{$key}->{$sample_name}=$value;
		}#3
		$counting_hash{raw_files}->{$sample_name}=split(',', $sample_info{$sample_name}->{'raw_files'});
		
		#combine and export files from %combined_file_tails
		my @file_heads=split(',', $sample_info{$sample_name}->{'R1_file_names_head'});
		foreach my $file_tail(@combined_file_tails){#3
			my @sub_files=map {$sample_dir.'/'.$_.'.'.$file_tail} @file_heads;
			my $in_sub_files_str=join (' ', @sub_files);
			my $out_file=$variables{ligations_dir}.'/'.$sample_name.'.'.$file_tail;
			system("cat $in_sub_files_str > $out_file"); #combine
		}#3
	}#2
	
	#statistics
	open my($OUT), ">", $variables{statistics_file} or die;
	print $OUT join("\t", 'item_names', @sample_names), "\n";
	foreach my $item_name(sort (keys %counting_hash) ){
		my @counting_num=map { $counting_hash{$item_name}->{$_} } @sample_names;
		print $OUT join("\t", $item_name, @counting_num), "\n";
	}
	close($OUT);
}

#######################################
#	print "###export counting table by viewpoints index.\n";
	#0=displayid, 1=coordinate, 2=size_bin_no, 3=enzyme_bin_no, 4=co_size_bin_no, 5=right_co_size_bin_no, 
	#6=co_enzyme_bin_no, 7=right_co_enzyme_bin_no,
	#8=co_ref, 9=col_offset, 10=col_seq, 11=tag, 12=sample_name
	#12:sample_name is not in the file but will be added in function S5_Co_localization_counting 
	#col_1 always sample names

############
sub HiC_S3_cisRC_counting{
	my ($variables_pointer, $ligation_file)=@_;
	my %variables=%$variables_pointer;
	
	#read statistics file
	my $statistics_pointer=read_statistics($variables{statistics_file});
	my %statistics=%$statistics_pointer; 
	
	print "\tRead $ligation_file\n";
	my (%size_counting, %enzyme_counting);
	open my($IN), "<", $ligation_file or die;
	while(my $line=<$IN>){
		chomp($line);
		my ($coordinate, $ref_1, $ref1_size_bin_name, $ref1_right_size_bin_name, $ref1_enzyme_bin_name, 
			$ref1_right_enzyme_bin_name, $seq_1, $ref_2, $ref2_size_bin_name, $ref2_right_size_bin_name, 
			$ref2_enzyme_bin_name, $ref2_right_enzyme_bin_name, $seq_2, $alignment_type)=split("\t", $line);
		#setup index of co-localizations
		my (@size_index, @enzyme_index, $index);
		foreach my $a($ref1_size_bin_name, $ref1_right_size_bin_name){
			foreach my $b($ref2_size_bin_name, $ref2_right_size_bin_name){
				$index=($a<=$b) ? join(',', $ref_1, $a, $b) : join(',', $ref_1, $b, $a);
				push(@size_index, $index) unless List::Util::first {$_ eq $index} @size_index;
			}
		}
		foreach my $a($ref1_enzyme_bin_name, $ref1_right_enzyme_bin_name){
			foreach my $b($ref2_enzyme_bin_name, $ref2_right_enzyme_bin_name){
				$index=($a<=$b) ? join(',', $ref_1, $a, $b) : join(',', $ref_1, $b, $a);
				push(@enzyme_index, $index) unless List::Util::first {$_ eq $index} @enzyme_index;
			}
		}
		#assign RC
		my $RC=int(10/@size_index)/10;
		foreach my $si(@size_index){
			my ($a, $b, $c)=split(',', $si);
			if (exists $size_counting{$a}->{$b}->{$c}){ 
				$size_counting{$a}->{$b}->{$c} += $RC;
			}
			else{ $size_counting{$a}->{$b}->{$c} = $RC; }
		}
		$RC=int(10/@enzyme_index)/10;
		foreach my $ei(@enzyme_index){
			my ($a, $b, $c)=split(',', $ei);
			if (exists $enzyme_counting{$a}->{$b}->{$c}){ 
				$enzyme_counting{$a}->{$b}->{$c} += $RC;
			}
			else{ $enzyme_counting{$a}->{$b}->{$c} = $RC; }
		}
	}
	close($IN);
	
	my @a=split('/', $ligation_file);
	my $name=$a[-1];
	$name=~s/\./_/;
	my $csv_header=join(',', 'chr', 'bin1', 'bin2', 'dist', 'RC');
	#export size_bin
	my $out_file=$variables{ligation_frequency_dir}.'/Size_bin_'.$name.'_RC.csv';
	print "\t\texport RCs into $out_file\n";
	open my($OUT1), ">", $out_file or die;
	print $OUT1 "$csv_header\n";
	foreach my $ref1(sort {$a<=>$b} (keys %size_counting)){
		my $pointer1=$size_counting{$ref1};
		my %hash1=%$pointer1;
		foreach my $pos1(sort {$a<=>$b} (keys %hash1)){
			my $pointer2=$hash1{$pos1};
			my %hash2=%$pointer2;
			foreach my $pos2(sort {$a<=>$b} (keys %hash2)){
				my $RC=$size_counting{$ref1}->{$pos1}->{$pos2};
				my $dist=$pos2-$pos1;
				print $OUT1 join(',', $ref1, $pos1, $pos2, $dist, $RC), "\n";
			}
		}
	}
	close($OUT1);
	#export enzyme_bin
	$out_file=$variables{ligation_frequency_dir}.'/Enzyme_bin_'.$name.'_RC.csv';
	print "\t\texport RCs into $out_file\n";
	open my($OUT2), ">", $out_file or die;
	print $OUT2 "$csv_header\n";
	foreach my $ref1(sort (keys %enzyme_counting)){
		my $pointer1=$enzyme_counting{$ref1};
		my %hash1=%$pointer1;
		foreach my $pos1(sort (keys %hash1)){
			my $pointer2=$hash1{$pos1};
			my %hash2=%$pointer2;
			foreach my $pos2(sort (keys %hash2)){
				my $RC=$enzyme_counting{$ref1}->{$pos1}->{$pos2};
				my $dist=$pos2-$pos1;
				print $OUT2 join(',', $ref1, $pos1, $pos2, $dist, $RC), "\n";
			}
		}
	}
	close($OUT2);
	
	#
}
########################
sub HiC_S3_transRC_counting{
	my ($variables_pointer, $ligation_file)=@_;
	my %variables=%$variables_pointer;

	#read statistics file
	my $statistics_pointer=read_statistics($variables{statistics_file});
	my %statistics=%$statistics_pointer; 
	
	print "\tRead $ligation_file\n";
	my (%size_counting, %enzyme_counting);
	open my($IN), "<", $ligation_file or die;
	while(my $line=<$IN>){
		chomp($line);
		my ($coordinate, $ref_1, $ref1_size_bin_name, $ref1_right_size_bin_name, $ref1_enzyme_bin_name, 
			$ref1_right_enzyme_bin_name, $seq_1, $ref_2, $ref2_size_bin_name, $ref2_right_size_bin_name, 
			$ref2_enzyme_bin_name, $ref2_right_enzyme_bin_name, $seq_2, $alignment_type)=split("\t", $line);
		#setup index of co-localizations
		my (@size_index, @enzyme_index, $index);
		if($ref_1<$ref_2){
			foreach my $a($ref1_size_bin_name, $ref1_right_size_bin_name){
				foreach my $b($ref2_size_bin_name, $ref2_right_size_bin_name){
					$index=join(',', $ref_1, $ref_2, $a, $b);
					push(@size_index, $index) unless List::Util::first {$_ eq $index} @size_index;
				}
			}
			foreach my $a($ref1_enzyme_bin_name, $ref1_right_enzyme_bin_name){
				foreach my $b($ref2_enzyme_bin_name, $ref2_right_enzyme_bin_name){
					$index=join(',', $ref_1, $ref_2, $a, $b);
					push(@enzyme_index, $index) unless List::Util::first {$_ eq $index} @enzyme_index;
				}
			}
		}
		else{
			foreach my $a($ref2_size_bin_name, $ref2_right_size_bin_name){
				foreach my $b($ref1_size_bin_name, $ref1_right_size_bin_name){
					$index=join(',', $ref_2, $ref_1, $a, $b);
					push(@size_index, $index) unless List::Util::first {$_ eq $index} @size_index;
				}
			}
			foreach my $a($ref2_enzyme_bin_name, $ref2_right_enzyme_bin_name){
				foreach my $b($ref1_enzyme_bin_name, $ref1_right_enzyme_bin_name){
					$index=join(',', $ref_2, $ref_1, $a, $b);
					push(@enzyme_index, $index) unless List::Util::first {$_ eq $index} @enzyme_index;
				}
			}
		}
		#assign RC
		my $RC=int(10/@size_index)/10;
		foreach my $si(@size_index){
			my ($a, $b, $c,$d)=split(',', $si);
			if (exists $size_counting{$a}->{$b}->{$c}->{$d}){ 
				$size_counting{$a}->{$b}->{$c}->{$d} += $RC;
			}
			else{ $size_counting{$a}->{$b}->{$c}->{$d} = $RC; }
		}
		$RC=int(10/@enzyme_index)/10;
		foreach my $ei(@enzyme_index){
			my ($a, $b, $c,$d)=split(',', $ei);
			if (exists $enzyme_counting{$a}->{$b}->{$c}->{$d}){ 
				$enzyme_counting{$a}->{$b}->{$c}->{$d} += $RC;
			}
			else{ $enzyme_counting{$a}->{$b}->{$c}->{$d} = $RC; }
		}
	}
	close($IN);
	
	my @a=split('/', $ligation_file);
	my $name=$a[-1];
	$name=~s/\./_/;
	my $csv_header=join(',', 'chr1', 'bin1', 'chr2', 'bin2', 'RC');
	#export size_bin
	my $out_file=$variables{ligation_frequency_dir}.'/Size_bin_'.$name.'_RC.csv';
	print "\t\texport RCs into $out_file\n";
	open my($OUT1), ">", $out_file or die;
	print $OUT1 "$csv_header\n";
	foreach my $ref1(sort {$a<=>$b} (keys %size_counting)){
		my $pointer1=$size_counting{$ref1};
		my %hash1=%$pointer1;
		foreach my $ref2(sort {$a<=>$b} (keys %hash1)){
			my $pointer2=$hash1{$ref2};
			my %hash2=%$pointer2;
			foreach my $pos1(sort {$a<=>$b} (keys %hash2)){
				my $pointer3=$hash2{$pos1};
				my %hash3=%$pointer3;
				foreach my $pos2(sort {$a<=>$b} (keys %hash3)){
					my $RC=$size_counting{$ref1}->{$ref2}->{$pos1}->{$pos2};
					print $OUT1 join(',', $ref1, $pos1, $ref2, $pos2, $RC), "\n";
				}
			}
		}
	}
	close($OUT1);
	#export enzyme_bin
	$out_file=$variables{ligation_frequency_dir}.'/Enzyme_bin_'.$name.'_RC.csv';
	print "\t\texport RCs into $out_file\n";
	open my($OUT2), ">", $out_file or die;
	print $OUT2 "$csv_header\n";
	foreach my $ref1(sort (keys %enzyme_counting)){
		my $pointer1=$enzyme_counting{$ref1};
		my %hash1=%$pointer1;
		foreach my $ref2(sort (keys %hash1)){
			my $pointer2=$hash1{$ref2};
			my %hash2=%$pointer2;
			foreach my $pos1(sort (keys %hash2)){
				my $pointer3=$hash2{$pos1};
				my %hash3=%$pointer3;
				foreach my $pos2(sort (keys %hash3)){
					my $RC=$enzyme_counting{$ref1}->{$ref2}->{$pos1}->{$pos2};
					print $OUT2 join(',', $ref1, $pos1, $ref2, $pos2, $RC), "\n";
				}
			}
		}
	}
	close($OUT2);
	
	#
}

################################################
sub HiC_S4_Tscore_calculation{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	
	print "###Tscore calculation.###\n";
	sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_ligation_frequency_dir" , $variables{ligation_frequency_dir});
	sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_size_bin_csv" , $variables{size_bin_csv});
	sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_enzyme_bin_csv" , $variables{enzyme_bin_csv});
	sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_p_level" , $variables{RC_pvalue_threshold});
	sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_RC_noise" , $variables{RC_noise_threshold});
	sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_permutation_num" , $variables{RC_permutation_num});
	sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_sampling_ratio" , $variables{RC_sampling_ratio});
	
	my $file_names_pointer=sub_3C::files_list($variables{ligation_frequency_dir}, 'file_name');
	my @file_names=@$file_names_pointer;
	my @RC_csv_names=grep(/_RC.csv$/, @file_names);
	#run
	foreach my $RC_csv_name(@RC_csv_names){#3
		sub_3C::refresh_R_script($variables{Tscore_R_file}, "R_RC_csv_name" , $RC_csv_name);
		#run R script
		print "\tRead $RC_csv_name, and run $variables{Tscore_R_file}\n";
		system("Rscript $variables{Tscore_R_file}");
	}#3


}
#########################################################
1; 
# make sure the file returns true or require will not succeed!#



