#! /usr/bin/perl -w
use strict;
use warnings;
#use threads;
use List::MoreUtils;
use List::Util;
use File::Find;


#the file constains all subroutines required for running the pipeline of 3C_analyzer

package sub_basic;

################################################################
#get date and time
sub get_time{
	my ($begin_time_str, $end_time_str, $type)=@_;
	$end_time_str='0,0,0,0,0,0,0,0,0' if $end_time_str eq '0';
	$type='duration' unless $type;
	
	my $my_time=0;
	my ($sec1,$min1,$hour1,$monthday1,$month1,$year1,$weekday1,$yearday1,$isdaylight1)=split(',', $begin_time_str);
	$year1 += 1900;
	my $begin_time=$year1.'/'.$month1.'/'.$monthday1.', '.$hour1.':'.$min1.':'.$sec1;	
	my ($sec2,$min2,$hour2,$monthday2,$month2,$year2,$weekday2,$yearday2,$isdaylight2)=split(',', $end_time_str);
	$year2 += 1900;
	my $end_time=$year2.'/'.$month2.'/'.$monthday2.', '.$hour2.':'.$min2.':'.$sec2;
	my $duration_m=($year2-$year1)*365*24*60 + ($yearday2-$yearday1)*24*60 + ($hour2-$hour1)*60 + ($min2-$min1);
	#
	if($type eq 'duration'){
		$my_time = $begin_time.'-----'.$end_time.'. Duration: '.$duration_m.'min' if $duration_m > 0;
	}
	elsif($type eq 'early'){
		$my_time= ($duration_m>=0) ? $begin_time_str : $end_time_str ;
	}
	elsif($type eq 'late'){
		$my_time= ($duration_m<=0) ? $begin_time_str : $end_time_str ;
	}
	elsif($type eq 'seconds'){
		$my_time = $duration_m*60 if $duration_m>0;
	}
	elsif($type eq 'readable'){
		$my_time = $hour1.':'.$min1.', '.($month1+1).'/'.$monthday1.'/'.$year1;
	}
	#
	return($my_time);
}


#############################################
#type= 'replace' or 'add' or 'read', the default is replace
sub refresh_log{
	my ($log_file, $r_name, $r_value, $type)=@_;
	$type='replace' unless $type;
	
	my %variables;
	#read old data
	if(-f $log_file){
		my $variables_pointer=file_to_hash($log_file, "=");
		%variables=%$variables_pointer;
	}
	
	#refresh new data
	if (exists $variables{$r_name} and $type eq 'add'){
		$variables{$r_name}=$variables{$r_name}+$r_value;
	}
	else{
		$variables{$r_name}=$r_value;
	}
	
	#export refreshed data
	hash_to_file(\%variables, $log_file, "=");
	
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

##############################
#
sub file_operation{
	my ($file, $type)=@_;
	
	#
	my @array=split("/", $file);
	my $file_name=$array[-1];
	my($file_name_head, $file_name_tail)=split(/\./, $file_name);
		
	my $out;
	if($type eq 'file_name'){# file name
		$out=$file_name;
	}
	elsif($type eq 'name_head'){# file name head
		$out=$file_name_head;
	}
	elsif($type eq 'name_tail'){# file name tail
		$out=$file_name_tail;
	}
	elsif($type eq 'file_size'){# file size
		my @stats=stat($file);
		$out=$stats[7];
	}
	return($out);
}
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


#########################
sub copy_files{
	my ($in_dir, $out_dir, $file_tail)=@_;
	$in_dir .= '/' unless $in_dir=~/\/$/;
	$out_dir .= '/' unless $out_dir=~/\/$/;
	mkdir($out_dir, 0755) unless -d $out_dir;
	
	my $file_names_pointer=files_list($in_dir,'file_name');
	my @file_names=@$file_names_pointer;
	my @sub_file_names=grep(/\.$file_tail$/, @file_names);
	#
	my @out_files;
	if(@sub_file_names>0){#2
		foreach my$file_name(@sub_file_names){#3
			my $in_file=$in_dir.$file_name;
			my $out_file=$out_dir.$file_name;
			system("cp $in_file $out_file"); 
			push(@out_files, $out_file);
		}#3
	}#2
	
	return(\@out_files);
}

######################
#######3
sub combine_files{
	my ($sample_dir, $file_tail, $out_file_head)=@_;
    
    #get files before combination
	my $files_pointer=files_list($sample_dir, 'file');
	my @files=@$files_pointer;
	my @sub_files=grep(/$file_tail$/, @files);
	
	#combine files
	my $in_sub_files_str=join (' ', @sub_files);
	my $out_file=$out_file_head.'.'.$file_tail;
	system("cat $in_sub_files_str > $out_file") if @sub_files>0; #combine
	
}
######################
#
sub combine_log_files{
	my ($variables_pointer)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(",", $variables{sample_names});
	
	my $total_parts_status=0;
	my $found=1; #judge if running time;
	#print "Fresh sample log files in $variables{sample_log_dir}\n";
	foreach my $sample_name(@sample_names){#2
		my $log_files_pointer=copy_files($variables{result_dir}.'/'.$sample_name, $variables{sample_log_dir}.'/'.$sample_name, 'log');
		my @log_files=@$log_files_pointer;
		
		#initiate hash
		my %hash;
		$hash{sample_name}=$sample_name;
		$hash{status}='no';
		
		#refresh hash
		my @log_status;
		my $parts_status=0;
		foreach my $log_file(@log_files){#3
			my $log_variables_pointer=file_to_hash($log_file, '=');
			my %log_variables=%$log_variables_pointer;
			$log_variables{status}='no' unless exists $log_variables{status};
			
			unless($log_variables{status} eq 'no'){#4
				foreach my $key(keys %log_variables){#5
					if(exists $hash{$key} and $key eq 'beginning_time'){
						$hash{$key}=get_time($hash{$key}, $log_variables{$key},'early');
					}
					elsif(exists $hash{$key} and $key eq 'ending_time'){
						$hash{$key}=get_time($hash{$key}, $log_variables{$key},'late');
					}
					elsif(exists $hash{$key} and $key=~/_num$/){
						$hash{$key} += $log_variables{$key};
					}
					elsif($key eq 'status'){
						push(@log_status, $log_variables{$key}) unless List::Util::first {$_ eq $log_variables{$key}} @log_status;
						$parts_status++ if $log_variables{$key} eq 'off';
					}
					else{
						$hash{$key}=$log_variables{$key};
					}
				}#5
			}#4
		}#3
		$hash{'supposed_time'}=$sample_info{$sample_name}->{'supposed_time'};
		$hash{'parallel_parts'}=$sample_info{$sample_name}->{'parallel_parts'};
		$hash{'parts_status'}=$parts_status;
		$total_parts_status += $parts_status;
		
		#status by samples
		if(@log_status==1 and $log_status[0] eq 'off'){
			$hash{status}='off';
			$found=0;
			$hash{'running_time'}=get_time($hash{'beginning_time'}, $hash{'ending_time'},'duration');
		}
		elsif(List::Util::first {$_ eq 'on'} @log_status){
			$hash{status}='on';
			$hash{'ending_time'}='NA';
			$hash{'running_time'}='NA';
		}
		
		#export $hash;
		my $sample_log=$variables{sample_log_dir}.'/'.$sample_name.'.log';
		hash_to_file(\%hash, $sample_log, '=');
	}#2
	#refresh the file Total.log
	sub_basic::refresh_log($variables{total_log_file}, 'parts_status', $total_parts_status);
	
	return($found);
}

#################################
#calculate free space of $raw_data_dir and $result_dir
sub free_space {
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
############################################################
sub check_bowtie_index{
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

############################################
#######
sub hash_to_file{
	my($hash_pointer, $out_file, $sep)=@_;
	my %hash=%$hash_pointer;
	$sep="=" unless $sep; #default split char is '='
	
	if ($out_file){
		open my($OUT), ">", $out_file or die "Err: Cannot write $out_file!!!\n";
		foreach my$key(sort (keys %hash)){
				print $OUT join("", $key, $sep, $hash{$key}), "\n";
		}
		close($OUT);
	}
	else{
		print "Err: The file name is null when hash_to_file!!!\n";
	}
	#
}
########################
#
sub file_to_hash{
	my ($file, $sep, $key_index, $value_index)=@_;
	$sep='=' unless $sep;
	$key_index=0 unless $key_index;
	$value_index=1 unless $value_index;
	
	my %hash;
	if(-f $file){
		open my ($INFO), "<", $file or die;
		while (<$INFO>) {#2
			chomp($_);
			my @items = split(/$sep/, $_); #split on the tabs
			my $key=$items[$key_index];
			my $value=$items[$value_index];
			#print"$key:$value\n";
			if ($_=~/$sep/){
				if(exists $hash{$key}){
					$hash{$key} .= ','.$value ;  
				}
				else{
					$hash{$key} =$value ;  
				}
			}
		}#2
		close($INFO);
	}
	
	return(\%hash);
}
##########################
#read file as the hash with two levels
sub file_to_hash2{
	my($file, $sep)=@_;
	$sep="\t" unless $sep;
	
	#
	my %hash2;
	open my($OUT), "<", $file or die;
	my $header=<$OUT>;
	chomp($header);
	#print "####$header###\n";
	my @col_names=split(/$sep/,$header);
	#
	while (<$OUT>){
		chomp($_);
		my @items=split(/$sep/, $_);
		my $row_name=$items[0];
		for(my $i=1;$i<@items; $i++){
			my $col_name=$col_names[$i];
			$hash2{$row_name}->{$col_name}=$items[$i];
			#print "$row_name:$col_name=$items[$i]\n";
		}
	}
	close($OUT);
	
	return(\%hash2);
}


#############################
#export the hash with two levels into text file
sub hash2_to_file{
	my($hash2_pointer, $statistics_file, $sep)=@_;
	my %hash2=%$hash2_pointer;
	$sep="\t" unless $sep;
	
	#get keys1 and key2
	my @row_names=keys %hash2;
	@row_names= sort @row_names;
	my @col_names=map { keys $hash2{$_} } @row_names;
	@col_names=List::MoreUtils::uniq @col_names;
	@col_names= sort @col_names;
	
	open my($OUT), ">", $statistics_file or die;
	print $OUT join($sep, 'names', @col_names), "\n";
	foreach my $row_name(@row_names){
		my @counting_num;
		foreach my $col_name(@col_names){
			$hash2{$row_name}->{$col_name}=0 unless exists $hash2{$row_name}->{$col_name};
			push(@counting_num, $hash2{$row_name}->{$col_name});
		}
		print $OUT join($sep, $row_name, @counting_num), "\n";
	}
	close($OUT);
	#
}
###################
sub combine_hash2{
	my($h1_pointer, $h2_pointer)=@_;
	my %h1=%$h1_pointer;
	my %h2=%$h2_pointer;
	
	foreach my $row(keys %h2){#2
		my $hash_pointer=$h2{$row};
		my %hash=%$hash_pointer;
		foreach my $col(keys %hash){
			if(exists $h1{$row}->{$col}){
				$h1{$row}->{$col} += $h2{$row}->{$col};
			}
			else{
				$h1{$row}->{$col}=$h2{$row}->{$col};
			}
		}
	}#2
	
	return(\%h1);
}
##################
#
sub print_hash{
	my($pointer)=@_;
	my %hash=%$pointer;
	
	foreach my $key (sort (keys %hash)){
		print "$key=$hash{$key}\n";
	
	}

}
###########################################
#initiate the monitor log file, the system log file, and the sample log files
sub initiate_log_files{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my @sample_names=split(',', $variables{sample_names});
	my $total_beginning_time= join(",", localtime(time) );
	
	#initiate time_monitor.log
	my $total_supposed_time=0;
	my $total_parallel_parts=0;
	open my($MON), ">", $variables{time_monitor_file} or die; 
	print $MON join("\t", 'sample_names', 'supposed_time', 'beginning_time', 'ending_time', 'running_time', 'parallel_parts', 'parts_status', 'status'), "\n";
	foreach my $sample_name(@sample_names) {
		print $MON join("\t", $sample_name, $sample_info{$sample_name}->{supposed_time}, 'NA', 'NA', 'NA', 
							$sample_info{$sample_name}->{parallel_parts}, 0, 'no'), "\n";
		$total_supposed_time += $sample_info{$sample_name}->{supposed_time};
		$total_parallel_parts += $sample_info{$sample_name}->{parallel_parts};
	}
	print $MON join("\t", 'Total', $total_supposed_time, $total_beginning_time, 'NA', 'NA', $total_parallel_parts, 0, 'no'), "\n";
	close($MON);
	
	#initial Total.log file
	sub_basic::refresh_log($variables{total_log_file}, 'beginning_time', $total_beginning_time);
	sub_basic::refresh_log($variables{total_log_file}, 'supposed_time', $total_supposed_time);
	sub_basic::refresh_log($variables{total_log_file}, 'parallel_parts', $total_parallel_parts);
	sub_basic::refresh_log($variables{total_log_file}, 'parts_status', 0);
	sub_basic::refresh_log($variables{total_log_file}, 'status', 'on');

	#initiate and clear the monitor.log
	open my($SYS), ">", $variables{system_monitor_file} or die;
	print $SYS join("\t", 'Time', 'Duration(s)', 'CPU_usage(%)', 'Memory_usage(%)'), "\n";
	close($SYS);
	
	#initiate sample log files
	foreach my$sample_name(@sample_names){
		my $sample_dir=$variables{result_dir}.'/'.$sample_name;
		mkdir($sample_dir, 0755) unless -d $sample_dir;
		my $files_pointer=files_list($sample_dir,'file');
		my @files=@$files_pointer;
		my @log_files=grep(/\.log$/, @files);
		foreach my $log_file(@log_files){
			#refresh_log($log_file, 'status', 'no');
		}
	}
	#
	return($total_beginning_time);
}
###################
sub refresh_monitor_log{
	my $variables_pointer=$_[0];
	my %variables=%$variables_pointer;
	my @sample_names=split(',', $variables{sample_names});
	my @monitor_items=('beginning_time', 'ending_time', 'running_time', 'parts_status', 'status');
	
	#read the file time_monitor.log
	my $time_pointer=sub_basic::file_to_hash2($variables{time_monitor_file});
	my %time=%$time_pointer;
	
	#refresh by reading sample log
	foreach my $log_name(keys %time) {#2
		my $log_pointer=file_to_hash($variables{sample_log_dir}.'/'.$log_name.'.log', '=');
		my %log=%$log_pointer;
		$log{status}='no' unless exists $log{status};
		if($log{status} eq 'on' or $log{status} eq 'off'){#3
			foreach my $name(@monitor_items){
				$time{$log_name}->{$name}=$log{$name} if exists $log{$name};
			}
		}#3
	}#2
	
	#update
	sub_basic::hash2_to_file(\%time, $variables{time_monitor_file});
	#
	return(\%time);
}
######################
#initiate start/end time and sample directory
sub initiate_starting_time{
	my ($variables_pointer, $raw_file)=@_;
	my %variables=%$variables_pointer;
	my $sample_info_pointer=$variables{sample_info_pointer};
	my %sample_info=%$sample_info_pointer;
	my $raw_to_sample_pointer=$variables{raw_to_sample_pointer};
	my %raw_to_sample=%$raw_to_sample_pointer;
	
	#refresh sample_dir
	$variables{sample_name}=$raw_to_sample{$raw_file}; #sample name
	$variables{sample_dir}=$variables{result_dir}.'/'.$variables{sample_name};
	mkdir($variables{sample_dir}, 0755) unless -d $variables{sample_dir};
	#Note: The result directory is temporarily changed to $variables{sample_dir} for $variables{sample_name}

	#refresh sample_log
	$variables{beginning_time}=join(',', localtime(time) );

	my $raw_file_name_head=sub_basic::file_operation($raw_file, 'name_head');
	$variables{sample_out_file}=$variables{sample_dir}.'/'.$raw_file_name_head;
	$variables{sample_log}=$variables{sample_out_file}.'.log';
	$variables{supposed_time}=int($sample_info{$variables{sample_name} }->{'supposed_time'}/$sample_info{$variables{sample_name} }->{'parallel_parts'}); #unit is second
	refresh_log($variables{sample_log}, "sample_name", $variables{sample_name});
	refresh_log($variables{sample_log}, "beginning_time" , $variables{beginning_time});
	refresh_log($variables{sample_log}, "supposed_time" , $variables{supposed_time});
	refresh_log($variables{sample_log}, "status" , 'on');
	
	return(\%variables);
}

##########
sub initiate_ending_time{
	my ($log_file, $beginning_time)=@_;
	
	my $ending_time=join(',', localtime(time) );
	my $running_time=get_time($beginning_time, $ending_time, 'duration');
	refresh_log($log_file, "ending_time" , $ending_time);
	refresh_log($log_file, "running_time", $running_time);
	refresh_log($log_file, "status", 'off');
	
}
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

##########################################
sub reverse_complement {
	my $dna_seq = $_[0];

	# reverse the DNA sequence
	my $revcom_seq = reverse($dna_seq);
	# complement the reversed DNA sequence
	$revcom_seq =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;

	return ($revcom_seq);
}
##################################################################
#truncate sequence of adapter from read
sub truncate_end3_seq{#1
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
################
#split hybrid sequence by enzyme site with known sequence
sub split_hybrid_seq{
	my($sequence, $known_seq, $string, $type)=@_;
	
	my $trunc_seq=$sequence;
	if($known_seq=~/$sequence/){#2
		$trunc_seq='NA';
	}#2
	elsif($type eq 'head'){#2
		my $pos = 0;
		for (my $pos = index($sequence, $string, $pos); $pos >= 0; $pos = index($sequence, $string, $pos) ) {
				my $head_seq=substr($sequence, 0, $pos+length($string));
				my $tail_seq=substr($sequence, $pos );
				#print "$pos:$head_seq\n$tail_seq\n\n";
				if($known_seq=~/$tail_seq/ and length($tail_seq)>10 ){
					$trunc_seq=$head_seq;
					last;
				}
				$pos += length $string;
		}
	}#2
	elsif($type eq 'tail'){#2
		my $pos = length($sequence)-1;
		for ( my $pos = rindex($sequence, $string, $pos); $pos >= 0; $pos = rindex($sequence, $string, $pos) ) {
			my $head_seq=substr($sequence, 0, $pos+length($string));
			my $tail_seq=substr($sequence, $pos );
			 if($known_seq=~/$head_seq/ and length($head_seq)>10){
				$trunc_seq=$tail_seq;
				#print "##$known_seq##\n";
				#print"$pos:$head_seq\n$tail_seq\n\n";
				last;
			}
			$pos -= length $string;
		}
	}#2
	
	return($trunc_seq);
}

#####################
#split sequence by enzyme sites
sub split_seq{
	my ($seq, $enzyme_site)=@_;
	
	my $head_seq=$seq;
	my $tail_seq='NA';
	my @a=split(/$enzyme_site/, $seq);
	if (@a>=2){
		$tail_seq=pop @a;
		$head_seq=join($enzyme_site, @a);
		$head_seq .= $enzyme_site;
	}
	
	return($head_seq);
}
##########################
#get genome sequences from fasta file
sub genome_sequences{
	my($genome_fasta_file)=@_;
	
	my %genome;
	open my($IN), "<",  $genome_fasta_file or die ;
	my($chr,$seq);
	while( defined ($chr=<$IN>) && defined ($seq=<$IN>) ){#2
		chomp($chr, $seq);
		$chr=~s/^>//;
		#my $seq_len=length($seq);
		$genome{$chr}=$seq;
		print "\t$chr......\n";
	}#2
	
	return(\%genome);
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

#################
sub system_loading{
	
	#return load average
	open(my $LOAD, "uptime |") or die;
		my $line=<$LOAD>;
		chomp($line);
		my @array=split(" ", $line);
		my $load=$array[-1];
	close($LOAD);
	
	#return cpu usage and memory usage
	my $tmp_cpu=0;
	my $tmp_memory=0;
	open(my $PROC, "ps u |") or die;
	while (my $line=<$PROC>){
		chomp($line);
		my @array=split(" ", $line);
		if ($array[2] cmp '%CPU' and $array[2]>20){
			$tmp_cpu +=$array[2];
			$tmp_memory +=$array[3];
		}
	}
	close($PROC);
	
	my $out=join(',', $load, $tmp_cpu, $tmp_memory);
	return($out);
}


###############
sub scaling_normalization{
	my($hash2_pointer, $raw_reads_num)=@_;
	my %hash2=%$hash2_pointer;
	
	my %norm_hash2;
	foreach my $row(keys %hash2){
		my $hash_pointer=$hash2{$row};
		my %hash=%$hash_pointer;
		foreach my $col(keys %hash){
			my $counts=$hash2{$row}->{$col};
			$norm_hash2{$row}->{$col}=int(($counts*1e6)/$raw_reads_num + 0.5);
		}
	}

	return(\%norm_hash2);
}
#########################################################
1; 
# make sure the file returns true or require will not succeed!#


