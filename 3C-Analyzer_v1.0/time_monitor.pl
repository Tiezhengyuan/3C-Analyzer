#!/usr/bin/perl -w
use strict;
use Glib qw(TRUE FALSE);
use Gtk2 -init;
use List::Util;

our $load=0;
our $cpu=0;
our $memory=0;
my $time_monitor_log_file=$ARGV[0];
my $system_monitor_log_file=$ARGV[1];

our (@sample_names, $total_beginning_seconds);
open my($IN), "<", $time_monitor_log_file or die;
while(my $line=<$IN>){
	chomp($line);
	my($key, $value)=split("=", $line);
	my($sample_name, $time_name)=split(":", $key);
	push(@sample_names, $sample_name) unless List::Util::first {$sample_name eq $_} @sample_names;
	$total_beginning_seconds=$value if $key eq 'Total:beginning_time';
}
close($IN);


######################################
#GUI of monitor
my $window = Gtk2::Window->new;
$window->set_default_size ('500', '500');
$window->signal_connect (delete_event => sub {Gtk2->main_quit});
$window->set_title('Progressing monitor');
$window->set_position('center');
	my $table = Gtk2::Table->new(10,6,TRUE);
	$table->set_border_width(5);
		my $frame=Gtk2::Frame->new();
		$frame->set_shadow_type('out');
		$frame->set_border_width(5);
			my $sw=Gtk2::ScrolledWindow->new();
			$sw->set_policy('automatic', 'automatic');
				my $liststore = Gtk2::ListStore->new (qw(Glib::String Glib::String Glib::String Glib::String Glib::String Glib::String));
				my $view = Gtk2::TreeView->new ($liststore);
					my $renderer = Gtk2::CellRendererText->new;
					my $col = Gtk2::TreeViewColumn->new;
					$col->pack_start ($renderer, TRUE);
					$col->add_attribute ($renderer, text => 0);
					$col->set_title ("Number");
				$view->append_column ($col);
					$renderer = Gtk2::CellRendererText->new;
					$col = Gtk2::TreeViewColumn->new;
					$col->pack_start ($renderer, TRUE);
					$col->add_attribute ($renderer, text => 1);
					$col->set_title ("Sample name");
				$view->append_column ($col);
					$renderer = Gtk2::CellRendererProgress->new;
					$col = Gtk2::TreeViewColumn->new;
					$col->pack_start ($renderer, TRUE);
					$col->add_attribute ($renderer, text => 2);
					$col->set_title ("Percentage");
				$view->append_column ($col);
					$renderer = Gtk2::CellRendererText->new;
					$col = Gtk2::TreeViewColumn->new;
					$col->pack_start ($renderer, TRUE);
					$col->add_attribute ($renderer, text => 3);
					$col->set_title ("Beginning time");
				$view->append_column ($col);
					$renderer = Gtk2::CellRendererText->new;
					$col = Gtk2::TreeViewColumn->new;
					$col->pack_start ($renderer, TRUE);
					$col->add_attribute ($renderer, text => 4);
					$col->set_title ("Running time");
				$view->append_column ($col);
					$renderer = Gtk2::CellRendererText->new;
					$col = Gtk2::TreeViewColumn->new;
					$col->pack_start ($renderer, TRUE);
					$col->add_attribute ($renderer, text => 5);
					$col->set_title ("Peak of memory usage %");
				$view->append_column ($col);
				#initiate values of the above five columns
				my $sample_num=1;
				foreach my $sample_name(@sample_names){
					my $iter = $liststore->append;
					$liststore->set ($iter, 0,$sample_num, 1,$sample_name, 2,0, 3,'NA', 4,'NA', 5,0); # start at 50%
					$sample_num++;
				}
				#refresh values
				Glib::Timeout->add (1e5, \&increase_progress_timeout ); #subroutine #refresh screen per 100 second
			$sw->add($view);
		$frame->add($sw);
	$table->attach_defaults($frame,0,6,0,9);
	
		$frame=Gtk2::Frame->new();
		$frame->set_label('Average system load  of the last 15 minutes');
		$frame->set_border_width(5);
			my $entry_1 = Gtk2::Entry->new();
			$entry_1->set_text(0);
		$frame->add($entry_1);
	$table->attach_defaults($frame,0,2,9,10);
		$frame=Gtk2::Frame->new();
		$frame->set_label('Peak of total CPU usage %');
		$frame->set_border_width(5);
			my $entry_2 = Gtk2::Entry->new();
			$entry_2->set_text(0);
		$frame->add($entry_2);
	$table->attach_defaults($frame,2,4,9,10);
		$frame=Gtk2::Frame->new();
		$frame->set_label('Peak of total memory usage %');
		$frame->set_border_width(5);
			my $entry_3 = Gtk2::Entry->new();
			$entry_3->set_text(0);
		$frame->add($entry_3);
	$table->attach_defaults($frame,4,6,9,10);
$window->add ($table);
$window->show_all();
Gtk2->main;


###########
sub increase_progress_timeout {
	my $renderer = shift;
	
	#return load average, cpu usage and memory usage
	open(my $LOAD, "uptime |") or die;
	while (my $line=<$LOAD>){
		chomp($line);
		my @array=split(" ", $line);
		$load=$array[-1];
		$entry_1->set_text($load);
	}
	close($LOAD);
	
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
	#update system monitor log
	my $tmp_time=join(",", localtime(time) );
	my $tmp_duration_seconds=time-$total_beginning_seconds;
	open my($SYS), ">>", $system_monitor_log_file or die;
		print $SYS join("\t", $tmp_time, $tmp_duration_seconds, $tmp_cpu, $tmp_memory), "\n";
	close($SYS);
	#update the values of the time monitor display
	$cpu=$tmp_cpu if $cpu<$tmp_cpu;
	$memory=$tmp_memory if $memory<$tmp_memory;
	$entry_2->set_text($cpu) ;
	$entry_3->set_text($memory) ;
	
	my $found=TRUE;
	my %time;
	open my($IN), "<", $time_monitor_log_file or die;
	while(my $line=<$IN>){
		chomp($line);
		my($key, $value)=split("=", $line);
		my($sample_name, $time_name)=split(":", $key);
		$time{$sample_name}->{$time_name}=$value;
		$found=FALSE if $line=~/running_time/ and $line eq 'NA';
	}
	close($IN);
	my $samples_num=keys %time;
	
	my $num=0;
	#print "#$num:\n";
  	$liststore->foreach(sub{
		my($model, $path, $iter)=@_;
		my $sample_name = $model->get ($iter, 1);
		$num++;
		my $perc=0;
		unless ($time{$sample_name}->{beginning_time} eq 'NA' or $time{$sample_name}->{supposed_time} eq 'NA'){
			my ($sec1,$min1,$hour1,$monthday1,$month1,$year1,$weekday1,$yearday1,$isdaylight1)=split(",",$time{$sample_name}->{beginning_time});
			my ($sec2,$min2,$hour2,$monthday2,$month2,$year2,$weekday2,$yearday2,$isdaylight2)=localtime(time);
			my $duration=($monthday2-$monthday1)*24*60*60 + ($hour2-$hour1)*60*60 + ($min2-$min1)*60 + ($sec2-$sec1); #second
			$perc=int($duration*100/$time{$sample_name}->{supposed_time});
			#print "$sample_name=$perc\n";
		}
		$perc=100 if $perc > 100 or $time{$sample_name}->{running_time} cmp 'NA';
		
		$liststore->set ($iter, 2, $perc.'%');
		unless ($time{$sample_name}->{beginning_time} eq 'NA'){
			my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst)= split(",", $time{$sample_name}->{beginning_time} );
			my $beginning_time=$hour.':'.$min.', '.($mon+1).'/'.$mday.'/'.($year+1900);
			$liststore->set ($iter, 3, $beginning_time) ;
		}
		$liststore->set ($iter, 4, $time{$sample_name}->{running_time} ) unless $time{$sample_name}->{running_time} eq 'NA';
		my $memory_usage = $model->get ($iter, 5);
		$liststore->set ($iter, 5, $memory) if $perc<100 and $perc>0 and $memory_usage<$memory;

		if($num==$samples_num){	return(TRUE);		}
		else{	return(FALSE);		}
	});
	#print "$num:", localtime(time), "\n\n\n";

	return $found; # Call again
}
