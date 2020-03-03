#!/usr/bin/perl -w
use warnings;
use strict;
use Cwd;
use List::Util;
use Gtk2::Pango;
use threads;
use threads::shared;
use Glib qw/TRUE FALSE/;
use Gtk2 qw/-init -threads-init/;
use Bio::SeqIO;
use Bio::Perl;


############################################################
#get the directory of perl scripts involved in Pscore
our $perl_dir=Cwd::getcwd();

#get subroutines 
require $perl_dir."/eRNA_subroutines.pm";

#get variables
my $var_file=$ARGV[0];
#my $var_file='/home/yuan/mysql_pre/eRNA/result/variables.txt';
our $variables_pointer=E_RNA::process_info($var_file);
our %variables=%$variables_pointer;
$variables{bowtie_dir}=$variables{alignment_dir};
$variables{bowtie_output_dir}=$variables{result_dir};

#setup shared hash to control thread
my %shash;
share(%shash); #will work for first level keys
$shash{'go'} = 0;
$shash{'data'} = '';
$shash{'work'} = '';
$shash{'die'} = 0;

our (@fasta_names, @selected_names);

#################################################
#GUI interface
#create window
my $window = Gtk2::Window->new('toplevel');
$window->set_border_width(10);
$window->set_position('center');
$window->set_size_request(1000,800);
$window ->signal_connect( 'destroy' => \&exit );
$window->set_title('Bowtie2 GUI');
	my $table=Gtk2::Table->new(10,6,TRUE);
			#frame: running mode
			my $frame=Gtk2::Frame->new('Running mode of bowtie');
				my $cb_entry_1 = Gtk2::ComboBox->new_text;
				$cb_entry_1->append_text('Index building only');
				$cb_entry_1->append_text('Sequence alignment');
				$cb_entry_1->set_active(1);
			$frame->add($cb_entry_1);
		$table->attach_defaults($frame, 0,2,0,1);
			#frame: output directory
			$frame=Gtk2::Frame->new('Alignment output directory');
				my $out_dir_chooser=Gtk2::Button->new ('Select a folder');
				$out_dir_chooser->signal_connect(clicked=>sub{
					my $file_chooser =Gtk2::FileChooserDialog->new ('Folder selection',  undef, 'select-folder',   'gtk-cancel' => 'cancel', 'gtk-ok' => 'ok');
					my $preview_widget = Gtk2::Label->new ('wheeee');
					$preview_widget->set_line_wrap (TRUE);
					$preview_widget->set_size_request (150, -1);
					$file_chooser->set (preview_widget => $preview_widget,  preview_widget_active => TRUE);
					$file_chooser->signal_connect (selection_changed => sub {
						my $filename = $file_chooser->get_preview_filename;
						my $active = defined $filename && not -d $filename;
						if ($active) {
							my $size = sprintf '%.1fK', (-s $filename) / 1024;
							my $desc = `file '$filename'`;
							$desc =~ s/^$filename:\s*//;
							$preview_widget->set_text ("$size\n$desc");
						}
						$file_chooser->set (preview_widget_active => $active);
					});
					if ('ok' eq $file_chooser->run) {	
						$variables{bowtie_output_dir}=$file_chooser->get_filename;
					}
					$file_chooser->destroy;
					$out_dir_chooser->set_label($variables{bowtie_output_dir});
				});
			$frame->add($out_dir_chooser);
		$table->attach_defaults($frame, 2,4,0,1);
		#frame: bowtie directory
		$frame=Gtk2::Frame->new('Bowtie directory');
			my $dir_chooser=Gtk2::Button->new ('Select a folder');
			$dir_chooser->signal_connect(clicked=>sub{
				my $file_chooser =Gtk2::FileChooserDialog->new ('Folder selection',  undef, 'select-folder',   'gtk-cancel' => 'cancel', 'gtk-ok' => 'ok');
					my $preview_widget = Gtk2::Label->new ('wheeee');
					$preview_widget->set_line_wrap (TRUE);
					$preview_widget->set_size_request (150, -1);
					$file_chooser->set (preview_widget => $preview_widget,  preview_widget_active => TRUE);
					$file_chooser->signal_connect (selection_changed => sub {
						my $filename = $file_chooser->get_preview_filename;
						my $active = defined $filename && not -d $filename;
						if ($active) {
							my $size = sprintf '%.1fK', (-s $filename) / 1024;
							my $desc = `file '$filename'`;
							$desc =~ s/^$filename:\s*//;
							$preview_widget->set_text ("$size\n$desc");
						}
						$file_chooser->set (preview_widget_active => $active);
					});
					if ('ok' eq $file_chooser->run) {	
						$variables{bowtie_dir}=$file_chooser->get_filename;
					}
					$file_chooser->destroy;
					$dir_chooser->set_label($variables{bowtie_dir});
				});
			$frame->add($dir_chooser);
		$table->attach_defaults($frame, 4,6,0,1);
			#frame: input query fasta/fastq file
			$frame=Gtk2::Frame->new('Query input file (*.fa or *.fq format)');
				my $input_chooser =Gtk2::FileChooserButton->new ('select a file' , 'open');
				$input_chooser->set_filename($variables{rawdata_dir});
			$frame->add($input_chooser);
		$table->attach_defaults($frame, 0,2,1,2);
			#frame: references fasta file
			$frame=Gtk2::Frame->new('Reference sequences (*.fa format)');
				my $ref_chooser =Gtk2::FileChooserButton->new ('select a file' , 'open');
				$ref_chooser->set_filename($variables{ref_seq_dir});
			$frame->add($ref_chooser);
		$table->attach_defaults($frame, 2,4,1,2);
		#initiate output name
		$frame=Gtk2::Frame->new('Name of bowtie output');
			my $entry_bowtie_output = Gtk2::Entry->new();
			$entry_bowtie_output->set_text('bowtie_output');
		$frame->add($entry_bowtie_output);
	$table->attach_defaults($frame, 4,6,1,2);
		
		#Frame: alignment mode --end-to-end/--local
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Alignment mode");
			my $cb_alignment_mode = Gtk2::ComboBox->new_text;
			$cb_alignment_mode->append_text('--end-to-end --very-fast');
			$cb_alignment_mode->append_text('--end-to-end --fast');
			$cb_alignment_mode->append_text('--end-to-end --sensitive');
			$cb_alignment_mode->append_text('--end-to-end --very-sensitive');
			$cb_alignment_mode->append_text('--local --very-fast');
			$cb_alignment_mode->append_text('--local --fast');
			$cb_alignment_mode->append_text('--local --sensitive');
			$cb_alignment_mode->append_text('--local --very-sensitive');
			$cb_alignment_mode->set_active(1);
		$frame->add($cb_alignment_mode);
	$table->attach_defaults($frame, 0,2,2,3);
	
		#frame: Other options
			$frame=Gtk2::Frame->new();
			my $label=Gtk2::Label->new('Other options');
		$frame->add($label);
	$table->attach_defaults($frame, 1,5,3,4);
		#Frame:--mp 
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of maximum and minimum mismatch penalties");
			my $cb_mismatch_penality = Gtk2::ComboBox->new_text;
			$cb_mismatch_penality->append_text('4,2');
			$cb_mismatch_penality->append_text('6,2');
			$cb_mismatch_penality->append_text('6,3');
			$cb_mismatch_penality->set_active(1);
		$frame->add($cb_mismatch_penality);
		$table->attach_defaults($frame, 0,2,4,5);
		#Frame:--ma
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Score of match bonus (--ma)");
			my $cb_match_bonus = Gtk2::ComboBox->new_text;
			foreach(1,2,3,4,5,6){
				$cb_match_bonus->append_text($_);
			}
			$cb_match_bonus->set_active(1);
		$frame->add($cb_match_bonus);
	$table->attach_defaults($frame, 2,4,4,5);
		#Frame: --np
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Score of ambiguous penality (--np)");
			my $cb_ambiguous_penalty= Gtk2::ComboBox->new_text;
			foreach (1..4){
				$cb_ambiguous_penalty->append_text($_);
			}
			$cb_ambiguous_penalty->set_active(0);
		$frame->add($cb_ambiguous_penalty);
		$table->attach_defaults($frame, 4,6,4,5);
		#Frame: --rdg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the read gap penalities (--rdg)");
			my $cb_read_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_read_gap_penalty->append_text('3,2');
			$cb_read_gap_penalty->append_text('5,3');
			$cb_read_gap_penalty->append_text('8,4');
			$cb_read_gap_penalty->set_active(1);
		$frame->add($cb_read_gap_penalty);
		$table->attach_defaults($frame, 0,2,5,6);
		#Frame: --rfg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the reference gap penalities (--rfg)");
			my $cb_ref_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_ref_gap_penalty->append_text('3,2');
			$cb_ref_gap_penalty->append_text('5,3');
			$cb_ref_gap_penalty->append_text('8,4');
			$cb_ref_gap_penalty->set_active(1);
		$frame->add($cb_ref_gap_penalty);
		$table->attach_defaults($frame, 2,4,5,6);
		#Frame:-k/-a
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Multiple alignments reporting (-k/-a)");
			my $cb_multiple_alignments = Gtk2::ComboBox->new_text;
			$cb_multiple_alignments->append_text('Distinctly valid alignment');
			$cb_multiple_alignments->append_text('-k 1');
			$cb_multiple_alignments->append_text('-k 2');
			$cb_multiple_alignments->append_text('-k 3');
			$cb_multiple_alignments->append_text('-k 4');
			$cb_multiple_alignments->append_text('-k 5');
			$cb_multiple_alignments->append_text('-k 10');
			$cb_multiple_alignments->append_text('-a');
			$cb_multiple_alignments->set_active(0);
		$frame->add($cb_multiple_alignments);
		$table->attach_defaults($frame, 4,6,5,6);
		#Frame: -D
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Attempts of consecutive seed extension (-D)");
			my $cb_seed_attempts= Gtk2::ComboBox->new_text;
			foreach (15,20,30,50){
				$cb_seed_attempts->append_text($_);
			}
			$cb_seed_attempts->set_active(0);
		$frame->add($cb_seed_attempts);
		$table->attach_defaults($frame, 0,2,6,7);
		#Frame: -R
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Times of re-seeding (-R)");
			my $cb_reseeding= Gtk2::ComboBox->new_text;
			foreach (2,4,6,8){
				$cb_reseeding->append_text($_);
			}
			$cb_reseeding->set_active(0);
		$frame->add($cb_reseeding);
		$table->attach_defaults($frame, 2,4,6,7);
		#frame: -p/--num-threads
		$frame=Gtk2::Frame->new('Number of multiple threads (-p)');
			my $cb_threads_num = Gtk2::ComboBox->new_text;
				foreach (1,2,4,8,16){
					$cb_threads_num->append_text($_);
				}
			$cb_threads_num->set_active(0);
		$frame->add($cb_threads_num);
	$table->attach_defaults($frame, 4,6,6,7);
		#Frame: -I
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Minimum insert for paired-end alignments (-I)");
			my $cb_min_insert = Gtk2::ComboBox->new_text;
			foreach (0,20,40,60,80,100){
				$cb_min_insert->append_text($_);
			}
			$cb_min_insert->set_active(0);
		$frame->add($cb_min_insert);
	$table->attach_defaults($frame, 0,2,7,8);
		#Frame: -X
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Maximum insert for paired-end alignments (-X)");
			my $cb_max_insert = Gtk2::ComboBox->new_text;
			foreach (100, 140,180,220,250,300){
				$cb_max_insert->append_text($_);
			}
			$cb_max_insert->set_active(4);
		$frame->add($cb_max_insert);
	$table->attach_defaults($frame, 2,4,7,8);
		#Frame: --sam-nohead/--sam-noseq
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Header exported in the sam file");
		my $check_button_samhead = Gtk2::CheckButton->new ("Suppress header lines");
		$check_button_samhead->set_active (FALSE);
		$frame->add($check_button_samhead);
		$table->attach_defaults($frame, 0,2,8,9);
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Header exported in the sam file");
		my $check_button_samSQ = Gtk2::CheckButton->new ("Suppress SQ header lines");
		$check_button_samSQ->set_active (FALSE);
		$frame->add($check_button_samSQ);
		$table->attach_defaults($frame, 2,4,8,9);
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Records exported in the sam file");
		my $check_button_unaligned_reads = Gtk2::CheckButton->new ("Suppress records of unaligned reads");
		$check_button_unaligned_reads->set_active (FALSE);
		$frame->add($check_button_unaligned_reads);
		$table->attach_defaults($frame, 4,6,8,9);
		
		

		
		my $hbox = Gtk2::HBox->new( FALSE, 5 );
			my $pbar = Gtk2::ProgressBar->new();
			$pbar->set_pulse_step(.1);
			$pbar->hide; #needs to be called after show_all
		$hbox->add($pbar);
			my $tbutton = Gtk2::Button->new_with_label('Apply and run');
			my $sconnect;
			my $lconnect = $tbutton->signal_connect( clicked => \&launch);
		$hbox->add($tbutton);
			my $ebutton = Gtk2::Button->new_from_stock('gtk-quit');
			$ebutton->signal_connect( clicked =>\&exit );
		$hbox->add($ebutton);
	$table->attach_defaults($hbox, 0,6,10,11);
$window->add($table);

$window->show_all();
Gtk2->main;

#######################################
sub exit{
	$shash{'die'} = 1;
	foreach my $thread(threads->list()){
		$thread->join;
	}
	Gtk2->main_quit;
	return FALSE;
} 
#######################
sub launch{
	
	#save variables
	$variables{bowtie_running_mode}=$cb_entry_1->get_active_text();
	$variables{ref_seq_file}=$ref_chooser->get_filename;
	my @array=split("/", $variables{ref_seq_file});
	$variables{bowtie_index_name}=$array[-1];
	$variables{bowtie_index_name}=~s/\.fasta$|\.fa$//;
	print "Index=$variables{bowtie_index_name}\n";
	$variables{query_input_file}=$input_chooser->get_filename; #default is fasta file
	$variables{bowtie_output_name}=$entry_bowtie_output->get_text();
	
	$variables{bowtie_options}=join(" ", '-x', $variables{bowtie_index_name}, 
				'-S', $variables{bowtie_output_dir}.'/'.$variables{bowtie_output_name}.'.sam',
				$cb_alignment_mode->get_active_text(), 
				'-p', $cb_threads_num->get_active_text(),
				'--mp', $cb_mismatch_penality->get_active_text(),
				'--np', $cb_ambiguous_penalty->get_active_text(),
				'--rdg', $cb_read_gap_penalty->get_active_text(),
				'--rfg', $cb_ref_gap_penalty->get_active_text(),
				'-D', $cb_seed_attempts->get_active_text(),
				'-R', $cb_reseeding->get_active_text(),
			);
	#output: sam file
	$variables{bowtie_options} .= " --no-hd " if $check_button_samhead->get_active();
	$variables{bowtie_options} .= " --no-sq " if $check_button_samSQ->get_active();
	$variables{bowtie_options} .= " --no-unal " if $check_button_unaligned_reads->get_active();
	
	# input file 
	if ($variables{query_input_file}=~/\.fastq$|\.fq$/ and $variables{query_input_file}=~/_R1_/){
		my $R2_file=$variables{query_input_file};
		$R2_file=~s/_R1_/_R2_/;
		if (-f $R2_file){
			$variables{bowtie_options}=join(" ", $variables{bowtie_options}, 
										'-I', $cb_min_insert->get_active_text(),
										'-X', $cb_max_insert->get_active_text(),
										'-q', '-1', $variables{query_input_file}, '-2', $R2_file, );
		}
		else{	$variables{bowtie_options} .= " -q $variables{query_input_file} ";		}
	}
	elsif ($variables{query_input_file}=~/\.fastq$|\.fq$/ and $variables{query_input_file}=~/_R2_/){
		my $R1_file=$variables{query_input_file};
		$R1_file=~s/_R2_/_R1_/;
		if (-f $R1_file){
			$variables{bowtie_options}=join(" ", $variables{bowtie_options}, 
										'-I', $cb_min_insert->get_active_text(),
										'-X', $cb_max_insert->get_active_text(),
										'-q', '-1', $R1_file, '-2', $variables{query_input_file} );
		}
		else{	$variables{bowtie_options} .= " -q $variables{query_input_file} ";		}
	}
	else {
		$variables{bowtie_options} .= " -f $variables{query_input_file} ";
	}

	$variables{bowtie_options} .= " $cb_multiple_alignments->get_active_text() " unless $cb_multiple_alignments->get_active_text()=~/alignment/;
	if ($cb_alignment_mode->get_active_text()=~/--local/){
		$variables{bowtie_options} .= " --ma $cb_match_bonus->get_active_text() " ;
	}
	#print "$variables{bowtie_dir}\t$variables{bowtie_output_dir}\n";
	#print "$variables{ref_seq_file}\t$variables{query_input_file}\n";
	#print "$variables{bowtie_options}\n";

		
		#create 1 sleeping thread
		my $thread_1 = threads->new(\&work);
			
		$shash{'go'} = 1;
		$pbar->show;
		$tbutton->set_label('Stop');
		$tbutton->signal_handler_block($lconnect);
		$sconnect = $tbutton->signal_connect( clicked => sub{ 	
			$shash{'go'} = 0;
			$tbutton->set_label('Run');
			$tbutton->signal_handler_block ($sconnect);
			$tbutton->signal_handler_unblock ($lconnect);
		});

		Glib::Timeout->add (100, sub {
			if($shash{'go'} == 1){
				$pbar->set_text('Running!');
				$pbar->pulse;
				return TRUE;
			}
			else{	
				$pbar->set_text('OK! It is done!');
				$tbutton->set_label('Run');
				return FALSE;
			}
		});

}

################################################## #######
sub work{
	$|++;
	while(1){
		return if $shash{'die'} == 1;

		if ( $shash{'go'} == 1 ){#2
			#current directory
			chdir $variables{bowtie_dir};
			
			unless (-f $variables{bowtie_index_name}.'.1.bt2' and -f $variables{bowtie_index_name}.'.rev.1.bt2'){
				print "Build bowtie index\n\n\n";
				system("$variables{bowtie_dir}/bowtie2-build $variables{ref_seq_file} $variables{bowtie_index_name} ");
			}
			
			print "Bowtie alignment begin:\n\n\n";
			print "\n$variables{bowtie_dir}/bowtie2 $variables{bowtie_options}\n\n";
			system("$variables{bowtie_dir}/bowtie2 $variables{bowtie_options}");
			
			last if $shash{'go'} == 0;
			return if $shash{'die'} == 1;
			$shash{'go'} = 0;   #turn off 
		}#2
		else{ sleep 1; }
	}
}

#####################################
#used for items' selection
sub items_selection{
	my ($index_names_pointer)=@_;
	my @index_names=@$index_names_pointer;

	
	my $table=Gtk2::Table->new(4,8,TRUE);
	$table->set_sensitive (TRUE);
	
	#create a scrolled window that will host the treeview
	my $sw_1 = Gtk2::ScrolledWindow->new (undef, undef);
	$sw_1->set_shadow_type ('etched-out');
	$sw_1->set_policy ('automatic', 'automatic');
	$sw_1->set_border_width(5);
	
		#this is one of the provided base Gtk2::TreeModel classes.
		my $tree_store_1 = Gtk2::TreeStore->new(qw/Glib::String/);
		#fill it with arbitry data
		foreach (@index_names) {
			my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_1->set ($iter,0 => $_);
		}
	
		#this will create a treeview, specify $tree_store as its model
		my $tree_view_1 = Gtk2::TreeView->new($tree_store_1);
		#create a Gtk2::TreeViewColumn to add to $tree_view
		my $tree_column_1 = Gtk2::TreeViewColumn->new();
		$tree_column_1->set_title ("Sort candidates");
			#create a renderer that will be used to display info in the model
			my $renderer_1 = Gtk2::CellRendererText->new;
		#add this renderer to $tree_column. This works like a Gtk2::Hbox. so you can add more than one renderer to $tree_column			
		$tree_column_1->pack_start ($renderer_1, FALSE);
		$tree_column_1->add_attribute($renderer_1, text => 0);
		$tree_column_1->set_sort_column_id(0);# Allow sorting on the column
	
		#add $tree_column to the treeview
		$tree_view_1->append_column ($tree_column_1);
		$tree_view_1->set_search_column(0);# make it searchable
		$tree_view_1->set_reorderable(TRUE);# Allow drag and drop reordering of rows
	
	$sw_1->add($tree_view_1);
	$table->attach_defaults($sw_1, 0,3,0,6);

	#create a scrolled window that will host the treeview
	my $sw_2 = Gtk2::ScrolledWindow->new (undef, undef);
	$sw_2->set_shadow_type ('etched-out');
	$sw_2->set_policy ('automatic', 'automatic');
	$sw_2->set_border_width(5);
	
		#this is one of the provided base Gtk2::TreeModel classes.
		my $tree_store_2 = Gtk2::TreeStore->new(qw/Glib::String/);
		my $tree_view_2 = Gtk2::TreeView->new($tree_store_2);

		#create a Gtk2::TreeViewColumn to add to $tree_view
		my $tree_column_2 = Gtk2::TreeViewColumn->new();
		$tree_column_2->set_title ("Sort selected items");
			
			#create a renderer that will be used to display info in the model
			my $renderer_2 = Gtk2::CellRendererText->new;
		#add this renderer to $tree_column. 
		$tree_column_2->pack_start ($renderer_2, FALSE);
		
		 # set the cell "text" attribute to column 0 - retrieve text from that column in treestore 
		$tree_column_2->add_attribute($renderer_2, text => 0);
		$tree_column_2->set_sort_column_id(0);# Allow sorting on the column
	
		#add $tree_column to the treeview
		$tree_view_2->append_column ($tree_column_2);
		$tree_view_2->set_search_column(0);# make it searchable
		$tree_view_2->set_reorderable(TRUE);# Allow drag and drop reordering of rows

	$sw_2->add($tree_view_2);
	$table->attach_defaults($sw_2, 4,7,0,6);

	#Widget of table: button
	my $button_1=Gtk2::Button->new('>');        ######>> button: one transfered
	$button_1->signal_connect(clicked=>sub{
		my $tree_selection_1=$tree_view_1->get_selection;
		my($tree_store_1, $iter_1)=$tree_selection_1->get_selected;
		my $value=$tree_store_1->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_2->append(undef);  #add the value into the second tree 
		$tree_store_2->set($iter,0=> $value);
		push(@selected_names, $value);
	} );
	$table->attach_defaults($button_1, 3,4,1,2);
	my $button_2=Gtk2::Button->new('>>');        ###### >> button: all transfered
	$button_2->signal_connect(clicked=>sub{
		$tree_store_1->clear;
		foreach (@index_names) {
			my $iter = $tree_store_2->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_2->set ($iter,0 => $_);
		}
		@selected_names=@index_names;
	} );
	$table->attach_defaults($button_2, 3,4,2,3);
	my $button_3=Gtk2::Button->new('<');      ######<< button
	$button_3->signal_connect(clicked=>sub{
		my $tree_selection_2=$tree_view_2->get_selection;
		my($tree_store_2, $iter_2)=$tree_selection_2->get_selected;
		my $value=$tree_store_2->get($iter_2, 0);
		return unless $iter_2;
		$tree_store_2->remove($iter_2);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		foreach(my $i=0; $i<@selected_names; $i++){
			delete $selected_names[$i] if $selected_names[$i] eq $value;
		}
		@selected_names=grep($_, @selected_names);
	} );
	$table->attach_defaults($button_3, 3,4,3,4);
	my $button_4=Gtk2::Button->new('<<');      ######<< button
	$button_4->signal_connect(clicked=>sub{
		$tree_store_2->clear;
		foreach (@index_names) {
			my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_1->set ($iter,0 => $_);
		}
		undef @selected_names;
	} );
	$table->attach_defaults($button_4, 3,4,4,5);
	my $button_5=Gtk2::Button->new('Up');        ######Up button
	$button_5->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_2->get_selection;
		my($tree_store, $iter)=$tree_selection->get_selected;
		my $path = $tree_store->get_path($iter);
		$path->prev;
		my $prev_iter = $tree_store->get_iter($path);
		$tree_store->move_before($iter,$prev_iter );
		
		my $value=$tree_store->get($iter, 0);
		my $pre_value=$tree_store->get($prev_iter, 0);
		foreach(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i>0){
				$selected_names[$i-1]=$value;
				$selected_names[$i]=$pre_value;
			}
		}
		@selected_names=grep($_, @selected_names);
	} );
	$table->attach_defaults($button_5, 7,8,1,2);
	my $button_6=Gtk2::Button->new('Down');   ######Down button
	$button_6->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_2->get_selection;
		my($tree_store, $iter)=$tree_selection->get_selected;
		my $next_iter = $tree_store->iter_next($iter);
		$tree_store->move_before($next_iter,$iter );
				
		my $value=$tree_store->get($iter, 0);
		my $next_value=$tree_store->get($next_iter, 0);
		foreach(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i<@selected_names-1){
				$selected_names[$i+1]=$value;
				$selected_names[$i]=$next_value;
			}
		}
		@selected_names=grep($_, @selected_names);
	} );
	$table->attach_defaults($button_6, 7,8,4,5);

	return ($table);
}

