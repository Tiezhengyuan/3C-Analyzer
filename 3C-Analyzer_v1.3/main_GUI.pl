#! /usr/bin/perl -w
use warnings;
use strict;
use Cwd;
use Gtk2 qw/-init -threads-init/;
use Glib qw/TRUE FALSE/; 
use Gtk2::Pango;
use Gtk2::SimpleList;
use threads; 
use threads::shared;
use List::Util;
use List::MoreUtils;
use File::Find;

#get the directory of perl scripts involved in Pscore
our $perl_dir=Cwd::getcwd();
print $perl_dir;
#get subroutines
require $perl_dir."/functions_basic.pm"; #sub_basic::
require $perl_dir."/functions_3C.pm";  #sub_3C::
require $perl_dir."/functions_gui.pm";  #sub_gui::

##########33
#get the parameters of 3C_analyzer 
#DpnII='GATC',  HindIII='AAGCTT'
our %variables=('probe_index_len'=>15, 'fragment_index_len'=>20, 'threads_num'=>8, 'size_bin_len'=>10000, 
			'enzyme_site'=>'GAATTC', 'enzyme_bin_enzymes'=>1, 'enzyme_bin_lower'=>1000, 'enzyme_bin_upper'=>20000,
			'raw_file_format'=>'FASTQ', 'sequencing_end'=>2, 'colocalization_range'=>200,
			'adapter_5'=>'AGATCGGAAGAGCGTCGTGT',
			'adapter_3'=>'AGATCGGAAGAGCACACGTCTGAACTCCAGT',
			'adapter_P5'=>'AATGATACGGCGACCACCGAACACTCTTTCCCTACACGACGCTCTTCCGATCT',
			'adapter_len'=>14, 'min_trim_len'=>10, 'RC_noise'=>1, 'query_length'=>12,
			'3C_trimming'=>'no', '3C_colocalization'=>'no', '3C_other_detection'=>'no', 
			'3C_statistics'=>'no', '3C_RC_counting'=>'no', '3C_Tscore'=>'no', 	);
#default directories
$variables{software_dir}=$perl_dir;
$variables{raw_data_dir}=$perl_dir.'/raw_data';
$variables{result_dir}=$perl_dir.'/result';
$variables{alignment_dir}=$perl_dir.'/bowtie2';
$variables{ref_dir}=$perl_dir.'/references';
#initiate variables in variables.txt
sub_basic::hash_to_file(\%variables, $variables{result_dir}.'/variables.txt', "=");

#parameters of 3C
$variables{var_file}=$variables{result_dir}.'/variables.txt';
#edited by users in adavance
$variables{site_info_file}=$variables{result_dir}.'/site_info.csv';
#log file of 3C_analyzer
$variables{log_file}=$variables{result_dir}.'/3C_analyzer.log'; 

#setup shared hash to control thread
my %shash;
share(%shash); #will work for first level keys
$shash{'go'} = 0;
$shash{'predict'}=0;
$shash{'work'} = '';
$shash{'die'} = 0;


################################


##############################################################################
#standard window creation, placement, and signal connecting
my $window=sub_gui::default_window("3C-analyzer: Discovery of Chromatin Interactions from 3C-seq", 900, 600);
	my $table = Gtk2::Table->new(9,6,TRUE);
	$table->set_border_width(10);

		#create a menu bar with two menu items, one will have $menu_edit as a submenu
		my $menu_bar = Gtk2::MenuBar->new;
			#add first menu item
			my $menu_item_edit= Gtk2::MenuItem->new('Outer tools');
			my $menu_edit=&menu_tools();
			$menu_item_edit->set_submenu ($menu_edit);#set its sub menu
		$menu_bar->append($menu_item_edit);
			#add third menu item
			my $menu_item_help = Gtk2::MenuItem->new('Help');
			$menu_item_help->set_sensitive(TRUE);
			$menu_item_help->set_right_justified(TRUE);#push it off to the right
			my $menu_help=&menu_help();
			$menu_item_help->set_submenu ($menu_help);#set its sub menu
		$menu_bar->append($menu_item_help);	
	$table->attach_defaults($menu_bar,0,6,0,1);

	#add the second widget vpaned of vbox
	my $vpaned = Gtk2::VPaned->new;
	$vpaned->set_position (300);

	#add widget of notebook under vpaned
	my $notebook = Gtk2::Notebook->new;
		############ notebok: pre-processing 
		my $nb_sw=&nb_self();
		my $nb_label = Gtk2::Label->new_with_mnemonic ('Pre-processing');
		$notebook->append_page ($nb_sw, $nb_label);
		#### notebook: 4C-seq
		$nb_sw=&nb_4C();
		$nb_label = Gtk2::Label->new_with_mnemonic ('One-to-All');
		$notebook->append_page ($nb_sw, $nb_label);
		#### notebook: 3C-MTS or Capture-C
		$nb_sw=&nb_3C();
		$nb_label = Gtk2::Label->new_with_mnemonic ('Many-to-All');
		$notebook->append_page ($nb_sw, $nb_label);
		#### notebook: Hi-C
		$nb_sw=&nb_HiC();
		$nb_label = Gtk2::Label->new_with_mnemonic ('All-to-All');
		$notebook->append_page ($nb_sw, $nb_label);
	$notebook->show_all();
	$vpaned->add1 ($notebook);

	#add widget of frame-scrollwindow-textview under vpaned
	my $frame=Gtk2::Frame->new();
	$frame->set_label('Output viewing');
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $sw = Gtk2::ScrolledWindow->new ();
		$sw->set_shadow_type ('etched-out');
		$sw->set_policy ('automatic', 'automatic');
		$sw->set_border_width(10);
			my $text_view=Gtk2::TextView->new();
			$text_view->can_focus(1);
			$text_view->set_editable(0);
			$text_view->set_left_margin(5);
			$text_view->set_right_margin(5);
			$text_view->set_wrap_mode('GTK_WRAP_WORD_CHAR');
			#our $text_buffer = Gtk2::TextBuffer->new();
			#our $iter=$text_buffer->get_iter_at_offset(0);
			#$text_buffer->insert($iter, "Self-testing begins before running:\n\n\n");
			#$text_view->set_buffer($text_buffer);
			 #refresh screen per 1 second
		$text_view->show_all;
		$sw->add($text_view);
	$frame->add($sw);
	$vpaned->add2 ($frame);
	$vpaned->show_all();
	$table->attach_defaults ($vpaned, 0,6,1,8);

	#progress 
	$frame=Gtk2::Frame->new();
	$frame->set_label('Progress bar');
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $pbar = Gtk2::ProgressBar->new();
		$pbar->set_pulse_step(0.2);
		$pbar->hide; #needs to be called after show_all
	$frame->add($pbar);
	$table->attach_defaults ($frame, 0,3,8,9);
	
	$frame=Gtk2::Frame->new();
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $button_1=Gtk2::Button->new('Refresh output');
		$button_1->signal_connect(clicked=>sub{
				#my $content=`tac $variables{log_file}`; #refresh view from the end
				my $content=`cat $variables{log_file}`; #refresh view from the top
				#print "$variables{log_file}";
				my $text_buffer = $text_view->get_buffer();
				$text_buffer->set_text($content);
		});
	$frame->add($button_1);
	$table->attach_defaults ($frame, 3,4,8,9);
	
	$frame=Gtk2::Frame->new();
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $button_2=Gtk2::Button->new('Stop and quit');
		$button_2->signal_connect(clicked=>sub{
			$shash{'go'} =0;
			$shash{'die'} = 1;
			#foreach my $thread(threads->list()){
				#print "$thread\n";
				#$thread->join();
			#}
			$window->destroy;
			threads->exit;
			return TRUE;
		});
	$frame->add($button_2);
	$table->attach_defaults ($frame, 4,6,8,9);

$window->add($table);
$window->show_all;
Gtk2->main(); #our main event-loop

#######################################################

###########
#subroutines
####################################################################
#used in widget of menu
sub menu_tools{
	#Start of with a Gtk2::Menu widget
	my $menu_edit = Gtk2::Menu->new();
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		my $menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('Bowtie1 GUI', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			my $perl_script=$perl_dir."/tool_bowtie1.pl";
			system("perl $perl_script $variables{var_file}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut

	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		$menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('Bowtie2 GUI', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			my $perl_script=$perl_dir."/tool_bowtie2.pl";
			system("perl $perl_script $variables{var_file}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
	
	return($menu_edit);
}


###############################################################################
sub menu_view{
	#Start of with a Gtk2::Menu widget
	my $menu_edit = Gtk2::Menu->new();
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		my $menu_item = Gtk2::ImageMenuItem->new_from_stock ('QS Viewer: quality scores', undef);#add an Image Menu item using stock
		$menu_item->signal_connect('activate' => sub { 
				my $perl_script=$perl_dir."/viewer_QS.pl";
				system("perl $perl_script $variables{var_file}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item);#add first a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		$menu_item = Gtk2::ImageMenuItem->new_from_stock ('IL Viewer: insert length distribution', undef);#add an Image Menu item using stock
		$menu_item->signal_connect('activate' => sub { 
				my $perl_script=$perl_dir."/viewer_IL.pl";
				system("perl $perl_script $variables{var_file}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item);#add first a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		$menu_item = Gtk2::ImageMenuItem->new_from_stock ('SD Viewer: sequencing depth', undef);#add an Image Menu item using stock
		$menu_item->signal_connect('activate' => sub { 
				my $perl_script=$perl_dir."/viewer_SD.pl";
				system("perl $perl_script $variables{var_file}");
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item);#add first a item named cut
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
	return($menu_edit);
}

#used in widget of menu
sub menu_help{
	#Start of with a Gtk2::Menu widget
	my $menu_edit = Gtk2::Menu->new();
	
	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
		my $menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('Help document', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect(activate => sub { 

		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add first a item named cut

	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add  a tearoff item
		$menu_item_cut = Gtk2::ImageMenuItem->new_from_stock ('About', undef);#add an Image Menu item using stock
		$menu_item_cut->signal_connect('activate' => sub {
			#standard window
			my $window = sub_gui::default_window('About 3C_analyzer', 400, 300);
				$frame = Gtk2::Frame->new();
				$frame->set_border_width(5);
					my $table=Gtk2::Table->new(6,6, TRUE);
						my $label = Gtk2::Label->new("Version: 1.1");
					$table->attach_defaults($label, 0,6,0,1);
						$label = Gtk2::Label->new("Developer: Tiezheng Yuan, tyuan\@mcw.edu");
					$table->attach_defaults($label, 0,6,1,2);
						$label = Gtk2::Label->new("Copyright: 2012.09-2014.09 TiezhengYuan\n  3C_analyzer can be freely applied for non-commercial use!");
					$table->attach_defaults($label, 0,6,2,5);
						my $button=Gtk2::Button->new('Close');
					$table->attach_defaults($button, 2,4,5,6);
				$frame->add($table);
			$window->add($frame);
			$window->show_all();
			Gtk2->main();
		});#connet to the activate signal to catch when this item is selected
	$menu_edit->append($menu_item_cut);#add second: a item named cut

	$menu_edit->append(Gtk2::TearoffMenuItem->new);#add a tearoff item
	
	return($menu_edit);
}


#########################################
sub nb_self_basic_setup{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Directory setup', 800, 600);
	my $table = Gtk2::Table->new(5,6,TRUE);
		
		#frame: raw data dir
		my $frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Directory of raw data storage');
			my $subtable=Gtk2::Table->new(1,4, TRUE);
				my $entry_1 = Gtk2::Entry->new();
				$entry_1->set_text($variables{raw_data_dir});
			$subtable->attach_defaults($entry_1, 0,3,0,1);
				my $button_1=Gtk2::Button->new('Select a folder');
				$button_1->signal_connect('clicked' => sub{ 
					my $raw_data_dir=show_chooser('File Chooser type select-folder','select-folder');
					$entry_1->set_text($raw_data_dir);
				});
			$subtable->attach_defaults($button_1, 3,4,0,1);
		$frame->add($subtable);
	$table->attach_defaults($frame,0,6,0,1);
	
		#frame: result dir
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Directory of result storage');
			$subtable=Gtk2::Table->new(1,4, TRUE);
				my $entry_2 = Gtk2::Entry->new();
				$entry_2->set_text($variables{result_dir});
			$subtable->attach_defaults($entry_2, 0,3,0,1);
				my $button_2=Gtk2::Button->new('Select a folder');
				$button_2->signal_connect('clicked' => sub{ 
					my $result_dir=show_chooser('File Chooser type select-folder','select-folder');
					$entry_2->set_text($result_dir);
				});
			$subtable->attach_defaults($button_2, 3,4,0,1);
		$frame->add($subtable);
	$table->attach_defaults($frame,0,6,1,2);
	
		#frame:multiple threads number
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Number of the multi-threads');
			my $cb_threads_num= Gtk2::ComboBox->new_text;
			foreach (1,4,8,12,16,24,32,48,64){
				$cb_threads_num->append_text($_);
			}
			$cb_threads_num->set_active(2);
		$frame->add($cb_threads_num);
	$table->attach_defaults($frame,0,3,2,3);
	
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Length of size bins for co-localization (Kb)");
			my $entry_size_bin_len = Gtk2::Entry->new();
			$entry_size_bin_len->set_text($variables{size_bin_len}/1000);
		$frame->add($entry_size_bin_len);
	$table->attach_defaults($frame,3,6,2,3);
	
		#entry enzyme 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Primary restriction enzyme');
			my $entry_enzyme_site = Gtk2::Entry->new();
			$entry_enzyme_site->set_text($variables{enzyme_site});
		$frame->add($entry_enzyme_site);
	$table->attach_defaults($frame, 0,3,3,4);
	
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Number of enzymes of enzyme bins for co-localization");
			my $entry_enzyme_bin_enzymes = Gtk2::Entry->new();
			$entry_enzyme_bin_enzymes->set_text($variables{enzyme_bin_enzymes});
		$frame->add($entry_enzyme_bin_enzymes);
	$table->attach_defaults($frame,3,6,3,4);
	
			#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Shortest length of enzyme bins (Kb)");
			my $entry_enzyme_bin_lower = Gtk2::Entry->new();
			$entry_enzyme_bin_lower->set_text($variables{enzyme_bin_lower}/1000);
		$frame->add($entry_enzyme_bin_lower);
	$table->attach_defaults($frame,0,3,4,5);
	
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Longest length of enzyme bins (Kb)");
			my $entry_enzyme_bin_upper = Gtk2::Entry->new();
			$entry_enzyme_bin_upper->set_text($variables{enzyme_bin_upper}/1000);
		$frame->add($entry_enzyme_bin_upper);
	$table->attach_defaults($frame,3,6,4,5);
	
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Minimum range of colocalizations (bp)");
			my $entry_colocalization_range = Gtk2::Entry->new();
			$entry_colocalization_range->set_text($variables{colocalization_range});
		$frame->add($entry_colocalization_range);
	$table->attach_defaults($frame,0,3,5,6);
	
		#frame: check information display
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Checking information');
			my $textview = Gtk2::TextView->new();
			$textview->set_wrap_mode('GTK_WRAP_WORD_CHAR');
		$frame->add($textview);
	$table->attach_defaults($frame,0,6,6,7);
		
		#frame: check button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Check');
			$button->set_border_width(10);
			$button->signal_connect(pressed=>sub{
				my $textbuffer=Gtk2::TextBuffer->new();
				my $iter=$textbuffer->get_iter_at_offset(0);
				my $dir=$entry_1->get_text();
				my $raw_data_dir_judging= ( -d $dir) ? "Great! $dir used for raw data storage  exists.\n" : "Error: Wrong raw data directory!\n";
				$textbuffer->insert($iter, $raw_data_dir_judging);
				my $file_name_pointer=sub_basic::files_list($dir, 'incrusive_file_name');
				my @file_name=@$file_name_pointer;
				my @fastq_names=grep(/\.fastq$|\fq$/, @file_name);
				foreach (@fastq_names){
					$textbuffer->insert($iter, "The names $_ should contain R1 or R2.\n") unless $_=~/_R1_|_R2_/;
				}
				$dir=$entry_2->get_text();
				my $result_dir_judging= ( -d $dir) ? "Great! $dir used for result storage exists.\n" : "Error: Wrong result directory!\n";
				$textbuffer->insert($iter, $result_dir_judging);
				$textview->set_buffer($textbuffer);
				$textview->show_all();
			});
		$frame->add($button);
	$table->attach_defaults($frame,0,2,7,8);
		
		#frame: save button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			$button=Gtk2::Button->new('Save and close');
			$button->set_border_width(10);
			$button->signal_connect('clicked' => sub { 
				#print "$entry_1->get_text()\n", "$entry_2->get_text()\n";
				$variables{raw_data_dir}=$entry_1->get_text();
				$variables{result_dir}=$entry_2->get_text();
					$variables{site_info_file}=$variables{result_dir}.'/site_info.csv';
					$variables{var_file}=$variables{result_dir}.'/variables.txt';
					$variables{log_file}=$variables{result_dir}.'/output.log'; 
				$variables{threads_num}=$cb_threads_num->get_active_text();
				#size_bin_len
				$variables{size_bin_len}=$entry_size_bin_len->get_text()*1000;
				#enzyme
				$variables{enzyme_site}=uc $entry_enzyme_site->get_text();
				$variables{enzyme_bin_enzymes}=$entry_enzyme_bin_enzymes->get_text();
				$variables{enzyme_bin_lower}=$entry_enzyme_bin_lower->get_text()*1000;
				$variables{enzyme_bin_upper}=$entry_enzyme_bin_upper->get_text()*1000;
				$variables{colocalization_range}=$entry_colocalization_range->get_text();
				#initiate variables in variables.txt
				my %new_variables=%variables;
				delete $new_variables{var_file};
				delete $new_variables{site_info_file};
				delete $new_variables{log_file};
				sub_basic::hash_to_file(\%new_variables, $variables{var_file}, "=");
				#
				$window->destroy;
			} );
		$frame->add($button);
	$table->attach_defaults($frame,2,6,7,8);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

############################################
#Pops up a standard file chooser. # Specify a header to be displayed Specify a type depending on your needs 
#Optionally add a filter to show only certain files-will return a path, if valid

sub show_chooser {
	my($heading, $type, $filter) =@_;
	#$type can be:* 'open', * 'save', * 'select-folder' or * 'create-folder' 
	my $file_chooser =  Gtk2::FileChooserDialog->new ( 
                            $heading,
                            undef,
                            $type,
                            'gtk-cancel' => 'cancel',
                            'gtk-ok' => 'ok'
                        );
    (defined $filter)&&($file_chooser->add_filter($filter));
    
    #if action = 'save' suggest a filename
    ($type eq 'save')&&($file_chooser->set_current_name("suggeste_this_file.name"));

    my $filename;
    if ('ok' eq $file_chooser->run){    
       $filename = $file_chooser->get_filename;
    }

    $file_chooser->destroy;

    if (defined $filename){
        if ((-f $filename)&&($type eq 'save')) {
            my $overwrite =show_message_dialog( $window,
                                                'question'
                                                ,'Overwrite existing file:'."<b>\n$filename</b>"
                                                ,'yes-no'
                                    );
            return  if ($overwrite eq 'no');
        }
        return $filename;
    }
    return;
}
###########
###########################
sub CellRenderer_window{
	my ($matrix_pointer)=@_;
	my @matrix=@$matrix_pointer;
	my $row_num=@matrix;
	my $first_line_pointer=shift @matrix;
	my @first_line=@$first_line_pointer;
	my $col_num=@first_line;
	
	my @model_attr;
	for(my $i=0; $i<$col_num; $i++){
		$model_attr[$i]='Glib::String';
	}
	my $model = Gtk2::ListStore->new (@model_attr);
	#display content
	for(my $m=0; $m<$row_num; $m++){
		if($col_num==1){	$model->set ($model->append, 0 => $matrix[$m][0]);	}
		elsif($col_num==2){	$model->set ($model->append, 0 => $matrix[$m][0], 1=>$matrix[$m][1]);	}
		else{	$model->set ($model->append, 0 => $matrix[$m][0], 1=>$matrix[$m][1], 2=>$matrix[$m][2]);	}
	}
	
	#add column
	my $view = Gtk2::TreeView->new ($model);
	for(my $i=0; $i<$col_num; $i++){
		my $renderer = Gtk2::CellRendererText->new;
		my $column = Gtk2::TreeViewColumn->new_with_attributes ($first_line[$i], $renderer, text => $i);
		$view->append_column ($column);
	}
	
	#scroller window
	my $scroller = Gtk2::ScrolledWindow->new;
	$scroller->set_policy (qw(automatic automatic));
	$scroller->add ($view);
	
	#checkbutton
	my $check = Gtk2::CheckButton->new ('resizable columns');
	$check->set_active (FALSE);
	$check->signal_connect (toggled => sub {
		map { $_->set_resizable ($check->get_active); } $view->get_columns;
	});

	my $box = Gtk2::VBox->new;
	$box->add ($scroller);
	$box->pack_start ($check, FALSE, FALSE, 0);

	return($box);
}

##################################
sub nb_self_sample_management{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Sample management', 700, 400);
	my $table = Gtk2::Table->new(7,4,TRUE);
	
	#get sample names and fastq files
	my %sample_hash;
	my $sample_names_pointer=sub_basic::auto_sample_names($variables{raw_data_dir}, $variables{raw_file_format} );
	my @sample_names=@$sample_names_pointer;
	my $files_pointer=sub_basic::files_list($variables{raw_data_dir}, 'incrusive_file');
	my @files=@$files_pointer;
	my @raw_files=grep(/\.fastq$|\.fq$/, @files);
	
		#the list of sample names and fastq names
		my $vbox = Gtk2::VBox->new;
			my $scroll = Gtk2::ScrolledWindow->new;
			$scroll->set_policy ('never', 'automatic');

			# create and load the model
			my $liststore = Gtk2::ListStore->new ('Glib::String', 'Glib::String');
			foreach (@sample_names) {
				my $iter = $liststore->append;
				$liststore->set ($iter, 0, $_, 1, $_);
			}
			# now a view
			my $treeview_1 = Gtk2::TreeView->new ($liststore);
			$treeview_1->set_rules_hint (TRUE);
			$treeview_1->set_reorderable (TRUE);
				# regular editable text column for column 0
				my $renderer = Gtk2::CellRendererText->new;
				$renderer->set (editable => FALSE);
				my $column = Gtk2::TreeViewColumn->new_with_attributes ('Fastq names', $renderer, text => 0);
				$renderer->signal_connect (edited => sub {
							my ($cell, $text_path, $new_text, $model) = @_;
							my $path = Gtk2::TreePath->new_from_string ($text_path);
							my $iter = $model->get_iter ($path);
							$model->set ($iter, 0, $new_text);
						}, $liststore);
			$treeview_1->append_column ($column);
				#renderer column 1
				$renderer = Gtk2::CellRendererText->new;
				$renderer->set (editable => TRUE);
				$column = Gtk2::TreeViewColumn->new_with_attributes ('Sample names', $renderer, text => 1);
				$renderer->signal_connect (edited => sub {
							my ($cell, $text_path, $new_text, $model) = @_;
							my $path = Gtk2::TreePath->new_from_string ($text_path);
							my $iter = $model->get_iter ($path);
							$model->set ($iter, 1, $new_text);
						}, $liststore);
			$treeview_1->append_column ($column);

			$scroll->add ($treeview_1);
		$vbox->pack_start ($scroll, TRUE, TRUE, 0);
	$table->attach_defaults($vbox,0,3,0,3);
		#frame: check sample_info
		$frame= Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Checking information');
			my $sw = Gtk2::ScrolledWindow->new ();
			$sw->set_shadow_type ('etched-out');
			$sw->set_policy ('automatic', 'automatic');
			$sw->set_border_width(10);
				my $textview_2 = Gtk2::TextView->new();
				$textview_2->set_wrap_mode('GTK_WRAP_WORD_CHAR');
			$sw->add($textview_2);
		$frame->add($sw);
	$table->attach_defaults($frame,0,4,3,7);		
		#button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_1=Gtk2::Button->new('Check');
			$button_1->set_border_width(5);
			$button_1->signal_connect(clicked=>sub{
						%sample_hash=(); #empty out hash
						$liststore->foreach(sub{
							my($model, $path, $iter)=@_;
							my $fastq=$model->get ($iter, 0);
							my $sample=$model->get ($iter, 1);
							$sample_hash{$sample}= $sample_hash{$sample} ? $sample_hash{$sample}.','.$fastq : $fastq;
							FALSE;
						});
						#connection between sample name and fastq files
						my $textbuffer=Gtk2::TextBuffer->new();
						my $iter=$textbuffer->get_iter_at_offset(0);
						my $n=1;
						foreach my $sample_name(keys %sample_hash){#4
							my @raw_names=split(",", $sample_hash{$sample_name});
							my @sub_raw_files;
							foreach my $raw_name(@raw_names){
								foreach my $raw_file(@raw_files){
									push(@sub_raw_files, $raw_file) if $raw_file=~/\/$raw_name/;
								}
							}
							my $files_num=@sub_raw_files;
							$textbuffer->insert($iter, $n.": $sample_name (number of fastq files=$files_num):\n");
							foreach(@sub_raw_files){
								$textbuffer->insert($iter, "\t$_\n");
							}
							$n++;
						}#4
					$textview_2->set_buffer($textbuffer);
					$textview_2->show_all();
				});
		$frame->add($button_1);
	$table->attach_defaults($frame,3,4,0,1);
		#button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_2=Gtk2::Button->new('Save');
			$button_2->set_border_width(5);
			$button_2->signal_connect(pressed=>sub{
							my $sample_info_file=$variables{result_dir}.'/sample_info.csv';
							sub_basic::hash_to_file(\%sample_hash, $sample_info_file, ",");
						});
		$frame->add($button_2);
	$table->attach_defaults($frame,3,4,1,2);
		#button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_3=Gtk2::Button->new('Close');
			$button_3->set_border_width(5);
			$button_3->signal_connect(pressed=>sub{	$window->destroy;		});
		$frame->add($button_3);
	$table->attach_defaults($frame,3,4,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

##################################
sub nb_3C_viewpoints_locking{
	#standard window creation, placement, and signal connecting
	my $title='Primary enzyme sites ( '.$variables{enzyme_site}.' ) of reference genome sequences';
	my $window = sub_gui::default_window($title, 1000, 600);
	my $table = Gtk2::Table->new(7,4,TRUE);
	
	#get sample names and fastq files
	my $sample_names_pointer=sub_basic::auto_sample_names($variables{raw_data_dir}, $variables{raw_file_format} );
	my @sample_names=@$sample_names_pointer;
	my $files_pointer=sub_basic::files_list($variables{raw_data_dir}, 'incrusive_file');
	my @files=@$files_pointer;
	my @fastq_files=grep(/\.fastq$|\.fq$/, @files);
	
		#the list of sample names and fastq names
		my $frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $scroll = Gtk2::ScrolledWindow->new;
			$scroll->set_policy ('never', 'automatic');

			## create and load the model
			#my @site_info_items=('Site name', 'Chromosome', 'Site position', 'Start territory', 'End territory', 'Upstream length', 'downstream length');
			#my $liststore = Gtk2::ListStore->new ('Glib::String', 'Glib::String', 'Glib::Int', 'Glib::Int','Glib::Int', 'Glib::Int', 'Glib::Int');
			#my $treeview = sub_gui::column_treeview($liststore, \@site_info_items);   #subroutine
			#$scroll->add ($treeview);
				my $slist = Gtk2::SimpleList->new ('Target sites'=> 'bool',
					'Site names'=> 'text',  'Site position'=>'int', 'Site seq'=>'text',
					'Upstream detection'=> 'bool', 'Downstream detection'=> 'bool',
				);

				$slist->get_selection->set_mode ('multiple');
				# simple way to make text columns editable
				$slist->set_column_editable (1, TRUE);
				$slist->set_rules_hint (TRUE);
				$slist->set_reorderable (TRUE);
				map { $_->set_resizable (TRUE) } $slist->get_columns;
			$scroll->add ($slist);
			$scroll->show_all;
		$frame->add($scroll);
	$table->attach_defaults($frame,0,3,0,10);
		
		#frame: enzyme sites file selection
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Genome sequences in FASTA');
			my $file_chooser_genome =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_genome->set_filename($variables{ref_dir});
			$file_chooser_genome->hide;
		$frame->add($file_chooser_genome);
	$table->attach_defaults($frame,3,4,0,1);
		#frame:
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Chromosome');
			my $entry_chr = Gtk2::Entry->new();
			$entry_chr->set_text('8');
			$frame->add($entry_chr);
	$table->attach_defaults($frame, 3,4,1,2);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Start of genomic region');
			my $entry_start_region = Gtk2::Entry->new();
			$entry_start_region->set_text(128309266);
		$frame->add($entry_start_region);
	$table->attach_defaults($frame, 3,4,2,3);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('End of genomic region');
			my $entry_end_region = Gtk2::Entry->new();
			$entry_end_region->set_text(128552222);
		$frame->add($entry_end_region);
	$table->attach_defaults($frame, 3,4,3,4);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Length of viewpoints');
			my $entry_region_length = Gtk2::Entry->new();
			$entry_region_length->set_text(500);
		$frame->add($entry_region_length);
	$table->attach_defaults($frame, 3,4,4,5);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Expanding of genomic region');
			my $entry_expanding_len = Gtk2::Entry->new();
			$entry_expanding_len->set_text(1000);
		$frame->add($entry_expanding_len);
	$table->attach_defaults($frame, 3,4,5,6);
	
		#button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_1=Gtk2::Button->new('Analyze chromosome');
			$button_1->set_border_width(5);
			$button_1->signal_connect(clicked=>sub{
				my $site_chr=$entry_chr->get_text();
				my $start_region=$entry_start_region->get_text();
				my $end_region=$entry_end_region->get_text();
				my $region_length=$entry_region_length->get_text();
				$variables{genome_fasta_file}=$file_chooser_genome->get_filename;
				
				@{$slist->{data}}=(); #clear slist
				my $judging=0;
				$judging=2 unless -f $variables{genome_fasta_file};
				$judging=3 if $start_region>=$end_region;
				if ($judging==0){ #judging circyle
					my $in_obj = Bio::SeqIO->new(-file => $variables{genome_fasta_file}, -format => 'fasta');
					while (my $seq_obj = $in_obj->next_seq() ) {
						my $chr=$seq_obj->display_id;
						if($chr eq $site_chr){
							my $chr_seq=$seq_obj->seq();
							my $chr_len=$seq_obj->length;
							my $territory_seq=$seq_obj->subseq($start_region, $end_region);
							my $enzyme_sites_str=sub_3C::enzyme_sites_counting($territory_seq, $variables{enzyme_site}, $start_region);#subrountine
							my @enzyme_sites=split(",", $enzyme_sites_str);
							my $enzyme_num=@enzyme_sites;
							for(my $i=0; $i<@enzyme_sites; $i++){
								my $site_name=$chr.'_E'.($i+1);
								my $site_position=$enzyme_sites[$i];
								my $site_seq=substr($chr_seq, $site_position-1, 20);
								my @enzyme_info=(TRUE, $site_name, $site_position, $site_seq, TRUE, TRUE);
								push @{$slist->{data}}, \@enzyme_info;
							}
							sub_gui::popup_window('Notice!', 'Enzyme sites are listed on the left scrolled window!');
							$judging=1;
							$judging=4 if $start_region>=$chr_len or $end_region>=$chr_len or $start_region>=$end_region;
							$judging=5 if @enzyme_sites==0;
							last;
						}
					}
				}#judging circyle
				
				sub_gui::popup_window('Error!', 'No a selected chromosome or wrong chromosome input!') if $judging ==0;
				sub_gui::popup_window('Error!', 'Please select genome fasta file!') if $judging==2;
				sub_gui::popup_window('Error!', "Wrong! Start territory should be beyond end territory!") if $judging ==3;
				sub_gui::popup_window('Error!', "Wrong territory input of chromosome $site_chr!") if $judging ==4;
				sub_gui::popup_window('Error!', "No enzyme ($variables{enzyme_site}) are deteced!") if $judging ==5;
			});
		$frame->add($button_1);
	$table->attach_defaults($frame,3,4,6,7);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_2=Gtk2::Button->new("Refresh enzyme sites");
			$button_2->set_border_width(5);
			$button_2->signal_connect(clicked=>sub{
						my $site_chr=$entry_chr->get_text();
						my $boundary_expanding_len=$entry_expanding_len->get_text();
						my $region_length=$entry_region_length->get_text();
						$variables{genome_fasta_file}=$file_chooser_genome->get_filename;
						
						#export
						my @selected_pos;
						foreach my $pointer(@{$slist->{data}}){
							my @items=@$pointer;
							push(@selected_pos, $items[2]) if $items[0]==1;  #only selected items
						}
						my $start_territory=$selected_pos[0]-$boundary_expanding_len;
						my $end_territory=$selected_pos[-1]+$boundary_expanding_len;
						open my($OUT), ">", $variables{site_info_file} or die;
						foreach my $pointer(@{$slist->{data}}){
							my @items=@$pointer;
							if($items[0]==1){ #only selected items
								my $site_name=$items[1];
								my $site_position=$items[2];
								my $upstream=($items[4] ==1) ? $region_length : 0;
								my $downstream=($items[5] ==1) ? $region_length : 0;
								print $OUT join(',', $site_name, $site_chr, $site_position, $start_territory, $end_territory, $upstream, $downstream), "\n";
							}
						}
						close($OUT);
						sub_gui::popup_window('Notice!', "Enzyme sites were exported into $variables{site_info_file}!");
				});
		$frame->add($button_2);
	$table->attach_defaults($frame,3,4,7,8);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_3=Gtk2::Button->new("Append enzyme sites");
			$button_3->set_border_width(5);
			$button_3->signal_connect(clicked=>sub{
						my $site_chr=$entry_chr->get_text();
						my $boundary_expanding_len=$entry_expanding_len->get_text();
						my $region_length=$entry_region_length->get_text();
						$variables{genome_fasta_file}=$file_chooser_genome->get_filename;
						
						#export
						my @selected_pos;
						foreach my $pointer(@{$slist->{data}}){
							my @items=@$pointer;
							push(@selected_pos, $items[2]) if $items[0]==1;  #only selected items
						}
						my $start_territory=$selected_pos[0]-$boundary_expanding_len;
						my $end_territory=$selected_pos[-1]+$boundary_expanding_len;
						open my($OUT), ">>", $variables{site_info_file} or die;
						foreach my $pointer(@{$slist->{data}}){
							my @items=@$pointer;
							if($items[0]==1){ #only selected items
								my $site_name=$items[1];
								my $site_position=$items[2];
								my $upstream=($items[4] ==1) ? $region_length : 0;
								my $downstream=($items[5] ==1) ? $region_length : 0;
								print $OUT join(',', $site_name, $site_chr, $site_position, $start_territory, $end_territory, $upstream, $downstream), "\n";
							}
						}
						close($OUT);
						sub_gui::popup_window('Notice', "Enzyme sites were appended into $variables{site_info_file}!");
				});
		$frame->add($button_3);
	$table->attach_defaults($frame,3,4,8,9);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_4=Gtk2::Button->new('Close');
			$button_4->set_border_width(5);
			$button_4->signal_connect(pressed=>sub{	$window->destroy;		});
		$frame->add($button_4);
	$table->attach_defaults($frame,3,4,9,10);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}


######################################
sub nb_self_package_dependency{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Package dependency checking',600,600);
	my $table = Gtk2::Table->new(6,3,TRUE);
	
		my $frame_1 = Gtk2::Frame->new();
		$frame_1->set_border_width(5);
		$frame_1->set_label('Checking information');
			my $sw = Gtk2::ScrolledWindow->new ();
			$sw->set_shadow_type ('etched-out');
			$sw->set_policy ('automatic', 'automatic');
			$sw->set_border_width(5);
				my $textview = Gtk2::TextView->new();
				$textview->set_left_margin(5);
				$textview->set_wrap_mode('GTK_WRAP_WORD_CHAR');
			$sw->add($textview);
		$frame_1->add($sw);
	$table->attach_defaults($frame_1,0,3,0,5);
			
		my $frame_2 = Gtk2::Frame->new();
		$frame_2->set_border_width(5);
			my $hbox = Gtk2::HBox->new(FALSE,5);
				my $button_1=Gtk2::Button->new('Check');
				$button_1->set_border_width(10);
				$button_1->signal_connect(pressed=>sub{
					my $textbuffer=Gtk2::TextBuffer->new();
					my $iter=$textbuffer->get_iter_at_offset(0);
					#Perl packages
					$textbuffer->insert($iter, "Perl package dependency:\n");
					my @ISA=qw(Cwd File::Find Glib Gtk2 List::Util List::MoreUtils threads threads::shared);
					foreach (sort @ISA){
						$textbuffer->insert($iter, (eval "require $_") ? "\tOK\t"."$_\n" : "\tNO\t"."$_\n");
					}
					#the third-party tools
					$textbuffer->insert($iter, "\nThe third-party tools:\n");
					$textbuffer->insert($iter, (-f $variables{alignment_dir}.'/bowtie2-align') ? "\tOK\t" : "\tNO\t");
					$textbuffer->insert($iter, "Bowtie2 (The aligner)\n");
					$textbuffer->insert($iter, (-f $variables{alignment_dir}.'/bowtie2-build') ? "\tOK\t" : "\tNO\t");
					$textbuffer->insert($iter, "Bowtie2-build (Index building for alignments)\n");
					$textbuffer->insert($iter, (-f $variables{alignment_dir}.'/samtools') ? "\tOK\t" : "\tNO\t");
					$textbuffer->insert($iter, "samtools (SAM/BAM format conversion)\n");
					
					$textview->set_buffer($textbuffer);
					$textview->show_all();
				});
			$hbox->pack_start($button_1,TRUE,TRUE,10);
				my $button_2=Gtk2::Button->new('Close');
				$button_2->set_border_width(10);
				$button_2->signal_connect(clicked => sub { $window->destroy;}  );
			$hbox->pack_start($button_2,TRUE,TRUE,10);
		$frame_2->add($hbox);
	$table->attach_defaults($frame_2,0,3,5,6);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}
#############################################
sub nb_self{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request ('200', '200');

	my $table_stock=Gtk2::Table->new(3, 5, TRUE);
	$table_stock->set_border_width(30);
	$table_stock->set_row_spacings(30);
	$table_stock->set_size_request (100, 100);
	#1
	$table_stock->attach_defaults(&sub_gui::nb_button('Basic setup', \&nb_self_basic_setup), 0,2,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 2,3,0,1);
	#2
	$table_stock->attach_defaults(&sub_gui::nb_button('Sample management', \&nb_self_sample_management), 3,5,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('down', 'etched-out'), 4,5,1,2);
	#3
	$table_stock->attach_defaults(&sub_gui::nb_button('Package dependency checking', \&nb_self_package_dependency), 3,5,2,3);
	#$table_stock->attach_defaults(Gtk2::Arrow->new('left', 'etched-out'), 2,3,2,3);
	
	$scroller->add_with_viewport($table_stock);

	return($scroller);
}


######################################################################
sub nb_3C{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request (200, 200);

	my $table_stock=Gtk2::Table->new(3, 8, TRUE);
	$table_stock->set_border_width(30);
	$table_stock->set_row_spacings(30);
	$table_stock->set_size_request (100, 100);
	#1
	$table_stock->attach_defaults(&sub_gui::nb_button('Viewpoints Locking', \&nb_3C_viewpoints_locking), 0,2,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 2,3,0,1); #add an arrow
	#2
	$table_stock->attach_defaults(&sub_gui::nb_button('FASTQ trimming', \&nb_3C_fastq_trimming), 3,5,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 5,6,0,1);#add an arrow
	#3
	$table_stock->attach_defaults(&sub_gui::nb_button('Co-localization detection', \&nb_3C_colocalization_detection), 6,8,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('down', 'etched-out'), 6,8,1,2);#add an arrow
	#4
	$table_stock->attach_defaults(&sub_gui::nb_button('Co-localization counting', \&nb_3C_colocalization_counting), 6,8,2,3);
	$table_stock->attach_defaults(Gtk2::Arrow->new('left', 'etched-out'), 5,6,2,3);#add an arrow
	my $tbutton=Gtk2::Button->new('Apply and Run');
	$tbutton->signal_connect(clicked=>sub{
		my $thread_monitor=threads->new( sub{
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir.'/time_monitor.pl';
					system("perl $perl_script $variables{var_file}");
				}
				else{	sleep 1;	}
			}
		});

		my $thread_main=threads->new( sub{
			#local $SIG{KILL} = sub { threads->exit };
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir."/main_3CMTS.pl";
					system("perl $perl_script $variables{var_file} > $variables{log_file} 2>&1");
					$shash{'go'} = 0; #turn off main running 
				}
				else{	sleep 1;	}
			}
		});

		if($shash{'go'} == 0){
			$shash{'go'} = 1;
			$pbar->show;
			Glib::Timeout->add (100, sub {
				if($shash{'go'} == 1){
					$pbar->set_text('Running!');
					$pbar->pulse;
					return TRUE;
				}
				else{	
					$pbar->set_text('OK! It is done!');
					return FALSE;
				}
			});
		}
	});  
	$table_stock->attach_defaults($tbutton, 2,5,2,3);
	$scroller->add_with_viewport($table_stock);
	return($scroller);
}


#############################
sub nb_3C_fastq_trimming{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options for adapter trimming', 500, 300);
	my $table = Gtk2::Table->new(3,4,TRUE);
		
		#frame: 3' adapter sequences
		$frame=Gtk2::Frame->new("3\' end adapter sequences");
		$frame->set_border_width(5);
			my $entry_adapter_3 = Gtk2::Entry->new();
			$entry_adapter_3->set_text($variables{adapter_3});
		$frame->add($entry_adapter_3);
	$table->attach_defaults($frame, 0,4,0,1);
		#frame: 5' adapter sequences
		$frame=Gtk2::Frame->new("5\' end adapter sequences");
		$frame->set_border_width(5);
			my $entry_adapter_5 = Gtk2::Entry->new();
			$entry_adapter_5->set_text($variables{adapter_5});
		$frame->add($entry_adapter_5);
	$table->attach_defaults($frame, 0,4,1,2);
		#frame: adapter length
		$frame=Gtk2::Frame->new('Minimum length of the adatper');
		$frame->set_border_width(5);
			my $cb_min_adapter = Gtk2::ComboBox->new_text;
			foreach (8..16){
				$cb_min_adapter->append_text($_);
			}
			$cb_min_adapter->set_active(6);
		$frame->add($cb_min_adapter);
	$table->attach_defaults($frame, 0,2,2,3);
			#frame: min read seq
		$frame=Gtk2::Frame->new('Minimum length left in trimming');
		$frame->set_border_width(5);
			my $cb_min_trim_len = Gtk2::ComboBox->new_text;
			foreach (8..20){
				$cb_min_trim_len->append_text($_);
			}
			$cb_min_trim_len->set_active(2);
		$frame->add($cb_min_trim_len);
	$table->attach_defaults($frame, 2,4,2,3);

		#frame: check button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->set_border_width(5);
			$button->signal_connect(clicked=>sub{
				$variables{adapter_3}=$entry_adapter_3->get_text();
				$variables{adapter_5}=$entry_adapter_5->get_text();
				$variables{adapter_len}=$cb_min_adapter->get_active_text;
				$variables{min_trim_len}=$cb_min_trim_len->get_active_text;
				$variables{'3C_trimming'}='yes';
				#refresh log file
				sub_basic::refresh_log($variables{var_file}, 'adapter_3', $variables{'adapter_3'});
				sub_basic::refresh_log($variables{var_file}, 'adapter_5', $variables{'adapter_5'});
				sub_basic::refresh_log($variables{var_file}, 'adapter_len', $variables{'adapter_len'});
				sub_basic::refresh_log($variables{var_file}, 'min_trim_len', $variables{'min_trim_len'});
				sub_basic::refresh_log($variables{var_file}, '3C_trimming', $variables{'3C_trimming'});
				$window->destroy;
			});
		$frame->add($button);
	$table->attach_defaults($frame,0,4,3,4);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

#######################
sub nb_3C_colocalization_detection{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of genome mapping (Bowtie 2)', 1000, 500);
	my $table=Gtk2::Table->new(6,6,TRUE);
	
		#genome_fasta
		my $file_names_pointer=sub_basic::files_list($variables{ref_dir}, 'file_name');
		my @file_names=@$file_names_pointer;
		my @fa_names=grep(/\.fa$|\.fasta$/, @file_names);
		#frame: select genome sequence (*.fa)
		$frame=Gtk2::Frame->new('Genome sequences (.fa file)');
		$frame->set_border_width(5);
			my $cb_fasta_file = Gtk2::ComboBox->new_text;
			foreach (@fa_names){
				$cb_fasta_file->append_text($_);
			}
			$cb_fasta_file->set_active(0);
		$frame->add($cb_fasta_file);
	$table->attach_defaults($frame, 0,3,0,1);
	
		#Frame: alignment mode --end-to-end/--local
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Alignment mode");
			my $cb_alignment_mode = Gtk2::ComboBox->new_text;
			$cb_alignment_mode->append_text('--end-to-end --very-fast');
			$cb_alignment_mode->append_text('--end-to-end --fast');
			$cb_alignment_mode->append_text('--end-to-end --sensitive');
			$cb_alignment_mode->append_text('--end-to-end --very-sensitive');
			$cb_alignment_mode->set_active(3);
		$frame->add($cb_alignment_mode);
	$table->attach_defaults($frame, 3,6,0,1);
		#Frame: mapping quality 
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Cutoff for mapping quality");
			my $cb_mapping_quality = Gtk2::ComboBox->new_text;
			foreach(1,10,13,20){
				$cb_mapping_quality->append_text($_);
			}
			$cb_mapping_quality->set_active(1);
		$frame->add($cb_mapping_quality);
	$table->attach_defaults($frame, 0,3,1,2);

		#Frame: --np
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Score of ambiguous penality (--np)");
			my $cb_ambiguous_penalty= Gtk2::ComboBox->new_text;
			foreach (1..4){
				$cb_ambiguous_penalty->append_text($_);
			}
			$cb_ambiguous_penalty->set_active(0);
		$frame->add($cb_ambiguous_penalty);
		$table->attach_defaults($frame, 3,6,1,2);
		#Frame: --rdg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the read gap penalities (--rdg)");
			my $cb_read_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_read_gap_penalty->append_text('3,2');
			$cb_read_gap_penalty->append_text('5,3');
			$cb_read_gap_penalty->append_text('8,4');
			$cb_read_gap_penalty->set_active(1);
		$frame->add($cb_read_gap_penalty);
		$table->attach_defaults($frame, 0,3,2,3);
		#Frame: --rfg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the reference gap penalities (--rfg)");
			my $cb_ref_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_ref_gap_penalty->append_text('3,2');
			$cb_ref_gap_penalty->append_text('5,3');
			$cb_ref_gap_penalty->append_text('8,4');
			$cb_ref_gap_penalty->set_active(1);
		$frame->add($cb_ref_gap_penalty);
		$table->attach_defaults($frame, 3,6,2,3);
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
		$table->attach_defaults($frame, 0,3,3,4);
		#Frame: -I
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Minimum fragment length for valid paired-end alignments (-I)");
			my $cb_min_insert = Gtk2::ComboBox->new_text;
			foreach (0,20,40,60,80,100){
				$cb_min_insert->append_text($_);
			}
			$cb_min_insert->set_active(0);
		$frame->add($cb_min_insert);
	$table->attach_defaults($frame,3,6,3,4);
		#Frame: -X
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Maximum fragment length for valid paired-end alignments (-X)");
			my $cb_max_insert = Gtk2::ComboBox->new_text;
			foreach (100, 200,300,400,500,600){
				$cb_max_insert->append_text($_);
			}
			$cb_max_insert->set_active(0);
		$frame->add($cb_max_insert);
	$table->attach_defaults($frame,0,3,4,5);
		#Frame:--mp 
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of maximum and minimum mismatch penalties (--mp)");
			my $cb_mismatch_penality = Gtk2::ComboBox->new_text;
			$cb_mismatch_penality->append_text('4,2');
			$cb_mismatch_penality->append_text('6,2');
			$cb_mismatch_penality->append_text('6,3');
			$cb_mismatch_penality->set_active(1);
		$frame->add($cb_mismatch_penality);
	$table->attach_defaults($frame, 3,6,4,5);
	
	#Frame
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	my $button=Gtk2::Button->new('Save and close');
	$button->signal_connect(clicked=>sub{
		my $genome_fasta_file=$cb_fasta_file->get_active_text;
		$variables{genome_fasta_file}=$variables{ref_dir}.'/'.$genome_fasta_file;
		$variables{genome_index_name}=$genome_fasta_file;
		$variables{genome_index_name}=~s/\.fa$|\.fasta$//;
		$variables{genome_index}=$variables{alignment_dir}.'/'.$variables{genome_index_name};
		$variables{mapping_quality}=$cb_mapping_quality->get_active_text();
		#refresh log file
		sub_basic::refresh_log($variables{var_file}, 'genome_fasta_file', $variables{genome_fasta_file}) ;
		sub_basic::refresh_log($variables{var_file}, 'genome_index', $variables{genome_index}) ;
		sub_basic::refresh_log($variables{var_file}, 'genome_index_name', $variables{genome_index_name});
		sub_basic::refresh_log($variables{var_file}, 'mapping_quality', $variables{mapping_quality}) ;
		#bowtie options
		if (exists $variables{genome_fasta_file} and -f $variables{genome_fasta_file}){
			$variables{bowtie_options}=join(" ", 
							$variables{alignment_dir}.'/bowtie2',  #bowtie executable script
							$cb_alignment_mode->get_active_text(), 
							'--mp', $cb_mismatch_penality->get_active_text(),
							'--np', $cb_ambiguous_penalty->get_active_text(),
							'--rdg', $cb_read_gap_penalty->get_active_text(),
							'--rfg', $cb_ref_gap_penalty->get_active_text(),
							'-I', $cb_min_insert->get_active_text(),
							'-X', $cb_max_insert->get_active_text(),
							'--no-hd', '-x', $variables{genome_index}, #bowtie index
						);
			$variables{'3C_colocalization'}='yes';
			sub_basic::refresh_log($variables{var_file}, '3C_colocalization', $variables{'3C_colocalization'} );
			sub_basic::refresh_log($variables{var_file}, '3C_other_detection', 'yes' );
			sub_basic::refresh_log($variables{var_file}, 'bowtie_options', $variables{bowtie_options} );
			$window->destroy;
		}
		else{
			sub_gui::popup_window('Error!', "Please select a fasta file of the reference genome in Setup!");
		}

	});
	$button->set_border_width(5);
	$frame->add($button);
	$table->attach_defaults($frame,1,5,5,6);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}

#######################
sub nb_3C_colocalization_counting{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of co-localization counting', 400, 200);
	my $table=Gtk2::Table->new(3,6,TRUE);
		
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Noise of read counts");
			my $cb_reads_noise = Gtk2::ComboBox->new_text;
			foreach (1,2,3,4,5,10,20){
				$cb_reads_noise->append_text($_);
			}
			$cb_reads_noise->set_active(0);
		$frame->add($cb_reads_noise);
	$table->attach_defaults($frame,0,6,0,1);
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Threshold of probability");
			my $cb_cumulative_p = Gtk2::ComboBox->new_text;
			foreach(0.05,0.01,0.001){
				$cb_cumulative_p->append_text($_);
			}
			$cb_cumulative_p->set_active(1);
		$frame->add($cb_cumulative_p);
	$table->attach_defaults($frame,0,6,1,2);

		#Frame
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->signal_connect(clicked=>sub{
				$variables{RC_noise}=$cb_reads_noise->get_active_text();
				$variables{RC_cumulative_p}=$cb_cumulative_p->get_active_text();
				$variables{'3C_statistics'}='yes';
				$variables{'3C_RC_counting'}='yes';
				$variables{'3C_Tscore'}='yes';
				
				sub_basic::refresh_log($variables{var_file}, 'RC_noise', $variables{RC_noise} );
				sub_basic::refresh_log($variables{var_file}, 'RC_cumulative_p', $variables{RC_cumulative_p} );
				sub_basic::refresh_log($variables{var_file}, '3C_statistics', $variables{'3C_statistics'});
				sub_basic::refresh_log($variables{var_file}, '3C_RC_counting', $variables{'3C_RC_counting'});
				sub_basic::refresh_log($variables{var_file}, '3C_Tscore', $variables{'3C_Tscore'});
				$window->destroy;
			});
		$frame->add($button);
	$table->attach_defaults($frame,1,5,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}

######################################################################
sub nb_4C{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request (200, 200);

	my $table_stock=Gtk2::Table->new(3, 8, TRUE);
	$table_stock->set_border_width(30);
	$table_stock->set_row_spacings(30);
	$table_stock->set_size_request (100, 100);
	#1
	$table_stock->attach_defaults(&sub_gui::nb_button('Viewpoints Locking', \&nb_4C_viewpoints_locking), 0,3,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 3,5,0,1); #add an arrow
	#2
	#$table_stock->attach_defaults(&sub_gui::nb_button('FASTQ trimming', \&nb_4C_fastq_trimming), 3,5,0,1);
	#$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 5,6,0,1);#add an arrow
	#3
	$table_stock->attach_defaults(&sub_gui::nb_button('Co-localization detection', \&nb_4C_colocalization_detection), 5,8,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('down', 'etched-out'), 6,8,1,2);#add an arrow
	#4
	$table_stock->attach_defaults(&sub_gui::nb_button('Co-localization counting', \&nb_4C_colocalization_counting), 6,8,2,3);
	$table_stock->attach_defaults(Gtk2::Arrow->new('left', 'etched-out'), 5,6,2,3);#add an arrow
	#run
	my $tbutton=Gtk2::Button->new('Apply and Run');
	$tbutton->signal_connect(clicked=>sub{
		my $thread_monitor=threads->new( sub{
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir.'/time_monitor.pl';
					system("perl $perl_script $variables{var_file}");
				}
				else{	sleep 1;	}
			}
		});

		my $thread_main=threads->new( sub{
			#local $SIG{KILL} = sub { threads->exit };
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir."/main_4C.pl";
					system("perl $perl_script $variables{var_file} > $variables{log_file} 2>&1");
					$shash{'go'} = 0; #turn off main running 
				}
				else{	sleep 1;	}
			}
		});

		if($shash{'go'} == 0){
			$shash{'go'} = 1;
			$pbar->show;
			Glib::Timeout->add (100, sub {
				if($shash{'go'} == 1){
					$pbar->set_text('Running!');
					$pbar->pulse;
					return TRUE;
				}
				else{	
					$pbar->set_text('OK! It is done!');
					return FALSE;
				}
			});
		}
	});  
	$table_stock->attach_defaults($tbutton, 2,5,2,3);
	$scroller->add_with_viewport($table_stock);
	return($scroller);
}

#############################
sub nb_4C_fastq_trimming{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options for adapter trimming', 500, 250);
	my $table = Gtk2::Table->new(3,4,TRUE);
		
		#frame: 3' adapter sequences
		$frame=Gtk2::Frame->new("P5 Illumina adapter sequences");
		$frame->set_border_width(5);
			my $entry_adapter_P5 = Gtk2::Entry->new();
			$entry_adapter_P5->set_text($variables{adapter_P5});
		$frame->add($entry_adapter_P5);
	$table->attach_defaults($frame, 0,4,0,1);
		#frame: adapter length
		$frame=Gtk2::Frame->new('Minimum length of the adatper');
		$frame->set_border_width(5);
			my $cb_min_adapter = Gtk2::ComboBox->new_text;
			foreach (8..16){
				$cb_min_adapter->append_text($_);
			}
			$cb_min_adapter->set_active(6);
		$frame->add($cb_min_adapter);
	$table->attach_defaults($frame, 0,2,1,2);
			#frame: min read seq
		$frame=Gtk2::Frame->new('Minimum length left in trimming');
		$frame->set_border_width(5);
			my $cb_min_trim_len = Gtk2::ComboBox->new_text;
			foreach (8..20){
				$cb_min_trim_len->append_text($_);
			}
			$cb_min_trim_len->set_active(2);
		$frame->add($cb_min_trim_len);
	$table->attach_defaults($frame, 2,4,1,2);

		#frame: check button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->set_border_width(5);
			$button->signal_connect(clicked=>sub{
				$variables{adapter_P5}=$entry_adapter_P5->get_text();
				$variables{adapter_len}=$cb_min_adapter->get_active_text;
				$variables{min_trim_len}=$cb_min_trim_len->get_active_text;
				#refresh log file
				sub_basic::refresh_log($variables{var_file}, 'adapter_P5', $variables{adapter_P5});
				sub_basic::refresh_log($variables{var_file}, 'adapter_len', $variables{adapter_len});
				sub_basic::refresh_log($variables{var_file}, 'min_trim_len', $variables{min_trim_len});
				$window->destroy;
			});
		$frame->add($button);
	$table->attach_defaults($frame,0,4,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

##################################
sub nb_4C_viewpoints_locking{
	#standard window creation, placement, and signal connecting
	my $title='Primary enzyme sites ( '.$variables{enzyme_site}.' ) of reference genome sequences';
	my $window = sub_gui::default_window($title, 1100, 600);
	my $table = Gtk2::Table->new(7,4,TRUE);
	
	#get sample names and fastq files
	my $sample_names_pointer=sub_basic::auto_sample_names($variables{raw_data_dir}, $variables{raw_file_format} );
	my @sample_names=@$sample_names_pointer;
	my $files_pointer=sub_basic::files_list($variables{raw_data_dir}, 'incrusive_file');
	my @files=@$files_pointer;
	my @fastq_files=grep(/\.fastq$|\.fq$/, @files);
	my %site_seq;
	
	#the list of sample names and fastq names
	my $frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
		my $scroll = Gtk2::ScrolledWindow->new;
		$scroll->set_policy ('never', 'automatic');
		## create and load the model
		#my @site_info_items=('Site name', 'Chromosome', 'Site position', 'Start territory', 'End territory', 'Upstream length', 'downstream length');
		#my $liststore = Gtk2::ListStore->new ('Glib::String', 'Glib::String', 'Glib::Int', 'Glib::Int','Glib::Int', 'Glib::Int', 'Glib::Int');
		#my $treeview = sub_gui::column_treeview($liststore, \@site_info_items);   #subroutine
		#$scroll->add ($treeview);
			my $slist = Gtk2::SimpleList->new ('Target site'=> 'bool', 'Site name'=> 'text',  'Site position'=>'int', 
								'Upstream sequence of viewpoint'=>'text', 'Downstream sequence of viewpoint'=>'text',		);
			$slist->get_selection->set_mode ('multiple');
			# simple way to make text columns editable
			$slist->set_column_editable (1, TRUE);
			$slist->set_rules_hint (TRUE);
			$slist->set_reorderable (TRUE);
			map { $_->set_resizable (TRUE) } $slist->get_columns;
		$scroll->add ($slist);
		$scroll->show_all;
	$frame->add($scroll);
	$table->attach_defaults($frame,0,4,0,6);
		
		#frame: enzyme sites file selection
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Genome sequences in FASTA');
			my $file_chooser_genome =Gtk2::FileChooserButton->new ('select a file' , 'open');
			$file_chooser_genome->set_filename($variables{ref_dir});
			$file_chooser_genome->hide;
		$frame->add($file_chooser_genome);
	$table->attach_defaults($frame,0,1,6,7);
		#frame: select 4C library
		$frame=Gtk2::Frame->new('4C library');
		$frame->set_border_width(5);
			my $cb_4C_library = Gtk2::ComboBox->new_text;
			foreach (@sample_names){
				$cb_4C_library->append_text($_);
			}
			$cb_4C_library->set_active(0);
		$frame->add($cb_4C_library);
	$table->attach_defaults($frame, 1,2,6,7);
		#frame
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Upstream/downstream');
			my $cb_up_down=Gtk2::ComboBox->new_text;
			$cb_up_down->append_text('Up');
			$cb_up_down->append_text('Down'); 
			$cb_up_down->set_active(0);
		$frame->add($cb_up_down);
	$table->attach_defaults($frame, 2,3,6,7);
		#frame:
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Chromosome');
			my $entry_chr = Gtk2::Entry->new();
			$entry_chr->set_text('11');
			$frame->add($entry_chr);
	$table->attach_defaults($frame, 0,1,7,8);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Start of genomic region');
			my $entry_start_region = Gtk2::Entry->new();
			$entry_start_region->set_text(32181947);
		$frame->add($entry_start_region);
	$table->attach_defaults($frame, 1,2,7,8);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('End of genomic region');
			my $entry_end_region = Gtk2::Entry->new();
			$entry_end_region->set_text(32185207);
		$frame->add($entry_end_region);
	$table->attach_defaults($frame, 2,3,7,8);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
		$frame->set_label('Length of viewpoints');
			my $entry_view_length = Gtk2::Entry->new();
			$entry_view_length->set_text(100);
		$frame->add($entry_view_length);
	$table->attach_defaults($frame, 3,4,7,8);
	
		#button: detect enzyme sites along chromosome
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_1=Gtk2::Button->new('Analyze chromosome');
			$button_1->set_border_width(5);
			$button_1->signal_connect(clicked=>sub{
				undef %site_seq;
				$variables{genome_fasta_file}=$file_chooser_genome->get_filename;
				my $sample_name=$cb_4C_library->get_active_text();
				$variables{'4C_up_down'}=$cb_up_down->get_active_text();
				my $site_chr=$entry_chr->get_text();
				my $start_region=$entry_start_region->get_text();
				my $end_region=$entry_end_region->get_text();
				my $view_len=$entry_view_length->get_text();
				$view_len=100 if $view_len<10;
						
				@{$slist->{data}}=(); #clear slist
				my $judging=0;
				$judging=2 unless -f $variables{genome_fasta_file};
				$judging=3 if $start_region>=$end_region;
				if ($judging==0){ #judging circyle
					my $in_obj = Bio::SeqIO->new(-file => $variables{genome_fasta_file}, -format => 'fasta');
					while (my $seq_obj = $in_obj->next_seq() ) {
						my $chr=$seq_obj->display_id;
						if($chr eq $site_chr){
							my $chr_seq=$seq_obj->seq();
							my $chr_len=$seq_obj->length;
							my $territory_seq=$seq_obj->subseq($start_region, $end_region);
							my $enzyme_sites_str=sub_3C::enzyme_sites_counting($territory_seq, $variables{enzyme_site}, $start_region);#subrountine
							my @enzyme_sites=split(",", $enzyme_sites_str);
							my $enzyme_num=@enzyme_sites;
							for(my $i=0; $i<@enzyme_sites; $i++){
								my $site_name=$chr.'_E'.($i+1);
								my $site_position=$enzyme_sites[$i];
								my $upstream_seq=substr($chr_seq, $site_position-1-20+length($variables{enzyme_site}), 20);
								my $downstream_seq=sub_basic::reverse_complement(substr($chr_seq, $site_position-1, 20));
								my @enzyme_info=(TRUE, $site_name, $site_position, $upstream_seq, $downstream_seq);
								push @{$slist->{data}}, \@enzyme_info;
								my $region_seq=substr($chr_seq, $site_position-1-$view_len, $view_len*2);
								if ($variables{'4C_up_down'} eq 'Up'){
									my $up_seq=substr($chr_seq, $site_position-1-$view_len+length($variables{enzyme_site}), $view_len);
									$site_seq{$site_name}->{view_seq}= $up_seq;
									$site_seq{$site_name}->{region_seq}=$region_seq;
								}
								else{
									my $down_seq=sub_basic::reverse_complement(substr($chr_seq, $site_position-1, $view_len));
									$site_seq{$site_name}->{view_seq}= $down_seq;
									$site_seq{$site_name}->{region_seq}=sub_basic::reverse_complement($region_seq);
								}
							}
							sub_gui::popup_window('Notice', 'Enzyme sites are listed on the left scrolled window!');
							$judging=1;
							$judging=4 if $start_region>=$chr_len or $end_region>=$chr_len or $start_region>=$end_region;
							$judging=5 if @enzyme_sites==0;
							last;
						}
					}
				}#judging circyle
				
				sub_gui::popup_window('Error!', 'No a selected chromosome or wrong chromosome input!') if $judging ==0;
				sub_gui::popup_window('Error!', 'Please select genome fasta file!') if $judging==2;
				sub_gui::popup_window('Error!', "Wrong! Start territory should be beyond end territory!") if $judging ==3;
				sub_gui::popup_window('Error!', "Wrong territory input of chromosome $site_chr!") if $judging ==4;
				sub_gui::popup_window('Error!', "No enzyme ($variables{enzyme_site}) are deteced!") if $judging ==5;
			});
		$frame->add($button_1);
	$table->attach_defaults($frame,0,1,8,9);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_3=Gtk2::Button->new("Append enzyme sites");
			$button_3->set_border_width(5);
			$button_3->signal_connect(clicked=>sub{
						$variables{genome_fasta_file}=$file_chooser_genome->get_filename;
						my $sample_name=$cb_4C_library->get_active_text();
						$variables{'4C_up_down'}=$cb_up_down->get_active_text();
						my $site_chr=$entry_chr->get_text();
						my $view_len=$entry_view_length->get_text();
						$view_len=100 if $view_len<10;
						
						#export
						open my($OUT), ">>", $variables{site_info_file} or die;
						my $n=0;
						foreach my $pointer(@{$slist->{data}}){
							my @items=@$pointer;
							if ($items[0]==1){
								#my $site_name=$items[2];
								my $site_position=$items[2];
								print $OUT join(',', $sample_name, $site_chr, $site_position, $variables{'4C_up_down'}, $view_len), "\n";
								$n++;
							}
						}
						close($OUT);
						sub_gui::popup_window('Notice', "A total of $n enzyme sites were appended into $variables{site_info_file}!");
					});
		$frame->add($button_3);
	$table->attach_defaults($frame,1,2,8,9);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_2=Gtk2::Button->new("Clear site_info.csv");
			$button_2->set_border_width(5);
			$button_2->signal_connect(clicked=>sub{
						unlink($variables{site_info_file}) if -f $variables{site_info_file};
						sub_gui::popup_window('Notice', "$variables{site_info_file} was reset!");
				});
		$frame->add($button_2);
	$table->attach_defaults($frame,2,3,8,9);
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button_4=Gtk2::Button->new('Close');
			$button_4->set_border_width(5);
			$button_4->signal_connect(pressed=>sub{	$window->destroy;		});
		$frame->add($button_4);
	$table->attach_defaults($frame,3,4,8,9);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}
#######################
sub nb_4C_colocalization_detection{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of genome mapping (Bowtie 2)', 900, 450);
	my $table=Gtk2::Table->new(5,6,TRUE);
	
		#genome_fasta
		my $file_names_pointer=sub_basic::files_list($variables{ref_dir}, 'file_name');
		my @file_names=@$file_names_pointer;
		my @fa_names=grep(/\.fa$|\.fasta$/, @file_names);
		#frame: select genome sequence (*.fa)
		$frame=Gtk2::Frame->new('Genome sequences (.fa file)');
		$frame->set_border_width(5);
			my $cb_fasta_file = Gtk2::ComboBox->new_text;
			foreach (@fa_names){
				$cb_fasta_file->append_text($_);
			}
			$cb_fasta_file->set_active(0);
		$frame->add($cb_fasta_file);
	$table->attach_defaults($frame, 0,3,0,1);
	
		#Frame: alignment mode --end-to-end/--local
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Alignment mode");
			my $cb_alignment_mode = Gtk2::ComboBox->new_text;
			$cb_alignment_mode->append_text('--end-to-end --very-fast');
			$cb_alignment_mode->append_text('--end-to-end --fast');
			$cb_alignment_mode->append_text('--end-to-end --sensitive');
			$cb_alignment_mode->append_text('--end-to-end --very-sensitive');
			$cb_alignment_mode->set_active(3);
		$frame->add($cb_alignment_mode);
	$table->attach_defaults($frame, 3,6,0,1);
		#Frame: mapping quality 
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Cutoff for mapping quality");
			my $cb_mapping_quality = Gtk2::ComboBox->new_text;
			foreach(1,10,13,20){
				$cb_mapping_quality->append_text($_);
			}
			$cb_mapping_quality->set_active(1);
		$frame->add($cb_mapping_quality);
	$table->attach_defaults($frame, 0,3,1,2);
	
		#Frame: --np
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Score of ambiguous penality (--np)");
			my $cb_ambiguous_penalty= Gtk2::ComboBox->new_text;
			foreach (1..4){
				$cb_ambiguous_penalty->append_text($_);
			}
			$cb_ambiguous_penalty->set_active(0);
		$frame->add($cb_ambiguous_penalty);
		$table->attach_defaults($frame, 3,6,1,2);
		#Frame: --rdg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the read gap penalities (--rdg)");
			my $cb_read_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_read_gap_penalty->append_text('3,2');
			$cb_read_gap_penalty->append_text('5,3');
			$cb_read_gap_penalty->append_text('8,4');
			$cb_read_gap_penalty->set_active(1);
		$frame->add($cb_read_gap_penalty);
		$table->attach_defaults($frame, 0,3,2,3);
		#Frame: --rfg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the reference gap penalities (--rfg)");
			my $cb_ref_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_ref_gap_penalty->append_text('3,2');
			$cb_ref_gap_penalty->append_text('5,3');
			$cb_ref_gap_penalty->append_text('8,4');
			$cb_ref_gap_penalty->set_active(1);
		$frame->add($cb_ref_gap_penalty);
		$table->attach_defaults($frame, 3,6,2,3);
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
		$table->attach_defaults($frame, 0,3,3,4);
		#Frame:--mp 
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of maximum and minimum mismatch penalties (--mp)");
			my $cb_mismatch_penality = Gtk2::ComboBox->new_text;
			$cb_mismatch_penality->append_text('4,2');
			$cb_mismatch_penality->append_text('6,2');
			$cb_mismatch_penality->append_text('6,3');
			$cb_mismatch_penality->set_active(1);
		$frame->add($cb_mismatch_penality);
	$table->attach_defaults($frame, 3,6,3,4);
	
	#Frame
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	my $button=Gtk2::Button->new('Save and close');
	$button->signal_connect(clicked=>sub{
		my $genome_fasta_file=$cb_fasta_file->get_active_text;
		$variables{genome_fasta_file}=$variables{ref_dir}.'/'.$genome_fasta_file;
		$variables{genome_index_name}=$genome_fasta_file;
		$variables{genome_index_name}=~s/\.fa$|\.fasta$//;
		$variables{genome_index}=$variables{alignment_dir}.'/'.$variables{genome_index_name};
		$variables{mapping_quality}=$cb_mapping_quality->get_active_text();
		#refresh log file
		sub_basic::refresh_log($variables{var_file}, 'genome_fasta_file', $variables{genome_fasta_file}) ;
		sub_basic::refresh_log($variables{var_file}, 'genome_index', $variables{genome_index}) ;
		sub_basic::refresh_log($variables{var_file}, 'genome_index_name', $variables{genome_index_name});
		sub_basic::refresh_log($variables{var_file}, 'mapping_quality', $variables{mapping_quality}) ;
		#bowtie options
		if (exists $variables{genome_fasta_file} and -f $variables{genome_fasta_file}){
			$variables{bowtie_options}=join(" ", 
							$variables{alignment_dir}.'/bowtie2',  #bowtie executable script
							$cb_alignment_mode->get_active_text(), 
							'--mp', $cb_mismatch_penality->get_active_text(),
							'--np', $cb_ambiguous_penalty->get_active_text(),
							'--rdg', $cb_read_gap_penalty->get_active_text(),
							'--rfg', $cb_ref_gap_penalty->get_active_text(),
							'--no-hd', '-x', $variables{genome_index}, #bowtie index
						);
			$variables{'3C_colocalization'}='yes';
			sub_basic::refresh_log($variables{var_file}, '3C_colocalization', $variables{'3C_colocalization'} );
			sub_basic::refresh_log($variables{var_file}, 'bowtie_options', $variables{bowtie_options} );
			$window->destroy;
		}
		else{
			sub_gui::popup_window('Error!', "Please select a fasta file of the reference genome in Setup!");
		}

	});
	$button->set_border_width(5);
	$frame->add($button);
	$table->attach_defaults($frame,1,5,4,5);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}

#######################
sub nb_4C_colocalization_counting{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of co-localization counting', 400, 200);
	my $table=Gtk2::Table->new(3,6,TRUE);
		
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Noise of read counts");
			my $cb_reads_noise = Gtk2::ComboBox->new_text;
			foreach (1,2,3,4,5,10,20){
				$cb_reads_noise->append_text($_);
			}
			$cb_reads_noise->set_active(0);
		$frame->add($cb_reads_noise);
	$table->attach_defaults($frame,0,6,0,1);
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Threshold of probability");
			my $cb_cumulative_p = Gtk2::ComboBox->new_text;
			foreach(0.05,0.01,0.001){
				$cb_cumulative_p->append_text($_);
			}
			$cb_cumulative_p->set_active(1);
		$frame->add($cb_cumulative_p);
	$table->attach_defaults($frame,0,6,1,2);

		#Frame
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->signal_connect(clicked=>sub{
				$variables{RC_noise}=$cb_reads_noise->get_active_text();
				$variables{RC_cumulative_p}=$cb_cumulative_p->get_active_text();
				$variables{'3C_statistics'}='yes';
				$variables{'3C_RC_counting'}='yes';
				$variables{'3C_Tscore'}='yes';
				
				sub_basic::refresh_log($variables{var_file}, 'RC_noise', $variables{RC_noise} );
				sub_basic::refresh_log($variables{var_file}, 'RC_cumulative_p', $variables{RC_cumulative_p} );
				sub_basic::refresh_log($variables{var_file}, '3C_statistics', $variables{'3C_statistics'});
				sub_basic::refresh_log($variables{var_file}, '3C_RC_counting', $variables{'3C_RC_counting'});
				sub_basic::refresh_log($variables{var_file}, '3C_Tscore', $variables{'3C_Tscore'});
				$window->destroy;
			});
		$frame->add($button);
	$table->attach_defaults($frame,1,5,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}

######################################################################
sub nb_HiC{
	my $scroller = Gtk2::ScrolledWindow->new();
	$scroller->set_shadow_type ('etched-out');
	$scroller->set_policy ('automatic', 'automatic');
	$scroller->set_size_request (200, 200);

	my $table_stock=Gtk2::Table->new(3, 8, TRUE);
	$table_stock->set_border_width(30);
	$table_stock->set_row_spacings(30);
	$table_stock->set_size_request (100, 100);
	#1
	#$table_stock->attach_defaults(&sub_gui::nb_button('FASTQ trimming', \&nb_HiC_fastq_trimming), 0,2,0,1);
	#$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 2,3,0,1);
	#2
	$table_stock->attach_defaults(&sub_gui::nb_button('Co-localization detection', \&nb_HiC_colocalization_detection), 0,3,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('right', 'etched-out'), 3,5,0,1); #add an arrow
	#3
	$table_stock->attach_defaults(&sub_gui::nb_button('Co-localization counting', \&nb_HiC_colocalization_counting), 5,8,0,1);
	$table_stock->attach_defaults(Gtk2::Arrow->new('down', 'etched-out'), 6,8,1,2); #add an arrow
	#run
	my $tbutton=Gtk2::Button->new('Apply and Run');
	$tbutton->signal_connect(clicked=>sub{
		my $thread_monitor=threads->new( sub{
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir.'/time_monitor.pl';
					system("perl $perl_script $variables{var_file}");
				}
				else{	sleep 1;	}
			}
		});

		my $thread_main=threads->new( sub{
			#local $SIG{KILL} = sub { threads->exit };
			while(1){
				return if $shash{'die'} == 1;
				if ( $shash{'go'} == 1 ){
					my $perl_script=$perl_dir."/main_HiC.pl";
					system("perl $perl_script $variables{var_file} > $variables{log_file} 2>&1");
					$shash{'go'} = 0; #turn off main running 
				}
				else{	sleep 1;	}
			}
		});

		if($shash{'go'} == 0){
			$shash{'go'} = 1;
			$pbar->show;
			Glib::Timeout->add (100, sub {
				if($shash{'go'} == 1){
					$pbar->set_text('Running!');
					$pbar->pulse;
					return TRUE;
				}
				else{	
					$pbar->set_text('OK! It is done!');
					return FALSE;
				}
			});
		}
	});  
	$table_stock->attach_defaults($tbutton, 5,8,2,3);
	$scroller->add_with_viewport($table_stock);
	return($scroller);
}

#############################
sub nb_HiC_fastq_trimming{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::window('Options for adapter trimming', 500, 250);
	my $table = Gtk2::Table->new(3,4,TRUE);
		
		#frame: 3' adapter sequences
		$frame=Gtk2::Frame->new("P5 Illumina adapter sequences");
		$frame->set_border_width(5);
			my $entry_adapter_P5 = Gtk2::Entry->new();
			$entry_adapter_P5->set_text($variables{adapter_P5});
		$frame->add($entry_adapter_P5);
	$table->attach_defaults($frame, 0,4,0,1);
		#frame: adapter length
		$frame=Gtk2::Frame->new('Minimum length of the adatper');
		$frame->set_border_width(5);
			my $cb_min_adapter = Gtk2::ComboBox->new_text;
			foreach (8..16){
				$cb_min_adapter->append_text($_);
			}
			$cb_min_adapter->set_active(6);
		$frame->add($cb_min_adapter);
	$table->attach_defaults($frame, 0,2,1,2);
			#frame: min read seq
		$frame=Gtk2::Frame->new('Minimum length left in trimming');
		$frame->set_border_width(5);
			my $cb_min_trim_len = Gtk2::ComboBox->new_text;
			foreach (8..20){
				$cb_min_trim_len->append_text($_);
			}
			$cb_min_trim_len->set_active(2);
		$frame->add($cb_min_trim_len);
	$table->attach_defaults($frame, 2,4,1,2);

		#frame: check button
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->set_border_width(5);
			$button->signal_connect(clicked=>sub{
				$variables{adapter_P5}=$entry_adapter_P5->get_text();
				$variables{adapter_len}=$cb_min_adapter->get_active_text;
				$variables{min_trim_len}=$cb_min_trim_len->get_active_text;
				#refresh log file
				sub_basic::refresh_log($variables{var_file}, 'adapter_P5', $variables{adapter_P5});
				sub_basic::refresh_log($variables{var_file}, 'adapter_len', $variables{adapter_len});
				sub_basic::refresh_log($variables{var_file}, 'min_trim_len', $variables{min_trim_len});
				sub_basic::refresh_log($variables{var_file}, '3C_trimming', 'yes');
				$window->destroy;
			});
		$frame->add($button);
	$table->attach_defaults($frame,0,4,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show();
	Gtk2->main();
}

#######################
sub nb_HiC_colocalization_detection{
	#standard window creation, placement, and signal connecting
	my$window=sub_gui::default_window('Options of genome mapping (Bowtie 2)', 1000, 500);
	my $table=Gtk2::Table->new(6,6,TRUE);
	
		#genome_fasta
		my $file_names_pointer=sub_basic::files_list($variables{ref_dir}, 'file_name');
		my @file_names=@$file_names_pointer;
		my @fa_names=grep(/\.fa$|\.fasta$/, @file_names);
		#frame: select genome sequence (*.fa)
		$frame=Gtk2::Frame->new('Genome sequences (.fa file)');
		$frame->set_border_width(5);
			my $cb_fasta_file = Gtk2::ComboBox->new_text;
			foreach (@fa_names){
				$cb_fasta_file->append_text($_);
			}
			$cb_fasta_file->set_active(0);
		$frame->add($cb_fasta_file);
	$table->attach_defaults($frame, 0,3,0,1);
	
		#Frame: alignment mode --end-to-end/--local
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Alignment mode");
			my $cb_alignment_mode = Gtk2::ComboBox->new_text;
			$cb_alignment_mode->append_text('--end-to-end --very-fast');
			$cb_alignment_mode->append_text('--end-to-end --fast');
			$cb_alignment_mode->append_text('--end-to-end --sensitive');
			$cb_alignment_mode->append_text('--end-to-end --very-sensitive');
			$cb_alignment_mode->set_active(3);
		$frame->add($cb_alignment_mode);
	$table->attach_defaults($frame, 3,6,0,1);
		#Frame: mapping quality 
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Cutoff for mapping quality");
			my $cb_mapping_quality = Gtk2::ComboBox->new_text;
			foreach(1,10,13,20){
				$cb_mapping_quality->append_text($_);
			}
			$cb_mapping_quality->set_active(1);
		$frame->add($cb_mapping_quality);
	$table->attach_defaults($frame, 0,3,1,2);

		#Frame: --np
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Score of ambiguous penality (--np)");
			my $cb_ambiguous_penalty= Gtk2::ComboBox->new_text;
			foreach (1..4){
				$cb_ambiguous_penalty->append_text($_);
			}
			$cb_ambiguous_penalty->set_active(0);
		$frame->add($cb_ambiguous_penalty);
		$table->attach_defaults($frame, 3,6,1,2);
		#Frame: --rdg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the read gap penalities (--rdg)");
			my $cb_read_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_read_gap_penalty->append_text('3,2');
			$cb_read_gap_penalty->append_text('5,3');
			$cb_read_gap_penalty->append_text('8,4');
			$cb_read_gap_penalty->set_active(1);
		$frame->add($cb_read_gap_penalty);
		$table->attach_defaults($frame, 0,3,2,3);
		#Frame: --rfg
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of the reference gap penalities (--rfg)");
			my $cb_ref_gap_penalty = Gtk2::ComboBox->new_text;
			$cb_ref_gap_penalty->append_text('3,2');
			$cb_ref_gap_penalty->append_text('5,3');
			$cb_ref_gap_penalty->append_text('8,4');
			$cb_ref_gap_penalty->set_active(1);
		$frame->add($cb_ref_gap_penalty);
		$table->attach_defaults($frame, 3,6,2,3);
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
		$table->attach_defaults($frame, 0,3,3,4);
		#Frame: -I
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Minimum fragment length for valid paired-end alignments (-I)");
			my $cb_min_insert = Gtk2::ComboBox->new_text;
			foreach (0,20,40,60,80,100){
				$cb_min_insert->append_text($_);
			}
			$cb_min_insert->set_active(0);
		$frame->add($cb_min_insert);
	$table->attach_defaults($frame,3,6,3,4);
		#Frame: -X
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Maximum fragment length for valid paired-end alignments (-X)");
			my $cb_max_insert = Gtk2::ComboBox->new_text;
			foreach (100, 200,300,400,500,600){
				$cb_max_insert->append_text($_);
			}
			$cb_max_insert->set_active(0);
		$frame->add($cb_max_insert);
	$table->attach_defaults($frame,0,3,4,5);
		#Frame:--mp 
		$frame = Gtk2::Frame->new();
		$frame->set_label ("Scores of maximum and minimum mismatch penalties (--mp)");
			my $cb_mismatch_penality = Gtk2::ComboBox->new_text;
			$cb_mismatch_penality->append_text('4,2');
			$cb_mismatch_penality->append_text('6,2');
			$cb_mismatch_penality->append_text('6,3');
			$cb_mismatch_penality->set_active(1);
		$frame->add($cb_mismatch_penality);
	$table->attach_defaults($frame, 3,6,4,5);
	
	#Frame
	$frame = Gtk2::Frame->new();
	$frame->set_border_width(5);
	my $button=Gtk2::Button->new('Save and close');
	$button->signal_connect(clicked=>sub{
		my $genome_fasta_file=$cb_fasta_file->get_active_text;
		$variables{genome_fasta_file}=$variables{ref_dir}.'/'.$genome_fasta_file;
		$variables{genome_index_name}=$genome_fasta_file;
		$variables{genome_index_name}=~s/\.fa$|\.fasta$//;
		$variables{genome_index}=$variables{alignment_dir}.'/'.$variables{genome_index_name};
		$variables{mapping_quality}=$cb_mapping_quality->get_active_text();
		#bowtie options
		if (exists $variables{genome_fasta_file} and -f $variables{genome_fasta_file}){
			$variables{bowtie_options}=join(" ", 
					$variables{alignment_dir}.'/bowtie2',  #bowtie executable script
					$cb_alignment_mode->get_active_text(), 
					'--mp', $cb_mismatch_penality->get_active_text(), '--np', $cb_ambiguous_penalty->get_active_text(),
					'--rdg', $cb_read_gap_penalty->get_active_text(), '--rfg', $cb_ref_gap_penalty->get_active_text(),
					'-I', $cb_min_insert->get_active_text(), '-X', $cb_max_insert->get_active_text(),
					'--no-hd', '-x', $variables{genome_index}, #bowtie index
				);
			$variables{'3C_colocalization'}='yes';
			$variables{'3C_one_detection'}='yes';
			$variables{'3C_no_detection'}='yes';
			#refresh log file
			sub_basic::refresh_log($variables{var_file}, 'genome_fasta_file', $variables{genome_fasta_file}) ;
			sub_basic::refresh_log($variables{var_file}, 'genome_index', $variables{genome_index}) ;
			sub_basic::refresh_log($variables{var_file}, 'genome_index_name', $variables{genome_index_name});
			sub_basic::refresh_log($variables{var_file}, 'mapping_quality', $variables{mapping_quality}) ;
			sub_basic::refresh_log($variables{var_file}, '3C_colocalization', $variables{'3C_colocalization'} );
			sub_basic::refresh_log($variables{var_file}, '3C_one_detection', $variables{'3C_one_detection'} );
			sub_basic::refresh_log($variables{var_file}, '3C_no_detection', $variables{'3C_no_detection'} );
			sub_basic::refresh_log($variables{var_file}, 'bowtie_options', $variables{bowtie_options} );
			$window->destroy;
		}
		else{
			sub_gui::popup_window('Error!', "Please select a fasta file of the reference genome in Setup!");
		}

	});
	$button->set_border_width(5);
	$frame->add($button);
	$table->attach_defaults($frame,1,5,5,6);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}

#######################
sub nb_HiC_colocalization_counting{
	#standard window creation, placement, and signal connecting
	my $window = sub_gui::default_window('Options of co-localization counting', 400, 200);
	my $table=Gtk2::Table->new(3,6,TRUE);
		
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Noise of read counts");
			my $cb_reads_noise = Gtk2::ComboBox->new_text;
			foreach (1,2,3,4,5,10,20){
				$cb_reads_noise->append_text($_);
			}
			$cb_reads_noise->set_active(1);
		$frame->add($cb_reads_noise);
	$table->attach_defaults($frame,0,6,0,1);
		#Frame: 
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(2);
		$frame->set_label ("Threshold of probability");
			my $cb_cumulative_p = Gtk2::ComboBox->new_text;
			foreach(0.05,0.01,0.001){
				$cb_cumulative_p->append_text($_);
			}
			$cb_cumulative_p->set_active(1);
		$frame->add($cb_cumulative_p);
	$table->attach_defaults($frame,0,6,1,2);

		#Frame
		$frame = Gtk2::Frame->new();
		$frame->set_border_width(5);
			my $button=Gtk2::Button->new('Save and close');
			$button->signal_connect(clicked=>sub{
				$variables{'RC_noise'}=$cb_reads_noise->get_active_text();
				$variables{'RC_cumulative_p'}=$cb_cumulative_p->get_active_text();
				$variables{'3C_statistics'}='yes';
				$variables{'3C_RC_counting'}='yes';
				$variables{'3C_Tscore'}='yes';
				
				sub_basic::refresh_log($variables{var_file}, 'RC_noise', $variables{RC_noise} );
				sub_basic::refresh_log($variables{var_file}, 'RC_cumulative_p', $variables{RC_cumulative_p} );
				sub_basic::refresh_log($variables{var_file}, '3C_statistics', $variables{'3C_statistics'});
				sub_basic::refresh_log($variables{var_file}, '3C_RC_counting', $variables{'3C_RC_counting'});
				sub_basic::refresh_log($variables{var_file}, '3C_Tscore', $variables{'3C_Tscore'});
				$window->destroy;
			});
		$frame->add($button);
	$table->attach_defaults($frame,1,5,2,3);
	
	$table->show_all();
	$window->add($table);
	$window->show_all;
	Gtk2->main();
}



####################################################
#used for items' selection
sub items_selection{
	my ($index_names_pointer, $variable_file, $variable_name)=@_;
	my @index_names=@$index_names_pointer;
	my @selected_names;
	
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
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
	} );
	$table->attach_defaults($button_1, 3,4,1,2);
	my $button_2=Gtk2::Button->new('>>');        ###### >> button: all transfered
	$button_2->signal_connect(clicked=>sub{
		$tree_store_1->clear;
		foreach (@index_names) {
			my $iter = $tree_store_2->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_2->set ($iter,0 => $_);
		}
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @index_names) );
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
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
	} );
	$table->attach_defaults($button_3, 3,4,3,4);
	my $button_4=Gtk2::Button->new('<<');      ######<< button
	$button_4->signal_connect(clicked=>sub{
		$tree_store_2->clear;
		foreach (@index_names) {
			my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
			$tree_store_1->set ($iter,0 => $_);
		}
		sub_basic::refresh_log($variable_file, $variable_name, 'NA' );
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
		for(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i>0){
				$selected_names[$i-1]=$value;
				$selected_names[$i]=$pre_value;
			}
		}
		@selected_names=grep($_, @selected_names);
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
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
		for(my $i=0; $i<@selected_names; $i++){
			if ( ($selected_names[$i] eq $value) and $i<@selected_names-1){
				$selected_names[$i+1]=$value;
				$selected_names[$i]=$next_value;
			}
		}
		@selected_names=grep($_, @selected_names);
		sub_basic::refresh_log($variable_file, $variable_name, join(",", @selected_names) );
	} );
	$table->attach_defaults($button_6, 7,8,4,5);

	return ($table);
}

#####################
#used for items selection
sub multiple_items_selection{
	my ($index_names_pointer, $variable_file, $variable_name_A, $variable_name_B)=@_;
	my @index_names=@$index_names_pointer;
	my (@selected_A, @selected_B);
	
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
		$tree_column_1->set_title ("Sort sample names");
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
	$table->attach_defaults($sw_1, 4,7,0,6);

	  #group A
	#create a scrolled window that will host the treeview 
	my $sw_2 = Gtk2::ScrolledWindow->new (undef, undef);
	$sw_2->set_shadow_type ('etched-out');
	$sw_2->set_policy ('automatic', 'automatic');
	$sw_2->set_border_width(5);
		my $tree_store_2 = Gtk2::TreeStore->new(qw/Glib::String/);
		my $tree_view_2 = Gtk2::TreeView->new($tree_store_2);
		my $tree_column_2 = Gtk2::TreeViewColumn->new();
		$tree_column_2->set_title ("Sort selected group A");
			my $renderer_2 = Gtk2::CellRendererText->new;
		$tree_column_2->pack_start ($renderer_2, FALSE);
		$tree_column_2->add_attribute($renderer_2, text => 0);
		$tree_column_2->set_sort_column_id(0);# Allow sorting on the column
		$tree_view_2->append_column ($tree_column_2);
		$tree_view_2->set_search_column(0);# make it searchable
		$tree_view_2->set_reorderable(TRUE);# Allow drag and drop reordering of rows
	$sw_2->add($tree_view_2);
	$table->attach_defaults($sw_2, 0,3,0,6);

	#Widget of table: buttons of group A
	my $button_3=Gtk2::Button->new('<');      ######<< button
	$button_3->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_1->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_2->append(undef);  #add the value into the second tree 
		$tree_store_2->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_A, sub_gui::tree_store_str($tree_store_2) );
	} );
	$table->attach_defaults($button_3, 3,4,1,2);
	my $button_1=Gtk2::Button->new('>');        ######>> button: one transfered
	$button_1->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_2->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store->remove($iter_1);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_A, sub_gui::tree_store_str($tree_store_2) );
	} );
	$table->attach_defaults($button_1, 3,4,3,4);
	my $button_2=Gtk2::Button->new('>>');        ###### >> button: all transfered
	$button_2->signal_connect(clicked=>sub{
		my @names=split(',', sub_gui::tree_store_str($tree_store_2) );
		foreach(@names){
			unless($_ eq 'NA'){
				my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
				$tree_store_1->set($iter,0=> $_);
			}
		}
		$tree_store_2->clear;
		sub_basic::refresh_log($variable_file, $variable_name_A, 'NA' );
	} );
	$table->attach_defaults($button_2, 3,4,4,5);

	
	  #group B
	#create a scrolled window that will host the treeview 
	my $sw_3 = Gtk2::ScrolledWindow->new (undef, undef);
	$sw_3->set_shadow_type ('etched-out');
	$sw_3->set_policy ('automatic', 'automatic');
	$sw_3->set_border_width(5);
		my $tree_store_3 = Gtk2::TreeStore->new(qw/Glib::String/);
		my $tree_view_3 = Gtk2::TreeView->new($tree_store_3);
		my $tree_column_3 = Gtk2::TreeViewColumn->new();
		$tree_column_3->set_title ("Sort selected group B");
			my $renderer_3 = Gtk2::CellRendererText->new;
		$tree_column_3->pack_start ($renderer_3, FALSE);
		$tree_column_3->add_attribute($renderer_3, text => 0);
		$tree_column_3->set_sort_column_id(0);# Allow sorting on the column
		$tree_view_3->append_column ($tree_column_3);
		$tree_view_3->set_search_column(0);# make it searchable
		$tree_view_3->set_reorderable(TRUE);# Allow drag and drop reordering of rows
	$sw_3->add($tree_view_3);
	$table->attach_defaults($sw_3, 8,11,0,6);
	
	#Widget of table: button
		my $button_6=Gtk2::Button->new('>');        ######>> button: one transfered
	$button_6->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_1->get_selection;
		my($model, $iter_1)=$tree_selection->get_selected;
		my $value=$model->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_1->remove($iter_1);
		my $iter=$tree_store_3->append(undef);  #add the value into the second tree 
		$tree_store_3->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_B, sub_gui::tree_store_str($tree_store_3) );
	} );
	$table->attach_defaults($button_6, 7,8,1,2);
	my $button_4=Gtk2::Button->new('<');      ######<< button
	$button_4->signal_connect(clicked=>sub{
		my $tree_selection=$tree_view_3->get_selection;
		my($tree_store, $iter_1)=$tree_selection->get_selected;
		my $value=$tree_store->get($iter_1, 0);
		return unless $iter_1;
		$tree_store_3->remove($iter_1);
		my $iter=$tree_store_1->append(undef);  #add the value into the second tree 
		$tree_store_1->set($iter,0=> $value);
		sub_basic::refresh_log($variable_file, $variable_name_B, tree_store($tree_store_3) );
	} );
	$table->attach_defaults($button_4, 7,8,3,4);

	my $button_5=Gtk2::Button->new('<<');        ###### >> button: all transfered
	$button_5->signal_connect(clicked=>sub{
		my @names=split(',', sub_gui::tree_store_str($tree_store_3));
		foreach (@names) {
			unless($_ eq 'NA'){
				my $iter = $tree_store_1->append(undef); #the iter is a pointer in the treestore. We use to add data.
				$tree_store_1->set ($iter,0 => $_);
			}
		}
		$tree_store_3->clear;
		sub_basic::refresh_log($variable_file, $variable_name_B, 'NA' );
	} );
	$table->attach_defaults($button_5, 7,8,4,5);


	$table->show_all();
	return ($table);
}




