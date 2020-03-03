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



#get the directory of perl scripts
our $perl_dir=Cwd::getcwd();

#get subroutines
require $perl_dir."/functions_basic.pm"; #sub_basic::
require $perl_dir."/functions_3C.pm";  #sub_3C::

#
package sub_gui;

#######################################################
#used in subroutine of notebook
sub nb_button{
	my ($label, $window)=@_;
	my $button=Gtk2::Button->new($label);
	$button->signal_connect(clicked => $window);
	return($button);
}

#################
sub default_window{
    my($title, $width,  $height)=@_;
    
	my $window = Gtk2::Window->new('toplevel');
	$window->set_title($title);
	$window->signal_connect('delete_event' => sub { Gtk2->main_quit; });
	$window->set_border_width(5);
	$window->set_position('center');
	$window->set_size_request($width,$height);
	
	return($window);
}

################
sub progress_bar{
	my($interval)=@_;
	
	my $frame=Gtk2::Frame->new();
	$frame->set_label('Progress bar');
	$frame->set_shadow_type('out');
	$frame->set_border_width(5);
		my $pbar = Gtk2::ProgressBar->new();
		$pbar->set_pulse_step($interval);
		$pbar->hide; #needs to be called after show_all
	$frame->add($pbar);
	
	return($frame)
}

############################
sub popup_window{
	my ($title, $str)=@_;
	
	my $window = Gtk2::Window->new();
	$window->set_title($title);
	$window->signal_connect('delete_event' => sub { Gtk2->main_quit; });
	$window->set_border_width(5);
	$window->set_position('center');
	$window->set_size_request('400','100');
		my $vbox=Gtk2::VBox->new(0, 5); #here the first 0 indicates FALSE
			my $label=Gtk2::Label->new($str);
		$vbox->add($label);
			my $button=Gtk2::Button->new_with_label('Close');
			$button->signal_connect(clicked=>sub{		$window->destroy;	 });
		$vbox->add($button);
	$window->add($vbox);
	$window->show_all;
	Gtk2->main();
}

##########

################
sub tree_store_str{
	my ($tree_store)=@_;
	
	my @items;
	$tree_store->foreach(sub{
			my($model, $path, $iter)=@_;
			my $value = $model->get ($iter, 0);
			push(@items, $value);
			return(0);
		});
	my $str=join(',', @items);
	return($str);
}

##############################
sub column_treeview{
	my ($liststore, $site_info_items_pointer)=@_;
	my @site_info_items=@$site_info_items_pointer;
	
	my $treeview = Gtk2::TreeView->new ($liststore);
	$treeview->set_rules_hint (1); #TRUE
	$treeview->set_reorderable (1);#TRUE
	my $num=0;
	foreach (@site_info_items){
		# renderer column 0:site_name
		my $renderer = Gtk2::CellRendererText->new;
		$renderer->set (editable => 1);#TRUE
		my $column = Gtk2::TreeViewColumn->new_with_attributes ($_, $renderer, text => $num);
		$renderer->signal_connect (edited => sub {
						my ($cell, $text_path, $new_text, $model) = @_;
						my $path = Gtk2::TreePath->new_from_string ($text_path);
						my $iter = $model->get_iter ($path);
						$model->set ($iter, $num, $new_text);
				}, $liststore);
		$treeview->append_column ($column);
		$num++;
	}
	
	return($treeview);
}


##########################################################
1; 
# make sure the file returns true or require will not succeed!#

