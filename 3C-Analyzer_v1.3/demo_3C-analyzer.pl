#! /usr/bin/perl -w
use strict;
use warnings;
use threads;
use Bio::SeqIO;
use Bio::Seq::Quality;
use List::MoreUtils;
use List::Util;
use File::Find;

##############################
sub prepare_alignment{
	my($download_dir, $alignment_dir)=@_;
	
	print "Download bowtie software into $download_dir\n";
	my $files_http='http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip';
	system("wget -cq -P $download_dir $files_http");

	#
	my $zip_file=$download_dir.'download';
	print "Uncompress $zip_file\n";
	system("unzip -oq $zip_file -d $download_dir");
	
	#
	my $tmp_bowtie_dir=$download_dir.'bowtie2-2.2.5/*';
	system("cp $tmp_bowtie_dir $alignment_dir");
	
}
##################
sub prepare_samtools{
	my($download_dir, $alignment_dir)=@_;
	
	print "Download samtools into $download_dir\n";
	my $files_http='http://sourceforge.net/projects/samtools/files/samtools/1.2/samtools-1.2.tar.bz2';
	system("wget -cq -P $download_dir $files_http");
	
	#
	my $zip_file=$download_dir.'samtools-1.2.tar.bz2';
	print "Uncompress $zip_file\n";
	system("tar -jxf $zip_file -C $download_dir");
}
################3

#get the directory of perl scripts
my $perl_dir=Cwd::getcwd();
my $download_dir=$perl_dir.'/download/';
mkdir($download_dir, 0755) unless -d $download_dir;
my $alignment_dir=$perl_dir.'/bowtie/';
mkdir($alignment_dir, 0755) unless -d $alignment_dir;

#bowtie 
#prepare_alignment($download_dir, $alignment_dir);

#samtools
#repare_samtools($download_dir, $alignment_dir);

print "\n\nDemo analysis of 3C-analyzer is done.!\n\n";
