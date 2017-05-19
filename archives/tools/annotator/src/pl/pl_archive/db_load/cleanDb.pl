#!/usr/bin/perl
########################################################################
## This script collects all the versions for a given organism
## and for organisms with more than two versions,the script will
#  archive then delete the oldest versions including the genome
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : January 2013
#
#  Input : a directory path where to store the archived data
#
## set default variables
$dbname="graber_transcriptdb";
$host="harlequin"; $user="lnh";$pass="lucie98";

getopts('ho:');
if(($opt_h)||(!$opt_o)){
    print <<HELP;

 This script collects all the versions for a given organism
 and for organisms with more than two versions,the script will
 archive then delete the oldest versions 

Usage:
   perl cleanDb.pl -o output_directory

Arguments:
   -h  displays this help message
   -o  Output directory
  

Examples:
cmd: ./cleanDb.pl -o /scratch/data/graber_transcriptdb/archives/

HELP
exit;
}
chdir($opt_o);
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);


