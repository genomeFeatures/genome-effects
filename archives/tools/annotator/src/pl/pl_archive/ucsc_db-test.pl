#!/usr/bin/perl
#**************************************************************************************
# This program collects all the current organism/db/build listing
# from ucsc browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start&clade=xxx)
# where xxx is the organism_group name (mammal,vertabrate,...)
# 
#
# Output : three config files - commas-separated files
#   a. ucsc_organism_version.txt) with the following fields
#       1.organism_group,
#       2.organims_name, 
#       3.ucsc_db, 
#       4.organism_version
#   b. ucsc_geneannotations.txt with the following fields:
#       1.organism_group,
#       2.organims_name, 
#       3.ucsc_db, 
#       4.annotation tables list separated by ":"
#   c. ucsc_organisms_w_chromosome_data.txt file
#      contains list of organisms version that have chromosome assembly. Fields include
#      1.ucsc_db
#      2. path to the chromosome data 
#
# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : October 2013
#
## Usage : perl program_name
#*******************************************************************************
use File::stat;
use LWP;

my $browser = LWP::UserAgent->new;
   $browser->timeout(10);
   $browser->env_proxy;
   $browser->agent('Mozilla/5.0');
 
$host="http://genome.ucsc.edu/cgi-bin/hgTables?command=start"; $user="genome";
$dbhost="genome-mysql.cse.ucsc.edu";$ftp_path="ftp://hgdownload.cse.ucsc.edu/goldenPath";
my $connection="mysql -h$dbhost -u$user -A";$tm = localtime;

my($directory)=@ARGV; chomp($directory);$directory=~s/\/$//;
my $filename="organismgroup".".txt";
my $response = $browser->get($host, ":content_file" => $filename); open(IN,"$filename");
  if(IN){@filecontent=<IN> ;}close(IN);$more=1; my @organismgroups=(); $i=0;
   #loadOrganismGroup();
  if(@filecontent>0){
     while($more){ $line=shift(@filecontent);chomp($line);
         if($line=~/<SELECT NAME="clade"/){ $end=1;
            while($end){
                  $line=shift(@filecontent);chomp($line);
                  if(($line=~/<OPTION VALUE="(.+)">(.+)<\/OPTION>/)||
                     ($line=~/<OPTION SELECTED VALUE="(.+)">(.+)<\/OPTION>/)){
                      $organism_group=$1;$organism_group=~s/ //g;
                      push(@organismgroups,$organism_group);
                    }
                   elsif($line=~/<\/SELECT>/){$end=0;$more=0;}
              }}
  } }
 my %organismmap;
$version_file="ucsc_organism_version.txt";$genome_file="ucsc_organisms_w_chromosome_data.txt";
$annot_file="ucsc_geneannotations.txt";
open(O,">$version_file");open(G,">$genome_file");open(A,">$annot_file");
while(@organismgroups){
      $organismg=shift(@organismgroups);
      $newurl="$host"."&clade=$organismg";my $filename="organismList".".txt";
      my $response = $browser->get($newurl, ":content_file" => $filename);
      open(IN,"$filename");if(IN){@filecontent=<IN>;}close(IN);
      $data="";$more=1;
      while($more){ $line=shift(@filecontent);chomp($line);
         if($line=~/<SELECT NAME="org"/i){
             while($more){$line=shift(@filecontent);chomp($line);
                   if($line=~/<\/SELECT>/i){$more=0;}else{$data.=$line;}
      }}}#now generate organism list
     @orglist=split(/<OPTION /i,$data);
     while(@orglist>0){ $line=shift(@orglist); 
          if($line=~/VALUE="(.+)">(.+)<\/OPTION>/){
             #now get the versions
             $organism=$1;$newurl="$host"."&clade=$organismg&org=$organism";
             my $response = $browser->get($newurl, ":content_file" => $filename);
             open(IN4,"$filename");$more=1;@filecontent=<IN4>;close(IN4);$data="";
              while($more){ $line2=shift(@filecontent);chomp($line2);
                   if($line2=~/<SELECT NAME="db"/){
                      while($more){$line2=shift(@filecontent);chomp($line2);
                           if($line2=~/<\/SELECT>/i){$more=0;}else{$data.=$line2;}}}
              }
             @dbs=split(/<OPTION /i,$data);print "$organismg,$organism\n";
             while(@dbs>0){$db=shift(@dbs);
                  if($db=~/VALUE="(.+)">(.+)<\/OPTION>/i){$db=$1;$desc=$2;
                     $data="";$more=1;
                     $commands=("$connection -e 'use $db; show tables' > $filename");
                     system($commands);open(IN4,"$filename");@filecontent=();
                     while(<IN4>){chomp($_);push(@filecontent,$_);}close(IN4);
                     @tablelist=();@tablelist=grep(/Gene\s*$/,@filecontent);
                     $newurl="$ftp_path/$db";
                     my $response = $browser->get($newurl, ":content_file" => $filename);
                     open(IN4,"$filename");@filecontent=();@filecontent=<IN4>;
                     @genome=grep(/chromosomes/,@filecontent);
                     print G "$organismg,$organism,$db,$ftp_path/$db/chromosomes\n" if(@genome>0);
                     print O "$organismg,$organism,$db,$desc\n";
                     if(@tablelist>0){
                        print A "$organismg,$organism,$db,".join(":",@tablelist)."\n"; 
                     }
                     
             }}
     }}
}


