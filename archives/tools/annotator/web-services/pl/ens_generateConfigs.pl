#!/usr/bin/perl
#**************************************************************************************
# This program collects all the current organism/db/build listing
# from Ensembl database browser (http://useast.ensembl.org/index.html)
# 
#
# Output :A commas-separated file CONFIG FILE
#   a. ens_organism_version.txt) with the following fields
#       1.ens_db,
#       2.organims_cname, 
#       4.organism_sname
# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : October 2013
#
## Usage : perl program_name  -- this program is ran on the server side via cron job
#  how to map ensembl to ucsc annotations?
#  
#*******************************************************************************
use File::stat;
use LWP;

my $browser = LWP::UserAgent->new;
   $browser->timeout(10);
   $browser->env_proxy;
   $browser->agent('Mozilla/5.0');

$host="http://useast.ensembl.org/index.html"; my($directory)=@ARGV; 
chomp($directory);$directory=~s/\/$//;
$version_file="ens_organism_version.txt"; @organisms=();$filename="org_list_temp.txt";
my $response = $browser->get($host, ":content_file" => $filename); open(IN,"$filename");
  if(IN){@filecontent=<IN> ;}close(IN);$more=1;
  if(@filecontent>0){
     while($more){ $line=shift(@filecontent);chomp($line);
         if($line=~/<h3>All genomes<\/h3>/){ $end=1;
            while($end){
                  $line=shift(@filecontent);chomp($line);
                  if(($line=~/<option value="(.+)\/Info\/Index">(.+)<\/option>/)){
                      $ens_db=lc($1)."_core_";$ens_db=~s/ //g;$comon_name=$2; $sc_name=$1;$sc_name=~s/_/ /;
                      $data="$ens_db,$comon_name,$sc_name";
                      push(@organisms,$data);
                    }
                   elsif($line=~/<\/select>/){$end=0;$more=0;}
              }}
  } }
if($directory){$version_file="$directory/ens_organism_version.txt";}
if(@organisms>20){
   open(O,">$version_file");my $filename="organismgroup".".txt";
   while(@organisms>0){ $line=shift(@organisms);print O "$line\n";}
   close(O);
}
print "Program complete\n";
exit(0);

