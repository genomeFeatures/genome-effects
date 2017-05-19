#!/usr/bin/perl

use DBI;
use Time::localtime;
use LWP 5.64;

#**************************************************************************************
# This program collects all the current organism/db/build listing
# from ucsc browser (http://genome.ucsc.edu/cgi-bin/hgTables?command=start&clade=xxx)
# where xxx is the organism_group name (mammal,vertabrate,...)
# 
#
# Output : three config files 
#   a. ucsc_organism_version_date.txt) with the following fields
#       1.organism_group,
#       2.organims_name, 
#       3.ucsc_db, 
#       4.organism_version
#   b. ucsc_geneannotations_date.txt with the following fields:
#       1.organism_group,
#       2.organims_name, 
#       3.ucsc_db, 
#       4.annotation tables list separated by ":"
#   c. Configuration file
#       contains list of config files to be used
#
# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : May 2010
#
## Usage : perl program_name
#*******************************************************************************
#
my $browser = LWP::UserAgent->new;
   $browser->timeout(10);
   $browser->env_proxy;
   $browser->agent('Mozilla/5.0');
 

q{
 ***************************************************************
 given a filename containing a list of ucsc tables for a given
 organism version, and a list of annotations we would like 
 to download from ucsc, this function generates a new list
 of annotations to download
 
 **************************************************************
};
sub getCurrentTablesList{
  open(IN,"$filename"); my @currentTables;
  @currentTables=<IN> if(IN); 
  close(IN); $newtablelist="";$i=0;$j=0;
  #now for each table in our selected annotation, gene
  while($i<@annot){
        $geneannot=$annot[$i];chomp($geneannot);@foo = grep(/$geneannot/, @currentTables);  
        if(@foo>0){
           while(@foo>0){
                 $token=shift(@foo);
                 if($token=~/^\s*$geneannot\s*$/){
                    $token=~s/^\s*//; $token=~s/\s*$//;
                   if($j==0){$newtablelist="$token"; ++$j;}
                   else{$newtablelist.=":$token";}
                 }
            }
        }
        ++$i;
   }
  return $newtablelist; # I know this is not necessaire
}

$config=shift||"/home/lnh/docs/tools/annotator/src/pl/txutils/db_svn/docs/Configuration";
chomp($config);
open(CONF,"$config");
if(!CONF){
   print"Usage: perl run_update_pipeline.pl config_file\n";
   exit(1);
}
my $main_path="/scratch/data/ucsc"; # path to gene annotations
my $load_base_path="/home/lnh/docs/tools/annotator/src/pl/txutils/docs";
$host="genome-mysql.cse.ucsc.edu"; $user="genome";$graber_transcriptdb;$gb_host;
$gb_user;$gb_pass;
$db_dev_host;
@filecontent=<CONF>; 
my $generated_config="$load_base_path/prog_generated_config.txt";
open(OUT3,">$generated_config") or die "Bad file $generated_config $!\n";
while(@filecontent>0){
   $line=shift(@filecontent);chomp($line);
    print OUT3 "$line\n";
   my($variable,$path)=split(",",$line);
   if($variable=~/ucsc_load_base/){$load_base_path=$path;}
   elsif($variable=~/ucsc_db_host/){$host=$path;}
   elsif($variable=~/ucsc_db_usr/){$user=$path;}
   elsif($variable=~/ucsc_org_base/){$main_path=$path;}
   elsif($variable=~/graber_db_dev/){$graber_transcriptdb=$path;}
   elsif($variable=~/graber_db_host/){$gb_host=$path;}
   elsif($variable=~/graber_db_usr/){$gb_user=$path;}
   elsif($variable=~/graber_db_pass/){$gb_pass=$path;}
   elsif($variable=~/graber_db_dev_host/){$gb_dev_host=$path;}
}
close(CONF);
my $connection="mysql -h$host -u$user -A";$tm = localtime;
my ($mday, $mon, $yday) = ($tm->mday, ($tm->mon)+1, ($tm->year)+1900);
# Set the path to configuration files
my $path="$load_base_path/docs";my $config_dir=$path;
my $annotationlist="$path/ucsc_annotations_list.txt";my $currentTablePath="$load_base_path/currenttables"; 
# Configuration file Current corganims_group/organisms/db/version list
$current_config_file="ucsc_organism_version"."_"."$mon-$mday-$yday".".txt"; 
my $configFile="$path/$current_config_file"; 

# Configuration file Current corganims_group/organisms/db/version/tables list
my $annotationsConfig="$path/ucsc_geneannotations"."_"."$mon-$mday-$yday".".txt"; 
my $processlog="$path/data_load_log.log"; open(LOG,">$processlog") or die "$processlog: $!\n";

if(LOG){
  print LOG "\n*************************************************************\n";
  print LOG "Starting getCurrentOrganisms.pl program : $mday/$mon/$yday\n";
  print LOG "\n*************************************************************\n";
   my(@filecontent,@annot);
  open(IN,"$annotationlist");
  if(IN){@annot=<IN>;}
  close(IN);# now connect to the browser and get the organism group list
  $url="http://genome.ucsc.edu/cgi-bin/hgTables?command=start";
  my $filename="organismgroup".".txt";
  my $response = $browser->get($url, ":content_file" => $filename);
  open(IN,"$filename");
  if(IN){@filecontent=<IN> ;}
  close(IN);$more=1; my @organismgroups; $i=0;
  if(@filecontent<=0){
     print LOG "Bad connection to $url\n";close(LOG); exit(1); 
   }
   #loadOrganismGroup();
  if(@filecontent>0){
     while($more){
         $line=shift(@filecontent);chomp($line);
         if($line=~/<SELECT NAME="clade"/){ $end=1;
            while($end){
                  $line=shift(@filecontent);chomp($line);
                  if(($line=~/<OPTION VALUE="(.+)">(.+)<\/OPTION>/)||
                     ($line=~/<OPTION SELECTED VALUE="(.+)">(.+)<\/OPTION>/)){
                      $organism_group=$1;$organism_group=~ s/ //g;
                      $organismgroups[$i]=$organism_group;++$i;
                    }
                   elsif($line=~/<\/SELECT>/){$end=0;$more=0;}
              }
          }
       }
   }
   if(!(-d $currentTablePath)){
       mkdir("$currentTablePath",0777)||print LOG "$!\n";
   }
   # now for each organism group, get the corresponding list of organisms and organism versions
   open(OUT,">$configFile");open(OUT2,">$annotationsConfig"); #generate the list of config files
   my %organismmap;
   while(@organismgroups){
         $organismg=shift(@organismgroups);
         print "Processing $organismg\n";
         $newurl="$url"."&clade=$organismg";
         my $filename="organismgroup".".txt";
         my $response = $browser->get($newurl, ":content_file" => $filename);
        open(IN,"$filename");
        if(IN){@filecontent=<IN>;}close(IN);$more=1;
        if(@filecontent<=0){
           print LOG "$organismg group not loaded\n";$more=0;
        }
        while($more){
              $line=shift(@filecontent);chomp($line);
              if($line=~/<SELECT NAME="org"/i){
                 $end=1;
                 while($end){
                     $line=shift(@filecontent);chomp($line);
                     if(($line=~/<OPTION VALUE="(.+)">(.+)<\/OPTION>/)||
                        ($line=~/<OPTION SELECTED VALUE="(.+)">(.+)<\/OPTION>/)){
                        $organism=$1;$newurl="$url"."&clade=$organismg&org=$organism";
                        my $response = $browser->get($newurl, ":content_file" => $filename);
                        # die "$newurl error: ",$response->status_line unless $response->is_success;
                         open(IN4,"$filename"); #or die "Can not open $filename\n";
                       # @filecontent2=<IN4>; #close(IN);i
                        $more2=1; #get this organism current version
                        $ct=0; 
                        # print "Processing $filename\n";
                         while($more2){
                               $line2=<IN4>;chomp($line2);
                               #print "$line2\n";
                               if($line2=~/<SELECT NAME="db"/){
                                  $end2=1; $db=""; $version="";
                                  $line2=<IN4>;chomp($line2);
                                  $prev=<IN4>; chomp($prev); 
                                  #print"$line2\n";
                                  while($end2){
                                        #$line2=<IN>;chomp($line2);
                                        if(($line2=~/<OPTION SELECTED VALUE="(.+)">(.+)<\/OPTION>/i)||
                                           ($line2=~/<OPTION VALUE="(.+)">(.+)<\/OPTION>/i)){
                                            $db=$1; $version=$2;$organism=~s/\.\s+/\./;
                                           ++$ct; 
                                           print OUT "$organismg,$organism,$db,$version\n"; #get current tables list for this
                                           #print "Processinf $organismg,$organism,$db,$version\n";
                                           $thisor=$organism; $thisor=~s/\s+/-/g;
                                           $thisfilename="$currentTablePath/$thisor"."_$db"."_currentTables_"."$mon-$mday-$yday"; 
                                           $commands=("$connection -e 'use $db; show tables' > $thisfilename");
                                           system($commands); my $tablelist="--";
                                           #$tablelist = getCurrentTablesList();
                                           open(IN2,"$thisfilename");
                                           my @currentTables;
                                           if(IN2){
                                              @currentTables=<IN2>; 
                                            }
                                           close(IN2);
                                            my $newtablelist="";$i=0;$j=0;
                                            #now for each table in our selected annotation, gene
                                           while($i<@annot){
                                                $geneannot=$annot[$i];chomp($geneannot);
                                               if($geneannot=~/refFlatSquishMulti/i){
                                                  if(($organism=~/mouse/i)&&($db=~/mm9/i)){ #use Jesse's annotations only for mm9
                                                     if($j==0){
                                                        $tablelist="$geneannot"; ++$j;
                                                      }
                                                      else{
                                                        $tablelist.=":$geneannot";
                                                      }
                                                   }
                                               }
                                               else{
                                                @foo = grep(/$geneannot/, @currentTables);  
                                                if(@foo>0){
                                                   while(@foo>0){
                                                         $token=shift(@foo);
                                                         if($token=~/^\s*$geneannot\s*$/){
                                                            $token=~s/^\s*//; $token=~s/\s*$//;
                                                            if($j==0){
                                                               $tablelist="$token"; ++$j;
                                                            }
                                                            else{
                                                              $tablelist.=":$token";
                                                            }
                                                          }
                                                     }
                                                  }
                                                      
                                               }
                                                ++$i;
                                               }
                                           $end2=0 if(($ct==2)||($prev=~/<\/select>/i));
                                           $line2=$prev; print OUT2 "$organismg,$organism,$db,$tablelist\n";
                                         }
                                    }
                                  $more2=0;
                               }
                         }
                     }
                    elsif($line=~/<\/SELECT>/){
                       $end=0;$more=0;
                     }
                 }
          }
    }
 }
}
$tm = localtime;
my ($mday, $mon, $yday) = ($tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "getCurrentOrganisms.pl program complete : $month/$mday/$yday\n";
close(OUT);
if(OUT3){print OUT3 "org_version,$configFile\nannot_config,$annotationsConfig\n";}
close(OUT3);
#now create the symbolic link to organism_version config file
chdir($config_dir);
if(-f "current_organism_version"){system("rm current_organism_version");}
system("ln -s $configFile current_organism_version"); 




