#!/usr/bin/perl

###########################################################################################################################
## insert_new_snps.pl
#  This script generates unique database id for every new SNP
#
# Input: a configuration file

#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#   Implimentation date : March 2012
#   Modification date : April 2012

#   Usage: perl insert_new_snps.pl -c <list_of_chromosome> -d <SNP config file>
#
#   where:
#    -h  displays this help message
#    -c  A file containing a list of chromosome names (one chromosome per line)  
#    -f  SNP config file
###########################################################################################################################
use DBI;
use Time::localtime;

use vars qw ($opt_h $opt_c $opt_f);
use Getopt::Std;

getopts('hc:f:');
if($opt_h||!$opt_f) {
    print <<HELP;

 This script rinserts new SNPs into the database

Usage:

   perl insert_new_snps.pl -c <list_of_chromosome> -f <SNP config file>

Arguments:
   -h  displays this help message
   -c  A file containing a list of chromosome names (one chromosome per line)  
   -f  SNP config file

Example: perl insert_new_snps.pl -c chromosome_list_file.txt -f snp_base_dir/snp_source.config
HELP
exit;
}

my $user ='lnh';
my $pwd  ='lucie98';

$dbname ='cgd_snpdb';
$host ='cgd-dev.jax.org';

my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pwd);
if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 

$get_source_id="select source_id from snp_source where source_name=?";
$qh_getSourceid=$dbh->prepare($get_source_id);
$get_strains="select strain_id from snp_strain_synonym where synonym_name=?";
$qh_getStrainId=$dbh->prepare($get_strains);
$get_strain="select strain_id from snp_strain where strain_name=?";
$qh_getStrain=$dbh->prepare($get_strain);
$get_max_strainid="select max(strain_id) from snp_strain";
$qh_getMaxStrainid=$dbh->prepare($get_max_strainid);
$get_strain_source="select* from snp_strain_by_source where strain_id=? and source_id=?";
$qh_get_strain_source=$dbh->prepare($get_strain_source);
my $insertStrain="insert ignore into snp_strain values(?,?)";
$qh_insertStrain=$dbh->prepare($insertStrain);
my $insertStrainSyn="insert ignore into snp_strain_synonym values(?,?)";
$qh_insertStrainSyn=$dbh->prepare($insertStrainSyn);
my $insertStrainBySource="insert ignore into snp_strain_by_source values(?,?)";
$qh_insertStrainBySource=$dbh->prepare($insertStrainBySource);

my $drop_snp_temp="drop temporary table if exists snp_position_temp";
my $qh_drop_snp_temp = $dbh->prepare($drop_snp_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $create_snp_temp="create temporary table snp_position_temp(snpid int unsigned default 0,
          chromosome_id tinyint default 0,bp_position int unsigned default 0,
         index(snpid),index(chromosome_id),index(chromosome_id,bp_position),index(bp_position))";
my $qh_create_snp_temp = $dbh->prepare($create_snp_temp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $load_this="load data local infile ? into table snp_position_temp";
my $qh_load_this = $dbh->prepare($load_this)or die "Couldn't prepare statement: " . $dbh->errstr;

my $update_this="update snp_position_temp t, snp_position p ";
   $update_this.=" set t.snpid=p.snpid where t.chromosome_id=p.chromosome_id and t.bp_position=p.bp_position";
my $qh_update_this = $dbh->prepare($update_this)or die "Couldn't prepare statement: " . $dbh->errstr;

my $insert_snp="insert into snp_position(chromosome_id,bp_position) ";
   $insert_snp.="select distinct chromosome_id,bp_position from snp_position_temp where snpid=0 ";
my $qh_insert_snp = $dbh->prepare($insert_snp)or die "Couldn't prepare statement: " . $dbh->errstr;

my $get_snp_rowcount="select count(*) as rowcount from snp_position_temp ";
my $qh_snp_rowcount = $dbh->prepare($get_snp_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;

my $get_snpid_rowcount="select count(*) as rowcount from snp_position_temp where snpid>0";
my $qh_snpid_rowcount = $dbh->prepare($get_snpid_rowcount)or die "Couldn't prepare statement: " . $dbh->errstr;

my $analyze="analyze table snp_position,snp_accession,snp_by_source";
$analyze="analyze table snp_accession";
my $qh_analyze = $dbh->prepare($analyze)or die "Couldn't prepare statement: " . $dbh->errstr;
 
########################################################################################
#$local_id\t$snpid\t$chr_id\t$bp_pos\t$id_type\t$source_id\is-found
$qh_delete_acc=$dbh->prepare("delete from snp_accession where source_id=?");
$qh_drop_acc_temp=$dbh->prepare("drop temporary table if exists snp_accession_temp");
my $create_acc_temp="create temporary table snp_accession_temp(snpid int unsigned default 0,
            accession_id char(20),chromosome_id tinyint default 0,
            bp_position int unsigned default 0,snpid_id_type tinyint default 0,
            source_id smallint default 0,is_found tinyint default 0,
            index(chromosome_id,bp_position),
           index(snpid),index(accession_id),index(snpid_id_type),index(source_id))";
my $qh_create_acc_temp = $dbh->prepare($create_acc_temp)or die "Couldn't prepare statement: ".$dbh->errstr;
my $load_acc_temp="load data local infile ? into table snp_accession_temp";
my $qh_load_acc = $dbh->prepare($load_acc_temp)or die "Couldn't prepare statement: ".$dbh->errstr;

my $qh_acc_rowcount=$dbh->prepare("select count(*) as rowcount from snp_accession_temp");
my $update_acc="update snp_accession_temp t, snp_position p set t.snpid=p.snpid";
   $update_acc.=" where t.chromosome_id=p.chromosome_id and t.bp_position=p.bp_position";
my $qh_update_acc_snp=$dbh->prepare($update_acc);

my $update_acc_temp="update snp_accession_temp t, snp_accession p set t.is_found=1 ";
   $update_acc_temp.=" where t.snpid=p.snpid and t.accession_id=p.accession_id and t.source_id=p.source_id";
my $qh_update_acc_temp = $dbh->prepare($update_acc_temp)or die "Couldn't prepare statement: ".$dbh->errstr;

my $insert_snp_acc="insert into snp_accession(accession_id,snpid_id_type,snpid,source_id) ";
   $insert_snp_acc.="select distinct accession_id,snpid_id_type,snpid,source_id from snp_accession_temp where is_found=0 ";
my $qh_insert_acc= $dbh->prepare($insert_snp_acc)or die "Couldn't prepare statement: ".$dbh->errstr;
### we also need to insert these accessions into the search engine cgd_items
#my $qh_create_cgd_items=$dbh->prepare("create temporary table if not exists cgd_items_temp(";
############################################################################################
#print SC "$local_id\t$chr_id\t$bp_pos\t$source_id\t0\n";
$qh_drop_snpbs_temp=$dbh->prepare("drop temporary table if exists snp_by_source_temp");
my $create_snpbs_temp="create temporary table snp_by_source_temp(snpid int unsigned default 0,
            chromosome_id tinyint default 0,bp_position int unsigned default 0,
            source_id smallint default 0,is_found tinyint default 0,
           index(snpid),index(chromosome_id,bp_position),index(source_id))";
my $qh_create_snpbs_temp = $dbh->prepare($create_snpbs_temp)or die "Couldn't prepare statement: ".$dbh->errstr;
my $load_snpbs_temp="load data local infile ? into table snp_by_source_temp";
my $qh_load_snpbs = $dbh->prepare($load_snpbs_temp)or die "Couldn't prepare statement: ".$dbh->errstr;

my $qh_snpbs_rowcount=$dbh->prepare("select count(*) as rowcount from snp_by_source_temp");
my $update_snpbs="update snp_by_source_temp t, snp_position p set t.snpid=p.snpid";
   $update_snpbs.=" where t.chromosome_id=p.chromosome_id and t.bp_position=p.bp_position";
my $qh_update_snpbs_snp=$dbh->prepare($update_snpbs);

my $update_snpbs_temp="update snp_by_source_temp t, snp_by_source p set t.is_found=1 ";
   $update_snpbs_temp.=" where t.snpid=p.snpid and t.source_id=p.source_id";
my $qh_update_snpbs_temp = $dbh->prepare($update_snpbs_temp)or die "Couldn't prepare statement: ".$dbh->errstr;

my $insert_snp_snpbs="insert into snp_by_source(snpid,source_id) ";
   $insert_snp_snpbs.="select distinct snpid,source_id from snp_by_source_temp where is_found=0 ";
my $qh_insert_snpbs= $dbh->prepare($insert_snp_snpbs)or die "Couldn't prepare statement: ".$dbh->errstr;


###############################
# get command line arguments 
$chr_list_file=$opt_c;$config_file=$opt_f;
#now validate the config file and get global variables
if(!(-f $config_file)){
  print STDERR "Bad configuration file name :$config_file\n";
  exit(1);
}
open(CONF,"$config_file");
if(!CONF){
  print STDERR "Could not open the configuration file $config_file:$!\n";
  exit(1);
}
@config=<CONF>;@variable=grep(/=/,@config);
my ($snp_dir,$source,$file_type,$output_dir,$source_id);my $count =0;$source_id=0;
$bp_pos_index=-1;$chrom_index=-1;$end_strain_index=-1;$start_strain_index=-1;
$file_prefix="";$file_sufix="";$all_in_one=0;$snp_file_name="";$source_abrev="";
$ssid_index=-1;
foreach my $line(@variable){# set global variables from the config file
  chomp($line); 
  if($line=~/SNP_SOURCE=(.+)$/){$source=$1;}
  if($line=~/SOURCE_ABREV=(.+)$/){$source_abrev=$1;}
  if($line=~/SNP_BASE_DIR=(.+)$/){$snp_dir=$1;}
  if($line=~/PIPELINE_DIR=(.+)$/){$output_dir=$1;}

  if($line=~/SNPID_INDEX=(.+)$/){$snpid_index=$1;}
  if($line=~/SNP_CHROM_INDEX=(.+)$/){$chrom_index=$1;}
  if($line=~/SNP_BP_POSITION_INDEX=(.+)$/){$bp_pos_index=$1;}
  if($line=~/SNP_FIRST_STRAIN_INDEX=(.+)$/){$start_strain_index=$1;}
  if($line=~/SNP_LAST_STRAIN_INDEX=(.+)$/){$end_strain_index=$1;}
  if($line=~/SNP_SSID_INDEX=(.+)$/){$ssid_index=$1;}
  if($line=~/SNP_RSID_INDEX=(.+)$/){$rsid_index=$1;}

  if($line=~/ALL_IN_ONE=(.+)$/){$all_in_one=$1;}
  if($line=~/SNP_FILE_NAME=(.+)$/){$snp_file_name=$1;}
  if($line=~/SNP_FILE_TYPE=(.+)$/){$file_type=$1;}
  if($line=~/SNP_FILE_PREFIX=(.+)$/){$file_prefix=$1;} 
  if($line=~/SNP_FILE_SUFIX=(.+)$/){$file_sufix=$1;} 
}
#Check config file settings
if(($bp_pos_index==-1)||($chrom_index==-1)){
  print STDERR "One of the following indice has a negative value in the config file:start_strain_index,bp_pos_index,chrom_index\n"; 
  exit(1);
}
if(!(-d $snp_dir)|| !(-d $output_dir)){
  print STDERR "One of the following is not a directory:$snp_dir or $output_dir.Check the config file settings\n"; 
  exit(1);
}
my($line,@linecontent,$chr,$chr_id);
if($snp_dir =~/(.*)\/$/){ 
        $snp_dir=$1; #format the directory name
}
if($output_dir =~/(.*)\/$/){ #format the directory name
        $output_dir=$1;
}
@chromosomes=(); 
if($all_in_one==0){
    open(IN,"$chr_list_file")or die "Can't open file:$!\n";@chromosomes=<IN>;close(IN);
 }
 else{$chromosomes[0]=$snp_file_name;}
open(LOG,">insert_log.log");
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Starting load process :  $mon/$mday/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
$qh_getSourceid->execute($source);
if($qh_getSourceid->rows<=0){ #source not found in our database
   print STDERR "The specified source was not found in our database\n"; 
   exit(1);
}
else{($source_id)=$qh_getSourceid->fetchrow_array();}
#delete the previous accessions of this source
$qh_delete_acc->execute($source_id);
while(@chromosomes>0){   # process file
   $chr=shift(@chromosomes);chomp($chr);print "Processing $chr;\n";
   if($all_in_one==0){
     if($chr=~/^0(\d)$/){$chr_id=$1;} 
     else{$chr_id=$chr;}  $chr_temp=$chr;
     $chr_id=20 if($chr =~/X/i);$chr_id=21 if($chr =~/Y/i);
     $chr_id=22 if($chr =~/M/i);$chr_id=23 if($chr =~/U/i);
     $snp_file_name="";$chr_temp="M" if(($chr =~/M/i)&&($source_id==20));
     $snp_file_name=$snp_dir."/$file_prefix".$chr_temp."$file_sufix";
  }
  else{$snp_file_name="$snp_dir/$snp_file_name";}
  if(!(-f $snp_file_name)){print "The file $snp_file_name does not exists\n";next;}
  open(SNP,"$snp_file_name") or die "Can't open $snp_file_name:$!";
  $line=<SNP>;chomp($line);  #get the file header line
  my(@linecontent);$delim="\t";# default file type 
  if($file_type==2){$delim=",";} # commas separated file
  @linecontent= split("$delim",$line); 
  if($end_strain_index==-1){$end_strain_index= @linecontent;--$end_strain_index;}
  if($end_strain_index<$start_strain_index){
     print "Bad strain index $end_strain_index<$start_strain_index\n";next;
  }
  $i=$start_strain_index;#now check if these strains already exists into our database
 # while($i<=$end_strain_index){ #I will restore this later
   while($i<$start_strain_index){
      $strain_name=$linecontent[$i];$strain_id=0;$strain_name=~s/^\s+//;$strain_name=~s/\s+$//;
      $strain_name=~s/,$//;$qh_getStrain->execute($strain_name);#check first in snp_strain table
      if($qh_getStrain->rows>0){
         ($strain_id)=$qh_getStrain->fetchrow_array();
         $qh_getStrainId->execute($strain_name); #check in synonym table
         if($qh_getStrainId->rows<=0){$qh_insertStrainSyn->execute($strain_id,$strain_name);}
      }
      else{ #new strain or synonymous
         $qh_getStrainId->execute($strain_name); #check first in snp_synonym
         if($qh_getStrainId->rows>0){($strain_id)=$qh_getStrainId->fetchrow_array();} #strain exists
         else{ #insert new strain
            $qh_getMaxStrainid->execute();
            if($qh_getMaxStrainid->rows>0){
               ($strain_id)=$qh_getMaxStrainid->fetchrow_array();++$strain_id;
            }
            $qh_insertStrain->execute($strain_id,$strain_name);$qh_getStrain->execute($strain_name);
            if($qh_getStrain->rows>0){
               ($strain_id)=$qh_getStrain->fetchrow_array();
               $qh_insertStrainSyn->execute($strain_id,$strain_name);
            }
         }
      }
      $qh_get_strain_source->execute($source_id,$strain_id);
      if($qh_get_strain_source->rows<=0){ ##insert this strain for this source if does not already exists 
         $qh_insertStrainBySource->execute($source_id,$strain_id);
      }
      $strains{$i}=$strain_id; #a dictionary mapping strain index to strain id   
     if($skip_one==0){++$i; }
     else{$i+=2;} 
   } 
  ########################################
  $count=0;@row_map=();$rows=0;@acc_map=();@snp_bys=();$main_rows=0;$acc_rows=0;
  $sbs_rows=0;$snp_main_file=$output_dir."/snp_position_temp_file.txt";
  $snp_acc_file=$output_dir."/snp_accession_temp_file.txt";
  $snp_bySource_file=$output_dir."/snp_bySource_temp_file.txt";
  while($line=<SNP>){
      chomp($line);$count +=1;$id="";$bp_pos=0;$local_id=0;
      @linecontent= split("$delim",$line);
      if(@linecontent>0){
         $chr_id=0;$local_b6allele;$local_snpallele;
         $bp_pos=$linecontent[$bp_pos_index];$chr=$linecontent[$chrom_index];$chr=~s/^chr//i;
         if($chr=~/^0(\d)$/){$chr_id=$1;}elsif($chr=~/(\d+)/){$chr_id=$1;}
         else{$chr_id=$chr;
            $chr_id=20 if($chr =~/X/i);$chr_id=21 if($chr =~/Y/i);
            $chr_id=22 if($chr =~/M/i);$chr_id=23 if($chr =~/U/i);
         }  
         if(($bp_pos=~/\d+/)&&($chr_id>0)){
             $bp_pos=~s/\s*//g;$local_id=0;
             #push(@snp_bys,"$chr_id-$bp_pos-$source_id");push(@row_map,"$chr_id-$bp_pos");
             if($ssid_index>=0){$ssid=$linecontent[$ssid_index];
               push(@acc_map,"$local_id=$ssid=$chr_id=$bp_pos=2=$source_id");
             }
            %rsidmap=();
             if($rsid_index>=0){
               $rsid=$linecontent[$rsid_index];
               push(@acc_map,"$local_id=$rsid=$chr_id=$bp_pos=1=$source_id");
               $rsidmap{$rsid}=1;
             }
             if(($snpid_index>=0)){ #generate accession ids
                $snpid=$linecontent[$snpid_index];@snpid=split(/\|/,$snpid);
                if(@snpid>0){
                   while(@snpid){$snpid=shift(@snpid); 
                      if($snpid=~/rs\d+/){
                         if(!exists($rsidmap{$snpid})){
                            push(@acc_map,"$local_id=$snpid=$chr_id=$bp_pos=1=$source_id");
                            $rsidmap{$snpid}=1;
                          }
                       }
                      else{ 
                        if($snpid ne ""){push(@acc_map,"$local_id=$snpid=$chr_id=$bp_pos=3=$source_id");}
                        else{push(@acc_map,"$local_id=$source_abrev$chr_id-$bp_pos=$chr_id=$bp_pos=3=$source_id");}
                      }
                    }
                }
                else{ 
                   if($snpid=~/rs\d+/){
                       if(!exists($rsidmap{$snpid})){
                          push(@acc_map, "$local_id=$snpid=$chr_id=$bp_pos=1=$source_id");
                          $rsidmap{$snpid}=1;
                       }
                   }
                   if(($snpid_index>=0)&&($snpid ne "")){
                       push(@acc_map,"$local_id=$snpid=$chr_id=$bp_pos=3=$source_id")if(!($snpid=~/rs\d+/));
                    }
                   else{push(@acc_map,"$local_id=$source_abrev$chr_id-$bp_pos=$chr_id=$bp_pos=3=$source_id");}       
                }
            }
        }
         ################
          if($count%100000==0){
            open(ACC,">$snp_acc_file");#open(SC,">$snp_bySource_file"); open(MAIN,">$snp_main_file");$main_rows=0;
            $acc_rows=0;$sbs_rows=0;
            if(MAIN){
                foreach my $key(@row_map){
                  ($chr_id,$pos)=split("-",$key);
                   print MAIN "0\t$chr_id\t$pos\n";$main_rows+=1;
                }
             }
            close(MAIN);
             if(-f $snp_main_file){
                $qh_drop_snp_temp->execute();$qh_create_snp_temp->execute();
                $qh_load_this->execute("$snp_main_file");
                $qh_snp_rowcount->execute();$id_rows=0;
                if($qh_snp_rowcount->rows>0){
                   ($load_rows)=$qh_snp_rowcount->fetchrow_array();
                   if($load_rows==$main_rows){ #load was success
                      $qh_update_this->execute();
                      $qh_snpid_rowcount->execute();($id_rows)=$qh_snpid_rowcount->fetchrow_array();
                      $qh_insert_snp->execute() ;$qh_update_this->execute();
                      print LOG "File snp_pos: $load_rows of $main_rows ";
                      print LOG " were loaded -$id_rows of these SNPs are already in our db\n";
                   }
                }
             
            }
            #load the accessions 
           if(ACC){
                foreach my $key(@acc_map){
                  ($local_id,$snpid,$chr_id,$bp_pos,$id_type,$source_id)=split("=",$key);
                   print ACC "$local_id\t$snpid\t$chr_id\t$bp_pos\t$id_type\t$source_id\t0\n";$acc_rows+=1;
                }
            }
            close(ACC);
           if(-f $snp_acc_file){
                $qh_drop_acc_temp->execute();$qh_create_acc_temp->execute();
                $qh_load_acc->execute("$snp_acc_file");
                $qh_acc_rowcount->execute();$rows=0;
                if($qh_acc_rowcount->rows>0){
                   ($load_rows)=$qh_acc_rowcount->fetchrow_array();
                   if($load_rows==$acc_rows){ #load was success
                      $qh_update_acc_snp->execute() or die "load fail:".mysql_error(); 
                      $qh_update_acc_temp->execute();
                      $qh_insert_acc->execute() or die "load fail:".mysql_error();
                     print LOG "File ACCESSION:  $load_rows OF $main_rows were loaded \n";
                   }
                }
              
            }
           ###load SNP by source
           #if(SC){ $local_id=0;$sbs_rows=0;
           #     foreach my $key(@snp_bys){
           #       ($chr_id,$bp_pos,$source_id)=split("-",$key);
           #        print SC "$local_id\t$chr_id\t$bp_pos\t$source_id\t0\n";$sbs_rows+=1;
           #     }
           # }
           #close(SC);
           #if(-f $snp_bySource_file){
           #     $qh_drop_snpbs_temp->execute();$qh_create_snpbs_temp->execute();
           #     $qh_load_snpbs->execute("$snp_bySource_file");$qh_snpbs_rowcount->execute();
           #     if($qh_snpbs_rowcount->rows>0){
           #        ($load_rows)=$qh_snpbs_rowcount->fetchrow_array();
           #        if($load_rows==$sbs_rows){ #load was success
           #           $qh_update_snpbs_snp->execute();$qh_update_snpbs_temp->execute();$qh_insert_snpbs->execute() ;
           #           print LOG "File SNP by source: $load_rows of $sbs_rows were loaded \n";
           #        }
           #     }
              
           # }
            print"$count processed \n";$main_rows=0;$acc_rows=0;$sbs_rows=0;@row_map=();@acc_map=();@snp_bys=();
        } #end of count%10000
         #################################
     } #end of @linecontent>0 check
   }#end of this chromosome SNP file loop
    #insert the last
 ############################
 q{ if(@row_map>0){
    open(MAIN,">$snp_main_file");$main_rows=0;
    if(MAIN){
       foreach my $key(@row_map){
           ($chr_id,$pos)=split("-",$key);
            print MAIN "0\t$chr_id\t$pos\n";$main_rows+=1;
        }
     }
     close(MAIN);
     if(-f $snp_main_file){ $id_rows=0;
          $qh_drop_snp_temp->execute();$qh_create_snp_temp->execute();
          $qh_load_this->execute("$snp_main_file") or die "load fail:".mysql_error();$qh_snp_rowcount->execute();
          if($qh_snp_rowcount->rows>0){
             ($load_rows)=$qh_snp_rowcount->fetchrow_array();
              if($load_rows==$main_rows){ #load was success
                 $qh_update_this->execute();
                 $qh_snpid_rowcount->execute();($id_rows)=$qh_snpid_rowcount->fetchrow_array();
                 $qh_insert_snp->execute();$qh_update_this->execute();
                 print LOG "File snp_pos: $load_rows of $main_rows ";
                 print LOG " were loaded -$id_rows of these SNPs are already in our db\n";
              }
          }
       }print "Inserting $main_rows of the main last segment\n";
   }};
   if(@acc_map>0){
      open(ACC,">$snp_acc_file");
      if(ACC){
         foreach my $key(@acc_map){
             ($local_id,$snpid,$chr_id,$bp_pos,$id_type,$source_id)=split("=",$key);
             print ACC "$local_id\t$snpid\t$chr_id\t$bp_pos\t$id_type\t$source_id\t0\n";$acc_rows+=1;
         }
       }
      close(ACC);
      if(-f $snp_acc_file){
          $qh_drop_acc_temp->execute();$qh_create_acc_temp->execute();
          $qh_load_acc->execute("$snp_acc_file");
          $qh_acc_rowcount->execute();
          if($qh_acc_rowcount->rows>0){
            ($load_rows)=$qh_acc_rowcount->fetchrow_array();
             if($load_rows==$acc_rows){ #load was success
                $qh_update_acc_snp->execute(); 
                $qh_update_acc_temp->execute();
                $qh_insert_acc->execute();
                print LOG "File ACCESSION:  $load_rows OF $main_rows were loaded \n";
              }
          }
         
      }print "Inserting $acc_rows of the accession last segment\n";
   }
   q{if(@snp_bys>0){
      open(SC,">$snp_bySource_file");
      if(SC){$local_id=0;$sbs_rows=0;
         foreach my $key(@snp_bys){
            ($chr_id,$bp_pos,$source_id)=split("-",$key);
            print SC "$local_id\t$chr_id\t$bp_pos\t$source_id\t0\n";$sbs_rows+=1;
          }
       }
       close(SC);
       if(-f $snp_bySource_file){
           $qh_drop_snpbs_temp->execute();$qh_create_snpbs_temp->execute();
           $qh_load_snpbs->execute("$snp_bySource_file");$qh_snpbs_rowcount->execute();
           if($qh_snpbs_rowcount->rows>0){
             ($load_rows)=$qh_snpbs_rowcount->fetchrow_array();
              if($load_rows==$sbs_rows){ #load was success
                 $qh_update_snpbs_snp->execute();$qh_update_snpbs_temp->execute();$qh_insert_snpbs->execute() ;
                 print LOG "File SNP by source: $load_rows of $sbs_rows were loaded \n";
              }
           }
       }
      print "Inserting $sbs_rows of the snp by source last segment\n";
   }};
   print "Done with chromosome $snp_file_name, $count lines in file,$source_id\n";  
 } #end of chromosome loop
$qh_analyze->execute();
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
print LOG "\n*************************************************************\n";
print LOG "Program Ends:  $mon/$mday/$yday @ $hour:$min:$sec \n";
print LOG "\n*************************************************************\n";
close(LOG);
system("rm $output_dir/*_temp_file.txt");
print"program completed\n";
 
exit(0);

