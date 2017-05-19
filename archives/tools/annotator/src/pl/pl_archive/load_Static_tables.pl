#!/usr/bin/perl

use DBI;
use Time::localtime;

#***********************************************************************************************
# This script loads new chromosomes, annotations,organisms, organism groups
# for a given assembly built into our database
# Tables affected: chromosome, chromosome_data, organism, gene_prediction,organism_group
#
# Input : three config files:
#   1. Config file 1: ucsc_organism_version_xxx.txt
#      a.organism_group
#      b.organism_name
#      c.ucsc_db_version
#      d.organism_build
# 
#   2.the file annotations_format.log contains the following fields:
#     a. annotation_name
#     b. organism
#     c. organism_group
#     d. a list of fieldCount,fields and the corresponding index
#
#   3.the file downloaded_annotations.log contains the list of all the annotations
#     we downloaded from ucsc for a given organism. One annotation per line
#     format:
#     a.organism_group
#     b.organism
#     c.gene_prediction name
#     d.generated file format :db_gene-prediction_date

# Author: Lucie N. Hutchins
#         Scientific Software Engineer
#
# Date : June 2010
#
#Usage : laod_geneannotations.pl ucsc_organism_version_xxx.txt annotations_format.log 
#        downloaded_annot.log
#********************************************************************************************
#

$config=shift||"/scratch/data/ucsc/ucsc_load_process/configuration/prog_generated_config.txt";
chomp($config);
open(CONF,"$config");
if(!CONF){
   print"Usage: perl load_Static_tables.pl config file\n";
   exit(1);
}
@filecontent=<CONF>; my($configpath,);
my $main_path="/scratch/data/ucsc"; # path to gene annotations
my $load_base_path="/scratch/data/ucsc/ucsc_load_process";
my ($orgv_configFile,$annoFortmatconfig,$downloadedAnnot_config);
while(@filecontent>0){
   $line=shift(@filecontent);chomp($line);
   my($variable,$base_path)=split(",",$line);
   if($variable=~/ucsc_load_base/){
      $load_base_path=$base_path;
   }
   elsif($variable=~/ucsc_org_base/){
      $main_path=$base_path;
   }
   elsif($variable=~/ucsc_table_list/){
      $downloadedAnnot_config=$base_path;
   }
   elsif($variable=~/org_version/){
      $orgv_configFile=$base_path;
   }
   elsif($variable=~/annotations_format/){
      $annoFortmatconfig=$base_path;
   }
  elsif($variable=~/graber_db_dev/){
       $dbname=$base_path;
  }
  elsif($variable=~/graber_db_hostdev/){
      $host=$base_path;
   }
   elsif($variable=~/graber_db_usr/){
      $user=$base_path;
   }
   elsif($variable=~/graber_db_pass/){
      $pwd=$base_path;
   }
}
close(CONF);
#$dbname="graber_transcriptdb";
#$host="demon.jax.org"; $user="lnh";$pwd="lucie";
# Set the path to configuration files
my $config_path="$load_base_path/configuration";
my $sql_path="$load_base_path/src/sql";
my $data_path="$load_base_path/data_load";
my %annotmap;
$tm = localtime;
my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
my $processlog="$config_path/data_load_log.log";

open(LOG,">>$processlog");
if(LOG){
  chomp($downloadedAnnot_config);chomp($orgv_configFile);chomp($annoFortmatconfig);
   print LOG "\n\nStarting laod_Static_tables.pl program:  $mday/$mon/$yday @ $hour:$min:$sec\n";
   #start data parsing process
   print LOG "Organism\tjob_type\tjob_starts\tjob_ends\n";
   #open the config fileS;
   open(OGV,"$orgv_configFile");open(DL,"$downloadedAnnot_config");
   if(!OGV ||!DL){ #if one of the config files has issues report and quit
      $log="Can not open the following file(s) :\n";
      $log.="$orgv_configFile," if(!OGV); $log.="$downloadedAnnot_config," if(!DL);
      print LOG "$log\n";
      close(LOG); exit(1);
   }
 my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user,$pwd);
 if(!$dbh){ print LOG "Could not connect to the database :$!\n";close(LOG); exit(1);} 
  my $getChr ="select chromosome_id,chromosome_name from chromosome";              # get current chromosomes list
  my $getOrgg ="select organism_group_id,organism_group_name from organism_group"; # get current organism_group list
  my $getOrg="select organism_id,organism_name from organism";                     #get current organisms list
  my $getOrgv ="select organism_version_id,organism_id,ucsc_db,coordinate_build from organism_version";
  my $getGannot="select gene_prediction_id,gene_prediction_name from gene_prediction";    #get current gene annotation list
  
  my $qh_chrlist    = $dbh->prepare($getChr)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_orgglist   = $dbh->prepare($getOrgg)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_orglist    = $dbh->prepare($getOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_orgvlist   = $dbh->prepare($getOrgv)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_gannotlist = $dbh->prepare($getGannot)or die "Couldn't prepare statement: " . $dbh->errstr;
  
  #get specific query 
  my $getThisChr="select chromosome_id from chromosome where chromosome_name=?";
  my $getThisOrgg="select organism_group_id from organism_group where organism_group_name=?";
  my $getThisOrg="select organism_id from organism where organism_name=?";
  my $getThisOrgv="select organism_version_id from organism_version where organism_id=?";
     $getThisOrgv.=" and ucsc_db=? and coordinate_build=?";
  my $getThisGannot="select gene_prediction_id from gene_prediction where gene_prediction_name=?";
  
  my $qh_gtchr      = $dbh->prepare($getThisChr)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_gtorgg     = $dbh->prepare($getThisOrgg)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_gtorg      = $dbh->prepare($getThisOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_gtorgv      = $dbh->prepare($getThisOrgv)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_gtgannot   = $dbh->prepare($getThisGannot)or die "Couldn't prepare statement: " . $dbh->errstr;
  
  #insert new items
  my $insertChr    ="insert ignore into chromosome(chromosome_name) values(?)";
  my $insertOrg    ="insert ignore into organism(organism_name,organism_group_id) values(?,?)";
  my $insertOrgg   ="insert ignore into organism_group(organism_group_name) values(?)";
  my $insertOrgv   ="insert ignore into organism_version(organism_id,ucsc_db,coordinate_build) values(?,?,?)";
  my $insertGannot ="insert ignore into  gene_prediction(gene_prediction_name) values(?)";
  
  my $qh_ichr   = $dbh->prepare($insertChr)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_iorg   = $dbh->prepare($insertOrg)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_iorgg  = $dbh->prepare($insertOrgg)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_iorgv  = $dbh->prepare($insertOrgv)or die "Couldn't prepare statement: " . $dbh->errstr;
  my $qh_igannot  = $dbh->prepare($insertGannot)or die "Couldn't prepare statement: " . $dbh->errstr;
 
## chromosome load
  $qh_drop_chrom_temp_table =$dbh->prepare("drop temporary table if exists chromosome_temp");
  $create_chrom_temp_table="create temporary table if not exists chromosome_temp(chromosome_id mediumint unsigned,";
  $create_chrom_temp_table.="chromosome_name varchar(255),index(chromosome_id),index(chromosome_name))";
  $qh_create_chrom_temp_table=$dbh->prepare($create_chrom_temp_table);
  $qh_load_chrom_temp_table=$dbh->prepare("LOAD DATA LOCAL INFILE ? INTO table chromosome_temp");
  $update_chrom_temp_table=" UPDATE chromosome_temp t, chromosome c ";
  $update_chrom_temp_table.=" set t.chromosome_id=c.chromosome_id where t.chromosome_id=0 ";
  $update_chrom_temp_table.=" and t.chromosome_name=c.chromosome_name";
  $qh_update_chrom_temp_table=$dbh->prepare($update_chrom_temp_table);
  $insert_new_chrom="INSERT INTO chromosome(chromosome_name) SELECT distinct chromosome_name";
  $insert_new_chrom.=" from chromosome_temp where chromosome_id=0";
  $qh_insert_new_chrom=$dbh->prepare($insert_new_chrom);

  $qh_chrom_row_count=$dbh->prepare("select count(*) from chromosome_temp");
  $qh_chrom_data_row_count=$dbh->prepare("select count(*) from chromosome_data_temp");

  $qh_drop_chrom_data_temp=$dbh->prepare("drop temporary table if exists chromosome_data_temp");
  $create_chrom_data_temp="create temporary table if not exists chromosome_data_temp(";
  $create_chrom_data_temp.="organism_version_id smallint,chromosome_id mediumint unsigned,";
  $create_chrom_data_temp.="chromosome_name varchar(255),fileName varchar(255),";
  $create_chrom_data_temp.="chrom_size int unsigned,is_found tinyint default 0,index(chromosome_id),";
  $create_chrom_data_temp.="INDEX(organism_version_id),index(organism_version_id,chromosome_id),index(chromosome_name))";
  $qh_create_chrom_data_temp=$dbh->prepare($create_chrom_data_temp);
  $qh_load_chrom_data_temp_table=$dbh->prepare("LOAD DATA LOCAL INFILE ? INTO table chromosome_data_temp");
  $update_chrom_data_temp_table="UPDATE chromosome_data_temp t, chromosome_temp c";
  $update_chrom_data_temp_table.=" set t.chromosome_id=c.chromosome_id where t.chromosome_name=c.chromosome_name";
  $qh_update_chrom_data_temp_table=$dbh->prepare($update_chrom_data_temp_table);

   $update_chrom_data_temp_table="UPDATE chromosome_data_temp t, chromosome_data c";
   $update_chrom_data_temp_table.=" set t.is_found=1 where t.organism_version_id=c.organism_version_id";
   $update_chrom_data_temp_table.=" and t.chromosome_id=c.chromosome_id";
   $qh_update_found_chrom_data_temp_table=$dbh->prepare($update_chrom_data_temp_table);
   $insert_chrom_data_temp_table=" INSERT INTO chromosome_data select organism_version_id,chromosome_id,";
   $insert_chrom_data_temp_table.=" fileName,chrom_size from chromosome_data_temp where is_found=0";
   $qh_insert_chrom_data_temp_table=$dbh->prepare($insert_chrom_data_temp_table);
  my(%chrmap,%orgmap,%orggmap,%orgvmap,%genpmap);
  # get the list of all the chromosomes
  $qh_chrlist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
  while(@row=$qh_chrlist->fetchrow_array()){
        $chr_id=$row[0]; $chr_name=lc($row[1]);$chr_name=~s/^\s+//;$chr_name=~s/\s+$//;
        $chrmap{"$chr_name"}=$chr_id;
   }
  #get the list of all organism_groups
  $qh_orgglist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
  while(@row=$qh_orgglist->fetchrow_array()){
        $orgg_id=$row[0]; $orgg_name=lc($row[1]);$orgg_name=~s/^\s+//;$orgg_name=~s/\s+$//;
        $orggmap{"$orgg_name"}=$orgg_id;
  }
  #get the list of all organisms
  $qh_orglist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
  while(@row=$qh_orglist->fetchrow_array()){
        $org_id=$row[0]; $org_name=lc($row[1]);$org_name=~s/^\s+//;$org_name=~s/\s+$//;
        $orgmap{"$org_name"}=$org_id;
  }
  #get the list of all organism versions
  $qh_orgvlist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
  while(@row=$qh_orgvlist->fetchrow_array()){
        $orgv_id=$row[0]; $org_id=$row[1];$ucsc_db=lc($row[2]);
        $coord_build=lc($row[3]);
        $ucsc_db=~s/^\s+//;$ucsc_db=~s/\s+$//;
        $coord_build=~s/^\s+//; $coord_build=~s/\s+$//;
        $orgmvap{$org_id}{$ucsc_db}{$coord_build}=$orgv_id;
  }
   #get the list of all gene_prediction
  $qh_gannotlist->execute() or die "Can't execute query: " . $dbh->errstr . "\n";
  while(@row=$qh_gannotlist->fetchrow_array()){
        $annot_id=$row[0]; $annot_name=lc($row[1]);
        $annot_name=~s/^\s+//;$annot_name=~s/\s+$//;
        $genpmap{$annot_name}=$annot_id;
  }
   # load the config files into memory
   my(%orgv_configFilemap,%annoFortmatconfigmap,%downloadedAnnot_configmap,$chromload,%orgvnewmap,%annot_map);
   #OGV |!DL |!ANO
   while($line=<OGV>){
       chomp($line);($organism_group,$organism,$db,$version)=split(",",$line);
       $orgv_configFilemap{"$organism"}{$db}{"orgg"}=$organism_group;
       $orgv_configFilemap{"$organism"}{$db}{"db"}=$db;
       $orgv_configFilemap{"$organism"}{$db}{"version"}=$version;
       $orgvnewmap{"$organism"}{$db}=1;
   }
   while($line=<DL>){
       chomp($line);($organism_group,$organism,$db,$annotation,$filename)=split(",",$line);
       $downloadedAnnot_configmap{lc($organism)}{$db}{$annotation}=$filename;
   }
    
 q{  while($line=<ANO>){ #ensGene:horse:mammal:equCab2:16,
       chomp($line);($annotation,$organism,$organism_group,$db,$fields)=split(":",$line);
       $annoFortmatconfigmap{$annotation}{lc($organism)}{$db}=$fields; #the first field is fieldcount
       $annot_map{lc($organism)}{$annotation}=1;
   }
  close(ANO);};
  close(DL);close(OGV);
  $chromload=0; #this flag is set if we are loading chromosome info
  #now for each organism, load organism/organism_group/organism_version tables ,load the chromInfo into chromosome
  for my $organism(sort keys(%orgv_configFilemap)){
      for my $db(%{$orgvnewmap{"$organism"}}){
          $organism_group=$orgv_configFilemap{"$organism"}{$db}{"orgg"};
          $ucsc_db=$orgv_configFilemap{"$organism"}{$db}{"db"};
          $coord_build=$orgv_configFilemap{"$organism"}{$db}{"version"};
          $dirname=lc($organism); $dirname=~s/^\s+//; $dirname=~s/\s+$//;$dirname=~s/\s+/-/;
          $organism_group=~s/^\s+//;$organism_group=~s/\s+$//; $orgg_id=0;$organism_id=0; $orgv_id=0;
          $organism=~s/^\s+//;$organism=~s/\s+$//;
          $ucsc_db=~s/^\s+//;$ucsc_db=~s/\s+$//;
          $coord_build=~s/^\s+//; $coord_build=~s/\s+$//;
          print LOG "Processing :$organism_group,$ucsc_db,$organism\n";
          # now check if this organism group is already in the database 
          $orgg_id=$orggmap{lc($organism_group)};        
          if(($orgg_id<=0)&&($organism_group=~/\w+/)){ #this is a new organism group
             $qh_iorgg->execute($organism_group);     #insert this new organism group then get the id
             $qh_gtorgg->execute($organism_group);    #now get the id of this organism_group
             if(@row=$qh_gtorgg->fetchrow_array()){
                $orgg_id=$row[0];
                $orggmap{lc($organism_group)}=$orgg_id;
             }
             else{
                print LOG "---\n$organism_group was not inserted\n";
             }
          }
         $org_id=$orgmap{lc($organism)};
         if(($org_id<=0)&&($organism=~/\w+/)){#this is a new organism
            $qh_iorg->execute($organism,$orgg_id); #insert this new organism then get the id
            $qh_gtorg->execute($organism);  
            if(@row=$qh_gtorg->fetchrow_array()){
               $org_id=$row[0];
               $orgmap{lc($organism)}=$org_id;
            } 
            else{
               print LOG "---\n$organism was not inserted\n";
            }
         }
        if($org_id>0){
            #now check the organism version
            $orgv_id=$orgmvap{$org_id}{lc($ucsc_db)}{lc($coord_build)};
            if(($orgv_id<=0)&&($ucsc_db=~/\w+/)){
                  $qh_iorgv->execute($org_id,$ucsc_db,$coord_build);
                  $qh_gtorgv->execute($org_id,$ucsc_db,$coord_build);  
                 if(@row=$qh_gtorgv->fetchrow_array()){
                    $orgv_id=$row[0];
                    $orgmvap{$org_id}{lc($ucsc_db)}{lc($coord_build)}=$orgv_id;
                 }
                 else{
                    print LOG "---\n$organism $ucsc_db,$coord_build was not inserted\n";
                 } 
             }
             #now process each annotation from the list for this organism
             $annot_id=0; $exist=0; $neworg=lc($organism);
             for my $annot_name(sort keys (%{$annot_map{$neworg}})){
                     next if($annot_name=~/chromInfo|author|cell|development|description|est_align|est_orientation/i);
                     next if($annot_name=~/gbCdnaInfo|geneName|library|mrnaClone|sex|sequence|source|tissue/i);
                     #remove leading and trailing space
                     $annot_name=~s/^\s+//;$annot_name=~s/\s+$//;
                     #check if this annotation exists in our db
                     $annot_id=$genpmap{lc($annot_name)};
                     if(($annot_id<=0)&&($annot_name=~/\w+/)){
                         $qh_igannot->execute($annot_name);#insert this new chrom then get the id
                         $qh_gtgannot->execute($annot_name);  
                         if(($annot_id)=$qh_gtgannot->fetchrow_array()){
                             $genpmap{lc($annot_name)}=$annot_id;
                         } 
                      }
                      if($annot_id<=0){
                         print LOG "---\n$organism:$annot_name was not inserted\n";
                      }
              } #end of foreach annotation
             #now load the chromInfo file of this organism 
             $chrominfo="chromInfo";
             $chrominfo_file=$downloadedAnnot_configmap{lc($organism)}{$ucsc_db}{$chrominfo};
             #$chrominfo_file="$chrominfofile";
             $chrom_file="$dirname"."_chrom.txt";$chrom_data_file="$dirname"."_chrom_data.txt";
             $chrom_data_file="$data_path/$chrom_data_file";
             $chrom_file="$data_path/$chrom_file";
             if(-f $chrominfo_file){       #we are loading the chromosome info for this organism
                 $chromload=1;$chromCount=0;$chromDataCount=0;
                 open(CHR,"$chrominfo_file") or die("Can't open $chrominfo_file : $!\n ");
                 if(CHR){
                     $header=<CHR>;chomp($header);$line="";@fields=split("\t",$header);
                     if(@fields>0){$index=0;
                        while(@fields>0){$field=shift(@fields);chomp($field);$line.=",$field>$index"; ++$index;}
                     }
                      @filecontent=<CHR>; #load chrom file in memory
                      #get this organism chrominfo file fields index
                      $chromindex=-1; $fileindex=-1; $sizeindex=-1;
                      if($line=~/,chrom>(\d+)/i){$chromindex=$1;}
                      if($line=~/,fileName>(\d+)/i){$fileindex=$1;}
                      if($line=~/,size>(\d+)/i){$sizeindex=$1;}
                      if($chromindex>=0 && $fileindex>=0){
                         open(CHRD,">$chrom_data_file");open(CHROM,">$chrom_file");
                         while($line=shift(@filecontent)){
                               chomp($line); next if($line=~/fileName/);
                               @linecontent=split("\t",$line);
                               $chrom=$linecontent[$chromindex] if($chromindex>=0 && $chromindex<@linecontent);
                               $fileName=$linecontent[$fileindex] if($fileindex>=0 && $fileindex<@linecontent);
                               $chromsize=$linecontent[$sizeindex] if($sizeindex>=0 && $sizeindex<@linecontent);
                               $chr_id=0;#check if this chromosome already exists in our db
                               if($chromindex>=0 && $chromindex<@linecontent){
                                  $chrom=~s/^\s+//;$chrom=~s/\s+$//;$chrom=~s/chr//i;                                             
                                  $chr_id=$chrmap{lc($chrom)};
                                  if($chr_id<=0){$chr_id=0;}
                                  if(CHROM){
                                     print CHROM "$chr_id\t$chrom\n";++$chromCount;
                                  }$exist=0;
                                  if(CHRD){ #generate the chromosome_data file
                                      if($orgv_id>0){
                                         print CHRD "$orgv_id\t$chr_id\t$chrom\t$fileName\t$chromsize\t$exist\n";
                                         ++$chromDataCount;
                                      }
                                  }
                                  
                              } #end of $chromindex<@linecontent
                         } #end of while loop
                        close(CHRD);close(CHROM);
                    }#end of $chromindex >=0 && $fileindex>=0
              }# end of if(CHR){
              #now load the chromosome data
             if((-f $chrom_data_file)||(-f $chrom_file)){
                   $qh_drop_chrom_temp_table->execute();$qh_create_chrom_temp_table->execute();
                   $qh_load_chrom_temp_table->execute($chrom_file);$qh_chrom_row_count->execute();
                   if($qh_chrom_row_count->rows>0){
                      ($rowcount)=$qh_chrom_row_count->fetchrow_array();
                      if($rowcount==$chromCount){
                         $qh_update_chrom_temp_table->execute();
                         $qh_insert_new_chrom->execute();$qh_update_chrom_temp_table->execute();
                         $qh_drop_chrom_data_temp->execute();$qh_create_chrom_data_temp->execute();
                         $qh_load_chrom_data_temp_table->execute($chrom_data_file);
                         $qh_chrom_data_row_count->execute();
                         if($qh_chrom_data_row_count->rows>0){
                            ($rowcount)=$qh_chrom_data_row_count->fetchrow_array();
                            if($rowcount==$chromDataCount){
                               $qh_update_chrom_data_temp_table->execute();
                               $qh_update_found_chrom_data_temp_table->execute();
                               $qh_insert_chrom_data_temp_table->execute();
                            }else{
                                 print LOG "Load Failed.$organism:$ucsc_db . ";
                                 print LOD "Only $rowcount of $chromDataCount lines were loaded for $chrom_data_file\n";}
                         }else{
                               print LOG "Load Failed.$organism:$ucsc_db. $chrom_data_file failed\n";}
                      }else{
                       print LOG "Load Failed.$organism:$ucsc_db . Only";
                       print LOG " $rowcount of $chromCount lines were loaded for $chrom_file\n";}
                   }else{
                       print LOG "Load Failed.$organism:$ucsc_db. $chrom_file failed\n";}
                   $qh_drop_chrom_temp_table->execute();$qh_drop_chrom_data_temp->execute();
                  print "$organism:$ucsc_db chromosome data loaded\n";
              }
        } #end of  if(-f $chrominfo_file){
     } #end of $org_id>0
    } #end of this version of the organism       
  } #end of loop organism
}
$tm = localtime;
 my ($sec,$min,$hour,$mday, $mon, $yday) = ($tm->sec,$tm->min,$tm->hour,$tm->mday, ($tm->mon)+1, ($tm->year)+1900);
   print LOG "End Static tables updates : $mon/$mday/$yday :$hour:$min:$sec\n";
close(LOG);

print "load_Static_tables.pl program complete\n";


