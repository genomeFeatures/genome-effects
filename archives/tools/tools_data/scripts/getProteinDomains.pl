#!/usr/bin/perl

########################################################################
## This script collects all the protein domains of every gene
## listed in the input gene file or comand line 
#
#   Author: Lucie Hutchins
#   Department of Research, Bioinformatics
#   Dr. Joel Graber's Lab
#   The Jackson Laboratory
#
#   Implimentation date : January 2012
#
#  Input : list of gene names,one gene name per line or a list of genes
#          separated by a commas character if a comand line
#  Output: a tab-delimited file with all the protein domains for 
#          for every gene on the list
#
# Data sources: http://www.uniprot.org/
#               http://pfam.sanger.ac.uk/
#**********************************************************************
## set default variables

use DBI;
use vars qw ($opt_h $opt_n $opt_o $opt_f $opt_F);
use Getopt::Std;
use LWP 5.64;
#use Time::localtime;
use POSIX;

$organism="Mouse"; # default organism
$inputfile="";$names="";$outputfile="";
@geneList=();
my $browser = LWP::UserAgent->new;
   $browser->timeout(10);$browser->env_proxy;
   $browser->agent('Mozilla/5.0');

getopts('ho:n:f:F:');

if($opt_o){$organism = $opt_o;} if($opt_n){$names=$opt_n;}
if($opt_f){$inputfile=$opt_f;} if($opt_F){$outputfile=$opt_F;}
if(($opt_h)||(!$opt_n && !$opt_f)){
    print <<HELP;

 This script collects and displays all the protein domains of the specified gene(s)
 from UniProt and Sanger databases

Usage:

   perl getProteinDomains.pl [-o <organism (Mouse,human,rat)>][-n <gene(s) name>][-f <input file name>][-F <output file name>]

Arguments:
   -h  displays this help message
   -o  organism common name (Mouse,Human,Rat,...)
   -n  gene name (if you have more than one gene, separate gene names with a commas)
       use this option if you want to specify gene name(s) on the command line
   -f  a file containing the list of gene names to query.One name per line
       use this option if you you have all your gene names in a file
   -F  the output file name (default is stdout)

Examples:
a. single gene name on the command line and result stored in a file called pax6_proteinDomains.txt
  cmd: ./getProteinDomains.pl -n pax6 -F pax6_proteinDomains.txt

b. multiple genes on the command line and result stored in a file called myGenes_proteinDomains.txt
  cmd: ./getProteinDomains.pl -n pax6,lamb3,kit -F myGenes_proteinDomains.txt

c. Input is a file myGenes.txt and result stored in a file called myGenes_proteinDomains.txt
  cmd: /getProteinDomains.pl -f myGenes.txt -F myGenes_proteinDomains.txt

d. single gene name on the command line and result displayed on stdout
  cmd: ./getProteinDomains.pl -n pax6


HELP
exit;
}
if($names ne ""){ #genes were specified on the comand line
   if($names=~/,/){
      @geneList=split(",",$names);
   }
   else{@geneList[0]=$names;}
}
else{ #gene names are listed in a file
  if(-f "$inputfile"){
    open(IN,"$inputfile");
    if(!IN){print "Bad input file $inputfile\n";exit;}
    @geneList=<IN>;
  }
  else{ print $HELP; exit;}
}
$uniprot_url_prefix="http://www.uniprot.org/uniprot/?query=";
$uniprot_url_sufix="&sort=score";
$sanger_url="http://pfam.sanger.ac.uk/protein/";
if($opt_F){
   open(OUT,">$outputfile");
   if(!OUT){print "Could not open the output file :$outputfile -$!\n";exit;}
   print OUT "Organism\tGene\tUniprot-ProteinID\tUniprot-ProteinName\tSange-DomainID\tSanger-DomainName\tstart_aa\tend_aa\n";
}
else{
   print  "Organism\tGene\tUniprot-ProteinID\tUniprot-ProteinName\tSanger-DomainID\tSanger-DomainName\tstart_aa\tend_aa\n";
}
foreach my $gene(@geneList){
      chomp($gene); $gene=~s/^\s+//;$gene=~s/\s+$//;
      $uniport_url="$uniprot_url_prefix$gene+$organism";
      my $filename="proteinList".".txt"; my @filecontent;
      if(-f $filename){system("rm $filename");}
      my $response = $browser->get($uniport_url, ":content_file" => $filename);
      open(IN,"$filename");if(IN){@filecontent=<IN> ;}close(IN);
      $data="";$result="";$pageindex=0;$totalhit=0;$pageSize=25;$numberPages=0;
      @results=grep(/results for <strong>/,@filecontent);next if(@results<=0);
      $result=shift(@results); chomp($result);if($result=~/results\s+for\s+<strong>/){$result=$`;}
      if($result=~/<p id="status">/){$result=$';}$result=~s/<strong>//g;$result=~s/<\/strong>//g;
      if($result=~/(\d+\s*-\s*\d+\s+of)?\s*(\d+)/){$total=$2;}$numberPages=ceil($total/$pageSize);$ofset=0;
      while($pageindex<$numberPages){ #get the number of pages then collects domains on each page
         if($pageindex==0){$ofset=$pageindex;}
         else{$ofset=$pageindex*$pageSize; }
         ++$pageindex; #print "Processing page $pageindex offset is $ofset\n";
         #get the domains listed on this page
         $uniport_url="$uniprot_url_prefix$gene+$organism&offset=$ofset";
         my $response = $browser->get($uniport_url, ":content_file" => $filename);
         open(IN,"$filename");if(IN){@filecontent=<IN> ;}close(IN);
         $data="";$rows="";
         while(@filecontent>0){
              $line=shift(@filecontent);chomp($line); 
             if($line=~/organismHeader/){
                $line=~s/^\s+//; $line=~s/\s+$//;$data=$line;$more=1;
                while(($more)&&(@filecontent>0)){
                       $line=shift(@filecontent);chomp($line); $line=~s/^\s+//; $line=~s/\s+$//;
                       if($line=~/<\/table>/){$data.=$`;last;}
                       else{$data.=$line;}
                  } 
                 last;
               }
         }
         if($data ne ""){ #get rows
           @protein_rows=split(/<tr>/,$data); #print "We have ".@protein_rows." protein(s) on this page\n";
           foreach $row(@protein_rows){
                $proteinid="";$protein_name="";
                next if(!($row=~/$organism/i));
                if($row=~/addOrAppendCart\('(.+)'\)/){$proteinid=$1;}
                if($row=~/class="protein_names"><div class="short">(.*)<\/div><div class="long"/){$protein_name=$1;}
                $protein_name=~s/<strong>//i;$protein_name=~s/<\/strong>//i;
                #now get Sanger domains name and cordinates <table class="resultTable details"
                next if($proteinid eq "");
                $sanger_domains_url="$sanger_url$proteinid"; $domain_temp="domains.txt";
                if(-f $domain_temp){system("rm $domain_temp");}
                $response = $browser->get($sanger_domains_url, ":content_file" => $domain_temp);
                open(IN,"$domain_temp");if(IN){@filecontent=<IN> ;}close(IN);$dom="";
                while(@filecontent>0){
                    $line=shift(@filecontent);chomp($line);
                    if($line=~/<table class="resultTable details"/){
                       $line=~s/^\s+//; $line=~s/\s+$//;$dom=$line;$more=1;
                       while(($more)&&(@filecontent>0)){
                              $line=shift(@filecontent);chomp($line); $line=~s/^\s+//; $line=~s/\s+$//;
                              if($line=~/<\/table>/){$dom.=$`;last;}
                              else{$dom.=$line;}
                        }
                       last;
                    }
                  }
                  if($dom ne ""){
                    @domain_rows=split(/<\/tr>/,$dom); #print "We have ".@domain_rows." domains for this protein $proteinid\n";
                    foreach $domain(@domain_rows){
                        next if($domain=~/class="inactive"/);
                        if($domain=~/class="(.*)">(.*)<\/td><td><a href=".*">(.*)<\/a><\/td><td>(\d+)<\/td><td>(\d+)<\/td>/){
                          $domain_id=$1; $domain_name=$3;$start_aa=$4;$end_aa=$5;if($domain_id=~/class="(.+)/){$domain_id=$1;}
                          $rows.="$organism\t$gene\t$proteinid\t$protein_name\t$domain_id\t$domain_name\t$start_aa\t$end_aa\n";}
                    }
                  }
                
            }
         }
        else{print "This page has no protein data\n";}
        if($opt_F) {print OUT "$rows";}
        else{print "$rows";}
      }
 }
if($opt_F){close(OUT);}
print "Program complete\n";
exit(0);

