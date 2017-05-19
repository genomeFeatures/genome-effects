 <?php
 
function processQueryCSV($datatype,$snp_loc,$genotypeChek,$selectedSources,$confScore,$selectedStrains,$rowcount,$page,$fields,$querystring){  
   /////////////////////////////////
   ///////////////////////////////// 
   global $SNP_HOME; global  $SNP_SEARCH_HOME;global $CGD_HOME;global $UNIT;global $CGDTYPE;
   if(empty($CGDTYPE))loadCgdType();
   $htmloutput="";$chr=$datatype{"chr"};$startpos=$datatype{"startpos"};
   $endpos=$datatype{"endpos"};$gene_id= $datatype{"geneid"};$query_single=$datatype{"query_single"};
   $transcript_id= $datatype{"transcriptid"}; $item_type_id=$datatype{"qtype"};
   $sourcedata=getSNPsource($selectedSources);$sourcelist=$sourcedata{"sourcelist"};
   $strainsList=getSelectedStrainList($selectedStrains);
   $is_imputed= $sourcedata{"imputedflag"};$qtype_name=$CGDTYPE{$item_type_id}{"type_name"};
   $chrStr=$chr;if($chr==20)$chrStr="X"; if($chr==21)$chrStr="Y"; if($chr==22)$chrStr="Mt";
  /////////////////
    $sources=split(":",$selectedSources);
   $strains=split(":",$selectedStrains);
  //////////////////
   $querySum="Query Summary:\015\012";
   $querySum.="SNP Source(s) : $sourcelist\015\012";
   $querySum.="$qtype_name:$query_single\015\012";
   $sources=split(":",$selectedSources);
   if($item_type_id>0)
       $querySum.="Genomic Location: Chr$chrStr: $startpos - $endpos bp\015\012";
   $snplocationdata=setSNPloc($snp_loc); // get the user specified genomic region and annotations
   $locCount=$snplocationdata[0];$snp_location=$snplocationdata[1];$location_name=$snplocationdata[2];
   $is_intergenic=$snplocationdata[3];
   $SNPlist= getSNPsIDs($con,$chr,$startpos,$endpos,$selectedSources,$snp_location,
                       $is_intergenic,$gene_id,$transcript_local_id,$locCount);
   $querySum.="Annotations :$location_name\015\012"; 
   $querySum.="SNPs Match : $SNPlist[0]\015\012";  
   $querySum.="---------------------------------------------------------------------------\015\012";
   if($SNPlist[0]<=0){ // no SNPs found within the specified genomic region
       print"$querySum \015\012 Zero result , our database does not 
              have SNPs within the specified genomic region for the selected source(s)\015\012";  
   }
   else{
     print"$querySum"; $strain_names="";$currenttotalSNP=$SNPlist[0];
     $header="\015\012chrom,bp_position";
      //now display strains header
     $j=3;$strainheader="";$strainIndex;
     if(is_array($strainsList)){
        while(list($strain_id,$strain_name)=each($strainsList)){
            $header.=",$strain_name";
            $strainIndex[$j]=$strain_id;
           ++$j;
        }
      }
      $header.=",SNPID,Alleles";
      if($fields{"v_type"}==1) $header.=",Variation_Type";
      if($fields{"g_name"}==1) $header.=",MGI_Gene";
      if($fields{"f_class"}==1) $header.=",function_class";
      if($fields{"cpg"}==1) $header.=",Is_CpG_site";
      if($fields{"m_allele"}==1) $header.=",MinorAllele%";
      if($fields{"m_genotype"}==1) $header.=",MissingGeno%";
      $header.="\015\012";
      print"$header";$startindex=0;
      //now get for every snp on the list, display the data summary
      displayRowsCSV($SNPlist,$datatype,$selectedStrains,$is_imputed,$selectedSources,$strainIndex,$snp_loc,$fields);
   }
  
 }
/*************************************************************** 
  This is the main engine
  Author: Lucie Hutchins, Scientific Software Engineer
  Date : December 16 2009
  Modified: July 2010
***************************************************************/
function processQuery($datatype,$snp_loc,$genotypeChek,$selectedSources,$confScore,$selectedStrains,$rowcount,$page,$fields,$querystring){     
   global $SNP_HOME; global  $SNP_SEARCH_HOME;global $CGD_HOME;global $UNIT;global $CGDTYPE;
   if(empty($CGDTYPE))loadCgdType();
   $htmloutput="";$chr=$datatype{"chr"};$startpos=$datatype{"startpos"};
   $endpos=$datatype{"endpos"};$gene_id= $datatype{"geneid"};$query_single=$datatype{"query_single"};
   $transcript_id= $datatype{"transcriptid"}; $item_type_id=$datatype{"qtype"};
   $sourcedata=getSNPsource($selectedSources);$sourcelist=$sourcedata{"sourcelist"};
   $strainsList=getSelectedStrainList($selectedStrains);
   $is_imputed= $sourcedata{"imputedflag"};$qtype_name=$CGDTYPE{$item_type_id}{"type_name"};
   $chrStr=$chr;if($chr==20)$chrStr="X"; if($chr==21)$chrStr="Y"; if($chr==22)$chrStr="Mt";
   $querySum="<table  cellspacing='2' cellpadding='2' style='background:#ffffff; width:94%;margin:auto;font-size:0.9em;padding-bottom:1em;'>
               <tr>
                 <td style='width:auto;background:#efffff;border-right:solid 1px  #B2B4B6;height:auto;'>
                    <div class='sr' style='width:auto;overflow:hidden;margin:auto;padding:0.em;height:auto;'>
                  <h2>Query Summary:</h2>";
   $querySum.="<ul style='height:10em;'><li><b>SNP Source(s) : </b>$sourcelist</li>";
   $querySum.=" <li><b>$qtype_name:</b> $query_single</nobr></li>";
   /************ set filter option and query summary ****************/
   $sg=$fields{"sg"};
   $filterlink="$SNP_SEARCH_HOME?query_single=$query_single&amp;sg=$sg&amp;genotype=$genotypeChek&amp;confScore=$confScore";
   $filterlink .="&amp;rcount=$rowcount&amp;source=$selectedSources";
   $sources=split(":",$selectedSources);
   foreach($sources as $source_id){$source_id=trim($source_id);$filterlink.="&amp;source=$source_id";}
   $strains=split(":",$selectedStrains);
   //if(!isset($_SESSION['genome'])) 
   //          loadQuerySNPsInSession($chr,$startpos,$endpos); 
   if(is_array($strains)&&!empty($strains[0]))
      foreach($strains as $strain_id)$filterlink .="&amp;st=$strain_id";
   if($item_type_id>0)
       $querySum.=" <li><b>Genomic Location:</b> Chr$chrStr: $startpos - $endpos bp</li>";
   $snplocationdata=setSNPloc($snp_loc); // get the user specified genomic region and annotations
   $locCount=$snplocationdata[0];$snp_location=$snplocationdata[1];$location_name=$snplocationdata[2];
   $is_intergenic=$snplocationdata[3];
   $SNPlist= getSNPsIDs($con,$chr,$startpos,$endpos,$selectedSources,$snp_location,
                       $is_intergenic,$gene_id,$transcript_local_id,$locCount);
   $querySum.="<li><b>Annotations :</b> $location_name</li>"; 
   $querySum.="<!-- <li><b>SNPs Match :</b> $SNPlist[0]</li>--></ul></div></td>";  
   // get the filters
   $filters="<td style='width:15em;background:#efffff;border-right:solid 1px #B2B4B6;'>";
          //<div class='sr' style='margin:auto;padding:0.em;'><h2>Filter Results:</h2>";
   $totalIntergenic=0;$totalIntronic=0;$total3primeutr=0;$total5primeutr=0;
   $totalallexons=0;$totalsyncoding=0;$totalnonsyncoding=0;
   $totalSNPS=$SNPlist[0];
   $options=""; $type=4;
   if($fields{"v_type"}==1) $options.="&vt=1";
   if($fields{"g_name"}==1) $options.="&gn=1";
   if($fields{"f_class"}==1) $options.="&fc=1";
   if($fields{"cpg"}==1) $options.="&cp=1";
   if($fields{"m_allele"}==1) $options.="&ma=1";
   if($fields{"m_genotype"}==1) $options.="&mg=1";
   $totalIntronic=getQuerySNPsCount($con,$chr,$startpos,$endpos,$selectedSources,0,$gene_id,$transcript_id,$sg);
   $totalallexons=getQuerySNPsCount($con,$chr,$startpos,$endpos,$selectedSources,10,$gene_id,$transcript_id,$sg);
   $total3primeutr=getQuerySNPsCount($con,$chr,$startpos,$endpos,$selectedSources,34,$gene_id,$transcript_id,$sg);
   $totalsyncoding=getQuerySNPsCount($con,$chr,$startpos,$endpos,$selectedSources,1,$gene_id,$transcript_id,$sg);
   $totalnonsyncoding=getQuerySNPsCount($con,$chr,$startpos,$endpos,$selectedSources,2,$gene_id,$transcript_id,$sg);
   $filters.="<ul><li style='margin-top:1.5em;'><a name='snp_loc=12' href='#' onclick=\"updatePageResult($type,this);\" ><b>SNPs Match: $totalSNPS</b></a></li>";
   $filters.="<li class='sr' style='margin:auto;padding:0.em;margin-top:1em;'><h2>Filter Results By:</h2></li>";
   if($totalIntronic>0)
      $filters.="<li><a onclick=\"updatePageResult($type,this);\" name='snp_loc=0' href='#'>Intronic ($totalIntronic)</a></li>";
  if($totalallexons>0)
      $filters.="<li><a onclick=\"updatePageResult($type,this);\" name='snp_loc=10' href='#'>SNPs on Exons ($totalallexons)</a></li>";
   if($total3primeutr>0)
      $filters.="<li><a onclick=\"updatePageResult($type,this);\" name='snp_loc=34' href='#'>UTR ($total3primeutr)</a></li>";
   if($totalsyncoding>0)
      $filters.="<li><a onclick=\"updatePageResult($type,this);\" name='snp_loc=1' href='#'>Coding Synonymous ($totalsyncoding)</a></li>";
   if($totalnonsyncoding>0)
      $filters.="<li><a onclick=\"updatePageResult($type,this);\" name='snp_loc=2' href='#'>Coding Non Synnonymous ($totalnonsyncoding)</a></li>";
   $filters .="</ul></td>";
   $filters.="  <td style='width:12em;background:#efffff;border-right:solid 1px  #B2B4B6;'>
                       <div class='sr'><h2>Imputed Genotype Legend</h2>
                         <table ><tr><th></th><th>A</th><th>T</th><th>C</th><th>G</th></tr>
                                <tr><td>High Conf</td><td class='ah'></td><td class='th'></td>
                                            <td class='ch'></td><td class='gh'></td></tr>
                                 <tr><td>Med. Conf</td><td class='am'></td><td class='tm'></td>
                                            <td class='cm'></td><td class='gm'></td></tr>
                                 <tr><td>Low Conf</td><td class='al'></td><td class='tl'></td>
                                            <td class='cl'></td><td class='gl'></td></tr>
                       </table></div></td>";
   $filters.=" <td style='background:#efffff;'><div style='width:100%; margin:auto;background:#333;padding-left:0.5em;padding-top:0.5em;
                   padding-bottom:0.5em;
                   font-weight:bold;color:#eee;'>Result Columns(Click to toggle)</div>
                   <div style='width:23em; height:10em;overflow:scroll;'> 
                     <table id='hfield' border='1' style='width:99%;margin:auto;'>";
   $filters.=getColumns($strainsList,$fields);
   $filters.="</table></div></td></tr></table>";
   if($totalSNPS<=0){ // no SNPs found within the specified genomic region
       print"$querySum<br/><hr/><div class='sr'> Zero result ,our database does not 
          have SNPs for your current filter settings. Please change your filters setting and try again</div>";  
   }
   else{
     $header.="$querySum $filters";
     print"$header";$strain_names;$currenttotalSNP=$SNPlist[0];
     print"<div style=' font-size:0.85em; margin:auto; height:2.5em;padding:1em; color:navy;'>
                      <p style='background:#fff;width:90%;margin:auto;text-align:left;'><strong>Hint</strong> :
                      You can rearrange columns order by drag and drop.You can also hide any
                      column by clicking the coresponding checkbox above
                     </p> </div>";
     $type=0; 
     $result.=displayResulttable($currenttotalSNP,$strainsList,$SNPlist,$rowcount,$page,$datatype,$selectedStrains,
                        $is_imputed,$selectedSources,$strainIndex,$snp_loc,$fields,$querystring,$filterlink,$type); 
   print"<div id='tresult' style='background:#eee;'><a href='$filterlink' id='hq'>hd</a>$result</div>";             
   }
  
 }

/********************************************************
 this function returns the list of columns in snptable
*********************************************************/
function getColumns($strainsList,$fields){
  $columnno=1;$tableid="snptable";$type=1;
  $columns="<tr><td><input type='checkbox' class='check' name='check' id='chrom-1'
                  value='t$columnno' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                 <label>chrom</label>";
  ++$columnno;
  $columns.="<input type='checkbox' class='check' name='check' id='position-1'
                 value='t$columnno' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                 <label>Position</label></td></tr>";
  ++$columnno;
  // now get the strains list
  if(is_array($strainsList)){
       $data=""; $i=0;
        while(list($strain_id,$strain_name)=each($strainsList)){
           // if($i%3==0)$data.="</br>";
            $columns.="<tr><td><input type='checkbox' class='check' name='check' 
                    value='t$columnno' id='st=$strain_id' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                    <label>$strain_name</label><t/d></tr>";
            
            ++$i;++$columnno;
         }
    //$columns.="$data";
   }
   $columns.="<tr><td><input type='checkbox' value='t$columnno' class='check' 
                name='check' id='snpid-1' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                <label>SNPID</label></td></tr>";
   ++$columnno;
   if($fields{"v_type"}==1){
         $field_name="Variation Type";
         $columns.="<tr><td><input type='checkbox' value='t$columnno' name='check' 
                       class='check' id='vt=1' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                       <label>$field_name</label></td></tr>"; ++$columnno;
  }
  if($fields{"g_name"}==1){
         $field_name="MGI Gene";
         $columns.="<tr><td><input type='checkbox' value='t$columnno' class='check' id='gn=1'
                        name='check' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                       <label>$field_name</label></td></tr>"; ++$columnno;
  }
  if($fields{"f_class"}==1){
         $field_name="Function class";
         $columns.="<tr><td><input type='checkbox' value='t$columnno' name='check' id='fc=1'
                       class='check' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                       <label>$field_name</label></td></tr>"; ++$columnno;
  }
  if($fields{"cpg"}==1){
     $field_name="CpG";
     $columns.="<tr><td><input type='checkbox' value='t$columnno' name='check' id='cp=1'
                       class='check' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                       <label>$field_name</label></td></tr>"; ++$columnno;
   }
   $field_name="Alleles";
   $columns.="<tr><td><input type='checkbox' value='t$columnno' name='check' id='allele-1'
                       class='check' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                       <label>$field_name</label></td></tr>"; ++$columnno;
   if($fields{"m_allele"}==1){
      $field_name="Minor Allele %";
      $columns.="<tr><td><input type='checkbox' value='t$columnno' name='check' id='ma=1'
                       class='check' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                       <label>$field_name</label></td></tr>"; ++$columnno;  
   }
   if($fields{"m_genotype"}==1){
        $field_name="Missing genotype %";
        $columns.="<tr><td><input type='checkbox' value='t$columnno' name='check' id='mg=1'
                       class='check' onclick=\"toggleColumn(this,'$tableid',$type);\" checked/>
                       <label>$field_name</label></td></tr>";
   }
   return $columns;
}
/***********************************************************************************
* This function displays the detail page of a given SNP ***************************
************************************************************************************/
function displaySnpDetailView($snpid,$selectedStrains){
  global  $ASSEMBLY_BUILD;global $NCBI_SNP_PAGE; global $GENOTYPEMAP;
  global $MGI_GENEDETAIL;global $ENTREZ_GENE;global  $ENSEMBL_CURRENT_PAGE;global $ENSEMBL_CURRENT_PEPTIDE;
  // get the accession id,chr,pos,
  $snpdata=getsnpdata($link,$snpid);
  if(!empty($snpdata)){
     $chromosome_id=$snpdata{"chr"};$bp_position=$snpdata{"pos"};
      $black6_allele=$snpdata{"b6allele"};$snp_allele=$snpdata{"snpallele"};
      $is_CpG_site=$snpdata{"iscpg"};$is_intergenic=$snpdata{"istergenic"};
      $cpg="no"; if($is_CpG_site)$cpg="yes"; $is_conflict=$snpdata{"is_conflict"};
      $mutation="Transversion";$mutation_type=$snpdata{"mutation_type"};
      if($mutation_type==1)$mutation="Transition";$impgenotypedata="";$genotypeData="";
      $chromosome_name=$chromosome_id;
      if($chromosome_id==20) $chromosome_name="X" ;
      if($chromosome_id==21) $chromosome_name="Y" ;
      if($chromosome_id==22) $chromosome_name="Mt" ;
      $chromosome_name="Chr$chromosome_name";
      $conflict_list="No problems were identified for this SNP";
      if($is_conflict==1){         
         $conflict_list="<b>This SNP has the following issue(s):</b><br/><ul>";
         // then get the list of conflict of this SNP
         $conflict_list.= getSNPConflictList($link,$snpid);
         $conflict_list .="</ul>";
       }
      
      $i=0; $all=1; $geneid=""; $intergenic=0;$gene_start=0;$gene_end=0;$selectedSources="";
      $snpaccession_ids=getSNPaccession($snpid,$selectedSources,$all);
      $accession_ids= split(",",$snpaccession_ids);$accession_id="";
      while($i<count($accession_ids)){
         $pos=preg_match('/rs\d+/',$accession_ids[$i]);
         if($pos!=0){
            if($i==0)$accession_id="<a href='$NCBI_SNP_PAGE$accession_ids[$i]'>$accession_ids[$i]</a>";
            else
               $accession_id.="||<a href='$NCBI_SNP_PAGE$accession_ids[$i]'>$accession_ids[$i]</a>";
         }
         else{
            if($i==0)$accession_id=$accession_ids[$i];
            else
                $accession_id .=" || $accession_ids[$i]";
         }
        
         ++$i;
      }
      $source_list= getSNPsources($con, $snpid);
      $sourcelist= $source_list{"list"};$is_imputed= $source_list{"type"};
      $selectedSources= $source_list{"ids"};$impgenoSum="";
      $sourcelists= split(":",$selectedSources);
      ///////////////////////
      // check id all the strains have the same genotype
      // get minor allele frequency :
       $impgenoSum=array();$genoSumarry=array();
       if(count($selectedSources)==1)
               $genoSumarry= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,$is_imputed);
       else{
               $genoSumarry= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,0);
               if($is_imputed==1)$impgenoSum= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,$is_imputed);   
        }
       $minorallele=""; $genomap="";
       $genomap= parseMinorAllele($genoSumarry,$impgenoSum,$snp_allele,$black6_allele);
       $minallelef=computeMAF($genomap,$black6_allele,$snp_allele,$sourcelists);
       $minal  =  $minallelef{"minal"};
       $minmaf =  $minallelef{"minmaf"};$minorallele="";
       $minalleles=split(":",$minal);
       $minmafs =split(":",$minmaf);$t=0;
       $allelecount=count($minalleles);
       while($t<$allelecount){
         
             if($t==0){$minorallele ="($minalleles[$t]) $minmafs[$t]";}
             else{$minorallele .=":($minalleles[$t]) $minmafs[$t]";}
             ++$t;
       }
       $nmap   =  $minallelef{"nmap"};
       $totalStrain = $minallelef{"totalStrain"};
       //////////////
      
       //$minorallele = "$minal - $minmaf "; 
      print "<table border='1' style='width:98%;margin:auto;'><tr><td class='sectionName'><em>SNP Summary</em></td><td colspan='3'><table class='dtable' border='1'><tr>
                        <td class='sectionName'>Chromosome</td><td class='sectionOptions'>$chromosome_name</td>
                        <td class='sectionName'><em> Source(s)</em></td><td class='sectionOptions'>$sourcelist</td>
                   </tr>";
       print"      <tr><td class='sectionName'><em>Bp_Position</em></td><td class='sectionOptions'>$bp_position ($ASSEMBLY_BUILD)</td>
                      <td class='sectionName'><em>Submitter_ID(s)</em></td><td class='sectionOptions'>$accession_id</td>
                    </tr>";
       print"       <tr><td class='sectionName'><em>Alleles</em></td><td class='sectionOptions'><font color='red'> $black6_allele</font>
                             /$snp_allele [C57BL/6J allele=$black6_allele]</td>
                   <td class='sectionName' width='150'><em>Strain Count</em></td><td class='sectionOptions' width='170'>$totalStrain</td>
                     </tr>";
       print "<tr><td class='sectionName'><em>Variation Type </em></td><td class='sectionOptions'>$mutation </td>
                  <td class='sectionName'><em>Minor Allele Frequency %</em></td><td class='sectionOptions'>$minorallele</td>              
              </tr>";
       print"<tr><td class='sectionName'><em>CpG island</em></td><td class='sectionOptions'>$cpg</td>
                <td class='sectionName'><em>Missing allele frequency (by source) %</em></td>
                <td class='sectionOptions'>$nmap</td>
             </tr>";
      print "</table></td></tr>";
      //////////////////////////////////
      print"<tr><td class='sectionName' style='background:#ddd;'><em>SNP Check</em></td>
              <td colspan='3'class='sectionOptions' style='background:#efefef;'> $conflict_list</td></tr>";
      // now get the sources of this snp
      /**/
       //now display strains header
      $strainsList = getStrainList($con,$snpid,$is_imputed,$selectedStrains);
      $j=0;$strainheader="";$strainIndex; $impgenotypedata=""; $genotypeData="";
      $geno_allele=""; $imp_allele=""; $selectedSources="";
     
      if($is_imputed==1)$impgenotypedata=getGenotypedata($snpid,$selectedStrains,$selectedSources,$is_imputed);
      $genotypeData =getGenotypedata($snpid,$selectedStrains,$selectedSources,0);
     if(is_array($strainsList)){
         $genoheader="<tr>"; $thisgeno="<tr>";
      while(list($strain_id,$strain_name)=each($strainsList)){
           $len=strlen($strain_name);
            $verticalName="<ul>";
            for($p=0;$p<=$len-1;++$p){
                if(!empty($strain_name[$p]))$verticalName.="<li>$strain_name[$p]</li>";
            }
            $verticalName.="</ul>";
            $genoheader.="<td valign=bottom class='strainlabel' width='4px'>$verticalName</td>";
            $conf="H";$char="N";$impchar="";$genochar=""; $html=1;
            $strainFoundinImputed=0; //sets flag for imputed genotype exists for this strain
            $strainFoundinGeno=0; //sets flag for  strain genotype exists for selected source(s)
            if(!empty($impgenotypedata{"$strain_id"})){ //get the imputed genotype
                $thisStraingeno=split(":",$impgenotypedata{"$strain_id"}); 
                list($has_conflict,$geno)=split(":",getThisGenotype($thisStraingeno,$html));
                $thisgeno.=$geno;
               ++$l;
             }
             else if(!empty($genotypeData{"$strain_id"})){
                   $thisStraingeno=split(":",$genotypeData{"$strain_id"}); 
                   list($has_conflict,$geno)=split(":",getThisGenotype($thisStraingeno,$html));
                   $thisgeno.=$geno;
                   ++$l;
             }
             else{  // this strain does not have a genotype data / strain not provided by this source
                 $conf="H";$char="N";
                 if($strain_id==7){
                    $thisgeno.= $GENOTYPEMAP{"$black6_allele"}{"$conf"};
                    $thisgeno .="$black6_allele</td>";++$l;
                  }
                 else{$thisgeno .="<td class='n'>N</td>";++$l;}      
             }
          ++$k;
       }
      $genoheader .="</tr>"; $thisgeno.="</tr>";
      }
      // now get the genotype for each strain,
      // get the genotype data from snp_genotype
      // if imputed, get the genotype data from imputed
      $snplocation="";
      if($is_intergenic==0){
         $snplocdata=getSNPGenelocation($snpid);
         $snplocation="<table border='1'  style='background:#fff;'>";
         while(list($gene_id,$transcriptdata)=each($snplocdata)){
               $mgigene=getMgiGene($gene_id,0);
               $ensembl_symbol=getEnsGene($gene_id);
               $gene="<tr class='sectionOptions' style='background:#fff;'><td><b> Ensembl Gene: $ensembl_symbol</b>";
               if(!empty($mgigene)){
                  $mgigene_symbol=$mgigene[1];$mgigene_id=$mgigene[0];
                  $mgigene="<b> MGI Symbol:</b> <a href='$MGI_GENEDETAIL/accession_report.cgi?id=$mgigene_id' 
                           target='_blank'> $mgigene_symbol</a>";
                  $gene.=" || $mgigene";
               }
               $entrezgene =getEntrezGene($gene_id);
              if(!empty($entrezgene)){
                  $entrezgeneid= $entrezgene[0];
                  $entrezgene="<b> Entrez Gene ID:</b> <a href='$ENTREZ_GENE$entrezgeneid' target='_blank'>$entrezgeneid</a>";
                  $gene.=" || $entrezgene";
               }
               $snplocation.="$gene</td></tr><tr><td>";
              $txtable="<table  class='innertranstable' border='1'>
                     <tr class='sectionOptions'><th>Transcript</th><th class='innertranstable'>Protein</th><th>CDS Start</th>
                            <th class='innertranstable'>CDS Stop</th><th>Type</th><th>T_start</th><th>T_end</th>";
              $txtable.="<th>Strand</th><th>Exon Count</th><th>SNP implication</th></tr>";
               while(list($lockey,$transcriptlist)=each($transcriptdata)){
                     $transcriptdatal=split(":",$transcriptlist);
                     foreach($transcriptdatal as $transc_id){
                        if($transc_id>0){$datum="";
                           $translation =getProteinData($transc_id);
                           $transc_data=getTranscriptInfo($transc_id);
                           $t_accession_id=$transc_data{"id"}; $tstart=$transc_data{"tstart"};
                           $location="Intronic";$loc="Intronic";
                           if($lockey==1){
                              $location=" <b> Coding Syn </b>";
                              //$location.= getSNPimplication($link,$snpid,$transc_id,$bp_position);
                            }
                           else if($lockey==2){
                              $location=" <b> Coding NonSyn </b>";
                              //$location.= getSNPimplication($link,$snpid,$transc_id,$bp_position);
                           }
                           else if($lockey==5){
                              $location="<b> Coding NonSyn - Stop Gained  </b>";
                              //$location.= getSNPimplication($link,$snpid,$transc_id,$bp_position);
                           }
                           else if($lockey==6){
                              $location=" <b> Coding NonSyn-- Stop Lost </b>";
                              //$location.= getSNPimplication($link,$snpid,$transc_id,$bp_position);
                           }
                           else if($lockey==7){
                              $location="<b> Coding NonSyn-- Start Gained</b>";
                             // $location.= getSNPimplication($link,$snpid,$transc_id,$bp_position);
                           }
                           else if($lockey==8){
                              $location=" <b> Coding NonSyn-- Start Lost </b>";
                             // $location.= getSNPimplication($link,$snpid,$transc_id,$bp_position);
                           }
                           else if($lockey==3)$location="3'UTR";
                           else if($lockey==4)$location="5'UTR";
                           else if($lockey==0){
                                $intron_rank=getIntronRank($transc_id,$bp_position);
                                $location="Intronic (rank: $intron_rank)";
                           }
                           $cds_start="N/A";$cds_end="N/A";$protein="N/A";
                           $geneData ="<td>$protein</td><td>$cds_start</td><td>$cds_end</td>";
                           if(!empty($translation)){
                              list($cds_start,$cds_end,$protein)=split(":",$translation);
                              //$geneData ="<td><a href='$ENSEMBL_CURRENT_PEPTIDE$t_accession_id' target='_blank'>
                              //$protein</a></td><td>$cds_start</td><td>$cds_end</td>";
                           }
                           $exon_count=getExonCount($transc_id);
                           $tend=$transc_data{"tend"}; $tstrand=$transc_data{"strand"};
                           if($lockey==1 ||$lockey==2||$lockey>4){
                              
                              $datum.="<tr class='sectionOptions' style='background:#fff;font-weight:bold;'><td>$t_accession_id  
                                       </td><td><a href='$ENSEMBL_CURRENT_PEPTIDE$t_accession_id' target='_blank'>$protein</a>
                                       </td><td>$cds_start </td><td>$cds_end </td><td>$T_type </td><td>$tstart
                                       </td><td>$tend</td><td>$tstrand</td><td>$exon_count</td><td>$location</td></tr>";
                              $datum.= "<tr class='sectionOptions'><td colspan='10'>";
                              $datum.=getSNPimplication($link,$snpid,$transc_id,$bp_position)."</td></tr>";
                           }
                           else{
                              $datum.="<tr class='sectionOptions' style='background:#fff;font-weight:bold;'>
                                           <td>$t_accession_id</td>$geneData
                                           <td>$T_type</td><td>$tstart</td><td>$tend</td>
                                            <td>$tstrand</td><td>$exon_count</td><td>$location</td></tr>";
                           }
                          //display this trans
                          $snplocation.="$txtable $datum</table>";
                         }
                     }
               }
            
               //$snplocation.="$datum";
         }
         print"<tr><td class='sectionName'><em>Gene Annotations </em></td>
               <td colspan='3'>$snplocation</table></td></tr>";  
      }
      else{  // display intergenic
          //SNP not located on transcript
             $snplocation="<em>SNP locate within Intergenic region </em>";//echo "Intergenic location";
             $data=getIntergenic($chromosome_id,$bp_position);
             print"<tr><td class='sectionName'><em>Gene Annotations </em></td>
                    <td  colspan='3' style='background:#fff;'>$snplocation $data</td></tr>";             
                       
      }  
     print"<tr><td class='sectionName'><em>Measured Variation</em></td><td colspan='3'class='sectionOptions' style='background:#fff;'>
            <table border='1'>$genoheader $thisgeno</table></td></tr></table>";    
    }  // end of snpid>0       
}

?>
