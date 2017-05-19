<?php
/* *****************************************************
  Author: Lucie Ndzana Hutchins
  Date : July 2008
  
  Note: This utility stores global variables and path
****************************************************** */

$ITEM_TYPES{"2"}="Gene(s):gene_by_annotation:gene_name";
$ITEM_TYPES{"3"}="Transcript(s):transcript_by_annotation:transcript_name";
$ITEM_TYPES{"4"}="EST(s):est_align:qName";
$LIB_DATA_PATH ="/data/www/htdocs/estLib/estLib/data";
//$ORGANISMESTCOUNT ="organism_estscount_102208.xml"; /* stores ESTs count per organism */
$ORGANISMESTCOUNT="organism_est_count.xml";
$LIBESTSUMTEXT ="all_lib_estsum.txt"; /* stores ESTs statistical summary for each library */
$LIBESTSUM ="libraries_est_summary.xml"; /* stores ESTs statistical summary for each library */
$ORGANISMESTSUM="organismsEstsSum.txt"; /* stores ESTs statistical summary for each organism */
$ORGANISFILE= "organism_list.xml";    /* stores the list of all organisms found in Libdat table
                                         fields stored: date,
                                          organism_name,Submitter Count, Libraries Count, 
                                          Strain Count, Organ Count, Tissue Count,Cell_type Count
                                       */
$ORGANIS_SUBMITTER_FILE= "organism_submitter_list.xml"; /* stores the list of all organisms found in Libdat table
                                          with corresponding submitter list
                                         fields stored: date,
                                          organism_name,Submitter */
$ORGANIS_ORGAN_FILE= "organism_organ_list.xml"; /* stores the list of all organs of each organism found in Libdat table
                                          with corresponding organs list
                                         fields stored: date,
                                          organism_name,organ */
                                          
$ORGANIS_TISSUE_FILE= "organism_tissue_list.xml"; /* stores the list of all organisms found in Libdat table
                                          with corresponding tissues list
                                         fields stored: date,
                                          organism_name,organ */

$SUBMITTERLIST ="libsource_list.xml" ; /* stores the list of all the institutions that provide
                                           ESTs Libraries:
                                           Fields stored: date,
                                             Institution ID, Name,Organism Count,
                                             Libraries count, Organ Count, Tissue Count, Cell_type Count
                                        */
$LIBRARY_STATS ="library_summary.xml";  /* stores a summary statistics of each library:
                                           Fields stored: date,
                                            Library_id, Library_name, ests_minLen,
                                           ests_maxLen,ests_avgLen,ests_count
                                         */
$LIBRARY_DATA ="library_data.xml";     /* stores supporting info for each library:
                                           Fields stored: date,
                                            Library_id, Institution,Organism,Organ,Tissue,Cell_type
                                         */
?>
