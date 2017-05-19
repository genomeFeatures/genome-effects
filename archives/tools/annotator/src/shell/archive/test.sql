dbhost=genome-mysql.cse.ucsc.edu
ucscusr=genome
proteindb=uniProt
ensdb=ensembldb.ensembl.org
ensusr=annonymous
TOP=`pwd`
taxon_file="organims_taxonomy.txt";

read -a taxons <<< $(mysql -h${dbhost} -u${ucscusr} -A ${proteindb} -se "SELECT concat(taxon,'=',replace(val,' ','::'))as rowdata FROM commonName");
#duplicate handler for stdout
exec 4<&1
exec 1> organims_taxonomy.txt
for row  in ${taxons[@]}; do
    taxon_id="0";taxon_name="";
    for token in `echo "$row" |grep -o -e "[^::]*"`;do
        if [ "$taxon_id" == "0" ];then taxon_id=${token}
        else taxon_name+=" ${token}"
        fi
    done
    #this is where I load 
    id=${taxon_id/=/$'\t'};echo "${id}${taxon_name}"
done
#restore stdout
exec 1<&4
