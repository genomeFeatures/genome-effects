dbhost=genome-mysql.cse.ucsc.edu
ucscusr=genome
proteindb=uniProt
ensdb=ensembldb.ensembl.org
ensusr=annonymous
TOP=`pwd`
#load taxonomy table from ucsc
read -a taxons <<< $(mysql -h${dbhost} -u${ucscusr} -A ${proteindb} -se "SELECT concat(taxon,'=',replace(val,' ','::')) as rowdata FROM commonName limit 1000");
#duplicate handler for stdout
exec 4<&1
exec 1>  uniProt.organims_taxonomy.txt
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
#load featureType table from ucsc
read -a features <<< $(mysql -h${dbhost} -u${ucscusr} -A ${proteindb} -se "SELECT concat(id,'=',replace(val,' ','::')) as rowdata FROM featureType limit 100");
#duplicate handler for stdout
exec 4<&1
exec 1> uniProt.featureType.txt
for row  in ${features[@]}; do
    feature_id="0";feature_name="";
    for token in `echo "$row" |grep -o -e "[^::]*"`;do
        if [ "$feature_id" == "0" ];then feature_id=${token}
        else feature_name+=" ${token}"
        fi
    done
    #this is where I load 
    id=${feature_id/=/$'\t'};echo "${id}${feature_name}"
done
#restore stdout
exec 1<&4

#load featureClass table from ucsc - uniProt.featureClass:  A class of protein feature
read -a features <<< $(mysql -h${dbhost} -u${ucscusr} -A ${proteindb} -se "SELECT concat(id,'=',replace(val,' ','::')) as rowdata FROM featureClass limit 100");
#duplicate handler for stdout
exec 4<&1
exec 1>  uniProt.featureClass.txt
for row  in ${features[@]}; do
    feature_id="0";feature_name="";
    for token in `echo "$row" |grep -o -e "[^::]*"`;do
        if [ "$feature_id" == "0" ];then feature_id=${token}
        else feature_name+=" ${token}"
        fi
    done
    #this is where I load 
    id=${feature_id/=/$'\t'};echo "${id}${feature_name}"
done
#restore stdout
exec 1<&4
echo "Program complete"
