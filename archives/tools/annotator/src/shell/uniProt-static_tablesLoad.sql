dbhost=genome-mysql.cse.ucsc.edu
ucscusr=genome
proteindb=uniProt
ensdb=ensembldb.ensembl.org
ensusr=annonymous
TOP=`pwd`
#load taxonomy commonName table from ucsc
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
#load keyword table from ucsc - uniProt.keyword: keyword definition
read -a features <<< $(mysql -h${dbhost} -u${ucscusr} -A ${proteindb} -se "SELECT concat(id,'=',replace(val,' ','::')) as rowdata FROM keyword limit 100");
#duplicate handler for stdout
exec 4<&1
exec 1> uniProt.keyword.txt
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

#load extDb table from ucsc - uniProt.extDb: external database local id assignment
read -a features <<< $(mysql -h${dbhost} -u${ucscusr} -A ${proteindb} -se "SELECT concat(id,'=',replace(val,' ','::')) as rowdata FROM extDb limit 100");
#duplicate handler for stdout
exec 4<&1
exec 1> uniProt.extDb.txt
for row  in ${features[@]}; do
    db_id="0";db_name="";
    for token in `echo "$row" |grep -o -e "[^::]*"`;do
        if [ "$db_id" == "0" ];then db_id=${token}
        else db_name+=" ${token}"
        fi
    done
    #this is where I load 
    id=${db_id/=/$'\t'};echo "${id}${db_name}"
done
#restore stdout
exec 1<&4

echo "Program complete"
