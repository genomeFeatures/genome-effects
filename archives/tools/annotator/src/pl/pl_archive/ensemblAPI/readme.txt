Some of the questions Joel suggested include:
What are we holding in our database that they aren't?  
What functionality do they have that we don't?  
What functionality do they have that we don't?  
How dependent is their code on also having a full ENSEMBL installation?  
How much effort would it take to (a) implement their package here, and (b) convert/load our data into their format?

lnh@edge:~/work/projects/graber_transcriptdb/src/pl/ensemblAPI> ./variationOnTranscript.pl 
species/group	homo_sapiens/core
database	homo_sapiens_core_62_37g
host:port	ensembldb.ensembl.org:5306

lnh@edge:~/work/projects/graber_transcriptdb/src/pl/ensemblAPI> ./variationOnTranscript.pl 
species/group	mus_musculus/core
database	mus_musculus_core_62_37o
host:port	ensembldb.ensembl.org:5306


The use of the registry ensures you will load the correct versions of the Ensembl databases for the software release it can find on a database instance. Using the registry object, you can then create any of number of database adaptors. Each of these adaptors is responsible for generating an object of one type. The Ensembl variation API uses a number of object types that relate to the data stored in the database. For example, in order to generate variation objects, you should first create a variation adaptor

#
#ensembl_website_68
#mus_musculus_variation_68_38
#
Issue:
You have to re install the modules for every new version of ensembl database.

Solution:
Create a script (perl) that checks for version updates for your organisms of intrest
The script could run one a week/month and alerts you if changes found.

Step:
get the content of ensembl browser and get the list of organisms and their verion and compare to
what the current API has if different then


