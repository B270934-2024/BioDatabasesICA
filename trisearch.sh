#!/bin/bash
if [[ $1 == "" ]] || [[ ! -e $1 ]] && [[ $1 != "-h" ]]; then
	echo "No file provided. Please type './trisearch -h' for help."
	exit 
elif [[ $2 == "" ]] || [[ $3 == "" ]] && [[ $1 != "-h" ]];then
	echo "Username/Password empty. Please type ./trisearch -h for help."
	exit
elif [[ $1 == "-h" ]] || [[ $2 == "-h" ]] || [[ $3 == "-h" ]]; then
	echo "Please provide a file containing Ensembl Gene IDs as the first argument. Subsequently, provide your username and password for mySQL (for default account, use s2761220 and B77u0ygQ respectively). This will then download data from Ensembl, NCBI and Uniprot. It will provide an output table called finalTable.tsv, which is stored as tab separated values. It also retains the data found in txt files without headers, named 'ensemblfiles.tsv', 'ncbifiles.tsv', and 'uniprotfiles.tsv' respectively. These are uploaded to mySQL which can then be interacted with to produce vertically separated values to the screen (to aid readability), and a horizontal table to finalTable.tsv." 
exit
#Check to see what our initial input is, which should be entered after the file name and their login details. If '-h' is entered, this will provide help to the user.
elif [[ -e $1 ]]; then
rm -f *.tmp
rm -f ensemblfiles.tsv
rm -f ncbifiles.tsv
rm -f uniprotfiles.tsv
here=$(pwd)
count=0
seqSearch=$(cat ${1}) #Open the given file, and read its contents into an array.
for sequence in ${seqSearch[@]};do
	if [[ ${sequence} =~ ENS([A-Z]{3}G|G)[0-9]{11} ]];then #check if sequence is a valid ID, either 3 characters followed by a G, or ENSG for humans, followed by 11 numbers.
	curl -X GET https://rest.ensembl.org/lookup/id/${sequence} > ${sequence}.json.tmp #Using the Ensembl rest API, download JSON files with the given IDs
		echo -e "\nSearched information on ${sequence} in Ensembl.\n"
	else 
		echo -e "\n${sequence} is not a valid ID. Skipping...\n"
		count=count+1 # increase count of failed sequences by one.
		if [[ ${count} -eq $(wc -l <${1}) ]];then #if the file contains all failed sequences...
			echo -e "No files containing Ensembl Gene IDs were detected. Exiting."
			exit #inform the user, then exit the programme.
		fi
	fi
done

seqjson=$(ls *.json.tmp)

for sequence in ${seqjson[@]};do
	touch ${sequence}.txt.tmp
	cat "${sequence}" | grep 'species' | awk '{print $2}'> ${sequence}.txt.tmp
	cat "${sequence}" | grep 'display_name'| awk '{print $2}'>> ${sequence}.txt.tmp
	asnam=$(cat "${sequence}" | grep '^assembly_name' |awk '{print $2}')
	btype=$(cat "${sequence}" | grep '^biotype'| awk '{print $2}')
	obtype=$(cat "${sequence}" | grep '^object_type'| awk '{print $2}')
	desc=$(cat "${sequence}" | grep '^description'|awk '{print $2}')
	disnam=$(cat "${sequence}" | grep '^display_name'|awk '{print $2}')
	strt=$(cat "${sequence}" | grep '^start'|awk '{print $2}')
        sqen=$(cat "${sequence}" | grep '^end'|awk '{print $2}')
	seqre=$(cat "${sequence}" | grep '^seq_region_name'|awk '{print $2}')
	spec=$(cat "${sequence}" | grep '^species'|awk '{print $2}')
	strn=$(cat "${sequence}" | grep '^strand:'|awk '{print $2}')
	echo -e "${asnam}\t${btype}\t${obtype}\t${desc}\t${disnam}\t${strt}\t${sqen}\t${seqre}\t${spec}\t${strn}" >>ensemblfiles.tsv
	echo Converted ${sequence} into TSV format.
done 
# Converts the files from JSON to TSV, selecting only these specific pieces of information from each file, then outputting it to results.json.results. This also takes out the species and protein code name, putting them into a different file, a .json.txt file. 

seqjsontxt=$(ls *.json.tmp.txt.tmp)

for sequence in ${seqjsontxt[@]};do
	echo -e "Reading file" ${sequence}
	seqthead=$(head -n1 ${sequence}) #species
	seqttail=$(tail -n1 ${sequence}) #identifier
	#Firstly, take the species and identifier from the file, and put it into the esearch query to look into NCBI
	esearch -db gene -query "${seqthead}[Organism] AND ${seqttail}[Gene] AND (alive[prop])" | efetch -format tabular > ${seqttail}.tsv.tmp
	echo -e "Searched NCBI for ${seqttail} from ${seqthead}"
	wget -O ${seqttail}uniprot.tsv.search.tmp "https://rest.uniprot.org/uniprotkb/search?fields=accession,gene_names,organism_name,sequence,length&query=gene:${seqttail}+organism_name:${seqthead}+reviewed:true&format=tsv"
	#Then, search using the uniprot rest API to find only the specific information about them. This downloads TSV files.
	echo -e Searched UNIPROT for ${seqttail} and ${seqthead}
	sleep 2
	#Then, rest 2 seconds. This is to prevent NCBI or UNIPROT from preventing access by this programme.
done

tsvs=$(ls *.tsv.tmp)

for file in ${tsvs[@]};do
	awk 'BEGIN{FS="\t"}{if ($1 != "tax_id") {print $0}}' ${file} >> ncbifiles.tsv
done
#This puts all the TSV files into a new file, headless.txt. This removes the headers from each file.
uniprots=$(ls *.search.tmp)

for file in ${uniprots[@]};do
	awk 'BEGIN{FS="\t"}{if ($1 != "Entry") {print $0}}' ${file} >> uniprotfiles.tsv
done
#This copies the same information as above, making sure there's no headers in the final output file.

echo "DROP TABLE IF EXISTS genes_from_ensembl; DROP TABLE IF EXISTS genes_from_ncbi; DROP TABLE IF EXISTS genes_from_uniprot; CREATE table if not exists genes_from_ensembl(assembly_name varchar(50),biotype varchar (50),object_type varchar(50), description varchar(255),display_name varchar(50),start int,end int,seq_region_name int,species varchar(50),strand int);create table if not exists genes_from_ncbi(tax_id varchar(20),Org_name varchar(50),GeneID varchar(50),CurrentID varchar(50),Status varchar(50),Symbol varchar(50),Aliases varchar(50),description varchar(255),other_designations varchar(50),map_location int,chromosome int, genomic_nucleotide_accession_version varchar(50),start_position_on_the_genomic_accession int,end_position_on_the_genomic_accession int,orientation varchar(50),exon_count int);create table if not exists genes_from_uniprot(Entry varchar(20),Gene varchar(50) PRIMARY KEY, Organism varchar(50),Sequence varchar(255),Length int);load data local infile '${here}/ensemblfiles.tsv' into table genes_from_ensembl;load data local infile '${here}/ncbifiles.tsv' into table genes_from_ncbi;load data local infile '${here}/uniprotfiles.tsv' into table genes_from_uniprot;alter table genes_from_uniprot ADD COLUMN sequence_start varchar(10);UPDATE genes_from_uniprot SET sequence_start=substring(genes_from_uniprot.sequence,1,10);alter table genes_from_ncbi ADD COLUMN gene_length int; UPDATE genes_from_ncbi SET gene_length=genes_from_ncbi.end_position_on_the_genomic_accession - genes_from_ncbi.start_position_on_the_genomic_accession;" > command.sql
#This large command is echoed to command.sql, which creates and loads data from the previous files. This also appends a very short versionof the sequence to the file.

echo "select genes_from_ensembl.display_name,genes_from_ensembl.species, genes_from_ensembl.biotype,genes_from_uniprot.length AS length_of_protein, genes_from_ncbi.description as long_description, genes_from_ensembl.description as summary, genes_from_ncbi.start_position_on_the_genomic_accession as gene_start,genes_from_ncbi.end_position_on_the_genomic_accession as gene_end,genes_from_ncbi.gene_length, genes_from_ncbi.chromosome,genes_from_uniprot.sequence_start,genes_from_ncbi.orientation,genes_from_ncbi.exon_count from genes_from_ensembl LEFT JOIN genes_from_ncbi ON genes_from_ensembl.display_name = genes_from_ncbi.symbol LEFT JOIN genes_from_uniprot ON genes_from_ensembl.display_name LIKE substring_index((genes_from_uniprot.Gene),' ',1)\G;" > searchverticaloutput.sql

echo "select genes_from_ensembl.display_name,genes_from_ensembl.species, genes_from_ensembl.biotype,genes_from_uniprot.length AS length_of_protein, genes_from_ncbi.description as long_description, genes_from_ensembl.description as summary, genes_from_ncbi.start_position_on_the_genomic_accession as gene_start,genes_from_ncbi.end_position_on_the_genomic_accession as gene_end,genes_from_ncbi.gene_length, genes_from_ncbi.chromosome,genes_from_uniprot.sequence_start,genes_from_ncbi.orientation,genes_from_ncbi.exon_count from genes_from_ensembl LEFT JOIN genes_from_ncbi ON genes_from_ensembl.display_name = genes_from_ncbi.symbol LEFT JOIN genes_from_uniprot ON genes_from_ensembl.display_name LIKE substring_index((genes_from_uniprot.Gene),' ',1);" > searchhorizontaloutput.sql

#This queries all 3 tables we just created, using LEFT JOIN to append the data from all of them together. This is because sometimes we are not searching for proteins, which will not have entries on Uniprot. The first outputs it to the screen in a clean vertical format, the second allows it to be saved as tab separated values, which are easier to work with.

mysql s2761220 -u${2} -p${3} < command.sql
mysql s2761220 -u${2} -p${3} < searchverticaloutput.sql
mysql s2761220 -u${2} -p${3} < searchhorizontaloutput.sql > finalTable.tsv
if [[ $(wc -l < finalTable.tsv) -ne 0 ]]; then
	echo -e "Table created! This table displays: the name, the type, the sequence length, the chromosome, the first 10 characters of the sequence, the orientation, and the exon count. All searched data is located in 3 files: uniprotfiles.tsv, ncbifiles.tsv, ensemblfiles.tsv. They are also on your mySQL server, and in the text file finalTable.tsv"
else
	echo Error! Table not created.
fi
rm -f command.sql
rm -f *.tmp
fi
#This removes the temporary files we create, to clean up the folder at the end.
