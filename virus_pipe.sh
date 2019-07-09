#!/usr/bin/env bash 

#Working as of 2019-07-09 (YYYY-MM-DD)


#####################USES STANDARD IN######################
# virus_pipe.sh YOUR_FILE.fastq genome_to_map_to #
#Change cores:                                            #
CPU=8                                                     #
#Memory ~ CPU*10                                          #
memory=64G                                                #
###########################################################

echo "Working with: $1"
date +"%m/%d/%Y %H:%M:%S $HOSTNAME" 

test_fastq="$(basename $1 .fastq)"
test_fastqgz="$(basename $1 .fastq.gz)"
test_fasta="$(basename $1 .fasta)"

if [[ $1 = $test_fastqgz".fastq.gz" ]];
then 
	newfile="$(basename $1 .fastq.gz)"
	gzip -d $1
elif [[ $1 = $test_fastq".fastq" ]];
then 
	newfile="$(basename $1 .fastq)"
elif [[ $1 = $test_fasta."fasta" ]];
then
	newfile="$(basename $1 .fasta)"
else
	echo "File extension not recognized, need fasta or fastq"
	exit
fi

mkdir $newfile
cd $newfile

echo $1 | tr '-' '\t' | awk '{print ">" $4"_""SE" "\n" "AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"$4"ATCTCGTATGCCGTCTTCTGCTTG"}' > ${newfile}-adapters.fa
echo $1 | tr '-' '\t' | awk '{print ">" $4"_""SE_rc" "\t" "AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"$4"ATCTCGTATGCCGTCTTCTGCTTG"}' | reverse.py >> ${newfile}-adapters.fa

if [ $2 =  "genome_to_map_to" ]; 
then 
	Dictionary="location"
else
	echo "Dictionary error"
	exit
fi

#Trimm and remove adapaters
java -jar /local/cluster/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads $CPU -phred33 \
../${newfile}.fastq ${newfile}_trimmed.fastq \
ILLUMINACLIP:${newfile}-adapters.fa:2:30:10 \
LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5

#FastQC
fastqc ${newfile}_trimmed.fastq

#Gzip orginal fasta file
gzip ../${newfile}.fastq

#Trinity assembly
Trinity \
	--single ${newfile}_trimmed.fastq \
	--seqType fq --max_memory $memory --CPU $CPU --bflyCalculateCPU --output ${newfile}_trinity

#BWA mem
bwa mem -t $CPU $Dictionary ${newfile}_trinity/Trinity.fasta > ${newfile}.sam

#Tar and gzip Trinity folder then remove orginal
tar -zcvf ${newfile}_trinity.tar.gz ${newfile}_trinity --remove-files

#Pull out unmapped reads
samtools view -h -f 4 ${newfile}.sam > ${newfile}_unmapped.sam

#Convert unmapped reads to fastq
java -jar /local/cluster/picard-tools-2.0.1/dist/picard.jar SamToFastq I=${newfile}_unmapped.sam F=${newfile}_unmapped.fastq VALIDATION_STRINGENCY=LENIENT

#Convert fastq to fasta
sed -n '1~4s/^@/>/p;2~4p' ${newfile}_unmapped.fastq > ${newfile}_unmapped.fasta

#Gzip fastq
gzip ${newfile}_unmapped.fastq

#Cap3
/local/cluster/bin/cap3 ${newfile}_unmapped.fasta > ${newfile}_cap.out
cat ${newfile}_unmapped.fasta.cap.singlets ${newfile}_unmapped.fasta.cap.contigs > ${newfile}_unmapped_cap3.fasta

#Gzip fasta
gzip ${newfile}_unmapped.fasta

#Convert fasta to tab
fasta_formatter -t -i ${newfile}_unmapped_cap3.fasta -o ${newfile}_unmapped_cap3.tab

#Title of all unmapped seq to a file
awk '{print $1}' ${newfile}_unmapped_cap3.tab > ${newfile}_unmapped_cap3.title

#Blastn
/local/cluster/bin/blastn \
		-query ${newfile}_unmapped_cap3.fasta \
		-db nt	\
		-out ${newfile}-blastn-report_e_1.tab \
		-evalue 1 \
		-num_threads $CPU \
		-num_alignments 1 \
		-outfmt "6 qseqid staxids sscinames scomnames sblastnames sskingdoms stitle salltitles qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop sstrand qcovs qcovhsp"
	
#Collect blastn seq titles and seqs
awk '{print $1}' ${newfile}-blastn-report_e_1.tab > ${newfile}_unmapped_cap3_blastn_hits.title
awk 'NR == FNR{a[$1];next} $1 in a' ${newfile}-blastn-report_e_1.tab ${newfile}_unmapped_cap3.tab | awk '{print ">"$1"\n"$2}' > ${newfile}_unmapped_cap3_blastn_hits.fasta

	
#Collect not blastn seq and seqs
awk 'NR==FNR{a[$1]=1;next}!($1 in a)' ${newfile}_unmapped_cap3_blastn_hits.title ${newfile}_unmapped_cap3.title > ${newfile}-no-blastn.title
awk 'NR==FNR{a[$1]=1;next}!($1 in a)' ${newfile}_unmapped_cap3_blastn_hits.title ${newfile}_unmapped_cap3.tab | awk '{print ">"$1"\n"$2}' > ${newfile}-no-blastn.fasta

#Blastx
	/local/cluster/bin/blastx \
		-query ${newfile}-no-blastn.fasta \
		-db nr \
		-out ${newfile}-blastx-report_e_10.tab \
		-evalue 10 \
		-num_threads $CPU \
		-num_alignments 1 \
		-outfmt "6 qseqid staxids sscinames scomnames sblastnames sskingdoms stitle salltitles qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop sstrand qcovs qcovhsp"

#Collect blastx seq titles and seqs
awk '{print $1}' ${newfile}-blastx-report_e_10.tab > ${newfile}_unmapped_cap3_blastx_hits.title
awk 'NR == FNR{a[$1];next} $1 in a' ${newfile}_unmapped_cap3_blastx_hits.title ${newfile}_unmapped_cap3.tab | awk '{print ">"$1"\n"$2}' > ${newfile}_unmapped_cap3_blastx_hits.fasta
 	
#Collect not blastx seq titles and seqs
awk 'NR==FNR{a[$1]=1;next}!($1 in a)' ${newfile}_unmapped_cap3_blastx_hits.title ${newfile}-no-blastn.title > ${newfile}-no-blastn-no-blastx.title
awk 'NR==FNR{a[$1]=1;next}($1 in a)' ${newfile}-no-blastn-no-blastx.title ${newfile}_unmapped_cap3.tab | awk '{print ">"$1"\n"$2}' > ${newfile}-no-blastn-no-blastx.fasta
	
#Collect virus blastn titles to new file
awk '/virus/||/viral/||/viroid/' ${newfile}-blastn-report_e_1.tab | virus_database.py > ${newfile}_unmapped_cap3_blastn_virus_hits.tab
awk '{print $1}' ${newfile}_unmapped_cap3_blastn_virus_hits.tab > ${newfile}_unmapped_cap3_blastn_virus_hits.title
awk 'NR == FNR{a[$1];next} $1 in a' ${newfile}_unmapped_cap3_blastn_virus_hits.title ${newfile}_unmapped_cap3.tab | awk '{print ">"$1"\n"$2}' > ${newfile}_unmapped_cap3_blastn_virus_hits.fasta

#Collect virus blastx titles to new file
awk '/virus/||/viral/||/viroid/' ${newfile}-blastx-report_e_10.tab | virus_database.py > ${newfile}_unmapped_cap3_blastx_virus_hits.tab
awk '{print $1}' ${newfile}_unmapped_cap3_blastx_virus_hits.tab > ${newfile}_unmapped_cap3_blastx_virus_hits.title
awk 'NR == FNR{a[$1];next} $1 in a' ${newfile}_unmapped_cap3_blastx_virus_hits.title ${newfile}_unmapped_cap3.tab | awk '{print ">"$1"\n"$2}' > ${newfile}_unmapped_cap3_blastx_virus_hits.fasta

#Make reports
cat ${newfile}_unmapped_cap3_blastn_virus_hits.tab | tr ' ' '@' | awk '{print $4}' | tr '@' ' ' | sort | uniq -c | sed -e 's/ *//' -e 's/ /\t/' | tr '@' ' ' | awk -v lane=$newfile -v species=$2 '{print $0 "\t" "blastn" "\t" lane "\t" species}' > ${newfile}_virus_report.txt
cat ${newfile}_unmapped_cap3_blastx_virus_hits.tab | tr ' ' '@' | awk '{print $4}' | tr '@' ' ' | sort | uniq -c | sed -e 's/ *//' -e 's/ /\t/' | tr '@' ' ' | awk -v lane=$newfile -v species=$2 '{print $0 "\t" "blastx" "\t" lane "\t" species}' >> ${newfile}_virus_report.txt

#Add title to report
sed -i 'Number of hits \t Virus name \t blast \t lane' ${newfile}_virus_report.txt

#Add titles blastn
sed -i '1iqseqid \t staxids \t sscinames \t scomnames \t sblastnames \t sskingdoms \t stitle \t salltitles \t qgi \t qacc \t qaccver \t qlen \t sseqid \t sallseqid \t sgi \t sallgi \t sacc \t saccver \t sallacc \t slen \t qstart \t qend \t sstart \t send \t qseq \t sseq \t evalue \t bitscore \t score \t length \t pident \t nident \t mismatch \t positive \t gapopen \t gaps \t ppos \t frames \t qframe \t sframe \t btop \t sstrand \t qcovs \t qcovhsp' ${newfile}-blastn-report_e_1.tab
sed -i '1iqseqid \t staxids \t sscinames \t scomnames \t sblastnames \t sskingdoms \t stitle \t salltitles \t qgi \t qacc \t qaccver \t qlen \t sseqid \t sallseqid \t sgi \t sallgi \t sacc \t saccver \t sallacc \t slen \t qstart \t qend \t sstart \t send \t qseq \t sseq \t evalue \t bitscore \t score \t length \t pident \t nident \t mismatch \t positive \t gapopen \t gaps \t ppos \t frames \t qframe \t sframe \t btop \t sstrand \t qcovs \t qcovhsp' ${newfile}_unmapped_cap3_blastn_virus_hits.tab

#Add titles blastx
sed -i '1iqseqid \t staxids \t sscinames \t scomnames \t sblastnames \t sskingdoms \t stitle \t salltitles \t qgi \t qacc \t qaccver \t qlen \t sseqid \t sallseqid \t sgi \t sallgi \t sacc \t saccver \t sallacc \t slen \t qstart \t qend \t sstart \t send \t qseq \t sseq \t evalue \t bitscore \t score \t length \t pident \t nident \t mismatch \t positive \t gapopen \t gaps \t ppos \t frames \t qframe \t sframe \t btop \t sstrand \t qcovs \t qcovhsp' ${newfile}-blastx-report_e_10.tab
sed -i '1iqseqid \t staxids \t sscinames \t scomnames \t sblastnames \t sskingdoms \t stitle \t salltitles \t qgi \t qacc \t qaccver \t qlen \t sseqid \t sallseqid \t sgi \t sallgi \t sacc \t saccver \t sallacc \t slen \t qstart \t qend \t sstart \t send \t qseq \t sseq \t evalue \t bitscore \t score \t length \t pident \t nident \t mismatch \t positive \t gapopen \t gaps \t ppos \t frames \t qframe \t sframe \t btop \t sstrand \t qcovs \t qcovhsp' ${newfile}_unmapped_cap3_blastx_virus_hits.tab

mkdir run_files
mv *.title *.sam *.readcount *.cap.* run_files/

#mailer your_email_address@hello_world.com "$newfile Finished" ${newfile}_virus_report.txt