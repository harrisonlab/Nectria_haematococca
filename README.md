# Nectria_haematococca
Commands for the analysis of Nectria haematococca (F. solani) ex. pea genomes.

#Data organisation

Data was copied from the raw_data repository to a local directory for assembly
and annotation.

```bash
    mkdir -p /home/groups/harrisonlab/project_files/N.haematococca
  	cd /home/groups/harrisonlab/project_files/N.haematococca
```

F. solani genomes were downloaded to the working directory


Two of the downloaded genomes contained metadata in the fasta header. This information
was remove using the following commands:

```bash
for Assembly in $(ls assembly/external_groups/F.solani/*/assembly/* | grep -e 'fna' | grep -v 'masked'); do
  Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
  Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
  echo "$Organism - $Strain"
  OutFile=$(echo $Assembly | sed 's/.fna/.fasta/g' | sed 's._//g')
  cat $Assembly | sed "s/ .*//" > $OutFile
done
```

# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
  for Assembly in $(ls assembly/external_groups/F.solani/*/assembly/* | grep -e 'fasta' | grep -v 'masked' | grep -e 'JS-169' -e 'IMV_00293'); do
    Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
    Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
    qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
    qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
  done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```

## Gene prediction


# Gene Prediction


Gene prediction followed three steps:
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.


#Gene prediction

Gene prediction was performed for Fusarium genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1

## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* RNAseq data from the Fusarium HAPI project was used
* qc of RNA seq data is detailed in the README file of this repository:


RNAseq data was downloaded from ncbi:
https://www.ncbi.nlm.nih.gov/sra/SRX1810229[accn]
https://www.ncbi.nlm.nih.gov/sra/SRX1810230[accn]
https://www.ncbi.nlm.nih.gov/sra/SRX1810231[accn]


```bash
# cd /home/groups/harrisonlab/project_files/N.haematococca
mkdir -p /data/scratch/armita/N.haematococca
cd /data/scratch/armita/N.haematococca
# FSOAMB
OutDir=raw_rna/F.solani/FMR4391/amphotericine_48h
mkdir -p $OutDir
fastq-dump --split-files -A SRR3609441 --gzip --outdir $OutDir
fastq-dump --split-files -A SRR3609442 --gzip --outdir $OutDir
fastq-dump --split-files -A SRR3609443 --gzip --outdir $OutDir
# FSOPSC
OutDir=raw_rna/F.solani/FMR4391/posaconazole_48h
mkdir -p $OutDir
fastq-dump --split-files -A SRR3609444 --gzip --outdir $OutDir
fastq-dump --split-files -A SRR3609445 --gzip --outdir $OutDir
fastq-dump --split-files -A SRR3609446 --gzip --outdir $OutDir
#FSONTC
OutDir=raw_rna/F.solani/FMR4391/dimethyl-sulfoxide_48h
mkdir -p $OutDir
fastq-dump --split-files -A SRR3609447 --gzip --outdir $OutDir
fastq-dump --split-files -A SRR3609448 --gzip --outdir $OutDir
fastq-dump --split-files -A SRR3609449 --gzip --outdir $OutDir
```



#### Aligning

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

<!-- ```bash
  for Reads in $(ls ../../../../../data/scratch/armita/N.haematococca/raw_rna/F.solani/FMR4391/*/*.fastq.gz); do
    Timepoint=$(echo $Reads| rev | cut -d '/' -f2 | rev)
    Strain=$(echo $Reads| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Reads | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain - $Timepoint"
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    echo $Reads
    OutDir=qc_rna/paired/$Strain/$Timepoint/interlevered
    qsub $ProgDir/rna_qc_fastq-mcf_unpaired.sh $Reads $IlluminaAdapters DNA $OutDir
  done
``` -->

```bash
for ReadsF in $(ls ../../../../../data/scratch/armita/N.haematococca/raw_rna/F.solani/FMR4391/*/*_1.fastq.gz); do    
Strain=$(echo $ReadsF| rev | cut -d '/' -f3 | rev)
Organism=$(echo $ReadsF | rev | cut -d '/' -f4 | rev)
TimePoint=$(echo $ReadsF | rev | cut -d '/' -f2 | rev)
echo "$Organism - $Strain - $TimePoint"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
ReadsR=$(echo "$ReadsF" | sed 's/_1.fastq.gz/_2.fastq.gz/g')
OutDir=../../../../../data/scratch/armita/N.haematococca/qc_rna/paired/$Organism/$Strain/$TimePoint
echo $ReadsF
echo $ReadsR
qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters RNA $OutDir
done
```
<!--
```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for File in $(ls -d ../N.haematococca/qc_rna/paired/*/*/interlevered/unpaired/*_trim.fq.gz); do
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        while [ $Jobs -gt 1 ]; do
        sleep 1m
        printf "."
        Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
        done
        printf "\n"
        Timepoint=$(echo $File | rev | cut -f4 -d '/' | rev)
        echo "$Timepoint"
        Prefix=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_trim.fq.gz//g')
        OutDir=../../../../../data/scratch/armita/N.haematococca/alignment/star/$Organism/$Strain/$Timepoint/$Prefix
        ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
        qsub $ProgDir/sub_star_unpaired.sh $Assembly $File $OutDir
      done
    done
  done
``` -->


```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for RNADir in $(ls -d ../../../../../data/scratch/armita/N.haematococca/qc_rna/paired/F.solani/FMR4391/*); do
FileNum=$(ls $RNADir/F/*_trim.fq.gz | wc -l)
for num in $(seq 1 $FileNum); do
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
done
printf "\n"
FileF=$(ls $RNADir/F/*_trim.fq.gz | head -n $num | tail -n1)
FileR=$(ls $RNADir/R/*_trim.fq.gz | head -n $num | tail -n1)
echo $FileF
echo $FileR
Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed "s/_1_trim.fq.gz//g")
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
echo "$Timepoint"
OutDir=../../../../../data/scratch/armita/N.haematococca/alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star.sh $Assembly $FileF $FileR $OutDir
done
done
done
```


Accepted hits .bam file were concatenated and indexed for use for gene model training:


```bash
qlogin -pe smp 8
cd /home/groups/harrisonlab/project_files/N.haematococca
for OutDir in $(ls -d ../../../../../data/scratch/armita/N.haematococca/alignment/star/*/*); do
  Strain=$(echo $OutDir | rev | cut -d '/' -f1 | rev)
  Organism=$(echo $OutDir | rev | cut -d '/' -f2 | rev)
  echo "$Organism - $Strain"
  # For all alignments
  BamFiles=$(ls $OutDir/*/*/*.sortedByCoord.out.bam | grep 'SRR'| tr -d '\n' | sed 's/.bam/.bam /g')
  mkdir -p $OutDir/concatenated
  samtools merge -@ 8 -f $OutDir/concatenated/concatenated.bam $BamFiles
done
logout
```


#### Braker prediction

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
AcceptedHits=$(ls ../../../../../data/scratch/armita/N.haematococca/alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

** Number of genes predicted:  **


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls ../../../../../data/scratch/armita/N.haematococca/alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```


Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa| grep -v 'IMV_00293'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain
    CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf)
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
  done
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

Note - Ensure that the "TPSI_appended.fa" assembly file is correct.



```bash
for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3 | grep -v 'IMV_00293'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
# -
# This section is edited
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
$ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
# -
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
done
```

In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


 ```bash
 for Gff in $(ls gene_pred/final_genes/*/*/final/final_genes_appended.gff3 | grep -w -e 'FON_63' -e 'Stocks4' | grep -w -e 'Stocks4'); do
 	Strain=$(echo $Gff | rev | cut -d '/' -f3 | rev)
 	Organism=$(echo $Gff | rev | cut -d '/' -f4 | rev)
 	echo "$Strain - $Organism"
 	cat $Gff | grep -w 'gene' | wc -l
 done
 ```

 ```
 Stocks4 - F.oxysporum_fsp_mathioli
 20361
 FON_63 - F.oxysporum_fsp_narcissi
 20545
 ```


 In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
  * no duplicate genes were identified


 ```bash
 for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3); do
 Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
 Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
 echo "$Organism - $Strain"
 FinalDir=gene_pred/final_genes/$Organism/$Strain/final
 mkdir -p $FinalDir
 ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
 # $ProgDir/remove_dup_features.py --inp_gff $GffAppended
 # $ProgDir/remove_dup_features.py --inp_gff $GffAppended | grep -A2 'Duplicate gene found' | tail -n1 | cut -f2 -d'=' > $FinalDir/filter_list.tmp
 GffFiltered=$FinalDir/filtered_duplicates.gff
 # cat $GffAppended | grep -v -w -f $FinalDir/filter_list.tmp > $GffFiltered
 # rm $FinalDir/filter_list.tmp
 $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
 GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
 LogFile=$FinalDir/final_genes_appended_renamed.log
 ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
 $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
 rm $GffFiltered

 Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
 $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed

 # The proteins fasta file contains * instead of Xs for stop codons, these should
 # be changed
 sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta
 done
 ```


## Assessing the Gene space in predicted transcriptomes:

 ```bash
 for Assembly in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.gene.fasta); do
 Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
 Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
 echo "$Organism - $Strain"
 ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
 # BuscoDB="Fungal"
 BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
 OutDir=gene_pred/busco/$Organism/$Strain/genes
 qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
 done
 ```

 ```bash
 	for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do  
 		echo $File;
 		cat $File | grep -e '(C)' -e 'Total';
 	done
 ```


#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteins $InterProRaw
done
```


## B) SwissProt


```bash
for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../uniprot/swissprot
SwissDbName=uniprot_sprot
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

## Small secreted proteins

Putative effectors identified within Augustus gene models using a number
of approaches:

* A) From Braker gene models - Signal peptide & small cystein rich protein


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
* SigP
* biopython
* TMHMM


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
 for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta); do
   SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
   ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   SplitDir=gene_pred/braker_split/$Organism/$Strain
   mkdir -p $SplitDir
   BaseName="$Organism""_$Strain"_braker_preds
   $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
   for File in $(ls $SplitDir/*_braker_preds_*); do
   Jobs=$(qstat | grep 'pred_sigP' | wc -l)
   while [ $Jobs -gt '20' ]; do
   sleep 10
   printf "."
   Jobs=$(qstat | grep 'pred_sigP' | wc -l)
   done
   printf "\n"
   echo $File
   qsub $ProgDir/pred_sigP.sh $File signalp-4.1
   done
 done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
 for SplitDir in $(ls -d gene_pred/braker_split/*/*); do
   Strain=$(echo $SplitDir | cut -d '/' -f4)
   Organism=$(echo $SplitDir | cut -d '/' -f3)
   InStringAA=''
   InStringNeg=''
   InStringTab=''
   InStringTxt=''
   SigpDir=braker_signalp-4.1
   echo "$Organism - $Strain"
   for GRP in $(ls -l $SplitDir/*_braker_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
     InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.aa";  
     InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp_neg.aa";  
     InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.tab";
     InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.txt";  
   done
   cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
   cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
   tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
   cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
 done
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
 for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta); do
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
   qsub $ProgDir/submit_TMHMM.sh $Proteome
 done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
# echo "$Organism - $Strain"
NonTmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $NonTmHeaders
SigP=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_aug_sp.aa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $NonTmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
# echo "Number of SigP proteins:"
TotalProts=$(cat $SigP | grep '>' | wc -l)
# echo "Number without transmembrane domains:"
SecProt=$(cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l)
# echo "Number of gene models:"
SecGene=$(cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l)
# A text file was also made containing headers of proteins testing +ve
PosFile=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
cat $PosFile | cut -f1 > $TmHeaders
printf "$Organism\t$Strain\t$TotalProts\t$SecProt\t$SecGene\n"
done
```

```
F.solani	IMV_00293	1682	1378	1376
F.solani	JS-169	1317	1055	1055
F.solani	Nacha2	1630	1333	1332
```



### C) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    BaseName="$Organism"_"$Strain"_EffectorP
    OutDir=analysis/effectorP/$Organism/$Strain
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
    qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

Note - this doesnt exclude proteins with TM domains or GPI anchors

```bash
  for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | grep -v 'Effector probability:' | cut -f1 > $Headers
    printf "EffectorP headers:\t"
    cat $Headers | wc -l
    Secretome=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    printf "Secreted effectorP headers:\t"
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/final_genes/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  done
```
```
F.solani - IMV_00293
EffectorP headers:	3152
Secreted effectorP headers:	304
F.solani - JS-169
EffectorP headers:	2626
Secreted effectorP headers:	219
F.solani - Nacha2
EffectorP headers:	2686
Secreted effectorP headers:	266
```

## SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash
for Secretome in $(ls gene_pred/braker_signalp-4.1/*/*/*_final_sp_no_trans_mem.aa); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
printf "Number of effectors predicted by EffectorP:\t"
EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
cat $EffectorP | wc -l
printf "Number of SSCPs predicted by both effectorP and this approach: \t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
echo ""
done
```

```
F.solani - IMV_00293
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	260
number of SSC-rich genes:	260
Number of effectors predicted by EffectorP:	304
Number of SSCPs predicted by both effectorP and this approach: 	186

F.solani - JS-169
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	166
number of SSC-rich genes:	166
Number of effectors predicted by EffectorP:	219
Number of SSCPs predicted by both effectorP and this approach: 	112

F.solani - Nacha2
% cysteine content threshold set to:	3
maximum length set to:	300
No. short-cysteine rich proteins in input fasta:	229
number of SSC-rich genes:	229
Number of effectors predicted by EffectorP:	266
Number of SSCPs predicted by both effectorP and this approach: 	168
```

### C) Identification of MIMP-flanking genes

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
GeneGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_appended_renamed.gff3)
OutDir=analysis/mimps/$Organism/$Strain
mkdir -p "$OutDir"
echo "$Organism - $Strain"
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
$ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
echo "The number of mimps identified:"
cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
echo "The following transcripts intersect mimps:"
MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
cat $MimpProtsTxt | wc -l
cat $MimpGenesTxt | wc -l
echo ""
done
```

Extraction of the B-tubulin region of IMV00293 revealed it to be a F. oxysporum
isolate and distinct from the other two isolates.
```
F.solani - IMV_00293
The number of mimps identified:
17
The following transcripts intersect mimps:
0
0

F.solani - JS-169
The number of mimps identified:
0
The following transcripts intersect mimps:
0
0

F.solani - Nacha2
The number of mimps identified:
0
The following transcripts intersect mimps:
0
0
```
<!--
Those genes that were predicted as secreted and within 2Kb of a MIMP
were identified:

```bash
for File in $(ls analysis/mimps/*/*/*_genes_in_2kb_mimp.txt | grep -w -e 'FON_63' -e 'Stocks4'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev | sed 's/_chromosomal//g')
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ProtsFile=$(echo $File | sed 's/genes/prots/g')
Secretome=$(ls gene_pred/final_genes_signalp-4.1/$Organism/$Strain/*_final_sp_no_trans_mem.aa)
OutFile=$(echo "$File" | sed 's/.gff/_secreted.gff/g')
SecretedHeaders=$(echo "$Secretome" | sed 's/.aa/_headers.txt/g')
cat $Secretome | grep '>' | tr -d '>' | sed 's/-p.//g' > $SecretedHeaders
SecretedMimps=$(echo "$File" | sed 's/.txt/_secreted_headers.txt/g')
cat $ProtsFile $SecretedHeaders | cut -f1 | sort | uniq -d > $SecretedMimps
cat $SecretedMimps | wc -l
cat $SecretedHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $File > tmp.txt
cat tmp.txt | wc -l
cat $SecretedHeaders | cut -f1 | cut -f1 -d '.' | sort | uniq | grep -f $File > tmp.txt
cat $SecretedHeaders | cut -f1 | sort | uniq | grep -f $File > tmp2.txt
MimpsGff=$(ls analysis/mimps/$Organism/$Strain/*_mimps.gff)
GenesIn2Kb=$(ls analysis/mimps/$Organism/$Strain/"$Strain"_genes_in_2kb_mimp.gff)
SecretedMimpsGff=$(echo $GenesIn2Kb | sed 's/.gff/_secreted.gff/g')

ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_gff_for_sigP_hits.pl tmp2.txt $GenesIn2Kb secreted_mimp ID > $SecretedMimpsGff
# cat $SecretedMimpsGff | grep -w 'mRNA' | wc -l
# cat $MimpsGff | grep -w -f $SecretedMimps > $SecretedMimpsGff
done
``` -->
