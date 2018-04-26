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
```


```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '125' | grep -v 'old'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for RNADir in $(ls -d ../fusarium/qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -e 'prelim' -e 'PDA' | grep 'PDA'); do
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
        Prefix=$(echo $FileF | rev | cut -f1 -d '/' | rev | sed "s/_R.*_trim.fq.gz//g")
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
  BamFiles=$(ls $OutDir/*/*/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
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
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
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
AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```


Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
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
for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3 | grep '1177'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/ncbi_edits_repmask/*_softmasked_repeatmasker_TPSI_appended.fa)
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

Codingquary was noted to predict a gene that went beyond the end of contig 47 in
isolate 1177.

As such this gene was removed manually:

```bash
GffAppended=$(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '1177')
cp $GffAppended tmp.gff
cat tmp.gff | grep -v 'CUFF_8208_1_74' > $GffAppended
```
