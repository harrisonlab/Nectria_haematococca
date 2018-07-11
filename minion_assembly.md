# Nectria_haematococca
Commands for the analysis of Nectria haematococca (F. solani) ex. pea genomes.

This document details the commands used to assemble and annotate minion sequence data.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/N.haematococca

The following is a summary of the work presented in this Readme.

The following processes were applied to Fusarium genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation


# 0. Building of directory structure

### Minion Data

```bash
	RawDatDir=/data/seq_data/minion/2018/20180222_AJ516/AJ516
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum_fsp_lactucae/AJ516

  RawDatDir=/data/seq_data/minion/2018/20180222_AJ520/AJ520
	# The scond run was output to a subfolder of the 20180222_AJ520 folder
  # RawDatDir=
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  mkdir -p $ProjectDir/raw_dna/minion/F.oxysporum_fsp_lactucae/AJ520
```


Data was basecalled again using Albacore 2.02 on the minion server:

```bash
screen -a
ssh nanopore@nanopore

# upgrade albacore
wget https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.2.7-cp34-cp34m-manylinux1_x86_64.whl
~/.local/bin/read_fast5_basecaller.py --version
pip3 install --user ont_albacore-2.2.7-cp34-cp34m-manylinux1_x86_64.whl --upgrade
~/.local/bin/read_fast5_basecaller.py --version

mkdir FoLatucae_26-04-18
cd FoLatucae_26-04-18

# Oxford nanopore AJ516 Run 1
Organism=F.oxysporum_fsp_lactucae
Strain=AJ516
Date=22-02-18
FlowCell="FLO-MIN106"
Kit="SQK-LSK108"
RawDatDir=/data/seq_data/minion/2018/20180222_AJ516/AJ516/GA10000/reads
OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date

mkdir -p $OutDir
cd $OutDir
~/.local/bin/read_fast5_basecaller.py \
  --flowcell $FlowCell \
  --kit $Kit \
  --input $RawDatDir \
  --recursive \
  --worker_threads 24 \
  --save_path albacore_v2.2.7 \
  --output_format fastq,fast5 \
  --reads_per_fastq_batch 4000

cat $OutDir/albacore_v2.2.7/workspace/pass/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz
tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
mkdir -p $FinalDir
mv $OutDir/${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
chmod +rw -R $FinalDir


	# OutDir=/data/scratch/nanopore_tmp_data/Fusarium/albacore_v2.2.7
  # mkdir -p $OutDir
  # chmod +w $OutDir
  # cp Alt_albacore_v2.10_barcode0*.fastq.gz $OutDir/.
  # chmod +rw $OutDir/Alt_albacore_v2.10_barcode0*.fastq.gz
  # scp Alt_albacore_v2.10_barcode*.fastq.gz armita@192.168.1.200:$OutDir/.
  # tar -cz -f Alt_albacore_v2.10_demultiplexed.tar.gz Alt_albacore_v2.10_demultiplexed
  # OutDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
  # mv Alt_albacore_v2.10_demultiplexed.tar.gz $OutDir/.
  # chmod +rw $OutDir/Alt_albacore_v2.10_demultiplexed.tar.gz

	# Oxford nanopore AJ516 Run 2
	Organism=F.oxysporum_fsp_lactucae
	Strain=AJ516
	Date=18-04-18
	FlowCell="FLO-MIN106"
	Kit="SQK-LSK108"
	RawDatDir=/data/seq_data/minion/2018/20180418_AJ16/AJ16/GA10000/reads
	OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date
	mkdir -p $OutDir

	mkdir -p $OutDir
	cd $OutDir
	~/.local/bin/read_fast5_basecaller.py \
	  --flowcell $FlowCell \
	  --kit $Kit \
	  --input $RawDatDir \
	  --recursive \
	  --worker_threads 24 \
	  --save_path albacore_v2.2.7 \
	  --output_format fastq,fast5 \
	  --reads_per_fastq_batch 4000

	cd ~/FoLatucae_26-04-18
	cat $OutDir/albacore_v2.2.7/workspace/pass/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz
	tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

	FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	mkdir -p $FinalDir
	mv ${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
	chmod +rw -R $FinalDir

	# Oxford nanopore AJ520 Run 1
	Organism=F.oxysporum_fsp_lactucae
	Strain=AJ520
	Date=22-02-18
	FlowCell="FLO-MIN106"
	Kit="SQK-RBK001"
	RawDatDir=/data/seq_data/minion/2018/20180222_AJ520/AJ520/GA20000/reads
	OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date
	mkdir -p $OutDir

	mkdir -p $OutDir
	cd $OutDir
	~/.local/bin/read_fast5_basecaller.py \
	  --flowcell $FlowCell \
	  --kit $Kit \
	  --input $RawDatDir \
	  --recursive \
	  --worker_threads 24 \
	  --save_path albacore_v2.2.7 \
	  --output_format fastq,fast5 \
	  --reads_per_fastq_batch 4000

	cat $OutDir/albacore_v2.2.7/workspace/pass/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz
	tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

	FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	mkdir -p $FinalDir
	mv ${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
	chmod +rw -R $FinalDir

	# Oxford nanopore AJ516 Run 2
	Organism=F.oxysporum_fsp_lactucae
	Strain=AJ520
	Date=18-04-18
	FlowCell="FLO-MIN106"
	Kit="SQK-LSK108"
	# The second run of AJ520 reads have not yet been copied over.
	# RawDatDir=/data/seq_data/minion/2018/20180222_AJ520/AJ520/GA30000/reads
	OutDir=~/FoLatucae_26-04-18/$Organism/$Strain/$Date
	mkdir -p $OutDir

	mkdir -p $OutDir
	cd $OutDir
	~/.local/bin/read_fast5_basecaller.py \
	  --flowcell $FlowCell \
	  --kit $Kit \
	  --input $RawDatDir \
	  --recursive \
	  --worker_threads 24 \
	  --save_path albacore_v2.2.7 \
	  --output_format fastq,fast5 \
	  --reads_per_fastq_batch 4000
		cat $OutDir/albacore_v2.2.7/workspace/pass/*.fastq | gzip -cf > ${Strain}_${Date}_albacore_v2.2.7.fastq.gz
		tar -cz -f ${Strain}_${Date}_albacore_v2.2.7.tar.gz $OutDir

		FinalDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
		mkdir -p $FinalDir
		mv ${Strain}_${Date}_albacore_v2.2.7.* $FinalDir/.
		chmod +rw -R $FinalDir
```

Sequence data was moved into the appropriate directories

```bash
	ProjectDir=/home/groups/harrisonlab/project_files/fusarium
	cd $ProjectDir

	Species="F.oxysporum_fsp_lactucae"
	Strain="AJ516"
	OutDir=$ProjectDir/raw_dna/minion/$Species/$Strain
	mkdir -p $OutDir
	RawDatDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	cp -s $RawDatDir/AJ516_18-04-18_albacore_v2.2.7.fastq.gz $OutDir/.
	cp -s $RawDatDir/AJ516_22-02-18_albacore_v2.2.7.fastq.gz $OutDir/.

	Species="F.oxysporum_fsp_lactucae"
	Strain="AJ520"
	OutDir=$ProjectDir/raw_dna/minion/$Species/$Strain
	mkdir -p $OutDir
	RawDatDir=/data/scratch/nanopore_tmp_data/Alternaria/albacore_v2.2.7
	cp $RawDatDir/AJ520_18-04-18_albacore_v2.2.7.fastq.gz $OutDir/.
	cp $RawDatDir/AJ520_22-02-18_albacore_v2.2.7.fastq.gz $OutDir/.
```

### MiSeq data

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180213_M04465_0067_000000000-BJ4DR/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/fusarium
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ516
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/AJ516_S1_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/AJ516_S1_L001_R2_001.fastq.gz .
  OutDir=$ProjectDir/raw_dna/paired/F.oxysporum_fsp_lactucae/AJ520
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/AJ520_S2_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/AJ520_S2_L001_R2_001.fastq.gz .
  cd $ProjectDir
```


## Assembly

### Removal of adapters

Splitting reads and trimming adapters using porechop
```bash
	for RawReads in $(ls  raw_dna/minion/*/*/*.fastq.gz | grep 'fsp_lactucae'); do
    Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
    Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
    echo "$Organism - $Strain"
  	OutDir=qc_dna/minion/$Organism/$Strain
  	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
  	qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
  done
```

#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz | grep 'lactucae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/* | grep 'lactucae'); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```

# Identify sequencing coverage

For Minion data:
```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz | grep 'lactucae'); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=60
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/* | grep 'Stocks4'); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```
MinION coverage was:
```
Stocks4	65.05
```

For Miseq data:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*q.gz | grep 'lactucae'); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
		qsub $ProgDir/run_fastqc.sh $RawData;
		GenomeSz=60
		OutDir=$(dirname $RawData)
		qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
	done
```

```bash
	for StrainDir in $(ls -d qc_dna/paired/*/* | grep 'lactucae'); do
		Strain=$(basename $StrainDir)
		printf "$Strain\t"
		for File in $(ls $StrainDir/*/*.txt); do
			echo $(basename $File);
			cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
		done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
	done
```

Miseq coverage was:
```
Stocks4	58.66
```


### Read correction using Canu

```bash
for TrimReads in $(ls qc_dna/minion/*/*/*q.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.6/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 58m $Strain $OutDir
done
```


### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.6/*/*/*.trimmedReads.fasta.gz | grep 'F.oxysporum_fsp_lactucae'); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```


Quast and busco were run to assess the effects of racon on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'F.oxysporum_fsp_lactucae'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
	echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
	OutDir=gene_pred/busco/$Organism/$Strain/assembly
	BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
	qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ520'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg | grep 'F.oxysporum_fsp_lactucae' | grep 'AJ516'); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon_2libs.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/*.fasta | grep 'round_10' | grep 'F.oxysporum_fsp_lactucae'); do
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta | grep 'F.oxysporum_fsp_lactucae'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/F.*/*/racon*/*.fasta | grep 'F.oxysporum_fsp_lactucae'); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/A*/*/assembly/*/short_summary_*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```