# Commands to build reference databases for Blast searches in the AHDB fusarium project

LS effectors need to be identified that show divergence between ff.spp or are
unique to a particular ff.spp. As such, blast databases of reference genomes need
to be constructed and searched against.


Key genomes on this list are:

* Fusarium oxysporum
* Fusarium oxysporum f. sp. cepae
* Fusarium oxysporum f. sp. narcissi
* Fusarium oxysporum f. sp. mathioli
* Fusarium oxysporum f. sp. lycopersici
* Fusarium proliferatum
* Fusarium redolens
* Fusarium avenaceum
* Fusarium culmorum
* Fusarium solani
* Fusarium gramminearum
* Fusarium poae
* Fusarium tricinctum
* Fusarium equiseti
* Fusarium lactis
* Fusarium cerealis
* Fusarium sambucinum
* Fusarium sacchari
* Fusarium acuatum
* Fusarium fujikuroi
* Fusarium coercileum
* Fusarium flocciferum

blast databases were made in the fusarium project and are found on the cluster
at ../analysis/AHDB_blast

## F.solani effectorP analysis

### Prepare queries JS-169 and Nacha2

Create a list of genes from FoM and FoN LS regions
```bash
for File in $(ls analysis/effectorP/F.solani/*/*_EffectorP_secreted_headers.txt | grep -v 'IMV_00293'); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  GeneFa=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended_renamed.cdna.fasta)
  cat $File | cut -f1 > tmp.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $GeneFa --headers tmp.txt \
    | sed "s/>/>${Strain}|/g"
done > analysis/AHDB_blast/F.solani_effP.fa
```

### Perform BLAST

Blast vs the reference genome databases:

```bash
for RefGenome in $(ls ../fusarium/analysis/metagenomics/reference_genomes/van_dam/renamed/*.fna); do
Prefix=$(basename $RefGenome .fna)
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 20s
printf "."
Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
done
printf "\n"
# Prefix=$(echo $RefGenome | cut -f5,6 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
OutDir=analysis/AHDB_blast/vs_ref_genomes/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/F.solani_effP.fa dna $RefGenome $OutDir
done
```

Blast against themselves

```bash
for RefGenome in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa); do
Prefix=$(basename $RefGenome _contigs_softmasked_repeatmasker_TPSI_appended.fa)
OutDir=analysis/AHDB_blast/vs_ref_genomes/$Prefix
# Prefix=$(echo $RefGenome | cut -f5,6 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/F.solani_effP.fa dna $RefGenome $OutDir
done
```

Blast vs sequenced genomes:

```bash
for RefGenome in $(ls ../fusarium/repeat_masked/F.*/*/ncbi_submission/*_contigs_unmasked.fa | grep -v 'old'); do
Prefix=$(echo $RefGenome | cut -f4,5 -d '/' --output-delimiter '_')
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/F.solani_effP.fa dna $RefGenome $OutDir
done
```

Blast vs work in progress genomes:

```bash
for RefGenome in $(ls ../fusarium/assembly/spades/F.*/*/filtered_contigs/contigs_min_500bp*.fasta | grep -w -e 'FON129' -e 'FON77' -e 'FON81' -e 'FON89' -e 'Straw465' -e 'F81' -e 'FOP1-EMR' -e 'R2' -e '15-074' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6' -e 'D2' -e 'PG8' -e 'L5' -e 'A1-2' -e 'HB6' | grep -v renamed | grep -v 'HB6'); do
Prefix=$(echo $RefGenome | cut -f5,6 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/F.solani_effP.fa dna $RefGenome $OutDir
done
for RefGenome in $(ls ../fusarium/assembly/spades_pacbio/F.*/*/filtered_contigs/contigs_min_500bp.fasta | grep -w -e '55'); do
Prefix=$(echo $RefGenome | cut -f5,6 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/AHDB_blast/vs_seq_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/AHDB_blast/F.solani_effP.fa dna $RefGenome $OutDir
done
```

## Summarise blast hists

```bash
CsvFiles=$(ls analysis/AHDB_blast/vs_*_genomes/*/*F.solani_effP.fa_hits.csv | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
Headers=$(echo $CsvFiles | sed 's&analysis/AHDB_blast/vs_ref_genomes/&&g' | sed 's&analysis/AHDB_blast/vs_seq_genomes/&&g' | sed -r "s&_F.solani_effP.fa_hits.csv&&g" | sed -r "s&/\S*&&g"  | sed 's&/van_dam&&g')
OutDir=analysis/AHDB_blast/vs_ref_genomes/extracted
mkdir -p $OutDir
Genomes=$(ls analysis/AHDB_blast/vs_*_genomes/*/*_genome.fa | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v 'GCA_000599445.1' | grep -e 'vs_ref_genomes' -e 'vs_seq_genomes' | grep -v 'F.oxysporum_/' | grep -v -e 'GCA_000599445.1' -e '77-13-4')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/blast_searches
$ProgDir/blast_parse_AHDB.py --blast_csv $CsvFiles --headers $Headers --genomes $Genomes --identity 0.50 --evalue 1e-30 --out_prefix $OutDir/F.solani_effP
```


A number of interesting genes (as present in JS-169) were identified which showed presence in F. solani, but absent in other Fusarium.

```
JS-169|g297.t1
JS-169|g1364.t1
JS-169|g1548.t1
JS-169|g1634.t1
JS-169|g2349.t1
JS-169|g2358.t1
JS-169|g3150.t1
JS-169|g3156.t1
JS-169|g3907.t1
JS-169|g4181.t1
JS-169|g4671.t1
JS-169|g5415.t1
JS-169|g5689.t1
JS-169|g7048.t1
JS-169|g8051.t1
JS-169|g8329.t1
JS-169|g8399.t1
JS-169|g8720.t1
JS-169|g9452.t1
JS-169|g9606.t1
JS-169|g10403.t1
JS-169|g10978.t1
JS-169|g11008.t1
JS-169|g11028.t1
JS-169|g11116.t1
```



These alignments were downloaded for further study.


```bash
ExtractDir=analysis/AHDB_blast/vs_ref_genomes/extracted
mkdir -p analysis/AHDB_blast/selected_loci
GeneList="JS-169|g297.t1 JS-169|g1364.t1 JS-169|g1548.t1 JS-169|g1634.t1 JS-169|g2349.t1 JS-169|g2358.t1 JS-169|g3150.t1 JS-169|g3156.t1 JS-169|g3907.t1 JS-169|g4181.t1 JS-169|g4671.t1 JS-169|g5415.t1 JS-169|g5689.t1 JS-169|g7048.t1 JS-169|g8051.t1 JS-169|g8329.t1 JS-169|g8399.t1 JS-169|g8720.t1 JS-169|g9452.t1 JS-169|g9606.t1 JS-169|g10403.t1 JS-169|g10978.t1 JS-169|g11008.t1 JS-169|g11028.t1 JS-169|g11116.t1"
x=$(echo $GeneList | sed 's/|/_/g')
for Name in $x; do
  ls $ExtractDir/F.solani_effP_${Name}_hits.fa
  cp $ExtractDir/F.solani_effP_${Name}_hits.fa analysis/AHDB_blast/selected_loci/.
done

InterPro=$(ls gene_pred/interproscan/F.solani/JS-169/JS-169_interproscan.tsv)
x=$(echo $GeneList | sed "s/JS-169|//g")
for Gene in $x; do
  cat $InterPro | grep $Gene
done >> analysis/AHDB_blast/selected_loci/interpro.tsv

```

```bash
mkdir -p analysis/AHDB_blast/all_targets
for File in $(ls analysis/AHDB_blast/vs_*_genomes/*/*_genome.fa); do
cp $File analysis/AHDB_blast/all_targets/.
done
tar -cz -f analysis/AHDB_blast/all_targets.tar.gz analysis/AHDB_blast/all_targets
rm -r analysis/AHDB_blast/all_targets
```
