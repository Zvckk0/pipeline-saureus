#!/usr/bin/env bash
set -euo pipefail
conda activate staph_pipeline

SAMPLE="ERR15680387"
REF="GCF_000010465.1_Newman.fasta"
THREADS=8

echo ">>> [0/7] Téléchargement des données brutes (ENA)"
mkdir -p 00_Raw_Data && cd 00_Raw_Data
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR156/087/${SAMPLE}/${SAMPLE}_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR156/087/${SAMPLE}/${SAMPLE}_2.fastq.gz
cd ..

echo ">>> [1/7] Contrôle qualité et trimming (fastp)"
mkdir -p 01_QC
fastp -i 00_Raw_Data/${SAMPLE}_1.fastq.gz -I 00_Raw_Data/${SAMPLE}_2.fastq.gz \
    -o 01_QC/${SAMPLE}_fastp_R1.fastq.gz -O 01_QC/${SAMPLE}_fastp_R2.fastq.gz \
    --detect_adapter_for_pe --html 01_QC/fastp.html --json 01_QC/fastp.json \
    --thread ${THREADS}
fastqc 01_QC/${SAMPLE}_fastp_R*.fastq.gz -o 01_QC/
multiqc 01_QC/ -o 01_QC/multiqc/

echo ">>> [2/7] Assemblage de novo (Shovill)"
mkdir -p 02_Assembly
shovill --R1 01_QC/${SAMPLE}_fastp_R1.fastq.gz --R2 01_QC/${SAMPLE}_fastp_R2.fastq.gz \
    --outdir 02_Assembly/shovill_out --cpus ${THREADS} --ram 16 --minlen 500 --force
quast.py 02_Assembly/shovill_out/contigs.fa -o 02_Assembly/quast_out --threads ${THREADS}
checkm lineage_wf -x fa 02_Assembly/shovill_out/ 02_Assembly/checkm_out/ --threads ${THREADS}

echo ">>> [3/7] Annotation génomique (DFAST + KofamScan)"
mkdir -p 03_Annotation
dfast --genome 02_Assembly/shovill_out/contigs.fa --out 03_Annotation/dfast_out \
    --minimum_length 200 --cpu ${THREADS}
exec_annotation -f mapper -o 03_Annotation/resultats_kegg.txt \
    -p profiles/ -k ko_list --cpu ${THREADS} 03_Annotation/dfast_out/protein.faa

echo ">>> [4/7] Profilage AMR, virulence et plasmides"
mkdir -p 04_AMR
abricate --db card 02_Assembly/shovill_out/contigs.fa > 04_AMR/resultats_card.tsv
abricate --db vfdb 02_Assembly/shovill_out/contigs.fa > 04_AMR/resultats_vfdb.tsv
abricate --db plasmidfinder 02_Assembly/shovill_out/contigs.fa > 04_AMR/resultats_plasmides.tsv
amrfinder -n 02_Assembly/shovill_out/contigs.fa --organism Staphylococcus_aureus \
    -o 04_AMR/resultats_amrfinder.tsv --threads ${THREADS}
mykrobe predict --sample ${SAMPLE} --species staph \
    --seq 01_QC/${SAMPLE}_fastp_R1.fastq.gz 01_QC/${SAMPLE}_fastp_R2.fastq.gz \
    --format csv --out 04_AMR/${SAMPLE}_mykrobe.csv --threads ${THREADS}
python3 arbitre_amr.py

echo ">>> [5/7] Typage moléculaire (MLST + spa)"
mkdir -p 05_Typage
mlst --scheme saureus 02_Assembly/shovill_out/contigs.fa > 05_Typage/resultats_mlst.tsv
spaTyper --fasta 02_Assembly/shovill_out/contigs.fa -o 05_Typage/resultats_spa.tsv

echo ">>> [6/7] SNP calling et phylogénie (Snippy + FastTree)"
mkdir -p 06_SNP_Phylo
snippy --outdir 06_SNP_Phylo/snippy_out --ref ${REF} \
    --R1 01_QC/${SAMPLE}_fastp_R1.fastq.gz --R2 01_QC/${SAMPLE}_fastp_R2.fastq.gz \
    --cpus ${THREADS}
snippy-core --ref ${REF} -o 06_SNP_Phylo/core 06_SNP_Phylo/snippy_out
fasttree -nt -gtr 06_SNP_Phylo/core.full.aln > 06_SNP_Phylo/mon_arbre.nwk

echo ">>> [7/7] Pipeline terminé. Résultats disponibles dans 00_Raw_Data -> 06_SNP_Phylo/
