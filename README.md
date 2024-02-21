# dragon-fruit-genome-annotation

# 1-building the database using repeatmodeler
    #!/bin/bash
    #
    #SBATCH --job-name=build_database
    #SBATCH --output=build_database.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --time=25:00:00
    #SBATCH --mem=200G
    
    module load repeatmodeler/2.0.1
    
    BuildDatabase -name condor_database condor_assembly.bp.p_ctg.fasta

# 2-making model of database using repeatmodeler
    #!/bin/bash
    #
    #SBATCH --job-name=model_database
    #SBATCH --output=model_database.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --time=25:00:00
    #SBATCH --mem=200G
    
    module load repeatmodeler/2.0.1
    
    RepeatModeler -database condor_database -pa 8

# 3-masking the genome using repeatmasker
    #!/bin/bash
    #
    #SBATCH --job-name=repeatmasker
    #SBATCH --output=repeatmasker.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --time=25:00:00
    #SBATCH --mem=200G
    
    module load repeatmasker/4.1.1
    
    GENOME="condor_assembly.bp.p_ctg.fasta"
    
    RepeatMasker \
    -pa 8 \
    -gff \
    -nolow \
    -xsmall \
    -dir soft_mask \
    -lib condor_database-families.fa \
    ${GENOME}

# RNA-Seq data analysis for evidence
    #1-making reference using the assembled genome of condor
    #!/bin/bash
    #
    #SBATCH --job-name=reference_index
    #SBATCH --output=reference_index.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --time=10:00:00
    #SBATCH --mem=200G
    
    module load hisat2/2.1.0
    
    hisat2-build condor_assembly.bp.p_ctg.fasta condor_index


    #2-mapping the NCBI downloaded RNA-Seq reads adn sorted with samtools
    #!/bin/bash
    #
    #SBATCH --job-name=SRR19863829
    #SBATCH --output=SRR19863829.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --time=04:00:00
    #SBATCH --mem=50G
    
    module load hisat2/2.1.0 samtools/1.16.1
    
    hisat2 -p 32 --dta -x /ibex/project/c2141/dragon-fruit/hifireads_condor/RNA-seq_data/rnaseq_ref/condor_index \
    -1 SRR19863829_R1_trim.fastq.gz \
    -2 SRR19863829_R2_trim.fastq.gz | \
    samtools sort -@ 32 -o SRR19863829.bam

    #3-extract the alignment summary from mapping results
    cat *.summary | grep -E '*overall alignment rate'

# Different plant protein sequence data as evidence
    #download the protein sequence data of different genome from NCBI GENOME
    #1-combine the all protein sequence data present in .fasta
    cat ~/path/*.fa > all_protein.fa

    #2-remove the duplicate sequences from all_protein.fa
    #!/bin/bash
    #
    #SBATCH --job-name=remove_duplicates
    #SBATCH --output=remove_duplicates.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --time=10:00:00
    #SBATCH --mem=100G
    
    module load seqkit/0.10.1
    
    seqkit rmdup -s < all_protein.fa > ghb_db_cg_os_le.fa
        
# 4-Run BRAKER for gneome annotation using the evidence of RNA-Seq mapped data in .bam files and protein sequence data of different genomes

    #!/bin/bash
    #
    #SBATCH --job-name=braker
    #SBATCH --output=braker.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --gres=gpu:2
    #SBATCH --time=80:00:00
    #SBATCH --mem=200G
    
    module load singularity/3.6 braker/2.1.4
    
    
    singularity exec -B $PWD,/ibex,/sw --nv braker3.sif braker.pl --genome=/ibex/scratch/tariqr/genome_annotation/Braker_annotation/repeatmodel/soft_mask/condor_assembly.bp.p_ctg.fasta.masked \
    --prot_seq=/ibex/scratch/tariqr/genome_annotation/protein_sequences/ghb_db_cg_os_le.fa \
    --bam=/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/DRR128251.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/DRR128253.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/DRR128254.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR13775188.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR13775193.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR13775195.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR13775200.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR13775202.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR13775203.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216397.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216398.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216399.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216400.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216401.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216403.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216404.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR15216405.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR19863825.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR19863829.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR2924904.bam,\
    /ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR3203780.bam,/ibex/scratch/tariqr/genome_annotation/RNA_Seq_files/SRR8327215.bam \
    --species=condor --threads=32 --useexisting --GENEMARK_PATH=/ibex/scratch/tariqr/genome_annotation/Braker_annotation/gmes_linux_64_4/ \
    --PROTHINT_PATH=/ibex/scratch/tariqr/genome_annotation/Braker_annotation/gmes_linux_64_4/ProtHint/bin

    

