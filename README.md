# diploid genome-annotation

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
    
    BuildDatabase -name genome_database genome.fasta

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
    
    RepeatModeler -database genome_database -pa 8

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
    
    GENOME="genome.fasta"
    
    RepeatMasker \
    -pa 8 \
    -gff \
    -nolow \
    -xsmall \
    -dir soft_mask \
    -lib genome_database-families.fa \
    ${GENOME}

# RNA-Seq data analysis for evidence
    #1-making reference using the assembled genome
    #!/bin/bash
    #
    #SBATCH --job-name=reference_index
    #SBATCH --output=reference_index.%j.out
    #SBATCH --partition=batch
    #SBATCH --cpus-per-task=32
    #SBATCH --time=10:00:00
    #SBATCH --mem=200G
    
    module load hisat2/2.1.0
    
    hisat2-build genome.fasta genome_index


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
    
    hisat2 -p 32 --dta -x ~/path/RNA-seq_data/rnaseq_ref/genome_index \
    -1 reads_R1_trim.fastq.gz \
    -2 reads_R2_trim.fastq.gz | \
    samtools sort -@ 32 -o reads.bam

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
    
    
    singularity exec -B $PWD,/ibex,/sw --nv braker3.sif braker.pl --genome=/ibex/scratch/tariqr/genome_annotation/Braker_annotation/repeatmodel/soft_mask/genome.fasta.masked \
    --prot_seq=/ibex/scratch/tariqr/genome_annotation/protein_sequences/ghb_db_cg_os_le.fa \
    --bam=~path/RNA_Seq_files/DRR128253.bam,\
    ~path/RNA_Seq_files/SRR13775188.bam,\
    ~path/RNA_Seq_files/SRR13775195.bam,\
    ~path/RNA_Seq_files/SRR13775202.bam,\
    --species=condor --threads=32 --useexisting --GENEMARK_PATH=/ibex/scratch/tariqr/genome_annotation/Braker_annotation/gmes_linux_64_4/ \
    --PROTHINT_PATH=/ibex/scratch/tariqr/genome_annotation/Braker_annotation/gmes_linux_64_4/ProtHint/bin

# 4-Run BRAKER3 for gneome annotation using the evidence of RNA-Seq mapped data in .bam files and protein sequence data of different genomes

	#!/bin/bash
	#
	#SBATCH --job-name=braker3
	#SBATCH --output=braker3.%j.out
	#SBATCH --partition=batch
	#SBATCH --cpus-per-task=32
	#SBATCH --time=300:00:00
	#SBATCH --mem=200G
	
	module load singularity/3.9.7 braker/3.0.3/singularity
	
	singularity exec -B $PWD,$HOME,/ibex,/sw /ibex/sw/rl9c/braker/3.0.3/singularity/braker3_latest.sif braker.pl \
	 --workingdir=~path/annotation/braker_annotation/braker1/ \
	 --genome=~path/genome.fasta.masked \
	 --prot_seq=/ibex/project/c2141/dragon-fruit/hifireads_condor/purge_haplotigs/annotation/protein_sequences/ghb_db_cg_os_le.fa \
	 --bam=~path/RNA-Seq_files/DRR128251.bam,\
	 ~path/RNA-Seq_files/DRR128253.bam,\
	 ~path/RNA-Seq_files/DRR128254.bam,\
	 ~path/RNA-Seq_files/SRR13775188.bam,\
	 --AUGUSTUS_CONFIG_PATH=~path/braker1/Augustus/config/ \
	 --crf --threads=32 --species=condor_purplepulp --busco_lineage=eudicots_odb10 --filterOutShort
	
