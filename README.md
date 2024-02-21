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


