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

