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

