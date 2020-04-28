def init_vars():

    global EXPERIMENT_SETTINGS
    EXPERIMENT_SETTINGS = {                             # global dictionary of experiment settings
        "Experiment name": "",                          # name of experiment
        "Experiment directory": "",                     # experiment directory
        "Initial condition name": "",                   # name of initial condition (day 0)
        "Final condition name": "",                     # name of final condition (day X)
        "Library reference file": "",                   # path to library reference file
        "FDR": "",                                      # false discovery rate percentage as decimal
        "CPU cores": ""                                 # number of CPU cores to use
    }

    global LIBRARY_READ_PATH, LIST_INITIAL_PATHS, LIST_FINAL_PATHS, LIST_CAS9_NEG_I_PATHS, LIST_CAS9_NEG_F_PATHS

    LIBRARY_READ_PATH = ""                              # path to library reads file
    LIST_INITIAL_PATHS = ""                             # list of initial condition sample files
    LIST_FINAL_PATHS = ""                               # list of final condition sample files
    LIST_CAS9_NEG_I_PATHS = ""                          # list initial Cas9  neg files
    LIST_CAS9_NEG_F_PATHS = ""                          # list final Cas9  neg files

    global FILE_FLAGS
    FILE_FLAGS = {
        "Bowtie2 index name": "bowtie2-index",          # name of bowtie2 index folder
        "Bowtie2 index full path": "",                  # path to bowtie2 index
        "Trimmed read file": "-trimmed.fastq",          # file suffix for trimmed reads
        "Fasta library file": "-fasta.fa",              # file suffix for bam file
        "Aligned bam file": ".bam",                     # file suffix for bam files
        "Alignments dir": "trimmed-reads",              # folder for alignments
        "Trimmed dir": "trimmed-reads",                 # folder for trimmed-reads
        "Readcounts dir": "readcounts",                 # folder for read counts
        "Cas9 pos tag": "-Cas9-pos",                    # tag to add to Cas9 positive files
        "Cas9 neg tag": "-Cas9-neg"                     # tag to add to Cas9 negative files
    }

    global CAS9_TYPE_NAMING
    CAS9_TYPE_NAMING = [
        "_Cas9-",                                       # tag to add to Cas9 positive files - used in pre-processing
        "_Cas9neg-"                                     # tag to add to Cas9 negative files - used in pre-processing
    ]

    global COLUMN_NAMES
    COLUMN_NAMES = {
        "Library": "Library",                           # name for library column (without the .ES) in the self.ES df
        "Initial": "Initial.ES",                        # name for initial column (with the .ES) in the self.ES df
        "Final": "Final.ES",                            # name for final column (with the .ES) in the self.ES df
        "Normalized": ".norm",                          # suffix for normalized values in intermediate dfs
        "log2FC": ".log2FC",                            # suffix for log2FC values in intermediate dfs
        "average": ".avg",                              # suffix for average values in intermediate dfs
        "sgw": ".sgw",                                  # suffix for weighted sgRNA values in intermediate dfs
        "rank": ".rank",                                # suffix for gene rank values in intermediate dfs
        "ES": ".ES",                                    # suffix for enrichment scores in final self.ES df
        "ER": "EnrichmentRatio",                        # name of Enrichment Ratio column in self.ES df
        "Cas9 pos": "_Cas9-",                           # suffix for Cas9-positive samples - used in TRACS algorithm
        "Cas9 neg": "_Cas9neg-"                         # suffix for Cas9-negative samples - used in TRACS algorithm
    }

