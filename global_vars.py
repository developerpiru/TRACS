def init_vars():

    global EXPERIMENT_SETTINGS
    EXPERIMENT_SETTINGS = {                             # global dictionary of experiment settings
        "Experiment name": "",                          # name of experiment
        "Experiment directory": "",                     # experiment directory
        "Library type": "",                             # type of CRISPR screen library used
        "Initial condition name": "",                   # name of initial condition (day 0)
        "Final condition name": "",                     # name of final condition (day X)
        "Cas9 Negative Control": "",                    # Cas9- negative control used; YES/NO
        "Library reference file": "",                   # path to library reference file
        "Library read file": ""                         # path to library read file
    }

    global LIBRARY_READ_PATH, LIST_INITIAL_PATHS, LIST_FINAL_PATHS, LIST_CAS9_NEG_I_PATHS, LIST_CAS9_NEG_F_PATHS

    LIBRARY_READ_PATH = ""                              # path to library reads file
    LIST_INITIAL_PATHS = ""                             # list of initial condition sample files
    LIST_FINAL_PATHS = ""                               # list of final condition sample files
    LIST_CAS9_NEG_I_PATHS = ""                          # list initial Cas9  neg files
    LIST_CAS9_NEG_F_PATHS = ""                          # list final Cas9  neg files

    # For testing:
    #LIST_INITIAL_PATHS = ["1", "2", "3"]
    #LIST_FINAL_PATHS = ["1", "2", "3"]
    #LIST_CAS9_NEG_I_PATHS = ["1", "2", "3"]
    #LIST_CAS9_NEG_F_PATHS = ["1", "2", "3"]

    global FILE_FLAGS
    FILE_FLAGS = {
        "Bowtie2 index name": "bowtie2-index",          # name of bowtie2 index folder
        "Bowtie2 index full path": "",                  # path to bowtie2 index
        "Trimmed read file": "-trimmed.fastq",          # file suffix for trimmed reads
        "Fasta library file": "-fasta.fa",              # file suffix for bam file
        "Aligned bam file": ".bam",                     # file suffix for bam files
        "Alignments dir": "alignments",                 # folder for alignments
        "Trimmed dir": "trimmed-reads",                 # folder for trimmed-reads
        "Readcounts dir": "readcounts",                 # folder for read counts
        "Cas9 pos tag": "-Cas9-pos",                    # tag to add to Cas9 positive files
        "Cas9 neg tag": "-Cas9-neg"                     # tag to add to Cas9 negative files
    }

    global REPLICATE_NAMING
    REPLICATE_NAMING = [
        "-1, -2, -3",
        "_1, _2, _3",
        "-A, -B, -C",
        "_A, _B, _C"
    ]

    global CAS9_TYPE_NAMING
    CAS9_TYPE_NAMING = [
        "-Cas9-pos",                                   # tag to add to Cas9 positive files
        "-Cas9-neg"                                    # tag to add to Cas9 negative files
    ]

    global REPLICATE_SELECTED
    REPLICATE_SELECTED = REPLICATE_NAMING[0]
