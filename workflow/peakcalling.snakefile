#####################################
## Snakefile for ChIP peak calling ##
#####################################

import os
#conda_prefix = str(os.environ["CONDA_PREFIX"])

import sys
#sys.path.insert(1, conda_prefix+"/lib/python"+str(sys.version_info[0])+"."+str(sys.version_info[1])+"/site-packages")

from typing import List
import pathlib
import re
import numpy
import pandas as pd
import math
from itertools import combinations


### Define the genomes (ZWART source)
genomes_table = pd.DataFrame({'genome_id': ['hg38', 'hg19'],
                              'fasta': ['/shared/data/Zwartlab/snakepipes_indices/hg38/BWAIndex/genome.fa',
                                        '/shared/data/Zwartlab/snakepipes_indices/hg19/BWAIndex/genome.fa'],
                              'blacklist': ['/shared/data/Zwartlab/snakepipes_indices/hg38/annotation/blacklist.bed',
                                            '/shared/data/Zwartlab/snakepipes_indices/hg19/annotation/blacklist.bed'],
                              'effective_genomeSize': [2945849067, 2861343702],
                              'ignore_for_normalization': ["KI270728.1 KI270727.1 KI270442.1 KI270729.1 GL000225.1 KI270743.1 GL000008.2 GL000009.2 KI270747.1 KI270722.1 GL000194.1 KI270742.1 GL000205.2 GL000195.1 KI270736.1 KI270733.1 GL000224.1 GL000219.1 KI270719.1 GL000216.2 KI270712.1 KI270706.1 KI270725.1 KI270744.1 KI270734.1 GL000213.1 GL000220.1 KI270715.1 GL000218.1 KI270749.1 KI270741.1 GL000221.1 KI270716.1 KI270731.1 KI270751.1 KI270750.1 KI270519.1 GL000214.1 KI270708.1 KI270730.1 KI270438.1 KI270737.1 KI270721.1 KI270738.1 KI270748.1 KI270435.1 GL000208.1 KI270538.1 KI270756.1 KI270739.1 KI270757.1 KI270709.1 KI270746.1 KI270753.1 KI270589.1 KI270726.1 KI270735.1 KI270711.1 KI270745.1 KI270714.1 KI270732.1 KI270713.1 KI270754.1 KI270710.1 KI270717.1 KI270724.1 KI270720.1 KI270723.1 KI270718.1 KI270317.1 KI270740.1 KI270755.1 KI270707.1 KI270579.1 KI270752.1 KI270512.1 KI270322.1 GL000226.1 KI270311.1 KI270366.1 KI270511.1 KI270448.1 KI270521.1 KI270581.1 KI270582.1 KI270515.1 KI270588.1 KI270591.1 KI270522.1 KI270507.1 KI270590.1 KI270584.1 KI270320.1 KI270382.1 KI270468.1 KI270467.1 KI270362.1 KI270517.1 KI270593.1 KI270528.1 KI270587.1 KI270364.1 KI270371.1 KI270333.1 KI270374.1 KI270411.1 KI270414.1 KI270510.1 KI270390.1 KI270375.1 KI270420.1 KI270509.1 KI270315.1 KI270302.1 KI270518.1 KI270530.1 KI270304.1 KI270418.1 KI270424.1 KI270417.1 KI270508.1 KI270303.1 KI270381.1 KI270529.1 KI270425.1 KI270396.1 KI270363.1 KI270386.1 KI270465.1 KI270383.1 KI270384.1 KI270330.1 KI270372.1 KI270548.1 KI270580.1 KI270387.1 KI270391.1 KI270305.1 KI270373.1 KI270422.1 KI270316.1 KI270338.1 KI270340.1 KI270583.1 KI270334.1 KI270429.1 KI270393.1 KI270516.1 KI270389.1 KI270466.1 KI270388.1 KI270544.1 KI270310.1 KI270412.1 KI270395.1 KI270376.1 KI270337.1 KI270335.1 KI270378.1 KI270379.1 KI270329.1 KI270419.1 KI270336.1 KI270312.1 KI270539.1 KI270385.1 KI270423.1 KI270392.1 KI270394.1 X Y MT",
                                                           "X Y MT",
                                                           "KI270728.1 KI270727.1 KI270442.1 KI270729.1 GL000225.1 KI270743.1 GL000008.2 GL000009.2 KI270747.1 KI270722.1 GL000194.1 KI270742.1 GL000205.2 GL000195.1 KI270736.1 KI270733.1 GL000224.1 GL000219.1 KI270719.1 GL000216.2 KI270712.1 KI270706.1 KI270725.1 KI270744.1 KI270734.1 GL000213.1 GL000220.1 KI270715.1 GL000218.1 KI270749.1 KI270741.1 GL000221.1 KI270716.1 KI270731.1 KI270751.1 KI270750.1 KI270519.1 GL000214.1 KI270708.1 KI270730.1 KI270438.1 KI270737.1 KI270721.1 KI270738.1 KI270748.1 KI270435.1 GL000208.1 KI270538.1 KI270756.1 KI270739.1 KI270757.1 KI270709.1 KI270746.1 KI270753.1 KI270589.1 KI270726.1 KI270735.1 KI270711.1 KI270745.1 KI270714.1 KI270732.1 KI270713.1 KI270754.1 KI270710.1 KI270717.1 KI270724.1 KI270720.1 KI270723.1 KI270718.1 KI270317.1 KI270740.1 KI270755.1 KI270707.1 KI270579.1 KI270752.1 KI270512.1 KI270322.1 GL000226.1 KI270311.1 KI270366.1 KI270511.1 KI270448.1 KI270521.1 KI270581.1 KI270582.1 KI270515.1 KI270588.1 KI270591.1 KI270522.1 KI270507.1 KI270590.1 KI270584.1 KI270320.1 KI270382.1 KI270468.1 KI270467.1 KI270362.1 KI270517.1 KI270593.1 KI270528.1 KI270587.1 KI270364.1 KI270371.1 KI270333.1 KI270374.1 KI270411.1 KI270414.1 KI270510.1 KI270390.1 KI270375.1 KI270420.1 KI270509.1 KI270315.1 KI270302.1 KI270518.1 KI270530.1 KI270304.1 KI270418.1 KI270424.1 KI270417.1 KI270508.1 KI270303.1 KI270381.1 KI270529.1 KI270425.1 KI270396.1 KI270363.1 KI270386.1 KI270465.1 KI270383.1 KI270384.1 KI270330.1 KI270372.1 KI270548.1 KI270580.1 KI270387.1 KI270391.1 KI270305.1 KI270373.1 KI270422.1 KI270316.1 KI270338.1 KI270340.1 KI270583.1 KI270334.1 KI270429.1 KI270393.1 KI270516.1 KI270389.1 KI270466.1 KI270388.1 KI270544.1 KI270310.1 KI270412.1 KI270395.1 KI270376.1 KI270337.1 KI270335.1 KI270378.1 KI270379.1 KI270329.1 KI270419.1 KI270336.1 KI270312.1 KI270539.1 KI270385.1 KI270423.1 KI270392.1 KI270394.1 X Y MT M"
                                                         ]})

# Define general variables
genome_used = (str(config["genome"])).lower()
blacklist = (genomes_table[genomes_table['genome_id']==genome_used]).blacklist.iloc[0]
genomeSize = (genomes_table[genomes_table['genome_id']==genome_used]).effective_genomeSize.iloc[0]
ignore_for_normalization = (genomes_table[genomes_table['genome_id']==genome_used]).ignore_for_normalization.iloc[0]

if ((eval(str(config["paired_end"])) == True)):
    read_extension = "--extendReads"
else:
    read_extension = "--extendReads "+str(config["fragment_length"])


if ((eval(str(config["use_macs3"])) == True)):
    macs_version = "macs3"
else:
    macs_version = "macs2"


### working directory
home_dir = os.path.join(config["output_directory"],"")
shell('mkdir -p {home_dir}')
workdir: home_dir


### get the unique samples names and other variables
# loading the sample table
sample_metadata = pd.read_csv(str(config["sample_config_table"]),  sep='\t+', engine='python')   # target_id | input_id | broad
sample_metadata = sample_metadata.iloc[:,0:3].set_axis(['target_id', 'input_id', 'broad'], axis=1, inplace=False)
TARGETNAMES = list(numpy.unique(list(sample_metadata.target_id)))
INPUTNAMES = list(numpy.unique(list(sample_metadata.input_id)))
SAMPLENAMES = list(numpy.unique(TARGETNAMES + INPUTNAMES))


# Get bam list
if not (os.path.exists(config["runs_directory"])):
    os.system("printf '\033[1;31m\\n!!! *runs_directory* does not exist !!!\\n\\n\033[0m'")
else:
    BAMS = next(os.walk(config["runs_directory"]))[2]
    RUNNAMES = numpy.unique([re.sub(rf"{config['bam_suffix']}$", "", i) for i in BAMS])



### MACS2 or MACS3?
if ((eval(str(config["use_macs3"])) == True)):
    PEAKCALLER = "macs3"
else:
    PEAKCALLER = "macs2"

### other generic variable
if ((eval(str(config["remove_duplicates"])) == True)):
    DUP = "dedup"
else:
    DUP = "mdup"



### Optional analysis outputs
if ((eval(str(config["perform_plotFingerprint"])) == True)):
    plotFingerprint_results = expand("05_Quality_controls_and_statistics/plotFingerprint/{target}_fingerPrinting_plot.pdf", target = TARGETNAMES)
else:
    plotFingerprint_results = []

if ((eval(str(config["perform_fragmentSizeDistribution"])) == True) & (eval(str(config["paired_end"])) == True)):
    fragmentSizeDistribution_results = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_metrics.txt"
else:
    fragmentSizeDistribution_results = []

if (len(SAMPLENAMES) > 2):
    PCA_wholeGenome_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome.pdf"
    PCA_wholeGenome_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome.pdf"
else:
    PCA_wholeGenome_12 = []
    PCA_wholeGenome_23 = []


if (len(TARGETNAMES) > 2):
    PCA_atPeaks_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.1-2_heatmap_atPeaks.pdf"
    PCA_atPeaks_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.2-3_heatmap_atPeaks.pdf"
else:
    PCA_atPeaks_12 = []
    PCA_atPeaks_23 = []


if (len(SAMPLENAMES) > 1):
    correlation_heatmap_wholeGenome_pearson = expand("05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome.{ext}", ext = ["pdf", "txt"])
    correlation_heatmap_wholeGenome_spearman = expand("05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome.{ext}", ext = ["pdf", "txt"])
else:
    correlation_heatmap_wholeGenome_pearson = []
    correlation_heatmap_wholeGenome_spearman = []


if (len(TARGETNAMES) > 1):
    correlation_heatmap_atPeaks_pearson = expand("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks.{ext}", ext = ["pdf", "txt"])
    correlation_heatmap_atPeaks_spearman = expand("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks.{ext}", ext = ["pdf", "txt"])
else:
    correlation_heatmap_atPeaks_pearson = []
    correlation_heatmap_atPeaks_spearman = []


if ((eval(str(config["paired_end"])) == True)):
    peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
else:
    peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES)


### Generation of global wildcard_constraints
# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
            raise ValueError("constraint_to(): Expected a list, got str instead")
    return "({})".format("|".join(values))

wildcard_constraints:
    SAMPLE = constraint_to(SAMPLENAMES),
    TARGET = constraint_to(TARGETNAMES),
    INPUT = constraint_to(INPUTNAMES)


ruleorder: fastQC_filtered_BAM > normalized_bigWig > raw_bigWig

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all funtions
if (set(SAMPLENAMES) <= set(RUNNAMES)):
    rule AAA_initialization:
        input:
            fastqc_bam_zip = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            plotFingerprint_results = plotFingerprint_results,
            fragmentSizeDistribution_results = fragmentSizeDistribution_results,
            normalized_bigWig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
            raw_bigWig = expand(os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_raw.coverage_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
            correlation_heatmap_wholeGenome_pearson = correlation_heatmap_wholeGenome_pearson,
            correlation_heatmap_wholeGenome_spearman = correlation_heatmap_wholeGenome_spearman,
            PCA_wholeGenome_12 = PCA_wholeGenome_12,
            PCA_wholeGenome_23 = PCA_wholeGenome_23,
            peaks = peaks,
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html",
            correlation_heatmap_atPeaks_pearson = correlation_heatmap_atPeaks_pearson,
            correlation_heatmap_atPeaks_spearman = correlation_heatmap_atPeaks_spearman,
            PCA_atPeaks_12 = PCA_atPeaks_12,
            PCA_atPeaks_23 = PCA_atPeaks_23,
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        shell:
            """
            printf '\033[1;36mPipeline ended!\\n\033[0m'
            """
else:
    missing_samples = '\\n  - '.join(list(set(SAMPLENAMES) - set(RUNNAMES)))
    os.system("printf '\033[1;31m\\n!!! Not all bam files are avalable in the input directory !!!\\n\\nPlease provide files for:\\n  - "+missing_samples+"\\n\\n\033[0m'")

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================

if (eval(str(config["skip_bam_filtering"])) == False):
    rule MAPQ_filter:
        input:
            source_bam = os.path.join(config["runs_directory"], ''.join(["{SAMPLE}", config["bam_suffix"]]))
        output:
            bam_mapq_only = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), ".bam"]))),
            bam_mapq_only_sorted = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"]))),
            bam_mapq_only_sorted_index = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam.bai"])))
        params:
            sample = "{SAMPLE}",
            MAPQ_threshold = config["MAPQ_threshold"]
        threads:
            max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
        shell:
            """
            printf '\033[1;36m{params.sample}: filtering MAPQ and re-indexing...\\n\033[0m'

            $CONDA_PREFIX/bin/samtools view -@ {threads} -h -q {params.MAPQ_threshold} {input.source_bam} -o {output.bam_mapq_only}

            $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted} {output.bam_mapq_only_sorted_index}
            """


    if ((eval(str(config["paired_end"])) == True) & (eval(str(config["umi_present"])) == True)):
        rule gatk4_markdups_umiAware:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
                umi_metrics = "01_BAM_filtered/umi_metrics/{SAMPLE}_UMI_metrics.txt",
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["remove_duplicates"])).lower(),
                sample = "{SAMPLE}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
            threads:
                workflow.cores
            shell:
                """
                printf '\033[1;36m{params.sample}: UMI-aware gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/umi_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

                $CONDA_PREFIX/bin/gatk UmiAwareMarkDuplicatesWithMateCigar \
                --INPUT {input.bam_mapq_only_sorted} \
                --OUTPUT {output.bam_mdup} \
                --REMOVE_DUPLICATES {params.remove_duplicates} \
                --MAX_EDIT_DISTANCE_TO_JOIN 1 \
                --UMI_METRICS_FILE {output.umi_metrics} \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --UMI_TAG_NAME RX \
                --CREATE_INDEX true \
                --VALIDATION_STRINGENCY STRICT \
                --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

                $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
                """
    else: # Single-end/no-UMI dedup
        rule gatk4_markdups:
            input:
                bam_mapq_only_sorted = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam"])),
                bam_mapq_only_sorted_index = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_sorted.bam.bai"]))
            output:
                bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
                bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
                dup_metrics = "01_BAM_filtered/MarkDuplicates_metrics/{SAMPLE}_MarkDuplicates_metrics.txt",
                flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
            params:
                remove_duplicates = (str(config["remove_duplicates"])).lower(),
                sample = "{SAMPLE}"
            log:
                out = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.out",
                err = "01_BAM_filtered/MarkDuplicates_logs/{SAMPLE}_MarkDuplicates.err"
            threads:
                workflow.cores
            shell:
                """
                printf '\033[1;36m{params.sample}: 'standard' gatk MarkDuplicates...\\n\033[0m'

                mkdir -p 01_BAM_filtered/MarkDuplicates_metrics
                mkdir -p 01_BAM_filtered/MarkDuplicates_logs
                mkdir -p 01_BAM_filtered/flagstat

                $CONDA_PREFIX/bin/gatk MarkDuplicates \
                --INPUT {input.bam_mapq_only_sorted} \
                --OUTPUT {output.bam_mdup} \
                --REMOVE_DUPLICATES {params.remove_duplicates} \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --CREATE_INDEX true \
                --VALIDATION_STRINGENCY LENIENT \
                --METRICS_FILE {output.dup_metrics} 2> {log.out} > {log.err}

                $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
                """
else:
    rule bam_link__skip_filtering:
        input:
            source_bam = os.path.join(config["runs_directory"], ''.join(["{SAMPLE}", config["bam_suffix"]]))
        output:
            bam_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            bai_mdup = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"]))
        params:
            sample = "{SAMPLE}"
        threads:
            max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
        shell:
            """
            printf '\033[1;36m{params.sample} (skip filtering): linking bam, indexing and computing flagstat...\\n\033[0m'

            mkdir -p 01_BAM_filtered/flagstat

            BAM_REAL=$(realpath {input.source_bam})
            ln -s $BAM_REAL {output.bam_mdup}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mdup} {output.bai_mdup}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """



rule fastQC_filtered_BAM:
    input:
        bam_mapq = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai_mapq = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        html = os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.html"])),
        zip = os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join("02_fastQC_on_BAM_filtered/"),
        sample = "{SAMPLE}"
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    shell:
        """
        mkdir -p 02_fastQC_on_BAM_filtered

        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir {params.fastQC_BAMs_outdir} {input.bam_mapq}
        """

# ------------------------------------------------------------------------------

rule plotFingerprint:
    input:
        target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        lorenz_curve_pdf = os.path.join("05_Quality_controls_and_statistics/plotFingerprint/{TARGET}_fingerPrinting_plot.pdf"),
        quality_metrics = os.path.join("05_Quality_controls_and_statistics/plotFingerprint/quality_metrics/{TARGET}_fingerPrinting_quality_metrics.txt")
    params:
        sample = "{TARGET}",
        sample_config_table = config["sample_config_table"],
        input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
        read_extension = read_extension,
        blacklist = blacklist
    threads:
        max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
    log:
        out = "05_Quality_controls_and_statistics/plotFingerprint/logs/{TARGET}_fingerPrinting_log.out",
        err = "05_Quality_controls_and_statistics/plotFingerprint/logs/{TARGET}_fingerPrinting_log.err"
    shell:
        """
        printf '\033[1;36m{params.sample}: plotting fingerprint...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/plotFingerprint/logs
        mkdir -p 05_Quality_controls_and_statistics/plotFingerprint/quality_metrics

        INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)

        $CONDA_PREFIX/bin/plotFingerprint \
        -b {input.target_bam} \
        01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
        --JSDsample 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
        -plot {output.lorenz_curve_pdf} \
        {params.read_extension} \
        --ignoreDuplicates \
        --outQualityMetrics {output.quality_metrics} \
        --labels {params.sample} ${{INPUT_ID}} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule fragmentSizeDistribution:
    input:
        all_bams = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), sample = SAMPLENAMES),
        all_bais = expand(os.path.join("01_BAM_filtered", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])), sample = SAMPLENAMES)
    output:
        fragment_distribution_plot = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_plot.pdf",
        fragmentSize_metrics = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_metrics.txt",
        fragmentSize_RawFragmentLengths = "05_Quality_controls_and_statistics/fragmentSize_distribution/fragmentSize_distribution_RawFragmentLengths.tab"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist
    threads:
        workflow.cores
    log:
        out = "05_Quality_controls_and_statistics/fragmentSize_distribution/logs/fragmentSize_distribution_log.out",
        err = "05_Quality_controls_and_statistics/fragmentSize_distribution/logs/fragmentSize_distribution_log.err"
    shell:
        """
        printf '\033[1;36mPlotting fragment size distribution...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/fragmentSize_distribution/logs

        $CONDA_PREFIX/bin/bamPEFragmentSize \
        --bamfiles {input.all_bams} \
        --binSize 1000000 \
        --blackListFileName {params.blacklist} \
        --samplesLabel {params.labels} \
        --histogram {output.fragment_distribution_plot} \
        --table {output.fragmentSize_metrics} \
        --outRawFragmentLengths {output.fragmentSize_RawFragmentLengths} \
        -p {threads} > {log.out} 2> {log.err}
        """

# ------------------------------------------------------------------------------

rule normalized_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        normalized_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])),
    params:
        sample = "{SAMPLE}",
        blacklist = blacklist,
        genomeSize = genomeSize,
        ignore_for_normalization = ignore_for_normalization,
        read_extension = read_extension,
        bw_binSize = config["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLE}_fragmentSize_distribution_log.out",
        err = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLE}_fragmentSize_distribution_log.err"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating RPGC normalized bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/RPGC_normalized/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.normalized_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing RPGC \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        --ignoreDuplicates \
        {params.read_extension} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule raw_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"]))
    output:
        raw_bigWig = os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{SAMPLE}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_raw.coverage_bs", str(config["bigWig_binSize"]), ".bw"])),
    params:
        sample = "{SAMPLE}",
        blacklist = blacklist,
        genomeSize = genomeSize,
        ignore_for_normalization = ignore_for_normalization,
        read_extension = read_extension,
        bw_binSize = config["bigWig_binSize"]
    threads:
        max(math.floor(workflow.cores/len(SAMPLENAMES)), 1)
    log:
        out = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLE}_fragmentSize_distribution_log.out",
        err = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLE}_fragmentSize_distribution_log.err"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating raw coverage bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/raw_coverage/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.raw_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing None \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        --ignoreDuplicates \
        {params.read_extension} \
        -p {threads} > {log.out} 2> {log.err}
        """

# ------------------------------------------------------------------------------

rule multiBigwigSummary_wholeGenome:
    input:
        all_norm_bigwig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), sample = SAMPLENAMES),
    output:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome.npz",
        multiBigWig_matrix_wholeGenome_raw = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome.txt"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads:
        workflow.cores
    log:
        out = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_log.out",
        err = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/multiBigWigSummary_matrix_wholeGenome_log.err"
    shell:
        """
        printf '\033[1;36mComputing multiBigwigSummary matrix (whole genome)...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs

        $CONDA_PREFIX/bin/multiBigwigSummary bins \
        -b {input.all_norm_bigwig} \
        -o {output.multiBigWig_matrix_wholeGenome} \
        --outRawCounts {output.multiBigWig_matrix_wholeGenome_raw} \
        --labels {params.labels} \
        --binSize 1000 \
        --chromosomesToSkip {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule correlations_wholeGenome:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome.npz"
    output:
        correlation_heatmap_wholeGenome_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome.pdf",
        correlation_heatmap_wholeGenome_pearson_tb = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_pearson.correlation_heatmap_wholeGenome.txt",
        correlation_heatmap_wholeGenome_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome.pdf",
        correlation_heatmap_wholeGenome_spearman_tb = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_spearman.correlation_heatmap_wholeGenome.txt"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_log.out",
        err_pearson = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_pearson.correlation_heatmap_wholeGenome_log.err",
        out_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_log.out",
        err_spearman = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_spearman.correlation_heatmap_wholeGenome_log.err"
    shell:
        """
        printf '\033[1;36mPlotting sample correlations (whole genome)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Pearson correlation whole genome RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_wholeGenome_pearson} \
        --outFileCorMatrix {output.correlation_heatmap_wholeGenome_pearson_tb} \
        --colorMap {params.heatmap_color} > {log.out_pearson} 2> {log.err_pearson}


        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Spearman correlation whole genome RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_wholeGenome_spearman} \
        --outFileCorMatrix {output.correlation_heatmap_wholeGenome_spearman_tb} \
        --colorMap {params.heatmap_color} > {log.out_spearman} 2> {log.err_spearman}
        """



rule PCA_wholeGenome:
    input:
        multiBigWig_matrix_wholeGenome = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/multiBigWigSummary_matrix_wholeGenome.npz"
    output:
        PCA_wholeGenome_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.1-2_heatmap_wholeGenome.pdf",
        PCA_wholeGenome_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/sample_correlation_PCA.2-3_heatmap_wholeGenome.pdf"
    params:
        labels = ' '.join(SAMPLENAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_log.out",
        err_12 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.1-2_heatmap_wholeGenome_log.err",
        out_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_log.out",
        err_23 = "05_Quality_controls_and_statistics/sample_comparisons_wholeGenome/logs/sample_correlation_PCA.2-3_heatmap_wholeGenome_log.err"
    shell:
        """
        printf '\033[1;36mPlotting PCA (whole genome)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 1 2 \
        --plotTitle 'PCA whole genome: PC1 vs PC2 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_12} > {log.out_12} 2> {log.err_12}

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_wholeGenome} \
        --labels {params.labels} \
        --PCs 2 3 \
        --plotTitle 'PCA whole genome: PC2 vs PC3 (RPGC normalized coverage)' \
        --plotFile {output.PCA_wholeGenome_23} > {log.out_23} 2> {log.err_23}
        """

# ------------------------------------------------------------------------------

if ((eval(str(config["paired_end"])) == True)):
    rule macs_callpeak_PE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES),
        output:
            peaksPE = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
        params:
            sample = "{TARGET}",
            macs_version = macs_version,
            sample_config_table = config["sample_config_table"],
            input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            read_extension = read_extension,
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        log:
            out = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.out",
            err = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAMPE_log.err"
        shell:
            """
            printf '\033[1;36m{params.sample}: calling peaks ({params.macs_version})...\\n\033[0m'

            mkdir -p 04_Called_peaks/logs

            INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)
            CALL_BROAD=$(grep -w {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                BROAD=""
            else
                BROAD="--broad"
            fi

            $CONDA_PREFIX/bin/{params.macs_version} callpeak \
            -t {input.target_bam} \
            -c 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
            -f BAMPE \
            -g {params.genomeSize} \
            -q {params.macs_qValue_cutoff} \
            --keep-dup all \
            --outdir 04_Called_peaks \
            --name {params.sample}.filtered.BAMPE ${{BROAD}} > {log.err} 2> {log.out}
            """
else:
    rule phantom_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES)
        output:
            phantom = '04_Called_peaks/phantom/{TARGET}.phantom.spp.out'
        log:
            out = '04_Called_peaks/phantom/logs/{TARGET}.phantom.log'
        params:
            sample = "{TARGET}",
            sample_config_table = config["sample_config_table"],
            input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36m{params.sample}: calculating phantom peak...\\n\033[0m'
            mkdir -p 04_Called_peaks/phantom/logs

            INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)

            ${{CONDA_PREFIX}}/bin/Rscript ${{CONDA_PREFIX}}/bin/run_spp.R -rf -c='{input.target_bam}' -i="01_BAM_filtered/${{INPUT_ID}}{params.input_suffix}" -savp -out={output.phantom} &> {log.out}
            """



    rule fragment_length:
        input:
            phantom = '04_Called_peaks/phantom/{TARGET}.phantom.spp.out'
        output:
            fragment_length_phanthom = temp('04_Called_peaks/phantom/{TARGET}.fragment_length')
        shell:
            """
            awk '{{print $3}}' < {input.phantom} | tr ',' '\\t' | awk '{{if($1!=0) print $1; else print $2}}' > {output.fragment_length_phanthom}
            """



    rule macs_callpeak_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
            input_bam_all = expand(os.path.join("01_BAM_filtered", ''.join(["{input}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])), input = INPUTNAMES),
            phantom = '04_Called_peaks/phantom/{TARGET}.fragment_length'
        output:
            peaksSE = "04_Called_peaks/{TARGET}.filtered.BAM_peaks.xls"
        params:
            sample = "{TARGET}",
            macs_version = macs_version,
            sample_config_table = config["sample_config_table"],
            input_suffix = "_mapq"+str(config["MAPQ_threshold"])+"_"+DUP+"_sorted.bam",
            genomeSize = genomeSize,
            macs_qValue_cutoff = config["macs_qValue_cutoff"],
            blacklist = blacklist
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        log:
            out = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAM_log.out",
            err = "04_Called_peaks/logs/{TARGET}_macs.callpeak.BAM_log.err"
        shell:
            """
            printf '\033[1;36m{params.sample}: calling peaks ({params.macs_version})...\\n\033[0m'

            mkdir -p 04_Called_peaks/logs

            INPUT_ID=$(grep -w {params.sample} {params.sample_config_table} | cut -f 2)
            CALL_BROAD=$(grep -w {params.sample} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                BROAD=""
            else
                BROAD="--broad"
            fi

            EXTSIZEPHANTOM=$(cat {input.phantom}) ${{BROAD}}

            if [ "$EXTSIZEPHANTOM" -lt 1 ]; then
              EXTSIZEPHANTOM=200
            fi

            $CONDA_PREFIX/bin/{params.macs_version} callpeak \
            -t {input.target_bam} \
            -c 01_BAM_filtered/${{INPUT_ID}}{params.input_suffix} \
            -f BAM \
            -g {params.genomeSize} \
            --nomodel \
            -q {params.macs_qValue_cutoff} \
            --outdir 04_Called_peaks \
            --name {params.sample}.filtered.BAM \
            --extsize $EXTSIZEPHANTOM > {log.err} 2> {log.out}
            """

# ------------------------------------------------------------------------------

if (eval(str(config["skip_bam_filtering"])) == False):
    picard_metrics_file = expand("01_BAM_filtered/MarkDuplicates_metrics/{sample}_MarkDuplicates_metrics.txt", sample = SAMPLENAMES)
    picard_metrics_dir = "01_BAM_filtered/MarkDuplicates_metrics"
else:
    picard_metrics_file = []
    picard_metrics_dir = []


if ((eval(str(config["paired_end"])) == True)):
    rule multiQC_PE:
        input:
            fastqc = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            picard_metrics = picard_metrics_file,
            flagstat = expand(os.path.join("01_BAM_filtered/flagstat/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"])), sample = SAMPLENAMES),
            peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
        output:
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html"
        params:
            out_directory = "05_Quality_controls_and_statistics/multiQC/",
            multiqc_report_name = "multiQC_report.html",
            picard_metrics_dir = picard_metrics_dir
        threads: 1
        log:
            out = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.out",
            err = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.err"
        shell:
            """
            mkdir -p 05_Quality_controls_and_statistics/multiQC/

            $CONDA_PREFIX/bin/multiqc -f \
            -o {params.out_directory} \
            -n {params.multiqc_report_name} \
            --dirs 02_fastQC_on_BAM_filtered {params.picard_metrics_dir} 01_BAM_filtered/flagstat 04_Called_peaks > {log.err} 2> {log.out}
            """
else:
    rule multiQC_SE:
        input:
            fastqc = expand(os.path.join("02_fastQC_on_BAM_filtered/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_fastqc.zip"])), sample = SAMPLENAMES),
            picard_metrics = picard_metrics_file,
            flagstat = expand(os.path.join("01_BAM_filtered/flagstat/", ''.join(["{sample}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted_flagstat.txt"])), sample = SAMPLENAMES),
            peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES),
            phanthom = expand('04_Called_peaks/phantom/{target}.phantom.spp.out', target = TARGETNAMES)
        output:
            multiqc_report = "05_Quality_controls_and_statistics/multiQC/multiQC_report.html"
        params:
            out_directory = "05_Quality_controls_and_statistics/multiQC/",
            multiqc_report_name = "multiQC_report.html",
            picard_metrics_dir = picard_metrics_dir
        threads: 1
        log:
            out = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.out",
            err = "05_Quality_controls_and_statistics/multiQC/multiQC_report_log.err"
        shell:
            """
            printf '\033[1;36mGenerating multiQC report...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/multiQC/

            $CONDA_PREFIX/bin/multiqc -f \
            -o {params.out_directory} \
            -n {params.multiqc_report_name} \
            --dirs 02_fastQC_on_BAM_filtered {params.picard_metrics_dir} 01_BAM_filtered/flagstat 04_Called_peaks > {log.err} 2> {log.out}
            """

# ------------------------------------------------------------------------------

if ((eval(str(config["paired_end"])) == True)):
    rule merge_all_peaks_PE:
        input:
            peaks = expand("04_Called_peaks/{target}.filtered.BAMPE_peaks.xls", target = TARGETNAMES)
        output:
            concat_peaks = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat.bed"),
            concat_peaks_sorted = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat_sorted.bed"),
            merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
        params:
            labels = ' '.join(SAMPLENAMES),
            blacklist = blacklist,
            ignore_for_normalization = ignore_for_normalization
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36mMerging all peaks...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/

            PEAK_LIST=$(ls 04_Called_peaks/*_peaks.*Peak | grep -v gapped)

            for i in ${{PEAK_LIST}}
            do
                cut -f 1-6 $i >> {output.concat_peaks}
            done

            sort -V -k1,1 -k2,2 {output.concat_peaks} > {output.concat_peaks_sorted}

            $CONDA_PREFIX/bin/bedtools merge -i {output.concat_peaks_sorted} | sort -V -k1,1 -k2,2 > {output.merged_peaks_sorted}
            """
else:
    rule merge_all_peaks_SE:
        input:
            peaks = expand("04_Called_peaks/{target}.filtered.BAM_peaks.xls", target = TARGETNAMES)
        output:
            concat_peaks = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat.bed"),
            concat_peaks_sorted = temp("05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_concat_sorted.bed"),
            merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
        params:
            labels = ' '.join(SAMPLENAMES),
            blacklist = blacklist,
            ignore_for_normalization = ignore_for_normalization
        threads:
            workflow.cores
        shell:
            """
            printf '\033[1;36mMerging all peaks...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/

            PEAK_LIST=$(ls 04_Called_peaks/*_peaks.*Peak | grep -v gapped)

            for i in ${{PEAK_LIST}}
            do
                cut -f 1-6 $i >> {output.concat_peaks}
            done

            sort -V -k1,1 -k2,2 {output.concat_peaks} > {output.concat_peaks_sorted}

            $CONDA_PREFIX/bin/bedtools merge -i {output.concat_peaks_sorted} | sort -V -k1,1 -k2,2 > {output.merged_peaks_sorted}
            """



rule multiBigwigSummary_atPeaks:
    input:
        all_norm_bigwig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{target}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_RPGC.normalized_bs", str(config["bigWig_binSize"]), ".bw"])), target = TARGETNAMES),
        merged_peaks_sorted = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/all_peaks_merged_sorted.bed"
    output:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz",
        multiBigWig_matrix_atPeaks_raw = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.txt"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads:
        workflow.cores
    log:
        out = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/multiBigWigSummary_matrix_atPeaks_log.out",
        err = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/multiBigWigSummary_matrix_atPeaks_log.err"
    shell:
        """
        printf '\033[1;36mComputing multiBigwigSummary matrix (at peaks)...\\n\033[0m'

        mkdir -p 05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs

        $CONDA_PREFIX/bin/multiBigwigSummary BED-file \
        --BED {input.merged_peaks_sorted} \
        -b {input.all_norm_bigwig} \
        -o {output.multiBigWig_matrix_atPeaks} \
        --outRawCounts {output.multiBigWig_matrix_atPeaks_raw} \
        --labels {params.labels} \
        --binSize 1000 \
        --chromosomesToSkip {params.ignore_for_normalization} \
        --blackListFileName {params.blacklist} \
        -p {threads} > {log.out} 2> {log.err}
        """



rule correlations_atPeaks:
    input:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    output:
        correlation_heatmap_atPeaks_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks.pdf",
        correlation_heatmap_atPeaks_pearson_tb = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_pearson.correlation_heatmap_atPeaks.txt",
        correlation_heatmap_atPeaks_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks.pdf",
        correlation_heatmap_atPeaks_spearman_tb = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_spearman.correlation_heatmap_atPeaks.txt"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization,
        heatmap_color = config["correlation_heatmap_colorMap"]
    threads: 1
    log:
        out_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_pearson.correlation_heatmap_atPeaks_log.out",
        err_pearson = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_pearson.correlation_heatmap_atPeaks_log.err",
        out_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_spearman.correlation_heatmap_atPeaks_log.out",
        err_spearman = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_spearman.correlation_heatmap_atPeaks_log.err"
    shell:
        """
        printf '\033[1;36mPlotting sample correlations (at peaks)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Pearson correlation at peaks RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_atPeaks_pearson} \
        --outFileCorMatrix {output.correlation_heatmap_atPeaks_pearson_tb} \
        --colorMap {params.heatmap_color} > {log.out_pearson} 2> {log.err_pearson}


        $CONDA_PREFIX/bin/plotCorrelation \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --corMethod spearman \
        --whatToPlot heatmap \
        --skipZeros \
        --plotNumbers \
        --removeOutliers \
        --plotTitle 'Spearman correlation at peaks RPGC normalized coverage' \
        --plotFile {output.correlation_heatmap_atPeaks_spearman} \
        --outFileCorMatrix {output.correlation_heatmap_atPeaks_spearman_tb} \
        --colorMap {params.heatmap_color} > {log.out_spearman} 2> {log.err_spearman}
        """



rule PCA_atPeaks:
    input:
        multiBigWig_matrix_atPeaks = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/multiBigWigSummary_matrix_atPeaks.npz"
    output:
        PCA_atPeaks_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.1-2_heatmap_atPeaks.pdf",
        PCA_atPeaks_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/sample_correlation_PCA.2-3_heatmap_atPeaks.pdf"
    params:
        labels = ' '.join(TARGETNAMES),
        blacklist = blacklist,
        ignore_for_normalization = ignore_for_normalization
    threads: 1
    log:
        out_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.1-2_heatmap_atPeaks_log.out",
        err_12 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.1-2_heatmap_atPeaks_log.err",
        out_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.2-3_heatmap_atPeaks_log.out",
        err_23 = "05_Quality_controls_and_statistics/sample_comparisons_atPeaks/logs/sample_correlation_PCA.2-3_heatmap_atPeaks_log.err"
    shell:
        """
        printf '\033[1;36mPlotting PCA (at peaks)...\\n\033[0m'

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --PCs 1 2 \
        --plotTitle 'PCA at peaks: PC1 vs PC2 (RPGC normalized coverage)' \
        --plotFile {output.PCA_atPeaks_12} > {log.out_12} 2> {log.err_12}

        $CONDA_PREFIX/bin/plotPCA \
        -in {input.multiBigWig_matrix_atPeaks} \
        --labels {params.labels} \
        --PCs 2 3 \
        --plotTitle 'PCA at peaks: PC2 vs PC3 (RPGC normalized coverage)' \
        --plotFile {output.PCA_atPeaks_23} > {log.out_23} 2> {log.err_23}
        """

# ------------------------------------------------------------------------------

if ((eval(str(config["paired_end"])) == True)):
    rule MACS_peak_QC_PE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
	        peaks = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks.xls"
        output:
	        qc = temp("05_Quality_controls_and_statistics/peaks_stats/{TARGET}.filtered.BAMPE_peaks.qc.txt")
        params:
            sample_config_table = config["sample_config_table"],
            peak_prefix = "04_Called_peaks/{TARGET}.filtered.BAMPE_peaks",
            blacklist = blacklist,
            target = "{TARGET}",
            genomeSize = genomeSize
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        shell:
            """
            printf '\033[1;36m{params.target}: computing peak stats...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/peaks_stats/

            # define peak file
            CALL_BROAD=$(grep -w {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                CALLING_MODE="narrow"
                PEAK="{params.peak_prefix}.narrowPeak"
                PEAK_CHR="{params.peak_prefix}_chr.narrowPeak"
            else
                CALLING_MODE="broad"
                PEAK="{params.peak_prefix}.broadPeak"
                PEAK_CHR="{params.peak_prefix}_chr.broadPeak"
            fi

            # get the number of peaks
            peak_count=$(wc -l < $PEAK)

            # get the number of mapped reads
            mapped_reads=$($CONDA_PREFIX/bin/samtools view -c -F 4 {input.target_bam})

            # calculate the number of alignments overlapping the peaks
            # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
            reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

            # calculate Fraction of Reads In Peaks
            frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

            # compute peak genome coverage
            peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
            genome_size={params.genomeSize}
            genomecov=$(bc -l <<< "$peak_len/$genome_size")

            # rounding fractions
            genomecov_round=$(printf "%.5f\n" "$genomecov")
            frip_round=$(printf "%.3f\n" "$frip")

            # write peak-based QC metrics to output file
            printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}

	        # add chr to peak files
            $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a $PEAK -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > $PEAK_CHR
	        """

    rule aggregate_FRiP_PE:
        input:
            qc = expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAMPE_peaks.qc.txt", target = TARGETNAMES)
        output:
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        params:
            all_qc = ' '.join(expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAMPE_peaks.qc.txt", target = TARGETNAMES))
        threads: 1
        shell:
            """
            # print header of FRiP report
            printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
            cat {params.all_qc} >> {output.aggregated_qc}
            """

#*******************************************************************
else:
    rule MACS_peak_QC_SE:
        input:
            target_bam = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bam"])),
            target_bai = os.path.join("01_BAM_filtered", ''.join(["{TARGET}_mapq", str(config["MAPQ_threshold"]), "_", DUP, "_sorted.bai"])),
	        peaks = "04_Called_peaks/{TARGET}.filtered.BAM_peaks.xls"
        output:
	        qc = temp("05_Quality_controls_and_statistics/peaks_stats/{TARGET}.filtered.BAM_peaks.qc.txt")
        params:
            sample_config_table = config["sample_config_table"],
            peak_prefix = "04_Called_peaks/{TARGET}.filtered.BAM_peaks",
            blacklist = blacklist,
            target = "{TARGET}",
            genomeSize = genomeSize
        threads:
            max(math.floor(workflow.cores/len(TARGETNAMES)), 1)
        shell:
            """
            printf '\033[1;36m{params.target}: computing peak stats...\\n\033[0m'

            mkdir -p 05_Quality_controls_and_statistics/peaks_stats/

            # define peak file
            CALL_BROAD=$(grep -w {params.target} {params.sample_config_table} | cut -f 3 | sed -e 's/\\(.*\\)/\\L\\1/')

            if [ $CALL_BROAD == "false" ]; then
                CALLING_MODE="narrow"
                PEAK="{params.peak_prefix}.narrowPeak"
                PEAK_CHR="{params.peak_prefix}_chr.narrowPeak"
            else
                CALLING_MODE="broad"
                PEAK="{params.peak_prefix}.broadPeak"
                PEAK_CHR="{params.peak_prefix}_chr.broadPeak"
            fi

            # get the number of peaks
            peak_count=$(wc -l < $PEAK)

            # get the number of mapped reads
            mapped_reads=$($CONDA_PREFIX/bin/samtools view -c -F 4 {input.target_bam})

            # calculate the number of alignments overlapping the peaks
            # exclude reads flagged as unmapped (unmapped reads will be reported when using -L)
            reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -F 4 -L $PEAK {input.target_bam})

            # calculate Fraction of Reads In Peaks
            frip=$(bc -l <<< "$reads_in_peaks/$mapped_reads")

            # compute peak genome coverage
            peak_len=$(awk '{{total+=$3-$2}}END{{print total}}' $PEAK)
            genome_size={params.genomeSize}
            genomecov=$(bc -l <<< "$peak_len/$genome_size")

            # rounding fractions
            genomecov_round=$(printf "%.5f\n" "$genomecov")
            frip_round=$(printf "%.3f\n" "$frip")

            # write peak-based QC metrics to output file
            printf '{params.target}\\t'$CALLING_MODE'\\t'$peak_count'\\t'$frip_round'\\t'$genomecov_round'\\n' > {output.qc}

	        # add chr to peak files
            $CONDA_PREFIX/bin/bedtools subtract -nonamecheck -a $PEAK -b {params.blacklist} | awk '{{if (length($1) <3 && $1 !="MT"){{print "chr"$0}} else {{print $0}} }}' > $PEAK_CHR
	        """

    rule aggregate_FRiP_SE:
        input:
            qc = expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAM_peaks.qc.txt", target = TARGETNAMES)
        output:
            aggregated_qc = "05_Quality_controls_and_statistics/peaks_stats/all_samples_FRiP_report.tsv"
        params:
            all_qc = ' '.join(expand("05_Quality_controls_and_statistics/peaks_stats/{target}.filtered.BAM_peaks.qc.txt", target = TARGETNAMES))
        threads: 1
        shell:
            """
            # print header of FRiP report
            printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
            cat {params.all_qc} >> {output.aggregated_qc}
            """


# ------------------------------------------------------------------------------
#                                 END pipeline
# ------------------------------------------------------------------------------

