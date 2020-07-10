# -*- coding: utf-8 -*-

# lunarc configuration file
# configuration file with sample list (yaml list)
# mm10 genome directory

configfile: "/projects/fs3/raquelgg/msc/trim28/src/config_files/config.yaml"

LOCATION = config["samples"]
LOCATION_STR = config["samples_str"]

INVIVO_ADULT_KO = config["invivo_adult_ko"]
INVIVO_ADULT_CTRL = config["invivo_adult_ctrl"]

INVIVO_BD_KO = config["invivo_bd_ko"]
INVIVO_BD_CTRL = config["invivo_bd_ctrl"]

INVITRO_CRISPR_KO_G3 = config["invitro_crispr_ko_g3"]
INVITRO_CRISPR_CTRL = config["invitro_crispr_ctrl"]

INVITRO_CRISPR_KO_G4 = config["invitro_crispr_ko_g4"]
INVITRO_CRISPR_CTRL = config["invitro_crispr_ctrl"]

INVITRO_CRISPR_KO_G13 = config["invitro_crispr_ko_g13"]
INVITRO_CRISPR_CTRL = config["invitro_crispr_ctrl"]

INVIVO_CRISPR_KO_G3 = config["invivo_crispr_ko_g3"]
INVIVO_CRISPR_CTRL = config["invivo_crispr_ctrl"]

INVIVO_CRISPR_KO_G4 = config["invivo_crispr_ko_g4"]
INVIVO_CRISPR_CTRL = config["invivo_crispr_ctrl"]

INVIVO_CRISPR_KO_G13 = config["invivo_crispr_ko_g13"]
INVIVO_CRISPR_CTRL = config["invivo_crispr_ctrl"]

INVIVO_CRISPR_KO_G3G4G13 = config["invivo_crispr_ko_g3g4g13"]
INVIVO_CRISPR_CTRL = config["invivo_crispr_ctrl"]
INVIVO_BD = config["invivo_bd"]
INVIVO_CRISPR_KO = config["invivo_crispr_ko"]
INVITRO_CRISPR_KO = config["invitro_crispr_ko"]

CHIPSEQ = config["chipseq"]
CHIPSEQ_LOC = config["chipseq_location"]
CUTNRUN = config["cut_n_run"]
CUTNRUN_LOC = config["cut_n_run_location"]
CHIPSEQ_LOC_WO_INP = config["cut_n_run_location_wo_inputs"]
CHIPSEQ_INP = config["cut_n_run_location_inputs"]
ADULT_CHIPSEQ  = config["adult_wt"]

# Run as:
# snakemake -j 5 --cluster-config /projects/fs3/raquelgg/msc/trim28/src/config_files/lunarc_config.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} --tasks-per-node {cluster.tasks-per-node}  -t {cluster.time} -o {cluster.o} -e {cluster.e} -J {cluster.J} -N {cluster.N}" --latency-wait 60

rule all:
    input:
        "../8_TEchipseq/MMERVK10C/plots/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.pdf"
        #"../8_TEchipseq/MMERVK10C/matrices/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.mat"

rule mapping:
    input:
        "/projects/fs3/raquelgg/msc/trim28/data/raw_samples/{location}_R1.fastq.gz",
        "/projects/fs3/raquelgg/msc/trim28/data/raw_samples/{location}_R2.fastq.gz",
        "/projects/fs1/common/genome/lunarc/indicies/star/mouse/GRCm38.p6/",
        "/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf"
    output:
        "../1_mapping/rnaseq/{location}Aligned.out.bam",
        "../1_mapping/rnaseq/{location}Aligned.sortedByCoord.out.bam"
    shell:
        """
        #echo Mapping reads from {wildcards.location} to mm10!
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml STAR/2.6.0c

        STAR --runThreadN 10 \
        --readFilesCommand gunzip -c \
        --outSAMattributes All \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir {input[2]} \
        --sjdbGTFfile {input[3]} \
        --outFileNamePrefix ../1_mapping/rnaseq/{wildcards.location} \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNoverLmax 0.03  \
        --readFilesIn  {input[0]} {input[1]}
        module purge
        """
rule indexing:
    input:
        "../1_mapping/rnaseq/{location}Aligned.sortedByCoord.out.bam"
    output:
        "../1_mapping/rnaseq/{location}Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml SAMtools/1.4

        samtools index -b {input}

        module purge
        """
rule bigwig:
    input:
        "../1_mapping/rnaseq/{location}Aligned.sortedByCoord.out.bam",
        "../1_mapping/rnaseq/{location}Aligned.sortedByCoord.out.bam.bai"
    output:
        "../1_mapping/rnaseq/{location}Aligned.sortedByCoord.out.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml Python/3.5.2

        bamCoverage --normalizeUsingRPKM -b {input[0]} -o {output}

        module purge
        """

rule featureCounts_gene:
    input:
        annotation="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        samples=expand("../1_mapping/{sample}Aligned.sortedByCoord.out.bam", sample=LOCATION)
    output:
        "../3_readcount/3_gene/gene_count_matrix_2.csv"
    shell:
        """
        ml GCC/7.3.0-2.30  OpenMPI/3.1.1
        ml Subread/1.6.3

        featureCounts -F GTF -s 2 -g gene_id -a {input.annotation} -o {output[1]} {input.samples}

        module purge
        """
rule featureCounts_TEs:
    input:
        annotation="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE_LINE_SINE_LTR.gtf",
        samples=expand("../1_mapping/{sample}Aligned.sortedByCoord.out.bam", sample=LOCATION)
    output:
        "../3_readcount/1_TEs/TE_count_matrix_2.csv"
    shell:
        """
        ml GCC/7.3.0-2.30  OpenMPI/3.1.1
        ml Subread/1.6.3

        featureCounts -s2 -F GTF -g transcript_id -a {input.annotation} -o {output} {input.samples}

        module purge
        """
rule featureCounts_ERVK:
    input:
        annotation="../../../annotations/mm10/repeatmasker/mm10_rmsk_TE.gtf",
        samples=expand("../1_mapping/{sample}Aligned.sortedByCoord.out.bam", sample=LOCATION)
    output:
        "../10_neighborgenes/mm10_rmsk_ERVK.gtf",
        "../3_readcount/1_TEs/ERVK_count_matrix_2.csv"
    shell:
        """
        ml GCC/7.3.0-2.30  OpenMPI/3.1.1
        ml Subread/1.6.3

        grep -w ERVK {input.annotation} > {output[0]}
        featureCounts -s2 -F GTF -g transcript_id -a {output[0]} -o {output[1]} {input.samples}

        module purge
        """
rule featureCounts_fullMMERVK10C:
    input:
        annotation="../data/mm10/repeatmasker/mm10_fullERVs_merge_repeatmasker_parsed.gtf",
        samples=expand("../1_mapping/{sample}Aligned.sortedByCoord.out.bam", sample=LOCATION)
    output:
        "../3_readcount/1_TEs/fullMMERVK10C_count_matrix_2.csv"
    shell:
        """
        ml GCC/7.3.0-2.30  OpenMPI/3.1.1
        ml Subread/1.6.3

        featureCounts -s2 -F GTF -g transcript_id -a {input.annotation} -o {output} {input.samples}

        module purge
        """

rule multimapping:
    input:
        "/projects/fs3/raquelgg/msc/trim28/data/raw_samples/{location}_R1.fastq.gz",
        "/projects/fs3/raquelgg/msc/trim28/data/raw_samples/{location}_R2.fastq.gz",
        "/projects/fs1/common/genome/lunarc/indicies/star/mouse/GRCm38.p6/",
        "/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf"
    output:
        "../5_multimapping/{location}Aligned.out.bam",
        "../5_multimapping/{location}Aligned.sortedByCoord.out.bam"
    shell:
        """
        echo Mapping reads from {wildcards.location} to GRCm38.p6!
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml STAR/2.6.0c

        STAR --runThreadN 10 \
        --readFilesCommand gunzip -c \
        --outSAMattributes All \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile {input[3]} \
        --genomeDir {input[2]} \
        --outFileNamePrefix ../5_multimapping/{wildcards.location} \
        --outFilterMultimapNmax 100 \
        --winAnchorMultimapNmax 200  \
        --readFilesIn  {input[0]} {input[1]}
        module purge
        """
rule multimap_indexing:
    input:
        "../5_multimapping/{location}Aligned.sortedByCoord.out.bam"
    output:
        "../5_multimapping/{location}Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml SAMtools/1.4

        samtools index -b {input}

        module purge
        """
rule multimap_bigwig:
    input:
        "../5_multimapping/{location}Aligned.sortedByCoord.out.bam",
        "../5_multimapping/{location}Aligned.sortedByCoord.out.bam.bai"
    output:
        "../5_multimapping/{location}Aligned.sortedByCoord.out.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml Python/3.5.2

        bamCoverage --normalizeUsingRPKM -b {input[0]} -o {output}

        module purge
        """
rule TEtranscripts_invivo_adult:
    input:
        ko=expand("../5_multimapping/{invivo_adult_ko}Aligned.out.bam", invivo_adult_ko=INVIVO_ADULT_KO),
        ctrl=expand("../5_multimapping/{invivo_adult_ctrl}Aligned.out.bam", invivo_adult_ctrl=INVIVO_ADULT_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invivo_adult/invivo_adult.cntTable",
        "../6_TEtranscripts/invivo_adult/invivo_adult_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invivo_adult/invivo_adult --stranded reverse

        module purge
        """
rule TEtranscripts_invivo_braindev:
    input:
        ko=expand("../5_multimapping/{invivo_bd_ko}Aligned.out.bam", invivo_bd_ko=INVIVO_BD_KO),
        ctrl=expand("../5_multimapping/{invivo_bd_ctrl}Aligned.out.bam", invivo_bd_ctrl=INVIVO_BD_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invivo_bd/invivo_bd.cntTable",
        "../6_TEtranscripts/invivo_bd/invivo_bd_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invivo_bd/invivo_bd --stranded reverse

        module purge
        """
rule TEtranscripts_invitro_crispr_g3:
    input:
        ko=expand("../5_multimapping/{invitro_crispr_ko_g3}Aligned.out.bam", invitro_crispr_ko_g3=INVITRO_CRISPR_KO_G3),
        ctrl=expand("../5_multimapping/{invitro_crispr_ctrl}Aligned.out.bam", invitro_crispr_ctrl=INVITRO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invitro_crispr_g3/invitro_crispr_g3.cntTable",
        "../6_TEtranscripts/invitro_crispr_g3/invitro_crispr_g3_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invitro_crispr_g3/invitro_crispr_g3 --stranded reverse

        module purge
        """
rule TEtranscripts_invitro_crispr_g4:
    input:
        ko=expand("../5_multimapping/{invitro_crispr_ko_g4}Aligned.out.bam", invitro_crispr_ko_g4=INVITRO_CRISPR_KO_G4),
        ctrl=expand("../5_multimapping/{invitro_crispr_ctrl}Aligned.out.bam", invitro_crispr_ctrl=INVITRO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invitro_crispr_g4/invitro_crispr_g4.cntTable",
        "../6_TEtranscripts/invitro_crispr_g4/invitro_crispr_g4_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invitro_crispr_g4/invitro_crispr_g4 --stranded reverse

        module purge
        """
rule TEtranscripts_invitro_crispr_g13:
    input:
        ko=expand("../5_multimapping/{invitro_crispr_ko_g13}Aligned.out.bam", invitro_crispr_ko_g13=INVITRO_CRISPR_KO_G13),
        ctrl=expand("../5_multimapping/{invitro_crispr_ctrl}Aligned.out.bam", invitro_crispr_ctrl=INVITRO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invitro_crispr_g13/invitro_crispr_g13.cntTable",
        "../6_TEtranscripts/invitro_crispr_g13/invitro_crispr_g13_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invitro_crispr_g13/invitro_crispr_g13 --stranded reverse

        module purge
        """
rule TEtranscripts_invivo_crispr_g3:
    input:
        ko=expand("../5_multimapping/{invivo_crispr_ko_g3}Aligned.out.bam", invivo_crispr_ko_g3=INVIVO_CRISPR_KO_G3),
        ctrl=expand("../5_multimapping/{invivo_crispr_ctrl}Aligned.out.bam", invivo_crispr_ctrl=INVIVO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invivo_crispr_g3/invivo_crispr_g3.cntTable",
        "../6_TEtranscripts/invivo_crispr_g3/invivo_crispr_g3_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invivo_crispr_g3/invivo_crispr_g3 --stranded reverse

        module purge
        """
rule TEtranscripts_invivo_crispr_g4:
    input:
        ko=expand("../5_multimapping/{invivo_crispr_ko_g4}Aligned.out.bam", invivo_crispr_ko_g4=INVIVO_CRISPR_KO_G4),
        ctrl=expand("../5_multimapping/{invivo_crispr_ctrl}Aligned.out.bam", invivo_crispr_ctrl=INVIVO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invivo_crispr_g4/invivo_crispr_g4.cntTable",
        "../6_TEtranscripts/invivo_crispr_g4/invivo_crispr_g4_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invivo_crispr_g4/invivo_crispr_g4 --stranded reverse

        module purge
        """
rule TEtranscripts_invivo_crispr_g13:
    input:
        ko=expand("../5_multimapping/{invivo_crispr_ko_g13}Aligned.out.bam", invivo_crispr_ko_g13=INVIVO_CRISPR_KO_G13),
        ctrl=expand("../5_multimapping/{invivo_crispr_ctrl}Aligned.out.bam", invivo_crispr_ctrl=INVIVO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invivo_crispr_g13/invivo_crispr_g13.cntTable",
        "../6_TEtranscripts/invivo_crispr_g13/invivo_crispr_g13_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invivo_crispr_g13/invivo_crispr_g13 --stranded reverse

        module purge
        """
rule TEtranscripts_invivo_crispr_g3g4g13:
    input:
        ko=expand("../5_multimapping/{invivo_crispr_ko_g3g4g13}Aligned.out.bam", invivo_crispr_ko_g3g4g13=INVIVO_CRISPR_KO_G3G4G13),
        ctrl=expand("../5_multimapping/{invivo_crispr_ctrl}Aligned.out.bam", invivo_crispr_ctrl=INVIVO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invivo_crispr_g3g4g13/invivo_crispr_g3g4g13.cntTable",
        "../6_TEtranscripts/invivo_crispr_g3g4g13/invivo_crispr_g3g4g13_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invivo_crispr_g3g4g13/invivo_crispr_g3g4g13 --stranded reverse

        module purge
        """
rule TEtranscripts_invivo_crispr:
    input:
        ko=expand("../5_multimapping/{invivo_crispr_ko}Aligned.out.bam", invivo_crispr_ko=INVIVO_CRISPR_KO),
        ctrl=expand("../5_multimapping/{invivo_crispr_ctrl}Aligned.out.bam", invivo_crispr_ctrl=INVIVO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invivo_crispr/invivo_crispr.cntTable",
        "../6_TEtranscripts/invivo_crispr/invivo_crispr_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invivo_crispr/invivo_crispr --stranded reverse

        module purge
        """
rule TEtranscripts_invitro_crispr:
    input:
        ko=expand("../5_multimapping/{invitro_crispr_ko}Aligned.out.bam", invitro_crispr_ko=INVITRO_CRISPR_KO),
        ctrl=expand("../5_multimapping/{invitro_crispr_ctrl}Aligned.out.bam", invitro_crispr_ctrl=INVITRO_CRISPR_CTRL),
        gtf="/projects/fs1/common/genome/lunarc/genomes/mouse/GRCm38.p6/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf",
        TEgtf="/projects/fs3/raquelgg/msc/trim28/data/mm10/repeatmasker/mm10_rmsk_TE.gtf"
    output:
        "../6_TEtranscripts/invitro_crispr/invitro_crispr.cntTable",
        "../6_TEtranscripts/invitro_crispr/invitro_crispr_DESeq2.R"
    shell:
        """
        ml icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132
        ml TEToolkit/2.0.3-Python-2.7.13

        TEtranscripts -t {input.ko} \
        -c {input.ctrl} --GTF {input.gtf} --TE {input.TEgtf} --format BAM --mode multi --project ../6_TEtranscripts/invitro_crispr/invitro_crispr --stranded reverse

        module purge
        """

### CUT AND RUN ###
rule mapping_cut_n_run:
    input:
        "/projects/fs3/raquelgg/cut_n_run/data/raw_samples/fastq_files/mouse/{sample}_R1_001.fastq.gz",
        "/projects/fs3/raquelgg/cut_n_run/data/raw_samples/fastq_files/mouse/{sample}_R2_001.fastq.gz"
    output:
        "../1_mapping/chipseq/cut_n_run/{sample}/{sample}_bt2.sens-loc.sam"
    shell:
        """
        #echo Mapping reads from {wildcards.sample} to mm10!
        ml GCC/7.3.0-2.30  OpenMPI/3.1.1
        ml Bowtie2/2.3.4.2

        bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
        -p 10 -x /projects/fs1/common/genome/lunarc/indicies/bowtie2/mouse/mm10/mm10 -1 {input[0]} -2 {input[1]} -S {output}

        module purge
        """
rule sam_filtering_cut_n_run:
    input:
        "../1_mapping/chipseq/cut_n_run/{sample}/{sample}_bt2.sens-loc.sam"
    output:
        sam="../1_mapping/chipseq/filtered/cut_n_run/{sample}/{sample}_bt2.sens-loc.mapq10.sam",
        statsU="../1_mapping/chipseq/filtered/cut_n_run/{sample}/{sample}_bt2.sens-loc.mapq10.stats",
        statsM="../1_mapping/chipseq/cut_n_run/{sample}/{sample}_bt2.sens-loc.orig.stats"
    shell:
        """
            ml GCC/5.4.0-2.26  OpenMPI/1.10.3
            ml SAMtools/1.4
            echo "> Filter MAPQ >10 to get unique reads:"
            echo ">> samtools view -q 10 -h {input} > {output.sam}"
            samtools view -q 10 -h {input} > {output.sam}
            echo "> Stats unique reads:"
            echo ">> samtools stats {output.sam} > {output.statsU}"
            samtools stats {output.sam} > {output.statsU}
            echo "> Stats all mapped reads:"
            echo ">> samtools stats {input} > {output.statsM}"
            samtools stats {input} > {output.statsM}
        """
rule samToBam_cut_n_run:
    input:
        samU="../1_mapping/chipseq/filtered/cut_n_run/{sample}/{sample}_bt2.sens-loc.mapq10.sam",
        samOrig="../1_mapping/chipseq/cut_n_run/{sample}/{sample}_bt2.sens-loc.sam"
    output:
        bamM="../1_mapping/chipseq/cut_n_run/{sample}/{sample}_bt2.sens-loc.bam",
        bamU="../1_mapping/chipseq/filtered/cut_n_run/{sample}/{sample}_bt2.sens-loc.mapq10.bam"
    shell:
        """
            ml GCC/5.4.0-2.26  OpenMPI/1.10.3
            ml SAMtools/1.4
            echo "> Convert to bam"
            echo "> samtools view -Sb {input.samOrig} > {output.bamM}"
            samtools view -Sb {input.samOrig} > {output.bamM}
            echo "> samtools view -Sb {input.samU} > {output.bamU}"
            samtools view -Sb {input.samU} > {output.bamU}
        """
rule BamToBw_cut_n_run:
    input:
        "../1_mapping/chipseq/filtered/cut_n_run/{sample}/{sample}_bt2.sens-loc.mapq10.bam"
    output:
        "../1_mapping/chipseq/filtered/cut_n_run/{sample}/{sample}_bt2.sens-loc.mapq10.bw"
    shell:
        """
            ml GCC/5.4.0-2.26  OpenMPI/1.10.3
            ml SAMtools/1.4
            ml foss/2016b
            ml Python/3.5.2
            basename=$(basename {input} .bam)
            dirname=$(dirname {input})
            sorted="$dirname/$basename.sorted.bam"
            samtools sort -o $sorted {input}
            samtools index -b $sorted
            bamCoverage --normalizeUsingRPKM -b $sorted -o {output}
        """

rule mapping_adult_chipseq:
    input:
        "/projects/fs3/jakobssonlab/h3k9me3_08.17_JiangY/{sample}.fastq.gz"
    output:
        "../1_mapping/chipseq/h3k9me3/{sample}/{sample}_bt2.sens-loc.sam"
    shell:
        """
        module purge
        mkdir -p ../1_mapping/chipseq/h3k9me3/{wildcards.sample}

        #echo Mapping reads from {wildcards.sample} to mm10!
        ml GCC/7.3.0-2.30  OpenMPI/3.1.1
        ml Bowtie2/2.3.4.2

        bowtie2 --sensitive-local -p 10 -x /projects/fs1/common/genome/lunarc/indicies/bowtie2/mouse/mm10/mm10 -U {input} -S {output}
        module purge
        """
rule sam_filtering_adult_chipseq:
    input:
        "../1_mapping/chipseq/h3k9me3/{sample}/{sample}_bt2.sens-loc.sam"
    output:
        sam="../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.sam",
        statsU="../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.stats",
        statsM="../1_mapping/chipseq/h3k9me3/{sample}/{sample}_bt2.sens-loc.orig.stats"
    shell:
        """
            ml GCC/5.4.0-2.26  OpenMPI/1.10.3
            ml SAMtools/1.4
            echo "> Filter MAPQ >10 to get unique reads:"
            echo ">> samtools view -q 10 -h {input} > {output.sam}"
            samtools view -q 10 -h {input} > {output.sam}
            echo "> Stats unique reads:"
            echo ">> samtools stats {output.sam} > {output.statsU}"
            samtools stats {output.sam} > {output.statsU}
            echo "> Stats all mapped reads:"
            echo ">> samtools stats {input} > {output.statsM}"
            samtools stats {input} > {output.statsM}
        """
rule samToBam_adult_chipseq:
    input:
        samU="../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.sam",
        samOrig="../1_mapping/chipseq/h3k9me3/{sample}/{sample}_bt2.sens-loc.sam"
    output:
        bamM="../1_mapping/chipseq/h3k9me3/{sample}/{sample}_bt2.sens-loc.bam",
        bamU="../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.bam"
    shell:
        """
            ml GCC/5.4.0-2.26  OpenMPI/1.10.3
            ml SAMtools/1.4
            echo "> Convert to bam"
            echo "> samtools view -Sb {input.samOrig} > {output.bamM}"
            samtools view -Sb {input.samOrig} > {output.bamM}
            echo "> samtools view -Sb {input.samU} > {output.bamU}"
            samtools view -Sb {input.samU} > {output.bamU}
        """
rule sorting_adult_chipseq:
    input:
        "../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.bam"
    output:
        "../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.sorted.bam"
    shell:
        """
        ml GCC/7.3.0-2.30
        ml SAMtools/1.9

        samtools sort -o {output} {input}

        module purge
        """
rule BamToBw_adult_chipseq:
    input:
        "../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.bam"
    output:
        "../1_mapping/chipseq/filtered/h3k9me3/{sample}/{sample}_bt2.sens-loc.mapq10.bw"
    shell:
        """
            ml GCC/5.4.0-2.26  OpenMPI/1.10.3
            ml SAMtools/1.4
            ml foss/2016b
            ml Python/3.5.2
            basename=$(basename {input} .bam)
            dirname=$(dirname {input})
            sorted="$dirname/$basename.sorted.bam"
            samtools sort -o $sorted {input}
            samtools index -b $sorted
            bamCoverage --normalizeUsingRPKM -b $sorted -o {output}
        """

rule bamCompare_mean_invitrog3:
    input:
        sample845="../1_mapping/rnaseq/invitro_crispr/ko/g3/845/845Aligned.sortedByCoord.out.bam",
        sample870="../1_mapping/rnaseq/invitro_crispr/ko/g3/870/870Aligned.sortedByCoord.out.bam"
    output:
        g3="../1_mapping/rnaseq/invitro_crispr/ko/g3/invitro_crispr_ko_g3.mean.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.sample845} -b2 {input.sample870} -o {output.g3} -of bigwig --sampleLength 2000 --ratio mean

        module purge
        """
rule bamCompare_mean_invitrog13:
    input:
        sample851="../1_mapping/rnaseq/invitro_crispr/ko/g13/851/851Aligned.sortedByCoord.out.bam",
        sample876="../1_mapping/rnaseq/invitro_crispr/ko/g13/876/876Aligned.sortedByCoord.out.bam"
    output:
        g13="../1_mapping/rnaseq/invitro_crispr/ko/g13/invitro_crispr_ko_g13.mean.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.sample851} -b2 {input.sample876} -o {output.g13} -of bigwig --sampleLength 2000 --ratio mean

        module purge
        """
rule bamCompare_mean_invitrog4:
    input:
        sample848="../1_mapping/rnaseq/invitro_crispr/ko/g4/848/848Aligned.sortedByCoord.out.bam",
        sample873="../1_mapping/rnaseq/invitro_crispr/ko/g4/873/873Aligned.sortedByCoord.out.bam"
    output:
        g4="../1_mapping/rnaseq/invitro_crispr/ko/g4/invitro_crispr_ko_g4.mean.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.sample848} -b2 {input.sample873} -o {output.g4} -of bigwig --sampleLength 2000 --ratio mean

        module purge
        """
rule bamCompare_mean_invitroctrl:
    input:
        sample842="../1_mapping/rnaseq/invitro_crispr/ctrl/ctrl/842/842Aligned.sortedByCoord.out.bam",
        sample867="../1_mapping/rnaseq/invitro_crispr/ctrl/ctrl/867/867Aligned.sortedByCoord.out.bam"
    output:
        ctrl="../1_mapping/rnaseq/invitro_crispr/ctrl/ctrl/invitro_crispr_ctrl.mean.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.sample842} -b2 {input.sample867} -o {output.ctrl} -of bigwig --sampleLength 2000 --ratio mean

        module purge
        """

rule bamCompare_mean_emxko:
    input:
        sample60ctx2="../1_mapping/rnaseq/invivo_bd/ko/60ctx2/60ctx2Aligned.sortedByCoord.out.bam",
        sample73ctx2="../1_mapping/rnaseq/invivo_bd/ko/73ctx2/73ctx2Aligned.sortedByCoord.out.bam"
    output:
        emxko="../1_mapping/rnaseq/invivo_bd/ko/60ctx2_73ctx2.mean.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.sample60ctx2} -b2 {input.sample73ctx2} -o {output.emxko} -of bigwig --sampleLength 2000 --ratio mean

        module purge
        """
rule bamCompare_mean_emxctrl:
    input:
        sample64ctx2="../1_mapping/rnaseq/invivo_bd/ctrl/64ctx2/64ctx2Aligned.sortedByCoord.out.bam",
        sample65ctx2="../1_mapping/rnaseq/invivo_bd/ctrl/65ctx2/65ctx2Aligned.sortedByCoord.out.bam"
    output:
        emxctrl="../1_mapping/rnaseq/invivo_bd/ctrl/64ctx2_65ctx2.mean.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.sample64ctx2} -b2 {input.sample65ctx2} -o {output.emxctrl} -of bigwig --sampleLength 2000 --ratio mean

        module purge
        """

rule bamCompare_mean_cas:
    input:
        g3="../1_mapping/rnaseq/invivo_crispr/ko/cas_g3/cas_g3Aligned.sortedByCoord.out.bw",
        g4="../1_mapping/rnaseq/invivo_crispr/ko/cas_g4/cas_g4Aligned.sortedByCoord.out.bw",
        g13="../1_mapping/rnaseq/invivo_crispr/ko/cas_g13/cas_g13Aligned.sortedByCoord.out.bw",
        g3g4g13="../1_mapping/rnaseq/invivo_crispr/ko/cas_g3g4g13/cas_g3g4g13Aligned.sortedByCoord.out.bw"
    output:
        g3="../1_mapping/rnaseq/invivo_crispr/ko/cas_g3/cas_g3_quater.bw",
        g4="../1_mapping/rnaseq/invivo_crispr/ko/cas_g4/cas_g4Aligned.sortedByCoord.out.bw",
        g13="../1_mapping/rnaseq/invivo_crispr/ko/cas_g13/cas_g13Aligned.sortedByCoord.out.bw",
        g3g4g13="../1_mapping/rnaseq/invivo_crispr/ko/cas_g3g4g13/cas_g3g4g13Aligned.sortedByCoord.out.bw",
        mean="../1_mapping/rnaseq/invivo_crispr/ko/cas_guides_mean.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.g3} -o {output.g3} -of bigwig --sampleLength 2000 --ratio mean

        module purge
        """

rule bamCompare_subtract_cutnrun:
    input:
        signal="../1_mapping/chipseq/filtered/cut_n_run/CR021_S10/CR021_S10_bt2.sens-loc.mapq10.sorted.bam",
        igg="../1_mapping/chipseq/filtered/cut_n_run/CR020_S9/CR020_S9_bt2.sens-loc.mapq10.sorted.bam"
    output:
        "../1_mapping/chipseq/filtered/cut_n_run/CR020_S9_CR021_10_bt2.sens-loc.mapq10.sorted.bw"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        bamCompare -b1 {input.signal} -b2 {input.igg} -o {output} -of bigwig --sampleLength 2000 --ratio subtract

        module purge
        """

rule fullMMERVK10C:
    input:
        fullERVs="../8_TEchipseq/mm10_fullERVs_score300.bed",
        MMERVK10Cint="../8_TEchipseq/MMERVK10C/MMERVK10C-int_mm10.bed" # From TEtranscripts' annotation file
    output:
        unstranded50="/projects/fs3/raquelgg/msc/trim28/8_TEchipseq/MMERVK10C/score300/mm10_fullMMERVK10C_unstranded50.bed"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml BEDTools/2.26.0

        bedtools intersect -u -f 0.5 -a {input.fullERVs} -b {input.MMERVK10Cint} > {output.unstranded50}

        module purge
        """

rule computeMatrix_MMERVK10C_invitro_crispr:
    input:
        cutnrun="../1_mapping/chipseq/filtered/cut_n_run/CR020_S9_CR021_10_bt2.sens-loc.mapq10.sorted.bw",
        invitrog3="../1_mapping/rnaseq/invitro_crispr/ko/g3/invitro_crispr_ko_g3.mean.bw",
        invitrog13="../1_mapping/rnaseq/invitro_crispr/ko/g13/invitro_crispr_ko_g13.mean.bw",
        invitrog4="../1_mapping/rnaseq/invitro_crispr/ko/g4/invitro_crispr_ko_g4.mean.bw",
        invitroctrl="../1_mapping/rnaseq/invitro_crispr/ctrl/ctrl/invitro_crispr_ctrl.mean.bw",
        stranded50="/projects/fs3/raquelgg/msc/trim28/8_TEchipseq/MMERVK10C/score300/mm10_fullMMERVK10C_stranded50.bed"
    output:
        combined_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.mat",
        trim28_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invitro_crispr/mm10_fullMMERVK10C_Trim28_stranded50.means.mat",
        cutnrun_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_stranded50.means.mat"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        computeMatrix scale-regions -S {input.cutnrun} {input.invitroctrl} {input.invitrog3} {input.invitrog4} {input.invitrog13}  -R {input.stranded50}  -b 2000 -a 2000 -out {output.combined_stranded50}
        computeMatrix scale-regions -S {input.invitrog3} {input.invitrog4} {input.invitrog13} {input.invitroctrl} -R {input.stranded50} -b 2000 -a 2000 -out {output.trim28_stranded50}
        computeMatrix scale-regions -S {input.cutnrun}  -R {input.stranded50}  -b 2000 -a 2000 -out {output.cutnrun_stranded50}

        module purge
        """
rule computeMatrix_MMERVK10C_invivo_bd:
    input:
        cutnrun="../1_mapping/chipseq/filtered/cut_n_run/CR045_S2_CR021_S10_bt2.sens-loc.mapq10.norm.bw",
        emxko="../1_mapping/invivo_bd/ko/60ctx2_73ctx2.mean.bw",
        emxctrl="../1_mapping/invivo_bd/ctrl/64ctx2_65ctx2.mean.bw",
        stranded50="/projects/fs3/raquelgg/msc/trim28/8_TEchipseq/MMERVK10C/score300/mm10_fullMMERVK10C_stranded50.bed"
    output:
        combined_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invivo_bd/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.mat"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        computeMatrix scale-regions -S {input.cutnrun} {input.emxctrl} {input.emxko} -R {input.stranded50} -b 2000 -a 2000 -out {output.combined_stranded50}

        module purge
        """
rule computeMatrix_MMERVK10C_invivo_crispr_cas:
    input:
        adultH3K9me3="../1_mapping/chipseq/filtered/h3k9me3/wt/GSM2643045/GSM2643045_bt2.sens-loc.mapq10.norm.bw",
        cas_g3="../1_mapping/invivo_crispr/ko/cas_g3/cas_g3Aligned.sortedByCoord.out.bw",
        cas_g4="../1_mapping/invivo_crispr/ko/cas_g4/cas_g4Aligned.sortedByCoord.out.bw",
        cas_g13="../1_mapping/invivo_crispr/ko/cas_g13/cas_g13Aligned.sortedByCoord.out.bw",
        cas_g3g4g13="../1_mapping/invivo_crispr/ko/cas_g3g4g13/cas_g3g4g13Aligned.sortedByCoord.out.bw",
        cas_ctrl="../1_mapping/invivo_crispr/ctrl/cas_ctrl/cas_ctrlAligned.sortedByCoord.out.bw",
        stranded50="/projects/fs3/raquelgg/msc/trim28/8_TEchipseq/MMERVK10C/score300/mm10_fullMMERVK10C_stranded50.bed"
    output:
        combined_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invivo_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50_cas.means.mat"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        computeMatrix scale-regions -S {input.adultH3K9me3} {input.cas_ctrl} {input.cas_g3} {input.cas_g4} {input.cas_g13} {input.cas_g3g4g13} -R {input.stranded50} -b 2000 -a 2000 -out {output.combined_stranded50}

        module purge
        """

rule plotHeatmap_invitro_crispr_MMERVK10Cint_combined_means:
    input:
        combined_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.mat",
        trim28_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invitro_crispr/mm10_fullMMERVK10C_Trim28_stranded50.means.mat",
        cutnrun_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_stranded50.means.mat"
    output:
        combined_stranded50="../8_TEchipseq/MMERVK10C/plots/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.pdf",
        combined_stranded50_tab="../8_TEchipseq/MMERVK10C/plots/score300/invitro_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.tab"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        plotHeatmap -m {input.combined_stranded50} -o {output.combined_stranded50} --zMin 0 --yMin 0 --yMax 2 \
        --outFileNameMatrix {output.combined_stranded50_tab} --plotFileFormat pdf  \
        --colorMap Greys Greys OrRd OrRd OrRd \
        --regionsLabel "MMERVK10C" \
        --samplesLabel "H3K9me3" "CRISPR Control" "CRISPR g3 KO" "CRISPR g4 KO" "CRISPR g13 KO" \
        --sortUsing max --sortUsingSamples 3

        """
rule plotHeatmap_invivo_bd_MMERVK10Cint_combined_means:
    input:
        combined_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invivo_bd/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.mat"
    output:
        combined_stranded50="../8_TEchipseq/MMERVK10C/plots/score300/invivo_bd/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50.means.pdf"
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        plotHeatmap -m {input.combined_stranded50} -o {output.combined_stranded50} --zMin 0 --yMin 0  --plotFileFormat pdf  \
        --colorMap Greys Greys OrRd \
        --regionsLabel "MMERVK10C" --samplesLabel "H3K9me3" "Emx - Control" "Emx - KO" \
        --sortUsing max --sortUsingSamples 3

        module purge
        """
rule plotHeatmap_invivo_crispr_MMERVK10Cint_combined_means_cas:
    input:
        combined_stranded50="../8_TEchipseq/MMERVK10C/matrices/score300/invivo_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50_cas.means.mat",
    output:
        combined_stranded50="../8_TEchipseq/MMERVK10C/plots/score300/invivo_crispr/mm10_fullMMERVK10C_H3K9me3_Trim28_stranded50_cas.means.pdf",
    shell:
        """
        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml foss/2016b
        ml Python/3.5.2

        plotHeatmap -m {input.combined_stranded50} -o {output.combined_stranded50} --zMin 0 \
        --yMin 0 --yMax 2  --plotFileFormat pdf \
         --colorMap Greys Greys OrRd OrRd OrRd OrRd --regionsLabel "MMERVK10C" --samplesLabel "H3K9me3" "Control" "Cas g3" "Cas g4" "Cas g13" "Cas g3g4g13" --sortUsing max --sortUsingSamples 3

        module purge
        """

rule neighborgene_MMERVK10C:
  input:
      "/projects/fs3/raquelgg/msc/trim28/8_TEchipseq/MMERVK10C/score300/mm10_fullMMERVK10C_unstranded50.bed",
      "/projects/fs3/raquelgg/annotations/mm10/gencode/gencode.vM20.annotation.gene.bed"
  output:
      "../10_neighborgenes/fullMMERVK10C_10kb_windows.bed",
      "../10_neighborgenes/fullMMERVK10C_25kb_windows.bed",
      "../10_neighborgenes/fullMMERVK10C_50kb_windows.bed",
      "../10_neighborgenes/fullMMERVK10C_10kb_windows_intersect_genes.bed",
      "../10_neighborgenes/fullMMERVK10C_25kb_windows_intersect_genes.bed",
      "../10_neighborgenes/fullMMERVK10C_50kb_windows_intersect_genes.bed"
  shell:
      """
      awk '{{if($2-10000 < 1){{print $1,1,$3+10000,$4,$5,$6,$7}}else{{print $1,$2-10000,$3+10000,$4,$5,$6,$7}}}}' OFS="\t" {input[0]} > {output[0]}
      awk '{{if($2-25000 < 1){{print $1,1,$3+25000,$4,$5,$6,$7}}else{{print $1,$2-25000,$3+25000,$4,$5,$6,$7}}}}' OFS="\t"  {input[0]} > {output[1]}
      awk '{{if($2-50000 < 1){{print $1,1,$3+50000,$4,$5,$6,$7}}else{{print $1,$2-50000,$3+50000,$4,$5,$6,$7}}}}' OFS="\t"  {input[0]} > {output[2]}

      ml GCC/5.4.0-2.26  OpenMPI/1.10.3
      ml BEDTools/2.26.0

      bedtools intersect -u -a {input[1]} -b {output[0]} > {output[3]}
      bedtools intersect -u -a {input[1]} -b {output[1]} > {output[4]}
      bedtools intersect -u -a {input[1]} -b {output[2]} > {output[5]}

      module purge
      """
rule neighborgene_ERVK_upreg:
    input:
        "../10_neighborgenes/upregulated/emx_upreg_ERVK.bed",
        "../10_neighborgenes/not_upregulated/emx_notupreg_ERVK.bed",
        "/projects/fs3/raquelgg/annotations/mm10/gencode/gencode.vM20.annotation.gene.bed"
    output:
        "../10_neighborgenes/upregulated/emx_upreg_ERVK_10kb_windows.bed",
        "../10_neighborgenes/upregulated/emx_upreg_ERVK_25kb_windows.bed",
        "../10_neighborgenes/upregulated/emx_upreg_ERVK_50kb_windows.bed",
        "../10_neighborgenes/not_upregulated/emx_notupreg_ERVK_10kb_windows.bed",
        "../10_neighborgenes/not_upregulated/emx_notupreg_ERVK_25kb_windows.bed",
        "../10_neighborgenes/not_upregulated/emx_notupreg_ERVK_50kb_windows.bed",
        "../10_neighborgenes/upregulated/emx_upreg_ERVK_10kb_windows_intersect_genes.bed",
        "../10_neighborgenes/upregulated/emx_upreg_ERVK_25kb_windows_intersect_genes.bed",
        "../10_neighborgenes/upregulated/emx_upreg_ERVK_50kb_windows_intersect_genes.bed",
        "../10_neighborgenes/not_upregulated/emx_notupreg_ERVK_10kb_windows_intersect_genes.bed",
        "../10_neighborgenes/not_upregulated/emx_notupreg_ERVK_25kb_windows_intersect_genes.bed",
        "../10_neighborgenes/not_upregulated/emx_notupreg_ERVK_50kb_windows_intersect_genes.bed"
    shell:
        """
        awk '{{if($2-10000 < 1){{print $1,1,$3+10000,$4,$5,$6,$7}}else{{print $1,$2-10000,$3+10000,$4,$5,$6,$7}}}}' OFS="\t" {input[0]} > {output[0]}
        awk '{{if($2-25000 < 1){{print $1,1,$3+25000,$4,$5,$6,$7}}else{{print $1,$2-25000,$3+25000,$4,$5,$6,$7}}}}' OFS="\t"  {input[0]} > {output[1]}
        awk '{{if($2-50000 < 1){{print $1,1,$3+50000,$4,$5,$6,$7}}else{{print $1,$2-50000,$3+50000,$4,$5,$6,$7}}}}' OFS="\t"  {input[0]} > {output[2]}

        awk '{{if($2-10000 < 1){{print $1,1,$3+10000,$4,$5,$6,$7}}else{{print $1,$2-10000,$3+10000,$4,$5,$6,$7}}}}' OFS="\t"  {input[1]} > {output[3]}
        awk '{{if($2-25000 < 1){{print $1,1,$3+25000,$4,$5,$6,$7}}else{{print $1,$2-25000,$3+25000,$4,$5,$6,$7}}}}' OFS="\t"  {input[1]} > {output[4]}
        awk '{{if($2-50000 < 1){{print $1,1,$3+50000,$4,$5,$6,$7}}else{{print $1,$2-50000,$3+50000,$4,$5,$6,$7}}}}' OFS="\t"  {input[1]} > {output[5]}

        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml BEDTools/2.26.0

        bedtools intersect -u -a {input[2]} -b {output[0]} > {output[6]}
        bedtools intersect -u -a {input[2]} -b {output[1]} > {output[7]}
        bedtools intersect -u -a {input[2]} -b {output[2]} > {output[8]}

        bedtools intersect -u -a {input[2]} -b {output[3]} > {output[9]}
        bedtools intersect -u -a {input[2]} -b {output[4]} > {output[10]}
        bedtools intersect -u -a {input[2]} -b {output[5]} > {output[11]}

        module purge
        """
rule neighborgene_ERVK:
    input:
        "../../../annotations/mm10/repeatmasker/mm10_rmsk_TE.bed",
        "/projects/fs3/raquelgg/annotations/mm10/gencode/gencode.vM20.annotation.gene.bed"
    output:
        "../10_neighborgenes/ERVK_10kb_windows.bed",
        "../10_neighborgenes/ERVK_25kb_windows.bed",
        "../10_neighborgenes/ERVK_50kb_windows.bed",
        "../10_neighborgenes/ERVK_10kb_windows_intersect_genes.bed",
        "../10_neighborgenes/ERVK_25kb_windows_intersect_genes.bed",
        "../10_neighborgenes/ERVK_50kb_windows_intersect_genes.bed"
    shell:
        """
        grep -w ERVK {input[0]} | awk '{{if($2-10000 < 1){{print $1,1,$3+10000,$4,$5,$6,$7}}else{{print $1,$2-10000,$3+10000,$4,$5,$6,$7}}}}' OFS="\t" > {output[0]}
        grep -w ERVK {input[0]} | awk '{{if($2-25000 < 1){{print $1,1,$3+25000,$4,$5,$6,$7}}else{{print $1,$2-25000,$3+25000,$4,$5,$6,$7}}}}' OFS="\t" > {output[1]}
        grep -w ERVK {input[0]} | awk '{{if($2-50000 < 1){{print $1,1,$3+50000,$4,$5,$6,$7}}else{{print $1,$2-50000,$3+50000,$4,$5,$6,$7}}}}' OFS="\t" > {output[2]}

        ml GCC/5.4.0-2.26  OpenMPI/1.10.3
        ml BEDTools/2.26.0

        bedtools intersect -u -a {input[1]} -b {output[0]} > {output[3]}
        bedtools intersect -u -a {input[1]} -b {output[1]} > {output[4]}
        bedtools intersect -u -a {input[1]} -b {output[2]} > {output[5]}

        module purge
        """
