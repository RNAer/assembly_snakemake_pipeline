import tempfile
import os
import glob
configfile: "config.yaml"
#IN = /home/snurk/pipeline_test/in
#OUT = /home/snurk/pipeline_test/out

DATA_ROOT = config["DATA_ROOT"]
PREPS = config["PREPS"]
SAMPLES = config["SAMPLES"]
SIZES = config['SIZES']
TMP_DIR_ROOT = config["TMP_DIR_ROOT"]

localrules: all, gather_assemblies

shell.prefix('source deactivate; ')

rule all:
    input:
        expand("quality/{sample}/{size}/metaquast/done.txt", sample=SAMPLES, size=SIZES),
        expand("quality/{sample}/{size}/quast/report.html", sample=SAMPLES, size=SIZES)

# combine the multi-lanes of the reads
rule combine:
    input:
        DATA_ROOT + '/{sample}/PE150_HiSeq4000/{prep}/{size}',
    output:
        R1 = temp('{sample}/PE150_HiSeq4000/{prep}/{size}/R1.fastq.gz'),
        R2 = temp('{sample}/PE150_HiSeq4000/{prep}/{size}/R2.fastq.gz')
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            tmp_r1 = temp_dir + '/R1.fastq.gz'
            tmp_r2 = temp_dir + '/R2.fastq.gz'
            r1 = ' '.join(sorted(glob.glob(input[0] + '/*R1*.fastq.gz')))
            print('zcat %s | gzip -c > %s && cp %s {output.R1}' % (r1, tmp_r1, tmp_r1))
            shell('zcat %s | gzip -c > %s && cp %s {output.R1}' % (r1, tmp_r1, tmp_r1))
            r2 = ' '.join(sorted(glob.glob(input[0] + '/*R2*.fastq.gz')))
            print('zcat %s | gzip -c > %s && cp %s {output.R2}' % (r2, tmp_r2, tmp_r2))
            shell('zcat %s | gzip -c > %s && cp %s {output.R2}' % (r2, tmp_r2, tmp_r2))

rule subsample:
    input:
        "{sample}/PE150_HiSeq4000/{prep}/{size}/{R}.fastq.gz"
    output:
        "{sample}/PE150_HiSeq4000/{prep}/{size}/{R}_%s.fastq.gz" % config['SUB_READ_CNT']
    log:
        'log/{sample}/{prep}/{size}/{R}_subsample.log'
    message: "subsampling {input} -> {output}"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            tmp_out = temp_dir + "/sub.fastq.gz"
            shell("%s sample -s 239 -2 {input} %s 2> {log} | gzip -c > %s && cp %s {output}"
                % (config["SEQTK"], config["SUB_READ_CNT"], tmp_out, tmp_out))

rule clean:
    input:
        R1="{sample}/PE150_HiSeq4000/{prep}/{size}/R1_%s.fastq.gz" % config['SUB_READ_CNT'],
        R2="{sample}/PE150_HiSeq4000/{prep}/{size}/R2_%s.fastq.gz" % config['SUB_READ_CNT']
    output:
        R1="{sample}/PE150_HiSeq4000/{prep}/{size}/R1_%s_clean.fastq.gz" % config['SUB_READ_CNT'],
        R2="{sample}/PE150_HiSeq4000/{prep}/{size}/R2_%s_clean.fastq.gz" % config['SUB_READ_CNT']
    log:
        'log/{sample}/{prep}/{size}/clean.log'
    message: "cleaning {input} -> {output}"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            tmp_R1 = temp_dir + "/R1.fastq.gz"
            tmp_R2 = temp_dir + "/R2.fastq.gz"
            shell("%s -m 36 -o %s -p %s {input.R1} {input.R2} &> {log} && cp %s {output.R1} && cp %s {output.R2}"
                % (config["CUTADAPT_PRX"], tmp_R1, tmp_R2, tmp_R1, tmp_R2))

rule metaspades:
    input:
        R1="{sample}/PE150_HiSeq4000/{prep}/{size}/R1_%s_clean.fastq.gz" % config['SUB_READ_CNT'],
        R2="{sample}/PE150_HiSeq4000/{prep}/{size}/R2_%s_clean.fastq.gz" % config['SUB_READ_CNT']
    output:
        "{sample}/PE150_HiSeq4000/{prep}/{size}/metaspades/scaffolds.fasta"
    log:
        'log/{sample}/{prep}/{size}/metaspades.log'
    message: "metaspades: {input} -> {output}"
    params:
        threads=12,
        mem=120,
        spades=config["SPADES"]
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("{params.spades} --meta -t {params.threads} -m {params.mem} -1 {input.R1} -2 {input.R2} -o %s &> {log} && "
            " cp -r %s/{{spades.log,*.fasta,corrected}} {wildcards.sample}/PE150_HiSeq4000/{wildcards.prep}/{wildcards.size}/metaspades/" % (temp_dir, temp_dir))
        #shell("{params.spades}/spades.py --meta -t {params.threads} -m {params.mem} -1 {input.R1} -2 {input.R2} -o {wildcards.path}/metaspades &> {log}")

rule gather_assemblies:
    input:
        "{sample}/PE150_HiSeq4000/{prep}/{size}/metaspades/scaffolds.fasta"
    output:
        "quality/{sample}/{size}/{prep}.fasta"
    message: "gather_assemblies: {input} -> {output}"
    shell:
        "rm -rf {output} && ln -s $(readlink -e {input}) {output}"
        #"ln -s $(readlink -e {input}) {output}"

rule quast:
    input:
        assemblies=expand("quality/{{sample}}/{{size}}/{prep}.fasta", prep=PREPS)
    output:
        "quality/{sample}/{size}/quast/report.html"
    log:
        "log/{sample}/{size}/quast.log"
    params:
        threads=8,
        quast_dir=config["QUAST_DIR"]
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("{params.quast_dir}/quast.py -t {params.threads} -o %s {input.assemblies} &> {log} &&"
                  " cp -r %s/* quality/{wildcards.sample}/{wildcards.size}/quast" % (temp_dir, temp_dir))

rule metaquast:
    input:
        assemblies=expand("quality/{{sample}}/{{size}}/{prep}.fasta", prep=PREPS)
    output:
        "quality/{sample}/{size}/metaquast/done.txt"
    log:
        "log/{sample}/{size}/metaquast.log"
    params:
        refs="refs/{sample}/",
        threads=8,
        quast_dir=config["QUAST_DIR"]
    run:
        #shell("{params.metaquast} -t {params.threads} -o quality/{wildcards.sample}/metaquast -R {input.refs} {input.assemblies} &> {log}")
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            ref = ','.join(glob.glob(params.refs + '*.fa'))
            # shell("if [ -d {params.refs} ] ; then {params.quast_dir}/metaquast.py -t {params.threads} -o %s -R %s {input.assemblies} &> {log} &&"
            #       " cp -r %s/* quality/{wildcards.sample}/metaquast ; fi; touch {output}" % (temp_dir, ref, temp_dir))
            shell("{params.quast_dir}/metaquast.py --gene-finding --meta -t {params.threads} -o %s -R %s {input.assemblies} &> {log} &&"
                  " cp -r %s/* quality/{wildcards.sample}/{wildcards.size}/metaquast && touch {output}" % (temp_dir, ref, temp_dir))
