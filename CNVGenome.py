from subprocess import call
import sys
import os

def lumpy_analysis(path_bam, path_file, bam_file, sample_name):
    '''
    Lumpy pipeline to call SV, according to Bcbio is great to detect variants between
    1-250 bp, also with results for 1000-25000 bp
    '''

    print "**** Lumpy analysis ****"
    # Making a directory to store the lumpy's archives
    call("mkdir %s/lumpy" % path_file, shell=True)
    print "-Lumpy file created"

    # Extracting the discordant paired-end alignments
    call("samtools view -b -F 1294 %s > %s/lumpy/%s.discordants.unsorted.bam" % \
    (path_bam, path_file, sample_name), shell=True)
    print "-Discordant alignments extracted"

    # Extract the split-read alignments
    call("samtools view -h %s | /home/bioinfo/lumpy-sv/scripts/extractSplitReads_BwaMem \
    -i stdin | samtools view -Sb > %s/lumpy/%s.splitters.unsorted.bam" % \
    (path_bam, path_file, sample_name), shell=True)
    print "-Split-read alignments extracted"

    # Sort both alignments
    call("samtools sort {path}/lumpy/{sample}.discordants.unsorted.bam -o \
    {path}/lumpy/{sample}.discordants.bam".format(path=path_file, sample=sample_name), shell=True)

    call("samtools sort {path}/lumpy/{sample}.splitters.unsorted.bam -o \
    {path}/lumpy/{sample}.splitters.bam".format(path=path_file, sample=sample_name), shell=True)
    print "-Alignments sorted"

    # Running LUMPY express
    print "lumpyexpress -v -B {bam} -S {path}/lumpy/{sample}.splitters.bam -D \
    {path}/lumpy/{sample}.discordants.bam -T {path}/lumpy/tmp -o \
    {path}/lumpy/{sample}.bam.vcf".format(bam=path_bam, path=path_file, sample=sample_name)
    call("lumpyexpress -v -B {bam} -S {path}/lumpy/{sample}.splitters.bam -D \
    {path}/lumpy/{sample}.discordants.bam -T {path}/lumpy/tmp -o \
    {path}/lumpy/{sample}.bam.vcf".format(bam=path_bam, path=path_file, sample=sample_name), shell=True)

    # Moving the vcf result to the file
    call("mv {bam_file}.vcf {path}/lumpy/{sample}".format(bam_file=bam_file, \
    path=path_file, sample=sample_name), shell=True)

    print "**** LUMPY express runned, bye! ****"

def delly_analysis(path_bam, path_file, bam_file, sample_name):
    '''
    Delly pipeline to detect SV variants, according to Bcbio is better to call
    variants between 250 bp and 1000 bp, also with good results for 1000-25000 bp.
    '''
    print "**** Delly analyse ****"
    # Making a directory to store the lumpy's archives
    call("mkdir %s/delly" % path_file, shell=True)
    print "-Delly file created"

    for process in ["DEL", "DUP", "INV", "TRA", "INS"]:
        print "-Starting delly analyse, type: %s" % process
        #germline SV calling
        call("delly call -t {process} -x /home/bioinfo/nextgen-pipeline/hg19.excl -g \
        /home/bioinfo/nextgen-pipeline/ucsc.hg19.fasta -o {path}/delly/{sample}.{process}.bcf \
        {path_bam}".format(process=process, path=path_file,sample=sample_name, path_bam=path_bam), shell=True)
        #bcf archive to vcf archive
        print "Bcftools covertion to VCF"
        call("bcftools view {path}/delly/{sample}.{process}.bcf > \
        {path}/delly/{sample}.{process}.vcf".format(path=path_file, sample=sample_name, process=process), \
        shell=True)
        print "bcftools view {path}/delly/{sample}.{process}.bcf > \
        {path}/delly/{sample}.{process}.vcf".format(path=path_file, sample=sample_name, process=process)

def manta_analysis(path_bam, path_file, bam_file, sample_name):
    '''
    Mant pipeline to call SV and indels
    '''
    print "**** Manta analysis ****"

    call("mkdir %s/manta" % path_file, shell=True)
    print "-Manta file created"

    #tumor-only analysis
    call("configManta.py --tumorBam %s --referenceFasta /home/bioinfo/nextgen-pipeline/ucsc.hg19.fasta \
    --runDir %s/manta" % (path_bam, path_file), shell=True)
    #execution on a single node
    call("python %s/manta/runWorkflow.py -m local -j 8" % path_file, shell=True)

def main():
    path_bam = os.path.abspath(sys.argv[1])
    path_file, bam_file = path_bam.rsplit("/", 1)
    path_file = path_file.rsplit("/", 1)[0]
    sample_name = bam_file.split("_", 1)[0]
    analysis = sys.argv[2]
    try:
        if analysis == "lumpy":
            print "--------Just running Lumpy pipeline--------"
            lumpy_analysis(path_bam, path_file, bam_file, sample_name)
        elif analysis == "delly":
            print "--------Just running Delly pipeline--------"
            delly_analysis(path_bam, path_file, bam_file, sample_name)
        elif analysis == "manta":
            print "--------Just running Manta pipeline--------"
            manta_analysis(path_bam, path_file, bam_file, sample_name)
        elif analysis == "all":
            print "--------Running all CNV pipelines--------"
            lumpy_analysis(path_bam, path_file, bam_file, sample_name)
            delly_analysis(path_bam, path_file, bam_file, sample_name)
            manta_analysis(path_bam, path_file, bam_file, sample_name)
    except:
        print "Try the parameters: lumpy, delly, manta and all"

main()
