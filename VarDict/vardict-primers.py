#coding = utf-8
import os

pathFile = os.listdir(os.path.dirname(os.path.realpath(__file__)))

for arq in pathFile:
    if arq.startswith("primer"):
        pipeline = arq.split(".")[0].split("-")[-1]
        with open("vardict-%s-primers.bed" %pipeline, "w") as vardict_bed:
            with open(arq, "r") as primers:
                for line in primers:
                    if line.startswith("Customer"):
                        continue
                    else:
                        try:
                            chromo = line.strip().split("\t")[1]
                            primer_start = line.strip().split("\t")[2]
                            primer_end = line.strip().split("\t")[3]
                            gene_info = pipeline
                            score = "0"
                            strand = "."
                            exon_start = int(primer_start) + len(line.strip().split("\t")[4])
                            exon_end = int(primer_end) - len(line.strip().split("\t")[5])
                            vardict_bed.write("chr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                            chromo, primer_start, primer_end, gene_info, score, strand,
                            exon_start, exon_end))
                        except:
                            print line.strip().split("\t")
