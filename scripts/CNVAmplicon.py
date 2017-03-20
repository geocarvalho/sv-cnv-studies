from subprocess import call
import os
import sys

def quandico_call(path_bam, path_file, bam_file, sample_name):

def main():
    path_bam = os.path.abspath(sys.argv[1])
    path_file, bam_file = path_bam.rsplit("/", 1)
    path_file = path_file.rsplit("/", 1)[0]
    sample_name = bam_file.split("_", 1)[0]
    analysis = sys.argv[2]
