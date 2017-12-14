#!/complgen/bin/anaconda/bin/python2.7
# Script to check in WES whether the gender as determined from the bam
# (by counting mappings to chrY normalized to chrX) matches the one specified in gentli
# wdecoster

import sys
import pysam
import os
import concurrent.futures as cfutures
import glob


class Wes_sample(object):
    def __init__(self, bam, gender):
        self.bam = bam
        self.gender = gender


def main():
    bams = [os.path.realpath(path) for path in glob.glob(sys.argv[1] + '/*.bam')]
    with cfutures.ProcessPoolExecutor() as executor:
        for i in executor.map(get_gender, bams):
            print("{}\t{}".format(i.bam, i.gender))


def valid_read(read):
    if read.mapping_quality >= 5 and read.reference_end and read.reference_start is not None:
        return True
    else:
        return False


def get_gender(bam):
    '''Determine the gender of a bam file.

    Based on the reads mapping between but not in the PAR regions
    of the Y chromosome normalized to the counts on chromosome X'''
    workfile = pysam.AlignmentFile(bam, "rb")
    yreadcount = sum([valid_read(read) for read in workfile.fetch(region='chrY:2781479-56887902')])
    refreadcount = sum([valid_read(read) for read in workfile.fetch(region='chrX')])
    countratio = yreadcount / float(refreadcount)
    if countratio <= 0.03:
        return Wes_sample(bam, 'f')
    elif 0.03 < countratio < 0.09:
        return Wes_sample(bam, 'u')
    else:
        return Wes_sample(bam, 'm')


if __name__ == "__main__":
    main()
