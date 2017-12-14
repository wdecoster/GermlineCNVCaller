#!/complgen/bin/anaconda/bin/python2.7
# Script to check in WES whether the gender as determined from the bam
# (by counting mappings to chrY normalized to chrX) matches the one specified in gentli
# wdecoster

import sys
import pysam
import os
import concurrent.futures as cfutures
import glob
import matplotlib.pyplot as plt


class Wes_sample(object):
    def __init__(self, bam, yreads, xreads, totalreads):
        self.bam = bam
        self.name = os.path.basename(bam).split('_')[0].lower().replace('map-rdsbwa-', '')
        self.yratio = yreads / totalreads
        self.xratio = xreads / totalreads
        self.gender = self.infer_gender()

    def infer_gender(self):
        self.countratio = self.yratio / self.xratio
        if self.countratio <= 0.03:
            return 'f'
        elif 0.03 < self.countratio < 0.09:
            return 'u'
        else:
            return 'm'


def main():
    bams = [os.path.realpath(path) for path in glob.glob(sys.argv[1] + '/*.bam')]
    with cfutures.ProcessPoolExecutor() as executor:
        res = list(executor.map(get_gender, bams))
        for i in res:
            print("{}\t{}".format(os.path.basename(i.bam).replace('.bam', ''), i.gender))
        plot_genders(res)


def valid_read(read):
    """Check if a read is properly mapped."""
    if read.mapping_quality >= 5 and read.reference_end and read.reference_start is not None:
        return True
    else:
        return False


def get_gender(bam):
    '''Determine the gender of a bam file.

    Based on the reads mapping between but not in the PAR regions
    of the Y chromosome normalized to the counts on chromosome X'''
    workfile = pysam.AlignmentFile(bam, "rb")
    return Wes_sample(
        bam=bam,
        yreads=sum([valid_read(read) for read in workfile.fetch(region='chrY:2649520-59034050')]),
        xreads=sum([valid_read(read) for read in workfile.fetch(region='chrX:2699520-154931044')]),
        totalreads=workfile.mapped)


def plot_genders(res):
    c_dict = {'f': 'pink', 'm': 'blue', 'u': 'grey'}
    plt.scatter(
        x=[i.xratio for i in res],
        y=[i.yratio for i in res],
        c=[c_dict[i.gender] for i in res],
        s=2)
    xmax = max([i.xratio for i in res])
    ymax = max([i.yratio for i in res])
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
    for i in res:
        if i.gender == "u":
            plt.annotate(i.name, xy=(i.xratio, i.yratio), size=6)
    plt.xlabel("Normalised read count on chrX")
    plt.ylabel("Normalised read count on chrY")
    plt.savefig("gender_plot.png")


if __name__ == "__main__":
    main()
