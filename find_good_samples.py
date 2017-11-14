import gentlisecret
import fdb
import os
import sys
import glob


class WES_experiment(object):
    def __init__(self, sid, individual, sample_id, experiment, description, project):
        self.sid = sid
        self.individual = individual
        self.sample_id = sample_id
        self.experiment = experiment
        self.description = description
        self.project = project
        base = "/storage/macchiato/complgen3/results/" + project + "/exomes/hg19s/"
        exp_folder = "exp" + str(experiment) + "-hg19s-" + description.replace(' ', '_') + "/"
        indiv_folder = "*_" + str(sid) + "_hg19s/"
        bam = "map-rdsbwa-*_" + str(sid) + "_hg19s.bam"
        bampath = glob.glob(base + exp_folder + "samples/" + indiv_folder + bam)
        if len(bampath) == 1:
            self.bampath = bampath[0]
            self.bamname = os.path.basename(self.bampath)
        else:
            sys.exit("Path to bam of {} incorrect:\n{}".format(individual, bampath))

    def make_link(self, destdir):
        os.symlink(self.bampath, os.path.join(destdir, self.bamname))
        os.symlink(self.bampath + '.bai', os.path.join(destdir, self.bamname + '.bai'))


def find_samples(kit, min_cov_frac):
    con = fdb.connect(dsn='molgenvz.cde.ua.ac.be:/home/firebird/gentli.fdb',
                      user='wdc',
                      password=gentlisecret.pw,
                      role='NBD_SC')
    cur = con.cursor()
    cur.execute('select "id", "sample_individual", "sample_sample", "experiment_id", \
               "experiment_description", "wes_project" from "NBD_SC:full_wes" \
               where \
               "wes_exome_capture_kit" = ? and \
               "ngs_seq_pct_target20x" > ?', (kit, min_cov_frac))
    return [WES_experiment(*i) for i in cur.fetchall()]


def main():
    for sample in find_samples("seqcapv3", 75):
        sample.make_link('bams')


if __name__ == '__main__':
    main()
