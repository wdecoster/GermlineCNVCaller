import gentlisecret
import fdb
import sys


def search_gender(indiv):
    sex_genotype = {"f": "XX", "m": "XY"}
    con = fdb.connect(dsn='molgenvz.cde.ua.ac.be:/home/firebird/gentli.fdb',
                      user='wdc',
                      password=gentlisecret.pw,
                      role='NBD_SC')
    cur = con.cursor()
    cur.execute('select "id", "gender" from "NBD_SC:individual" \
               where "id" = ?', (indiv, ))
    return sex_genotype[cur.fetchall()[0][1]]


def main():

    print("SAMPLE_NAME\tSEX_GENOTYPE")
    for sample in sys.stdin.readlines():
        print("{}\t{}".format(sample.strip(), search_gender(sample.split('_')[0].lower())))


if __name__ == '__main__':
    main()