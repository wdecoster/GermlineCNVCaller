import gentlisecret
import fdb
import sys


def search_gender(indiv):
    sex_genotype = {"f": "chrXchrX", "m": "chrXchrY"}
    con = fdb.connect(dsn='molgenvz.cde.ua.ac.be:/home/firebird/gentli.fdb',
                      user='wdc',
                      password=gentlisecret.pw,
                      role='NBD_SC')
    cur = con.cursor()
    cur.execute('select "id", "gender" from "NBD_SC:individual" \
               where "id" = ?', (indiv, ))
    try:
        return sex_genotype[cur.fetchall()[0][1]]
    except KeyError:
        sys.stderr.write("No gender found for {}\n".format(indiv))


def main():

    print("SAMPLE_NAME\tSEX_GENOTYPE")
    for sample in sys.stdin.readlines():
        try:
            print("{}\t{}".format(sample.strip(), search_gender(sample.split('_')[0].lower())))
        except IndexError:
            sys.stderr.write("Problem parsing {}\n".format(sample.strip()))


if __name__ == '__main__':
    main()
