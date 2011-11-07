__author__ = 'will'
import argparse
import csv
import urllib
import multiprocessing
import re
from operator import itemgetter
import sys
from itertools import izip, repeat, tee, imap
import json


def make_picr_request(intup):
    accession, output_id = intup
    url = 'http://www.ebi.ac.uk/Tools/picr/rest/'
    url += 'getUPIForAccession?accession=%(acc)s&database=%(db)s'
    return urllib.urlopen(url % {'acc':accession, 'db':output_id}).read()


def extract_picr_ids(xmlin):
    return re.findall('\<accession\>(.*?)\</accession\>', xmlin)


def check_picr_database(database):

    available = set(x.strip() for x in open('picr_databases.txt'))
    return database in available

def convert_using_picr(id_iter, output_id, num_processes=4):

    assert check_picr_database(output_id)

    pool = multiprocessing.Pool(num_processes)

    input_iter, output_iter = tee(id_iter, 2)
    iterable = izip(input_iter,repeat(output_id))
    res_iter = pool.imap(make_picr_request, iterable, chunksize = 10)
    outconv = []

    for xmlout, input_id in izip(res_iter, output_iter):
        outconv.append((input_id, set(extract_picr_ids(xmlout))))

    return outconv


def make_mygene_request(accession):
    #print 'http://mygene.info/query?q='+accession
    return urllib.urlopen('http://mygene.info/query?q='+accession).read()

def extract_mygene_ids(injson):
    #print injson
    obj = json.loads(injson)
    return set(row['id'] for row in obj['rows'] if row['id'].isdigit())

def convert_using_mygene(id_iter, num_processes = 4):

    pool = multiprocessing.Pool(num_processes)

    input_iter, output_iter = tee(id_iter, 2)
    res_iter = pool.imap(make_mygene_request, input_iter, chunksize = 10)
    outconv = []

    for xmlout, input_id in izip(res_iter, output_iter):
        outconv.append((input_id, set(extract_mygene_ids(xmlout))))

    return outconv





if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Convert IDS.')
    parser.add_argument('-o','--output-DB', type = str,
                        default = 'SWISSPROT', dest = 'database',
                        help = 'The dabase to convert too.')
    parser.add_argument('-n', type = int, default = 4, dest = 'numprocesses',
                        help = 'Number of processes to bombard the webservice with. DONT BE A DICK!')

    args = parser.parse_args()
    indata = imap(itemgetter(0), csv.reader(sys.stdin, delimiter = '\t'))

    if args.database.lower() == 'entrez':
        output = convert_using_mygene(indata, num_processes = args.numprocesses)
    elif check_picr_database(args.database):
        output = convert_using_picr(indata, args.database, num_processes = args.numprocesses)

    else:
        raise KeyError, 'Could not understand the output database: %s' % args.database

    for idval, newids in output:
        print '%s\t%s' % (idval, '|'.join(newids))