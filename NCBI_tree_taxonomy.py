#!/usr/bin/env python3
import argparse
import re
import sys
import os.path
import sqlite3
from Bio import Phylo
from Bio import Entrez


def assembly_info(assembly):
    result = dict()
    assembly_id = Entrez.read(Entrez.esearch(db='assembly', term=assembly))['IdList']
    if len(assembly_id) == 0:
        return None
    summary = Entrez.read(
        Entrez.esummary(db='assembly', id=assembly_id[0]), validate=False
    )['DocumentSummarySet']['DocumentSummary'][0]
    '''Summary example:
    {'DocumentSummarySet': DictElement({'DocumentSummary': [DictElement({'RsUid': '255868'
    'Taxid': '456320', 
    'Organism': 'Methanococcus voltae A3 (euryarchaeotes)',
    'SpeciesTaxid': '2188',
    'SpeciesName': 'Methanococcus voltae'
    '''

    if summary:
        result['id'] = assembly_id[0]
        result['SpeciesTaxid'] = summary['SpeciesTaxid']
        result['SpeciesName'] = summary['SpeciesName']
        result['Taxid'] = summary['Taxid']
        result['Organism'] = summary['Organism']
        return result
    else:
        print(f'Error with {assembly}')
        return None


def main():
    args = parse_args()
    Entrez.email = 'email@domain.edu'
    conn = sqlite3.connect(args.database)
    c = conn.cursor()
    create_cmd = '''
        CREATE TABLE IF NOT EXISTS assemblies
        (acc text, id text, taxid text, organism text,
        speciestaxid text, species text)
    '''
    conn.execute(create_cmd)

    insert_cmd = '''
        INSERT INTO assemblies (acc, id, taxid, organism, speciestaxid, species)
        VALUES (?, ?, ?, ?, ?, ?)
    '''
    
    regex = re.compile(r'(GC[AF]_\d+\.\d)')
    tree = Phylo.read(args.input, 'newick')

    for leaf in tree.get_terminals():
        leaf_spec = ''
        m = regex.search(leaf.name)
        if m:
            leaf_acc = m.group()
        else:
            continue
        #leaf_acc = regex.search(leaf.name).group(0)

        # Check if accession is already in the DB
        c.execute('SELECT EXISTS(SELECT 1 FROM assemblies where acc=?)', (leaf_acc,))
        if c.fetchone()[0] == 1:
            res = c.execute('SELECT * FROM assemblies where acc=?', (leaf_acc,))
            tax_array = res.fetchone()
            leaf_spec = tax_array[5]
            if args.print:
                print('\t'.join([tax_array[0], tax_array[4], tax_array[5]]))
        else:
            r = assembly_info(leaf_acc)
            if r:
                if args.print:
                    print('\t'.join([
                        leaf_acc,
                        r['SpeciesTaxid'],
                        r['SpeciesName']]))
                c.execute(insert_cmd, (
                    leaf_acc, r['id'], r['Taxid'], r['Organism'],
                    r['SpeciesTaxid'], r['SpeciesName']
                ))
                conn.commit()
            leaf_spec = r['SpeciesName']

        if len(leaf_spec) > 0:
            leaf.name = f'{leaf_spec} {leaf_acc}'

    Phylo.write(tree, args.output, 'newick')


def parse_args():
    parser = argparse.ArgumentParser(
        description='Find taxonomy of NCBI Assembly accessions in phylogenetic trees.'
        )
    parser.add_argument('-d', '--database', required=True,
        help='SQLite3 database to store the data.')    
    parser.add_argument('-i', '--input',
        help='Input tree.')
    parser.add_argument('-o', '--output',
        help='Output tree.')
    parser.add_argument('-p', '--print', action='store_true',
        help='Print found taxa to stdout.')
    
    return parser.parse_args()


if __name__ == '__main__':
    main()
