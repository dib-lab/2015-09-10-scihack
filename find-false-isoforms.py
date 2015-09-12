#!/usr/bin/env python

import argparse
import os
import sys

from khmer import khmer_args
import screed

import sbt

def iterkmers(seq, K):
    for start in xrange(len(seq) - K + 1):
        yield seq[start:start+K]

def search_sequence(node, seq, threshold):
    presence = [node.graph.get(kmer) for kmer in iterkmers(seq, node.graph.ksize())]
    if sum(presence) >= int(threshold * (len(seq) - node.graph.ksize() + 1)):
        return 1
    return 0

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--sbt')
    parser.add_argument('--assembly', nargs='+')
    parser.add_argument('--threshold', type=float, default=0.9)
    parser.add_argument('--print-tree', action='store_true', default=False)
    args = parser.parse_args()

    # Load the Sequence Bloom Tree
    tree = sbt.load_sbt(args.sbt_file)

    for fn in args.assembly:
        with open(fn + '.good', 'wb') as good_fp, open(fn + '.bad', 'wb') as bad_fp:
            print('*** Searching', fn)
            n_good = 0
            n_bad = 0
            for record in screed.open(fn):
                results = tree.find(search_sequence, record.sequence, args.threshold)
                if results:
                    n_good += 1
                    good_fp.write('>{name}:[{support}]\n{seq}'.format(
                        name=record.name, support=','.join([x.name for x in results]),
                        seq=record.sequence))
                else:
                    n_bad += 1
                    bad_fp.write('>{name}\n{seq}'.format(name=record.name, seq=record.sequence))
            print(n_good, 'good', n_bad, 'bad transcripts')

    if args.print_tree:
        sbt.print_sbt(tree)

if __name__ == '__main__':
    main()
