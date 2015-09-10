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
    parser.add_argument('--sbt-file')
    parser.add_argument('--search', nargs='+')
    parser.add_argument('--threshold', type=float, default=0.9)
    parser.add_argument('--print-tree', action='store_true', default=False)
    args = parser.parse_args()

    tree = sbt.load_sbt(args.sbt_file)

    if args.search:
        for fn in args.search:
            print '*** Searching', fn
            for record in screed.open(fn):
                print '---\n', record.name
                print [x.metadata for x in tree.find(search_sequence, record.sequence, args.threshold)]

    if args.print_tree:
        sbt.print_sbt(tree)

if __name__ == '__main__':
    main()
