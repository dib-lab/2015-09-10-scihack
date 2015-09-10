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

class NodegraphFactory(object):

    def __init__(self, args):
        self.args = args

    def create_nodegraph(self):
        return khmer_args.create_nodegraph(self.args)

def main():

    parser = khmer_args.build_nodegraph_args()
    parser.add_argument('--samples', nargs='+')
    parser.add_argument('--save-prefix')
    parser.add_argument('--print-tree', action='store_true', default=False)
    args = parser.parse_args()

    factory = NodegraphFactory(args)
    root = sbt.Node(factory)

    for sample_fn in args.samples:
        print '*** Build node for', sample_fn
        leaf = sbt.Leaf(sample_fn, factory.create_nodegraph())
        print '--- Consuming file...'
        leaf.graph.consume_fasta(sample_fn)
        print '--- Adding node to SBT...'
        root.add_node(leaf)
        print '--- Done with', sample_fn

    print '\n*** Saving to disk'
    fn = sbt.save_sbt(root, args.save_prefix)
    print '--- Save to', fn

    if args.print_tree:
        sbt.print_sbt(root)

if __name__ == '__main__':
    main()
