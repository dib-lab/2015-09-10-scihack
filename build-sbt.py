#!/usr/bin/env python

from __future__ import print_function

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

    if not args.save_prefix:
        print('No save prefix specified! Exiting...', file=sys.stderr)
        sys.exit(1)

    factory = NodegraphFactory(args)
    root = sbt.Node(factory)

    for sample_fn in args.samples:
        print('*** Build node for', sample_fn)
        leaf = sbt.Leaf(os.path.basename(sample_fn),
                        os.path.basename(sample_fn),
                        factory.create_nodegraph())
        fname = os.path.join('.sbt.' + args.save_prefix,
                             ".".join([args.save_prefix, os.path.basename(sample_fn), 'sbt']))
        if os.path.exists(fname):
            print('--- Loading existing file...')
            leaf.graph.load(fname)
        else:
            print('--- Consuming file...')
            leaf.graph.consume_fasta(sample_fn)
        print('--- Adding node to SBT...')
        root.add_node(leaf)
        print('--- Done with', sample_fn)

    if args.print_tree:
        sbt.print_sbt(root)

    print('\n*** Saving to disk')
    fn = sbt.save_sbt(root, args.save_prefix)
    print('--- Save to', fn)

if __name__ == '__main__':
    main()
