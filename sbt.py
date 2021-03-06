#!/usr/bin/env python

"""
A trial implementation of sequence bloom trees, Solomon & Kingsford, 2015.

This is a simple in-memory version where all of the graphs are in
memory at once; to move it onto disk, the graphs would need to be
dynamically loaded for each query.

To try it out, do::

    factory = GraphFactory(ksize, tablesizes)
    root = Node(factory)

    graph1 = factory.create_nodegraph()
    # ... add stuff to graph1 ...
    leaf1 = Leaf("a", graph1)
    root.add_node(leaf1)

For example, ::

    # filenames: list of fa/fq files
    # ksize: k-mer size
    # tablesizes: Bloom filter table sizes

    factory = GraphFactory(ksize, tablesizes)
    root = Node(factory)

    for filename in filenames:
        graph = factory.create_nodegraph()
        graph.consume_fasta(filename)
        leaf = Leaf(filename, graph)
        root.add_node(leaf)

then define a search function, ::

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.graph.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * len(seq)):
            return 1
        return 0
"""

from __future__ import print_function, unicode_literals

from collections import namedtuple
import hashlib
import json
import math
import os
import random
import shutil
from tempfile import NamedTemporaryFile

import khmer
from khmer import khmer_args
from random import randint
from numpy import array


NodePos = namedtuple("NodePos", ["pos", "node"])


class GraphFactory(object):
    "Build new nodegraphs (Bloom filters) of a specific (fixed) size."

    def __init__(self, ksize, starting_size, n_tables):
        self.ksize = ksize
        self.starting_size = starting_size
        self.n_tables = n_tables

    def create_nodegraph(self):
        return khmer.Nodegraph(self.ksize, self.starting_size, self.n_tables)


class SBT(object):

    def __init__(self, factory):
        self.factory = factory
        self.nodes = [None]

    def add_node(self, node):
        try:
            pos = self.nodes.index(None)
        except ValueError:
            # There aren't any empty positions left.
            # Extend array
            current_size = len(self.nodes)
            self.nodes += [None] * (current_size + 1)
            pos = current_size

        if pos == 0:  # empty tree
            self.nodes[0] = node
            return

        p = self.parent(pos)
        c1, c2 = self.children(p.pos)
        if isinstance(p.node, Leaf):
            # Create a new internal node
            # node and parent are children of new internal node
            n = Node(self.factory, name="internal." + str(p.pos))
            self.nodes[p.pos] = n

            self.nodes[c1.pos] = p.node
            self.nodes[c2.pos] = node

            n.graph.update(p.node.graph)
            n.graph.update(node.graph)

            # update all parents!
            p = self.parent(p.pos)
            while p:
                p.node.graph.update(node.graph)
                p = self.parent(p.pos)

            return

        # Doing this this way always return that parent node is a Leaf,
        # so this piece of code is not necessary
#        if isinstance(p.node, Node):
#            # Parent node can accommodate a new child.
#            p.node.graph.update(node.graph)
#            if self.nodes[c1.pos] is None:
#                self.nodes[c1.pos] = node
#            elif self.nodes[c2.pos] is None:
#                self.nodes[c2.pos] = node

    def find(self, search_fn, *args):
        matches = []
        visited, queue = set(), [0]
        while queue:
            node_p = queue.pop(0)
            node_g = self.nodes[node_p]
            if node_p not in visited and node_g is not None:
                visited.add(node_p)
                if search_fn(node_g, *args):
                    if isinstance(node_g, Leaf):
                        matches.append(node_g)
                    elif isinstance(node_g, Node):
                        queue.extend(c.pos for c in self.children(node_p))
        return matches

    def parent(self, pos):
        if pos == 0:
            return None
        p = int(math.floor((pos - 1) / 2))
        return NodePos(p, self.nodes[p])

    def children(self, pos):
        c1 = 2 * pos + 1
        c2 = 2 * (pos + 1)
        return (NodePos(c1, self.nodes[c1]),
                NodePos(c2, self.nodes[c2]))

    def save(self, tag):
        dirname = '.sbt.' + tag

        if not os.path.exists(dirname):
            os.makedirs(dirname)

        structure = []
        for node in self.nodes:
            if node is None:
                structure.append(None)
                continue

            data = {
                'filename': os.path.join('.sbt.' + tag,
                                         '.'.join([tag, node.name, 'sbt'])),
                'name': node.name
            }
            if isinstance(node, Leaf):
                data['metadata'] = node.metadata

            node.graph.save(data['filename'])
            structure.append(data)

        fn = tag + '.sbt.json'
        with open(fn, 'w') as fp:
            json.dump(structure, fp)

        return fn

    @staticmethod
    def load(sbt_fn):
        with open(sbt_fn) as fp:
            nodes = json.load(fp)

        if nodes[0] is None:
            # TODO error!
            raise ValueError("Empty tree!")

        sbt_nodes = []

        ksize, tablesize, ntables, _, _, _ = khmer.extract_nodegraph_info(nodes[0]['filename'])
        factory = GraphFactory(ksize, tablesize, ntables)

        for node in nodes:
            if node is None:
                sbt_nodes.append(None)
                continue

            graph = khmer.load_nodegraph(node['filename'])

            if 'metadata' in node:
                # only Leaf nodes have metadata
                l = Leaf(node['metadata'], graph)
                sbt_nodes.append(l)
            else:
                n = Node(factory, name=node['name'])
                n.graph = graph
                sbt_nodes.append(n)

        tree = SBT(factory)
        tree.nodes = sbt_nodes

        return tree

    def print_dot(self):
        print("""
        digraph G {
        nodesep=0.3;
        ranksep=0.2;
        margin=0.1;
        node [shape=circle];
        edge [arrowsize=0.8];
        """)

        for i, node in enumerate(self.nodes):
            if node is None:
                continue

            p = self.parent(i)
            if p is not None:
                if isinstance(node, Leaf):
                    print('"', p.pos, '"', '->', '"', node.name, '";')
                else:
                    print('"', p.pos, '"', '->', '"', i, '";')
        print("}")

    def print(self):
        visited, stack = set(), [0]
        while stack:
            node_p = stack.pop()
            node_g = self.nodes[node_p]
            if node_p not in visited and node_g is not None:
                visited.add(node_p)
                depth = int(math.floor(math.log(node_p + 1, 2)))
                print(" " * 4 * depth, node_g)
                if isinstance(node_g, Node):
                    stack.extend(c.pos for c in self.children(node_p)
                                       if c.pos not in visited)

class Node(object):
    "Internal node of SBT; has 0, 1, or 2 children."

    def __init__(self, factory, name=None):
        self.factory = factory
        self.graph = factory.create_nodegraph()

        self.name = name

    def __str__(self):
        return '*Node:{name} [occupied: {nb}, fpr: {fpr:.2}]'.format(
                name=self.name, nb=self.graph.n_occupied(),
                fpr=khmer.calc_expected_collisions(self.graph, True, 1.1))


class Leaf(object):
    def __init__(self, metadata, nodegraph, name=None):
        self.metadata = metadata
        if name is None:
            name = metadata
        self.name = name
        self.graph = nodegraph

    def __str__(self):
        return '**Leaf:{name} [occupied: {nb}, fpr: {fpr:.2}] -> {metadata}'.format(
                name=self.name, metadata=self.metadata,
                nb=self.graph.n_occupied(),
                fpr=khmer.calc_expected_collisions(self.graph, True, 1.1))


def filter_distance( filter_a, filter_b, n=1000 ) :
    """
    Compute a heuristic distance per bit between two Bloom
    filters.

    filter_a : First filter
    filter_b : Second filter
    n        : Number of positions to compare (in groups of 8)
    """
    A = filter_a.graph.get_raw_tables()
    B = filter_b.graph.get_raw_tables()
    distance = 0
    for q,p in zip( A, B ) :
        a = array( q, copy=False )
        b = array( p, copy=False )
        for i in map( lambda x : randint( 0, len(a) ), range(n) ) :
            distance += sum( map( int, [ not bool((a[i]>>j)&1)
                                           ^ bool((b[i]>>j)&1)
                                         for j in range(8) ] ) )
    return distance / ( 8.0 * len(A) * n )


def test_simple():
    factory = GraphFactory(5, 100, 3)
    root = SBT(factory)

    leaf1 = Leaf("a", factory.create_nodegraph())
    leaf1.graph.count('AAAAA')
    leaf1.graph.count('AAAAT')
    leaf1.graph.count('AAAAC')

    leaf2 = Leaf("b", factory.create_nodegraph())
    leaf2.graph.count('AAAAA')
    leaf2.graph.count('AAAAT')
    leaf2.graph.count('AAAAG')

    leaf3 = Leaf("c", factory.create_nodegraph())
    leaf3.graph.count('AAAAA')
    leaf3.graph.count('AAAAT')
    leaf3.graph.count('CAAAA')

    leaf4 = Leaf("d", factory.create_nodegraph())
    leaf4.graph.count('AAAAA')
    leaf4.graph.count('CAAAA')
    leaf4.graph.count('GAAAA')

    leaf5 = Leaf("e", factory.create_nodegraph())
    leaf5.graph.count('AAAAA')
    leaf5.graph.count('AAAAT')
    leaf5.graph.count('GAAAA')

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def search_kmer(obj, seq):
        return obj.graph.get(seq)

    leaves = [leaf1, leaf2, leaf3, leaf4, leaf5 ]
    kmers = [ "AAAAA", "AAAAT", "AAAAG", "CAAAA", "GAAAA" ]

    def search_kmer_in_list(kmer):
        x = []
        for l in leaves:
            if l.graph.get(kmer):
                x.append(l)

        return set(x)

    for kmer in kmers:
        assert set(root.find(search_kmer, kmer)) == search_kmer_in_list(kmer)

    print('-----')
    print([ x.metadata for x in root.find(search_kmer, "AAAAA") ])
    print([ x.metadata for x in root.find(search_kmer, "AAAAT") ])
    print([ x.metadata for x in root.find(search_kmer, "AAAAG") ])
    print([ x.metadata for x in root.find(search_kmer, "CAAAA") ])
    print([ x.metadata for x in root.find(search_kmer, "GAAAA") ])

def test_longer_search():
    ksize = 5
    factory = GraphFactory(ksize, 100, 3)
    root = SBT(factory)

    leaf1 = Leaf("a", factory.create_nodegraph())
    leaf1.graph.count('AAAAA')
    leaf1.graph.count('AAAAT')
    leaf1.graph.count('AAAAC')

    leaf2 = Leaf("b", factory.create_nodegraph())
    leaf2.graph.count('AAAAA')
    leaf2.graph.count('AAAAT')
    leaf2.graph.count('AAAAG')

    leaf3 = Leaf("c", factory.create_nodegraph())
    leaf3.graph.count('AAAAA')
    leaf3.graph.count('AAAAT')
    leaf3.graph.count('CAAAA')

    leaf4 = Leaf("d", factory.create_nodegraph())
    leaf4.graph.count('AAAAA')
    leaf4.graph.count('CAAAA')
    leaf4.graph.count('GAAAA')

    leaf5 = Leaf("e", factory.create_nodegraph())
    leaf5.graph.count('AAAAA')
    leaf5.graph.count('AAAAT')
    leaf5.graph.count('GAAAA')

    root.add_node(leaf1)
    root.add_node(leaf2)
    root.add_node(leaf3)
    root.add_node(leaf4)
    root.add_node(leaf5)

    def kmers(k, seq):
        for start in range(len(seq) - k + 1):
            yield seq[start:start + k]

    def search_transcript(node, seq, threshold):
        presence = [ node.graph.get(kmer) for kmer in kmers(ksize, seq) ]
        if sum(presence) >= int(threshold * (len(seq) - ksize + 1)):
            return 1
        return 0

    try1 = [ x.metadata for x in root.find(search_transcript, "AAAAT", 1.0) ]
    assert set(try1) == set([ 'a', 'b', 'c', 'e' ]), try1 # no 'd'

    try2 = [ x.metadata for x in root.find(search_transcript, "GAAAAAT", 0.6) ]
    assert set(try2) == set([ 'a', 'b', 'c', 'd', 'e' ])

    try3 = [ x.metadata for x in root.find(search_transcript, "GAAAA", 1.0) ]
    assert set(try3) == set([ 'd', 'e' ]), try3
