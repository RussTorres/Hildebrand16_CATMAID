#!/usr/bin/env python
'''
Skeleton algorithms for morphology analysis on
    networkx representations of CATMAID skeletons
'''
import logging

import networkx


def name(sk):
    name = sk['neuron']['neuronname']
    if ' ' in name:
        tokens = name.split()
        if len(tokens) == 2 and tokens[0] == 'neuron':
            return tokens[1]
    return name


def nodes(sk):
    vertices = {}
    for v in sk['vertices']:
        if sk['vertices'][v]['type'] == 'skeleton':
            vertices[v] = sk['vertices'][v]
    return vertices


def soma(sk):
    soma = None
    for v in sk['vertices']:
        for l in sk['vertices'][v]['labels']:
            if l == 'soma':
                if soma is not None:
                    logging.critical(
                        "Found 2 somas [%s, %s] in neuron %s",
                        soma, v, name(sk))

                    raise ValueError(
                        "Found 2 somas [{}, {}] in neuron {}".format(
                            soma, v, name(sk)))
                # soma = sk['vertices'][v]
                soma = v
                break
    return soma


def dedges(sk):
    # don't include edges to missing vertices
    dedges = {}
    conns = sk['connectivity']
    verts = sk['vertices']
    for cid in conns:
        if cid not in verts:
            continue
        for pid in conns[cid]:
            if pid not in verts:
                continue
            if conns[cid][pid]['type'] != 'neurite':
                continue
            dedges[cid] = dedges.get(cid, []) + [pid, ]
    return dedges


def redges(neuron):
    redges = {}
    dedges = neuron.dedges
    for cid in dedges:
        for pid in dedges[cid]:
            redges[pid] = redges.get(pid, []) + [cid, ]
    return redges


def edges(neuron):
    edges = {}
    for cid in neuron.dedges:
        for pid in neuron.dedges[cid]:
            edges[cid] = edges.get(cid, []) + [pid, ]
            edges[pid] = edges.get(pid, []) + [cid, ]
    return edges


def tags(sk):
    all_tags = {}
    for v in sk['vertices']:
        for l in sk['vertices'][v]['labels']:
            all_tags[l] = all_tags.get(l, []) + [v, ]
    return all_tags


def root(neuron):
    sg = networkx.topological_sort(neuron.dgraph)
    if not len(sg):
        if len(neuron.nodes) == 1:
            return neuron.nodes.keys()[0]
        else:
            logging.warning(
                "Attempt to get root for neuron[%s] with %s nodes",
                neuron.name, len(neuron.nodes))
            return None
    return sg[0]


def leaves(neuron):
    return [n for n in neuron.dgraph if neuron.dgraph.out_degree(n) == 0]


def get_id(skeleton):
    return skeleton.get('id', None)


def annotation(skeleton):
    return skeleton['neuron'].get('annotations', [])


def bifurcations(neuron):
    return [n for n in neuron.dgraph if neuron.dgraph.out_degree(n) > 1]
