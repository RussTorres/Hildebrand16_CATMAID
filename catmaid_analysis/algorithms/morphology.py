#!/usr/bin/env python

import networkx

try:
    import numpy
    has_numpy = True
except ImportError:
    has_numpy = False


try:
    xrange
except NameError as E:
    xrange = range


def node_array(neuron, node_list=None, include_nid=False):
    '''
    node_array(neuron, node_list=None)
    Creates an array of cartesian coordinates (x, y and z) from an input Neuron
    Input:
     - Neuron    (from Neuron class)
     - node_list (list of nodes with xyz coordinates)
    Output:
     - Array of all nodes (with xyz coordinates) that are associated
       with the input Neuron
    '''
    if node_list is None:
        node_list = neuron.nodes.keys()
    nodes = []
    for nid in node_list:
        node = neuron.nodes[nid]
        if include_nid:
            nodes.append((nid, node['x'], node['y'], node['z']))
        else:
            nodes.append((node['x'], node['y'], node['z']))
    if has_numpy:
        return numpy.array(nodes)
    return nodes


def node_position(node):
    '''returns position of node or connector'''
    position = (node['x'], node['y'], node['z'])
    if has_numpy:
        return numpy.array(position)
    return position


def distance(neuron, v0, v1):
    """
    distance(neuron, v0, v1)
    Returns the distance between v0 and v1
    Input:
     - Neuron (from Neuron class)
     - v0 (string representing node ID from Neuron)
     - v1 (string representing node ID from Neuron)
    Output:
     - distance between v0 and v1 in 3D space
    """
    return sum([
        (neuron.skeleton['vertices'][v0][k] -
         neuron.skeleton['vertices'][v1][k]) ** 2.
        for k in ('x', 'y', 'z')]) ** 0.5


def find_path(neuron, v0, v1):
    """
    find_path(neuron, v0, v1)
    Finds the shortest path between v0 and v1 using networkx functions
    Input:
     - Neuron (from Neuron class)
     - v0 (string representing node ID from Neuron)
     - v1 (string representing node ID from Neuron)
    Output:
     - An array of node IDs, or a Non-directed networkx graph representing the
       shortest path between v0 and v1
    """
    return networkx.shortest_path(neuron.graph, v0, v1)


def path_length(neuron, v0, v1):
    """
    path_length(neuron, v0, v1)
    Finds the length of the shortest path between v0 and v1 using
    distance function
    Input:
     - Neuron (from Neuron class)
     - v0 (string representing node ID from Neuron)
     - v1 (string representing node ID from Neuron)
    Output:
     - sum of distances between each node in the shortest path between v0 and
       v1 in 3D space
    """
    path = find_path(neuron, v0, v1)
    return sum([
        distance(neuron, path[i], path[i+1]) for i in xrange(len(path) - 1)])


def branch_order(neuron, v, base=None):
    """if base is None, default to soma"""
    if base is None:
        base = neuron.soma
        if neuron.soma is None:
            return None
    path = networkx.shortest_path(neuron.graph, base, v)
    order = 0
    for v in path:
        if len(neuron.edges[v]) > 2:
            order += 1
    return order


def unique_neurites(neu, base=None):
    '''
    This function generates lists of unique neurites based off branching
        structure of neuron object
    '''
    if base is None:
        base = neu.root
    neurites = []
    for leaf in neu.leaves:
        path = networkx.shortest_path(neu.graph, base, leaf)
        branchinpath = set(path).intersection(
            set(neu.bifurcations).union(set([base])))
        bpargs = [path.index(bp) for bp in branchinpath]
        for bparg in sorted(bpargs, reverse=True):
            neurites.append(path[bparg:])
            path = path[:(bparg + 1)]
    return list(set([tuple(neurite) for neurite in neurites]))


def total_pathlength(neuron):
    return sum([sum([distance(neuron, path[i], path[i+1])
                     for i in xrange(len(path) - 1)])
                for path in unique_neurites(neuron)])
