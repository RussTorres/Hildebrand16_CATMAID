#!/usr/bin/env python
import networkx


def dgraph(neuron):
    """
    Creates a networkx Graph out of the skeleton edges.
    Note: not a network graph of all skeletons. Is a graph representation
    of a single neuron.
    """
    dgraph = networkx.DiGraph()
    dedges = neuron.dedges
    for cid in dedges:
        for pid in dedges[cid]:
            dgraph.add_edge(pid, cid)
    return dgraph


def graph(neuron):
    """
    Creates a networkx Graph out of a skeleton edges. Instead of
    Directed it is Undirected.
    """
    return neuron.dgraph.to_undirected()
