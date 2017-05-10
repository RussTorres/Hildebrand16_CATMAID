#!/usr/bin/env python
"""
Neuron class
lazy attributes, so algorithms can get attributes without requiring
lots of pre-computation on startup
Basic attributes:
    - skeleton [loaded from json on init]
        - neuron name
        - vertices
    - meta attributes [lazy]
        - soma
        - axons
        - tags
        - root
"""
import logging
import json

from . import algorithms


try:
    xrange
except NameError as E:
    xrange = range

try:
    unicode
except NameError as E:
    unicode = str


def load_skeleton(skeleton):
    '''
    load_skeleton
    loads skeleton from file. Parameter skeleton must be a file name
    returns dictionary skeleton object.
    '''
    if isinstance(skeleton, (str, unicode)):
        with open(skeleton, 'r') as skf:
            sk = json.load(skf)
        skid = skeleton.strip('skel.json')
        if skid.isdigit():
            sk['id'] = int(skid)
        return sk
    elif isinstance(skeleton, dict):
        return skeleton
    elif isinstance(skeleton, list):
        skeleton = algorithms.skeleton_json_new_to_old.convert_new_to_old(
            skeleton)
        logging.info("Converting skeleton to dictionary.")
        return skeleton
    raise Exception("Failed to load skeleton: {}".format(skeleton))


def lazy(f):
    """Make a lazy method
    The resulting method will evaluate once, the first time called.
    """
    attr = '_{}'.format(f.__name__)

    def wrapped(self):
        if not hasattr(self, attr):
            setattr(self, attr, f(self))
        return getattr(self, attr)

    wrapped.__name__ = f.__name__
    wrapped.__doc__ = f.__doc__
    return wrapped


def lazyproperty(f):
    """Shortcut to make a lazy property"""
    return property(lazy(f))


class Neuron(object):
    """
    Neuron Object holds a skeleton and preforms actions on the skeleton.
    Parameters
    ----------
    skeleton: either a skeleton(dictonary) or a filename with a skeleton in it.
    """
    def __init__(self, skeleton):
        self.skeleton = load_skeleton(skeleton)

    def __repr__(self):
        return "<%s.%s at %s: neuron_id: %s, skeleton_id: %s>" % (
            self.__module__, self.__class__.__name__, hex(id(self)),
            self.name, self.skeleton_id)

    @lazyproperty
    def nodes(self):
        return algorithms.skeleton.nodes(self.skeleton)

    @lazyproperty
    def edges(self):
        """ [child][parent] """
        return algorithms.skeleton.edges(self)

    @lazyproperty
    def dedges(self):
        """ [child][parent] """
        return algorithms.skeleton.dedges(self.skeleton)

    @lazyproperty
    def redges(self):
        """ [parent][child] """
        return algorithms.skeleton.redges(self)

    @lazyproperty
    def skeleton_id(self):
        return algorithms.skeleton.get_id(self.skeleton)

    @lazyproperty
    def name(self):
        return algorithms.skeleton.name(self.skeleton)

    @lazyproperty
    def dgraph(self):
        return algorithms.graph.dgraph(self)

    @lazyproperty
    def graph(self):
        return algorithms.graph.graph(self)

    @lazyproperty
    def soma(self):
        return algorithms.skeleton.soma(self.skeleton)

    @lazyproperty
    def tags(self):
        return algorithms.skeleton.tags(self.skeleton)

    @lazyproperty
    def root(self):
        return algorithms.skeleton.root(self)

    @lazyproperty
    def leaves(self):
        return algorithms.skeleton.leaves(self)

    # TODO annotation support via json api dropped in recent version
    # @lazyproperty
    # def annotations(self):
    #     return algorithms.skeleton.annotation(self.skeleton)
    @lazyproperty
    def bifurcations(self):
        return algorithms.skeleton.bifurcations(self)
