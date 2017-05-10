#!/usr/bin/env python
'''
Simple reading of skeleton jsons from CATMAID.

'''

import json
import logging
import requests
from . import neuron
from .algorithms import utils

CATMAID_SERVER = 'http://hildebrand16.neurodata.io/catmaid/'
PROJECT_ID = 6
SKELETON_SOURCE = '{server}{project}/'.format(
    server=CATMAID_SERVER, project=PROJECT_ID)
DEFAULT_TIMEOUT = None


class ServerSkeletonSource:
    def __init__(self, skel_source, cache=True):
        self._cache = ({} if cache else None)
        self._skel_source = skel_source

    def get_skeleton_ids(self):
        self.skeleton_ids = requests.get(skel_source).json()

    def clear_cache(self):
        self._cache = ({} if self._cache is not None else None)

    def _fetchskeleton(self, skeleton_id):
        '''
        retrieves skeleton from current API and
            returns a more human-readable legacy format
        '''
        skeletonurl = '{s}/skeleton/{sk}/json'.format(
            s=self._skel_source, sk=skeleton_id)
        newskel = requests.get(skeletonurl).json()
        return utils.get_legacy_skeleton(newskel)

    def get_skeleton(self, skeleton_id):
        if self._cache is None:
            skeleton = self._fetchskeleton(skeleton_id)
        elif skeleton_id in self._cache:
            logging.debug("fetching skeleton {}".format(skeleton_id))
            return self._cache[skeleton_id].skeleton
        else:
            skeleton = self._fetchskeleton(skeleton_id)
        if self._cache is not None:
            logging.debug("caching skeleton {}".format(skeleton_id))
            self._cache[skeleton_id] = neuron.Neuron(skeleton)
        return skeleton

    def get_neuron(self, skeleton_id):
        n = neuron.Neuron(self.get_skeleton(skeleton_id))
        n.skeleton['id'] = (skeleton_id if n.skeleton_id is None
                            else n.skeleton_id)
        return n

    def all_neurons_iter(self):
        return (self.get_neuron(sk) for sk in self.get_skeleton_ids())

    def all_skeletons_iter(self):
        return (self.get_skeleton(sk) for sk in self.get_skeleton_ids())
