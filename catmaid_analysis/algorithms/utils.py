#!/usr/bin/env python
'''
'''
import logging


def get_legacy_skeleton(skeleton_json):
    """
    -- old --
    vert types are always 'skeleton' or 'connector'
    conn types are always 'presynaptic_to', 'postsynaptic_to', or 'neurite'
    for skeleton verts
        labels, radius, type, x, y, z
    for connector verts
        labels, reviewer_id, type, x, y, z
    for all conns: no sub-keys
    also ['neuron']['neuronname']
    """
    if len(skeleton_json) == 5:
        name, nodes, tags, connectors, _ = skeleton_json
        skid = None
        neuron_id = None
        annotations = None
    else:
        name, nodes, tags, connectors, _, skid, neuron_id, annotations = skeleton_json
    verts = {}
    conns = {}
    # add verticies and neurites
    for nid, pid, uid, x, y, z, r, conf, ct, et in nodes:
        if pid is not None:
            conns[str(nid)] = {str(pid): {'type': 'neurite'}}
        verts[str(nid)] = {
            'x': x, 'y': y, 'z': z, 'radius': r,
            'type': 'skeleton', 'labels': [],
            'confidence': conf,
        }
    # add connectors
    for tid, cid, post, x, y, z, ct in connectors:
        # 'connector' verts
        scid = str(cid)
        # neuron to conn relation: e.g. neuron is presynaptic_to the conn
        stype = 'postsynaptic_to' if post else 'presynaptic_to'
        if scid not in verts:
            # NOTE review information is not being copied over
            verts[str(cid)] = {
                'x': x, 'y': y, 'z': z, 'type': 'connector', 'labels': []
            }
        # 'conns'
        stid = str(tid)
        if stid in conns:
            conns[stid][str(cid)] = {'type': stype}
        else:
            conns[stid] = {str(cid): {'type': stype}}
    # copy over tags
    for t in tags:
        for nid in tags[t]:
            verts[str(nid)]['labels'].append(t)
    return {
        'neuron': {
            'neuronname': name,
            'id': neuron_id,
            'annotations': annotations},
        'vertices': verts,
        'connectivity': conns,
        'id': skid}
