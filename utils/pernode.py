"""Pure per-node option distribution used by the MPI code paths.

Kept free of DDFacet / casacore / mpi4py imports so it can be unit-tested on a
laptop and in CI. ``scripts/pipeline.py`` calls it from ``full_clearcache_mpi``
and ``redo_dppp_di_mpi``.
"""
import copy


def per_node_config(base_options, dico_mslist, dico_fullmslist):
    """Build one option dict per node that has a full mslist.

    A node can appear in ``dico_fullmslist`` without a per-node entry in
    ``dico_mslist`` (more MPI nodes than measurement sets). Such a node falls
    back to the global ``base_options["mslist"]`` instead of raising KeyError
    (regression: commit de84e2b).

    Returns a list of ``(node, options)`` pairs, where ``options`` is a deep
    copy of ``base_options`` with ``mslist`` and ``full_mslist`` set.
    """
    configs = []
    for node in dico_fullmslist:
        options = copy.deepcopy(base_options)
        options["mslist"] = dico_mslist.get(node, base_options["mslist"])
        options["full_mslist"] = dico_fullmslist[node]
        configs.append((node, options))
    return configs
