"""Regression test for the node/mslist distribution (commit de84e2b).

On a run with more MPI nodes than measurement sets, a node can have an entry in
the full mslist but no per-node mslist. The distribution must fall back to the
global mslist instead of raising KeyError.
"""
import pytest

from utils.pernode import per_node_config, node_mslist_dict


def test_node_missing_from_partial_mslist_falls_back():
    base = {"mslist": "all_mslist.txt", "full_mslist": "all_full.txt", "other": 1}
    dico_mslist = {"nodeA": "local_nodeA_mslist.txt"}
    dico_full = {"nodeA": "local_nodeA_full.txt", "nodeB": "local_nodeB_full.txt"}

    configs = dict(per_node_config(base, dico_mslist, dico_full))

    # Node A keeps its own list.
    assert configs["nodeA"]["mslist"] == "local_nodeA_mslist.txt"
    # Node B has no partial list: it must fall back, not raise.
    assert configs["nodeB"]["mslist"] == "all_mslist.txt"
    assert configs["nodeB"]["full_mslist"] == "local_nodeB_full.txt"


def test_every_fullmslist_node_gets_a_config():
    base = {"mslist": "all_mslist.txt"}
    dico_full = {"n1": "f1", "n2": "f2", "n3": "f3"}
    configs = per_node_config(base, {}, dico_full)
    assert [node for node, _ in configs] == ["n1", "n2", "n3"]


def test_options_are_deep_copied():
    base = {"mslist": "all_mslist.txt", "nested": {"k": 1}}
    configs = per_node_config(base, {}, {"n1": "f1", "n2": "f2"})
    configs[0][1]["nested"]["k"] = 99
    assert base["nested"]["k"] == 1


def test_full_mslist_uses_the_full_dictionary():
    # At 12 nodes the sim has 6 MSs in mslist.txt and 12 in big-mslist.txt, so the
    # two dicts cover different nodes. A full_mslist step must use the full dict, or
    # it silently skips the 6 nodes that only have a full-list entry.
    dico_mslist = {"n%d" % i: "local_n%d_mslist.txt" % i for i in range(6)}
    dico_full = {"n%d" % i: "local_n%d_full_mslist.txt" % i for i in range(12)}
    assert node_mslist_dict(dico_mslist, dico_full, "mslist") is dico_mslist
    assert node_mslist_dict(dico_mslist, dico_full, "full_mslist") is dico_full
    assert len(node_mslist_dict(dico_mslist, dico_full, "full_mslist")) == 12
    with pytest.raises(ValueError):
        node_mslist_dict(dico_mslist, dico_full, "bogus")
