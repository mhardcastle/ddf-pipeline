from utils import parset


def test_option_list_entries_are_well_formed():
    for entry in parset.option_list:
        assert len(entry) >= 4
        section, name, otype, _default = entry[:4]
        assert isinstance(section, str) and section
        assert isinstance(name, str) and name
        assert callable(otype)


def test_core_options_present():
    names = {entry[1] for entry in parset.option_list}
    for expected in ("mslist", "ndir", "imsize"):
        assert expected in names
