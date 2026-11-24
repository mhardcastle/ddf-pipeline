from utils.options import options


def test_reads_values_from_config(tmp_path):
    cfg = tmp_path / "t.cfg"
    cfg.write_text("[data]\nmslist = foo.txt\n[image]\nimsize = 1234\n")
    option_list = (
        ("data", "mslist", str, None),
        ("image", "imsize", int, 20000),
    )
    result = options([str(cfg)], option_list)
    assert result["mslist"] == "foo.txt"
    assert result["imsize"] == 1234


def test_falls_back_to_default_when_absent(tmp_path):
    cfg = tmp_path / "empty.cfg"
    cfg.write_text("")
    option_list = (("image", "imsize", int, 20000),)
    result = options([str(cfg)], option_list)
    assert result["imsize"] == 20000
