from utils.termsize import get_terminal_size_linux


def test_returns_width_height_pair():
    size = get_terminal_size_linux()
    assert isinstance(size, tuple)
    assert len(size) == 2
