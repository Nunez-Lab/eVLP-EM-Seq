from combine_dss import Runner

import mock_dss as md
from random import seed


def input_mock(s):
    return iter([x + "\n" for x in s.splitlines()])


class OutputMock:
    def __init__(self):
        self.s = ""

    def write(self, x):
        self.s += x

    def get(self):
        return self.s


def test_string(t, b, o):
    om = OutputMock()
    r = Runner(md.Maker.CHRS, input_mock(t), input_mock(b), om)
    r.run()
    res = om.get()
    assert res == o, "\n" + md.diagnostic(t, b, o, res)


if __name__ == "__main__":
    for i in range(100):
        seed(i)
        t, b, o = md.make_triple(10000)
        test_string(t, b, o)
