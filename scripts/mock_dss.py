from random import randint


class Maker:
    CHRS = ["chr1", "chr2", "chr3", "puc_chr", "lambda"]

    def __init__(self):
        self.chri = 0
        self.pos = 0

    def make_row(self):
        chrom = self.CHRS[self.chri]
        pos = self.pos
        n = randint(1, 10)
        x = randint(1, n)

        if self.chri < len(self.CHRS) - 1 and randint(0, 10) > -1:
            self.chri += 1
        self.pos += randint(1, 10)

        return f"{chrom}\t{pos}\t{n}\t{x}\n"

    def make(self, length):
        o = "chr\tpos\tn\tx\n"
        for i in range(length):
            o += self.make_row()
        return o

    @staticmethod
    def split_row(row):
        chrom, pos, n, x = row.strip().split("\t")
        pos = int(pos)
        n = int(n)
        x = int(x)
        x1 = randint(0, x)
        x2 = x - x1
        n1 = randint(x1, n - x2)
        n2 = n - n1
        assert x1 <= n1
        assert x2 <= n2
        assert x1 + x2 == x
        assert n1 + n2 == n
        row1 = ""
        row2 = ""
        if n1 > 0:
            row1 = f"{chrom}\t{pos}\t{n1}\t{x1}\n"
        if n2 > 0:
            row2 = f"{chrom}\t{pos+1}\t{n2}\t{x2}\n"
        return row1, row2

    @staticmethod
    def split(s):
        xs = s.splitlines()
        t = xs[0] + "\n"
        b = xs[0] + "\n"
        for line in xs[1:]:
            left, right = Maker.split_row(line)
            t += left
            b += right
        return t, b


def make_triple(length):
    while True:
        o = Maker().make(length)
        t, b = Maker.split(o)
        if t.splitlines()[1].split("\t")[0] != b.splitlines()[1].split("\t")[0]:
            continue
        return t, b, o


def diagnostic(t, b, o, res):
    return f"{t}-----\n{b}-----\n{o}-----\n{res}"
