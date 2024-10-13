#!/usr/bin/env python3

import sys


class Runner:
    def __init__(self, chrom_order, tf, bf, of):
        self.chrom_order = chrom_order
        self.tf = tf
        self.bf = bf
        self.of = of

        self.cur_row = None

    def _stage(self, row):
        if self.cur_row is None:
            self.cur_row = row
        elif (
            row["chrom"] == self.cur_row["chrom"]
            and row["pos"] == self.cur_row["pos"]
        ):
            self.cur_row["n"] += row["n"]
            self.cur_row["x"] += row["x"]
        else:
            self._commit()
            self.cur_row = row

    def _commit(self):
        self.of.write(
            f"{self.cur_row['chrom']}\t{self.cur_row['pos']}\t{self.cur_row['n']}\t{self.cur_row['x']}\n"
        )

    def _lte(self, left, right):
        # None is maximal
        if right is None:
            return True
        if left is None:
            return False
        left_idx = self.chrom_order.index(left["chrom"])
        right_idx = self.chrom_order.index(right["chrom"])
        if left_idx < right_idx:
            return True
        if left_idx > right_idx:
            return False
        return left["pos"] <= right["pos"]

    def _parse(self, f, pos_offset):
        try:
            row = next(f)
        except StopIteration:
            return None
        chrom, pos, n, x = row.strip().split("\t")
        # if chrom not in self.chrom_order:
        #     self.chrom_order.append(chrom)
        return {
            "chrom": chrom,
            "pos": int(pos) + pos_offset,
            "n": int(n),
            "x": int(x),
        }

    def run(self):
        header = next(self.tf)
        next(self.bf)
        self.of.write(header)

        top = self._parse(self.tf, 0)
        bot = self._parse(self.bf, -1)

        while True:
            if top is None and bot is None:
                break
            elif self._lte(top, bot):
                self._stage(top)
                top = self._parse(self.tf, 0)
            else:
                self._stage(bot)
                bot = self._parse(self.bf, -1)

        self._commit()


def usage():
    print(
        "usage: {sys.argv[0]} <path to OT dss input> <path to OB dss input> <desired path to combined input>"
    )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        usage()
        sys.exit(0)

    if len(sys.argv) != 4:
        print("error: incorrect number of arguments")
        usage()
        sys.exit(1)

    TOP_PATH = sys.argv[1]
    BOTTOM_PATH = sys.argv[2]
    OUTPUT_PATH = sys.argv[3]

    chrom_order = []
    with open(TOP_PATH, "r") as tf:
        next(tf)
        for line in tf:
            chrom = line.split("\t")[0]
            if chrom not in chrom_order:
                chrom_order.append(chrom)

    with open(TOP_PATH, "r") as tf:
        with open(BOTTOM_PATH, "r") as bf:
            with open(OUTPUT_PATH, "w") as of:
                r = Runner(chrom_order, tf, bf, of)
                r.run()
