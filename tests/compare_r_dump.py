#!/usr/bin/env python3
"""把 Python 实现的中间量与原 R 流程的 dump 逐项比对。

用法:
    python3 tests/compare_r_dump.py <dump_dir> <ab1> <chr> <pos> <ref> <alt> <fasta>

dump_dir 由 tests/make_r_dump.r 生成。
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sanger_plot import abif, align, sangerseq as ss  # noqa: E402
from sanger_plot.reference import Fasta  # noqa: E402
from sanger_plot.variant import adjust_position  # noqa: E402
from sanger_plot.pipeline import (  # noqa: E402
    locate_anchor,
    locate_and_trim,
    parse_strand,
    strand_name,
)

FAIL = []

# forward/reverse 的全部可接受写法（大小写不敏感、允许前后空白）+ 期望的归一化结果
_STRAND_ALIAS_EXPECTED = (
    ("forward", "forward"), ("F", "forward"), ("f", "forward"), ("fwd", "forward"),
    ("+", "forward"), ("reverse", "reverse"), ("R", "reverse"), ("r", "reverse"),
    ("rev", "reverse"), ("-", "reverse"), ("REVERSE", "reverse"),
    ("  Reverse  ", "reverse"), (True, "reverse"), (False, "forward"),
)
# 这些必须被 ValueError 拒绝
_STRAND_BAD_VALUES = ("sideways", "", "fw", "rr", "0")


def _strand_alias_table():
    """把别名表跑成字符串结果，顺带确认非法值确实抛 ValueError。"""
    got = ["%r->%s" % (alias, strand_name(parse_strand(alias))) for alias, _ in _STRAND_ALIAS_EXPECTED]
    for bad in _STRAND_BAD_VALUES:
        try:
            parse_strand(bad)
        except ValueError:
            got.append("%r->ValueError" % (bad,))
        else:
            got.append("%r->NO-RAISE" % (bad,))
    return got


def _strand_alias_want():
    want = ["%r->%s" % (alias, name) for alias, name in _STRAND_ALIAS_EXPECTED]
    want += ["%r->ValueError" % (bad,) for bad in _STRAND_BAD_VALUES]
    return want


def check(name, got, want, tol=0.0):
    ok = False
    if isinstance(want, str) and isinstance(got, str):
        ok = got == want
    else:
        try:
            ok = np.allclose(np.asarray(got, dtype=float), np.asarray(want, dtype=float),
                             rtol=tol, atol=tol, equal_nan=True)
        except Exception:
            ok = got == want
    print("  %-28s %s" % (name, "OK" if ok else "FAIL"))
    if not ok:
        FAIL.append(name)
        g = np.asarray(got)
        w = np.asarray(want)
        if g.shape == w.shape and g.size <= 12:
            print("      got :", g)
            print("      want:", w)
        else:
            print("      got shape", g.shape, "want shape", w.shape)
            try:
                if g.shape == w.shape:
                    d = np.nonzero(~np.isclose(g.astype(float), w.astype(float), equal_nan=True))[0]
                    print("      first diffs at", d[:10], "n_diff=", len(d))
                    if len(d):
                        i = d[0]
                        print("      got[i]=", g[i], "want[i]=", w[i])
            except (TypeError, ValueError):
                # 字符串列表等无法转数值的情形：直接把两边打出来
                print("      got  :", list(g)[:20])
                print("      want :", list(w)[:20])


def rd(dump_dir, fname):
    with open(os.path.join(dump_dir, fname), "r") as fh:
        lines = fh.read().splitlines()
    return lines


def rd_table(dump_dir, fname):
    lines = rd(dump_dir, fname)
    header = lines[0].split("\t")
    rows = [l.split("\t") for l in lines[1:] if l]
    return header, rows


def rd_num_table(dump_dir, fname):
    """读成浮点矩阵，R 的 NA 转成 NaN。"""
    _, rows = rd_table(dump_dir, fname)
    return np.array([[np.nan if x == "NA" else float(x) for x in r] for r in rows])


def rd_scalar(dump_dir, fname):
    lines = rd(dump_dir, fname)
    if lines and lines[0] == "v":
        lines = lines[1:]
    lines = [l for l in lines if l != ""]
    return lines[-1] if lines else ""


def main():
    dump_dir, ab1, chrom, pos, ref, alt, fasta_path = sys.argv[1:8]
    pos = int(pos)

    print("== 1. ABIF 解析 ==")
    raw = abif.read_abif(ab1)
    check("ab1.trace", raw.trace, rd_num_table(dump_dir, "ab1_trace.tsv"))
    check("ab1.peak_pos", raw.peak_pos, rd_num_table(dump_dir, "ab1_peakpos.tsv"))

    check("ab1.primary", raw.primary_seq, rd_scalar(dump_dir, "ab1_primary.txt"))
    check("ab1.secondary", raw.secondary_seq, rd_scalar(dump_dir, "ab1_secondary.txt"))

    print("== 2. makeBaseCalls(ratio=0.3) ==")
    mbc = ss.make_base_calls(raw, ratio=0.3)  # read_abif 已直接返回 Sangerseq
    check("mbc.peak_amp", mbc.peak_amp, rd_num_table(dump_dir, "mbc_peakamp.tsv"))
    check("mbc.peak_pos", mbc.peak_pos, rd_num_table(dump_dir, "mbc_peakpos.tsv"))
    check("mbc.primary", mbc.primary_seq, rd_scalar(dump_dir, "mbc_primary.txt"))
    check("mbc.secondary", mbc.secondary_seq, rd_scalar(dump_dir, "mbc_secondary.txt"))

    print("== 3. 参考序列 (原为 makeSCF+read.scf 的产物) ==")
    fa = Fasta(fasta_path)
    region = fa.fetch(chrom, pos - 300, pos + 300)
    ref500 = region[29:529]
    check("ref500", ref500, rd_scalar(dump_dir, "ref500.txt"))

    final = dict(l.split("\t") for l in rd(dump_dir, "final.txt") if "\t" in l)

    print("== 4. setAllelePhase ==")
    phased, table, pa, refstrand = ss.set_allele_phase(mbc, ref500, trim5=29, trim3=30)
    pa_txt = dict(l.split("\t") for l in rd(dump_dir, "pa.txt") if "\t" in l)
    check("seedseq", "".join(table.consensus[29:len(table.consensus) - 30]),
          rd_scalar(dump_dir, "seedseq.txt"))
    if float(pa_txt["score_rc"]) > float(pa_txt["score_fwd"]):
        want_start = int(pa_txt["subRC_start"])
    else:
        want_start = int(pa_txt["sub_start"])
    check("pa.subject_start", pa.subject_start, want_start)
    check("phased.primary", phased.primary_seq, final["primarySeq"])
    check("phased.secondary", phased.secondary_seq, final["secondarySeq"])

    print("== 5. swalign 等价（复用 pipeline.locate_anchor）==")
    p2 = locate_anchor(phased.primary_seq, ref500, 21)
    want2 = final["posref2"].split("|")
    check("posref2.r_start", p2.r_pos, int(want2[4]))
    check("posref2.r_end", p2.r_end, int(want2[5]))
    check("posref2.strand", "-" if p2.rc else "+", want2[6])
    check("posref2.cigar", p2.cigar_str, want2[7])

    print("== 6. 断点位置校正（INDEL 重复序列逻辑）==")
    dlen = 21
    dlen2 = 10
    # ref/alt 为空即"只给位点"；原脚本用字面量 "ref"/"alt" 当哨兵，
    # 这里用 None 表达同一语义（adjust_position 两者都接受）。
    my_altposs = adjust_position(
        chrom, pos, ref or None, alt or None,
        "-" if p2.rc else "+", fa,
    )
    check("altposs_final", my_altposs, int(final["altposs_final"]))
    altposs = int(final["altposs_final"])

    print("== 7. 最终画图参数 ==")
    want1 = final["posref"].split("|")
    p1, trim5len, trim3len = locate_and_trim(
        phased.primary_seq, chrom, altposs, dlen, p2.rc, fa
    )
    check("posref.r_start", p1.r_pos, int(want1[4]))
    check("posref.r_end", p1.r_end, int(want1[5]))
    check("posref.strand", "-" if p1.rc else "+", want1[6])
    check("posref.cigar", p1.cigar_str, want1[7])

    check("trim5len", trim5len, int(final["trim5len"]))
    check("trim3len", trim3len, int(final["trim3len"]))
    check("len", len(phased.secondary_seq), int(final["len"]))

    refseq1 = fa.fetch(chrom, altposs - dlen2, altposs + dlen2)
    if want1[6] == "-":  # 读序比对到参考反向互补（数据驱动，保留）
        refseq1 = align._revcomp(refseq1)
    # 反向测序方向：原 R 脚本由文件名正则推断，Python 侧改为 --strand 入参。
    # 这里用 R dump 的 Rdir_match 反推该传什么，既复现 R 行为又顺带验证 parse_strand。
    want_strand = "reverse" if final["Rdir_match"] == "TRUE" else "forward"
    check("strand(归一化)", strand_name(parse_strand(want_strand)), want_strand)
    if parse_strand(want_strand):
        refseq1 = align._revcomp(refseq1)
    check("refseq1", refseq1, final["refseq1"])

    print("== 8. parse_strand 别名表 ==")
    check("别名归一化", _strand_alias_table(), _strand_alias_want())

    print()
    if FAIL:
        print("FAILED: %d 项 -> %s" % (len(FAIL), ", ".join(FAIL)))
        return 1
    print("ALL OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
