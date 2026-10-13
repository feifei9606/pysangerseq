"""变异描述解析 + INDEL 断点位置校正。

输入格式（与原脚本一致，两种分隔符都支持）：
    12-14796923-G-A        12:14796923:G:A
    12-14796923            （只给位点，ref/alt 为 None）
    6-26091000-CTTGGT-CTT  （INDEL，ref/alt 长度不等）

`adjust_position` 负责把"变异位点"挪到实际应该标注的位置：在重复序列
（poly-A、短串联重复等）里同一段插入/缺失有多个等价写法，原脚本用
"把重复单元从参考序列里剥掉、数剥掉的长度" 的办法定位真正的断点。
逻辑照搬原 R 实现，未做语义改动。

> ⚠️ 唯一一处修正：R 在 MNV 分支里有 `refsequence != "ref" && altsequence != "alt"`
> 的哨兵守卫；第一版移植漏了它，导致"只给位点 + 反链读序"时 `altposs` 被平移
> `len("alt") == 3`，图上表现为变异位点偏离正中、且红色箭头消失。
> 现由 `adjust_position` 开头的守卫显式处理。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: Optional[str]
    alt: Optional[str]

    @property
    def is_indel(self) -> bool:
        if self.ref is None or self.alt is None:
            return False
        return len(self.ref) != len(self.alt)

    def describe(self) -> str:
        if self.ref is None or self.alt is None:
            return "%s:%d" % (self.chrom, self.pos)
        return "%s:%d-%s>%s" % (self.chrom, self.pos, self.ref, self.alt)


def parse_variant(spec: str) -> Variant:
    """解析 `chr-pos-ref-alt` / `chr:pos:ref:alt`。"""
    spec = spec.strip()
    if ":" in spec:
        parts = spec.split(":")
    elif "-" in spec:
        parts = spec.split("-")
    else:
        raise ValueError("变异描述格式不对: %r" % spec)
    if len(parts) < 2:
        raise ValueError("变异描述至少需要 chr 和 pos: %r" % spec)
    chrom = parts[0]
    try:
        pos = int(parts[1])
    except ValueError:
        raise ValueError("变异位置不是整数: %r" % parts[1])

    ref = parts[2] if len(parts) > 2 and parts[2] != "" else None
    alt = parts[3] if len(parts) > 3 and parts[3] != "" else None
    if ref is None or alt is None:
        ref = alt = None
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt)


def _strip_repeat_unit(seq: str) -> str:
    """求 seq 的最小重复单元。

    对应原 R 的 for(i in 1:ceiling(inslen/2)) + while 剥除循环：
    能把 seq 完整剥空的最小前缀即为重复单元，否则返回 seq 本身。
    """
    if not seq:
        return seq
    for i in range(1, (len(seq) + 1) // 2 + 1):
        cand = seq[:i]
        rest = seq
        while rest.startswith(cand):
            rest = rest[len(cand) :]
        if len(rest) == 0:
            return cand
    return seq


def _posins(anchor: str, repeat_unit: str) -> int:
    """anchor 中剥掉重复单元后，返回"被剥掉的总长度"。

    对应原 R 中 nchar(insrefseqerf) - nchar(insrefseq1)。
    """
    s = anchor
    while s.startswith(repeat_unit):
        s = s[len(repeat_unit) :]
    for i in range(1, len(repeat_unit) + 1):
        if s[i - 1 : i] != repeat_unit[i - 1 : i]:
            if i != 1:
                s = s[i - 1 :]
            break
    return len(anchor) - len(s)


def adjust_position(
    chrom: str,
    altposs: int,
    refsequence: Optional[str],
    altsequence: Optional[str],
    strand: str,
    fasta,
) -> int:
    """把变异位点校正到重复序列中的正确断点位置。

    strand 为第一次 swalign 得到的链方向（'+'/'-'），即原脚本的 posref2[7]。

    refsequence / altsequence 为 None 表示"只给了位点、没有 ref/alt"。
    原脚本用字面量 `"ref"` / `"alt"` 当哨兵，并在 MNV 分支里**显式排除**它们
    （`nchar(refsequence) > 1 && refsequence != "ref" && ...`）—— 也就是
    只给位点时不做任何位置校正。两个哨兵字面量这里也一并认（便于按 R 风格调用）。
    """
    if not refsequence or not altsequence or refsequence == "ref" or altsequence == "alt":
        return altposs

    minus = strand == "-"
    ref_len = len(refsequence)

    if ref_len < len(altsequence):
        # ------------------------------ 插入
        if altsequence.startswith(refsequence):
            inserted = altsequence[ref_len:]
            unit = _strip_repeat_unit(inserted)
            start = altposs + ref_len
            anchor = fasta.fetch(chrom, start, start + 599)  # faidx 取 601bp，只拼前 600
            posins = _posins(anchor, unit)
            if not minus:
                altposs = altposs + posins + ref_len
            # minus 链位置不动
        else:
            if minus:
                altposs = altposs + ref_len - 1

    elif ref_len > len(altsequence):
        # ------------------------------ 缺失
        alt_len = len(altsequence)
        if refsequence.startswith(altsequence):
            deleted = refsequence[alt_len:]
            unit = _strip_repeat_unit(deleted)
            start = altposs + ref_len
            anchor = fasta.fetch(chrom, start, start + 599)
            posins = _posins(anchor, unit)
            if not minus:
                altposs = altposs + posins + alt_len
            else:
                altposs = altposs + len(deleted)
        else:
            if minus:
                altposs = altposs + ref_len - 1

    elif ref_len > 1 and len(altsequence) > 1:
        # ------------------------------ 等长多碱基（MNV）
        # 注意：只有"真·多碱基"才进来。1bp SNV（ref/alt 各 1 字符）不走这里；
        # "只给位点"的哨兵已在函数开头 return，也不会走到这里。
        if minus:
            altposs = altposs + len(altsequence)

    return altposs
