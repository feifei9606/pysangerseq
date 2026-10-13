"""basecalling 与等位基因相位校正。

等价于 sangerseqR 的：
    makeBaseCalls(obj, ratio)    → 由 trace 重新做碱基判读，给出主/次两条读序
    setAllelePhase(obj, refseq, trim5, trim3)
                                 → 用参考序列作骨架，把主读序校正为"参考等位"，
                                   次读序校正为"非参考等位"

术语
----
primary   参考等位（对应图上"目标序列"里与 ref 相同的那条）
secondary 非参考等位
两者相等即纯合/野生，不等即杂合。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from .align import IUPAC_CODE_MAP, LocalAlignmentResult, _revcomp, local_align, merge_iupac

_BASES = ("A", "C", "G", "T")


class SangerseqError(RuntimeError):
    pass


@dataclass
class Sangerseq:
    """对应 sangerseqR 的 sangerseq 对象。"""

    trace: np.ndarray  # (n, 4) 列序 A,C,G,T
    peak_pos: np.ndarray  # (n, 4) 峰位（trace 下标 + 1）
    peak_amp: np.ndarray  # (n, 4) 峰高
    primary_seq: str
    secondary_seq: str
    primary_id: str = ""
    secondary_id: str = ""


# --------------------------------------------------------------------------
# getpeaks
# --------------------------------------------------------------------------
def getpeaks(trace: np.ndarray) -> np.ndarray:
    """等价 sangerseqR::getpeaks —— 返回局部极大值 (索引(1-based), 值)。

    注意原实现用 `rep(..., times = rle$lengths)`，若出现"平台"
    （相邻相等值）则该平台上的每个点都会被算作峰；这里刻意保持同样行为。
    """
    values: List[float] = []
    lengths: List[int] = []
    for v in trace:
        if values and v == values[-1]:
            lengths[-1] += 1
        else:
            values.append(v)
            lengths.append(1)
    if not values:
        return np.zeros((0, 2))

    ext = [float("-inf")] + list(values) + [float("-inf")]
    d1 = np.diff(np.asarray(ext, dtype=float))
    d2 = np.diff(np.sign(d1))
    peak_runs = d2 == -2
    expanded = np.repeat(peak_runs, lengths)
    idx = np.nonzero(expanded)[0] + 1  # 1-based
    return np.column_stack([idx.astype(float), np.asarray(trace, dtype=float)[idx - 1]])


def _peakvalues(peaks: np.ndarray, pstart: float, pstop: float) -> Tuple[float, float]:
    """窗口内最大峰：(值, 位置)；窗口为空返回 (0, nan)。"""
    if peaks.shape[0] == 0:
        return 0.0, float("nan")
    mask = (peaks[:, 0] > pstart) & (peaks[:, 0] < pstop)
    region = peaks[mask]
    if region.shape[0] == 0:
        return 0.0, float("nan")
    k = int(np.argmax(region[:, 1]))  # 与 which.max 一致：取第一个最大值
    return float(region[k, 1]), float(region[k, 0])


# --------------------------------------------------------------------------
# makeBaseCalls
# --------------------------------------------------------------------------
def make_base_calls(obj: Sangerseq, ratio: float = 0.33) -> Sangerseq:
    """等价 sangerseqR::makeBaseCalls。"""
    peaks = [getpeaks(obj.trace[:, c]) for c in range(4)]
    primarypeaks = np.asarray(obj.peak_pos[:, 0], dtype=float)
    n = len(primarypeaks)
    if n == 0:
        raise SangerseqError("peakPosMatrix 为空，无法 basecall")

    diffs = np.diff(np.concatenate([[0.0], primarypeaks]))
    starts = primarypeaks - 0.5 * diffs
    stops = np.concatenate([primarypeaks[:-1] + 0.5 * diffs[1:], [primarypeaks[-1] + 0.5 * diffs[-1]]])

    temp_pos = np.full((n, 4), np.nan)
    temp_amp = np.full((n, 4), np.nan)
    primary: List[str] = []
    secondary: List[str] = []

    for i in range(n):
        amps = np.zeros(4)
        poss = np.full(4, np.nan)
        empty = True
        for c in range(4):
            a, p = _peakvalues(peaks[c], starts[i], stops[i])
            amps[c] = a
            poss[c] = p
            if not np.isnan(p):
                empty = False
        if empty:
            # 与 R 的 `next` 一致：不留 basecall，该行保持 NA 并最终被丢弃
            continue
        temp_amp[i] = amps
        temp_pos[i] = poss

        mx = amps.max()
        if mx > 0:
            keep = (amps / mx) >= ratio
        else:
            keep = np.zeros(4, dtype=bool)
        letters: List[Optional[str]] = [b if k else None for b, k in zip(_BASES, keep)]
        # order(signals, decreasing=TRUE)：稳定排序，并列时保持 A,C,G,T 顺序
        order = np.argsort(-amps, kind="stable")
        letters = [letters[k] for k in order]

        n_ok = sum(1 for x in letters if x is not None)
        if n_ok == 4 or n_ok == 0:
            primary.append("N")
            secondary.append("N")
        elif n_ok > 1:
            primary.append(letters[0])
            rest = sorted(x for x in letters[1:4] if x is not None)
            secondary.append(merge_iupac(["".join(rest)])[0])
        else:
            primary.append(letters[0])
            secondary.append(letters[0])

    keep_rows = ~np.all(np.isnan(temp_pos), axis=1)
    return Sangerseq(
        trace=obj.trace,
        peak_pos=temp_pos[keep_rows],
        peak_amp=temp_amp[keep_rows],
        primary_seq="".join(primary),
        secondary_seq="".join(secondary),
        primary_id="sangerseq package primary basecalls",
        secondary_id="sangerseq package secondary basecalls",
    )


# --------------------------------------------------------------------------
# basecalldf
# --------------------------------------------------------------------------
@dataclass
class BasecallTable:
    primary: List[Optional[str]]
    secondary: List[Optional[str]]
    consensus: List[str]
    possibilities: List[Optional[str]]
    ref: List[str] = field(default_factory=list)
    noref: List[str] = field(default_factory=list)
    newprimary: List[str] = field(default_factory=list)
    newsecondary: List[str] = field(default_factory=list)


def basecalldf(obj: Sangerseq) -> BasecallTable:
    """等价 sangerseqR:::basecalldf（含 data.frame 的长度循环规则）。"""
    p = list(obj.primary_seq.upper())
    s = list(obj.secondary_seq.upper())
    n = max(len(p), len(s))
    if len(p) == 0 or len(s) == 0:
        raise SangerseqError("primary/secondary 序列为空")
    # R 的 data.frame 会把短的那个循环补齐
    p = [p[i % len(p)] for i in range(n)]
    s = [s[i % len(s)] for i in range(n)]

    primary = [IUPAC_CODE_MAP.get(ch) for ch in p]
    secondary = [IUPAC_CODE_MAP.get(ch) for ch in s]
    if any(x is None for x in primary) or any(x is None for x in secondary):
        bad = [ch for ch, x in zip(p + s, primary + secondary) if x is None]
        raise SangerseqError("序列中含非 IUPAC 字母: %s" % sorted(set(bad)))

    consensus = list(primary)
    for i in range(n):
        if primary[i] != secondary[i] or len(consensus[i]) > 1:
            consensus[i] = "N"
    possibilities = list(consensus)
    for i in range(n):
        if possibilities[i] == "N":
            possibilities[i] = primary[i] + secondary[i]
    return BasecallTable(primary=primary, secondary=secondary, consensus=consensus, possibilities=possibilities)


# --------------------------------------------------------------------------
# setAllelePhase
# --------------------------------------------------------------------------
def set_allele_phase(
    obj: Sangerseq, refseq: str, trim5: int = 0, trim3: int = 0
) -> Tuple[Sangerseq, BasecallTable, LocalAlignmentResult, str]:
    """等价 sangerseqR::setAllelePhase。

    返回 (校正后的 Sangerseq, basecall 表, 比对结果, 参考链方向)。
    """
    table = basecalldf(obj)
    n = len(table.consensus)
    if trim5 + trim3 >= n:
        raise SangerseqError("trim5+trim3 超过读长")
    seedseq = "".join(table.consensus[trim5 : n - trim3])

    # 与原 R 一致：pairwiseAlignment(seedseq, ref) 与 (seedseq, revcomp(ref)) 取分高者。
    # 注意 R 选中 revcomp 后会把 refseq 整体换成反向互补，坐标也就落在新的 refseq 上，
    # 这里必须保持同样的语义（不要映射回原坐标）。
    pa = local_align(seedseq, refseq)
    pa_rc = local_align(seedseq, _revcomp(refseq))
    if pa_rc.score > pa.score:
        chosen = pa_rc
        refstrand = "Reference (revcomp)"
        refvector = list(_revcomp(refseq))
    else:
        chosen = pa
        refstrand = "Reference"
        refvector = list(refseq.upper())

    if chosen.query_width < 10:
        raise SangerseqError(
            "种子序列太短：至少需要 10bp 的良好匹配序列（匹配到 %d bp）" % chosen.query_width
        )

    seqstart = chosen.query_start
    refstart = chosen.subject_start
    startoffset = refstart - seqstart - trim5
    if startoffset < 0:
        refvector = ["N"] * abs(startoffset) + refvector
        startoffset = 0
    end = n + startoffset
    if end > len(refvector):
        refvector = refvector + ["N"] * (end - len(refvector))

    table.ref = refvector[startoffset:end]

    # noref：把 ref 字符（作为正则字符组）从 possibilities 里删掉
    table.noref = []
    for i in range(n):
        r = table.ref[i]
        pat = re.compile("[%s]" % re.escape(r))
        table.noref.append(pat.sub("", table.possibilities[i]))

    table.newprimary = [""] * n
    table.newsecondary = [""] * n
    for i in range(n):
        ref_not_option = len(table.noref[i]) == len(table.possibilities[i])
        if ref_not_option:
            table.newprimary[i] = table.possibilities[i]
            table.newsecondary[i] = table.possibilities[i]
        else:
            table.newprimary[i] = table.ref[i]
            table.newsecondary[i] = table.noref[i]
        if len(table.newsecondary[i]) == 0:
            table.newsecondary[i] = table.ref[i]

    primary_seq = "".join(merge_iupac(table.newprimary))
    secondary_seq = "".join(merge_iupac(table.newsecondary))

    out = Sangerseq(
        trace=obj.trace,
        peak_pos=obj.peak_pos,
        peak_amp=obj.peak_amp,
        primary_seq=primary_seq,
        secondary_seq=secondary_seq,
        primary_id=refstrand,
        secondary_id="NonReference",
    )
    return out, table, chosen, refstrand
