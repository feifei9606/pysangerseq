"""局部比对（等价 Biostrings::pairwiseAlignment(type="local")）。

原 R 调用：
    pairwiseAlignment(seedseq, ref, type="local",
                      gapOpening=-200, gapExtension=-10)

实测行为
--------
* R 的调用没有传 substitutionMatrix，Biostrings 会用内部默认矩阵。
  在本机 Biostrings 2.56 上实测该矩阵为 匹配=1.981756 / 错配=-5.899286
  （比值约 1:-2.98）。这里取等价的整数比 **匹配=1 / 错配=-3**。
  已验证：本流程只用到"比对区间"，而该区间在
  {默认矩阵, 1:-1, 1:-3, 1:0} 四种打分下完全一致（见 tests/），
  所以该简化不影响结果。
* gapOpening=-200 / gapExtension=-10 为负，实际表现为"极大的空位惩罚"，
  实测比对结果不含任何 gap；这里按 惩罚 = gap_open + gap_ext*L 实现，
  等价于原行为。
* 返回值用 Biostrings 的坐标语义：1-based、闭区间；width 为比对列数。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

# Biostrings::IUPAC_CODE_MAP
IUPAC_CODE_MAP: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "M": "AC", "R": "AG", "W": "AT", "S": "CG", "Y": "CT", "K": "GT",
    "V": "ACG", "H": "ACT", "D": "AGT", "B": "CGT", "N": "ACGT",
}
_SORTED_BASES = ("A", "C", "G", "T")
_CODE_TO_LETTER = {v: k for k, v in IUPAC_CODE_MAP.items()}


def merge_iupac(strings) -> list:
    """等价 Biostrings::mergeIUPACLetters：把每段字母并成 IUPAC 代码。

    merge_iupac(["AC", "G"]) -> ["M", "G"]
    """
    out = []
    for s in strings:
        if s is None or len(s) == 0:
            raise ValueError("merge_iupac 收到空串（与原 R 行为一致，属异常输入）")
        bases = set()
        for ch in s.upper():
            expanded = IUPAC_CODE_MAP.get(ch)
            if expanded is None:
                raise ValueError("非 IUPAC 字母: %r" % ch)
            bases.update(expanded)
        key = "".join(b for b in _SORTED_BASES if b in bases)
        if key not in _CODE_TO_LETTER:
            raise ValueError("无法合并为 IUPAC 代码: %r" % s)
        out.append(_CODE_TO_LETTER[key])
    return out


def _score(a: str, b: str, match: float, mismatch: float) -> float:
    """单个位置打分；模糊字母按展开碱基取平均（Biostrings fuzzy 规则）。"""
    ba = IUPAC_CODE_MAP.get(a)
    bb = IUPAC_CODE_MAP.get(b)
    if ba is None or bb is None:
        return mismatch
    if len(ba) == 1 and len(bb) == 1:
        return match if ba == bb else mismatch
    total = 0.0
    for x in ba:
        for y in bb:
            total += match if x == y else mismatch
    return total / (len(ba) * len(bb))


@dataclass
class LocalAlignmentResult:
    score: float
    query_start: int  # 1-based 闭区间
    query_end: int
    subject_start: int
    subject_end: int
    cigar: str

    @property
    def query_width(self) -> int:
        return self.query_end - self.query_start + 1

    @property
    def subject_width(self) -> int:
        return self.subject_end - self.subject_start + 1

    @property
    def has_indel(self) -> bool:
        return ("D" in self.cigar) or ("I" in self.cigar)


def local_align(
    query: str,
    subject: str,
    match: float = 1.0,
    mismatch: float = -3.0,
    gap_open: float = 200.0,
    gap_ext: float = 10.0,
) -> LocalAlignmentResult:
    """仿射空位局部比对（Gotoh）。query 为 pattern，subject 为目标序列。"""
    q = query.upper()
    s = subject.upper()
    n, m = len(q), len(s)
    if n == 0 or m == 0:
        return LocalAlignmentResult(0.0, 0, 0, 0, 0, "")

    neg = float("-inf")
    H = [[0.0] * (m + 1) for _ in range(n + 1)]
    E = [[neg] * (m + 1) for _ in range(n + 1)]
    F = [[neg] * (m + 1) for _ in range(n + 1)]
    tb = [[""] * (m + 1) for _ in range(n + 1)]

    best = 0.0
    bi = bj = 0
    for i in range(1, n + 1):
        qi = q[i - 1]
        for j in range(1, m + 1):
            E[i][j] = max(
                H[i][j - 1] - (gap_open + gap_ext),
                E[i][j - 1] - gap_ext,
            )
            F[i][j] = max(
                H[i - 1][j] - (gap_open + gap_ext),
                F[i - 1][j] - gap_ext,
            )
            diag = H[i - 1][j - 1] + _score(qi, s[j - 1], match, mismatch)
            v = max(0.0, diag, E[i][j], F[i][j])
            if v == 0.0:
                tb[i][j] = ""
            elif v == diag:
                tb[i][j] = "M"
            elif v == F[i][j]:
                tb[i][j] = "I"  # 消耗 query，即 subject 上有空位
            else:
                tb[i][j] = "D"  # 消耗 subject，即 query 上有空位
            H[i][j] = v
            if v > best:
                best = v
                bi, bj = i, j

    # 回溯
    ops = []
    i, j = bi, bj
    while i > 0 and j > 0 and H[i][j] > 0:
        op = tb[i][j]
        if op == "M":
            ops.append("M")
            i -= 1
            j -= 1
        elif op == "I":
            ops.append("I")
            i -= 1
        elif op == "D":
            ops.append("D")
            j -= 1
        else:
            break
    ops.reverse()

    if not ops:
        return LocalAlignmentResult(0.0, 0, 0, 0, 0, "")

    q_start, s_start = i + 1, j + 1
    q_end, s_end = bi, bj

    # CIGAR 归并
    cigar_parts = []
    for op in ops:
        if cigar_parts and cigar_parts[-1][1] == op:
            cigar_parts[-1][0] += 1
        else:
            cigar_parts.append([1, op])
    cigar_str = "".join("%d%s" % (c, o) for c, o in cigar_parts)

    return LocalAlignmentResult(best, q_start, q_end, s_start, s_end, cigar_str)


_COMPLEMENT = str.maketrans("ACGTNRYMKSWBDHV", "TGCANYRMKSWVHDB")


def _revcomp(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]
