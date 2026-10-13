"""Smith-Waterman 局部比对 —— 移植自 `swalign` 库（原脚本的外部命令依赖）。

原 R 脚本通过
    python2 swalign -summary /dev/stdout <ref> <query>
调用 swalign 命令行，输出 8 列：
    query  ref  query_start  query_end  ref_start  ref_end  strand  cigar

这里把该库的 DP / 回溯 / CIGAR 归并逐行移植到 Python 3，保证坐标语义、
打分矩阵（IdentityScoringMatrix：相同给 match，不同给 mismatch）、
gap 计分与并列取后者（`val >= max_val`）等细节完全一致。

`summary()` 对应命令行里 "正负链各跑一遍、取分高者" 的行为。

坐标约定（与 swalign 相同）：0-based、半开区间 [start, end)。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

_COMPLEMENT = str.maketrans("ACGTNMRWSYKVHDXB", "TGCANKYWSRMBDHVX")


def revcomp(seq: str) -> str:
    """反向互补（保留大小写之外统一转大写）。"""
    return seq.upper().translate(_COMPLEMENT)[::-1]


class IdentityScoringMatrix:
    """swalign.NucleotideScoringMatrix —— 简单的相同/不同打分。"""

    def __init__(self, match: int = 1, mismatch: int = -1):
        self.match = match
        self.mismatch = mismatch

    def score(self, one: str, two: str, wildcard=None) -> int:
        if wildcard and (one in wildcard or two in wildcard):
            return self.match
        return self.match if one == two else self.mismatch


def _reduce_cigar(operations: List[str]) -> List[Tuple[int, str]]:
    count = 1
    last = None
    ret: List[Tuple[int, str]] = []
    for op in operations:
        if last and op == last:
            count += 1
        elif last:
            ret.append((count, last.upper()))
            count = 1
        last = op
    if last:
        ret.append((count, last.upper()))
    return ret


def _cigar_str(cigar: List[Tuple[int, str]]) -> str:
    return "".join("%d%s" % (num, op) for num, op in cigar)


@dataclass
class Alignment:
    """swalign.Alignment 的等价体（只保留本流程用得到的字段）。"""

    query: str
    ref: str
    q_pos: int
    r_pos: int
    cigar: List[Tuple[int, str]]
    score: float
    rc: bool = False

    def __post_init__(self) -> None:
        self.q_len = 0
        self.r_len = 0
        self.matches = 0
        self.mismatches = 0
        i, j = self.r_pos, self.q_pos
        for count, op in self.cigar:
            if op == "M":
                self.q_len += count
                self.r_len += count
                for _ in range(count):
                    if self.query[j] == self.ref[i]:
                        self.matches += 1
                    else:
                        self.mismatches += 1
                    i += 1
                    j += 1
            elif op == "I":
                self.q_len += count
                j += count
                self.mismatches += count
            elif op == "D":
                self.r_len += count
                i += count
                self.mismatches += count
        self.q_end = self.q_pos + self.q_len
        self.r_end = self.r_pos + self.r_len

    @property
    def cigar_str(self) -> str:
        return _cigar_str(self.cigar)


def local_align(
    ref: str,
    query: str,
    scoring: IdentityScoringMatrix,
    gap_penalty: float = -1,
    gap_extension_penalty: float = -1,
    rc: bool = False,
) -> Alignment:
    """单方向局部比对。ref/query 的 0-based 索引分别对应 DP 的列/行。"""
    ref = ref.upper()
    query = query.upper()
    nrow = len(query) + 1
    ncol = len(ref) + 1

    # cell = (val, op, runlen)
    score_m = [[0.0] * ncol for _ in range(nrow)]
    op_m = [[" "] * ncol for _ in range(nrow)]
    run_m = [[0] * ncol for _ in range(nrow)]
    for row in range(1, nrow):
        op_m[row][0] = "i"
    for col in range(1, ncol):
        op_m[0][col] = "d"

    max_val = 0.0
    max_row = 0
    max_col = 0

    for row in range(1, nrow):
        for col in range(1, ncol):
            mm_val = score_m[row - 1][col - 1] + scoring.score(query[row - 1], ref[col - 1])

            if op_m[row - 1][col] == "i":
                ins_run = run_m[row - 1][col]
                if score_m[row - 1][col] == 0:
                    ins_val = 0.0  # 起始不罚
                else:
                    ins_val = score_m[row - 1][col] + gap_extension_penalty
            else:
                ins_run = 0
                ins_val = score_m[row - 1][col] + gap_penalty

            if op_m[row][col - 1] == "d":
                del_run = run_m[row][col - 1]
                if score_m[row][col - 1] == 0:
                    del_val = 0.0
                else:
                    del_val = score_m[row][col - 1] + gap_extension_penalty
            else:
                del_run = 0
                del_val = score_m[row][col - 1] + gap_penalty

            cell_val = max(mm_val, del_val, ins_val, 0.0)

            if del_run and cell_val == del_val:
                val = (cell_val, "d", del_run + 1)
            elif ins_run and cell_val == ins_val:
                val = (cell_val, "i", ins_run + 1)
            elif cell_val == mm_val:
                val = (cell_val, "m", 0)
            elif cell_val == del_val:
                val = (cell_val, "d", 1)
            elif cell_val == ins_val:
                val = (cell_val, "i", 1)
            else:
                val = (0.0, "x", 0)

            # 注意：>= 使得并列时取更靠后的格子（与 swalign 一致）
            if val[0] >= max_val:
                max_val = val[0]
                max_row = row
                max_col = col

            score_m[row][col] = val[0]
            op_m[row][col] = val[1]
            run_m[row][col] = val[2]

    row, col = max_row, max_col
    aln: List[str] = []
    while True:
        val = score_m[row][col]
        if val <= 0:
            break
        op = op_m[row][col]
        aln.append(op)
        if op == "m":
            row -= 1
            col -= 1
        elif op == "i":
            row -= 1
        elif op == "d":
            col -= 1
        else:
            break

    aln.reverse()
    cigar = _reduce_cigar(aln)
    return Alignment(
        query=query, ref=ref, q_pos=row, r_pos=col, cigar=cigar, score=max_val, rc=rc
    )


def summary(
    ref: str,
    query: str,
    match: int = 2,
    mismatch: int = -1,
    gap_penalty: float = -1,
    gap_extension_penalty: float = -1,
) -> Alignment:
    """等价于 `swalign -summary`：正负链各比对一次，取分高者（并列取正链）。"""
    scoring = IdentityScoringMatrix(match, mismatch)
    best = None
    for strand in "+-":
        if strand == "-":
            aln = local_align(
                ref,
                revcomp(query),
                scoring,
                gap_penalty,
                gap_extension_penalty,
                rc=True,
            )
        else:
            aln = local_align(ref, query, scoring, gap_penalty, gap_extension_penalty, rc=False)
        if best is None or aln.score > best.score:
            best = aln
    assert best is not None
    return best
