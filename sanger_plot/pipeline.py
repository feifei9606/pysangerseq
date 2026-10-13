"""端到端流程编排 —— 对应原 R 脚本从 args 解析到出图的全过程。

与原脚本的差异（都是有意的）
----------------------------
* 参考子串直接从 FASTA 取，不再走 `samtools faidx` + `makeSCF` + `read.scf`
  的"伪造 SCF"绕路（那只是为了在 R 里拿到子串）。
* `swalign` 调用换成内置实现（swalign_local），不再需要 python2。
* 原脚本写作 `|sed -n 2p` 取 swalign 输出，实测取到的是 dump 的第 2 行
  "Ref  : ..."（只有 5 个字段，`posref2[7]` 为 NA 直接报错），
  实际应为最后一行数据行；这里按语义取数据行。
* **测序方向改为显式入参**：原脚本用文件名正则（`R_xxx.ab1` / `-R...`）猜正反向，
  现在由 `run(strand=...)` / `--strand` 决定，默认 forward。反向样本必须显式传
  `reverse`（或 `R`）才能得到与旧流程一致的镜像图。
  注意与 `posref.rc` 那处翻转是两回事：后者是"读序比对到参考反向互补"，数据驱动，
  与原脚本一致，保持不动。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Union

from .abif import read_abif
from .align import _revcomp
from .reference import Fasta
from .render import RenderOptions, RenderResult, chromatogram
from .sangerseq import make_base_calls, set_allele_phase
from .swalign_local import Alignment, summary as sw_summary
from .variant import Variant, adjust_position, parse_variant

# 与原脚本一致的默认参考基因组
DEFAULT_FASTA = {
    "hg19": "/sg2/19.yuli/02WES/Sangerseq/human_g1k_v37.fasta",
    "hg38": "/sg2/19.yuli/02WES/Sangerseq/Homo_sapiens_assembly38.fasta",
}

# 测序方向。原脚本是从文件名正则（`R_xxx.ab1` / `-R...`）猜的，
# 这里改成显式入参：CLI / 库调用 / 测试三处都走 parse_strand 归一化。
STRAND_FORWARD_ALIASES = frozenset({"forward", "f", "fwd", "+"})
STRAND_REVERSE_ALIASES = frozenset({"reverse", "r", "rev", "-"})
STRAND_META = "forward|reverse|F|R"


class PipelineError(RuntimeError):
    pass


def parse_strand(value: Union[str, bool]) -> bool:
    """归一化测序方向，返回 True 表示反向（需要镜像）。

    接受 forward/f/fwd/+ 与 reverse/r/rev/-（大小写不敏感，前后空白忽略）；
    也接受 bool（True 视为 reverse），方便库调用。
    无法识别时抛 ValueError，消息里列出全部可接受取值。
    """
    if isinstance(value, bool):
        return value
    if not isinstance(value, str):
        raise ValueError(
            "测序方向必须是字符串或布尔值，收到 %r（可接受: %s）" % (value, STRAND_META)
        )
    key = value.strip().lower()
    if key in STRAND_FORWARD_ALIASES:
        return False
    if key in STRAND_REVERSE_ALIASES:
        return True
    raise ValueError(
        "无法识别的测序方向 %r（可接受: %s，大小写不敏感）" % (value, STRAND_META)
    )


def strand_name(reverse: bool) -> str:
    """把布尔方向还原成规范名，用于日志与测试断言。"""
    return "reverse" if reverse else "forward"


@dataclass
class RunResult:
    outfile: Optional[str]
    variant: Variant
    altposs: int
    read_length: int
    trim5len: int
    trim3len: int
    refseq1: str
    reverse: bool
    primary_seq: str
    secondary_seq: str
    posref2: Alignment
    posref: Alignment
    render: RenderResult
    # 有差异、但差异位点不在展示窗口正中（第 dlen2+1 位）→ 可疑
    center_mismatch: bool

    @property
    def strand(self) -> str:
        """归一化后的方向名（"forward" / "reverse"）。"""
        return strand_name(self.reverse)

    @property
    def center_pos_1based(self) -> int:
        """展示窗口正中那一位的 1-based 序号（dlen=21 时为 11）。"""
        return len(self.render.lab3) // 2 + 1

    def center_warning(self) -> Optional[str]:
        """差异位点不在正中时返回一条提示文案，否则 None。

        正中位点是"断点"的定义位置（`altposs`），箭头也只画在那里；
        一旦偏移，通常意味着坐标/ref-alt/`--strand` 有一样不对，
        或者锚定比对落到了错的位置。
        """
        if not self.center_mismatch:
            return None
        return (
            "差异位点在第 %s 位，不在展示窗口正中（第 %d 位）——"
            "请核对变异坐标 / ref-alt / --strand；"
            "若样本在该位点本就是野生型可忽略此提示"
            % ("、".join(str(k + 1) for k in self.render.diffindex),
               self.center_pos_1based)
        )


def run(
    ab1_path: str,
    variant_spec: str,
    genome: str = "hg38",
    fasta_path: Optional[str] = None,
    outfile: Optional[str] = None,
    dlen: int = 21,
    ratio: float = 0.3,
    phase_trim5: int = 29,
    phase_trim3: int = 30,
    strand: Union[str, bool] = "forward",
    opts: Optional[RenderOptions] = None,
    write: bool = True,
) -> RunResult:
    """跑完整流程并出图。

    ratio / phase_trim5 / phase_trim3 对应原脚本里写死的
    makeBaseCalls(ratio = 0.3) 与 setAllelePhase(trim5 = 29, trim3 = 30)。

    strand 为测序方向（默认 "forward"），决定出图是否镜像。原脚本是从文件名
    正则推断的，现在必须显式给出：反向样本要传 "reverse"（或 "R"）。
    """
    if genome not in DEFAULT_FASTA and fasta_path is None:
        raise PipelineError("未知基因组 %r，且未提供 fasta" % genome)
    fa_path = fasta_path or DEFAULT_FASTA[genome]

    try:
        reverse = parse_strand(strand)
    except ValueError as exc:
        raise PipelineError(str(exc))

    if dlen < 1:
        raise PipelineError("dlen 必须为正整数，收到 %r" % dlen)
    if dlen % 2 == 0:
        dlen += 1
    dlen2 = (dlen - 1) // 2

    variant = parse_variant(variant_spec)

    base = read_abif(ab1_path)

    fa = Fasta(fa_path)
    try:
        chrom = _resolve_chrom(fa, variant.chrom)
        altposs = variant.pos
        chrom_len = fa.chrom_length(chrom)
        if altposs < 1 or altposs > chrom_len:
            raise PipelineError(
                "变异位点 %s:%d 超出染色体范围（%s 全长 %d）"
                % (chrom, altposs, chrom, chrom_len)
            )

        # ---- 1. 参考子串：原为 samtools faidx(±300) → makeSCF → read.scf → subseq(30,500)
        region = fa.fetch(chrom, altposs - 300, altposs + 300)
        ref500 = region[29:529]
        if len(ref500) != 500:
            raise PipelineError(
                "参考骨架只取到 %d bp（需要 500 bp，来自 %s:%d-%d）；"
                "位点过于靠近染色体/contig 边界"
                % (len(ref500), chrom, altposs - 300, altposs + 300)
            )

        # ---- 2. basecall + 相位校正
        called = make_base_calls(base, ratio=ratio)
        phased, _table, _pa, _strand = set_allele_phase(
            called, ref500, trim5=phase_trim5, trim3=phase_trim3
        )
        read_seq = phased.primary_seq
        alt_seq = phased.secondary_seq
        read_len = len(alt_seq)

        # ---- 3. 第一次 swalign：判断链方向 + 探测断点附近是否有 INDEL
        posref2 = locate_anchor(read_seq, ref500, dlen)

        # ---- 4. 变异位点校正
        # variant.ref / variant.alt 为 None 时（只给位点）adjust_position 不做校正，
        # 与原脚本的 "ref"/"alt" 哨兵语义一致。
        altposs = adjust_position(
            chrom, altposs, variant.ref, variant.alt,
            "-" if posref2.rc else "+", fa,
        )

        # ---- 5. 第二次 swalign：定位读序，算 trim
        posref, trim5len, trim3len = locate_and_trim(
            read_seq, chrom, altposs, dlen, posref2.rc, fa
        )

        # ---- 6. 展示用参考窗口
        # 注意这里有两处互相独立的翻转，别混：
        #   a) posref.rc —— 读序比对到参考的反向互补，数据驱动，必须保留；
        #   b) reverse   —— 反向测序（以前由文件名正则推断，现在由 --strand 显式给出）。
        refseq1 = fa.fetch(chrom, altposs - dlen2, altposs + dlen2)
        if len(refseq1) != dlen:
            raise PipelineError(
                "展示参考窗口只取到 %d bp（需要 %d bp，来自 %s:%d-%d）"
                % (len(refseq1), dlen, chrom, altposs - dlen2, altposs + dlen2)
            )
        if posref.rc:
            refseq1 = _revcomp(refseq1)
        if reverse:
            refseq1 = _revcomp(refseq1)

        if trim5len < 0 or trim3len < 0:
            raise PipelineError(
                "trim5len/trim3len 出现负值 (%d/%d)，说明断点定位失败"
                "（原脚本同样会异常）；请检查变异位点是否与读序匹配"
                % (trim5len, trim3len)
            )

        if outfile is None and write:
            outfile = ab1_path + ".png"

        rr = chromatogram(
            phased, refseq1,
            trim5=trim5len, trim3=trim3len,
            outfile=outfile if write else None,
            dlen=dlen,
            reverse=reverse,
            opts=opts or RenderOptions(),
        )

        # 有差异但不在正中 → 可疑（正中位点即 altposs，箭头也只画在那里）
        center_mismatch = bool(rr.diffindex) and (dlen2 not in rr.diffindex)

        return RunResult(
            outfile=outfile if write else None,
            variant=variant,
            altposs=altposs,
            read_length=read_len,
            trim5len=trim5len,
            trim3len=trim3len,
            refseq1=refseq1,
            reverse=reverse,
            primary_seq=phased.primary_seq,
            secondary_seq=phased.secondary_seq,
            posref2=posref2,
            posref=posref,
            render=rr,
            center_mismatch=center_mismatch,
        )
    finally:
        fa.close()


def _resolve_chrom(fa: Fasta, chrom: str) -> str:
    """染色体名解析：带/不带 `chr` 前缀都容忍，找不到时给出可读报错。"""
    if fa.has_chrom(chrom):
        return chrom
    alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    if fa.has_chrom(alt):
        return alt
    sample = ", ".join(fa.chroms()[:5])
    raise PipelineError(
        "参考基因组中没有染色体 %r（也尝试过 %r）；该文件前几个序列名为: %s ..."
        % (chrom, alt, sample)
    )


def _has_indel(cigar) -> bool:
    """对应原脚本的 grepl("[D|I]", posref[8])。

    注意 swalign 的 Alignment.cigar 是 [(count, op), ...] 这样的列表，
    而 CIGAR 字符串在 .cigar_str 里 —— 两种都接受，避免调用点写错。
    """
    if isinstance(cigar, str):
        return ("D" in cigar) or ("I" in cigar)
    return any(op in ("D", "I") for _, op in cigar)


def locate_anchor(read_seq: str, ref500: str, dlen: int) -> Alignment:
    """第一次 swalign：用 500bp 参考骨架上的锚定片段判读序链方向 + 探测 INDEL。

    对应原脚本第 104-108 行：先取 ref 的第 251..271 位（21bp）当锚定序列；
    若比对结果含 `D/I`（说明这一段正好压在变异上），就换成第 273..293 位重来一次。

    该函数由 pipeline 与回归测试共用，避免两处各写一份而漂移。
    """
    if len(ref500) < 293:
        raise PipelineError(
            "参考骨架只有 %d bp，无法取到锚定片段（需要 >= 293 bp）" % len(ref500)
        )
    posref2 = sw_summary(read_seq, ref500[250:271])
    if _has_indel(posref2.cigar):
        posref2 = sw_summary(read_seq, ref500[272:293])
    return posref2


def locate_and_trim(
    read_seq: str,
    chrom: str,
    altposs: int,
    dlen: int,
    rev_strand: bool,
    fasta: Fasta,
):
    """用断点旁 2*dlen2 bp 的锚定片段定位读序，算出 trim5/trim3。

    对应原脚本第 243-289 行：默认取断点"上游"的锚定片段（负链时取下游），
    若比对结果含 INDEL（说明锚定片段正好压在变异上），就换到另一侧再算一次，
    并把 trim5len 的公式从 `+dlen2` 换成 `-dlen2-1`。

    返回 (posref, trim5len, trim3len)。
    """
    dlen2 = (dlen - 1) // 2
    read_len = len(read_seq)

    def chunk(minus_side: bool):
        if minus_side:
            return altposs + 1, altposs + dlen2 + dlen2
        return altposs - dlen2 - dlen2, altposs - 1

    posref = sw_summary(read_seq, fasta.fetch(chrom, *chunk(rev_strand)))
    trim5len = min(posref.r_pos, posref.r_end) + dlen2
    trim3len = read_len - trim5len - dlen
    if _has_indel(posref.cigar):
        posref = sw_summary(read_seq, fasta.fetch(chrom, *chunk(not rev_strand)))
        trim5len = min(posref.r_pos, posref.r_end) - dlen2 - 1
        trim3len = read_len - trim5len - dlen
    return posref, trim5len, trim3len

