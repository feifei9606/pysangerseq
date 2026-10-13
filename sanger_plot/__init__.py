"""一代测序（Sanger）自动化峰图绘制 —— 纯 Python 实现。

本包是 `sangerseq-refalt_AIO.r` 的完整重构，去掉了 R / sangerseqR /
Biostrings / Cairo / makeSCF / swalign(python2) 全部外部依赖，
只保留 Python + numpy + matplotlib。

模块划分
--------
abif               ABIF(.ab1) 文件解析
reference          FASTA(.fai) 随机访问取序列
sangerseq          basecalling / 等位基因相位校正（对应 sangerseqR）
swalign_local      Smith-Waterman 局部比对（对应 swalign 命令行）
align              Biostrings pairwiseAlignment 等价实现
variant            变异描述解析 + INDEL 断点位置校正
render             峰图绘制
pipeline           端到端流程编排
"""

from .abif import read_abif, AbifError
from .reference import Fasta
from .sangerseq import Sangerseq, make_base_calls, set_allele_phase
from .variant import Variant, parse_variant, adjust_position
from .render import chromatogram, RenderOptions
from .pipeline import (
    DEFAULT_FASTA,
    STRAND_META,
    PipelineError,
    RunResult,
    locate_anchor,
    locate_and_trim,
    parse_strand,
    run,
    strand_name,
)

__all__ = [
    "read_abif",
    "AbifError",
    "Fasta",
    "Sangerseq",
    "make_base_calls",
    "set_allele_phase",
    "Variant",
    "parse_variant",
    "adjust_position",
    "chromatogram",
    "RenderOptions",
    "DEFAULT_FASTA",
    "STRAND_META",
    "PipelineError",
    "RunResult",
    "locate_anchor",
    "locate_and_trim",
    "parse_strand",
    "run",
    "strand_name",
]

__version__ = "1.1.0"
