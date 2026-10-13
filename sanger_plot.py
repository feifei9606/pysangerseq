#!/usr/bin/env python3
"""一代测序（Sanger）自动化峰图绘制 —— 命令行入口。

用法（与原 R 脚本的前两个参数顺序保持一致）：
    ./sanger_plot.py 12-14796923-G-A C7197-...-F_....ab1
    ./sanger_plot.py 12-14796923-G-A sample.ab1 --genome hg19 -o out.png
    ./sanger_plot.py 12-14796923-G-A sample-R.ab1 --strand reverse   # 反向测序
    ./sanger_plot.py 6-26091000-CTTGGT-CTT sample.ab1 --fasta /path/ref.fa

默认输出为 `<ab1 文件>.png`（与原脚本一致）。

注意：正向/反向**不再从文件名推断**。以前文件名里带 `-R...` 会被自动镜像，
现在必须显式传 `--strand reverse`（或 `--strand R`）。
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sanger_plot.pipeline import (  # noqa: E402
    DEFAULT_FASTA,
    STRAND_META,
    PipelineError,
    parse_strand,
    run,
)
from sanger_plot.render import DEFAULT_CJK_FONT, RenderOptions  # noqa: E402


def _strand_arg(value: str) -> str:
    """argparse 的校验钩子：非法方向直接在参数解析阶段报错，并原样保留入参用于回显。"""
    try:
        parse_strand(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc))
    return value


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sanger_plot.py",
        description="一代测序峰图绘制（Python 重构版，无 R / swalign / makeSCF 依赖）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "变异描述格式：chr-pos-ref-alt 或 chr:pos:ref:alt（如 12-14796923-G-A）\n"
            "--strand 取值：%s（大小写不敏感；默认 forward）\n"
            "注意：正向/反向不再从文件名推断，反向样本必须显式传 --strand reverse" % STRAND_META
        ),
    )
    p.add_argument("variant", help="变异描述，如 12-14796923-G-A")
    p.add_argument("ab1", help=".ab1 文件路径")
    p.add_argument("--genome", default="hg38", choices=sorted(DEFAULT_FASTA),
                   help="参考基因组（默认 hg38）")
    p.add_argument("--fasta", default=None, help="自定义参考基因组 fasta（优先级高于 --genome）")
    p.add_argument("-o", "--out", default=None, help="输出 png（默认 <ab1>.png）")
    p.add_argument("--strand", default="forward", metavar=STRAND_META, type=_strand_arg,
                   help="测序方向，默认 forward（反向测序需显式指定 reverse 或 R）")
    p.add_argument("--dlen", type=int, default=21, help="展示窗口宽度（奇数，默认 21）")
    p.add_argument("--ratio", type=float, default=0.3, help="basecall 次峰比例阈值（默认 0.3）")
    p.add_argument("--phase-trim5", type=int, default=29, help="setAllelePhase trim5（默认 29）")
    p.add_argument("--phase-trim3", type=int, default=30, help="setAllelePhase trim3（默认 30）")

    p.add_argument("--snp-max-diffs", type=int, default=3,
                   help="差异位点 <= 该值判为 SNP 并逐个加框，否则按 INDEL 展示（默认 3）")
    p.add_argument("--cex-base", type=float, default=0.8, help="碱基字母相对字号（默认 0.8）")
    p.add_argument("--cex-mtext", type=float, default=0.6,
                   help='"参考序列/目标序列" 相对字号（默认 0.6）')
    p.add_argument("--letter-scale", type=float, default=1.13,
                   help="碱基字母高度标定系数（换字体后需重新标定；严格照搬 R 用 1.0，默认 1.13）")
    p.add_argument("--showhets", action="store_true",
                   help="用浅蓝底色标出杂合峰区间（原脚本 showhets 选项，默认关）")
    p.add_argument("--showtrim", action="store_true",
                   help="用红色阴影标出被 trim 掉的区间（原脚本 showtrim 选项，默认关）")

    p.add_argument("--width-px", type=int, default=2000, help="画布宽（像素，默认 2000）")
    p.add_argument("--height-px", type=int, default=1000, help="画布高（像素，默认 1000）")
    p.add_argument("--dpi", type=int, default=300, help="dpi（默认 300）")
    p.add_argument("--base-font", default=None, help="碱基字母字体文件")
    p.add_argument("--cjk-font", default=DEFAULT_CJK_FONT, help="中文字体文件")
    p.add_argument("--dry-run", action="store_true", help="只跑流程不落盘，打印关键参数")
    p.add_argument("-q", "--quiet", action="store_true", help="不打印结果摘要")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    opts = RenderOptions(
        snp_max_diffs=args.snp_max_diffs,
        cex_base=args.cex_base,
        cex_mtext=args.cex_mtext,
        letter_scale=args.letter_scale,
        showhets=args.showhets,
        showtrim=args.showtrim,
        width_px=args.width_px,
        height_px=args.height_px,
        dpi=args.dpi,
        base_font=args.base_font,
        cjk_font=args.cjk_font,
    )
    try:
        res = run(
            ab1_path=args.ab1,
            variant_spec=args.variant,
            genome=args.genome,
            fasta_path=args.fasta,
            outfile=args.out,
            dlen=args.dlen,
            ratio=args.ratio,
            phase_trim5=args.phase_trim5,
            phase_trim3=args.phase_trim3,
            strand=args.strand,
            opts=opts,
            write=not args.dry_run,
        )
    except PipelineError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2
    except ValueError as exc:  # 例如 ab1 不是 ABIF、变异描述格式不对
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2

    # 差异位点不在正中 → 打提示（不改退出码，便于批量跑时在日志里发现问题）
    warn = res.center_warning()
    if warn:
        print("WARN: %s" % warn, file=sys.stderr)

    if args.quiet:
        return 0

    r = res.render
    print("变异            : %s" % res.variant.describe())
    print("校正后断点      : %s:%d" % (res.variant.chrom, res.altposs))
    print("读长            : %d" % res.read_length)
    print("trim5 / trim3   : %d / %d" % (res.trim5len, res.trim3len))
    print("展示参考窗口    : %s" % res.refseq1)
    print("读序链方向      : %s" % ("-" if res.posref.rc else "+"))
    print("测序方向        : %s（--strand 入参）" % res.strand)
    print("面板数          : %d" % r.n_plots)
    print("差异位点(1-based): %s" % [k + 1 for k in r.diffindex])
    print("加框位点        : %s" % [k + 1 for k in r.box_indices])
    print("箭头 x          : %s" % (None if r.arrow_pos is None else round(r.arrow_pos, 2)))
    print("输出            : %s" % (res.outfile or "(dry-run，未落盘)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
