"""峰图绘制 —— 等价于原 R 脚本里的 chromatogramref()。

版面完全按原图反推（2000x1000px @300dpi 的参考图实测）：
  * 绘图区(plt) = 左 6% 下 12% 上 88% 右 100%，即 mar = c(2,2,2,0)
    （1 行 = 0.2in = 60px）
  * R 的 plot() 默认 yaxs="r"，会上下各留 4% 空隙：
    y 轴实际范围 = (-0.04*yl, 1.04*yl)，所以 y=0 落在绘图区内侧
  * 参考序列字母贴在绘图区上沿外 ~0.5 行，目标序列字母贴在
    下沿外 ~0.5 行；方框 -> y ∈ [-0.14*yl, -0.03*yl]
  * 字母位置按等间距铺开（原脚本刻意这样做，见 pos <- 1:n*peakwidth - peakwidth/2）
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from .sangerseq import Sangerseq  # noqa: E402

# R 的 "green"/"blue"/"black"/"red"（正链 A,C,G,T 顺序），负链颜色整体反转
_FWD_COLORS = ("green", "blue", "black", "red")
_REV_COLORS = ("red", "black", "blue", "green")
_RGB = {"green": "#00FF00", "blue": "#0000FF", "black": "#000000", "red": "#FF0000"}

_BASE_LETTERS = ("A", "C", "G", "T")

# 默认字体（可用环境变量覆盖）
DEFAULT_CJK_FONT = "/sg2/7.hudahui/CNV/cnv_master/cnv/src/utility/fonts/SourceHanSansSC-Regular.otf"


@dataclass
class RenderOptions:
    snp_max_diffs: int = 3
    cex_mtext: float = 0.6
    cex_base: float = 0.8
    letter_scale: float = 1.13
    ylim: float = 3.0
    showhets: bool = False
    showtrim: bool = False
    width_px: int = 2000
    height_px: int = 1000
    dpi: int = 300
    pointsize: float = 12.0
    base_font: Optional[str] = None
    cjk_font: Optional[str] = None


@dataclass
class RenderResult:
    """绘制过程中的关键中间量，便于自动化比对。"""

    n_plots: int
    tracewidth: float
    ylims: tuple
    aveposition: np.ndarray
    pos: np.ndarray
    lab1: List[str]
    lab2: List[str]
    lab3: List[str]
    diffindex: List[int]
    box_indices: List[int]
    arrow_pos: Optional[float]
    plotrange: tuple
    breaks: List[float]


def _row_means_nan(m: np.ndarray) -> np.ndarray:
    """等价 R 的 rowMeans(m, na.rm=TRUE)：全是 NA 的行返回 NaN。"""
    cnt = np.sum(~np.isnan(m), axis=1)
    s = np.nansum(m, axis=1)
    out = np.full(m.shape[0], np.nan)
    nz = cnt > 0
    out[nz] = s[nz] / cnt[nz]
    return out


def chromatogram(
    obj: Sangerseq,
    refseq1: str,
    trim5: int = 0,
    trim3: int = 0,
    outfile: Optional[str] = None,
    dlen: int = 21,
    reverse: bool = False,
    opts: Optional[RenderOptions] = None,
) -> RenderResult:
    """绘制峰图。outfile 为 None 时不落盘（仅计算）。

    reverse=True 表示反向测序：曲线镜像、碱基字母反向互补、A/C/G/T 配色整体反转。
    该值由调用方显式给出（来自 `--strand`），渲染层不再从文件名推断。
    """
    opts = opts or RenderOptions()

    traces = np.asarray(obj.trace, dtype=float)
    bc1 = list(obj.primary_seq)
    bc2 = list(obj.secondary_seq)
    ref = [c.upper() for c in refseq1]

    cols = ["A", "C", "G", "T"]
    peak_pos = np.asarray(obj.peak_pos, dtype=float)

    # ---- aveposition（原脚本：行均值，再优先取 basecalls2 对应通道的峰位） ----
    ave = _row_means_nan(peak_pos)
    for i in range(len(ave)):
        if i >= len(bc2):
            break
        c = bc2[i]
        if c in cols:
            v = peak_pos[i, cols.index(c)]
            if not np.isnan(v):
                ave[i] = v

    n0 = len(ave)
    bc1 = bc1[:n0]
    bc2 = bc2[:n0]

    # ---- 按 trim 截取 ----
    if not opts.showtrim:
        if trim5 + trim3 > len(bc1):
            bc1, bc2, ave = [], [], np.array([])
        else:
            end = len(bc1) - trim3
            bc1 = bc1[trim5:end]
            bc2 = bc2[trim5:end]
            ave = ave[trim5:end]
    bc1_arr = np.array(bc1, dtype=object)
    bc2_arr = np.array(bc2, dtype=object)
    if len(ave) == 0:
        raise ValueError("trim 之后没有任何碱基，无法绘图")

    # ---- trace 裁剪与坐标平移 ----
    if trim3 is not None:
        n2 = int(min(np.nanmax(ave) + 10, traces.shape[0]))
        traces = traces[: max(n2, 1)]
    if trim5 is not None:
        offset = int(max(1, ave[0] - 10))
        traces = traces[offset - 1 :]
        ave = ave - (offset + 3)

    maxsignal = traces.max(axis=1)
    q75 = float(np.percentile(maxsignal, 75))
    q25 = float(np.percentile(maxsignal, 25))
    yl = q75 + opts.ylim * (q75 - q25)
    ylims = (0.0, yl)

    ntr = traces.shape[0]
    p = np.concatenate([[0.0], ave, [float(ntr)]])
    midp = np.diff(p) / 2.0
    starts = ave - midp[:-1]
    ends = ave + midp[1:]

    het_lab = bc1_arr == bc2_arr
    starthets = np.where(het_lab, np.nan, starts)
    endhets = np.where(het_lab, np.nan, ends)

    idx = np.arange(1, len(bc1) + 1)
    trimmed = (idx <= trim5) | (idx > (len(bc1) - trim3))
    starttrims = np.where(trimmed, np.nan, starts)
    endtrims = np.where(trimmed, np.nan, ends)

    valuesperbase = ntr / len(bc1)
    tracewidth = dlen * valuesperbase

    breaks = [1.0]
    while breaks[-1] + tracewidth <= ntr:
        breaks.append(breaks[-1] + tracewidth)
    numplots = len(breaks)

    # ---- 版面 ----
    dpi = opts.dpi
    fig_w_in = opts.width_px / dpi
    fig_h_in = opts.height_px / dpi
    left_m = 0.4 / fig_w_in
    right_m = 0.0
    band_h = 1.0 / numplots
    # R 里 mar=c(2,2,2,0) 是固定 0.4in 的上下边距；面板一多（例如 showtrim=True 时
    # 不裁剪读序，panel 数会激增）就会把绘图区挤成负高度，R 自己也报
    # "figure margins too large"。这里对边距加一个"不超过面板高度 30%"的上限：
    # numplots=1 时 30% > 12%，取值与 R 完全一致；面板多时退化为可按比例绘制。
    top_m = min(0.4 / fig_h_in, 0.30 * band_h)
    bottom_m = min(0.4 / fig_h_in, 0.30 * band_h)
    line_h_px = opts.pointsize * 1.2 / 72.0 * dpi  # 1 行 = 60px @300dpi/12pt

    font_cjk = font_manager.FontProperties(
        fname=opts.cjk_font or DEFAULT_CJK_FONT
    )
    font_base = font_manager.FontProperties(fname=opts.base_font) if opts.base_font else font_cjk

    fig = plt.figure(figsize=(fig_w_in, fig_h_in), dpi=dpi, facecolor="white")

    usr_x0 = 1.0 - 0.04 * (tracewidth - 1.0)
    usr_x1 = tracewidth + 0.04 * (tracewidth - 1.0)
    usr_y0 = -0.04 * yl
    usr_y1 = 1.04 * yl

    result = None

    for j, i in enumerate(breaks, start=1):
        band_top = 1.0 - (j - 1) * band_h
        band_bottom = 1.0 - j * band_h
        ax_rect = [
            left_m,
            band_bottom + bottom_m,
            1.0 - left_m - right_m,
            band_h - top_m - bottom_m,
        ]
        if ax_rect[2] <= 0 or ax_rect[3] <= 0:
            raise ValueError(
                "绘图区尺寸非正（面板数=%d, 宽=%.4f, 高=%.4f）；"
                "请加大 --width-px/--height-px，或改用更小的 --dlen"
                % (numplots, ax_rect[2], ax_rect[3])
            )
        ax = fig.add_axes(ax_rect, facecolor="white", xlim=(usr_x0, usr_x1), ylim=(usr_y0, usr_y1))
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.patch.set_alpha(0.0)

        # 该面板的碱基
        rng = (ave >= i) & (ave < (i + tracewidth))
        starthet = starthets[rng] - tracewidth * (j - 1)
        starthet = np.clip(starthet, 0, None)
        endhet = endhets[rng] - tracewidth * (j - 1)
        endhet = np.clip(endhet, None, tracewidth)
        lab1 = [x for x, k in zip(bc1, rng) if k]
        lab2 = [x for x, k in zip(bc2, rng) if k]
        lab3 = [x for x, k in zip(ref, rng) if k]
        starttrim = starttrims[rng] - tracewidth * (j - 1)
        endtrim = endtrims[rng] - tracewidth * (j - 1)

        maxvalues = traces.max(axis=1)
        numtop = float(maxvalues[:10].min())
        numtail = float(maxvalues[-10:].min())
        left = int(np.nonzero(maxvalues == numtop)[0][0])
        tail = np.nonzero(maxvalues == numtail)[0]
        right = len(maxvalues) - (int(tail[-1]) + 1)

        lo = int(i) + left
        hi = int(min(i + tracewidth, ntr)) - right
        pr = np.arange(lo, hi + 1)

        # 原脚本此处先算 pos <- aveposition[range] - tracewidth*(j-1) 再用它取长度，
        # 紧接着就把 pos 覆盖成等间距序列，所以实际只需要"本面板碱基个数"。
        n_bases = int(rng.sum())
        peakwidth = (tracewidth - left - right) / n_bases
        pos = np.arange(1, n_bases + 1) * peakwidth - peakwidth / 2.0

        is_rev = reverse

        plot_idx = pr[::-1] if is_rev else pr
        plot_idx = plot_idx[(plot_idx >= 1) & (plot_idx <= ntr)]
        x = np.arange(1, len(plot_idx) + 1)

        for c in range(4):
            color = (_REV_COLORS if is_rev else _FWD_COLORS)[c]
            ax.plot(x, traces[plot_idx - 1, c], color=_RGB[color], lw=0.8, solid_joinstyle="round")

        if opts.showhets:
            for a, b in zip(starthet, endhet):
                _add_span(ax, a, b, usr_y0, usr_y1, facecolor="#D5E3F7")
        if opts.showtrim:
            for a, b in zip(starttrim, endtrim):
                _add_span(ax, a, b, usr_y0, usr_y1, facecolor="red", alpha=0.3)

        # R 向：把碱基字母反向互补后再展示，保证正反向出图一致
        if is_rev:
            lab2 = _revcomp_chars(lab2)
            lab1 = _revcomp_chars(lab1)

        diffindex = [k for k in range(len(lab2)) if k < len(lab3) and lab2[k] != lab3[k]]

        # 1 行在 y 轴上的高度（用户坐标）。plt 高度 = 1.08*yl 对应 ax_rect[3]*height_px 像素
        line_h = line_h_px * yl / (ax_rect[3] * opts.height_px / 1.08)
        # 基线位置按参考图实测（都写成"plt 边缘 ± 若干行"，换画布尺寸时同比缩放）
        y_ref_base = usr_y1 + 0.193 * line_h
        y_btm_base = usr_y0 - 0.808 * line_h
        y_ref_label = usr_y1 + 0.191 * line_h
        y_btm_label = usr_y0 - 0.809 * line_h
        # 字体高度标定：SimSun 的拉丁大写约为 0.875em，SourceHanSansSC 约 0.73em，
        # 乘 1.2 后视觉字号与原图一致；要严格照搬 R 的 pt 值就把它设成 1.0。
        letter_pt = opts.pointsize * opts.cex_base * opts.letter_scale

        box_indices: List[int] = []
        # 原脚本：diff <= 3bp 视为 SNP（逐个加框），否则按 INDEL 只展示次等位字母
        snp_max_diffs = opts.snp_max_diffs
        is_snp = len(diffindex) <= snp_max_diffs
        for k in range(len(lab1)):
            if k >= len(lab3):
                break
            l1, l2, l3 = lab1[k], lab2[k], lab3[k]
            c1 = "black" if l1 in _BASE_LETTERS else "purple"
            c2 = "black" if l2 in _BASE_LETTERS else "purple"
            l1n, l2n = _resolve_iupac_pair(l1, l2, l3)

            ax.text(pos[k], y_ref_base, l3, fontproperties=font_base,
                    fontsize=letter_pt, ha="center", va="baseline",
                    color=c1 if l3 in _BASE_LETTERS else "purple", clip_on=False)

            if is_snp:
                if l1n == l2n and l3 == l2n:  # 野生型
                    ax.text(pos[k], y_btm_base, l2n, fontproperties=font_base,
                            fontsize=letter_pt, ha="center", va="baseline",
                            color=c1, clip_on=False)
                elif l1n == l2n and l3 != l2n:  # 纯合突变
                    ax.text(pos[k], y_btm_base, l2n, fontproperties=font_base,
                            fontsize=letter_pt, ha="center", va="baseline",
                            color=c1, clip_on=False)
                    _draw_box(ax, pos[k], peakwidth, yl)
                    box_indices.append(k)
                else:  # 杂合突变
                    ax.text(pos[k], y_btm_base, l1n + l2n, fontproperties=font_base,
                            fontsize=letter_pt, ha="center", va="baseline",
                            color=c1, clip_on=False)
                    _draw_box(ax, pos[k], peakwidth, yl)
                    box_indices.append(k)
            else:
                ax.text(pos[k], y_btm_base, l2n, fontproperties=font_base,
                        fontsize=letter_pt, ha="center", va="baseline",
                        color=c2, clip_on=False)

        # ---- 箭头 ----
        # 注意：箭头只画在正中位点（pos[dlen2]）；"SNP 未检出"时 R 干脆不画，
        # 此时 arrow_pos 必须是 None，否则调用方（CLI / 测试）会误以为图上有箭头。
        arrow_pos = None
        if len(pos) > 10:
            cand = float(pos[10])
            draw_arrow = True
            if len(diffindex) > snp_max_diffs:
                if abs(diffindex[-1] - 10) < abs(diffindex[0] - 10):
                    cand += peakwidth / 2.0
                else:
                    cand -= peakwidth / 2.0
            elif len(diffindex) <= snp_max_diffs and len(lab1) > 10 and len(lab3) > 10 \
                    and lab1[10] == lab2[10] and lab3[10] == lab2[10]:
                draw_arrow = False  # SNP 未检出
            if draw_arrow:
                arrow_pos = cand
                # R 的 arrows(code=1) 箭头在 (x0,y0) 端，即朝下指向峰；shrink 置 0
                # 才能让箭头端点严格落在 0.90yl / 0.65yl 上（默认 2pt 会缩掉约 8px）
                ax.annotate(
                    "",
                    xy=(arrow_pos, 0.65 * yl),
                    xytext=(arrow_pos, 0.90 * yl),
                    arrowprops=dict(arrowstyle="-|>", color="red", lw=1.0,
                                    mutation_scale=10, shrinkA=0, shrinkB=0),
                    annotation_clip=False,
                )

        # ---- 左右两侧中文说明 ----
        x_label = usr_x0 + (usr_x1 - usr_x0) * 0.0202
        ax.text(x_label, y_ref_label, "参考序列", fontproperties=font_cjk,
                fontsize=opts.pointsize * opts.cex_mtext, ha="right", va="baseline",
                clip_on=False)
        ax.text(x_label, y_btm_label, "目标序列", fontproperties=font_cjk,
                fontsize=opts.pointsize * opts.cex_mtext, ha="right", va="baseline",
                clip_on=False)

        result = RenderResult(
            n_plots=numplots,
            tracewidth=tracewidth,
            ylims=ylims,
            aveposition=ave.copy(),
            pos=pos.copy(),
            lab1=lab1,
            lab2=lab2,
            lab3=lab3,
            diffindex=diffindex,
            box_indices=box_indices,
            arrow_pos=arrow_pos,
            plotrange=(int(pr[0]), int(pr[-1])),
            breaks=list(breaks),
        )

    if outfile:
        os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
        fig.savefig(outfile, dpi=dpi, facecolor="white")
    plt.close(fig)
    assert result is not None
    return result


def _add_span(ax, a: float, b: float, y0: float, y1: float, **kwargs) -> None:
    """画一条横跨面板的底色/阴影带。

    R 的 `rect()` 允许 xleft > xright（会自动交换），但 matplotlib 要求宽高非负，
    所以这里统一取 min/abs；并用 `clip_on=False` 对齐 R 的 `xpd=NA`（trim 区间通常
    就在展示窗口之外的边缘，裁剪掉就看不见了）。
    """
    if np.isnan(a) or np.isnan(b):
        return
    ax.add_patch(
        Rectangle(
            (min(a, b), y0),
            abs(b - a),
            y1 - y0,
            edgecolor="none",
            zorder=0,
            clip_on=False,
            **kwargs,
        )
    )


def _draw_box(ax, x_center: float, peakwidth: float, yl: float) -> None:
    """原样照搬：rect 的 xleft 比 xright 大（R 不在意顺序）。"""
    ax.add_patch(
        Rectangle(
            (x_center - peakwidth / 2.0, -0.14 * yl),
            peakwidth,
            0.11 * yl,
            fill=False,
            edgecolor="black",
            linewidth=1.0,
            clip_on=False,
        )
    )


_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "M": "K", "K": "M", "R": "Y", "Y": "R", "W": "W", "S": "S",
    "V": "B", "B": "V", "H": "D", "D": "H", "N": "N",
}


def _revcomp_chars(chars: Sequence[str]) -> List[str]:
    out = []
    for ch in reversed(list(chars)):
        out.append(_COMPLEMENT.get(ch, ch))
    return out


def _resolve_iupac_pair(l1: str, l2: str, l3: str):
    """原脚本对 R/Y/M/K/S/W 这类"双等位写成同一模糊码"的情况做拆分。"""
    if l1 == l2 and l1 in ("R", "Y", "M", "K", "S", "W"):
        pair = {
            "R": ("A", "G"), "Y": ("C", "T"), "M": ("A", "C"),
            "K": ("T", "G"), "S": ("C", "G"), "W": ("A", "T"),
        }[l1]
        if l3 == pair[0]:
            return pair[0], pair[1]
        if l3 == pair[1]:
            return pair[1], pair[0]
    return l1, l2
