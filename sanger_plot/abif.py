"""ABIF(.ab1) 文件解析。

等价于 sangerseqR::read.abif() + sangerseq(abif) 方法，逐字节对齐其
偏移量/字节序/取数方式，因此结果与原 R 流程完全一致。

要点（与原 R 实现保持一致，勿"顺手优化"）
------------------------------------------
* 目录项偏移：dataoffset(header) + i*elementsize，字段依次为
  name(4) tagnumber(4) elementtype(2) elementsize(2)
  numelements(4) datasize(4) dataoffset(4)。
* 数据起始：datasize > 4 时取该目录项的 dataoffset；否则取目录项内
  偏移 20 处的 4 字节内联值。R 里数组下标是 1-based 且切片是闭区间，
  所以数据区实际会多读 1 个字节——这里刻意保留（对字符串类型有影响，
  调用方会再截断）。
* 整数一律 **大端**。sangerseqR 用 SInt16（有符号）读 DATA/PLOC/P1AM，
  未归一化的原始值 >32767 会变成负数，这里同样保留该行为。
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .sangerseq import Sangerseq

# Biostrings::DNA_ALPHABET
DNA_ALPHABET = set("ACGTMRWSYKVHDBN-+.")


class AbifError(ValueError):
    """ABIF 文件格式错误。"""


# --------------------------------------------------------------------------
# 标量解码
# --------------------------------------------------------------------------
def _sint16(buf: bytes, off: int) -> int:
    return struct.unpack_from(">h", buf, off)[0]


def _sint32(buf: bytes, off: int) -> int:
    return struct.unpack_from(">i", buf, off)[0]


def _decode(etype: int, raw: bytes, numelements: int) -> Any:
    """按 elementtype 解码，语义与 sangerseqR 中的 *Int*/f32/f64/RTC 一致。"""
    if etype == 1:  # UInt8
        return list(raw[:numelements])
    if etype == 2:  # char (rawToChar)
        return raw.decode("latin-1")
    if etype == 3:  # UInt16
        return list(struct.unpack_from(">%dH" % numelements, raw, 0))
    if etype == 4:  # SInt16
        return list(struct.unpack_from(">%dh" % numelements, raw, 0))
    if etype == 5:  # SInt32
        return list(struct.unpack_from(">%di" % numelements, raw, 0))
    if etype == 7:  # float32
        return list(struct.unpack_from(">%df" % numelements, raw, 0))
    if etype == 8:  # float64
        return list(struct.unpack_from(">%dd" % numelements, raw, 0))
    if etype == 10:  # date
        return {
            "year": _sint16(raw, 0),
            "month": raw[2],
            "day": raw[3],
        }
    if etype == 11:  # time
        return {"hour": raw[0], "minute": raw[1], "second": raw[2], "hsecond": raw[3]}
    if etype == 18:  # pString
        n = struct.unpack_from(">b", raw, 0)[0]
        return raw[1 : 1 + n].decode("latin-1")
    if etype == 19:  # cString
        return raw[: max(len(raw) - 1, 0)].decode("latin-1")
    # 1024+ 为私有类型，原样返回字节
    return raw


@dataclass
class AbifFile:
    """解析后的 ABIF 内容。"""

    version: int
    elementsize: int
    numelements: int
    directory: List[Tuple[str, int, int, int, int]]  # (name, tagnumber, etype, size, offset)
    data: Dict[str, Any]  # key = "<name>.<tagnumber>"

    def get(self, key: str) -> Any:
        if key not in self.data:
            raise AbifError("ABIF tag not found: %s" % key)
        return self.data[key]

    def has(self, key: str) -> bool:
        return key in self.data


def parse_abif(path: str) -> AbifFile:
    """解析 ABIF 文件的头部、目录与各 tag 数据。"""
    with open(path, "rb") as fh:
        buf = fh.read()

    if len(buf) < 30 or buf[0:4] != b"ABIF":
        raise AbifError("file not in ABIF format: %s" % path)

    version = _sint16(buf, 4)
    elementsize = _sint16(buf, 16)
    numelements = _sint32(buf, 18)
    dataoffset = _sint32(buf, 26)

    directory: List[Tuple[str, int, int, int, int]] = []
    entries: List[Tuple[str, int, int, int, int, int, int]] = []
    for i in range(numelements):
        d = dataoffset + i * elementsize
        if d + 24 > len(buf):
            break
        name = buf[d : d + 4].decode("latin-1")
        tagnumber = _sint32(buf, d + 4)
        etype = _sint16(buf, d + 8)
        esize = _sint16(buf, d + 10)
        nel = _sint32(buf, d + 12)
        dsize = _sint32(buf, d + 16)
        doff = _sint32(buf, d + 20)
        directory.append((name, tagnumber, etype, esize, nel))
        entries.append((name, tagnumber, etype, esize, nel, dsize, doff))

    data: Dict[str, Any] = {}
    for i, (name, tagnumber, etype, esize, nel, dsize, doff) in enumerate(entries):
        # datasize > 4：数据在文件别处（用目录项记录的 dataoffset）
        # 否则：4 字节内联在目录项偏移 20 处
        start = doff if dsize > 4 else (dataoffset + i * elementsize + 20)
        # 与 R 的闭区间切片一致：多读 1 字节
        raw = buf[start : start + nel * esize + 1]
        data["%s.%d" % (name, tagnumber)] = _decode(etype, raw, nel)

    return AbifFile(
        version=version,
        elementsize=elementsize,
        numelements=numelements,
        directory=directory,
        data=data,
    )


# --------------------------------------------------------------------------
# sangerseq(abif)
# --------------------------------------------------------------------------
def _clean_basecalls(text: str, limit: int) -> str:
    """对应 R：剔除非 DNA_ALPHABET 字符后截取前 limit 个。"""
    kept = "".join(c for c in text if c in DNA_ALPHABET)
    return kept[:limit]


def sangerseq_from_abif(obj: AbifFile) -> Sangerseq:
    """由 ABIF 对象构造 trace/峰位/碱基序列，等价 sangerseq(abif)。

    返回公共的 `sangerseq.Sangerseq`（与 make_base_calls / set_allele_phase 同一个类型），
    调用方不需要再手工包一层。
    """
    trace_cols = [obj.get(k) for k in ("DATA.9", "DATA.10", "DATA.11", "DATA.12")]
    n = min(len(c) for c in trace_cols)
    tracematrix = np.column_stack([np.asarray(c[:n], dtype=float) for c in trace_cols])

    fwo = obj.get("FWO_.1")
    order = []
    for base in "ACGT":
        idx = fwo.find(base)  # regexpr 取首次出现位置，1-based
        if idx < 0:
            raise AbifError("FWO_ tag %r 中找不到碱基 %s" % (fwo, base))
        order.append(idx)
    trace = tracematrix[:, order]

    ploc = np.asarray(obj.get("PLOC.2"), dtype=float)
    nloc = len(ploc)
    primary = _clean_basecalls(obj.get("PBAS.2"), nloc)
    if len(primary) != nloc:
        # R 只在字符被剔除时给 warning，不中断
        pass

    if obj.has("P2BA.1"):
        secondary = _clean_basecalls(obj.get("P2BA.1"), nloc)
    else:
        secondary = ""

    if obj.has("P1AM.1"):
        p1 = np.asarray(obj.get("P1AM.1"), dtype=float)
    else:
        p1 = np.full(nloc, np.nan)
    if obj.has("P2AM.1"):
        p2 = np.asarray(obj.get("P2AM.1"), dtype=float)
    else:
        p2 = np.full(nloc, np.nan)

    nan_col = np.full((nloc, 1), np.nan)
    peak_pos = np.hstack([(ploc + 1).reshape(-1, 1), (ploc + 1).reshape(-1, 1), nan_col, nan_col])
    peak_amp = np.hstack([p1.reshape(-1, 1), p2.reshape(-1, 1), nan_col, nan_col])

    return Sangerseq(
        trace=trace,
        peak_pos=peak_pos,
        peak_amp=peak_amp,
        primary_seq=primary,
        secondary_seq=secondary,
    )


def read_abif(path: str) -> Sangerseq:
    """一步到位：读取 .ab1 并返回 `Sangerseq`（trace / 峰位 / 碱基序列）。"""
    return sangerseq_from_abif(parse_abif(path))
