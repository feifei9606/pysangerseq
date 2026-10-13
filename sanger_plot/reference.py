"""FASTA 随机访问 —— 替代 `samtools faidx`。

优先读同目录的 `.fai`；没有就现场扫描建立内存索引（不改动原文件，
因为参考基因组往往在别的账号目录下、无写权限）。

坐标一律使用 **1-based 闭区间**，与 samtools faidx 命令行一致，
方便和原 R 脚本的输出逐项对照。
"""

from __future__ import annotations

import os
from typing import Dict, List, Tuple


class Fasta:
    """按 `.fai` 索引随机取子序列。"""

    def __init__(self, path: str):
        if not os.path.exists(path):
            raise FileNotFoundError("参考基因组不存在: %s" % path)
        self.path = path
        self.fai_path = path + ".fai"
        self._fh = open(path, "rb")
        # name -> (length, offset_of_first_base, bases_per_line, bytes_per_line)
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        if os.path.exists(self.fai_path):
            self._load_fai()
        else:
            self._build_index()

    # ---------------------------------------------------------------- 索引
    def _load_fai(self) -> None:
        with open(self.fai_path, "r") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                name, length, offset, lb, lw = parts[:5]
                self.index[name] = (int(length), int(offset), int(lb), int(lw))

    def _build_index(self) -> None:
        """没有 .fai 时扫描一遍文件建索引（只在内存中，不落盘）。"""
        name = None
        length = 0
        first_base = 0
        bases_per_line = 0
        bytes_per_line = 0
        with open(self.path, "rb") as fh:
            pos = 0
            while True:
                line = fh.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    if name is not None:
                        self.index[name] = (length, first_base, bases_per_line, bytes_per_line)
                    name = line[1:].split()[0].decode("latin-1")
                    length = 0
                    first_base = pos + len(line)
                    bases_per_line = 0
                    bytes_per_line = 0
                else:
                    stripped = line.rstrip(b"\r\n")
                    if bases_per_line == 0:
                        bases_per_line = len(stripped)
                        bytes_per_line = len(line)
                    length += len(stripped)
                pos += len(line)
        if name is not None:
            self.index[name] = (length, first_base, bases_per_line, bytes_per_line)

    # ---------------------------------------------------------------- 取序列
    def chrom_length(self, chrom: str) -> int:
        if chrom not in self.index:
            raise KeyError("参考基因组中没有序列: %s" % chrom)
        return self.index[chrom][0]

    def has_chrom(self, chrom: str) -> bool:
        return chrom in self.index

    def chroms(self) -> List[str]:
        return list(self.index)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """取 `chrom:start-end`（1-based 闭区间），越界自动夹紧，返回大写。

        与 samtools faidx 一致：超出部分被裁掉。
        """
        if chrom not in self.index:
            raise KeyError("参考基因组中没有序列: %s" % chrom)
        length, offset, lb, lw = self.index[chrom]
        if start < 1:
            start = 1
        if end > length:
            end = length
        if start > end:
            return ""
        out: List[bytes] = []
        # 只读首行所在的那个"行", 然后按固定步长跳读
        p = start - 1
        remaining = end - start + 1
        while remaining > 0:
            line_no, col = divmod(p, lb)
            take = min(lb - col, remaining)
            self._fh.seek(offset + line_no * lw + col)
            out.append(self._fh.read(take).rstrip(b"\r\n"))
            p += take
            remaining -= take
        return b"".join(out).decode("latin-1").upper()

    def close(self) -> None:
        self._fh.close()

    def __enter__(self) -> "Fasta":
        return self

    def __exit__(self, *exc) -> None:
        self.close()
