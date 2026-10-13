#!/usr/bin/env bash
# 端到端回归：对每个 (ab1, 变异) 组合，先跑原始 R 流程 dump 中间量，再和 Python 实现逐项比对。
#
# 覆盖 SNV / 插入 / 缺失 / 重复序列（INDEL 断点校正路径）。
#
# 用法:  bash tests/run_all.sh
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
WORK="${WORK:-/tmp/sanger_regress}"
RSCRIPT="${RSCRIPT:-/bi/software/R-4.0.2/bin/Rscript}"
PYTHON="${PYTHON:-/bi/software/Anaconda3/bin/python}"
ORIG_R="${ORIG_R:-/sg2/7.hudahui/tmp/20260819/script.r}"
AB1_DIR="${AB1_DIR:-/sg2/7.hudahui/tmp/20260819}"
SANGERSEQ_DATA="${SANGERSEQ_DATA:-$AB1_DIR}"
FASTA_HG38="${FASTA_HG38:-/sg2/19.yuli/02WES/Sangerseq/Homo_sapiens_assembly38.fasta}"
FASTA_HG19="${FASTA_HG19:-/sg2/19.yuli/02WES/Sangerseq/human_g1k_v37.fasta}"
SAMTOOLS="${SAMTOOLS:-/bi/5.liangyuanyuan/sw/anaconda3/bin/samtools}"
MAKESCF="${MAKESCF:-/bi/5.liangyuanyuan/sw/io_lib-1.14.8/progs/makeSCF}"
PY2="${PY2:-/bi/software/Python-2.7.18/bin/python2}"
SWALIGN="${SWALIGN:-/bi/5.liangyuanyuan/sw/swalign/bin/swalign}"

export R_LIBS="${R_LIBS:-/bi/software/R-4.0.2/lib64/R/library}"
export LD_LIBRARY_PATH="/bi/7.hudahui/soft/conda/lib:/bi/software/R-4.0.2/lib64/R/lib:/usr/lib64:/bi/software/staden-2.0.0b11-2016-linux-x86_64/lib:/bi/5.liangyuanyuan/sw/io_lib-1.14.8/lib:/bi/software/Python-2.7.18/lib:${LD_LIBRARY_PATH:-}"

rm -rf "$WORK"; mkdir -p "$WORK"
cp "$AB1_DIR"/*.ab1 "$WORK"/ 2>/dev/null || true

F_FILE=$(basename "$(ls "$WORK"/*-F_*.ab1 2>/dev/null | head -1)")
R_FILE=$(basename "$(ls "$WORK"/*-R_*.ab1 2>/dev/null | head -1)")
if [ -z "$F_FILE" ] || [ -z "$R_FILE" ]; then
    echo "在 $AB1_DIR 找不到预期的 F/R ab1 文件" >&2
    exit 2
fi

# 变异列表：SNV / 插入 / 缺失 / poly-A 重复 / 只给位点
CASES=(
    "12-14796923-G-A"
    "12-14796923-G-GTT"
    "12-14796923-GT-G"
    "12-14796923-G-GA"
    "12-14796923"      # 只给位点、无 ref/alt —— 原脚本用字面量 "ref"/"alt" 当哨兵，
                       # 若不显式排除，反链读序会被平移 len(alt)=3（曾真实发生的 bug）
)

pass=0; fail=0; skip=0
for VAR in "${CASES[@]}"; do
    for AB1 in "$F_FILE" "$R_FILE"; do
        TAG="$(echo "$VAR" | tr -c 'A-Za-z0-9' '_')__${AB1: -34}"
        DUMP="$WORK/dump_${TAG}"
        # 必须在 ab1 所在目录跑（原脚本要求），否则 makeSCF 的临时文件会散落
        ( cd "$WORK" && "$RSCRIPT" "$HERE/make_r_dump.r" "$ORIG_R" "$DUMP" "$VAR" "$AB1" \
              hg38 21 "$SAMTOOLS" "$MAKESCF" "$PY2" "$SWALIGN" "$FASTA_HG19" "$FASTA_HG38" \
              >/dev/null 2>&1 )
        if [ ! -f "$DUMP/final.txt" ]; then
            echo "SKIP  $VAR  $AB1   (原始 R 流程未产出结果)"
            skip=$((skip+1)); continue
        fi
        IFS='-' read -r CHR POS REF ALT <<< "$VAR"
        if "$PYTHON" "$HERE/compare_r_dump.py" "$DUMP" "$WORK/$AB1" "$CHR" "$POS" "$REF" "$ALT" "$FASTA_HG38" > "$DUMP/result.log" 2>&1; then
            echo "PASS  $VAR  $AB1"
            pass=$((pass+1))
        else
            echo "FAIL  $VAR  $AB1   -> $DUMP/result.log"
            sed -n '1,200p' "$DUMP/result.log" | grep -E "FAIL|got|want" | head -20
            fail=$((fail+1))
        fi
    done
done

echo
echo "== 汇总: PASS=$pass FAIL=$fail SKIP=$skip =="
[ "$fail" -eq 0 ]
