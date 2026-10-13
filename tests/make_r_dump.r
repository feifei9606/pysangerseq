#!/usr/bin/env Rscript

# 生成"原始 R 流程"的中间量 dump，供 compare_r_dump.py 做逐项比对。
#
# 原理：把原脚本 script.r 从 chromatogramref() 定义处截断（跳过画图），
# 再追加一段 dump 代码。因为不画图，所以不需要 Cairo；比对完全在数值层面。
#
# 用法:
#   Rscript make_r_dump.r <script.r 路径> <dump 输出目录> <postion> <ab1> <hg19|hg38> \
#            <dlen> <samtools> <makeSCF> <python2> <swalign> <hg19.fa> <hg38.fa>
#
# 两个必须的修补（原脚本在当前环境跑不通，原因见 README）：
#   1) 去掉 library(Cairo)：本机缺 libEGL，Cairo 包加载失败（此处本来也不画图）
#   2) `sed -n 2p` -> `tail -1`：swalign 的 -summary 数据行在 stdout 最后一行；
#      原写法取到的是 dump 的第 2 行 "Ref  : ..."，只有 5 个字段，会报 NA 错误

argv <- commandArgs(TRUE)
script_r <- argv[1]
outdir   <- argv[2]
rest     <- argv[-(1:2)]

src <- paste(readLines(script_r, warn = FALSE), collapse = "\n")
cut <- regexpr("chromatogramref <- function\\(", src)
if (cut < 0) stop("在 ", script_r, " 中找不到 chromatogramref 定义")
head <- substr(src, 1, cut - 1)
head <- gsub("library\\(Cairo, quietly=TRUE\\)", "", head)
head <- gsub("\\|sed -n 2p", "|tail -1", head)

dump <- '
######################## FULL DUMP ########################
od <- OUTDIR
dir.create(od, showWarnings = FALSE)

wt  <- function(df, f) write.table(df, file.path(od, f), sep="\t", quote=FALSE,
                                   col.names=TRUE, row.names=FALSE)
wtc <- function(x, f) write.table(data.frame(v=as.character(x)), file.path(od, f),
                                  sep="\t", quote=FALSE, row.names=FALSE)

wt(as.data.frame(hetsangerseq@traceMatrix),   "ab1_trace.tsv")
wt(as.data.frame(hetsangerseq@peakPosMatrix), "ab1_peakpos.tsv")
wtc(as.character(hetsangerseq@primarySeq),    "ab1_primary.txt")
wtc(as.character(hetsangerseq@secondarySeq),  "ab1_secondary.txt")

wt(as.data.frame(hetcalls@peakPosMatrix), "mbc_peakpos.tsv")
wt(as.data.frame(hetcalls@peakAmpMatrix), "mbc_peakamp.tsv")
wtc(as.character(hetcalls@primarySeq),    "mbc_primary.txt")
wtc(as.character(hetcalls@secondarySeq),  "mbc_secondary.txt")

wtc(as.character(ref), "ref500.txt")

bc <- sangerseqR:::basecalldf(hetcalls)
wt(bc, "basecalldf.tsv")

n <- nrow(bc)
seedseq <- paste(bc$consensus[(29+1):(n-30)], collapse="")
wtc(seedseq, "seedseq.txt")

pa   <- pairwiseAlignment(seedseq, ref, type="local", gapOpening=-200, gapExtension=-10)
paRC <- pairwiseAlignment(seedseq, reverseComplement(DNAString(ref)), type="local",
                          gapOpening=-200, gapExtension=-10)
writeLines(c(
  paste0("score_fwd\t", score(pa)),
  paste0("score_rc\t",  score(paRC)),
  paste0("pat_start\t", start(pattern(pa))),      paste0("pat_end\t", end(pattern(pa))),
  paste0("sub_start\t", start(subject(pa))),      paste0("sub_end\t", end(subject(pa))),
  paste0("patRC_start\t", start(pattern(paRC))),  paste0("patRC_end\t", end(pattern(paRC))),
  paste0("subRC_start\t", start(subject(paRC))),  paste0("subRC_end\t", end(subject(paRC)))
), file.path(od, "pa.txt"))

writeLines(c(
  paste0("infile\t",            infile),
  paste0("name\t",              name),
  paste0("chr\t",               chr),
  paste0("input_pos\t",         ainfo[2]),
  paste0("input_ref\t",         refsequence),
  paste0("input_alt\t",         altsequence),
  paste0("altposs_final\t",     altposs),
  paste0("dlen\t",              dlen),
  paste0("dlen2\t",             dlen2),
  paste0("trim5len\t",          trim5len),
  paste0("trim3len\t",          trim3len),
  paste0("len\t",               len),
  paste0("refseq1\t",           as.character(refseq1)),
  paste0("refseq1_len\t",       nchar(as.character(refseq1))),
  paste0("posref2\t",           paste(posref2, collapse="|")),
  paste0("posref\t",            paste(posref,  collapse="|")),
  paste0("primarySeq\t",        as.character(hetseqalleles@primarySeq)),
  paste0("secondarySeq\t",      as.character(hetseqalleles@secondarySeq)),
  paste0("Rdir_match\t",        grepl("R_[a-zA-Z0-9]+.ab1", name, perl=TRUE) ||
                                grepl("-R[0-9S]*_[a-zA-Z0-9-]+_[a-zA-Z0-9]+_[a-zA-Z0-9.]+.ab1", name, perl=TRUE))
), file.path(od, "final.txt"))
cat("DUMP OK ->", od, "\n")
'
dump <- gsub("OUTDIR", sprintf('"%s"', outdir), dump, fixed = TRUE)

script <- paste0(head, dump, "\n")

tmp <- tempfile(fileext = ".r")
writeLines(script, tmp)

# 后缀参数（与 script.r 的位置参数完全一致）
args <- c(rest[1], rest[2], rest[3], rest[4], rest[5], rest[6],
          rest[7], rest[8], rest[9], rest[10])
system2(file.path(R.home("bin"), "Rscript"), c(tmp, args))
