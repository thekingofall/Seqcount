#!/usr/bin/env python3


import os
import glob
import gzip
import argparse
import sys
import subprocess
from tqdm import tqdm
import shutil

def compute_complement(seq: str) -> str:
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(comp_dict.get(base, base) for base in seq)

def is_match(seq: str, pattern: str, tol: int) -> bool:
    """
    如果 tol == 0，则直接进行子串匹配；
    否则采用滑动窗口对每个窗口计算汉明距离，
    只要有一个窗口的距离 <= tol，即认为匹配成功。
    """
    m = len(pattern)
    if len(seq) < m:
        return False
    if tol == 0:
        return pattern in seq
    else:
        for i in range(len(seq) - m + 1):
            window = seq[i:i+m]
            dist = sum(1 for a, b in zip(window, pattern) if a != b)
            if dist <= tol:
                return True
        return False

def process_single_file(file: str, pattern_fw: str, pattern_comp: str,
                        pattern_rev: str, pattern_revcomp: str,
                        max_reads: int, tolerance: int) -> str:
    """
    处理单个 FASTQ 文件，统计匹配计数，返回一行制表分隔的统计结果：
    Filename, Sampled_Reads, Count_Forward, Pct_Forward, Count_Complement, Pct_Complement,
    Count_Reverse, Pct_Reverse, Count_RevComp, Pct_RevComp, Total_Matches, Pct_Total
    """
    count_fw = 0
    count_comp = 0
    count_rev = 0
    count_revcomp = 0
    total_reads = 0

    try:
        with gzip.open(file, "rt") as fh:
            # 判断是否显示进度条
            # 新增 --show-pbar 或 -P 参数时，强制显示进度条；否则看命令行是否包含 --noheader、-n 或 --parallel（默认并行时不显示）
            force_show = any(arg in sys.argv for arg in ("--show-pbar", "-P"))
            if force_show:
                show_pbar = True
            else:
                show_pbar = True
                for arg in sys.argv:
                    if arg in ("--noheader", "-n", "--parallel", "-p"):
                        show_pbar = False
                        break

            pbar = tqdm(total=max_reads, desc=f"Processing {file}") if show_pbar else None

            while total_reads < max_reads:
                header = fh.readline()
                if not header:
                    break
                seq_line = fh.readline().strip()
                fh.readline()  # 略过 '+'
                fh.readline()  # 略过质量值行
                total_reads += 1
                if pbar:
                    pbar.update(1)
                seq_lower = seq_line.lower()
                if is_match(seq_lower, pattern_fw, tolerance):
                    count_fw += 1
                if is_match(seq_lower, pattern_comp, tolerance):
                    count_comp += 1
                if is_match(seq_lower, pattern_rev, tolerance):
                    count_rev += 1
                if is_match(seq_lower, pattern_revcomp, tolerance):
                    count_revcomp += 1

            if pbar:
                pbar.close()
    except Exception as e:
        print(f"处理文件 {file} 时出错：{e}", file=sys.stderr)
        return ""
    sample_total = total_reads if total_reads > 0 else max_reads
    pct_fw = (count_fw / sample_total * 100) if sample_total else 0
    pct_comp = (count_comp / sample_total * 100) if sample_total else 0
    pct_rev = (count_rev / sample_total * 100) if sample_total else 0
    pct_revcomp = (count_revcomp / sample_total * 100) if sample_total else 0
    total_matches = count_fw + count_comp + count_rev + count_revcomp
    pct_total = (total_matches / sample_total * 100) if sample_total else 0

    result_line = (f"{file}\t{sample_total}\t{count_fw}\t{pct_fw:.4f}\t"
                   f"{count_comp}\t{pct_comp:.4f}\t{count_rev}\t{pct_rev:.4f}\t"
                   f"{count_revcomp}\t{pct_revcomp:.4f}\t{total_matches}\t{pct_total:.4f}")
    return result_line

def run_parallel(args, pattern_fw: str, pattern_comp: str,
                 pattern_rev: str, pattern_revcomp: str):
    """
    并行模式：生成命令文件，每行指令调用本脚本 '--fqfile' 模式处理单个文件，
    然后利用 Parafly 执行并行任务，最后合并各任务输出为最终结果文件。
    """
    files = sorted(glob.glob("*.fq.gz"))
    if not files:
        print("当前目录下未找到 .fq.gz 文件。", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmpdir, exist_ok=True)
    cmd_file = os.path.join(args.tmpdir, "commands.txt")
    with open(cmd_file, "w") as fcmd:
        for file in files:
            base = os.path.basename(file)
            out_file = os.path.join(args.tmpdir, base + ".result.txt")
            script_path = os.path.abspath(sys.argv[0])
            # 如果想在并行任务中也显示进度条，则同时传入 --show-pbar (-P)
            cmd = (f"python3 {script_path} --fqfile {file} -f {args.forward} "
                   f"-m {args.max_reads} -t {args.tolerance} -o {out_file} -n -P")
            fcmd.write(cmd + "\n")

    parafly_cmd = ["ParaFly", "-c", cmd_file, "-CPU", str(args.threads)]
    print("开始使用 Parafly 并行处理...")
    result = subprocess.run(parafly_cmd)
    if result.returncode != 0:
        print("Parafly 执行失败！", file=sys.stderr)
        sys.exit(1)

    result_lines = []
    for file in files:
        base = os.path.basename(file)
        out_file = os.path.join(args.tmpdir, base + ".result.txt")
        if os.path.exists(out_file):
            with open(out_file, "r") as fin:
                line = fin.read().strip()
                if line:
                    result_lines.append(line)
        else:
            print(f"警告：未找到输出文件 {out_file}", file=sys.stderr)

    result_lines.sort()
    try:
        with open(args.outfile, "w") as fout:
            header_line = ("Filename\tSampled_Reads\tCount_Forward\tPct_Forward\t"
                           "Count_Complement\tPct_Complement\tCount_Reverse\tPct_Reverse\t"
                           "Count_RevComp\tPct_RevComp\tTotal_Matches\tPct_Total\n")
            fout.write(header_line)
            for line in result_lines:
                fout.write(line + "\n")
    except Exception as e:
        print(f"写入输出文件 {args.outfile} 失败：{e}", file=sys.stderr)
        sys.exit(1)
    print(f"\n结果已成功输出到 {args.outfile}")
    print(f"临时文件保存在 {args.tmpdir}，如不需要请手动删除。")
    # 如希望自动删除临时文件夹，可取消下面这行注释
    # shutil.rmtree(args.tmpdir)

def main():
    parser = argparse.ArgumentParser(
        description="对 FASTQ 文件进行序列匹配统计，支持并行处理（默认启用 Parafly 并行）"
    )
    parser.add_argument("--forward", "-f",
                        default="CTATAGCGAAACATCGGCCACCTATA",
                        help="正向 DNA 序列，默认为 CTATAGCGAAACATCGGCCACCTATA")
    parser.add_argument("--max_reads", "-m", type=int, default=10000,
                        help="处理的最大 FASTQ 记录数，默认 10000 条")
    parser.add_argument("--tolerance", "-t", type=int, default=0,
                        help="匹配时允许的最大错误碱基数，默认 0（精确匹配）")
    parser.add_argument("--outfile", "-o", default="result.txt",
                        help="输出结果的 txt 文件，默认 result.txt")
    parser.add_argument("--fqfile", "-q",
                        help="指定单个 FASTQ 文件进行处理", default=None)
    parser.add_argument("--noheader", "-n", action="store_true",
                        help="在输出中不包含表头（适用于并行子任务）")
    # 新增一个参数，用于强制显示进度条
    parser.add_argument("--show-pbar", "-P", action="store_true",
                        help="强制显示进度条（即使在并行模式下也显示，可能导致输出混乱）")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-p", "--parallel", dest="parallel", action="store_true",
                       help="使用 Parafly 并行处理 (默认开启)")
    group.add_argument("--no-parallel", dest="parallel", action="store_false",
                       help="不使用 Parafly 并行处理")
    parser.set_defaults(parallel=True)
    parser.add_argument("--threads", "-T", type=int, default=4,
                        help="并行处理的线程数，默认 4")
    parser.add_argument("--tmpdir", "-d", default="tmp_parafly",
                        help="存放并行任务临时文件的目录，默认 tmp_parafly")
    args = parser.parse_args()

    forward = args.forward.strip()
    complement = compute_complement(forward)
    reverse = forward[::-1]
    reverse_complement = compute_complement(reverse)
    pattern_fw = forward.lower()
    pattern_comp = complement.lower()
    pattern_rev = reverse.lower()
    pattern_revcomp = reverse_complement.lower()

    if args.fqfile:
        result_line = process_single_file(args.fqfile, pattern_fw, pattern_comp,
                                            pattern_rev, pattern_revcomp,
                                            args.max_reads, args.tolerance)
        if args.outfile:
            try:
                with open(args.outfile, "w") as f:
                    if not args.noheader:
                        header_line = ("Filename\tSampled_Reads\tCount_Forward\tPct_Forward\t"
                                       "Count_Complement\tPct_Complement\tCount_Reverse\tPct_Reverse\t"
                                       "Count_RevComp\tPct_RevComp\tTotal_Matches\tPct_Total\n")
                        f.write(header_line)
                    f.write(result_line + "\n")
            except Exception as e:
                print(f"写入输出文件 {args.outfile} 失败：{e}", file=sys.stderr)
                sys.exit(1)
        else:
            if not args.noheader:
                header_line = ("Filename\tSampled_Reads\tCount_Forward\tPct_Forward\t"
                               "Count_Complement\tPct_Complement\tCount_Reverse\tPct_Reverse\t"
                               "Count_RevComp\tPct_RevComp\tTotal_Matches\tPct_Total")
                print(header_line)
            print(result_line)
        sys.exit(0)

    if args.parallel:
        run_parallel(args, pattern_fw, pattern_comp, pattern_rev, pattern_revcomp)
    else:
        files = sorted(glob.glob("*.fq.gz"))
        if not files:
            print("当前目录下未找到 .fq.gz 文件。", file=sys.stderr)
            sys.exit(1)
        try:
            with open(args.outfile, "w") as outf:
                header_line = ("Filename\tSampled_Reads\tCount_Forward\tPct_Forward\t"
                               "Count_Complement\tPct_Complement\tCount_Reverse\tPct_Reverse\t"
                               "Count_RevComp\tPct_RevComp\tTotal_Matches\tPct_Total\n")
                outf.write(header_line)
                print(header_line, end="")
                for f in files:
                    result_line = process_single_file(f, pattern_fw, pattern_comp,
                                                      pattern_rev, pattern_revcomp,
                                                      args.max_reads, args.tolerance)
                    outf.write(result_line + "\n")
                    print(result_line)
            print(f"\n结果已成功输出到 {args.outfile}")
        except Exception as e:
            print(f"写入输出文件 {args.outfile} 失败：{e}", file=sys.stderr)
            sys.exit(1)

if __name__ == '__main__':
    main()
