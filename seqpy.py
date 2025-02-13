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

def extract_matching_reads(file_r1: str, file_r2: str, pattern_fw: str, pattern_comp: str,
                         pattern_rev: str, pattern_revcomp: str, max_reads: int, 
                         tolerance: int, output_dir: str) -> tuple:
    """
    提取匹配的reads对到新的文件中
    返回：提取的reads数量
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 构建输出文件名
    base = os.path.basename(file_r1).replace('_1.fq.gz', '')
    out_r1 = os.path.join(output_dir, f"{base}_matched_R1.fq.gz")
    out_r2 = os.path.join(output_dir, f"{base}_matched_R2.fq.gz")
    
    count = 0
    total_reads = 0
    
    # 判断是否显示进度条
    force_show = any(arg in sys.argv for arg in ("--show-pbar", "-P"))
    if force_show:
        show_pbar = True
    else:
        show_pbar = True
        for arg in sys.argv:
            if arg in ("--noheader", "-n", "--parallel", "-p"):
                show_pbar = False
                break

    pbar = tqdm(total=max_reads, desc=f"Extracting from {base}") if show_pbar else None
    
    with gzip.open(file_r1, 'rt') as f1, \
         gzip.open(file_r2, 'rt') as f2, \
         gzip.open(out_r1, 'wt') as o1, \
         gzip.open(out_r2, 'wt') as o2:
        
        while total_reads < max_reads:
            # 读取R1
            header1 = f1.readline()
            if not header1:
                break
            seq1 = f1.readline()
            plus1 = f1.readline()
            qual1 = f1.readline()
            
            # 读取R2
            header2 = f2.readline()
            seq2 = f2.readline()
            plus2 = f2.readline()
            qual2 = f2.readline()
            
            total_reads += 1
            
            # 检查是否匹配
            seq_lower = seq1.lower().strip()
            if (is_match(seq_lower, pattern_fw, tolerance) or
                is_match(seq_lower, pattern_comp, tolerance) or
                is_match(seq_lower, pattern_rev, tolerance) or
                is_match(seq_lower, pattern_revcomp, tolerance)):
                
                # 写入匹配的reads对
                o1.write(header1)
                o1.write(seq1)
                o1.write(plus1)
                o1.write(qual1)
                
                o2.write(header2)
                o2.write(seq2)
                o2.write(plus2)
                o2.write(qual2)
                
                count += 1
                
            if pbar:
                pbar.update(1)
    
    if pbar:
        pbar.close()
    
    return count, total_reads

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
    每个文件的输出文件名将为 "base_最大记录数.result.txt"
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
            out_file = os.path.join(args.tmpdir, f"{base}_{args.max_reads}.result.txt")
            script_path = os.path.abspath(sys.argv[0])
            cmd = (f"python3 {script_path} --fqfile {file} -f {args.forward} "
                   f"-m {args.max_reads} -t {args.tolerance} -o {out_file} -n")
            if args.show_pbar:
                cmd += " -P"
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
        out_file = os.path.join(args.tmpdir, f"{base}_{args.max_reads}.result.txt")
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
    parser.add_argument("--outfile", "-o", default=None,
                        help="输出结果的 txt 文件，默认 result_{max_reads}_{tolerance}.txt")
    parser.add_argument("--fqfile", "-q",
                        help="指定单个 FASTQ 文件进行处理", default=None)
    parser.add_argument("--noheader", "-n", action="store_true",
                        help="在输出中不包含表头（适用于并行子任务）")
    parser.add_argument("--show-pbar", "-P", action="store_true",
                        help="强制显示进度条（即使在并行模式下也显示，可能导致输出混乱）")
    parser.add_argument("--extract", "-e", action="store_true",
                        help="提取匹配的reads对到新文件")
    parser.add_argument("--outdir", "-O", default="matched_reads",
                        help="提取的reads输出目录，默认为matched_reads")
    
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

    # 若用户未指定 outfile，则自动命名
    if args.outfile is None:
        args.outfile = f"result_seq{args.forward}_{args.max_reads}reads_tol{args.tolerance}.txt"

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

    # 处理所有文件
    if args.parallel:
        run_parallel(args, pattern_fw, pattern_comp, pattern_rev, pattern_revcomp)
    else:
        files_r1 = sorted(glob.glob("*_1.fq.gz"))
        if not files_r1:
            print("当前目录下未找到 _1.fq.gz 文件。", file=sys.stderr)
            sys.exit(1)

        try:
            with open(args.outfile, "w") as outf:
                header_line = ("Filename\tSampled_Reads\tCount_Forward\tPct_Forward\t"
                               "Count_Complement\tPct_Complement\tCount_Reverse\tPct_Reverse\t"
                               "Count_RevComp\tPct_RevComp\tTotal_Matches\tPct_Total\n")
                outf.write(header_line)
                print(header_line, end="")
                
                for file_r1 in files_r1:
                    # 处理统计
                    result_line = process_single_file(file_r1, pattern_fw, pattern_comp,
                                                      pattern_rev, pattern_revcomp,
                                                      args.max_reads, args.tolerance)
                    outf.write(result_line + "\n")
                    print(result_line)
                    
                    # 如果需要提取reads
                    if args.extract:
                        file_r2 = file_r1.replace('_1.fq.gz', '_2.fq.gz')
                        if os.path.exists(file_r2):
                            print(f"\n提取匹配reads: {file_r1} 和 {file_r2}")
                            extracted_count, total_processed = extract_matching_reads(
                                file_r1, file_r2,
                                pattern_fw, pattern_comp, pattern_rev, pattern_revcomp,
                                args.max_reads, args.tolerance, args.outdir
                            )
                            print(f"在处理的 {total_processed} 条reads中提取到 {extracted_count} 对匹配的reads")
                        else:
                            print(f"警告：未找到配对文件 {file_r2}")
                            
            print(f"\n结果已成功输出到 {args.outfile}")
            if args.extract:
                print(f"匹配的reads已保存到 {args.outdir} 目录")
                
        except Exception as e:
            print(f"处理过程出错：{e}", file=sys.stderr)
            sys.exit(1)

if __name__ == '__main__':
    main()
