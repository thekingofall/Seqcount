
# FASTQ 序列匹配统计工具

## 功能描述
本工具用于在 `.fq.gz` 格式的 FASTQ 文件中统计特定 DNA 序列及其变体（互补链/反向链/反向互补链）的出现频率，支持以下特性：
- 多模式匹配（精确匹配/容错匹配）
- 并行处理加速
- 进度条显示
- 结果百分比统计

## 环境要求
- Python 3.6+
- 必需 Python 库：
  ```bash
  pip install argparse tqdm
  ```
- 并行处理需要 [ParaFly](https://github.com/ParaFly/ParaFly)（已预装时可跳过）

## 快速开始
```bash
# 单文件处理模式
python fastq_matcher.py -q sample.fq.gz -f "CTATAGCGAAACATCGGCCACCTATA" -m 10000 -t 1 -o result.txt

# 并行处理模式（默认）
python fastq_matcher.py -f "YOUR_SEQUENCE" --parallel -T 8 -d tmp_dir -o final_result.txt
```

## 参数说明
| 参数 | 缩写 | 说明 |
|------|------|------|
| `--forward` | `-f` | 设置目标正向序列（默认：CTATAGCGAAACATCGGCCACCTATA）|
| `--max_reads` | `-m` | 最大处理reads数（默认：10000）|
| `--tolerance` | `-t` | 允许的错配碱基数（默认：0）|
| `--outfile` | `-o` | 输出文件名（默认：result.txt）|
| `--threads` | `-T` | 并行线程数（默认：4）|
| `--tmpdir` | `-d` | 临时文件目录（默认：tmp_parafly）|
| `--no-parallel` |  | 禁用并行处理 |
| `--show-pbar` | `-P` | 强制显示进度条 |

## 工作流程
1. **序列预处理**  
   自动生成四种变体模式：
   - 正向序列（原序列小写）
   - 互补序列
   - 反向序列
   - 反向互补序列

2. **匹配规则**  
   - `tolerance=0` 时进行精确子串匹配
   - `tolerance>0` 时使用滑动窗口计算汉明距离

3. **输出格式**  
   制表符分隔的结果文件包含以下列：
   ```
   Filename  Sampled_Reads  [四种模式的计数/百分比]  Total_Matches  Pct_Total
   ```

## 典型应用场景
1. **引物有效性验证**  
   ```bash
   python fastq_matcher.py -f "YOUR_PRIMER_SEQ" -t 2 --parallel
   ```

2. **接头污染检测**  
   ```bash
   python fastq_matcher.py -f "ILLUMINA_ADAPTER_SEQ" -m 50000 -T 12
   ```

3. **质量控制监控**  
   ```bash
   python fastq_matcher.py -f "CONTROL_SEQUENCE" --no-parallel -P
   ```

## 注意事项
1. 输入文件需符合命名规范 `*.fq.gz`
2. 并行处理时临时目录会自动创建，任务完成后需手动清理
3. 百分比计算结果保留四位小数
4. 当使用高容错阈值时（tolerance > 3），建议增加采样量

> 遇到问题时，可添加 `-P` 参数查看实时处理进度，或检查临时目录中的中间结果文件。

