#!/bin/bash

# 默认值
input_pattern="*gz"
output_file_arg=""
patterns_string=""
default_output_dir="CountFold"
lines_to_process_param=""
default_lines_to_process="100000"
num_parallel_jobs_param="" # 从 -j 选项获取
sequence_description_param="" # 从 -d 选项获取
default_sequence_description="未说明序列名字" # 新增列的默认值

# 函数：获取DNA序列的反向互补序列
get_reverse_complement() {
    local seq=$1
    local complemented_seq=$(echo "$seq" | tr 'ATCGatcg' 'TAGCTAGC')
    local reverse_complemented_seq=$(echo "$complemented_seq" | rev)
    echo "$reverse_complemented_seq"
}

# 函数：处理单个gz文件（将在后台运行）
process_single_file_in_background() {
    local gz_file="$1"
    local current_sample_name
    local current_cmd_prefix_for_grep 
    local current_denominator_for_percentage
    local -a current_pattern_array 
    local current_all_fwd_count current_all_rc_count
    local original_pattern_iter rc_pattern_iter 
    local cmd_all_fwd_chain_parts cmd_all_rc_chain_parts 
    local local_has_patterns_fwd local_has_patterns_rc
    local cmd_all_fwd_count_full cmd_all_rc_count_full
    local fwd_count_output rc_count_output
    local all_fwd_count_val all_rc_count_val
    local current_all_fwd_percentage_val current_all_rc_percentage_val
    local total_lines_output 

    current_sample_name=$(basename "$gz_file" .gz)
    
    current_cmd_prefix_for_grep="zcat \"$gz_file\""
    if [ ! -z "$effective_lines_to_process" ] && [ "$process_all_lines" -eq 0 ]; then
        current_cmd_prefix_for_grep+=" | head -n $effective_lines_to_process"
    fi

    total_lines_output=$(eval "$current_cmd_prefix_for_grep | wc -l" 2>/dev/null)
    current_denominator_for_percentage=$(echo "$total_lines_output" | awk '{print $1}')
    if ! [[ "$current_denominator_for_percentage" =~ ^[0-9]+$ ]]; then
        current_denominator_for_percentage=0
    fi

    IFS=',' read -r -a current_pattern_array <<< "$patterns_string"
    
    current_all_fwd_count=0
    current_all_rc_count=0

    cmd_all_fwd_chain_parts="$current_cmd_prefix_for_grep"
    local_has_patterns_fwd=0
    for original_pattern_iter in "${current_pattern_array[@]}"; do
        if [ ! -z "$original_pattern_iter" ]; then
            cmd_all_fwd_chain_parts+=" | grep '$original_pattern_iter'"
            local_has_patterns_fwd=1
        fi
    done

    if [ "$local_has_patterns_fwd" -eq 1 ]; then
        cmd_all_fwd_count_full="$cmd_all_fwd_chain_parts | wc -l"
        fwd_count_output=$(eval "$cmd_all_fwd_count_full" 2>/dev/null)
        all_fwd_count_val=$(echo "$fwd_count_output" | awk '{print $1}')
        if [[ "$all_fwd_count_val" =~ ^[0-9]+$ ]]; then
            current_all_fwd_count=$all_fwd_count_val
        fi
    fi

    cmd_all_rc_chain_parts="$current_cmd_prefix_for_grep"
    local_has_patterns_rc=0
    for original_pattern_iter in "${current_pattern_array[@]}"; do
        if [ ! -z "$original_pattern_iter" ]; then
            rc_pattern_iter=$(get_reverse_complement "$original_pattern_iter") 
            cmd_all_rc_chain_parts+=" | grep '$rc_pattern_iter'"
            local_has_patterns_rc=1
        fi
    done

    if [ "$local_has_patterns_rc" -eq 1 ]; then
        cmd_all_rc_count_full="$cmd_all_rc_chain_parts | wc -l"
        rc_count_output=$(eval "$cmd_all_rc_count_full" 2>/dev/null)
        all_rc_count_val=$(echo "$rc_count_output" | awk '{print $1}')
        if [[ "$all_rc_count_val" =~ ^[0-9]+$ ]]; then
            current_all_rc_count=$all_rc_count_val
        fi
    fi
            
    current_all_fwd_percentage_val="N/A"
    current_all_rc_percentage_val="N/A"

    if [ "$current_denominator_for_percentage" -gt 0 ]; then
        current_all_fwd_percentage_val=$(awk -v c="$current_all_fwd_count" -v t="$current_denominator_for_percentage" 'BEGIN {printf "%.2f%%", (c/t)*100}')
        current_all_rc_percentage_val=$(awk -v c="$current_all_rc_count" -v t="$current_denominator_for_percentage" 'BEGIN {printf "%.2f%%", (c/t)*100}')
    fi
            
    # sequence_description 是从父脚本继承的全局变量
    echo -e "$current_sample_name\t$sequence_description\t$patterns_string\t$current_all_fwd_count\t$current_all_rc_count\t$current_denominator_for_percentage\t$current_all_fwd_percentage_val\t$current_all_rc_percentage_val"
}


# 使用帮助函数
usage() {
    echo "用法: $0 -p \"序列1,序列2,...\" [-d \"序列描述\"] [-i \"输入文件匹配模式\"] [-o \"输出文件名或路径\"] [-N <行数>] [-j <任务数>]"
    echo "选项:"
    echo "  -p <序列列表>   : (必需) 逗号分隔的多个原始搜索序列。"
    echo "                    脚本将计算两种共现计数："
    echo "                    1. 全正向共现: 同时匹配所有原始序列的行数。"
    echo "                    2. 全反向互补共现: 同时匹配所有原始序列各自的反向互补序列的行数。"
    echo "  -d <描述>       : (可选) 对当前查询的序列组合的描述性名称。默认为 \"$default_sequence_description\"。"
    echo "  -i <输入模式>   : (可选) 输入 .gz 文件的匹配模式。 默认: \"*gz\""
    echo "  -o <输出文件>   : (可选) 输出 TSV 文件的路径和名称."
    echo "                    如果只提供文件名, 文件将保存在 '$default_output_dir/文件名'."
    echo "                    如果提供路径, 将使用该路径."
    echo "                    如果省略, 结果将默认保存到 '$default_output_dir/YYYYMMDD_HHMMSS_序列摘要[_m行数].tsv'."
    echo "  -N <行数>       : (可选) 只处理每个文件的前 N 行。默认为 $default_lines_to_process 行。如果无效，则处理所有行。"
    echo "  -j <任务数>     : (可选) 并行处理文件的最大任务数。默认尝试使用CPU核心数，否则为4。"
    echo "  -h              : 显示此帮助信息。"
    echo "                    所有表格输出也会同时显示在屏幕上。"
    exit 1
}

# 解析命令行参数
while getopts "p:d:i:o:N:j:h" opt; do # 添加 d: 到 getopts
    case ${opt} in
        p ) patterns_string=$OPTARG ;;
        d ) sequence_description_param=$OPTARG ;; # 处理 -d 选项
        i ) input_pattern=$OPTARG ;;
        o ) output_file_arg=$OPTARG ;;
        N ) lines_to_process_param=$OPTARG ;;
        j ) num_parallel_jobs_param=$OPTARG ;; 
        h ) usage ;;
        \? ) echo "无效的选项: -$OPTARG" 1>&2; usage ;;
        : ) echo "选项 -$OPTARG 需要一个参数。" 1>&2; usage ;;
    esac
done

# 检查必需的参数 -p
if [ -z "$patterns_string" ]; then
    echo "错误: 必须通过 -p 选项提供搜索序列。" 1>&2
    usage
fi

# 确定序列描述 (全局变量)
sequence_description="$default_sequence_description"
if [ ! -z "$sequence_description_param" ]; then
    sequence_description="$sequence_description_param"
fi

# 确定要处理的行数 (全局变量)
effective_lines_to_process=""
process_all_lines=0 
head_param_for_filename="" 

if [ -z "$lines_to_process_param" ]; then 
    effective_lines_to_process="$default_lines_to_process"
    head_param_for_filename="_m${default_lines_to_process}"
elif [[ "$lines_to_process_param" =~ ^[1-9][0-9]*$ ]]; then 
    effective_lines_to_process="$lines_to_process_param"
    head_param_for_filename="_m${lines_to_process_param}"
else 
    echo "警告: -N 的值 '$lines_to_process_param' 不是一个有效的正整数。将对每个文件处理所有行。" >&2
    effective_lines_to_process="" 
    process_all_lines=1 
    head_param_for_filename="" 
fi

# 确定并行任务数 (全局变量)
DEFAULT_MAX_JOBS=4
if [ ! -z "$num_parallel_jobs_param" ]; then
    if [[ "$num_parallel_jobs_param" =~ ^[1-9][0-9]*$ ]]; then
        num_parallel_jobs="$num_parallel_jobs_param"
    else
        echo "警告: -j 的值 '$num_parallel_jobs_param' 不是一个有效的正整数。将使用默认值 $DEFAULT_MAX_JOBS。" >&2
        num_parallel_jobs="$DEFAULT_MAX_JOBS"
    fi
else
    if command -v nproc &> /dev/null; then
        num_parallel_jobs=$(nproc)
    else
        num_parallel_jobs="$DEFAULT_MAX_JOBS"
    fi
fi
echo "信息: 将使用 $num_parallel_jobs 个并行任务处理文件。" >&2


# 确定最终的输出文件路径
final_output_path=""
if [ -z "$output_file_arg" ]; then
    timestamp=$(date +%Y%m%d_%H%M%S)
    sanitized_sequences_part=$(echo "$patterns_string" | tr ',' '_' | sed 's/[^a-zA-Z0-9_-]//g')
    max_seq_part_len=40 
    if [ ${#sanitized_sequences_part} -gt $max_seq_part_len ]; then
        filename_seq_part="${sanitized_sequences_part:0:$max_seq_part_len}_etc"
    else
        filename_seq_part="$sanitized_sequences_part"
    fi
    if [ -z "$filename_seq_part" ] || ([ "$filename_seq_part" == "_etc" ] && [ ${#sanitized_sequences_part} -le ${#filename_seq_part} ]) ; then
        filename_seq_part="patterns"
    fi
    if [ $process_all_lines -eq 1 ]; then 
        head_param_for_filename_actual=""
    else 
        head_param_for_filename_actual=$head_param_for_filename
    fi
    dynamic_default_filename="${timestamp}_${filename_seq_part}${head_param_for_filename_actual}.tsv"
    final_output_path="$default_output_dir/$dynamic_default_filename"
else
    if [[ "$output_file_arg" == */* ]]; then
        final_output_path="$output_file_arg"
    else
        final_output_path="$default_output_dir/$output_file_arg"
    fi
fi

output_directory=$(dirname "$final_output_path")
if [ "$output_directory" != "." ] && [ ! -d "$output_directory" ]; then
    if [ ! -z "$output_directory" ]; then
        mkdir -p "$output_directory" || { echo "错误: 无法创建目录 '$output_directory'." >&2; exit 1; }
    fi
fi

shopt -s nullglob 

# 定义表头
header=$(echo -e "样本\t序列描述\t查询序列组合\t全正向共现计数\t全反向互补共现计数\t总处理行数\t全正向比例(%)\t全反向互补比例(%)")

# 收集所有后台作业的输出
body_unsorted=$(
    active_jobs=0
    for file_to_process in $input_pattern; do
        if [ ! -f "$file_to_process" ]; then
            echo "警告: '$file_to_process' 不是一个有效的文件，已跳过。" >&2 
            continue
        fi

        process_single_file_in_background "$file_to_process" &

        active_jobs=$((active_jobs + 1))
        if [ "$active_jobs" -ge "$num_parallel_jobs" ]; then
            wait -n 
            active_jobs=$((active_jobs - 1))
        fi
    done
    wait 
)

# 将表头和排序后的主体组合输出到 tee
(
    echo "$header"
    if [ ! -z "$body_unsorted" ]; then # 只有在有数据行时才排序和打印
        echo "$body_unsorted" | sort -t $'\t' -k1,1
    fi
) | tee "$final_output_path"

exit 0
