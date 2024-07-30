#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @FileName  :calc_specific_len.py
# @Time      :2024/07/21 01:21:56
# @Author    :Yuchen@rlab

import pysam
import pandas as pd
import numpy as np
import os
import glob
from tqdm import tqdm

def extract_reads_by_length(bam_file, read_lengths):
    """
    提取指定长度的 reads 并返回包含这些 reads 的字典，同时返回原始 BAM 文件中的 mapped reads 数量。
    
    :param bam_file: 输入的 BAM 文件路径
    :param read_lengths: 指定的 reads 长度列表
    :return: 包含指定长度的 reads 的字典，以及原始 BAM 文件中的 mapped reads 数量
    """
    filtered_reads = {length: {'+': [], '-': []} for length in read_lengths}
    filtered_reads['<21'] = {'+': [], '-': []}
    filtered_reads['>24'] = {'+': [], '-': []}

    with pysam.AlignmentFile(bam_file, "rb",threads=12) as bam_in:
        total_mapped_reads = bam_in.mapped
        for read in bam_in:
            if not read.is_unmapped:
                strand = '+' if not read.is_reverse else '-'
                if read.query_length < 21:
                    filtered_reads['<21'][strand].append(read)
                elif read.query_length > 24:
                    filtered_reads['>24'][strand].append(read)
                elif read.query_length in read_lengths:
                    filtered_reads[read.query_length][strand].append(read)
    
    return filtered_reads, total_mapped_reads


def calculate_base_coverage(filtered_reads, total_mapped_reads):
    """
    计算 BAM 文件中的碱基覆盖度，并生成包含深度和 RPM 的 DataFrame。
    
    :param filtered_reads: 包含指定长度的 reads 的字典
    :param total_mapped_reads: 原始 BAM 文件中的 mapped reads 数量
    :return: 包含碱基覆盖度的 DataFrame
    """
    dataframes = {}
    
    for length, strands in filtered_reads.items():
        coverage_data = []
        coverage_dict = {}
        
        for strand, reads in strands.items():
            for read in reads:
                contig = read.reference_name
                if contig not in coverage_dict:
                    coverage_dict[contig] = {'+': {}, '-': {}}
                
                for position in range(read.reference_start, read.reference_end):
                    if position not in coverage_dict[contig][strand]:
                        coverage_dict[contig][strand][position] = 0
                    coverage_dict[contig][strand][position] += 1
            
            for contig, strand_coverage in coverage_dict.items():
                for strand, coverage in strand_coverage.items():
                    for position, depth in coverage.items():
                        rpm = (depth / total_mapped_reads) * 1e6
                        coverage_data.append([contig, position + 1, strand, depth, rpm])
        
        df = pd.DataFrame(coverage_data, columns=['Ref', 'Pos', 'Strand', 'Depth', 'RPM'])
        dataframes[length] = df
    
    return dataframes

def merge_coverage_dfs(coverage_dfs, bam_name):
    """
    横向合并多个覆盖度 DataFrame 成一个 DataFrame，按照指定的表头。
    
    :param coverage_dfs: 一个包含不同长度 reads 的覆盖度 DataFrame 字典
    :param bam_name: BAM 文件的名称，用于填充 "bamName" 列
    :return: 横向合并后的 DataFrame
    """
    df_list = []
    lengths = ['<21', 21, 22, 23, 24, '>24']
    columns = ['Chromosome', 'Position', 'Strand', 'bamName', 'RNA Size', 'RPM']

    for length in lengths:
        if length in coverage_dfs:
            df = coverage_dfs[length]
            df['bamName'] = bam_name
            df['RNA Size'] = length
            df = df.rename(columns={'Ref': 'Chromosome', 'Pos': 'Position'})
            df = df[['Chromosome', 'Position', 'Strand', 'bamName', 'RNA Size', 'RPM']]
            df_list.append(df)
    
    merged_df = pd.concat(df_list, ignore_index=True)
    return merged_df

def calculate_log2_ratio(df_21nt, df_22nt):
    """
    计算两个 DataFrame 对应位置的 log2(21nt/22nt)。
    
    :param df_21nt: 21nt 的 DataFrame
    :param df_22nt: 22nt 的 DataFrame
    :return: 包含 log2(21nt/22nt) 的 DataFrame
    """
    merged_df = pd.merge(df_21nt, df_22nt, on=['Ref', 'Pos'], suffixes=('_21nt', '_22nt'), how='outer')
    merged_df.fillna(0, inplace=True)  # 将NaN值填充为0
    merged_df['log2_21nt_22nt'] = np.log2((merged_df['RPM_21nt'] + 1) / (merged_df['RPM_22nt'] + 1))
    return merged_df[['Ref', 'Pos', 'log2_21nt_22nt']]

if __name__ == "__main__":
    BAM_FILES = glob.glob("/home/chenyc/bioinfo/Project/Cooperative_Project/xulab_tomato_cmv_srna/5_map2cmv/*.bam")
    READ_LENS = [21, 22, 23, 24]
    for bam_file in tqdm(BAM_FILES):
        # 提取指定长度的 reads
        filtered_reads,total_mapped_reads = extract_reads_by_length(bam_file, READ_LENS)

        # 计算碱基覆盖度，并生成包含深度和 RPM 的 DataFrame
        coverage_dfs = calculate_base_coverage(filtered_reads, total_mapped_reads)

        output_file = bam_file.replace('.bam', '_specific_len.xlsx')
        merged_df = merge_coverage_dfs(coverage_dfs, os.path.basename(bam_file).replace('_R1_aligned.bam', ''))
        merged_df.to_csv(output_file.replace('.xlsx', '.csv'), index=False)

        # # 分别获取21nt和22nt的 DataFrame
        # df_21nt = coverage_dfs[21]
        # df_22nt = coverage_dfs[22]

        # # 计算 log2(21nt/22nt)
        # log2_ratio_df = calculate_log2_ratio(df_21nt, df_22nt)

        # # 输出到同一个 Excel 文件的不同 sheet
        # with pd.ExcelWriter(output_file) as writer:
        #     df_21nt.to_excel(writer, sheet_name='21nt_Coverage', index=False)
        #     df_22nt.to_excel(writer, sheet_name='22nt_Coverage', index=False)
        #     log2_ratio_df.to_excel(writer, sheet_name='log2_21nt_22nt', index=False)

        print(f"Results have been written to {output_file}")