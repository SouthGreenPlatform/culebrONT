#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
from pathlib import Path
from collections import OrderedDict
from pprint import pprint as pp
version = "0.0.1"
pd.set_option("display.precision", 2)


class AutoVivification(OrderedDict):
    """
    Implementation of perl's autovivification feature.

    Example:

    >>> a = AutoVivification()
    >>> a[1][2][3] = 4
    >>> a[1][3][3] = 5
    >>> a[1][2]['test'] = 6
    >>> print a
    >>> {1: {2: {'test': 6, 3: 4}, 3: {3: 5}}}

    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def main():
    dico_step_files = {}
    for key, list_files in snakemake.params.dico.items():
        dico_step_files[key] = [path.replace("{fastq}",snakemake.params.sample) for path in list_files]
    stat_time = snakemake.output.stat_time
    sample = snakemake.params.sample
    output_dir = Path(snakemake.params.out_dir)
    dico_benchmark_time = AutoVivification()

    for bench_file in dico_step_files["assembly_list"]:
        bench_file_path = Path(bench_file)
        assembler = bench_file_path.stem
        df = pd.read_csv(bench_file, sep="\t")
        dico_benchmark_time["STEP-ASSEMBLERS"][f"ASSEMBLERS"][assembler] = f'{df["s"][0]:.2f}'
     
        # check if MINIASM to load MINIASM_MINIPOLISH file
        if assembler in ("MINIASM"):
            df = pd.read_csv(f"{output_dir}/BENCHMARK/ASSEMBLER/MINIASM_MINIPOLISH.txt", sep="\t")
            df2 = pd.read_csv(f"{output_dir}/BENCHMARK/ASSEMBLER/MINIMAP4MINIASM.txt", sep="\t")
            dico_benchmark_time["STEP-ASSEMBLERS"][f"MINIMAP"][assembler] = f'{df2["s"][0]:.2f}'
            dico_benchmark_time["STEP-ASSEMBLERS"][f"MINIPOLISH"][assembler] = f'{df["s"][0]:.2f}'
        
        # check if SMARTDENOVO to load FASTQ2FASTA file
        if assembler in ("SMARTDENOVO", "SHASTA"):
            df = pd.read_csv(f"{output_dir}/BENCHMARK/ASSEMBLER/FASTQ2FASTA.txt", sep="\t")
            dico_benchmark_time["STEP-ASSEMBLERS"][f"FASTQ2FASTA"][assembler] = f'{df["s"][0]:.2f}'

    if "polishers_list" in dico_step_files:
        for bench_file in dico_step_files["polishers_list"]:
            bench_file_path = Path(bench_file)
            assembler, polisher = bench_file_path.stem.split("_")
            df = pd.read_csv(bench_file, sep="\t")
            dico_benchmark_time["STEP-POLISHER"][f"{polisher}"][assembler] = f'{df["s"][0]:.2f}'

    if "correction_list" in dico_step_files:
        for bench_file in dico_step_files["correction_list"]:
            bench_file_path = Path(bench_file)
            assembler = bench_file_path.stem
            corrector = bench_file_path.parent.stem          
            # check if NANOPOLISH or MEDAKA to load all step file
            if corrector in ("NANOPOLISH", "MEDAKA"):
                summ_step = 0
                summ_variant = 0
                
                variant_list = list(bench_file_path.parent.glob(f"{assembler}*-VARIANTS.txt"))
                # print(f"\n\n{[elm.name for elm in variant_list]}{[elm.parent for elm in variant_list][0]}")
                for bench_file_glob in variant_list:
                    df = pd.read_csv(bench_file_glob, sep="\t")
                    summ_variant += df["s"][0]                
                mean_time = summ_variant/len(variant_list)
                for bench_file_glob in bench_file_path.parent.glob(f"{assembler}*.txt"):
                    if bench_file not in variant_list:
                        df = pd.read_csv(bench_file_glob, sep="\t")
                        summ_step += df["s"][0]
                dico_benchmark_time["STEP-CORRECTION"][f"{corrector}"][assembler] = f'{mean_time+summ_step:.2f}'

    df=pd.DataFrame.from_dict(dico_benchmark_time)
    dataframe_benchmark= df.T.stack().apply(pd.Series)
    with open(stat_time, "w") as benchmark_file:
        with pd.option_context('display.float_format', '{:0.2f}'.format):
        # print(f"dico_benchmark_time:\n{dataframe_benchmark}\n")
            dataframe_benchmark.to_csv(benchmark_file, index=True)


if __name__ == '__main__':
    main()
