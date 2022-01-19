import pandas as pd
import argparse
from pathlib import Path
from collections import OrderedDict
from pprint import pprint as pp
import re
version = "0.0.1"

def replace_all(repls, str):
    """
    Function that take a dictionary and text variable and return text variable with replace 'Key' from dictionary with 'Value'.

    :param repls: a python dictionary
    :type repls: dict()
    :param str: a string where remplace some words
    :type str: str()
    :rtype: str()
    :return: - txt with replace 'Key' of dictionary with 'Value' in the input txt

    Example:
        >>> text =  "i like apples, but pears scare me"
        >>> print(replace_all({"apple": "pear", "pear": "apple"}, text))
        i like pears, but apples scare me
    """
    import re
    return re.sub('|'.join(re.escape(key) for key in repls.keys()), lambda k: repls[k.group(0)], str)


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
    
    assembly_list = snakemake.params.assembly_list
    quality_list = snakemake.params.quality_list
    out_repository = snakemake.params.out_dir
    out_stats = snakemake.output.stat
    summary_flagstats = snakemake.input.summary
    
    dico_flagstats_time = AutoVivification()


    for flagstats_file in sorted(summary_flagstats):
        flagstats_file_path = Path(flagstats_file)
        assembler, step, _ = flagstats_file_path.stem.split("-")
        
        step = replace_all({"STEP_":"","_STARTFIXED":""}, step)
        with open(flagstats_file, 'r') as flagstat:
                for line in flagstat:
                    line = line.strip()
                    if not line.startswith('#') or not line.startswith('\n'):
                        # print(line)
                        if 'in total' in line:
                            total = str(re.search("([0-9]+ [\+] [0-9]+) in total", line).groups()).strip('(').strip(')').strip('\',')
                        elif 'secondary' in line:
                            secondary = str(re.search("([0-9]+ [\+] [0-9]+) secondary", line).groups()).strip('(').strip(')').strip('\',')
                        elif 'supplementary' in line:
                            supplementary = str(re.search("([0-9]+ [\+] [0-9]+) supplementary", line).groups()).strip('(').strip(')').strip('\',')
                        elif 'duplicates' in line:
                            duplicates = str(re.search("([0-9]+ [\+] [0-9]+) duplicates", line).groups()).strip('(').strip(')').strip('\',')
                        elif 'mapped (' in line:
                            mapped = str(re.search("([0-9]+ [\+] [0-9]+) mapped", line).groups()).strip('(').strip(')').strip('\',')
                        elif 'mapped (' in line:
                            mapped = str(re.search("([0-9]+ [\+] [0-9]+) mapped", line).groups()).strip('(').strip(')').strip('\',')
                        #elif 'paired' in line:
                        #    paired = str(re.search("([0-9]+ [\+] [0-9]+) paired in sequencing", line).groups()).strip('(').strip(')').strip('\',')


                            #paired = re.search("([0-9]+ [\+] [0-9]+) paired in sequencing", line).groups()
                            #read1 = re.search("([0-9]+ [\+] [0-9]+) read1", line).groups()
                            #read2 = re.search("([0-9]+ [\+] [0-9]+) read2", line).groups()
                            #properly = re.search("([0-9]+ [\+] [0-9]+) properly paired", line).groups()
                            #mate = re.search("([0-9]+ [\+] [0-9]+) with itself and mate mapped", line).groups()
                            #singletons = re.search("([0-9]+ [\+] [0-9]+) singletons", line).groups()
                            #mate_diff_chr = re.search("([0-9]+ [\+] [0-9]+) with mate mapped to a different chr$", line).groups()
                            #mapq5 = re.search("([0-9]+ [\+] [0-9]+) with mate mapped to a different chr (mapQ>=5)", line).groups()


        dico_flagstats_time[assembler][step]['Total'] = f'{total}'
        dico_flagstats_time[assembler][step]['Secondary'] = f'{secondary}'
        dico_flagstats_time[assembler][step]['Supplementary'] = f'{supplementary}'
        dico_flagstats_time[assembler][step]['Duplicates'] = f'{duplicates}'
        dico_flagstats_time[assembler][step]['Mapped'] = f'{mapped}'
        #dico_flagstats_time[assembler][step]['Paired'] = f'{paired}'

    df = pd.DataFrame.from_dict(dico_flagstats_time)
    #df.columns.values[0] = "Assembler"
    #df.columns.values[1] = "Step"
    dataframe_flag = df.T.stack().apply(pd.Series)
    with open(out_stats, "w") as flagstats_out:
        #columns = ['Assembler', "Step", "Total", "Secondary", "Supplementary", "Duplicates","Mapped"]
        #print(f"dataframe_flag:\n{dataframe_flag}\n")
        dataframe_flag.to_csv(flagstats_out, index=True)

if __name__ == '__main__':
    main()
