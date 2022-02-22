import pandas as pd
import argparse
from pathlib import Path
from collections import OrderedDict
from pprint import pprint as pp
import re
version = "0.0.1"
pd.set_option("display.precision", 2)


def replace_all(repls, str):
    """
    Function that take a dictionnary and text variable and return text variable with replace 'Key' from dictionnary with 'Value'.

    :param repls: a python dictionary
    :type repls: dict()
    :param str: a string where remplace some words
    :type str: str()
    :rtype: str()
    :return: - txt with replace 'Key' of dictionnary with 'Value' in the input txt

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
    summary_busco = snakemake.input.summary
    
    dico_busco_time = AutoVivification()

    # f"{output_dir}{{fastq}}/ASSEMBLERS/{{assemblers}}/QUALITY/{{quality_step}}/BUSCO/{{assemblers}}_{{quality_step}}_short-summary-BUSCO.txt",

    for busco_file in sorted(summary_busco):
        busco_file_path = Path(busco_file)
        assembler, step, _ = busco_file_path.stem.split("-")
        
        step = replace_all({"STEP_":"","_STARTFIXED":""}, step)
        with open(busco_file, 'r') as busco:
                for line in busco:
                    line = line.strip()
                    if not line.startswith('#') or not line.startswith('\n'):
                        # print(line)
                        if line.startswith('C'):
                            complete, single, duplicate, fragmented, missing, total = re.search(r"C:(\d+\.+\d+)%\[S:(\d+\.+\d+)%,D:(\d+\.+\d+)%\],F:(\d+\.+\d+)%,M:(\d+\.+\d+)%,n:(\d+)",line).groups()

        dico_busco_time[assembler][step]['Complete'] = f'{complete}%'
        dico_busco_time[assembler][step]['Single'] = f'{single}%'
        dico_busco_time[assembler][step]['Duplicate'] = f'{duplicate}%'
        dico_busco_time[assembler][step]['Fragmented'] = f'{fragmented}%'
        dico_busco_time[assembler][step]['Missing'] = f'{missing}%'
        dico_busco_time[assembler][step]['Total Genes'] = f'{total}'
     

    df=pd.DataFrame.from_dict(dico_busco_time)
    dataframe_busco= df.T.stack().apply(pd.Series)
    with open(out_stats, "w") as busco_file:
        with pd.option_context('display.float_format', '{:0.4f}'.format):
        #print(f"dataframe_busco:\n{dataframe_busco}\n")
            dataframe_busco.to_csv(busco_file, index=True)


if __name__ == '__main__':
    main()
