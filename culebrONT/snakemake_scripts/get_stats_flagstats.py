import pandas as pd
from collections import defaultdict, OrderedDict
from pprint import pprint as pp
import re
from pathlib import Path
version = "0.0.1"
pd.set_option("display.precision", 2)


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
    return re.sub('|'.join(re.escape(key) for key in repls.keys()), lambda k: repls[k.group(0)], str)


def parse_flagstat_multi_report(file):
    """
    Take a list of files, parse the data assuming it's a flagstat file
    Returns csv file
    """
    flagstat_regexes = {
        "total": r"(\d+) \+ (\d+) in total \(QC-passed reads \+ QC-failed reads\)",
        "secondary": r"(\d+) \+ (\d+) secondary",
        # "supplementary": r"(\d+) \+ (\d+) supplementary",
        "duplicates": r"(\d+) \+ (\d+) duplicates",
        # "mapped": r"(\d+) \+ (\d+) mapped \((.+):(.+)\)",
        "paired in sequencing": r"(\d+) \+ (\d+) paired in sequencing",
        # "read1": r"(\d+) \+ (\d+) read1",
        # "read2": r"(\d+) \+ (\d+) read2",
        "properly paired": r"(\d+) \+ (\d+) properly paired \((.+):(.+)\)",
        # "with itself and mate mapped": r"(\d+) \+ (\d+) with itself and mate mapped",
        "singletons": r"(\d+) \+ (\d+) singletons \((.+):(.+)\)",
        "mate mapped to a diff chr": r"(\d+) \+ (\d+) with mate mapped to a different chr",
        "mate mapped to a diff chr (mapQ >= 5)": r"(\d+) \+ (\d+) with mate mapped to a different chr mapQ>=5\)",
    }
    parsed_data = {}
    # re_groups = ["passed", "failed", "passed_pct", "failed_pct"]
    re_groups = ["passed", "failed"]

    file_str = Path(file).open("r").read()
    for k, r in flagstat_regexes.items():
        r_search = re.search(r, file_str, re.MULTILINE)
        if r_search:
            for i, j in enumerate(re_groups):
                try:
                    key = "{}_{}".format(k, j)
                    val = r_search.group(i + 1).strip("% ")
                    parsed_data[key] = float(val) if ("." in val) else int(val)
                except IndexError:
                    pass  # Not all regexes have percentages
                except ValueError:
                    parsed_data[key] = float("nan")
    # Work out the total read count
    try:
        parsed_data["(%) mapped_pass"] = f'{(parsed_data["total_passed"]/parsed_data["total"])*100:.2f}%'
        parsed_data["(%) properly_paired"] = f'{(parsed_data["properly paired_passed"]/parsed_data["total"])*100:.2f}%'
    except KeyError as e:
        print(e)
        pass
    return parsed_data


def main():
    assembly_list = snakemake.params.assembly_list
    quality_list = snakemake.params.quality_list
    out_repository = snakemake.params.out_dir
    out_stats = snakemake.output.stat
    summary_flagstats = snakemake.input.summary

    dico_flagstats_time = defaultdict(OrderedDict)

    for flagstats_file in sorted(summary_flagstats):
        flagstats_file_path = Path(flagstats_file)
        assembler, step, _ = flagstats_file_path.stem.split("-")

        step = replace_all({"STEP_":"","_STARTFIXED":""}, step)
        dico_flagstats_file = parse_flagstat_multi_report(flagstats_file)

        dico_flagstats_time[assembler][step] = dico_flagstats_file

    dataframe_mapping_stats = pd.DataFrame.from_dict(dico_flagstats_time)
    dataframe_mapping_stats = dataframe_mapping_stats.T.stack().apply(pd.Series)
    dataframe_mapping_stats.drop(dataframe_mapping_stats.filter(regex='_failed').columns, axis=1, inplace=True)
    dataframe_mapping_stats.columns = [name.replace("_passed", "") for name in dataframe_mapping_stats.columns]
    dataframe_mapping_stats.reset_index(inplace=True)
    try:
        dataframe_mapping_stats.rename({"level_0":'Assembler',"level_1":'Steps' }, axis='columns', inplace=True, errors="raise")
    except KeyError as e:
        print(dataframe_mapping_stats)

    with open(out_stats, "w") as out_csv_file:
        print(f"Library size:\n{dataframe_mapping_stats}\n")
        dataframe_mapping_stats.to_csv(out_csv_file, index=False, header=True, sep="\t")


if __name__ == '__main__':
    main()
