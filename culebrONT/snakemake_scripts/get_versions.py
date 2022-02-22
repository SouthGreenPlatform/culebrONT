import pandas as pd
import argparse
from pathlib import Path
from collections import OrderedDict
from pprint import pprint as pp
import re
version = "0.0.1"

import pandas as pd
from pathlib import Path


def main():
    dir = snakemake.params.dir
    #print(dir)
    output = snakemake.output.csv
    rep = Path(dir)
    rows = []
    paths = Path(rep).glob('*-version.txt')
    for path in sorted(paths):
        file_name = path.stem.split('-')
        tool = file_name[0].split('_')[0]
        dep = file_name[0].split('_')[1] if len(file_name[0].split('_')) > 1 else '-'
        filev = str(path)
        with open(filev, "r") as fd:
            rows.append([tool, dep, fd.readline().rstrip()])
    df = pd.DataFrame(rows, columns=["Tool", "Secondary", "Version"])
    #print(df)
    df.to_csv(output, index=True)


if __name__ == '__main__':
    main()
