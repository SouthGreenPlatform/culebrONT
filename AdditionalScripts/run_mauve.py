#!/usr/bin/env python3
"""
Takes several final polished and fixstarted assembly fasta files,
create corresponding symbolic link in output directory and
run progressive Mauve to create a multiple alignment of these genomes
"""

import os
import re
import subprocess

fasta_files = snakemake.input.liste
out_dir = snakemake.params.out_dir
xmfa = snakemake.output.xmfa
dir_aggregated = snakemake.params.dir_aggregated
log_e = snakemake.log.error
log_o = snakemake.log.output

# find common path of fasta files
common_prefix = os.path.commonprefix(fasta_files)
reduced_paths = [re.sub(f"({common_prefix})", "", path) for path in fasta_files]  # gsub common path to remove it from all inputFiles

dest_paths = []
# remove symbolic links from another run
for path in reduced_paths:
    fields = path.split("/")[0:-1]  # split and remove common leave (startfixed_asm.fasta)
    fields.append("asm.fasta")
    label = '-'.join(fields)
    path = os.path.join(out_dir, label)
    dest_paths.append(path)

# recovery path of renamed symbolic files
for src, dst in zip(fasta_files, dest_paths):
    if os.path.isfile(dst):
        os.remove(dst)
    os.symlink(src, dst)
stringList = ' '.join([str(item) for item in dest_paths])

# convert list on string
cmd = f"progressiveMauve --output={xmfa} {stringList}"

#launching MSA
#os.system(cmd)
subprocess.run(cmd, shell=True, capture_output=False, stdout=open(log_o, 'w'), stderr=open(log_e, 'w'))

#rename fasta on xmfa
subst = f's\'|{dir_aggregated}||\'ig'
cmd = f"sed -i {subst} {xmfa}"
#os.system(cmd)
subprocess.run(cmd, shell=True, capture_output=False)
