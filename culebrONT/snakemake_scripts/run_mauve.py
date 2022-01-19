#!/usr/bin/env python3
"""
Takes several final polished and fixstarted assembly fasta files,
create corresponding symbolic link in output directory and
run progressive Mauve to create a multiple alignment of these genomes
"""

import os
import subprocess

fasta_files = snakemake.input.liste
out_dir = snakemake.params.out_dir
dir = snakemake.params.dir
xmfa = snakemake.output.xmfa
dir_aggregated = snakemake.params.dir_aggregated
log_e = snakemake.log.error
log_o = snakemake.log.output

dest_paths = []
for src in fasta_files:
    dst = os.path.join(out_dir, os.path.basename(src))
    if os.path.exists(dst):
        os.remove(dst)
    os.symlink(src, dst)
    dest_paths.append(dst)
stringList = ' '.join([str(item) for item in dest_paths])

#recovery version of mauve
version_file = f"{dir}/versions/MAUVE-version.txt"
cmd = f"progressiveMauve --version "
subprocess.run(cmd, shell=True, capture_output=False, stderr=open(version_file, 'w'))

# convert list on string
cmd = f"progressiveMauve --output={xmfa} {stringList}"

#launching MSA
subprocess.run(cmd, shell=True, capture_output=False, stdout=open(log_o, 'w'), stderr=open(log_e, 'w'))

#rename fasta on xmfa
subst = f's\'|{dir_aggregated}||\'ig'
cmd = f"sed -i {subst} {xmfa}"
subprocess.run(cmd, shell=True, capture_output=False)
