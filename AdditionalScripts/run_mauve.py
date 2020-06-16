#!/usr/bin/env python3
"""
Takes several final polished and fixstarted assembly fasta files,
create corresponding symbolic link in output directory and
run progressive Mauve to create a multiple alignment of these genomes
"""

import os
import re
import subprocess
import fileinput

fasta_files = snakemake.input.liste
out_dir = snakemake.params.out_dir
xmfa = snakemake.output.xmfa
circular = snakemake.params.circular
dir_quast = snakemake.params.dir_quast
log_e = snakemake.log.error
log_o = snakemake.log.output

if circular:
    common_prefix = os.path.commonprefix(fasta_files)  # find common path
    reduced_paths = [re.sub(f"({common_prefix})", "", path) for path in fasta_files]  # gsub common path to remove it from all inputFiles
    dest_paths = []
    for path in reduced_paths:
        fields = path.split("/")[0:-1]  # split and remove common leave (startfixed_asm.fasta)
        fields.append("asm.fasta")
        label = '-'.join(fields)
        path = os.path.join(out_dir, label)
        dest_paths.append(path)

    for src, dst in zip(fasta_files, dest_paths):  # symlink the fasta files in outDir to have them if need be
        if os.path.isfile(dst):
            os.remove(dst)
        os.symlink(src, dst)

    cmd = ["progressiveMauve", f"--output={xmfa}"] + dest_paths
else:
    cmd = ["progressiveMauve", f"--output={xmfa}"] + fasta_files

#launching MSA
mauve_out = subprocess.run(cmd, check=True, capture_output=True, text=True)

#rename fasta on xmfa
subst = f's\'|{dir_quast}||\'ig' if not circular else f's\'|{out_dir}||\'ig'
print('**SB**', circular, subst)
cmd = f"sed -i {subst} {xmfa}"
os.system(cmd)
print('**CMD**',cmd)

#cmd = ["sed", "-i", subst, xmfa]
#mauve_renaming = subprocess.run(cmd, check=True, capture_output=True, text=True)
#print(mauve_renaming.stdout)
#print(mauve_renaming.stderr)

print(mauve_out.stdout)
print(mauve_out.stderr)