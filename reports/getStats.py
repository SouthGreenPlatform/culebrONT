import re
import pandas as pd

assemblers = snakemake.params.liste_assemblers
busco_steps = snakemake.params.liste_busco
out_repository = snakemake.params.out_dir
summary = snakemake.input.summary

results = {}
id = 1
for assembler in assemblers:
    for step in busco_steps:
        with open(f"{out_repository}/ASSEMBLERS/{assembler}/QUALITY/{step}/BUSCO_RESULTS/short_summary_BUSCO_RESULTS.txt", 'r') as busco:
            for line in busco:
                if line.startswith('#'):
                    continue
                if line.startswith('\n'):
                    continue
                col = re.split('\s', line)
                if col[1].startswith("C"):
                    percentage = re.split(",", col[1])
                    CompletePercentage = re.split(":", re.split("%", percentage[0])[0])[1]
                    FragmentedPercentage = re.split(":", re.split("%", percentage[2])[0])[1]
                    MissingPercentage = re.split(":", re.split("%", percentage[3])[0])[1]
                    NumberOfGenes = re.split(":", re.split("%", percentage[4])[0])[1]
                    results[id] = [assembler,step,CompletePercentage, FragmentedPercentage, MissingPercentage, NumberOfGenes]
                    id = id + 1

df = pd.DataFrame(results)
df.index = ['Assembler', 'Step', 'Complete', 'Fragmented', 'Missing', 'Total Genes']
df.to_csv(snakemake.output.stat)
