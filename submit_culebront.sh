#!/bin/bash
#SBATCH --job-name culebrONT
#SBATCH --output slurm-%x_%j.log
#SBATCH --error slurm-%x_%j.log

config="None"
# module help
function help
{
	printf "\033[36m####################################################################\n";
	printf "#                     submit_culebrONT                             #\n";
	printf "####################################################################\n";
	printf "
 Input:
	Configuration file for snakemake

 Exemple Usage: ./submit_culebrONT.sh -c config.yaml

 Usage: ./submit_culebrONT.sh -c {file}
	options:
		-c {file} = Configuration file for run CulebrONT

		-h = see help\n\n"
	exit 0
}


##################################################
## Parse command line options.
while getopts c:h: OPT;
	do case $OPT in
		c)	config=$OPTARG;;
		h)	help;;
		\?)	help;;
	esac
done

if [ $# -eq 0 ]; then
	help
fi

##################################################
## Main code
##################################################

if [ $config != "None" ] && [ -e $config ]; then

	config=`readlink -m $config`
	module load conda
	env=${HOME}/.conda/envs/snakemake
	[ ! -d $env ] && echo -e "## [$(date) - culebrONT]\t Creating conda environment for snakemake" && conda env create -f envs/environment.yaml -n snakemake

	source activate snakemake

	module load singularity
	module load python/3.7
	module load graphviz/2.40.1
	#snakemake --unlock
	echo "CONFIG FILE IS $config"

	# SLURM JOBS WITHOUT PROFILES
	snakemake --nolock --use-conda --use-singularity --cores -p --verbose -s $CULEBRONT/Snakefile \
	--latency-wait 60 --keep-going --restart-times 1 --rerun-incomplete  \
	--configfile $config \
	--cluster "python3 $CULEBRONT/slurm_wrapper.py $config $CULEBRONT/cluster_config.yaml" \
	--cluster-config $CULEBRONT/cluster_config.yaml \
	--cluster-status "python3 $CULEBRONT/slurm_status.py"


	# USING PROFILES
	#snakemake --nolock --use-singularity --use-conda --cores -p --verbose -s Snakefile --configfile config.yaml \
	#--latency-wait 60 --keep-going --restart-times 1 --rerun-incomplete --cluster-config cluster_config.yaml --profile slurm-culebrONT
# if arguments empty
else
	echo "configfile = "$config
	printf "\033[31m \n\n You must add a valid config file !!!!!!!!!!!! \n\n"
	help
fi
