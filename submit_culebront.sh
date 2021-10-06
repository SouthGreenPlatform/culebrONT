#!/bin/bash

# module help
function help
{
    printf "\033[36m####################################################################\n";
    printf "#                     submit_culebrONT                             #\n";
    printf "####################################################################\n";
    printf "
 Usage: ./submit_culebrONT.sh -c {file} -k {file} -a {additional} -p {profile}
    options:
        -c {file} = Configuration file for run CulebrONT
        -k {file} = Configuration file cluster for run CulebrONT
        -p {path} = Cluster profile
        -a {string} =  additional snakemake arguments
        -h = see help

 Example usage LOCAL:
       ./submit_culebrONT.sh -c config.yaml
       ./submit_culebrONT.sh -c config.yaml -a \"--dry-run --cores 1\"
 Example usage CLUSTER:
       ./submit_culebrONT.sh -c config.yaml -p CulebrONT_SLURM
       ./submit_culebrONT.sh -c config.yaml -p CulebrONT_SLURM -k cluster_config.yaml
       ./submit_culebrONT.sh -c config.yaml -p CulebrONT_SLURM -k cluster_config.yaml -a \"--dry-run --cores 1\"

 At least a config.yaml file is needed to run CulebrONT*
 If you are on LOCAL MODE, -p and -k arguments are NOT needed
 If you are on CLUSTER MODE, please provide the profile path generated according of your favorite scheduler (-p ).
    * You can modify cluster resources giving a cluster_config.yaml (-k) (overwrite profile cluster configuration)\n"
    exit 0
}


##################################################
## Parse command line options.
while getopts c:k:h:p:a: OPT;
    do case $OPT in
        c)    config=$OPTARG;;
        k)    cluster_config=$OPTARG;;
        p)    profile=$OPTARG;;
        a)    additional=$OPTARG;;
        h)    help;;
        \?)    help;;
    esac
done

if [ $# -eq 0 ]; then
    help
fi

##################################################
## main
##################################################

# check config
if [ -n "$config" ] && [ -e "$config" ]; then
  config=$(realpath "$config")
  echo "CONFIG FILE IS $config"
else
  printf "\033[31m \n\n ERROR : you need to provide a valid CONFIG FILE ! \n\n"
  echo "CONFIG FILE IS $config"
  exit;
fi

# cluster_config T et profile T
if [[ -d "${profile}" ]] && [[ -n "${cluster_config}" ]] && [[ -e "$cluster_config" ]]; then
  profile=$(realpath "$profile")
  cluster_config=$(realpath "$cluster_config")
  echo "PROFILE DIR IS" "$profile"
  echo "CLUSTER CONFIG IS" "$cluster_config"
  echo "snakemake -p -s $CULEBRONT/Snakefile \
      --configfile $config \
      --cluster-config $cluster_config \
      --profile $profile \
      $additional"
  snakemake -p -s ${CULEBRONT}/Snakefile \
      --configfile ${config} \
      --cluster-config ${cluster_config} \
      --profile ${profile} \
      ${additional}

# cluster_config F et profile T
elif [[ -d "${profile}" ]] && [[ -z ${cluster_config} ]] && [[ ! -e $cluster_config ]]; then
  profile=$(realpath "$profile")
  echo "PROFILE DIR IS $profile"
  echo "CLUSTER CONFIG IN PROFILE IS ${profile}/cluster_config.yaml"
  echo "snakemake -p -s $CULEBRONT/Snakefile \
      --configfile $config \
      --profile $profile \
      $additional"
  snakemake -p -s ${CULEBRONT}/Snakefile \
      --configfile ${config} \
      --profile ${profile} \
      ${additional}

# cluster_config F, profile F
elif [[ ! $profile ]] && [[ ! $cluster_config ]]; then
  echo "You are running on $(uname -a)"
  echo "snakemake -p -s $CULEBRONT/Snakefile \
  --configfile $config \
  --use-singularity --singularity-args --bind ${test_dir} \
  $additional"
  snakemake -p -s ${CULEBRONT}/Snakefile \
  --configfile ${config} \
  --use-singularity --singularity-args "--bind ${test_dir}" \
  ${additional}

# cluster_config T, profile F
elif [[ ! $profile ]] && [[ $cluster_config ]]; then
  echo "You need verify CulebrONT arguments compatiblity. Please provide path to CLUSTER profile and cluster_config.yaml"
  help
  exit;
fi
