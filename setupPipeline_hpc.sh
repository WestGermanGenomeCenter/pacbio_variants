
# setup the smk8 compatible workflow
# sometimes you need to 'module load condaforge' before the next command works
cd /gpfs/project/projects/bmfz_gtl/software/pb_variants
conda env create -n smk9 -f smk9.yaml
conda activate smk9
mkdir -p ~/.config
mkdir -p ~/.config/snakemake

user="$(whoami)" # create clusterLogs_$user, edit the profile to dump logs there aswell

mkdir /gpfs/project/projects/bmfz_gtl/software/pb_variants/clusterlogs_$user -p 




cp -r /gpfs/project/projects/bmfz_gtl/software/snakemake_profile/pbs_pacbio ~/.config/snakemake/. # this only works for HPC HILBERT users

sed -i "s/clusterLogs/clusterlogs_$user/g" ~/.config/snakemake/pbs_pacbio/config.yaml # so that for each user there is a specific dir with clusterlogs