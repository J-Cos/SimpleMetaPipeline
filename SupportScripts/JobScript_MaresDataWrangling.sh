

#PBS -lwalltime=48:00:00
#PBS -lselect=1:ncpus=100:mem=1200gb

module load anaconda3/personal 

cd /rds/general/user/jcw120/ephemeral/mares_bar_25Jan2021/mares_bar_taxonomy

ncbitax2lin --nodes-file nodes.dmp --names-file names.dmp


