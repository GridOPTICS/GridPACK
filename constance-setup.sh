#unset default python environment 

unset PYTHONPATH

unset PYTHONHOME

 

#java is required for the current simulator

module load java

#set conda path

source /share/apps/python/anaconda3.2019.3/etc/profile.d/conda.sh

#set conda and ray path

export PATH="/share/apps/python/ray/bin:$PATH"

#set up python-3.4 and ray env in conda

source activate /share/apps/python/ray

