# Tutorial for Lerner Research Institute High performance Linux cluster at Cleveland Clinic
Research who currently work at Cleveland Clinic and have a Linux account in LRI HPC node can freely use the latest version of 'divine'

1. SSH to HPC's login node.
	```bash
	ssh 10.88.8.18
	```

1. Link python binary and pip into your PATH
	```bash
	mkdir -p $HOME/bin
	ln -s /cm/shared/apps/python/2.7.13/bin/pip $HOME/bin
	ln -s /cm/shared/apps/python/2.7.13/bin/python2.7 $HOME/bin/python
	```

1. Open your bash configuration to define Divine path and environment
	```bash
	vim $HOME/.bashrc
	# add the following lines at the end of the file,
	
	# divine setup
	export PATH=$HOME/.local/bin:$PATH #if not defined
	export PATH=$HOME/bin:$PATH #if not defined

	export DIVINE=/mnt/isilon/data/w_QHS/hwangt-share/apps/divine
	export GCN=$DIVINE/gcn
	export GCN_DATA_DIR=$DIVINE/gcndata
	export GCN_DB_DIR=$DIVINE/gcndb
	export GCN_LOGFILE=$DIVINE/gcn/logs/gcn.log
	export PYTHONPATH=$DIVINE/python_libs/lib/python2.7/site-packages
	export PYTHONPATH=$DIVINE:$PYTHONPATH
	#save the file and exit (e.g., :wq!)
	```

1. Apply divine environment and install some required python modules
	```bash
	source ~/.bashrc
	pip install --user -r $DIVINE/requirements.txt
	python $DIVINE/gcn/bin/prioritize/divine.py --help
	```
	
1. Copy Divine example directory to your Divine project directory to work with
	```bash
	mkdir -p $HOME/lustre/projects/divine
	cd $HOME/lustre/projects/divine
	cp -r $DIVINE/gcn/bin/prioritize/examples ./
	cd examples
	
	# run the existing example scripts as described in 
	# https://github.com/hwanglab/divine/blob/master/documents/tutorial/divine_tutorial.md
	# for example,

	bash ./runme_pfeisffer.sh
	
	# check the result in the directory
	cd ./Pfeiffer
	```
1. Modify the example bash scripts for your sbatch script