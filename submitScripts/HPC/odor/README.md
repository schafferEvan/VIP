# How to run post-processing steps on a High Performance Computing cluster

## Copy the up-to-date VIP repository from your local computer to the appropriate folder on your cluster

```bash
scp -r ~/VIP/ <username>@<address-to-cluster-submit-node>:~/code/
```

## Navigate to this folder on your cluster

```bash
cd ~/code/VIP/submitScripts/odor/
```

## Submit the jobs in the order of the numbers in the shell scripts. Follow the steps below to make any user-specific changes.

### `1_sumImage.sh`

Notes: 
- `1_sumImage.sh` calls upon `compile_sumImage.m`
- Update `experimentFolder` (line 4) in `compile_sumImage.m` to the path on your cluster that contains your experimental files

### `2_runWatershedSegmentation.sh`

Notes:
- `2_runWatershedSegmentation.sh` calls upon `run_watershed_segmentation.m`
- Update `experimentFolder` (line 4) in `run_watershed_segmentation.m` to the directory on your cluster that contains your experimental files

### `3_extractF.sh`

Notes:
- `3_extractF` calls upon `extract_F_from_conComp.m`
- Update `experimentFolder` (line 4) in `extract_F_from_conComp.m` to the directory on your cluster that contains your experimental files

### `4_align.sh`

Notes:
- `4_align.sh` calls upon `alignImagingAndBehaviorMultiTrial.m`
- Update the paths in lines 9-10 and 12 in `4_align.sh` to point to the correct directories
- Update `traceFolder` (line 4) in `alignImagingAndBehaviorMultiTrial.m` to the `Yproj` directory

### 5_compileRaw.sh

Notes:
- `5_compileRaw.sh` calls upon `compileRaw.py`
- Set up a Python virtual environment for this repository

#### To set up a virtual environment on the cluster, follow the following steps:

##### Creating a new virtual environment:

If you do not already have a directory in your `$HOME` directory where your virtual environments live, first make one.

```bash
mkdir ~/my_envs
```

Then launch an interactive session and load Python. On my cluster, it would look like this:

```bash
srun --pty -t 60:00 -c 4 --mem-per-cpu=8gb -A axs /bin/bash
module load anaconda/3-5.3.1
```

Afterwards, I will create an environment with a reasonably descriptive name. In this case, I will call it `VIP`.

```bash
pyvenv install ~/my_envs/VIP
```

##### To activate this virtual environment:

I run this command.

```bash
source ~/my_envs/VIP/bin/activate
```

When working in this virtual environment, the virtual environment name will be shown in brackets in front of the `user-host-prompt` string.

```bash
(VIP) user@host:~$
```

##### To install modules on the virtual environment:

After activating the environment, installing modules is the same as always. The difference is that the modules will be stored in `/path/to/virtev/lib`. An easy way of install modules is to use `pip`.

Before install modules, the first thing to do is to update pip:

```python
pip install --upgrade pip
```

Then you can install modules as you like, for example numpy:

```python
pip install -Iv numpy
```

##### To deactivate a virtual environment:

Quitting a virtual environment can be done using the command `deactivate`, which was loaded using the `source` command upon activating the virtual environment.

```python
deactivate
```

Source: [https://wiki.anunna.wur.nl/index.php/Virtual_environment_Python_3.4_or_higher](https://wiki.anunna.wur.nl/index.php/Virtual_environment_Python_3.4_or_higher)

##### The correct version of sklearn, skimage, numpy and scipy to be using are:

```python
numpy._version_
'1.17.0'
>>> sklearn._version_
'0.21.2'
>>> skimage._version_
'0.15.0'
>>> scipy._version_
'1.1.0'
```

To install these, you can run the following commands.
```python
pip install -Iv numpy==1.17.0
pip install -Iv scikit-learn==0.21.2
pip install -Iv scikit-image==0.15.0
pip install -Iv scipy==1.1.0
```

- Update `expID` (line 21) to the correct data directory and `expDir` (line 22) to the correct path in `compileRaw.py`.

### 6_postProcessOne.sh

Notes: 

- `6_postProcessOne.sh` calls upon `postProcessOne.py`
- Update `expID` (line 47) to the correct data directory and `rootFolder` (line 48) to the correct path in `postProcessOne.py`. 

### 7_GMMreg.sh

Notes:

- `7_GMMreg.sh` calls upon `GMMreg_toCommonCoords.m`
- Update `experimentFolder` (line 4) and `figureFolder` (line 5) to the correct paths in `GMMreg_toCommonCoords.m`.
- To compile the files need for the Gaussian transformation in `GMMreg_toCommonCoords.m`, run the following steps. 
    
    ```shell
    srun --pty -t 60:00 -c 4 --mem-per-cpu=8gb -A axs /bin/bash
    module load matlab/2019b
    matlab -nodisplay -nodesktop
    mex mex_GaussTransform.c GaussTransform.c -output mex_GaussTransform
    exit
    ```

### 8_appendAlignedCentroids.sh

Notes:

- `8_appendAlignedCentroids.sh` calls upon `append_aligned_centroids.py`
- Update `expID` (line 20), expDir (line 21) and template_file (line 44) to the correct paths in `append_aligned_centroids.py`.

