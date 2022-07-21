- Note from Jonathan
This project is archived at the state it was in upon graduation.
Since graduating, I have been happy to continue development on this workflow for while working on the NASA GeneLab team.
You can reach out to me at Jonathan.d.oribello@nasa.gov for the most recent development on this workflow.


# NASA Pipeline: GL-DPPD-7101-C

This is a nextflow implementation of a NASA bioinformatics pipeline for Jonathan Oribello's Bioinformatics Master's Project.

## Installation

Ensure Conda is installed and available (tested with Conda 4.8.3):
<https://docs.anaconda.com/anaconda/install/>

Clone repository, depth is set to only get the a shallow copy.

```bash
git clone --depth=1 https://github.com/J-81/masterProject.git && cd masterProject
```

Create conda environment using package main environment file and activate
```bash
conda env create -f envs/main.yml && conda activate main
```



## Usage

Running with Tower monitoring (Optional, Recommended):
- Login and setup here: [NextflowTower](https://tower.nf)
- **Ensure Nextflow environment variables are set**

```bash
nextflow run main.nf -c config/default.config -with-tower
```

Running with without Tower monitoring (No Extra Setup):

```bash
nextflow run main.nf -c config/default.config
```

## Contributing
Outside contribution is not welcome at this time.

## License
[MIT](https://choosealicense.com/licenses/mit/)
