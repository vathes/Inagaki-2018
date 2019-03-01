# Inagaki et al., 2018, 2019

Data pipeline for Inagaki et al., 2018, 2019 from Svoboda Lab.

This project presents a DataJoint pipeline design for the data accompanying the papers:
>Hidehiko K. Inagaki, Lorenzo Fontolan, Sandro Romani & Karel Svoboda. "Discrete Attractor Dynamics Underlies Persistent Activity in the Frontal Cortex" (2019) Nature (https://doi.org/10.1038/s41586-019-0919-7)

>Hidehiko K. Inagaki, Miho Inagaki, Sandro Romani and Karel Svoboda. "Low-Dimensional and Monotonic Preparatory Activity in Mouse Anterior Lateral Motor Cortex" (2018) Jneurosci (https://doi.org/10.1523/JNEUROSCI.3152-17.2018)

The data in original MATLAB format (.mat) have been ingested into a DataJoint data pipeline. 

The data: (Not available)

## Design DataJoint data pipeline 
This repository contains the **Python 3.7** code of the DataJoint data pipeline design for this dataset, as well as scripts for data ingestions and visualization.
 
![Pipeline diagram of intracellular and extracellular](images/erd_from_sess.png)

## Conversion to NWB 2.0
This repository contains the **Python 3.7** code to convert the DataJoint pipeline into NWB 2.0 format (See https://neurodatawithoutborders.github.io/)
Each NWB file represents one recording session. The conversion script can be found [here](scripts/datajoint_to_nwb.py)

## Demonstration of the data pipeline
Data queries and usages are demonstrated in this [Jupyter Notebook](notebooks/Inagaki-2018-examples.ipynb), where several figures from the paper are reproduced. 







