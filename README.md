# adrian_striatum_analysis
*This repository contains the Matlab code that provides an interface to the database of chronic rat Neuropixels recordings performed by Adrian Bondy in various striatal subregions*

### What dataset are we talking about?
The dataset consists of several dozens of Neuropixels recordings from probes implanted chronically in rats performing the Poisson Clicks task. The recordings were performed to target a variety of sites in the striatum.

### How to find the data
The actual data resides in `\\bucket.pni.princeton.edu\brody\abondy\adrian_striatum_analysis`. To work locally for faster database access, copy this folder to some local directory `my_path` and then edit the file `get_parameters.m` such that the variable `pc_data_path` is set to `my_path`:<br>
i.e. replace <br>
`P.pc_data_path = fullfile('X:','abondy','adrian_striatum_analysis');` <br>
with<br>
`P.pc_data_path = my_path;`

There are two ways you can view an index of the data:
  1. Run `paths = get_data_paths()` to get a cell array of strings listing the paths to the datafiles ("Cells" files) corresponding to individual sessions. These files include "all" the data, include the spike times and associated behavioral event times. See next section for a detailed description of a "Cells" file.
  2. Run `T = load_cells_table()` to load into the workspace a table with a row for each individual cell in the database, and columns describing features of those cells, including:
    - d

  This table does not contain spike times but is useful for providing an overview of the properties of the cells in the database, for example the 3-D anatomical location of the entire database of recorded cells, the distribution of waveform properties across the dataset, or the contribution to the database of particular rats.

### "Cells" files
A "Cells" file is the structure that will loaded into your workspace when loading one of the datafiles:<br>
`Cells = load(paths{i})`

These contain a large number of fields providing information about the recording session. The key ones are:
 -
 -
 -

### Dependencies
There are a number of dependencies in this repository on other Brody Lab repositories, namely `npx-utils` and `labwide_pbups_analysis`.
