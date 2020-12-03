This code is mainly here for review. If you want to run the model please follow the instructions below.
Many simulations, even with multiprocessing, run for a very long time. As a rule of thumb on my machine simulation is about 3-4x slower than real time (i.e. 3-4s for a 1000ms simulation). This was created to run in a 64bit linux environment. I cannot really help with windows or mac but just contact me if you have questions, maybe I can help.

# Installation Instructions

Here are installation instructions that should allow testing of the smallexc
model as used for the J Neurophysiol manuscript (Koert & Kuenzel, "Small dendritic 
synapses enhance temporal coding in a model of cochlear nucleus bushy cells").
This is for linux systems, specifically ubuntu 18.04.

1) unzip the smallexc repository wherever you keep your model code, in this example it
could be in:

```~/Documents/Python/KoertKuenzelSBC```

2) create a python3 virtual environment:

```python3 -m venv smallexc_venv```

3) activate the venv:

```source smallexc_venv/bin/activate```

4) install packages that are required (NEURON is included in this!)

```pip3 install -r requirements.txt```

For some reason I have to do this twice, possibly due to dependency issues. I am still investigating.

5) change to the smallexc dir

```cd ~/Documents/Python/KoertKuenzelSBC```

6) Run the test stimulation. The actual code is in the subdirectory "smallexc". 
Simulations "expect" to be run from the root level of the repository:

```python smallexc/smallexc_test.py```

After a few seconds a figure should appear with some traces. Also, the test will check
whether the subdirectories figs and results are present, that most of the model programs
expect to exist.
If this succeeds you should be ready to reproduce all figures of the paper.

# Instructions for specific simulations
Here I will list the commands (and some explanations) necessary to reproduce any figure of the paper.
All "expect" to be started from the root directory of the repository. Most produce data-output in the result directory. Figures are saved as pdf in the figs subdirectory.

Do not forget:

* activate the virtual environment
* cd to the root of the repository

## Figure1
This figure combines results from two simulations:
#### Figure1 panels B-D
```python smallexc/Figure1.py```
#### Figure1 panels E-F
So for a quick glance use:

```python smallexc/impedance_testing_stylized.py 1 3 4```

For the paper we used:

```python smallexc/impedance_testing_stylized.py 1 9 4```
## Figure2
#### Figure2 panels A-E
```python smallexc/Figure2.py```

Calculating the 2D-Figures will take a lot of time if done in the same resolution as we did for the paper. A conservative estimate is that the simulation runs ~10% realtime. Thus one condition of 5sec will take one core about 30-60s. Please adjust number of conditions, duration of individual runs and dt in lines 328-332 of the code, if you just want to get a first glance. Also, adjust the cores to the number of logical cores you have. In parallel simulations runtime (using python multiprocessing) scales linear with the number of cores, but running more jobs than cores seems to be very unefficient.

To get a first impression try the settings below (lines 328ff of Figure2.py):

```Python3
conditions = 8
conditions_2d = 5
cores = 4
dur = 1000.0
dt = 0.025
```

#### Figure2 panels F-G
This is part of the result of "impedance_testing_stylized.py" described for Figure1, panels E-F above.

## Figure3
```python smallexc/Figure3.py```

See description of Figure 2, panels A-E above.

## Figure4
#### Figure4 panels A-E
```python smallexc/Figure4.py```

See description of Figure 2, panels A-E above.
#### Figure2 panel F
```python smallexc/Figure4F.py```

Adjust settings in line 160-164 if you want a first impression.

## Figure5
#### Figure5 panels A-G
```python smallexc/Figure5.py```

See description of Figure 2, panels A-E above.

#### Figure5 panels H-J
This is part of the result of "impedance_testing_stylized.py" described for Figure1, panels E-F above.

## Figure6
#### Figure6 panels B-E
```python smallexc/Figure6.py```

Calculates repetitions of one frequency/level/cf combination. This is quite quick. Adjust stimulus and simulation settings in lines 31-37.

#### Figure6 panels F-G
```python smallexc/Figure6_partB.py```

Despite running 2D data this should be rather fast. Still, you can adjust the settings of the simulation in lines 275-279.

#### Figure6 panels G5&G6
```python smallexc/Figure6_partC.py```

Should take about the time as Figure6_partB.py

Adjust settings in lines 193-197.

## Figure7
#### Figure7, panels A-I
```python smallexc/Figure7.py```

This is controlled by command line parameters passed to the program. The first parameter is the cf of the unit to simulate, the second parameter (0 or 1) controls whether the synapses should be on the dendrite (0, default) or the soma (1). For the data in the figure (this is the default...) you pass:

```python smallexc/Figure7.py 1500.1 0```

Per default 25 frequencies and 25 levels are presented in the grid, running on 20 cores. This can be adjusted in lines 511-513.

#### Figure7, panels J-K
```python smallexc/Figure7.py 500.1 0```

```python smallexc/Figure7.py 1000.1 0```

```python smallexc/Figure7.py 1500.1 0```

```python smallexc/Figure7.py 2000.1 0```

```python smallexc/Figure7.py 2500.1 0```

See above for explanations.

#### Figure7, panels J-K, entrainment indices
The entrainment data was added as an afterthought upon reviewer comments. It is generated for the stylized AND one 3D model for one CF at one level (but many frequencies) by:

```python smallexc/entrainmenttest.py 0 25 1500.1 77```

The parameters to pass are: load data (0 or 1), number of frequencies, cf, db SPL. For the data in Figure 7 and 11 we used:

```python smallexc/entrainmenttest.py 0 25 500.1 77```

```python smallexc/entrainmenttest.py 0 25 1000.1 77```

```python smallexc/entrainmenttest.py 0 25 1500.1 77```

```python smallexc/entrainmenttest.py 0 25 2000.1 77```

```python smallexc/entrainmenttest.py 0 25 2500.1 77```


## Figure8
#### Figure8, panels A
```python smallexc/Figure8a.py 0 25 20 1389.0 1501.0 65.0 1```

The parameters to pass, in order, are: weload, nconds, ncores, freq, cf, db, examples

* weload: whether to load existing data (1) or force rerun (0)
* nconds: number of conditions
* ncores: number of cores
* freq: stimulus frequency in Hz
* cf: cf of unit to simulate in Hz
* db: stimulus level in dB SPL
* examples: whether to generate (and save) example traces (0 or 1)

This program appends results to a textfile "results/Figure.txt". We (manually) entered these results for a range of different stimulus conditions into table ("results/Figure8a_data.csv"). This data is read and analyzed by:

```python smallexc/Figure8a_extradata.py```

The outcome of this was only presented in text-form in the paper.

#### Figure8, panels B
```python smallexc/Figure8b.py```

Simulation parameters (number of conditions, cores etc) can be adjusted in lines 289-292.


#### Figure8, panels C
```python smallexc/Figure8c.py```

Simulation parameters (number of conditions, cores etc) can be adjusted in lines 242-245.

#### Figure8, panel D
This part of the data generated by simulations in:

```python smallexc/impedance_testing_balancing.py 1 25 20 3000.0 0.01 1```

Control by input parameters: load, conds, cores, dur, dt, mp

* load: 1 loads existing data, 0 reruns no matter what
* conds: how many nxn conditions to test
* cores: how many jobs to spawn
* dur: duration of chirp
* dt: sampling interval
* mp: 1 uses parallel processing

## Figure9
#### Figure9, panels A
```python smallexc/Figure9a.py```

Simulation parameters (number of conditions, cores etc) can be adjusted in lines 303-306.

#### Figure9, panels B
```python smallexc/Figure9b.py```

Simulation parameters (number of conditions, cores etc) can be adjusted in lines 214-217.

## Figure10
#### Figure10, panels D/E/G/H
```python smallexc/impedance_testing_3d.py 0 9 20 1```

```python smallexc/impedance_testing_3d.py 0 9 20 2```

```python smallexc/impedance_testing_3d.py 0 9 20 3```

You can pass arguments to this function as in:

```python3 impedance_testing_3d.py weload nconds ncores cell```

* weload: 1 = load data instead of calculating new (use 0 or anything to recalc)
* nconds: Nr. of conditions varied for Ra (processing time)
* ncores: how many processes you'd like to spaw
* cell: which cell (only 1,2 or 3 at the moment (see semodels.BS_3D))

#### Figure10, panels F
```python smallexc/impedance_testing_3d_summary.py```

After you ran the above for all three example cells (the function will just fail if you did not!) this will read the result data and generate the summary data presented in these panels.

## Figure11
#### Figure11, panels A-I
```python smallexc/Figure11.py 0 25 20 1 1500.1 0 0.064```

This is controlled by command line parameters passed to the program. In order they are:

* load: 1 loads existing data, 0 reruns no matter what
* conds: how many nxn conditions to test
* cores: how many jobs to spawn
* cell: which cell (only 1,2 or 3 at the moment (see semodels.BS_3D))
* cf: cf of the unit to simulate
* synapsemodel: 0 places synapses randomly in the dendrite, 1 places all synapses on the soma.
* gsyn: total conductance of the dendritic synapses

#### Figure11, panels J-K
```python smallexc/Figure11.py 0 25 20 1 500.1 0 0.064```

```python smallexc/Figure11.py 0 25 20 1 1000.1 0 0.064```

```python smallexc/Figure11.py 0 25 20 1 1500.1 0 0.064```

```python smallexc/Figure11.py 0 25 20 1 2000.1 0 0.064```

```python smallexc/Figure11.py 0 25 20 1 1500.1 0 0.064```

```python smallexc/Figure11.py 0 25 20 2 1500.1 0 0.064```

```python smallexc/Figure11.py 0 25 20 3 1500.1 0 0.064```

```python smallexc/Figure11.py 0 25 20 1 1500.1 1 0.064```

See above for explanations.

#### Figure11, panels J-K, entrainment indices
```python smallexc/entrainmenttest.py 0 25 1500.1 77```

See Figure 7 for details.
