# sbcnrnpy
NEURON/Python code for SBC simulations

This code is mainly here for review. If you want to run the model I recommend getting the Example.py to work. Then go from there. Many simulations, even with multiprocessing, run for a very long time. Performance is fantastically poor. Sorry. If you just want to try, use lower dt, lower number of conditions etc. Experiments are usually defined in the __main__ part, for example after line 291 in Figure2.py. Here, the variables "conditions", "conditions_2d", "dur" and "dt" are defined and have strong impact on the duration of the simulation. As a rule of thumb on my machine simulation is about 3-4x slower than real time (i.e. 3-4s for a 1000ms simulation).
This was created to run in a 64bit linux environment. I cannot help with windows or mac.

You need to be able to import neuron in python for this to run. I compiled from source, so here (https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3881#p16729) for some tips. See here (https://neuron.yale.edu/neuron/docs/scripting-neuron-basics) for python + NEURON basics.

You also need the cochlea package (https://github.com/mrkrd/cochlea), best via pip install cochlea. You probably need pandas and other dependencies for cochlea.
