# Jets in *pp* at the LHC

## Overview:

One uses high energy jets to study quantum chromodynamics. This code consists of 2 parts and for each there is a separate code. The first part is generating a few events using *Pythia8* and storing the information about the different particles in a *.root* file. The second part is using *FastJet3* to find jets in events and plot them in a correlation plot of rapidity and azimuthal angle. 

Disclaimer: This code assumes that the user is working with operating system linux and has installed *ROOT*, *Pythia8* and *FastJet3*. For those, who don't, follow instructions in the following links:

https://root.cern/install/

https://pythia.org/

https://fastjet.fr/

## Pythia
For generation of a few events (100), one uses the *jetGenerator.cc* macro with *jetGenerator.cmnd* as the command file for easier manipulation. One needs to first compile the file with command:

<pre><code> g++ -I/path/to/Pythia/include `root-config --cflags` jetGenerator.cc -o jetGenerator -lpythia8 -L/path/to/Pythia/lib `root-config --libs` </pre></code>

where one needs to specify the paths to include and lib directories. This compiles the *jetGenerator.cc* and creates an executable file *jetGenerator*, but one can still change the input settings in *jetGenerator.cmnd*. Sometimes, there is a problem with finding the missing shared libraries, that problem can be solved with command:

<pre><code> export LD_LIBRARY_PATH=/path/to/Pythia/lib:$LD_LIBRARY_PATH </pre></code>

after that, one should be able to run the generation of events by running the command below.

<pre><code> ./jetGenerator </pre></code>

The output is the *jetGen.root* file, which contains a TTree and a correlation plot of rapidity and azimuthal angle. 

## FastJet