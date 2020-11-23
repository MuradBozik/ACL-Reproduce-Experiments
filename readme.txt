Author: Yasha Pushak
Created: 2018-08-21

This directory contains data and scripts similar to those we used in the 
2018 PPSN paper: Algorithm Configuration Landscapes: More Benign than 
Expected?

The running time data is exactly the same data that was used for this paper.

The scripts have been refactored, with excess computations and unused 
statistics removed. Due to a time-crunch before the deadline, we massively
paralellized our code in order to get the results to run fast enough, but
this was done in quite a messy way that would make the scripts overly-
complicated to share. Hence we have done the refactoring and removing of
unused code.

Fortunately, this will make the code easier for you to understand and use
in your own applications. Unfortunately, it means the exact random seeds
used to analyze the data presented in the original paper are unavailable.
However, running this refactored code should produce qualitatively identical
results, or at least very nearly so. We have found that in some runs the
bootstrap confidence intervals for the Regions200 response for CPLEX's
mip_strategy_rinsheur are sometimes slightly larger than those in the run
used in the paper. This small increase in size meant that when we performed
our test to see if this parameter response was interesting, it was deemed
unteresting, unlike in the results in the original paper. However, visual
inspection of the parameter response (see supplementary material PDF) shows
that the response for this parameter is clearly very close to being deemed
un-interesting because of how flat it is.

Apart from small deviations like this, which may in some cases cause a 
derived, binary feature such as "interesting" to slightly pass their 
threshold, the results from running this code should be qualitatively 
identical to those in the original paper. The data included here is that
from one such replicate run, which you can either further analyze yourself, 
or use to verify that you are able to correctly run the code yourself. 

----------------------------------------------------------------------------
----------------------------Running the Code--------------------------------

I recommend that you start by running the file exampleEAX.py, which contains
an example of how to run our code on the EAX tsp-rue-1000-3000 scenario.

If you want to reproduce our results, you can run exampleAll.py, which
will (sequentially) run our code on all of the scenarios. Running this
script will take a while. On our machines, it takes just under 43 hours
to process the data in all of the scenarios sequentially.

ACL.py does the bulk of the computations. We have also included helper.py, 
which contains a collection of helper functions that are used throughout our
code. 

If you do run exampleAll.py, you can generate a summary of the results by 
running gatherSummary.py after it is done.

In order to run the code, you will need numpy and scipy installed. 
Otherwise the code only requires standard python libraries.

You can re-create the figure containing the four example parameter response
slices in the paper by running createSlicePlot.plt with gnuplot.
A pdf (which we have already created for you using the replicate data)
will be created called "Parameter-Response-Slices.pdf". 

This code will likely only run on Unix-based systems, it hasn't been tested
on any other operation systems. However, it is also likely that any and all
modifications that may be required to run the code on other operation 
systems are contained in "helper.py", were you should be able to implement
new versions of the "os.system(*)" commands that are used with relative 
ease.

----------------------------------------------------------------------------
------------------------Running your own Experiments------------------------

If you want to run your own, similar experiments, you will need to first 
collect running time data for a bunch of parameter response slices, and
then format it in the same way we have here.

In particular, you will need to create 'scenario' and 'phase' directories.
The phase directory must be a sub-directory of the scenario directory.

Inside of these directories, you will need to create a configs.csv file:
[scenario]/[phase]/configs.csv

It is a comma-separated file containing key-value pairs, one per line.
The key should be a unique ID that identifies a parameter configuration for
your target algorithm. The value, should be the configuration itself.

Our scripts require that you use a particular format for the ID that encodes
additional information about that configuration. In particular, you must
name it as follows:
p[PARAMETER NAME]-[value number]
where 'PARAMETER NAME' is the name of the parameter (response slice) 
corresponding to that configuration, it must be in all caps; and where 
'value number' is an integer (indexed from 0) that indicates the position
of that parameter value in its slice.

For example, imagine that we wanted to create a simple parameter response
slice for EAX, which has two parameter Npop and Nch. Let the parameter
response slices be centered around the point Npop=100 and Nch=30. We will
use 5 parameter values in each slice, for Npop we will use: 
30, 70, 90, 100, 110, 130, 170
and for Nch we will use:
16, 24, 28, 30, 32, 36, 44
Then, the corresponding configs.csv file would be as follows:

pNCH-0, -Nch '16' -Npop '100'
pNCH-1, -Nch '24' -Npop '100'
pNCH-2, -Nch '28' -Npop '100'
pNCH-3, -Nch '30' -Npop '100'
pNCH-4, -Nch '32' -Npop '100'
pNCH-5, -Nch '36' -Npop '100'
pNCH-6, -Nch '44' -Npop '100'
pNPOP-0, -Nch '30' -Npop '30'
pNPOP-1, -Nch '30' -Npop '70'
pNPOP-2, -Nch '30' -Npop '90'
pNPOP-3, -Nch '30' -Npop '100'
pNPOP-4, -Nch '30' -Npop '110'
pNPOP-5, -Nch '30' -Npop '130'
pNPOP-6, -Nch '30' -Npop '170'

You will also need to provide the running time data corresponding to each
of these configurations. The code will expect this to be contained in the
directory [scenario]/[phase]/data in several csv files. These files must
also follow a particular naming convention, that is:

runtimes-[configuration id].csv or .log

e.g., runtimes-pNCH-0.csv 

These files must contain all of the running time data, with one independent
target algorithm run stored on each line. The syntax is:

instance name, seed, instance size, running time, miscellaneous data

The actual values of the seed themselves are not too important, provided 
that they are unique. If you do not have the original values, you may use
0,1,2,3,...

The values of the instance sizes are also not important. They are only
required by our scripts for legacy purposes. However, our code does use
them to sort the instances in the outputed files. If you have a desired
instance order, feel free to use any set of integers to order them here.

The running times should be floats, and you can use 'inf' to represent 
censored runs.

The miscellaneous data is optional. You may either include it or omit it. 
Either way, our code will ignore it. 

See the data provided with this code for example running time files.

Once you have all of those, you're all set! Please refer to exampleEAX.py
for a demo of how to run our code on your data. 

----------------------------------------------------------------------------
---------------------Creating Parameter Response Slices---------------------

We have also included in ACL.py the code that we used to sample parameter
values along the parameter response slices. We expect that you will have a
different means of running and collecting target algorithm running times,
so you're on your own for devising a means to collect this data, and for
setting up where you want to centre your slices. However, you may use our
code for doing the sampling itself. 

In particular, there are two functions for you to use: realRange(), which 
can be used to sample a slice along a real-valued parameter, and
integerRange(), which is for integer-valued parameters.

They both are both used the same way. They require a minimum of 4 parameters
and have an additional 4 parameters that can control the properties of the
outputed slice sample and the verbosity of the code. You four required 
parameters are: 
centre - This point is guaranteed to be in the sample, and it will be at
         the centre of the sampled points.
minVal - The minimum allowed parameter value (inclusive)
maxVal - The maximum allowed parameter value (inclusive)
n      - The number of parameter values included in the sample. Must be odd.

You can generate a sample for EAX's Npop parameter using, e.g.:

>>> integerRange(100,20,1000,15)
[40, 65, 82, 93, 100, 107, 118, 135, 160, 197, 253, 337, 462, 651, 933]

where each integer is the number of points included in the sample. Keep in
mind that this is a randomized procedure, so the results will vary.

There are three other parameters that can control the outputed samples. For
example, there is a parameter p that can be used to control the exponent
used when sampling. Our code uses a default value of p=1.5, but if you
increase p, then your points will be more tightly clustered around the 
centre, and vice versa if you decrease p (note tha p should be greater than 
1).

E.g.,:

>>> ACL.integerRange(100,20,1000,15,p=3)
[43, 81, 94, 97, 98, 99, 100, 101, 102, 103, 106, 119, 157, 271, 614] 
with many values near 100

and 
>>> ACL.integerRange(100,20,1000,15,p=1.1)
[38, 61, 81, 100, 119, 139, 162, 186, 214, 244, 277, 313, 353, 397, 445]
where the numbers are more evenly spaced apart

You can also choose a different "balance" i.e., how the maximum fraction of
points that it attempts to place on one side of the centre (note that if, 
for example, the centre is either the min or max value, then this cannot be
satisfied). We used balance=0.75, but if you want the points more strictly
balanced you can use balance=0.5, e.g.,

>>> ACL.integerRange(100,20,1000,15,balance=0.5)
[30, 55, 71, 82, 90, 95, 98, 100, 102, 105, 110, 118, 129, 145, 170]
where you can see that 100 is exactly in the middle. However, we found
this to be less desirable, since you will notice that the sampled points 
do not obtain as much coverage of the parameter's allowed range.

You can also change the function used to sample the points. We used
and exponential function (which takes in the parameter p). You can define
your own function if you wish. However, we recommend against it. We first
experimented using a polynomial function as well, but we found that
floating point errors caused many of the points near to the centre to in 
fact be the centre. 

Finally, you can also increase the verbosity of the code, if you wish.

Have fun! :) 
