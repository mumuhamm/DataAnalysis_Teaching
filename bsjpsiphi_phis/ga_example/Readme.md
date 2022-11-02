The cut optimization uses the Genetic Algorithm (GA) from the TMVA library. This script runs the GA in a standalone way.
To build the executable use "make"

In MyGA.cc you have to set the variables for the cuts, allowed ranges for each variable and the input trees.
The input is made by a signal Tree (MC) and a "background" tree wich is usually taken by the data Side Bands. 
Weights can be assigned to both signal and background events to represent the actual number of events expected in the full data sample.

The optimization is based on the "classical" significance (s/sqrt(s+b)), however it is straightforward to change to a custom function.
To run the optimization simply issue ./MyGa

The execution time can vary by few minutes to few hours (even days) depending on the number of variables used and how wide is the allowed range.
For the record the optimization on the JpsiMu sample 2017-2018 can run in ~20 mins, while for the JpsiTkTk, which has much more events, the rnning time can be ~1hour, using 8 variables