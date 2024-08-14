# TE-analysis-for-multi-pathway-interaction

This readme file will help replicate the results in the paper titled: "Disentangling coexisting sensory pathways of interaction in collective dynamics.”
The directory “code” contains all the codes and the directory "data" contains the sample datasets required to replicate the results of this work. 

The programs and functions are detailed below:
## “Analysis.m”:
This MATLAB code processes the data and calculates transfer entropy (TE) and conditional TE for different delays between time series and plots all the relevant results of TE analysis for an airfoil-flag system coupled via two pathways of interaction (hydrodynamics + electromechanical). 
## “Analysis_noMech.m”:
This MATLAB code processes the data and calculates TE and conditional TE for different delays between time series and plots all the relevant results of TE analysis for an airfoil-flag system coupled via a single interaction pathway (hydrodynamics).
## “Seasonally_adjust.py”:
This Python code conducts the seasonal adjustment to remove seasonal trends from the relevant timeseries.
## “symbolize_data.m”:
This MATLAB function calculates a symbolized timeseries based on the input timeseries and embedding dimension (only for m=2,3,4).
## “permutate.m”:
This MATLAB function randomly shuffles one time series (x) while maintaining the dynamics between two other time series (y,z) intact.
## “compute_entropy.m”:
This MATLAB function computes the Shannon entropy of a set of random variables.
## “transfer_entropy_delay.m”:
This MATLAB function calculates the transfer entropy (TE) from a source Y to a target X for a given time delay (del).
## “cond_transfer_entropy_delay.m”:
This MATLAB function calculates the conditional TE from a source Y to a target X for a given time delay (del), accounting for a third time series Z with a time delay (del2) with respect to target X.
## “downsample.m”:
This MATLAB function downsamples a timeseries (u) recorded at time instances (t) to udown at time instances tdown.
