# T10_Chudnovskiy_Simulation

## Introduction

The [attached code](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/t10Model_nondim.py) and [Latex writeup](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/latexWriteup/main.pdf) contain my reproduction of [Chudnovskiy's 2003 tearing mode simulation on the T-10 tokamak](http://iopscience.iop.org/article/10.1088/0029-5515/43/8/307/meta) and also makes use of a [followup T-10 arictle by Ivanov](https://aip.scitation.org/doi/10.1063/1.4897174) to fill in some of the details.  The original papers on this simulation provided only the original equations and only a few of the associated initial conditions and assumptions.  Many details, including any details on the numerical scheme, were not provided, and I made many interpretive guesses to fill in the blanks.  

## Status

Note that as of the present (Dec. 5, 2018), the code appears to be mostly working.  I've been able to mostly reproduce the feedforward and feedback portions of the simulation.  There are a few differences that are likely attributable to small differences in initial conditions and possibly small bugs.  I've tried altering the q-profile to create the 5 Guass amplitudes in the 2003 paper but have had limited success.


## Initial conditions:

![Initial conditions](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/latexWriteup/images/jAndQ.png)

## Feedforward results:

![Feedforward results](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/latexWriteup/images/feedforwardResults.png)
![Feedforward animation](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/feedforwardAnimation.gif)


## Feedback results:

#### Mode suppression

![Feedforward results, mode suppression](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/latexWriteup/images/feedbackResults.png)


#### Mode amplification

![Feedforward results, mode amplification](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/latexWriteup/images/feedbackResults2.png)

## Todo

TODO(John) Experiment with different ICs to make the results agree better

TODO(John) Cleanup code to make it more user-friendly.  

TODO(John) Figure out why the sign on the feedback gain is wrong.

TODO(John) Double check the non-dimensionalized equations and then implement them

