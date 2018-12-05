# T10_Chudnovskiy_Simulation
The [attached code](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/t10Model.py) and [Latex writeup](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/latexWriteup/main.pdf) contains my reproduction of [Chudnovskiy's 2003 tearing mode simulation on the T-10 tokamak](http://iopscience.iop.org/article/10.1088/0029-5515/43/8/307/meta) and also makes use of a [followup T-10 arictle by Ivanov](https://aip.scitation.org/doi/10.1063/1.4897174) to fill in some of the details.   

Note that as of the present (Nov. 27, 2018), the code appears to be mostly working.  I've been able to mostly reproduce the feedforward and feedback portions of the simulation.  There are a few differences that are likely attributable to small differences in initial conditions and possibly small bugs.  I've tried altering the q-profile to create the 5 Guass amplitudes in the 2003 paper but have had limited success.

Feedforward results:

![Feedforward results](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/feedforwardResults.png)
![Feedforward animation](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/feedforwardAnimation.gif)


Feedback results:

![Feedforward results](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/feedbackResults.png)



TODO(John) Experiment with different ICs to make the results agree better

TODO(John) Cleanup code to make it more user-friendly.  

