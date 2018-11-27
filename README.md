# T10_Chudnovskiy_Simulation
The attached code and Latex writeup contains my reproduction of [Chudnovskiy's 2003 tearing mode simulation on the T-10 tokamak](http://iopscience.iop.org/article/10.1088/0029-5515/43/8/307/meta).  
My code also makes use of a [followup arictle by Ivanov](https://aip.scitation.org/doi/10.1063/1.4897174) to fill in some of the details.   

Note that as of the present (Nov. 27, 2018), the code appears to be mostly working.  I've been able to mostly reproduce the feedforward portion of the simulation.  There are a few differences that are likely attributable to small differences in initial conditions and possibly small bugs.

Feedforward results:
![Feedforward results](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/feedforwardResults.png)
![Feedforward animation](https://github.com/jwbrooks0/T10_Chudnovskiy_Simulation/blob/master/feedforwardAnimation.gif)

TODO(John) Does the slow time scale need a 2pi ?

TODO(John) Implement feedback portion of code
