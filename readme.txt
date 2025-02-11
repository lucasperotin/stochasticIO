


############## USAGE ##################
- "g++ simu.cc"
- "./a.out args.txt"



############## MORE DETAILS BELOW ##################

#### argument file structure #####

Each line of the argument file corresponds to one experiment and should be in the form
"sampleFile totBand winLen heuristic typeWin outfile"
- sampleFile is the location of your instance
- totBand is the total bandwidth of the system
- winLen is the size of the window. It should be a number between 0 and 1. Set to 0.5, it means that the window corresponds to 50% of the length of the smaller application if it was alone in the platform.
- heuristic supports "equal", "FIFOCom", "greedyYield", "lookAheadGreedyYield", "Set10", "nextEvLexMin"
- typeWin can take 2 values : "v1", in which case we use the chosen heuristic from start, and the utilization is computed on the whole window
			      "v2", in which case we use the heuristic "equal" and until we switch to our chosen heuristic at the middle of the window. 
					The utilization is computed only on the second half in this case
- outfile is the location of the file where the result will be writen down. Note that it is in append mode so if you run it twice you'll get twice the result.
As of now it simply writes the minimum yield at the end of the window

example: see in args.txt


#### sample file structure ####

Each line of the sample file corresponds to one application and should be in the form
"name bandwidth age progress agePhase i wIter InitialType x_1 x_2 ... x_i"
- name is the name of the application (i.e. an identifier)
- bandwidth is the maximum bandwidth that we can allocate to the application
- age = -release as the window starts at 0
- progress is the time that would have been required if the application was solo. Therefore the initial yield is progress/age and we must have progress<=age
- agePhase is the amount of time that has passed since the current phase was released
- i is the number of phases
- wIter is the average length of an iteration for gopi's generation, put 0 otherwise
- InitialType is 'C' if the first phase is communication and 'W' otherwise
- x_1, ... x_i are the volume of each phase. This is not normalized, for instance for a communication phase, if we have a maximum bandwidth of 0.5, it will take time at least 2 x_k to process

Example : see files/sample.txt

#### Use debug parameter ####
At the begininning of the main, there is a debuga parameter.
It can be usefull to show the state of the applications at every decision point, for debug or for understanding the code. Check out the code for more info.

#####
