# lilgp_acgp_kernel
ACGP kernel modified for counting functional bloat and improved crossover

In this code we modify the Crossover and Mutation functions to track the count of introns and propogate the count to print to INT file. 
Code also supports multiple subpopulations. We run experiments with two apps 'multiplexer' and 'ant'. 
With the analysis of the introns, we modify the Crossover function to restrict the source and destination nodes to be an intron. 

We have three schemes:

a) Both introns.

b) One of the source/destination to be intron.

c) Both non-introns.

This project led us to conclude that scheme (b) was more optimized with fitness and the tree size of the individuals in the population. Work was done under supervision of Dr. Cezary Janikow for the Course CS6340 at University of Missouri - St. Louis.
