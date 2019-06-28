Project for including binary stars in clusters (both open and globular)

Uses cosmic (Breivik) to generate binaries and evolve them for each cluster

Need to create a script to get the velocity dispersions for clusters based on their masses, radii, etc.
Will use these velocity dispersion (sigma values) to find hard-soft period cutoff 

HS Boundary: period at which the binding energy of a binary is higher than the average kinetic energy of stars in a cluster (Geller 2015). 'soft binaries' have longer periods and will be disrupted more easily by another objects gravity. Globular clusters should have shorter cutoff periods (more densely packed -- more interactions) while open clusters should have a longer cutoff period

Will sample binaries below this period and evolve them using cosmic, then pass them through our EBLSST code to find period recovery.
