# Coriolis-under-ice-positioning
## Coriolis under-ice positioning algorithm :

We implemented the “Terrain-following interpolation for under-ice floats” method presented by Kaihe Yamazaki during ADMT 22.

The method is described in Appendix A of Kaihe Yamazaki et al. article ([https://doi.org/10.1029/2019JC015406](https://doi.org/10.1029/2019JC015406)).

The “Terrain-following” method explained in this paper should be understood before reading the documentation and using the code.

We also tried to improve the method by considering in situ data measured by the float to constrain the algorithm : 
- The average drift depth
- The max depth of the profile 
- The grounded flag 

## Docummentation
[How to use the code](https://github.com/euroargodev/Coriolis-under-ice-positioning/blob/main/estimate_profile_locations_V1.0_20220825.pdf)
