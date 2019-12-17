# Propeller Blade Optimizer
This project outlines the airfoil selection, turbine design considerations, and modeling process for a self-starting micro wind turbine.

## Set constraints
Set constraints included a maximum rotor diameter of 30cm and chord lengths between 1.5-5 cm. The design operational wind speed was 10 m/s, and design tip speed ratio was 5.

## Airfoil selection and blade design
The E214 airfoil profile was selected for the blades based on thickness to chord ratio, which was the lowest acceptable value based on an initial stress calculation. Lift and drag coefficient data was extrapolated for post-stall modeling using Viterna’s equations. The airfoil’s performance was modeled using Blade Element Momentum (BEM) Theory. The Glauert empirical relationship, commonly used for the correction of turbulent wake states (a > 0.5), was not used due to abnormal results. Thus, thrust and torque values were skeptically analyzed when reviewing performance predictions. A three-blade turbine, with a chord that decreases linearly from 45mm to 15mm (from root to tip) was chosen for the design.

## Angle of attack and Coefficient of power
Using a relatively low angle of attack distribution to drive better performance at lower tip speed ratios, and the angle of twist was calculated using BEM theory. Though the coefficient of power found for a tip speed ratio of 5 exceeded that of Betz limit, suggesting modeling error, the trends of the plots produced still provided insight as to the design’s performance. Modeling the performance of the turbine over a range of tip speed ratios provided a power coefficient curve that was indicative of an optimal operation at tip speed ratios of 3 to 5.