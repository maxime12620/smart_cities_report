# Dynamic Models for Building Energy Management

The notebooks can be run interactively on MyBinder.com by clicking on the button below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/cghiaus/dm4bem/HEAD)

## Tutorials
- t01: Read weather data and solar radiation from EnergyPlus files.
- t02: Model and dynamc simulation of a simple wall with capacities in all nodes.
- t03: Model and dynamic simulation of a cubic building with proportional feed-back control.
- t04: Assembling thermal circuits.
- t05: Heating and cooling with dead-band.
- t06: Control by changing the inputs during the itegration in time.
- t07: Free-cooling by switching between models and changing the inputs during the integration.

## Sessions
1.	**Model**
-	Plan of a simple building.
-	Hypothesis for boundary conditions, windows, doors, and wall composition.
-	Write down the adjancy matrix **A**, the conductance matrix **G** and the capacity matrix **C**.
2.	**Pyhton implementation**
-	Implement matrices **A, G,** and **C**.
-	Calculate the solar flows.
-	Write the input vectors b and f in time.
-	Simulate a step response
-	Simulate the response to weather data.
3.	**Debugging and optimization**
-	Debug.
-	Complex controllers (dead-band, model predictive control).
4.	**Written report and publication**
-	Write report in Jupyter notebook.
-	Publish the report on GitHub and MyBinder

**References**

1. C. Ghiaus (2013). Causality issue in the heat balance method for calculating the design heating and cooling load. *Energy* 50: 292-301
[DOI 10.1016/j.energy.2012.10.024](http://dx.doi.org/10.1016/j.energy.2012.10.024), [HAL 03605823]( https://hal.archives-ouvertes.fr/hal-03605823/document)

2. C. Ghiaus, N. Ahmad (2020). Thermal circuits assembling and state-space extraction for modelling heat transfer in buildings, *Energy*, 195:117019
[DOI 10.1016/j.energy.2020.117019](https://doi.org/10.1016/j.energy.2020.117019), [HAL 03600778](https://hal.archives-ouvertes.fr/hal-03600778/document)

3. C. Ghiaus (2021). Dynamic Models for Energy Control of Smart Homes, in *S. Ploix M. Amayri, N. Bouguila (eds.) Towards Energy Smart Homes*, Online ISBN: 978-3-030-76477-7, Print ISBN: 978-3-030-76476-0, Springer, pp. 163-198 (ref.)
[DOI 10.1007/978-3-030-76477-7_5](https://doi.org/10.1007/978-3-030-76477-7_5), [HAL 03578578](https://hal.archives-ouvertes.fr/hal-03578578/document)

4. J. Kneifel (2013). Annual Whole Building Energy Simulation of the NIST Net Zero Energy Residential Test Facility Design, *NIST Technical Note 1767*, [DOI 10.6028/NIST.TN.1767](https://doi.org/10.6028/NIST.TN.1767)

[Teaching support:](https://filesender.renater.fr/?s=download&token=f1a3d994-0efc-4abe-bf3e-ea8085a1b22e) (link valid till 13 April 2022)

# Questions for exam
1. Defintion of science.
2. Reproducibility crises: definition and how to overcome it.
3. Definition of physical and computational causality.
4. Conservation lwas: 2 examples.
5. Relation between conservation laws and symmetry in physics.
6. Constitutive laws: one example of a univseral law and one of a phenomenological law.
7. Explain why id the SI system of units there are only seven fundamantal units.
8. What is the difference between the classical SI system of units and the newly addopted one? Why is this differnce important?
9. What is the relationship between energy and temperature?
10. Draw thhe basic netword for heat trasfer modelling. Explain each element of the network: temperaure nodes, flow branches, conductances, capacities, temperature sources, flow sources.
11. Draw the framework for obtaining the *difussion equation*.
12. Show the analogy between:
    - heat transfer
    - mass transfer
    - electrical conduction
13. Define the modes of heat transfer and give the expression of conductance for:
    - conduction
    - convection
    - radiation
    - advection.
14. Conservation of energy in steady-state and in dynamics.
15. Definition of sensible heat.
16. Surface phenomena and volume phenomena in energy balance equation.
17. Costitituve law for:
    - heat conduction,
    - heat convection,
    - heat advection,
    - radiative heat exchnage.
18. Linear form of heat transfer by radiation.
19. Draw a wall and a window. Make a thermal network model of this system.
20. Difference between *Differenctial Algebraic Equations* model ans *state-space* representation.
