# Riemann-Solver

This program solves the 1D hydrodynamical Riemann problem given any inital state. For more infomation on the Riemann problem refer to the [wikipedia](https://en.wikipedia.org/wiki/Riemann_problem#:~:text=A%20Riemann%20problem%2C%20named%20after,in%20the%20domain%20of%20interest) page. To use this code the makefile should be ran and then the program can be started in the terminal as follows 
```
make all

./solver p_L P_L v_L p_R P_R v_R
```
where L signifies the left state and R signifies the right state. Also p = density, P = pressure, and v = velocity. It's important to note that this code is preset to deal with an adiabatic gas (gamma = 5/3). 

Below is an example of the programs output: 
<br/>
Inputs:
<br/>
&nbsp; Left State: ( p_L, P_L, v_L ) = ( 2.0, 1.0. -2.0 )
<br/>
&nbsp; Right State: ( p_R, P_R, v_R ) = ( 1.0, 0.4, -0.2 )
<br/>
Outputs:
<br/>
&nbsp; Intermediate State: ( p_1, P_1, v_1, p_2 ) = ( 0.424, 0.0755, -0.895, 0.368 ) 
&nbsp; NOTE: P_1 = P_2 and v_1 = v_2 this is always the case
<br/>
&nbsp; Shock Speeds: No shocks 
<br/>
&nbsp; Graphical solution:
  
  ![Sol 4](https://user-images.githubusercontent.com/60577496/114378234-e4556e80-9b5d-11eb-85c5-c2df8a6030a8.png)
  
  This program utilizes GNUplot inorder to generate the plots.
