# Modelling of a Thermal Energy Storage using Schumann Model

We model a packed bed thermal energy storage(TES) in MATLAB using a Schumman model. It is a one-dimensional model i.e. the temperature gradient along the radius of the packed bed is assumed to be constant and only temperature gradient along the lateral axis is considered. 
Schumann model consists of two equations, one for the solid storage material and one for the fluid flowing through the packed bed. The two equations are listed below:

<p align="center">
<img src="https://latex.codecogs.com/svg.image?\large&space;\frac{\partial&space;T_s}{\partial&space;t}&space;=\frac{h_v}{\rho_s\&space;c_s\&space;(1-\epsilon)}&space;&space;(T_a-T_s)&plus;&space;\frac{k_{s,eff}}{\rho_s&space;c_s&space;(1-\epsilon)}&space;&space;\frac{\partial^2T}{\partial&space;x^2}" title="https://latex.codecogs.com/svg.image?\large \frac{\partial T_s}{\partial t} =\frac{h_v}{\rho_s\ c_s\ (1-\epsilon)} (T_a-T_s)+ \frac{k_{s,eff}}{\rho_s c_s (1-\epsilon)} \frac{\partial^2T}{\partial x^2}" />
</p>

<p align = "center">

<img src="https://latex.codecogs.com/svg.image?\large&space;\frac{\partial&space;T_a}{\partial&space;t}&space;&plus;&space;\frac{G}{\rho_a&space;\epsilon}&space;\frac{\partial&space;T_a}{\partial&space;t}&space;=&space;\&space;\frac{h_v}{\rho_s&space;c_s\epsilon}&space;(T_s&space;-&space;T_a)&space;&plus;&space;\frac{UD\pi}{\rho_a&space;c_a&space;A\epsilon}(T_{inf}&space;-T_a)" title="https://latex.codecogs.com/svg.image?\large \frac{\partial T_a}{\partial t} + \frac{G}{\rho_a \epsilon} \frac{\partial T_a}{\partial t} = \ \frac{h_v}{\rho_s c_s\epsilon} (T_s - T_a) + \frac{UD\pi}{\rho_a c_a A\epsilon}(T_{inf} -T_a)" />

</p>


The equations mentioned above are solved in MATLAB using built-in pdepe solver and the temperature gradient is obtained at each time step for both solid and fluid


