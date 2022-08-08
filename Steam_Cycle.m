%%% This code tries to calcuate thermodynamic properties at all points of
%%% steam cycle.
%%% Point 1 represents inlet of pump
%%% Point 2 reperesents inlet of heat exchanger
%%% Point 3 represent inlet of turbine
%%% Point 4 represent inlet of condenser

%%% All property values are obtained from CoolProp from its online
%%% repository (https://ibell.pythonanywhere.com/)

%%% The steam turbine chosen is D-R RLH series (0-1.865 MW)
%%% The inlet steam pressure and temperature for turbine is upto 63 bar and 482 C
%%% The back pressure (assumed as outlet pressure) is upto 21 bar
%%% https://www.directindustry.com/prod/siemens-power-genereration/product-23116-2019858.html

    
%%% Point 3 (Turbine Inlet)
%%% The temperature and pressure are given so other properties can be found
P3 = 63e5;          %Pressure in Pascal
T3 = 482 + 273;     %Temperature in Kelvin
s3 = py.CoolProp.CoolProp.PropsSI('S','T',T3,'P',P3,'Water'); %Entropy (J/kg/K)
h3 = py.CoolProp.CoolProp.PropsSI('H','T',T3,'P',P3,'Water'); %Enthalpy (J/kg)


%%% Point 4 (Turbine Outlet)
P4 = 0.5e5;
s4_s = s3;          %Assuming Isentropic Expansion
T4_s = py.CoolProp.CoolProp.PropsSI('T','S',s4_s,'P',P4,'Water');
h4_s = py.CoolProp.CoolProp.PropsSI('H','S',s4_s,'P',P4,'Water');
eta_turbine = 0.9;  %Isentropic Efficiency
h4 = h3 - eta_turbine*(h3-h4_s);
P_req = 1.4e6;      %Power required from turbine
eta_shaft = 0.95;   %Mechanical Efficiency of generator
eta_gen = 0.95;     %Electrical efficiency of generator
m_steam =  P_req/((h3 - h4)*eta_shaft*eta_gen); %Mass Flow Rate of Steam

%%% Point 1 (Inlet of Pump)
%%% Assuming that there is no pressure loss across condenser

P1 = P4;
x = 0; %The output of condenser is a saturated liquid
T1    = py.CoolProp.CoolProp.PropsSI('T','Q',0,'P',P4,'Water');
rho_1 = py.CoolProp.CoolProp.PropsSI('D','Q',0,'P',P4,'Water'); %Density in kg/m3
v1 = 1/rho_1; %Specific Volume
h1 = py.CoolProp.CoolProp.PropsSI('H','Q',0,'P',P4,'Water');
s1 = py.CoolProp.CoolProp.PropsSI('S','Q',0,'P',P4,'Water');


%%% Point 2 (HEX Inlet)
%%% We assume a pressure drop of 10% across heat exchanger
P2 = 1.1*P3;
%%% We assume isentropic compression first with an efficiency of 90%
eta_compressor = 0.9;
s2_s = s1;
T2_s = py.CoolProp.CoolProp.PropsSI('T','P',P2,'S',s2_s,'Water');
h2_s = py.CoolProp.CoolProp.PropsSI('H','P',P2,'S',s2_s,'Water');
h2 = h1 + v1*(P2-P1)/eta_compressor;

P_fromHEX = m_steam*(h3-h2);
efficiency_cycle = 100*P_req/P_fromHEX


