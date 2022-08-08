function [E_stored_fluid, E_lost_time, Ts_output, Tf_output, P_pump] = Discharger(discharge_power, Ts_fromprevious, Tf_fromprevious)

global pa ca ka H A G U T_inf ps cs ks e h D  xmesh dx Tf_input Ts_input ks_eff

Tf_input = Tf_fromprevious;
Ts_input = Ts_fromprevious;

T = 273.15+ (750+9.8)/2; 
%https://www.siemensgamesa.com/en-int/-/media/siemensgamesa/downloads/en
%/products-and-services/hybrid-power-and-storage/
%etes/siemens-gamesa-etes-ad-teaser-industrial-decarbonization.pdf

%https://en.climate-data.org/europe/germany/hamburg/hamburg-69/
T_inf = 9.8+273;
%T_in = 750+273;

%For solid at T
%ps = 2083;              
%cs = 827.69+(0.3339*T);
%ks =(807/(350 +T - 273.15) + 0.64);

ps = 2680;
cs = 1068;
ks = 2.5;

pa = PropsSI('D','T',T,'P',101325,'Air');
ca = PropsSI('C','T',T,'P',101325,'Air');
ka = PropsSI('L','T',T,'P',101325,'Air'); 
ua  = PropsSI('V','T',T,'P',101325,'Air');
e = 0.4;

ks_eff = 1/(e/ka + (1-e)/ks);

%For air at T


dp = 0.02;
rho_c_eff = e*pa*ca+ps*cs*(1-e);
eff0 =  1;
E_des = 30*3.6e9;
DT_TES = (750-9.8);

V = (E_des)/(rho_c_eff*DT_TES*eff0) % [m3]
asp_ratio = 1.25; 
D = ((4*V)/(pi*asp_ratio))^(1/3);%[m]
H= asp_ratio*D ;%[m];;
A=(pi*D^2)/4 %[m2]=pi*D^2/4 storage cross section area


% Calculations done manually
%stcycle_eff = 0.93; 

%power_HEXout = stcycle_eff*power_fromturbine;
%Typical value from internet
%HEX_eff = 0.90;
%power_HEXin = HEX_eff*power_HEXout;



m_flow =  discharge_power/(ca*(Tf_input(1) - T_inf));
G = m_flow/A;

alpha_p = (700/(6*(1-e)))*(G^0.76)*(dp^0.24);
h = 700*(G/dp)^0.76;
Bi = h*(2*dp/3)/(ks)
inv = 1/h + dp/(ks*10);
h_eff = 1/inv
h = h_eff


Re = G*dp/ua
Prin = (ca*ua)/ka;
%Tav = (T+T_inf)/2
Tav = T;
Uin = (ka/dp)*((2.58*Re^(1/3)*Prin^(1/3))+(0.094*Re^0.8*Prin^0.4));
Rhoa_out = (6.75*10^-18*(Tav^6))-(2.429*10^-14*(Tav^5))+(3.561*10^-11*(Tav^4))-(2.799*10^-8*(Tav^3))+(1.343*10^-5*(Tav^2))-(0.004509*Tav)+1.274;
Visc_out = (8.118*10^-15*Tav^3)-(2.243*10^-11*Tav^2)+(4.76*10^-8*Tav)+(1.743*10^-5);
Re_out = (D*5*Rhoa_out)/Visc_out;
%  Prandtl number %%%%
Ca_out = (2.42*10^-10*(Tav^4))-(7.131*10^-7*(Tav^3))+(0.0006581*(Tav^2))-(0.008615*Tav)+1006;
Ka_out = (9.381*10^-12*(Tav^3))-(2.592*10^-8*(Tav^2))+(7.298*10^-5*(Tav))+0.02477;
Prout = (Ca_out*Visc_out)/Ka_out;
Uout = (Ka_out/D)*(0.037*Re_out^(4/5)-871)*Prout^(1/3)
rin = D/2; ins1 = 0.05; steel = 0.02; ins2 = 0.25; 
R2 = rin+ins1; R3 = R2+steel; R4 = R3+ins2;
Kins1 = 1.25; Kins2 = 0.14; Ksteel = 20;
U_ext = 1/((1/Uin)+((rin/R4)*1/Uout)+(rin*((Kins1^-1*log(R2/rin))+(Ksteel^-1*log(R3/R2))+(Kins2^-1*log(R4/R3)))))

U = U_ext;


disp('it starts')
e
G
pa
dp
ua
A
a = 1.75*(1-e)/e^3;
b = G^2/(pa*dp);
c = (150*(1-e)^2)/e^3;
d = (G*ua)/(pa*dp^2);
delP = H * ((a*b + c*d))
P_pump = delP*A*G/pa
h;
ks;

%U of an insulated wall
U = 0.678;


m = 0;
dt = 60;  %60 sec
n_mesh = 500 + 1;    
t_max = 120;
n_step = t_max/dt+1;
xmesh = linspace(0,H,n_mesh);
dx = H/(n_mesh-1);
tstep = linspace(0,t_max,n_step);
options = odeset('RelTol',1e-1,'AbsTol',1e-2);
%sol = pdepe(m,@Charge_eq,@ini_cd,@boundary,xmesh,tstep,options);
sol = pdepe(m,@Discharge_eq,@ini_cd_discharge,@boundary_discharge,xmesh,tstep,options);

Ts = sol(:,:,2);
Tf = sol(:,:,1);
Ts_output = Ts(end,:);
Tf_output = Tf(end,:);
%size(Ts)
% h = surf(xmesh,tstep/3600,u)
% set(h,'linestyle','none');
% colorbar()
% xlabel('x')
% ylabel('t')
% zlabel('u')
%save('temp_high.mat','Ts')

%
%figure(1)
t_hr = tstep/3600;
t_min = tstep/60;
% plot(t_hr, transpose(Ts(:,1:25:end)))
% %plot(t_min, (Ts(:,1:20:end)))
% legend(sprintf('At Tank Height'))
% legend(string(round(xmesh(1:25:end),2)))
% 
% [leg,att] = legend('show');
% title(leg,'Tank Depth (m)')
% leg.Title.Visible = 'on'
% xlabel('Time (hr)')
% ylabel('Temperature')
% 
% 
% figure(2)
% plot_t = 60;
% size(Ts)
% plot(xmesh,transpose((Ts(1:plot_t:end,:))),'-o')
% 
% hold on
% plot(xmesh,transpose((Tf(1:plot_t:end,:))))
% 
% legend('Solid 1200', 'Solid 2400', 'Solid 3600', 'Solid 4800' ...
% ,'Fluid 1200', 'Fluid 2400', 'Fluid 3600', 'Fluid 4800')
% hold off
% legend(string(round(tstep(1:plot_t:end)/60,2)));
% [leg,att] = legend('show');
% title(leg,'Time (min)')
% leg.Title.Visible = 'on'
% xlabel('Tank Height')
% ylabel('Temperature')
% 

%dx = H/n_mesh;
%This is for energy contained in the fluid
E_stored_fluid = [0; m_flow*dt*ca*cumsum((Tf(1:end-1,1)-T_inf))];
E_stored_fluid = E_stored_fluid';

%E_stored_time =[0; mass/n_mesh*cs*sum(temp_diff,2)];
%E_stored_solid = cumsum(E_stored_time);
E_lost_time = [];
%E_lost_time(1) = abs(sum(ps*cs*A*(1-e)* (Ts(2,:) - Ts(1,:))*dx));
for i = 1: length(tstep)
    %disp('run')
    E_lost_time(i) = abs(sum(ps*cs*A*(1-e)* (Ts(1,:) - Ts(i,:))*dx));
end
%size(E_lost_time)
%size(t_hr)
%E_stored_solid = ps*A*cs*(1-e)*trapz(xmesh,Ts(
%figure(3)
%yyaxis left
E_stored_MWh = E_stored_fluid/3.6e9;
 
%E_stored_solidMWh = E_stored_solid/3.6e9;
E_lost_solidMWh1 = E_lost_time/3.6e9;

a = 1.75*(1-e)/e^3;
b = G^2/(pa*dp);
c = (150*(1-e)^2)/e^3;
d = (G*ua)/(pa*dp^2);
delP = H * ((a*b + c*d));
P_pump = delP*A*G/pa;


%plot(t_hr, E_stored_MWh, t_hr, E_stored_solidMWh,t_hr, E_stored_solidMWh1)
% plot(t_hr, E_stored_MWh,t_hr, E_lost_solidMWh1);
% hold on
SOC = 1 - E_lost_time/E_des;

%yyaxis right
%plot(t_hr,SOC)

%legend('Fluid', 'Solid', 'SOC')
%legend(string(round(tstep(1:144:end)/3600,2)))
%[leg,att] = legend('show');
%title(leg,'Time (hr)')
%leg.Title.Visible = 'on'
%xlabel('Time')
%hold off
end

