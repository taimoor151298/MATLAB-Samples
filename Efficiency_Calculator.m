close all
clear all
power = 5.4e6;
dt_sec = 60;
initial_temp = [9.8 + 273;9.8+273];
size_battery = 30*3.6e9;
[SOC,final_temp,~, P_pump] = Charger(power,initial_temp);
temp_matrix = final_temp;

%for i = 1:9
%    %disp('hello ji')
%    %disp(temp_matrix(i,:))
%    [SOC_matrix, out_temp] = Charger(power, temp_matrix(:,i));
%    SOC = [SOC, SOC_matrix(2:end)+ SOC(end)];
%    temp_matrix = [temp_matrix, out_temp];
%    
%end
E_stored = SOC*size_battery;
efficiency_withoutloss = diff(E_stored,1)/(power*dt_sec);
efficiency_withloss = (diff(E_stored,1) - P_pump*dt_sec)/(power*dt_sec);
%figure()
%hold on
%plot(SOC(1:end-1),efficiency_withoutloss)

plot(SOC(1:end-1), efficiency_withloss, 'LineWidth', 2)
xlabel('SOC')
ylabel('Charging Efficiency \eta_{ch}')
%hold off
%legend('Without Loss', 'With Loss')
figure()
plot(SOC)
%to_csv = [SOC(1:end-1)' , efficiency'];
%writematrix( to_csv, 'efficiency.csv');
%SOC = linspace(0,SOC(end),540);


[curve1, goodness1, output1] = fit(SOC(106:253)',efficiency_withloss(106:253)','poly1');
[curve2, goodness2, output2] = fit(SOC(254:401)',efficiency_withloss(254:401)','poly1');
[curve3, goodness3, output3] = fit(SOC(1:401)',efficiency_withloss(1:401)','poly1');
[curve4, goodness4, output4] = fit(SOC(1:401)',efficiency_withloss(1:401)','poly2');
goodness1
curve1
goodness2
curve2
%p1 = plot(curve1,SOC(106:253),efficiency_withloss(106:253), 'b-');
hold on

%p2 = plot(curve2,SOC(254:401),efficiency_withloss(254:401), 'b-');
temp_a = curve1(SOC(106:253));
temp_b = curve2(SOC(254:401));
temp_c = [temp_a; temp_b];
figure(6)
p2 = plot(SOC(106:401), efficiency_withloss(106:401), 'LineWidth', 2)
hold on
p1 = plot(SOC(106:401), temp_c, 'LineWidth',2)


%p2 = plot(curve3,SOC(1:401),efficiency_withloss(1:401));
set(findall(gca, 'Type', 'Line'),'LineWidth',2)
xlabel('SOC')
ylabel('Charge Efficiency \eta_{ch}')
legend('Actual Efficiency','Linearized Efficiency')
%plot(SOC(106:end-1), efficiency_withloss(106:end) )
grid
hold off

average_1 = mean(efficiency_withloss(106:401));
%efficiency = linspace(0.4342,0.6531,6000);
%Pbat = linspace(0,5.4,6000);
%surf(efficiency, Pbat, efficiency'*Pbat)

