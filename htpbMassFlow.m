function htpbMassFlow(portDiameter)
clc
constants;
fuelLength = 36;
fuelLength = convlength(fuelLength, 'in', 'm');
portDiameter = convlength(portDiameter, 'in', 'm');
%burn rate calulation
Rr = convlength(Rr, 'in', 'm');

%burn time calulation
Hu = convlength(Hu, 'in', 'm');

Hd = convlength(Hd, 'in', 'm');

Hb = (Hd - portDiameter - Hu) / 2;    %HTPB burned
burnTime = Hb / Rr;     %time it takes to burn the given amount of HTPB

midPortDiam = ((Hd-Hu)-portDiameter)/2  + portDiameter;

a=[1:burnTime
    1:burnTime
    1:burnTime];

for i=1:burnTime
    Rd = 2 * ((portDiameter/2) + Rr);
    Vo = pi*((portDiameter/2)^2)*fuelLength;
    Vi = pi*((Rd/2)^2)*fuelLength;
    Vr = Vi - Vo;                           %volume of HTPB burend in one second in m^3
    a(1, i) = convlength(portDiameter, 'm', 'in');
    a(2, i) = Vr * Dh;                      %weight of HTPB burned in one second in kg
    
    portDiameter = portDiameter + (Rr * 2);
end

figure(1)
plot(a(1,:),a(2,:));
xlabel('Port Diameter (inches)')
ylabel('HTPB Mass Flow Rate (kg/s)')
title('HTPB Mass Flow Rate Fluctuation')

Rd = 2 * ((midPortDiam/2) + Rr);
Vo = pi*((midPortDiam/2)^2)*fuelLength;
Vi = pi*((Rd/2)^2)*fuelLength;
Vr = Vi - Vo;                           %volume of HTPB burend in one second in m^3
midBurn = Vr * Dh;
nitrousFlow = midBurn * OF;

for i=1:burnTime
    a(3, i) = nitrousFlow / a(2, i);
end

for i=1:burnTime
    a(4, i) = i;
end

for i=1:burnTime
    a(5, i) = nitrousFlow + a(2,i);
end

figure(2)
plot(a(1,:),a(3,:));
title('O/F Fluctuation')
xlabel('Port Diameter (inches)')
ylabel('O/F Ratio')

figure(3)
plot(a(4,:),a(5,:));
title('Total Mass Flow Rate Fluctuation')
xlabel('Time (seconds)')
ylabel('Total Mass Flow Rate (kg/s)')
t=0:1:10;
startBurn = a(2,1);
m_dot = startBurn + nitrousFlow;
thrust_new=0;
%% Below is the data that was found using CEA. Using different OF ratios, new values for mass flow rate, specific heat constant and molecular mass were calculated
%% Follow are those values
m_dots= [0,6.70909611081257,6.73671818995132,6.76434026909007,6.79196234822883,6.81958442736758,6.84720650650633,6.87482858564509,6.90245066478384,6.93007274392259,6.95769482306134,6.98531690220009];
k_ray=[1.5546 1.5628 1.5708 1.5787 1.5865 1.5941 1.6016 1.6090 1.6162 1.6233 1.6302]; 
Mm_val = [26.845    26.682 26.518 26.354 26.188 26.023 25.856 25.689 25.522 25.354 25.185];
for i=1:1:11
    thrust_new(i)=Thrust_output(m_dots(i),k_ray(i),Mm_val(i));%% Here is the thrust matrix
end
thrust_new=convforce(thrust_new,'N','lbf')
plot(t,thrust_new)
axis([0 12 0 4000])
title('Thrust curve')
xlabel('time')
ylabel('Thrust')

[Dt_inches,De_inches,aRatio,Ve,M,thrust,poundsForce] = Thrust(m_dot)

thrust_new=0;
m_dots= [0,6.70909611081257,6.73671818995132,6.76434026909007,6.79196234822883,6.81958442736758,6.84720650650633,6.87482858564509,6.90245066478384,6.93007274392259,6.95769482306134,6.98531690220009];
k=[1.5546 1.5628 1.5708 1.5787 1.5865 1.5941 1.6016 1.6090 1.6162 1.6233 1.6302]; 
Mm_val = [26.845    26.682 26.518 26.354 26.188 26.023 25.856 25.689 25.522 25.354 25.185];
for i=1:1:11
    thrust_new(i)=Thrust_output(m_dots(i),k(i),Mm_val(i));
end

thrust_new=convforce(thrust_new,'N','lbf')
plot(t,thrust_new)
axis([0 12 0 4000])
nitrousMass = nitrousFlow * burnTime;
Ti = convlength(Ti, 'in', 'm');
To = convlength(To, 'in', 'm');

volumeNitrous = nitrousMass / Dn;       %volume the tank will have to be in m^3
nitrousTankLength = volumeNitrous / (pi*((Ti/2)^2));      %Total nitrous tank length in m
nitrousTankLengthFt = nitrousTankLength * 3.28084

%weight of Nitrous Tank Body
Vv = pi*((To/2)^2)*nitrousTankLength;
Vz = pi*((Ti/2)^2)*nitrousTankLength;
Va = Vv - Vz;   %volume of aluminum in Nitrous Tank Body in m^3
nitrousTankMass = Va * Al;   %mass of Nitrous Tank Body in kg

%weight of Combustion Chamber Paper Core
Ci = convlength(Ci, 'in', 'm');
Vl = pi*((Ci/2)^2)*fuelLength;
Vk = pi*((Hd/2)^2)*fuelLength;
Vp = Vl - Vk;   %volume of Combustion Chamber Paper Core in m^3
paperCoreMass = Pr * Vp;   %mass of Combustion Chamber Paper Core in kg

%weight of Combustion Chamber Body
Bm = convlength(Bm, 'in', 'm');
Nm = convlength(Nm, 'in', 'm');
Pl = convlength(Pl, 'in', 'm');
Pe = convlength(Pe, 'in', 'm');
It = convlength(It, 'in', 'm');
Co = convlength(Co, 'in', 'm');

Cl = Bm + Nm + Pl + Pe + It + fuelLength;   %Combustion Chamber body length in m
Vw = pi*((Co/2)^2)*Cl;
Vq = pi*((Ci/2)^2)*Cl;
Vc = Vw - Vq;
combustionChamberMass = Vc * Al;   %weight of combustion chamber body in kg

Vf = pi*((Hd/2)^2)*fuelLength;
Vd = pi*((portDiameter/2)^2)*fuelLength;
Vh = Vf - Vd;     %Volume of HTPB grain in m^3
htpbMass = Vh * Dh;       %weight of HTPB grain in kg

%initial mass of rocket motor
rocketMass = Nc+Bh+Ni+Pq+Gr+Gi+Cr+Cc+Sp+Si+Sc+nitrousTankMass+paperCoreMass+combustionChamberMass+htpbMass+nitrousMass;
%a(6,1) = rocketMass;

forceGravity = rocketMass * -g;
F = thrust(1) + forceGravity;
e = F / rocketMass;
d = (.5 * e);       %assumed 1 second interval
b(1,1) = d;

for i=1:burnTime
    rocketMass = rocketMass - a(2,i) - nitrousFlow;
    a(6, i) = rocketMass;
end

for i=2:burnTime
forceGravity = a(6, i-1) * -g;
F = thrust(1) + forceGravity;
e = F / a(6, i-1);
d = (.5 * e) * (i^2);
b(1,i) = d + b(1,i-1);
end

figure(4)
plot(a(4,:),a(6,:));
title('Rocket Mass Fluctuation')
xlabel('Time (seconds)')
ylabel('Rocket Mass (kg)')

% figure(5)
% plot(a(4,:),b(1,:));
% title('Rocket Flight Height')
% xlabel('Time (seconds)')
% ylabel('Height (m)')


for i=1:burnTime
    a(7, i) = Thrust(a(5,i));
end

for i=1:burnTime
    a(8, i) = a(7,i) * .224808943;
end

%figure(6)
%plot(a(4,:),a(7,:));
%title('Total Thrust')
%xlabel('Time (seconds)')
%ylabel('Thrust (N)')

% figure(6)
% plot(a(4,:),a(8,:));
% title('Total Thrust')
% xlabel('Time (seconds)')
% ylabel('Thrust (lbf)')

%disp()












