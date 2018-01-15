function [thrust]=Thrust_output(m_dot,k,Mm_val)
constants;
%Pa = Pa * .4;
Tt = (2/(k+1)) * Tc;    %temp of throat
Pt = (2/(k+1))^(k/(k-1)) * Pc;      %pressure of throat
At = (m_dot/Pt) * sqrt((R * Tt)/(Mm_val * k));
r = sqrt(At/pi);
Dt_meters = r * 2;
Dt_inches = Dt_meters * 39.3701;

M = sqrt((2/(k-1))*((Pc/Pa)^((k-1)/k)-1));

Ae = (At/M)*((1+((k-1)/2)*M^2)/((k+1)/2))^((k+1)/(2*(k-1)));
re = sqrt(Ae/pi);
De_meters = re * 2;
De_inches = De_meters * 39.3701;
aRatio = Ae/At;
Te = ((1 + ((k-1)/2) * M^2)^-1) * Tc;
Ve = M * sqrt((R/Mm_val)*k*Te);



thrust = Ve * m_dot;
poundsForce = .224808943 * thrust;