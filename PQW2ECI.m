function y=PQW2ECI(arg_prg,inc_angle,RAAN)
arg_prg=input('Enter the argument of perigee : ');
inc_angle=input('Enter the Inclination angle : ');
RAAN=input('Enter the RAAN : ');

R1=[cos(RAAN) sin(RAAN) 0;-sin(RAAN) cos(RAAN) 0; 0 0 1];
R2=[1 0 0; 0 cos(inc_angle) sin(inc_angle); 0 -sin(inc_angle) cos(inc_angle)];
R3=[cos(arg_prg) sin(arg_prg) 0; -sin(arg_prg) cos(arg_prg) 0; 0 0 1];
y=R1*R2*R3;