
t=input(' Enter the time [YYYY,MM,DD,hh,mm,ss]  ');
ENU=input(' Enter the ENU position, n x 3 matrix ') ;
el_mask = input(' Enter the Elevation Mask (deg) : ');

for i = 1 : size(ENU,1)
    RE(i) = ENU(i,1);
    RN(i) = ENU(i,2);
    RU(i) = ENU(i,3);
    Rrel(i) = sqrt(RE(i)^2+RN(i)^2+RU(i)^2);
    Az_rad(i) = acos(RN(i)/sqrt(RE(i)^2+RN(i)^2));
    El_rad(i) = asin(RU/Rrel);
end

% 1. ECI to ECEF DCM

function DCM=ECI2ECEF_DCM(time)

time=datetime(t); % Current Time 
time_zero = datetime(2000,01,01,12,00,00); % Universal Time
jd=juliandate(time_zero);
theta_g0=siderealTime(jd); 

delta_t=seconds(time - time_zero); % Duration

w_earth = 7.27*10^(-5); % rad/s

theta_g_rad = w_earth * delta_t + theta_g0;
theta_g = theta_g_rad * 180 / pi ;

DCM=[cos(theta_g) sin(theta_g), 0; -sin(theta_g) cos(theta_g), 0; 0 0 1 ] ;
disp(' ECI to ECEF DCM : ');
disp(DCM)
end


% 2. ENU to Azimuth Angle

function az = azimuth(ENU)

Az=Az_rad * 180 / pi ; 

disp(' Azimuth Angle ( deg ) : ');
disp(Az);
end

% 3. ENU to Elevation Angle

function el = elevation(ENU, el_mask);

El = El_rad * 180 / pi ;
el_mask_deg = el_mask * 180 / pi ;

if El > el_mask_deg
    disp(' Elevation Angle (deg) : ');
    disp(El);
else
    disp(' NaN ');
end

end