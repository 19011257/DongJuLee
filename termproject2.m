

% Term Project 2

clear all
clc

load('nav.mat');

t=input(' Enter the time [YYYY,MM,DD,hh,mm,ss]  ');
t=datetime(t);
t_end = t + days(1);

el_mask = input('Enter the Elevation Mask Angle : ');

lon_r = input(' Enter the Radar Position Longitude (deg) ') ;
lat_r = input(' Enter the Radar Position Latitude (deg) ') ;


p_gps = nav.GPS.a * (1 - (nav.GPS.e)^2 );
p_qzss = nav.QZSS.a * (1 - (nav.QZSS.e)^2 );
p_bds = nav.BDS.a * (1 - (nav.BDS.e)^2 );

mu= 3.986004418 * 10^5;

n_gps = sqrt(mu/(nav.GPS.a)^3);     % n = mean motion
n_qzss = sqrt(mu/(nav.QZSS.a)^3);
n_bds = sqrt(mu/(nav.BDS.a)^3);     

lon_GPS=[];
lat_GPS=[];
lon_QZSS=[];
lat_QZSS=[];
lon_BDS=[];
lat_BDS=[];

El_GPS=[];
Az_GPS=[];
El_QZSS=[];
Az_QZSS=[];
El_BDS=[];
Az_BDS=[];

while t < t_end 
    % Calculate Eccentric Anomaly
    M_gps = n_gps * seconds(t - datetime(nav.GPS.toc)) + nav.GPS.M0 ;
    M_qzss = n_qzss * seconds(t - datetime(nav.QZSS.toc)) + nav.QZSS.M0 ;
    M_bds = n_bds * seconds(t - datetime(nav.BDS.toc)) + nav.BDS.M0 ;

    fun_g = @(E) E - (nav.GPS.e)*sin(E) - M_gps ;
    fun_q = @(E) E - (nav.QZSS.e)*sin(E) - M_qzss ;
    fun_b = @(E) E - (nav.BDS.e)*sin(E) - M_bds ;
    
    E_gps = fzero(fun_g, M_gps);
    E_qzss = fzero(fun_q, M_qzss);  
    E_bds = fzero(fun_b, M_bds);
    
    % Calculate the True anomaly
    r2d = 180/pi ; 
    nu_g = atan2(((sqrt(1-nav.GPS.e^2)*sin(E_gps))/(1-(nav.GPS.e)*cos(E_gps))),((cos(E_gps)-nav.GPS.e)/(1-(nav.GPS.e)*cos(E_gps))));
    nu_q = atan2(((sqrt(1-nav.QZSS.e^2)*sin(E_qzss))/(1-(nav.QZSS.e)*cos(E_qzss))),((cos(E_qzss)-nav.QZSS.e)/(1-(nav.QZSS.e)*cos(E_qzss))));
    nu_b = atan2(((sqrt(1-nav.BDS.e^2)*sin(E_bds))/(1-(nav.BDS.e)*cos(E_bds))),((cos(E_bds)-nav.BDS.e)/(1-(nav.BDS.e)*cos(E_bds))));
    
    r_g_s = p_gps / (1 + nav.GPS.e * cos(nu_g));
    r_q_s = p_qzss / (1 + nav.QZSS.e * cos(nu_q));  % distance r for scalar
    r_b_s = p_bds / (1 + nav.BDS.e * cos(nu_b));
    
    r_g = [r_g_s * cos(nu_g); r_g_s * sin(nu_g); 0];
    r_q = [r_q_s * cos(nu_q); r_q_s * sin(nu_q); 0];   % vector r on Perifocal frame
    r_b = [r_b_s * cos(nu_b); r_b_s * sin(nu_b); 0];
    
    % Rotation Matrix : PQW to ECI
    R1_g=[cos(nav.GPS.OMEGA) -sin(nav.GPS.OMEGA) 0;sin(nav.GPS.OMEGA) cos(nav.GPS.OMEGA) 0; 0 0 1];
    R2_g=[1 0 0; 0 cos(nav.GPS.i) -sin(nav.GPS.i); 0 sin(nav.GPS.i) cos(nav.GPS.i)];
    R3_g=[cos(nav.GPS.omega) -sin(nav.GPS.omega) 0; sin(nav.GPS.omega) cos(nav.GPS.omega) 0; 0 0 1];
    
    r_g_eci = R1_g * R2_g * R3_g * r_g ;
    
    R1_q=[cos(-nav.QZSS.OMEGA) sin(-nav.QZSS.OMEGA) 0;-sin(-nav.QZSS.OMEGA) cos(-nav.QZSS.OMEGA) 0; 0 0 1];
    R2_q=[1 0 0; 0 cos(-nav.QZSS.i) sin(-nav.QZSS.i); 0 -sin(-nav.QZSS.i) cos(-nav.QZSS.i)];
    R3_q=[cos(-nav.QZSS.omega) sin(-nav.QZSS.omega) 0; -sin(-nav.QZSS.omega) cos(-nav.QZSS.omega) 0; 0 0 1];
    
    r_q_eci = R1_q * R2_q * R3_q * r_q ;
    
    R1_b=[cos(-nav.BDS.OMEGA) sin(-nav.BDS.OMEGA) 0;-sin(-nav.BDS.OMEGA) cos(-nav.BDS.OMEGA) 0; 0 0 1];
    R2_b=[1 0 0; 0 cos(-nav.BDS.i) sin(-nav.BDS.i); 0 -sin(-nav.BDS.i) cos(-nav.BDS.i)];
    R3_b=[cos(-nav.BDS.omega) sin(-nav.BDS.omega) 0; -sin(-nav.BDS.omega) cos(-nav.BDS.omega) 0; 0 0 1];
    
    r_b_eci = R1_b * R2_b * R3_b * r_b ;
    
    % Calculate Sidereal time and 
    w_earth = 7.292116 * 10^(-5) ;  % rad/sec
    
    t_ust = datetime(2000,01,01,12,00,00); 
    jd=juliandate(t_ust);
    theta_g0=siderealTime(jd); 
    delta_t=seconds(t - t_ust); 
    theta_g = mod((w_earth * delta_t)*r2d + theta_g0,360);
    
    % DCM : ECI to ECEF frame
    C_ie=[cosd(theta_g) sind(theta_g), 0; -sind(theta_g) cosd(theta_g), 0; 0 0 1 ] ;
    
    r_g_ecef = C_ie * r_g_eci ;
    r_q_ecef = C_ie * r_q_eci ;
    r_b_ecef = C_ie * r_b_eci ;
    
    wgs84 = wgs84Ellipsoid('kilometer');
    [lat_g,lon_g,h_g] = ecef2geodetic(wgs84,r_g_ecef(1),r_g_ecef(2),r_g_ecef(3),"degrees");

    [lat_q,lon_q,h_q] = ecef2geodetic(wgs84,r_q_ecef(1),r_q_ecef(2),r_q_ecef(3),"degrees");

    [lat_b,lon_b,h_b] = ecef2geodetic(wgs84,r_b_ecef(1),r_b_ecef(2),r_b_ecef(3),"degrees");


    % ECEF frame to ENU frame
    gst=theta_g*pi/180 ;
    R_earth = 6371000; 
    R_ecef = [R_earth*cos(lon_r)*cos(lat_r); R_earth*cos(lat_r)*sin(lon_r); R_earth*sin(lat_r)];
    C_en=[-sin(lon_r+gst), cos(lon_r+gst), 0 ; -cos(lon_r+gst)*sin(lat_r), -sin(lon_r+gst)*sin(lat_r), cos(lat_r) ; cos(lon_r+gst)*cos(lat_r), cos(lat_r)*sin(lon_r+gst), sin(lat_r)];
    
    r_g_enu = C_en * (r_g_ecef - R_ecef);
    r_q_enu = C_en * (r_q_ecef - R_ecef);
    r_b_enu = C_en * (r_b_ecef - R_ecef);

    El_g=r2d*asin(r_g_enu(3)/sqrt(r_g_enu(1)^2+r_g_enu(2)^2+r_g_enu(3)^2));
    Az_g=r2d*acos(r_g_enu(2)/sqrt(r_g_enu(1)^2+r_g_enu(2)^2));
    
    if El_g <= el_mask
        El_g = 0;
    end
    
    if Az_g <0
        Az_g = mod(Az_g,360);
    end
    if r_g_enu(1) < 0
        Az_g = 360 - Az_g ;
    else
        Az_g = Az_g ;
    end


    El_q=r2d*asin(r_q_enu(3)/sqrt(r_q_enu(1)^2+r_q_enu(2)^2+r_q_enu(3)^2));
    Az_q=r2d*acos(r_q_enu(2)/sqrt(r_q_enu(1)^2+r_q_enu(2)^2));
    if El_q <= el_mask
        El_q = 0;
    end
   
    if Az_q <0
        Az_q = mod(Az_q,360);
    end
    if r_q_enu(1) < 0
        Az_q = 360 - Az_q ;
    else
        Az_q = Az_q ;
    end

    El_b=r2d*asin(r_b_enu(3)/sqrt(r_b_enu(1)^2+r_b_enu(2)^2+r_b_enu(3)^2));
    Az_b=r2d*acos(r_b_enu(2)/sqrt(r_b_enu(1)^2+r_b_enu(2)^2));
    if El_b <= el_mask
        El_b = 0;
    end
    
    if Az_b <0
        Az_b = mod(Az_b,360);
    end
    if r_b_enu(1) < 0
        Az_b = 360 - Az_b ;
    else
        Az_b = Az_b ;
    end

    lon_GPS=[lon_GPS ; lon_g];
    lat_GPS=[lat_GPS ; lat_g];
    lon_QZSS=[lon_QZSS; lon_q];
    lat_QZSS=[lat_QZSS; lat_q];
    lon_BDS=[lon_BDS; lon_b];
    lat_BDS=[lat_BDS; lat_b];

    El_GPS=[El_GPS; El_g];
    Az_GPS=[Az_GPS; Az_g];
    El_QZSS=[El_QZSS; El_q];
    Az_QZSS=[Az_QZSS; Az_q];
    El_BDS=[El_BDS; El_b];
    Az_BDS=[Az_BDS; Az_b];

    t = t + minutes(1);

end

%legend('GPS','QZSS','BDS')