semimajor_axis = input('Enter the semimajor axis : ');
eccentricity = input('Enter the eccentricity : ');
true_anomaly = input('Enter the true anomaly : ');

disp(solveRangeInPerifocalFrame);
disp(solveVelocityInPerifocalFrame);

function rangeInPQW = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly)
    true_anomaly_rad = deg2rad(true_anomaly);

    PQW_r = semimajor_axis * (1 - eccentricity^2) / (1 + eccentricity * cos(true_anomaly_rad));

    rangeInPQW = [PQW_r * cos(true_anomaly_rad);
                  PQW_r * sin(true_anomaly_rad);
                  0];
end


function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly)
    true_anomaly_rad = deg2rad(true_anomaly);

    % Gravitational constant (μ) in km^3/s^2
    u = 3.986004418*10^5;

    % Calculate the range (position vector) in the perifocal frame
    p = semimajor_axis * (1 - eccentricity^2);

    % Calculate the magnitude of the velocity in the perifocal frame
    PQW_v = sqrt(u / p);

    % Construct the velocity vector in the perifocal frame
    velocityInPQW = [-PQW_v * sin(true_anomaly_rad);
                     PQW_v * (eccentricity + cos(true_anomaly_rad));
                     0];
end
