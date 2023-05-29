semimajor_axis = input('Enter the semimajor axis : ');
eccentricity = input('Enter the eccentricity : ');
true_anomaly = input('Enter the true anomaly : ');

result_r = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly);
disp('rangeInPQW : ');
disp(result_r);
result_v = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly);
disp('velocityInPQW =');
disp(result_v);

function rangeInPQW = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly)
    true_anomaly_rad = true_anomaly*pi/180;

    PQW_r = semimajor_axis * (1 - eccentricity^2) / (1 + eccentricity * cos(true_anomaly_rad));

    rangeInPQW = [PQW_r * cos(true_anomaly_rad);
                  PQW_r * sin(true_anomaly_rad);
                  0];
end


function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly)
    true_anomaly_rad = true_anomaly*pi/180;

    u = 3.986004418*10^5;

    p = semimajor_axis * (1 - eccentricity^2);

    PQW_v = sqrt(u / p);

    velocityInPQW = [-PQW_v * sin(true_anomaly_rad);
                     PQW_v * (eccentricity + cos(true_anomaly_rad));
                     0];
end
