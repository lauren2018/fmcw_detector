function y = outlier_test(distance, distance_est, tolerance)
    if (abs(distance-distance_est) > tolerance)
        y = 1;
    else
        y = 0;
    end
end