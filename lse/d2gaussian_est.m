function F = d2gaussian_est(x,xdata)
 F = x(6) + x(1) * exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) ));
 
 %% x: parameters
 % x1: total number of photons
 % x2: x axis mean
 % x3: x axis stddev
 % x4: y axis mean
 % x5: y axis stddev
 % x6: background photons/pixel
 
 
