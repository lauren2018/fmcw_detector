function F = pnoise_model(fnoise,snoise, x,xdata)

    omegac = 2*pi*fnoise;
    F = 10*log10(2*(snoise + x(2) ./ ( omegac^2 + (2*pi * ( xdata - x(1)*(xdata(2)-xdata(1)) ) ).^2 ) + x(2) ./ ( omegac^2 + (2*pi * ( xdata + x(1)*(xdata(2)-xdata(1)) ) ).^2 ))) ;

 %F = x(6) + x(1) * exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) ));
 
 %% x: parameters
 % x1: center location
 % x2: tauc
 % x3: shot noise floor
 
end