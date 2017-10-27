function F = pnoise_model(x,xdata)

    F = x(3) + x(2) ./ ( 1 + (pi * (x(2)) * ( xdata - x(1)*(xdata(2)-xdata(1)) ) ).^2 ) ;

 %F = x(6) + x(1) * exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) ));
 
 %% x: parameters
 % x1: center location
 % x2: tauc
 % x3: shot noise floor
 
end