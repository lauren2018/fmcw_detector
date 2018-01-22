function F = pnoise_model_data(x,xdata)

    tauc = 0.4e-7;
%     F =10^(-79/10)+ x(2) ./ ( 1 + (pi * (0.5e-6) * ( xdata - x(1)*(xdata(2)-xdata(1)) ) ).^2 ) + x(2) ./ ( 1 + (pi * (0.5e-6) * ( xdata + x(1)*(xdata(2)-xdata(1)) ) ).^2 ) ;
    F = (2*10^-9.12+(x(2) ./ ( 1 + (pi * (tauc) * ( xdata -xdata(1) - x(1)*(xdata(2)-xdata(1)) ) ).^2 ) + x(2) ./ ( 1 + (pi * (tauc) * ( xdata+xdata(1) + x(1)*(xdata(2)-xdata(1)) ) ).^2 )));
%     F = 10*(log10(x(2) ./ ( 1 + (pi * (tauc) * ( xdata -xdata(1) - x(1)*(xdata(2)-xdata(1)) ) ).^2 ))) ;
%     F = 10*log10(x(2) ./ ( 1 + (pi * (tauc) * ( xdata - x(1)*solu(xdata(2)-xdata(1)) ) ).^2 ) + x(2) ./ ( 1 + (pi * (tauc) * ( xdata+ x(1)*(xdata(2)-xdata(1)) ) ).^2 )) ;

    %     F =0.514e-9+ x(2) ./ ( 1 + (pi * (x(3)) * ( xdata - x(1)*(xdata(2)-xdata(1)) ) ).^2 ) + x(2) ./ ( 1 + (pi * (x(3)) * ( xdata + x(1)*(xdata(2)-xdata(1)) ) ).^2 ) ;

 %F = x(6) + x(1) * exp(-((xdata(:,:,1)-x(2)).^2/(2*x(3)^2) + (xdata(:,:,2)-x(4)).^2/(2*x(5)^2) ));
 
 %% x: parameters
 % x1: center location
 % x2: tauc
 % x3: shot noise floo6
 
end