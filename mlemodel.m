function E = mlemodel (xdata,params)
% The expected counts per pixel. (SOM Eq. 12)
    E = params(3) + (params(2) ./ ( 1 + (pi * (params(2)) * ( xdata - params(1)*(xdata(2)-xdata(1)) ) ).^2 ) ) ;
%E = params (1) * twoDGauss (x,y,params (2),params (3),params (4),params(5))+params (6);
end