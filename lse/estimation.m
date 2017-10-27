clc; clear; close all;

% N = 200*linspace(1,10,10);
% N = [1000 2000];
N = 1000;
NA = 1.4;
nm = 1.5;
z0 = 0e-9;
lambda = 600e-9;
pixelsize = 100e-9;
numofpixels = 10;
backgroundphotons = 100;
rep = 100;

pixelindex = linspace(0, (numofpixels/2)*pixelsize, numofpixels/2+1);

qraw = @(x, y) uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2 + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2;

qarray = zeros(1, (numofpixels/2)^2);
psfmodel = zeros(numofpixels/2, numofpixels/2);
coordinatesx = zeros(numofpixels/2, numofpixels/2);
coordinatesy = zeros(numofpixels/2, numofpixels/2);
coordinates = zeros(numofpixels, numofpixels,2);

iter = ((numofpixels^2/4));

for i=1:iter

    k = mod(i, numofpixels/2);
    if (k==0) 
        k=numofpixels/2;
    end
    j = ((i-k)/(numofpixels/2))+1;

    qarray(i) = integral2(qraw,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    
end

Az0 = sum(qarray)*4;

muarray = qarray / Az0;

for i=1:numofpixels/2
    for j=1:numofpixels/2
%         psfmodel(i,j) = muarray((i-1)*numofpixels/2 + j);
        coordinatesx(i,j) = pixelsize * j - pixelsize * 0.5;
        coordinatesy(i,j) = pixelsize * i - pixelsize * 0.5;
    end
end

psfmodel = [rot90(psfmodel,2) flipud(psfmodel) ; fliplr(psfmodel) psfmodel];
coordinatesx = [-1*rot90(coordinatesx,2) flipud(coordinatesx) ; -1*fliplr(coordinatesx) coordinatesx];
coordinatesy = [rot90(coordinatesy,2) flipud(coordinatesy) ; -1*fliplr(coordinatesy) -1*coordinatesy];
coordinates(:,:,1) = coordinatesx*1e9;
coordinates(:,:,2) = coordinatesy*1e9;

background = backgroundphotons * ones(size(psfmodel));

expected_size = 0.61e9*lambda/NA/2;

middlepatch = integral2(qraw,-1*pixelsize/2,1*pixelsize/2,-1*pixelsize/2,1*pixelsize/2)/Az0;

estimated_size_ls = zeros(length(N), rep);
estimated_size_mle = zeros(length(N), rep);
error_ls = zeros(1,length(N));
error_mle = zeros(1,length(N));

% OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
% UBOUND = [2*middlepatch, 1e9*pixelindex(end), 2*expected_size, 1e9*pixelindex(end), 2*expected_size, 2*backgroundphotons];
% LBOUND = [0, -1e9*pixelindex(end), 0, -1e9*pixelindex(end), 0, 0];

for h=1:length(N)
%     expectedsignal = zeros(numofpixels, numofpixels, rep);
    params = [N(h) * middlepatch, 0, expected_size, 0, expected_size, backgroundphotons];
    for i=1:rep
        expectedsignal = poissrnd(N(h)*psfmodel) + poissrnd(background);
        solparz_pixel_ls = lsqcurvefit(@d2gaussian_est,params,coordinates,expectedsignal);
        estimated_size_ls(h,i) = 2*NA*sqrt(solparz_pixel_ls(3)^2+solparz_pixel_ls(5)^2)/0.61;
        solparz_pixel_mle = MLEwG(coordinates, expectedsignal, params);
        estimated_size_mle(h,i) = 2*NA*sqrt(solparz_pixel_mle(3)^2+solparz_pixel_mle(5)^2)/0.61;
    end
    error_ls(h) = std(estimated_size_ls(h,:));
    error_mle(h) = std(estimated_size_mle(h,:));
end




