function solution = pnoise_fit(perin,freq,params)
    
    datafun = @(params) sum((perin - pnoise_model(params, freq)).^2);
    options = optimset ('MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-9);
%     options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
        
% 'StepTolerance',1e-10,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,...
%         'MaxFunctionEvaluations', 400,'MaxIterations',400);    
%     lb = [];
%     ub = [];
    solution = fminsearchbnd(datafun, params, [0 1e-12 0], [length(freq) Inf Inf], options);
%     solution = lsqcurvefit(@pnoise_model, params, freq, perin, [0 0 0], [], options);
end
