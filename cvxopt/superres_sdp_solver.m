function f_est = superres_sdp_solver(y, delta_noise, n)

    Tc = (n-1)/2;
    cvx_solver sdpt3
    cvx_precision best
    cvx_begin sdp quiet
        variable X(n+1,n+1) hermitian;
        X >= 0;
        X(n+1,n+1) == 1;
        trace(X) == 2;
        for j = 1:n-1,
            sum(diag(X,j)) == X(n+1-j,n+1);
        end
        maximize(real(X(1:n,n+1)'*y)-norm(X(1:n,n+1))*delta_noise)
    cvx_end

    dual_op=cvx_optval;
    u = X(1:n,n+1);

    aux_u =- conv(u,flipud(conj(u)));
    aux_u(2*Tc+1)=1+aux_u(2*Tc+1);
    roots_pol = roots(flipud(aux_u));

    % Isolate roots on the unit circle
    roots_detected = roots_pol(abs(1-abs(roots_pol)) < 1e-4);
    [auxsort,ind_sort]=sort(real(roots_detected));
    roots_detected = roots_detected(ind_sort);

    % Roots are double so take 1 out of 2 and compute argument
    f_rec_roots = angle(roots_detected(1:2:end))/2/pi;

    % Argument is between -1/2 and 1/2 so convert angles
    f_rec_roots(f_rec_roots < 0)= f_rec_roots(f_rec_roots < 0) + 1;  
    f_est=sort(f_rec_roots);
    
end