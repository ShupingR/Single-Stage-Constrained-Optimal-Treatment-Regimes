% est values of a regime using analytic form
function [c,ceq, DC, DCeq] = analy_constraint_grad(theta, X, z0, z2, n, c)
  %  tic;
    theta = theta/norm(theta);
    z1 = X * theta;
    h1 =  n^(-2/6)*var(z1);
    c = mean( z0 + z2 .* (1- 2*normcdf(-1* z1/h1)) ) - c;
    ceq = [];
    % ceq = (theta' * theta - 1)^2;
    if nargout > 2
        DC= mean( 2*z2/h1 * normpdf(-1 * z1/h1) * X , 1);
        % DCeq = 4 * ( theta' * theta -1) * theta' ;
        DCeq = [];
    end
    % toc;
end