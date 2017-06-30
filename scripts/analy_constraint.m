% est values of a regime using analytic form
function [c,ceq] = analy_constraint(theta, X, z0, z2, n, c)
  %  tic;
    theta = theta/norm(theta);
    z1 = X * theta;
    hSquare = [ n^(-2/6)*var(z1), n^(-2/6)*var(z2) ];
    c = mean( z0 + z2 .* (1- 2*normcdf(-1* z1/hSquare(1))) ) - c;
    ceq = [];
    % ceq = (theta' * theta - 1)^2;
    % toc;
end