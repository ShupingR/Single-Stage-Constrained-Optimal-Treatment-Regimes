% estimate objective value of a regime using analytic form
function val = analy_objective(theta, X, z0, z2, n, sign)
  %  tic;
    theta = theta/norm(theta);
    z1 = X * theta;
   % hSquare = [ n^(-2/6)*var(z1), n^(-2/6)*var(z2) ];
    h1 = n^(-2/6)*var(z1);
    val = sign * mean( z0 + z2 .* (1- 2*normcdf(-1* z1/h1)) );
   % toc;
end