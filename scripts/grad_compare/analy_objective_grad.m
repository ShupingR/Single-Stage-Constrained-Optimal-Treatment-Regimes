% estimate objective value of a regime using analytic form
function [f, g] = analy_objective_grad(theta, X, z0, z2, n)
  %  tic;
   theta = theta/norm(theta);
    z1 = X * theta;
   % hSquare = [ n^(-2/6)*var(z1), n^(-2/6)*var(z2) ];
    h1 = n^(-2/6)*var(z1);
    f = -1 * mean( z0 + z2 .* (1- 2*normcdf(-1* z1/h1)) );
    if nargout > 1 % gradient required
        g = -1 * mean( 2*z2/h1 * normpdf(-1 * z1/h1) * X , 1);
    end
   % toc;
end