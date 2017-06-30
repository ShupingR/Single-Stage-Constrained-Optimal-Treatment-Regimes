function val = kde_val(theta, X, z0, z2, n)
    % X: 2 dimension , it does not include intercept
    % n: dataset sample size
    % nkde: bootstrap sample size
    % column vector for theta, indexing a regime
    tic;
    nkde = 10000;
    thetaN = theta/norm(theta);
    z1 = X * thetaN;
    z = [z1, z2];
    % smooth bootstrap generates samples equals to KDE 
    % silverman rule of thumbe for bandzwidth, noise std
    % product kde assume independency in two dim
    hSquare = [ n^(-2/6)*var(z1), n^(-2/6)*var(z2) ];
    rng(1000,'twister');
    zRs = datasample(z, nkde, 'replace', true) + mvnrnd( [0 0], hSquare, nkde);    
    z1Rs = zRs(:, 1);
    z2Rs = zRs(:, 2);
    val =  mean(z0) + mean(sign(z1Rs) .* z2Rs);
    toc;
end
