function  [ betaSol, object_val, exitflag ] = min_loss(targetProb, targetRatio)
    n = 100000;
    muX = [0, 0];
    sigmaX = [1,  0 ; 0, 1];
    X = mvnrnd(muX, sigmaX, n);
    Xint = [ones(n,1) X];
    beta0 = [0.5; -0.5; 0.5; 0.5; -0.5; 0.5];   % Starting point
    lb = -5*ones(1,6);
    ub = 5*ones(1,6);
    options = optimset( 'Algorithm','interior-point', ...
                               'PlotFcns',@optimplotfval,'Display','iter' ); 
    my_objective = @(beta)total_loss (beta, targetProb, targetRatio, Xint);
    % [ betaSol, object_val, exitflag ] = simulannealbnd(my_objective, beta0,lb,ub);
    problem = createOptimProblem('fmincon', 'objective', my_objective, ...
                                                'x0', beta0, 'lb',lb,'ub', ub, 'options', options);
    ms = MultiStart('StartPointsToRun', 'all', 'Display','on');
    [ betaSol, object_val, exitflag ] = run(ms, problem, 40);   
end