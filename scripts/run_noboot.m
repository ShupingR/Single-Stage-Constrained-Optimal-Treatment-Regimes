clc
clear
rng(2017,'twister');

num_parfor = 44;
parpool(num_parfor)

num_mc = 200;
% num_boot
nrs = 5;

n =1000;
a1 = [0.5; -0.5; 1];
a2 = [-0.5; 0.5; -1];

sim_parm = load('sim_parm');
num_set = size(sim_parm, 1);

for setting = 1:num_set
  tic;
  % read in simulation setting
  b1 = sim_parm(setting, 1:3)';
  b2 = sim_parm(setting, 4:6)';

  % solver setting up
  options = optimset('Algorithm','interior-point','Display','off');
  theta0 = [-0.5; 0.5]; lb = [-1; -1]; ub = [1; 1];

  % find the cMin and cMax for all MC iterations
  [X, alpha1Hat, beta1Hat, alpha2Hat, beta2Hat, cMin, cMax] = gen_data(50000, a1, b1, a2, b2);
  Xint = [ones(50000,1) X]; z0 = Xint*alpha1Hat;  z2_1 = Xint*beta1Hat; z2_2 = Xint*beta2Hat; % z1 = Xint*theta;
  my_constraint_min = @(theta) analy_objective(theta, X, z0, z2_2, n, 1);
  problem = createOptimProblem('fmincon', 'objective', my_constraint_min, 'x0', theta0, 'lb',lb,'ub', ub, 'options', options);
  ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
  [ ~, constraint_min, ~ ] = run(ms, problem, nrs);

  my_constraint_max = @(theta) analy_objective(theta, X, z0, z2_2, n, -1);
  problem = createOptimProblem('fmincon', 'objective', my_constraint_max, 'x0', theta0, 'lb',lb,'ub', ub, 'options', options);
  ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
  [ ~, constraint_max, ~ ] = run(ms, problem, nrs);
  constraint_max = -1*constraint_max;

  cList = linspace( constraint_min , constraint_max, 20 );

  % MC loop
  parfor mc = 1:num_mc
    % open fileID for this setting
      % write out file header
    fileID = fopen(strcat('output_setting_', num2str(setting)', '_mc_', num2str(mc), '.txt'), 'a');
    fprintf(fileID,'bound, object_val, constraint_val, exitlflag, thetaSol1, thetaSol2 \r\n');
    rng( (2017+mc), 'twister' );

    % generate MC data
    [X, alpha1Hat, beta1Hat, alpha2Hat, beta2Hat, ~, ~] = gen_data(n, a1, b1, a2, b2);
    Xint = [ones(n,1) X];
    z0 = Xint*alpha1Hat;
    % z1 = Xint*theta;
    z2_1 = Xint*beta1Hat;
    z2_2 = Xint*beta2Hat;

    %% options for fmincon
    options = optimset('Algorithm','interior-point','Display','off');
    theta0 = [-0.5; 0.5]; lb = [-1; -1]; ub = [1; 1];
    %                      'PlotFcns',@optimplotfval,'Display','off' );

    %% the maximum objective function value
    my_objective = @(theta) analy_objective(theta, X, z0, z2_1, n, -1);
    problem = createOptimProblem('fmincon', 'objective', my_objective,  'x0', theta0, 'lb',lb,'ub', ub,'options', options);
    ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
    [ thetaSol, object_val, exitflag ] = run(ms, problem, nrs);
    constraint_val = analy_objective(thetaSol, X, z0, z2_2, n, 1);
    object_val = -1 * object_val;
    A1 = [222, object_val, nan, exitflag, thetaSol(1)/norm(thetaSol), thetaSol(2)/norm(thetaSol)];
    fprintf(fileID,'%4.4f, %4.4f, %4.4f, %d, %4.4f, %4.4f \r\n',A1);

    % the minimum constraint function value
    my_constraint = @(theta) analy_objective(theta, X, z0, z2_2, n, 1);
    problem = createOptimProblem('fmincon', 'objective', my_constraint, 'x0', theta0, 'lb',lb,'ub', ub, 'options', options);
    ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
    [ thetaSol, object_val, exitflag ] = run(ms, problem, nrs);
    min_constraint = object_val;
    A2 = [-222, nan, min_constraint, exitflag, thetaSol(1)/norm(thetaSol), thetaSol(2)/norm(thetaSol)];
    fprintf(fileID,'%4.4f, %4.4f, %4.4f, %d, %4.4f, %4.4f \r\n',A2);

    %% estimate constrained optimal regime index paramters
    ic = 1;
    for c = cList
      if (c > min_constraint)
        my_objective = @(theta) analy_objective(theta, X, z0, z2_1, n, -1);
        my_constraint = @(theta) analy_constraint(theta, X, z0, z2_2, n, c);
        problem = createOptimProblem('fmincon', 'objective', my_objective, 'x0', theta0, 'lb',lb,'ub', ub,'nonlcon', my_constraint, 'options', options);
        ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
        [ thetaSol, object_val, exitflag ] = run(ms, problem, nrs);
        object_val = -1 * object_val;
        constraint_val = analy_objective(thetaSol, X, z0, z2_2, n, 1);
        A3 = [c, object_val, constraint_val, exitflag, thetaSol(1)/norm(thetaSol), thetaSol(2)/norm(thetaSol)];
        fprintf(fileID,'%4.4f, %4.4f, %4.4f, %d, %4.4f, %4.4f \r\n',A3);
        ic = ic + 1;
      else
        A3 = nan(1,6);
        fprintf(fileID,'%4.4f, %4.4f, %4.4f, %d, %4.4f, %4.4f \r\n',A3);
      end
    end
    fclose(fileID);
  end
  time = toc;
  dlmwrite('time', time, '-append');
end

delete(gcp('nocreate'))
