clc
clear
rng(2017,'twister');

num_parfor = 40;
parpool(num_parfor)

num_mc = 80;
num_boot = 20;
nrs = 5;
n =1000;

a1 = [0.5; -0.5; 1];
a2 = [-0.5; 0.5; -1];

sim_parm = load('sim_parm');
num_set = size(sim_parm, 1);

tic;
for set = 1:num_set
  b1 = sim_parm(set, 1:3)';
  b2 = sim_parm(set, 4:6)';

  % solver set up
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

  cList = linspace( constraint_min , constraint_max, 15 );

  % MC loop
  parfor mc = 1:num_mc
    rng( (2017+mc), 'twister' );
    [X, alpha1Hat, beta1Hat, alpha2Hat, beta2Hat, ~, ~] = gen_data(n, a1, b1, a2, b2);
    for boot = 1:num_boot
      fileID = fopen(strcat('output_setting_', num2str(set), '_mc_', num2str(num_mc), '_boot_', num2str(boot), '_nrs_', num2str(nrs), '.txt'),'a');
      if boot ==1
        fprintf(fileID,'bound, object_val, constraint_val, exitlflag, thetaSol1, thetaSol2 \r\n');
      end
      rng( (2017+boot), 'twister' );
      X_boot = datasample(X, n);
      Xint_boot = [ones(n,1) X_boot];
      z0 = Xint_boot*alpha1Hat;
      % z1 = Xint*theta;
      z2_1 = Xint_boot*beta1Hat;
      z2_2 = Xint_boot*beta2Hat;

      %% options for fmincon
      options = optimset('Algorithm','interior-point','Display','off');
      theta0 = [-0.5; 0.5]; lb = [-1; -1]; ub = [1; 1];
      %                      'PlotFcns',@optimplotfval,'Display','off' );

      %% the maximum objective function value
      my_objective = @(theta) analy_objective(theta, X_boot, z0, z2_1, n, -1);
      problem = createOptimProblem('fmincon', 'objective', my_objective,  'x0', theta0, 'lb',lb,'ub', ub,'options', options);
      ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
      [ thetaSol, object_val, exitflag ] = run(ms, problem, nrs);
      constraint_val = analy_objective(thetaSol, X_boot, z0, z2_2, n, 1);
      object_val = -1 * object_val;
      A1 = [222, object_val, nan, exitflag, thetaSol(1)/norm(thetaSol), thetaSol(2)/norm(thetaSol)];
      fprintf(fileID,'%4.4f, %4.4f, %4.4f, %d, %4.4f, %4.4f \r\n',A1);

      % the minimum constraint function value
      my_constraint = @(theta) analy_objective(theta, X_boot, z0, z2_2, n, 1);
      problem = createOptimProblem('fmincon', 'objective', my_constraint, 'x0', theta0, 'lb',lb,'ub', ub, 'options', options);
      ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
      [ thetaSol, object_val, exitflag ] = run(ms, problem, nrs);
      min_constraint = object_val;
      A2 = [-222, nan, min_constraint, exitflag, thetaSol(1)/norm(thetaSol), thetaSol(2)/norm(thetaSol)];
      fprintf(fileID,'%4.4f, %4.4f, %4.4f, %d, %4.4f, %4.4f \r\n',A2);

      %% estimate constrained optimal regime index paramters
%       objValList = nan(length(cList), 1);
%       constValList = nan(length(cList), 1);
%       exitflagList = nan(length(cList), 1);
      ic = 1;
      for c = cList
        if (c > min_constraint)
          my_objective = @(theta) analy_objective(theta, X_boot, z0, z2_1, n, -1);
          my_constraint = @(theta) analy_constraint(theta, X_boot, z0, z2_2, n, c);
          problem = createOptimProblem('fmincon', 'objective', my_objective, 'x0', theta0, 'lb',lb,'ub', ub,'nonlcon', my_constraint, 'options', options);
          ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
          [ thetaSol, object_val, exitflag ] = run(ms, problem, nrs);
          object_val = -1 * object_val;
          constraint_val = analy_objective(thetaSol, X_boot, z0, z2_2, n, 1);
          A3 = [c, object_val, constraint_val, exitflag, thetaSol(1)/norm(thetaSol), thetaSol(2)/norm(thetaSol)];
          fprintf(fileID,'%4.4f, %4.4f, %4.4f, %d, %4.4f, %4.4f \r\n',A3);
          ic = ic + 1;
        else
          % do nothing
        end
      end
      fclose(fileID);
    end
  end
end
toc;
delete(gcp('nocreate'))
% %fileID = fopen('output.txt','a');
% %fprintf(fileID,'k, kappa, fval, exitlfag, tau1, tau2, tau3, tau4, tau5, tau6 \r\n');
% parpool(npar)
%
% parfor k = 1:nk
%     mytime = cputime;
%     kappa = kappa_list(k)
%     my_objective = @(tau) objective_function_two( tau, sample);
%     my_constraint = @(tau) constraint_function_two( tau, sample, kappa );
%     problem = createOptimProblem('fmincon', 'objective', my_objective, ...
%                                              'x0', tau0,  'nonlcon', my_constraint, 'options', options);
%
%     ms = MultiStart('StartPointsToRun', 'all', 'Display','on');
% %
%     [tauSol, fval, exitflag] = run(ms, problem, ns);
%     objective_val = -1 * fval;
%     A =[ k, kappa, objective_val, exitflag,vec2mat(tauSol, length(tauSol)) ] ;
%     fileID = fopen('output_3.txt','a');
%     fprintf(fileID,'%d, %4.4f, %4.4f, %d, %4.4f, %4.4f,%4.4f, %4.4f, %4.4f, %4.4f \r\n',A);
% end
% etime =
% timeID = fopen('output_time.txt','a');
% fprintf(timeID,' time: %10.2f',etime);
% delete(gcp('nocreate'))
% %profile viewer
% %profsave
