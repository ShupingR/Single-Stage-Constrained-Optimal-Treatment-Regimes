targetProbList=  [0.3, 0.5, 0.7];
targetRatioList = [ 0.5, 1, 1.5];
betaSolList = nan( length(targetProbList)*length(targetRatioList), 12);
rsquare = 0.7;
rconst =  ( 1 - rsquare^2 ) / rsquare^2;
alpha1 = [0.5; -0.5; 1];
alpha2 = [-0.5; 0.5; 1];
i = 1;
for targetProb = targetProbList
    for targetRatio = targetRatioList
        [ betaSol, object_val, exitflag ] = min_loss(targetProb, targetRatio);
        beta = betaSol(:);
        nb = length(beta);
        beta1 = beta(1:nb/2);
        beta2 = beta(nb/2 +1: end);
        var1 = rconst * ( alpha1(2:end)'*alpha1(2:end) + beta1'*beta1 );
        var2 = rconst * ( alpha1(2:end)'*alpha1(2:end) + beta2'*beta2 );
        betaSolList(i, :) = ...
        [betaSol', var1 , var2, object_val, exitflag, targetProb, targetRatio];
        i = i + 1;
    end
end

dataset({betaSolList  'b11','b12','b13','b21','b22','b23', 'var1', 'var2',...
            'object_val', 'exitflag', 'targetProb', 'targetRatio'})
dlmwrite('beta_sim_design_40', betaSolList);
