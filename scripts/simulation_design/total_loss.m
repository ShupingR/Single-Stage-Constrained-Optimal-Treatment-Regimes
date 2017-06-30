function val = total_loss (beta, targetProb, targetRatio, Xint)
    L = length(beta);
    beta1 = beta(1:L/2);
    beta2 = beta(L/2+1:L);
    beta1 = beta1/norm(beta1);
    beta2 = beta2/norm(beta2);
    val = relative_effect_ratio_loss(beta1, beta2, targetRatio, Xint) + ...
            prob_disagree_loss (beta1, beta2, targetProb, Xint);
end