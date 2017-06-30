function loss = prob_disagree_loss (beta1, beta2, targetProb, Xint)
  % function to calculate the square of the difference between
  % targetprob and numerical prob , probability of disagreement with two strategies
  prob = mean( (Xint*beta1) .* (Xint*beta2) > 0 ) ;
  loss = (prob-targetProb)^2;
end