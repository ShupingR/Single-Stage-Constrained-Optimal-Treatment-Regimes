function loss = relative_effect_ratio_loss (beta1, beta2, targetRatio, Xint)
  % relative effect ratio optimal regimes disagree 
  disagreeIndicator = (Xint*beta1) .* (Xint*beta2) > 0;
  part1 = mean( abs(Xint * beta1) .* disagreeIndicator );
  part2 = mean( abs(Xint * beta2) .* disagreeIndicator );
  relativePart1 = part1 / mean(abs(Xint * beta1));
  relativePart2 = part2 / mean(abs(Xint * beta2));
  ratio = relativePart1/relativePart2;
  loss = (ratio - targetRatio)^2;
end