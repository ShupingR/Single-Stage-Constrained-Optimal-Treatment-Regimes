function [X, a1Hat, b1Hat, a2Hat, b2Hat, cMin, cMax] = gen_data(n, alpha1, beta1, alpha2, beta2)
    % Generate trainset
    muX = [0, 0];
    sigmaX = [1,  0 ; 0, 1];
    X = mvnrnd(muX, sigmaX, n);
    Xint = [ones(n,1) X];
    % SIGMA = [1 0.8;0.8 1];
    % X = mvtrnd(SIGMA,3,n);
    % Generate error term in linear regressions for Y1 and Y2
    muE1 = 0;
    sigmaE1 = 1.2;
    e1 = normrnd(muE1, sigmaE1, n, 1);

    muE2 = 0;
    sigmaE2 = 1.5;
    e2 = normrnd(muE2, sigmaE2, n, 1);

    % Generate A treatment assignment (-1, 1)
    a = randi(2, n, 1);
    a(a==2) = -1;
    A = repmat(a, 1, size(Xint, 2));

    % Generate Y1 and Y2
    Y1 = Xint*alpha1 + a.*(Xint*beta1) + e1;
    Y2 = Xint*alpha2 + a.*(Xint*beta2) + e2;
    Xall = [Xint, A.*Xint];
    cMin = min(Y2);
    cMax =max(Y2);
    % least square estimate
    mdlY1 = fitlm(Xall,Y1,'Intercept',false);
    mdlY2 = fitlm(Xall,Y2,'Intercept',false);

    parm1Hat = mdlY1.Coefficients.Estimate;
    a1Hat = parm1Hat(1:3);
    b1Hat = parm1Hat(4:6);

    parm2Hat = mdlY2.Coefficients.Estimate;
    a2Hat = parm2Hat(1:3);
    b2Hat = parm2Hat(4:6);
end
