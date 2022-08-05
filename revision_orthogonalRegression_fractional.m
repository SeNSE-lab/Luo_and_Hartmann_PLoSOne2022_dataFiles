function [J, fitWhisker, rotWhisker, coefficients] = revision_orthogonalRegression_fractional(P, whisker, angle_min, angle_max)
%revision_orthogonalRegression_fractional takes a whisker (Mx2) and fit it to y=Ax^2,
%by model ii regression (or reduced major axis regression). The function
%returns the optimal coefficients of the model for the regression, as well
%as the fitted whisker. The optimization method is MATLAB fmincon.
%
%       [J, x, y, coef] = orthogonalRegression_fractional(whisker)
%
% By Yifu
% 2021/12/15


% model y=ax^P
% generate initial conditions 1:
A0 = whisker(end,2)/whisker(end,1)^P;
% generate initial conditions 2:
fit_model = @(a,x) a*x.^P;
g = fittype(fit_model, 'coefficients', {'a'});
fitRes = fit(whisker(:,1), whisker(:,2), g, 'StartPoint', A0);
A0 = fitRes.a;


coefficients = fmincon(@(coef)generateWhiskerCostFunc(coef, P, whisker), ...
        [0, A0],...
        [],[],[],[],[angle_min, 0],[angle_max, Inf],[],...
        optimset('MaxIter',10000,'display','off'));

% J is the optimized residual sum of squres for this whisker
[J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coefficients, P, whisker);




% take actual whisker, coefficients to generate fitted model, calcualte
% sum of squares
function [J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coef, P, whisker)
    numPts = size(whisker, 1);
    rotWhisker = whisker*rot2(-coef(1),'deg');
    
    [fitX, fitY] = acoeff2cart_fractional(coef(2), P, whiskerLength(whisker), numPts);
    fitWhisker = [fitX; fitY]';
    
    % option 1: find J with minimum distance 
    %           from each point on "whisker" to any point on "fitWhisker"
    x_diff = rotWhisker(:,1) - fitX;
    y_diff = rotWhisker(:,2) - fitY;
    J = sum(min(x_diff.^2 + y_diff.^2, [], 2));
    
    
    % option 2: find J with distance
    %           from each point on "whisker" to corresponding point on "fitWhisker"
%     residuals_x = fitWhisker(:,1) - whisker(:,1);
%     residuals_y = fitWhisker(:,2) - whisker(:,2);
%     J = sum(residuals_x.^2+residuals_y.^2);      % RSS
    

end


end
