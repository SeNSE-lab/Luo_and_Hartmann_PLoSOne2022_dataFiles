function [J, fitWhisker, rotWhisker, coefficients] = revision_orthogonalRegression_quadratic(whisker, angle_min, angle_max)
%
% By Yifu
% 2022/06/14


% model y=ax^2
% generate initial conditions 1:
A0 = whisker(end,2)/whisker(end,1)^2;
% generate initial conditions 2:
fit_model = @(a,x) a*x.^2;
g = fittype(fit_model, 'coefficients', {'a'});
fitRes = fit(whisker(:,1), whisker(:,2), g, 'StartPoint', A0);
A0 = max([fitRes.a,1e-5]);
initial = [0, A0];

coefficients = fmincon(@(coef)generateWhiskerCostFunc(coef, whisker),...
        initial, ...% initial value
        [],[],[],[],[angle_min,0],[angle_max,Inf],[],...
        optimset('MaxIter',4000,'display','off'));

% J is the optimized residual sum of squres for this whisker
[J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coefficients, whisker);





% take actual whisker, coefficients to generate fitted model, calcualte
% sum of squares
function [J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coef, whisker)
    numPts = size(whisker, 1);
    rotWhisker = whisker*rot2(-coef(1),'deg');
    
    [fitX, fitY] = acoeff2cart_quadratic(coef(2), whiskerLength(rotWhisker), numPts);
    fitWhisker = [fitX', fitY'];
    
    % option 1: find J with minimum distance 
    %           from each point on "whisker" to any point on "fitWhisker"
    %            colVec   - rowVec
    x_diff = rotWhisker(:,1) - fitX;
    y_diff = rotWhisker(:,2) - fitY;
    J = sum(min(x_diff.^2 + y_diff.^2, [], 2));
    
    % option 2: find J with distance
    %           from each point on "whisker" to corresponding point on
    %           "fitWhisker"
%     residuals_x = fitWhisker(:,1) - whisker(:,1);
%     residuals_y = fitWhisker(:,2) - whisker(:,2);
%     J = sum(residuals_x.^2+residuals_y.^2);      % RSS
    

end


end
