function [J, fitWhisker, rotWhisker, coefficients] = revision_orthogonalRegression_ellipse(whisker, angle_min, angle_max)
%
% By Yifu
% 2022/06/15


% model ellipse
% generate initial conditions:
fit_model = @(a,b,x) b * (1-sqrt(1-x.^2/a^2));
g = fittype(fit_model, 'coefficients', {'a','b'});
fitRes = fit(whisker(:,1), whisker(:,2), g, 'StartPoint', [200, 5],...
                'Lower', [whisker(end,1), whisker(end,2)],...
                'Upper', [5*whisker(end,1), abs(5*whisker(end,2))]...
                );
A0 = fitRes.a;
B0 = fitRes.b;
initial = [0, A0, B0];

coefficients = fmincon(@(coef)generateWhiskerCostFunc(coef, whisker),...
        initial, [],[],[],[],...
        [angle_min,whisker(end,1),whisker(end,2)],...
        [angle_max,5*whisker(end,1),abs(5*whisker(end,2))],[],...
        optimset('MaxIter',4000,'display','off'));

% J is the optimized sum of squres for this whisker
[J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coefficients, whisker);



% take actual whisker, coefficients to generate fitted model, calcualte
% sum of distances
function [J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coef, whisker)
    numPts = size(whisker, 1);
    rotWhisker = whisker*rot2(-coef(1),'deg');
    
    [fitX, fitY] = acoeff2cart_ellipse(coef(2), coef(3), rotWhisker(end,1), numPts);
    fitWhisker = [fitX', fitY'];
    
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
