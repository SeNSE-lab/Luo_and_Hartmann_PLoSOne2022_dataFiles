function [J, fitWhisker, rotWhisker, coefficients] = revision_orthogonalRegression_poly5(whisker, angle_min, angle_max)
%
% By Yifu
% 2022/06/15


% model y=ax^5+bx^4+cx^3+dx^2
% generate initial conditions 1:
A0 = whisker(end,2)/whisker(end,1)^5/4;
B0 = whisker(end,2)/whisker(end,1)^4/4;
C0 = whisker(end,2)/whisker(end,1)^3/4;
D0 = whisker(end,2)/whisker(end,1)^2/4;
% generate initial conditions 2:
fit_model = @(a,b,c,d,x) a*x.^5 + b*x.^4 + c*x.^3 + d*x.^2;
g = fittype(fit_model, 'coefficients', {'a','b','c','d'});
fitRes = fit(whisker(:,1), whisker(:,2), g, 'StartPoint', [A0, B0, C0, D0]);
A0 = fitRes.a;
B0 = fitRes.b;
C0 = fitRes.c;
D0 = fitRes.d;
initial = [0,A0, B0, C0, D0];

coefficients = fmincon(@(coef)generateWhiskerCostFunc(coef, whisker), ...
        initial,...
        [],[],[],[],[angle_min,-Inf,-Inf,-Inf,-Inf],[angle_max,Inf,Inf,Inf,Inf],[],...
        optimset('MaxIter',4000,'display','off'));
    
% J is the optimized residual sum of squres for this whisker
[J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coefficients, whisker);



% take actual whisker, coefficients to generate fitted model, calcualte
% sum of distances
function [J, fitWhisker, rotWhisker] = generateWhiskerCostFunc(coef, whisker)
    numPts = size(whisker, 1);
    rotWhisker = whisker*rot2(-coef(1),'deg');
    
    [fitX, fitY] = acoeff2cart_poly5(coef(2), coef(3), coef(4), coef(5), whiskerLength(rotWhisker), numPts);
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
