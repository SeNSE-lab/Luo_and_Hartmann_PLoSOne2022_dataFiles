function [x,y] = acoeff2cart_poly5(a, b, c, d, s, numPts)
xMax0=s;

xMax=fminsearch(@LOCAL_FitS,xMax0, optimset('Algorithm','active-set','TolFun',1e-2,'display','off'),a,b,c,d,s,numPts);
x=0:xMax/(numPts-1):xMax;
y=a.*(x.^5)+b.*(x.^4)+c.*(x.^3)+d.*(x.^2);
if ~isempty(x) && ~isempty(y)
    [x,y]=equidist(x,y,numPts);
else
    x = NaN; y = NaN;
end
end


function e=LOCAL_FitS(x, a, b, c, d, s, numPts)
% local fit function
    xVals=0:x/(numPts-1):x;
    yVals=a.*(xVals.^5) + b.*(xVals.^4) + c.*(xVals.^3) + d.*(xVals.^2);

    sGuess=arclength(xVals,yVals);
    e=abs(s-sGuess);

end