function [x,y] = acoeff2cart_fractional(a,P,s,numPts)
xMax0=s;

xMax=fminsearch(@LOCAL_FitS,xMax0, optimset('Algorithm','active-set','TolFun',1e-2,'display','off'),a,P,s,numPts);
x=0:xMax/(numPts-1):xMax;
y=a.*(x.^P);
if ~isempty(x) && ~isempty(y)
    [x,y]=equidist(x,y,numPts);
else
    x = NaN; y = NaN;
end
end


function e=LOCAL_FitS(x,a,P,s,numPts)
% local fit function
    xVals=0:x/(numPts-1):x;
    yVals=a.*(xVals.^P);

    sGuess=arclength(xVals,yVals);
    e=abs(s-sGuess);

end