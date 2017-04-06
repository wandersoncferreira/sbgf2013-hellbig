function y = costapp(n,m)
if m >= 1,
    y=[0.5*(1 - cos(pi/m.*[0:m-1]')); ones(n-2*m,1); 0.5*(1-cos(pi/m.*[m-1:-1:0]'))];
else
    y = ones(n,1);
    
end