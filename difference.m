function [SOLi, SOLd] = difference(Xi,Yi,Xd,Yd,SOLi,SOLd,a,ii,dd)


Si = abs( Yi - Xi   );
Sd = abs( Yd - Xd   );
 




pi = Si < 0.5*a; 
pd = Sd < 0.5*a; 


SOLi(pi) = SOLi(pi) + 1;
SOLd(pd) = SOLd(pd) + 1;





