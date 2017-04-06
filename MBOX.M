function g=mbox(x,y,z,xb,yb,zb,Im,Dm,Ir,Dr);
% T=mbox(x,y,z,xb,yb,zb,Im,Dm,Ir,Dr)
%Total intensity magnetic anomaly of a rectangular prism
%described by xb, yb, zb computed at x, y, z.
%
%Im,Dm : Inclination & declination of magnetization (degrees)
%Ir,Dr : Inclination & declination of reference field (degrees)

[n,m]=size(z);
if n*m==1, z=z*ones(size(x));end

f=pi/180;Im=Im*f;Ir=Ir*f;Dm=(90-Dm)*f;Dr=(90-Dr)*f;
L=cos(Im)*cos(Dm);
M=cos(Im)*sin(Dm);
N=sin(Im);
l=cos(Ir)*cos(Dr);
m=cos(Ir)*sin(Dr);
n=sin(Ir);
flx=(m*N+n*M)/2;fly=(l*N+n*L)/2;flz=(l*M+m*L)/2;
ftx=-l*L;fty=-m*M;ftz=-n*N;

g=zeros(size(x));
for i=1:2
      xx=xb(i)-x;
      for j=1:2
        yy=yb(j)-y;
        for k=1:2
          zz=zb(k)-z;
          sgn=(-1)^(i+j+k);
          R=sqrt(xx.^2+yy.^2+zz.^2);
          g=g+sgn*( flx*log((R+xx)./(R-xx))...
                   +fly*log((R+yy)./(R-yy))...
                   +flz*log((R+zz)./(R-zz))...
                   +ftx*atan(yy.*zz./xx./R)...
                   +fty*atan(zz.*xx./yy./R)...
                   +ftz*atan(xx.*yy./zz./R));
        end
      end
end
