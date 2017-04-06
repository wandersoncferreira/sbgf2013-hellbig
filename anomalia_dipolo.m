function anom = anomalia_dipolo(Ir,Dr,Im,Dm,fgrid)
%
%ANOMALIA_DIPOLO(Ir,Dr,Im,Dm,fgrid)
%
%Calcula anomalia de intensidade total de dipolo
%sobre malha regular. Grava resultado no formato  
%Grid Golden Software (Surfer)  
%  
%  Ir: Inclinação do campo de referência (graus)
%  Dr: Declinação do campo de referência (graus)
%  Im: Inclinação do momento magnético (graus)
%  Dm: Declinação do momento magnético (graus)
%  fgrid: Nome do arquivo de saida (Grid Surfer)

close all

fout=fname_ext(fgrid,'grd');

nr=101;
nc=81;
h=10;
mult=5;

p=max(nr,nc);

delta=h*mult/(p/2-0.5);

ye=((1:nc)-(nc+1)/2)*delta;
xn=((1:nr)-(nr+1)/2)*delta;

[yye,xxn]=meshgrid(ye,xn);
zz=zeros(size(xxn));

x0=30;
y0=-20;
z0=h;

x1 = -30;
y1 = -20;

x2 = 0;
y2 = 30;

Im1 = Im + 20;
Dm1 = Dm + 20;

Im2 = Im - 20;
Dm2 = Dm - 20;

anom1=1e7*anom_dip_r2(x0,y0,z0,xxn,yye,zz,Im,Dm,Ir,Dr);
anom2=1e7*anom_dip_r2(x1,y1,z0,xxn,yye,zz,Im1,Dm1,Ir,Dr);
anom3=1e7*anom_dip_r2(x2,y2,z0,xxn,yye,zz,Im2,Dm2,Ir,Dr);

anom = anom1 + anom2 + anom3;

contour(ye,xn,anom,30);grid on;
axis equal
colorbar
save_grd(fout,anom,ye,xn);

fprintf('\nCÁLCULO DE ANOMALIA DE DIPOLO\n===============\n\n')
fprintf('Direção de referência\nIr: %.1f\nDr: %.1f\n\n',Ir,Dr);
fprintf('Direção do momento magnético\nIm: %.1f\nDm: %.1f\n\n',Im,Dm);
fprintf('Grid gravado em: %s\n',fout);
%-------------------------------------------
function T=anom_dip_r2(x0,y0,z0,x,y,z,Im,Dm,Ir,Dr);

%Dipolo localizado em x0,y0,z0
%Anomalia calculada em x,y,z
%Inclinações e declinações em graus

if any(size(x)~=size(y)) | any(size(x)~=size(z))
	error('*** ANOM_DIP_R2 *** x,y,z not all same size')
end

x=x-x0;
y=y-y0;
z=z-z0;

f=pi/180;

Im=Im*f;
Dm=Dm*f;
Ir=Ir*f;
Dr=Dr*f;

%a=[cos(Im)*sin(Dm);cos(Im)*cos(Dm);sin(Im)];
%b=[cos(Ir)*sin(Dr);cos(Ir)*cos(Dr);sin(Ir)];

a=[cos(Im)*cos(Dm);cos(Im)*sin(Dm);sin(Im)];
b=[cos(Ir)*cos(Dr);cos(Ir)*sin(Dr);sin(Ir)];

r=sqrt(x.^2+y.^2+z.^2);

t3=r.^(-3);
t5=r.^(-5);

Uxx=3*x.^2.*t5-t3;
Uyy=3*y.^2.*t5-t3;
Uzz=3*(z.^2).*t5-t3;
Uxy= 3*x.*y.*t5;
Uxz=3*x.*z.*t5;
Uyz=3*y.*z.*t5;

%Uxx	Uxy	Uxz
%Uxy	Uyy	Uyz
%Uxz	Uyz	Uzz

T= (a(1)*Uxx+a(2)*Uxy+a(3)*Uxz)*b(1)+...
   (a(1)*Uxy+a(2)*Uyy+a(3)*Uyz)*b(2)+...
   (a(1)*Uxz+a(2)*Uyz+a(3)*Uzz)*b(3);
%-------------------------------------------
function save_grd(fname,grd,xc,yr)
%Function dummy=SAVE_GRD(fname,grid,xc,yr)
%
%	Saves GoldenSoftware grid file (binary format).
 
 
fid=fopen(fname,'w');
zmin=min(min(grd(isfinite(grd))));
zmax=max(max(grd(isfinite(grd))));
gr=grd';
i=find(isnan(gr));
gr=single(gr);
gr(i)=ones(size(i))*1.70141e+38;
[nc,nr]=size(gr);
 
id=['DSBB']';
lx=length(xc);
ly=length(yr);
fwrite(fid,id,'uchar');
fwrite(fid,nc,'short');
fwrite(fid,nr,'short');
fwrite(fid,[xc(1);xc(lx);yr(1);yr(ly);zmin;zmax],'double');
fwrite(fid,gr,'float');
fclose(fid);
%-------------------------------------------
function rslt=fname_ext(fin,ext)
%Checks and appends extension to filename

if ~isstr(fin);error('Fname must be a string');end
if ~isstr(ext);error('Extension must be a string');end
ip=findstr(fin,'.');

if ~isempty(ip)
	act=fin(ip(end)+1:end);
	if strcmp(upper(act),upper(ext))
		rslt=fin;
	else
		rslt=[fin,'.',ext];
	end
else
	rslt=[fin,'.',ext];
end
%-------------------------------------------
%-------------------------------------------

