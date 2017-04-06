function [I6, I7, I8, I9,Mag, ir, dm,gridAS] = hellbig(grid_gxf,Ir,Dr,d)


%Lendo o grid que será trabalhado.
[grid,X,Y] = load_gxf(grid_gxf); 


gridAS = AS(grid,X,Y);

%Gerando as componentes do campo magnético.
[Grid_X,Grid_Y,Grid_Z] = H2Hi(grid,Ir,Dr,X,Y);


%Função que gera os pesos em X e em Y. Definindo o zero em relação ao
%centro do grid.


%Gerando os valores das integrais.
[I1,I2,I3,I4,I5,I6, I7, I8, I9] = integration(Grid_X, Grid_Y, Grid_Z,d);


%Calculando o módulo da magnetização total

Mag = sqrt(I8.^2 + I9.^2 + 0.25*(I6 + I7).^2);

%Calculando a inclinação

ir = asin(((I6 + I7)./(2*Mag)));
ir = ir.*(180/pi);

%Calculando a declinação

dm = atan2(I8,I9);
dm = dm.*(180/pi);


function [HxR, HyR, HzR] = H2Hi(grid,Ir,Dr,X,Y)

%Lendo o grid que será utilizado posteriormente.
[n,m] = size(grid);


grid = boundary(grid,max(floor(0.05*n),10));
%Aplicando a Transformada de Fourier no grid
[gridFFT] = fft_gxf(grid,X,Y);

%Aplicando os filtros de transformação de componentes.
Hxx = gridFFT.*H2Hx(Ir,Dr);
Hyy = gridFFT.*H2Hy(Ir,Dr);
Hzz = gridFFT.*H2Hz(Ir,Dr);


%Aplicando a transformada Inversa de Fourier;
Hx = ifft2(Hxx);
Hy = ifft2(Hyy);
Hz = ifft2(Hzz);


%Desfazendo o flipper.
Hx = unbound(Hx,max(floor(0.05*n),10));
Hy = unbound(Hy,max(floor(0.05*n),10));
Hz = unbound(Hz,max(floor(0.05*n),10));

%Pegando as dimensões corretas do grid original, após os calculos da FFT
HxR = real(Hx(1:n,1:m));
HyR = real(Hy(1:n,1:m));
HzR = real(Hz(1:n,1:m));

function HX = H2Hx(Ir,Dr)
global G_uu G_vv G_rr
    
f = pi/180;
Ir = Ir*f;
Dr = Dr*f;


v1 = sin(Ir);
v2 = cos(Ir)*sin(Dr);
v3 = cos(Ir)*cos(Dr);

G_uu(1,1) = 1.0; G_vv(1,1) = 1.0; G_rr(1,1) = 1.0;
HX = 1i*G_uu./(G_rr*v1 + 1i*(G_uu*v2 + G_vv*v3));


function HY = H2Hy(Ir,Dr)
global G_uu G_vv G_rr

f = pi/180;
Ir = Ir*f;
Dr = Dr*f;

v1 = sin(Ir);
v2 = cos(Ir)*sin(Dr);
v3 = cos(Ir)*cos(Dr);

G_uu(1,1) = 1.0; G_vv(1,1) = 1.0; G_rr(1,1) = 1.0;
HY = 1i*G_vv./(G_rr*v1 + 1i*(G_uu*v2 + G_vv*v3));

function HZ = H2Hz(Ir,Dr)
global G_uu G_vv G_rr

f = pi/180;
Ir = Ir*f;
Dr = Dr*f;

v1 = sin(Ir);
v2 = cos(Ir)*sin(Dr);
v3 = cos(Ir)*cos(Dr);

G_uu(1,1) = 1.0; G_vv(1,1) = 1.0; G_rr(1,1) = 1.0;
HZ = G_rr./(G_rr*v1 + 1i*(G_uu*v2 +  G_vv*v3));


function [H1,H2,H3,H4,H5,H6,H7,H8,H9] = integration(Grid_X,Grid_Y,Grid_Z,d)
[n m] = size(Grid_X);

C = (-1/(2*pi));

%Peso da Integração de Simpson's, vetor coluna
peso = 2*ones(d,1);
peso(2:2:end-1)=4;
peso(1)=1;
peso(end)=1;


%Vetor que define os pesos de cada janela em relação a x e y.
vWindow = (-floor(d/2):1:floor(d/2));

%Separação entre as coordenadas

hx = 100;
hy = 100;


%Define os grids resultado de cada integral de Helbig

H1 = NaN(n,m);
H2 = NaN(n,m);
H3 = NaN(n,m);
H4 = NaN(n,m);
H5 = NaN(n,m);
H6 = NaN(n,m);
H7 = NaN(n,m);
H8 = NaN(n,m);
H9 = NaN(n,m);



for i=1:n-(d-1)
    for j=1:m-(d-1)
   
        %Grids não corrigidos
    gX = Grid_X(i:(d-1+i),j:(d-1+j));
    gY = Grid_Y(i:(d-1+i),j:(d-1+j));
    gZ = Grid_Z(i:(d-1+i),j:(d-1+j));
    
    
    
    
    %Calculo das integrais não corrigidas
    I1 = hx*hy*peso'*(gX*peso)/9;
    I2 = hx*hy*peso'*(gY*peso)/9;
    I3 = hx*hy*peso'*(gZ*peso)/9;

    
    %grid de uns para correção de I1 a I3, ou seja, retirar um termo
    %constante;
    u = ones(d,d);
    Iu = hx*hy*peso'*(u*peso)/9;
    
    %Correção de I1, I2 e I3
    gXc = gX - I1/Iu*u;
    gYc = gY - I2/Iu*u;
    gZc = gZ - I3/Iu*u;
    
    
    %Grids corrigidos com I1, I2, I3, mas não corrigidos por I4 e I5;
    xgYc = gYc*diag(vWindow);
    ygXc = (gXc'*diag(vWindow(:)))';
    
    
    %Calculo de I4 e I5 utilizando a correção anterior, a fim de corrigir a
    %tendencia 
    
    I4 = hx*hy*peso'*(ygXc*peso)/9;
    I5 = hx*hy*peso'*(xgYc*peso)/9;
 
    
    
    %Correção de I4.
    I4l = hx*hy*peso'*((ones(d)'*diag(vWindow(:).*vWindow(:)))'*peso)/9;
    gX_final = gXc - (ones(d)'*diag(vWindow(:)*I4/I4l))';
    
    
    %Correção de I5
    I5l = hx*hy*peso'*((ones(d)*diag(vWindow.*vWindow))*peso)/9;
    gY_final = gYc - ones(d)*diag(vWindow*I5/I5l);
    
    
    
    
    %Criando conjunto de Grids finais Corrigidos com I1,I2, I3, I4 e I5.
    
    xgX_final = gX_final*diag(vWindow);
    %xgY_final = gY_final*diag(vWindow);
    xgZ_final = gZc*diag(vWindow);
    %ygX_final = (gX_final'*diag(vWindow(:)))';
    ygY_final = (gY_final'*diag(vWindow(:)))';
    ygZ_final = (gZc'*diag(vWindow(:)))';
    
    
    
    %Integrais de 1 a 6 deveriam ser próximas de zero.
    %H1(floor(d/2)+i,floor(d/2)+j) = hx*hy*peso'*(gXc*peso)/9;
    
    %H2(floor(d/2)+i,floor(d/2)+j) = hx*hy*peso'*(gYc*peso)/9;
    
    %H3(floor(d/2)+i,floor(d/2)+j) = hx*hy*peso'*(gZc*peso)/9;
    
    %H4(floor(d/2)+i,floor(d/2)+j) = hx*hy*peso'*(ygX_final*peso)/9;
    
    %H5(floor(d/2)+i,floor(d/2)+j) = hx*hy*peso'*(xgY_final*peso)/9;
    
    %Integrais de 6 a 9.
    
    H6(floor(d/2)+i,floor(d/2)+j) = C*hx*hy*peso'*(xgX_final*peso)/9;
    
    H7(floor(d/2)+i,floor(d/2)+j) = C*hx*hy*peso'*(ygY_final*peso)/9;
    
    H8(floor(d/2)+i,floor(d/2)+j) = C*hx*hy*peso'*(xgZ_final*peso)/9;
    
    H9(floor(d/2)+i,floor(d/2)+j) = C*hx*hy*peso'*(ygZ_final*peso)/9;
    
    

    end
    
end

function gridAS = AS(grid,X,Y)
global G_uu G_vv G_rr

[n,m] = size(grid);

gridFFT = fft_gxf(grid,X,Y);

FdelX = (1i*G_uu).*gridFFT;
FdelY = (1i*G_vv).*gridFFT;
FdelZ = (G_rr).*gridFFT;

delX = ifft2(FdelX);
delY = ifft2(FdelY);
delZ = ifft2(FdelZ);

delX = real(delX(1:n,1:m));
delY = real(delY(1:n,1:m));
delZ = real(delZ(1:n,1:m));

gridAS = sqrt((delX.^2  +  delY.^2  + delZ.^2));