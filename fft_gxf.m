function gridFFT = fft_gxf(grid,x,y)
global G_uu G_vv G_rr  
%Aplicando o filtro de Redução ao Pólo.
%Im = 38.91; Dm = 53.36; Ir = - 29.71; Dr = -20.58;

   
%Separando as dimensões do grid inserido.
[n,m] = size(grid);


%Criando uma janela para minimizar os efeitos de borda (Leakage);
wn = costapp(n,max(floor(0.01*n),5));
wm = costapp(m,max(floor(0.01*m),5))';
window = wn(:,ones(1,m)).*wm(ones(n,1),:);


%Calculando os incrementos.
dx = mean(diff(x));
dy = mean(diff(y));


%Calculando o número de amostrar. Usa-se a potencia de 2 mais próxima do
%valor real, pois assim a DFT acelera seus calculos se tornando a FFT.
n2 = 2^nextpow2(n);
m2 = 2^nextpow2(m);

%Calcula o incremento da frequencia.
du = 2*pi/(m2*dx); 
dv = 2*pi/(n2*dy);

%Criando a matriz grid que será utilizada para calcular a FFT.
grid= grid.*window; %Aplicando a Janela ao grid de dados.
media = mean(mean(grid)); 
gr(1:n,1:m) = grid - media;


%Aplicando a FFT no grid após reduzir sua média..
gridFFT = fft2(gr,n2,m2);

%Centralizar o espectro em torno do zero.
u = fftshift((-m2/2: m2/2 -1)*du);
v = fftshift((-n2/2:n2/2 -1)'*dv);


UU =u(ones(n2,1),:);
VV =v(:,ones(1,m2));
RR = sqrt(UU.^2 + VV.^2);


%Definindo valores para variaveis globais que serão utilizadas por outras
%funções.
G_uu = UU;
G_vv = VV;
G_rr = RR;