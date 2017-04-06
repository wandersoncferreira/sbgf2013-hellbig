function RTP

global    X Y UU VV


y = inputdlg({'Input GRID GXF','Output GRID GXF','Total Inclination','Total Declination','IGRF Inclination','IGRF Declination'},'Reduction to the Pole Operator',[1 38; 1 38; 1 38; 1 38; 1 38; 1 38]);
    
    if mean(size(y)) == 0
           return
    end
    
    GRD_Input = y{1};
    GRD_Output = y{2};
    Im = str2double(y{3});
    Dm = str2double(y{4});
    Ir = str2double(y{5});
    Dr = str2double(y{6});

% [1] Carregando o grid e suas coordenadas

[GR,X,Y] = load_gxf(GRD_Input);


% [2] Aplicando a FFT 
gridFFT = fft_gxf(GR);

% [3] Aplicando o filtro de Redução ao Pólo.

ff = redpole(Im,Dm,Ir,Dr);
grid_rpole = gridFFT.*ff;
grid_rpole(1,1) = 0.0;

% [4] Tomando a transformada Inversa de Fourier.

grid_final = ifft2(grid_rpole);

% [5] Visualizando os resultados
[n,m] = size(GR);
grid_final = real(grid_final(1:n,1:m));


flag = true;
while flag
    c = listdlg('PromptString','Qual gráfico você gostaria de visualizar?',...
                'SelectionMode','single', 'ListString',{'Grid.gxf','Grid_Reduzido'},...
                'Name','Visualizando os mapas','ListSize',[230 130]);
            
                if c == 1
                    figure
                    title('ORIGINAL ANOMALY');
                    pcolor(GR);
                    colorbar('EastOutside'); colormap('jet');
                end
                                            
                if c == 3
                    figure
                    title('REDUCED ANOMALY');
                    pcolor(grid_final);
                    colorbar('EastOutside'); colormap('jet');
                end
                
                if mean(size(c)) == 0
                    flag = false;
                end
end


%Salvando o grid reduzido para ser visualizado no GEOSOFT
save_gxf(GRD_Input,grid_final,GRD_Output);
    


%Aplicando o filtro de Redução ao Pólo.
%Im = 38.91; Dm = 53.36; Ir = - 29.71; Dr = -20.58;

function [gridFFT] = fft_gxf(grid)

global G_uu G_rr G_vv  X Y UU VV

%Lendo o grid e separando suas dimensões.
[n,m] = size(grid);


%Criando uma janela para minimizar os efeitos de borda (Leakage);
wn = costapp(n,max(floor(0.05*n),5));
wm = costapp(m,max(floor(0.05*m),5))';
window = wn(:,ones(1,m)).*wm(ones(n,1),:);

%Calculando os incrementos.
dx = mean(diff(X));
dy = mean(diff(Y));


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

%Transforma as componentes u e v calculadas acima em duas matrizes de
%tamanhos n2 e m2 respectivamente.
UU =u(ones(n2,1),:);
VV =v(:,ones(1,m2));

%Podemos chamar rr de módulo das componentes uu e vv.
rr = sqrt(UU.^2 + VV.^2);

%Definindo valores para variaveis globais que serão utilizadas por outras
%funções.
G_uu = UU;
G_vv = VV;
G_rr = rr;

function y = costapp(n,m)

if m >= 1,
    y=[0.5*(1 - cos(pi/m.*[0:m-1]')); ones(n-2*m,1); 0.5*(1-cos(pi/m.*[m-1:-1:0]'))];
else
    y = ones(n,1);
    
end

function ff = redpole(Im,Dm, Ir, Dr)

global G_uu G_vv G_rr 

f = pi/180;

%Transformar os valores dos ângulos para radianos.
Im = Im*f;
Dm = Dm*f;
Ir = Ir*f;
Dr = Dr*f;

%Definindo os termos da redução ao pólo para Declinação e Inclinação
%da Magnetização Remanescente.
a1 = cos(Im)*sin(Dm);
a2 = cos(Im)*cos(Dm);
a3 = sin(Im);

%Definindo os termos da redução ao pólo para Declinação e Inclinação 
%da Magnetização Referencia (IGRF).
b1 = cos(Ir)*sin(Dr);
b2 = cos(Ir)*cos(Dr);
b3 = sin(Ir);

G_uu(1,1) = 1.0; G_vv(1,1) = 1.0; G_rr(1,1) = 1.0;

ff = G_rr.^2./((1i*b1*G_uu + 1i*b2*G_vv + b3*G_rr).*(1i*a1*G_uu + 1i*a2*G_vv + a3*G_rr));

function [grid_gxf,X,Y] = load_gxf(grid_name)

%A função recebe como parâmetro o nome do grid em GXF e retorna uma matriz
%com os dados de intensidade e duas matrizes X,Y referentes as coordenadas espaciais X,Y de cada ponto.
% 28/03/2013 Foram adicionados controles de existencia da entrada grid_name
% e variavel que calcula o tempo de processamento.

if ~exist(grid_name,'file')
    error(['***Arquivo não encontrado: ', grid_name]);
end

%Obtendo parâmetros do arquivo.
leitura= textread(grid_name,'%s','whitespace',' \b\t');
t = size(leitura,1);

for i = 1:t  
    switch leitura{i}
        case '#POINTS'
            POINTS = str2double( leitura{i+1} );
        case '#ROWS'
            ROWS = str2double( leitura{i+1} );
        case '#XORIGIN'
            XORIGIN = str2double( leitura{i+1} );
        case '#YORIGIN'
            YORIGIN = str2double( leitura{i+1} );
        case '#PTSEPARATION'
            PTSEPARATION = str2double( leitura{i+1} );
        case '#RWSEPARATION'
            RWSEPARATION = str2double( leitura{i+1} );
        case '#DUMMY'
            DUMMY = str2double( leitura{i+1} );  
            break;
    end
end
    

%Gerando as coordenadas X,Y do grid.
i=1:1:POINTS;
    X = XORIGIN + (i -1)*PTSEPARATION;

i = 1:1:ROWS;
    Y = YORIGIN + (i-1)*RWSEPARATION;

X = X';
Y = Y';

%Limpando variaveis que não serão usadas mais.
clear leitura t i a b PTSEPARATION RWSEPARATION XORIGIN YORIGIN 

%Mantem os dados do cabeçalho e pega os dados do grid do arquivo gxf
Grid = textread(grid_name,'%s','headerlines',49);


%Transformar o Grid lido em uma Matriz ('#POINTS')Linhas x ('#ROWS')Colunas
M = size(Grid,1);

j = 1;
w = 1;
gr = ones(POINTS,ROWS);


for i = 1:M
    
    gr(j,w) = str2double(Grid{i});
    
    if gr(j,w) == DUMMY
        gr(j,w) = NaN;
    end
    
    if j == POINTS
        w = w + 1;
        j = 0;
        
        if w == ROWS + 1
            break
        end
    end
    
    j = j + 1;
end

grid_gxf = gr';

function save_gxf(grid_name,grid,save_name)


%A funçao recebe como parametros o nome do grid original em GXF, e a matriz do grid
%após suas alterações, ou seja, já pronta para ser salva. A saida da função será um
%novo grid em GXF chamado GRID_MATLAB.gxf


fid = fopen(grid_name,'r');
fad = fopen(save_name,'w+');

%Determinar o numero de linhas e colunas e o valor DUMMY.

leitura= textread(grid_name,'%s','whitespace',' \b\t');
t = size(leitura,1);


for i = 1:t
    switch leitura{i}
        case '#POINTS'
            POINTS = str2double( leitura{i+1} );
        case '#ROWS'
            ROWS = str2double( leitura{i+1} );
        case '#DUMMY'
            DUMMY = str2double( leitura{i+1} );
            break;
    end           
end
    

%Limpando variaveis que não serão mais utilizadas.

clear leitura M N t a

%Obtendo o grid que será salvo no novo arquivo.

gr = grid';
gr(isnan(gr)) = DUMMY;

%Salvar o grid após realizar operações sobre ele.Tomar cuidado pois os
%espaços são importantes na comparação de string abaixo.

i = 1;
w = 1;
flag = true;

while 1
    tline = fgetl(fid);
    fprintf(fad,tline);
    fprintf(fad,'\n');
    
    if strcmp(tline,'#GRID                                                                           ')  
        
        A = fgetl(fid);
        B = find(isspace(A));
        number = size(B,2) + 1;
        
        
        
        while flag
            
            %Imprimindo o DUMMY separadamente.
            
            if gr(i,w) == DUMMY
                fprintf(fad,'%.0e ', DUMMY);
            else
                fprintf(fad,'%.5f ', gr(i,w));
            end
            
            
            if mod(i,number) == 0
                fprintf(fad,'\n');
            end
        
            if i == POINTS
                fprintf(fad,'\n');
                i = 0;
                w = w + 1;
                if w == ROWS + 1
                    break
                end
            end
            i = i + 1;
        end       
        break
    end
end

fclose(fad);
fclose(fid);




