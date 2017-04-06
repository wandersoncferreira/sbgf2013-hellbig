function [Grid,X,Y] = load_gxf(grid_name)


%A função recebe como parâmetro o nome do grid em GXF e retorna uma matriz
%com os dados de intensidade e duas matrizes X,Y referentes as coordenadas espaciais X,Y de cada ponto.
% 28/03/2013 Foram adicionados controles de existencia da entrada grid_name
% e variavel que calcula o tempo de processamento.

if ~exist(grid_name,'file')
    error(['***Arquivo não encontrado: ', grid_name]);
end

%Obtendo parâmetros do arquivo.
fid = fopen(grid_name);
leitura= textscan(fid,'%s');

fclose(fid);

leitura = leitura{1};

t = size(leitura,1);

for i =1:t
    switch leitura{i}
        case '#POINTS'
            POINTS = sscanf(leitura{i+1}, '%f');
        case '#ROWS'
            ROWS = sscanf(leitura{i+1}, '%f');
        case '#XORIGIN'
            XORIGIN = sscanf(leitura{i+1}, '%f');
        case '#YORIGIN'
            YORIGIN = sscanf(leitura{i+1}, '%f');
        case '#PTSEPARATION'
            PTSEPARATION = sscanf(leitura{i+1}, '%f');
        case '#RWSEPARATION'
            RWSEPARATION = sscanf(leitura{i+1},'%f');
        case '#DUMMY'
            DUMMY = sscanf(leitura{i+1},'%f');
        case '#GRID'
            Grid = leitura(i+1:end);
    end
    
    
end
    

%Gerando as coordenadas X,Y do grid.
i=1:1:POINTS;
    X = XORIGIN + (i -1)*PTSEPARATION;

i = 1:1:ROWS;
    Y = YORIGIN + (i-1)*RWSEPARATION;
    

%Limpando variaveis que não serão usadas mais.
clear t i PTSEPARATION RWSEPARATION XORIGIN YORIGIN 



%Converter grid de cell para double. A função str2double é extremanente
%lenta.
S = sprintf('%s*', Grid{:});
Grid = sscanf(S, '%f*');

Grid(Grid==DUMMY) = NaN;
Grid = reshape(Grid,POINTS,ROWS)';


