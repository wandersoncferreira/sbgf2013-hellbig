function [x,y,sum_m] = centro_massa(imagem)
% Fun��o que calcula a massa e o centro de massa numa imagem binarizada,
% considerando unit�ria a massa de cada pixel preto.
% 
% A origem do referencial de uma imagem � o canto superior esquerdo, os 
% eixos s�o crescentes para baixo e para a esquerda, as unidades s�o em 
% pixel o que equivaler� a uma unidade dos ind�ces da matriz onde a imagem
% est� guardada.

[linhas, colunas] = size(imagem);
% Inicializar os numeradores dos centros de massa horizontal e vertical e 
% a massa total
num_x = 0;
num_y = 0;
sum_m = 0;
% Ciclo para percorrer toda a imagem

for i = 1:linhas,
    for j = 1:colunas,
% Para cada pixel preto encontrado actualizar os numeradores 
% dos centros de massa horizontal e vertical e a massa total 
if imagem(i,j) > 0.4*max(imagem(:));
    num_x = num_x + j*imagem(i,j);
    num_y = num_y + i*imagem(i,j);
    sum_m = sum_m + imagem(i,j);  
end
    end
end

% Se foi encontrado algum pixel preto na imagem calcular os centros de 
% massa horizontal e vertical, caso contr�rio atribuir o valor 0
if sum_m ~= 0
x = num_x/sum_m;
y = num_y/sum_m;
else
x = 0;
y = 0;
end