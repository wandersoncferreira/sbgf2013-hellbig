function imagem_bw = binariza(imagem, limiar)
% Fun��o que binariza uma imagem em n�veis de cinzento guardada numa
% matriz cujos elementos s�o da classe uint8.
%
% A fun��o recebe um valor de limiar acima do qual o pixel ser� tornado
% completamente branco, e abaixo do qual sera tornado completamente preto.
%
% A origem do referencial de uma imagem � o canto superior esquerdo, os 
% eixos s�o crescentes para baixo e para a esquerda, as unidades s�o em 
% pixel o que equivaler� a uma unidade dos ind�ces da matriz onde a imagem
% est� guardada.

[linhas, colunas] = size(imagem);
% Pr�-aloca��o de mem�ria
imagem_bw(linhas, colunas) = 0;
% Ciclo para percorrer toda a imagem
for i = 1:linhas,
for j = 1:colunas,
if imagem(i,j) > limiar
% Se o valor do pixel actual for maior que o valor de limiar 
% torn�-lo completamente branco
imagem_bw(i,j) = 0;
else
% Caso contr�rio torn�-lo completamente preto
imagem_bw(i,j) = 255;
end
end
end
