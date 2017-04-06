function imagem_bw = binariza(imagem, limiar)
% Função que binariza uma imagem em níveis de cinzento guardada numa
% matriz cujos elementos são da classe uint8.
%
% A função recebe um valor de limiar acima do qual o pixel será tornado
% completamente branco, e abaixo do qual sera tornado completamente preto.
%
% A origem do referencial de uma imagem é o canto superior esquerdo, os 
% eixos são crescentes para baixo e para a esquerda, as unidades são em 
% pixel o que equivalerá a uma unidade dos indíces da matriz onde a imagem
% está guardada.

[linhas, colunas] = size(imagem);
% Pré-alocação de memória
imagem_bw(linhas, colunas) = 0;
% Ciclo para percorrer toda a imagem
for i = 1:linhas,
for j = 1:colunas,
if imagem(i,j) > limiar
% Se o valor do pixel actual for maior que o valor de limiar 
% torná-lo completamente branco
imagem_bw(i,j) = 0;
else
% Caso contrário torná-lo completamente preto
imagem_bw(i,j) = 255;
end
end
end
