function [grid_p,k] = boundary(grid,k)

%Onde k é o tamanho da "continuação" da anomalia. Devemos observar que a
%"continuação" se faz por meio de matrizes quadradas e dessa forma, como um
%grid não é usualmente quadrado, não é possivel adicionar uma réplicas
%exatas do grid com suas dimensões.


[n m] = size(grid);

if (k > m || k > n)
    error(['***Dimensões erradas, a "continuação" não pode exceder a seguinte dimensão: ', num2str(min(n,m))]);   
end


%Matrizes das Laterais

M1 = grid(1:n, 1:k);
M2 = grid(1:n, m -(k-1):m);

V1 = grid(1:k, 1:m);
V2 = grid(n - (k-1):n, 1:m);

%Matrizes dos vertices

Q1 = grid(1:k, 1:k);
Q2 = grid(1:k, m -(k-1):m);

Q3 = grid(n -(k-1):n, 1:k);
Q4 = grid(n -(k-1):n, m- (k-1):m);

[p,q] = size(Q3);
[v w] = size(Q2);
[x z] = size(Q4);


%Invertendo as matrizes antes de uní-las

M1 = fliplr(M1);
M2 = fliplr(M2);

V1 = flipud(V1);
V2 = flipud(V2);

Q1 = Q1(k:-1:1,k:-1:1);
Q3 = Q3(p:-1:1,q:-1:1);

Q2 = Q2(v:-1:1,w:-1:1);
Q4 = Q4(x:-1:1,z:-1:1);


%Unindo a nova matriz.
grid_p = [ Q1 V1 Q2; M1 grid M2; Q3 V2 Q4];