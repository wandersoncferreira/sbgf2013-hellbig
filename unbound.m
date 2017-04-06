function gridR = unbound(gridp,k)

[n2 m2] = size(gridp);

%Retirando as linhas colocadas pela função flipper.

gridp(n2-(k-1):n2, :) = [];
gridp(:, m2-(k-1):m2) = [];
gridp(1:k,:)= [];
gridp(:,1:k) = [];

gridR = gridp;
