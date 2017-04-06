
grid = 'JUQUIA_MAG.gxf';
Ir = -26.77;
Dr = -15.88;

%Valor que deveria retornar: Im: 43.20, Dm: -175.2

j = 3:2:25;

for i=1:11
    [I6,I7,I8,I9,M{1,i}, im{1,i}, dm{1,i},gridAS{1,i}] = hellbig(grid,Ir,Dr,j(i));
end


%Extended Direct Method.
[n m ] = size(M{1,1});
SOLi = zeros(n,m);
SOLd = zeros(n,m);

for k=1:n
    for v=1:m

        for i=1:11
         s(i) = im{1,i}(k,v);
        end
    
    p = s(1);

        for j = 1:10
            if abs(s(j+1) - p) < 1*j
            SOLi(k,v) = SOLi(k,v) + 1;
            p = s(j+1);
            end
        end
    end
end


for k=1:n
    for v=1:m

        for i=1:11
         s(i) = dm{1,i}(k,v);
        end
    
    p = s(1);

        for j = 1:10
            if abs(s(j+1) - p) < 1*j
            SOLd(k,v) = SOLd(k,v) + 1;
            p = s(j+1);
            end
        end
    end
end








%figure
%contourf(dm1,10)
%title('Momento - janela 3x3');
%colorbar;



