 clear all; close all;


grid = 'norte.gxf';
Ir = -12;
Dr = -14;

%Valor que deveria retornar: Im: -27.6, Dm: -18.6

j = 3:2:25;

for i=1:12
    [I6,I7,I8,I9,M{1,i}, im{1,i}, dm{1,i}, gridAS{1,i}] = hellbig(grid,Ir,Dr,j(i));
end

%Extended Direct Method.
[n m ] = size(M{1,1});
SOLi = zeros(n,m);
SOLd = zeros(n,m);

SOLi = SOLi(:);
SOLd = SOLd(:);


 ii = zeros(n,m);
 dd = zeros(n,m);
 
Teste = zeros(n,m);

set(0,'DefaultFigureWindowStyle','docked');
figure; imagesc(gridAS{1,1}); axis 'xy'; colorbar
handle1 = str2double(input('Quantas anomalias foram localizadas?: ','s'));
fprintf('Selecione a anomalia da melhor forma possível.\n');


for i=1:handle1

 k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    k = waitforbuttonpress;
    point2 = get(gca,'CurrentPoint');    % button up detected
    finalRect = rbbox;
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    
    
    x = round(x);
    y = round(y);
    

    p{i} = [y(1), y(3), x(1), x(3)];
    
    hold on
    
    Teste(p{i}(1):p{i}(2),p{i}(3):p{i}(4)) = gridAS{1,1}(p{i}(1):p{i}(2),p{i}(3):p{i}(4));
    [px,py,massa] = centro_massa(Teste);
    px = round(px);
    py = round(py);
    plot(px,py,'+r','MarkerSize',14,'LineWidth',4,'MarkerEdgeColor','k');
    
    Teste=zeros(n,m);
   
    
    axis manual
    v = plot(x,y,'--r','linewidth',3);
    text(x(1)-2,y(3)+ 10,['Anomalia', num2str(i)],'FontSize', 12, 'FontWeight', 'bold','BackgroundColor',[.7 .9 .7]); 

    hold off
    

end


fprintf('A cruz na imagem, indica o CENTRO DE MASSA da anomalia selecionada.\n');
fprintf('Selecione uma região ao redor do CENTRO DE MASSA\n.');




NN=11;

for kk=1:NN
	for i=1:NN-kk+1
    [SOLi, SOLd] = difference(im{1,i},im{1,i+kk},dm{1,i},dm{1,i+kk},SOLi,SOLd,kk);
	end
end


SOLi = reshape(SOLi,n,m);
SOLd = reshape(SOLd,n,m);

Soli = zeros(n,m);
Sold = zeros(n,m);



M = size(SOLi,1);
idx = sub2ind(size(SOLi),px,py);
offset = [-1, M, 1, -M,(M+1),(-M+1),(-M-1),(M-1)];
neighb = bsxfun(@plus,idx,offset);
solucoesI = SOLi(neighb);




for i=1:handle1
    
    fprintf('ANOMALIA %.0f: \n \n Inclinação (graus): %.2f +- %.2f\n Declinação (graus): %.2f +- %.2f\n\n',i,si(i),dvi(i),sd(i),dvd(i));
    
end





figure;
imagesc(SOLi);
axis 'xy'
title('Solução para a Inclinação');
colorbar;

figure;
imagesc(SOLd);
axis 'xy'
title('Solução para a Declinação');
colorbar;



%figure
%contourf(dm1,10)
%title('Momento - janela 3x3');
%colorbar;



