close all

load '../code/serial/listcentres.out'
load '../code/serial/numcentres.out'

numcentres = [0; numcentres];
listcentres = listcentres + 1;

for i = 1:S
    figure(i); hold all;
    plot(x(:,1),x(:,2),'.b','markersize',12); 
    plot(x(listcentres(numcentres(i)+1:numcentres(i+1)),1),...
        x(listcentres(numcentres(i)+1:numcentres(i+1)),2),...
        'xr','markersize',12,'linewidth',5); 
    axis equal tight;
    D = pdist2(x(listcentres(numcentres(i)+1:numcentres(i+1)),:),x(listcentres(numcentres(i)+1:numcentres(i+1)),:));
    min(min(D(D>0)))
end