N = 100;
d = 2;
h1 = 0.3;
beta = 0.5;
S = 3;
proc = [5];

rng(10);

x = rand(N,d);
x = bsxfun(@minus,x,min(x));
x = bsxfun(@rdivide,x,max(x));

% figure;
% plot(x(:,1),x(:,2),'.','markersize',12); 
% axis equal tight;

fileID = fopen('data.txt','w');
dlmwrite('data.txt',[N d],'delimiter',' ','-append')
dlmwrite('data.txt',proc,'delimiter',' ','-append')
dlmwrite('data.txt',[h1 beta S],'delimiter',' ','-append')
dlmwrite('data.txt',x,'delimiter',' ','-append')
fclose(fileID);