function [f] = myfunc(x,n,m,o,net,inputs,targets)  
k=0;
for i=1:n
     for j=1:m
        k=k+1;
        xi(i,j)=x(k);
     end
end
for i=1:n
     k=k+1;
     xl(i)=x(k);
     xb1(i,1)=x(k+n);
end
for i=1:o
     k=k+1;
     xb2(i,1)=x(k);
end
net.iw{1,1}=xi;
net.lw{2,1}=xl;
net.b{1,1}=xb1;
net.b{2,1}=xb2;
f=sum((net(inputs)-targets).^2)/length(inputs); 