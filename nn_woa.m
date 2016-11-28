
clc
tic
close all 
clear all 
rng default   
filename = 'wdbc1.xlsx'; 
sheetname1 = 'Sheet1'; 
sheetname2 = 'Sheet2'; 
gbest = 0.0; 
input = xlsread(filename,sheetname1,'A1:Z10000'); 
target = xlsread(filename,sheetname2,'A1:Z10000');   
inputs=input'; 
targets=target';   
m=length(inputs(:,1)); 
o=length(targets(:,1));   
n=10; net=feedforwardnet(n); 
net=configure(net,inputs,targets); 
kk=m*n+n+n+o;   
for j=1:kk     
  LB(1,j)=-1.5;     
  UB(1,j)=1.5; 
end 
pop=10; 
for i=1:pop     
    for j=1:kk         
      xx(i,j)=LB(1,j)+rand*(UB(1,j)-LB(1,j));     
    end 
end   
maxrun=1; 
for run=1:maxrun     
    fun=@(x) myfunc(x,n,m,o,net,inputs,targets);
    x0=xx; 
 
% pso initialization----------------------------------------------start   
    x=x0;       % initial population     
    %v=0.1*x0;   % initial velocity     
    for i=1:pop         
        f0(i,1)=fun(x0(i,:));
    end     
    [fmin0,index0]=min(f0);     
    pbest=x0;               % initial pbest     
    gbest=x0(index0,:);     % initial gbest     
    % pso initialization------------------------------------------------end 
    % pso algorithm---------------------------------------------------start     
    %c1=1.5; 
    %c2=2.5;
    
    ite=1;
    maxite=100;
    tolerance=1;
    while ite<=maxite %&& tolerance>10^-8
          %w=0.1+rand*0.4;         % pso velocity updates 
          
          A1=2-ite*((2)/maxite);
          A2=-1+ite*((-1)/maxite);
          
        %  for i=1:pop          
        %      for j=1:kk
        %         v(i,j)=w*v(i,j)+c1*rand*(pbest(i,j)-x(i,j))...  
        %                  +c2*rand*(gbest(1,j)-x(i,j)); 
        %      end  
        %  end
          % pso position update
          for i=1:pop
              
                r1=rand(); % r1 is a random number in [0,1]
                r2=rand(); % r2 is a random number in [0,1]
        
                A=2*A1*r1-A1;  % Eq. (2.3) in the paper
                C=2*r2;      % Eq. (2.4) in the paper
        
        
                 b=1;               %  parameters in Eq. (2.5)
                 l=(A2-1)*rand+1;   %  parameters in Eq. (2.5)
        
                 p = rand();        % p in Eq. (2.6)
              
              for j=1:kk
                  
                if p<0.5   
                    if abs(A)>=1
                        rand_leader_index = floor(pop*rand()+1);
                        X_rand = x(rand_leader_index, :);
                        D_X_rand=abs(C*X_rand(j)-x(i,j));       % Eq. (2.7)
                        x(i,j)=X_rand(j)-A*D_X_rand;            % Eq. (2.8)
                    
                    elseif abs(A)<1
                        D_Leader=abs(C*pbest(j)-x(i,j));   % Eq. (2.1)
                        x(i,j)=pbest(j)-A*D_Leader;        % Eq. (2.2)
                    end
                
                elseif p>=0.5
              
                    distance2Leader=abs(pbest(j)-x(i,j));  % Eq. (2.5)
                    x(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+pbest(j);
                
                end
                  
                  
              %    x(i,j)=x(i,j)+v(i,j);   
              end  
          end
          % handling boundary violations
          for i=1:pop      
              for j=1:kk 
                if x(i,j)<LB(j)
                  x(i,j)=LB(j);      
                elseif x(i,j)>UB(j)
                  x(i,j)=UB(j); 
                end          
              end        
          end
          
          % evaluating fitness
          for i=1:pop 
             f(i,1)=fun(x(i,:));    
          end
          
          % updating pbest and fitness
          for i=1:pop      
            if f(i,1)<f0(i,1)    
               pbest(i,:)=x(i,:);
               f0(i,1)=f(i,1);        
            end      
          end     
           
          [fmin,index]=min(f0);   % finding out the best particle 
          ffmin(ite,run)=fmin;    % storing best fitness  
          ffite(run)=ite;         % storing iteration count  
          
          % updating gbest and best fitness  
          
          if fmin<fmin0    
             gbest=pbest(index,:);
             fmin0=fmin;      
          end       
           
          % calculating tolerance 
          
         % if ite>100;       
         %     tolerance=abs(ffmin(ite-100,run)-fmin0);
         % end
          
          
          % displaying iterative results    
          if ite==1      
              disp(sprintf('Iteration    Best particle    Objective fun')); 
          end   
          disp(sprintf('%8g  %8g          %8.4f',ite,index,fmin0));   
          ite=ite+1;     
    end   
    % pso algorithm-----------------------------------------------------end    
    xo=gbest;
    fval=fun(xo);
    xbest(run,:)=xo;
    ybest(run,1)=fun(xo);
    disp(sprintf('****************************************'));
    disp(sprintf('    RUN   fval       ObFuVa'));
    disp(sprintf('%6g %6g %8.4f %8.4f',run,fval,ybest(run,1)));
end
toc   
% Final neural network model 
disp('Final nn model is net_f') 
net_f = feedforwardnet(n); 
net_f=configure(net_f,inputs,targets); 
[a b]=min(ybest); xo=xbest(b,:);
k=0; 
for i=1:n     
  for j=1:m         
    k=k+1;      
    xi(i,j)=xo(k);
  end 
end 
for i=1:n 
    k=k+1;   
    xl(i)=xo(k); 
    xb1(i,1)=xo(k+n);
end 
for i=1:o   
    k=k+1;
    xb2(i,1)=xo(k);
end 
net_f.iw{1,1}=xi; 
net_f.lw{2,1}=xl;
net_f.b{1,1}=xb1; 
net_f.b{2,1}=xb2;   
%Calculation of MSE 
err = sum((net_f(inputs)-targets).^2)/length(net_f(inputs)) 
%Regression plot 
scatter(targets,net_f(inputs)) 
 
disp('Trained ANN net_f is ready for the use'); 
%Trained ANN net_f is ready for the use
plot(ybest(run,1));

