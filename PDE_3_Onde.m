clear
clc

%%DISCRETIZZAZIONE%%
t0=0;
h=0.1;     %spacesteps
k=0.1;     %timesteps
tn=300;
t=t0:k:tn;
Nt= length(t);
xn=5;
x1=-xn;
x=x1:h:xn;
Nx =length(x);

%%STABILITA%%
lambda=k/h;
alpha=0.50;
c=alpha/lambda;
disp(['velocità=' num2str(c)]);
disp(['alpha=' num2str(alpha)]);
%%INIZIALIZZAZIONE%%
u=zeros(Nx,Nt);
v=zeros(Nx,Nt);
u0=zeros(Nx,1);
udot=zeros(Nx,1);  

%%PARAMETRI PER ERRORE%%
T=(xn-x1)/c;  %tempo per tornare al centro
tau=fix((tn-t0)/T);  %non prendo la condizione iniziale 
errore=zeros(tau,1);
n_star=zeros(tau,1);

%%CONDIZIONI INIZIALI%%
initial=input('1-gaussiana; 2-u0 nullo; 3-una componente ');
   switch(initial)
   case 1
            for i=1:Nx
                    u0(i)=exp(-x(i)^2);
                    udot(i)=0;
            end

            %ANALITICA%
            for j=1:Nt
                for i=1:Nx
                    u(i,j)=0.5*(exp(-(x(i)+c*t(j))^2)+exp(-(x(i)-c*t(j))^2));
                end
            end
   case 2 
            for i=1:Nx
                    u0(i)=0;
                    udot(i)=-2*c*x(i)*exp(-x(i)^2);
            end

            %ANALITICA%
            for j=1:Nt
                for i=1:Nx
                    u(i,j)=0.5*(exp(-(x(i)+c*t(j))^2)-exp(-(x(i)-c*t(j))^2));
                end
            end
            
   case 3 
            for i=1:Nx
                    udot(i)=2*c*x(i)*exp(-x(i)^2);
                    u0(i)=exp(-x(i)^2);
            end

            %ANALITICA%
            for j=1:Nt
                for i=1:Nx
                    u(i,j)=exp(-(x(i)-c*t(j))^2);
                end
            end
    
   end 
   
v(:,1)=u0(:);

%%SOLUZIONE NUMERICA + CONDIZIONI AL CONTORNO%%
boundary=input('1-trasparente; 2-Dirichlet; 3-Von Neumann');
for j=2:Nt
    if j==2
        v(:,j)=u0(:)+k*udot(:);
    else
        for i=2:Nx-1 %SOLUZIONE NUMERICA%
            v(i,j)=2*v(i,j-1)-v(i,j-2)+(c*lambda)*(c*lambda)*(v(i+1,j-1)-2*v(i,j-1)+v(i-1,j-1));
        end
    end
    
    switch(boundary)
        case 1 %trasparenza
            v(1,j)=v(1,j-1)+c*lambda*(v(2,j-1)-v(1,j-1));
            v(Nx,j)=v(Nx,j-1)-c*lambda*(v(Nx,j-1)-v(Nx-1,j-1));
        case 2
            v(1,j)=0;
            v(Nx,j)=0;
        case 3
            v(1,j)=v(2,j);
            v(Nx,j)=v(Nx-1,j);
    end
end

%%ERRORE%%
switch (boundary)
case 1 %%TRASPARENZA%% 
        for j=1:Nt
                AO=v(:,j).^2 ;
                sum_O=sum(AO,1);
                errore(j)=sqrt(sum_O)/Nt; 
        end
case {2,3} 
    for s=1:tau
        n_star(s)=fix(s*T/h +1); %numero istanti di ritorno al centro
        switch (initial)
        case 2 %%caso u0=0%% 
            AO=v(:,n_star(s)).^2 ;
            sum_O=sum(AO,1);
            errore(s)=sqrt(sum_O);
        case {1,3}
            AO=(u0(:)-abs(v(:,n_star(s)))).^2;
            B=u0(:).^2;
            sum_B =sum(B,1);
            sum_O=sum(AO,1);
            errore(s)=sqrt(sum_O./sum_B);
        end
    end
end

%PLOT%%
% figure(1)
% for j=1:Nt/15
%     plot(x,u0,x,v(:,j),'-o')
%     legend ("Condizione Iniziale","Numerica")
%     ylim([-1.2 1.2])
%     xlabel("x")
%     ylabel("u(x)")
%     title("t=", num2str(t(j)))
%     pause(0.01);
% end

%PLOT 3D%%
figure(2)
surf(x,t,v','FaceAlpha',0.5,"LineStyle", "None")
title('Soluzione Numerica')
grid on
colorbar
xlabel('X ')
ylabel('t ')

%PLOT ERRORE%%
figure (3)
switch (boundary)
    case 1
        semilogy(t,errore)
    case {2,3}
        semilogy(t(n_star),errore,'-*')
end
title("RMSE")
xlabel('t')
ylabel('RMSE(t)')
xlim ([0,150])
grid on

