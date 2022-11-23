clear
clc

%%%%%%%%%%%%%%%%%%%%%%
%DIFFUSIONE COSENO e SENO
%%%%%%%%%%%%%%%%%%%%%%

%parametri%
L=2;    %larghezza
U=2;   %ampiezza

%vettore dei punti discreti
t0=0;
h=0.05;     %spacesteps
k=0.01;     %timesteps
t=t0:k:10;
Nt= length(t);
xn=L*0.5;
x1=-xn;
x=x1:h:xn;
Nx =length(x);

sigma=k/(h^2); 
alpha=0.2;
b=sqrt(alpha/sigma);    %radice coefficiente diffusione
disp(['Coefficiente di diffusione: b=' num2str(b)]);

u0=zeros(Nx,1);
u=zeros(Nx,Nt);

%%PERTURBAZIONE INIZIALE SENO O COSENO%%
a=input('1-coseno; 2-seno');
switch (a)
    case 1
        for i=1:Nx
            u0(i)=U*cos(pi*x(i)/L); 
        end
        for j= 1:Nt
            for i=1:Nx
                EXP=-t(j)*(pi*b/L)^2;
                u(i,j)=U*exp(EXP)*cos(pi*x(i)/L); 
            end
        end
    case 2
        for i=1:Nx
            u0(i)=U*sin(2*pi*x(i)/L); 
        end  
        for j= 1:Nt
            for i=1:Nx
                EXP=-4*pi*pi*b*b*t(j)/(L*L);
                u(i,j)=U*exp(EXP)*sin(2*pi*x(i)/L); 
            end
        end
end

%%EULERO%%
%stabile se alpha<=0.5 
%sottostima se 1/6<alpha<0.5
vEU(:,1)=u0;
vEU(1,:)=0;
vEU(Nx,:)=0;
for j= 2:Nt
        
        for i= 2:Nx-1
            
               vEU(i,j)= vEU(i,j-1)+(b*b*sigma)*(vEU(i+1,j-1)-...
                   2*vEU(i,j-1)+vEU(i-1,j-1));
        end
end

%%LEAP FROG%% INSTABILE%%
%differenze prime bilanciate per d/dt al primo membro della PDE%
vLF(:,1)=u0;
vLF(1,:)=0;
vLF(Nx,:)=0;
vLF(:,2)=vEU(:,2);
 for j= 3:Nt
        for i= 2:Nx-1

                  vLF(i,j)= vLF(i,j-2)+2*(b*b*sigma)*(vLF(i+1,j-1)...
                      -2*vLF(i,j-1)+ vLF(i-1,j-1));
        end
 end
 
%%BACKWARD EULERO%%
%incondizionatamente stabile, terzo ordine, matrici con BC nulle%
v1 = zeros(Nx,1); 
v2 = zeros(Nx,1); %sotto
v3 = zeros(Nx,1); %sopra
v1(:)=1+2*b*b*sigma;
v2(:)=-b*b*sigma;  %sotto
v3(:)=-b*b*sigma; %sopra
v1(1)=1;
v1(Nx)=1;
v3(2)=0;
v3(Nx)=0;
v2(Nx-1)=0;

A = spdiags([v2 v1 v3],-1:1,Nx,Nx);
matEulero = full(A);

B_BE=eye(Nx,Nx);
B_BE(1,1)=0;
B_BE(Nx,Nx)=0;
%%RISOLVO%%
vBE=zeros(Nx,Nt);
vBE(:,1)= u0;
D=matEulero\B_BE;
for j= 2:Nt 
       
               vBE(:,j)= D*vBE(:,j-1);    
end

%%CRANK-NICOLSON%%
%media di forward and backward Eulero, sempre stabile, terzo ordine%
g1 = zeros(Nx,1); 
g2 = zeros(Nx,1); %sotto
g3 = zeros(Nx,1); %sopra
g1(:)=1+alpha;
g2(:)=-alpha*0.5;  %sotto
g3(:)=-alpha*0.5; %sopra
g1(1)=1;
g1(Nx)=1;
g3(2)=0;
g3(Nx)=0;
g2(Nx-1)=0;

G = spdiags([g2 g1 g3],-1:1,Nx,Nx);
matCrank = full(G);

o1 = zeros(Nx,1); 
o2 = zeros(Nx,1); %sotto
o3 = zeros(Nx,1); %sopra
o1(:)=1-alpha;
o2(:)=alpha*0.5;  %sotto
o3(:)=alpha*0.5; %sopra
o1(1)=1;
o1(Nx)=1;
o3(2)=0;
o3(Nx)=0;
o2(Nx-1)=0;

O = spdiags([o2 o1 o3],-1:1,Nx,Nx);
B_CN = full(O);
%%RISOLVO%%
vCN=zeros(Nx,Nt);
vCN(:,1)= u0;
Q=matCrank\B_CN;
for j= 2:Nt 
               vCN(:,j)= Q*vCN(:,j-1);    
end
    
%%ERRORE%%
errore_EU=zeros(Nt,1); %migliore
errore_BE=zeros(Nt,1); %peggiore
errore_LF=zeros(Nt,1);
errore_CN=zeros(Nt,1); %medio
 for j=3:Nt
    AEU=(u(:,j)-vEU(:,j)).^2 ; %differenze quadratiche
    ABE=(u(:,j)-vBE(:,j)).^2 ;
    ALF=(u(:,j)-vLF(:,j)).^2 ;
    ACN=(u(:,j)-vCN(:,j)).^2 ;
    
    B=u(:,j).^2;   %denominatore
    sum_B =sum(B,1) ;
    
    sum_EU=sum(AEU,1); %calcolo somma delle differenze
    sum_BE=sum(ABE,1);
    sum_LF=sum(ALF,1);
    sum_CN=sum(ACN,1);
    
    errore_EU(j)=sqrt(sum_EU./sum_B);
    errore_BE(j)=sqrt(sum_BE./sum_B);
    errore_LF(j)=sqrt(sum_LF./sum_B);
    errore_CN(j)=sqrt(sum_CN./sum_B);   
 end    
 
%PLOT%
for j=1:Nt/15
     plot (x, vEU(:,j),x, vBE(:,j),'-.',x, vCN(:,j),'--',x,u0,'-g',x, vLF(:,j), '--');
     legend ("Eulero", "Backward Eulero", "Crank-Nicolson", "Condizione Iniziale", "Leap Frog")
     ylim ([-2.3 2.3]);
     title("t=", num2str(t(j)))
     pause(0.1);
end  
 
%PLOT DELL'ERRORE% 
%eulero è il migliore, backward eulero il peggiore
figure (1)
semilogy(t,errore_EU,t,errore_CN,t, errore_BE, t, errore_LF)
legend ("Eulero", "Crank Nicolson", "Backward Eulero", "Leap Frog")
title("RMSE")
xlabel('t')
ylabel('RMSE(t)')
ylim([0,1])
grid on

%PLOT DIFFERENZA ANALITICA E NUMERICA%% 
figure (2)
surf(t,x,u-vBE,'FaceAlpha',0.5,"LineStyle", "None")
xlabel('t')
ylabel('x')
colorbar
title("ERRORE BACKWARD EULERO 3D")
figure (3)
surf(t,x,u-vCN,'FaceAlpha',0.5,"LineStyle", "None")
xlabel('t')
ylabel('x')
colorbar
title("ERRORE CRANK NICOLSON 3D")
figure (4)
surf(t,x,u-vEU,'FaceAlpha',0.5,"LineStyle", "None")
xlabel('t')
ylabel('x')
colorbar
title("ERRORE EULERO 3D")
figure (5)
surf(t,x,u-vLF,'FaceAlpha',0.5,"LineStyle", "None")
xlabel('t')
ylabel('x')
colorbar
title("ERRORE LEAP FROG 3D")



