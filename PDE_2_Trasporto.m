clear
clc

%%%%%%%%%%%%%%%%%%%%%%
%TRASPORTO IMPLICITE
%%%%%%%%%%%%%%%%%%%%%%

%%DISCRETIZZAZIONE%%
t0=0;
h=0.1;     %spacesteps
k=0.1;     %timesteps
t=t0:k:10;
Nt= length(t);
xn=8;
x1=-xn;
x=x1:h:xn;
Nx =length(x);

%%PARAMETRI SOLUZIONE%%
c=-0.5;    %velocità analitica
L=2;      %larghezza
U=2;      %ampiezza

lam=k/h; 
alpha=c*lam;
disp(['Cond di stabilità: alpha=' num2str(abs(c)*lam)]);

%%INIZIALIZZAZIONE%%
vUP=zeros(Nx,Nt);
vEU=zeros(Nx,Nt);
vLW=zeros(Nx,Nt);
vLF=zeros(Nx,Nt);
u0=zeros(Nx,1);
u=zeros(Nx,Nt);
v_BE=zeros(Nx,Nt);
v_CN=zeros(Nx,Nt);

%%CONDIZIONE INIZIALE%%
a=input('1-unipolare; 2-bipolare ');
switch (a)
    case 1 %iniziale unipolare
        for i=1:Nx
            if abs(x(i))<= L/2
                u0(i)=U/2*(1+cos(2*pi*x(i)/L));
            else 
                u0(i)= 0;
            end 
        end
        for j= 1:Nt %analitica unipolare
            for i=1:Nx
                X=x(i)+c*t(j);
                if abs(X)<= L/2  
                    u(i,j)=(U/2)*(1+cos(2*pi*X/L));
                else
                    u(i,j)= 0;
                end
            end
        end

    case 2 %iniziale bipolare
        for i=1:Nx
            if abs(x(i))<= L/2
                u0(i)=U*(sin(2*pi*x(i)/L));
            else 
                u0(i)= 0;
            end
        end
        for j= 1:Nt %analitica bipolare
            for i=1:Nx
                X=x(i)+c*t(j);
              if abs(x(i))<= L/2  
                 u(i,j)=U*(sin(2*pi*X/L));
              else 
                  u(i,j)= 0;
              end
            end
        end
end

%%INSERISCO COND INIZIALE PER OGNI METODO NUMERICO%
vLF(:,1)=u0;
vLW(:,1)=u0;
vEU(:,1)=u0;
vUP(:,1)=u0;
v_BE(:,1)= u0;
v_CN(:,1)= u0;

%%UPWIND and REVERSE UPWIND I ordine%%
%%reverse upwind non è mai stabile, upwind è stabile condizionatamente per alpha<1%%
%%sottostima fase th per alpha < 0.5%%
%reverse usa diff prime all'indietro, upwind in avanti%
if c>0
    for j= 2:Nt
            for i= 1:Nx-1
                vUP(i,j)= vUP(i,j-1)+c*lam*(vUP(i+1,j-1)-vUP(i,j-1));
            end
    vUP(Nx,j)=0;  %estremo superiore del dominio spaziale
    end
else %REVERSE UPWIND% INSTABILE
        for j= 2:Nt
            for i= 2:Nx
                   vUP(i,j)= vUP(i,j-1)+c*lam*(vUP(i,j-1)-vUP(i-1,j-1));
            end
        vUP(1,j)=0;  %estremo inferiore del dominio spaziale
        end 
end

%%FORWARD EULERO%% INSTABILE
%%instabile, usa diff bilanciate a tn-> cfl errore%%
for j= 2:Nt
        for i= 2:Nx-1         
               vEU(i,j)= vEU(i,j-1)+0.5*c*lam*(vEU(i+1,j-1)-vEU(i-1,j-1));
        end
            if c>0   
                vEU(1,j)= vEU(1,j-1)+c*lam*(vEU(2,j-1)-vEU(1, j-1));
                vEU(Nx,j)=0;  %estremo superiore del dominio spaziale
            else
                vEU(1,j)=0;   %estremo inferiore del dominio spaziale
                vEU(Nx,j)=vEU(Nx,j-1)+c*lam*(vEU(Nx,j-1)-vEU(Nx-1,j-1));  
            end 
end

%%LAX WENDROFF II ordine%% SOVRASTIMA FASE ONDA MONOCROMATICA
%differenze prime con effetto diffusivo, considero la media fra due passi intermedi%
%stabile con alpha<1%
for j= 2:Nt
        for i= 2:Nx-1     
          vLW(i,j)= vLW(i,j-1)+c*lam*0.5*(vLW(i+1,j-1)-vLW(i-1,j-1))...
              +((c*lam)^2)*0.5*(vLW(i+1,j-1)-2*vLW(i,j-1)+vLW(i-1,j-1));
        end
            if c>0  
                vLW(1,j)= vLW(1,j-1)+c*lam*(vLW(2,j-1)-vLW(1, j-1));
                vLW(Nx,j)=0;  %estremo superiore del dominio spaziale 
            else
                vLW(1,j)=0;   %estremo inferiore del dominio spaziale
                vLW(Nx,j)=vLW(Nx,j-1)+c*lam*(vLW(Nx,j-1)-vLW(Nx-1,j-1));  
            end 
end

%%LEAP FROG (III ordine, instabile)%% SOVRASTIMA PER V+ E VICEVERSA PER V-
%differenze prime bilanciate, ha due fattori d'evoluzione%
%se discriminante >0 ho complex coniug di modulo 1, condizionatamente stabile%
%se discriminante <0 ho numero reale-> non è mai stabile% 
vLF(:,2)=vUP(:,2);  %riempio seconda entrata con Upwind
 for j= 3:Nt
        for i= 2:Nx-1
                  vLF(i,j)= vLF(i,j-2)+c*lam*(vLF(i+1,j-1)- vLF(i-1,j-1));
        end
            if c>0    
                vLF(1,j)= vLF(1,j-1)+c*lam*(vLF(2,j-1)-vLF(1, j-1));
                vLF(Nx,j)=0;  %estremo superiore del dominio spaziale
            else
                vLF(1,j)=0;   %estremo inferiore del dominio spaziale
                vLF(Nx,j)=vLF(Nx,j-1)+c*lam*(vLF(Nx,j-1)-vLF(Nx-1,j-1));
            end 
 end

%%ERRORE ESPLICITI%%
errore_EU=zeros(Nt,1);
errore_UP=zeros(Nt,1); %PEGGIORE
errore_LF=zeros(Nt,1); %MEDIO
errore_LW=zeros(Nt,1); %MIGLIORE

 for j=2:Nt
    AEU=(u(:,j)-vEU(:,j)).^2 ;
    ALF=(u(:,j)-vLF(:,j)).^2 ;
    ALW=(u(:,j)-vLW(:,j)).^2 ;
    AUP=(u(:,j)-vUP(:,j)).^2 ;
    
    B=u(:,j).^2;
    sum_B =sum(B,1) ;
    
    sum_EU=sum(AEU,1);
    sum_UP=sum(AUP,1);
    sum_LF=sum(ALF,1);
    sum_LW=sum(ALW,1);
       
    errore_EU(j)=sqrt(sum_EU./sum_B);
    errore_UP(j)=sqrt(sum_UP./sum_B);
    errore_LF(j)=sqrt(sum_LF./sum_B);
    errore_LW(j)=sqrt(sum_LW./sum_B);

 end   

%le matrici contengono le condizioni al contorno con metodo upwind%

%%BACKWARD EULERO I ORDINE%% SOVRASTIMA FASE ONDA MONOCROMATICA
%diff prime bilanciate a tn+1, sempre stabile, primo ordine%
%%creo la matrice A che andrò a invertire%%
v1 = ones(Nx,1); 
v2 = zeros(Nx,1);
v3 = zeros(Nx,1);
v2(:)=c*lam/2;  %sotto
v3(:)=-c*lam/2; %sopra
v3(2)=0;
v3(Nx)=0;
v2(Nx-1)=0;

A = spdiags([v2 v1 v3],-1:1,Nx,Nx);
matEulero = full(A);

%matrice B%
b1 = zeros(Nx,1); 
b2 = zeros(Nx,1);
b3 = zeros(Nx,1);

if c>0
    b1(:)=1;
    b2(:)=0; %sotto
    b3(:)=0; %sopra
    b3(2)=c*lam;
    b1(1)=1-c*lam;
    b1(Nx)=0;
    S = spdiags([b2 b1 b3],-1:1,Nx,Nx);
    B = full(S);
else
    b1(:)=1;
    b2(:)=0; %sotto
    b3(:)=0; %sopra
    b2(Nx-1)=-c*lam;
    b1(1)=0;
    b1(Nx)=1+c*lam;
    S = spdiags([b2 b1 b3],-1:1,Nx,Nx);
    B = full(S);
end

D=matEulero\B; %fa anche l'inverso
for j= 2:Nt 
               v_BE(:,j)= D*v_BE(:,j-1);
end

%%CRANK NICOLSON%% SECONDO ORDINE, SOVRASTIMA FASE 
%media di forward and backward Eulero, sempre stabile,%
CN1 = ones(Nx,1); 
CN2 = zeros(Nx,1);
CN3 = zeros(Nx,1);
CN2(:)=alpha/4;  %sotto
CN3(:)=-alpha/4; %sopra
CN3(2)=0;
CN3(Nx)=0;

K = spdiags([CN2 CN1 CN3],-1:1,Nx,Nx);
matCN = full(K);

bb1 = ones(Nx,1); 
bb2 = zeros(Nx,1);
bb3 = zeros(Nx,1);
bb2(:)=-alpha/4; %sotto
bb3(:)=alpha/4; %sopra
bb1(1)=1-alpha;
bb2(2)=alpha;
bb1(Nx)=0;
bb3(Nx)=0;

J = spdiags([bb2 bb1 bb3],-1:1,Nx,Nx);
B_CN = full(J);

E=(matCN)\B_CN;
for j= 2:Nt 
               v_CN(:,j)= E*v_CN(:,j-1);
               if c>0
                   v_CN(1,j)=v_CN(1,j-1)+c*lam*(v_CN(2,j-1)-v_CN(1,j-1));
                   v_CN(Nx,j)=0;
               else
                   v_CN(Nx,j)=v_CN(Nx,j-1)+c*lam*(v_CN(Nx,j-1)-v_CN(Nx-1,j-1));
                   v_CN(1,j)=0;
               end
end

%%ERRORE IMPLICITI%%
errore_BE=zeros(Nt,1);
errore_CN=zeros(Nt,1); %MIGLIORE 

 for j=1:Nt
    ABE=(u(:,j)-v_BE(:,j)).^2 ;
    ACN=(u(:,j)-v_CN(:,j)).^2 ;
    
    B=u(:,j).^2;
    sum_B =sum(B,1) ;
    
    sum_BE=sum(ABE,1);
    sum_CN=sum(ACN,1);
       
    errore_BE(j)=sqrt(sum_BE./sum_B);
    errore_CN(j)=sqrt(sum_CN./sum_B);
 end 
 
%%PLOT%%
%%upwind and backward Eulero hanno ampiezza minore: sono i meno accurati
%%cran nicolson e leap frog hanno code pronunciate 
%%upwind e backward eulero perdono d'ampiezza (poco accurati)
figure(1)
for j=1:Nt/2
     plot (x, u0,x, v_BE(:,j),"-.",x,v_CN(:,j),'--',x,vUP(:,j),x,vLF(:,j),x, vLW(:,j));
     legend ('Condizione iniziale',"Backward Eulero","Crank Nicolson","Upwind","Leap Frog", "Lax Wendroff")
     grid on
     ylim([-2.5 2.5])
     title("t=", num2str(t(j)))
     xlabel('t')
     ylabel('u(x,t)')
     pause(0.1); 
     
end

%%PLOT ERRORI%% 
%Lax e Wendroff è migliore,poi Leap Frog & Crank Nicolson
figure (2)
semilogy(t,errore_UP,t,errore_LF,t,errore_LW,t,errore_BE,t,errore_CN)
legend ("Upwind","Leap Frog", "Lax Wendroff",'Backward Eulero','Crank-Nicolson')
title("RMSE")
xlabel('t')
ylabel('RMSE(t)')
xlim ([0,5])
grid on




