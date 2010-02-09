clear all; clf; 
%===================================================
% Cas 1D Lineaire                                 %
% Calcul et Tracer des Fonctions de Forme EF      %
%===================================================
% Nombre de Particules : Discretisation
N=10; h=1/N; xp = [0.0:h:1.0]; 
nnodes = length(xp);
% =====================
% Tracer des fonctions de Formes
% =====================
clear xe; he=h/10;
xe = [0.0:he:1.0];
neval=length(xe);

%choix de l'approche (EF, MLS)
app='MLS';

  % Fonction de Formes
% ==================
switch app;
case 'EF';

Forme=zeros(nnodes,neval);
DForme=zeros(nnodes,neval);
for j = 1:neval
  xg  = xe(j); 
   [phi,dphi] = fEF(xg,xp,he);
   for i=1:nnodes 
       Forme(i,j)=phi(i);
       DForme(i,j)=dphi(i); 
   end;
end
y=0.01*ones(1,nnodes);
%plot2d(xp,y,style=-1);
for i=1:length(xp)
%  plot2d(xe,Forme(i,:),style=i);
end
%plot(xe,Forme(2,:));
%plot(xe,Forme(5,:));

case 'MLS';
  Forme=zeros(nnodes,neval);
  DForme=zeros(nnodes,neval);
for j=1:neval
    xg=xe(j); 
    %degré du polynôme d'approximation
    mp=2;
    
    %définition du support de la fonction poids
    dm=3.1;
    %type de fonction poids (constante, gaussienne, harmonique, spline quadratique)
    tpefct='spline quadratique';
    
   [phi,dphi] = fMLS2(xg,xp,h,mp,dm,tpefct);
   %disp(phi)

   for i=1:nnodes 
       %disp(phi)
       Forme(i,j)=phi(i);
       DForme(i,j)=dphi(i);
   end;
end
clear phi; clear dphi;
y=0.01*ones(1,nnodes);
%plot2d(xp,y,style=-1);

for i=1:length(xp)
%  plot2d(xe,Forme(i,:),style=i);
end

%plot2d(xe,Forme(2,:),style=2);
%plot2d(xe,DForme(2,:),style=3);
hold on;
plot(xe,Forme(5,:),'r');
hold on;
plot(xe,DForme(5,:),'b');

  
end

%
% Construction de la solution u(x)=Sum_I u_I phi_I
% ================================================
if (N==10) 
	u=[0. 0.4 0.5 1.4 1.5 1.45 1.4 1.2 1. 0.5 0.45];
%  u=[0. .12  0.28 .39 0.35 0.05  0.48 .52 0.64 0.65 0.5];
	sol=zeros(1,neval);
	for j=1:neval
	sol(j)=0.;
	for i=1:nnodes
	sol(j)=sol(j)+u(i)*Forme(i,j);
	end
    end
    hold on;
	plot(xe,sol);
    hold on;
	plot(xp,u);
end

