clear all;  
close all;
%===================================================
% Cas 1D Lineaire                                 %
% Calcul et Tracer des Fonctions de Forme EF      %
%===================================================
% Nombre de Particules : Discretisation
N=10; h=1/N; xp = [0.0:h:1.0]; 
nnodes = length(xp);
r = 2.2*h ;      %D�finition du rayon d'influence

% =====================
% Tracer des fonctions de Formes
% =====================
clear xe; he=h/10;
xe = [0.0:he:1.0];
neval=length(xe)
% Fonction de Formes
% ==================
Forme=zeros(nnodes,neval);
DForme=zeros(nnodes,neval);
for j = 1:neval
   xg  = xe(j); 
   [phi,dphi] = fMLS(xg,xp,he,r);
   Forme(:,j)=phi;
   DForme(:,j)=dphi; 
end
y=0.01*ones(1,nnodes);
% figure 
% hold
% plot(xp,y);
% for i=1:length(xp)
%   plot(xe,Forme(i,:));
% end

figure
hold
plot(xe,Forme(2,:));
plot(xe,Forme(5,:));
plot(xe,DForme(2,:));
plot(xe,DForme(5,:));
title 'Fonction de forme'
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
    figure
    hold
	plot(xe,sol,'r');
	plot(xp,u,'o');
end

