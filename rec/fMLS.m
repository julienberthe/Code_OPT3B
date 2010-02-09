% ============================================================
% Fonction de Forme Elements Finis 1D Lineaire
% ============================================================
% ============================================================
function [phiMLS,dphiMLS] = fMLS(xg,xp,he,r)
  % Paramètres de la fonction de depart :
  % xg : point courant de l'evalutation
  % xp : intervalle [0;1] discrétisé
  % he : dixieme du pas en espace
  % r  : nouvel argument, rayon d'influence de la fMLS consideree
h = xp(2)-xp(1);                   % parce que le pas en espace n'est pas passé en argument
rinf = r ;
phiMLS  = zeros(1,length(xp));     % declaration des fonctions phi
dphiMLS = zeros(1,length(xp));     % derive des fonctions phi
[p,pX] = build_P(xg);                   % Construit le chargement    
m = length(p);                             % Taille de la base d'interpolation 
P = zeros(length(xp),m);           % declaration de la matrice P

omega = zeros(length(xp),length(xp));
domega = zeros(length(xp),length(xp));
%D�finition des fonctions de forme


for i=1:length(xp)
    dist= abs((xg-xp(i))/rinf);
    if dist <=1
       omega(i,i)=1-3*dist^2+2*dist^3;
       % derivee
	domega(i,i) = -6*dist+6*dist^2;
    end
    P(i,:) = build_P(xp(i));
end
A = P'*omega*P;
Ax = P'*domega*P;
AX = -A^(-1)*Ax*A^(-1) ;
B = P'*omega ;
BX = P'*domega ;
phiMLS = p*(A^-1)*P'*omega;
dphiMLS = pX * A^(-1) * B + p * AX * B + p * A^(-1) * BX;
% Calculs pour la déricée de Phi
 




%  for i = 1:length(phiMLS)
%      if i==1
%          dphiMLS(1,i) = (phiMLS(2)-phiMLS(1))/h;
%      else
%          if i==length(phiMLS)
%              dphiMLS(1,i) = (phiMLS(length(phiMLS))-phiMLS(length(phiMLS)-1))/h;
%          else
%              dphiMLS(1,i) = (phiMLS(i+1)-phiMLS(i-1))/2/h;
%          end
%      end
%  end
%  end

