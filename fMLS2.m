function [phiMLS,dphiMLS]=fMLS2(xg,xp,h,mp,dm,tpefct)


%%construction de la matrice A
A=zeros(mp+1);
Ad=zeros(mp+1);
for ii=1:length(xp)
    %calcul de la matrice d'approximation
    [vecp,matpp]=ppt(xp(ii),mp);
    %Evalution de la fonction poids
    [pds,dpds]=poids(xg-xp(ii),dm,h,tpefct);
    %Construction de A et Ad
    A=A+pds*matpp;
    Ad=Ad+dpds*matpp;
end
    
%%construction de la matrice B
B=zeros(mp+1,length(xp));
for ii=1:length(xp)
   %calcul du vecteur d'approximation
   [vecp,matpp]=ppt(xp(ii),mp);
   %evaluation de la fonction poids
   [pds,dpds]=poids(xg-xp(ii),dm,h,tpefct);
   %construction de B et Bd
   
   B(:,ii)=pds*vecp;
   
   Bd(:,ii)=dpds*vecp;
end

%%évaluation de la fonction polynomiale aux points de Gauss et de sa
%%dérivée
[vecp,matpp]=ppt(xg,mp);
[vecpd,matppd]=dppt(xg,mp);


%%construction des phi et dphi
%vecp
phiMLS=vecp'*inv(A)*B;


invAd=-inv(A)*Ad*inv(A);
dphiMLS=vecpd'*inv(A)*B+vecp'*inv(A)*Bd+vecp'*invAd*B;

end