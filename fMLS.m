%============================================================
% Fonction de Forme MLS
% ============================================================
% ============================================================
function [phiMLS,dphiMLS] = fMLS(xg,xp,h,mm,dm,tpefct,enri)
  %degré du polynome d'approximation
  Mp=mm;  
 nbfct=0;
  %initialisation
  if(enri==0)
    A=zeros(Mp+1);
    Ad=zeros(Mp+1);
     B=zeros(Mp+1,length(xp));
  Bd=zeros(Mp+1,length(xp));
  else
      nbfct=1;  %nb de fonctions ajoutées dans l'enrichissement
        A=zeros(Mp+1+nbfct);
    Ad=zeros(Mp+1+nbfct);
     B=zeros(Mp+1+nbfct,length(xp));
  Bd=zeros(Mp+1+nbfct,length(xp));
  end
  
  %construction de l'opérateur A
  for i=1:length(xp)
    [vecp1,matpp1]=ppt(xp(i),Mp,enri,nbfct);
    %disp(matpp1)
    [momo1,mama1]=poids(xg-xp(i),dm,h,tpefct);
    %disp(xg-xp(i))
    %disp(momo1);
    A=A+momo1*matpp1;
    Ad=Ad+mama1*matpp1;
    
  end
  
  %disp(A)
  
  %construction de l'opérteur B
 
  for i=1:length(xp)
    [vecp2,matpp2]=ppt(xp(i),Mp,enri,nbfct);
    [momo2,mama2]=poids(xg-xp(i),dm,h,tpefct);
    %disp('vecp')
    %disp(momo2)
    B(:,i)=momo2*vecp2;
    %disp(B)
    Bd(:,i)=mama2*vecp2;

end

  %disp('ici')
  %construction des fonctions de forme
  [vecp3,matpp3]=ppt(xg,Mp,enri,nbfct);
  [vecpd3,matppd3]=dppt(xg,Mp,enri,nbfct);

  phiMLS=vecp3'*inv(A)*B;
  Add=-inv(A)*Ad*inv(A);

  dphiMLS=vecpd3'*inv(A)*B+vecp3'*inv(A)*Bd+vecp3'*Add*B;%vecpd3'*inv(A)*B+


end
