% ============================================================
% Fonction de Forme MLS
% ============================================================
% ============================================================
function [phiMLS,dphiMLS] = fMLS(xg,xp,h,mm,dm,tpefct)
  %degré du polynome d'approximation
  Mp=mm;  
 
  %construction de l'opérateur A
  A=zeros(Mp+1);
  Ad=zeros(Mp+1);
  for i=1:length(xp)
    [vecp1,matpp1]=ppt(xp(i),Mp);
    %disp(matpp1)
    [momo1,mama1]=poids(xg-xp(i),dm,h,tpefct);
    %disp(xg-xp(i))
    %disp(momo1);
    A=A+momo1*matpp1;
    Ad=Ad+mama1*matpp1;
    
  end
  
  
  %construction de l'opérteur B
  B=zeros(Mp+1,length(xp));
  Bd=zeros(Mp+1,length(xp));
  for i=1:length(xp)
    [vecp2,matpp2]=ppt(xp(i),Mp);
    [momo2,mama2]=poids(xg-xp(i),dm,h,tpefct);
    %disp('vecp')
    %disp(momo2)
    B(:,i)=momo2*vecp2;
    %disp(B)
    Bd(:,i)=mama2*vecp2;
   
    %disp(vecp2)
    %disp(xg-xp(i))
    %disp('et la')
end

  %disp('ici')
  %construction des fonctions de forme
  [vecp3,matpp3]=ppt(xg,Mp);
  [vecpd3,matppd3]=dppt(xg,Mp);
    %vecpd3
  %vecp3
  %xg
  phiMLS=vecp3'*inv(A)*B;
  Add=-inv(A)*Ad*inv(A);
  dphiMLS=vecpd3'*inv(A)*B+vecp3'*inv(A)*Bd+vecp3'*Add*B;
  %dphiMLS=zeros(length(xp),1);
   %Bd
   %B
   


end
