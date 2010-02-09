% ============================================================
% Fonction de Forme MLS
% ============================================================
% ============================================================
function [phiMLS,dphiMLS] = fMLS(xg,xp,h,mm,dm,tpefct)
  %degrÃ© du polynome d'approximation
  Mp=mm;  
 
  %construction de l'opÃ©rateur A
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
  
  %disp(A)
  
  %construction de l'opÃ©rteur B
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
  %disp('premier terme')
  %vecpd3'*inv(A)*B
  %phiMLS
  %vecpd3'
  %A
  %inv(A)
  %disp('deuxième terme')
  %vecp3'*inv(A)*Bd
  %Bd
  %disp('troisième terme')
  %vecp3'*Add*B
  %vecp3'*Add
  %vecpd3'*inv(A)
  %inv(A)*Ad*inv(A)
  %B
%   vecpd3'
%   inv(A)*B
%   B
%   inv(A)
%   cond(A)
% if(mm>1)
%   for j=2:(length(phiMLS)-1)
%       %dphiMLS(1,j)=(phiMLS(1,j+1)+phiMLS(1,j-1)-2*phiMLS(1,j))/(2*h);
%     dphiMLS(1,j)=(phiMLS(1,j+1)-phiMLS(1,j))/(h);
%   end
%   %dphiMLS(1,1)=(phiMLS(1,2)-phiMLS(1,1))/h;
%   dphiMLS(1,length(phiMLS))=(phiMLS(1,length(phiMLS))-phiMLS(1,length(phiMLS)-1))/h;
% else
  dphiMLS=vecpd3'*inv(A)*B+vecp3'*inv(A)*Bd+vecp3'*Add*B;%vecpd3'*inv(A)*B+
  %dphiMLS=zeros(length(xp),1);
   %Bd
   %B
%end
%dphiMLS
%h

end
