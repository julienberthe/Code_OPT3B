%calcul de la matrice PPt

function [vecpd,matppd]=dppt(xx,degm,enri,nbfct)
if(enri==0)
  vecpd=zeros(degm+1,1);
else
        vecpd=zeros(degm+nbfct+1,1);
end
  for iii=0:degm
    if(iii==0);
      vecpd(iii+1)=0;
    else
     vecpd(iii+1)=xx^(iii-1);
    end   
  end
  if(enri~=0)
      vecpd(degm+2)=-sin(1/0.04*xx);
      %vecpd(degm+3)=cos(1/0.04*xx);
  end
      
 matppd=vecpd*vecpd';
 
end
    
    
