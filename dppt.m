%calcul de la matrice PPt

function [vecpd,matppd]=dppt(xx,degm)
  vecpd=zeros(degm+1,1);
  for iii=0:degm
    if(iii==0);
      vecpd(iii+1)=0;
    else
     vecpd(iii+1)=xx^(iii-1);
   end
   
  end
 matppd=vecpd*vecpd';
 
end
    
    
