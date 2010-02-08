%calcul de la matrice PPt

function [vecpd,matppd]=dppt(xx,degm)
  vecpd=zeros(degm+1,1);
  for i=0:degm
    if(i==0);
      vecpd(i+1)=0;
    else
     vecpd(i+1)=xx^(i-1);
   end
   
  end
 matppd=vecpd*vecpd';
 
end
    
    
