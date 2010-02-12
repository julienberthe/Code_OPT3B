%calcul de la matrice PPt

function [vecp,matpp]=ppt(xx,degm,enri,nbfct)
if(enri==0)
 vecp=zeros(degm+1,1);
else
      vecp=zeros(degm+nbfct+1,1);
end
  for jj=0:degm;
    vecp(jj+1)=xx^jj;
  end
  
  if(enri~=0)
      vecp(degm+2)=cos(1/0.04*xx);
      %vecp(degm+3)=sin(1/0.04*xx);
  end
  
 matpp=vecp*vecp';
end
    
    
