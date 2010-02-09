%calcul de la matrice PPt

function [vecp,matpp]=ppt(xx,degm)
  vecp=zeros(degm+1,1);
  for jj=0:degm;
    vecp(jj+1)=xx^jj;
  end
  
 matpp=vecp*vecp';
 
end
    
    
