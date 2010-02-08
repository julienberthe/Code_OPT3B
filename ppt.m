%calcul de la matrice PPt

function [vecp,matpp]=ppt(xx,degm)
  vecp=zeros(degm+1,1);
  for i=0:degm;
    vecp(i+1)=xx^i;
  end
  
 matpp=vecp*vecp';
end
    
    
