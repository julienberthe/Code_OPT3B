//calcul de la matrice PPt

function [vecp,matpp]=ppt(x,degm)
  vecp=zeros(degm+1,1);
  for i=0:degm;
    vecp(i+1)=x^i;
  end
 matpp=vecp*vecp';
 //disp(vecp)
 //disp(vecp')
//disp(matpp) 
 endfunction
    
    
