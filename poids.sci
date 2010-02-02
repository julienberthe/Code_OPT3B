//Fonction poids
function [pds,dpds]=poids(x,hmin,hmax)
  
  app='constante';
  
  select app;
    
  case 'constante';
    //disp('constante')
    if((x<=-hmin) & (x>=-hmax))
        pds=1;
        dpds=0; 
    else
      pds=0;
      dpds=0;
    end
    
//  case 'gaussienne'
//    if((abs(x)/h)<=1)
//        pds=(1/(%pi*(h^2)))*exp((-1)*(x^2));
//        dpds=0; 
//    else
//      pds=0;
//      dpds=0;
//    end
    
   case 'gaussienne' 
    if((x<=-hmin) & (x>=-hmax))
        pds=(1/(%pi*((hmax-hmin)^2)))*exp((-1)*(x^2));
        dpds=0; 
    else
      pds=0;
      dpds=0;
    end
     
  
case 'harmonique';
   if((x<=-hmin) & (x>=-hmax))
        pds=cos(%pi*x);
        dpds=-sin(%pi*x); 
    else
      pds=0;
      dpds=0;
    end
end


  
endfunction
