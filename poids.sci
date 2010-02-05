//Fonction poids
function [pds,dpds]=poids(x,hmin,hmax,tpefct)
  //type de fonction poids
  app=tpefct;
  
 // disp('hmin')
    //disp(hmin)
   // disp('hmax')
   // disp(hmax)
   // disp('x')
    //disp(x)
    
  select app;
    
  case 'constante';
    //disp('constante')
    if((x>=hmin) & (x<=hmax))
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

    case 'spline quadratique'
    //disp('spline quadratique')
    if((x>=hmin) & (x<=hmax))
      hh=abs(hmax-hmin);
      xx=x/hh;
        pds=1-6*xx^2+8*xx^3-3*xx^4;
        dpds=-12*xx+24*xx^2-12*xx^3; 
    else
      pds=0;
      dpds=0;
    end
    
   case 'gaussienne';
    //disp('gaussienne')
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
