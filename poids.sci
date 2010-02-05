//Fonction poids
function [pds,dpds]=poids(x,dm,h,tpefct)
  //type de fonction poids
  app=tpefct;
  
  //support de la fonction poids
  hmls=dm*h;
  
  //distance relative
  s=abs(x)/hmls;
  
 // disp('hmin')
    //disp(hmin)
   // disp('hmax')
   // disp(hmax)
   // disp('x')
    //disp(x)
    
  select app;
    
  case 'constante';
    //disp('constante')
    if(s<=1)
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
    if(s<=1)      
        pds=1-6*s^2+8*s^3-3*s^4;
        dpds=-12*s+24*s^2-12*s^3; 
    else
      pds=0;
      dpds=0;
    end
    
   case 'gaussienne';
    //disp('gaussienne')
    if(s<=1)
        pds=(1/(%pi*(hmls^2)))*exp((-1)*(s^2));
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
