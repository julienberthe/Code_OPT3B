//Fonction poids
function [pds,dpds]=poids(x,dm,h,tpefct)
  //type de fonction poids
  app=tpefct;
  
  //support de la fonction poids
  hmls=dm*h;
  
  //distance relative
  s=abs(x)/hmls;
  ss=x/hmls;

    
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
    
    case 'spline quadratique';
    //disp('spline quadratique')
    if(s<=1)      
        pds=1-6*s^2+8*s^3-3*s^4;
        if (ss>=0)
          dpds=-12*ss+24*ss^2-12*ss^3; 
        else
          dpds=-12*ss-24*ss^2-12*ss^3;
        end       
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
  disp('harmonique')
   if(s<=1)
        pds=cos(%pi*s);
        dpds=-sin(%pi*s); 
    else
      pds=0;
      dpds=0;
    end
end


  
endfunction
