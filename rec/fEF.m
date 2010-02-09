% ============================================================
% Fonction de Forme Elements Finis 1D Lineaire
% ============================================================
% ============================================================
function [phiEF,dphiEF] = fEF(xg,xp,he)
h = xp(2)-xp(1);
phiEF  = zeros(1,length(xp));
dphiEF = zeros(1,length(xp));
for i=1:length(xp)
    dist=(xg-xp(i))/h;
		if (abs(dist)==0)
		   phiEF(i)=1;
    else
    	 if ((abs(dist))<(1-he/10))
    	    if (dist>0)
    			    dphiEF(i)=-1/h;
	 	   	   phiEF(i)=1-dist;
    	    else
    			    dphiEF(i)= 1/h;
	 	   	   phiEF(i)=1+dist;
   		    end
		 end
		end
end

