clear all;  
getf("fGauss.sci"); getf("fEF.sci");getf("fMLS.sci");
//===============================================================
// Cas 1D :                                                    //
// Resolution de l'equation :   u,xx + f = 0                   //
// ==========================  -k0*u(0)= -u,x(0) // u(0) = 0   //
//                              u,x(1)= T                      //
//                                                             //
// f= 6*d*x                                                    //                                                        //
// u = f0 + e*x - d*x^3	                                       //
//===============================================================
// Parametres de l'equation
// ========================
T    = -2;  d = 5;
k0   = 10;           k00=1; if (k0==0) k00=0; k0=1; end
// Nombre de Particules : Discretisation
N=19; h=1/N; xp = [0.0:h:1.0]; nnodes = length(xp); ncells = nnodes-1;
// Choix de la Methode: Methode Elements Finis
// ====================
MEF=0;
MLSType='constante';
mp=1;
np=2;

// Points de Gauss
// ===============
[gg,weight,jac] = fGauss(h,ncells); hhg=gg(2)-gg(1);
// ==========================
// Initialistion des Matrices : Mise a Zero
// ==========================
k  = zeros(nnodes) ; f  = zeros(nnodes,1);  GG = zeros(nnodes,1);
// Boucle sur les points de Gauss
// ==============================
for j = 1:length(gg)
   xg = gg(j);
   weight1=weight(j);
   // Calcul Phi(xg), dPhi(xg)
   if (MEF==1) [phi,dphi] = fEF(xg,xp,hhg); end;
   if (MEF==0) [phi,dphi] = fMLS(xg,xp,hhg,mp,np,MLSType); end;
   // Calcul Matrice de rigidite : k et Second Membre : f
   if j == 1
    GG(1:3,1) = -phi(1:3)';
    k = k+k00*k0*phi'*phi;
   else
	 	if j<length(gg)
	    	 k = k+(weight1*jac)*(dphi'*dphi);
        fbody=6*d*xg ;
	     	f = f+(weight1*fbody*jac)*phi';
	   end
	   if j==length(gg)
	     f= f+T*phi';
	   end
  end
end
// ==========
// Resolution
// ==========
q=[0];
if (k00*k0==0)
  Encastrement = k0*k00;
  mat = [k GG; GG' zeros(1)];
  depl = mat\[f' q]'; 
  nnodesT=length(depl)-1
else
  Ressort = k0*k00;
  mat = [k]; depl = mat\[f];
  nnodesT=length(depl)
end;
u = depl(1:nnodesT); uxp = depl(1:nnodes);
// =====================
// Tracer de la solution
// =====================
clear xe; clear sol; he=h/10;
xe = [0.0:he:1.0];
// Fonction de Formes
for j = 1:length(xe)
   xg  = xe(j); [phi,dphi] = fEF(xg,xp,he);
   for i=1:nnodesT Forme(j,i)=phi(i); end;
end
// Construction de la solution u=Sum_i u_i Forme_i
sol=zeros(length(xe));
for j=1:length(xe)
	sol(j)=0.; for i=1:nnodesT sol(j)=sol(j)+u(i)*Forme(j,i); end
end
//
plot2d(xe,sol,style=1);
plot2d(xp,uxp,style=-1);
