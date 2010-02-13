clear all;clf;close all;
%getf("fGauss.sci"); getf("fEF.sci");getf("fMLS.sci");
%===============================================================
% Cas 1D :                                                    %
% Resolution de l'equation :   u,xx + f = 0                   %
% ==========================  -k0*u(0)= -u,x(0) % u(0) = 0   %
%                              u,x(1)= T                      %
%                                                             %
% f= 6*d*x                                                    %                                                        %
% u = f0 + e*x - d*x^3	                                       %
%===============================================================
% Parametres de l'equation
% ========================
T    = -2;   d = 30; L=1; b=3;a=-1;c=1/0.04;
k0   = 10;           k00=1; if (k0==0) k00=0; k0=1; end
% Nombre de Particules : Discretisation
N=50; h=1/N; xp = [0.0:h:1.0]; nnodes = length(xp); ncells = nnodes-1;
he=h/10;
% Choix de la Methode: Methode Elements Finis
% ====================
MEFvar=0;
MLSType='spline quadratique';
DER=0; %calcul des dérivées exactes (1, seulement pour mp=1) ou par différences centrées (0)
enri=1;  %enrichissement (1) ou non (0)
PUM=1;  %PUM si 1 rien si 0
mp=3;
dm=5.1;

% Points de Gauss
% ===============
[gg,weight,jac] = fGauss(h,ncells); hhg=gg(2)-gg(1);
% ==========================
% Initialistion des Matrices : Mise a Zero
% ==========================
if(PUM==0)
k  = zeros(nnodes) ; f  = zeros(nnodes,1);  GG = zeros(nnodes,1);
else
k  = zeros(2*nnodes) ; f  = zeros(2*nnodes,1);  GG = zeros(2*nnodes,1);
end
% Boucle sur les points de Gauss
% ==============================
if(MEFvar==1 || DER==1)    
for j = 1:length(gg)
   xg = gg(j);
   weight1=weight(j);
   % Calcul Phi(xg), dPhi(xg)
   if (MEFvar==1) [phi,dphi] = fEF(xg,xp,hhg); end;
   if (MEFvar==0&&DER==1) [phi,dphi] = fMLS(xg,xp,h,mp,dm,MLSType,enri); end;
   for i=1:nnodes 
            %disp(phi)
            Forme(i,j)=phi(i);
        end;
   % Calcul Matrice de rigidite : k et Second Membre : f
   if(PUM==0)
   if j == 1
    GG(1:3,1) = -phi(1:3)';
    k(1:nnodes,1:nnodes) = k(1:nnodes,1:nnodes)+k00*k0*phi'*phi;
   else
	 	if j<length(gg)
	    	 k(1:nnodes,1:nnodes) = k(1:nnodes,1:nnodes)+(weight1*jac)*(dphi'*dphi);
        %fbody=6*d*xg ;
        fbody= 6*d*xg+b*c*c*cos(c*(xg-a));
	     	f(1:nnodes,1) = f(1:nnodes,1)+(weight1*fbody*jac)*phi';
	   end
	   if j==length(gg)
	     f(1:nnodes,1)= f(1:nnodes,1)+T*phi';
	   end
   end
   end
end
end

if(MEFvar==0&&DER==0)  %Calcul des fonctions de forme par différences centrées
    disp('Différences centrées')
    for j = 1:length(gg)
        xg = gg(j);
        %weight1=weight(j);
        [phi,dphi] = fMLS(xg,xp,h,mp,dm,MLSType,enri);
        for i=1:nnodes 
            %disp(phi)
            Forme(i,j)=phi(i);
        end;
    end
    for i=1:nnodes
        for j=1:length(gg)
            if(j<length(gg))
            DForme(i,j)=(Forme(i,j+1)-Forme(i,j))/(gg(j+1)-gg(j));
            %else
            %DForme(i,:)=dphi(i);
            end
            %%dérivées des fonctions de forme calculées par différence
            %%finie
        DForme(i,length(gg))=(Forme(i,length(gg))-Forme(i,length(gg)-1))/(gg(length(gg))-gg(length(gg)-1));
        end
    end;
    if(PUM==0)
    for j = 1:length(gg)
        xg = gg(j);
        weight1=weight(j);
        if j == 1
            GG(1:3,1) = -Forme(1:3,j);
            k(1:nnodes,1:nnodes) = k(1:nnodes,1:nnodes)+k00*k0*Forme(:,j)*Forme(:,j)';
        else
            if j<length(gg)
                k(1:nnodes,1:nnodes) = k(1:nnodes,1:nnodes)+(weight1*jac)*(DForme(:,j)*DForme(:,j)');
                %fbody=6*d*xg ;
                fbody= 6*d*xg+b*c*c*cos(c*(xg-a));
                f(1:nnodes,1) = f(1:nnodes,1)+(weight1*fbody*jac)*Forme(:,j);
            end
            if j==length(gg)
            f(1:nnodes,1)= f(1:nnodes,1)+T*Forme(:,j);
            end
        end
    end
    end
end

if(PUM==1)
    for i = 1:nnodes
        for j=1:length(gg)
            Forme(i+nnodes,j)=Forme(i,j)*(cos(c*(gg(j)-a)));%(1/c)
        end;
    end;
    for i=1:2*nnodes
        for j=1:length(gg)
            if(j<length(gg))
            DForme(i,j)=(Forme(i,j+1)-Forme(i,j))/(gg(j+1)-gg(j));
            %else
            %DForme(i,:)=dphi(i);
            end
            %%dérivées des fonctions de forme calculées par différence
            %%finie
        DForme(i,length(gg))=(Forme(i,length(gg))-Forme(i,length(gg)-1))/(gg(length(gg))-gg(length(gg)-1));
        end
    end;
    for j = 1:length(gg)
        xg = gg(j);
        weight1=weight(j);
        if j == 1
            GG(1:3,1) = -Forme(1:3,j);
            GG(nnodes+1:nnodes+3,1) = -Forme(nnodes+1:nnodes+3,j);
            k = k+k00*k0*Forme(:,j)*Forme(:,j)';
        else
            if j<length(gg)
                k = k+(weight1*jac)*(DForme(:,j)*DForme(:,j)');
                %fbody=6*d*xg ;
                fbody= 6*d*xg+b*c*c*cos(c*(xg-a));
                f = f+(weight1*fbody*jac)*Forme(:,j);
            end
            if j==length(gg)
            f= f+T*Forme(:,j);
            end
        end
    end
    
        
    
end
    


% ==========
% Resolution
% ==========
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
% =====================
% Tracer de la solution
% =====================
%clear Forme; %clear DForme;
clear xe; clear sol; he=h/10;
xe = [0.0:he:1.0];
% Fonction de Formes
for j = 1:length(xe)
   xg  = xe(j); 
   if (MEFvar==1) [phi,dphi] = fEF(xg,xp,hhg); end;
   if (MEFvar==0) [phi,dphi] = fMLS(xg,xp,h,mp,dm,MLSType,enri); end;
   if(PUM==0)
   for i=1:(nnodesT) Forme2(j,i)=phi(i); end;
   else
       for i=1:(nnodesT/2) Forme2(j,i)=phi(i); end;
       for i=1:(nnodesT/2) Forme2(j,i+(nnodesT/2))=Forme2(j,i)*(cos(c*(xg-a))); end;
   end
end
% Construction de la solution u=Sum_i u_i Forme_i
sol=zeros(length(xe),1);
solreel=zeros(length(xe),1);
for j=1:length(xe)
	sol(j)=0.; for i=1:nnodesT sol(j)=sol(j)+u(i)*Forme2(j,i); end
    solreel(j)=-30*xe(j)^3+3*cos(25*xe(j)+25)+(88+75*sin(50))*xe(j)-3*cos(25)-(15/2)*sin(25)+44/5+(15/2)*sin(50);
end
%
figure;
plot(xe,sol,'g');
hold on;
plot(xe,solreel,'b');
% hold on;
% plot(xe,Forme2(:,5),'r');
% hold on;
% plot(xe,Forme2(:,16),'m');


%%Erreur sur la solution
erreur=zeros(size(sol,1),1);
for ii=1:size(sol,1)
    erreur(ii)=abs(solreel(ii,1)-sol(ii,1));%/solreel(ii,1);
end

figure;
hold on
plot(xe,erreur,'r')

%moyenne de l'erreur
moy=mean(erreur);
disp(moy)
maxx=max(erreur);
disp(maxx)
somme=sum(erreur)/L;
disp(somme)

%calcul de l'énergie potentielle


%hold on 
%plot(xe,log(erreur),'r')

% hold on;
% plot(gg,DForme(5,:),'m');
% hold on;
% plot(gg,Forme(5,:),'r');
