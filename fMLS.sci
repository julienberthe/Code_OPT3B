// ============================================================
// Fonction de Forme MLS
// ============================================================
// ============================================================
function [phiMLS,dphiMLS] = fMLS(xg,xp,he,mm,nn,tpefct)
  //degré du polynome d'approximation
  Mp=mm;  
  //Nombre de particules dans la zone d'influences:
  N=nn;
  
  // Construction des bornes du support des fonctions de formes
  compt=0;
  toto=modulo(100*xg,1000*he);
  if(toto>9.99)
      toto=0;
    end
  if(toto<0.01)
    compt=1;
  end
  //disp(compt)
  minus=xg;
  maxus=xg;
  hmin=0;
  hmax=0;
  while(compt<N)
    //disp('max')
    maxus=maxus+he;
    //  disp(maxus)
    caca=modulo(100*maxus,1000*he);
    //disp(caca)
    if(caca>9.99)
      caca=0;
    end
    if(caca<0.01 & maxus<=1.001)
      //disp('je suis la')
      hmax=maxus-xg+0.001;
      compt=compt+1;
    end
    if(compt<N)
    //disp('min')
    minus=minus-he;
    //disp(minus)
    coco=modulo(100*minus,1000*he);
    if(coco>9.99)
      coco=0;
    end
    //disp(coco)
    if(coco<0.01 & minus>=-0.001)
      //disp('je suis ici')
      hmin=minus-xg-0.001;
      compt=compt+1;
    end
    end
  end
  //disp(xg)
  //disp(hmin)
  //disp(hmax)
  //disp(compt)
  //Support des fonctions de formes
  
  //h=(N/2)*(xp(2)-xp(1));
//  disp('etude')
//  disp(xg)
//  disp('hmin')
//  disp(hmin)
//  disp('hmax')
//  disp(hmax)
  
  //construction de l'opérateur A
  A=zeros(Mp+1);
  Ad=zeros(Mp+1);
  for i=1:length(xp)
    [vecp1,matpp1]=ppt(xp(i),Mp);
    // disp(matpp1)
    [momo1,mama1]=poids(xp(i)-xg,hmin,hmax,tpefct);//,xg
    //disp(xg-xp(i))
    //disp(momo1);
    A=A+momo1*matpp1;
    Ad=Ad+mama1*matpp1;
    
  end
  
  //disp(A)
  
  //construction de l'opérteur B
  B=zeros(Mp+1,length(xp));
  Bd=zeros(Mp+1,length(xp));
  for i=1:length(xp)
    [vecp2,matpp2]=ppt(xp(i),Mp);
    [momo2,mama2]=poids(xp(i)-xg,hmin,hmax,tpefct);
    disp('vecp')
    disp(vecp2)
    B(:,i)=momo2*vecp2;
    //disp(B)
    Bd(:,i)=mama2*vecp2;
    //disp(xg-xp(i))
  end
  
  //construction des fonctions de forme
  [vecp3,matpp3]=ppt(xg,Mp);
  [vecpd3,matppd3]=dppt(xg,Mp);
  phiMLS=vecp3'*inv(A)*B;
  //disp(A)
  //disp(inv(A))
  //dphiMLS=vecpd3'*inv(Ad)*Bd;
  dphiMLS=zeros(length(xp),1);
   

endfunction