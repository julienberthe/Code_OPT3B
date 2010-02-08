// ============================================================
// Fonction de Forme MLS
// ============================================================
// ============================================================
function [phiMLS,dphiMLS] = fMLS(xg,xp,h,mm,dm,tpefct)
  //degré du polynome d'approximation
  Mp=mm;  
 
  //construction de l'opérateur A
  A=zeros(Mp+1);
  Ad=zeros(Mp+1);
  for i=1:length(xp)
    [vecp1,matpp1]=ppt(xp(i),Mp);
    //disp(matpp1)
    [momo1,mama1]=poids(xg-xp(i),dm,h,tpefct);
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
    [momo2,mama2]=poids(xg-xp(i),dm,h,tpefct);
    //disp('vecp')
    //disp(momo2)
    B(:,i)=momo2*vecp2;
    //disp(B)
    Bd(:,i)=mama2*vecp2;
    //disp(vecp2)
    //disp(xg-xp(i))
    //disp('et la')
end

  //disp('ici')
  //construction des fonctions de forme
  [vecp3,matpp3]=ppt(xg,Mp);
  [vecpd3,matppd3]=dppt(xg,Mp);
  //disp(vecpd3)
  //disp(vecp3)
  phiMLS=vecp3'*inv(A)*B;
  Add=-inv(A)*Ad*inv(A);
  dphiMLS=vecpd3'*inv(A)*B+vecp3'*inv(A)*Bd+vecp3'*Add*B;
  //dphiMLS=zeros(length(xp),1);
   
   
   
 //    // Construction des bornes du support des fonctions de formes
//  compt=0;
//  toto=modulo(100*xg,1000*he);
//  if(toto>9.99)
//      toto=0;
//    end
//  if(toto<0.01)
//    compt=1;
//  end
//  //disp(compt)
//  minus=xg;
//  maxus=xg;
//  hmin=0;
//  hmax=0;
   
//     while(compt<N)
//      //disp('max')
//      maxus=maxus+he;
//      //  disp(maxus)
//      caca=modulo(100*maxus,1000*he);
//      //disp(caca)
//      if(caca>9.99)
//        caca=0;
//      end
//      if(caca<0.01 & maxus<=1.001)
//        //disp('je suis la')
//        hmax=maxus-xg+0.001;
//        compt=compt+1;
//      end
//      if(compt<N)
//      //disp('min')
//      minus=minus-he;
//      //disp(minus)
//      coco=modulo(100*minus,1000*he);
//      if(coco>9.99)
//        coco=0;
//      end
//      //disp(coco)
//      if(coco<0.01 & minus>=-0.001)
//        //disp('je suis ici')
//        hmin=minus-xg-0.001;
//        compt=compt+1;
//      end
//      end
//    end

endfunction
