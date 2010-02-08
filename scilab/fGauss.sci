// ============================================================
// ============================================================
// Integration : Points de Gauss et Poids
// ============================================================
// ============================================================
function [gg,weight,jac] = fGauss(h,ncells)
// gg = Points Extremes + 4 Points par segments
// Sur Element de reference [-1, 1] :
h2=h/2;
jac=h2;
gg0     = [-0.8611363115; -0.3399810435; 0.3399810435; 0.8611363115];
weight0 = [ 0.3478548451;  0.6521451548; 0.6521451548; 0.3478548451];
weight=[0;weight0];
gg1=h2*(gg0+ones(4,1));
hhg=gg1(1);
gg=[0;gg1];
for ii=1:(ncells-1) weight=[weight;weight0]; end;
for ii=1:(ncells-1) gg1=gg1+h*ones(4,1) ; gg=[gg;gg1]; end;
gg=[gg ; 1.];
weight=[weight;0];
endfunction