function v = ptvalue(x,alpha,beta, lambda)
%prospect theory value fuction 
if (x<0)
   v = lambda.*abs(x).^beta;
else
    v = x.^alpha; 
end;
