function [Q, R, Z, I] = CHASEChoiceMatrices(d, pstay, theta, tau)
%step size
%Delta = alpha.*sigma.*sqrt(tau);

%total number of states
m = 2.*theta+1; %2.*round(theta./Delta) + 1;
%states - exclusive the two absorbing boundaries
x = -(m-1)/2:(m-1)/2;
%transition probabilities
p1 = repmat((1-pstay)./2 .* (1-d),1,length(x));
p2 = repmat((1-pstay)./2 .* (1+d),1,length(x)); 
p3 = repmat(pstay,1,length(x)); %stay
%constructing the transition probability matrix
tm = zeros(m-2,m);
i = (1:m-2)';
j = i+1;
pi = sub2ind(size(tm),i,i);%lower
ri = sub2ind(size(tm),i,j);%middle
qi = sub2ind(size(tm),i,j+1);%upper
tm(pi) = p1(2:m-1);
tm(ri) = p3(2:m-1);
tm(qi) = p2(2:m-1);
Q = sparse(tm(:,2:m-1));
%Q = tm(:,2:m-1);
R = sparse([tm(:,1) tm(:,m)]);
%R = [tm(:,1) tm(:,m)];
%initial disribution,
Z =  disnormalpdf(0, tau, x(2:(length(x)-1))); %dblLaplacePDFTrunc(x(2:(length(x)-1)),tau,theta);
%Z1 = zeros(1,m-2);
%Z1(1, (m-3)/2+1) = 1;
%numZ = round(z_range.*(m-2));
%pz = 1./numZ;
%bottomZ = round(((m-2)-numZ)./2);
%Z1(bottomZ:numZ) = pz; 
%identity matrix
I = speye(m-2);
%I = eye(m-2);