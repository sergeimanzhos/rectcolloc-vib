
global xNN nd Nterms xnNNs xminp xmaxp xmint xmaxt 
[xNN, nd, Nterms, xnNNs, xminp, xmaxp, xmint, xmaxt] = PESload('fitN100NNs10.dat');


function v = Cl2COpes(Rij)
% received Rij in Bohr, matrix Npts x D
% order of Rij: CO, CCl1, CCl2, OCl1, OCl2, Cl1Cl2
% returns cm-1
Req = [1.1816  1.7518  1.7518  2.6047  2.6047  2.9021]*1.8897261254535; % Req = [COeq CCl1eq CCl2eq OCl1eq OCl2eq Cl1Cl2eq];
[Npts, ndim] = size(Rij);
Reqmat = ones(Npts,1)*Req;
xpes = exp(-5.0*(Rij-Reqmat)./Reqmat);
v = xNNperm(xpes); % in cm-1
end

% **** PES **************
function [xNN, nd, Nterms, xnNNs, xminp, xmaxp, xmint, xmaxt]=PESload(filename)
  PESdata = dlmread(filename);
  nd = PESdata(1,1);
  Nterms = PESdata(1,2);
  xnNNs = PESdata(1,3);
  for i=1:nd,
      xminp(i) = PESdata(2,i);
      xmaxp(i) = PESdata(3,i);
  end;    
  xmint = PESdata(2,nd+1);
  xmaxt = PESdata(2,nd+2);
  for i=1:Nterms,
      for j=1:nd+3,
          xNN(i,j) = PESdata(i+3,j);
      end;
  end;
end

function [out] = xneuron(x)
   out = 2d0./(1d0+exp(-2d0*x))-1d0; 
end

function out = xNNperm(x)
% this is symmetrized PES
    out = 0.5*(xNNimport(x)+xNNimport(x(:,[1 3 2 5 4 6])));
end

function out = xNNimport(x) 
% this is the NN PES before symmetrization. 
%What is called in the main program is it averaged over symmetry equivalent points (xNNperm).
global xNN nd Nterms xnNNs xminp xmaxp xmint xmaxt

   [Ndata, ndim] = size(x);
   % scale input x pn = 2*(p-minp)/(maxp-minp) - 1;
   xn = 2*(x-ones(Ndata,1)*xminp)./(ones(Ndata,1)*(xmaxp-xminp))-1;

   E=zeros(Ndata,1);
   for j=1:Ndata,
       for i=1:Nterms,
           xdot = dot(xNN(i,1:nd),xn(j,1:ndim));
           E(j) = E(j)+xneuron(xdot+xNN(i,ndim+2))*xNN(i,ndim+1)+xNN(i,ndim+3);
       end;
   E(j)=E(j)/xnNNs;    
   out(j) = xmint+0.5d0*(E(j)+1d0)*(xmaxt-xmint);
   end;
   
   out = out'; % to return a column
end