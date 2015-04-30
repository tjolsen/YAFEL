function [ x, w ] = GaussLegQuadRule( polyOrder )
%GaussLegQuadRule(polyOrder) computes the Gauss-Legendre quadrature rule
%over [-1, 1] (nodes and weights) for a given polynomial degree (and corresponding
%number of quadrature points) following the procedure outlined at:
%
% http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
%

if(nargin == 0)
    polyOrder = 5;
end



%Newton-Raphson iteration to find roots of N-th order legrange polynomial
ii = (1:polyOrder)';
x = cos(pi*(ii-0.25)/(polyOrder+.5));
P = LegPoly(polyOrder, x);


TOLERANCE = 1.0e-14;
while(any(abs(P(:,end))>TOLERANCE))
    dP = (polyOrder./(x.^2 - 1)).*(x.*P(:,end) - P(:,end-1));
    
    x = x - P(:,end)./dP;
    P = LegPoly(polyOrder, x);
end

dP = (polyOrder./(x.^2 - 1)).*(x.*P(:,end) - P(:,end-1));
w = 2./( (1-x.^2).*(dP.^2) );

end

function P = LegPoly(N,x)
%LegPoly: Computes Legendre polynomials of order N (P0...PN) at given "x"
%values

%ensure that x is a column vector... much easier to work with
if(size(x,1) < size(x,2))
    x = x';
end

P = zeros(length(x), N+1);

P(:,1) = 1;
P(:,2) = x;

for n = 2:N
    P(:,n+1) = ( (2*n-1)*x.*P(:,n) - (n-1)*P(:,n-1) )/n;
end


end