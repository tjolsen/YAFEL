x = [0 0 0;
     1 0 0;
     1.1 .5 1;
     0.7 1 0;
     2 0 0.25];
 
 tets = [1 2 3 4;
     2 3 4 5];
 
 face = [2,4,3];
 
 g = [2,.8,1]';
 
 phi = x*g;
 
 xa = x(2,:)';
 xb = x(3,:)';
 xc = x(4,:)';
 xi = mean(x(tets(1,:),:))';
 xj = mean(x(tets(2,:),:))';
 xface = mean(x(face,:))';
 
 phi_a = dot(xa,g);
 phi_b = dot(xb,g);
 phi_c = dot(xc,g);
 phi_i = dot(xi,g);
 phi_j = dot(xj,g);
 
 e1 = xc-xa; e1 = e1/norm(e1);
 e2 = xb-xa; e2 = e2/norm(e2);
 e3 = xj-xi; e3 = e3/norm(e3);
 n = cross(xc-xa, xb-xa); n = n/norm(n);
  
 Ecov = [e1,e2,e3];
 V = dot(e1,cross(e2,e3));
f1 = cross(e2,e3)/V;
f2 = cross(e3,e1)/V;
f3 = cross(e1,e2)/V;

Econtra = [f1,f2,f3];


dj = ((xj-xface)'*n);
di = -((xi-xface)'*n);

xij_j = xj - e3*dj/(n'*e3);
xij_i = xi + e3*di/(n'*e3);

epsilon = xface-xij_j;


dpdxcov = [(phi_c-phi_a)/norm(xc-xa);
    (phi_b-phi_a)/norm(xb-xa);
    (phi_j-phi_i)/norm(xj-xi)];



error = Econtra*dpdxcov - g

%------------------------------------------------
figure(1);
hold on;
for ti = 1:size(tets,1)
    t = tets(ti,:);
    plot3(x(t([1,2]), 1),x(t([1,2]), 2),x(t([1,2]), 3),'k.-');
    plot3(x(t([1,3]), 1),x(t([1,3]), 2),x(t([1,3]), 3),'k.-');
    plot3(x(t([1,4]), 1),x(t([1,4]), 2),x(t([1,4]), 3),'k.-');
    plot3(x(t([2,3]), 1),x(t([2,3]), 2),x(t([2,3]), 3),'k.-');
    plot3(x(t([2,4]), 1),x(t([2,4]), 2),x(t([2,4]), 3),'k.-');
    plot3(x(t([3,4]), 1),x(t([3,4]), 2),x(t([3,4]), 3),'k.-');
end
plot3([xi(1) xj(1)],[xi(2) xj(2)],[xi(3) xj(3)],'bo-');
plot3(xface(1),xface(2),xface(3),'rx');
plot3(xij_j(1),xij_j(2),xij_j(3),'ro');

quiver3(kron([1,1,1]',xij_j(1)), kron([1,1,1]',xij_j(2)), kron([1,1,1]',xij_j(3)), ...
    [e1(1),e2(1),e3(1)]',...
    [e1(2),e2(2),e3(2)]',...
    [e1(3),e2(3),e3(3)]',...
    .25);

quiver3(kron([1,1,1]',xij_j(1)), kron([1,1,1]',xij_j(2)), kron([1,1,1]',xij_j(3)), ...
    [f1(1),f2(1),f3(1)]',...
    [f1(2),f2(2),f3(2)]',...
    [f1(3),f2(3),f3(3)]',...
    .25);
axis equal

hold off;