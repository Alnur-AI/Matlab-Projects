function dydt = odefcn(t,z,k1,vb,alpha,theta,k,S,m)
  dydt = zeros(2,1);
  dydt(1) = z(2);
  dydt(2) = ((-k1)*(z(2)^2) + ( vb*sin(alpha)-z(2)*sin(theta) ) *sin(theta)*k*S )/m;
end
