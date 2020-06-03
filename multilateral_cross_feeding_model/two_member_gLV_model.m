function dxdt = two_member_gLV_model(t,x,c,beta,k)

dxdt = zeros(2,1);
dxdt(1) = c(1)*x(2)*(x(1)/(x(1)+beta))*(1-sum(x)/k);
dxdt(2) = c(2)*x(1)*(x(2)/(x(2)+beta))*(1-sum(x)/k);

end