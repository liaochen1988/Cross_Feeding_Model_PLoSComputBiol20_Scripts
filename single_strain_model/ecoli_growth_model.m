function dxdt = ecoli_growth_model(t,x,S_e,k_m,m_e1,m_e2,k_e1,k_e2,Km_s,Km_m,Ki_r,k_r,d_r,alpha,Ym,K_ia,antibiotic)

Mhat = x(1);
E2 = x(2);
R = x(3);
E1 = (alpha-m_e2*E2)/m_e1;

J_in = k_e1*E1*S_e/(S_e+Km_s);
J_con = k_e2*E2*Mhat/(Km_m+Mhat)*K_ia/(K_ia+antibiotic);
J_grow = J_con*Ym;


dxdt = zeros(3,1);
dxdt(1) = J_in - J_con - J_grow*Mhat - k_m*Mhat;
dxdt(2) = alpha*J_grow/m_e2*Ki_r/(Ki_r+R)-J_grow*E2;
dxdt(3) = k_r*Km_m/(Km_m+Mhat) - d_r*R;

end

