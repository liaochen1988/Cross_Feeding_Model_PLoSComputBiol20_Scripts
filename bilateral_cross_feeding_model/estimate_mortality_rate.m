% By Chen Liao, Memorial Sloan Kettering Cancer Center, June, 2 2019

% This script fits an exponential model to cfu data to estimate mortality
% rate of the lysine and leucine auxotroph

figure();
hold on;

% lysine auxotroph
t = [0.021, 10.062, 17.957, 25.020, 34.086, 46.020, 60.030, 68.030];
cfu = [8.620e6, 2.338e6, 1.831e6, 1.073e6, 3.450e5, 1.502e5, 9.215e4, 5.275e4];
plot(t,cfu,'o')
set(gca,'YScale','log');
tf = [0:0.1:100];
p = polyfit(t,log(cfu),1);
plot(tf,exp(p(2)+p(1)*tf),'-');

% leucine auxotroph
t = [0,10.009,24.956,33.959,45.923,59.976,67.995,78.960];
cfu = [5.418e6, 3.142e6, 6.029e6, 5.978e6, 5.374e6, 5.565e6, 5.234e6, 3.771e6];
plot(t,cfu,'o')
set(gca,'YScale','log');
tf = [0:0.1:100];
p = polyfit(t,log(cfu),1);
plot(tf,exp(p(2)+p(1)*tf),'-');
