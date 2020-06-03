function err = residual_function(params,y_obs)

%% unpack params
k_Mee2014       = params(1); % carrying capactiy
beta_Mee2014    = params(2); % Michaelis constant
C12_Mee2014     = params(3:2:end);
C21_Mee2014     = params(4:2:end);
assert (length(C12_Mee2014) == length(C21_Mee2014));

%% Simulate Mee's Lotka Volterra model
tol = 1e-6;
option = odeset('RelTol',tol,'AbsTol',tol*ones(1,2),'NonNegative',ones(1,2));

y_sim = zeros(length(C12_Mee2014),2);
for k=1:length(C12_Mee2014)
    [~,y] = ode15s(@two_member_gLV_model, [0,84], [1e7,1e7], option, [C12_Mee2014(k),C21_Mee2014(k)], beta_Mee2014, k_Mee2014);
    y_sim(k,1) = y(end,1)/1e7;
    y_sim(k,2) = y(end,2)/1e7;
end

%% compute difference between observed and simulated data

% SMAPE: symmetric mean absolute percentage error
err = sum(abs(y_obs(:)-y_sim(:))./(abs(y_obs(:))+abs(y_sim(:))))/length(y_obs(:));

% MSE: mean squared error
% err = sum((y_obs(:)-y_sim(:)).^2)/length(y_obs(:));

end

