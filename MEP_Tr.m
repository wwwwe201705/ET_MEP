function [ EMEP, HMEP, GMEP, I0 ] = MEP_Tr(Rn,TA,qs,Is,z)
%Type 1 Dense Canopy
%Type 0 grass, cropfield
% ratio = 1; % z_unstalbe / z_stable 
z_unstable = z;
ratio=1;
type=1;
% z_stable = z/10;
z_stable = z/ratio;

%MEP function
Lv = 2.5E6;
Rv = 461;
T0 = 273;
KA = T0 + TA;
alpha = 1; beta = 5; gamma2 = 9; kappa = 0.41; cp = 1000; g = 9.81;rho = 1.18;

% MEP
N = length(Rn);
HMEP = zeros(N,1);
I0 = ones(N,1).*1000;
sigma = Lv .^ 2 ./ cp ./ Rv .* qs ./ KA .^2; 
B = 6 .* (sqrt(1+11./36.* sigma)-1);
% ef_n = find(~isnan(Rn));

if type
    HMEP = Rn ./ (1+B);
else
    for i =1:N
        if Rn(i)>0
            C1 = sqrt(3) / alpha;
            C2 = gamma2 / 2;
            I0(i) = rho * cp * sqrt(C1 * kappa * z_unstable) * ( C2 * kappa * z_unstable * g / (rho * cp * KA(i)))^(1/6);
        else
            C1 = 1/(1+2*alpha);
            C2 = 2*beta;    
            I0(i) = rho * cp * sqrt(C1 * kappa * z_stable) * ( C2 * kappa * z_stable * g / (rho * cp * KA(i)))^(1/6);
        end
        
        if isnan(Rn(i))
            HMEP(i) = NaN;
        else
%             Rn(i)
            HMEP(i) = fzero(@(x) B(i)*x*abs(x)^(1/6) + x*abs(x)^(1/6)+ B(i)/sigma(i)*Is/I0(i)*x -abs(x)^(1/6)*Rn(i), Rn(i)*0.8);
            
        end
    end
    HMEP = HMEP + Rn - Rn; %Make NaN
end
    EMEP = B.*HMEP;
    GMEP = Rn - EMEP - HMEP;
end

