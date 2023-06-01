function[BIAS,MAB,RMSE,MAPD,RMSEs,RMSEu,KGE,nse] = f_rmse(y,yhat)
bias = yhat - y;
BIAS = nanmean(bias);
ABSbias = abs(yhat - y);
MAB = mean(ABSbias);
sqerr = bias.^2;
meansqerr = mean(sqerr);
RMSE = sqrt(meansqerr);
percentDEV = abs((y - yhat)./y); percentDEV(percentDEV>100) = NaN;
MAPD = 100.*nanmean(percentDEV);

%% SYSTEMATIC AND UNSYSTEMATIC RMSE
X = [y ones(size(y))];
Y = yhat;
m = X\Y;
%m(1) = m; m(2) = c;

PRED  = m(1).*y + m(2);
RMSEs = sqrt(sum((PRED - y).^2)./length(y));
RMSEu = sqrt(sum((yhat - PRED).^2)./length(y));


%% Kling-Gupta efficiency
r           = corrcoef(y,yhat);
Rpearson    = r(1,2);
SDyhat      = std(yhat);
SDy         = std(y);

MUyhat      = std(yhat);
MUy         = std(y);


KGE = 1 - sqrt((Rpearson - 1).^2 + ((SDyhat./SDy) - 1).^2 + ((MUyhat./MUy) - 1).^2);
%%NSE
sim=y;
obs=yhat;
ave_obs = sum(obs(:))/numel(obs);   %实测数据平均数
Numerator = sum(power(obs-sim,2));  %分子
Denominator = sum(power(obs-ave_obs,2));%分母
nse = 1 - Numerator/Denominator;
end
