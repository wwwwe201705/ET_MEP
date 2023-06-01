function [ET,E_obs,T] = MEP_ET(name)
%UNTITLED 此处提供此函数的摘要
%   此处提供详细说明

load(name)
T=data(:,1);
TA=data(:,2);
Rh=data(:,3);
Rn=data(:,4);
Ts=data(:,5);
E_obs=data(:,7);
s=size(data);
model=0;
if s(2)==9
fvc=data(:,9);
else
    fvc=1;
   
end
%%
z=2.5;
Is=600;


[ e_s ] = e_sat( TA);


[qs] = q_surf( Rh./100.*e_s);


[ GMEP,EMEP,HMEP ] = MEP_Ev(Rn, Ts, qs, Is, z);    %%计算裸土蒸发
Ev=EMEP;
[EMEP1,HMEP1,GMEP1 ] =MEP_Tr(Rn,TA,qs,Is,z);    %%计算植被蒸散

Tr=EMEP1;

ET=(1-fvc).*Ev+fvc.*Tr;


E_obs=data(:,7);
% nse=NSE(ET,E_obs);
x=E_obs;
y=ET;                                   
mdl = fitlm(x,y);
R2=mdl.Rsquared.Adjusted;

LE=ET;
LE_stic=E_obs;
[BIAS,MAB,RMSE,MAPD,RMSEs,RMSEu,KGE,nse]= f_rmse(LE(~isnan(LE) & ~isnan(LE_stic)),LE_stic(~isnan(LE) & ~isnan(LE_stic)));
data(:,10)=ET;
data(1:4,11)=[nse,R2,RMSE,KGE];
name1=strcat(name,'.mat');


 save(name1,'data')  %%保存计算结果及评价指标
end