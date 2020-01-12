%% Asset Pricing 2, Homework 2
%% Question 1
% load data
clear;
p30=xlsread('StockPortfolios.csv','B5:AE604');
p2=xlsread('FFMktFct.csv','B5:C604');

% construct data
T=604-4;
N=30;
rf=p2(:,1);
Rf=p2(:,1)+1;
ell=ones(T,1);
F=p2(:,2)-rf;
z=p30./100-repmat(rf,1,30);

% (a) Assumption IID
muf1=mean(F);
b1=var(F)^(-1)*muf1;
stdf1=std(F);
M1=1-(F-muf1)*b1;
alpha1=z'*M1./T;
e1=zeros(T,N);
for i=1:30
    temp=cov(z(:,i),F)/var(F);
    temp=temp(1,2);
    e1(:,i)=z(:,i)-F*(temp);
end
test1=T.*alpha1'*inv(cov(e1))*alpha1./(1+muf1'*inv(var(F))*muf1);
HJ1=alpha1'*inv(cov(e1))*alpha1;

% (b)
% GMM Model
clear b muf;
syms b muf;
vars=[b muf];
M=1-(F-muf)*b;                            % symbol function
moment=[z'*M./T;F'*M./T;(F-muf)'*ell./T]; % symbol function
% Construct d
beta=(z-repmat(mean(z),T,1))'*F./(T*var(F));
Sigmaf=cov(F);
muft=mean(F); % temporary proxy for muf used in d
muz=mean(z)'; % used in d
d=[-beta*Sigmaf, muz*b; -Sigmaf, muft*b; 0, -1];

% First Stage Estimation
eqns=d'*moment;
[b_f, muf_f]=vpasolve(eqns,vars,[b1,muf1]);
Mf=subs(M,[b,muf],[b_f,muf_f]);
alpha_f=z'*Mf./T;
test_f=T.*alpha_f'*inv(cov(e1))*alpha_f./(1+muf_f'*inv(var(F))*muf_f);
HJ_f=alpha_f'*inv(cov(e1))*alpha_f;

% Second Stage Estimation
moment_vec=[z.*Mf,F.*Mf,F-muf_f];
S=cov(moment_vec);
eqns2=d'*inv(S)*moment;
[b_s, muf_s]=vpasolve(eqns2,vars,[b1,muf1]);
Ms=subs(M,[b,muf],[b_s,muf_s]);    % substitute estimates from second stage
std(Ms)                            % standard errors of Ms
alpha_s=z'*Ms./T;                   % pricing errors
test_s=T.*alpha_s'*inv(cov(e1))*alpha_s./(1+muf_s'*inv(var(F))*muf_s);
HJ_s=alpha_s'*inv(cov(e1))*alpha_s;