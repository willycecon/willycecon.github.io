is ,%% Asset Pricing 2, Homework 2
%% Question 2
% load data
clear;
p30=xlsread('StockPortfolios.csv','B5:AE604');
p2=xlsread('FFMktFct.csv','B5:C604');

% construct data
T=600;
N=30;
rf=p2(:,1);
Rf=p2(:,1)+1;
ell=ones(T,1);
z=p30./100-repmat(rf,1,30);

% construct PCs
Sigmaz=cov(z);
[V,D]=eig(Sigmaz);
[eigvalue,ind]=sort(diag(D),'descend');
Ds=D(ind,ind);
Vs=V(:,ind);
F=z*Vs(:,1:5);

% (a) Assumption IID
mufa=mean(F)';
ba=inv(Ds(1:5,1:5))*mufa;

% (b)
% GMM Model
syms b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5;
b=[b1 b2 b3 b4 b5]';
muf=[muf1 muf2 muf3 muf4 muf5]';
vars=[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5];
M=1-(F-repmat(muf',T,1))*b;                            % symbol function
moment=[z'*M./T;F'*M./T;(F-repmat(muf',T,1))'*ell./T]; % symbol function

% Construct d
muft=mean(F)';         % temporary proxy for muf used in d
d=[-Vs(:,1:5)*Ds(1:5,1:5) Vs(:,1:5)*muft*b'; -Ds(1:5,1:5) muft*b'; zeros(5,5), -eye(5)];

% First Stage Estimation
eqns_f=d'*moment;
[b1_f,b2_f,b3_f,b4_f,b5_f,muf1_f,muf2_f,muf3_f,muf4_f,muf5_f]=vpasolve(eqns_f,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[ba',mufa']);
Mf=subs(M,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[b1_f,b2_f,b3_f,b4_f,b5_f,muf1_f,muf2_f,muf3_f,muf4_f,muf5_f]);
muf_f=[muf1_f,muf2_f,muf3_f,muf4_f,muf5_f];
b_f=[b1_f,b2_f,b3_f,b4_f,b5_f];
% Second Stage Estimation
moment_vec=[z.*Mf,F.*Mf,F-muf_f];
S=cov(moment_vec);
eqns_s=d'*pinv(S)*moment;   % use pseudo-inverse as weighting matrix
[b1_s,b2_s,b3_s,b4_s,b5_s,muf1_s,muf2_s,muf3_s,muf4_s,muf5_s]=vpasolve(eqns_s,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[ba',mufa']);
Ms=subs(M,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[b1_s,b2_s,b3_s,b4_s,b5_s,muf1_s,muf2_s,muf3_s,muf4_s,muf5_s]);    % substitute estimates from second stage
muf_s=[muf1_s,muf2_s,muf3_s,muf4_s,muf5_s];
b_s=[b1_s,b2_s,b3_s,b4_s,b5_s];

% (c) exclude the factor pricing conditions
% only N+K moments
momentc=[z'*M./T;(F-repmat(muf',T,1))'*ell./T];       % symbol function
% construct d
dc=[-Vs(:,1:5)*Ds(1:5,1:5) Vs(:,1:5)*muft*b'; zeros(5,5), -eye(5)];
% First Stage
eqns_cf=dc'*momentc;
[b1_cf,b2_cf,b3_cf,b4_cf,b5_cf,muf1_cf,muf2_cf,muf3_cf,muf4_cf,muf5_cf]=vpasolve(eqns_cf,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[ba',mufa']);
M_cf=subs(M,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[b1_cf,b2_cf,b3_cf,b4_cf,b5_cf,muf1_cf,muf2_cf,muf3_cf,muf4_cf,muf5_cf]);
muf_cf=[muf1_cf,muf2_cf,muf3_cf,muf4_cf,muf5_cf];
b_cf=[b1_cf,b2_cf,b3_cf,b4_cf,b5_cf];
% Second Stage
momentc_vec=[z.*M_cf,F-muf_cf];
Sc=cov(momentc_vec);
eqns_cs=dc'*inv(Sc)*momentc;   
[b1_cs,b2_cs,b3_cs,b4_cs,b5_cs,muf1_cs,muf2_cs,muf3_cs,muf4_cs,muf5_cs]=vpasolve(eqns_cs,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[ba',mufa']);
Ms=subs(M,[b1 b2 b3 b4 b5 muf1 muf2 muf3 muf4 muf5],[b1_cs,b2_cs,b3_cs,b4_cs,b5_cs,muf1_cs,muf2_cs,muf3_cs,muf4_cs,muf5_cs]);    % substitute estimates from second stage
muf_cs=[muf1_cs,muf2_cs,muf3_cs,muf4_cs,muf5_cs];
b_cs=[b1_cs,b2_cs,b3_cs,b4_cs,b5_cs];