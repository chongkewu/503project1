clear
close all
[snr,n1,var,lam,n2,N,Pd,Pfa] = Initial();
Roc_Monte(snr,n1,var,lam,n2,N,Pd,Pfa);
Roc_Ana(snr,n1,var,lam,n2,N,Pd,Pfa);


function [snr,n1,var,lam,n2,N,Pd,Pfa] = Initial()
    snr = -1:5;
    n1 = length(snr);
    var = 1./(10.^(snr./10));
    lam = 0:0.1:1;
    n2 = length(lam);
    N = 3*10^4;
    Pd = zeros(n1,n2);
    Pfa = zeros(n1,n2);
end
function Roc_Monte(snr,n1,var,lam,n2,N,Pd,Pfa)
for i=1:n1 % different variance
    for j = 1:n2 %different lambda
        X_H0 = randn(N,1)*sqrt(var(i));
        X_H1 = 1+randn(N,1)*sqrt(var(i));
        n_P_d = sum(X_H1 >= lam(j));
        n_P_fa = sum(X_H0 >= lam(j));
        if n_P_d <= n_P_fa %test1
            fprintf(" Numbers of Pd smaller than Pfa\n");
        end
        Pd(i,j) = n_P_d/N;
        Pfa(i,j) = n_P_fa/N;
    end
end
plt1(Pfa,Pd,n1);
end

function plt1(Pfa,Pd,n1)
figure
for i = 1:n1
    plot(Pfa(i,:),Pd(i,:)); 
    hold on
end
legend('\sigma = 1.2589','\sigma = 1',...
    '\sigma = 0.7943','\sigma = 0.6310',...
    '\sigma = 0.5012','\sigma = 0.3981',...
    '\sigma = 0.3162',...
    'Location','northwest')
title('ROC (Monte Carlo)')
axis([ 0 0.55 0.45 1 ])
xlabel('P_{fa}')
ylabel('P_d')
end

function Roc_Ana(snr,n1,var,lam,n2,N,Pd,Pfa)
mu = 0;
for i = 1:n1
    pd0 = makedist('Normal',mu,sqrt(var(i)));
    pd1= makedist('Normal',mu+1,sqrt(var(i)));
    Pd(i,:) = 1 - cdf(pd1,lam);
    Pfa(i,:) = 1- cdf(pd0,lam);
end
plt1(Pfa,Pd,n1);title('ROC (Analytical)')
end

