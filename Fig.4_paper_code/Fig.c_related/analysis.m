
%% -delta_sigma / sigma / Temperature change (/K) 
% [nucleus]
% 7.0E-3 COS0_nuc
% 8.5E-3 COS1_nuc
% 3.4E-3 COS2_nuc
% -9.4E-3 COS3_nuc
% 2.1E-3 COS4_nuc
% -8.6E-3 COS5_nuc
% 4.1E-3 COS6_nuc
% 7.1E-3 COS7_nuc
% 8.1E-3 COS8_nuc
% 1.40E-2 COS9_nuc
% 1.22E-2 COS10_nuc
% 5.7E-3 COS11_nuc

nuc=[5.7E-3, 1.22E-2, 1.40E-2, 8.1E-3, 7.1E-3, 4.1E-3, -8.6E-3, 2.1E-3, -9.4E-3, 3.4E-3, 8.5E-3, 7.0E-3]

%% cyto
% -8.1E-3 COS1_cyto
% -6.2E-3 COS2_cyto
% -1.2E-2 COS3_cyto
% -3E-2 COS4_cyto
% -5.4E-2 COS5_cyto
% -3.3E-2 COS6_cyto
% -1.67E-2 COS7_cyto
% -8.0E-3 COS8_cyto
% -1.45E-2 COS9_cyto
% -8.0E-3 COS10_cyto
% -6.0E-3 COS11_cyto

cyto=[-6.0E-3,-8.0E-3,-1.45E-2,-8.0E-3,-1.67E-2,-3.3E-2,-5.4E-2,-3E-2,-1.2E-2,-8.1E-3,-6.2E-3]

c=[1 1 1 1 1 1 1 1 1 1 1 1]
R=[2 2 2 2 2 2 2 2 2 2 2]

x=[1,2]
err=[std(nuc) std(cyto)]
av=[mean(nuc) mean(cyto)]
figure;plot(c,nuc,'o')
hold on
plot(R,cyto,'o')
bar(x,av)
errorbar(x,av,err)






