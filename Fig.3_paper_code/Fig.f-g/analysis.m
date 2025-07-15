
%% raw
change_control=[4.4 4.4 3.3 4.83 1.69 2.38 4.04 3.22 1.56 3.99 2.26 3.79 2.23 4.41]*10^-3*0.207*7/0.2*1000;%fg/um^2/K
temperature_control=[6.5 5.3 4.4 4.5 3.2 3.2 5.01 4.91 5.03 4.53 3.20 4.56 4.13  5.26];
cp_control=change_control./temperature_control

drymass_control=[1.40 0.49 0.58 0.86 0.46 0.68 0.67 1.01 0.51 1.19 0.75 0.72 0.81 1.10]*1000%fg/um^2
cp_control=[4.9 XX XX 7.7 XX XX 5.7 XX XX 6.4 XX XX XX 6.0]%fg/um^2

c=[1 0.95 1 1 1 1 1 1 1 1 1 1 1 1]

change_RNA=[3.5 3.6 3.8 1.1 2.84 1.46 2.06 2.72 2.77  2.88 2.14 2.12]*10^-3*0.207*7/0.2*1000;
temperature_RNA=[4.4 4.7 5.9 4.2 5.3 3.6 3.98 4.16 4.41 3.51 3.35 3.83];
cp_RNA=change_RNA./temperature_RNA;
drymass_RNA=[0.118 0.0969 0.112 0.0708 0.0840 0.0684 0.0693 0.110 0.105 0.1759 0.0793 0.0666]*0.207*7/0.2*1000;
R=[2 2 2 2 2.05 1.95 2 2 2 2 2 2]

x=[1,2]
err=[std(cp_control) std(cp_RNA)]
av=[mean(cp_control) mean(cp_RNA)]
figure;plot(c,cp_control,'o')
hold on
plot(R,cp_RNA,'o')
bar(x,av)
errorbar(x,av,err)

av = [mean(cp_control) mean(cp_RNA)] %drymass 18.8% down
mean(cp_RNA)/mean(cp_control)
av = [mean(drymass_control) mean(drymass_RNA)] %drymass 13.0% down
mean(drymass_RNA)/mean(drymass_control)

x=[1,2]
err=[std(drymass_control) std(drymass_RNA)]
av=[mean(drymass_control) mean(drymass_RNA)]
figure;plot(c,drymass_control,'o')
hold on
plot(R,drymass_RNA,'o')
bar(x,av)
errorbar(x,av,err)