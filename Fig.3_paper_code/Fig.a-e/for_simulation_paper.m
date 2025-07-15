%% Temporal plot of experimental result(fg/um^2)
xcenter=168;ycenter=168;zcenter=29;
load('cell1_pt_zint')
nr=cell1_pt_zint;
shift=120;
num_time=30;
countnr=zeros(shift+1+num_time-1,1);
for k=1:shift+1+num_time-1
    count_num=0;
for x=xcenter-5:xcenter+5
    for y=ycenter-5:ycenter+5
      countnr(k,1)=countnr(k,1)+nr(y,x,k); 
      count_num=count_num+1;
    end
end
countnr(k,1)=countnr(k,1)/count_num;
end
figure;plot(([1:150]-5)*0.02,-countnr(1:shift+1+num_time-1),'o')

countnr(6:125)=countnr(6:125)-countnr(6);
figure;plot(([1:150]-5)*0.02,-countnr(1:shift+1+num_time-1),'o')
temporal_evolution(:,1)=([1:150]-5)*0.02;
temporal_evolution(:,2)=countnr(1:shift+1+num_time-1);
temporal_evolution_drymass(:,1) = temporal_evolution(:,1);
temporal_evolution_drymass(:,2) = temporal_evolution(:,2)*0.207*7/0.2*1000;
figure;plot(temporal_evolution_drymass(:,1), -temporal_evolution_drymass(:,2))
csvwrite('temporal_evolution_drymass.csv',temporal_evolution_drymass)

%% Mask generation
mask=zeros(351,351);
mask(115:215,136:203)=1;

load('cell_off_RI')
a=mask.*cell_off_RI;
figure;imagesc(mask.*cell_off_RI);daspect([1 1 1]);caxis([0 0.012])

mask(203:215,136:150)=0;
mask(205:215,189:203)=0;
mask(210:215,176:188)=0;
mask(115:129,144:151)=0;
mask(115:134,185:203)=0;
mask(135:144,192:203)=0;
mask(189:202,136:142)=0;
mask(115:145,136:143)=0;
mask(115:120,152:184)=0;

mask=mask(51:300,51:300);
cell_off_RI_crop=cell_off_RI(51:300,51:300);
figure;imagesc(cell_off_RI_crop.*mask);daspect([1 1 1]);colormap gray

mask=mask(118-70:118+70,118-70:118+70);
figure;imagesc(mask);daspect([1 1 1]);
csvwrite('mask_crop.csv',mask)

%% initial concentration of drymass (fg/um^2)
load('cell_off_zint')
off_zint= cell_off_zint(51:300,51:300)*0.207*7/0.2*1000;%zint
off_zint=off_zint(118-70:118+70,118-70:118+70);
figure;imagesc(off_zint);
csvwrite('initial_con_drymass.csv',off_zint) 

%% temperature (K)
load('cell1_pt_RI')
pt_0_02s=cell1_pt_RI(51:300,51:300,6);%use RI image
temperature = pt_0_02s(118-70:118+70,118-70:118+70)*10^4;
figure;imagesc(temperature);daspect([1 1 1]);colorbar;colormap jet;%caxis([-6.5E-4 10E-4])
csvwrite('temperature.csv',temperature) 
