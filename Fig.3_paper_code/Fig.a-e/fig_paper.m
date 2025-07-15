%% temporal plot (Fig.3a)
load('cell1_pt_zint')
load('cell1_pt_RI')
load('cell_off_zint')
load('cell_off_RI')

% Center o fthe heating spot (left graph)
xcenter=168;
ycenter=168;
zcenter=29;

% Periphery of the heating spot (right graph)
% xcenter=174;
% ycenter=193;
% zcenter=29;

figure;imagesc(cell1_pt_RI(:,:,125)-cell1_pt_RI(:,:,6));daspect([1 1 1])
nr=cell1_pt_RI;
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

%% temperature image (Fig.2b)
pt_0_02s=cell1_pt_RI(51:300,51:300,6);%use RI image
off=cell_off_RI(51:300,51:300);%use RI image
figure;imagesc(pt_0_02s*10^4);daspect([1 1 1]);colorbar;colormap jet;%caxis([-6.5E-4 10E-4])
figure;imagesc(off);daspect([1 1 1]);colorbar;colormap jet;

temperature = pt_0_02s*10^4;
input =temperature;
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/6.4*255;
GrayIndex = uint8(GrayIndex);
Map       =gray(255);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='temperature.bmp';
imwrite(RGB3,filename);
daspect([1 1 1])

input =off-(-9E-4);
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/(0.022-(-9E-4))*255;
GrayIndex = uint8(GrayIndex);
Map       =gray(255);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='off.bmp';
imwrite(RGB3,filename);
daspect([1 1 1])

%% drymass change (Fig.3c)
pt_2_4s=-(((cell1_pt_zint(51:300,51:300,125)-cell1_pt_zint(51:300,51:300,6))))*0.207*7/0.2*1000; %[fg/um^2]

save pt_2_4s.mat pt_2_4s
figure;imagesc(pt_2_4s);daspect([1 1 1]);colorbar;colormap jet;caxis([-39.64 20.89])

% colormap
(0-(-39.64))/(20.89-(-39.64))*255
mymap1(1:167,1)=(1:167)/167;
mymap1(1:167,2)=(1:167)/167;
mymap1(1:167,3)=1;
mymap1(168:256,1)=1;
mymap1(168:256,2)=1-(1:89)/89;
mymap1(168:256,3)=1-(1:89)/89;
mymap = colormap(mymap1);
save('MyColormap','mymap')

% Merge plus and minus images by image J to generate Figure
% off_zint=cell_off_zint(51:300,51:300);%off: zint
% pt_2_4s_plus=pt_2_4s.*(off_zint>0.027); %off: zint
% pt_2_4s_plus=pt_2_4s_plus.*(pt_2_4s_plus>0);
% input =pt_2_4s_plus;
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(20.89)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_2_4s_plus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])
% 
% pt_2_4s_minus=pt_2_4s.*(off_zint>0.027);%off: zint
% pt_2_4s_minus=-pt_2_4s_minus.*(pt_2_4s_minus<0);
% input =pt_2_4s_minus;
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/39.64*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_2_4s_minus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])

%% Experiment time sequence images (Fig.2d)
% Experiment
pt_2_4s=-((cell1_pt_zint(51:300,51:300,125)-cell1_pt_zint(51:300,51:300,6)))*0.207*7/0.2*10^3; %ƒ¢‚š–ƒ°ƒ¢‚Ž/ƒ¿–‚P‚O‚O‚O[fg/um^2];
pt_1_0s=-((cell1_pt_zint(51:300,51:300,55)-cell1_pt_zint(51:300,51:300,6)))*0.207*7/0.2*10^3; %ƒ¢‚š–ƒ°ƒ¢‚Ž/ƒ¿–‚P‚O‚O‚O[fg/um^2];
pt_0_1s=-((cell1_pt_zint(51:300,51:300,10)-cell1_pt_zint(51:300,51:300,6)))*0.207*7/0.2*10^3; %ƒ¢‚š–ƒ°ƒ¢‚Ž/ƒ¿–‚P‚O‚O‚O[fg/um^2];
off_zint=cell_off_zint(51:300,51:300); %off: zint
figure;imagesc(pt_1_0s);caxis([-5.66*7 2.98*7]);daspect([1 1 1])
figure;imagesc(pt_2_4s);caxis([-5.66*7 2.98*7]);daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% pt_2_4s_plus=pt_2_4s.*(off_zint>0.027); %off: zint
% input =pt_2_4s_plus.*(pt_2_4s_plus>0);
% input=input(68-15:158+15,90-15:153+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(2.98*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_2_4s_plus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])
% 
% pt_2_4s_minus=pt_2_4s.*(off_zint>0.027);%off: zint
% input =pt_2_4s_minus.*(pt_2_4s_minus<0);
% input=-input(68-15:158+15,90-15:153+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(5.66*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_2_4s_minus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% pt_1_0s_plus=pt_1_0s.*(off_zint>0.027);%off: zint
% input =pt_1_0s_plus.*(pt_1_0s_plus>0);
% input=input(68-15:158+15,90-15:153+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(2.98*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_1_0s_plus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])
% 
% pt_1_0s_minus=pt_1_0s.*(off_zint>0.027);%off: zint
% input =pt_1_0s_minus.*(pt_1_0s_minus<0);
% input=-input(68-15:158+15,90-15:153+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(5.66*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_1_0s_minus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% pt_0_1s_plus=pt_0_1s.*(off_zint>0.027);%off: zint
% input =pt_0_1s_plus.*(pt_0_1s_plus>0);
% input=input(68-15:158+15,90-15:153+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(2.98*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_0_1s_plus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])
% 
% pt_0_1s_minus=pt_0_1s.*(off_zint>0.027);%off: zint
% input =pt_0_1s_minus.*(pt_0_1s_minus<0);
% input=-input(68-15:158+15,90-15:153+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(5.66*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='pt_0_1s_minus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])

% Overlay off image by image J to generate Figure
% off_RI = cell_off_RI(51:300,51:300);
% input =off_RI(68-15:158+15,90-15:153+15)-(-9E-4);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(0.022-(-9E-4))*255;
% GrayIndex = uint8(GrayIndex);
% Map       =gray(255);
% RGB3      = ind2rgb(GrayIndex, Map);
% figure;imshow(RGB3);
% filename='off_crop.bmp';
% imwrite(RGB3,filename);
% daspect([1 1 1])

%% Simulation time sequence images
% Simulation
sim_0_1s=double(readmatrix('sim_1s.csv'));
sim_1s=double(readmatrix('sim_10s.csv'));
sim_2_4s=double(readmatrix('sim_24s.csv'));
figure;imagesc(sim_1s);caxis([-5.66*7 2.98*7]);daspect([1 1 1])
figure;imagesc(sim_2_4s);caxis([-5.66*7 2.98*7]);daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% input =sim_0_1s.*(sim_0_1s>0);
% input=input(21-15:111+15,43-15:106+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(2.98*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='sim_0_1s_plus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])
% 
% input =-sim_0_1s.*(sim_0_1s<0);
% input=input(21-15:111+15,43-15:106+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(5.66*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='sim_0_1s_minus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% input =sim_1s.*(sim_1s>0);
% input=input(21-15:111+15,43-15:106+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(2.98*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='sim_1s_plus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])
% 
% input =sim_1s.*(sim_1s<0);
% input=-input(21-15:111+15,43-15:106+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(5.66*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='sim_1s_minus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% input =sim_2_4s.*(sim_2_4s>0);
% input=input(21-15:111+15,43-15:106+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(2.98*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='sim_2_4s_plus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])
% 
% input =-sim_2_4s.*(sim_2_4s<0);
% input=input(21-15:111+15,43-15:106+15);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(5.66*7)*255;
% GrayIndex = uint8(GrayIndex);
% figure;imshow(GrayIndex);
% filename='sim_2_4s_minus.bmp';
% imwrite(GrayIndex,filename);
% daspect([1 1 1])

%% profile comparison between sim and exp
load('pt_2_4s.mat')
figure;imagesc(pt_2_4s);caxis([-5.66*7 2.98*7])
exp_2_4s= pt_2_4s(118-70:118+70,118-70:118+70);
sim_2_4s=double(readmatrix('sim_24s.csv'));

% x-axis
figure;plot(([71-32:71+32+2]-72)*0.207,exp_2_4s(71,71-32:71+32+2),'o')
hold on
plot(([71-32:71+32+2]-72)*0.207,sim_2_4s(71,71-32:71+32+2))

%% temporal evolution (Fig.2e)
sim_temp=double(readmatrix('sim_temporal.csv'));
temporal=double(readmatrix('temporal_evolution_drymass.csv'));
figure;plot(temporal(1:150,1),-temporal(1:150,2),'o')
hold on
plot(0:0.1:2.4,sim_temp)