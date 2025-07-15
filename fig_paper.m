%% Generation of Fig.5 for each time point 
cd 'G:\soret_timelapse\20231021\cell2\time1' %dataset you want to analyze
%% RI image (Fig.5a)
load('cell_off_RI.mat')
off=cell_off_RI(51:300,51:300);
figure;imagesc(off);daspect([1 1 1]);colorbar;colormap gray;caxis([-9E-4 0.03])

input =cell_off(51:300,51:300)-(-9E-4);
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/(0.03-(-9E-4))*255;
GrayIndex = uint8(GrayIndex);
Map       =gray(255);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='time1.bmp'; %change the name for proper time
imwrite(RGB3,filename);
daspect([1 1 1])

%% Delta_sigma/sigma (%) (Fig.5a)
load('cell1_pt_zint.mat')
load('cell_off_zint.mat')
pt_2_4s=(cell1_pt_zint(51:300,51:300,125)-cell1_pt_zint(51:300,51:300,6))./cell_off_zint(51:300,51:300);
off_zint=cell_off_zint(51:300,51:300);
figure;imagesc(off_zint);daspect([1 1 1]);colorbar;colormap jet;caxis([0 0.15])
figure;imagesc(-pt_2_4s.*(off_zint>0.07));daspect([1 1 1]);colorbar;caxis([-0.034 0.025]);colormap mymap

% ratio colormap
(0-(-0.034))/(0.025-(-0.034))*255
mymap1(1:147,1)=(1:147)/147;
mymap1(1:147,2)=(1:147)/147;
mymap1(1:147,3)=1;
mymap1(148:256,1)=1;
mymap1(148:256,2)=1-(1:109)/109;
mymap1(148:256,3)=1-(1:109)/109;
ax = gca;
mymap = colormap(mymap1);
save('MyColormap','mymap')

% input =-(pt_2_4s).*(off_zint>0.07)-(-0.034);
% max(max(input))
% min(min(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/(0.025-(-0.034))*255;
% GrayIndex = uint8(GrayIndex);
% Map       =colormap(mymap);
% RGB3      = ind2rgb(GrayIndex, Map);
% filename='pt_2_4s.bmp';
% imwrite(RGB3,filename);
% daspect([1 1 1])

%% time evolution (Fig.5b)
xcenter=172;
ycenter=196;
nr=cell1_pt_zint./cell_off_zint;
countnr=zeros(150,1);
for k=1:150
    count_num=0;
for x=xcenter-1:xcenter+1
    for y=ycenter-5:ycenter+5
      countnr(k,1)=countnr(k,1)+nr(y,x,k); 
      count_num=count_num+1;
    end
end
countnr(k,1)=countnr(k,1)/count_num;
end
a=countnr(1:150)-mean(countnr(1:5));
a(6:125)=a(6:125)-a(6);
figure;plot(([1:150]-5)*0.02,a,'o')

%% time-lapse change of Delta_sigma/sigma (%) (Fig.5c)

% Cell1 
time1=[0 35/60 1+18/60 1+53/60 2+47/60 3+52/60 5+27/60 6+27/60];% hour
change1=[3.02 3.07 3.10 3.26 3.40 4.03 0.13 0.016];% (%)

% Cell2 
time2=[0 34/60 1+11/60 1+53/60 2+31/60 3];% hour
change2=[3.18 2.88 2.16 1.84 0.79 0.18];% (%)

% Cell3
time3=[0 30/60 1+10/60 2+23/60]% hour
change3=[2.91 3.15 2.27 1.58]% (%)

% Cell4 
time4=[0 32/60 1+10/60 1+42/60 2+14/60]% hour
change4=[3.89 2.20 0.84 0.86 0.03]% (%)

figure;plot(time1,change1,'o')
hold on
plot(time2,change2,'o')
plot(time3,change3,'o')
plot(time4,change4,'o')