%% image at 2s, temperature 
xcenter=168;
ycenter=168;
zcenter=27;%cell6 nuc

load('cell1_pt_RI.mat')
load('cell1_pt_zint.mat')
load('cell_off_RI.mat')
load('cell_off_zint.mat')
off_RI = cell_off_RI(50:250,80:280);
off_zint = cell_off_zint(50:250,80:280);
temperature = cell1_pt_RI(50:250,80:280,6)*10^4;
pt_2s_ratio = -(cell1_pt_zint(50:250,80:280,105)-cell1_pt_zint(50:250,80:280,6))./off_zint;

figure;imagesc(temperature);daspect([1 1 1]);colorbar;caxis([0 6.9])
figure;imagesc(off_RI);daspect([1 1 1]);colorbar;colormap gray;caxis([-9E-4 0.017])
figure;imagesc(pt_2s_ratio.*(off_zint>0.027));daspect([1 1 1]);colorbar;colormap(mymap);caxis([-0.066 0.035])

% ratio colormap
(0-(-0.035))/(0.066-(-0.035))*255
mymap1(1:168,1)=(1:168)/168;
mymap1(1:168,2)=(1:168)/168;
mymap1(1:168,3)=1;
mymap1(169:256,1)=1;
mymap1(169:256,2)=1-(1:88)/88;
mymap1(169:256,3)=1-(1:88)/88;
ax = gca;
mymap = colormap(mymap1);
save('MyColormap','mymap')

input =off_RI-(-9E-4);
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/(0.017-(-9E-4))*255;
GrayIndex = uint8(GrayIndex);
Map       =gray(255);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='off.bmp';
imwrite(RGB3,filename);
daspect([1 1 1])


pt_2s_ratio=pt_2s_ratio.*(off_zint>0.027);
input =pt_2s_ratio-(-0.066);
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/(0.066-(-0.035))*255;
GrayIndex = uint8(GrayIndex);
Map       =colormap(mymap);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='pt_2s_ratio.bmp';
imwrite(RGB3,filename);
daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% input =-pt_2s_ratio.*(pt_2s_ratio<0);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/0.066*255;
% GrayIndex = uint8(GrayIndex);
% Map       =gray(255);
% RGB3      = ind2rgb(GrayIndex, Map);
% figure;imshow(RGB3);
% filename='pt_2s_minus_ratio.bmp';
% imwrite(RGB3,filename);
% daspect([1 1 1])
% 
% input =pt_2s_ratio.*(pt_2s_ratio>0);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/0.035*255;
% GrayIndex = uint8(GrayIndex);
% Map       =gray(255);
% RGB3      = ind2rgb(GrayIndex, Map);
% figure;imshow(RGB3);
% filename='pt_2s_plus_ratio.bmp';
% imwrite(RGB3,filename);
% daspect([1 1 1])

%% temporal evolution
nr=-cell1_pt_zint./cell_off_zint;
countnr=zeros(150,1);
for k=1:150
    count_num=0;
for x=xcenter-5:xcenter+5
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
