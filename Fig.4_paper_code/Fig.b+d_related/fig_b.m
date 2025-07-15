%% image at 2s, temperature 
xcenter=168; 
ycenter=168; 
zcenter=27;%cell6 cyto

load('cell1_pt_RI.mat')
load('cell1_pt_zint.mat')
load('cell_off_RI.mat')
load('cell_off_zint.mat')
off_RI = cell_off_RI(80:280,40:240);
off_zint = cell_off_zint(80:280,40:240);
temperature = cell1_pt_RI(80:280,40:240,6)*10^4;
pt_2s_ratio = -(cell1_pt_zint(80:280,40:240,105)-cell1_pt_zint(80:280,40:240,6))./off_zint;

figure;imagesc(temperature);daspect([1 1 1]);colorbar;caxis([0 5.5])
figure;imagesc(off_RI);daspect([1 1 1]);colorbar;colormap gray;caxis([-9E-4 0.017])
figure;imagesc(pt_2s_ratio.*(off_zint>0.027));daspect([1 1 1]);colorbar;colormap(mymap);caxis([-0.037 0.053])

%ratio
(0-(-0.053))/(0.037-(-0.053))*255
mymap1(1:106,1)=(1:106)/106;
mymap1(1:106,2)=(1:106)/106;
mymap1(1:106,3)=1;
mymap1(106:256,1)=1;
mymap1(106:256,2)=1-(1:151)/151;
mymap1(106:256,3)=1-(1:151)/151;
ax = gca;
mymap = colormap(mymap1);
save('MyColormap','mymap')

input =off-(-9E-4);
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

pt_2s_ratio=-pt_2s_ratio.*(off>0.027);
input =pt_2s_ratio-(-0.037);
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/(0.037-(-0.053))*255;
GrayIndex = uint8(GrayIndex);
Map       =colormap(mymap);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='pt_2s_ratio.bmp';
imwrite(RGB3,filename);
daspect([1 1 1])

% Merge plus and minus images by image J to generate Figure
% pt_2s_ratio=pt_2s_ratio.*(off>0.027);
% input =pt_2s_ratio.*(pt_2s_ratio>0);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/0.053*255;
% GrayIndex = uint8(GrayIndex);
% Map       =gray(255);
% RGB3      = ind2rgb(GrayIndex, Map);
% figure;imshow(RGB3);
% filename='pt_2s_plus_ratio.bmp';
% imwrite(RGB3,filename);
% daspect([1 1 1])
% 
% input =-pt_2s_ratio.*(pt_2s_ratio<0);
% max(max(input))
% GrayIndex=input;
% GrayIndex=GrayIndex/0.037*255;
% GrayIndex = uint8(GrayIndex);
% Map       =gray(255);
% RGB3      = ind2rgb(GrayIndex, Map);
% figure;imshow(RGB3);
% filename='pt_2s_minus_ratio.bmp';
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
