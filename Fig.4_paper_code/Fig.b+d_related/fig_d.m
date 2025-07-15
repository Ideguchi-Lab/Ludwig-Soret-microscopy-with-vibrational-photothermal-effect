xcenter=168;
ycenter=168;
zcenter=27;%cell6 cyto

load('cell1_pt_RI.mat')
load('cell1_pt_zint.mat')
load('cell_off_RI.mat')
load('cell_off_zint.mat')
off_zint = cell_off_zint(80:280,40:240);
pt_0_04s_ratio=-(cell1_pt_zint(80:280,40:240,7)-cell1_pt_zint(80:280,40:240,6))./off_zint;
pt_0_08s_ratio=-(cell1_pt_zint(80:280,40:240,9)-cell1_pt_zint(80:280,40:240,6))./off_zint;
figure; imagesc(pt_0_04s_ratio.*(off_zint>0.027));daspect([1 1 1]);colorbar;colormap(mymap);caxis([-4.5E-3 4.5E-3])
figure; imagesc(pt_0_08s_ratio.*(off_zint>0.027));daspect([1 1 1]);colorbar;colormap(mymap);caxis([-4.5E-3 4.5E-3])

%ratio
(0-(-4.5E-3))/(4.5E-3-(-4.5E-3))*255
mymap1(1:128,1)=(1:128)/128;
mymap1(1:128,2)=(1:128)/128;
mymap1(1:128,3)=1;
mymap1(129:256,1)=1;
mymap1(129:256,2)=1-(1:128)/128;
mymap1(129:256,3)=1-(1:128)/128;
ax = gca;
mymap = colormap(mymap1);
save('MyColormap','mymap')

pt_0_04s_ratio=pt_0_04s_ratio.*(off_zint>0.027);
input =pt_0_04s_ratio-(-4.5E-3);
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/(4.5E-3-(-4.5E-3))*255;
GrayIndex = uint8(GrayIndex);
Map       =colormap(mymap);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='pt_0_04s_ratio.bmp';
imwrite(RGB3,filename);
daspect([1 1 1])

pt_0_08s_ratio=pt_0_08s_ratio.*(off_zint>0.027);
input =pt_0_08s_ratio-(-4.5E-3);
max(max(input))
GrayIndex=input;
GrayIndex=GrayIndex/(4.5E-3-(-4.5E-3))*255;
GrayIndex = uint8(GrayIndex);
Map       =colormap(mymap);
RGB3      = ind2rgb(GrayIndex, Map);
figure;imshow(RGB3);
filename='pt_0_08s_ratio.bmp';
imwrite(RGB3,filename);
daspect([1 1 1])