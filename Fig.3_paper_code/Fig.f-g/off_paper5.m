%% measurement&frequnecy_av
clear all
close all
path=('I:\soret_paper\RNA\control\control14\');
path_bg=('I:\soret_paper\RNA\control\control14\bg\');

%% Set the coordinates for the region of interest (center of heating spot in this case)
xcenter=160; ycenter=193; zcenter=25;%RNA10
% xcenter=161; ycenter=193; zcenter=28;%RNA11
% xcenter=162; ycenter=193; zcenter=29;%RNA12
xcenter=163; ycenter=193; zcenter=27;%control14

%% Set background regions used for normalization (manually defined)
y1=100;y2=200;x1=30;x2=50;%RNA10
% y1=210;y2=250;x1=210;x2=260;%RNA11
% y1=100;y2=200;x1=50;x2=100;%RNA12
y1=200;y2=300;x1=10;x2=50;%control14

%% get file names
rootpath_on=path_bg; 
rootpath1=strcat(path,'1\');
files1=dir([rootpath1,'*.pgm']);
rootpath2=strcat(path,'2\');
files2=dir([rootpath2,'*.pgm']);
rootpath3=strcat(path,'3\');
files3=dir([rootpath3,'*.pgm']);

files_on = dir([rootpath_on,'*.pgm']);
files1 = dir([rootpath1,'*.pgm']);
files2 = dir([rootpath2,'*.pgm']);
files3 = dir([rootpath3,'*.pgm']);
%% physical parameters
num_rot=10; % Number of illumination angles
wav=705*10^-9;% Illumination wavelength in meters
NA=0.85;% Numerical aperture
pixelsize=12*10^-6/(200/1.8*150/100);% Effective pixel size in the object space (meters)
ROI=[1,1,1004,1004];

%% calculate necessary physical parameters
dim=ROI(3)+1;% Size of the ROI
freq_per_pixel=1/(pixelsize*dim);% Frequency resolution
aperturesize=2*round((2*NA/wav)/freq_per_pixel/2)+1;% Aperture size in frequency domain

%% Detect center frequency component for normalization
t1=(double(imread(strcat(rootpath_on,files_on(1).name))));
t1=t1(720-502:720+502,720-502:720+502);
Meas_t1=fftshift(fft2(t1));
figure;imagesc(log(abs(Meas_t1)))
a=zeros(1005,1005);
for k=0:9
t1=(double(imread(strcat(rootpath_on,files_on(k+1).name))));
t1=t1(720-502:720+502,720-502:720+502);
Meas_t1=fftshift(fft2(t1));
a=a+Meas_t1;
figure(1);imagesc(log(abs(a)));daspect([1 1 1]);drawnow
Meas_t1=Meas_t1(1:400,1:1000);

[v,l]=max(abs(Meas_t1(:)));
[ii, jj] = ind2sub(size(Meas_t1), l);
coor_stack(:,:,k+1)=[jj,ii];%ii:y, jj+700:x
ii
jj
end

% Determine central frequency point based on average of selected frames
center_x=round((coor_stack(1,1,1)+coor_stack(1,1,6))/2);%821;
center_y=round((coor_stack(1,2,1)+coor_stack(1,2,6))/2);%198;

%% Create a circular mask in frequency domain
lengthx=aperturesize;
lengthy=lengthx;
mask=zeros(lengthx,lengthy);
for k=1:lengthx
    for l=1:lengthy
        if (k-(lengthx-1)/2-1)^2+(l-(lengthy-1)/2-1)^2 <= ((lengthx-1)/2)^2
            mask(k,l)=1;
        end
    end
end

%% select the center for each angle
k=1
for l=1:num_rot
    if(l<=4)    
    temp_on=double(imread(strcat(rootpath1,files1(160*(l-1)+k+2).name)));
    temp_on=temp_on(720-502:720+502,720-502:720+502);
    
    elseif(l<=8)
    temp_on=double(imread(strcat(rootpath2,files2(160*(l-5)+k+2).name)));
    temp_on=temp_on(720-502:720+502,720-502:720+502);

    elseif(l<=10)
    temp_on=double(imread(strcat(rootpath3,files3(160*(l-9)+k+2).name)));
    temp_on=temp_on(720-502:720+502,720-502:720+502);
    end

    temp_off=(double(imread(strcat(rootpath_on,files_on(l).name) )));
    temp_off=temp_off(720-502:720+502,720-502:720+502);
      
    Meas_on=fftshift(fft2(temp_on));
    temp=coor_stack(:,:,l);
    Meas_on=Meas_on/Meas_on(temp(2),temp(1));%normalized by DC intensity
    Meas_on=Meas_on(center_y-(lengthy-1)/2:center_y+(lengthy-1)/2,center_x-(lengthx-1)/2:center_x+(lengthx-1)/2).*mask;
    Meas_on=padarray(Meas_on,[round(aperturesize/2) round(aperturesize/2)],0);
     
    Meas_off=fftshift(fft2(temp_off));
    temp=coor_stack(:,:,l);
    Meas_off=Meas_off/Meas_off(temp(2),temp(1));
    Meas_off=Meas_off(center_y-(lengthy-1)/2:center_y+(lengthy-1)/2,center_x-(lengthx-1)/2:center_x+(lengthx-1)/2).*mask;
    Meas_off=padarray(Meas_off,[round(aperturesize/2) round(aperturesize/2)],0);   
    figure(3);imagesc(log(abs(Meas_off)));daspect([1 1 1]);drawnow
    
    temp_on=ifft2(ifftshift(Meas_on));
    temp_off=ifft2(ifftshift(Meas_off)); 
    temp_diff=temp_on./temp_off;
   
    temp_diff=temp_diff/mean2(temp_diff(y1:y2,x1:x2));
    figure;imagesc(phase_unwrap(angle(temp_diff)));
    
    %% determine shift vector and shifted mask
    [mymax,temp_idx]=max(log(abs(Meas_on)));
    [mymax,temp_idx_x]=max(mymax);
    temp_idx_y=temp_idx(temp_idx_x);
    idx_x(l,k)=temp_idx_x;
    idx_y(l,k)=temp_idx_y;      
    [x,y]=size(Meas_on);
    shift_x(l,k)=round(x/2)-temp_idx_x;
    shift_y(l,k)=round(x/2)-temp_idx_y;
    temp_mask=padarray(mask,[round(aperturesize/2) round(aperturesize/2)],0);
    temp_mask=circshift(temp_mask,[shift_y(l,k), shift_x(l,k)]);

    %% calculate the output
    diff_stack(:,:,l,k)=temp_diff;
    freq_stack(:,:,l,k)=fftshift(fft2(temp_diff)).*temp_mask;
    mask_stack(:,:,l,k)=temp_mask;
    figure(1);imagesc(log(abs(freq_stack(:,:,l,k))));daspect([1 1 1]);drawnow
l
end

%% prepare parameters for ODT reconstruction
temp_img=diff_stack(:,:,1,1);
[dim dim]=size(temp_img);
dim_img=dim;
dim_freq=floor(dim/2);
dimx=dim_img;
dimy=dim_img;
dimz=dim_img;

%% define physical parameter
n_m=1.330;
pixelsize=wav/(4*NA);
freq_pixel_size=1/(dim_img*pixelsize); %size of one pixel in the simulated frequency space
f_lambda=1/wav;

%% create mask
mask=zeros(dim_freq,dim_freq);
for k=1:dim_freq
    for l=1:dim_freq
        if (k-(dim_freq-1)/2-1)^2+(l-(dim_freq-1)/2-1)^2 <= ((dim_freq-1)/2)^2
        mask(k,l)=1;
        end
    end
end
mask=padarray(mask,[(dim_img-dim_freq)/2,(dim_img-dim_freq)/2]);

%% calculate fz sphere
fx=linspace(-0.5/pixelsize,0.5/pixelsize,dimx);
fy=linspace(-0.5/pixelsize,0.5/pixelsize,dimy); 
dfxy=abs(fx(1)-fx(2));
[fx,fy]=meshgrid(fx,fy);
fz=(n_m*f_lambda)^2-fx.^2-fy.^2;
fz=sqrt(fz).*mask;
filt_size=aperturesize*2-9;
filter=zeros(filt_size,filt_size)+1;
filter=padarray(filter,[(dimx-filt_size)/2,(dimx-filt_size)/2]);
filter=imgaussfilt(filter,3);

%% centerize all the measured freq spectra
for k=1
    for l=1:num_rot
        
    %% calculate bg-subtracted Rytov field
    temp_field=diff_stack(:,:,l,k);
    temp_field(:,1)=temp_field(:,2);
    temp_field(1,:)=temp_field(2,:);
    temp_amp=abs(temp_field);
    temp_phase=phase_unwrap(angle(temp_field*exp(1i*-2.5)));
    temp_phase=temp_phase-mean2(temp_phase(y1:y2,x1:x2));
    Rytov_field=log(temp_amp)+1i*(temp_phase);
    Rytov_freq=fftshift(fft2(Rytov_field))*pixelsize^2;

    %% estimate the illumination wavevector  
    fx_illum=shift_x(l,k)*freq_pixel_size;
    fy_illum=shift_y(l,k)*freq_pixel_size;
    fz_illum=sqrt((n_m*f_lambda)^2-fy_illum^2-fx_illum^2);
    Rytov_freq=Rytov_freq.*mask_stack(:,:,l,k);

    %% stack necessary parameters
    freq_Rytov_stack(:,:,l,k)=Rytov_freq;
    fz_stack(:,:,l,k)=circshift(fz,[shift_y(l,k),shift_x(l,k)]);%k_z
    fz_illum_stack(:,:,l,k)=fz_illum;%ki_z
    fx_illum_stack(l)=fx_illum;
    fy_illum_stack(l)=fy_illum;
    temp_NA=sqrt(fx_illum^2+fy_illum^2)/f_lambda    
    end
end
mask_stack=mask_stack>0;

%% perform diffraction tomography
clear fz_coor_ref fz_corr_ref_matrix
[x y]=size(mask_stack(:,:,1));
factor=round(wav/2/(n_m-sqrt(n_m^2-0.56^2))/pixelsize/2);%7;
div=round(503*(74/104)/factor)/351;
fz_passband=1/pixelsize/factor; %bandwidth of the reconstructed frequency spectra in z direction
fz_frequency_per_pixel=fz_passband/(round(dimz*div)-1);%-1
fz_coor_ref_matrix=zeros(x,y,round(dimz*div));
temp_fz_matrix=zeros(x,y,round(dimz*div));

for s=1:dimz*div
    fz_coor_ref(s)=fz_passband/2-fz_frequency_per_pixel*(s-1); %coordinate of the kz in the reconstructed frequency space
    fz_coor_ref_matrix(:,:,s)=fz_coor_ref(s);
end

K=1
%map 2D scattering spectrum to 3D scattering spectrum
freq_3d=zeros(dimx,dimy,round(dimz*div)); %matrix to which the frequency spectrum will be stored
freq_3d_mask=zeros(dimx,dimy,round(dimz*div)); %matrix to count how many times each pixel of the 3D frqeuenecy space is averaged
temp_fz_matrix=zeros(x,y,round(dimz*div));

for l=1:num_rot     
    F=2*i*2*pi*fz_stack(:,:,l,K).*freq_Rytov_stack(:,:,l,K); %relate the 2D scattering spectrum to the 3D scattering potetial F          
    temp_mask=mask_stack(:,:,l,K);%freqeucny bassband for the current loop, defined by the collection objective

    for s=1:dimz*div
        temp_fz_matrix(:,:,s)=fz_stack(:,:,l,K)-fz_illum_stack(:,:,l,K); %kz coordinate of the current loop
    end
    
    [mymin idx]=min(abs(fz_coor_ref_matrix-temp_fz_matrix),[],3);
    idx=idx.*temp_mask;
    for u=1:dimx
        for v=1:dimy            
            if idx(u,v) > 0
                freq_3d(u,v,idx(u,v))=freq_3d(u,v,idx(u,v))+F(u,v);
                freq_3d_mask(u,v,idx(u,v))=freq_3d_mask(u,v,idx(u,v))+1;
            end   
        end
    end
    K
end

%average mult-counted frequency components
freq_3d=freq_3d./((freq_3d_mask+(freq_3d_mask==0)));

%visualize the calculated frequency in ky-kz domain
[x y z]=size(freq_3d);
temp=freq_3d(round(x/2),:,:);   
temp=squeeze(temp);

%come back to the real 3D space and show refractive index
n_m=1.33;
reconst_3d=ifftn(ifftshift(freq_3d))/(pixelsize)^3/factor; %normalize after ifft to retrieve correct r-index value
rindex=sqrt(n_m^2-reconst_3d/(f_lambda)^2/(2*pi)^2); %calculated r-index value
rindex=fftshift(rindex,3);

close all
temp_rindex=rindex;
cell_off_zint=sum((squeeze(real(temp_rindex(:,:,zcenter-10:zcenter+10))))-n_m,3);
% cell_off_zint=cell_off_zint-(mean2(cell_off_zint(185:263,68:110))+mean2(cell_off_zint(54:103,225:257)))/2+0.01;%RNA10
% cell_off_zint=cell_off_zint-(mean2(cell_off_zint(210:250,210:260))+mean2(cell_off_zint(210:250,210:260)))/2;%RNA11
% cell_off_zint=cell_off_zint-(mean2(cell_off_zint(125:155,98:117))+mean2(cell_off_zint(220:290,290:324)))/2;%RNA12
cell_off_zint=cell_off_zint-mean2(cell_off_zint(230:285,48:100))+0.01;%control14

drymass=mean2(cell_off_zint(ycenter-5:ycenter+5,xcenter-5:xcenter+5))*0.207*7/0.2*1000;
disp(['drymass (fg/um^2) = ', num2str(drymass)])
figure;imagesc(cell_off_zint);daspect([1 1 1]);colorbar;drawnow;colormap gray
save cell_off_zint.mat cell_off_zint

