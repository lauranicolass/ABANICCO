clear;clc;close all
currentFolder = pwd;
currletter=currentFolder(1:2);
directory=fullfile(currentFolder,'\data');
directory2=fullfile(currentFolder,'\example_images');
addpath(directory,directory2)
load(fullfile(directory,'Regions_Angles.mat'))

load(fullfile(directory,'Regions_Angles1.mat'))


an_im=imread(fullfile(directory2,'totoro.jpg'));

wantfigure=1; % change to 0 if uninterested in obtaining figures
numberofshades=15; % number of colors detected per category
morphology=10; 
% This is the standard option to fill small holes and remove small objects of a 500x500 image.
% Larger images will benefit from higher number. 
% Smaller images should use a smaller option.

%% This is the short version:

[Results,Regions_Angles1a,subTablesa] = AnalyzeColorImage_final2_short(an_im, Regions_Angles,Regions_Angles1,numberofshades,morpho,wantfigure);


%% This has the full shade analysis: the number of shades is limited to 8 so that they fit in the final display

% [Results,Regions_Angles1a,subTablesa] = AnalyzeColorImage_final2(an_im, Regions_Angles,Regions_Angles1,morpho,wantfigure);

%% Retrieve colors:

blue=Results(10).ResultingRGB;
figure;imshow(blue,[]);title('The blue','FontSize',15,FontName='Arial')

%% Create your own clusters:

% Define the regions present in cluster:
colors_forblue={'Blue','Ultramarine','Teal'};

mask=zeros(size(Results(1).Resultingmask));
for i=1:length(colors_forblue)
colorw=colors_forblue{i};
ind_colorw=find(strcmp({Results.Name}.',colorw));
maskind=Results(ind_colorw).Resultingmask;
mask=logical(mask+maskind);
end

blue_mask=mask;clear mask
New_blue=bsxfun(@times, an_im, cast(blue_mask,class(an_im)));
figure;imshow(New_blue,[]);title('The blues cluster','FontSize',15,FontName='Arial')

%% Re-define boundaries:
% re-convert to lab:
%conver to lab
imlab=rgb2lab(an_im);
a=imlab(:,:,2); A=a(:);
b=imlab(:,:,3); B=b(:);
points=[A B];
AB=[A B];

[thetas,rhos] = cart2pol(AB(:,2),AB(:,1));

thetas2=wrapTo360(rad2deg( thetas(:)));

ABpolar=[thetas rhos];

% Using the polar description decide the new angles and radious of the
% region:
newtheta1=90;
newtheta2=180;
rhonew1=30;


indbk1=find(thetas2>newtheta1);
indbk2=find(thetas2<newtheta2);
indbk3=find(rhos>rhonew1);

ABC_inter=intersect(intersect(indbk1,indbk2,'stable'),indbk3,'stable');

keepbk=zeros(size(A));
keepbk(ABC_inter)=1;
maskbk=reshape(keepbk,size(a));
bk_new2 = bsxfun(@times, an_im, cast(maskbk,class(an_im)));
figure;imshow(bk_new2);title('The new segmentation','FontSize',15,FontName='Arial')