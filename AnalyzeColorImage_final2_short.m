function [Regions_Angles,Regions_Angles1,subTables] = AnalyzeColorImage_final2_short(an_im, Regions_Angles,Regions_Angles1,numshades,morpho,figs)
LAB=rgb2lab((an_im));

%% Classify pixels according to where they lie in AB histogram
l=LAB(:,:,1);
a=LAB(:,:,2);
b=LAB(:,:,3);
L=l(:);
A=a(:);
B=b(:);
points=[A B];
AB=[A B];
rb=floor(AB+127);

totalpix_anim=size(a,1)*size(a,2);

figure('WindowState','maximized','Color',[1 1 1]);
tiledlayout('flow','TileSpacing','Compact');
nexttile
imshow(an_im);
title('New Image')
for i=1:length(Regions_Angles)
    ud=Regions_Angles(i).Points;
    [C] = intersect(rb,ud,'rows');
    indices = find(ismember(rb,C,'rows'));
    keepupperred=zeros(size(A));
    keepupperred(indices)=1;
    blue=reshape(keepupperred,size(a));
    %%% REMOVE SMALL OBJECTS AND FILL SMALL HOLES
    BW2 = bwareaopen(blue,morpho);
    binaryImage2 = ~bwareaopen(~BW2, 2*morpho);
    segmented_image = bsxfun(@times, an_im, cast(binaryImage2,class(an_im)));
    %%%
    nexttile
    imshow(segmented_image)
    title(strcat(Regions_Angles(i).Name,'-Region'))

    Regions_Angles(i).Resultingmask=binaryImage2;
    Regions_Angles(i).ResultingRGB=segmented_image;
    Regions_Angles(i).Resultingmask_percentage=(sum(binaryImage2(:))/totalpix_anim)*100;
    clear ud udx udy
end

%% Calculate Shades
Th=[];
R=[];
Per=[];
Lvalues=[];
Colors=[];
Colors2=[];
Names=[];
Base=[];
Deg=[];
Vectorplot=[];


for i=1:length(Regions_Angles)
    colorcode=Regions_Angles(i).Code;
    colorname={Regions_Angles(i).Name};
    base=Regions_Angles(i).base;
    im=Regions_Angles(i).ResultingRGB;
    im_mask=Regions_Angles(i).Resultingmask;
    [X_no_dither,map] = rgb2ind(im,numshades,'nodither');
    maskper=Regions_Angles(i).Resultingmask_percentage;
    if maskper==100
        X_no_dither=ones(size(X_no_dither));
    end
    xnological=logical(X_no_dither);
    totalpix=sum(xnological(:));
    
    if height(map)>1
     map(1,:)=[];
    end

    Regions_Angles(i).X_no_dither=X_no_dither;
    Regions_Angles(i).map=map;

    for j=1:numshades
        maski=X_no_dither==j;
        %          figure; imshow(maski)
        percent0=(sum(maski(:))/totalpix)*100;
        
        if percent0>5
            %             figure; imshow(maski)
            shades(j).mask=maski;
            shades(j).shade_RGB=(map(j,:));
            shades(j).shade_LAB=rgb2lab(map(j,:));
            shades(j).shade_L=shades(j).shade_LAB(1);
            [shades(j).theta,shades(j).rho] = cart2pol(shades(j).shade_LAB(3),shades(j).shade_LAB(2));
            shades(j).thetadeg=wrapTo360(rad2deg( shades(j).theta));
            percent=(sum(maski(:))/totalpix_anim)*100;
            shades(j).percent=percent;
        end
    end
    if exist('shades')
        shades = shades(all(~cellfun(@isempty,struct2cell(shades))));
        Regions_Angles(i).shades=shades;
        Th=[Th;[shades.theta].'];
        R=[R;[shades.rho].'];
        Per=[Per;[shades.percent].'];
        Lvalues=[Lvalues;[shades.shade_L].'];
        Deg=[Deg;[shades.thetadeg].'];
        Colors=[Colors;cell2mat({shades.shade_RGB}.')];
        Colors2=[Colors2;repmat(colorcode,[length(shades) 1])];
        Names=[Names;repmat(colorname,[length(shades) 1])];
        Base=[Base;repmat(base,[length(shades) 1])];
        Vectorplot=[Vectorplot ((-1).^(0:length(shades)-1).*2)];

        clear shades
    end

end

Vectorplot=[Vectorplot].';

if length(Th)<2
tb0 = table(Th,R,Lvalues,Deg,Per,Colors,cellstr(Colors2),Names,Base,Vectorplot,'VariableNames', {'Th', 'R','Lvalues','Deg','Per','Colors','Colors2','Names','Base','Vectorplot'});

else

tb0 = table(Th,R,Lvalues,Deg,Per,Colors,Colors2,Names,Base,Vectorplot,'VariableNames', {'Th', 'R','Lvalues','Deg','Per','Colors','Colors2','Names','Base','Vectorplot'});

end
%% Organize the Shades
tb0 = sortrows(tb0,'Colors2','ascend');
[sortedcolor,sortIdx]=sortrows(cellstr( Colors2 ),'ascend');
[G,TID] = findgroups(sortedcolor);
[B, C]= unique(Names(sortIdx),'stable');

[G2,TID2] = findgroups(cellstr( Names(sortIdx)));
T_split = splitapply( @(varargin) varargin, tb0 , G);
subTables = cell(size(T_split, 1));
% for i = 1:size(T_split, 1)
%     if i==IndexLight
%         subTables{i} = tb0(IndexLight_true,1:6);
%     else
%         subTables{i} = table(T_split{i, :}, 'VariableNames', ...
%             tb0.Properties.VariableNames);
%     end
% end
Names2cell = table2cell(tb0(:,8));
for i = 1:size(T_split, 1)
    ttsplit_length=length(T_split{i,8});
    if ttsplit_length==1
        namesearch=T_split{i,8};
        idnamesearch=find(strcmp(Names2cell, namesearch{1, 1}));

        subTables{i} = tb0(idnamesearch,1:end);
    else
        subTables{i} = table(T_split{i, :}, 'VariableNames', ...
            tb0.Properties.VariableNames);
    end
end


lines=Regions_Angles(1).lines;
lines2=Regions_Angles(1).lines2;

%% Display results:
figure('WindowState','maximized','Color',[1 1 1]);
for p=1:length(Regions_Angles)
    t1=Regions_Angles(p).Theta1;
    t2=Regions_Angles(p).Theta2;
    r1=Regions_Angles(p).Rho1;
    r2=Regions_Angles(p).Rho2;
    [xmesh,ymesh] = meshgrid(min(t1,t2):0.01:max(t1,t2),r1:0.1:r2);
    xmesh = xmesh(:);
    ymesh = ymesh(:);
    polarscatter(xmesh,ymesh,'MarkerFaceColor',Regions_Angles(p).Code,'MarkerEdgeColor',Regions_Angles(p).Code,'LineWidth',0.5);
    hold on
end
legend({Regions_Angles.Name},'AutoUpdate','off','FontSize',14)
th1 = linspace(0,2*pi,50);
th2 = linspace(Regions_Angles(length(Regions_Angles)-1).Theta1,Regions_Angles(length(Regions_Angles)-1).Theta2);
r1 = 3;
r2 = 35;
polarplot(th1,r1+zeros(size(th1)),'k')
hold on
polarplot(th2,r2+zeros(size(th2)),'k')
hold on
hold on
polarplot(lines, lines2*98,'k','LineWidth',0.5)
hold on
for i=1:length(subTables)
    tb1=subTables{i, 1};
    s = polarscatter(tb1,'Th','R','filled','ColorVariable','Colors');
    s.SizeData = 100;
    hold on
end
set(gcf,'color','w');

%%
figure('WindowState','maximized','Color',[1 1 1]);
for p=1:length(Regions_Angles)
    t1=Regions_Angles(p).Theta1;
    t2=Regions_Angles(p).Theta2;
    r1=Regions_Angles(p).Rho1;
    r2=Regions_Angles(p).Rho2;
    [xmesh,ymesh] = meshgrid(min(t1,t2):0.01:max(t1,t2),r1:0.1:r2);
    xmesh = xmesh(:);
    ymesh = ymesh(:);
    polarscatter(xmesh,ymesh,'MarkerFaceColor',Regions_Angles(p).Code,'MarkerEdgeColor',Regions_Angles(p).Code,'LineWidth',0.5);
    hold on
end
legend({Regions_Angles.Name},'AutoUpdate','off','FontSize',14)
polarplot(lines, lines2*98,'k','LineWidth',0.5)
hold on
for i=1:length(subTables)
    tb1=subTables{i, 1};
    s = polarbubblechart(tb1,'Th','R','Per','Colors','MarkerFaceAlpha',1);
    hold on
end
th1 = linspace(0,2*pi,50);
th2 = linspace(Regions_Angles(length(Regions_Angles)-1).Theta1,Regions_Angles(length(Regions_Angles)-1).Theta2);
r1 = 3;
r2 = 35;
polarplot(th1,r1+zeros(size(th1)),'k')
hold on
polarplot(th2,r2+zeros(size(th2)),'k')
hold on
bubblelegend('Percentage of Shade in Image','FontSize',14)
hold on
polarscatter([Regions_Angles.base],[Regions_Angles.Rho1],'k*');
set(gcf,'color','w');

%% Show mixing
Intervals=Regions_Angles1(1).AllIntervals;
Intervals_rad=Regions_Angles1(1).AllIntervals_Rad;
Intervals_rad1 = [Intervals_rad.x1; Intervals_rad.x2].';
[Bist,~,Y]=unique(Intervals_rad1,'rows','stable');
[C,X] = hist(Y,unique(Y));
Z = find(ismember(Y,X(C>1)));

for i=1:length(subTables)
    tb1=subTables{i, 1};
    for j=1:height(tb1)
        iswithbrown=0;
        th=wrapTo360(rad2deg( tb1.Th(j)));
        R=tb1.R(j);
        belongs=[];
        for k=1:length(Intervals)
            if th<Intervals(k).x2 && th>Intervals(k).x1
                belongs(k)=1;
            else
                belongs(k)=0;
            end
        end
        idxs=find(belongs==1).';
        idxs=sort(idxs,'descend');
        bcolors1=Intervals(idxs(1)).Color;
        if strcmp(bcolors1,tb1.Names(j))
            idxs=idxs;
        else
            idxs=sort(idxs,'ascend');
            bcolors1=Intervals(idxs(1)).Color;
        end
        val1=th*Intervals(idxs(1)).p1+Intervals(idxs(1)).p2;
        tb1.Val1(j)=val1;
        tb1.Val1_name(j)={bcolors1};
        if length(idxs)>1
            bcolors2=Intervals(idxs(2)).Color;
            totcolors=2;
            val2=th*Intervals(idxs(2)).p1+Intervals(idxs(2)).p2;
            tb1.Val2(j)=val2;
            tb1.Val2_name(j)={bcolors2};
            tb1.totcolors(j)=totcolors;
        else
            totcolors=1;
            tb1.totcolors(j)=totcolors;
        end
        % rad:
        if th<90 && th>0
            belongs_rad=[];
            % check only if it is in between brown and pure
            if R<Intervals_rad(Z(1)).x2 && R>Intervals_rad(Z(1)).x1
                iswithbrown=1;
                tb1.totcolors(j)=totcolors;
                bcolors1_rad=Intervals_rad(Z(1)).Color;
                bcolors2_rad=Intervals_rad(Z(2)).Color;
                val1_rad=R*Intervals_rad(Z(1)).p1+Intervals_rad(Z(1)).p2;
                val2_rad=R*Intervals_rad(Z(2)).p1+Intervals_rad(Z(2)).p2;
                if strcmp(bcolors1,'Pure')
                    tb1.Val_pure(j)=val1_rad;
                    tb1.Val_brown(j)=val2_rad;
                else
                    tb1.Val_pure(j)=val2_rad;
                    tb1.Val_brown(j)=val1_rad;
                end

            end

        end
        if iswithbrown==0
            if totcolors==1
                colorname=strcat(num2str(val1*100),'% of',{' '},bcolors1);
            elseif totcolors==2
                colorname=strcat(num2str(val1*100),'% of',{' '},Intervals(idxs(1)).Color,' and ',{' '}, num2str(val2*100),'% of',{' '},Intervals(idxs(2)).Color);
            end
        elseif iswithbrown==1
            if totcolors==1
                colorname=strcat(num2str(tb1.Val_brown(j)*100),'% of',{' '},' Brown, and',{' '}, num2str(tb1.Val_pure(j)*val1*100),'% of',{' '},Intervals(idxs(1)).Color);

            else
                colorname=strcat(num2str(tb1.Val_brown(j)*100),'% of',{' '},' Brown, ',{' '}, num2str(tb1.Val_pure(j)*val1*100),'% of',{' '},Intervals(idxs(1)).Color,', and ',{' '}, num2str(tb1.Val_pure(j)*val2*100),'% of',{' '},Intervals(idxs(2)).Color);
            end
        end
        tb1.ShadeName(j)=colorname;
        tb1.Angle(j)=th;

    end
    subTables{i, 1}=tb1;
end


figure('WindowState','maximized','Color',[1 1 1]);
t = tiledlayout('flow','TileSpacing','Compact');
nexttile
imshow(an_im,[])
for i=1:length(subTables)
    tb1=subTables{i, 1};
    nexttile
    s = bubblechart(tb1,'R','Angle','Per','Colors','MarkerFaceAlpha',1);
    title(strcat('The',{' '},tb1.Names(1),'s'),'FontSize',15,FontName='Arial')
    xlabel('Radius',FontSize=15,FontName='Arial')
    ylabel('Angle (º)',FontSize=15,FontName='Arial')
    % bubblelegend('Percentage of Shade in Image')
    nexttile
    imshow(Regions_Angles((find(strcmp({Regions_Angles.Name}.',tb1.Names(1))))).ResultingRGB)
    title(strcat('The Segmented',{' '},tb1.Names(1),'s'),'FontSize',12,FontName='Arial')

end



for i=1:length(subTables)
coloris=subTables{i, 1}.Names{1};
subTables{i, 2}=coloris;
end


end