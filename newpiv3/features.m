tic
clear all
close all


K=[525 0 319.5;
    0 525 239.5;
    0 0 1];

depth=dir('depth*.mat');
images=dir('*.png');

for im_pairs=1:length(depth)-2

load(depth(im_pairs).name)
depth1=depth_array;
im1=imread(images(im_pairs).name);

load(depth(im_pairs+1).name)
depth2=depth_array;
im2=imread(images(im_pairs+1).name);

cl1=reshape(im1,480*640,3);
cl2=reshape(im2,480*640,3);

xyz1=zeros(3, length(depth1(:)));
xyz2=zeros(3, length(depth2(:)));

u = repmat(1:640,480,1);
u = u(:)-K(1,3);
v = repmat((1:480)',640,1);
v=v(:)-K(2,3);
xyz1=zeros(length(u),3);
xyz2=zeros(length(u),3);


xyz1(:,3) = double(depth1(:))*0.001; % Convert to meters
xyz1(:,1) = (xyz1(:,3)/K(1,1)) .* u ;
xyz1(:,2) = (xyz1(:,3)/K(2,2)) .* v;
xyz2(:,3) = double(depth2(:))*0.001; % Convert to meters
xyz2(:,1) = (xyz2(:,3)/K(1,1)) .* u ;
xyz2(:,2) = (xyz2(:,3)/K(2,2)) .* v;


 im1g = single(rgb2gray(im1));
 im2g = single(rgb2gray(im2));
 [fa, da] = vl_sift(im1g) ;
 [fb, db] = vl_sift(im2g) ;
 [matches, scores] = vl_ubcmatch(da, db) ;
 
%% Ransac
u1=zeros(1,length(matches));
v1=zeros(1,length(matches));
u2=zeros(1,length(matches));
v2=zeros(1,length(matches));

for i=1:length(matches)
    u1(i)=fa(1,matches(1,i));
    v1(i)=fa(2,matches(1,i));
    u2(i)=fb(1,matches(2,i));
    v2(i)=fb(2,matches(2,i));    
end

pontos1=zeros(length(matches),3);
pontos2=zeros(length(matches),3);
for i=1:length(matches)
        indice1(i)=480*floor(u1(i))+floor(v1(i));
        indice2(i)=480*floor(u2(i))+floor(v2(i));

        pontos1(i,:)=xyz1(indice1(i),:);
        pontos2(i,:)=xyz2(indice2(i),:)';
end

p1=pointCloud(xyz1,'Color',cl1);
p2=pointCloud(xyz2,'Color',cl2);

inlier_max=0;
iterations=0;
while iterations<100
    iterations=iterations+1;
    random_index=randsample(length(matches),4);
    random_set1=zeros(4,3);
    random_set2=zeros(4,3);
    for i=1:4
       random_set1(i,:)=pontos1(random_index(i),:);
       random_set2(i,:)=pontos2(random_index(i),:);
    end
    
    centroid1=mean(random_set1);
    cdistancia1=random_set1-repmat(centroid1,4,1);
    centroid2=mean(random_set2);
    cdistancia2=random_set2-repmat(centroid2,4,1);

    Covariance=cdistancia2'*cdistancia1;

    [U,S,V] = svd(Covariance);

    R=V*U';
    T=centroid1'-R*centroid2';
    
    error=pontos1'-R*pontos2'-repmat(T,1,length(matches));
    error = sqrt(sum(abs(error).^2));  
    classifier=error<0.015;  
    n_inliers=sum(classifier);
    if n_inliers>inlier_max
       inlier_max=n_inliers;
       max_class=classifier;
    end
end
inlier_max
classifier=max_class;
i=1;
while i<=length(classifier)
    if classifier(i)==0
        classifier(i)=[];
        pontos1(i,:)=[];
        pontos2(i,:)=[];
    else
        i=i+1;
    end
end

%% 

    
    centroid1=mean(pontos1);
    cdistancia1=pontos1-repmat(centroid1,length(pontos1),1);
    centroid2=mean(pontos2);
    cdistancia2=pontos2-repmat(centroid2,length(pontos2),1);

    Covariance=cdistancia2'*cdistancia1;

    [U,S,V] = svd(Covariance);

    R=V*U';
    T=centroid1'-R*centroid2';
    
    if im_pairs==1
        transforms(1)=struct('R',eye(3),'T',zeros(3,1));
        transforms(im_pairs+1)=struct('R',R,'T',T);
    else
        transforms(im_pairs+1)=struct('R',transforms(im_pairs).R*R,'T',transforms(im_pairs).R*T+transforms(im_pairs).T);
    end
    
    xyz2=transforms(im_pairs+1).R*xyz2'+repmat(transforms(im_pairs+1).T,1,480*640);
    xyz2=xyz2';
cl1=reshape(im1,480*640,3);
cl2=reshape(im2,480*640,3);
 
p1=pointCloud(xyz1,'Color',cl1);
p2=pointCloud(xyz2,'Color',cl2);

figure(1);
if im_pairs==1
    %showPointCloud(p1);
    cloud{1}=p1
    XYZ=xyz1;
    CL=cl1;
end

 cloud{im_pairs+1}=p2
[tform,cloud{im_pairs+1}]=pcregrigid(cloud{im_pairs+1},cloud{im_pairs},'InlierRatio',inlier_max/307200);
transforms(im_pairs+1).R=tform.T(1:3,1:3)'*transforms(im_pairs+1).R;
transforms(im_pairs+1).T=tform.T(1:3,1:3)'*transforms(im_pairs+1).T+tform.T(4,1:3)';
XYZ=[XYZ;p2.Location];
CL=[CL;cl2];
% showPointCloud(cloud{im_pairs+1});
[XYZ,iold,inew]=unique(fix(XYZ*1000),'rows');
XYZ=XYZ/1000;

CL=CL(iold,:);
end
p=pointCloud(XYZ,'Color',CL);
showPointCloud(p);
toc