function [pcloud, transforms]=reconstruction(image_name, depth_cam, rgb_cam, Rdtrgb,Tdtrgb)
    
    image_name=image_name(1:length(image_name));
    
    for im_pairs=1:length(image_name)-1
        
        %indicates in the command window what pair of images being analised
        display([image_name(im_pairs).rgb,' and ',image_name(im_pairs+1).rgb])
        
        %load depth and rgb of image 1 of the pair
        load(image_name(im_pairs).depth);
        depth1=depth_array;
        im1=imread(image_name(im_pairs).rgb);
        
        %load depth and rgb of image 2 of the pair
        load(image_name(im_pairs+1).depth);
        depth2=depth_array;
        im2=imread(image_name(im_pairs+1).rgb);
        
        %gets xyz from each images depths
        xyz1=get_xyzasus(depth1(:),[480 640],1:640*480,depth_cam.K,1,0);
        xyz2=get_xyzasus(depth2(:),[480 640],1:640*480,depth_cam.K,1,0);
        
        %gets color info of respective xyz points from each images 
        rgbd1=get_rgbd(xyz1,im1,Rdtrgb,Tdtrgb,rgb_cam.K);
        rgbd2=get_rgbd(xyz2,im2,Rdtrgb,Tdtrgb,rgb_cam.K);
        cl1=reshape(rgbd1,480*640,3);
        cl2=reshape(rgbd2,480*640,3);
        
        %grayscale of the 2 images
        im1g = single(rgb2gray(im1));
        im2g = single(rgb2gray(im2));
        
        %find the features of each image using SIFT
        [fa, da] = vl_sift(im1g) ;
        [fb, db] = vl_sift(im2g) ;
        
        %matching of each image features. Some matchings are correct (inliers) and
        %others aren't (outliers). There are features that are not matched
        %at all
        [matches, scores] = vl_ubcmatch(da, db) ;

        %% Ransac-to distinguish inliers from outliers
        %(u1,v1)-rgb coordinates of the 1st image features with match
        u1=fa(1,matches(1,:));
        v1=fa(2,matches(1,:));
        
        %(u2,v2)-rgb coordinates of 2nd image features whose index is equal
        %to index of its 1st image match
        u2=fb(1,matches(2,:));
        v2=fb(2,matches(2,:));
                        
        %from the (u,v), it gets the xyz index 
        indice1=480*floor(u1)+floor(v1);
        indice2=480*floor(u2)+floor(v2);
        
        %xyz of the matched features (inliers and outliers)
        points1=xyz1(indice1,:);
        points2=xyz2(indice2,:);
            
        %max number of inliers detected by RANSAC (initialization)
        inlier_max=0;
        iterations=0;
        
        %try 100 models
        while iterations<100
            iterations=iterations+1;
            %gets 4 random indexes(1 to length of matches)
            random_index=randsample(length(matches),4);
            
            %gets the points with the respective random indexes
            random_set1=points1(random_index,:);
            random_set2=points2(random_index,:);
            
            %calculates the centroid of each 4 points (1 per image) and the
            %vector distance between those points and their centroid
            centroid1=mean(random_set1);
            cdistance1=random_set1-repmat(centroid1,4,1);
            centroid2=mean(random_set2);
            cdistance2=random_set2-repmat(centroid2,4,1);
            
            %calculates covariance matrix for SVD
            Covariance=cdistance2'*cdistance1;

            [U,S,V] = svd(Covariance);
            
            %best R,T transformation with pairs of four points
            R=V*U';
            T=centroid1'-R*centroid2';
            
            %calculates the square error por every point (inlier and
            %outlier
            error=points1'-R*points2'-repmat(T,1,length(matches));
            error = sum(abs(error).^2);
            
            %it is boolean, if error is less than a threshold, it's
            %concidered an inlier for the suggested model(1). if not, it's an outlier (0)
            classifier=error<0.002;
            
            %number of the suggested model's inliers found
            n_inliers=sum(classifier);
            
            %if this model found more inliers than the previous ones
            if n_inliers>inlier_max
                inlier_max=n_inliers;
                
                %registers the boolean map with more inliers
                max_class=classifier;
            end
        end
      
        i=1;
        while i<=length(max_class)
            
            %if points in index i are outliers
            if max_class(i)==0
                
                %erase them from the matrix of points that are going to the
                %be used to find the R,T between images
                max_class(i)=[];
                points1(i,:)=[];
                points2(i,:)=[];
            else
                i=i+1;
            end
        end

        %% Calculate R,T with the inliers
        
        %calculates the centroids of the inliers and the vector distance
        %between the inliers and their respective centroid
        centroid1=mean(points1);
        cdistance1=points1-repmat(centroid1,length(points1),1);
        centroid2=mean(points2);
        cdistance2=points2-repmat(centroid2,length(points2),1);
        
        %Coveariance matrix for SVD
        Covariance=cdistance2'*cdistance1;

        [U,S,V] = svd(Covariance);
        
        %Rigid Body Transformation from image 2 to image 1
        R=V*U';
        T=centroid1'-R*centroid2'; 

        %If the 1st pair of images are being analised
        %NOTE: transforms(i) has transformations of image i to the origin image (1st of all images)
        if im_pairs==1
            
            % Transformation of the origin image to origin image is R=I and T=0
            transforms(1)=struct('R',eye(3),'T',zeros(3,1));
            
            %Transformation of the 2nd (of all) image to the origin is R=R and T=T 
            transforms(im_pairs+1)=struct('R',R,'T',T);
        else
            %Chain: R=previous_R*R ; T=previous_R*T+previous_T
            transforms(im_pairs+1)=struct('R',transforms(im_pairs).R*R,'T',transforms(im_pairs).R*T+transforms(im_pairs).T);
        end
        
        %apply transformation to xyz of image 2
        xyz2=transforms(im_pairs+1).R*xyz2'+repmat(transforms(im_pairs+1).T,1,480*640);
        xyz2=xyz2';
        
        %if it's the 1st pair of image being analised
        if im_pairs==1
            
            %XYZ will have the xyz of all images
            XYZ=xyz1;
            
            %CL will have the xyz of all images
            CL=cl1;
        end
        
        xyz1=transforms(im_pairs).R*xyz1'+repmat(transforms(im_pairs).T,1,480*640);
        xyz1=xyz1';
        p1=pointCloud(xyz1,'Color',cl1);
        p1=pcdownsample(p1,'gridAverage',0.01);
        %point cloud of image 2 transformed
        p2=pointCloud(xyz2,'Color',cl2);
        p2=pcdownsample(p2,'gridAverage',0.01);
        %ICP to adjust p2 to p1
        [tform]=pcregrigid(p2,p1,'InlierRatio',inlier_max/307200);
        
        %Update the transforms with the new transformation
        transforms(im_pairs+1).R=tform.T(1:3,1:3)'*transforms(im_pairs+1).R;
        transforms(im_pairs+1).T=tform.T(1:3,1:3)'*transforms(im_pairs+1).T+tform.T(4,1:3)';
        
        xyz2=tform.T(1:3,1:3)*xyz2'+repmat(tform.T(4,1:3)',1,480*640);
        %XYZ and CL have the updated information of the xyz of the image 2 
        XYZ=[XYZ;xyz2'];
        CL=[CL;cl2];
        
        %filters the xyz points. If there are points that are less than a
        %milimeter apart from each other, than those points are considered
        %the same point
        [XYZ,iold,inew]=unique(fix(XYZ*1000),'rows');
        XYZ=XYZ/1000;
        CL=CL(iold,:);
        
    end
    
    %shows the final point cloud
    p=pointCloud(XYZ,'Color',CL);
    showPointCloud(p);
    pcloud=[double(XYZ),double(CL)];
   
end