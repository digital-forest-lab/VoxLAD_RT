%this code is used to combine LAD and WAD arrays into one single array
%formated for FLiESvox

clear

voxel_size=0.3;
min_dist=voxel_size-(voxel_size/2);

plot='SERC'
file_LAD=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot,'/',plot,'_5mgrid_LAD_no_zeros.asc_occlusion_corr_thresh_30with_omega_v.asc'); %this is the file created by LeWoS containing the point coordinates and probabilities of wood class
file_WAD=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot,'/',plot,'_5mgrid_WAD_no_zeros.asc_occlusion_corr_thresh_15.asc');
%file_LAD=strcat('F:\ReadRXP\voxTrace_',plot,'\',plot,'_5mgrid_LAD_no_zeros.asc_occlusion_corr_thresh_30with_omega_v.asc'); %this is the file created by LeWoS containing the point coordinates and probabilities of wood class
%file_WAD=strcat('F:\ReadRXP\voxTrace_',plot,'\',plot,'_5mgrid_WAD_no_zeros.asc_occlusion_corr_thresh_15.asc');
fid_LAD=fopen(file_LAD);
fid_WAD=fopen(file_WAD);

% path=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where bin files are
% path_res=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where result files are written
%path=strcat('E:\',plot,'\'); %this is where bin files are
%path_res=strcat('G:\',plot,'\'); %this is where result files are written

% start_pos=1;
% end_pos=242;

file_res=strcat(file_LAD,'_combined_WAD.asc');
fidres=fopen(file_res,'w');
% file_res_wl_sep_sample=strcat(path_res,'wood_leaf_sample.asc');
% fidres_wl_sample=fopen(file_res_wl_sep_sample,'w');
% distance_around_center_to_sample=5; %this is the area that will be output to vizualise the wood leaf separation, in meters

%wood_probability_threshold=0.9;



%fid_lewos=fopen(file_lewos);

disp('Reading LAD array and creating KDTree_object');
tic

%Do a first loop to put x y z coord from LAD in a KDtree
i=1;
while ~feof(fid_LAD)
    xtemp=fscanf(fid_LAD, '%f', [1]);
    if isempty(xtemp), break, end
    xmin(i)=xtemp;
    ymin(i)=fscanf(fid_LAD, '%f', [1]);
    zmin(i)=fscanf(fid_LAD, '%f', [1]);
    
    xmax=fscanf(fid_LAD, '%f', [1]);
    ymax=fscanf(fid_LAD, '%f', [1]);
    zmax=fscanf(fid_LAD, '%f', [1]);
    
    lad(i)=fscanf(fid_LAD, '%f', [1]);
    omega_v(i)=fscanf(fid_LAD, '%f', [1]);
    g_code(i)=fscanf(fid_LAD, '%d', [1]);
    i=i+1;
    
end
coords(:,1)=xmin;
coords(:,2)=ymin;
coords(:,3)=zmin;
coords_LAD_KD=KDTreeSearcher(coords);
toc

%Do a second loop to go through WAD and print
disp('Going through WAD array and outputing: wood and leaf, and wood only');
tic
clear xmin ymin zmin
i=1;
while ~feof(fid_WAD)
    xtemp=fscanf(fid_WAD, '%f', [1]);
    if isempty(xtemp), break, end
    xmin(i)=xtemp;
    ymin(i)=fscanf(fid_WAD, '%f', [1]);
    zmin(i)=fscanf(fid_WAD, '%f', [1]);
    
    xmax=fscanf(fid_WAD, '%f', [1]);
    ymax=fscanf(fid_WAD, '%f', [1]);
    zmax=fscanf(fid_WAD, '%f', [1]);
    
    wad=fscanf(fid_WAD, '%f', [1]);
    coords_WAD(1)=xmin(i);
    coords_WAD(2)=ymin(i);
    coords_WAD(3)=zmin(i);
    idx=knnsearch(coords_LAD_KD,coords_WAD);
    
    %if the coords are the same, then it is the same voxel which has both leaf and wood,
    %otherwise, this is a voxel which has wood and no leaf
    distance=pdist2 (coords_WAD,coords(idx,:));
   
    if distance<min_dist
%     if coords_WAD==coords(idx,:)
%      
%         %voxel has both wood and leaves
        if lad(idx)>0.0005 || wad >0.0005, fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %.3f 0 0 0 %.3f\n',xmin(i),ymin(i),zmin(i),xmax,ymax,zmax,lad(idx),g_code(idx),wad,omega_v(idx));end
    else
        
        %voxel has only wood (lad=0)
        if wad >0.0005,fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f 0.0 1 %.3f 0 0 0 1.0\n',xmin(i),ymin(i),zmin(i),xmax,ymax,zmax,wad);end
        
    end
    i=i+1;
end
toc
tic
disp('Creating WAD KDTree_object');
coords_WAD_(:,1)=xmin;
coords_WAD_(:,2)=ymin;
coords_WAD_(:,3)=zmin;
coords_WAD_KD=KDTreeSearcher(coords_WAD_);

toc

%Do a third loop through LAD to print those voxels that have no WAD
disp('Going through LAD array and outputing: leaf only');
tic
clear lad omega_v g_code
frewind(fid_LAD);

while ~feof(fid_LAD)
    xtemp=fscanf(fid_LAD, '%f', [1]);
    if isempty(xtemp), break, end
    coords_LAD_(1)=xtemp;
    coords_LAD_(2)=fscanf(fid_LAD, '%f', [1]);
    coords_LAD_(3)=fscanf(fid_LAD, '%f', [1]);
    
    xmax=fscanf(fid_LAD, '%f', [1]);
    ymax=fscanf(fid_LAD, '%f', [1]);
    zmax=fscanf(fid_LAD, '%f', [1]);
    
    
    lad=fscanf(fid_LAD, '%f', [1]);
    omega_v=fscanf(fid_LAD, '%f', [1]);
    g_code=fscanf(fid_LAD, '%d', [1]);
    %     coords_LAD_(1)=xmin;
    %     coords_LAD_(2)=ymin;
    %     coords_LAD_(3)=zmin;
    idx=knnsearch(coords_WAD_KD,coords_LAD_);
    
     distance=pdist2 (coords_LAD_,coords_WAD_(idx,:));
      %if the coords are the same, then it is the same voxel which has both leaf and wood,
    %otherwise, this is a voxel which has leaf and no wood
    
    %if coords_LAD_~=coords_WAD_(idx,:)
    if distance>min_dist
        %voxel has leaves and no wood inside (WAD =0)
        
        if lad>0.0005, fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f %.3f %d 0.0 0 0 0 %.3f\n',coords_LAD_(1),coords_LAD_(2),coords_LAD_(3),xmax,ymax,zmax,lad,g_code,omega_v);end
  
    end
end


toc
fclose all;