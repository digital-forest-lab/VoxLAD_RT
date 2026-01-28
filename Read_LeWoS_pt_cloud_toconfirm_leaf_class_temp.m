%this code is used to get the number of points in a LeWoS classification to
%assign the array size in WooTrace to save processing time

% because this code creates a wood point cloud, the work cannot be split
% manually between workers, unless the wood point cloud files are stiched
% together afterwards

clear 


plot='Bartlett'
file_lewos=strcat('F:\ReadRXP\voxTrace_',plot,'\',plot,'_Pointcloud_for_lewos.txt_LeWoS_results.bin'); %this is the file created by LeWoS containing the point coordinates and probabilities of wood class

% path=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where bin files are
% path_res=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where result files are written
path=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where bin files are
path_res=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where result files are written

% start_pos=1;
% end_pos=98;
worker=8;

%***************

switch worker
    case 1
        start_pos=1;
        end_pos=12;
    case 2
        start_pos=13;
        end_pos=25;
    case 3
        start_pos=26;
        end_pos=38;
    case 4
        start_pos=39;
        %end_pos=124;
        end_pos=51; %for 7 x 7 grid at Bartlett
    case 5
        start_pos=52;
        end_pos=64;
    case 6
        start_pos=65;
        end_pos=77;
    case 7
        start_pos=78;
        end_pos=90;
    case 8
        start_pos=91;
        end_pos=98;
end


%file_res_wood=strcat(path_res,'wood_point_cloud.bin');
%fidres_wood=fopen(file_res_wood,'w');
file_res_wl_sep_sample=strcat(path_res,'wood_leaf_sample_worker_',num2str(worker),'.asc');
fidres_wl_sample=fopen(file_res_wl_sep_sample,'w');
distance_around_center_to_sample=5; %this is the area that will be output to vizualise the wood leaf separation, in meters

wood_probability_threshold=0.9;
wood_probability_threshold_lower=0.1;

tic

fid_lewos=fopen(file_lewos);

disp('reading lewos and creating KDTree_object');
tic
num_pts= fread(fid_lewos,1,'uint32')
lewos_pts = fread(fid_lewos,[num_pts 3],'double');
lewos_prob=fread(fid_lewos,[num_pts 1],'double');

lewos_pts_=KDTreeSearcher(lewos_pts);
toc

%get an estimate of the number of leaves in a point cloud to allocate
%memory, use a scan in the plot center
pos_center=floor((end_pos-start_pos)/2);
if rem(pos_center,2)==0, pos_center=pos_center-1;end
file_rxp=strcat(path,'ScanPos',num2str(pos_center),'.bin');
fid_rxp=fopen(file_rxp);
bidon=fread(fid_rxp,6,'double');
 num_leaves=0;
 num_wood=0;
    disp('getting number of leaves in plot center');
    tic
    while ~feof(fid_rxp)
        
        %event_a=fscanf(fid_rxp_ascii, '%d',[1]);
        event=fread(fid_rxp,1,'int');
        if isempty(event), break, end
        switch event
            case 0
                bidon=fread(fid_rxp,6,'double');
            case 1 %if event==1
                num_leaves=num_leaves+1;
                %leaves(leaf_num,:)=fread(fid_rxp,3,'double');
                xyz=fread(fid_rxp,3,'double');
                %             if rem(num_leaves,100000)==0
                %                 num_leaves
                %                 toc
                %                 tic
                %             end
            case 2 % elseif event>1   %point is wood or noise 
                num_wood=num_wood+1;
                xyz=fread(fid_rxp,3,'double');
            case 3
                xyz=fread(fid_rxp,3,'double');
        end
        
    end
    toc
    fclose(fid_rxp);
    num_leaves
    leaves=zeros(num_leaves*10,3); %allocating 10 times the amount of leaves found in the plot center scan, if this is not the maximum number of leaves found throughout all scans there will be a significant loss in efficiency
 num_wood
    woods=zeros(num_wood*10,3); %allocating 10 times the amount of wood found in the plot center scan, if this is not the maximum number of wood found throughout all scans there will be a significant loss in efficiency

for pos=start_pos:end_pos
    
    pos
    tic
    if pos==206 && plot=="SERC", continue,end
    
    file_rxp=strcat(path,'ScanPos',num2str(pos),'.bin'); %this is the binary file created by C++ code containing pulse and hit events
    file_rxp_res=strcat(path_res,'ScanPos',num2str(pos),'lewos.bin'); %This is the output of this code containing the same info as above but with modified event number for some of the leaf hits
    %file_rxp_ascii=strcat(path,'test_ascii');
    
    fid_rxp=fopen(file_rxp);
    fidres_rxp=fopen(file_rxp_res,'w');
    %fid_rxp_ascii=fopen(file_rxp_ascii);
    
    disp('storing leaves');
    
    bbox_x_min=fread(fid_rxp,1,'double');
    bbox_y_min=fread(fid_rxp,1,'double');
    bbox_z_min=fread(fid_rxp,1,'double');
    bbox_x_max=fread(fid_rxp,1,'double');
    bbox_y_max=fread(fid_rxp,1,'double');
    bbox_z_max=fread(fid_rxp,1,'double');
    
    plot_center=[bbox_x_min+(bbox_x_max-bbox_x_min)/2 bbox_y_min+(bbox_y_max-bbox_y_min)/2];
    
    %count to number of leaves to allocate array space and avoid having
    %zeros in array which would be considered in knnsearch
   
    %this next loop is to record the point xyz for all leaves in current
    %classification based on reflectance
    leaf_num=0;
    wood_num=0;
    num_lines=0;
    tic
    while ~feof(fid_rxp)
        
        %event_a=fscanf(fid_rxp_ascii, '%d',[1]);
        event=fread(fid_rxp,1,'int');
        if isempty(event), break, end
        num_lines=num_lines+1;
        switch event %if event==0
            case 0
                bidon=fread(fid_rxp,6,'double');
            case 1
                leaf_num=leaf_num+1;
                leaves(leaf_num,:)=fread(fid_rxp,3,'double');
                %             if rem(leaf_num,100000)==0
                %                 leaf_num
                %                 toc
                %                 tic
                %             end
            case 2  %elseif event>1 %point is wood or noise
                wood_num=wood_num+1;
                woods(wood_num,:)=fread(fid_rxp,3,'double');
            case 3
                xyz=fread(fid_rxp,3,'double');
        end
    end
    toc
    
    %searching for array indexes in lewos for points corresponding to leaf
    %in relfetance based classif
    disp('looking for lewos match');
    tic
    %leaves=gpuArray(leaves); %using large array on gpu does not work, a
    %limitation is set at size(uint32)
    idx=knnsearch(lewos_pts_,leaves(1:leaf_num,:)); %idx in position 1 contains the index of the point in lewos_pts that is closest to the firest leaf
     idx_wood=knnsearch(lewos_pts_,woods(1:wood_num,:)); %idx in position 1 contains the index of the point in lewos_pts that is closest to the firest leaf
  
    toc
    %creating a new binary file with modified classes
    frewind(fid_rxp);
    bidon=fread(fid_rxp,6,'double');
    fwrite(fidres_rxp,bidon,'double');
    disp('writing new bin file');
    tic
    leaf_num=0;
    wood_num=0;
    %while ~feof(fid_rxp)
    for i=1:num_lines
        
        event=fread(fid_rxp,1,'int');
        %if isempty(event), break, end
        
        switch event
            case 0  %if event==0
                bidon=fread(fid_rxp,6,'double');
                fwrite(fidres_rxp,event,'int');
                fwrite(fidres_rxp,bidon,'double');
            case 1  %elseif event==1
                xyz=fread(fid_rxp,3,'double');
                leaf_num=leaf_num+1;
                if lewos_prob(idx(leaf_num))>wood_probability_threshold
                    %point is wood in lewos, revert class to wood
                    fwrite(fidres_rxp,2,'int');
                    %fwrite(fidres_wood,xyz,'double');
                    if sqrt((xyz(1)-plot_center(1))^2+(xyz(2)-plot_center(2))^2)<distance_around_center_to_sample
                        fprintf(fidres_wl_sample,'%.3f %.3f %.3f 255 0 0\n',xyz(1),xyz(2),xyz(3));
                    end
                else %leaf class is maintained
                    fwrite(fidres_rxp,event,'int');
                    if sqrt((xyz(1)-plot_center(1))^2+(xyz(2)-plot_center(2))^2)<distance_around_center_to_sample
                        fprintf(fidres_wl_sample,'%.3f %.3f %.3f 0 255 0\n',xyz(1),xyz(2),xyz(3));
                    end
                end
                fwrite(fidres_rxp,xyz,'double');
            case 2  %elseif event==2
                xyz=fread(fid_rxp,3,'double');
                wood_num=wood_num+1;
                if lewos_prob(idx_wood(wood_num))<wood_probability_threshold_lower
                    %point is wood in intensity, but lower than threshold
                    %in Lewos, revert to noise since class is too uncertain
                    fwrite(fidres_rxp,3,'int');
                    %fwrite(fidres_wood,xyz,'double');
                    if sqrt((xyz(1)-plot_center(1))^2+(xyz(2)-plot_center(2))^2)<distance_around_center_to_sample
                        fprintf(fidres_wl_sample,'%.3f %.3f %.3f 0 0 255\n',xyz(1),xyz(2),xyz(3));
                    end
                else %wood class is maintained
                    fwrite(fidres_rxp,event,'int');
                     if sqrt((xyz(1)-plot_center(1))^2+(xyz(2)-plot_center(2))^2)<distance_around_center_to_sample
                        fprintf(fidres_wl_sample,'%.3f %.3f %.3f 0 0 0\n',xyz(1),xyz(2),xyz(3));
                    end
                end
                fwrite(fidres_rxp,xyz,'double');
            case 3  %else %if event==3
                xyz=fread(fid_rxp,3,'double');
                fwrite(fidres_rxp,event,'int');
                fwrite(fidres_rxp,xyz,'double');
        end
        
    end
    toc
end
toc
fclose all;