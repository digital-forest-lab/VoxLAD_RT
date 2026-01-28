%this code reads the point cloud and launches LeWoS
%Because this code requires loading the point cloud in memory, I may need
%to sprob_mbt the array in multiple parts and run on the 128Gb RAM Desktop
%the Computer Vision System Toolbox is requried, 29usd

clear
tic
%file_ptcloud='E:\temp riscan project for octree\Pasoh.RiSCAN\Pasoh_nofloor_octree2cm_1streturn.txt';
file_ptcloud='F:\ReadRXP\voxTrace_Bartlett\Bartlett_Pointcloud_for_lewos.txt';
%file_bin=strcat(file_ptcloud,'.bin');
%file_ptcloud='D:\Wang leaf wood separation model\LeWoS-master\sample_input_file_from_Moorthy.txt';

fid_ptcloud=fopen(file_ptcloud);
%fid_bin=fopen(file_bin);

output_results_pointstream=1; %in ascii pointstream format using color wheel 
output_results_probability=1; %in binary format

ft_threshold=0.125;
paral=1;
plot=0;

%__________end of input section

if output_results_pointstream
  %  file_ptstream_wl=strcat(file_ptcloud,'_wl.asc');
    file_ptstream_prob=strcat(file_ptcloud,'_prob.asc');
    
 %   fid_res_wl=fopen(file_ptstream_wl,'w');
    fid_res_prob=fopen(file_ptstream_prob,'w');
end


if output_results_probability
 file_res=strcat(file_ptcloud,'_LeWoS_results.bin');
 fid_res=fopen(file_res,'w');
end

num_pts=0;
while ~feof(fid_ptcloud)
    num_pts=num_pts+1;
    %xyz=fscanf(fid_bin, '%f', [3]);
    fgetl(fid_ptcloud);
    
    %if isempty(bidon), continue, end
end
num_pts
points=zeros(num_pts,3);
pt_num=0;
frewind(fid_ptcloud);
%while ~feof(fid_ptcloud)
for pt_num=1:num_pts
    points(pt_num,:)=fscanf(fid_ptcloud, '%f', [3]);
%     bidon=fscanf(fid_ptcloud, '%f', [1]);
%     if isempty(bidon), break, end
%     pt_num=pt_num+1;
%     points(pt_num,1)=bidon;
%     
%     points(pt_num,2)=fscanf(fid_ptcloud, '%f', [1]);
%     points(pt_num,3)=fscanf(fid_ptcloud, '%f', [1]);
    
end
toc
tic
message="starting LeWoS"
[prob_mb] = RecursiveSegmentation_release(points, ft_threshold, paral, plot);
toc
tic
message="writing results to file"

% if output_results_probability
%      for pt=1:size(prob_mb)
%         fprintf(fid_res, '%f %f %f %f\n',points(pt,1),points(pt,2),points(pt,3),prob_mb(pt));
%      end
% end
num_pts
size_points=size(points)
size_prob_mb=size(prob_mb)
if output_results_probability
    fwrite(fid_res,size(points,1),'uint32');
    fwrite(fid_res,points,'double');
    fwrite(fid_res,prob_mb,'double');
end
    


% write results in file
if output_results_pointstream
    for pt=1:size(prob_mb)
        
%         if BiLabel_Regu(pt)==1
%             fprintf(fid_res_wl, '%f %f %f 0 0 0\n',points(pt,1),points(pt,2),points(pt,3));
%         else
%             fprintf(fid_res_wl, '%f %f %f 128 128 128\n',points(pt,1),points(pt,2),points(pt,3));
%         end
        if prob_mb(pt)==0
            fprintf(fid_res_prob, '%f %f %f 255 0 255\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.1
            fprintf(fid_res_prob, '%f %f %f 125 0 255\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.15
            fprintf(fid_res_prob, '%f %f %f 0 0 255\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.2
            fprintf(fid_res_prob, '%f %f %f 0 125 255\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.3
            fprintf(fid_res_prob, '%f %f %f 0 255 255\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.4
            fprintf(fid_res_prob, '%f %f %f 0 255 125\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.5
            fprintf(fid_res_prob, '%f %f %f 0 255 0\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.6
            fprintf(fid_res_prob, '%f %f %f 125 255 0\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.7
            fprintf(fid_res_prob, '%f %f %f 255 255 0\n',points(pt,1),points(pt,2),points(pt,3));
        elseif prob_mb(pt)<0.8
            fprintf(fid_res_prob, '%f %f %f 255 125 0\n',points(pt,1),points(pt,2),points(pt,3));
        else
            fprintf(fid_res_prob, '%f %f %f 255 0 0\n',points(pt,1),points(pt,2),points(pt,3));
        end
    end
end
toc
fclose all;