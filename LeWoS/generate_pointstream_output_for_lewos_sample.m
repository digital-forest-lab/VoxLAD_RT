%this code reads a lewos classification and creates a pointstream file for
%a sample to vizualise probabilities using a color wheel

clear

plot='MorganMonroe'
path=strcat('E:\ReadRXP\voxTrace_',plot,'\');
file_lewos=strcat(path,'Morgan_nofloor_octree2cm_1streturn.txt_results.bin'); %this is the file created by LeWoS containing the point coordinates and probabilities of wood class
fid_lewos=fopen(file_lewos);

file_res_wl_sep_sample=strcat(path,'lewos_prob_sample.asc');
fid_res_prob=fopen(file_res_wl_sep_sample,'w');
distance_around_center_to_sample=7;

file_rxp=strcat(path,'ScanPos1.bin'); %this is the binary file created by C++ code containing pulse and hit events
fid_rxp=fopen(file_rxp);

bbox_x_min=fread(fid_rxp,1,'double');
bbox_y_min=fread(fid_rxp,1,'double');
bbox_z_min=fread(fid_rxp,1,'double');
bbox_x_max=fread(fid_rxp,1,'double');
bbox_y_max=fread(fid_rxp,1,'double');
bbox_z_max=fread(fid_rxp,1,'double');

plot_center=[bbox_x_min+(bbox_x_max-bbox_x_min)/2 bbox_y_min+(bbox_y_max-bbox_y_min)/2];

num_pts= fread(fid_lewos,1,'uint32')
points = fread(fid_lewos,[num_pts 3],'double');
prob_mb=fread(fid_lewos,[num_pts 1],'double');

for pt=1:num_pts
    
    if sqrt((points(pt,1)-plot_center(1))^2+(points(pt,2)-plot_center(2))^2)<distance_around_center_to_sample
        
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

fclose all;