%this code generates a grid of origin points and directions for needles originating
%above a voxel array towards the ground

%-------input
%plot="EMS";
grid_spacing=0.3; %meters
distance_above_canopy=10; %meters
%plot_list=["EMS" "SERC" "Pasoh" "MorganMonroe"];
plot_list=["Hubbard"];
zen_angle_list=[15 30 45 60];
%zen_angle_list=[0];
buffer_on_plot_edges=1; % in meters
adjustement_for_unit_vector=0.001;

%-------end of input
for plot=plot_list
    for zenith_angle=zen_angle_list %degrees
        %path=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot','/');
        
        if plot=="Random"
            file_grid=strcat('H:\random array relative variance\needle_grid_theta_',num2str(zenith_angle),'_2.txt');
            bbox_x_min=0;
            bbox_y_min=0;
            bbox_z_min=0;
            bbox_x_max=10.2;
            bbox_y_max=10.2;
            bbox_z_max=10.2;
        else
            
            %path=strcat('H:\ReadRXP\voxTrace_',plot','\');
            path=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot','/');
           
            
            file_grid=strcat(path,'Relative Variance/needle_grid_theta_',num2str(zenith_angle),'.txt');
            %get bbox
            
            file_bbox=strcat(path,plot,'_readRXP_batch.bat');
            fid_bbox=fopen(file_bbox);
            fscanf(fid_bbox, '%s', [12]);
            
            bbox_x_min=fscanf(fid_bbox, '%d', [1]);
            bbox_y_min=fscanf(fid_bbox, '%d', [1]);
            bbox_z_min=fscanf(fid_bbox, '%d', [1]);
            bbox_x_max=fscanf(fid_bbox, '%d', [1]);
            bbox_y_max=fscanf(fid_bbox, '%d', [1]);
            bbox_z_max=fscanf(fid_bbox, '%d', [1]);
            
            fclose(fid_bbox);
        end
        
        fid_grid=fopen(file_grid,'w');
        zenith_angle=deg2rad(zenith_angle);
        
        grid_height=bbox_z_max+distance_above_canopy;
        x_unit_vector=0;
        y_unit_vector=sin(zenith_angle);
        z_unit_vector=-cos(zenith_angle);
        
        grid_shift_along_y_axis=-(grid_height-bbox_z_min)*tan(zenith_angle);
        
        test_last_line=rem((bbox_x_max-buffer_on_plot_edges)-(bbox_x_min+buffer_on_plot_edges),grid_spacing)==0;
        number_of_needles_x=ceil(((bbox_x_max-buffer_on_plot_edges)-(bbox_x_min+buffer_on_plot_edges))/grid_spacing)+test_last_line;
        test_last_line=rem((bbox_y_max-buffer_on_plot_edges+grid_shift_along_y_axis)-(bbox_y_min+buffer_on_plot_edges+grid_shift_along_y_axis),grid_spacing)==0;
        number_of_needles_y=ceil(((bbox_y_max-buffer_on_plot_edges+grid_shift_along_y_axis)-(bbox_y_min+buffer_on_plot_edges+grid_shift_along_y_axis))/grid_spacing)+test_last_line;
        number_of_needles=number_of_needles_x*number_of_needles_y;
        fprintf(fid_grid, '%d\n',number_of_needles);
        
        fprintf(fid_grid, '%.6f %.6f %.6f\n', x_unit_vector+adjustement_for_unit_vector,...
            y_unit_vector+adjustement_for_unit_vector,z_unit_vector+adjustement_for_unit_vector);
        i=0;
        for x=bbox_x_min+buffer_on_plot_edges:grid_spacing:bbox_x_max-buffer_on_plot_edges
            for y=bbox_y_min+buffer_on_plot_edges+grid_shift_along_y_axis:grid_spacing:bbox_y_max-buffer_on_plot_edges+grid_shift_along_y_axis
                fprintf(fid_grid, '%.2f %.2f %.2f\n', x,y,grid_height);
                i=i+1;
            end
        end
        if i~=number_of_needles, message='error on num of needles, revise text file first line',end
    end
end
fclose all;