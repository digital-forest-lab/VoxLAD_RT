% Compile the mex function
%run: mex wooRaytrace.cpp

%****This code only works for square floor matrices, there appears to be a
%problem with wooRatrace indexing the voxels when floor is not square.
%23 May 2020: this apparently was solve using a comment in
%https://www.mathworks.com/matlabcentral/fileexchange/56527-fast-raytracing-through-a-3d-grid
%There is a bug in the sub2ind function in the wooRayTrace.cpp
%Change: return (Y + (X - 1)*gridSize[0] + (Z - 1)*gridSize[1]*gridSize[0]);
%To: return (Y + (X - 1)*gridSize[1] + (Z - 1)*gridSize[1]*gridSize[0]);
%I though I had tried this before, but I may have understood wrongly that
%the issue was arising when the voxels are not square
%the wooraytrace.cpp was recompiled

% This code reads a grid of starting points from which needles are
% virtually inserted in the canopy. The file also defines the direction in
%which the needles are inserted. The canopy is defined by LAD arrays
% which are split into very thin vertical layers.

%the code now requires the  Statistics and Machine Learning Toolbox because
%the Nilson equation uses the randsample function

clear

fig=0;

%-----Inputs
%plot="EMS";
voxel_size=0.3;
voxel_size_z=voxel_size %0.01;
number_of_runs=1;
zen_ang=[30]; %[0 15 30 45 60]
%zen_ang=[15 30 45 60];
layer_thickness=1; %0 is original voxel size layer thickness to compute vertical profile, 1 is 1 m layer thickness
plot_list=["Hubbard"];
omega_v_present=1;
%plot_list=["EMS" "SERC" "Pasoh" "MorganMonroe"];
max_n=30; %this is the maximum number of hits to consider in the Nilson equation within a voxel
%-----End of inputs


%for plot_name=["EMS" "SERC" "Pasoh" "MorganMonroe"]
for plot_name=plot_list
    
    tic
    if plot_name=="Tonzi"
        bbox_x_min=0;
        bbox_y_min=0;
        bbox_z_min=0;
        bbox_x_max=100;
        bbox_y_max=100;
        bbox_z_max=12;
        path=strcat('H:\Tonzi relative variance\');
        %path=strcat('/Volumes/LaCie SSD/Tonzi relative variance/');
        
        file_lad=strcat(path,'landscape_v2.asc copy.txt');
        fid_lad=fopen(file_lad);
    elseif plot_name=="Random"
        bbox_x_min=0;
        bbox_y_min=0;
        bbox_z_min=0;
        bbox_x_max=100;
        bbox_y_max=100;
        bbox_z_max=12;
        path=strcat('H:\random array relative variance\');
        %path=strcat('/Volumes/LaCie SSD/Tonzi relative variance/');
        
        file_lad=strcat(path,'random_array_2.txt');
        fid_lad=fopen(file_lad);
    else
        path=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot_name,'/');
        %path=strcat('H:\ReadRXP\voxTrace_',plot_name,'\');
        %path=strcat('E:\',plot_name,'\');
        %get bbox
        
        file_bbox=strcat(path,plot_name,'_readRXP_batch.bat');
        fid_bbox=fopen(file_bbox);
        fscanf(fid_bbox, '%s', [12]);
        
        bbox_x_min=fscanf(fid_bbox, '%d', [1]);
        bbox_y_min=fscanf(fid_bbox, '%d', [1]);
        bbox_z_min=fscanf(fid_bbox, '%d', [1]);
        bbox_x_max=fscanf(fid_bbox, '%d', [1]);
        bbox_y_max=fscanf(fid_bbox, '%d', [1]);
        bbox_z_max=fscanf(fid_bbox, '%d', [1]);
        
        fclose(fid_bbox);
        if omega_v_present
            file_lad=strcat(path,plot_name,'_5mgrid_LAD_no_zeros.asc_occlusion_corr_thresh_30with_omega_v.asc');
        else
            file_lad=strcat(path,plot_name,'_5mgrid_LAD_no_zeros.asc_occlusion_corr_thresh_30.asc');
        end
        fid_lad=fopen(file_lad);
    end
    
    lenght_y=bbox_y_max-bbox_y_min; %this is the plot size in y direction in meters
    bbox_y_min=bbox_y_min-lenght_y;
    
    %for testing:
    % bbox_x_min=10;
    % bbox_x_max=30;
    % bbox_y_min=10;
    % bbox_y_max=30;
    % bbox_z_min=324;
    % bbox_z_max=360;
    
    
    if rem((bbox_x_max-bbox_x_min),voxel_size)~=0||rem((bbox_y_max-bbox_y_min),voxel_size)~=0||rem((bbox_z_max-bbox_z_min),voxel_size)~=0, message='warning: voxel array size does not allow good fit with voxel size. THIS WILL RESULT IN INCORRECT STATISTICS FOR VOXELS BOUNDING THE BBOX',end
    
    dimX=double(round((bbox_x_max-bbox_x_min)/voxel_size));
    dimY=double(round((bbox_y_max-bbox_y_min)/voxel_size));
    dimZ=double(round((bbox_z_max-bbox_z_min)/voxel_size_z));
    
    gridSize = [dimX dimY dimZ];
    gridBounds = [bbox_x_min bbox_y_min bbox_z_min bbox_x_max bbox_y_max bbox_z_max];
    bbox_min=gridBounds([1 2 3]);
    bbox_max=gridBounds([4 5 6]);
    
    %create an array to reference the voxel indices provided by wooRaytrace.cpp and
    %the coordinates of the voxel (Minimum vertex). NO I don't because the
    %following code from wootracetest.m provides them as per the following
    %lines:
    % for nVoxel = 1:size(indexes,1)
    %             [dy, dx, dz] = ind2sub([dimX dimY dimZ],indexes(nVoxel,1));
    %             vx = bbox_x_min + [(dx-1)*voxel_size (dx-1)*voxel_size + voxel_size];
    %             vy = bbox_y_min + [(dy-1)*voxel_size (dy-1)*voxel_size + voxel_size];
    %             vz = bbox_z_min + [(dz-1)*voxel_size (dz-1)*voxel_size + voxel_size];
    
    
    
    
    % file_res=strcat(path,'test_results.txt');
    % fid_res=fopen(file_res,'w');
    
    %event='Pulse';
    %count=0;
    
    %pathlenghts_hits=zeros(dimX*dimY*dimZ,1);
    
    %load LAD matrix
    
    
    
    dimZ_LAD=double(round((bbox_z_max-bbox_z_min)/voxel_size));
    
    LAD_G=zeros(dimX*dimY*dimZ_LAD,3);
    
    while ~feof(fid_lad)
        
        xmin=fscanf(fid_lad, '%f', [1]);
        if isempty(xmin), break, end
        ymin=fscanf(fid_lad, '%f', [1]);
        zmin=fscanf(fid_lad, '%f', [1]);
        xmax=fscanf(fid_lad, '%f', [1]);
        ymax=fscanf(fid_lad, '%f', [1]);
        zmax=fscanf(fid_lad, '%f', [1]);
        lad=fscanf(fid_lad, '%f', [1]);
        if omega_v_present
            omega=fscanf(fid_lad, '%f', [1]);
        else
            omega=1;
        end%***********new
        G_code=fscanf(fid_lad, '%f', [1]);
        
        
        dx=(xmin/voxel_size)+1;
        dy=(ymin/voxel_size)+1;
        dz=(zmin/voxel_size)+1;
        %         I want to calculate a unique index for each voxel
        %         used previously:
        %          dz=ceil(idx/(dimX*dimY));
        %             remain=idx-(dz-1)*dimX*dimY;
        %             dx=ceil(remain/dimY);
        %             dy=remain-(dx-1)*dimY;
        %             vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
        %             vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
        %             vz = bbox_z_min + [(dz-1)*voxel_size (dz)*voxel_size];
        %             I want to reverse this
        index=int32(dy+(dx-1)*dimY+(dz-1)*dimX*dimY);
        LAD_G(index,1)=lad;
        LAD_G(index,2)=omega;
        LAD_G(index,3)=G_code;
        
        %populating the duplicate array
        dy=((ymin+lenght_y)/voxel_size)+1;
        index=int32(dy+(dx-1)*dimY+(dz-1)*dimX*dimY);
        LAD_G(index,1)=lad;
        LAD_G(index,2)=omega;
        LAD_G(index,3)=G_code;
        
        
    end
    toc
    
    %for zenith_angle=[0 15 30 45 60]
    for zenith_angle=zen_ang
        tic
        pathlenghts=zeros(dimX*dimY*dimZ,1);
        
        if plot_name=="Tonzi" || plot_name=="Random"
            
            %file_grid=strcat('/Volumes/LaCie SSD/Tonzi relative variance/needle_grid_theta_',num2str(zenith_angle),'.txt');
            file_grid=strcat(path,'needle_grid_theta_',num2str(zenith_angle),'.txt');
        else
            %file_grid=strcat(path,'Relative Variance\needle_grid_theta_',num2str(zenith_angle),'_10cm_spacing.txt');
            file_grid=strcat(path,'Relative Variance/needle_grid_theta_',num2str(zenith_angle),'.txt');
        end
        
        fid_grid=fopen(file_grid);
        
        number_of_needles=fscanf(fid_grid, '%d', [1]);
        
        
        
        
        
        erecto=[0.001813, 0.014469, 0.038908, 0.072257, 0.11037, 0.148837, 0.182804, 0.208301, 0.222241];
        plano=[0.217993772,  0.205566574, 0.181747777, 0.149126381, 0.111758486, 0.07430429,  0.040992495, 0.016051798, 0.002458428];
        spheri=[0.013816774, 0.044002464, 0.072956086, 0.099621579, 0.1232069, 0.143184018, 0.158760891, 0.169409487, 0.175041802];
        plagio=[0.007084426, 0.053135816, 0.126252094, 0.19247555,  0.220395847, 0.197361602, 0.133842925, 0.059766886, 0.009684853];
        uniform=[0.111111111,  0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111];
        
        %Fliesvox LUT used 3 G functions:
        % 1 = spherical leaf angle distribution
        % 2 = planophile leaf angle distribution
        % 3 = erectrophile leaf angle distribution
        
        theta_l=deg2rad([5, 15, 25, 35, 45, 55, 65, 75, 85]);
        
        %read first line to get needle zenith angle, all needles are sent in the
        %same direction
        %         OrigX=fscanf(fid_grid, '%f', [1]);
        %         OrigY=fscanf(fid_grid, '%f', [1]);
        %         OrigZ=fscanf(fid_grid, '%f', [1]);
        dirX=fscanf(fid_grid, '%f', [1]);
        dirY=fscanf(fid_grid, '%f', [1]);
        dirZ=fscanf(fid_grid, '%f', [1]);
        %frewind(fid_grid);
        
        theta=asin(dirY)
        %if theta>pi/2,theta=pi-theta;end . this should not appen
        
        if plot_name=="Tonzi" || plot_name=="Random"
            
            %file_res=strcat('/Volumes/LaCie SSD/Tonzi relative variance/results_theta_',num2str(round(rad2deg(theta))),'.txt');
            file_res=strcat(path,'results_theta_',num2str(round(rad2deg(theta))),'.txt');
        else
            if omega_v_present
                %file_res=strcat(path,'Relative Variance/results_theta_',num2str(round(rad2deg(theta))),'_10cm_spacing.txt');
                if layer_thickness
                    file_res=strcat(path,'Relative Variance/results_theta_',num2str(round(rad2deg(theta))),'.txt');
                else
                    file_res=strcat(path,'Relative Variance/results_theta_',num2str(round(rad2deg(theta))),'_30cmlayers.txt');
                end
            else
                if layer_thickness
                    file_res=strcat(path,'Relative Variance/results_theta_',num2str(round(rad2deg(theta))),'no_omega.txt');
                else
                    file_res=strcat(path,'Relative Variance/results_theta_',num2str(round(rad2deg(theta))),'no_omega_30cmlayers.txt');
                end
            end
        end
        
        fid_res=fopen(file_res,'w');
        %There are 3 possible G values, compute them once and store in
        %G(1,2 or 3)
        G=zeros(4,1);
        
        for dist=1:4
            switch dist
                case 1
                    %distr is spherical
                    g=spheri;
                case 2
                    %distr is planophile
                    g=plano;
                case 3
                    %distr is erectophile
                    g=erecto;
                case 4
                    %distr is uniform
                    g=uniform;
            end
            for i=1:9
                x_=acos(cot(theta)*cot(theta_l(i)));
                if theta <= pi/2-theta_l(i), S=cos(theta)*cos(theta_l(i));
                else S=cos(theta)*cos(theta_l(i))*(1+(2*(tan(x_)-x_)/pi));
                end
                G(dist)=G(dist)+g(i)*S;
            end
        end
        
        
        
        
        
        for run=1:number_of_runs
            if layer_thickness
                vertical_hits=zeros(number_of_needles,bbox_z_max-bbox_z_min+1);
            else
                vertical_hits=zeros(number_of_needles,dimZ+1);
            end
            hits=zeros(number_of_needles,1);
            
            needle_number=1;
            while ~feof(fid_grid)
                %hits(needle_number)=0;
                
                
                OrigX=fscanf(fid_grid, '%f', [1]);
                if isempty(OrigX), break, end
                OrigY=fscanf(fid_grid, '%f', [1]);
                OrigZ=fscanf(fid_grid, '%f', [1]);
                %                 dirX=fscanf(fid_grid, '%f', [1]);
                %                 dirY=fscanf(fid_grid, '%f', [1]);
                %                 dirZ=fscanf(fid_grid, '%f', [1]);
                
                %     position = ftell(fid_rxp);
                %     event=fscanf(fid_rxp, '%d',[1]);
                %     if isempty(event), event=0;end
                
                
                %    if event
                %         x=fscanf(fid_rxp, '%f', [1]);
                %         y=fscanf(fid_rxp, '%f', [1]);
                %         z=fscanf(fid_rxp, '%f', [1]);
                %
                %         lineCoords = [OrigX OrigY OrigZ x y z]; %line is set by pulse origin and first hit
                %  else
                %     fseek(fid_rxp,position,-1);
                %there was no hit, need to find the intersection of pulse with bbox
                %I have 6 faces to check
                origin    = [OrigX, OrigY, OrigZ];
                direction = -[dirX, dirY, dirZ];
                vmin      = [bbox_x_min, bbox_y_min, bbox_z_min];  % vertex min
                vmax      = [ bbox_x_max,  bbox_y_max,  bbox_z_max];        % vertex max
                
                [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
                if ~flag,message='error on flag in rayboxinteresction for intersection bbox',end
                intersection_bbox = origin + tmin*direction;
                %these next lines are to put the intersection exactly on the bbox
                %limit, as a problem with the call to wooRaytrace was observed when
                %a point on linecoords is outside (even slightly, but how slightly
                %is not clear) the gridbounds
                intersection_bbox(intersection_bbox>bbox_max)=bbox_max(intersection_bbox>bbox_max);
                intersection_bbox(intersection_bbox<bbox_min)=bbox_min(intersection_bbox<bbox_min);
                
                %         figure;
                %         hold on;
                %         grid on;
                %
                %         % box (voxel)
                %         vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
                %         faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
                %         h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
                %         set(h,'FaceAlpha',0.5);
                %
                %         % origin
                %         text(origin(1), origin(2), origin(3), 'origin');
                %         plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);
                %
                %         % direction
                %         quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);
                %
                %         % intersection
                %         plot3(intersection_bbox(1), intersection_bbox(2), intersection_bbox(3), 'r.', 'MarkerSize', 15);
                %
                %         view(60,30);
                %         axis tight;
                %         xlabel('x');
                %         ylabel('y');
                %         zlabel('z');
                
                lineCoords = [OrigX OrigY OrigZ intersection_bbox]; %line is set by pulse origin and intersection with bbox
                
                % end
                
                
                
                
                if fig
                    gridImage = zeros(dimX,dimY,dimZ);
                    figure;
                    
                    ha3 = subplot(1,1,1);
                    axes(ha3);
                    hold all;
                    xlabel('X');
                    ylabel('Y');
                    zlabel('Z');
                    view(45,45);
                    box on;
                    grid on;
                    axis equal;
                end
                
                
                
                %fprintf('Line coordinates:%s\n',num2str(lineCoords));
                
                
                indexes = wooRaytrace(gridSize,gridBounds,lineCoords);
                
                if fig
                    axes(ha3);
                    for nVoxel = 1:size(indexes,1)
                        [dy, dx, dz] = ind2sub([dimX dimY dimZ],indexes(nVoxel,1));
                        vx = bbox_x_min + [(dx-1)*voxel_size (dx-1)*voxel_size + voxel_size];
                        vy = bbox_y_min + [(dy-1)*voxel_size (dy-1)*voxel_size + voxel_size];
                        vz = bbox_z_min + [(dz-1)*voxel_size_z (dz-1)*voxel_size_z + voxel_size_z];
                        fv.vertices = [[vx(1) vy(1) vz(1)];[vx(2) vy(1) vz(1)];[vx(2) vy(2) vz(1)];[vx(1) vy(2) vz(1)]; ...
                            [vx(1) vy(1) vz(2)];[vx(2) vy(1) vz(2)];[vx(2) vy(2) vz(2)];[vx(1) vy(2) vz(2)]];
                        fv.faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
                        h = patch(fv, 'FaceColor', 'blue', 'EdgeColor', 'white');
                        set(ha3,'XLim',[bbox_x_min bbox_x_max],'YLim',[bbox_y_min bbox_y_max],'ZLim',[bbox_z_min bbox_z_max],'XTick',bbox_x_min: bbox_x_max,'YTick',bbox_y_min :bbox_y_max,'ZTick',bbox_z_min :bbox_z_max);
                        h.FaceAlpha = 0.2;
                    end
                    plot3([lineCoords(1) lineCoords(4)]',[lineCoords(2) lineCoords(5)]',[lineCoords(3) lineCoords(6)]','-r*');
                end
                
                %for the traversed voxels containing the hit, compute omega I statistics
                %if event
                
                %may 2020: Here I am only looking at the voxel whithin which the hot occured,
                %that is why the code in this section is not vectorized
                
                %         [dy, dx, dz] = ind2sub([dimX dimY dimZ],indexes(size(indexes,1),1));
                %         vx = bbox_x_min + [(dx-1)*voxel_size (dx-1)*voxel_size + voxel_size];
                %         vy = bbox_y_min + [(dy-1)*voxel_size (dy-1)*voxel_size + voxel_size];
                %         vz = bbox_z_min + [(dz-1)*voxel_size (dz-1)*voxel_size + voxel_size];
                
                %second version:
                %         idx=double(indexes(size(indexes,1),1));
                %         dz=ceil(idx/(dimX*dimY));
                %         remain=idx-(dz-1)*dimX*dimY;
                %         dx=ceil(remain/dimY);
                %         dy=remain-(dx-1)*dimY;
                %         vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
                %         vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
                %         vz = bbox_z_min + [(dz-1)*voxel_size (dz)*voxel_size];
                %
                %         origin= [x y z];
                %         direction = [dirX, dirY, dirZ];
                %         vmin      = [vx(1), vy(1), vz(1)];  % vertex min
                %         vmax      = [ vx(2), vy(2), vz(2)];        % vertex max
                %
                %         [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
                %         intersection1 = origin + tmin*direction;
                %
                %         if fig
                %             vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
                %             faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
                %             h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
                %             set(h,'FaceAlpha',0.5);
                %
                %             % origin
                %             text(origin(1), origin(2), origin(3), 'origin');
                %             plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);
                %
                %             % direction
                %             quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);
                %
                %             % intersection
                %             plot3(intersection1(1), intersection1(2), intersection1(3), 'r.', 'MarkerSize', 15);
                %         end
                %
                %         pathlenghts_hits(indexes(size(indexes,1)))=pathlenghts_hits(indexes(size(indexes,1)))+sqrt((origin(1)-intersection1(1))^2+(origin(2)-intersection1(2))^2+(origin(3)-intersection1(3))^2);
                %
                %         %count=count+1;
                %         %fprintf(fid_res,'pathlenght hit: %f\n', pathlenght_hit);
                %     end
                
                
                %go through traversed voxels to compute GAP omega statistics, I will use the
                %rayBoxIntersection code
                %for nVoxel = 1:size(indexes,1)-event %is this is a pulse that ended with a hit, I need to exclude the last voxel indice
                %TODO: vertorize this for loop:
                
                
                
                %         [dy, dx, dz] = ind2sub([dimX dimY dimZ],indexes(nVoxel,1));
                %         vx = bbox_x_min + [(dx-1)*voxel_size (dx-1)*voxel_size + voxel_size];
                %         vy = bbox_y_min + [(dy-1)*voxel_size (dy-1)*voxel_size + voxel_size];
                %         vz = bbox_z_min + [(dz-1)*voxel_size (dz-1)*voxel_size + voxel_size];
                %TODO: replace Hit0 and Pulse with 0 and 1 to avoid strcmp call and
                %speed up code (modif in cpp code)..DONE
                
                %second version (avoiding the call to in2sub:
                %         idx=double(indexes(nVoxel,1));
                %         dz=ceil(idx/(dimX*dimY));
                %         remain=idx-(dz-1)*dimX*dimY;
                %         dx=ceil(remain/dimY);
                %         dy=remain-(dx-1)*dimY;
                
                %vectorized:
                %idx=double(indexes(nVoxel,1));
                dz=ceil(double(indexes)/(dimX*dimY));
                remain=double(indexes)-(dz-1)*dimX*dimY;
                dx=ceil(remain/dimY);
                dy=remain-(dx-1)*dimY;
                
                dz_larger=ceil(dz/(voxel_size/voxel_size_z)); %this is the index of the voxel in the original
                %matrix format where voxels are cubes, it will be used to exctract
                %the LAD and G function for the voxel being processed
                
                vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
                vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
                vz = bbox_z_min + [(dz-1)*voxel_size_z (dz)*voxel_size_z];
                
                origin    = [OrigX, OrigY, OrigZ];
                direction = [dirX, dirY, dirZ];
                vmin      = [vx(:,1), vy(:,1), vz(:,1)];  % vertex min
                vmax      = [vx(:,2), vy(:,2), vz(:,2)];        % vertex max
                
                
                
                [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
                if ~flag,message='error on flag in rayboxinteresction for intersection 1',end
                intersection1 = origin + tmin*direction;
                
                if fig
                    %figure;
                    hold on;
                    grid on;
                    
                    % box (voxel)
                    %             vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
                    %             faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
                    %             h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
                    %             set(h,'FaceAlpha',0.5);
                    %
                    %             % origin
                    %             text(origin(1), origin(2), origin(3), 'origin');
                    %             plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);
                    %
                    %             % direction
                    %             quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);
                    %
                    % intersection
                    for i=1:size(intersection1,1)
                        plot3(intersection1(i,1), intersection1(i,2), intersection1(i,3), 'r.', 'MarkerSize', 15);
                    end
                    
                    view(60,30);
                    %axis tight;
                    xlabel('x');
                    ylabel('y');
                    zlabel('z');
                end
                
                %looking for exit point,
                %if strcmp(event,'Hit0'),origin= [x y z], else origin = [intersection_bbox];end
                direction = -[dirX, dirY, dirZ];
                % box (voxel)
                [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
                if ~flag,message='error on flag in rayboxinteresction for intersection 2',end
                intersection2 = origin + tmin*direction;
                
                if fig
                    %             vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
                    %             faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
                    %             h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
                    %             set(h,'FaceAlpha',0.5);
                    %
                    %             % origin
                    %             text(origin(1), origin(2), origin(3), 'origin');
                    %             plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);
                    %
                    %             % direction
                    %             quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);
                    %
                    % intersection
                    for i=1:size(intersection2,1)
                        plot3(intersection2(i,1), intersection2(i,2), intersection2(i,3), 'g.', 'MarkerSize', 15);
                    end
                    %plot3(intersection2(1), intersection2(2), intersection2(3), 'r.', 'MarkerSize', 15);
                end
                
                pathlenghts(indexes(:,1))=sqrt((intersection2(:,1)-intersection1(:,1)).^2+(intersection2(:,2)-intersection1(:,2)).^2+(intersection2(:,3)-intersection1(:,3)).^2);
                
                %determining if there is a contact with leaf by needle
                
                for i=1:size(indexes,1)
                    
                    index=int32(dy(i)+(dx(i)-1)*dimY+(dz_larger(i)-1)*dimX*dimY);
                    %G=0.85; %for now, need to compute it based on zentih angle of needle
                    
                    lad=LAD_G(index,1);
                    if lad==0, continue,end
                    path_length=pathlenghts(indexes(i));
                    omega=LAD_G(index,2);;
                    
                    %x=rand; %this is a random number between 0 and 1. I will use it to generate a 0 or 1 (contact)
                    %based on the probability of a contact, if x is lower that the
                    %probability, then there is a contact
                    
                    %P1=1-exp(-lad*path_length*G(LAD_G(index,2)));
                    
                    for n=1:max_n
                        summation=0;
                        for ii=0:n-1
                            %NOTE: the omega were removed from these
                            %equations (only one omega in case of
                            %summation) because the lad value is not the
                            %actual lad but the effective lad, it already
                            %is equal to lad*omega.
                            %effective_lad=lad*omega
                            iteration=(factorial(n-1)/(factorial(n-1-ii)*factorial(ii)))*((omega*lad*path_length*G(LAD_G(index,3)))^(ii+1)/factorial(ii+1))*(1-omega)^(n-ii-1);
                            %iteration=(factorial(n-1)/(factorial(n-1-ii)*factorial(ii)))*((omega^2*lad*path_length*G(LAD_G(index,3)))^(ii+1)/factorial(ii+1))*(1-omega)^(n-ii-1);
                            summation=summation+iteration;
                        end
                        P(n+1)=exp(-lad*path_length*G(LAD_G(index,3)))*summation;
                        %P(n+1)=exp(-omega*lad*path_length*G(LAD_G(index,3)))*summation; %Probabilities of n hits are stores with n+1 index because index 0 is not allowed
                    end
                    P(1)=exp(-lad*path_length*G(LAD_G(index,3)));
                    
                    num_hits = randsample(max_n+1,1,true,P);
                    
                    hits(needle_number)=hits(needle_number)+(num_hits-1); %there is a hit
                    if layer_thickness
                        height=int16((vz(i,1)+vz(i,2))/2)-bbox_z_min+1; %original line for 1m thick layers
                    else
                        height=int16((vz(i,2)-bbox_z_min+1)/voxel_size_z); %modified line for original voxel size layer thickness
                    end
                    vertical_hits(needle_number,height)=vertical_hits(needle_number,height)+(num_hits-1);

                end
                needle_number=needle_number+1;
                
                %count=count+1;
                %fprintf(fid_res,'pathlenght gap: %f\n', pathlenght);
                %end
                
            end
            frewind(fid_grid);
            bidon=fscanf(fid_grid, '%d', [1]);
            bidon=fscanf(fid_grid, '%f', [3]);
            variance(run)=var(hits);
            average(run)=mean(hits);
            relative_variance(run)=var(hits)/mean(hits);
            fprintf(fid_res, 'variance: %.3f  Mean: %.3f  Relative variance: %.3f\n',var(hits), mean(hits), var(hits)/mean(hits));
        end
        %image the last run
        needle=1;
        horizontal_view=zeros(int32(sqrt(size(hits,1))),int32(sqrt(size(hits,1))));
        for i=1:int32(sqrt(size(hits,1)))
            for j=1:int32(sqrt(size(hits,1)))
                horizontal_view(i,j)=hits(needle);
                %if horizontal_view(i,j)>20,horizontal_view(i,j)=1;end %to get better contrast in image
                needle=needle+1;
            end
        end
        figure
        image(horizontal_view);
        figure
        %averaged_vertical_hits=mean(vertical_hits,3);
        plot(sum(vertical_hits));
        average_Mean=mean(average)
        fprintf(fid_res, 'average variance: %.3f  average Mean: %.3f  average Relative variance: %.3f\n',mean(variance), mean(average), mean(relative_variance));
        fprintf(fid_res, 'omega: %.3f\n\n',2/(mean(relative_variance)+1));
        fprintf(fid_res, 'height - omega - variance - mean - relative variance\n');
        for h=1:size(vertical_hits,2)
            rel_var_of_layer=var(vertical_hits(:,h))/mean(vertical_hits(:,h));
            fprintf(fid_res, '%d %.3f %.3f %.3f %.3f\n',h,2/(rel_var_of_layer+1),var(vertical_hits(:,h)),mean(vertical_hits(:,h)),rel_var_of_layer);
        end
        
        fclose all;
        toc
    end
end