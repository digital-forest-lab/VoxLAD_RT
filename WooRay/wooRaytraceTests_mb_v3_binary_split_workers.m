% Compile the mex function
%run: mex wooRaytrace.cpp

%****This code only works for square floor matrices, there appears to be a
%problem with wooRatrace indexing the voxels when floor is not square
clear
%***********Inputs
fig=0;

voxel_size=0.3;
plot='Hubbard'

path=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where bin files are
%path=strcat('E:\',plot,'\'); %this is where bin files are
%path_res=strcat('F:\ReadRXP\voxTrace_',plot,'\'); %this is where statistic files are written
path_res=path; %strcat('E:\',plot,'\'); %this is where statistic files are written
%path=strcat('E:\test_RXP\');

worker=8;

%***************

switch worker
    case 1
        start_pos=1;
        end_pos=12; %31
    case 2
        start_pos=13; %32
        end_pos=24; %62
    case 3
        start_pos=25; %63
        end_pos=37; %changed from 93 
    case 4
        start_pos=38; %94
        end_pos=50; %124
    case 5
        start_pos=51; %125
        end_pos=63;%155
    case 6
        start_pos=64;%156
        end_pos=76;%186
    case 7
        start_pos=77;%187
        end_pos=89;%217
    case 8
        start_pos=90;%218
        end_pos=98;%242
    case 999 %does not split between workers
        start_pos=1;
        end_pos=98;
end
% start_pos=217; %this can be split between workers
% end_pos=242;


for pos=start_pos:end_pos
    pos
    if pos==206 && plot=="SERC", continue,end
    
    tic
    
    file_name=strcat('ScanPos',num2str(pos),'lewos.bin');
    file_rxp=strcat(path,file_name);
    file_rxp_res=strcat(path_res,file_name);
    %file_rxp_ascii=strcat(path,'test_ascii');
    
    fid_rxp=fopen(file_rxp);
    %fid_rxp_ascii=fopen(file_rxp_ascii);
    
    file_res_gap=strcat(file_rxp_res,'_statistics_gaps.txt');
    fidres_gap=fopen(file_res_gap,'w');
    file_res_hit=strcat(file_rxp_res,'_statistics_hits.txt');
    fidres_hit=fopen(file_res_hit,'w');
    
    bbox_x_min=fread(fid_rxp,1,'double');
    bbox_y_min=fread(fid_rxp,1,'double');
    bbox_z_min=fread(fid_rxp,1,'double');
    bbox_x_max=fread(fid_rxp,1,'double');
    bbox_y_max=fread(fid_rxp,1,'double');
    bbox_z_max=fread(fid_rxp,1,'double');
    
    % bbox_x_min=-30;
    % bbox_x_max=30;
    % bbox_y_min=-30;
    % bbox_y_max=30;
    % bbox_z_min=50;
    % bbox_z_max=80;
    
    if rem((bbox_x_max-bbox_x_min),voxel_size)~=0||rem((bbox_y_max-bbox_y_min),voxel_size)~=0||rem((bbox_z_max-bbox_z_min),voxel_size)~=0, message='warning: voxel array size does not allow good fit with voxel size. THIS WILL RESULT IN INCORRECT STATISTICS FOR VOXELS BOUNDING THE BBOX',end
    
    
    dimX=double(round((bbox_x_max-bbox_x_min)/voxel_size));
    dimY=double(round((bbox_y_max-bbox_y_min)/voxel_size));
    dimZ=double(round((bbox_z_max-bbox_z_min)/voxel_size));
    
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
    
    
    
    %event='Pulse';
    %count=0;
    pathlenghts_gaps=zeros(dimX*dimY*dimZ,1);
    pathlenghts_hits=zeros(dimX*dimY*dimZ,1);
    Ni_leaf=zeros(dimX*dimY*dimZ,1);
    
    while ~feof(fid_rxp)
        
        %event_a=fscanf(fid_rxp_ascii, '%d',[1]);
        event=fread(fid_rxp,1,'int');
        if isempty(event), break, end
        
        %OrigX=fscanf(fid_rxp, '%f', [1]);
        OrigX=fread(fid_rxp,1,'double');
        OrigY=fread(fid_rxp,1,'double');
        OrigZ=fread(fid_rxp,1,'double');
        dirX=fread(fid_rxp,1,'double');
        dirY=fread(fid_rxp,1,'double');
        dirZ=fread(fid_rxp,1,'double');
        
        %     OrigX_a=fscanf(fid_rxp_ascii, '%f', [1]);
        %     OrigY_a=fscanf(fid_rxp_ascii, '%f', [1]);
        %     OrigZ_a=fscanf(fid_rxp_ascii, '%f', [1]);
        %     dirX_a=fscanf(fid_rxp_ascii, '%f', [1]);
        %     dirY_a=fscanf(fid_rxp_ascii, '%f', [1]);
        %     dirZ_a=fscanf(fid_rxp_ascii, '%f', [1]);
        %     if round(OrigX*100000)~=round(OrigX_a*100000)
        %         message='pause herE';
        %     end
        %count=count+1;
        position = ftell(fid_rxp);
        %     position_a = ftell(fid_rxp_ascii);
        %event=fscanf(fid_rxp, '%d',[1]);
        event=fread(fid_rxp,1,'int');
        %     event_a=fscanf(fid_rxp_ascii, '%d',[1]);
        if isempty(event), event=0;end
        
        
        if event
            x=fread(fid_rxp,1,'double');
            y=fread(fid_rxp,1,'double');
            z=fread(fid_rxp,1,'double');
            %         x_a=fscanf(fid_rxp_ascii, '%f', [1]);
            %         y_a=fscanf(fid_rxp_ascii, '%f', [1]);
            %         z_a=fscanf(fid_rxp_ascii, '%f', [1]);
            %count=count+1;
            lineCoords = [OrigX OrigY OrigZ x y z]; %line is set by pulse origin and first hit
        else
            fseek(fid_rxp,position,-1);
            %         fseek(fid_rxp_ascii,position_a,-1);
            %there was no hit, need to find the intersection of pulse with bbox
            %I have 6 faces to check
            origin    = [OrigX, OrigY, OrigZ];
            direction = -[dirX, dirY, dirZ];
            vmin      = [bbox_x_min, bbox_y_min, bbox_z_min];  % vertex min
            vmax      = [ bbox_x_max,  bbox_y_max,  bbox_z_max];        % vertex max
            
            [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
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
            
        end
        
        
        
        
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
        if isempty(indexes),continue,end %the tls is outside the bbox and the pulse did not enter the bbox, or the hit is outside the bbox, which should never occur because of the bounds set in readRXP
        
        if fig
            axes(ha3);
            for nVoxel = 1:size(indexes,1)
                [dy, dx, dz] = ind2sub([dimX dimY dimZ],indexes(nVoxel,1));
                vx = bbox_x_min + [(dx-1)*voxel_size (dx-1)*voxel_size + voxel_size];
                vy = bbox_y_min + [(dy-1)*voxel_size (dy-1)*voxel_size + voxel_size];
                vz = bbox_z_min + [(dz-1)*voxel_size (dz-1)*voxel_size + voxel_size];
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
        if event==1 %hit is a leaf
            %         [dy, dx, dz] = ind2sub([dimX dimY dimZ],indexes(size(indexes,1),1));
            %         vx = bbox_x_min + [(dx-1)*voxel_size (dx-1)*voxel_size + voxel_size];
            %         vy = bbox_y_min + [(dy-1)*voxel_size (dy-1)*voxel_size + voxel_size];
            %         vz = bbox_z_min + [(dz-1)*voxel_size (dz-1)*voxel_size + voxel_size];
            
            %second version:
            idx=double(indexes(size(indexes,1),1));
            dz=ceil(idx/(dimX*dimY));
            remain=idx-(dz-1)*dimX*dimY;
            dx=ceil(remain/dimY);
            dy=remain-(dx-1)*dimY;
            vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
            vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
            vz = bbox_z_min + [(dz-1)*voxel_size (dz)*voxel_size];
            
            origin= [x y z];
            direction = [dirX, dirY, dirZ];
            vmin      = [vx(1), vy(1), vz(1)];  % vertex min
            vmax      = [ vx(2), vy(2), vz(2)];        % vertex max
            
            [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
            
            intersection1 = origin + tmin*direction;
            
            if fig
                vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
                faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
                h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
                set(h,'FaceAlpha',0.5);
                
                % origin
                text(origin(1), origin(2), origin(3), 'origin');
                plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);
                
                % direction
                quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);
                
                % intersection
                plot3(intersection1(1), intersection1(2), intersection1(3), 'r.', 'MarkerSize', 15);
            end
            if sqrt((origin(1)-intersection1(1))^2+(origin(2)-intersection1(2))^2+(origin(3)-intersection1(3))^2) > voxel_size * sqrt(3)
                message='error hits pathlenght too long, check bbox to be integer factor of voxel size'
                pathlong=sqrt((origin(1)-intersection1(1))^2+(origin(2)-intersection1(2))^2+(origin(3)-intersection1(3))^2)
                origin
            end
            pathlenghts_hits(indexes(size(indexes,1)))=pathlenghts_hits(indexes(size(indexes,1)))+sqrt((origin(1)-intersection1(1))^2+(origin(2)-intersection1(2))^2+(origin(3)-intersection1(3))^2);
            Ni_leaf(indexes(size(indexes,1)))=Ni_leaf(indexes(size(indexes,1)))+1;
            %count=count+1;
            %fprintf(fid_res,'pathlenght hit: %f\n', pathlenght_hit);
        end
        
        
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
        %PROBLEM: I am also including a pathelenght_gap on the voxel containing the hit, I want to
        %remove the last indexes value: indexes(size(indexes,1),1). SOLVED on next line:
        if event~=0,indexes(size(indexes,1):end)=[];end
        if isempty(indexes),continue,end
        dz=ceil(double(indexes)/(dimX*dimY));
        remain=double(indexes)-(dz-1)*dimX*dimY;
        dx=ceil(remain/dimY);
        dy=remain-(dx-1)*dimY;
        
        vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
        vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
        vz = bbox_z_min + [(dz-1)*voxel_size (dz)*voxel_size];
        
        origin    = [OrigX, OrigY, OrigZ];
        direction = [dirX, dirY, dirZ];
        vmin      = [vx(:,1), vy(:,1), vz(:,1)];  % vertex min
        vmax      = [vx(:,2), vy(:,2), vz(:,2)];        % vertex max
        
        
        
        [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
        intersection1 = origin + tmin*direction;
        
        if fig
            %figure;
            hold on;
            grid on;
            
            % box (voxel)
            vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
            faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
            h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
            set(h,'FaceAlpha',0.5);
            
            % origin
            text(origin(1), origin(2), origin(3), 'origin');
            plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);
            
            % direction
            quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);
            
            % intersection
            plot3(intersection1(1), intersection1(2), intersection1(3), 'r.', 'MarkerSize', 15);
            
            view(60,30);
            axis tight;
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end
        
        %looking for exit point,
        %if strcmp(event,'Hit0'),origin= [x y z], else origin = [intersection_bbox];end
        direction = -[dirX, dirY, dirZ];
        % box (voxel)
        [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax);
        intersection2 = origin + tmin*direction;
        
        if fig
            vertices = [vmax(1) vmin(2) vmin(3); vmax(1) vmax(2) vmin(3); vmin(1) vmax(2) vmin(3); vmin(1) vmax(2) vmax(3); vmin(1) vmin(2) vmax(3); vmax(1) vmin(2) vmax(3); vmin; vmax ];
            faces = [1 2 3 7; 1 2 8 6; 1 6 5 7; 7 5 4 3; 2 8 4 3; 8 6 5 4];
            h= patch('Vertices',vertices,'Faces',faces,'FaceColor','green');
            set(h,'FaceAlpha',0.5);
            
            % origin
            text(origin(1), origin(2), origin(3), 'origin');
            plot3(origin(1), origin(2), origin(3), 'k.', 'MarkerSize', 10);
            
            % direction
            quiver3(origin(1), origin(2), origin(3), direction(1), direction(2), direction(3), 15);
            
            % intersection
            plot3(intersection2(1), intersection2(2), intersection2(3), 'r.', 'MarkerSize', 15);
        end
        if sqrt((intersection2(:,1)-intersection1(:,1)).^2+(intersection2(:,2)-intersection1(:,2)).^2+(intersection2(:,3)-intersection1(:,3)).^2) > voxel_size * sqrt(3)
            message='error gaps pathlenght too long, check bbox to be integer factor of voxel size'
            pathlong=sqrt((intersection2(:,1)-intersection1(:,1)).^2+(intersection2(:,2)-intersection1(:,2)).^2+(intersection2(:,3)-intersection1(:,3)).^2)
            origin
        end
        pathlenghts_gaps(indexes(:,1))=pathlenghts_gaps(indexes(:,1))+sqrt((intersection2(:,1)-intersection1(:,1)).^2+(intersection2(:,2)-intersection1(:,2)).^2+(intersection2(:,3)-intersection1(:,3)).^2);
        
        %count=count+1;
        %fprintf(fid_res,'pathlenght gap: %f\n', pathlenght);
        %end
        
    end
    
    %Write results to files (one for gap pathlenghts, one for hits pathlenghts
    %and Number of hits on leaves
    
    fprintf(fidres_gap,'%.2f %.2f %.2f %.2f %.2f %.2f\n',bbox_x_min,bbox_y_min,bbox_z_min,bbox_x_max,bbox_y_max,bbox_z_max);
    fprintf(fidres_gap,'%.2f\n',voxel_size);
    
    index=find(pathlenghts_gaps);
    for i=1:size(index,1)
        
        dz=ceil(double(index(i))/(dimX*dimY));
        remain=double(index(i))-(dz-1)*dimX*dimY;
        dx=ceil(remain/dimY);
        dy=remain-(dx-1)*dimY;
        
        vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
        vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
        vz = bbox_z_min + [(dz-1)*voxel_size (dz)*voxel_size];
        
        fprintf(fidres_gap,'%.2f %.2f %.2f %.2f %.2f %.2f %.3f\n',vx(1),vy(1),vz(1),vx(2),vy(2),vz(2),pathlenghts_gaps(index(i)));
        
    end
    
    index=find(pathlenghts_hits);
    for i=1:size(index,1)
        dz=ceil(double(index(i))/(dimX*dimY));
        remain=double(index(i))-(dz-1)*dimX*dimY;
        dx=ceil(remain/dimY);
        dy=remain-(dx-1)*dimY;
        
        vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
        vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
        vz = bbox_z_min + [(dz-1)*voxel_size (dz)*voxel_size];
        
        fprintf(fidres_hit,'%.2f %.2f %.2f %.2f %.2f %.2f %.3f %d\n',vx(1),vy(1),vz(1),vx(2),vy(2),vz(2),pathlenghts_hits(index(i)),Ni_leaf(index(i)));
        
    end
    toc
end
fclose all;
