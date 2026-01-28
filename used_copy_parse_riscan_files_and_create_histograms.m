%this code reads riscan files, crops the plot area ignores the deviation
%above 25 points,  targets the first and single hits and corrects the
%reflectance for those points, and creates 2 histograms, one for singles
%only, one for singles and first. It can also create the lvox input files
clear

%_____inputs
plot_name='Bartlett';


path=strcat('C:\Temporary\for histograms\'); %this is where the batch, in and lvox input files are writen
path_rpy=strcat('F:\ReadRXP\voxTrace_Bartlett\'); %this is where RPY and ascii files are

fileres_first=strcat(path,'histogram_',plot_name,'_first.csv');
fileres_single_and_first=strcat(path,'histogram_',plot_name,'_single_and_first.csv');

start_pos=1;
end_pos=98;
%positions=[1 49 97 145 193 241] %diagonal
positions=[1 17 37 49 61 81 85] %diagonal

%number of scans done per day(used to see effect of rain events on return
%intensity)
num_scans_day1=114;
num_scans_day2=num_scans_day1+128;
num_scans_day3=num_scans_day2+999;

flag=1; %flag 1= generate histograms (does not write to files, only generate histogram matrix), 2:
%prepare parsed files (uses a refl threshold to separate l\w)

only_create_batch=0; %this is set to 1 to only create the .in and batch files for runing lvox (in case the flag-2 run was interupted), otherwise set to 0

minimum_dist=3.5; %this is the minimum distance from scanner to be included in the leaf wood refl threshold, this is because closer points are not accurate for intensity and are numerous
maximum_dist=100; %above this distance a point is ignored for histogram
distance_above_bbox_z_to_ignore=4; %EMS at 5, Morgan 11, Pasoh 9, serc 4

%: these are not known or used if flag=1
leaf_noise_threshold=8;
leaf_wood_threshold=36;

res=0.3; %resolution at which to run lvox (voxel size)

res_angu=0.04; %angle between consecutive laser shots

if strcmp(plot_name,'Bartlett') %updated April 2024
bbox_x_min=-47;
bbox_y_min=-8;
bbox_z_min=250;
bbox_x_max=3;
bbox_y_max=42;
bbox_z_max=280;
end


if strcmp(plot_name,'MorganMonroe') %updated oct 2020
bbox_x_min=-51;
bbox_y_min=-6;
bbox_z_min=239;
bbox_x_max=9;
bbox_y_max=54;
bbox_z_max=290;
end

if strcmp(plot_name,'SERC') %updated oct 2020
 
bbox_x_min=-53;
bbox_y_min=-2;
bbox_z_min=-15;
bbox_x_max=7;
bbox_y_max=58;
bbox_z_max=48;
end

if strcmp(plot_name,'EMS') %updated oct 2020
   
bbox_x_min=-6;
bbox_y_min=-59;
bbox_z_min=324;
bbox_x_max=54;
bbox_y_max=1;
bbox_z_max=360;
end

if strcmp(plot_name,'Pasoh') %updated oct 2020
   
bbox_x_min=-62;
bbox_y_min=-1;
bbox_z_min=127;
bbox_x_max=-2;
bbox_y_max=59;
bbox_z_max=193;
end

%_____end of inputs

xmin=bbox_x_min-10; %oct 2018: this does not define the bbox, it appears to be just a limit beyond which points are not written to file, not sure why I don't just use bbox x-1, for example, and why point beyond the bbox are included at all, perhaps it is to include those points in the point clouds in case I want to run lvox with larger bbox latter on. Also, a certain buffer is needed because some scans may be slightly outside the bbox, and those points close to it may not be counted as Nbefore
xmax=bbox_x_max+10;
ymin=bbox_y_min-10;
ymax=bbox_y_max+10;

if flag==2
    file=strcat(path,plot_name,'_batch.bat');
    fidres_batch=fopen(file,'w');
    
    
    % different angles for vert and hor positions are set
    %mod(pos,2) gives 0 for even (hor scans) and 1 for odd (vetical scans)
    %so I should use index+1, 1 will refer to hor scans, 2 to vert scans
%     min_hor_angle(1)=325;
%     max_hor_angle(1)=70; %was 395
%     min_vert_angle(1)=55;
%     max_vert_angle(1)=70; %was 125
%      min_hor_angle(2)=0;
%     max_hor_angle(2)=360;
%     min_vert_angle(2)=30;
%     max_vert_angle(2)=100; %this is not really max angle, but angle of opening, so if the angle entered here is 100, and the min vert angle is 30, then the max vert angle is 130, this is because JFT changed the angle input format in lvox..
   %March 16 2020 modifications to FOV:
    min_hor_angle(1)=331;
    max_hor_angle(1)=58; %was 395
    min_vert_angle(1)=61;
    max_vert_angle(1)=58; %was 125
     min_hor_angle(2)=0;
    max_hor_angle(2)=360;
    min_vert_angle(2)=32;
    max_vert_angle(2)=98; %this is not really max angle, but angle of opening, so if the angle entered here is 100, and the min vert angle is 30, then the max vert angle is 130, this is because JFT changed the angle input format in lvox..
   
    
    
    %___ end of input section (oct 2018: this is an old end mark)
    
    file=strcat(path,plot_name,'.in');
    fidres_all=fopen(file,'w');
    fprintf(fidres_all,'%d %d %d %d %d %d\n', bbox_x_min, bbox_y_min, bbox_z_min, bbox_x_max, bbox_y_max, bbox_z_max);
end

if flag==1, refl_hist=zeros(101,4,2);end


for pos=positions 
%for pos=start_pos:end_pos
    %         for pos=51:51
    pos
    
    if pos<10
        file_rpy=strcat(path_rpy,'ScanPos00',num2str(pos),'.RPY');
        %file= strcat(path,'ScanPos00',num2str(pos),'.ascii');
        file= strcat(path_rpy,'ScanPos00',num2str(pos),'.ascii');
        
    elseif pos<100
        file_rpy=strcat(path_rpy,'ScanPos0',num2str(pos),'.RPY');
        %file= strcat(path,'ScanPos0',num2str(pos),'.ascii');
        file= strcat(path_rpy,'ScanPos0',num2str(pos),'.ascii');
    else
        file_rpy=strcat(path_rpy,'ScanPos',num2str(pos),'.RPY');
        %file= strcat(path,'ScanPos',num2str(pos),'.ascii');
        file= strcat(path_rpy,'ScanPos',num2str(pos),'.ascii');
        %
    end
    
    if flag==2 && only_create_batch==0
        %         if pos<10
        %             fileres= strcat(path,'ScanPos00',num2str(pos),'_lvox_input.txt');
        %         elseif pos<100
        %             fileres= strcat(path,'ScanPos0',num2str(pos),'_lvox_input.txt');
        %         else
        %             fileres= strcat(path,'ScanPos',num2str(pos),'_lvox_input.txt');
        %             %
        %         end
        fileres= strcat(path,plot_name,'_pos',num2str(pos),'_lvox_input.txt');
        fidres=fopen(fileres,'w');
    end
    fid_rpy=fopen(file_rpy);
    fid=fopen(file);
   
    
    for j=1:2,fgets(fid_rpy);end  %skip two first lines
    bidon=fscanf(fid_rpy,'%c',[4]); %skip 4 first characters of line
    if rem(pos, 2) ~= 0
        yaw=fscanf(fid_rpy, '%f', [1]); %these lines are to assign the yaw to both odd and even scans, as per March 2020 results verification
    else
        bidon=fscanf(fid_rpy, '%f', [1]);
    end
    fgets(fid_rpy);
    bidon=fscanf(fid_rpy,'%c',[2]); %skip 2 first characters of line
    x_tls=fscanf(fid_rpy, '%f', [1]);
    fgets(fid_rpy);
    bidon=fscanf(fid_rpy,'%c',[2]); %skip 2 first characters of line
    y_tls=fscanf(fid_rpy, '%f', [1]);
    fgets(fid_rpy);
    bidon=fscanf(fid_rpy,'%c',[2]); %skip 2 first characters of line
    z_tls=fscanf(fid_rpy, '%f', [1]);
    
    
    if only_create_batch==0
        while ~feof(fid) %begin parsing through point cloud
            
            
            x=fscanf(fid, '%f', [1]);
            if isempty(x), break, end
            
            y=fscanf(fid, '%f', [1]);
            
            z=fscanf(fid, '%f', [1]);
            deviation= fscanf(fid, '%d', [1]);
            
            r=fscanf(fid, '%f', [1]);
            hit_num=fscanf(fid, '%f', [1]);
            num_hits=fscanf(fid, '%f', [1]);
            if x<xmin | x>xmax, continue,end
            if y<ymin | y>ymax, continue,end
            
            if hit_num>1, continue, end
            refl=10^((r+20)/10);
            d=sqrt((x-x_tls)^2+(y-y_tls)^2+(z-z_tls)^2);
            if d<36
                slope=0.7155+0.04562*d-0.003741*d^2+0.00007449*d^3; %this equ from prism graphpad
                
                refl=refl*slope;
            end
            
            if flag==2
                %I need to integrate generate lvox input code here. TODO
                
                if refl<leaf_noise_threshold | deviation>25
                    
                    fprintf(fidres,'%f %f %f 155 155 155\n',x,y,z);
                    
                    %point is a leaf
                elseif  refl<leaf_wood_threshold
                    fprintf(fidres,'%f %f %f 1 1 1\n',x,y,z);
                    
                    %point is wood
                else %refl>leaf_wood_threshold
                    
                    fprintf(fidres,'%f %f %f 255 255 255\n',x,y,z);
                    
                end
                
            end
            
            if flag==1
                if deviation>25 || d<minimum_dist || d>maximum_dist || (z-bbox_z_min)< distance_above_bbox_z_to_ignore, continue,end %differnces in floor level go to 11 m at Morgan between lower and higher
                
                refl=round(refl);
                
                if pos<=num_scans_day1 %day 1
                    day=1;
                elseif pos<=num_scans_day2 %day 2
                    day=2;
                elseif pos<num_scans_day3%day 3
                    day=3;
                else
                    day=4;
                end
                
                if num_hits>1, tag=2;
                    %2: point is first of numerous hits
                    if refl > 0 && refl <= 100, refl_hist(refl,day,2)=refl_hist(refl,day,2)+1;end
                else tag=1; %1: it is a single hit
                    if refl > 0 && refl <= 100, refl_hist(refl,day,2)=refl_hist(refl,day,2)+1;end
                    if refl > 0 && refl <= 100, refl_hist(refl,day,1)=refl_hist(refl,day,1)+1;end
                end
                
            end
            
            
        end %end of parsing through point cloud
    end
    
    fclose(fid);
    fclose(fid_rpy);
    if flag==2
        
        line=horzcat(plot_name,'_pos',num2str(pos),'_lvox_input.txt Riegl ',num2str(x_tls),' ',num2str(y_tls),' ',num2str(z_tls),' ',num2str(yaw),' ', num2str(min_hor_angle(mod(pos,2)+1)),' ', num2str(max_hor_angle(mod(pos,2)+1)),' ', num2str(min_vert_angle(mod(pos,2)+1)),' ', num2str(max_vert_angle(mod(pos,2)+1)),' ',num2str(res_angu),' ', num2str(mod(pos,2)));
        %Create batch and .in files for lvox
        if only_create_batch==0
            fclose(fidres);
        end
            file=strcat(path,plot_name,'_pos',num2str(pos),'.in');
            fidres_in=fopen(file,'w');
            %
            fprintf(fidres_in,'%d %d %d %d %d %d\n', bbox_x_min, bbox_y_min, bbox_z_min, bbox_x_max, bbox_y_max, bbox_z_max);
            
           
            fprintf(fidres_in,line);
            fclose(fidres_in);
        %end %previous end for if only_create_batch==0
        
        fprintf(fidres_all,line);
        fprintf(fidres_all,'\n');
        line=horzcat('lvox -r',num2str(res),' -i ',plot_name,'_pos',num2str(pos),'.in -o C1.asc -h ',plot_name,'_60m_pos',num2str(pos),'_',num2str(res),'m.hit -I100 --operator -1 -C 100 200 -q1');
        fprintf(fidres_batch,line);
        fprintf(fidres_batch,'\n');
        
        
    end
    
    
    
end

writematrix(refl_hist(:,:,2),fileres_single_and_first);
writematrix(refl_hist(:,:,1),fileres_first);

fclose all;



