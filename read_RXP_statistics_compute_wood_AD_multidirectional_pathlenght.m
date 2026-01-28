%2 April 2021: this is a modified version of the LAD computation code used
%to calculate the wood area density for FliesVox. G value of 0.5 will be
%used

%April 15, 2020------------------
% this code is a modified version of the code used to read
%lvox hits to now read statistics computed from the
%wooRaytraceTests_mb_v3_binary.m code which uses RXP outputs from the
%VZ-400 converted to binary files using the readRXP_MB_binary_class.exe
%code

%there are two files containing statistics for each scna position produced by
%wooRaytraceTests_mb_v3_binary.m ,the format of these files is
%for the gaps:
%first line is bbox x min, y min, z min, x max, y max, z max
% second line is voxel size
% from third line on, one line per voxel: xmin y min zmin, xmax ymax zmax, total pathlenght of all pulses that traversed this voxel
%for the hits, the file has no header and each line is associated with one
%voxel statistics:
% xmin y min zmin, xmax ymax zmax, total pathlenght of all pulses up to the
% hit location, number of hits inside the voxel (only leaf hits as
% classified by readRXP_MB_binary_class.exe as 1: leaf, 2: wood or noise
% (these are excluded from the hits statistics, but the pulses that hit wood or resulted in a noise event are
% considered in the gaps statistics up to the entry point in the current
% voxel))
%
%So what this code will do is
%1. for all scans - read each line assiated with a voxel, compute the G funciton for the voxel as seen by the scna position, generate a
%3(xyz min) x 4(gap pathlenght, hit pathlenght, number of hits/G(theta) function, LAD (to compute next step)) x 242 (number of scans) matrix
%containing the summation of all statistics from the different scan positions

%2. For all voxels - do the summation as before:
%N_numerator=N_numerator+ni_l(scan)*(1/G(scan));
%                   N_denominator=N_denominator+probe_lenght(scan);
%probe_lenght(scan) was: ng(scan)*delta_g(scan)+ni_l(scan)*delta_i_l(scan);
%this will be a fourth entry

%3. write to files LAD matrix, pathelenght (required to process occlusion in next code to use), statistics,
%average contribution per scan and LAD stabilization in same format as before
%------------

%There are three sets of coordinates:
% April 2020: this is because in the output matrices, the x,y,z coordinates all start at 0,0,0, and not at bbox min
%So we need to bring the coordinates currently in the statistics files to
%what is called MMC, where 0,0,0 is at the lower corner of the matrix array
%The MSC (sequencial) coordinate system is also needed because when data
%about the array is stored in a matrix the indices need to be integers, so
%to store values in the matrix the MSC need to be used
%TLS metric coordinates (TMC): point cloud
%Matrix metric coordinates (MMC): from bbox borders
%matrix sequencial coordinates (MSC): voxel number in sequencial order in
%x,y,z

clear
tic
%************Begin input section

start_pos=1; %I do not think these can be split between workers
end_pos=98;%242;

plot='Bartlett'
height_instrument=1.7;
grid=5 %this is the distance between scasn positions, used to simulate different field protocols. Values 5 , 10 , 15 or 20
occlusion_threshold=15; %April 2021: this is only to display result, no incidence on output
%plot_size=60;

path_rpy=strcat('F:\ReadRXP\voxTrace_',plot,'\');
%path_rpy=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot,'/');

%path_stats=strcat('F:\ReadRXP\voxTrace_',plot,'\');
path_stats=strcat('F:\ReadRXP\voxTrace_',plot,'\');
%path_stats=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot,'/');
%path_in_files=strcat('F:\lvox_',plot,'\');
%path=strcat('/Volumes/LaCie SSD/lvox_',plot,'/');

%************end of input section

grid_size=strcat('_',num2str(grid),'mgrid');

file_stats_gaps=strcat(path_stats,'ScanPos',num2str(start_pos),'lewos.bin_statistics_gaps.txt');

%nscans=end_pos-start_pos+1

%define leaf angle distribution function (asign selection from list below to g variable)
erecto=[0.001813, 0.014469, 0.038908, 0.072257, 0.11037, 0.148837, 0.182804, 0.208301, 0.222241];
plano=[0.217993772,  0.205566574, 0.181747777, 0.149126381, 0.111758486, 0.07430429,  0.040992495, 0.016051798, 0.002458428];
spheri=[0.013816774, 0.044002464, 0.072956086, 0.099621579, 0.1232069, 0.143184018, 0.158760891, 0.169409487, 0.175041802];
plagio=[0.007084426, 0.053135816, 0.126252094, 0.19247555,  0.220395847, 0.197361602, 0.133842925, 0.059766886, 0.009684853];
uniform=[0.111111111,  0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111, 0.111111111];
% % 
% G2m=[0.040540541,0.02027027,0.087837838,0.060810811,0.108108108,0.155405405,0.128378378,0.175675676,0.222972973];
% G4m=[0.032258065,0,0.048387097,0.129032258,0.107526882,0.080645161,0.161290323,0.166666667,0.274193548];
% G6m=[0.018348624,0.033639144,0.03058104,0.064220183,0.088685015,0.137614679,0.128440367,0.220183486,0.278287462];
% G8m=[0.017182131,0.017182131,0.020618557,0.06185567,0.113402062,0.12371134,0.164948454,0.213058419,0.268041237];
% G10m=[0.025,0.01,0.03,0.065,0.09,0.16,0.155,0.23,0.235];


theta_l=deg2rad([5, 15, 25, 35, 45, 55, 65, 75, 85]);

fid_stats_gaps=fopen(file_stats_gaps);


bbox_xmin=fscanf(fid_stats_gaps, '%f', [1]);
bbox_ymin=fscanf(fid_stats_gaps, '%f', [1]);
bbox_zmin=fscanf(fid_stats_gaps, '%f', [1]);
bbox_xmax=fscanf(fid_stats_gaps, '%f', [1]);
bbox_ymax=fscanf(fid_stats_gaps, '%f', [1]);
bbox_zmax=fscanf(fid_stats_gaps, '%f', [1]);

voxel_size=fscanf(fid_stats_gaps, '%f', [1]);
DimX=double(round((bbox_xmax-bbox_xmin)/voxel_size));
DimY=double(round((bbox_ymax-bbox_ymin)/voxel_size));
DimZ=double(round((bbox_zmax-bbox_zmin)/voxel_size));

plot_size=bbox_xmax-bbox_xmin;

fclose(fid_stats_gaps);

switch grid
    case 5
        for i=start_pos:end_pos
            scan_pos(i)=i;
        end
    case 10
        j=1;
        for i=[1,5,9,13,17,21,45,49,53,57,61,65,89,93,97,101,105,109,133,137,141,145,149,153,177,181,185,189,193,197,221,225,229,233,237,241]
            for above_scan=[0,1]  %use these two lines if grid sizes 10,15 or 20 m is used
                %i=i+above_scan;
                scan_pos(j)=i+above_scan;
                j=j+1;
            end
        end
    case 15
        j=1;
        for i=[1,7,13,19,69,75,81,87,133,139,145,151,201,207,213,219]
            for above_scan=[0,1]  %use these two lines if grid sizes 10,15 or 20 m is used
                %i=i+above_scan;
                scan_pos(j)=i+above_scan;
                j=j+1;
            end
        end
    case 20
        j=1;
        for i=[25,33,41,113,121,129,201,209,217]
            for above_scan=[0,1]  %use these two lines if grid sizes 10,15 or 20 m is used
                %i=i+above_scan;
                scan_pos(j)=i+above_scan;
                j=j+1;
            end
        end
end

% for pos=[1,5,9,13,17,21,45,49,53,57,61,65,89,93,97,101,105,109,133,137,141,145,149,153,177,181,185,189,193,197,221,225,229,233,237,241] %10m grid
% for pos=[1,7,13,19,69,75,81,87,133,139,145,151,201,207,213,219]  %15 m grid
% for pos=[25,33,41,113,121,129,201,209,217]  %20 m grid

%     for above_scan=[0,1]  %use these two lines if grid sizes 10,15 or 20 m is used
%     pos=pos+above_scan;

%use below for single trees
% file_in=strcat(path,tree,'_pos.in'); %to get scan positions

%Use below for clusters:
%file_in=strcat(path_in_files,plot,'.in'); %to get scan positions


% for size=[10 20 30 50 100 200]

% if size==10
%     filename='tram_20m_pos';
% else

% end
%size

%this below reads the hit files, one for each scan
% for i=start_num:nscans-1+start_num
%     %     filename=strcat(plot,'_pos',num2str(i),'_',num2str(size/100),'m');
%     %     %use above for single trees
%
%     %use below for clusters
%     file=strcat(path,plot,'_60m_pos',num2str(i),'_0.3m.hit');
%     % filename=strcat('Q3_20m_pos',num2str(i));
%     %     filename=strcat('tram_20m_pos',num2str(i));
%     %file=strcat(path,filename,'.hit');
%     fid(i)=fopen(file);
%     bbox_xmin=fscanf(fid(i), '%f', [1]);
%     bbox_ymin=fscanf(fid(i), '%f', [1]);
%     bbox_zmin=fscanf(fid(i), '%f', [1]);
%     bbox_xmax=fscanf(fid(i), '%f', [1]);
%     bbox_ymax=fscanf(fid(i), '%f', [1]);
%     bbox_zmax=fscanf(fid(i), '%f', [1]);
%
%     vox_size=fscanf(fid(i), '%f', [1]);
%
%     x_dim=fscanf(fid(i), '%d', [1]);
%     y_dim=fscanf(fid(i), '%d', [1]);
%     z_dim=fscanf(fid(i), '%d', [1]);
%
%     nscans_no_good=fscanf(fid(i), '%d', [1]);
% end


%define different g functions for different heights (maximum 3). The height entered refers to the upper limit of the stratum. Keep in
%mind that z=0 is not the ground level, but the start of the bounding box.
%For example, to define a first stratum going from the gound
%to 4 m above the ground, when the bbox starts 0.7 m above ground, use a
%height of (4-0.7)=3.3
%if start of bbox is 0m in TMC, 4m above ground is 4-1.7=2.3
%height values are in MMC system
%FliesVox codes:
% 1: spherical
% 2: planophile
% 3: erectophile
% 4: uniform. NEW from 2022

%ground is at -1.7 in TMC

%Oct 2018: this below appears incorect to me, one G function is used to
%derive LAD values, and a different G function is used in Fliesvox to
%compute the radiative transfer, this is physically incoherent, both should
%be the same, hence if FliesVox only takes predefined functions, I should
%use the closest one to compute LAD values.

%For wood:
%G=0.5;  nov 2024 now using G functions derived from QSM for EMS and Pasoh
% 
% if strcmp(plot,'MorganMonroe')
%     nG=2; %number of G functions defined
%     
%     height(1)=30; % this is in meters above the floor of the voxel matrix
%     g_height(1,:)=plano;
%     FV_code(1)=2;
%     
%     height(2)=9999;
%     g_height(2,:)=uniform;
%     FV_code(2)=1; %uniform not available in Fliesvox!! april 2021: it is now : code 4
%     
%     leaf_angle='_planoanduniform';
% end
% 
% if strcmp(plot,'SERC')
%     nG=2; %number of G functions defined
%     
%     height(1)=30; % this is in meters above the floor of the voxel matrix
%     g_height(1,:)=plano;
%     FV_code(1)=2;
%     
%     height(2)=9999;
%     g_height(2,:)=erecto;
%     FV_code(2)=3;
%     leaf_angle='_planoanderecto';
% end
% 
if strcmp(plot,'EMS')
    nG=2; %number of G functions defined
    
    height(1)=5; % this is in meters above the floor of the voxel matrix
    g_height(1,:)=erecto;
    FV_code(1)=3;

     height(2)=9999; % this is in meters above the floor of the voxel matrix
    g_height(2,:)=uniform;
    FV_code(2)=4;
    %leaf_angle='_plano';
end
% 
if strcmp(plot,'Pasoh')
    nG=3; %number of G functions defined
    
    height(1)=5; % this is in meters above the floor of the voxel matrix
    g_height(1,:)=erecto;
    FV_code(1)=3;

    height(2)=35; % this is in meters above the floor of the voxel matrix
    g_height(2,:)=uniform;
    FV_code(2)=4;

    height(3)=9999; % this is in meters above the floor of the voxel matrix
    g_height(3,:)=plano;
    FV_code(3)=2;
    %leaf_angle='_uniform';
end

if strcmp(plot,'Bartlett')
    nG=1; %number of G functions defined
    
    height(1)=9999; % this is in meters above the floor of the voxel matrix
    g_height(1,:)=spheri;
    FV_code(1)=1;
    %leaf_angle='_uniform';
end
% %
% height(3)=7;
% g_height(3,:)=G6m;
% FV_code(3)=3;
%
% height(4)=9;
% g_height(4,:)=G8m;
% FV_code(4)=3;
%
% height(5)=99;
% g_height(5,:)=G10m;
% FV_code(5)=3;

%height=height-bbox_zmin-height_instrument; %to convert TMC coord to MMC
%***This is wrong, zmin starts at zero see below it is the variable
%compared to height, so the height is defined from. bbox has no business
%here, and instrument height is not very relevant

%this next line to output results in a different folder not to erase
%current results
%path=strcat(path,'test\')

%I need to get all the voxel size and plot size out of the file snames, this info will be in the statistics file

%filename=strcat(plot,'_',num2str(size),'cm');
file=strcat(path_stats,plot,grid_size,'_WAD_no_zeros.asc');
fidres_lad=fopen(file,'w');
% file=strcat(path,filename,'test2LAD.asc')
% file=strcat('F:\Protected\lvox_tram\tram_20m_10cm_LAD.asc'); %************CHANGED


file=strcat(path_stats,plot,grid_size,'_pathlenght_wood.asc');
% file=strcat(path,filename,'test2fVol.asc')
% file=strcat('F:\Protected\lvox_tram\tram_20m_10cm_fVol.asc');%************CHANGED
fidres_pathlenght=fopen(file,'w');

% file=strcat(path,filename,'_',num2str(plot_size),'m',leaf_angle,grid_size,'_ng_deltag.asc')
% fidres_ng_deltag=fopen(file,'w');

file=strcat(path_stats,plot,grid_size,'_statistics_wood.asc');
fidres_stats=fopen(file,'w');

% file=strcat(path,filename,'_',num2str(plot_size),'m',leaf_angle,grid_size,'_debug.asc')
% fidres_debug=fopen(file,'w');

file=strcat(path_stats,plot,grid_size,'_WAD_stabilization.csv');
fidres_LADstab=fopen(file,'w');

% file=strcat(path,filename,'_',num2str(plot_size),'m',leaf_angle,grid_size,'_average_contribution_per_scan.csv')
% fidres_LADcontrib=fopen(file,'w');

% fid_sc_pos=fopen(file_in);
% fgets(fid_sc_pos);


%***Warning: below may need to be modified if another scanner type is used,
%which changes the number of parameters passed on to lvox. Checlk the .in
%file conmtaining scan positions if getting an error
% for i=start_num:nscans-1+start_num
%     bidon=fscanf(fid_sc_pos,'%s',[1]);
%     bidon=fscanf(fid_sc_pos,'%s',[1]);
%     %tls positions are converted from TLS metric coordinates (TMC) to
%     %Matrix metric coordinates (MMC)
%     x_tls(i)=fscanf(fid_sc_pos,'%f',[1])-bbox_xmin;
%     y_tls(i)=fscanf(fid_sc_pos,'%f',[1])-bbox_ymin;
%     z_tls(i)=fscanf(fid_sc_pos,'%f',[1])-bbox_zmin;
%     bidon_=fscanf(fid_sc_pos,'%f',[7]);
%
% end

count=0;
% count_all=0;
count_occ=0;
% count_occ_=0;
% count_occ_all=0;
%sum_fVol=0;
sum_LA=0;
% sum_G=0;
% count_G=0;

max_probe_total_lenght=500; %in meters
max_num_voxel_LAD_estimates_per_probe_bin=100000;
max_num_voxels=DimX*DimY*DimZ; %3000000;
LAD_deviations=NaN(max_num_voxel_LAD_estimates_per_probe_bin,max_probe_total_lenght);
%probe_lenght_vector=NaN(1,max_num_voxels);
min_total_probe_lenght=2000; %I may want to make this longer depending on the 5m results at SERC runing now
%tot_probe_lenght_per_voxel=NaN(DimX*DimY*DimZ,1);

LAD_values=NaN(DimX*DimY*DimZ,max_probe_total_lenght);

%LAD_contrib_single_scan=zeros(nscans,1);
%number_of_contributions=zeros(nscans,1);

Master=zeros(DimX*DimY*DimZ,4);

g_code=zeros(DimX*DimY*DimZ,1);

%------new code



for pos=scan_pos(1,:)  %begin loop through scan position
    pos
    if pos==206 && plot=="SERC", continue,end
    
    file_stats_gaps=strcat(path_stats,'ScanPos',num2str(pos),'lewos.bin_statistics_gaps.txt');
    file_stats_hits=strcat(path_stats,'ScanPos',num2str(pos),'lewos.bin_statistics_wood_hits.txt');
    fid_stats_gaps=fopen(file_stats_gaps);
    fid_stats_hits=fopen(file_stats_hits);
    
    fgets(fid_stats_gaps);
    fgets(fid_stats_gaps);
    
    
    if pos<10
        file_rpy=strcat(path_rpy,'ScanPos00',num2str(pos),'.RPY');
        
    elseif pos<100
        file_rpy=strcat(path_rpy,'ScanPos0',num2str(pos),'.RPY');
        
    else
        file_rpy=strcat(path_rpy,'ScanPos',num2str(pos),'.RPY');
        
    end
    
    fid_rpy=fopen(file_rpy);
    
    %RPY example format :
    %     Roll=0.207358
    % Pitch=-0.804565
    % Yaw=0.976073
    % X=-50.743374
    % Y=30.720732
    % Z=132.709891
    fgets(fid_rpy);
    fgets(fid_rpy);
    fgets(fid_rpy);
    bidon=fscanf(fid_rpy,'%c',[2]); %skip 4 first characters of line
    x_tls=fscanf(fid_rpy, '%f', [1])-bbox_xmin;
    %=temp;
    fgets(fid_rpy);
    bidon=fscanf(fid_rpy,'%c',[2]); %skip 4 first characters of line
    y_tls=fscanf(fid_rpy, '%f', [1])-bbox_ymin;
    fgets(fid_rpy);
    bidon=fscanf(fid_rpy,'%c',[2]); %skip 4 first characters of line
    z_tls=fscanf(fid_rpy, '%f', [1])-bbox_zmin;
    
    
    %get x,y,z of voxel, convert to 0,0,0 coordinate system, for that voxel compute G
    
    while ~feof(fid_stats_gaps)
        xmin=fscanf(fid_stats_gaps, '%f', [1])-bbox_xmin;
        if isempty(xmin), break, end
        ymin=fscanf(fid_stats_gaps, '%f', [1])-bbox_ymin;
        zmin=fscanf(fid_stats_gaps, '%f', [1])-bbox_zmin;
        xmax=fscanf(fid_stats_gaps, '%f', [1])-bbox_xmin;
        ymax=fscanf(fid_stats_gaps, '%f', [1])-bbox_ymin;
        zmax=fscanf(fid_stats_gaps, '%f', [1])-bbox_zmin;
        
        pathlenght_gaps=fscanf(fid_stats_gaps, '%f', [1]);
        
%         need to reverse this: xmin = (x-1)*vox_size; to get sequencial coord
%         populate matrix,
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
        index=int32(dy+(dx-1)*DimY+(dz-1)*DimX*DimY);
        %Master: 1)gap pathlenght 2) hit pathlenght 3) number of hits/G(theta) 4) LAD
        Master(index,1)=Master(index,1)+pathlenght_gaps;
        
        
        
    end   %end of loop going through all gap pathlenghts for each pulse for the current scan position
    
    
    %assessing what to record for the hits
    %at this stage we have recorded all gap lenghts for all voxels for the
    %given scan, so additional hit pathlenghts will allow producing the
    %total path lenght for the given scan position and any given voxel
    
    %LAD_values is not declared, inefficient code, what is the maximum size of
    %LAD_values? I cannot know, it has lenght of total pathlenght, which is a
    %function of voxel size and scan reosluton.. wrong, it is declared
    %within each voxel, so it's value is reset as the code goes thourgh
    %voxels. This is because i nthe old code verison, we go thourgh each
    %voxel, and then through each scan, so as the code loops thourgh the
    %scan positions it stores the changing LAD estimate. In the new code
    %verison I go through all scan positions and then SOME of the voxels
    %(only those that were seen from that scan position).
    %failure in logic: I loop through all scan positions, and throuhg all
    %voxels contaiing leaves, so I can record. Wait the issue is that the
    %voxels change so I need to record the info for each one to see how it
    %evolves as new scans come in
    %So in order to
    %get the same data as before to derive the graph I use in the paper, I
    %will save the LAD_values for each scan positon and each voxel. perhaps
    %this is not needed, let's look at what is done with LAD_values:
    %for each voxel:
    
    
    %      first_non_NaN_index = find(isnan(probe_lenght_vector), 1);
    %probe_lenght_vector is a vector containing the total probe
    %lenght for each voxel, so its maximum size should be dimx*dimy*dimz
    %this variable is only used at the end to write results.
    %not very informative stat
    %probe_lenght_vector(first_non_NaN_index)=total_probe_lenght;
    %this next one tot_probe_lenght_per_voxel is also only used at the end for writing
    %results, in the new code I have this info in the matrix by adding indice 1 and 2. this is useful
    %for deriving the % of voxels below a given probe lenhgt
    %                 tot_probe_lenght_per_voxel(x,y,z)=total_probe_lenght;
    %
    %                 if total_probe_lenght>min_total_probe_lenght %I am including a voxel only if the total lenght of the probe is above a minimum value to exclude computing statistics on voxel that were not well observed (final LAD estimate not reliable)
    %                     for probe_lenght_=1:max_probe_total_lenght
    %                         if ~isnan(LAD_values(probe_lenght_))
    %                             %find the index of the first NaN in the first
    %                             %index
    %                             first_non_NaN_index = find(isnan(LAD_deviations(:,probe_lenght_)), 1);
    %                             LAD_deviations(first_non_NaN_index,probe_lenght_)=LAD_values(probe_lenght_)-LAD;
    %                         end
    %                     end
    %                 end
    
    
    while ~feof(fid_stats_hits)  %begin of loop going through each pulse with a hit for the current scan position
        xmin=fscanf(fid_stats_hits, '%f', [1])-bbox_xmin;
        if isempty(xmin), break, end
        ymin=fscanf(fid_stats_hits, '%f', [1])-bbox_ymin;
        zmin=fscanf(fid_stats_hits, '%f', [1])-bbox_zmin;
        xmax=fscanf(fid_stats_hits, '%f', [1])-bbox_xmin;
        ymax=fscanf(fid_stats_hits, '%f', [1])-bbox_ymin;
        zmax=fscanf(fid_stats_hits, '%f', [1])-bbox_zmin;
        
        pathlenght_hits=fscanf(fid_stats_hits, '%f', [1]);
        Ni_leaf=fscanf(fid_stats_hits, '%d', [1]);
        
        
        dx=(xmin/voxel_size)+1;
        dy=(ymin/voxel_size)+1;
        dz=(zmin/voxel_size)+1;
        index=int32(dy+(dx-1)*DimY+(dz-1)*DimX*DimY);
        vox_center_x=xmin+voxel_size/2;
        vox_center_y=ymin+voxel_size/2;
        vox_center_z=zmin+voxel_size/2;
        
        for i=1:nG
            if zmin<height(i)
                g=g_height(i,:);
                g_code(index)=FV_code(i);
                break
            end
        end
        theta=acos((vox_center_z-z_tls)/sqrt((x_tls-vox_center_x)^2+(y_tls-vox_center_y)^2+(z_tls-vox_center_z)^2));
        if theta>pi/2,theta=pi-theta;end
        
        G=0.0;
        for i=1:9
            x_=acos(cot(theta)*cot(theta_l(i)));
            if theta <= pi/2-theta_l(i), S=cos(theta)*cos(theta_l(i));
            else S=cos(theta)*cos(theta_l(i))*(1+(2*(tan(x_)-x_)/pi));
            end
            G=G+g(i)*S;
        end
        %Master: 1)gap pathlenght 2) hit pathlenght 3) number of hits/G(theta) 4) LAD
        %I don't think LAD is needed in the Master matrix, it can be
        %computed from all the other indices. MAster is now 3 indexes
        Master(index,2)=Master(index,2)+pathlenght_hits;
        Master(index,3)=Master(index,3)+Ni_leaf/G;
        %careful for confusion on total_probe_lenght, it increases with
        %each encrement through scan positions. This goes also for two
        %previous Master entries
        total_probe_lenght=Master(index,1)+Master(index,2);
        %LAD_values is storing the information about only the given
        %voxel
        if total_probe_lenght<max_probe_total_lenght&&pathlenght_hits>0, LAD_values(index,ceil(total_probe_lenght))=Master(index,3)/total_probe_lenght;end %this is a matrix containing the LAD values changing as new information is added from different scan posisitons
        
    end  %end of loop going through each pulse with a hit for the current scan position
    
    
end  %end of loop going through the scan positions

%At this stage I have a matrix containing total pathlenghts and Ni/G


%End of code new version

% x_dim
% for x=1:x_dim
%     x
%     for y=1:y_dim
%         for z=1:z_dim
%             xmin = (x-1)*vox_size;
%             ymin = (y-1)*vox_size;
%             zmin = (z-1)*vox_size;
%             xmax = xmin+vox_size;
%             ymax = ymin+vox_size;
%             zmax = zmin+vox_size;
%             vox_center_x=xmin+(vox_size/2);
%             vox_center_y=ymin+(vox_size/2);
%             vox_center_z=zmin+(vox_size/2);
%
%             for i=1:nG
%                 if zmin<height(i)
%                     g=g_height(i,:);
%                     g_code=FV_code(i);
%                     break
%                 end
%             end


%             N_numerator=0;
%             N_denominator=0;
%             fVol=0;
%             LAD_values=NaN(1,max_probe_total_lenght);
%             total_probe_lenght=0;
%             total_ng=0;
%             total_deltag=0;
%             contrib_flag=zeros(nscans,1);


%             for scan=start_num:nscans-1+start_num %modifications below are to simulate lower density scaning patterns
%             for scan=[1,5,9,13,17,21,45,49,53,57,61,65,89,93,97,101,105,109,133,137,141,145,149,153,177,181,185,189,193,197,221,225,229,233,237,241] %10m grid
%for scan=[1,7,13,19,69,75,81,87,133,139,145,151,201,207,213,219]  %15 m grid
%for scan=[25,33,41,113,121,129,201,209,217]  %20 m grid

%                      for above_scan=[0,1]  %use these two lines if grid sizes 10,15 or 20 m is used
%                          scan=scan+above_scan;


%                 ni(scan)=fscanf(fid(scan), '%u', [1]);
%                 nt(scan)=fscanf(fid(scan), '%u', [1]);
%                 nb(scan)=fscanf(fid(scan), '%u', [1]);
%
%                 d_bar(scan)=fscanf(fid(scan), '%f', [1]);
%                 delta_b(scan)=fscanf(fid(scan), '%f', [1]);
%                 ni_l(scan)=fscanf(fid(scan), '%u', [1]);
%                 delta_i_l(scan)=fscanf(fid(scan), '%f', [1]);
%                 delta_i_l_bar(scan)=fscanf(fid(scan), '%f', [1]);
%                 ni_w(scan)=fscanf(fid(scan), '%u', [1]);
%                 delta_i_w(scan)=fscanf(fid(scan), '%f', [1]);
%                 delta_i_w_bar(scan)=fscanf(fid(scan), '%f', [1]);
%                 ni_hg(scan)=fscanf(fid(scan), '%u', [1]);
%                 delta_i_hg(scan)=fscanf(fid(scan), '%f', [1]);
%                 delta_i_hg_bar(scan)=fscanf(fid(scan), '%f', [1]);
%                 sel_scan_bidon=fscanf(fid(scan), '%u', [1]);
%
%                 %conditions set to exclude unreliable statistics from
%                 %lvox hit files:
%                 %if ni(scan)==0,continue,end  %this is to avoid increasing probe lenght if nothing is found inside the voxel (lvox statistics have been observed to be unreliable when nothing is found inside)
%                 %if ni(scan)==0&&ni_w(scan)==0&&ni_hg(scan)==0,continue,end  %voxel is certainly inside FOV if these conditions are false. If lvox were working correctly, using the rule nt=0 should allow to determine if the voxel is within the instrumment field of view
%                 if nt(scan)==0, continue, end %this is now determining if the voxel is within the field of view, and replacing the above line
%                 %if ni(scan)==0&&nb(scan)==0,continue,end %this indicates a voxel outside the field of view of the scanner, some statistics for those cases can still increase (or decrease) the probe_lenght variable
%                 if d_bar(scan)==-1,continue,end %this arises when Nt<Nb, it indicates an impossibility which lvox appears to flag with d_bar=-1


%                 ng(scan)=nt(scan)-nb(scan)-ni_l(scan)-ni_w(scan)-ni_hg(scan);
%                 if ng(scan)<0,ng(scan)=0;end %this is to exclude incoherent statistics from lvox
%                 total_ng=total_ng+ng(scan);
%
%                 if ng(scan)~=0
%                     delta_g(scan)=(d_bar(scan)*nt(scan)-delta_b(scan)*nb(scan)-ni_w(scan)*(delta_i_w(scan)+delta_i_w_bar(scan))-ni_l(scan)*(delta_i_l(scan)+delta_i_l_bar(scan))-ni_hg(scan)*(delta_i_hg(scan)+delta_i_hg_bar(scan)))/ng(scan);
%                 else
%                     delta_g(scan)=0;
%                 end
%                 if delta_g(scan)>sqrt(3*vox_size^2),delta_g(scan)=d_bar(scan);end %this is to exclude incoherent statistics from lvox, the choice of assigning d_bar is arbitrary, setting a maximum at vox_size is the theoretical maximum value for a path lenght throuhg a cube
%                 if delta_g(scan)<0,delta_g(scan)=d_bar(scan);end %this has been observed and may happen if Ng is very low (but observed as high as 46 on a small sample), again due to incoherent lvox statistics vs effective number of returns from instrument
%                 total_deltag=total_deltag+ delta_g(scan);
%
%used to investigate calc of probe_lenght
%                     if delta_g(scan)>vox_size
%                         delta_g(scan)
%                         ngscan=ng(scan)
%                         dbar_=d_bar(scan)
%                         nt_=nt(scan)
%                         deltab=delta_b(scan)
%                         nb_=nb(scan)
%                         ni_wood=ni_w(scan)
%                         delta_wood=delta_i_w(scan)+delta_i_w_bar(scan)
%                         ni_leaf=ni_l(scan)
%                         delta_leaf=delta_i_l(scan)+delta_i_l_bar(scan)
%                         ni_noise=ni_hg(scan)
%                         delta_noise=delta_i_hg(scan)+delta_i_hg_bar(scan)
%                     end

%                 sumg=0.0;
%                 theta=acos((vox_center_z-z_tls(scan))/sqrt((x_tls(scan)-vox_center_x)^2+(y_tls(scan)-vox_center_y)^2+(z_tls(scan)-vox_center_z)^2));
%                 if theta>pi/2,theta=pi-theta;end
%
%                 G(scan)=0.0;
%                 for i=1:9
%                     x_=acos(cot(theta)*cot(theta_l(i)));
%                     if theta <= pi/2-theta_l(i), S=cos(theta)*cos(theta_l(i));
%                     else S=cos(theta)*cos(theta_l(i))*(1+(2*(tan(x_)-x_)/pi));
%                     end
%                     G(scan)=G(scan)+g(i)*S;
%
%                 end
%
%                 probe_lenght(scan)=ng(scan)*delta_g(scan)+ni_l(scan)*delta_i_l(scan);
%                 total_probe_lenght=total_probe_lenght+probe_lenght(scan);
%
%                 contrib_flag(scan)=1;
%                 N_numerator=N_numerator+ni_l(scan)*(1/G(scan));
%                 N_denominator=N_denominator+probe_lenght(scan);

%LAD_values is a vector with LAD estimates after each
%contributing scan position. To add this to new code, I
%can record the same vector, whcih is independant of a
%particular voxel or scan position
%                 if probe_lenght(scan)>0&&ni_l(scan)>0,LAD_values(ceil(total_probe_lenght))=N_numerator/N_denominator;end
%                 add=(ng(scan)*delta_g(scan)+ni_l(scan)*delta_i_l(scan)+ni_w(scan)*delta_i_w(scan)+ni_hg(scan)*delta_i_hg(scan))/(nt(scan)*d_bar(scan));
%                 if add>0,fVol=fVol+add;end
%                 count_G=count_G+1;
%                 sum_G=sum_G+G(scan);

%end  %this is the end of above scan loop for 10, 15 and 20 m grids
%             end
%             sel_scan=start_num;
%             max_num_entering=0;
%             for scan=start_num:nscans-1+start_num
%                 if ni(scan)==0,continue,end
%                 if nt(scan)-nb(scan)>max_num_entering, sel_scan=scan;end
%             end

%             N_numerator=0;
%             N_denominator=0;
%
%             fVol=0;
%
%             %oct 2019: it is not clear to me why I need to restart another
%             %loop here, I think this can be fused with above loop. done
%
% %             LAD_values=zeros(1,nscans);
% %             fVol_values=zeros(1,nscans);
%
% LAD_values=NaN(1,max_probe_total_lenght);
%
%                 total_probe_lenght=0;
%probe_lenght_bins=zeros(1,10000);
%             views_count=0;
%             fVol_=0;
%for scan=start_num:nscans-1+start_num   %modifications below are to simulate lower density scaning patterns
%for scan=[1,5,9,13,17,21,45,49,53,57,61,65,89,93,97,101,105,109,133,137,141,145,149,153,177,181,185,189,193,197,221,225,229,233,237,241] %10 m grid
% for scan=[1,7,13,19,69,75,81,87,133,139,145,151,201,207,213,219];  %15 m grid
%for scan=[25,33,41,113,121,129,201,209,217]  %20 m grid

%                 for above_scan=[0,1];
%                    scan=scan+above_scan;

%                     add_=(ng(scan)*delta_g(scan)+ni_l(scan)*delta_i_l(scan)+ni_w(scan)*delta_i_w(scan)+ni_hg(scan)*delta_i_hg(scan))/(nt(scan)*d_bar(scan))
%                     if add_>0,fVol_=fVol_+add_;end

%remeber on ly process if ni gt 1
%      if ni(scan)==0&&ni_w(scan)==0&&ni_hg(scan)==0,continue,end  %voxel is certainly inside FOV if these conditions are false
%if ni(scan)==0,continue,end  %WOW, If I interpret this
%correctly, it means I ignore information from scans
%that did not observed any sort of material inside the
%voxel when calculating LAD, correcting this to
%consider pulses that went through withour hiting
%anything as valid information would decrease LAD
%estimates. This is replace with the following
%conditions (same as used above):
%if ni(scan)==0&&nb(scan)==0,continue,end %this indicates a voxel outside the field of view of the scanner, some statistics for those cases can still increase (or decrease) the probe_lenght variable
%                     if d_bar(scan)==-1,continue,end %this arises when Nt<Nb, it indicates an impossibility which lvox appears to flag with d_bar=-1
%                      total_probe_lenght=total_probe_lenght+probe_lenght(scan);
%                     N_numerator=N_numerator+ni_l(scan)*(1/G(scan));
%                     N_denominator=N_denominator+probe_lenght(scan);
%            views_count=views_count+1;
%        if probe_lenght(scan)>0&&ni_l(scan)>0,LAD_values(ceil(total_probe_lenght))=N_numerator/N_denominator;end



%                     add=(ng(scan)*delta_g(scan)+ni_l(scan)*delta_i_l(scan)+ni_w(scan)*delta_i_w(scan)+ni_hg(scan)*delta_i_hg(scan))/(nt(scan)*d_bar(scan));
%                     if add>0,fVol=fVol+add;end
%                     fVol_values(views_count)=fVol;
%                     count_G=count_G+1;
%                     sum_G=sum_G+G(scan);
%fprintf(fidres_debug, '%.1f %.1f %.1f %.1f %.1f %.1f xtls:%.1f ytls:%.1f ztls:%.1f scan:%d G:%.3f probe lenght:%.3f Ng:%d delta g:%.3f nt:%d d bar:%.3f nb:%d delta b:%.3f ni_l:%d Numerator:%.3f Deno:%.3f LAD:%.3f fVol:%.3f \n',xmin,ymin,zmin,xmax,ymax,zmax,x_tls(scan),y_tls(scan),z_tls(scan),scan,G(scan),probe_lenght(scan), ng(scan),delta_g(scan),nt(scan), d_bar(scan),nb(scan), delta_b(scan),ni_l(scan),N_numerator, N_denominator,N_numerator/N_denominator,fVol);
%                  end
%             end
%
%             LAD=N_numerator/N_denominator;




%             N=ni_l(sel_scan)/(ng*delta_g+ni_l(sel_scan)*delta_i_l(sel_scan));

%             fVol=(ng*delta_g+ni_l(sel_scan)*delta_i_l(sel_scan)+ni_w(sel_scan)*delta_i_w(sel_scan)+ni_hg(sel_scan)*delta_i_hg(sel_scan))/(nt(sel_scan)*d_bar(sel_scan));
%             fVol=(ng*delta_g+ni_l*delta_i_l+ni_hg*delta_i_hg)/(nt(sel_scan)*d_bar);
%fVol=(ng*delta_g+ni_l*delta_i_l)/(nt(sel_scan)*d_bar);

%             if fVol_<0.15, count_occ_all=count_occ_all+1;end
%             count_all=count_all+1;

%             if LAD>0 %*************CHANGED

%loop to record average contribution of each scan to all
%voxel LADs
%                 if total_probe_lenght>occlusion_threshold
%                     for scan=start_num:nscans
%                         if contrib_flag(scan)==0, continue,end
%                         contribution=(N_numerator-ni_l(scan)*(1/G(scan)))/(N_denominator-probe_lenght(scan))-LAD;
%                         if isnan(contribution), contribution=LAD;, end %this means that N_denominator-probe_lenght(scan) is zero, and that
%                         %the given scan is the only contributor to the LAD of the
%                         %given voxel
%                         LAD_contrib_single_scan(scan)=LAD_contrib_single_scan(scan)+contribution;
%                         number_of_contributions(scan)=number_of_contributions(scan)+1;
%                     end
%                 end

%                 first_non_NaN_index = find(isnan(probe_lenght_vector), 1);
%                 %probe_lenght_vector is a vector containing the total probe
%                 %lenght for each voxel, so its maximum size should be dimx*dimy*dimz
%                 %this variable is only used at the end to write results.
%                 %not very informative stat
%                 probe_lenght_vector(first_non_NaN_index)=total_probe_lenght;
%this next one tot_probe_lenght_per_voxel is also only used at the end for writing
%results, in the new code I have this info. this is useful
%for deriving the % of voxels below a given probe lenhgt
%                 tot_probe_lenght_per_voxel(x,y,z)=total_probe_lenght;
%
%                 if total_probe_lenght>min_total_probe_lenght %I am including a voxel only if the total lenght of the probe is above a minimum value to exclude computing statistics on voxel that were not well observed (final LAD estimate not reliable)
%                     for probe_lenght_=1:max_probe_total_lenght
%                         if ~isnan(LAD_values(probe_lenght_))
%                             %find the index of the first NaN in the first
%                             %index
%
%                             first_non_NaN_index = find(isnan(LAD_deviations(:,probe_lenght_)), 1);
%
%                             LAD_deviations(first_non_NaN_index,probe_lenght_)=LAD_values(probe_lenght_)-LAD;
%
%
%
%                         end
%                     end
%                 end
%                 if views_count>40
%                     fprintf(fidres_LADstab, '\n\ndebut, LAD:%.3f\n',LAD);
%                     for i=1:views_count
%                         fprintf(fidres_LADstab, '%.3f ',fVol_values(i));
%                         fprintf(fidres_LADstab, '%.3f\n',(LAD-LAD_values(i)));
%                     end
%                 end


%                 if fVol<0
%                     fVol=0;
%                     %                     ng
%                     %                     delta_g
%                     %                     N
%                 end
%                 %                 if N>0
%                 count=count+1;
%                 sum_fVol=sum_fVol+fVol;
%                 end

%the following values (xmin, xmax, etc) are in meters from the bbox border
%(matrix metric coordinates, MMC)
%x,y,z are matrix sequencial coordinates (MSC)




%                 z_tls(sel_scan); ????


%                 LAD=N/G;
%if LAD>10,LAD,end
%                 if fVol_<0.15, count_occ_=count_occ_+1;end

%November 2019: replacing the rule below used to determine if voxel is
%occluded, now using number of pulses comoputed as total path lenght/voxel
%size. I first use 50 pulses (15m/0.3m), I am trying to end up in the same
%ballpark as when I used fvol for the 4 sites
%                 if fVol<0.15
%                     count_occ=count_occ+1;
% %                     LAD=0;
%                     LAD=999;
%
%                 end


%                 if total_probe_lenght<occlusion_threshold
%                     count_occ=count_occ+1;
%                     %    LAD=999; %Occlusion is no longer detemined here, but in generate_vertical_profile...
%
%                 end
%
%                 sum_LA=sum_LA+LAD*vox_size^3
%                 fprintf(fidres_lad_10, '%.3f %.3f %.3f %.3f %.3f %.3f %f %d\n',xmin,ymin,zmin,xmax,ymax,zmax,LAD,g_code);
%                 fprintf(fidres_pathlenght_10, '%.3f %.3f %.3f %.3f %.3f %.3f %.2f\n',xmin,ymin,zmin,xmax,ymax,zmax,total_probe_lenght);
%                 fprintf(fidres_ng_deltag, '%.3f %.3f %.3f %.3f %.3f %.3f %d %.2f\n',xmin,ymin,zmin,xmax,ymax,zmax,total_ng, total_deltag);
%                 %             else
%                 %                 fprintf(fidres_lad_10, '%.3f %.3f %.3f %.3f %.3f %.3f 0 %d\n',xmin,ymin,zmin,xmax,ymax,zmax,g_code);
%             end
%         end
%     end
% end

%Need to include in processing of final matrix Master
%Master matrix now contains, for each voxel, the sum of gaps pathlenghts,
%hit pathlenghts, and Ni/G

%  if total_probe_lenght<occlusion_threshold
%                     count_occ=count_occ+1;
%                     %    LAD=999; %Occlusion is no longer detemined here, but in generate_vertical_profile...
%
%  end

tot_probe_lenght_per_voxel=(Master(:,1)+Master(:,2));

voxels_containing_material=find(Master(:,3)>0&tot_probe_lenght_per_voxel(:,1)>0); %this is because pathlenght hits is recorded with 3 decimals (mm) and it is possible some total pathlenght values are below one milimeter, resulting in a division by zero when computing LAD
count=size(voxels_containing_material,1)




occ=find(tot_probe_lenght_per_voxel(voxels_containing_material,1)<occlusion_threshold);
%occ=find(Master(occ,3)>0);  %this is because the definition of an occluded voxel is a voxel within which material was detected and is below the pathlenght threshold
count_occ=size(occ,1)
proportion_of_occ=count_occ/count

%To get voxel coordinates:
%used previously:
%          dz=ceil(idx/(dimX*dimY));
%             remain=idx-(dz-1)*dimX*dimY;
%             dx=ceil(remain/dimY);
%             dy=remain-(dx-1)*dimY;
%             vx = bbox_x_min + [(dx-1)*voxel_size (dx)*voxel_size];
%             vy = bbox_y_min + [(dy-1)*voxel_size (dy)*voxel_size];
%             vz = bbox_z_min + [(dz-1)*voxel_size (dz)*voxel_size];

% for Ni_leaf_over_G=nonzeros(Master(:,3))
%     count=size(Ni_leaf_over_G,1);
%     %              sum_LA=sum(Ni_leaf_over_G*voxel_size^3)
% end


Master(voxels_containing_material,4)=Master(voxels_containing_material,3)./tot_probe_lenght_per_voxel(voxels_containing_material,1); %this is the LAD of each voxel
sum_LA=sum(Master(:,4));

for idx=(voxels_containing_material).'  %this last .' is to convert the vector from a column to a line, so that the for loop considers each element of idx seperately instead of passing idx as one vector
    dz=ceil(idx/(DimX*DimY));
    remain=idx-(dz-1)*DimX*DimY;
    dx=ceil(remain/DimY);
    dy=remain-(dx-1)*DimY;
    vx = [(dx-1)*voxel_size (dx)*voxel_size];
    vy = [(dy-1)*voxel_size (dy)*voxel_size];
    vz = [(dz-1)*voxel_size (dz)*voxel_size];
    fprintf(fidres_lad, '%.3f %.3f %.3f %.3f %.3f %.3f %f\n',vx(1),vy(1),vz(1),vx(2),vy(2),vz(2),Master(idx,4));
    fprintf(fidres_pathlenght, '%.3f %.3f %.3f %.3f %.3f %.3f %.2f\n',vx(1),vy(1),vz(1),vx(2),vy(2),vz(2),tot_probe_lenght_per_voxel(idx));
    if tot_probe_lenght_per_voxel(idx)>min_total_probe_lenght
        for probe_lenght_=1:max_probe_total_lenght
            if ~isnan(LAD_values(idx,probe_lenght_))
                %find the index of the first NaN in the first index
                first_non_NaN_index = find(isnan(LAD_deviations(:,probe_lenght_)), 1);
                LAD_deviations(first_non_NaN_index,probe_lenght_)=LAD_values(idx,probe_lenght_)-Master(idx,4);
            end
        end
    end
end
%             sum_LA=sum_LA+LAD*voxel_size^3

% fprintf(fidres_lad, '%.3f %.3f %.3f %.3f %.3f %.3f %f %d\n',xmin,ymin,zmin,xmax,ymax,zmax,LAD,g_code);
% fprintf(fidres_pathlenght, '%.3f %.3f %.3f %.3f %.3f %.3f %.2f\n',xmin,ymin,zmin,xmax,ymax,zmax,total_probe_lenght);
% fprintf(fidres_ng_deltag, '%.3f %.3f %.3f %.3f %.3f %.3f %d %.2f\n',xmin,ymin,zmin,xmax,ymax,zmax,total_ng, total_deltag);

% fprintf(fidres_stats, 'Plot: %s voxel size: %.1f grid size: %s number of G functions:%d\n',plot, voxel_size,grid_size,nG);
% for i=1:nG
%     fprintf(fidres_stats, 'G function %d: %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f until %d meters\n',i,g_height(i,:), height(i));
% end
fprintf(fidres_stats, 'proportion of occlusion %.4f\n',count_occ/count);
fprintf(fidres_stats, 'number of occluded voxels %d\n',count_occ);
fprintf(fidres_stats, 'number of voxels containing material (detected) %d\n',count);
% fprintf(fidres_stats, 'Average fVol %.4f\n',sum_fVol/count);

% fprintf(fidres_stats, 'Average probe lenght per voxel %.4f\n',mean(probe_lenght_vector,'omitnan'));
% fprintf(fidres_stats, 'Standard deviation in probe lenght per voxel %.4f\n',std(probe_lenght_vector,'omitnan'));
fprintf(fidres_stats, 'total leaf area %.1f\n',sum_LA);
fprintf(fidres_stats, 'LAI %.3f\n',sum_LA/plot_size^2);
% fprintf(fidres_stats, 'volume occluded %d\n',count_occ*voxel_size^3);
% fprintf(fidres_stats, 'count, sum and average G %d, %d, %d\n',count_G,sum_G, sum_G/count_G);
%next 3 lines not very useful
% average_total_probe_lenght=mean(tot_probe_lenght_per_voxel,'all','omitnan');
% st_dev_total_probe_lenght=std(tot_probe_lenght_per_voxel(:),'omitnan');
% fprintf(fidres_stats, 'average_total_probe_lenght per voxel: %.2f  st_dev_total_probe_lenght per voxel: %.2f\n',average_total_probe_lenght,st_dev_total_probe_lenght);

%this next line, use matrix addition indices 1 and 2 instead of tot_probe_lenght_per_voxel
%need to interpret this line in new code, calc from new LAD_values matrix:
%LAD_value is now stored as : LAD_values(ceil(total_probe_lenght),index)=Master(index,3)/total_probe_lenght;



for i=[5,10,15,20,25,50,75,100,200,300,400,500],fprintf(fidres_stats, 'Percent of voxels having a total probe lenght lower than %d: %.1f\n',i,size(find(tot_probe_lenght_per_voxel(voxels_containing_material,1)<i),1)/count*100);end

averages=mean(LAD_deviations,'omitnan');
stand_dev=std(LAD_deviations,'omitnan');
for probe_lenght_=1:max_probe_total_lenght,fprintf(fidres_LADstab, '%.1f,%.3f,%.3f\n',probe_lenght_-0.5,averages(probe_lenght_),stand_dev(probe_lenght_));end

% fprintf(fidres_LADcontrib,'scan, total contribution, number of contributions, average contribution\n');
% for scan=start_num:nscans
%     fprintf(fidres_LADcontrib,'%d, %.3f, %.0f, %.3f\n',scan,LAD_contrib_single_scan(scan),number_of_contributions(scan),LAD_contrib_single_scan(scan)/number_of_contributions(scan));
% end

%all good from here



% proportion_of_occ_considering_data_from_no_hit_scans=count_occ_/count
% proportion_of_occ_all=count_occ_all/count_all
% average_fVol=sum_fVol/count
total_leaf_area=sum_LA
LAI=total_leaf_area/(plot_size^2)
%  LAI_uncorrected=sum_LA/(x_dim*y_dim*vox_size^2)
% volume_occluded=count_occ*vox_size^3
% count_G_=count_G
% sum_G_=sum_G



fclose all;
%end
toc
beep
