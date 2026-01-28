%the program will write four files: one for points, one for topology, one
%for cell types (repetition of the number 12) , and one for voxel values
%once all 4 files are written, they will be combined into one single file
%in the VTK format

%-------------------------INPUT



clear
%if shea tree is used, change lvox_tram to lvox bellow:
% file=strcat('C:\Protected\Postdoc\Data Analysis\Variance in contact numbers\',name);

plot='whatever';
% % mean_lad=0.39;
%
% file=strcat('C:\lvox_',tree,'\',tree,'_10cm_LAD_no_zeros.asc');

% file='F:Protected\lvox_Q4\Q4_10cm_LAD_no_zeros.asc_occlusion_corrected_.asc_recentered.asc';
file='/Users/mbeland/Documents/Research/Professor/FLiESvox outputs/Bartlett/Bartlett_5mgrid_WAD_no_zeros.asc';
%file='/Volumes/LaCie SSD/lvox_SERC/SERC_30cm_60m_planoanderecto_5mgrid_LAD_no_zeros.asc_occlusion_corr_thresh_100.asc';

g_code_present_in_file=0;
lad_present_in_file=1;
header_in_file=0; %flag for header present
bound=0; %flag to use boundary values
%these are there to exclude points outside these boundaries:
minx=0;
miny=0;
maxx=48;
maxy=48;


%-------------------------

%fileres='C:\Protected\Fliesvox\ellipoid.vtk';
% file=strcat('F:\Protected\lvox_Q3\',name);
fid=fopen(file);

%IF header is present:
if header_in_file
    bbox_xmin=fscanf(fid, '%f', [1]);
    bbox_ymin=fscanf(fid, '%f', [1]);
    bbox_zmin=fscanf(fid, '%f', [1]);
    bbox_xmax=fscanf(fid, '%f', [1]);
    bbox_ymax=fscanf(fid, '%f', [1]);
    bbox_zmax=fscanf(fid, '%f', [1]);
    
    vox_size=fscanf(fid, '%f', [1]);
end



fileres_pt=strcat('C:\Temporary\temp_pt_',plot,'.asc');
fileres_topo=strcat('C:\Temporary\temp_topo_',plot,'.asc');
fileres_lad=strcat('C:\Temporary\temp_lad_',plot,'.asc');
fileres=strcat(file,'.vtk');


fidres_pt=fopen(fileres_pt,'w');
fidres_topo=fopen(fileres_topo,'w');
fidres_lad=fopen(fileres_lad,'w');
fidres=fopen(fileres,'w');


num_pts=0;
num_cells=0;
lad=1;

while ~feof(fid)
    xmin=fscanf(fid, '%f', [1]);
    %if xmin>5, break,end %**************************************************************
    if isempty(xmin), break, end
    ymin=fscanf(fid, '%f', [1]);
    zmin=fscanf(fid, '%f', [1]);
    xmax=fscanf(fid, '%f', [1]);
    ymax=fscanf(fid, '%f', [1]);
    zmax=fscanf(fid, '%f', [1]);
    if lad_present_in_file, lad=fscanf(fid, '%f', [1]);else lad=1;end
    if g_code_present_in_file,g_code=fscanf(fid, '%d', [1]);end
    if lad<=0 | isinf(lad), continue,end
    if bound &&(xmin<minx || ymin<miny || xmax>maxx || ymax>maxy), continue, end
    %             if lad==999, lad=mean_lad;end
    
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmin,ymin,zmin);
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmax,ymin,zmin);
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmin,ymax,zmin);
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmax,ymax,zmin);
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmin,ymin,zmax);
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmax,ymin,zmax);
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmin,ymax,zmax);
    fprintf(fidres_pt, '%.3f %.3f %.3f\n',xmax,ymax,zmax);
    %8 points were written in the pt file
    fprintf(fidres_topo, '8 ');
    for i=0:7
        fprintf(fidres_topo, '%lu ',num_pts+i);
    end
    fprintf(fidres_topo, '\n');
    
    num_pts=num_pts+8;;
    num_cells=num_cells+1;
    
    %fprintf(fidres_types, '11\n');
    fprintf(fidres_lad, '%.4f\n',lad);
    
    
    
end

% num_pts
% num_cells
fclose all;
fidres_pt=fopen(fileres_pt);
fidres_topo=fopen(fileres_topo);
%fidres_types=fopen(fileres_types);
fidres_lad=fopen(fileres_lad);
fidres=fopen(fileres,'w');


fprintf(fidres, '# vtk DataFile Version 2.0\n');
fprintf(fidres, 'Voxel visualization_ tree #2\n');
fprintf(fidres, 'ASCII\n');
fprintf(fidres, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fidres, 'POINTS %u float\n',num_pts);
for i=1:num_pts
    temp=fgets(fidres_pt);
    fprintf(fidres, '%s',temp);
end

fprintf(fidres, 'CELLS %lu %lu\n',num_cells,num_cells*9);
for i=1:num_cells
    temp=fgets(fidres_topo);
    fprintf(fidres, '%s',temp);
end

fprintf(fidres, 'CELL_TYPES %lu\n',num_cells);
for i=1:num_cells
    
    fprintf(fidres, '11\n');
end

fprintf(fidres, 'CELL_DATA %lu\n',num_cells);
fprintf(fidres, 'SCALARS LAD float 1\n');
fprintf(fidres, 'LOOKUP_TABLE default\n');
for i=1:num_cells
    temp=fgets(fidres_lad);
    fprintf(fidres, '%s',temp);
end


fclose all;


beep