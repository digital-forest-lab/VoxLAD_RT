%this program reads an SOP matrix produced by RiScan for transforming
%vertical to horizontal positions. Then reads all SOP matrices
%corresponding to odd scans and applies a combination of two
%transformations to be applied to all even scans (after automatic
%registration is performed)

%Martin Beland


%_______ input
path='F:\ReadRXP\voxTrace_Bartlett\Bartlett.RiSCAN\';
name='ScanPos002.dat';

start_scan=1;
num_scans=98;

%---- end of input

file=strcat(path,name);
fid_vert_to_hor=fopen(file);

M1=zeros(4,4);
M2=zeros(4,4);

for i=1:4
         M1(i,:)= fscanf(fid_vert_to_hor, '%f', [4]);
end
    

for scan=start_scan:2:num_scans
    if scan<10, name=strcat('ScanPos00',num2str(scan),'.DAT'); elseif scan<100 name=strcat('ScanPos0',num2str(scan),'.DAT'); else name=strcat('ScanPos',num2str(scan),'.DAT');end
    
    file=strcat(path,name);
    fid_scan=fopen(file);
    
    for i=1:4
         M2(i,:)= fscanf(fid_scan, '%f', [4]);
    end
    if scan<10-1, name=strcat('ScanPos00',num2str(scan+1),'.DAT'); elseif scan<100-1 name=strcat('ScanPos0',num2str(scan+1),'.DAT'); else name=strcat('ScanPos',num2str(scan+1),'.DAT');end
    file_res=strcat(path,name);
    fid_res=fopen(file_res,'w');
    M3=M2*M1;
    for i=1:4
        for j=1:4
            fprintf(fid_res, '%.15f ',M3(i,j));
        end
        fprintf(fid_res, '\n');
    end
end
fclose all;
