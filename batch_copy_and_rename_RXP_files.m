%this code gets all RXP files for a site and copies them to a unique
%directory and renames them for processing by readRXP

%SERC issue:
% Riscan - correct position	
%rule: if pos(in riscan project) >=206 and pos <=230, asign to pos+1 in correct grid, if pos = 231, ignore this
%position as it was not acquired, so do not include in RXP export, and add
%an empty file for statistics results on gaps and hits (see excel file in
%lvox_SERC folder on 4tb drive) april 14, 2020


plot='Bartlett';
path=strcat('F:\ReadRXP\voxTrace_',plot,'\',plot,'.RiSCAN\SCANS\');

plot_target='Bartlett';
dest=strcat('F:\ReadRXP\voxTrace_',plot_target,'\');

start_pos=1;
end_pos=98;

for pos=start_pos:end_pos
   
    if pos==231 && plot=="SERC", continue,end
   
       
    if pos<10
        path2=strcat(path,'ScanPos00',num2str(pos),'\SINGLESCANS\');
    elseif pos<100
        path2=strcat(path,'ScanPos0',num2str(pos),'\SINGLESCANS\');
    else
        path2=strcat(path,'ScanPos',num2str(pos),'\SINGLESCANS\');
    end
    cd (path2)
    
    if pos >=206 && pos <=230 && plot=="SERC"
       %deal with scan pos sequence issue
       pos=pos+1
    end
   
    files = dir('*.rxp'); 
    for id = 1:length(files)
        % Get the file name (minus the extension)
        [~, f] = fileparts(files(id).name);
        
        % Convert to number
        %num = str2double(f(end-1));
        if ~isletter(f(end))
            % If numeric, rename and copy
            
            copyfile(files(id).name, sprintf('%sScanPos%d.rxp', dest,pos));
            
        end
    end
    
end
  
