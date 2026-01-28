%this code prepares the batch file for executing readRXP on all files for a
%given site

%readRXP_MB_binary_class -input EMS_pos2.rxp -trans ScanPos002_EMS.DAT -output EMS.bin -minRefl 8 -maxRefl 29 -bounds -6 -59 324 54 1 374.1

plot='Bartlett';
path=strcat('F:\ReadRXP\voxTrace_',plot,'\');

fileres=strcat(path,plot,'_readRXP_batch.bat');
fidres=fopen(fileres,'w');

start_pos=1;
end_pos=98;

switch plot
    case 'EMS'
        bbox='-6 -59 323 54 1 362';
        minRefl=10;
        maxRefl=40;
    case 'SERC'
        bbox='-53 -2 -16 7 58 47';
        minRefl=10.5;
        maxRefl=40;
    case 'Pasoh'
        bbox='-62 -1 127 -2 59 193';
        minRefl=11;
        maxRefl=45;
    case 'MorganMonroe'
        bbox='-51 -6 239 9 54 299';
        minRefl=10.5;
        maxRefl=40;
    case 'Bartlett'
        bbox='-46 -7 250 2 41 280';
        minRefl=10;
        maxRefl=40;
end
for pos=start_pos:end_pos
       if pos==206 && plot=="SERC", continue,end
    if pos<10
        DAT_name=strcat('ScanPos00',num2str(pos),'.DAT');
    elseif pos<100
        DAT_name=strcat('ScanPos0',num2str(pos),'.DAT');
    else
        DAT_name=strcat('ScanPos',num2str(pos),'.DAT');
    end
    
    input_name=strcat('ScanPos',num2str(pos),'.rxp');
    
    output_name=strcat('ScanPos',num2str(pos),'.bin');
    
    line=strcat("readRXP_MB_binary_class -input ",input_name," -trans ",DAT_name," -output ",output_name," -minRefl ",num2str(minRefl)," -maxRefl ",num2str(maxRefl)," -bounds ",bbox);
    
    fprintf(fidres,'%s\n',line);
end
fclose all;