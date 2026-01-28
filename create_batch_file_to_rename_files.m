fileres=strcat('C:\Temporary\for histograms\batch_rename.bat');
fidres=fopen(fileres,'w');

for i=1:98
    if i<10
fprintf(fidres, 'rename ScanPos00%d* ScanPos00%d.ascii\n',i,i);
    elseif i<100
        fprintf(fidres, 'rename ScanPos0%d* ScanPos0%d.ascii\n',i,i);
    else
        fprintf(fidres, 'rename ScanPos%d* ScanPos%d.ascii\n',i,i);
    end
end
fclose all;