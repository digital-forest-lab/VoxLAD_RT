%generate vertical foliage density profile


clear

%inputs:
%------------------------


%setting a value for maximum allowable lad value:
max_lad=56; %this corresponds to about 20 leaves per 10x10x10 cm voxel, area of 1 leaf is about 0.0028m2 , 20 leaves is 0.056 m2

voxel_size=0.3;

max_tree_height=120; %this is only for setting a max size on a vector

occlusion_threshold=30;

number_of_voxels_above_ground_to_exclude=4; %this aims to exclude the ground in the leaf area matrix

ground_tile_size=3; %this is a gathering of 3x3 voxel tiles to find the ground, it serves to avoid holes in the ground coverage, i.e. voxel columns within which the ground was not hit by any point. the number needs to be divisible by voxel size!

%layer_thickness_for_voxel_clumping=0; %meters. this is the thickness of the upper layer on which to apply within voxel clumping
%within_voxel_clumping_factor=0.4; %only for upper layer defined above,
%these two above are now defined using an equation

file_will_be_resampled=0; %using 0 does not include lad=0 points, using 1 does include them. Nov 2019: why would I want to include them is not clear. April 2024: this is likely to produce the LAI profile on the first run to output the profile variable which will be used in Prism to fit a function, the code is then run a second time with the appropriate fitting parameters. No this is wrong, I don't remember what this is for

plot_size=48; %was 60 %[m] this is used only for the cluster, irrelevant for individual trees

%path='/Volumes/LaCie SSD/lvox_MorganMonroe/';
%for plot=["EMS" "MorganMonroe" "SERC" "Pasoh"];
for plot=["Bartlett"];
    %for grid=["5" "10" "15" "20"];
    for grid=["5"];
        
        path=strcat('F:\ReadRXP\voxTrace_',plot,'\');
        %path=strcat('E:\',plot,'\');
        %path=strcat('/Volumes/LaCie SSD/ReadRXP/voxTrace_',plot,'/');
        
        %------------------------------
        switch plot
             case {"Bartlett"}
                poisson=0.85;
%                 plot_max_tree_height=28;
%                 min_omega=0.38; %this is the within voxel clumping at the canopy top
                %plateau=14.05;
                plateau=1.442;
                y0=-2.65;
                k=9.619;
                x0=0.3396;
            case {"SERC"} 
                poisson=0.9; % this is the correction to apply to LAD based on the leaf size vs voxel size relation as per Beland et al 2014 AFM. it is 1/the overestimation provided by eq 5 in the AFM paper
                %plateau=-1;
                plateau=1.442;
                y0=-2.65;
                k=9.619;
                %x0=0.3396;
                x0=0.4; %modified nov 2022 to get LAI closer to leaf litter value
            case {"Pasoh"}
                 poisson=0.9;
                 
                 plateau=1.442;
                y0=-2.65;
                k=9.619;
                x0=0.3396;
            case {"EMS"}
                poisson=0.85;
%                 plot_max_tree_height=28;
%                 min_omega=0.38; %this is the within voxel clumping at the canopy top
                %plateau=14.05;
                plateau=1.442;
                y0=-2.65;
                k=9.619;
                x0=0.3396;
            case {"MorganMonroe"}
                poisson=0.85;
                %                 plot_max_tree_height=36;
                %                 min_omega=0.32; %this is the within voxel clumping at the canopy top
                %                 plateau=22.45;
                plateau=1.442;
                y0=-2.65;
                k=9.619;
                x0=0.3396;
        end
        %K=3.258; %this was for first attempt liniking height with omega
        %voxel
        
        file=strcat(path,plot,'_',grid,'mgrid_LAD_no_zeros.asc')  %for Mac
        fid=fopen(file);
        file_pathlenght=strcat(path,plot,'_',grid,'mgrid_pathlenght.asc')  %for Mac
        fid_pathlenght=fopen(file_pathlenght);
        
        
        
        % ground_level=-9999; %Leaf area below this value will be removed from the matrix and will not contribute to canopy LAI
        
        %Nov 2019: editing this code to add data on height vs mean LAD per voxel
        %and percent of occluded voxels per height level
        
        %march 14 2020: with current 3 x 3m tiles to map the ground, there appears
        %to be some ground included in the leaf area sums, this could be because of
        %a slope within the 3x3 m area? increasing the height above ground not to
        %consider
        
        % file=strcat('F:Protected\lvox_M7\M7_10cm_LAD_no_zeros.asc')
        %file=strcat('d:\lvox_MorganMonroe\MorganMonroe_30cm_LAD_no_zeros.asc');
        
        
        
        %to view the profile: the variable profile contains the profile before
        %correction, and profile_ after correction. Copy the variable content into
        %prism to view and create plot
        
        %to compute LAI: (in meters)
        
        ground_level=zeros(plot_size/ground_tile_size+1,plot_size/ground_tile_size+1);
        ground_level(:)=20;% this is the z level above which no ground is expected , changed from 10 to 20 because of elevation variability at some sites
        ceiling_level=zeros(plot_size/ground_tile_size+1,plot_size/ground_tile_size+1);
        canopy_height_model=zeros(plot_size/voxel_size+1,plot_size/voxel_size+1);
        
        lad_above_max=int16.empty(0,2);  %index 1 is the lad value, 2 is the path lenght value
        above_idx=1;
        
        for vox_size=[voxel_size];
            vox_size
            
            
            % if vox_size==.1
            %     name='tram_20m_pos';
            % else
            %     name=strcat('Q3_20m_',num2str(vox_size*100),'cm');
            % end
            
            % vox_size=1;
            
            %--------------------------%____________end inputs section
            % file=strcat('F:Protected\lvox_Q3\',name,'_LAD.asc');
            
            
            
            % file=strcat('F:Protected\lvox_Q3\',name,'_occlusion_corrected_',num2str(file_will_be_resampled),'.asc');
            fileres=strcat(file,'_occlusion_corr_thresh_',num2str(occlusion_threshold),'.asc'); %I should not write to LAD file the real density including clumping, because in radiative transfer terms, it is the effective leaf area density that should be used
            %******************************** When using within voxel clumping, use below not to change file to use for RT modeling (which cannot consider within voxel clumping at the moment):
            %fileres=strcat(file,'test',num2str(occlusion_threshold),'.asc');
            fidres=fopen(fileres,'w');
            fileres_omega=strcat(file,'_occlusion_corr_thresh_',num2str(occlusion_threshold),'with_omega_v.asc');
            
            file_ground=strcat(file,'_ground_level_.asc');
            fid_ground=fopen(file_ground,'w');
            file_height=strcat(file,'_canopy_height_model_.asc');
            fid_height_model=fopen(file_height,'w');
            
            
            
            sum_la=zeros(max_tree_height*10,1); %times 10 means that bins are every decimeters
            sum_lad=zeros(max_tree_height*10,1);  % this is the sum of LAD values for each vertical bin
            count=zeros(max_tree_height*10,1);  %this is the number of voxels for which LAD>0, careful this EXCLUDES the occluded voxels
            profile=zeros(max_tree_height*10,1);  %this is the mean LAD per vertical layer
            occluded=zeros(max_tree_height*10,1);  %this is the number of occluded voxels per vertical layer
            
            %zmin_ variable will become the height above the lowest z value for the
            %given column (set by ground_level)
            
            while ~feof(fid)
                xmin=fscanf(fid, '%f', [1]);
                if isempty(xmin), break, end
                ymin=fscanf(fid, '%f', [1]);
                zmin=fscanf(fid, '%f', [1]);
                xmin_=round(xmin/ground_tile_size)+1; %these are used to find the minimum z value for each x,y column
                ymin_=round(ymin/ground_tile_size)+1;
                if zmin<ground_level(xmin_,ymin_), ground_level(xmin_,ymin_)=zmin;end
                
                xmax=fscanf(fid, '%f', [1]);
                ymax=fscanf(fid, '%f', [1]);
                zmax=fscanf(fid, '%f', [1]);
                if zmax>ceiling_level(xmin_,ymin_), ceiling_level(xmin_,ymin_)=zmax;end
                if zmax>canopy_height_model(round(xmin/voxel_size)+1,round(ymin/voxel_size)+1), canopy_height_model(round(xmin/voxel_size)+1,round(ymin/voxel_size)+1)=zmax;end
                lad=fscanf(fid, '%f', [1]);
                g_code=fscanf(fid, '%d', [1]);
                
            end
            frewind(fid);
            
            for x=0:plot_size/ground_tile_size-1
                for y =0:plot_size/ground_tile_size-1
                    xmin=x*ground_tile_size;
                    ymin=y*ground_tile_size;
                    fprintf(fid_ground, '%.3f %.3f %.3f %.3f %.3f %.3f 1\n',xmin,ymin,ground_level(x+1,y+1),xmin+ground_tile_size,ymin+ground_tile_size,ground_level(x+1,y+1)+0.1);
                end
            end
            for x=0:plot_size/voxel_size-1
                for y =0:plot_size/voxel_size-1
                    xmin=x*voxel_size;
                    ymin=y*voxel_size;
                    fprintf(fid_height_model, '%.3f %.3f %.3f\n',xmin+voxel_size/2,ymin+voxel_size/2,canopy_height_model(x+1,y+1)-voxel_size/2);
                end
            end
            
            
            %14 march 2020: problem with this algorithm: it assumes I see the ground
            %everywhere, which is not the case, result is that for some columns the
            %ground may be 4-5 m up. Possible solution: degrade the resolution of the
            %ground matrix to 1 x 1 m tiles. was degraded to 3 x 3 m on the basis of
            %tests on Morgan at 15 m scan pattern
            
            while ~feof(fid)
                xmin=fscanf(fid, '%f', [1]);
                if isempty(xmin), break, end
                ymin=fscanf(fid, '%f', [1]);
                zmin=fscanf(fid, '%f', [1]);
                xmin_=round(xmin/ground_tile_size)+1; %these are used to find the minimum z value for each x,y column
                ymin_=round(ymin/ground_tile_size)+1;
                zmin_=round((zmin-ground_level(xmin_,ymin_))*10)+1;%this operation is made to create a vector with z values every decimeter, 0.0m becomes 1, 0.5m becomes 6, 1.5m becomes 16
                %             zmin=round(fscanf(fid, '%f', [1])*10)+1;
                
                
                xmax=fscanf(fid, '%f', [1]);
                ymax=fscanf(fid, '%f', [1]);
                zmax=fscanf(fid, '%f', [1]);
                lad=fscanf(fid, '%f', [1]);
                g_code=fscanf(fid, '%d', [1]);
                
                %get path lenght value for each voxel
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                total_path_lenght=fscanf(fid_pathlenght, '%f', [1]);
                
                
                %             if zmin<ground_level, continue, end %voxel is at or below ground level
                if total_path_lenght<occlusion_threshold, occluded(zmin_)=occluded(zmin_)+1;end
                if total_path_lenght<occlusion_threshold || lad <= 0,continue,end
                
                lad=lad*poisson;
                if lad > max_lad
                    lad_above_max(above_idx,1)=lad;
                    lad_above_max(above_idx,2)=total_path_lenght;
                    above_idx=above_idx+1;
                    lad=max_lad;
                end
                
                la=lad*vox_size^3;
                sum_la(zmin_)=sum_la(zmin_)+la;  %DONE 14 march 2020, TODO : in order to assign corect height above ground, I would need to do a first loop only to find the ground values, then do this loop to record the leaf area density values, then the next loop to assign the average LAd values to occluded voxels
                sum_lad(zmin_)=sum_lad(zmin_)+lad;
                count(zmin_)=count(zmin_)+1;
            end
            profile=sum_lad./count;
            total_leaf_area_before_correction=sum(sum_la)
            frewind(fid);
            frewind(fid_pathlenght);
            
            %this next pass writes the file corrected
            while ~feof(fid)
                xmin=fscanf(fid, '%f', [1]);
                if isempty(xmin), break, end
                ymin=fscanf(fid, '%f', [1]);
                zmin=fscanf(fid, '%f', [1]);
                xmin_=round(xmin/ground_tile_size)+1; %these are used to find the minimum z value for each x,y column
                ymin_=round(ymin/ground_tile_size)+1;
                %zmin_=round(zmin*10)+1;
                zmin_=round((zmin-ground_level(xmin_,ymin_))*10)+1;
                
                xmax=fscanf(fid, '%f', [1]);
                ymax=fscanf(fid, '%f', [1]);
                zmax=fscanf(fid, '%f', [1]);
                lad=fscanf(fid, '%f', [1]);
                g_code=fscanf(fid, '%d', [1]);
                
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                bidon=fscanf(fid_pathlenght, '%f', [1]);
                total_path_lenght=fscanf(fid_pathlenght, '%f', [1]);
                
                if zmin<ground_level(xmin_,ymin_)+voxel_size*number_of_voxels_above_ground_to_exclude, continue, end %voxel is at or below ground level
                %             if zmin<ground_level, continue, end %voxel is at or below ground level
                if total_path_lenght<occlusion_threshold
                    lad=profile(zmin_);
                else
                    lad=lad*poisson;
                end
                
                %height_above_ground=zmax-ground_level(xmin_,ymin_);
                
                %omega_v=(-log((height_above_ground-plateau)/(plot_max_tree_height-plateau))/K)+min_omega;
                %no longer used, now omega is detemined from vertical
                %cumulative LAI below in next loop
%                 if height_above_ground>plateau
%                     omega_v=(-log((height_above_ground-plateau)/(plot_max_tree_height-plateau))/K)+min_omega;
%                     if omega_v>1,omega_v=1;end
%                     if omega_v<min_omega,omega_v=min_omega;end
%                 else
%                     omega_v=1;
%                 end
                
                %if zmax>ceiling_level(xmin_,ymin_)-layer_thickness_for_voxel_clumping,lad=lad/within_voxel_clumping_factor;end %max ceiling was used to get the max tree height locally, then assign voxel clumping at a distance down from there, the max tree height is now fixed for a plot
                if lad > max_lad
                    
                    lad=max_lad;
                end
                if isnan(lad)==1
                    lad=0;
                    meassage='error, lad is NaN, caused by absence of unoccluded voxels on the horizontal layer'
                end
                switch file_will_be_resampled
                    %                     case 0
                    %                         if lad>0, fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f %f %f %d\n',xmin,ymin,zmin,xmax,ymax,zmax,lad,omega_v,g_code);end
                    %                     case 1
                    %                         fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f %f %f %d\n',xmin,ymin,zmin,xmax,ymax,zmax,lad,omega_v,g_code);
                    case 0
                        if lad>0, fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f %f %d\n',xmin,ymin,zmin,xmax,ymax,zmax,lad,g_code);end
                    case 1
                        fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f %f %d\n',xmin,ymin,zmin,xmax,ymax,zmax,lad,g_code);
                end
                
                
            end
            fclose all;
            
            %this loop to generate the vertical profile of leaf area,
            %correted for occlusion, to use as a basis for asigning omega
            %voxel (function of vertical cumulative LAI)
             fid=fopen(fileres); %

            
            sum_la_=zeros(max_tree_height*10,1);
            %sum_la_omega=zeros(max_tree_height*10,1);
            sum_lad_=zeros(max_tree_height*10,1);
            %sum_lad_omega=zeros(max_tree_height*10,1);
            count_=zeros(max_tree_height*10,1);
            
            % this next pass is there to create a vertical LAD profiles with the corrected LAD values for
            % occlusion
            while ~feof(fid)
                xmin=fscanf(fid, '%f', [1]);
                if isempty(xmin), break, end
                ymin=fscanf(fid, '%f', [1]);
                zmin=fscanf(fid, '%f', [1]);
                xmin_=round(xmin/ground_tile_size)+1; %these are used to find the minimum z value for each x,y column
                ymin_=round(ymin/ground_tile_size)+1;
                %zmin_=round(zmin*10)+1;
                zmin_=round((zmin-ground_level(xmin_,ymin_))*10)+1;
                
                xmax=fscanf(fid, '%f', [1]);
                ymax=fscanf(fid, '%f', [1]);
                zmax=fscanf(fid, '%f', [1]);
                lad=fscanf(fid, '%f', [1]);
                %omega_v=fscanf(fid, '%f', [1]);
                g_code=fscanf(fid, '%d', [1]);
                if zmin<ground_level(xmin_,ymin_)+voxel_size*number_of_voxels_above_ground_to_exclude, continue, end %voxel is at or below ground level
                la=lad*vox_size^3;
                %lad_omega=lad/omega_v;
                %la_omega=lad_omega*vox_size^3;
                
                sum_la_(zmin_)=sum_la_(zmin_)+la;
                %sum_la_omega(zmin_)=sum_la_omega(zmin_)+la_omega;
                sum_lad_(zmin_)=sum_lad_(zmin_)+lad;
                %sum_lad_omega(zmin_)=sum_lad_omega(zmin_)+lad_omega;
                
                count_(zmin_)=count_(zmin_)+1;
            end
            LA_profile_for_LAI=(sum_la_);  %this _ symbol refers to the corrected values for occlusion
            
            fclose all;
            
            %create profile for corrected file:
            % file=strcat('F:Protected\lvox_Q3\',name,'_occlusion_corrected_',num2str(file_will_be_resampled),'.asc');
            fid=fopen(fileres); %
            
            fileres2=strcat(file,'_statistics.asc');
            %******************************** TESTING:
            %fileres2=strcat(file,'_tests.asc');
            fidres_stats=fopen(fileres2,'w');
            fidres_omega=fopen(fileres_omega,'w'); %don't understand, this appears not used, doesn't it erase the file?
            
            sum_la_=zeros(max_tree_height*10,1);
            sum_la_omega=zeros(max_tree_height*10,1);
            sum_lad_=zeros(max_tree_height*10,1);
            sum_lad_omega=zeros(max_tree_height*10,1);
            count_=zeros(max_tree_height*10,1);
            
            % this next pass is there to create a vertical LAD profiles with the corrected LAD values for
            % occlusion
            %and alos add the omega values to the corresponding result file 
            while ~feof(fid)
                xmin=fscanf(fid, '%f', [1]);
                if isempty(xmin), break, end
                ymin=fscanf(fid, '%f', [1]);
                zmin=fscanf(fid, '%f', [1]);
                xmin_=round(xmin/ground_tile_size)+1; %these are used to find the minimum z value for each x,y column
                ymin_=round(ymin/ground_tile_size)+1;
                %zmin_=round(zmin*10)+1;
                zmin_=round((zmin-ground_level(xmin_,ymin_))*10)+1;
                
                xmax=fscanf(fid, '%f', [1]);
                ymax=fscanf(fid, '%f', [1]);
                zmax=fscanf(fid, '%f', [1]);
                lad=fscanf(fid, '%f', [1]);
                %omega_v=fscanf(fid, '%f', [1]);
                g_code=fscanf(fid, '%d', [1]);
                if zmin<ground_level(xmin_,ymin_)+voxel_size*number_of_voxels_above_ground_to_exclude, continue, end %voxel is at or below ground level
                la=lad*vox_size^3;
                %compute omega_v
                cumulative_LAI=sum(LA_profile_for_LAI(zmin_:size(LA_profile_for_LAI)),'omitnan')/plot_size^2;
                if cumulative_LAI<plateau
                    omega_v=(-log((cumulative_LAI-plateau)/(y0-plateau))/k)+x0;
                else
                    omega_v=1;
                end
                if lad>0, fprintf(fidres, '%.3f %.3f %.3f %.3f %.3f %.3f %f %f %d\n',xmin,ymin,zmin,xmax,ymax,zmax,lad,omega_v,g_code);end
                lad_omega=lad/omega_v;
                la_omega=lad_omega*vox_size^3;
                
                sum_la_(zmin_)=sum_la_(zmin_)+la;
                sum_la_omega(zmin_)=sum_la_omega(zmin_)+la_omega;
                sum_lad_(zmin_)=sum_lad_(zmin_)+lad;
                sum_lad_omega(zmin_)=sum_lad_omega(zmin_)+lad_omega;
                
                count_(zmin_)=count_(zmin_)+1;
            end
            
            
            profile_=sum_lad_./count_;  %this _ symbol refers to the corrected values for occlusion
            profile_omega=sum_lad_omega./count_; %this is corrected for within voxel clumping
            total_leaf_area=sum(sum_la_)
            total_leaf_area_omega_v=sum(sum_la_omega)
            LAI_occ_corrected=total_leaf_area/plot_size^2
            LAI_occ_corrected_omega=total_leaf_area_omega_v/plot_size^2
            
            fprintf(fidres_stats, 'Poisson correction factor used: %.2f \n', poisson);
            fprintf(fidres_stats, 'Occlusion threshold used: %.1f meters\n', occlusion_threshold);
            fprintf(fidres_stats, 'Proportion of occluded voxels: %.2f\n', sum(occluded)/sum(count_));
            %fprintf(fidres_stats, 'Layer thickness used for within voxel clumping: %.1f meters\n', layer_thickness_for_voxel_clumping);
            %fprintf(fidres_stats, 'Within voxel clumping factor used: %.2f\n',within_voxel_clumping_factor );
            fprintf(fidres_stats, 'Total leaf area before occlusion correction: %.0f m2\nTotal leaf area after occlusion correction: %.0f m2\nLAI with occlusion corrected: %.3f m2/m2\n', total_leaf_area_before_correction, total_leaf_area, LAI_occ_corrected);
            fprintf(fidres_stats, 'Total leaf area after considering within voxel clumping: %.0f m2\nLAI considering within voxel clumping : %.3f m2/m2\n',total_leaf_area_omega_v, LAI_occ_corrected_omega);
            fprintf(fidres_stats, 'Number of voxels having LAD above limit: %.d \n', size(lad_above_max,1));
            fprintf(fidres_stats, 'height (m)/leaf area per layer/leaf area per layer (considering within voxel clumping)/Average leaf area density per voxel per layer/Average leaf area density per voxel per layer(considering within voxel clumping)/number of occluded voxels/number of voxels containing material/proportion of occluded voxels per layer\n');
            for h=1:max_tree_height*10
                if sum_la_(h)>0,fprintf(fidres_stats, '%.2f %.1f %.1f %.3f %.3f %d %d %.4f \n',(h-1)/10+voxel_size/2,sum_la_(h),sum_la_omega(h),profile_(h),profile_omega(h),occluded(h),count_(h),occluded(h)/count_(h));end
            end
            
            
            
        end
        fclose all;
    end
end
beep