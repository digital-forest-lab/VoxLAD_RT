function [flag ,tmin] = rayBoxIntersection(origin, direction, vmin, vmax)
%  Ray/box intersection using the Smits' algorithm
%
% Input:
%    origin.
%    direction.
%    box = (vmin,vmax)
% Output:
%    flag: (0) Reject, (1) Intersect.
%    tmin: distance from the ray origin.
% Author: 
%    Jesus Mena

    if (direction(1) >= 0) 
    	tmin = (vmin(:,1) - origin(1)) / direction(1);
    	tmax = (vmax(:,1) - origin(1)) / direction(1);
    else
    	tmin = (vmax(:,1) - origin(1)) / direction(1);
    	tmax = (vmin(:,1) - origin(1)) / direction(1);
    end
  
    if (direction(2) >= 0) 
        tymin = (vmin(:,2) - origin(2)) / direction(2);
        tymax = (vmax(:,2) - origin(2)) / direction(2);
    else
    	tymin = (vmax(:,2) - origin(2)) / direction(2);
    	tymax = (vmin(:,2) - origin(2)) / direction(2);
    end

%     if ( (tmin > tymax) || (tymin > tmax) )
%         flag = 0;
%         tmin = -1;
%     	return;
%     end
    flag((tmin > tymax) | (tymin > tmax))=0;
    tmin((tmin > tymax) | (tymin > tmax))=0;
       
%     if (tymin > tmin)
%         tmin = tymin;
%     end
    tmin(tymin > tmin)=tymin(tymin > tmin);
    
% 	if (tymax < tmax)
%         tmax = tymax;
%     end
    tmax(tymax < tmax)=tymax(tymax < tmax);
    
	if (direction(3) >= 0)
       tzmin = (vmin(:,3) - origin(3)) / direction(3);
       tzmax = (vmax(:,3) - origin(3)) / direction(3);
    else
       tzmin = (vmax(:,3) - origin(3)) / direction(3);
       tzmax = (vmin(:,3) - origin(3)) / direction(3);
    end


%     if ((tmin > tzmax) || (tzmin > tmax))
%         flag = 0;
%         tmin = -1;
%        return;
%     end
     flag((tmin > tzmax) | (tzmin > tmax))=0;
    tmin((tmin > tzmax) | (tzmin > tmax))=0;
%     
%     if (tzmin > tmin)
%         tmin = tzmin;
%     end
    tmin(tzmin > tmin)=tzmin(tzmin > tmin);
   
%     if (tzmax < tmax)
%         tmax = tzmax;
%     end
    tmax(tzmax < tmax)=tzmax(tzmax < tmax);
    
  % if( (tmin < t1) && (tmax > t0) )
      flag = 1;
  % else
  %    flag = 0;
  %    tmin = -1;
  % end;
end
