function [ x_nodes, t_nodes ] = load_nodes 


%--------------------------------------------------------------
%
% Michelle Koutnik (mkoutnik@ess.washington.edu)
% last updated: January 2009
%
% Define space and time nodes for values on calculation grid
% for accumulation rate, b_dot(x,t)
%
%---------------------------------------------------------------


global min_search_E min_search_fs min_search_bed
global lower_resolution


%   x- positions of nodes (meters)
%   ------------------------------
    
     x_nodes = (0:500:1.42165e+05);
 
   
if ( (lower_resolution == 1) || (min_search_E == 1) || (min_search_bed == 1) || (min_search_fs == 1) )
      
     x_nodes = (0:1000:1.42165e+05);

end
     
     
%   t_nodes must include the present day (t=0), as the last value!
%   t_nodes needs to be a column vector.
%   ------------------------------------
    t_nodes = [-20000:100:0]'; 
  

    
%   % For min search:  
%        t_nodes = [-20000:100:0]'; 
%        x_nodes = (17000:5000:300000);
%   
%  
 
     