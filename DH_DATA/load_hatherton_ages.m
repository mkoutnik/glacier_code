% load_hatherton_ages


% -------------------------------------------------------------------------
% From Trevor 10 July 2017:
% Here are the LGM limit elevations, ages, and distances along the flowline:
% 
% Danum Platform: 1500 m, 7.5 - 8.6 kyr BP, 45 km along Hatherton flowline
% Magnis Valley: 1330 m, 8.9 - 9.2 kyr BP, 30 km along Hatherton flowline
% Lake Wellman: 1120 m, 10 - 8.3 kyr BP, 12 km along Hatherton flowline
% Diamond Hill: no LGM limit, start of Darwin flowline
% 
% Note that these are the elevations of the drift limit, not of the glacier 
% centerline. At Lake Wellman, for instance, the drift limit is almost 10 km 
% from the center of Hatherton Glacier. This means that these are lower limits for the  
% centerline surface elevation. The Danum Platform and Magnis Valley elevations are 
% from closer to the glacier, so those should be closer to the true centerline elevation, 
% but still a lower estimate. I'm not sure how we're going to handle that. 
% Currently, the center of Hatherton Glacier is 100 m higher than the edges, but 
% that would probably not hold for a wider glacier. Some detailed grappling with 
% the data at each location might help, but I won't be able to do that for a while!
% 
% I just added DH_ages.mat to the Git repository under DH_DATA. 
% It contains tables of data from each location. Here is a brief description of each table:
% 
% BH - The Brown Hills are on the far side of Diamond Hill from Darwin Glacier. 
%      Don't bother with these data at this time. Mostly pre-exposed, and may not mean much for the main trunk of the glacier.
% 
% Brit2 - The Britannia 2 deposit in Dubris and Bibra Valleys (top of Hatherton) 
%         are from MIS 6 and are ~137,000 years old. Don't worry about these for now either. 
%         45 km along Hatherton flowline
% 
% DAN - Danum Platform (between Dubris and Bibra Valleys at the head of Hatherton Glacier). 
%       I keep this separate from the ages in the valley because of the large elevation 
%       differences between the two, even though they deglaciated at the same time. 
%       45 km along Hatherton flowline
% 
% DH - Diamond Hill erratics. These will be used to directly force the elevation boundary condition. 
%      Start of Darwin flowline (glacier mouth).
% 
% DVBV - Dubris and Bibra Valleys. Same deglaciation history as DAN, but lower elevations. 
%        45 km along Hatherton flowline
% 
% MVfloor - From the floor of Magnis Valley on Hatherton Glacier. 
%           These ages are more consistent than the ages from the valley walls. 
%           30 km along Hatherton flowline.
% 
% MVwalls - see above. Some of these are likely pre-exposed. 30 km along Hatherton flowline.
% 
% LW - Erratics from Lake Wellman. Almost all pre-exposed, with the exception of 
%      13-HAT-010-LW, 13-HAT-018-LW, and 13-HAT-020-LW. For now, better to 
%      use the C-14 algae history. 12 km along Hatherton flowline.
% 
% LWC14 - Algae C-14 ages from Lake Wellman give a much better history than 
%         the exposure ages in LW. 12 km along Hatherton flowline.

% -------------------------------------------------------------------------


% load DH_ages.mat
load LGM_values_use.mat   % Run script through Matlab 2015...!





