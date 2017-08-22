% extract ages from tables
% MATLAB 2012 doesn't read tables correctly...


load DH_ages.mat


% Diamond Hill -- near modern Darwin GL
DH_samples = DH{:,1};
DH_elevations = DH{:,2};
DH_Be10_ages = DH{:,3};
DH_Be10_int_error = DH{:,4};
DH_Be10_ext_error = DH{:,5};
DH_Al26_ages = DH{:,6};
DH_Al26_int_error = DH{:,7};
DH_Al26_ext_error = DH{:,8};

% Danum Platform -- 45 km along Hatherton
DAN_samples = DAN{:,1};
DAN_elevations = DAN{:,2};
DAN_Be10_ages = DAN{:,3};
DAN_Be10_int_error = DAN{:,4};
DAN_Be10_ext_error = DAN{:,5};
DAN_Al26_ages = DAN{:,6};
DAN_Al26_int_error = DAN{:,7};
DAN_Al26_ext_error = DAN{:,8};

% Magnis Valley -- 30 km along Hatherton
MVfloor_samples = MVfloor{:,1};
MVfloor_elevations = MVfloor{:,2};
MVfloor_Be10_ages = MVfloor{:,3};
MVfloor_Be10_int_error = MVfloor{:,4};
MVfloor_Be10_ext_error = MVfloor{:,5};
MVfloor_Al26_ages = MVfloor{:,6};
MVfloor_Al26_int_error = MVfloor{:,7};
MVfloor_Al26_ext_error = MVfloor{:,8};

% Lake Wellman C14 -- 12 km along Hatherton
LWC14_elevations = LWC14{:,1};
LWC14_ages = LWC14{:,2};
LWC14__error = LWC14{:,3};



% Prune out LGM values...
% Diamond Hill: no LGM limit
% Danum Platform: 1500 m, 7.5-8.6 kyr BP
% Magnis Valley: 1330 m, 8.9-9.2 kyr BP
% Lake Wellman: 1120 m, 11-8.3 kyr BP   

disp('Is position given value, or +10km ?? ')
disp(' ' )

DAN_index_use = find((DAN_Be10_ages >= 7500) & (DAN_Be10_ages <= 8600));
DAN_ages_use = DAN_Be10_ages(DAN_index_use);
DAN_elevations_use = DAN_elevations(DAN_index_use);
DAN_position = 45000;

MVfloor_index_use = find((MVfloor_Be10_ages >= 8900) & (MVfloor_Be10_ages <= 9200));
MVfloor_ages_use = MVfloor_Be10_ages(MVfloor_index_use);
MVfloor_elevations_use = MVfloor_elevations(MVfloor_index_use);
MVfloor_position = 30000;

LWC14_index_use = find((LWC14_ages >= 8300) & (LWC14_ages <= 10000));
LWC14_ages_use = LWC14_ages(LWC14_index_use);
LWC14_elevations_use = LWC14_elevations(LWC14_index_use);
LWC14_position = 12000;

save LGM_values_use.mat DAN_ages_use DAN_elevations_use DAN_position ...
                        MVfloor_ages_use MVfloor_elevations_use MVfloor_position ...
                        LWC14_ages_use LWC14_elevations_use LWC14_position
                    
                    
                    


