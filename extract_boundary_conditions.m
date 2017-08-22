% extract boundary conditions needed from different minimization runs


for ii = 1:8
    
eval(['load run_min_search' int2str(ii) 'a.mat'])

x_P_min = x_P; x_w_min = x_w; x_e_min = x_e;
E_P_min = E_P; E_w_min = E_w; E_e_min = E_e;
fs_P_min = fs_P; fs_w_min = fs_w; fs_e_min = fs_e;
x_P2_min = x_P2; x_w2_min = x_w2; x_e2_min = x_e2;
E_P2_min = E_P2; E_w2_min = E_w2; E_e2_min = E_e2;
fs_P2_min = fs_P2; fs_w2_min = fs_w2; fs_e2_min = fs_e2;


eval(['save min' int2str(ii) '.mat x_P_min x_w_min x_e_min x_P2_min x_w2_min x_e2_min E_P_min E_w_min E_e_min fs_P_min fs_w_min fs_e_min E_P2_min E_w2_min E_e2_min fs_P2_min fs_w2_min fs_e2_min deformation_only deformation_plus_sliding sliding_only'])

           
end



% create a test case 9 that uses ... 
ii = 9;

load run_min_search4a.mat
  x_P_min = x_P; x_w_min = x_w; x_e_min = x_e;
  E_P_min = E_P; E_w_min = E_w; E_e_min = E_e;
  fs_P_min = fs_P; fs_w_min = fs_w; fs_e_min = fs_e;

load run_min_search6a.mat
  x_P2_min = x_P2; x_w2_min = x_w2; x_e2_min = x_e2;
  E_P2_min = E_P2; E_w2_min = E_w2; E_e2_min = E_e2;
  fs_P2_min = fs_P2; fs_w2_min = fs_w2; fs_e2_min = fs_e2;

eval(['save min' int2str(ii) '.mat x_P_min x_w_min x_e_min x_P2_min x_w2_min x_e2_min E_P_min E_w_min E_e_min fs_P_min fs_w_min fs_e_min E_P2_min E_w2_min E_e2_min fs_P2_min fs_w2_min fs_e2_min deformation_only deformation_plus_sliding sliding_only'])
  