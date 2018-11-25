function [med_lambda,med_gamma,age_pc_infctd,age_first_infctn,dur_imm_rspnse]=calc_smmry_stats(str)

data=readtable(str,'Format','%q %q %f %q %f %f %f %f %f %f %f %f %f %f %f');
idx=strcmp(data.Test,'DAT');
med_lambda=median(data.b0(idx));
med_gamma=median(data.gamma(idx));
a=0:90;
med_prev=med_lambda/(med_lambda+med_gamma)*(1-exp(-(med_lambda+med_gamma)*a));
figure; plot(a,med_prev);
xlabel('Age (yrs)'); ylabel('Proportion DAT positive')
pc_infctd=0.1:0.1:0.4;
age_pc_infctd=-log(1-pc_infctd)/med_lambda;
age_first_infctn=1/med_lambda;
dur_imm_rspnse=1/med_gamma;