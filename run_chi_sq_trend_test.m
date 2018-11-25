function data=run_chi_sq_trend_test(data,VLidx,AgeLow,AgeHigh,StudyPops,VLCases)
%RUN_CHI_SQ_TREND_TEST Run chi-square trend (Cochrane-Armitage) test
data.chisqp=NaN(size(data,1),1);

for i=1:numel(VLidx)
   j=VLidx(i);
   y=data(j,[AgeLow,AgeHigh,StudyPops,VLCases]);
   y=reshape(table2array(y),numel(AgeLow),4);
   y(isnan(y(:,1))|y(:,1)<14,:)=[];
   [p,stats]=cochran_arm(y(:,3:4));
   data.chisqp(j)=p;
end



    