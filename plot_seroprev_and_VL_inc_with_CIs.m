function plot_seroprev_and_VL_inc_with_CIs(data,AgeLow,AgeHigh,DATidx,rK39idx,LSTidx,PCRidx,Mltidx,SeroPosPrevs,CIL,CIU,VLidx,VLIncs,VLCIL,VLCIU)
%PLOT_SEROPREV_AND_VL_INC_WITH_CIS Plot infection prevalences and VL
%incidences against age

%% PREVALENCE

% Matrices of bounds for age groups
AgeLowMat=table2array(data(:,AgeLow))';
AgeHighMat=table2array(data(:,AgeHigh))';

% Matrix of mid points of age groups
AgeMidptMat=(AgeLowMat+AgeHighMat+1)/2;

% Values for horizontal bars showing age groups
AgeL=AgeMidptMat-AgeLowMat;
AgeH=AgeHighMat+1-AgeMidptMat;

% Matrix of seropositive prevalences
SeroPosPerMat=table2array(data(:,SeroPosPrevs))';

% Values for error bars for prevalences
errL=SeroPosPerMat-table2array(data(:,CIL))';
errU=table2array(data(:,CIU))'-SeroPosPerMat;

% Set plot colours for unique studies 
Cmbidx=[DATidx;rK39idx;PCRidx;LSTidx;VLidx];
[AuthStrtDt,~,idx]=unique(data(Cmbidx,{'Author','StartDate'}),'rows');
clrs=num2cell(distinguishable_colors(size(AuthStrtDt,1)),2);
clrs{strcmp(AuthStrtDt.Author,'Bern et al, 2006')}=clrs{strcmp(AuthStrtDt.Author,'Bern et al, 2007')&strcmp(AuthStrtDt.StartDate,'2002')};
clrs{strcmp(AuthStrtDt.Author,'Hasker et al, 2013')&strcmp(AuthStrtDt.StartDate,'01-Mar-2007')}=clrs{strcmp(AuthStrtDt.Author,'Hasker et al, 2013')&strcmp(AuthStrtDt.StartDate,'2009')};
m={'+','o','*','.','x','s','d','v','^','<','>','p','h'}';
mrkrs=repmat(m,ceil(size(AuthStrtDt,1)/numel(m)),1);
mrkrs{strcmp(AuthStrtDt.Author,'Bern et al, 2006')}=mrkrs{strcmp(AuthStrtDt.Author,'Bern et al, 2007')&strcmp(AuthStrtDt.StartDate,'2002')};
mrkrs{strcmp(AuthStrtDt.Author,'Hasker et al, 2013')&strcmp(AuthStrtDt.StartDate,'01-Mar-2007')}=mrkrs{strcmp(AuthStrtDt.Author,'Hasker et al, 2013')&strcmp(AuthStrtDt.StartDate,'2009')};
l={'-.','-','--'}';
lstys=repmat(l,ceil(size(AuthStrtDt,1)/numel(l)),1);
lstys{strcmp(AuthStrtDt.Author,'Bern et al, 2006')}=lstys{strcmp(AuthStrtDt.Author,'Bern et al, 2007')&strcmp(AuthStrtDt.StartDate,'2002')};
lstys{strcmp(AuthStrtDt.Author,'Hasker et al, 2013')&strcmp(AuthStrtDt.StartDate,'01-Mar-2007')}=lstys{strcmp(AuthStrtDt.Author,'Hasker et al, 2013')&strcmp(AuthStrtDt.StartDate,'2009')};

% Age groups to plot
AgeGps=1:8;

InclActvExclPastVL=repmat({''},size(data,1),1);
ActvVL=strcmp(data.SeroTestInclClinVL,'Y');
PastVL=strcmp(data.SeroTestInclPastVL,'Y');
MxdTest=strcmp(data.Author,'Kaushal et al, 2017');
InclActvExclPastVL(ActvVL&PastVL)={'*'};
InclActvExclPastVL(~ActvVL&~PastVL)={'**'};
InclActvExclPastVL(ActvVL&~PastVL)={'*,**'};
InclActvExclPastVL(MxdTest)={'***'};

% Plot prevalence by age group for each diagnostic
% DAT
figure; h=errorbar(AgeMidptMat(AgeGps,DATidx),SeroPosPerMat(AgeGps,DATidx),errL(AgeGps,DATidx),errU(AgeGps,DATidx),'LineWidth',2);%,AgeL(AgeGps,DATidx),AgeH(AgeGps,DATidx))
% % Version w/o connecting lines
% errorbar(AgeMidptMat(AgeGps,DATidx),SeroPosPerMat(AgeGps,DATidx),errL(AgeGps,DATidx),errU(AgeGps,DATidx),AgeL(AgeGps,DATidx),AgeH(AgeGps,DATidx),'o')%'LineWidth',2)%
set(gca,'FontSize',14)
i=idx(ismember(Cmbidx,DATidx));
set(h,{'Color'},clrs(i),{'LineStyle'},lstys(i),{'Marker'},mrkrs(i),'MarkerSize',10)
xlabel('Age (yrs)'); ylabel('Proportion DAT positive')
l=legend(strcat(data.Author(DATidx),{' ('},data.ASMDefDAT(DATidx),{')'},{' ['},data.StartDate(DATidx),{']'},InclActvExclPastVL(DATidx)),'Location','northwest');
set(l,'FontSize',10)
ylim([0 0.4])
saveas(gcf,'DAT.eps','epsc')

% rK39
data.rK39test=repmat({''},size(data,1),1);
ELISAidx=strcmp(data.rK39ELISA,'Y');
RDTidx=strcmp(data.rK39RDT,'Y');
data.rK39test(RDTidx)={'RDT'};
data.rK39test(ELISAidx)={'ELISA'};
figure; h1=errorbar(AgeMidptMat(AgeGps,rK39idx),SeroPosPerMat(AgeGps,rK39idx),errL(AgeGps,rK39idx),errU(AgeGps,rK39idx),'Linewidth',2);%,AgeL(AgeGps,rK39idx),AgeH(AgeGps,rK39idx))
set(gca,'FontSize',14)
i=idx(ismember(Cmbidx,rK39idx));
set(h1,{'Color'},clrs(i),{'LineStyle'},lstys(i),{'Marker'},mrkrs(i),'MarkerSize',10)
xlabel('Age (yrs)'); ylabel('Proportion rK39 positive')
l1=legend(strcat(data.Author(rK39idx),{' ['},data.StartDate(rK39idx),{'] ('},data.rK39test(rK39idx),{')'},InclActvExclPastVL(rK39idx)),'Location','northwest');
set(l1,'FontSize',10)
xlim([0 90])
saveas(gcf,'rK39.eps','epsc')

% LST
figure; h2=errorbar(AgeMidptMat(AgeGps,LSTidx),SeroPosPerMat(AgeGps,LSTidx),errL(AgeGps,LSTidx),errU(AgeGps,LSTidx),'Linewidth',2); hold on
% figure; h2=errorbar(AgeMidptMat(AgeGps,LSTidx),SeroPosPerMat(AgeGps,LSTidx),errL(AgeGps,LSTidx),errU(AgeGps,LSTidx),AgeL(AgeGps,LSTidx),AgeH(AgeGps,LSTidx),'Linewidth',2); hold on
% a=0:100;
% plot(a,0.5752*(1-exp(-0.0496*a)),'k--','LineWidth',2)
set(gca,'FontSize',14)
i=idx(ismember(Cmbidx,LSTidx));
set(h2,{'Color'},clrs(i),{'LineStyle'},lstys(i),{'Marker'},mrkrs(i),'MarkerSize',10)
xlabel('Age (yrs)'); ylabel('Proportion LST positive')
l2=legend(strcat(data.Author(LSTidx),{' ['},data.StartDate(LSTidx),{']'},InclActvExclPastVL(LSTidx)),'Location','northeast');
set(l2,'FontSize',10)
% legend([strcat(data.Author(LSTidx),{' ['},data.StartDate(LSTidx),{']'},InclActvExclPastVL(LSTidx));'Catalytic model'],'Location','northeast')
xlim([0 90])
hold off
saveas(gcf,'LST.eps','epsc')

% PCR
figure; h3=errorbar(AgeMidptMat(AgeGps,PCRidx),SeroPosPerMat(AgeGps,PCRidx),errL(AgeGps,PCRidx),errU(AgeGps,PCRidx),'Linewidth',2);%,AgeL(AgeGps,PCRidx),AgeH(AgeGps,PCRidx))
i=idx(ismember(Cmbidx,PCRidx));
set(gca,'FontSize',14)
set(h3,{'Color'},clrs(i),{'LineStyle'},lstys(i),{'Marker'},mrkrs(i),'MarkerSize',10)
xlabel('Age (yrs)'); ylabel('Proportion PCR positive')
l3=legend(strcat(data.Author(PCRidx),{' ['},data.StartDate(PCRidx),{']'},InclActvExclPastVL(PCRidx)),'Location','northwest');
set(l3,'FontSize',10)
xlim([0 80])
saveas(gcf,'PCR.eps','epsc')

% SAME STUDY DIFFERENT TESTS
Type=cellfun(@(x)strrep(x,'PREVALENCE (',''),data.Type,'UniformOutput',false);
Type=cellfun(@(x)strrep(x,')',''),Type,'UniformOutput',false);
Type=cellfun(@(x)strrep(x,'(2002',''),Type,'UniformOutput',false);
for i=1:numel(Mltidx)
figure; errorbar(AgeMidptMat(AgeGps,Mltidx{i}),SeroPosPerMat(AgeGps,Mltidx{i}),errL(AgeGps,Mltidx{i}),errU(AgeGps,Mltidx{i}),'Linewidth',2);%,AgeL(AgeGps,PCRidx),AgeH(AgeGps,PCRidx))
set(gca,'FontSize',14)
xlabel('Age (yrs)'); ylabel('Proportion positive')
l4=legend(strcat(Type(Mltidx{i}),InclActvExclPastVL(Mltidx{i})),'Location','northwest');
set(l4,'FontSize',12)
xlim([0 90])
ylim([0 0.7])
saveas(gcf,['Mltidx' num2str(i) '.eps'],'epsc')
end

%% VL INCIDENCE

% Matrix of VL case incidences
VLCaseIncMat=table2array(data(:,VLIncs))';
% Values for error bars for incidences
VLerrL=VLCaseIncMat-table2array(data(:,VLCIL))';
VLerrU=table2array(data(:,VLCIU))'-VLCaseIncMat;

% Plot VL incidence by age group
figure; h5=errorbar(AgeMidptMat(AgeGps,VLidx),VLCaseIncMat(AgeGps,VLidx),VLerrL(AgeGps,VLidx),VLerrU(AgeGps,VLidx),'Linewidth',2);
set(gca,'FontSize',14)
i=idx(ismember(Cmbidx,VLidx));
set(h5,{'Color'},clrs(i),{'LineStyle'},lstys(i),{'Marker'},mrkrs(i),'MarkerSize',10)
xlabel('Age (yrs)'); ylabel('VL incidence (cases/1000/yr)')
l5=legend(strcat(data.Author(VLidx),{' ['},data.StartDate(VLidx),{' - '},data.EndDate(VLidx),{']'}));
set(l5,'FontSize',9)
saveas(gcf,'VLinc.eps','epsc')