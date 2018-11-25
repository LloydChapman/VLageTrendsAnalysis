function extract_raw_data(data,AgeLow,AgeHigh,SeroTests,SeroPoss)

%% Matrices of bounds for age groups and numbers tested and positive
AgeLowMat=table2array(data(:,AgeLow));
AgeHighMat=table2array(data(:,AgeHigh));
TestMat=table2array(data(:,SeroTests));
PosMat=table2array(data(:,SeroPoss));

%% Construct matrix of results for each study
for i=1:size(data,1)
    m=find(~isnan(PosMat(i,:)),1,'last');
    l=find(~isnan(TestMat(i,:)),1,'last');
    if ~isempty(m)&&~isempty(l)
        x=NaN(4,m);
        x(1,:)=AgeLowMat(i,1:m);
        x(2,:)=AgeHighMat(i,1:m);
        x(3,:)=TestMat(i,1:m);
        x(4,:)=PosMat(i,1:m);
        save(['data' num2str(i) '.mat'],'x')
    end
end