function print_table(data,AgeLow,AgeHigh,SeroTests,SeroPoss,SeroPosPrevs,CIL,CIU,OR,ORCIL,ORCIU,p,StudyPops,VLCases,VLIncs,VLCIL,VLCIU,VLRR,VLRRCIL,VLRRCIU,VLRRp)
%PRINT_TABLE Print table of data from each study to csv file

% Remove previously saved output
delete *.csv

for i=1:size(data,1)
    % Get data from ith study
    datai=data(i,:);
    fid=fopen([data.Author{i},'_' num2str(i) '.csv'],'a');
    % Stack data into long table format
    if (strcmp(data.Type(i),'INCIDENCE')||strcmp(data.Type(i),'INCIDENCE (RAW DATA)'))
        S=stack(datai,{AgeLow,AgeHigh,StudyPops,VLCases,VLIncs,VLCIL,VLCIU,VLRR,VLRRCIL,VLRRCIU,VLRRp},'NewDataVariableName',{'AgeLow','AgeHigh','StudyPop','NumVLCases','VLInc','CIL','CIU','RR','RRCIL','RRCIU','p'},'ConstantVariables',[]);
        fprintf(fid,[repmat('%s,',1,6) '\n'],'Age group (yrs)','n','No. cases','Incidence/1000/yr (95% CI)','RR (95% CI)','p');
        fmt='%.2f (%.2f-%.2f),';
    else
        S=stack(datai,{AgeLow,AgeHigh,SeroTests,SeroPoss,SeroPosPrevs,CIL,CIU,OR,ORCIL,ORCIU,p},'NewDataVariableName',{'AgeLow','AgeHigh','NumSeroTest','NumSeroPos','SeroPosPrev','CIL','CIU','OR','ORCIL','ORCIU','p'},'ConstantVariables',[]);
        fprintf(fid,[repmat('%s,',1,6) '\n'],'Age group (yrs)','n','No. positive','Prevalence (95% CI)','OR (95% CI)','p');
        fmt='%.3f (%.3f-%.3f),';
    end
    % Remove indicator variable from table
    S.Indicator=[];
    % Remove empty rows
    S=S(~all(isnan(table2array(S)),2),:);
    % Write output to file with name 'Author, YYYY.csv'
    str=sprintf(['%.0f-%.0f,' repmat('%.0f,',1,2) fmt '%.2f (%.2f-%.2f),%.3f,\n'], table2array(S)');
    str=strrep(str,'NaN','');
    str=strrep(str,' (-)','Ref.');
    fprintf(fid,'%s',str);
    fclose(fid);
end