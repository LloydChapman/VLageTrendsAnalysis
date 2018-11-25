function print_par_ests(x,fname)
    fid=fopen([fname '.txt'],'w');
    fprintf(fid,[repmat('%s & ',1,11) '%s \\\\ \n'],'{\bf Study}','{\bf Location}','{\bf Year}','{\bf Test}','$b_{0,i}$','{\bf 95\% CI}','$b_{1,i}$','{\bf 95\% CI}','$\gamma_i$','{\bf 95\% CI}','$-\log L$','AIC');
    y=table2cell(x);
    str=[];
    for i=1:size(y,1)
        str=[str sprintf([repmat('%s & ',1,2) '%.0f & %s & ' repmat('%.2g & (%.2g--%.2g) & ',1,3) '%.1f & %.1f \\\\ \n'],y{i,:})];
    end
    str=strrep(str,'NaN','');
    str=strrep(str,'(--)','');    
    fprintf(fid,'%s',str);
    fclose(fid);
end