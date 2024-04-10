function save_sys_to_file(system)
    fileID = fopen('systemTemp.m', 'w'); 
    fclose(fileID);
    if fileID == -1 
        error('File is not opened'); 
    end
    n = size(system, 1);
    fileID = fopen('systemTemp.m', 'a'); 
    fileHeader = ['function [ dx ] = systemTemp(t, x)\n\tdx = zeros(' num2str(n) ', 1);\n'];
    fprintf(fileID, fileHeader);
    
    for i = 1:n
        strTmp = system(i, 1);
        for j = 1:n
            toReplace = ['x' num2str(n + 1 - j)];
            replacement = ['x(' num2str(n + 1 - j) ')'];
            strTmp = strrep(strTmp, toReplace, replacement); 
        end
        str = ['\tdx(' num2str(i) ') = ' char(strTmp) ';\n'];
        fprintf(fileID, str);
    end
    fprintf(fileID, 'end');
    fclose(fileID);
end
