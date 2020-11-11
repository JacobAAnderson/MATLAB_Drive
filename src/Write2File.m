


function Write2File(filePath, format, data)

fileID = fopen(filePath,'a+'); 
fprintf(fileID, format, data');
fclose(fileID);

end



