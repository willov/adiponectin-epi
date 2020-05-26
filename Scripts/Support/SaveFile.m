function [] = SaveFile(fileName,variable, variableName)
%SAVEFILE Saves a file with "fileName" with the content of a "variable" with
%name "variableName"

S.(variableName)=variable;
save(fileName,'-struct', 'S')

end

