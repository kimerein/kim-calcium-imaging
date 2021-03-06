function comp=chooseComponents(C_df)

% Input from user to choose components
% Components refer to output of CNMF and should refer to the rows in C_df
% Default is all components, excluding last row of C_df, which is background

prompt='Which components to use in average? Enter their numbers.';
options.WindowStyle='normal';
c=inputdlg(prompt,'Choose components',1,{num2str(1:size(C_df,1)-1)},options);
try
    comp=str2num(c{1});
catch
    prompt='Please enter as a list of integers separated by spaces.';
    c=inputdlg(prompt,'Choose components',1,{num2str(1:size(C_df,1)-1)},options);
    try
        comp=str2num(c{1});
    catch
        comp=1:size(C_df,1)-1;
        disp('Improper input format. Using all components except background.');
    end
end     

if any(comp>size(C_df,1)) || any(comp<1)
    disp(['Components should be integers between 1 and ' num2str(size(C_df,1)) '.']);
    prompt=['Which components to use? Enter integers between 1 and ' num2str(size(C_df,1)) '.'];
    c=inputdlg(prompt,'Choose components',1,{num2str(1:size(C_df,1)-1)},options);
    try
        comp=str2num(c{1});
    catch
        comp=1:size(C_df,1)-1;
        disp('Improper input format. Using all components except background.');
    end
end

if any(comp>size(C_df,1)) || any(comp<1)
    comp=1:size(C_df,1)-1;
    disp('Get user-defined components FAILED. Using all components except background.');
end