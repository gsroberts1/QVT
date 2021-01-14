function v = load_dat(name, res)
%LOAD_DAT: Loads dat files in current directory.
%   
%   Used by: loadpcvipr.m
%   Dependencies: NONE

[fid,errmsg]= fopen(name,'r');
if fid < 0  % If name does not exist in directory
    disp(['Error Opening Data : ',errmsg]);
end

% Reads in as short, reshapes by image resolution (i.e. 320x320x320)
v = reshape(fread(fid,'short=>single'),res);
fclose(fid);

return