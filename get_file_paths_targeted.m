function fpaths = get_file_paths_targeted(folderpath, varargin)
% returns files on 'folderpath' that contain all the strings in varargin

%unwrap input if double wrapped
if iscell(varargin{1})
    varargin = varargin{1};
end

% all paths under umbrella
fpaths = get_file_paths_all(folderpath);

% constrain via strings
for istr = 1:length(varargin)
    fpaths = fpaths(contains(fpaths, varargin{istr}));
end


