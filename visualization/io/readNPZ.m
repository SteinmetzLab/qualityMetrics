function s = readNPZ(filename)
% readNPZ  Read a NumPy .npz archive into a struct (one field per array).
%
%   s = readNPZ(filename)
%
% A .npz file is just a (usually uncompressed) zip of .npy members. This
% extracts them to a temp folder and reads each with readNPY, returning a
% struct whose field names are the member names (e.g. s.amp_uV, s.peak_depth).
% Requires readNPY on the path.

if exist('readNPY', 'file') ~= 2
    error('readNPZ:dep', 'readNPY not found on the path.');
end

tmp = tempname;
mkdir(tmp);
cleaner = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>

names = unzip(filename, tmp);
s = struct();
for i = 1:numel(names)
    [~, nm, ext] = fileparts(names{i});
    if strcmpi(ext, '.npy')
        s.(matlab.lang.makeValidName(nm)) = readNPY(names{i});
    end
end
end
