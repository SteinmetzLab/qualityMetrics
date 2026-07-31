function data = readNPY(filename)
% readNPY  Minimal reader for NumPy .npy files (numeric arrays).
%
%   data = readNPY(filename)
%
% Handles the common case: .npy v1.0/2.0, numeric dtype, C- or Fortran-order.
% Returns the array with NumPy's shape preserved (row-major reshaped to MATLAB
% column-major correctly). Enough for the SpikeInterface property files
% (location, gain_to_uV, group, inter_sample_shift, offset_to_uV).

fid = fopen(filename, 'r', 'l');
if fid < 0, error('readNPY:open', 'Cannot open %s', filename); end
cleaner = onCleanup(@() fclose(fid));

% Read the magic as raw bytes — '*char' would locale-decode 0x93 and mangle it.
magic = fread(fid, 6, 'uint8')';
if ~isequal(magic, [147 double('NUMPY')])
    error('readNPY:magic', 'Not a .npy file: %s', filename);
end
major = fread(fid, 1, 'uint8');
fread(fid, 1, 'uint8');                       % minor
if major == 1
    hlen = fread(fid, 1, 'uint16');
else
    hlen = fread(fid, 1, 'uint32');
end
header = char(fread(fid, hlen, 'uint8')');    % ASCII header

% --- parse the header dict ---------------------------------------------
descr = regexp(header, '''descr''\s*:\s*''([^'']+)''', 'tokens', 'once');
descr = descr{1};
fortran = ~isempty(regexp(header, '''fortran_order''\s*:\s*True', 'once'));
shp = regexp(header, '''shape''\s*:\s*\(([^)]*)\)', 'tokens', 'once');
shape = str2num(['[' strrep(shp{1}, ',', ' ') ']']); %#ok<ST2NM>
if isempty(shape), shape = 1; end
if isscalar(shape), shape = [shape 1]; end

% --- map dtype ----------------------------------------------------------
endian = descr(1); key = descr(2:end);
map = struct('f8','double','f4','single','i8','int64','i4','int32', ...
             'i2','int16','i1','int8','u8','uint64','u4','uint32', ...
             'u2','uint16','u1','uint8','b1','uint8');
if ~isfield(map, key)
    error('readNPY:dtype', 'Unsupported dtype "%s" in %s', descr, filename);
end
cls = map.(key);

% reopen with correct endianness if big-endian
if endian == '>'
    pos = ftell(fid); clear cleaner; fclose(fid);
    fid = fopen(filename, 'r', 'b'); fseek(fid, pos, 'bof');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
end

raw = fread(fid, prod(shape), ['*' cls]);

% --- reshape honouring NumPy order --------------------------------------
if fortran
    data = reshape(raw, shape);
else
    data = reshape(raw, fliplr(shape));
    data = permute(data, numel(shape):-1:1);
end
end
