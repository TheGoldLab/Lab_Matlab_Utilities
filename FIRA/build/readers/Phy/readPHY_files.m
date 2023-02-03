function error_ = readPHY_files(filename)
% function error_ = readPHY_files(filename)
%
% Reads sorted spike data from a folder of Phy Numpy files. The
% passed-in filename should be the "params.py" file used by Phy, which sits
% in in the same folder as various .npy and .tsv files.
%
% This requires npy-matlab, to read the Numpy .npy files.
% https://github.com/kwikteam/npy-matlab
%
% This only looks for spike data, it doesn't look for any task or
% behavioral data.  The idea is to load an experiment from another data
% source beforehand, then load the spikes useing this function.
%
% This will label spikes by "cluster" id, as opposed to recording channel
% or unit ids.  It will assign all clusters to channel 0, and use each
% cluster id as the unit id.  As a result, spikes loaded from here should
% have ids in the range 0-999.  It might be useful to set FIRA.spm.spikes
% to a set of values in this range, so as to ignore spikes from other
% sources, like Plexon.
%
% putting the data into the global
%  FIRA data structure, as follows:
%
%    FIRA.header.filename ... filename
%    FIRA.header.filetype ... 'phy'
%    FIRA.header.paradigm ... 'xxx'
%
%   raw spike data is put into FIRA.raw.spikes.
%   spike channel ids are initialized into FIRA.spikes.
%
% Arguments:
%   filename  ... typically includes the full path
%
% Returns:
%   error_    ... duh
%   also fills appropriate fields of the global FIRA
%


% parse arguments
if nargin < 1 || isempty(filename)
    error_ = 1;
    return
end

% try locating the folder that contains "params.py"
if isfile(filename)
    phyDir = fileparts(filename);
elseif isfile([filename '.py'])
    phyDir = fileparts([filename '.py']);
else
    error_ = 2;
    return
end

%%%
% FILL IN FIRA HEADER INFO
%%%
global FIRA

error_ = 0;

if isempty(FIRA.header.filename)
    FIRA.header.filename = {filename};
else
    FIRA.header.filename{end+1} = filename;
end

if isempty(FIRA.header.filetype)
    FIRA.header.filetype = 'phy';
else
    FIRA.header.filetype = [FIRA.header.filetype, ', phy'];
end

FIRA.header.paradigm = 'xxx';


%% Read data from Phy.
phy.params = readParams(filename);

npyFiles = dir(fullfile(phyDir, '*.npy'));
for ii = 1:numel(npyFiles)
    npyFile = npyFiles(ii);
    data = readNPY(fullfile(npyFile.folder, npyFile.name));
    [~, fieldName] = fileparts(npyFile.name);
    phy.(fieldName) = data;
end

tsvFiles = dir(fullfile(phyDir, '*.tsv'));
for ii = 1:numel(tsvFiles)
    tsvFile = tsvFiles(ii);
    data = readtable(fullfile(tsvFile.folder, tsvFile.name), 'FileType', 'delimitedtext');
    [~, fieldName] = fileparts(tsvFile.name);
    phy.(fieldName) = data;
end


%% Convert Phy data into FIRA raw format.
keep_spikes = isfield(FIRA.raw, 'spikes');
if keep_spikes
    % Phy reports spikes in "clusters".
    % Assigning clusters as the "units" of "channel 0" results in FIRA
    % spike ids that are the same as the clusters, in the range 0-999.
    spikeClusterIds = double(phy.spike_clusters);
    uniqueClusterIds = unique(spikeClusterIds);
    id = [zeros(size(uniqueClusterIds)), uniqueClusterIds];

    % The verify() step filters spike ids in the data against a list passed
    % in by the caller.  It also sets up FIRA.spikes headers / bookkeeping.
    verify(FIRA.spm.spikes, id);

    % Phy spike times are in units of sample numbers.
    % Convert to milliseconds.
    spikeTimes = 1000 * double(phy.spike_times) ./ phy.params.sample_rate;
    clusterSpikeTimes = [spikeTimes, spikeClusterIds];
    FIRA.raw.spikes = cat(1, FIRA.raw.spikes, clusterSpikeTimes);
end


% Try to read params that were saved as Python code in params.py.
% This scrapes out simple variable assignments and saves the results in
% Matlab.  It's not a real Python evaluator!
function params = readParams(paramsFile)

params = struct();

lines = readlines(paramsFile, 'EmptyLineRule', 'skip');
for ii = 1:numel(lines)
    line = lines{ii};
    rawSplits = split(line, '=');
    key = strip(rawSplits{1});
    value = strip(rawSplits{2});
    switch value
        case 'True'
            params.(key) = true;
        case 'False'
            params.(key) = true;
        otherwise
            try
                params.(key) = eval(rawSplits{2});
            catch e
                warning('Error parsing value <%s> for key <%s>: %s', ...
                    value, key, e.message);
            end
    end
end
