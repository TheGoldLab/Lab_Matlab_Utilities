function dataInSessionDataFormat = convertFIRAToDataSession(dataInFIRAFormat, ...
    options)
% function dataInSessionDataFormat = convertFIRAToDataSession(dataInFIRAFormat, ...
%     options)

arguments
    dataInFIRAFormat
    options.trialSelectionArray = true(size(dataInFIRAFormat.ecodes.data,1),1);
    options.valueFields = 'all'; % continuous values
    options.idFields = 'all';    % discrete values
    options.timeFields = 'all';  % time values
    options.signalChannels = 'all';
    options.spikeIDs = 'all'; % 'all' or array of ids (dataInFIRAFormat.spikes.id)
end

%% Add header
%
dataInSessionDataFormat.header = dataInFIRAFormat.header;

%% Collect data from FIRA and put into a struct
%
% Use selectionArray to choose trials
numTrials = sum(options.trialSelectionArray);

%% Collect ecode values, IDs, and times into separate tables
%
fields = {...
    'values',   options.valueFields; ...
    'ids',      options.idFields; ...
    'times',    options.timeFields};
for ii = 1:size(fields, 1)

    % Check for entries
    if ~isempty(fields{ii,2})

        % Set up the table
        dataInSessionDataFormat.(fields{ii,1}) = table;

        if ischar(fields{ii,2})
            if strcmp(fields{ii,2}, 'all')
                Ltype = strcmp(fields{ii,1}(1:end-1), dataInFIRAFormat.ecodes.type);
                names = dataInFIRAFormat.ecodes.name(Ltype);
            else
                names = {fields{ii,2}};
            end
        else % if iscell(fields{ii,2})
            names = fields{ii,2};
        end

        % Add the data, per column
        for vv = 1:length(names)

            % Get the data
            inds = strcmp(names{vv}, dataInFIRAFormat.ecodes.name);
            if any(inds)
                ecodes = dataInFIRAFormat.ecodes.data(options.trialSelectionArray, inds);

                % Add it if anything found, otherwise just add nans
                if ~isempty(ecodes)
                    dataInSessionDataFormat.(fields{ii,1}).(names{vv}) = ecodes;
                else
                    dataInSessionDataFormat.(fields{ii,1}).(names{vv}) = nan(numTrials, 1);
                end
            end
        end
    end
end

%% Conditionally add spike data
%
dataInSessionDataFormat.spikes = [];
if ~isempty(options.spikeIDs) && isfield(dataInFIRAFormat, 'spikes')

    % check spikeIDs
    if ischar(options.spikeIDs) && strcmp(options.spikeIDs, 'all')
        options.spikeIDs = dataInFIRAFormat.spikes.id;
        Lspikes = true(size(options.spikeIDs));
    else % if isnumeric(spikeIDs)
        Lspikes = ismember(dataInFIRAFormat.spikes.id, options.spikeIDs);
    end

    if any(Lspikes)
        dataInSessionDataFormat.spikes = cell2table(...
            dataInFIRAFormat.spikes.data(options.trialSelectionArray, Lspikes), ...
            'VariableNames', cellstr(num2str(options.spikeIDs'))');
    end
end

%% Conditionally add signal (analog) data
% 
dataInSessionDataFormat.signals = [];
if ~isempty(options.signalChannels) && isfield(dataInFIRAFormat, 'analog')
    if strcmp(options.signalChannels, 'all')
        channelIndices = 1:length(dataInFIRAFormat.analog.name);
    else
        channelIndices = find(ismember(dataInFIRAFormat.analog.name, options.signalChannels));
    end
    trialIndices = find(options.trialSelectionArray);
    dataInSessionDataFormat.signals.sampleRates = ...
        dataInFIRAFormat.analog.store_rate(channelIndices);
    dataInSessionDataFormat.signals.data = cell(length(trialIndices), length(channelIndices));
    for tt = 1:length(trialIndices)
        for cc = 1:length(channelIndices)
            dataInSessionDataFormat.signals.data{tt, cc} = ...
                dataInFIRAFormat.analog.data(trialIndices(tt), channelIndices(cc)).values;
        end
    end
end