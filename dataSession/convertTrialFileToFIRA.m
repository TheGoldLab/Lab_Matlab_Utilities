function dataInFIRAFormat = convertTrialFileToFIRA(dataInTrialFileFormat, filename)
% function dataInFIRAFormat = convertTrialFileToFIRA(dataInTrialFileFormat, filename)
%
% Convert trialFile to FIRA
%

if nargin < 2 || isempty('filename')
    filename = 'Converted from trialFile';
end

% Get numTrials -- remember first and last are dummies
numTrials = size(dataInTrialFileFormat, 2)-2;

%% Set up FIRA
%
dataInFIRAFormat = struct( ...
    'header', struct(...
    'filename', filename, ...
    'filetype', 'Converted from trialFile', ...
    'paradigm', 'xxx', ...
    'subject', [], ...
    'session', [], ...
    'spmF', [], ...
    'numTrials', numTrials, ...
    'date', date, ...
    'messageH', [], ...
    'flags', []), ...
    'spm', 'Converted from trialFile', ...
    'ecodes', struct( ...
    'name', {{}}, ...
    'type', {{}}, ...
    'data', []), ...
    'spikes', struct( ...
    'channel', [], ...
    'unit', [], ...
    'id', [], ...
    'data', {{}}), ...
    'analog', struct( ...
    'name', {{}}, ...
    'acquire_rate', [], ...
    'store_rate', [], ...
    'gain', [], ...
    'offset', [], ...
    'zero', {{}}, ...
    'error', {{}}, ...
    'data', []), ...
    'dio', {{}});

% Return if no data given
if numTrials < 1
    return
end

%% Add ecodes
%
% Start with four "standard" fields
dataInFIRAFormat.ecodes.name = {'trial_num', 'trial_begin', 'trial_end', 'trial_wrt'};
dataInFIRAFormat.ecodes.type = {'value', 'time', 'time', 'time'};
dataInFIRAFormat.ecodes.data = nan(numTrials, 4);
dataInFIRAFormat.ecodes.data(:,1) = 1:numTrials;
dataInFIRAFormat.ecodes.data(:,2) = [dataInTrialFileFormat(2:end-1).start_time]';
dataInFIRAFormat.ecodes.data(:,3) = [dataInTrialFileFormat(2:end-1).end_time]';
dataInFIRAFormat.ecodes.data(:,4) = [dataInTrialFileFormat(2:end-1).wrt_time]';

% Now add all the "enhancements", by category.
% Loop through the trials
for tt = 1:numTrials

    % Loop through the enhancement categories present in each trial
    categories = fieldnames(dataInTrialFileFormat(tt+1).enhancement_categories);
    for cc = 1:length(categories)

        % Get list of enhancements of this type
        enhancements = dataInTrialFileFormat(tt+1).enhancement_categories.(categories{cc});
        numEnhancements = length(enhancements);
            
        % Loop through the enhancements
        for ee = 1:numEnhancements
            value = dataInTrialFileFormat(tt+1).enhancements.(enhancements{ee});
            if isstruct(value)

                % Loop through 
                doubleSecretEnhancements = fieldnames(value);
                numDoubleSecretEnhancements = length(doubleSecretEnhancements);

                for dd = 1:numDoubleSecretEnhancements

                    % Save the name/type (always "value", just cuz)
                    name = [enhancements{ee} '_' doubleSecretEnhancements{dd}];
                    Lcolumn = strcmp(name, dataInFIRAFormat.ecodes.name);
                    if ~any(Lcolumn)
                        dataInFIRAFormat.ecodes.name(end+1) = {name};
                        dataInFIRAFormat.ecodes.type(end+1) = {'value'};
                        dataInFIRAFormat.ecodes.data(:,end+1) = nan;
                        Lcolumn(end+1) = true;
                    end

                    % Save the data
                    value = dataInTrialFileFormat(tt+1).enhancements.(enhancements{ee}).(doubleSecretEnhancements{dd});
                    if isempty(value)
                        value = nan;
                    elseif length(value) > 1
                        fprintf('Warning: found >1 values for %s on trial %d\n', ...
                            name, tt+1)
                        value = value(1);
                    end
                    dataInFIRAFormat.ecodes.data(tt,Lcolumn) = value;
                end

            else

                % Save the name/type (always "value", just cuz)
                Lcolumn = strcmp(enhancements{ee}, dataInFIRAFormat.ecodes.name);
                if ~any(Lcolumn)
                    dataInFIRAFormat.ecodes.name(end+1) = enhancements(ee);
                    dataInFIRAFormat.ecodes.type(end+1) = categories(cc);
                    dataInFIRAFormat.ecodes.data(:,end+1) = nan;
                    Lcolumn(end+1) = true;
                end

                % Save the data
                if isempty(value)
                    value = nan;
                elseif length(value) > 1
                    fprintf('Warning: found >1 values for %s on trial %d\n', ...
                        enhancements{ee}, tt+1)
                    value = value(1);
                end
                dataInFIRAFormat.ecodes.data(tt,Lcolumn) = value;
            end
        end
    end
end

% convention is names are row vector
dataInFIRAFormat.ecodes.name = dataInFIRAFormat.ecodes.name(:)';

%% Add analog signals
%
signalNames = fieldnames(dataInTrialFileFormat(2).signals);
numSignals = length(signalNames);
dataInFIRAFormat.analog.data = struct( ...
    'start_time', cell(numTrials, numSignals), ...
    'length', [], ...
    'values', []);

% Loop through each signal type
for ss = 1:numSignals

    % Add name, sample rate
    dataInFIRAFormat.analog.name = ...
        cat(2, dataInFIRAFormat.analog.name, signalNames{ss});
    dataInFIRAFormat.analog.acquire_rate = ...
        cat(2, dataInFIRAFormat.analog.acquire_rate, ...
        dataInTrialFileFormat(2).signals.(signalNames{ss}).sample_frequency);
    dataInFIRAFormat.analog.store_rate = dataInFIRAFormat.analog.acquire_rate;

    % Add data
    for tt = 1:numTrials
        dataInFIRAFormat.analog.data(tt,ss).start_time = ...
            dataInTrialFileFormat(tt+1).signals.(signalNames{ss}).first_sample_time;
        dataInFIRAFormat.analog.data(tt,ss).values = ...
            dataInTrialFileFormat(tt+1).signals.(signalNames{ss}).signal_data;
        dataInFIRAFormat.analog.data(tt,ss).length = ...
            length(dataInFIRAFormat.analog.data(tt,ss).values);
    end
end

%% Add spikes
%
spikeChannelNames = setdiff(fieldnames(dataInTrialFileFormat(2).numeric_events), {'ecodes'});
numSpikeChannels = length(spikeChannelNames);

% Start the data celery
dataInFIRAFormat.spikes.data = {};

% Loop through the channels to count units
for ss = 2:numSpikeChannels

    % Parse channel number
    channelNumber = str2double(spikeChannelNames{ss}(end-2:end));

    % Keep track of units
    channelUnits = [];

    % Keep track of where we're storing the data
    channelRowStart = 0;

    % Loops through the trials and collect the data
    for tt = 1:numTrials

        % Get the data
        trialSpikes = dataInTrialFileFormat(tt+1).numeric_events.(spikeChannelNames{ss});

        if ~isempty(trialSpikes)

            % Get the units
            units = unique(trialSpikes(:,3));
            numUnits = length(units);

            % Loop through the units
            for uu = 1:numUnits

                % Check if we need to add the unit
                if ~any(channelUnits == units(uu))

                    % Add to lists of channels/units/ids
                    dataInFIRAFormat.spikes.channel = cat(2, dataInFIRAFormat.spikes.channel, channelNumber);
                    dataInFIRAFormat.spikes.unit = cat(2, dataInFIRAFormat.spikes.unit, units(uu));
                    dataInFIRAFormat.spikes.id = cat(2, dataInFIRAFormat.spikes.id, channelNumber*1000+units(uu));

                    % Add to list
                    channelUnits = cat(2, channelUnits, units(uu));

                    % Add data column
                    dataInFIRAFormat.spikes.data = ...
                        cat(2, dataInFIRAFormat.spikes.data, cell(numTrials, 1));
                end

                % Now add the spike data
                Lunit = trialSpikes(:,3) == units(uu);
                dataInFIRAFormat.spikes.data{tt,channelRowStart+uu} = trialSpikes(Lunit,1);
            end
        end
    end

    % Update the row counter
    channelRowStart = size(dataInFIRAFormat.spikes.data, 2);
end
