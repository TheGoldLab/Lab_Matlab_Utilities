function dataInFIRAFormat = convertTrialFileToFIRA(dataInTrialFileFormat, filename, spikeCategories)
% function dataInFIRAFormat = convertTrialFileToFIRA(dataInTrialFileFormat)
%
% Convert trialFile to FIRA
%

if nargin < 2 || isempty(filename)
    filename = 'Converted from trialFile';
end

if nargin < 3 || isempty(spikeCategories)
    spikeCategories = {'spikes'};
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

%% Add ecodes from trial enhancements.
%
% Start with four "standard" fields
dataInFIRAFormat.ecodes.name = {'trial_num', 'trial_begin', 'trial_end', 'trial_wrt'};
dataInFIRAFormat.ecodes.type = {'value', 'time', 'time', 'time'};
dataInFIRAFormat.ecodes.data = nan(numTrials, 4);
dataInFIRAFormat.ecodes.data(:,1) = 1:numTrials;
dataInFIRAFormat.ecodes.data(:,2) = [dataInTrialFileFormat(2:end-1).start_time]';
dataInFIRAFormat.ecodes.data(:,3) = [dataInTrialFileFormat(2:end-1).end_time]';
dataInFIRAFormat.ecodes.data(:,4) = [dataInTrialFileFormat(2:end-1).wrt_time]';

% Use trial "categories" to indicate which trial fields,
% including events and enhancements, contain ecodes.
% Categories have names like "id", "value", and "time".
% Event and enhancement fields have names like "fp_on", and "t1_angle".

% Loop through the trials
for tt = 1:numTrials

    % Loop through the categories present in each trial.
    categories = dataInTrialFileFormat(tt+1).categories;
    if ~isempty(categories)

        % Treat categories that aren't for spikes as ecodes.
        ecodeCategoryNames = setdiff(fieldnames(categories), spikeCategories);
        for cc = 1:length(ecodeCategoryNames)

            % Get names of trial enhancements in this category.
            categoryName = ecodeCategoryNames{cc};
            fieldsInCategory = categories.(categoryName );
            enhancements = dataInTrialFileFormat(tt+1).enhancements;
            enhancementNames = {};
            if ~isempty(enhancements)
                enhancementNames = fieldnames(enhancements);
            end
            enhancementsInCategory = intersect(fieldsInCategory, enhancementNames);
            numEnhancements = length(enhancementsInCategory);

            % Loop through the enhancements in this category.
            for ee = 1:numEnhancements
                enhancementName = enhancementsInCategory{ee};
                value = dataInTrialFileFormat(tt+1).enhancements.(enhancementName);
                if isstruct(value)

                    % Loop through
                    doubleSecretEnhancements = fieldnames(value);
                    numDoubleSecretEnhancements = length(doubleSecretEnhancements);

                    for dd = 1:numDoubleSecretEnhancements

                        % Save the name/type (always "value", just cuz)
                        name = [enhancementName '_' doubleSecretEnhancements{dd}];
                        Lcolumn = strcmp(name, dataInFIRAFormat.ecodes.name);
                        if ~any(Lcolumn)
                            dataInFIRAFormat.ecodes.name(end+1) = {name};
                            dataInFIRAFormat.ecodes.type(end+1) = {'value'};
                            dataInFIRAFormat.ecodes.data(:,end+1) = nan;
                            Lcolumn(end+1) = true;
                        end

                        % Save the data
                        value = dataInTrialFileFormat(tt+1).enhancements.(enhancementName).(doubleSecretEnhancements{dd});
                        if isempty(value)
                            value = nan;
                        elseif length(value) > 1
                            fprintf('Warning: found >1 values for %s at trial file index %d, FIRA index %d\n', name, tt+1, tt)
                            value = value(1);
                        end
                        dataInFIRAFormat.ecodes.data(tt,Lcolumn) = value;
                    end

                else

                    % Save the name/type (always "value", just cuz)
                    Lcolumn = strcmp(enhancementName, dataInFIRAFormat.ecodes.name);
                    if ~any(Lcolumn)
                        dataInFIRAFormat.ecodes.name{end+1} = enhancementName;
                        dataInFIRAFormat.ecodes.type{end+1} = categoryName;
                        dataInFIRAFormat.ecodes.data(:,end+1) = nan;
                        Lcolumn(end+1) = true;
                    end

                    % Save the data
                    if isempty(value)
                        value = nan;
                    elseif length(value) > 1
                        fprintf('Warning: found >1 values for %s at trial file index %d, FIRA index %d\n', enhancementName, tt+1, tt)
                        value = value(1);
                    end
                    dataInFIRAFormat.ecodes.data(tt,Lcolumn) = value;
                end
            end

            % Get names of trial numeric event buffers in this category.
            numeric = dataInTrialFileFormat(tt+1).numeric_events;
            numericNames = {};
            if ~isempty(numeric)
                numericNames = fieldnames(numeric);
            end
            numericInCatecory = intersect(fieldsInCategory, numericNames);
            numNumeric = length(numericInCatecory);

            % Loop through the numeric events in this category.
            for nn = 1:numNumeric
                numericName = numericInCatecory{nn};
                eventData = dataInTrialFileFormat(tt+1).numeric_events.(numericName);

                if isempty(eventData)
                    ecodeValue = nan;
                else
                    if size(eventData, 1) > 1
                        fprintf('Warning: found >1 values for %s at trial file index %d, FIRA index %d\n', numericName, tt+1, tt)
                    end
                    % Each event can have multiple values, like:
                    %   [timestamp]
                    %   [timestamp, value]
                    %   [timestamp, value1, value2, ...]
                    % Take the last one available.
                    ecodeValue = eventData(1,end);
                end

                % Create the new ecode and save the data.
                Lcolumn = strcmp(numericName, dataInFIRAFormat.ecodes.name);
                if ~any(Lcolumn)
                    dataInFIRAFormat.ecodes.name{end+1} = numericName;
                    dataInFIRAFormat.ecodes.type{end+1} = categoryName;
                    dataInFIRAFormat.ecodes.data(:,end+1) = nan;
                    Lcolumn(end+1) = true;
                end
                dataInFIRAFormat.ecodes.data(tt,Lcolumn) = ecodeValue;
            end

            % Get names of trial text event buffers in this category.
            text = dataInTrialFileFormat(tt+1).text_events;
            textNames = {};
            if ~isempty(text)
                textNames = fieldnames(text);
            end
            textInCategory = intersect(fieldsInCategory, textNames);
            numText = length(textInCategory);

            % Loop through the text events in this category.
            for tx = 1:numText
                textName = textInCategory{tx};
                eventData = dataInTrialFileFormat(tt+1).text_events.(textName);

                if isempty(eventData)
                    ecodeValue = nan;
                else
                    if size(eventData.timestamp_data, 1) > 1
                        fprintf('Warning: found >1 values for %s at trial file index %d, FIRA index %d\n', textName, tt+1, tt)
                    end
                    % Text events have a timestamp and a text string.
                    % Take the first timestamp.
                    ecodeValue = eventData.timestamp_data(1);
                end

                % Create the new ecode and save the data.
                Lcolumn = strcmp(textName, dataInFIRAFormat.ecodes.name);
                if ~any(Lcolumn)
                    dataInFIRAFormat.ecodes.name{end+1} = textName;
                    dataInFIRAFormat.ecodes.type{end+1} = categoryName;
                    dataInFIRAFormat.ecodes.data(:,end+1) = nan;
                    Lcolumn(end+1) = true;
                end
                dataInFIRAFormat.ecodes.data(tt,Lcolumn) = ecodeValue;
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

%% Add spikes from trial numeric events.
%

% Start the data celery
% Mmmm, celery...
dataInFIRAFormat.spikes.data = {};

% Loops through the trials and collect the data
for tt = 1:numTrials

    % Get spikes from trial numeric event lists.
    numeric = dataInTrialFileFormat(tt+1).numeric_events;
    numericNames = {};
    if ~isempty(numeric)
        numericNames = fieldnames(numeric);
    end

    % But only the numeric event lists that are in spike categories.
    spikeNames = {};
    categories = dataInTrialFileFormat(tt+1).categories;
    if ~isempty(categories)
        for sc = 1:length(spikeCategories)
            spikeCategory = spikeCategories{sc};
            if isfield(categories, spikeCategory)
                spikeNames = cat(1, spikeNames, categories.(spikeCategory));
            end
        end
    end
    spikeChannelNames = intersect(numericNames, spikeNames);

    % Loop through the channels to count units
    numSpikeChannels = length(spikeChannelNames);
    for ss = 1:numSpikeChannels

        % Parse channel number
        channelNumber = str2double(spikeChannelNames{ss}(end-2:end));

        % Keep track of units
        channelUnits = [];

        % Keep track of where we're storing the data
        channelRowStart = 0;

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
