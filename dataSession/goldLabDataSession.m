classdef goldLabDataSession
    % classdef goldLabDataSession
    %
    % Superclass for loading/converting/saving data files in the Gold Lab
    % Usage:
    %
    %  Load data from a converted file:
    %   data = goldLabDataSession.loadSession(filename, <optional property/value pairs>)
    %
    % Convert data from a single session:
    %   data = goldLabDataSession.convertSession(filename, <optional property/value pairs>)
    %
    %  Batch convert data from many sessions listed in a csv file:
    %   goldLabDataSession.convertSessions(csvFilename, tag, <optional property/value pairs>)
    %
    %  Save (converted) data from a single session in the standard place:
    %   goldLabDataSession.saveSession(data, filename, <optional property/value pairs>)
    %
    %  Assumes files are stored according to a particular directory
    %  structure.
    %  For raw, unsorted data:
    %       <baseDataDirectory>/<tag>/Data/<Monkey>/Raw/Unsorted/<filename>
    %  For raw, sorted data:
    %       <baseDataDirectory>/<tag>/Data/<Monkey>/Raw/Sorted/<filename>
    %  For converted (via FIRA), unsorted data:
    %       <baseDataDirectory>/<tag>/Data/<Monkey>/Converted/Unsorted/FIRA/<filename>
    %  For converted (via Pyramid), unsorted data:
    %       <baseDataDirectory>/<tag>/Data/<Monkey>/Converted/Unsorted/Pyramid/<filename>
    %  For converted (via FIRA), sorted data:
    %       <baseDataDirectory>/<tag>/Data/<Monkey>/Converted/Sorted/FIRA/<filename>
    %  For converted (via Pyramid), sorted data:
    %       <baseDataDirectory>/<tag>/Data/<Monkey>/Converted/Sorted/Pyramid/<filename>

    properties

        % Name of raw file (can include path and/or extension, or both can be
        %   set using defaults)
        rawFilename = [];

        % Name of converted file (can include path and/or extension, or both can be
        %   set using defaults)
        convertedFilename = [];

        % Where the data are found
        baseDataDirectory = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/';

        % Experiment 'tag' string used in data directory
        tag = 'AODR';

        % Monkey -- typically can be parsed from filename
        monkey = '';

        % Data are 'Sorted' or 'Unsorted'
        sortType = 'Sorted';

        % Define how to parse the data. NOTE: THIS IS THE "SPM" FILE FOR
        %   FIRA, EXPERIMENT_FILE FOR PYRAMID, JUST GIVING IT A STANDARD NAME
        convertSpecs = '';

        % Flag to force conversion even if converted file already exists
        forceConvert = false;

        % Flag indicating where to look for indexed filename
        indexedFilenameIsRaw = true;

        % Matlab format: 'trialFile' (only for Pyramid), 'FIRA', or 'dataSession'
        matlabFormat = 'dataSession';
    end

    methods

        %% Constructor method that can set properties
        %   e.g., call with:
        %   obj = goldLabDataSession('tag', <tag>);
        %
        function obj = goldLabDataSession(varargin)

            % Update local pathnames from file
            if exist('goldLabDataSessionLocalPathnames', 'file')
                obj.baseDataDirectory = goldLabDataSessionLocalPathnames('baseDataDirectory');
            end

            % Assume arguments are property/value pairs
            % Once I upgrade to R2023 I can use setProperties
            for ii = 1:2:nargin
                obj.(varargin{ii}) = varargin{ii+1};
            end
        end

        % Subclass must redefine convert methods
        function [data, obj] = convert(obj)
            data = [];
        end

        %% Load data from file
        %
        function [data, obj] = load(obj)

            % Check that file exists
            if exist(obj.convertedFilename, 'file')

                % Loop through the contents (should be 'data' and 'obj')
                for cc = who('-file', obj.convertedFilename)'
                    load(obj.convertedFilename, cc{:});
                end
            else
                data = [];
            end
        end

        %% Save converted data to file
        %
        function save(obj, data)

            % Save the data and the object using the correct name/ext
            save(obj.convertedFilename, 'data', 'obj');
        end

        function obj = updateFilenames(obj, filename)

            % Default file name/extensions
            fname = '*';
            fextConverted = '.mat';
            fextRaw = '';

            % Parse given filename
            if ischar(filename)

                % Get parts
                [fpath, fname, fext] = fileparts(filename);

                % Parse extension
                if ~isempty(fext) && ~strcmp(fext, '.mat')
                    fextRaw = fext;
                end
                fextConverted = '.mat';

                % If path given, put raw and converted in the same place
                if ~isempty(fpath)
                    obj.rawFilename = fullfile(fpath, [fname fextRaw]);
                    obj.convertedFilename = fullfile(fpath, [fname fextConverted]);
                    return
                end
            end

            % Parse monkey name
            if isempty(obj.monkey)
                if ischar(filename)
                    obj.monkey = upper(fname(1:2));
                else
                    error('Unknown Monkey')
                end
            end

            % Make default raw filename with path
            obj.rawFilename = fullfile( ...
                obj.baseDataDirectory, ...
                obj.tag, ...
                'Data', ...
                obj.monkey, ...
                'Raw', ...
                obj.sortType, ...
                [fname fextRaw]);

            % Make default converted filename with path
            className = class(obj);
            obj.convertedFilename = fullfile( ...
                obj.baseDataDirectory, ...
                obj.tag, ...
                'Data', ...
                obj.monkey, ...
                'Converted', ...
                obj.sortType, ...
                className(length('goldLabDataSessionConvertWith')+1:end), ...
                [fname fextConverted]);

            % Check for file index
            if isscalar(filename)

                % indexedFilenameIsRaw determines where to check for indexed file
                if obj.indexedFilenameIsRaw
                    files = dir(obj.rawFilename);
                else
                    files = dir(obj.convertedFilename);
                end

                % get filename from index
                if ~isempty(files)

                    % Case on sign of index
                    if filename > 0  % positive implies indexed from start
                        [~, fname, ext] = fileparts(files(filename).name);
                    elseif filename < 0 % negative implies indexed from end
                        [~, fname, ext] = fileparts(files(end+filename+1).name);
                    end

                    % Check if we found anything
                    if strcmp(fname, '*')
                        disp('No indexed file found')
                        return
                    end

                    % Put Humpty Dumpty back together again
                    [rpath, ~, fextRaw] = fileparts(obj.rawFilename);
                    [cpath, ~, fextConverted] = fileparts(obj.convertedFilename);
                    if ~isempty(ext)
                        if obj.indexedFilenameIsRaw
                            fextRaw = ext;
                        else
                            fextConverted = ext;
                        end
                    end
                    obj.rawFilename = fullfile(rpath, [fname fextRaw]);
                    obj.convertedFilename = fullfile(cpath, [fname fextConverted]);
                end
            end
        end
    end

    methods (Static)

        %% Convenience routines for loading/converting/saving
        %
        % Use these functions whenever possible.
        % varargin is optional list of property/value pairs

        %% Load converted session data
        %
        function [data, obj] = loadSession(filename, varargin)

            % Check arguments
            if nargin < 1 || isempty(filename)
                error('No file given');
            end

            % Get the appropriate object. Add load-specific default
            %   for indexedFilenameIsRaw property
            obj = goldLabDataSession.getObject(filename, varargin{:}, ...
                'indexedFilenameIsRaw', false);

            % Use that object to load the data
            [data, obj] = obj.load();

            % Check if we got something -- if not, try to convert
            if isempty(data)
                data = obj.convert();
            end

            % Convert to appropriate format
            if ~isempty(data)
                data = goldLabDataSession.convertToMatlabFormat(data, obj.matlabFormat);
            end
        end

        %% Convert session data
        %
        function [data, obj] = convertSession(filename, varargin)

            % Check arguments
            if nargin < 1 || isempty(filename)
                error('No file given');
            end

            % Get the appropriate object
            obj = goldLabDataSession.getObject(filename, varargin{:});

            % If not forcing convert, check if it already exists
            if ~obj.forceConvert

                % Try to load it
                [data, newObj] = obj.load();

                % Found it
                if ~isempty(data)

                    % Use the found data and object
                    newObj.matlabFormat = obj.matlabFormat;
                    obj = newObj;

                    % Give feedback
                    fprintf('Convert session: found existing <%s>\n', obj.convertedFilename)
                end
            else
                data = [];
            end

            % Use the object to convert the data
            if isempty(data)
                [data, obj] = obj.convert();
            end

            % Convert to appropriate format
            data = goldLabDataSession.convertToMatlabFormat(data, ...
                obj.matlabFormat, obj.rawFilename);

            % Save it
            obj.save(data)
        end

        %% Convert a batch of files from a csv spreadsheet
        %
        function convertSessions(csvFilename, tag, varargin)

            % Load the file into a table
            T = readtable(csvFilename, 'Delimiter', ',');

            % Loop through the entries
            for ii = 1:size(T, 1)

                % Be nice and check args (if blank use default)
                args = {};

                % Converter
                if any(strcmp(T.Properties.VariableNames, 'converter')) && ~isempty(T.converter{ii})
                    args = cat(2, args, 'converter', T.converter{ii});
                end

                % Sort Type
                if any(strcmp(T.Properties.VariableNames, 'sortType')) && ~isempty(T.sortType{ii})
                    args = cat(2, args, 'sortType', T.sortType{ii});
                end

                % Convert Specs
                if any(strcmp(T.Properties.VariableNames, 'convertSpecs')) && ~isempty(T.convertSpecs{ii})
                    args = cat(2, args, 'convertSpecs', T.convertSpecs{ii});
                end

                % Matlab format
                if any(strcmp(T.Properties.VariableNames, 'matlabFormat')) && ~isempty(T.matlabFormat{ii})
                    args = cat(2, args, 'matlabFormat', T.matlabFormat{ii});
                end

                % Give feedback
                fprintf('****Convert sessions: Converting <%s>\n', T.filename{ii})

                % Call convertSeession for each file given
                goldLabDataSession.convertSession( ...
                    T.filename{ii}, 'tag', tag, args{:}, varargin{:});
            end
        end

        %% Save session data
        %
        function saveSession(data, filename, varargin)

            % Check arguments
            if nargin < 2 || isempty(data) || isempty(filename)
                error('No data and/or filename given');
            end

            % Get the appropriate object
            obj = goldLabDataSession.getObject(filename, varargin{:});

            % Use that object to convert the data
            obj.save(data);
        end

        %% Get appropriate filetype
        %
        function dataOut = convertToMatlabFormat(dataIn, format, rawFilename)

            % check args
            if nargin < 1 || isempty(dataIn)
                dataOut = [];
                return
            end

            if nargin < 2 || isempty(format)
                format = 'dataSession';
            end

            if nargin < 3 || isempty(rawFilename)
                rawFilename = 'converted file';
            end

            % Default return unconverted
            dataOut = dataIn;

            % check structure fieldnames for keys
            fnames = fieldnames(dataIn);

            % Is in trialFile format
            if any(strcmp('enhancements', fnames))

                % Convert?
                if ~strcmp(format, 'trialFile')

                    % Convert to FIRA
                    dataOut = convertTrialFileToFIRA(dataIn, rawFilename);

                    % Possibly convert to sessionData
                    if strcmp(format, 'dataSession')
                        dataOut = convertFIRAToDataSession(dataOut);
                    end
                end

                % Is in FIRA format
            elseif any(strcmp('ecodes', fnames))

                % Can convert only to dataSession
                if strcmp(format, 'dataSession')
                    dataOut = convertFIRAToDataSession(dataIn);
                elseif strcmp(format, 'trialFile')
                    disp('Data in FIRA format, cannot convert to trialFile')
                end

                % Is in dataSession format already -- can't convert
            elseif any(strcmp('ids', fnames))
                if ~strcmp(format, 'dataSession')
                    fprintf('Data in dataSession format, cannot convert to <%s>\n', ...
                        format);
                end
            end
        end

        %% Get the appropriate object
        %
        function obj = getObject(filename, varargin)

            % Check varargin for converter, otherwise use default
            Lconverter = strcmp('converter', varargin);
            if any(Lconverter)
                ci = find(Lconverter);
                converter = varargin{ci+1};
                varargin(ci:ci+1) =[];
            else
                converter = 'Pyramid';
            end

            % Get object based on converter
            obj = feval(['goldLabDataSessionConvertWith' converter], ...
                varargin{:});

            % Update filenames
            obj = obj.updateFilenames(filename);
        end
    end
end