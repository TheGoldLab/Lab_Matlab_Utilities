classdef goldLabDataSessionConvertWithPyramid < goldLabDataSession
% classdef goldLabDataSessionConvertWithPyramid < goldLabDataSession
%
% Subclass of goldLabDataSession for converting raw data files using Pyramid
%
% NOTE: set machine/user-specific python search paths below

    properties

        % Python search path: string of colon-separated path names
        pythonSearchPath = '/Users/jigold/miniconda3/bin:/Users/jigold/miniconda3/condabin:';

        % Define where pyramid looks for files (single directory name)
        pyramidSearchPath = '/Users/jigold/GoldWorks/Mirror_jigold/Manuscripts/2022_dlPFC/mfiles/dataSession/pyramid/';

        % Should we also save the converted binary (hdf5) trialFile?        
        saveTrialFile = true;
    end

    methods

        % Constructor takes variable argument list (property/value pairs)
        function obj = goldLabDataSessionConvertWithPyramid(varargin)
            obj = obj@goldLabDataSession(varargin{:});

            % Update local pathnames from file (note that in
            % goldLabDataSession varargin overrides the file, but here
            % the file overrides varargin. it's just easier. Whatever.
            if exist('goldLabDataSessionLocalPathnames', 'file')
                obj.pythonSearchPath = goldLabDataSessionLocalPathnames('pythonSearchPath');
                obj.pyramidSearchPath = goldLabDataSessionLocalPathnames('pyramidSearchPath');
            end
        end

        % Convert file
        function [data, obj] = convert(obj)

            % Check for the raw file (possibly missing extension)
            [rpath, rname, rext] = fileparts(obj.rawFilename);
            if isempty(rext)
                checkRawFile = fullfile(rpath, [rname '.*']);
            else
                checkRawFile = obj.rawFilename;
            end
            checkedRawFile = dir(checkRawFile);
            if ~isempty(checkedRawFile)
                obj.rawFilename = fullfile(checkedRawFile(1).folder, ...
                    checkedRawFile(1).name);
            else
                error('Pyramid: raw file <%s> not found', obj.rawFilename);
            end

            % Check if we are saving the hdf5 file
            if obj.saveTrialFile
                [cpath, cname, ~] = fileparts(obj.convertedFilename);
                trialFilename = fullfile(cpath, [cname '.hdf5']);
            else
                trialFilename = './tmp.hdf5';
            end

            % Check that we have an experiment file
            if isempty(obj.convertSpecs)
                error('Pyramid conversion requires convertSpecs')
            else
               [cpath, cname, ~] = fileparts(obj.convertSpecs);
               obj.convertSpecs = fullfile(cpath, [cname '.yaml']);
            end

            % Make the trial file using pyramid
            % Set Python search paths
            currentPath = getenv('PATH');
            if ~contains(currentPath, obj.pythonSearchPath)
                setenv('PATH', cat(2, obj.pythonSearchPath, currentPath));
            end
            system(sprintf('pyramid convert --search-path %s --experiment %s --readers plexon_reader.plx_file=%s --trial-file %s', ...
                obj.pyramidSearchPath, obj.convertSpecs, obj.rawFilename, trialFilename));

            % Make the matlab struct in trialFile format
            trialFile = TrialFile(trialFilename);
            data = trialFile.read();
        end
    end
end
