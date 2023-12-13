classdef goldLabDataSessionConvertWithFIRA < goldLabDataSession
% classdef goldLabDataSessionConvertWithFIRA < goldLabDataSession

    properties

        % These are all arguments to 'bFile'

        % Chooses the "file reader" to read the data file
        %    'rex'             for Rex
        %    'nex'             for Nex
        %    'plx' or 'plexon' for Plexon
        fileType = [];

        % List of spike channels/units or keyword 'all'
        spikeList = 'all';

        % List of signal channels to keep
        signalList = 49:51;

        % Flag to store matlab comments
        keepMatlabCommands = false;

        % Flag to store dio commands
        keepDIO = true;

        % Any extra flags (typically used by the spm file)
        flags = false;
    end

    methods

        % Constructor takes variable argument list (property/value pairs)
        function obj = goldLabDataSessionConvertWithFIRA(varargin)
            obj = obj@goldLabDataSession(varargin{:});
        end

        % Convert file
        function [data, obj] = convert(obj)

            global FIRA

            % give some feedback
            fprintf('FIRA conversion\n')
            fprintf('raw file: <%s>\n', obj.rawFilename)
            fprintf('converted file: <%s>\n', obj.convertedFilename)
            fprintf('spm: <%s>\n', obj.convertSpecs')

            % Call this horrible function to do the work
            bFile( ...
                obj.rawFilename, ...
                obj.fileType, ...
                obj.convertSpecs, ...
                obj.convertedFilename, ...
                obj.spikeList, ...
                obj.signalList, ...
                obj.keepMatlabCommands, ...
                obj.keepDIO, ...
                obj.flags, ...
                []);

            % return a copy
            data = FIRA;
            clear FIRA
        end
    end
end