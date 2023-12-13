function [sacs_, bf_] = getFIRA_saccadesPerTrial(trial, varargin)
%function [sacs_, bf_] = getFIRA_saccades(trial, varargin)
%
% Convenience function recalibrate eye position to fixation, and return 
%   saccade parameters from the current trial in FIRA.spm
%
% Arguments:
%   trial        ... trial number
%   options:
%       num_saccades ... argument to findSaccades
%       recal        ... flag, whether to recalibrate
%
% Returns:
%   sacs_ ... list of saccade parameters computed by parseSaccades
%               ([] if no saccade found)
%   bf_   ... flag for (true if) broken fixation 
%               (eye record stops before fpoff)

% 6/2/07 jig added gain
%
% Copyright 2005 by Joshua I. Gold
%   University of Pennsylvania

if nargin < 1 || isempty(trial)
    trial = 1;
end

% Default options
options = struct( ...
    'num_saccades', 2, ...
    'recal',        true, ...
    'heye',         'horiz_eye', ...
    'veye',         'vert_eye', ...
    'fpoff',        'fp_off', ...
    'fpx',          'fp_x', ...
    'fpy',          'fp_x', ...
    'max_time',     1000, ... % ms
    'hgain',        1, ...
    'vgain',        1, ...
    'saccadeParser',@findSaccadesA);

for ii = 1:2:nargin-1
    options.(varargin{ii}) = varargin{ii+1};
end

global FIRA
sacs_ = [];

if isfield(FIRA, 'analog')

    % get indices of horizontal, vertical eye position data
    eyeX = strcmp(options.heye, FIRA.analog.name);
    eyeY = strcmp(options.veye, FIRA.analog.name);

    % check that we have appropriate analog data
    if any(eyeX) && ~isempty(FIRA.analog.data(trial, eyeX).values) && ...
            any(eyeY) && ~isempty(FIRA.analog.data(trial, eyeY).values)

        % check that store rate is the same for both channels
        if FIRA.analog.store_rate(eyeX) ~= FIRA.analog.store_rate(eyeY)
            disp('getFIRA_saccades: bad store rates')
        end
        
        % get the sample number in the analog data array corresponding
        % to the time of fp offset
        fpoff_samp = ceil((FIRA.ecodes.data(trial, strcmp(options.fpoff, FIRA.ecodes.name)) - ...
            FIRA.analog.data(trial, 1).start_time) * ...
            FIRA.analog.store_rate(eyeX)/1000);

        % check for enough data
        if isnan(fpoff_samp) || ...
                FIRA.analog.data(trial, eyeX).length <= fpoff_samp
           
            % return, with flag for broken fixation set
            bf_= true;
            return
        end            

        % possibly use gain to re-scale traces
        if options.hgain ~= 1
            FIRA.analog.data(trial, eyeX).values = options.hgain*...
                FIRA.analog.data(trial, eyeX).values;
        end
        if options.vgain ~= 1
            FIRA.analog.data(trial, eyeY).values = options.vgain*...
                FIRA.analog.data(trial, eyeY).values;
        end

        if options.recal

            % re-calibrate eye position to last 10 samples of fixation
            FIRA.analog.data(trial, eyeX).values = ...
                FIRA.analog.data(trial, eyeX).values - ...
                mean(FIRA.analog.data(trial, eyeX).values(fpoff_samp-10:fpoff_samp)) + ...
                FIRA.ecodes.data(trial, strcmp(options.fpx, FIRA.ecodes.name));

            FIRA.analog.data(trial, eyeY).values = ...
                FIRA.analog.data(trial,eyeY).values - ...
                mean(FIRA.analog.data(trial,eyeY).values(fpoff_samp-10:fpoff_samp)) + ...
                FIRA.ecodes.data(trial, strcmp(options.fpy, FIRA.ecodes.name));
        end

        % get saccade parameters
        if options.num_saccades > 0

            % look for at most max_time time (ms)
            end_ind = fpoff_samp + round(FIRA.analog.store_rate(1)/1000*options.max_time);
            if end_ind > FIRA.analog.data(trial, eyeX).length
                end_ind = FIRA.analog.data(trial, eyeX).length;
            end
            
            % call findSaccades to do the work
            sacs_ = feval(options.saccadeParser, ...
                FIRA.analog.data(trial, eyeX).values(fpoff_samp:end_ind), ...
                FIRA.analog.data(trial, eyeY).values(fpoff_samp:end_ind), ...
                FIRA.analog.store_rate(1), options.num_saccades, false);

            %             cla reset; hold on;
            %             xvals = FIRA.analog.data(trial, eyeX).values;
            %             yvals = FIRA.analog.data(trial, eyeY).values;
            %             tax = cumsum(ones(size(xvals)))-fpoff_samp;
            %             plot(tax, xvals, 'r-')
            %             plot(tax, yvals, 'b-')
            %             plot([0 0], [-30 30], 'k--')
            %             axis([-2500 2500 -30 30])
            %             title(sprintf('fp off = %.2f, length=%d', ...
            %                 FIRA.ecodes.data(trial, strcmp(options.fpoff, FIRA.ecodes.name)), ...
            %                 length(xvals)))
            %             r = input('next')
        end        
    end
end

% no broken fixation
bf_ = false;

