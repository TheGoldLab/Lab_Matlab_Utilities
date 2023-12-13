function continue_ = parse(t)
% function continue_ = parse(t)
%
% Parse method for class trial, which describes how
%   trial data is added to FIRA. Reads FIRA.raw and fills
%   FIRA.raw.trial. Intended to be called by buildFIRA_parse
%   (see that file for details)
%
% Input:
%   t     ... the trial object
%
% Output:
%   continue_ ... false when no more FIRA.raw.startcds

% Copyright 2005 by Joshua I. Gold
%   University of Pennsylvania

% fills FIRA.raw.trial
global FIRA

% check for startcds
if isempty(FIRA.raw.trial.num_startcds) || FIRA.raw.trial.num_startcds == 0
    continue_ = false;
    disp('trial.parse: no trials')
    return
end

% return flag, true if there's still startcds to parse after we're done
continue_ = true;

% loop until we find a good trial
while FIRA.raw.trial.startcd_index <= FIRA.raw.trial.num_startcds

    % get raw ecodes from current trial
    % JIG AND SJ added offsets 11/9/12
    ecodes = FIRA.raw.ecodes( ...
        max(1, (FIRA.raw.trial.startcds(FIRA.raw.trial.startcd_index)+ ...
        t.startcd_offset)) : ...
        FIRA.raw.trial.startcds(FIRA.raw.trial.startcd_index+1)+ ...
        t.endcd_offset,:);
    
    % update startcds index
    FIRA.raw.trial.startcd_index = FIRA.raw.trial.startcd_index + 1;

    % update total trial counter
    FIRA.raw.trial.total_count = FIRA.raw.trial.total_count + 1;
    
    % check any, all, no code flags to determine good trial
    if(     (isempty(t.anycd) ||  any(ismember(t.anycd, ecodes))) && ...
            (isempty(t.allcd) ||  all(ismember(t.allcd, ecodes))) && ...
            (isempty(t.nocd)  || ~any(ismember(t.nocd,  ecodes))))

        % save the trial ecodes
        FIRA.raw.trial.ecodes = ecodes;
            
        % update good trial counter
        FIRA.raw.trial.good_count = FIRA.raw.trial.good_count + 1;

        % get wrt, start, end times for the current trial and store
        % in FIRA.raw.trial ... use t_.tminmax for bounds.
        % wrt determines reference time for all events -- typically
        %  FPONCD, which can occur multiple times in trials with
        %  "flashing" fp -- so select the last one
        wrt = ecodes(ecodes(:,2) == t.wrt, 1);
        if isempty(wrt)
            wrt = ecodes(1,1);
        else
            wrt = wrt(end);
        end
        FIRA.raw.trial.wrt = wrt;
        
        % start_time is time of first ecodes, 
        %   bounded by t.tminmax(1)
        if (ecodes(1,1) - wrt) < t.tminmax(1)
            FIRA.raw.trial.start_time = t.tminmax(1) + wrt;
        else
            FIRA.raw.trial.start_time = ecodes(1,1);
        end
        
        % end_time is time of last ecodes, 
        %   bounded by t.tminmax(2)
        ee = ecodes(ecodes(:,2)>0, 1);
        if (ee(end,1) - wrt) > t.tminmax(2)
            FIRA.raw.trial.end_time = t.tminmax(2) + wrt;
        else
            FIRA.raw.trial.end_time = ee(end,1);
        end
        
        % store start time of NEXT trial ecode, in case
        %   we're waiting for low-priority messages that
        %   might come after the "official" end of this trial
        ind = FIRA.raw.trial.startcds(FIRA.raw.trial.startcd_index);
        if ind > size(FIRA.raw.ecodes, 1)
            FIRA.raw.trial.next_time = FIRA.raw.ecodes(end, 1);
        else
            FIRA.raw.trial.next_time = FIRA.raw.ecodes(ind, 1);
        end

        %         disp(' ')
        %
        %         disp(sprintf('Trial parse: begin=%.4f,end=%.4f,len=%.4f,next=%.4f', ...
        %             FIRA.raw.trial.start_time - FIRA.raw.trial.wrt, ...
        %             FIRA.raw.trial.end_time   - FIRA.raw.trial.wrt, ...
        %             FIRA.raw.trial.end_time   - FIRA.raw.trial.start_time, ...
        %             FIRA.raw.trial.next_time  - FIRA.raw.trial.wrt))
        
        % possibly add extra time
        if ~isempty(t.textra)
            if t.textra > 0
                FIRA.raw.trial.end_time = ...
                    min(FIRA.raw.trial.end_time + t.textra, FIRA.raw.trial.next_time);
            elseif t.textra < 0
                FIRA.raw.trial.end_time = ...
                    max(FIRA.raw.trial.end_time + t.textra, FIRA.raw.trial.start_time);
            end
        end
        
%         % some feedback
%         disp(sprintf('Trial extra: begin=%.4f,end=%.4f,len=%.4f,next=%.4f', ...
%             FIRA.raw.trial.start_time - FIRA.raw.trial.wrt, ...
%             FIRA.raw.trial.end_time   - FIRA.raw.trial.wrt, ...
%             FIRA.raw.trial.end_time   - FIRA.raw.trial.start_time, ...
%             FIRA.raw.trial.next_time  - FIRA.raw.trial.wrt))
            
        % found good trial, outta
        return
        
    % on the last trial, if it's bad, save it (except wrt)
    elseif FIRA.raw.trial.startcd_index > FIRA.raw.trial.num_startcds
               FIRA.raw.trial.ecodes = ecodes;
            
        % update good trial counter
        FIRA.raw.trial.good_count = FIRA.raw.trial.good_count + 1;

        % get wrt, start, end times for the current trial and store
        % in FIRA.raw.trial ... use t_.tminmax for bounds.
        % wrt determines reference time for all events
        wrt = ecodes(ecodes(:,2) == t.wrt, 1);
        if isempty(wrt)
            wrt = ecodes(1,1);
        end
        FIRA.raw.trial.wrt = []; % SPECIAL FLAG
        
        % start_time is bounded by t.tminmax(1)
        if (ecodes(1,1) - wrt) < t.tminmax(1)
            FIRA.raw.trial.start_time = t.tminmax(1) + wrt;
        else
            FIRA.raw.trial.start_time = ecodes(1,1);
        end
        
        % end_time is bounded by t.tminmax(2)
        if (ecodes(end,1) - wrt) > t.tminmax(2)
            FIRA.raw.trial.end_time = t.tminmax(2) + wrt;
        else
            FIRA.raw.trial.end_time = ecodes(end,1);
        end
        
        % possibly add extra time
        if ~isempty(t.textra)
            if t.textra > 0
                FIRA.raw.trial.end_time = ...
                    min(FIRA.raw.trial.end_time + t.textra, FIRA.raw.trial.next_time);
            elseif t.textra < 0
                FIRA.raw.trial.end_time = ...
                    max(FIRA.raw.trial.end_time + t.textra, FIRA.raw.trial.start_time);
            end
        end
    else
        % disp(sprintf('TRIAL PARSE: BAD TRIAL!!!!'))
    end    
end

% no more startcds
continue_ = false;
