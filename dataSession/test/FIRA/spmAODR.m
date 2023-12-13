function spmAODR(func)
% function spmAODR(func)
%
% File describing how to interpret data for FIRA
% Use input argument "func" to switch functions:
%   'init'      ... fills FIRA.spm
%   'trial'     ... parse current trial (in FIRA.trial)
%   'cleanup'   ... whatever you need to do at the end
%
% Returns:
%   nada, but can change FIRA
%

% created by jig 10/21/04
% modified from Sharath's spm724rg.m
% LD 2-21-2007
% LD 9-17-08 added returnsac, omitrew Omitrew only for data after 9/15/08

% FOR DOTS TASK, CHOICE 1 IS DEFINED AS UP/LEFT; CHOICE2 IS DEFINED AS
% DOWN/RIGHT. NOTE THAT FOR REWARD ASYMMETRY TASKS, DEFINITION OF REWCONT=0
% MIGHT NOT MEAN CHOICE1 IS PAIRED WITH SMALL REWARD.

% edited by Yunshu 5-17-2016
% change how rewsize is computed given more than one beep
% now the rewsize is a little shorter than (rewoff-rewon), because rew on is
% the first time the tti pulse for valve opening (1st in the four pulses) is sent, and rew off is
% the last time the tti pulse for valve closing (3rd in the four pulses) is sent, which means the
% difference includes the time of interpulse intervel.

% Yunshu 5-23-2016
% Before 20150202 ntti should be 2n, after that, ntti should be 4n. If
% ntti is other strange number, set rewsize to 0 and rew off = rew on.

global FIRA
declareEC_ecodes_AODR;

%% if called with empty FIRA.spm, fill it and return
if strcmp(func, 'init')
    
    % useful ecode markers for tagged ecodes
    cb = 7000;
    cm = 6000;
    cx = 8000;
    
    % FIRA.spm is a struct with fields corresponding to "dataTypes"
    % that have values that are argument lists to the given dataType
    % constructor method (empty means use defaults).
    % (or, if you prefer, an nx2 cell array, with the first column
    % corresponding to the fieldnames, the second to the arglist.)
    %
    % Thus, to create FIRA with dataType "ecodes" including a list of
    % computed values called "one" and "two", use:
    % spm.ecodes = {'compute', {'one' 'value'; 'two' 'value}}
    %
    % buildFIRA_init then converts these descriptors into actual dataType
    % objects used to parse data. Note that the arglists to buildFIRA_init
    % can override the choices of dataTypes here (so you can use the same
    % base spm file to create different FIRAs)
    %
    % see FIRA/build/dataTypes for possible dataTypes
    %     FIRA/build/buildFIRA_init.m for details of how the spm struct is
    %       used
    FIRA.spm = struct( ...
        'trial',  {{ ...
        'startcd', 1005,                             ...
        'anycd',   [EC.CORRECTCD EC.WRONGCD EC.NOCHCD EC.BRFIXCD], ...
        'allcd',   [] ...
        }}, ...
        'ecodes', {{  ...
        'extract', {   ...  % values to extract & save in FIRA{2}: <name> <type> <method> <args>
        
        % Codes about Timing of task:
        'fp_on'             'time'      @getEC_time     {EC.FPONCD      1}; ...  %FP on
        'ring_on'           'time'      @getEC_time     {EC.RINGON      1};...
        'sample_on'         'time'      @getEC_time     {EC.TRGC1CD     1}; ...   %sample on
        'targ1_off'         'time'      @getEC_time     {EC.TRGC2CD     1}; ...   % code dropped when MSAC target turns off
        'mean_on'           'time'      @getEC_time     {EC.TARGONCD    1}; ...
        'target_off'        'time'      @getEC_time     {EC.TARGOFFCD   1}; ...
        'fp_off'            'time'      @getEC_time     {EC.FPOFFCD     1}; ...
        'choices_on'        'time'      @getEC_time     {EC.TRGC3CD     1};...
        'sac_on'            'time'      @getEC_time     {EC.SACMAD      1}; ...
        'targ_acq'          'time'      @getEC_time     {EC.TRGACQUIRECD 1}; ...
        'fdbk_on'           'time'      @getEC_time     {EC.FDBKONCD    1}; ...
        'fdbk_off'          'time'      @getEC_time     {EC.FDBKOFFCD   1}; ...
        'all_off'           'time'      @getEC_time     {EC.ALLOFFCD    1}; ...
        
        %Codes about correct or not, both correct target and sample are
        %considered 'correct' here'
        'online_correct'    'time'      @getEC_time     {EC.CORRECTCD    1}; ...
        'online_error'      'time'      @getEC_time     {EC.WRONGCD    1}; ...
        'online_ncerr'      'time'      @getEC_time     {EC.NOCHCD    1}; ...
        'online_brfix'      'time'      @getEC_time     {EC.BRFIXCD    1}; ...
        %Location, id, and environment info
        'fp_x'              'value'     @getEC_tagged	{EC.I_FIXXCD    cb    cm    cx  0.1};   ...
        'fp_y'              'value'     @getEC_tagged	{EC.I_FIXYCD    cb    cm    cx  0.1};   ...
        't1_x'              'value'     @getEC_tagged	{EC.I_TRG1XCD   cb    cm    cx  0.1};   ...
        't1_y'              'value'     @getEC_tagged   {EC.I_TRG1YCD   cb    cm    cx  0.1};   ...
        'sample_x'          'value'     @getEC_tagged   {EC.I_TRG3XCD   cb    cm    cx  0.1};   ...
        'sample_y'          'value'     @getEC_tagged   {EC.I_TRG3YCD   cb    cm    cx  0.1};   ...
        't2_x'              'value'     @getEC_tagged   {EC.I_TRG2XCD   cb    cm    cx  0.1};   ...
        't2_y'              'value'     @getEC_tagged   {EC.I_TRG2YCD   cb    cm    cx  0.1};   ...
        'task_id'           'id'        @getEC_tagged   {EC.I_TASKIDCD  cb    cm    cx  1};     ...
        'trial_id'          'id'        @getEC_tagged   {EC.I_TRIALIDCD cb    cm    cx  1};     ...
        'seed_base'         'value'     @getEC_tagged   {EC.I_DTVARCD   cb    cm    cx  1};     ...
        'hazard'            'value'     @getEC_tagged   {EC.I_T1IDH     cb    cm    cx  0.001};     ...
        'noise'             'value'     @getEC_tagged   {EC.I_T1SIGMA   cb    cm    cx  0.1};     ...
        
        }, ...
        'compute', { ...                % names/types of fields to compute & save later in FIRA{2}
        'correct_target'        'id';       ... % 1 or 2
        'sample_id'             'id';       ... % neg are close to T1, pos are close to T2
        't1_angle'              'id';       ...
        't2_angle'              'id';       ...
        'sample_angle'          'value';    ...
        'llr'                   'value';    ... % evidence for T1 (-) vs T2 (+)
        'target_r'              'value';    ...
        'sac_angle'             'value';    ...
        'score'                 'id';       ... % offline score: 1=correct, 0=error, -1=nc, -2=brfix,-3=sample
        'online_score'          'value';    ... % online score: 1=correct, 0=error, -1=nc, -2=brfix
        'choice'                'id';       ... % 0=neither, 1=T1, 2=T2
        'score_match'           'value';    ...
        'RT'                    'value';    ...
        'tacp'                  'id';       ... % trials after change-point
        'llr_for_switch'        'value';    ...
        'saccades_x_end'         'value';    ...
        'saccades_y_end'         'value';    ...
        'saccades_duration'      'value';    ...
        'saccades_vmax'          'value';    ...
        'saccades_vector_distance' 'value';    ...
        }, ...
        'tmp', { ...  % values to extract & use temporarily: <name> <type> <method> <args>
        }}}, ...
        'spikes',  [], ...
        'analog',  [], ...
        'dio',     []);
    
    return
    
    %% parse the trial ... check whether anything is given in FIRA.trial
elseif strcmp(func, 'trial')
    
    % note:
    %   for recordings before 4/22/09, All_OffCD is the feedback for both correct and error trials.
    %     Reward follows All_OffCD in correct trials
    %   for recordings afterward,
    %     in error trials, FDBKONCD is dropped 400 ms after TRGACQUIRECD, All_OffCD dropped 400 ms later.
    %     in correct trials, REW DIO dropped 800 ms after TRGACQUIRECD
    %             if (getFIRA_ec(tti, 'taskid')==3)
    %                 setFIRA_ec(tti, {'fdbkon'}, getFIRA_ec(tti, {'All_Off'}));
    %             end
    % get this trial index
    tti = FIRA.raw.trial.good_count;
    task_id = getFIRA_ec(tti, 'task_id');
    
    % Check
    if task_id>=1
        sfixoff='fp_off';
    else
        sfixoff='';
    end
    
    % parse saccades
    if ~isempty(sfixoff)
        
        % inline functions to compute angles, distances in deg from xy
        ang_deg  = @(x,y) mod(atan2(y,x)*180/pi+360,360);
        ang_diff = @(x,y) abs(rad2deg(angdiff(deg2rad(x),deg2rad(y))));
        
        % Set t1/t2/sample angles
        setFIRA_ec(tti, 't1_angle', ang_deg(getFIRA_ec(tti, 't1_x'), getFIRA_ec(tti, 't1_y')));
        setFIRA_ec(tti, 't2_angle', ang_deg(getFIRA_ec(tti, 't2_x'), getFIRA_ec(tti, 't2_y')));
        
        if task_id == 1
            
            % For MSAC, set target
            correctTargetAngle = getFIRA_ec(tti, 't1_angle');
            setFIRA_ec(tti, 'correct_target', 1);

        elseif any(task_id == 2:5)
            
            % For AODR, figure out sample angle, correct/error target, LLR
            % set sample angle
            setFIRA_ec(tti, 'sample_angle', ang_deg(getFIRA_ec(tti, 'sample_x'), getFIRA_ec(tti, 'sample_y')));
            
            % parse trial id
            trial_id = getFIRA_ec(tti, 'trial_id')-100*task_id;
            
            % Parse LLR
            % task_adaptiveODR3.c "Task Info" menu has P1-P9, which
            %   corresponds to the probability of showing the cue
            %   at locations far from (P1) or close to (P9) the true
            %   target
            llr_id = mod(trial_id, 9); % 0-8 for T1/T2, used below
            [~,N] = fileparts(FIRA.header.filename{1});
            if strncmp(N, 'Ci', 2) % Cicero
                if task_id == 2
                    % ORDER: P1->P9
                    ps = [0.0 0.05 0.10 0.10 0.15 0.15 0.20 0.15 0.10];
                else
                    ps = [0.0 0.05 0.10 0.10 0.15 0.15 0.20 0.15 0.10];
                end
            else % Default MrM
                ps = [0.0 0.0 0.0 0.10 0.15 0.30 0.15 0.15 0.15];
            end
            
            if trial_id < 9 % T1 is correct
                correctTargetAngle = getFIRA_ec(tti, 't1_angle');
                errorTargetAngle = getFIRA_ec(tti, 't2_angle');
                setFIRA_ec(tti, 'correct_target', 1);
                setFIRA_ec(tti, 'sample_id', llr_id-4);
                setFIRA_ec(tti, 'llr', log10(ps(llr_id+1)) - log10(ps(end-llr_id)));
            else % T2 is correct
                correctTargetAngle = getFIRA_ec(tti, 't2_angle');
                errorTargetAngle = getFIRA_ec(tti, 't1_angle');
                setFIRA_ec(tti, 'correct_target', 2);
                setFIRA_ec(tti, 'sample_id', -(llr_id-4));
                setFIRA_ec(tti, 'llr', log10(ps(end-llr_id)) - log10(ps(llr_id+1)));
            end        
        end

        % Get and parse saccades
        % sacs columns are:
        %   1. latency; 2. duration; 3. max velocity; 4. peak velocity
        %   5. end x; 6. end y; 7. raw amplitude; 8. vector amplitude
        [sacs, bf] = getFIRA_saccadesPerTrial(tti, 'fpoff', 'fp_off', 'saccadeParser', @findSaccadesAODR);
        
        if bf
            % fixbreak or other errors
            setFIRA_ec(tti, 'score', -2);
        elseif isempty(sacs) || ~isfinite(sacs(1,1))
            % ncerr
            setFIRA_ec(tti, 'score', -1);
        else
            % Loop through the sacs, find the most appropriate choice
            % Set some defaults for checking saccade
            MIN_saccades_AMPLITUDE = 5;
            MAX_saccades_AMPLITUDE = 18;
            MAX_ANGULAR_DISTANCE  = 25;
            sacs = sacs(isfinite(sacs(:,1)),:);
            if ~isempty(sacs)
                targAcq = getFIRA_ec(tti, 'targ_acq') - getFIRA_ec(tti, 'fp_off');
                sacs(sacs(:,1)>targAcq,:) = [];
            end
            savedSac = [];
            score = -1; % default: no choice
            while score < 0 && ~isempty(sacs)
                
                % Get amp, angle of current saccade
                amp = sqrt(sacs(1,5).^2 + sacs(1,6).^2);
                ang = ang_deg(sacs(1,5), sacs(1,6));
                
                if task_id==1
                    
                    % Memory task
                    % Report correct trial
                    if (amp > MIN_saccades_AMPLITUDE && amp < MAX_saccades_AMPLITUDE) && ...
                            (ang_diff(ang, getFIRA_ec(tti, 't1_angle')) < MAX_ANGULAR_DISTANCE)
                        score = 1;
                    end
                else
                    
                    % AODR, parse outcome
                    %fprintf('task id=%d, trial id=%d, llr=%.2f\n', ...
                    %    task_id, getFIRA_ec(tti, 'trial_id'), getFIRA_ec(tti, 'llr'))
                    % Check for good saccade based on length and angle
                    if amp > MIN_saccades_AMPLITUDE && amp < MAX_saccades_AMPLITUDE
                        if ang_diff(ang, correctTargetAngle) < MAX_ANGULAR_DISTANCE
                            score = 1; % Correct!
                        elseif ang_diff(ang, errorTargetAngle) < MAX_ANGULAR_DISTANCE
                            score = 0; % Error!
                        elseif ang_diff(ang, getFIRA_ec(tti, 'sample_angle')) < MAX_ANGULAR_DISTANCE
                            score = -3; % Cue!
                        end
                    end
                end
                
                % Save sac
                if size(sacs,1)==1 || ismember(score, [-3 0 1])
                    savedSac = [sacs(1, [1 2 3 5 6]) ang amp];
                end
                sacs(1,:) = [];
            end
            
            % Store saccade stats
            if ~isempty(savedSac)
                setFIRA_ec(tti, 'RT',        savedSac(1));
                setFIRA_ec(tti, 'saccades_duration',   savedSac(2));
                setFIRA_ec(tti, 'saccades_vmax',  savedSac(3));
                setFIRA_ec(tti, 'saccades_x_end',  savedSac(4));
                setFIRA_ec(tti, 'saccades_y_end',  savedSac(5));
                setFIRA_ec(tti, 'sac_angle', savedSac(6));
                setFIRA_ec(tti, 'saccades_vector_distance',   savedSac(7));
            end
            
            % Save the score and choice
            setFIRA_ec(tti, 'score', score);
            if task_id > 1
                if score < 0 || score > 1
                    setFIRA_ec(tti, 'choice', 0);
                elseif (score == 1 && correctTargetAngle==getFIRA_ec(tti, 't1_angle')) || ...
                        (score == 0 && errorTargetAngle==getFIRA_ec(tti, 't1_angle'))
                    setFIRA_ec(tti, 'choice', 1)
                else
                    setFIRA_ec(tti, 'choice', 2)
                end
            
                % PLOTZ -- set to true to see trial-by-trial eye
                % positions
                if false % getFIRA_ec(tti, 'score') == -3
                    
                    if correctTargetAngle == getFIRA_ec(tti, 't1_angle')
                        TcX = getFIRA_ec(tti, 't1_x');
                        TcY = getFIRA_ec(tti, 't1_y');
                        TeX = getFIRA_ec(tti, 't2_x');
                        TeY = getFIRA_ec(tti, 't2_y');
                    else
                        TcX = getFIRA_ec(tti, 't2_x');
                        TcY = getFIRA_ec(tti, 't2_y');
                        TeX = getFIRA_ec(tti, 't1_x');
                        TeY = getFIRA_ec(tti, 't1_y');
                    end
                    axsz = 15;
                    xsr = FIRA.analog.data(tti, 2).values;
                    ysr = FIRA.analog.data(tti, 3).values;                    
                    xs = nanrunmean(FIRA.analog.data(tti, 2).values, 40);
                    ys = nanrunmean(FIRA.analog.data(tti, 3).values, 40);
                    ts = (0:length(xs)-1) + FIRA.analog.data(tti, 2).start_time - getFIRA_ec(tti, sfixoff);
                    score = getFIRA_ec(tti, 'score');
                    co = {'m' 'k' 'c' 'r' 'g'}; 
                    
                    % Top plot is vs time
                    subplot(2,1,1); cla reset; hold on;
                    plot([-2500 1500], [0 0], 'k:');
                    plot([0 0], [-axsz axsz], 'k-');
                    plot(ts, xs, 'b-');
                    plot(ts, ys, 'c-');
                    plot(ts, xsr, 'b:');
                    plot(ts, ysr, 'c:');
                    plot(targAcq.*[1 1], [-axsz axsz], 'k-', 'LineWidth', 4);
                    plot(getFIRA_ec(tti, 'RT').*[1 1], [-axsz axsz], '-', ...
                        'Color', co{score+4}, 'LineWidth', 2);
                    %for ss = 1:size(sacs, 1)
                    %    plot(sacs(ss,[1 1]), [-axsz axsz], '-', 'Color', co{ss}, 'LineWidth', 2);
                    %end
                    axis([-2500 1500 -axsz axsz]);
                    
                    % Second is x vs y
                    subplot(2,1,2); cla reset; hold on;
                    plot(TcX, TcY, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 15);
                    plot(TeX, TeY, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 15);
                    plot(getFIRA_ec(tti, 'sac_endx'), getFIRA_ec(tti, 'sac_endy'), ...
                        'kd', 'MarkerSize', 15, 'MarkerFaceColor', co{score+4});

                    %for ss = 1:size(sacs, 1)
                    %    plot(sacs(ss, 5), sacs(ss, 6), 'd', 'MarkerSize', 15, 'Color', co{ss}, 'MarkerFaceColor', co{ss});
                    %end
                    plot(xs(ts>0), ys(ts>0), '-', 'Color', [0.5 0.5 0.5]);
                    plot(xsr(ts>0), ysr(ts>0), ':', 'Color', [0.5 0.5 0.5]);
                    plot([-10 10], [0 0], 'k:');
                    plot([0 0], [-axsz axsz], 'k:');
                    axis([-axsz axsz -axsz axsz]);
                    r = input('next')
                end
            end
        end
        
        % Compare to on-line score
        Lscore = isfinite(getFIRA_ec(tti, ...
            {'online_brfix', 'online_ncerr', 'online_error', 'online_correct'}));
        if any(Lscore)
            setFIRA_ec(tti, 'online_score', find(Lscore,1)-3); % convert to -2 -> 1
            setFIRA_ec(tti, 'score_match', getFIRA_ec(tti, 'online_score') == getFIRA_ec(tti, 'score'));
        end
    end

else % cleanup

    % Add tacp, llr_for_switch
    numTrials = size(FIRA.ecodes.data,1);
    tacp = nans(numTrials, 1);
    correctTarget = getFIRA_ec([],'correct_target');
    noChangeTrials = [false; diff(correctTarget)==0];
    counter = 0;
    for ii = 1:numTrials
        if noChangeTrials(ii) > 0
            counter = counter + 1;
        elseif noChangeTrials(ii) == 0
            counter = 0;
        else
            counter = nan;
        end
        tacp(ii) = counter;
    end
    setFIRA_ec([],'tacp',tacp);

    % values.llr_for_switch: LLR for switch
    % Get LLRs encoded as favoring T1 (-) vs T2 (+)
    llr_for_switch = getFIRA_ec([],'llr');

    % Convert all to favoring current T (+) vs not (-)
    Lct = correctTarget==1;
    llr_for_switch(Lct) = -llr_for_switch(Lct);

    % Flip signs on non-cp trials
    Lncp = tacp~=0;
    llr_for_switch(Lncp) = -llr_for_switch(Lncp);
    setFIRA_ec([],'llr_for_switch',llr_for_switch);
end
