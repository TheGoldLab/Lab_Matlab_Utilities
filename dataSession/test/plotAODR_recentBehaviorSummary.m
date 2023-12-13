function plotAODR_recentBehaviorSummary(monkey, sessions, behaviorOnly)
% function plotAODR_recentBehaviorSummary(monkey, sessions, behaviorOnly)
%

% Check arguments
arguments
    monkey = 'MM';
    sessions = 5;
    behaviorOnly = true;
end

% sessions determines which sessions to plot -- either:
%   scalar: n most recent sessions
%   array: indices of sessions from file list
%   cell array: list of file names
if isempty(sessions)
    sessions = 1:5;
elseif isscalar(sessions)
    sessions = 1:sessions;
end
numSessions = length(sessions);

% Open the figure
wid     = 16; % total width
hts     = 16/numSessions.*0.9;
cols    = num2cell(ones(1,numSessions).*3);
[axs,~] = getPLOT_axes([], wid, hts, cols, 2, 1.5, [], 'ADPODR', true);
set(axs,'Units','normalized');

% Loop through the sessions
for ii = 1:numSessions
    
    % Call plot function with appropriate axes
    % If sessions is cell array of string names
    if iscell(sessions)
        plotADPODR_sessionBehavior(sessions{ii}, ...
            monkey, behaviorOnly, axs((ii-1)*3+(1:3)));
        
    else
        % Otherwise sessions is array of indices
        plotADPODR_sessionBehavior(sessions(ii), ...
            monkey, behaviorOnly, axs((ii-1)*3+(1:3)));
    end
end


