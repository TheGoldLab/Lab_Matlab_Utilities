function [fits, data] = plotAODR_sessionBehavior(options)
% function [fits, data] = plotAODR_sessionBehavior(options)
%

% Parse optional arguments
arguments
    options.filename     = -1; % tag for latest file
    options.monkey       = 'MM';
    options.sortType     = 'Sorted';
    options.converter    = 'Pyramid';
    options.convertSpecs = 'AODR_experiment';
    options.axs          = [];
end

% Get the data
[data, obj] = goldLabDataSession.convertSession( ...
    options.filename, ...
    'tag',          'AODR', ...
    'monkey',       options.monkey, ...
    'sortType',     options.sortType, ...
    'converter',    options.converter, ...
    'convertSpecs', options.convertSpecs);

% For printing on fig
[~,filename,~] = fileparts(obj.convertedFilename);
filename = strrep(filename, '_', '-');

%% Collect relevant data
% PMF
% Independent variable is "signed cues" -- which is cue location
%   signed by corresponding "LLR for switch"
hazards    = nonanunique(data.values.hazard);
numHazards = length(hazards);
sCues      = sign(data.values.llr_for_switch).*abs(data.ids.sample_id);
cues       = nonanunique(sCues);
numCues    = length(cues);

% Dependent variable is choice switch re: previous TRUE STATE
%   (note that this is not the same as a choice switch, because 
%   the previous choice might have been in error).
prevState = [nan; data.ids.correct_target(1:end-1)];
choice_switch = zeros(size(prevState));
choice_switch(data.ids.choice~=prevState)=1;
choice_switch(~ismember(data.ids.choice, [1 2]) | ~ismember(prevState,[1 2])) = nan;
Lswitch = choice_switch==1;
Lstay   = choice_switch==0;

% collect data in pdat
%   dim 1 is LLR
%   dim 2 is pmf, cmf T1/T2 choices
%   dim 3 is hazard
%   dim 4 is for data, ns, sem
pdat  = nans(numCues, 3, numHazards, 3); 
Lc1   = data.ids.choice==1;
Lc2   = data.ids.choice==2;
Ltask = data.ids.task_id>=2;
Lgood = Ltask & data.ids.score>=0;
if ~any(Lgood)
    if ~isempty(options.axs) && options.axs(1) ~= -1
        axes(options.axs(1));
        title(sprintf('%s: NO TRIALS', filename))
        fprintf('%s: NO TRIALS\n', filename)
    end
    return
end

% For each cue location
for ll = 1:numCues
    Lcue = Lgood & sCues==cues(ll);
    
    % For each hazard
    for hh = 1:numHazards
        Lh = Lcue & data.values.hazard == hazards(hh);
        
        % collect data
        pdat(ll,:,hh,1) = [ ...
            sum(Lh&Lswitch)./sum(Lh&(Lswitch|Lstay)).*100, ...
            mean(data.values.RT(Lh&Lc1), 'omitnan'), ...
            mean(data.values.RT(Lh&Lc2), 'omitnan')];
        
        % collect ns
        pdat(ll,:,hh,2) = [ ...
            sum(Lh&(Lswitch|Lstay)), sum(Lh&Lc1), sum(Lh&Lc2)];

        % collect sems
        pdat(ll,2:3,hh,3) = [ ...
            nanse(data.values.RT(Lh&Lc1)), ...
            nanse(data.values.RT(Lh&Lc2))];
    end
end

% Get logistic fit
%   input matrix is h1_bias, h2_bias, llr for switch, choice switch
ldat = cat(2, zeros(sum(Lgood), 2), sCues(Lgood), choice_switch(Lgood));
for hh = 1:numHazards
    ldat(data.values.hazard(Lgood)==hazards(hh), hh) = 1;
end
ldat = ldat(isfinite(ldat(:,4)),:);
fits = logist_fit(ldat, 'lu1');

%% Plotz
%
% Set up figure
if ~(length(options.axs)==1 && options.axs(1)==-1)

    axs = options.axs;
    if isempty(axs)
        wid     = 15; % total width
        hts     = 4.5.*ones(1,4);
        cols    = {1,2,1,2};
        [axs,~] = getPLOT_axes(1, wid, hts, cols, 2.2, [], [], 'ADPODR', true);
        set(options.axs,'Units','normalized');
    end
    fontSize = 12;
    co = {[4 94 167]./255, [194 0 77]./255};

    % 1. Show Logistic
    % Get max n to scale points
    if isfinite(axs(1))
        axes(axs(1)); cla reset; hold on;
        plot([-4 4], [50 50], 'k:');
        plot([0 0], [0 100], 'k:');
        xax = (-4:0.01:4)';
        lx  = cat(2, ones(size(xax)), xax);
        maxn = max(reshape(pdat(:,1,:,2),[],1));
        % For each hazard
        for hh = 1:numHazards
            % Plot raw data per LLR
            for ll = 1:numCues
                if pdat(ll,1,hh,2) > 0
                    plot(cues(ll), pdat(ll,1,hh,1), '.', 'MarkerSize', pdat(ll,1,hh,2)/maxn*50, ...
                        'Color', co{hh});
                end
            end
            plot(xax, logist_valLU1(fits([hh 3 4]), lx).*100, '-', 'Color', co{hh}, 'LineWidth', 2);

            % bias is just log((1-H)/H) -- see Glaze et al 2015
            h=text(-3.5, 75+(hh-1)*15, sprintf('H=%0.2f (%.2f)', ...
                1./(1+exp(-fits(hh))), hazards(hh)));
            set(h, 'Color', co{hh}, 'FontSize', fontSize)
        end
        axis([-4 4 0 100]);
        xlabel('Cue location(-=stay, +=switch)')
        ylabel('Prob(switch)');
        title(sprintf('%s: %d trials', filename, sum(Lgood)))
        set(gca,'FontSize', fontSize);
    end

    % 2. Show performance over trials
    %
    % Compute percent correct and bias relative to expected (from logistic
    % fits)
    if length(axs)>2 && isfinite(axs(2)) && isfinite(axs(3))

        % Collect some useful data
        windowSize = 25; % number of trials for running average
        choiceArr  = choice_switch(Lgood);
        cueArr     = sCues(Lgood);
        hazardArr  = data.values.hazard(Lgood);

        % Compute correct/error with respect to posterior
        postArr    = data.values.llr_for_switch(Lgood) + ...
            log(hazardArr./(1-hazardArr));
        correctArr = double(double(postArr>0)==choiceArr);
        correctArr(postArr==0) = nan;

        % For storing data to plot
        rmdat = nans(sum(Lgood), 4);

        % For computing logistic predictions
        lx = ones(windowSize,2);

        % Compute in hazard-specific blocks
        hBlocks = [1; find(diff(hazardArr)); sum(Lgood)];
        for bb = 1:length(hBlocks)-1
            blockLen = hBlocks(bb+1)-hBlocks(bb);
            if blockLen>=windowSize
                hfits = fits([find(hazardArr(hBlocks(bb))==hazards) 3 4]);
                for counter = (hBlocks(bb)+windowSize-1):(hBlocks(bb+1)-1)
                    inds = counter-windowSize+1:counter;
                    cArr = choiceArr(inds);
                    lx(:,2) = cueArr(inds);
                    ps = logist_valLU1(hfits, lx);
                    % bias real/expected
                    rmdat(counter,1:2) = [mean(cArr) mean(ps)];
                    % pct correct real/expected
                    ps(postArr(inds)<0) = 1-ps(postArr(inds)<0);
                    ps(postArr(inds)==0) = nan;
                    rmdat(counter,3:4) = ...
                        [mean(correctArr(inds), 'omitnan') ...
                        mean(ps, 'omitnan')];
                end
            end
        end

        % Bias, then pct cor
        trials = data.values.trial_num(Lgood);
        labels = {'p(switch)', 'p(cor)'};
        for vv = 1:2
            axes(axs(1+vv)); cla reset; hold on;
            plot(trials([1 end]), [0.5 0.5], 'k:');
            % For each hazard
            for hh = 1:numHazards
                Lh = hazardArr==hazards(hh);
                % points are actual, lines are expected
                % plot([trials(Lh) trials(Lh)]', rmdat(Lh,(vv-1)*2+(1:2))', '-', 'Color', co{hh});
                plot(trials(Lh), rmdat(Lh,(vv-1)*2+1), '.', 'Color', co{hh});
                plot(trials(Lh), rmdat(Lh,(vv-1)*2+2), '-', 'Color', co{hh});
            end
            set(gca, 'FontSize', fontSize);
            axis([trials([1 end])' 0 1]);
            ylabel(labels{vv});
            xlabel('Trial number');
        end
    end

    % 3. Show CMF T1 & T2 choices
    %
    if length(axs)>=4 && isfinite(axs(4))
        MIN_N = 5;
        lst = {'o-', 'o:'};
        hs  = nans(4,1);
        leg = cell(4,1);
        % Set up plot
        axes(axs(4)); cla reset; hold on;
        plot([0 0], [0 1000], 'k:');
        for cc = 1:2
            % For each hazard
            for hh = 1:numHazards
                % Plot error bars, means
                Lg = pdat(:,cc+1,hh,2)>=MIN_N;
                if any(Lg)
                    xs = cues(Lg);
                    ys = pdat(Lg,cc+1,hh,1);
                    es = pdat(Lg,cc+1,hh,3);
                    plot([xs xs]', [ys-es ys+es]', '-', 'Color', co{hh});
                    hs((cc-1)*2+hh) = plot(xs, ys, lst{cc}, 'MarkerSize', 10, ...
                        'Color', co{hh}, 'MarkerFaceColor', co{hh});
                    leg{(cc-1)*2+hh} = sprintf('T%d, H=%.2f', cc, hazards(hh));
                end
            end
        end
        legend(hs(isfinite(hs)), leg(isfinite(hs)));
        vals = reshape(pdat(Lg,2:3,:,1),[],1);
        axis([-4 4 min(vals)-50, max(vals)+50]);
        set(gca,'FontSize', fontSize);
        ylabel('RT (ms)')
        title('RT by choice and hazard')
        xlabel('Cue location(-=stay, +=switch)')
    end

    % 4. Show Eye positions (MGS and AODR)
    % Use different symbols for different target configurations
    if length(axs)>=6 && isfinite(axs(5)) && isfinite(axs(6))

        % Get all target configurations
        Ts = table2array(data.values(:, {'t1_x', 't1_y', 't2_x', 't2_y'}));
        gr = 0.9.*ones(1,3);
        tc = [153 216 201]./255;
        axlm = 18;
        %    colors = {[153 216 201]./255, [44 162 95]./255, 'm', 'k'};

        % MEMSAC
        axes(axs(5)); cla reset; hold on;
        plot([-15 15], [0 0], 'k:');
        plot([0 0], [-15 15], 'k:');
        Ltr = data.ids.task_id==1 & data.ids.score==1;
        % plot targets
        umTs = unique(Ts(Ltr,1:2), 'rows');
        plot(umTs(:,1), umTs(:,2), 'go', 'MarkerSize', 15, 'MarkerFaceColor', gr);
        % Plot eye positions
        plot(data.values.saccades_x_end(Ltr), data.values.saccades_y_end(Ltr), 'k.');
        axis([-axlm axlm -axlm axlm]);
        set(gca, 'FontSize', fontSize);
        xlabel('Horizontal position (dva)');
        ylabel('Vertical position (dva)');

        % AODR
        axes(axs(6)); cla reset; hold on;
        plot([-50 50], [0 0], 'k:');
        plot([0 0], [-50 50], 'k:');
        % Separate by choices: T1, T2, chose_cue, NC
        Lchc = [data.ids.choice==1 data.ids.choice==2, ...
            data.ids.score==-3, data.ids.score==-1];
        colors = {'k' 'k' 'm' 'c'};
        % Loop through each Target configuration
        uTs = unique(Ts(Lgood,:), 'rows');
        uTs = uTs(isfinite(uTs(:,4)),:);
        % For each set of target positions
        for rr = 1:size(uTs,1)

            % Get trials
            Ltrials = all(Ts==uTs(rr,:), 2);

            % For each choice
            for cc = 1:4

                % plot target
                if cc <= 2
                    plot(uTs(rr,(cc-1)*2+1), uTs(rr,(cc-1)*2+2), 'o', ...
                        'MarkerSize', 15, 'color', tc, ...
                        'MarkerFaceColor', gr);
                end

                % Plot eye positions
                plot(data.values.saccades_x_end(Ltrials&Lchc(:,cc)), ...
                    data.values.saccades_y_end(Ltrials&Lchc(:,cc)), ...
                    '.', 'Color', colors{cc});
            end
        end
        axis([-axlm axlm -axlm axlm]);
        set(gca, 'FontSize', fontSize);
    end
end
