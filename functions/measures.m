function out = measures(truth, estimates)

% Evaluation measures for the LOCATA challenge
%
% Inputs:
%   truth
%   estimates
%
% Outputs:
%   out
%
% Author: Christine Evers, c.evers@imperial.ac.uk
%
% Reference: 
%       [1]	A. A. Gorji, R. T. F. "Performance measures for multiple target tracking problems," FUSION, 2011
%       [2] C. Evers, H. Loellmann, H. Mellmann, A. Schmidt, H. Barfuss, P.
%       Naylor, W. Kellermann, "The LOCATA Challenge: Acoustic Source
%       Localization and Tracking," submitted to IEEE/ACM Transactions on
%       Audio, Speech, and Language Processing, 2019.
%
% Notice: This programm is part of the LOCATA evaluation release. 
%         Please report problems and bugs to info@locata-challenge.org.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE WORK (AS DEFINED BELOW) IS PROVIDED UNDER THE TERMS OF OPEN DATA
% COMMONS ATTRIBUTION LICENSE (ODC-BY) v1.0, WHICH CAN BE FOUND AT
% http://opendatacommons.org/licenses/by/1.0/.
% THE WORK IS PROTECTED BY COPYRIGHT AND/OR OTHER APPLICABLE LAW. ANY USE
% OF THE WORK OTHER THAN AS AUTHORIZED UNDER THIS LICENSE OR COPYRIGHT LAW
% IS PROHIBITED.
%
% BY EXERCISING ANY RIGHTS TO THE WORK PROVIDED HERE, YOU ACCEPT AND AGREE
% TO BE BOUND BY THE TERMS OF THIS LICENSE. TO THE EXTENT THIS LICENSE MAY
% BE CONSIDERED TO BE A CONTRACT, THE LICENSOR GRANTS YOU THE RIGHTS
% CONTAINED HERE IN CONSIDERATION OF YOUR ACCEPTANCE OF SUCH TERMS AND
% CONDITIONS.
%
% -------------------------------------------------------------------------
%
% Representations, Warranties and Disclaimer
%
% UNLESS OTHERWISE MUTUALLY AGREED TO BY THE PARTIES IN WRITING, LICENSOR
% OFFERS THE WORK AS-IS AND MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY
% KIND CONCERNING THE WORK, EXPRESS, IMPLIED, STATUTORY OR OTHERWISE,
% INCLUDING, WITHOUT LIMITATION, WARRANTIES OF TITLE, MERCHANTIBILITY,
% FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF
% LATENT OR OTHER DEFECTS, ACCURACY, OR THE PRESENCE OF ABSENCE OF ERRORS,
% WHETHER OR NOT DISCOVERABLE. SOME JURISDICTIONS DO NOT ALLOW THE
% EXCLUSION OF IMPLIED WARRANTIES, SO SUCH EXCLUSION MAY NOT APPLY TO YOU.
%
% Limitation on Liability.
%
% EXCEPT TO THE EXTENT REQUIRED BY APPLICABLE LAW, IN NO EVENT WILL
% LICENSOR BE LIABLE TO YOU ON ANY LEGAL THEORY FOR ANY SPECIAL,
% INCIDENTAL, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES ARISING OUT OF
% THIS LICENSE OR THE USE OF THE WORK, EVEN IF LICENSOR HAS BEEN ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

az_thresh = deg2rad(30);
time_acc = 1e-12;
% OSPA_cutoff = 45;   % in deg!
OSPA_cutoff = rad2deg(az_thresh);
OSPA_p = [1., 1.5, 2., 5.];

% True sources in recording:
src_names = fieldnames(truth.source);
num_trks = length([estimates.source]);
num_srcs = length(src_names);
all_src_idx = 1:num_srcs;

% Check whether an elevation is provided at all
% Note: Some participants provide no estimates for some source(s) which
% requires a check for all sources
flag_elev = 0;
for i = 1:num_trks
    if ~isempty(estimates.source(i).elevation)
        flag_elev = 1;
        break;
    end
end

if isempty(estimates.source(1).azimuth)
    warning('No azimuth estimates for first track!')
end

% Init:
N_false = zeros(length(truth.timestamps),1);
N_valid = zeros(length(truth.timestamps),1);
N_miss = zeros(length(truth.timestamps),num_srcs);
N_swap = zeros(length(truth.timestamps),num_srcs);
az_error = nan(length(truth.timestamps),num_srcs);
el_error = nan(length(truth.timestamps),num_srcs);
error_az_ss = nan(length(truth.timestamps), num_srcs, num_trks);
error_el_ss = nan(length(truth.timestamps), num_srcs, num_trks);
OSPA_dist = nan(length(truth.timestamps), length(OSPA_p));

% Save source-to-track and track-to-source assignments
assoc_src_mat = zeros(length(truth.timestamps), num_srcs);
assoc_trk_mat = zeros(length(truth.timestamps), num_trks);

trk_active_flag = zeros(length(truth.timestamps),num_trks);
prev_nonzero_assign = zeros(1,num_srcs);
for t = 1 : length(truth.timestamps)
    % Number of active sources and number of tracks at this time stamp:
    active_src_idx = [];
    for src_idx = 1 : num_srcs
        if truth.source.(src_names{src_idx}).VAD.activity(t)
            active_src_idx = [active_src_idx, src_idx];
        end
    end
    num_active_srcs = length(active_src_idx);
    
    active_trk_idx = [];
    for trk_idx = 1 : num_trks
        match_idx = find(abs(estimates.source(trk_idx).timestamps - truth.timestamps(t)) < time_acc, 1, 'first');
        if ~isempty(match_idx) && ~isnan(estimates.source(trk_idx).azimuth(match_idx))
            active_trk_idx = [active_trk_idx, trk_idx];
            
            % Track activity (for track swaps):
            trk_active_flag(t,trk_idx) = 1;
        end
    end
    num_active_trks = length(active_trk_idx);
    active_trk_idx = find(trk_active_flag(t,:));
    
    if num_active_srcs  == 0 && num_active_trks > 0
        %% No sources
        %
        % Eval num false tracks:
        N_false(t) = num_active_trks;
        
        %% OSPA
        
        for p_idx = 1 : length(OSPA_p)            
            OSPA_dist(t,p_idx) = ( 1/num_active_trks*( OSPA_cutoff^OSPA_p(p_idx)*abs(num_active_trks) ) )^(1/OSPA_p(p_idx));            
        end
    elseif num_active_srcs > 0 && num_active_trks == 0
        %% No tracks
        %
        % Eval num missing sources
        N_miss(t,active_src_idx) = 1;
        
         %% OSPA
        
        for p_idx = 1 : length(OSPA_p)
            OSPA_dist(t,p_idx) = ( 1/num_active_srcs*( OSPA_cutoff^OSPA_p(p_idx)*abs(num_active_srcs) ) )^(1/OSPA_p(p_idx));            
        end
    elseif num_active_srcs > 0 && num_active_trks > 0
        %% Active sources and active tracks
        %
        % 1) Distance matrix:
        dist_mat_az = inf(num_srcs,num_trks);
        for trk_idx = active_trk_idx
            for src_idx = all_src_idx
                this_trk_data = estimates.source(trk_idx);
                this_trk_timestamp = find(abs(this_trk_data.timestamps - truth.timestamps(t)) < time_acc, 1, 'first');
                dist_mat_az(src_idx,trk_idx) = abs(angdiff( truth.source.(src_names{src_idx}).azimuth(t), this_trk_data.azimuth(this_trk_timestamp)));
                
                % Equation (12a) in paper: Sanity check that angdiff
                % provides same result:
                paper_dist_mat_az = abs(mod(truth.source.(src_names{src_idx}).azimuth(t) - this_trk_data.azimuth(this_trk_timestamp) + pi, 2*pi) - pi);
                if (dist_mat_az(src_idx, trk_idx) - paper_dist_mat_az) > 1e-5
                    keyboard
                end
                error_az_ss(t,src_idx,trk_idx) = dist_mat_az(src_idx,trk_idx);
            end
        end
        
        % Elevation:
        if flag_elev
            dist_mat_el = inf(num_srcs,num_trks);
            for trk_idx = active_trk_idx
                for src_idx = all_src_idx
                    this_trk_data = estimates.source(trk_idx);
                    this_trk_timestamp = find(abs(this_trk_data.timestamps - truth.timestamps(t)) < time_acc, 1, 'first');
%                     dist_mat_el(src_idx,trk_idx) = abs(angdiff( truth.source.(src_names{src_idx}).elevation(t), this_trk_data.elevation(this_trk_timestamp)));
                    dist_mat_el(src_idx,trk_idx) = abs(truth.source.(src_names{src_idx}).elevation(t) - this_trk_data.elevation(this_trk_timestamp));
                    error_el_ss(t,src_idx,trk_idx) = dist_mat_az(src_idx,trk_idx);
                end
            end
        end
        
        % Save copy of distance matrix for OSPA before gating:
        dist_mat_4OSPA = dist_mat_az;
        
        % 2) Gating
        invalid_assigns = dist_mat_az > az_thresh;
        dist_mat_az(invalid_assigns) = inf;
        
        % 3) One-to-one mapping between tracks and sources:
        % Returns 1 x num_srcs matrix, where matrix elements correspond to
        % the index of the assigned source
        assignment = munkres(dist_mat_az);
        
        %% OSPA
        
        for p_idx = 1 : length(OSPA_p)
            OSPA_dist_mat = min(OSPA_cutoff,rad2deg(dist_mat_4OSPA(:))).^OSPA_p(p_idx);
            [~,cost]= Hungarian(OSPA_dist_mat);
                        
            %Calculate final distance
            OSPA_dist(t,p_idx) = ( 1/max(num_active_trks,num_active_srcs)*( OSPA_cutoff^OSPA_p(p_idx)*abs(num_active_trks-num_active_srcs)+ cost ) )^(1/OSPA_p(p_idx));            
        end
        
        %% Track accuracy
        
        % Azimuth:
        for src_idx = active_src_idx
            if assignment(src_idx) > 0
                az_error(t,src_idx) = dist_mat_az(src_idx,assignment(src_idx));
            end
        end
        
        % Inclination
        if flag_elev
            for src_idx = active_src_idx
                if assignment(src_idx) > 0
                    el_error(t,src_idx) = dist_mat_el(src_idx,assignment(src_idx));
                end
            end
        end
        
        %% Valid and false tracks:
        unassoc_trks = nonzeros(setdiff(active_trk_idx,assignment));
        valid_trks = intersect(active_trk_idx,assignment);
        assoc_trk_mat(t,valid_trks) = all_src_idx(assignment~=0);
        N_false(t) = N_false(t) + length(unassoc_trks);
        N_valid(t) = N_valid(t) + length(valid_trks);
        
        %% Missing source detections:
        unassoc_srcs = find(assignment == 0);
        assoc_src_mat(t,all_src_idx) = assignment;
        N_miss(t,unassoc_srcs) = N_miss(t,unassoc_srcs) + 1;
        
        %% Swapped tracks
        if any(prev_nonzero_assign ~= assignment)
            changed_src_idx = find(prev_nonzero_assign ~= assignment);
            reused_trk_IDs = intersect(nonzeros(prev_nonzero_assign(changed_src_idx)), nonzeros(assignment(changed_src_idx)));
            
            for src_idx = all_src_idx
                if prev_nonzero_assign(src_idx) == 0
                    % Source newborn - Ignore
                elseif assignment(src_idx) == 0
                    % Source currently either missed or died - Ignore
                elseif assignment(src_idx) ~= prev_nonzero_assign
                    N_swap(t,src_idx) = N_swap(t,src_idx) + 1;
                    
                    % Debug:
                    disp(' ')
                    disp(['Previous nonzero assignment matrix: ', mat2str(prev_nonzero_assign)])
                    disp(['Current assignment matrix: ', mat2str(assignment)])
                    disp(['Number of missing source detections at t: ', mat2str(N_miss(t,:))])
                    disp(['Number of track swaps at t: ', mat2str(N_swap(t,:))])
                end
            end
        end
        
        % Save assignment for eval of track swaps at next time step:
        for src_idx = all_src_idx
            if assignment(src_idx) ~= 0
                prev_nonzero_assign(src_idx) = assignment(src_idx);
            end
        end
    else
        % No sources, no tracks - do nothing
    end
end

%% MOTA

MOTA = 1/num_srcs * sum( N_false + N_miss + N_swap);

%% Average number of swaps in tracks (ANST) [1]

N_swaps_per_src = nan(num_srcs,1);
for src_idx = all_src_idx
    act_N_swap = N_swap(find(truth.source.(src_names{src_idx}).VAD.activity),src_idx);
    N_swaps_per_src(src_idx) = sum(act_N_swap);
end

%% Track latency [1] & Probability of detection

track_latency = zeros(num_srcs,1);
pd = cell(num_srcs,1);
for src_idx = all_src_idx
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.source.(src_names{src_idx}).VAD.activity);
    end
    
    num_miss_srt = zeros( length(VAD_srt_idx),1);
    miss_duration = zeros(length(VAD_srt_idx),1);
    for period_idx = 1 : length(VAD_srt_idx)
        this_miss = N_miss(VAD_srt_idx(period_idx):VAD_end_idx(period_idx),src_idx);
        % Check that at the track is detected for at least 1 sample
        % AND the track delay corresponds to the number of elements in
        % N_miss of this period that are 1.
        if ~isempty(find(this_miss, 1, 'first'))
            num_miss_srt(period_idx) = find(this_miss, 1, 'first')-1;
            miss_duration(period_idx) = truth.timestamps(num_miss_srt(period_idx)+1);
        end
        pd{src_idx}(period_idx) = sum(1-this_miss)/length(this_miss);
    end
    %     track_latency(src_idx) = mean(num_miss_srt);
    track_latency(src_idx) = mean(miss_duration);
end

%% Accuracy - apply VAD results

% Azimuth error:
mean_az_error = cell(num_srcs,1);
mean_az_error_ss = cell(num_srcs,1);
for src_idx = all_src_idx
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.source.(src_names{src_idx}).VAD.activity);
    end
    
    for period_idx = 1 : length(VAD_srt_idx)
        period_az_error = az_error(VAD_srt_idx(period_idx):VAD_end_idx(period_idx),src_idx);
        nan_idx = isnan(period_az_error);
        period_az_error(nan_idx) = [];
        mean_az_error{src_idx}(period_idx) = mean(period_az_error);
        
        % Evaluate azimuth error without gating:
        if num_srcs == 1
            period_az_error_ss = error_az_ss(VAD_srt_idx(period_idx):VAD_end_idx(period_idx),src_idx,:);
            nan_idx = isnan(period_az_error_ss);
            period_az_error_ss(nan_idx) = [];
            mean_az_error_ss{src_idx}(period_idx) = mean(period_az_error_ss(:));
        end
    end
end

if flag_elev
    % Elevation error:
    mean_el_error = cell(num_srcs,1);
    mean_el_error_ss = cell(num_srcs,1);
    for src_idx = all_src_idx
        VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
        VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
        if length(VAD_end_idx) < length(VAD_srt_idx)
            VAD_end_idx(end+1) = length(truth.source.(src_names{src_idx}).VAD.activity);
        end
        
        for period_idx = 1 : length(VAD_srt_idx)
            period_el_error = el_error(VAD_srt_idx(period_idx):VAD_end_idx(period_idx),src_idx);
            nan_idx = isnan(period_el_error);
            period_el_error(nan_idx) = [];
            mean_el_error{src_idx}(period_idx) = mean(period_el_error);
            
            % Evaluate elevation error without gating:
            period_el_error_ss = error_el_ss(VAD_srt_idx(period_idx):VAD_end_idx(period_idx),src_idx,:);
            nan_idx = isnan(period_el_error_ss);
            period_el_error_ss(nan_idx) = [];
            mean_el_error_ss{src_idx}(period_idx) = mean(period_el_error_ss(:));
        end
    end
end

%% Rate of track fragmentation [1] over voice-active periods

% Track swap rate
track_swap_rate = nan(num_srcs,1);
for src_idx = all_src_idx
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.source.(src_names{src_idx}).VAD.activity);
    end
    
    TSR{src_idx} = zeros(length(VAD_srt_idx),1);
    for period_idx = 1 : length(VAD_srt_idx)
        period_duration = truth.timestamps(VAD_end_idx(period_idx)) - truth.timestamps(VAD_srt_idx(period_idx));
        TSR{src_idx}(period_idx) = sum(N_swap(VAD_srt_idx(period_idx):VAD_end_idx(period_idx),src_idx))/period_duration;
    end
    track_swap_rate(src_idx) = mean(TSR{src_idx});
end

if any(isnan( track_swap_rate ) ) > 0
    warning('NaN in TFR detected!')
end

% Track fragmentation rate
track_frag_rate = nan(num_srcs,1);
for src_idx = all_src_idx
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.source.(src_names{src_idx}).VAD.activity);
    end
    
    % Rate of live tracks being interrupted by missing detections or
    % affected by track swaps:
    TFR{src_idx} = zeros(length(VAD_srt_idx),1);
    for period_idx = 1 : length(VAD_srt_idx)
        this_period = VAD_srt_idx(period_idx):VAD_end_idx(period_idx);
        period_duration = truth.timestamps(VAD_end_idx(period_idx)) - truth.timestamps(VAD_srt_idx(period_idx));
        N_fragment = length(find(N_miss(this_period(1:end-1),src_idx) == 1 & N_miss(this_period(2:end),src_idx) == 0)) + sum(N_swap(VAD_srt_idx(period_idx):VAD_end_idx(period_idx),src_idx));
        TFR{src_idx}(period_idx) = sum(N_fragment)/period_duration;
    end
    track_frag_rate(src_idx) = mean(TFR{src_idx});
end

if any(isnan( track_frag_rate ) ) > 0
    warning('NaN in TFR detected!')
end

%% False alarm rate [1] over recording duration

rec_duration = truth.timestamps(end)-truth.timestamps(1);
FAR = sum(N_false)/rec_duration;
if any(isnan(FAR))
    keyboard
end

%% False alarm rate only within voice-active periods [new in TASL paper]

VAP = zeros(1,length(truth.timestamps));
for src_idx = all_src_idx
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.source.(src_names{src_idx}).VAD.activity);
    end
    
    for period_idx = 1 : length(VAD_srt_idx)
        VAP(VAD_srt_idx(period_idx):VAD_end_idx(period_idx)) = 1;
    end
end
% N_false_VAP = N_false(find(VAP));
VAP_srt_idx = find(VAP(2:end) & ~VAP(1:end-1)) + 1;
VAP_end_idx = find(VAP(1:end-1) & ~VAP(2:end));
if length(VAP_end_idx) < length(VAP_srt_idx)
    VAP_end_idx(end+1) = length(truth.timestamps);
end
for period_idx = 1 : length(VAP_srt_idx)
    VAP_duration = truth.timestamps(VAP_end_idx(period_idx))-truth.timestamps(VAP_srt_idx(period_idx));
    FAR_VAP(period_idx) = sum(N_false(VAP_srt_idx(period_idx):VAP_end_idx(period_idx)))/VAP_duration;
end
FAR_VAP = mean(FAR_VAP);

%% Outptus

out.N_swaps_per_src = N_swaps_per_src;
out.N_valid = N_valid;
out.MOTA = MOTA;
out.track_swap_rate = track_swap_rate;
out.track_frag_rate = track_frag_rate;
out.track_latency = track_latency;
out.FAR = FAR;
out.FAR_VAP = FAR_VAP;
out.pd = pd;
out.mean_az_error = mean_az_error;
out.mean_az_error_ss = mean_az_error_ss;
if flag_elev
    out.mean_el_error = mean_el_error;
    out.mean_el_error_ss = mean_el_error_ss;
end
out.OSPA = OSPA_dist;

%% Plots for sanity check

scrn = get(0,'Screensize');
fig_size = scrn;
fig_size(3) = fig_size(3)/2;
f_ref = figure('Position', fig_size);
num_cols = 2 + 3 + 1; % VAD activity, azimuth against ground truth, N_miss, N_flag, N_swap + 1 for table with measures

subplot(num_cols,1,1); hold on;
for src_idx = 1 : num_srcs
    plot( truth.timestamps, truth.source.(src_names{src_idx}).VAD.activity )
end
xlabel('Time [s]')
ylabel('VAD')
grid on;
v = axis; v(2) = truth.timestamps(end); axis(v);

subplot(num_cols,1,2); hold on;
for src_idx = 1 : num_srcs
    plot( truth.timestamps, rad2deg(truth.source.(src_names{src_idx}).azimuth) )
end
for trk_idx = 1 : num_trks
    plot( estimates.source(trk_idx).timestamps, rad2deg(estimates.source(trk_idx).azimuth), 'x')
end
xlabel('Time [s]')
ylabel('Azimuth [deg]')
grid on;
v = axis; v(2) = truth.timestamps(end); v(3) = -180; v(4) = 180;  axis(v);

subplot(num_cols,1,3); hold on;
for src_idx = 1 : num_srcs
    plot(truth.timestamps, N_miss(:,src_idx))
end
xlabel('Time [s]')
ylabel('Missing track flag')
grid on;
v = axis; v(2) = truth.timestamps(end); axis(v);

subplot(num_cols,1,4); hold on;
plot(truth.timestamps, N_false)
xlabel('Time [s]')
ylabel('False tracks')
grid on;
v = axis; v(2) = truth.timestamps(end); axis(v);

subplot(num_cols,1,5); hold on;
for src_idx = 1 : num_srcs
    plot(truth.timestamps, N_swap(:,src_idx))
end
xlabel('Time [s]')
ylabel('Swapped tracks')
grid on;
v = axis; v(2) = truth.timestamps(end); axis(v);

% Add table of measures in final subplto:
subplot(num_cols,1,6); plot(3);
pos = get(subplot(num_cols,1,6),'position');
delete(subplot(num_cols,1,6))
% Create uitable - NOTE: uitable uses variable names as column names. Write
% into new variables with appropriate names:
RowNames = cell(num_srcs+1,1);
for src_idx = 1 : num_srcs
    RowNames{src_idx,1} = src_names{src_idx};
end
RowNames{end,1} = 'Full scenario';
tmp_FAR = FAR; FAR = nan(num_srcs+1,1); FAR(end) = tmp_FAR;
avg_TL = nan(num_srcs+1,1); avg_TL(1:num_srcs,1) = track_latency;
TSR = nan(num_srcs+1,1); TSR(1:num_srcs,1) = track_swap_rate;
avg_pd = nan(num_srcs+1,1); avg_pd(1:num_srcs,1) = cellfun(@mean,pd)';
num_valid_trks = nan(num_srcs+1,1); num_valid_trks(end,1) = sum(N_valid,1);
num_miss_trks = sum(N_miss,1)'; num_miss_trks(end+1,1) = nan;
num_false_trks = nan(num_srcs+1,1); num_false_trks(end) = sum(N_false);
num_swap_trks = sum(N_swap,1)'; num_swap_trks(end+1,1) = nan;
mat_tab = table(num_valid_trks, num_miss_trks, num_false_trks, num_swap_trks, avg_pd,FAR,avg_TL,TSR,'RowNames',RowNames);
ui_tab = uitable('Data',mat_tab{:,:},'ColumnName',mat_tab.Properties.VariableNames,...
    'RowName',mat_tab.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
set(ui_tab,'units','normalized')
set(ui_tab,'position',pos)

out.f_ref = f_ref;

%% Plots for presentation

truth_col = [9,49,69;
    16,120,150;
    130,147,86;
    188,161,54;
    194,87,26;
    154,38,23]/256;
est_col = flipud([60,100,120;
    67,171,201;
    181,198,137;
    239,212,105;
    245,139,76;
    205,89,74]/256);
LineWidth = 2;
FontSize = 22;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');
pos = [100,100,1250,350]; % position and size

out.f_presentation = [];
out.f_label = {};

% VAD
f_ref = figure('Position', pos); hold on;
leg = cell(num_srcs,1);
h_VAD = zeros(num_srcs,1);
for src_idx = 1 : num_srcs
    h_VAD(src_idx) = plot( truth.timestamps, truth.source.(src_names{src_idx}).VAD.activity, '-', 'color', truth_col(src_idx,:), 'LineWidth', LineWidth );
    leg{src_idx} = ['Source ', num2str(src_idx)];
end
v = axis; v(2) = truth.timestamps(end); v(3) = 0; v(4) = 1.1;  axis(v);
xlabel('Time [s]')
ylabel('Voice activity')
grid on;
legend(h_VAD, leg)
set(gca,'YTick',[0,1])
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
out.f_presentation(end+1) = f_ref;
out.f_label{end+1} = 'VAD';

% Azimuth
f_ref = figure('Position', pos); hold on;
leg = cell(num_srcs+num_trks,1);
h_az = zeros(num_srcs+num_trks,1);
for src_idx = 1 : num_srcs
    h_az(src_idx) = plot( truth.timestamps, rad2deg(truth.source.(src_names{src_idx}).azimuth), '-', 'color', truth_col(src_idx,:), 'LineWidth', LineWidth );
    leg{src_idx} = ['Source ', num2str(src_idx)];
end
for trk_idx = 1 : num_trks
    plot( estimates.source(trk_idx).timestamps, rad2deg(estimates.source(trk_idx).azimuth), '.', 'color', est_col(mod(trk_idx,size(est_col,1))+1,:), 'LineWidth', LineWidth)
    h_az(trk_idx+num_srcs) = plot( estimates.source(trk_idx).timestamps(1:10:end), rad2deg(estimates.source(trk_idx).azimuth(1:10:end)), 'o', 'color', est_col(mod(trk_idx,size(est_col,1))+1,:), 'LineWidth', LineWidth);
    leg{num_srcs+trk_idx} = ['Track ', num2str(trk_idx)];
end
v = axis; v(2) = truth.timestamps(end); v(3) = -180; v(4) = 180;  axis(v);
ylimits = get(gca, 'YLim');
% Plot voice activity:
for src_idx = 1 : num_srcs
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.timestamps);
    end
    x_VAD = [truth.timestamps(1), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(1)];
    y_VAD = [ylimits(1); ylimits(1)*ones(length(VAD_srt_idx(:)),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(1)*ones(length(VAD_end_idx),1); ylimits(1)];
    [sorted_x_VAD,sort_idx] = sort(x_VAD, 'ascend');
    sorted_y_VAD = y_VAD(sort_idx);
    h_az(num_srcs+num_trks+src_idx) = patch(sorted_x_VAD, sorted_y_VAD, truth_col(src_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    leg{num_srcs+num_trks+src_idx} = ['VAD Source ', num2str(src_idx)];
end
xlabel('Time [s]')
ylabel('Azimuth [deg]')
grid on;
legend(h_az, leg)
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
out.f_presentation(end+1) = f_ref;
out.f_label{end+1} = 'az';

% Elevation    
if sum(cellfun(@isempty,{estimates.source(:).elevation})) == 0
    f_ref = figure('Position', pos); hold on;
    leg_el = cell(num_srcs+num_trks,1);
    h_el = zeros(num_srcs+num_trks,1);
    for src_idx = 1 : num_srcs
        h_el(src_idx) = plot( truth.timestamps, rad2deg(truth.source.(src_names{src_idx}).elevation), '-', 'color', truth_col(src_idx,:), 'LineWidth', LineWidth );
        leg_el{src_idx} = ['Source ', num2str(src_idx)];
    end
    for trk_idx = 1 : num_trks
        plot( estimates.source(trk_idx).timestamps, rad2deg(estimates.source(trk_idx).elevation), '.', 'color', est_col(mod(trk_idx,size(est_col,1))+1,:), 'LineWidth', LineWidth)
        h_el(trk_idx+num_srcs) = plot( estimates.source(trk_idx).timestamps(1:10:end), rad2deg(estimates.source(trk_idx).elevation(1:10:end)), 'o', 'color', est_col(mod(trk_idx,size(est_col,1))+1,:), 'LineWidth', LineWidth);
        leg_el{num_srcs+trk_idx} = ['Track ', num2str(trk_idx)];
    end
    v = axis; v(2) = truth.timestamps(end); v(3) = -180; v(4) = 180;  axis(v);
    ylimits = get(gca, 'YLim');
    % Plot voice activity:
    for src_idx = 1 : num_srcs
        VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
        VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
        if length(VAD_end_idx) < length(VAD_srt_idx)
            VAD_end_idx(end+1) = length(truth.timestamps);
        end
        x_VAD = [truth.timestamps(1), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(1)];
        y_VAD = [ylimits(1); ylimits(1)*ones(length(VAD_srt_idx(:)),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(1)*ones(length(VAD_end_idx),1); ylimits(1)];
        [sorted_x_VAD,sort_idx] = sort(x_VAD, 'ascend');
        sorted_y_VAD = y_VAD(sort_idx);
        h_el(num_srcs+num_trks+src_idx) = patch(sorted_x_VAD, sorted_y_VAD, truth_col(src_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        leg_el{num_srcs+num_trks+src_idx} = ['VAD Source ', num2str(src_idx)];
    end
    xlabel('Time [s]')
    ylabel('Elevation [deg]')
    grid on;
    legend(h_el, leg_el)
    set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
    out.f_presentation(end+1) = f_ref;
    out.f_label{end+1} = 'el';
end

% Range
f_ref = figure('Position', pos); hold on;
if sum(cellfun(@isempty,{estimates.source(:).range})) == 0
    leg_r = cell(num_srcs+num_trks,1);
    h_r = zeros(num_srcs+num_trks,1);
else
    leg_r = cell(num_srcs,1);
    h_r = zeros(num_srcs,1);
end
for src_idx = 1 : num_srcs
    r = zeros(length(truth.timestamps),1);
    for t = 1 : length(truth.timestamps)
        % Apply rotation / translation of array to source
        R = squeeze(truth.array.rotation(:,t,:));
        p = truth.array.position(:,t);
        h = truth.source.(src_names{src_idx}).position(:,t);
        v = R'*(h - p);
        
        % Spherical position
        % NOTE: azimuth = 0 deg is defined along the positive y-axis,
        %       elevation = 0 deg is defined along the positive z-axis
        [tmp_az,tmp_el,r(t)] = mycart2sph(v(1), v(2), v(3));
        if tmp_az ~= truth.source.(src_names{src_idx}).azimuth
            keyboard
        end
        if tmp_el ~= truth.source.(src_names{src_idx}).elevation
            keyboard
        end
    end
    h_r(src_idx) = plot( truth.timestamps, r, '-', 'color', truth_col(src_idx,:), 'LineWidth', LineWidth );
    leg_r{src_idx} = ['Source ', num2str(src_idx)];
end
if sum(cellfun(@isempty,{estimates.source(:).range})) == 0
    for trk_idx = 1 : num_trks
        plot( estimates.source(trk_idx).timestamps, estimates.source(trk_idx).range, '.', 'color', est_col(mod(trk_idx,size(est_col,1))+1,:), 'LineWidth', LineWidth)
        h_r(trk_idx+num_srcs) = plot( estimates.source(trk_idx).timestamps(1:10:end), estimates.source(trk_idx).range(1:10:end), 'o', 'color', est_col(mod(trk_idx,size(est_col,1))+1,:), 'LineWidth', LineWidth);
        leg_r{num_srcs+trk_idx} = ['Track ', num2str(trk_idx)];
    end
end
% v = axis; v(2) = truth.timestamps(end); v(3) = -180; v(4) = 180;  axis(v);
v = axis; v(2) = truth.timestamps(end); axis(v); v = axis;
ylimits = get(gca, 'YLim');
% Plot voice activity:
for src_idx = 1 : num_srcs
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.timestamps);
    end
    x_VAD = [truth.timestamps(1), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(1)];
    y_VAD = [ylimits(1); ylimits(1)*ones(length(VAD_srt_idx(:)),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(1)*ones(length(VAD_end_idx),1); ylimits(1)];
    [sorted_x_VAD,sort_idx] = sort(x_VAD, 'ascend');
    sorted_y_VAD = y_VAD(sort_idx);
    h_r(end+1) = patch(sorted_x_VAD, sorted_y_VAD, truth_col(src_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    leg_r{end+1} = ['VAD Source ', num2str(src_idx)];
end
xlabel('Time [s]')
ylabel('Range [m]')
grid on;
legend(h_r, leg_r)
axis(v)
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
out.f_presentation(end+1) = f_ref;
out.f_label{end+1} = 'range';

% Missed detections
f_ref = figure('Position', pos); hold on;
leg_miss = cell(num_srcs,1);
for src_idx = 1 : num_srcs
    h_miss(src_idx) = plot(truth.timestamps, N_miss(:,src_idx), '-', 'color', est_col(src_idx,:), 'LineWidth', LineWidth);
    leg_miss{src_idx} = ['Source ', num2str(src_idx)];
end
ylimits = get(gca, 'YLim');
% Plot voice activity:
for src_idx = 1 : num_srcs
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.timestamps);
    end
    x_VAD = [truth.timestamps(1), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(1)];
    y_VAD = [ylimits(1); ylimits(1)*ones(length(VAD_srt_idx(:)),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(1)*ones(length(VAD_end_idx),1); ylimits(1)];
    [sorted_x_VAD,sort_idx] = sort(x_VAD, 'ascend');
    sorted_y_VAD = y_VAD(sort_idx);
    h_miss(num_srcs+src_idx) = patch(sorted_x_VAD, sorted_y_VAD, truth_col(src_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    leg_miss{num_srcs+src_idx} = ['VAD Source ', num2str(src_idx)];
end
xlabel('Time [s]')
ylabel('Number of missed estimates')
grid on;
v = axis; v(2) = truth.timestamps(end); axis(v);
legend(h_miss, leg_miss)
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
out.f_presentation(end+1) = f_ref;
out.f_label{end+1} = 'N_miss';

% Number of false alarms
f_ref = figure('Position', pos); hold on;
h_false = plot(truth.timestamps, N_false, 'k-', 'LineWidth', LineWidth);
leg_false{1} = ['Number of false alarms'];
ylimits = get(gca, 'YLim');
% Plot voice activity:
for src_idx = 1 : num_srcs
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.timestamps);
    end
    x_VAD = [truth.timestamps(1), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(1)];
    y_VAD = [ylimits(1); ylimits(1)*ones(length(VAD_srt_idx(:)),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(1)*ones(length(VAD_end_idx),1); ylimits(1)];
    [sorted_x_VAD,sort_idx] = sort(x_VAD, 'ascend');
    sorted_y_VAD = y_VAD(sort_idx);
    h_false(1+src_idx) = patch(sorted_x_VAD, sorted_y_VAD, truth_col(src_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    leg_false{1+src_idx} = ['VAD Source ', num2str(src_idx)];
end
xlabel('Time [s]')
ylabel('Number of false estimates')
grid on;
legend(h_false, leg_false);
v = axis; v(2) = truth.timestamps(end); axis(v);
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
out.f_presentation(end+1) = f_ref;
out.f_label{end+1} = 'N_false';

% Number of track swaps
f_ref = figure('Position', pos); hold on;
for src_idx = 1 : num_srcs
    h_swap(src_idx) = plot(truth.timestamps, N_swap(:,src_idx), '-', 'color', est_col(src_idx,:), 'LineWidth', LineWidth);
    leg_swap{src_idx} = ['Source ', num2str(src_idx)];
end
ylimits = get(gca, 'YLim');
% Plot voice activity:
for src_idx = 1 : num_srcs
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.timestamps);
    end
    x_VAD = [truth.timestamps(1), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(1)];
    y_VAD = [ylimits(1); ylimits(1)*ones(length(VAD_srt_idx(:)),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(1)*ones(length(VAD_end_idx),1); ylimits(1)];
    [sorted_x_VAD,sort_idx] = sort(x_VAD, 'ascend');
    sorted_y_VAD = y_VAD(sort_idx);
    h_swap(num_srcs+src_idx) = patch(sorted_x_VAD, sorted_y_VAD, truth_col(src_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    leg_swap{num_srcs+src_idx} = ['VAD Source ', num2str(src_idx)];
end
xlabel('Time [s]')
ylabel('Swapped tracks')
legend(h_swap, leg_swap);
grid on;
v = axis; v(2) = truth.timestamps(end); axis(v);
legend(h_swap,leg_swap);
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
out.f_presentation(end+1) = f_ref;
out.f_label{end+1} = 'N_swap';

% OSPA
f_ref = figure('Position', pos); hold on;
cmap = parula(length(OSPA_p));
for p_idx = 1 : length(OSPA_p)
    h_OSPA(p_idx) = plot(truth.timestamps, OSPA_dist(:,p_idx), '-', 'color', cmap(p_idx,:), 'LineWidth', LineWidth);
    leg_OSPA{p_idx} = ['OSPA $p=$', num2str(OSPA_p(p_idx))];
end
v = axis; v(2) = truth.timestamps(end); axis(v);
ylimits = get(gca, 'YLim');
% Plot voice activity:
for src_idx = 1 : num_srcs
    VAD_srt_idx = find(truth.source.(src_names{src_idx}).VAD.activity(2:end) & ~truth.source.(src_names{src_idx}).VAD.activity(1:end-1)) + 1;
    VAD_end_idx = find(truth.source.(src_names{src_idx}).VAD.activity(1:end-1) & ~truth.source.(src_names{src_idx}).VAD.activity(2:end));
    if length(VAD_end_idx) < length(VAD_srt_idx)
        VAD_end_idx(end+1) = length(truth.timestamps);
    end
    x_VAD = [truth.timestamps(1), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_srt_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(VAD_end_idx(:)), truth.timestamps(1)];
    y_VAD = [ylimits(1); ylimits(1)*ones(length(VAD_srt_idx(:)),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(2)*ones(length(VAD_srt_idx),1); ylimits(1)*ones(length(VAD_end_idx),1); ylimits(1)];
    [sorted_x_VAD,sort_idx] = sort(x_VAD, 'ascend');
    sorted_y_VAD = y_VAD(sort_idx);
    h_OSPA(length(OSPA_p)+src_idx) = patch(sorted_x_VAD, sorted_y_VAD, truth_col(src_idx,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    leg_OSPA{length(OSPA_p)+src_idx} = ['VAD Source ', num2str(src_idx)];
end
xlabel('Time [s]')
ylabel('OSPA [deg]')
legend(h_OSPA, leg_OSPA);
grid on;
legend(h_OSPA,leg_OSPA);
set(findall(gcf,'-property','FontSize'),'FontSize',FontSize)
out.f_presentation(end+1) = f_ref;
out.f_label{end+1} = 'OSPA';

end