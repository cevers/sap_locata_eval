function main_evaluation()

% Main function for LOCATA evaluation. This script produces (among many others) all figures and
% tables provided in [1].
%
% Inputs: Please adjust l.84-113 as appropriate (see comments below)
%
% Outputs: Creates directory as specified in results_dir containing:
%   - averaged_measures_[date].mat - structure containing:
%       * 'all_measures' - Evaluation measures per submission, task, array,
%         and recording.
%       * 'avg_measures' - Evaluation measures per submission, task, and
%          array, averaged over all corresponding recordings
%       * 'p_meta_data' - Meta data of submissions (author, ID, etc)
%   - One sub-dir per submission, task, recording, array
%     (results_dir/[submission_ID]/task[task_ID]/recording[recording_ID]/[array_ID]), containing:
%       * measures_presentation_[eval_measure].fig + measures_presentation_[eval_measure].png - time series of eval
%         measure
%       * measures.fig + measure.png - Subplots of VAD, az error, missing
%         track, FAR, swapped tracks, + counters used for eval measures
%   - Sub-dir containing statists of evaluation results
%       (results_dir/stats), containing:
%       * Bar-plots for each eval measure, per task
%         (results_dir/stats/[eval_measure]_task[task_ID].fig + results_dir/stats/[eval_measure]_task[task_ID].png)
%       * LaTeX tables for azimuth errors w/ and w/o gating
%         (results_dir/stats/az_error.tex + results_dir/stats/az_error_diff_gating.tex)
%
% Authors: Christine Evers, c.evers@imperial.ac.uk
%          Heiner Loellmann, Loellmann@LNT.de
%
% Reference: 
%       [1] C. Evers, H. Loellmann, H. Mellmann, A. Schmidt, H. Barfuss, P.
%       Naylor, W. Kellermann, "The LOCATA Challenge: Acoustic Source
%       Localization and Tracking," submitted to IEEE/ACM Transactions on
%       Audio, Speech, and Language Processing, 2019.
%
% Notice: This programm is part of the LOCATA corpus release. 
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

clear all; close all; clc;

% addpath(genpath('./'))
addpath(genpath('./utils/'))
addpath(genpath('./functions/'))  % main functions

format compact

%% Main settings

load_data_flag = 1;         % if 1, load mat-file with ground-truth and participants' data,
                            % if 0, this file is generated from the submissions and eval database
flag_save_data = 1;         % if 1, the results are saved
flag_baseline = 1;          % if 1, baseline evaluation results are appended to submissions
add_my_results_flag = 0;    % if 1, new results are included for evlauation, that were not submitted for the LOCATA challenge.

%% Settings for reading data from databases
%%%%%%%%%%%%%% Needs to be adapted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if load_data_flag
    % Mat file containing all participants submissions + ground truth data
    load_data_fname = '../data/p_data_and_GT_2020-01-25_130935.mat';
else
    % Folder with the submission results of the participants
    p_dir = '../../LOCATA_databases/Submissions/';
    
    % Folder with final Eval database with all ground-truth data
%     eval_dir = '../../LOCATA_databases/Eval_with_GT_corrected/';
    eval_dir = '~/Downloads/LOCATA_Renamed/datasets_anonymous/eval/';
    
    % Folder with ground-truth vad for eval database
    eval_vad_dir = eval_dir;
end

% Directory containing any new results to be included for evaluation, that were not submitted to the
% LOCATA challenge:
if add_my_results_flag
    n_dir = '../my_results/';
end

%% Settings for writing results

% Folder for generated ground-truth and participants data
% This is the folder that contains load_data_fname!
if ~load_data_flag
    data_dir = '../data';
    if flag_save_data
        if ~exist( data_dir, 'dir')
            mkdir( data_dir )
        end
    end
end

% Directory for generated evaluation results
results_dir = '../results_paper_v10_elcorrected/';
log_fname_root = [results_dir, filesep, 'log'];
disp(['Writing results to ', results_dir])
disp(['Log files in ', log_fname_root])

%% Load data

% Loading of submissions / baseline results - do not change the following settings:
num_submissions = 16;
valid_submission_ids = [1:4, 6:12, 15, 16];
if flag_baseline
    baseline_ID = 17;
end

opts = init();

% Read participants' submissions + ground truth data:
if load_data_flag
    disp(['Loading data from previously generated mat file: ', load_data_fname])
    load( load_data_fname )
else
    disp('Reading data from original submissions and database.')
    disp(['Participants submissions: ', p_dir])
    disp(['Evaluation database: ', eval_dir])
    
    % Read participants' submissions:
    [ estimates, p_meta_data ] = read_participants_results( p_dir, eval_dir, valid_submission_ids, opts ); % participants' results and meta data
    if num_submissions ~= length(p_meta_data) || num_submissions ~= size(estimates,1)
        error(['Valid IDs of LOCATA submissions are 1, ..., 16. You specified: ', valid_submission_ids])
    end
    
    % Read ground truth data:
    truth = read_ground_truth( eval_dir, eval_vad_dir ); % ground-truth data
    
    % Read baseline results to submissions:
    if flag_baseline
        % Append baseline ID to submission IDs:
        if baseline_ID ~= (num_submissions+1)
            error(['Baseline ID should correspond to num_submissions + 1. There are ', num_submissions, ', but you specified the baseline ID as: ', baseline_ID])
        end
        ids = [valid_submission_ids, baseline_ID];
        
        % Load baseline results from file:
        baseline_fname = [data_dir,filesep,'baseline_results_merged_2018-09-18_034612.mat'];
        disp(['Loading baseline results from file: ', baseline_fname])
        load(baseline_fname) % loads variable merged_DOA_estimator
        
        p_meta_data{baseline_ID} = participants_data( baseline_ID );
        estimates = append_baseline_results( estimates, merged_DOA_estimator ); % append baseline results
    else
        ids = valid_submission_ids;
    end
    
    % Compress and save participants' results + ground truth data to file:
    if flag_save_data
        save_fname = [ data_dir filesep 'p_data_and_GT_', datestr(now, 'yyyy-mm-dd_HHMMSS'), '.mat'];
        disp(['Compressing and saving submissions + ground truth to mat file: ', save_fname])
        save(save_fname, 'truth', 'estimates', 'p_meta_data', 'ids' );
    end
end

if add_my_results_flag
    disp(['Adding new results to evaluation that were not submitted for the LOCATA challenge.'])
    disp(['Location of new results: ', n_dir])
    
    % Read participants' submissions:
    [ estimates, p_meta_data ] = read_participants_results( p_dir, eval_dir, valid_submission_ids, opts ); % participants' results and meta data
    
    disp('Not completed / tested yet!')
    keyboard
end

%% Get evaluation measures for each submission, task, recording, array

% Create dir for saving eval results:
if flag_save_data && ~exist( results_dir, 'dir')
    mkdir( results_dir )
end

current_date = datestr(now, 'yyyy-mm-dd_HHMMSS');

for ID = ids %#ok - ids saved and loaded from load_data_fname  (see l.115)
    
    tasks = p_meta_data{ID}.tasks; %#ok - p_meta_data saved and loaded from load_data_fname  (see l.115)
    array_id = p_meta_data{ID}.array_id;
    name =  p_meta_data{ID}.name; % name for submission ID
    
    fprintf('\n Processing %s : ID %d\n', name, ID);
    
    array_names = {'Robot head', 'DICIT', 'Hearing aids', 'Eigenmike'}; % names for the plots
    
    % Runs over the tasks available for each submission
    % NOTE: This does not necessarily run from 1, ..., 6 as most
    % participants did not submit to all tasks.
    for task_cnt = 1:length(tasks)
        % Create directory for this task in results directory:
        task_idx = tasks{task_cnt};
        
        if flag_save_data 
            results_task_dir = [results_dir, filesep, name,'_ID_',num2str(ID), filesep ,'task', num2str(task_idx)];
            if ~exist(results_task_dir, 'dir')
                mkdir(results_task_dir)
            end
        end
        
        % Open log file:
        fid = fopen([log_fname_root, '_task', num2str(task_idx), '.txt'], 'a');
        fprintf(fid, '\nSubmission ID: %s\n', num2str(ID));
        
        recordings = p_meta_data{ID}.recordings{task_cnt};
        
        for rec_cnt = 1:length(recordings)
            % Create directory for this recording in results directory:
            this_recording = recordings(rec_cnt);
            rec_idx = this_recording;
            
            if flag_save_data
                if ~exist([results_task_dir, filesep, 'recording', num2str(this_recording)], 'dir')
                    mkdir([results_task_dir, filesep, 'recording', num2str(this_recording)])
                end
            end
            
            for arr_cnt = 1:length(array_id)
                arr_idx = array_id( arr_cnt ); % NOT necessarily equal to arr_cnt!!!
                arr_idx_tr = arr_idx;
                    
                if ~isempty( estimates{ID, arr_idx, rec_idx, task_idx} ) && ~isempty(truth{ arr_idx_tr, rec_idx, task_idx})     %#ok - estimates & truth saved and loaded from load_data_fname  (see l.115)
                    disp(['Task ', num2str(task_idx), ', Recording ', num2str(rec_idx), ', Array ', array_names{arr_idx}, ' Name: ', name])
                    
                    % Create directory for this array in results directory:
                    this_array = array_names{arr_idx};
                    if ~isempty(regexp(this_array, 'Robot'))
                        this_array = 'Robot_head';
                    end
                    if flag_save_data
                        this_save_dir = [results_task_dir, filesep, 'recording', num2str(this_recording), filesep, this_array, filesep];
                        if ~exist(this_save_dir, 'dir')
                            mkdir(this_save_dir)
                        end
                    end
                    if flag_save_data
                        this_save_dir = [results_task_dir, filesep, 'recording', num2str(this_recording), filesep, this_array, filesep];
                        if ~exist(this_save_dir, 'dir')
                            mkdir(this_save_dir)
                        end
                    end
                    
                    % Run measures:
                    all_measures{ID,arr_idx,rec_idx,task_idx} = measures(truth{ arr_idx_tr, rec_idx, task_idx}, squeeze([estimates{ID, arr_idx, rec_idx, task_idx}]));
                    
                    %% Save results of measures to log:
                    for this_src_idx = 1 : length(all_measures{ID,arr_idx,rec_idx,task_idx}.mean_az_error)
                        az_err = rad2deg(all_measures{ID,arr_idx,rec_idx,task_idx}.mean_az_error{this_src_idx});
                        fprintf(fid, 'Azimuth error source %s: %f\n', num2str(this_src_idx), az_err);
                        
                        if isfield( all_measures{ID,arr_idx,rec_idx,task_idx}, 'mean_el_error' )
                            el_err = rad2deg(all_measures{ID,arr_idx,rec_idx,task_idx}.mean_el_error{1});
                            fprintf(fid, 'Elevation error source %s: %f\n', num2str(this_src_idx), el_err);
                        end
                    end
                    
                    N_swaps_per_src = all_measures{ID,arr_idx,rec_idx,task_idx}.N_swaps_per_src;
                    fprintf(fid, 'Number of swaps: %f\n', N_swaps_per_src);
                    
                    FAR = all_measures{ID,arr_idx,rec_idx,task_idx}.FAR;
                    fprintf(fid, 'FAR: %f\n', FAR);
                    
                    pd = all_measures{ID,arr_idx,rec_idx,task_idx}.pd{1:end};
                    fprintf(fid, 'pd: %f\n', pd);
                    
                    TFR = all_measures{ID,arr_idx,rec_idx,task_idx}.track_frag_rate;
                    fprintf(fid, 'TFR: %f\n', TFR);
                    
                    TL_ins = all_measures{ID,arr_idx,rec_idx,task_idx}.track_latency;
                    fprintf(fid, 'Track latency [s]: %f\n', TL_ins);
                    
                    if task_idx == 1 || task_idx == 3 || task_idx == 5
                        az_error_total =  rad2deg(all_measures{ID,arr_idx,rec_idx,task_idx}.mean_az_error_ss{:});
                        fprintf(fid, 'Single-source azimuth error: %f\n', az_error_total);
                        
                        if isfield( all_measures{ID,arr_idx,rec_idx,task_idx}, 'mean_el_error_ss')
                            el_error_total = rad2deg(all_measures{ID,arr_idx,rec_idx,task_idx}.mean_el_error_ss{:});
                            fprintf(fid, 'Single-source elevation error: %f\n', el_error_total);
                        end
                    end
                    
                    figure(all_measures{ID,arr_idx,rec_idx,task_idx}.f_ref); suptitle(['Task ', num2str(task_idx), ', Rec ', num2str(rec_idx), ', ', array_names{arr_idx}, ', ', name])
                    if flag_save_data
                        saveas(all_measures{ID,arr_idx,rec_idx,task_idx}.f_ref, [this_save_dir filesep 'measures.fig'], 'fig');
                        saveas(all_measures{ID,arr_idx,rec_idx,task_idx}.f_ref, [this_save_dir filesep 'measures.png'], 'png');
                        
                        for fig_idx = 1 : length(all_measures{ID,arr_idx,rec_idx,task_idx}.f_presentation)
                            saveas(all_measures{ID,arr_idx,rec_idx,task_idx}.f_presentation(fig_idx), [this_save_dir filesep 'measures_presentation_', all_measures{ID,arr_idx,rec_idx,task_idx}.f_label{fig_idx}, '.fig'], 'fig');
                            saveas(all_measures{ID,arr_idx,rec_idx,task_idx}.f_presentation(fig_idx), [this_save_dir filesep 'measures_presentation_', all_measures{ID,arr_idx,rec_idx,task_idx}.f_label{fig_idx}, '.png'], 'png');
                        end
                    end
                    
                    fprintf(fid, '*******************************\n');
                elseif ids ~= baseline_ID
                    % should not occur for participants results
                    warning('Empty cell detected')
                    
                end
                % pause
                close all;
            end
        end
    end
    % close log file:
    fclose(fid);
    
    % Save the results - overwrite on every iteration
    if flag_save_data
        save( [results_dir, filesep, 'measures_', current_date, '.mat'], 'all_measures');
    end
end % eof ID loop

%% Evaluate average measures per submission, task and array

% Measures - Extract statistics over multiple recordings
avg_measures = measures_statistics(all_measures);           %#ok - Variable used for saving to mat file
save([ results_dir , filesep, 'averaged_measures_', current_date, '.mat'], 'avg_measures', 'p_meta_data', 'all_measures' );

%% Plots

% Plot figures
figures_for_paper(results_dir, current_date);

% Print numbers on screen:
avg_measures = stats_for_paper(results_dir, current_date, ids);

%% Edit plot for Task 3, Recording 2, Submissions 6 + 7:
fig_327_em = openfig([results_dir, '/Salvati_ID_7/task3/recording2/Eigenmike/measures_presentation_az.fig']);
obj_327_em_markers = fig_327_em.Children(2).Children(2);
obj_327_em_dots = fig_327_em.Children(2).Children(3);
pause(.1)
fig_326_em = openfig([results_dir, '/Lebarbenchon_ID_6/task3/recording2/Eigenmike/measures_presentation_az.fig']);
hold on; 
plot(obj_327_em_markers.XData, obj_327_em_markers.YData, 'b', 'Marker', 'o', 'MarkerSize', 6, 'linestyle', 'none', 'linewidth', 2 )
plot(obj_327_em_dots.XData, obj_327_em_dots.YData, 'b', 'Marker', '.', 'MarkerSize', 6', 'linestyle', 'none', 'linewidth', 2)
fig_326_em.Children(2).Children(2).DisplayName = 'Submission 7';
fig_326_em.Children(2).Children(4).DisplayName = 'Submission 6';
fig_326_em.Children(2).Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
saveas(fig_326_em, [results_dir filesep 'az_task3_rec2_sub6-7_em.fig'], 'fig');
saveas(fig_326_em, [results_dir filesep 'az_task3_rec2_sub6-7_em.png'], 'png');

%% Edit plot for Task 3, Recording 2, Submissions 6 + 7:
fig_347_robot = openfig([results_dir, '/Salvati_ID_7/task3/recording4/Robot_head/measures_presentation_az.fig']);
fig_346_em = openfig([results_dir, '/Lebarbenchon_ID_6/task3/recording4/Eigenmike/measures_presentation_az.fig']);
fig_343_dicit = openfig([results_dir, '/Qian_ID_3/task3/recording4/DICIT/measures_presentation_az.fig']);

obj_347_robot_truth = fig_347_robot.Children(2).Children(4);
obj_347_robot_marker = fig_347_robot.Children(2).Children(2);
obj_347_robot_dots = fig_347_robot.Children(2).Children(3);
obj_346_em_truth = fig_346_em.Children(2).Children(4);
obj_346_em_marker = fig_346_em.Children(2).Children(2);
obj_346_em_dots = fig_346_em.Children(2).Children(3);

figure(fig_343_dicit); hold on;
fig_343_dicit.Children(2).Children(2).DisplayName = 'Submission 3, DICIT';
fig_343_dicit.Children(2).Children(4).DisplayName = 'Ground truth, DICIT';
fig_343_dicit.Children(2).Children(4).Color = fig_343_dicit.Children(2).Children(2).Color;

plot(obj_347_robot_marker.XData, obj_347_robot_marker.YData, 'b', 'Marker', 'o', 'MarkerSize', 6, 'linestyle', 'none', 'linewidth', 2, 'linejoin', 'round', 'DisplayName', 'Submission 7, Robot Head' )
plot(obj_347_robot_dots.XData, obj_347_robot_dots.YData, 'b', 'Marker', '.', 'MarkerSize', 6', 'linestyle', 'none', 'linewidth', 2)
plot(obj_347_robot_truth.XData, obj_347_robot_truth.YData, 'b', 'linewidth', 2, 'DisplayName', 'Ground truth, Robot Head')

plot(obj_346_em_marker.XData, obj_346_em_marker.YData, 'color', [119 172 48]/255, 'Marker', 'o', 'MarkerSize', 6, 'linestyle', 'none', 'linewidth', 2, 'linejoin', 'round', 'DisplayName', 'Submission 6, Eigenmike' )
plot(obj_346_em_dots.XData, obj_346_em_dots.YData, 'color', [119 172 48]/255, 'Marker', '.', 'MarkerSize', 6', 'linestyle', 'none', 'linewidth', 2)
plot(obj_346_em_truth.XData, obj_346_em_truth.YData, 'color', [119 172 48]/255, 'linewidth', 2, 'DisplayName', 'Ground truth, Eigenmike')

fig_343_dicit.Children(2).YLim = [-50 50];
fig_343_dicit.OuterPosition(4) = 1.75 * fig_343_dicit.OuterPosition(4);
fig_343_dicit.Children(2).Children(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
fig_343_dicit.Children(2).Children(5).Annotation.LegendInformation.IconDisplayStyle = 'off';

saveas(fig_343_dicit, [results_dir filesep 'az_task3_rec4_sub3-6-7.fig'], 'fig');
saveas(fig_343_dicit, [results_dir filesep 'az_task3_rec4_sub3-6-7.png'], 'png');

%% Edit plot for Fig 10 & Fig 11 (remove VAD from azimuth plots for readability)

fname = 'Ban_ID_4/task2/recording5/Robot_head/measures_presentation_az';
rm_VAD_legend(results_dir, fname)

fname = 'Ban_ID_4/task2/recording5/Robot_head/measures_presentation_VAD';
rm_VAD_legend(results_dir, fname, 1)

fname = '/Ban_ID_4/task4/recording4/Robot_head/measures_presentation_az';
rm_VAD_legend(results_dir, fname)

fname = 'Ban_ID_4/task4/recording4/Robot_head/measures_presentation_VAD';
rm_VAD_legend(results_dir, fname, 1)

fname = 'Ban_ID_4/task6/recording2/Robot_head/measures_presentation_az';
rm_VAD_legend(results_dir, fname)

fname = 'Ban_ID_4/task6/recording2/Robot_head/measures_presentation_VAD';
rm_VAD_legend(results_dir, fname, 1)

% Fig 11:
fname = 'Ban_ID_4/task4/recording1/Robot_head/measures_presentation_az';
rm_VAD_legend(results_dir, fname)

fname = 'Liu_ID_2/task4/recording1/Eigenmike/measures_presentation_az';
rm_VAD_legend(results_dir, fname)

fname = 'Ban_ID_4/task4/recording1/Robot_head/measures_presentation_OSPA';
rm_VAD_legend(results_dir, fname)

fname = 'Liu_ID_2/task4/recording1/Eigenmike/measures_presentation_OSPA';
rm_VAD_legend(results_dir, fname)

fname = 'Liu_ID_2/task4/recording1/Robot_head/measures_presentation_az';
rm_VAD_legend(results_dir, fname)

fname = 'Kitic_ID_10/task4/recording1/Eigenmike/measures_presentation_az';
rm_VAD_legend(results_dir, fname)

fname = 'Liu_ID_2/task4/recording1/Robot_head/measures_presentation_OSPA';
rm_VAD_legend(results_dir, fname)

fname = 'Kitic_ID_10/task4/recording1/Eigenmike/measures_presentation_OSPA';
rm_VAD_legend(results_dir, fname)

fname = 'Ban_ID_4/task4/recording1/Robot_head/measures_presentation_VAD';
rm_VAD_legend(results_dir, fname, 1)

fname = 'Liu_ID_2/task4/recording1/Eigenmike/measures_presentation_VAD';
rm_VAD_legend(results_dir, fname, 1)

return



