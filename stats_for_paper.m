function avg_measures = stats_for_paper(results_dir, results_date, ids)

load([results_dir, filesep, 'averaged_measures_', results_date, '.mat'])

save_dir = [results_dir, '/stats/'];
if ~exist(save_dir)
    mkdir(save_dir)
    disp(['Created ', save_dir])
end

num_tasks = 6;

algs_135 = find(sum(squeeze(any(avg_measures.az_error([1,3,5],:,:),2)),1)==3);
algs_135(find(algs_135==17)) = [];

algs_246 = find(sum(squeeze(any(avg_measures.az_error([2,4,6],:,:),2)),1)==3);
algs_246(find(algs_246==17)) = [];

algs_24 = find(sum(squeeze(any(avg_measures.az_error([2,4],:,:),2)),1)==2);
algs_24(find(algs_24==17)) = [];

algs_multisource = find(sum(squeeze(any(avg_measures.az_error([2,4,6],:,:),2)),1)>=1);

disp(['Algorithms in $S_{1,3,5}$: ', num2str(algs_135)])
disp(['Algorithms in $S_{2,4,6}$: ', num2str(algs_246)])

%% Tasks 1, 3, 5

fprintf('******************************************\n')
fprintf('Averages of submissions to all three tasks 1, 3, 5\n')
fprintf(['Results evaluated for submission IDs: ', num2str(algs_135),'\n\n'])

%% Averages across all submissions submitted to specified tasks:

% Azimuth for single-source tasks 1, 3, 5:
mean_across_submissions( rad2deg(avg_measures.az_error_ss), algs_135, 1:2:num_tasks, 'azimuth error', 'deg' )

% Elevation for single-source tasks 1, 3, 5:
mean_across_submissions( rad2deg(avg_measures.el_error), algs_135, 1:2:num_tasks, 'elevation error', 'deg' )

% Probability of detection:
mean_across_submissions( avg_measures.pd*100, algs_135, 1:2:num_tasks, 'probability of detection', 'percent' )

% False alarm rate (across recording):
mean_across_submissions( avg_measures.FAR, algs_135, 1:2:num_tasks, 'FAR throughout recording', '' )

% False alarm rate (voice-activity only):
mean_across_submissions( avg_measures.FAR_VAP, algs_135, 1:2:num_tasks, 'FAR during voice-activity only', '' )

%% FAR for tasks 1, 3 and for the eigenmike

FAR_em_task13 = squeeze(avg_measures.FAR_VAP([1,3],4,:));
% Keep only submissions that were submitted for both tasks 1 and 3:
FAR_em_task13 = FAR_em_task13(:,~isnan(sum(FAR_em_task13,1)));
mean_FAR_em_task13 = mean(FAR_em_task13,2);
disp(['Average FAR during VAPs using eigenmike for task 1: ', num2str(mean_FAR_em_task13(1)), ', task 3: ', num2str(mean_FAR_em_task13(2))])

%% Average FAR per task:

mean_TL = zeros(size(avg_measures,1),1);
for task_idx = 1:2:size(avg_measures.track_latency,1)
    TL = squeeze(avg_measures.track_latency(task_idx,:,:));
    TL = TL(:,algs_135);
    mean_TL(task_idx) = mean(TL(~isnan(TL)));
end
fprintf(['Average TL over VAPs for tasks 1, 3, 5: ', mat2str(round(mean_TL([1,3,5]),1)), '.\n'])

%% Tables:

disp('******************************************')
disp(' ')

% Sort tables rows by single- (1,3,5) and multi-source (2,4,6) tasks:
task_135_246_idx = [1, 3, 5, 2, 4, 6];
task_135 = [1,3,5];
task_246 = [2,4,6];

% Azimuth error
save_latex_table(rad2deg(avg_measures.az_error(task_135_246_idx,:,ids)), [save_dir, 'az_error.tex'], 1);

% FAR during VAPs:
save_latex_table(avg_measures.FAR_VAP(task_135_246_idx,:,ids), [save_dir, 'FAR_voiceActivity.tex'], 0);

% FAR throughout recording:
save_latex_table(avg_measures.FAR(task_135_246_idx,:,ids), [save_dir, 'FAR_recording.tex'], 0);

% Pd:
save_latex_table(avg_measures.pd(task_135_246_idx,:,ids)*100, [save_dir, 'pd.tex'], 0);

% Track swap rate:
save_latex_table(avg_measures.track_swap_rate(task_135_246_idx,:,ids), [save_dir, 'track_swap_rate.tex'], 0);

% Track fragmentation rate:
save_latex_table(avg_measures.track_frag_rate(task_135_246_idx,:,ids), [save_dir, 'track_frag_rate.tex'], 0);

% Difference in azimuth error resulting from gating for single-source
% tasks:
diff = avg_measures.az_error - avg_measures.az_error_ss(:,:,:);
for task_idx = task_135
    offset_gating = squeeze(diff(task_idx,:,:));
    offset_inds = find(~isnan(offset_gating));
    nonoffset_inds = find(isnan(offset_gating));
    diff(task_idx,nonoffset_inds) = nan;
    [arr_inds, submission_inds] = ind2sub(size(offset_gating), offset_inds);
end
% diff = diff(1:2:end,:,:);
save_latex_table(abs(rad2deg(diff(task_135,:,ids))), [save_dir, 'az_error_diff_gating.tex'], 2);

% OSPA:
OSPA_246 = cell(size(avg_measures.OSPA));
for p_idx = 1 : length(OSPA_246)
    OSPA_246{p_idx} = avg_measures.OSPA{p_idx}(task_246,:,algs_multisource);
end
save_OSPA_table(OSPA_246([1,4]), [save_dir, 'OSPA.tex'], algs_multisource);

return

function save_latex_table( data, fname, az_table_flag )

if az_table_flag ==2
    caption_str = '\caption{Difference in average azimuth errors with and without gating, evaluated for single-source tasks 1, 3, 5 for all submissions and the baseline (BL). Submissions unaffected by gating, and hence outliers, are highlighted in bold font. Results marked by $\ast$ are inconclusive.}';
    label_str = '\label{table:azimuth_errors_gating_offset}';
elseif az_table_flag == 1
    caption_str = '\caption{Average azimuth errors during \ac{VAP}. Submissions corresponding to minimum average errors are highlighted in bold font. Column colour indicates type of algorithm, where white indicates frameworks involving only \acs{DoA} estimation (Submission IDs 1, 6, 9, 11, 12, 15, 16 and the baseline (BL)), and grey indicates frameworks that combine \ac{DoA} estimation with source tracking (Submission IDs 2, 3, 4, 7, 8, 10).}';
    label_str = '\label{table:azimuth_errors}';
elseif az_table_flag == 3
    caption_str = '\caption{Average \acs{OSPA} results. Submissions corresponding to minimum average errors are highlighted in bold font. Column colour indicates type of algorithm, where white indicates frameworks involving only \acs{DoA} estimation (Submission IDs 1, 6, 9, 11, 12, 15, 16 and the baseline (BL)), and grey indicates frameworks that combine \ac{DoA} estimation with source tracking (Submission IDs 2, 3, 4, 7, 8, 10).}';
    label_str = '\label{table:OSPA}';
else
    caption_str = '\caption{MyTableCaption}';
    label_str = '\label{table:MyTableLabel}';
end
    
[ntasks,narrays,nsubs] = size(data);
ntasks_ss = 3;
arrays = {'Robot Head', 'DICIT', 'Hearing Aids', 'Eigenmike'};
tasks = [1, 3, 5, 2, 4, 6];

% Add column (single/multi-source) and column 2 (array type) to data:
label_singlesource = '\multirow{12}{*}{\rotatebox[origin=c]{90}{Single Source}}';
label_multisource = '\multirow{12}{*}{\rotatebox[origin=c]{90}{Multiple Sources}}';
label_empty = ' & ';

table_in.data = permute(data,[2 1 3]);
table_in.data = round(reshape(table_in.data, size(table_in.data,1)*size(table_in.data,2),size(table_in.data,3)),1);
table_in.dataFormat = {'%.1f',size(data,3)};%{'%d',1,'%.3f',3,'%.0f',2,'%.3f',1}; % 1 col - int, following cols: float, then 2 cols int, last col float

% convert to latex table:
latex = latexTable(table_in);

% Edit table headers for journal:
[nrows,ncols] = size(latex);
idx = find(~cellfun(@isempty, regexp(latex, 'begin{tabular}')));
latex{idx} = '\begin{tabular}{|c|c|c||c|>{\columncolor{darkgrey}}c|>{\columncolor{darkgrey}}c|>{\columncolor{darkgrey}}c|c|>{\columncolor{darkgrey}}c|>{\columncolor{darkgrey}}c|c|>{\columncolor{darkgrey}}c|c|c|c|c|c|}';
latex_tmp = cell(size(latex,1)+2,1);
for row = 1 : nrows
    if row <= idx + 1
        latex_tmp{row} = latex{row};
    else
        latex_tmp{row+2} = latex{row};
    end
end
latex_tmp{idx+2} = '\multicolumn{2}{|c|}{\multirow{2}{*}{Task}} & \multirow{2}{*}{Array} & \multicolumn{14}{c|}{Submission ID}\\\cline{4-17}';
latex_tmp{idx+3} = '\multicolumn{2}{|c|}{} & & 1 & 2 & 3 & 4 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 15 & 16 & BL\\\hline\hline';
[nrows,ncols] = size(latex_tmp);

% Add two columns at beginning of each row (for single/multi-source label +
% array type)
end_idx = find(~cellfun(@isempty, regexp(latex_tmp, 'end{tabular}')));
array_ct = 1;
task_ct = 1;
for row = idx + 4 : end_idx-1
    if isempty(regexp(latex_tmp{row}, 'hline'))     % skip rows that correspond to 'hline'
        if array_ct > 1
            latex_tmp{row} = [label_empty, label_empty, arrays{array_ct}, label_empty, latex_tmp{row}];
        else
            latex_tmp{row} = [label_empty, ['\multirow{4}{*}{', num2str(tasks(task_ct)),'}'], label_empty, arrays{array_ct}, label_empty, latex_tmp{row}];
        end
        
        % Identify minimum error and bold the font:
        min_val = min(data(task_ct, array_ct,:));
        if mod(min_val, 1.0) == 0
            min_val = [num2str(min_val), '.0'];
        else
            min_val = round(min_val,1);
            if mod(min_val, 1.0) == 0
                min_val = [num2str(min_val), '.0'];
            else
                min_val = num2str(min_val);
            end
        end
        [min_srt, min_end] = regexp(latex_tmp{row}, min_val);
        bold_srt = '\textbf{';
        bold_end = '}';
        for min_idx = 1 : length(min_srt)
            latex_tmp{row} = [latex_tmp{row}(1:min_srt(min_idx)-1), bold_srt, latex_tmp{row}(min_srt(min_idx):min_end(min_idx)), bold_end, latex_tmp{row}(min_end(min_idx)+1:end)];
            
            % Compensate indices for addition of 'textbf'
            len_add_str = length([bold_srt, bold_end]);
            if length(min_srt) > 1
                min_srt = min_srt + len_add_str;
                min_end = min_end + len_add_str;
            end
            
        end
        
        % counter for array type:
        if array_ct < 4
            array_ct = array_ct + 1;
        else 
            array_ct = 1;
            task_ct = task_ct + 1;
        end
    else
        if array_ct > 1
            latex_tmp{row} = '\cline{3-17}';
        else
            latex_tmp{row} = '\cline{2-17}';
        end
    end
end
latex_tmp{row} = '\hline';

% Add label: single/multi-source:
latex_tmp{idx+4} = ['\multirow{12}{*}{\rotatebox[origin=c]{90}{Single Source}}', latex_tmp{idx+4}];
if ntasks == 6
    latex_tmp{idx+4 + narrays*ntasks_ss*2 - 1} = '\hline\hline';
    latex_tmp{idx+4 + narrays*ntasks_ss*2} = ['\multirow{12}{*}{\rotatebox[origin=c]{90}{Multiple Sources}}', latex_tmp{idx+4 + narrays*ntasks_ss*2}];
end


latex_tmp = latex_tmp(3:end-3);
latex = cell(length(latex_tmp)+6,1);
latex{1} = '\begin{table*}[tb]';
latex{2} = '\centering';
latex{3} = caption_str;
latex{4} = label_str;
latex(5 + (1:length(latex_tmp))) = latex_tmp;
latex{end} = '\end{table*}';

fid=fopen(fname,'w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);

return

function save_OSPA_table( data, fname, ids )

caption_str = '\caption{Average \acs{OSPA} results. Column colour indicates type of algorithm, where white indicates frameworks involving only \acs{DoA} estimation (Submission IDs 11, 12, 16 and the baseline (BL)), and grey indicates frameworks that combine \ac{DoA} estimation with source tracking (Submission IDs 2, 4, 10).}';
label_str = '\label{table:OSPA}';

[ntasks,narrays,nsubs] = size(data{1});
ntasks_ss = 3;
arrays = {'Robot Head', 'DICIT', 'Hearing Aids', 'Eigenmike'};
tasks = [2, 4, 6];
nps = length(data);

% Add column (single/multi-source) and column 2 (array type) to data:
% label_singlesource = '\multirow{12}{*}{\rotatebox[origin=c]{90}{Single Source}}';
label_multisource = '\multirow{12}{*}{\rotatebox[origin=c]{90}{Multiple Sources}}';
label_empty = ' & ';

% Permute data to fit into table:
table_in.data = nan(size(data{1},1), size(data{1},2), length(data)*size(data{1},3));
for p_idx = 1 : length(data)
    for task_idx = 1 : size(data{p_idx}, 1)
        for array_idx = 1 : size(data{p_idx}, 2)
            for submission_idx = 1 : size(data{p_idx}, 3)
                table_in.data(task_idx,array_idx, (submission_idx-1)*nps + p_idx) = data{p_idx}(task_idx,array_idx,submission_idx);
            end
        end
    end
end
table_in.data = permute(table_in.data,[2 1 3]);
table_in.data = round(reshape(table_in.data, size(table_in.data,1)*size(table_in.data,2),size(table_in.data,3)),1);
table_in.dataFormat = {'%.1f',size(table_in.data,2)};

% Convert to LaTeX table:
latex = latexTable(table_in);

%% Format columns of table

[nrows,ncols] = size(latex);
idx = find(~cellfun(@isempty, regexp(latex, 'begin{tabular}')));
latex{idx} = '\begin{tabular}{|c|c||>{\columncolor{darkgrey}}c|>{\columncolor{darkgrey}}c||>{\columncolor{darkgrey}}c|>{\columncolor{darkgrey}}c||>{\columncolor{darkgrey}}c|>{\columncolor{darkgrey}}c||c|c||c|c||c|c||c|c|}';
latex_tmp = cell(size(latex,1)+2,1);
for row = 1 : nrows
    if row <= idx + 1
        latex_tmp{row} = latex{row};
    else
        latex_tmp{row+2} = latex{row};
    end
end
latex_tmp{idx+2} = ['\multirow{4}{*}{\rotatebox[origin=c]{90}{Task}} & \multirow{4}{*}{Array} & \multicolumn{', num2str(nsubs*nps), '}{c|}{Submission ID}\\\cline{3-', num2str(nsubs*nps+2), '}'];
latex_tmp{idx+3} = [' & & \multicolumn{2}{|c|}{2} & \multicolumn{2}{|c|}{4} & \multicolumn{2}{|c|}{10} & \multicolumn{2}{|c|}{11} & \multicolumn{2}{|c|}{12} & \multicolumn{2}{|c|}{16} & \multicolumn{2}{|c|}{BL}\\\cline{3-', num2str(nsubs*nps+2), '}'];
[nrows,ncols] = size(latex_tmp);

% Add two columns at beginning of each row (for single/multi-source label +
% array type)
end_idx = find(~cellfun(@isempty, regexp(latex_tmp, 'end{tabular}')));
array_ct = 1;
task_ct = 1;
for row = idx + 4 : end_idx-1
    if isempty(regexp(latex_tmp{row}, 'hline'))     % skip rows that correspond to 'hline'
        if array_ct > 1
            latex_tmp{row} = [label_empty, arrays{array_ct}, label_empty, latex_tmp{row}];
        else
            latex_tmp{row} = [['\multirow{4}{*}{', num2str(tasks(task_ct)),'}'], label_empty, arrays{array_ct}, label_empty, latex_tmp{row}];
        end

        % counter for array type:
        if array_ct < 4
            array_ct = array_ct + 1;
        else 
            array_ct = 1;
            task_ct = task_ct + 1;
        end
    else
        if array_ct > 1
            latex_tmp{row} = ['\cline{2-', num2str(nsubs*nps+2), '}'];
        else
            latex_tmp{row} = ['\cline{1-', num2str(nsubs*nps+2), '}'];
        end
    end
end
latex_tmp{row} = '\hline';


%% Add additional rows

% New table conrtaining appropriate headers:
% Reshuffle to include each p_idx as a separate row in one "central" table:
latex_main = cell(size(latex_tmp,1)+1,1);
latex_tmp = latex_tmp(3:end-3);

latex_main{1} = '\begin{table*}[tb]';
latex_main{2} = '\centering';
latex_main{3} = caption_str;
latex_main{4} = label_str;
latex_main(5:8) = latex_tmp(1:4);
latex_main{9} = ['& & \multicolumn{2}{c|}{$p$} & \multicolumn{2}{c|}{$p$} & \multicolumn{2}{c|}{$p$} & \multicolumn{2}{c|}{$p$} & \multicolumn{2}{c|}{$p$} & \multicolumn{2}{c|}{$p$} & \multicolumn{2}{c|}{$p$} \\\cline{3-', num2str(nsubs*nps+2), '}'];
latex_main{10} = '& & $1$ & $5$ & $1$ & $5$ & $1$ & $5$ & $1$ & $5$ & $1$ & $5$ & $1$ & $5$ & $1$ & $5$ \\\hline\hline';
latex_main(11+(1:length(latex_tmp)-4)) = latex_tmp(5:end);
latex_main{end+1} = '\end{table*}';

latex_main

fid=fopen(fname,'w');
[nrows,ncols] = size(latex_main);
for row = 1:nrows
    fprintf(fid,'%s\n',latex_main{row,:});
end
fclose(fid);

return

function mean_across_submissions( data, subm_inds, task_inds, data_type_str, unit )

for task_idx = task_inds
    squeeze_data = squeeze(data(task_idx,:,:));
    squeeze_data = squeeze_data(:,subm_inds);
    mean_data(task_idx) = mean(squeeze_data(~isnan(squeeze_data)));
end
fprintf(['Average ', data_type_str, ' for tasks ', num2str(task_inds), ': ', mat2str(round(mean_data(task_inds),1)), ' ', unit, '.\n'])

return