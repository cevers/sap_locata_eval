function figures_for_paper(results_dir, results_date)

% Script to display LOCATA results for the presentation
%
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

close all;

load( [results_dir, filesep, 'averaged_measures_', results_date, '.mat'])

figure_dir = [results_dir, filesep, 'stats/'];
if ~exist(figure_dir, 'dir')
    mkdir(figure_dir)
    disp(['Created ', figure_dir])
end

flag_final = 1;
fsize = 20;   % font size
pos = [100,100,1250,350]; % position and size
flag_title = 0; % if 1, title for each figure
flag_save = 1;  % if 1, figures are saved
if flag_save && ~exist( figure_dir, 'dir' )
    mkdir( figure_dir )
end
array_names = {'Robot Head', 'DICIT array', 'Hearing Aids', 'Eigenmike'}; 
% same ordering as for array index!
if ~flag_final
    flag_title = 0;
    fsize = 28;   % font size
end

az_limits = rad2deg([min(avg_measures.az_error(:)) 1.1*max(avg_measures.az_error(:))]);
az_limits_noGating = rad2deg([min(avg_measures.az_error_ss(:)) 1.1*max(avg_measures.az_error_ss(:))]);
el_limits = rad2deg([min(avg_measures.el_error(:)) 1.1*max(avg_measures.el_error(:))]);
pd_limits = [0 110];
TFR_limits = [min(avg_measures.track_frag_rate(:)) 1.1*max(avg_measures.track_frag_rate(:))];
TSR_limits = [min(avg_measures.track_swap_rate(:)) 1.1*max(avg_measures.track_swap_rate(:))];
TL_limits = [min(avg_measures.track_latency(:)) 1.1*max(avg_measures.track_latency(:))];
FAR_limits = [min(avg_measures.FAR(:)) 1.1*max(avg_measures.FAR(:))];
FAR_VAP_limits = [min(avg_measures.FAR_VAP(:)) 1.1*max(avg_measures.FAR_VAP(:))];

min_OSPA = inf;
max_OSPA = -inf;
for p_idx = 1 : length(avg_measures.OSPA)
    if min_OSPA > min(avg_measures.OSPA{p_idx}(:))
        min_OSPA = min(avg_measures.OSPA{p_idx}(:));
    end
    if max_OSPA < max(avg_measures.OSPA{p_idx}(:))
        max_OSPA = max(avg_measures.OSPA{p_idx}(:));
    end
end     
OSPA_limits = [min_OSPA max_OSPA];

for this_task = 1:6
    %% Azimuth
    
    % With gating
    generate_figure(this_task, rad2deg(avg_measures.az_error), rad2deg(sqrt(avg_measures.az_error_var)), 'Azimuth error [deg]', az_limits, array_names, fsize, pos, figure_dir, 'azimuth', flag_final, flag_title, flag_save)
    
    % Without gating (only for single-source tasks 1, 3, 5)

    if this_task==1 || this_task==3 || this_task==5
        generate_figure(this_task, rad2deg(avg_measures.az_error_ss), rad2deg(sqrt(avg_measures.az_error_ss_var)), 'Azimuth error [deg]', az_limits_noGating, array_names, fsize, pos, figure_dir, 'azimuth_noGating', flag_final, flag_title, flag_save)
    end
    
    %% Elevation
    
    generate_figure(this_task, rad2deg(avg_measures.el_error), rad2deg(sqrt(avg_measures.el_error_var)), 'Elevation error [deg]', el_limits, array_names, fsize, pos, figure_dir, 'elevation', flag_final, flag_title, flag_save)
    
    %% Pd
    
    generate_figure(this_task, avg_measures.pd*100, avg_measures.pd_var*100, '$p_d$ [\%]', pd_limits, array_names, fsize, pos, figure_dir, 'pd', flag_final, flag_title, flag_save)
    
    %% FAR
    
    % FAR during recording:
    generate_figure(this_task, avg_measures.FAR, avg_measures.FAR_var, 'False Alarm Rate', FAR_limits, array_names, fsize, pos, figure_dir, 'FAR', flag_final, flag_title, flag_save)

    % FAR during voice-activity periods only (added for T-ASL paper)
    generate_figure(this_task, avg_measures.FAR_VAP, avg_measures.FAR_VAP_var, 'False Alarm Rate', FAR_VAP_limits, array_names, fsize, pos, figure_dir, 'FAR_VAP', flag_final, flag_title, flag_save)
    
    %% TFR
    
    generate_figure(this_task, avg_measures.track_frag_rate, avg_measures.track_frag_rate_var, 'Track Fragm. Rate', TFR_limits, array_names, fsize, pos, figure_dir, 'TFR', flag_final, flag_title, flag_save)
    
    %% TSR
    
    generate_figure(this_task, avg_measures.track_swap_rate, avg_measures.track_swap_rate_var, 'Track Swap Rate', TSR_limits, array_names, fsize, pos, figure_dir, 'TSR', flag_final, flag_title, flag_save)
        
    %% Track Latency
    
    generate_figure(this_task, avg_measures.track_latency, avg_measures.track_latency_var, 'Track Latency [s]', TL_limits, array_names, fsize, pos, figure_dir, 'TL', flag_final, flag_title, flag_save)

    %% OSPA
    
    for p_idx = 1:length(avg_measures.OSPA)
        generate_figure(this_task, avg_measures.OSPA{p_idx}, avg_measures.OSPA_var{p_idx}, 'OSPA', OSPA_limits, array_names, fsize, pos, figure_dir, ['OSPA_pidx',num2str(p_idx)], flag_final, flag_title, flag_save)
    end

    close all;
end % eof task loop

function generate_figure(this_task, this_measure, this_measure_var, this_ylabel, ylimits, array_names, fsize, fpos, figure_dir, save_str, flag_final, flag_title, flag_save)

num_arrays = size(this_measure,2);
num_submissions = size(this_measure,3);

if flag_final
    cf = figure; hold on;
    clf
    fig_title = strcat( 'Task ', num2str(this_task),', All arrays');
else
    subplot(N,2,4)
end
set(gcf, 'position', fpos )
set( cf, 'Name', fig_title )

% Submissions:
submissions_for_task = sum(~isnan(squeeze(this_measure(this_task,:,:))),1) > 0;

% Bar plot
hBar=bar(squeeze(this_measure(this_task,:,:)),1,'grouped');      % basic grouped bar plot, keep handle
hFg=gcf; hAx=gca;             % handles to current figure, axes
ylim( ylimits )
hAx.XTickLabel=array_names; % label x axis
% hAx.YColor=get(hFg,'color'); 
%box off  % hide y axis/labels, outline

% Change colour map:
cmap = jet(num_submissions-1);
cmap(end+1,:) = [0 0 0];
for bar_idx = 1 : length(hBar)
    set(hBar(bar_idx), 'FaceColor', cmap(bar_idx,:))
end

% Whiskers / error bars:
x=repmat([1:num_arrays].',1,num_submissions);  % bar x nominals
xOff=bsxfun(@plus,x,[hBar.XOffset]); % bar midpoints x
hold all                             % keep axis limits, etc.
hEB=errorbar(xOff,squeeze(this_measure(this_task,:,:)),squeeze(this_measure_var(this_task,:,:)),'.');   % add errorbars

% Axes labels
ylabel(this_ylabel, 'interpreter', 'latex', 'FontSize',fsize)
ylim( ylimits )
xticklabels( array_names )
set(gca,'TickLabelInterpreter','latex', 'FontSize', fsize)

% Label submissions as text:
hT=[];              % placeholder for text object handles
for i=1:length(hBar)  % iterate over number of bar objects
    if i == 17
        hT = [hT,text(hBar(i).XData+hBar(i).XOffset, 1.05*squeeze(this_measure(this_task,:,i)), 'BL', ...
            'VerticalAlignment','bottom','horizontalalign','center', 'interpreter', 'latex', 'rotation', 0, 'fontsize', fsize-6)];
    else
        hT = [hT,text(hBar(i).XData+hBar(i).XOffset, 1.05*squeeze(this_measure(this_task,:,i)), num2str(i), ...
            'VerticalAlignment','bottom','horizontalalign','center', 'interpreter', 'latex', 'rotation', 0, 'fontsize', fsize-6)];
    end
end

% Draw lines between array groups:
x_baseline_perArray = hBar(num_submissions).XData+hBar(num_submissions).XOffset;
x_subm1_perArray = hBar(1).XData+hBar(1).XOffset;
dist_betwArrays = x_subm1_perArray(2:end)-x_baseline_perArray(1:3);
for array_idx = 1 : (num_arrays-1)
    line(x_baseline_perArray(array_idx)+dist_betwArrays(array_idx)/2*ones(1,2), ylimits, 'color', [0.5 0.5 0.5]);
end

% Horizontal grid only:
ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

% Legend of submission IDs in this task:
legend_labels = num2str(find(submissions_for_task)');
for legend_els = 1 : length(legend_labels)
    if legend_labels(legend_els,:) == '17'
        legend_labels(legend_els,:) = 'BL';
        break
    end
end
hLg=legend(hBar(submissions_for_task), legend_labels, ...
           'orientation','horizontal', ...
           'location','northoutside'); legend('boxoff') % add legend
set(hLg,'interpreter', 'latex')
       
if flag_title
    title( strcat( 'Task ', num2str(this_task)), 'interpreter', 'latex' )
end

if flag_save
    saveas(gcf, [figure_dir, filesep,  save_str, '_task',num2str(this_task),'.fig'], 'fig');
    saveas(gcf, [figure_dir, filesep, save_str, '_task',num2str(this_task),'.png'], 'png');
end
    
return