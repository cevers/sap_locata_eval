function out = measures_statistics(all_measures)

% function measures(truth, estimates)
% Evaluation measures for the LOCATA challenge
%
% Inputs:
%   all_measures: cell array of measures obtained for each algorithm,
%   array, recording and task
%
% Outputs: 
%   out: Structure of average statistics (mean and variance) over all
%   measures (az & el error, pd, FAR, track latency, track fragmentation
%   rate)
%   out.az_error      : azimuth error
%   out.az_error_var  : azimuth variance
%   out.el_error      : elevation error
%   out.el_error_var  : elevation variance
%   out.pd            : probablility of detection
%   out.pd_var        : PD variance
%   out.track_latency : track latency
%   out.track_latency_var : variance of track latency
%   out.track_frag_rate   : track fragmentation rate
%   out.track_frag_rate_var : variance of track fragmentation rate
%   Each field contains a cell array of size no. of tasks (6) X no. of
%   arrays (4) X no. of submissions (16) for the fields
%
% Author: Christine Evers, c.evers@imperial.ac.uk
%
% Reference: LOCATA documentation for participants (v3)
%            www.locata-challenge.org
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

out.az_error = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.az_error_var = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.el_error = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.el_error_var = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.az_error_ss = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.az_error_ss_var = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.pd = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.pd_var = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.track_latency = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.track_latency_var = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.track_frag_rate = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.track_frag_rate_var = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.track_swap_rate = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.track_swap_rate_var = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.FAR = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
out.FAR_VAP = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));

% Number of OSPA p-values:
% Submission 4, robot head, recording 4, task 6 definitely contains OSPA estimates. Extract number of p-hypotheses from there
num_OSPA_p = size(all_measures{4,1,4,6}.OSPA,2);        
out.OSPA = cell(num_OSPA_p,1);
out.OSPA_var = cell(num_OSPA_p,1);
for p_idx = 1 : num_OSPA_p
    out.OSPA{p_idx} = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
    out.OSPA_var{p_idx} = nan(size(all_measures,4), size(all_measures,2), size(all_measures,1));
end

for task_idx = 1:size(all_measures,4)
    for arr_idx = 1:size(all_measures,2)
        for alg_idx = 1 : size(all_measures,1)
            az_error = [];
            el_error = [];
            az_error_ss = [];
            pd = [];
            FAR = [];
            FAR_VAP = [];
            track_latency = [];
            track_frag_rate = [];
            track_swap_rate = [];
            OSPA = cell(num_OSPA_p,1);
            for rec_idx = 1:size(all_measures,3)
                if ~isempty( all_measures{alg_idx,arr_idx,rec_idx,task_idx} )
                    % FAR is independent of source index:
                    FAR = [FAR all_measures{alg_idx,arr_idx,rec_idx,task_idx}.FAR];
                    FAR_VAP = [FAR_VAP all_measures{alg_idx,arr_idx,rec_idx,task_idx}.FAR_VAP];
                        
                    % OSPA:
                    for p_idx = 1 : num_OSPA_p
                        this_OSPA = all_measures{alg_idx,arr_idx,rec_idx,task_idx}.OSPA(:,p_idx)';
                        OSPA{p_idx} = [OSPA{p_idx} this_OSPA(~isnan(this_OSPA))];
                    end
                    
                    % All measures that depend on source index:
                    for src_idx = 1 : length(all_measures{alg_idx,arr_idx,rec_idx,task_idx}.mean_az_error)
                        az_error = [az_error all_measures{alg_idx,arr_idx,rec_idx,task_idx}.mean_az_error{src_idx}(~isnan(all_measures{alg_idx,arr_idx,rec_idx,task_idx}.mean_az_error{src_idx}))];
                        if sum(strcmp(fields(all_measures{alg_idx,arr_idx,rec_idx,task_idx}), 'mean_el_error')) > 0
                            el_error = [el_error all_measures{alg_idx,arr_idx,rec_idx,task_idx}.mean_el_error{src_idx}(~isnan(all_measures{alg_idx,arr_idx,rec_idx,task_idx}.mean_el_error{src_idx}))];
                        end
                        az_error_ss = [az_error_ss all_measures{alg_idx,arr_idx,rec_idx,task_idx}.mean_az_error_ss{src_idx}(~isnan(all_measures{alg_idx,arr_idx,rec_idx,task_idx}.mean_az_error_ss{src_idx}))];
                        pd = [pd all_measures{alg_idx,arr_idx,rec_idx,task_idx}.pd{src_idx}(~isnan(all_measures{alg_idx,arr_idx,rec_idx,task_idx}.pd{src_idx}))];
                        track_latency = [track_latency all_measures{alg_idx,arr_idx,rec_idx,task_idx}.track_latency(src_idx)];
                        track_frag_rate = [track_frag_rate all_measures{alg_idx,arr_idx,rec_idx,task_idx}.track_frag_rate(src_idx)];
                        track_swap_rate = [track_swap_rate all_measures{alg_idx,arr_idx,rec_idx,task_idx}.track_swap_rate(src_idx)];
                    end
                end
            end
            
            if ~isnan(az_error)
                out.az_error(task_idx, arr_idx, alg_idx) = mean(az_error);
                out.az_error_var(task_idx,arr_idx,alg_idx) = var(az_error);
                if ~isempty(el_error)
                    out.el_error(task_idx, arr_idx, alg_idx) = mean(el_error);
                    out.el_error_var(task_idx,arr_idx,alg_idx) = var(el_error);
                end
                out.az_error_ss(task_idx, arr_idx, alg_idx) = mean(az_error_ss);
                out.az_error_ss_var(task_idx,arr_idx,alg_idx) = var(az_error_ss);
                out.pd(task_idx, arr_idx, alg_idx) = mean(pd);
                out.pd_var(task_idx, arr_idx, alg_idx) = var(pd);
                out.FAR(task_idx, arr_idx, alg_idx) = mean(FAR);
                out.FAR_var(task_idx, arr_idx, alg_idx) = var(FAR);
                out.FAR_VAP(task_idx, arr_idx, alg_idx) = mean(FAR_VAP);
                out.FAR_VAP_var(task_idx, arr_idx, alg_idx) = var(FAR_VAP);
                out.track_latency(task_idx, arr_idx, alg_idx) = mean(track_latency);
                out.track_latency_var(task_idx, arr_idx, alg_idx) = var(track_latency);
                out.track_frag_rate(task_idx, arr_idx, alg_idx) = mean(track_frag_rate);
                out.track_frag_rate_var(task_idx, arr_idx, alg_idx) = var(track_frag_rate);
                out.track_swap_rate(task_idx, arr_idx, alg_idx) = mean(track_swap_rate);
                out.track_swap_rate_var(task_idx, arr_idx, alg_idx) = var(track_swap_rate);
                
                for p_idx = 1:num_OSPA_p
                    out.OSPA{p_idx}(task_idx,arr_idx,alg_idx) = mean(OSPA{p_idx});
                    out.OSPA_var{p_idx}(task_idx,arr_idx,alg_idx) = var(OSPA{p_idx});
                end
            end
        end
    end
end
