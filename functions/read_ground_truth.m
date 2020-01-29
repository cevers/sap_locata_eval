function truth = read_ground_truth( data_dir, vad_dir )
% function truth = read_ground_truth( data_dir, vad_dir )
% provides ground-truth data for all tasks and recordings of a given LOCATA database
%
% INPUT:
% data_dir: dev or eval database with ground-truth data
% vad_dir:  corresponding database with ground-truth vad
%
% OUTPUT:
% truth: cell array with ground-truth data
%        truth{arr_idx, rec_idx, task_idx}.timestamps
%        truth{arr_idx, rec_idx, task_idx}.time
%        truth{arr_idx,rec_idx,task_idx}.source.(src_in_rec{src_idx}).VAD.activity
%        truth{arr_idx,rec_idx,task_idx}.source.(src_in_rec{src_idx}).VAD.num_periods
%
% authors: Heiner Loellmann (LMS, FAU) and Christine Evers (Imperial)
%
% Remark: based load_participants_results.m by C.Evers
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

% Initialize settings
opts = init();
fs = 48e3;   % audio sampling frequency
tasks = 1:6; % all tasks

for task_idx = 1:length(tasks)
    this_task = tasks(task_idx);
    task_dir = [data_dir, filesep, 'task', num2str(this_task)];
    finpath = dir(task_dir);
    
    % Read all recording IDs available for this task:
    recordings = nan( length(finpath),1);
    for f_idx = 1 : length(finpath)
        if finpath(f_idx).isdir
            srt_idx = regexp(finpath(f_idx).name, 'recording', 'end');
            recordings(f_idx) = str2double(finpath(f_idx).name((srt_idx+1):end));
        end
    end
    recordings = unique(recordings(~isnan(recordings)));
    
    % Parse through all recordings within this task:
    for rec_idx = 1:length(recordings)
        this_recording = recordings(rec_idx);
        rec_dir = [task_dir, filesep, 'recording', num2str(this_recording)];
        finpath = dir(rec_dir);
        
        % Read all recording IDs available for this task:
        array_names = {};
        for f_idx = 1 : length(finpath)
            if finpath(f_idx).isdir
                array_names{end+1} = finpath(f_idx).name;
            end
        end
        
        % Only need ground truth for one of the arrays in each recording needed:
        % arr_idx = 1;
        array_names = unique(intersect(array_names, opts.valid_arrays));
        
        if isstring(array_names)
            array_names = {array_names}; % convert string to cell array
        end
        
        for arr_cnt = 1:length( array_names )
            
            arr_idx = array_index( array_names(arr_cnt) ); % uniqe index for cell array with ground-truth needed (hence NOT equal to arr_cnt!!!)
            
            this_array = array_names{arr_cnt};  % current array name
            array_dir = [rec_dir filesep array_names{arr_cnt}];
            
            % Load data from csv / wav files in database:
            fprintf('\n Loading data for task %d, recording %d, %s \n', this_task, this_recording, this_array )
            % [ audio_array, audio_source, position_array, position_source, required_time] = load_data(array_dir);
            [ ~ , ~, position_array, position_source, required_time] = load_data(array_dir);
            truth{arr_idx, rec_idx, task_idx} = get_truth(this_array, position_array, position_source, required_time);
            
            % Load timestamps of recording:
            truth{arr_idx, rec_idx, task_idx}.timestamps = elapsed_time(required_time.time);
            truth{arr_idx, rec_idx, task_idx}.timestamps = truth{arr_idx, rec_idx, task_idx}.timestamps(find(required_time.valid_flag));
            truth{arr_idx, rec_idx, task_idx}.time = required_time.time(:,find(required_time.valid_flag));
            
            % Get names of sources that are active in this recording:
            src_in_rec = fieldnames(truth{arr_idx,rec_idx,task_idx}.source);
            NoS = length(src_in_rec); % No. of sources
            
            % Ground-truth VAD
            fprintf('Extract ground-truth VAD ... ')
            dir_vad = [ vad_dir, filesep, 'task', num2str(this_task), filesep 'recording', num2str(this_recording), filesep, this_array ];
            
            for src_idx = 1 : NoS
                
                fprintf('\n Source %s', src_in_rec{src_idx} );
                vad_fname = [ dir_vad, filesep, 'VAD_', this_array, '_', src_in_rec{src_idx}, '.txt'];
                
                % Load VAD vector
                vad_table = readtable( vad_fname );
                vad_vector = table2array( vad_table);
                L_audio = length( vad_vector );
                
                t_stamps_audio = ( 0: L_audio - 1 )/fs; % 48kHz
                t_stamps_opti = truth{arr_idx, rec_idx, task_idx}.timestamps; % 120Hz
                
                VAD = zeros( length(  t_stamps_opti ), 1 );
                
                % VAD values @48kHz matched with timestamps @120Hz
                cnt = 1;
                for i = 2 : L_audio
                    
                    if t_stamps_audio(i)>= t_stamps_opti( cnt )
                        VAD(cnt) = vad_vector(i-1);
                        cnt = cnt+1;
                    end
                    
                    if cnt > length(VAD)
                        break
                    end
                end
                
                if cnt <= length(VAD)
                    
                    VAD(cnt:end) = vad_vector(end);
                    % Mismatch of up to 2 values to be expected
                    if cnt < length(VAD)-1
                        warning('VAD values do not match')
                        delta_theo = round(L_audio/400) - length(VAD)  % theoretical difference
                        delta_act  = length(VAD) - cnt  % actual difference
                    end
                end
                
                % Calculate number of active periods
                vad_cnt = 0;
                for i=2:length(VAD)
                    
                    if VAD(i-1)==0 && VAD(i)==1 % speech onset
                        vad_cnt = vad_cnt+1;
                    end
                end
                
                if sum(VAD)==0
                    error('No VAD values found!')
                end
                
                % Store VAD ground-truth in structure
                truth{arr_idx,rec_idx,task_idx}.source.(src_in_rec{src_idx}).VAD.num_periods =  max( [ 1 vad_cnt ] ); % equal to one if no VAD is performed
                truth{arr_idx,rec_idx,task_idx}.source.(src_in_rec{src_idx}).VAD.activity = VAD; % VAD decisions for each timestamp
                
            end % eof sources
            
            if isempty( truth{arr_idx,rec_idx,task_idx} )
                error('Empty truth variable')
            end
            
        end % eof array loop
        
        fprintf('\n .... finished.\n')
        
    end % eof recordings
    
end % eof tasks

