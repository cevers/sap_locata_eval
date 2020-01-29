function rm_VAD_legend(results_dir, fname, VAD_flag)

if nargin < 3
    VAD_flag = 0;
end

fig = openfig([results_dir filesep, fname, '.fig']);

nrs = 0;
for obj_idx = 1 : length(fig.Children(2).Children)
    if ~isempty(regexp(fig.Children(2).Children(obj_idx).DisplayName, 'VAD'))
        fig.Children(2).Children(obj_idx).Annotation.LegendInformation.IconDisplayStyle = 'off';
        fig.Children(2).Children(obj_idx).Visible = 'off';
    elseif ~isempty(regexp(fig.Children(2).Children(obj_idx).DisplayName, 'data'))
        fig.Children(2).Children(obj_idx).Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    legend off
    legend
end

if VAD_flag
    nsrcs = sum(~cellfun(@isempty, (regexp({fig.Children(2).Children.DisplayName}, 'Source'))));
    cols = hsv(nsrcs+1);
    src_idx = 1;
    for obj_idx = 1 : length(fig.Children(2).Children)
        if ~isempty(regexp(fig.Children(2).Children(obj_idx).DisplayName, 'Source'))
            fig.Children(2).Children(obj_idx).Color = cols(src_idx,:);
            src_idx = src_idx + 1;
        end
    end
end

saveas(fig, [results_dir filesep, fname, '.fig']);
saveas(fig, [results_dir filesep, fname, '.png']);

return
