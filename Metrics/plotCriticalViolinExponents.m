function plotCriticalViolinExponents(stats_list)

%states_list = ["sws", "other", "rem", "all"];
states_list = ["ISR", "other", "rem", "sws", "sws_non_ISR"];
event = "slownr";
event2 = "exceptslownr";
event_list = [event, "all", "all", "all", event2];

stats_names = {'\alpha', '\beta', '\gamma_{exp}', '\gamma_{area}', '\gamma_{shape}', '\chi', '\gamma_{area} - \gamma_{exp}'};
n_stats = 7;
n_conditions = 5; % sws, other, rem, slownr

%[fig, axs] = makeFigure('Critical Exponents', 'Critical Exponents', [2, 3]);
[fig,axs] = makeFigure('Critical Exponent','',[2, 4],'size',[35,20],'format','poster');

for stat = 1:n_stats
    s = stats_list{stat};
    
    % construire Y et group comme vecteurs
    Y = [];
    group = [];
    for c = 1:n_conditions
        if isfield(s, states_list(c))
            vals = s.(states_list(c))(:);
            Y = [Y; vals];
            group = [group; repmat(states_list(c), length(vals), 1)];
        end
    end

    [p,h] = StatTests(s);
    
    axes(axs(stat));
    colors = [0.20, 0.40, 0.70;   
          0.85, 0.45, 0.10;   
          0.15, 0.65, 0.35;  
          0.55, 0.25, 0.70;
          0.30, 0.85, 0.15];  
    xg = categorical(group);
    x_pos = 1.225 + (-0.5:0.1848:1.2392).';
    %v = violinplot(Y, 'GroupByColor',group);
    v = violinplot(ones(size(Y)),Y, 'GroupByColor',group, 'DensityWidth',1);
    for i = 1 : numel(v)
        v(i).FaceColor = colors(i,:);
        v(i).HandleVisibility = 'off';
    end
    xticks(x_pos);
    %xlim([0.5 numel(x_pos)+0.5])
    xticklabels(["ISR", "other", "rem", "sws", "sws non ISR"]);
    
   
  
    %clear xlim;
    display(xlim)
    %display(diff(xlim))
    %cla
    %plot(x_pos, zeros(size(x_pos)), 'LineStyle','none')  % pas de ligne ni de point  % fake plot numérique
    %xlim([0.5 numel(x_pos)+0.5])
    pBar(p.p1, x_pos, 0.05, draw=[true, true, false, false]);
% [n.s., *, **, ***] -> n'affiche que ** et ***
    title(axs(stat), stats_names{stat});
    ylabel(axs(stat), stats_names{stat});
    % if stat==7
    %     for c = 1:numel(v)
    %         y_pos = median(v(c).YData, 'omitnan');
    %         text(axs(stat), c, y_pos, ...
    %             sprintf('p=%.3f\nh=%d', p.p0(c), h.h0(c)), ...
    %             'HorizontalAlignment', 'center', ...
    %             'VerticalAlignment', 'middle', ...
    %             'FontSize', 8, ...
    %             'BackgroundColor', 'white');
    %     end
    % else
    %     for c = 1:numel(v)
    %         pairs = h.h1(h.h1(:,1)==c | h.h1(:,2)==c, :);
    %         p_pairs = p.p1(p.p1(:,1)==c | p.p1(:,2)==c, :);
    %         %display(pairs);
    %         %display(p_pairs);
    % 
    %         str = '';
    %         for k = 1:size(pairs,1)
    %             other = pairs(k,1) + pairs(k,2) - c;
    %             if c==4
    %                 display(other);
    %             end
    %             stars = repmat('*', 1, pairs(k,3));
    %             if isempty(stars), stars = 'ns'; end
    %             str = sprintf('%svs %s: p=%.3f %s\n', str, states_list(other), p_pairs(k,3), stars);
    %         end
    % 
    %         y_pos = median(v(c).YData, 'omitnan');
    %         text(axs(stat), c, y_pos, str, ...
    %             'HorizontalAlignment', 'center', ...
    %             'VerticalAlignment', 'middle', ...
    %             'FontSize', 7, ...
    %             'BackgroundColor', 'white');
    %     end
    % end
end

saveFig(gcf,"/mnt/hubel-data-103/Guillaume/NRAvalanche/Figures/svg_Files/ViolinPlots_pBars",'svg') % file_path: 'aaa/vvv/file_name'