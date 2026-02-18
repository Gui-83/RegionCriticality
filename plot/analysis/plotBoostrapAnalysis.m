function plotBoostrapAnalysis(Ealpha, Ebeta, Egamma1, Egamma2, Egamma3)
data = {Ealpha, Ebeta, Egamma1, Egamma2, Egamma3};
    labels = {'alpha', 'beta', 'gamma Area fit', 'Egamma scaling law', 'Egamma shape collapse'};
    
    % Create a figure
    figure;
    for i = 1:length(data)
        subplot(2, 3, i); % Create a 2x3 grid of subplots
        y = zscore(data{i});
        [~, p] = kstest(y);
        histogram(data{i}, 'Normalization', 'pdf', 'FaceAlpha', 0.7);
        legend(sprintf('%s, p=%.2f', labels{i}, p));
        xlabel('Value');
        ylabel('Probability Density');
        
        % Add legend for each subplot
        legend;
    end

    % Concatenate the data into a single array
    allData = cell2mat(data');
    
    subplot(2, 3, 6);
    boxplot(allData');
    title('Box Plot of Ealpha, Ebeta, Egamma1, Egamma2, and Egamma3');
    xlabel('Groups');
    ylabel('Values');
end

