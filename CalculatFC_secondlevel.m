clear all; close all;

%%%%%%%%%%%% make the FC mask based on the second-level analysis
load('ROI_infant_age.mat');
network_key = readtable('Network_key_Conn.csv'); 
% to index element in this table, use network_key{1,1}{1,1},and network_key{1,3}.

numROI = size(ROI, 2);
corr_name = ROI(1).names(1:numROI);

corr_h = [];   % Beta value
corr_F = [];   % Statistic value
corr_p = [];   % One-tailed p value
corr_z = [];   % load the individual data - have checked, consistent to the first-level output

for i = 1:numROI
	corr_h = [ corr_h ; ROI(i).h(1:numROI) ];
	corr_F = [ corr_F ; ROI(i).F(1:numROI) ];
	corr_p = [ corr_p ; ROI(i).p(1:numROI) ];
    corr_z = cat(3, corr_z, ROI(i).y);
end
numSub = size(corr_z,1);

% fetal_z = corr_z(1:129,:,:);
% infant_z = corr_z(130:188,:,:);
% figure
% histogram(mean(fetal_z,1));
% hold on
% histogram(mean(infant_z,1));
% legend('fetal FC','infant FC')


anti_corr_net = 1;
% Case two or several anti-correlated networks (ex: Default and Dorsal Attention)
% In our case, we need to see negative effect of age
if anti_corr_net
    tail = 'two-sided';
    corr_p = 2*min(corr_p, 1-corr_p);

% Case single network
else
    tail  = 'one-sided';

end

% Compute Bonferroni and FDR correction
% conn_fdr function is in conn main folder
alpha_bonf = 0.05 / ((numROI)*(numROI-1)/2);
corr_p(corr_p == 0)= 0.001;
vector_fdr = nonzeros(triu(corr_p)');
vector_fdr(isnan(vector_fdr)) = [];
corr_p_fdr = conn_fdr(vector_fdr);

% WRITE OUTPUT
fprintf('\nANALYSIS INFO');
fprintf('\n--------------------------------------');
fprintf('\nSTATISTICS');
fprintf('\n--------------------------------------');
fprintf([ '\n' num2str(numROI) ' x ' num2str(numROI-1) ' ROIs matrix ; ' tail]);
fprintf([ '\np-uncorrected:\t\t\t\t\t ' num2str(numel(corr_p(corr_p <= 0.05))/2)  ]);
fprintf([ '\np-bonferroni (alpha = ' num2str(round(alpha_bonf, 5)) '):\t ' num2str(numel(corr_p(corr_p <= alpha_bonf))/2) ]);
fprintf([ '\np-FDR corrected:\t\t\t\t ' num2str(numel(corr_p_fdr(corr_p_fdr <= 0.05))) ]);
fprintf('\n--------------------------------------\n');
% Outputs: p-FDR corrected 998

corr_p_fdr_mask = corr_p_fdr <= 0.05; % Fdr mask

corr_z_vec_thr = [];

for i = 1:numSub 
    z_mat = squeeze(corr_z(i, :, :));
    z_vec = nonzeros(triu(z_mat, 1)');
    z_vec_thr = z_vec .* corr_p_fdr_mask;
    corr_z_vec_thr = [corr_z_vec_thr, z_vec_thr];
end

F_vec = nonzeros(triu(corr_F, 1)');
F_vec_thr = F_vec .* corr_p_fdr_mask;

pos_ind = F_vec_thr > 0;
pos_ind_mat = repmat(pos_ind, 1, numSub);
neg_ind = F_vec_thr < 0;
neg_ind_mat = repmat(neg_ind, 1, numSub);

mean_Pos_fc = sum(corr_z_vec_thr .* pos_ind_mat) ./ sum(pos_ind_mat);
mean_Pos_fc = mean_Pos_fc';

mean_Neg_fc = sum(corr_z_vec_thr .* neg_ind_mat) ./ sum(neg_ind_mat);
mean_Neg_fc = mean_Neg_fc';

% load('dist_mat.mat')
% dist_vec = nonzeros(triu(cord_mat, 1)');
% dist_pos = dist_vec.* pos_ind;
% dist_neg = dist_vec.* neg_ind;
% 
% figure
% histogram(nonzeros(dist_pos));
% hold on
% histogram(nonzeros(dist_neg));
% legend('pos','neg')

%%%%%%%%%%%%%%%%%% Make network (8x8) correlation matrix
network_sum = zeros(8,8,numSub);
network_count = zeros(8,8,numSub);

for i = 1:numSub 
    z_mat = squeeze(corr_z(i, :, :));
    z_vec = nonzeros(triu(z_mat, 1)');
    z_vec_thr = z_vec .* corr_p_fdr_mask;
    z_mat_thr = nan(numROI);
    z_mat_thr(tril(true(numROI), -1)) = z_vec_thr;
    
    for m = 1:numROI
        netID_x = network_key{m,2};
        for n = 1:m-1
            netID_y = network_key{n,2};
            network_sum(netID_x, netID_y, i)= network_sum(netID_x, netID_y, i) + z_mat_thr(m, n);
            if network_sum(netID_x, netID_y, i)~=0
                network_count(netID_x, netID_y, i)= network_count(netID_x, netID_y, i) + 1;
            end
        end
    end
end

network_avg = zeros(36,numSub);
network_avg_mat = zeros(8,8,numSub);
for i = 1:numSub
    sub_sum = squeeze(network_sum(:,:,i));
    sub_cnt = squeeze(network_count(:,:,i));
    sub_avg = zeros(8,8);
    for m = 1:8
        for n = 1:m
            sub_avg(m,n) = (sub_sum(m,n) + sub_sum(n,m))/(sub_cnt(m,n) + sub_cnt(n,m));
        end
    end
    
    lowerTriangleVector = sub_avg(tril(true(size(sub_avg))));
    network_avg(:,i) = lowerTriangleVector;
    network_avg_mat(:,:,i) = sub_avg;
end

%plot correlation with age
demo = xlsread('Data_for_R_n184_final.xls','Data_for_R_n184_final','D2:D185');
age = demo(2:185,4);
preterm = demo(:,7);

for row = 1:8
    for col = 1:row
        % Calculate the current plot index
        plotIndex = (row - 1) * 8 + col;
        
        % Create subplot
        subplot(8, 8, plotIndex);
  
        % Calculate correlation coefficient
        y  = squeeze(network_avg_mat(row, col, :));
        [r,p] = corr(y, age,'rows', 'complete');
        
        % Create correlation plot
        scatter(age, y, 'filled');
        if p < 0.001
            title(['r= ', num2str(r), '***']);
        elseif p > 0.001 && p < 0.01
            title(['r= ', num2str(r), '**']);
        elseif p > 0.01 && p <0.05
            title(['r= ', num2str(r), '*']);
        else
            title(['r= ', num2str(r)]);
        end
            
       
        % Compute the coefficients of the trendline
        validIndices = ~isnan(age) & ~isnan(y);
        coefficients = polyfit(age(validIndices), y(validIndices), 1);
        trendline = coefficients(1) * age + coefficients(2);
        
        hold on; % Hold the current plot
        plot(age, trendline, 'r', 'LineWidth', 2);
        hold off; % Release the hold on the plot
        
        % Remove x and y labels from all plots except the bottom and left ones
        if row ~= 8
            set(gca, 'XTickLabel', []);
        end
        if col ~= 1
            set(gca, 'YTickLabel', []);
        end
    end
end

% Adjust layout
sgtitle('Grid of Correlation Plots', 'FontSize', 12);

% Adjust figure size if needed
set(gcf, 'Position', [100, 100, 1200, 1000]);


%%%%%%%%%%%% plot reterm
for row = 1:8
    for col = 1:row
        % Calculate the current plot index
        plotIndex = (row - 1) * 8 + col;
        
        % Create subplot
        subplot(8, 8, plotIndex);
  
        % Calculate correlation coefficient
        y  = squeeze(network_avg_mat(row, col, :));
        
        y1 = y(preterm == 0);
        x1 = age(preterm == 0);
        [r1,p1] = corr(x1, y1,'rows', 'complete');
        
        y2 = y(preterm == 1);
        x2 = age(preterm == 1);
        [r2,p2] = corr(x2, y2,'rows', 'complete');
        
        % Create correlation plot
        scatter(x1, y1);
        hold on
        scatter(x2, y2);
        
        if p1 < 0.001
            title(['r= ', num2str(r1), '***']);
        elseif p1 > 0.001 && p < 0.01
            title(['r= ', num2str(r1), '**']);
        elseif p1 > 0.01 && p <0.05
            title(['r= ', num2str(r1), '*']);
        else
            title(['r= ', num2str(r1)]);
        end
            
       
        % Compute the coefficients of the trendline
        coefficients1 = polyfit(x1, y1, 1);
        trendline1 = coefficients1(1) * x1 + coefficients1(2);
        
        coefficients2 = polyfit(x2, y2, 1);
        trendline2 = coefficients2(1) * x2 + coefficients2(2);
        
        hold on; % Hold the current plot
        plot(x1, trendline1, 'Color',[0 1 0], 'LineWidth', 2);
        plot(x2, trendline2, 'r', 'LineWidth', 2);
        
        hold off; % Release the hold on the plot
        
        % Remove x and y labels from all plots except the bottom and left ones
        if row ~= 8
            set(gca, 'XTickLabel', []);
        end
        if col ~= 1
            set(gca, 'YTickLabel', []);
        end
    end
end

% Adjust layout
sgtitle('Grid of Correlation Plots', 'FontSize', 12);

% Adjust figure size if needed
set(gcf, 'Position', [100, 100, 1200, 1000]);

