% 主函数：分析不同alpha下的beta值并拟合关系，输出结果表格
function [result_table, fit_result] = beta_alpha_analysis(alpha_folders,num1,num2)
    % 基础路径
    base_dir = 'E:\LungFEM_Study\Results\RawData';
    
    % 定义不同alpha对应的文件夹
    % alpha_folders = {
    %     'B01_FZRH_E03E13_phi59_alpha1_tau25\Relax', ...  % alpha = 0.1
    %     'B02_FZRH_E03E13_phi59_alpha3_tau25\Relax', ...  % alpha = 0.3
    %     'B03_FZRH_E03E13_phi59_alpha5_tau25\Relax', ...  % alpha = 0.5
    %     'B04_FZRH_E03E13_phi59_alpha7_tau25\Relax', ...  % alpha = 0.7
    %     'B05_FZRH_E03E13_phi59_alpha9_tau25\Relax'       % alpha = 0.9
    % };
    
    % 对应的alpha值
    alpha_values = 0.1:0.2:0.9;
    % num_alpha = length(alpha_values);
    
    % 预存储所有结果
    all_data = []; % 存储所有数据点
    avg_beta = zeros(1, num2-num1); % 存储平均beta值
    
    % 循环处理每个alpha文件夹
    for i = 1:length(alpha_values)
        folder_path = fullfile(base_dir, alpha_folders{i+num1-1});
        num_files = 5;  % 每个文件夹有5个应变文件
        
        % 获取该alpha下的所有beta值
        [beta_vec, strains] = get_beta_from_folder(folder_path, num_files);
        
        % 存储结果
        avg_beta(i) = mean(beta_vec);
        
        % 为当前alpha创建数据块
        alpha_col = repmat(alpha_values(i), num_files, 1);
        strain_col = strains(:) * 100; % 转换为百分比
        beta_col = beta_vec(:);
        
        % 添加到总数据集
        all_data = [all_data; alpha_col, strain_col, beta_col];
    end
    
    % 创建结果表格
    result_table = array2table(all_data, ...
        'VariableNames', {'Alpha', 'Strain_%', 'Beta'});
    
    % 创建平均beta值表格
    avg_table = table(alpha_values', avg_beta', ...
        'VariableNames', {'Alpha', 'Average_Beta'});
    
    % 显示结果表格
    disp('所有数据点结果:');
    disp(result_table);
    disp(' ');
    disp('平均Beta值:');
    disp(avg_table);
    
    % 保存结果为CSV文件
    writetable(result_table, 'beta_results_all.csv');
    writetable(avg_table, 'beta_results_avg.csv');
    disp('结果已保存为CSV文件: beta_results_all.csv 和 beta_results_avg.csv');
    
    % 线性拟合 (使用平均beta值)
    p = polyfit(alpha_values, avg_beta, 1);
    fit_line = polyval(p, alpha_values);
    
    % 创建拟合结果表格
    fit_result = table(p(1), p(2), ...
        'VariableNames', {'Slope', 'Intercept'});
    
    % 显示拟合结果
    disp(' ');
    disp('Beta与Alpha关系拟合结果:');
    disp(fit_result);
    fprintf('拟合公式: β = %.4fα + %.4f\n', p(1), p(2));
    
    % 绘制beta-alpha关系图
    figure('Position', [100, 100, 800, 600], 'Name', 'Beta与Alpha关系分析');
    hold on; grid on; box on;
    
    % 绘制所有数据点
    scatter(result_table.Alpha, result_table.Beta, 80, 'filled', ...
            'DisplayName', '所有数据点');
    
    % 绘制平均beta值
    plot(alpha_values, avg_beta, 'ko-', 'LineWidth', 2, 'MarkerSize', 10, ...
         'DisplayName', '平均beta值');
    
    % 绘制拟合线
    plot(alpha_values, fit_line, 'r--', 'LineWidth', 2, ...
         'DisplayName', sprintf('拟合: β = %.4fα + %.4f', p(1), p(2)));
    
    % 设置图形属性
    xlabel('Alpha (α)', 'FontSize', 14);
    ylabel('指数 Beta (β)', 'FontSize', 14);
    title('Beta与Alpha关系分析', 'FontSize', 16);
    legend('Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    
    % 添加公式标注
    text(0.5, min(avg_beta)*0.95, sprintf('拟合公式: β = %.4fα + %.4f', p(1), p(2)), ...
         'FontSize', 12, 'BackgroundColor', 'white');
    
    hold off;
    
    % 保存图形
    saveas(gcf, 'beta_alpha_relationship.png');
    disp('图形已保存为: beta_alpha_relationship.png');
end
