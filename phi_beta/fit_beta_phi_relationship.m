function fit_beta_phi_relationship(phi,beta)
    % 如果未提供数据，使用默认值
    if nargin < 1
        phi = [0.1, 1, 25, 50, 100];
        beta = [0.054, 0.0996, 0.1138, 0.1032, 0.0897];
    end
    
    % 确保数据为行向量
    phi = phi(:)';
    beta = beta(:)';
    
    % 创建数据表格
    data_table = table(phi', beta', 'VariableNames', {'Tau', 'Beta'});
    disp('原始数据:');
    disp(data_table);
    
    % 保存数据
    writetable(data_table, 'beta_tau_data.csv');
    
    % 尝试多种非线性拟合
    fit_models = {
        % 模型名称, 模型函数, 初始参数
        '高斯模型', @(b,x) b(1)*exp(-((x-b(2))/b(3)).^2), [max(beta), median(phi), range(phi)/4];
        '洛伦兹模型', @(b,x) b(1)./((x-b(2)).^2 + b(3)), [max(beta), median(phi), range(phi)/2];
        '二次多项式', @(b,x) b(1)*x.^2 + b(2)*x + b(3), [0.0001, 0.01, min(beta)];
        '指数衰减', @(b,x) b(1)*exp(-b(2)*x) + b(3), [max(beta)-min(beta), 0.01, min(beta)];
        '对数正态分布', @(b,x) b(1)*exp(-(log(x)-b(2)).^2/(2*b(3)^2)), [max(beta), log(median(phi)), 0.7];
        'S型曲线', @(b,x) b(1)./(1 + exp(-b(2)*(x-b(3)))), [max(beta), 0.1, median(phi)];
        '双指数', @(b,x) b(1)*exp(-b(2)*x) + b(3)*exp(-b(4)*x), [max(beta)/2, 0.01, max(beta)/2, 0.001];
        '线性模型', @(b,x) b(1)*x + b(2), [(max(beta)-min(beta))/range(phi), min(beta)];
    };
    
    % 准备拟合
    best_fit = [];
    best_rsquared = -Inf;
    fit_results = cell(size(fit_models, 1), 1);
    
    figure('Position', [100, 100, 1200, 800], 'Name', 'Beta与Tau关系拟合比较');
    hold on; grid on; box on;
    scatter(phi, beta, 100, 'filled', 'k', 'DisplayName', '原始数据');
    
    % 创建颜色循环
    colors = lines(size(fit_models, 1));
    
    % 对每个模型进行拟合
    for i = 1:size(fit_models, 1)
        model_name = fit_models{i, 1};
        model_func = fit_models{i, 2};
        init_params = fit_models{i, 3};
        
        % 执行非线性拟合
        try
            % 设置优化选项
            opts = statset('nlinfit');
            opts.MaxIter = 1000;
            opts.Robust = 'on';
            
            % 拟合模型
            mdl = fitnlm(phi, beta, model_func, init_params, 'Options', opts);
            
            % 计算R²
            y_pred = predict(mdl, phi');
            ss_res = sum((beta - y_pred').^2);
            ss_tot = sum((beta - mean(beta)).^2);
            rsquared = 1 - (ss_res / ss_tot);
            
            % 存储结果
            fit_results{i} = struct(...
                'Model', model_name, ...
                'Parameters', mdl.Coefficients.Estimate, ...
                'SE', mdl.Coefficients.SE, ...
                'R_squared', rsquared, ...
                'ModelObj', mdl);
            
            % 绘制拟合曲线
            tau_min = min(phi(phi>0));
            tau_max = max(phi);
            tau_fit = logspace(log10(tau_min*0.5), log10(tau_max*2), 300);
            beta_fit = predict(mdl, tau_fit');
            
            plot(tau_fit, beta_fit, 'LineWidth', 2, 'Color', colors(i,:), ...
                 'DisplayName', sprintf('%s (R²=%.4f)', model_name, rsquared));
            
            % 检查是否是最佳拟合
            if rsquared > best_rsquared
                best_rsquared = rsquared;
                best_fit = fit_results{i};
            end
        catch ME
            fprintf('模型 %s 拟合失败: %s\n', model_name, ME.message);
            fit_results{i} = [];
        end
    end
    
    % 设置图形属性
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'linear');
    xlabel('Tau (τ)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Beta (β)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Beta-phi', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 优化图例位置
    lgd = legend('Location', 'bestoutside');
    lgd.FontSize = 10;
    set(gca, 'FontSize', 12);
    
    % 添加网格
    grid on;
    grid minor;
    
    % 保存图形
    saveas(gcf, 'beta_tau_fits_comparison.png');
    
    % 显示最佳拟合结果
    if ~isempty(best_fit)
        fprintf('\n最佳拟合模型: %s (R² = %.4f)\n', best_fit.Model, best_fit.R_squared);
        disp('拟合参数:');
        
        % 创建参数表格
        param_table = table();
        param_table.Parameter = cell(numel(best_fit.Parameters), 1);
        param_table.Value = best_fit.Parameters(:);
        param_table.SE = best_fit.SE(:);
        disp(param_table);
        
        % 创建所有模型拟合结果表格
        fit_table = table();
        for i = 1:size(fit_models, 1)
            if ~isempty(fit_results{i})
                fit_table.Model{i} = fit_results{i}.Model;
                fit_table.R_squared(i) = fit_results{i}.R_squared;
                
                % 格式化参数字符串
                param_str = '';
                for p = 1:numel(fit_results{i}.Parameters)
                    param_str = [param_str, sprintf('%.4f ± %.4f, ', ...
                        fit_results{i}.Parameters(p), fit_results{i}.SE(p))];
                end
                fit_table.Parameters{i} = param_str(1:end-2);
            else
                fit_table.Model{i} = fit_models{i,1};
                fit_table.R_squared(i) = NaN;
                fit_table.Parameters{i} = '拟合失败';
            end
        end
        
        % 按R²排序
        fit_table = sortrows(fit_table, 'R_squared', 'descend');
        
        disp('所有模型拟合结果:');
        disp(fit_table);
        
        % 保存拟合结果
        writetable(fit_table, 'beta_tau_fit_results.csv');
        
        % 使用最佳模型预测更多点
        tau_pred = logspace(log10(tau_min*0.3), log10(tau_max*3), 100);
        beta_pred = predict(best_fit.ModelObj, tau_pred');
        
        % 创建预测结果表格
        pred_table = table(tau_pred', beta_pred, 'VariableNames', {'Tau', 'Beta_Predicted'});
        writetable(pred_table, 'beta_tau_predictions.csv');
        
        % 绘制最佳拟合曲线
        figure('Position', [100, 100, 900, 700], 'Name', '最佳拟合模型');
        hold on; grid on; box on;
        
        % 原始数据
        scatter(phi, beta, 120, 'filled', 'k', 'DisplayName', '原始数据');
        
        % 最佳拟合曲线
        plot(tau_pred, beta_pred, 'r-', 'LineWidth', 3, ...
             'DisplayName', sprintf('%s (R²=%.4f)', best_fit.Model, best_fit.R_squared));
        
        % 设置图形属性
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'linear');
        xlabel('Tau (τ)', 'FontSize', 16, 'FontWeight', 'bold');
        ylabel('Beta (β)', 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('最佳拟合模型: %s', best_fit.Model), 'FontSize', 18, 'FontWeight', 'bold');
        
        % 添加公式标注
        params = best_fit.Parameters;
        se = best_fit.SE;
        switch best_fit.Model
            case '高斯模型'
                eq_str = sprintf('β = (%.4f±%.4f) * exp(-((τ - (%.4f±%.4f))/(%.4f±%.4f))²)', ...
                    params(1), se(1), params(2), se(2), params(3), se(3));
            case '洛伦兹模型'
                eq_str = sprintf('β = (%.4f±%.4f) / ((τ - (%.4f±%.4f))² + (%.4f±%.4f))', ...
                    params(1), se(1), params(2), se(2), params(3), se(3));
            case '二次多项式'
                eq_str = sprintf('β = (%.6f±%.6f)τ² + (%.4f±%.4f)τ + (%.4f±%.4f)', ...
                    params(1), se(1), params(2), se(2), params(3), se(3));
            case '指数衰减'
                eq_str = sprintf('β = (%.4f±%.4f)*exp(-(%.4f±%.4f)τ) + (%.4f±%.4f)', ...
                    params(1), se(1), params(2), se(2), params(3), se(3));
            case '对数正态分布'
                eq_str = sprintf('β = (%.4f±%.4f) * exp(-(ln(τ)-(%.4f±%.4f))²/(2*(%.4f±%.4f)²))', ...
                    params(1), se(1), params(2), se(2), params(3), se(3));
            case 'S型曲线'
                eq_str = sprintf('β = (%.4f±%.4f) / (1 + exp(-(%.4f±%.4f)(τ-(%.4f±%.4f)))', ...
                    params(1), se(1), params(2), se(2), params(3), se(3));
            case '双指数'
                eq_str = sprintf('β = (%.4f±%.4f)*exp(-(%.4f±%.4f)τ) + (%.4f±%.4f)*exp(-(%.4f±%.4f)τ)', ...
                    params(1), se(1), params(2), se(2), params(3), se(3), params(4), se(4));
            case '线性模型'
                eq_str = sprintf('β = (%.4f±%.4f)τ + (%.4f±%.4f)', ...
                    params(1), se(1), params(2), se(2));
        end
        
        % 添加公式框
        annotation('textbox', [0.15, 0.15, 0.7, 0.1], 'String', eq_str, ...
                   'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
                   'FontSize', 12, 'EdgeColor', 'black', 'LineWidth', 1);
        
        % 添加图例
        legend('Location', 'best', 'FontSize', 12);
        set(gca, 'FontSize', 14);
        
        % 添加网格
        grid on;
        grid minor;
        
        % 保存最佳拟合图形
        saveas(gcf, 'best_beta_tau_fit.png');
        
        % 保存高分辨率图片
        print('best_beta_tau_fit_highres', '-dpng', '-r300');
    else
        error('没有模型成功拟合数据。');
    end
    
    % 创建模型性能比较图
    figure('Position', [100, 100, 800, 600], 'Name', '模型性能比较');
    hold on; grid on; box on;
    
    % 收集所有R²值
    r2_values = [];
    model_names = {};
    for i = 1:size(fit_models, 1)
        if ~isempty(fit_results{i})
            r2_values(end+1) = fit_results{i}.R_squared;
            model_names{end+1} = fit_models{i,1};
        end
    end
    
    % 绘制条形图
    [sorted_r2, sort_idx] = sort(r2_values, 'descend');
    sorted_names = model_names(sort_idx);
    
    barh(sorted_r2, 'FaceColor', [0.2, 0.6, 0.8]);
    set(gca, 'YTick', 1:length(sorted_names));
    set(gca, 'YTickLabel', sorted_names);
    xlabel('R²值', 'FontSize', 14, 'FontWeight', 'bold');
    title('模型拟合性能比较', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 添加数值标签
    for i = 1:length(sorted_r2)
        text(sorted_r2(i) + 0.01, i, sprintf('%.4f', sorted_r2(i)), ...
             'VerticalAlignment', 'middle', 'FontSize', 10);
    end
    
    % 设置图形属性
    set(gca, 'FontSize', 12);
    xlim([0, min(1.1, max(sorted_r2)*1.1)]);
    
    % 保存图形
    saveas(gcf, 'model_performance_comparison.png');
end