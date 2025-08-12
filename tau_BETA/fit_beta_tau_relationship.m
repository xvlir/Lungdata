function fit_beta_tau_relationship(tau,beta)
    % 已知数据
    
    % 创建数据表格
    data_table = table(tau', beta', 'VariableNames', {'Tau', 'Beta'});
    disp('原始数据:');
    disp(data_table);
    
    % 保存数据
    writetable(data_table, 'beta_tau_data.csv');
    
    % 尝试多种非线性拟合
    fit_models = {
        % 模型名称, 模型函数, 初始参数
        '高斯模型', @(b,x) b(1)*exp(-((x-b(2))/b(3)).^2), [0.12, 25, 25];
        '洛伦兹模型', @(b,x) b(1)./((x-b(2)).^2 + b(3)), [0.12, 25, 100];
        '二次多项式', @(b,x) b(1)*x.^2 + b(2)*x + b(3), [0.0001, 0.01, 0.05];
        '指数衰减', @(b,x) b(1)*exp(-b(2)*x) + b(3), [0.06, 0.01, 0.05];
        '对数正态分布', @(b,x) b(1)*exp(-(log(x)-b(2)).^2/(2*b(3)^2)), [0.12, 3, 1];
        'S型曲线', @(b,x) b(1)./(1 + exp(-b(2)*(x-b(3)))), [0.12, 0.1, 25];
        '双指数', @(b,x) b(1)*exp(-b(2)*x) + b(3)*exp(-b(4)*x), [0.06, 0.01, 0.05, 0.001];
    };
    
    % 准备拟合
    best_fit = [];
    best_rsquared = -Inf;
    fit_results = cell(size(fit_models, 1), 1);
    
    figure('Position', [100, 100, 1000, 600]);
    hold on; grid on; box on;
    scatter(tau, beta, 100, 'filled', 'k', 'DisplayName', '原始数据');
    
    % 对每个模型进行拟合
    for i = 1:size(fit_models, 1)
        model_name = fit_models{i, 1};
        model_func = fit_models{i, 2};
        init_params = fit_models{i, 3};
        
        % 执行非线性拟合
        try
            mdl = fitnlm(tau, beta, model_func, init_params);
            
            % 计算R²
            y_pred = predict(mdl, tau');
            ss_res = sum((beta - y_pred').^2);
            ss_tot = sum((beta - mean(beta)).^2);
            rsquared = 1 - (ss_res / ss_tot);
            
            % 存储结果
            fit_results{i} = struct(...
                'Model', model_name, ...
                'Parameters', mdl.Coefficients.Estimate, ...
                'R_squared', rsquared, ...
                'ModelObj', mdl);
            
            % 绘制拟合曲线
            tau_fit = logspace(log10(0.05), log10(150), 200);
            beta_fit = predict(mdl, tau_fit');
            plot(tau_fit, beta_fit, 'LineWidth', 2, 'DisplayName', ...
                sprintf('%s (R²=%.4f)', model_name, rsquared));
            
            % 检查是否是最佳拟合
            if rsquared > best_rsquared
                best_rsquared = rsquared;
                best_fit = fit_results{i};
            end
        catch
            fprintf('模型 %s 拟合失败\n', model_name);
        end
    end
    
    % 设置图形属性
    set(gca, 'XScale', 'log');
    xlabel('Tau (τ)', 'FontSize', 14);
    ylabel('Beta (β)', 'FontSize', 14);
    title('Beta与Tau关系拟合', 'FontSize', 16);
    legend('Location', 'best', 'FontSize', 10);
    set(gca, 'FontSize', 12);
    
    % 保存图形
    saveas(gcf, 'beta_tau_fits.png');
    
    % 显示最佳拟合结果
    fprintf('\n最佳拟合模型: %s (R² = %.4f)\n', best_fit.Model, best_fit.R_squared);
    disp('拟合参数:');
    disp(best_fit.Parameters);
    
    % 创建拟合结果表格
    fit_table = table();
    for i = 1:size(fit_models, 1)
        if ~isempty(fit_results{i})
            fit_table.Model{i} = fit_results{i}.Model;
            fit_table.R_squared(i) = fit_results{i}.R_squared;
            param_str = sprintf('%.4f, ', fit_results{i}.Parameters);
            fit_table.Parameters{i} = param_str(1:end-2);
        end
    end
    
    % 按R²排序
    fit_table = sortrows(fit_table, 'R_squared', 'descend');
    
    disp('所有模型拟合结果:');
    disp(fit_table);
    
    % 保存拟合结果
    writetable(fit_table, 'beta_tau_fit_results.csv');
    
    % 使用最佳模型预测更多点
    tau_pred = logspace(log10(0.05), log10(200), 50);
    beta_pred = predict(best_fit.ModelObj, tau_pred');
    
    % 创建预测结果表格
    pred_table = table(tau_pred', beta_pred, 'VariableNames', {'Tau', 'Beta_Predicted'});
    writetable(pred_table, 'beta_tau_predictions.csv');
    
    % 绘制最佳拟合曲线
    figure('Position', [100, 100, 800, 600]);
    hold on; grid on; box on;
    
    % 原始数据
    scatter(tau, beta, 100, 'filled', 'k', 'DisplayName', '原始数据');
    
    % 最佳拟合曲线
    plot(tau_pred, beta_pred, 'r-', 'LineWidth', 2, ...
         'DisplayName', sprintf('%s (R²=%.4f)', best_fit.Model, best_fit.R_squared));
    
    % 设置图形属性
    set(gca, 'XScale', 'log');
    xlabel('Tau (τ)', 'FontSize', 14);
    ylabel('Beta (β)', 'FontSize', 14);
    title(sprintf('最佳拟合: %s', best_fit.Model), 'FontSize', 16);
    legend('Location', 'best', 'FontSize', 12);
    set(gca, 'FontSize', 12);
    
    % 添加公式标注
    params = best_fit.Parameters;
    switch best_fit.Model
        case '高斯模型'
            eq_str = sprintf('β = %.4f * exp(-((τ - %.4f)/%.4f)^2)', params(1), params(2), params(3));
        case '洛伦兹模型'
            eq_str = sprintf('β = %.4f / ((τ - %.4f)^2 + %.4f)', params(1), params(2), params(3));
        case '二次多项式'
            eq_str = sprintf('β = %.6fτ² + %.4fτ + %.4f', params(1), params(2), params(3));
        case '指数衰减'
            eq_str = sprintf('β = %.4f*exp(-%.4fτ) + %.4f', params(1), params(2), params(3));
        case '对数正态分布'
            eq_str = sprintf('β = %.4f * exp(-(ln(τ)-%.4f)^2/(2*%.4f^2))', params(1), params(2), params(3));
    end
    text(0.5, 0.1, eq_str, 'FontSize', 12, 'Units', 'normalized', ...
         'BackgroundColor', 'white');
    
    % 保存最佳拟合图形
    saveas(gcf, 'best_beta_tau_fit.png');
end