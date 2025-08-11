function [t_data, stress_data, t_fit, stress_fit, A, beta, A_err, beta_err, max_strain, max_stress] = ...
    fit_stress_relaxation(time, stress, filename)
    
    % 确定最大应变值（从文件名解析）
    strain_str = regexp(filename, 'Strain(\d+)', 'tokens');
    if ~isempty(strain_str)
        strain_num = str2double(strain_str{1}{1});
        max_strain = strain_num / 10;  % 01表示10%，02表示20%，依此类推
    else
        max_strain = NaN;
        warning('无法从文件名解析应变值');
    end
    
    % 找到应力最大值位置和值
    [max_stress, max_idx] = max(stress);
    
    % 提取松弛阶段数据
    t_data = time(max_idx:end);
    stress_data = stress(max_idx:end);
    
    % 准备拟合数据（对数变换）
    valid_idx = stress_data > 0;  % 确保正值
    t_data = t_data(valid_idx);
    stress_data = stress_data(valid_idx);
    
    % 检查是否有足够数据点
    if numel(t_data) < 2
        error('有效数据点不足（小于2个），无法进行拟合');
    end
    
    log_t = log(t_data);
    log_stress = log(stress_data);
    
    % 执行线性拟合
    mdl = fitlm(log_t, log_stress);
    coeff = mdl.Coefficients.Estimate;
    se = mdl.Coefficients.SE;
    
    % 提取拟合参数
    A = exp(coeff(1));
    beta = -coeff(2);
    A_err = A * se(1);  % 误差传递
    beta_err = se(2);
    
    % 生成拟合曲线
    t_fit = logspace(log10(min(t_data)), log10(max(t_data)), 100);
    stress_fit = A * t_fit.^(-beta);
end