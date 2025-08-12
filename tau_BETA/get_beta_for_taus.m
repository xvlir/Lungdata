function beta_values = get_beta_for_taus(tau_values, base_dir)
    % tau_values: array of tau values
    % base_dir: base directory where results are stored
    % Output: beta_values corresponding to each tau

    % 预分配
    beta_values = nan(size(tau_values));
    
    % 获取所有文件夹
    dir_list = dir(base_dir);
    dir_list = dir_list([dir_list.isdir]); % 只保留文件夹
    dir_names = {dir_list.name};
    dir_names = dir_names(~ismember(dir_names, {'.', '..'})); % 去掉 . 和 ..
    
    % 遍历 tau
    for i = 1:length(tau_values)
        tau_target = tau_values(i);
        matched_folder = '';
        
        % 在所有文件夹中搜索 tau
        for j = 1:length(dir_names)
            name = dir_names{j};
            % 用正则匹配 _tau 后面的数字
            tokens = regexp(name, '_tau([\d\.]+)$', 'tokens');
            if ~isempty(tokens)
                tau_in_folder = str2double(tokens{1}{1});
                % 允许浮点误差比较
                if abs(tau_in_folder - tau_target) < 1e-6
                    matched_folder = name;
                    break;
                end
            end
        end
        
        if isempty(matched_folder)
            warning('No folder found for tau = %g', tau_target);
            continue;
        end
        
        % 拼接 CSV 路径
        csv_path = fullfile(base_dir, matched_folder, 'Relax', 'Strain04.csv');
        
        if ~isfile(csv_path)
            warning('File not found: %s', csv_path);
            continue;
        end
        
        % 读取 CSV
        data = readmatrix(csv_path);
        time_data = data(:,1);
        stress_data = data(:,3);
        
        % 调用拟合函数
        [~, ~, ~, ~, ~, beta, ~, ~, ~, ~] = ...
            fit_stress_relaxation(time_data, stress_data, csv_path);
        
        beta_values(i) = beta;
    end
    
    % tau-beta 拟合
    fit_beta_tau_relationship(tau_values, beta_values);
end
