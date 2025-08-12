function [phi_values, beta_values] = get_beta_for_phis(base_dir, alpha_folders)
    % 预分配
    n = length(alpha_folders);
    phi_values = nan(1, n);
    beta_values = nan(1, n);
    
    for i = 1:n
        folder_rel = alpha_folders{i};
        csv_path = fullfile(base_dir, folder_rel, 'Strain04.csv');
        
        if ~isfile(csv_path)
            warning('File not found: %s', csv_path);
            continue;
        end
        
        % 从文件夹名中提取 phi 值
        tokens = regexp(folder_rel, 'phi(\d+)', 'tokens');
        if isempty(tokens)
            warning('No phi found in folder: %s', folder_rel);
            continue;
        end
        phi_values(i) = str2double(tokens{1}{1});
        
        % 读取数据
        data = readmatrix(csv_path);
        time_data = data(:,1);
        stress_data = data(:,3);
        
        % 调用拟合函数
        [~, ~, ~, ~, ~, beta, ~, ~, ~, ~] = ...
            fit_stress_relaxation(time_data, stress_data, csv_path);
        
        beta_values(i) = beta;
    end
end
