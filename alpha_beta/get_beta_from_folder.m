function [beta_vec, strains] = get_beta_from_folder(folder_path, num_files)
    % 预存储结果
    beta_vec = zeros(1, num_files);
    strains = zeros(1, num_files);
    
    % 循环处理每个文件
    for i = 1:num_files
        % 构建文件名
        filename = fullfile(folder_path, sprintf('Strain0%d.csv', i));
        
        % 读取数据
        try
            data = readmatrix(filename);
            time = data(:,1);
            stress = data(:,3);
        catch
            error('无法读取文件: %s', filename);
        end
        
        % 拟合数据并获取beta值
        [~, ~, ~, ~, ~, beta, ~, ~, max_strain] = fit_stress_relaxation(time, stress, filename);
        
        % 存储结果
        beta_vec(i) = beta;
        strains(i) = max_strain;
    end
end