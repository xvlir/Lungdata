function fit_beta_phi_linear(phi, beta, save_csv_path)
    % 线性拟合 beta = m * phi + c
    % phi: phi 数组
    % beta: beta 数组
    % save_csv_path: CSV 文件保存路径

    % 确保列向量
    phi = phi(:);
    beta = beta(:);

    % 线性拟合
    p = polyfit(phi, beta, 1); % p(1)=m, p(2)=c
    beta_fit = polyval(p, phi);

    % 计算 R²
    SS_res = sum((beta - beta_fit).^2);
    SS_tot = sum((beta - mean(beta)).^2);
    R2 = 1 - SS_res/SS_tot;

    % 保存数据到 CSV
    data_out = [phi, beta];
    writematrix(data_out, save_csv_path);

    % 绘图
    figure;
    scatter(phi, beta, 60, 'filled', 'DisplayName', 'Data Points');
    hold on;
    phi_range = linspace(min(phi), max(phi), 200);
    beta_range = polyval(p, phi_range);
    plot(phi_range, beta_range, 'r-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

    % 拟合公式字符串
    fit_eq = sprintf('\\beta = %.4f \\times \\phi + %.4f,  R^2 = %.4f', p(1), p(2), R2);
    text(mean(phi), max(beta), fit_eq, 'FontSize', 10, 'HorizontalAlignment', 'center');

    % 图形细节
    xlabel('\phi', 'FontSize', 12);
    ylabel('\beta', 'FontSize', 12);
    title('\beta-\phi Linear Fit', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on;

    % 输出拟合参数
    fprintf('Linear Fit: beta = %.6f * phi + %.6f\n', p(1), p(2));
    fprintf('R-squared = %.6f\n', R2);
end
