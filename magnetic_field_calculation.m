clear; clc; close all;
%% ===================== 物理常数与几何参数 =====================
MU0 = 4*pi*1e-7;          % 真空磁导率 [H/m]
J0 = 2e7;                 % 体电流密度 [A/m²]
EPS = 1e-4;              % 数值稳定性常数
COEF_K = MU0*J0/(4*pi);   % 磁场计算核心系数 [T]

% 线圈几何参数
R_INNER = 2;              % 内侧直线段距中心距离 [m]
L_STRAIGHT = 3;           % 直线段长度 [m]
Z_PLANE = 0;              % z平面坐标 [m]

% 离散化参数
N_BASE_POINTS = 100;      % 线圈基础离散点数
N_ROT_ANGLES = 18;        % 绕Y轴旋转分段数（360°均分）
TOTAL_ROT_ANGLE = 360;    % 总旋转角度 [deg]
N_ARC_SEGMENTS = 4;       % 每个圆弧段的细分段数

% 线圈截面尺寸（矩形）
C_HALF = 0.25;            % 半宽度 [m]
D_HALF = 0.25;            % 半高度 [m]
%% ===================== 圆弧几何参数定义 =====================
arc_radius = zeros(1, 6);
arc_angle_deg = [70, 60, 50, 50, 60, 70];  % 6段圆弧角度
arc_angle_rad = deg2rad(arc_angle_deg);
arc_center_x = zeros(1, 6);
arc_center_y = zeros(1, 6);
arc_start_angle = zeros(1, 6);
arc_end_angle = zeros(1, 6);

% 第1段圆弧（右上）
arc_radius(1) = 0.6;
arc_center_x(1) = arc_radius(1) + R_INNER;
arc_center_y(1) = L_STRAIGHT / 2;
arc_start_angle(1) = 0;
arc_end_angle(1) = arc_angle_deg(1);

% 第2段圆弧（右中）
arc_radius(2) = 1.0;
arc_center_x(2) = arc_center_x(1) + arc_radius(1)*cos(arc_angle_rad(1))*(arc_radius(2)-arc_radius(1))/arc_radius(1);
arc_center_y(2) = arc_center_y(1) - arc_radius(1)*sin(arc_angle_rad(1))*(arc_radius(2)-arc_radius(1))/arc_radius(1);
arc_start_angle(2) = arc_end_angle(1);
arc_end_angle(2) = arc_start_angle(2) + arc_angle_deg(2);

% 第3段圆弧（右下→左下，上半）
arc_radius(3) = arc_center_y(2)/sin(arc_angle_rad(3)) + arc_radius(2);
arc_center_x(3) = arc_center_x(2) - cos(arc_angle_rad(3))*(arc_radius(3) - arc_radius(2));
arc_center_y(3) = 0;
arc_start_angle(3) = arc_end_angle(2);
arc_end_angle(3) = arc_start_angle(3) + arc_angle_deg(3);

% 第4段圆弧（左下→右下，下半）
arc_radius(4) = arc_radius(3);
arc_center_x(4) = arc_center_x(3);
arc_center_y(4) = 0;
arc_start_angle(4) = arc_end_angle(3);
arc_end_angle(4) = arc_start_angle(4) + arc_angle_deg(4);

% 第5段圆弧（左中）
arc_radius(5) = arc_radius(2);
arc_center_x(5) = arc_center_x(2);
arc_center_y(5) = -arc_center_y(2);
arc_start_angle(5) = arc_end_angle(4);
arc_end_angle(5) = arc_start_angle(5) + arc_angle_deg(5);

% 第6段圆弧（左上）
arc_radius(6) = arc_radius(1);
arc_center_x(6) = arc_center_x(1);
arc_center_y(6) = -arc_center_y(1);
arc_start_angle(6) = arc_end_angle(5);
arc_end_angle(6) = arc_start_angle(6) + arc_angle_deg(6);
%% ===================== 线圈长度计算与基础点生成 =====================
len_straight = L_STRAIGHT;
len_arcs = arc_radius .* arc_angle_rad;
total_coil_length_upper = (len_straight + sum(len_arcs))/2;
avg_step_length = total_coil_length_upper / (N_BASE_POINTS - 1);

% 生成线圈基础点
coil_base_points = generateCoilBasePoints(R_INNER, L_STRAIGHT, arc_radius, arc_center_x, arc_center_y, ...
                                          arc_angle_rad, arc_start_angle, arc_end_angle, ...
                                          N_BASE_POINTS, avg_step_length, Z_PLANE);
coil_base_points(:,1) = coil_base_points(:,1) - 0.08;  % 内偏8cm，根据你的R_INNER调整
%% ===================== 绕Y轴旋转角度定义 =====================
rot_angles_y = linspace(0, TOTAL_ROT_ANGLE - TOTAL_ROT_ANGLE/N_ROT_ANGLES, N_ROT_ANGLES);
rot_angles_y_rad = deg2rad(rot_angles_y);
%% ===================== 局部坐标系原点与角度 =====================
local_frame_origins = generateLocalFrameOrigins(R_INNER, arc_radius, arc_center_x, arc_center_y, ...
                                                arc_start_angle, arc_end_angle, N_ARC_SEGMENTS, Z_PLANE, N_ROT_ANGLES, rot_angles_y_rad, L_STRAIGHT);

local_frame_angles = generateLocalFrameAngles(arc_angle_deg, N_ARC_SEGMENTS);
local_frame_angles_rad = deg2rad(local_frame_angles);
%% ===================== 坐标变换到局部坐标系 =====================
points_local_frame = rotateCoil(coil_base_points, N_ROT_ANGLES, rot_angles_y_rad, local_frame_origins, local_frame_angles_rad);
%% ===================== 梯形棱柱几何参数 =====================
[trap_bottom_width, trap_alpha, trap_beta] = generateTrapezoidParameters(L_STRAIGHT, arc_radius, arc_angle_rad, ...
                                                                          N_ARC_SEGMENTS, D_HALF);
%% ===================== 磁场计算（并行） =====================
B_local_frame = calculateMagneticFieldParallel(points_local_frame, trap_bottom_width, C_HALF, D_HALF, ...
                                               trap_alpha, trap_beta, COEF_K, EPS);
%% ===================== 完全反变换到全局坐标系 =====================
% 核心修改：所有磁场矢量完全反变换到全局坐标系后再求和
B_global_complete = fullyInverseTransformMagneticField(B_local_frame, local_frame_angles_rad, rot_angles_y_rad);
%% ===================== 全局坐标系下求和 =====================
% 现在所有磁场矢量都在全局坐标系下，可以直接叠加
B_total_coil = squeeze(sum(sum(B_global_complete, 1), 2));
% 先对25个局部坐标系求和（维度1），再对18个旋转角度求和（维度2）
%% ===================== 结果可视化 =====================
plotMagneticFieldResults(B_total_coil, total_coil_length_upper, N_BASE_POINTS);
%% =========================================================================
%  子函数定义
%% =========================================================================
function base_points = generateCoilBasePoints(r_inner, L_straight, arc_r, arc_cx, arc_cy, ...
                                              arc_theta, arc_start, arc_end, n_points, step_len, z_plane)
    base_points = zeros(n_points, 3);
    idx = 1;
    % 直线段
    n_line = max(1, round(L_straight/2 / step_len));
    y_line = linspace(0, L_straight/2, n_line);
    for i = 1:n_line
        if idx > n_points, break; end
        base_points(idx, :) = [r_inner, y_line(i), z_plane];
        idx = idx + 1;
    end
    % 上半部分圆弧
    for arc_idx = 1:3
        if idx > n_points, break; end
        remaining = n_points - idx + 1;
        arc_len = arc_r(arc_idx) * arc_theta(arc_idx);
        n_arc = max(1, min(round(arc_len / step_len), remaining));
        theta_seq = linspace(deg2rad(arc_start(arc_idx)), deg2rad(arc_end(arc_idx)), n_arc);
        for i = 1:n_arc
            if idx > n_points, break; end
            x = arc_cx(arc_idx) + arc_r(arc_idx) * cos(pi - theta_seq(i));
            y = arc_cy(arc_idx) + arc_r(arc_idx) * sin(pi - theta_seq(i));
            base_points(idx, :) = [x, y, z_plane];
            idx = idx + 1;
        end
    end
    while idx <= n_points
        base_points(idx, :) = base_points(idx-1, :);
        idx = idx + 1;
    end
end

function origins_all = generateLocalFrameOrigins(r_inner, arc_r, arc_cx, arc_cy, ...
                                             arc_start, arc_end, n_segments, z_plane ,n_angles ,rot_angles_y, L_STRAIGHT)
    n_arcs = 6;
    total_frames = n_segments + n_arcs * n_segments;
    origins = zeros(total_frames, 3);
    % ─────────────── 直线段：均匀分成4份 ───────────────
    y_positions = linspace(-L_STRAIGHT/2, L_STRAIGHT/2, n_segments + 1);  % 5个点 → 4个中点
    for seg = 1:n_segments
        y_mid = (y_positions(seg) + y_positions(seg+1)) / 2;
        origins(seg, :) = [r_inner, y_mid, z_plane];
    end
    current_idx = 1 + n_segments;
    for arc_idx = 1:n_arcs
        theta_split = linspace(deg2rad(arc_start(arc_idx)), deg2rad(arc_end(arc_idx)), n_segments + 1);
        split_pts = zeros(n_segments + 1, 3);
        for i = 1:length(theta_split)
            theta = theta_split(i);
            split_pts(i, 1) = arc_cx(arc_idx) + arc_r(arc_idx) * cos(pi - theta);
            split_pts(i, 2) = arc_cy(arc_idx) + arc_r(arc_idx) * sin(pi - theta);
            split_pts(i, 3) = z_plane;
        end
        for seg_idx = 1:n_segments
            origins(current_idx, :) = (split_pts(seg_idx, :) + split_pts(seg_idx + 1, :)) / 2;
            current_idx = current_idx + 1;
        end
    end
    origins_all = zeros(total_frames, n_angles, 3);
    for  i= 1:n_angles
       R_y = [cos(rot_angles_y(i)), 0, -sin(rot_angles_y(i));
               0,          1, 0;
               sin(rot_angles_y(i)), 0, cos(rot_angles_y(i))];
       for m = 1:total_frames
           origins_all(m, i, :) = R_y * (origins(m, :))';
       end

    end
end

function angles = generateLocalFrameAngles(arc_angle_deg, n_segments)
    n_arcs = 6;
    total_frames = n_segments + n_arcs * n_segments;
    angles = zeros(1, total_frames);
    % 直线段全部设为 0
    angles(1:n_segments) = 0;
    current_idx = n_segments + 1;
    cumulative_angle = 0;
    for arc_idx = 1:n_arcs
        angle_deg = arc_angle_deg(arc_idx);
        base_angle = cumulative_angle + angle_deg/8;
        for seg_idx = 1:n_segments
            angles(current_idx) = base_angle + angle_deg/4 * (seg_idx - 1);
            current_idx = current_idx + 1;
        end
        cumulative_angle = cumulative_angle + angle_deg;
    end
end

function points_local_frame = rotateCoil(base_points, n_angles, rot_angles_y, origins, rot_angles_z)
    [n_points, ~] = size(base_points);
    n_frames = length(rot_angles_z);
    points_local_frame = zeros(n_frames, n_angles, n_points, 3);
    for i = 1:n_angles
        R_y = [cos(rot_angles_y(i)), 0, -sin(rot_angles_y(i));
               0,          1, 0;
               sin(rot_angles_y(i)), 0, cos(rot_angles_y(i))];
        for k = 1:n_frames
            cos_t = cos(rot_angles_z(k));
            sin_t = sin(rot_angles_z(k));
            R_z = [cos_t,  sin_t, 0;
                  -sin_t,  cos_t, 0;
                   0,      0,     1];
            for j = 1:n_points
                translated = base_points(j, :) - squeeze(origins(k, i, :))';
                % 修正：旋转顺序应为 R_z * R_y' 而非 R_y * R_z
                % 场点在基平面（未旋转），源段已绕Y轴旋转
                % 正确的局部坐标变换：先用 R_y' 将场点带入源段参考系，再用 R_z 旋转到局部坐标系
                points_local_frame(k, i, j, :) = (R_z * R_y' * translated')';
            end
        end
    end
end

function [b_vals, alpha_vals, beta_vals] = generateTrapezoidParameters(L_straight, arc_r, arc_angle_rad, n_segments, d_half)
    n_arcs = 6;
    total_frames = n_segments + n_arcs * n_segments;
    b_vals = zeros(total_frames, 1);
    alpha_vals = zeros(total_frames, 1);
    beta_vals = zeros(total_frames, 1);
    % 直线段：复制4份相同参数
    b_vals(1:n_segments)     = L_straight/4;
    alpha_vals(1:n_segments) = 0;
    beta_vals(1:n_segments)  = 0;
    % 计算圆弧段的b参数
    b_arcs = zeros(1, n_arcs);
    for n = 1:n_arcs
        angle_rad = arc_angle_rad(n);
        r = arc_r(n);
        % b = (r*cos(θ/8) - d) * tan(θ/8)
        b_arcs(n) = (r*cos(angle_rad/8) - d_half) * sin(angle_rad/8)/cos(angle_rad/8);
    end
    current_idx = n_segments + 1;
    for arc_idx = 1:n_arcs
        angle_rad = arc_angle_rad(arc_idx);
        for seg_idx = 1:n_segments
            b_vals(current_idx) = b_arcs(arc_idx);
            alpha_vals(current_idx) = angle_rad/8;
            beta_vals(current_idx) = angle_rad/8;
            current_idx = current_idx + 1;
        end
    end
end

function B_field = calculateMagneticFieldParallel(local_coords, b_vals, c, d, alpha_vals, beta_vals, coef_k, eps)
    [n_frames, n_angles, n_points, ~] = size(local_coords);
    B_field = zeros(n_frames, n_angles, n_points, 3);
    parfor k = 1:n_frames
        B_frame = zeros(n_angles, n_points, 3);
        for i = 1:n_angles
            for j = 1:n_points
                P = squeeze(local_coords(k, i, j, :));
                [Bx, By, Bz] = trapezoidalPrismMagneticField(P, b_vals(k), c, d, ...
                                                   alpha_vals(k), beta_vals(k), coef_k, eps);
                B_frame(i, j, 1) = Bx;
                B_frame(i, j, 2) = By;
                B_frame(i, j, 3) = Bz;
            end
        end
        B_field(k, :, :, :) = B_frame;
    end
end

function [Bx, By, Bz] = trapezoidalPrismMagneticField(P, b, c, d, alpha, beta, coef_k, eps)
    x = -P(1); y = P(2); z = P(3);
    L1 = -d - x;    L2 = d - x;
    Q1 = -c - z;    Q2 = c - z;
    R1 = (d + x)*tan(alpha) + b - y;
    R2 = (d + x)*tan(beta) + b + y;
    By = 0;
    % 计算Bx
    Sx = 0;
    for n = 1:4
        [theta, R, Q] = getSxParameters(n, alpha, beta, R1, R2, Q2, Q1);
        Sxn = calculateSxn(theta, R, Q, L1, L2, eps);
        Sx = Sx + (-1)^(n-1) * Sxn;
    end
    Bx = Sx * coef_k;
    % 计算Bz
    Sz = 0;
    for n = 1:4
        [theta, R, L] = getSzParameters(n, alpha, beta, R1, R2, L2, L1);
        Szn = calculateSzn(theta, R, L, Q1, Q2, eps);
        Sz = Sz + (-1)^(n-1) * Szn;
    end
    Bz = Sz * coef_k;
end

function [theta, R, Q] = getSxParameters(n, alpha, beta, R1, R2, Q2, Q1)
    switch n
        case 1; theta = alpha; R = R1; Q = Q2;
        case 2; theta = alpha; R = R1; Q = Q1;
        case 3; theta = beta;  R = R2; Q = Q2;
        case 4; theta = beta;  R = R2; Q = Q1;
    end
end

function [theta, R, L] = getSzParameters(n, alpha, beta, R1, R2, L2, L1)
    switch n
        case 1; theta = alpha; R = R1; L = L2;
        case 2; theta = alpha; R = R1; L = L1;
        case 3; theta = beta;  R = R2; L = L2;
        case 4; theta = beta;  R = R2; L = L1;
    end
end

function Sxn = calculateSxn(theta, R, Q, L1, L2, eps)
    st = sin(theta); ct = cos(theta);
    Sxn = sxnIntegrand(L2, st, ct, R, Q, eps) - sxnIntegrand(L1, st, ct, R, Q, eps);
end

function s = sxnIntegrand(t, st, ct, R, Q, eps)
    term1 = t * asinh(eps +(t*st + R*ct) / (ct * sqrt(t^2 + Q^2+ eps))+ eps);
    term2 = R*ct * asinh((t + R*st*ct) / (ct * sqrt(R^2*ct^2 + Q^2 + eps) + eps)+ eps);
    term3 = Q * atan((Q^2*st - t*R*ct) / (Q*sqrt(t^2 + 2*R*t*st*ct + (R^2+Q^2)*ct^2 + eps) + eps)+ eps);
    s = term1 + term2 + term3;
end

function Szn = calculateSzn(theta, R, L, Q1, Q2, eps)
    st = sin(theta); ct = cos(theta);
    Szn = sznIntegrand(Q2, st, ct, R, L, eps) - sznIntegrand(Q1, st, ct, R, L, eps);
end

function s = sznIntegrand(t, st, ct, R, L, eps)
    term1 = t*st * asinh((L + R*st*ct) / (ct * sqrt(t^2 + R^2*ct^2 + eps) + eps)+ eps);
    term2 = -t * asinh((L*st + R*ct) / (ct * sqrt(t^2 + L^2 + eps) + eps)+ eps);
    term3 = -R*ct^2 * asinh((t*ct) / sqrt(L^2 + 2*L*R*st*ct + R^2*ct^2 + eps) + eps);
    term4 = -R*st*ct * atan((t*(L + R*st*ct)) / (R*ct * sqrt(t^2*ct^2 + L^2 + 2*L*R*st*ct + R^2*ct^2 + eps) + eps)+ eps);
    term5 = L * atan((t*(L*st + R*ct)) / (L * sqrt(t^2*ct^2 + L^2 + 2*L*R*st*ct + R^2*ct^2 + eps) + eps)+ eps);
    s = term1 + term2 + term3 + term4 + term5;
end

function B_global = fullyInverseTransformMagneticField(B_local, local_angles_rad, rot_angles_y_rad)
    % 完全反变换：将局部坐标系中的磁场矢量完全反变换到全局坐标系
    % 正变换链：全局 → 中间（R_y'）→ 局部（R_z）
    % 逆变换链：局部 → 中间（R_z_inv）→ 全局（R_y）
    [n_local_frames, n_rot_angles, n_points, ~] = size(B_local);
    B_global = zeros(size(B_local));
    for k = 1:n_local_frames
        cos_z = cos(local_angles_rad(k));
        sin_z = sin(local_angles_rad(k));
        R_z_inv = [cos_z,  -sin_z, 0;   % 绕Z轴旋转的逆矩阵（转置）
                  sin_z,  cos_z, 0;
                   0,      0,     1];
        for i = 1:n_rot_angles
            theta_y = rot_angles_y_rad(i);
            cos_y = cos(theta_y);
            sin_y = sin(theta_y);
            % 修正：使用 R_y 正向旋转矩阵而非 R_y 的逆矩阵
            % 逆变换应为 B_global = R_y * R_z_inv * B_local
            R_y_fwd = [cos_y,  0, -sin_y;
                       0,     1,  0;
                       sin_y, 0,  cos_y];
            for j = 1:n_points
                B_local_vec = squeeze(B_local(k, i, j, :));
                % 修正：先用 R_z_inv 撤销局部旋转，再用 R_y（正向）旋转回全局
                B_intermediate = R_z_inv * B_local_vec;
                B_global_vec = R_y_fwd * B_intermediate;
                B_global(k, i, j, :) = B_global_vec';
            end
        end
    end
end

function plotMagneticFieldResults(B_total, total_length, n_points)
    cumulative_length = linspace(0, total_length, n_points);
    Bz = B_total(:, 3);
    figure('Color', 'white', 'Position', [100, 100, 900, 600]);
    plot(cumulative_length, Bz, 'b-', 'LineWidth', 1.5, ...
         'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    xlabel('累计曲线长度 (m)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('B_z 磁场强度 (T)', 'FontSize', 12, 'FontWeight', 'bold');
    title('B_z磁场强度沿线圈曲线累计长度的变化（完全反变换后求和）', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor; box on;
    set(gca, 'FontSize', 10);
    legend('B_z 分量', 'Location', 'best');
    [Bz_max, idx_max] = max(Bz);
    [Bz_min, idx_min] = min(Bz);
    text(cumulative_length(idx_max), Bz_max, ...
        sprintf('  Max: %.3e T\n  at %.2f m', Bz_max, cumulative_length(idx_max)), ...
        'FontSize', 9, 'Color', 'red', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
    text(cumulative_length(idx_min), Bz_min, ...
        sprintf('  Min: %.3e T\n  at %.2f m', Bz_min, cumulative_length(idx_min)), ...
        'FontSize', 9, 'Color', 'blue', 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    fprintf('\n========== 磁场计算结果统计（完全反变换后求和）==========\n');
    fprintf('B_z 最大值: %.4e T (位置: %.2f m)\n', Bz_max, cumulative_length(idx_max));
    fprintf('B_z 最小值: %.4e T (位置: %.2f m)\n', Bz_min, cumulative_length(idx_min));
    fprintf('B_z 平均值: %.4e T\n', mean(Bz));
    fprintf('B_z 标准差: %.4e T\n', std(Bz));
    % 输出Bx和By的统计信息（验证反变换正确性）
    Bx = B_total(:, 1);
    By = B_total(:, 2);
    fprintf('\nBx 平均值: %.4e T（理论上应接近0）\n', mean(Bx));
    fprintf('By 平均值: %.4e T（理论上应接近0）\n', mean(By));
    fprintf('Bx 最大值: %.4e T\n', max(abs(Bx)));
    fprintf('By 最大值: %.4e T\n', max(abs(By)));
end
