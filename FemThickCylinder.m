%=========================================================================%
% CHƯƠNG TRÌNH MATLAB GIẢI BÀI TOÁN ỐNG TRỤ DÀY (PLANE STRAIN)
clear;clc;close all; 

%% ================== BLOCK 1: THÔNG SỐ ĐẦU VÀO (GIAO DIỆN UI) ===========================

% --- Tạo hộp thoại nhập liệu (Input Dialog) ---
prompt = {'Mô đun đàn hồi E (Pa):', 'Hệ số Poisson \nu:', 'Bán kính trong Ri (m):', ...
          'Bán kính ngoài Ro (m):', 'Áp suất trong p_i (Pa):', 'Áp suất ngoài p_o (Pa):',...
          'Kích thước phần tử (Mesh Sizing - m):', ...
          'Loại phần tử muốn vẽ (1: Tứ giác Q4, 2: Tam giác T3):'};
dlgtitle = 'Nhap thong so - 360 do & 1/4 Ong tru'; % Tiêu đề hộp thoại
dims = 1; % Số dòng cho mỗi ô nhập liệu

% Đặt các giá trị mặc định có sẵn vào hộp thoại (Người dùng có thể sửa)
definput = {'200e9', '0.3', '0.5', '1.0', '2e7', '1e7', '0.05', '1'}; 
answer = inputdlg(prompt, dlgtitle, dims, definput); % Hiển thị UI và hứng kết quả vào biến answer

% Nếu người dùng bấm "Cancel" (hủy) hoặc đóng cửa sổ, dừng chương trình
if isempty(answer), disp('Da huy qua trinh nhap lieu!'); return; end

% --- Chuyển đổi kiểu dữ liệu từ chuỗi (string) sang số (double) ---
E = str2double(answer{1});       % Mô đun đàn hồi E (Vật liệu thép khoảng 200 GPa)
nu = str2double(answer{2});      % Hệ số Poisson
Ri = str2double(answer{3});      % Bán kính trong
Ro = str2double(answer{4});      % Bán kính ngoài
pi_p = str2double(answer{5});    % Áp suất tác dụng lên mặt trong
po_p = str2double(answer{6});    % Áp suất tác dụng lên mặt ngoài
sizing = str2double(answer{7});  % Kích thước mục tiêu của một phần tử lưới
plot_choice = round(str2double(answer{8})); % Lựa chọn 1 (Q4) hoặc 2 (T3)

% --- Xây dựng Ma trận vật liệu đàn hồi [D] ---
% Bài toán ống trụ dài chịu áp suất là bài toán Biến dạng phẳng (Plane Strain)
% Công thức ma trận [D] đặc trưng cho liên hệ Ứng suất - Biến dạng trong Plane Strain
D = (E/((1+nu)*(1-2*nu))) * [1-nu, nu, 0; 
                             nu, 1-nu, 0; 
                             0, 0, (1-2*nu)/2];
                             
% --- Chuẩn bị dữ liệu cho Khảo sát hội tụ (Mesh Convergence) ---
% Tính số lượng phần tử N nằm trên bề dày (từ Ri đến Ro) dựa vào kích thước sizing
N_user = max(2, round((Ro - Ri) / sizing)); 

% Tạo một mảng Ns chứa các cấp độ chia lưới khác nhau (từ thô đến mịn) để vẽ đồ thị hội tụ
Ns = unique([max(2, round(N_user/4)), max(4, round(N_user/2)), N_user, N_user*2]); 

% Khởi tạo mảng chứa sai số bằng 0
errQ = zeros(1, length(Ns)); % Sai số của phần tử Tứ giác Q4
errT = zeros(1, length(Ns)); % Sai số của phần tử Tam giác T3
%% ================== BLOCK 2: KHẢO SÁT HỘI TỤ (GIẢI MÔ HÌNH 1/4) =======================
fprintf('Dang tinh toan mo hinh 1/4 ong tru...\n');
for idx = 1:length(Ns) % Duyệt qua từng cấp độ chia lưới trong mảng Ns
    N = Ns(idx); % Số phần tử trên bề dày ở bước lặp hiện tại
    
    % 1. Sinh lưới 1/4 (Tham số cuối là 0). Trả về tọa độ Nút (nodes) và kết nối Nút (eQ, eT)
    [nodes, eQ, eT] = mesh_gen(N, Ri, Ro, 0); 
    
    % 2. Gọi hàm giải bài toán FEM cho Q4 và T3
    [UQ, SQ] = solveFEM(nodes, eQ, 'Q', D, Ri, Ro, pi_p, po_p, 0); % Giải hệ Q4
    [UT, ST] = solveFEM(nodes, eT, 'T', D, Ri, Ro, pi_p, po_p, 0); % Giải hệ T3
    
    % 3. So sánh với nghiệm giải tích để tìm sai số
    % Tìm index (số thứ tự) của các nút nằm trên biên mặt trong (Bán kính = Ri)
    in_n = find(abs(sqrt(nodes(:,1).^2+nodes(:,2).^2)-Ri)<1e-5); 
    
    % Tính góc theta của các nút này (dùng để chiếu hệ tọa độ Decartes sang hệ tọa độ Cực)
    th = atan2(nodes(in_n,2), nodes(in_n,1));                                  
    
    % Tính chuyển vị hướng kính ur = ux*cos(theta) + uy*sin(theta)
    urQ = UQ(2*in_n-1).*cos(th) + UQ(2*in_n).*sin(th); % Của FEM Q4
    urT = UT(2*in_n-1).*cos(th) + UT(2*in_n).*sin(th); % Của FEM T3
    
    % Lấy nghiệm chuyển vị đúng (chính xác) từ phương trình giải tích Lame
    [ur_ex, ~, ~] = exact(Ri, Ri, Ro, pi_p, po_p, E, nu);
    
    % Tính phần trăm sai số trung bình (Sai số tương đối)
    errQ(idx) = mean(abs(urQ - ur_ex)/abs(ur_ex))*100;
    errT(idx) = mean(abs(urT - ur_ex)/abs(ur_ex))*100;
    
    % 4. Khi vòng lặp chạy đến độ mịn mà người dùng yêu cầu (N_user), lưu trữ lại kết quả để vẽ
    if N == N_user
        n_14=nodes; eQ_14=eQ; eT_14=eT;     % Lưu cấu trúc lưới
        UQ_14=UQ; SQ_14=SQ; UT_14=UT; ST_14=ST; % Lưu Chuyển vị (U) và Ứng suất (S)
    end 
end

%% ================== BLOCK 3: TÍNH TOÁN THÊM MÔ HÌNH FULL 360 ĐỘ =====================
fprintf('Dang tinh toan them mo hinh 360 do...\n');
% Sinh lưới và giải bài toán cho cấu trúc 360 độ (Tham số cuối = 1)
[n_360, eQ_360, eT_360] = mesh_gen(N_user, Ri, Ro, 1); 
[UQ_360, SQ_360] = solveFEM(n_360, eQ_360, 'Q', D, Ri, Ro, pi_p, po_p, 1);
[UT_360, ST_360] = solveFEM(n_360, eT_360, 'T', D, Ri, Ro, pi_p, po_p, 1);
fprintf('Hoan thanh!\n');

%% ================== BLOCK 4: VẼ CÁC ĐỒ THỊ VÀ HÌNH ẢNH YÊU CẦU ======================

% Gán biến tương ứng tùy thuộc người dùng muốn xem Q4 hay T3
is_T3 = (plot_choice == 2); % Trả về logic True (1) hoặc False (0)
if is_T3
    elem_360 = eT_360; U_360 = UT_360; S_360 = ST_360;
    elem_14 = eT_14; U_14 = UT_14; S_14 = ST_14;
    elem_name = 'T3 (Tam giac)';
else
    elem_360 = eQ_360; U_360 = UQ_360; S_360 = SQ_360;
    elem_14 = eQ_14; U_14 = UQ_14; S_14 = SQ_14;
    elem_name = 'Q4 (Tu giac)';
end

% Hình 1: Đồ thị khảo sát hội tụ lưới
figure('Name','Khao sat hoi tu');
plot(Ns, errQ, '-ob', Ns, errT, '-sr', 'LineWidth', 2); % Vẽ 2 đường sai số
legend('Q4 (Tu giac)','T3 (Tam giac)'); 
title('Sai so chuyen vi u_r (%)'); xlabel('So phan tu (N)'); grid on;

% Hình 2: Vẽ minh họa hình học ống và cấu trúc lưới phân tử đã chia
figure('Name','Mo hinh & Luoi (Full 360)');
subplot(1,2,1); hold on; th_plot = linspace(0, 2*pi, 200); % Tạo mảng góc 0->360 độ
plot(Ri*cos(th_plot), Ri*sin(th_plot), 'b', Ro*cos(th_plot), Ro*sin(th_plot), 'r', 'LineWidth', 2); % Vẽ vòng tròn trong và ngoài
axis equal; title('Mo hinh 360 do'); grid on;
subplot(1,2,2); 
% Dùng lệnh patch để vẽ các đa giác liên kết tạo thành lưới
patch('Faces',elem_360,'Vertices',n_360,'FaceColor','w', 'EdgeColor', [0.2 0.2 0.2]); 
axis equal; title(['Luoi phan tu ', elem_name]);

% Hình 3: Gọi hàm vẽ 6 biểu đồ nhiệt (contour) cho ống FULL 360 độ
draw_6_plots(n_360, elem_360, U_360, S_360, nu, ['Truong Ung suat/Chuyen vi 360 do - ', elem_name], is_T3);

% Hình 4: Gọi hàm vẽ 6 biểu đồ nhiệt cho 1/4 Ống trụ
draw_6_plots(n_14, elem_14, U_14, S_14, nu, ['Truong Ung suat/Chuyen vi 1/4 ong - ', elem_name], is_T3);

% Hình 5: Đồ thị so sánh kết quả dọc theo phương bán kính giữa FEM và Giải tích
% Tìm các nút nằm trên trục x (có tung độ y = 0) của mô hình 1/4
n_y0 = find(n_14(:,2)<1e-5); 
[r_y0, idx] = sort(n_14(n_y0,1)); n_y0 = n_y0(idx); % Sắp xếp thứ tự từ trong ra ngoài

% Lấy kết quả giải tích chuẩn Lame tại các bán kính này
[ur_ex, sr_ex, st_ex] = exact(r_y0, Ri, Ro, pi_p, po_p, E, nu);

% Trích xuất chuyển vị của FEM (chuyển sang hệ tọa độ cực)
th_nodes = atan2(n_14(:,2), n_14(:,1)); 
ux = U_14(1:2:end); uy = U_14(2:2:end); 
ur14 = ux.*cos(th_nodes) + uy.*sin(th_nodes); 

% Trích xuất ứng suất dọc mặt cắt trục ngang để vẽ so sánh
if is_T3
    % Đối với T3, ứng suất được lưu ở TRỌNG TÂM phần tử, nên ta tìm tọa độ tâm (cx, cy)
    cx = (n_14(elem_14(:,1),1) + n_14(elem_14(:,2),1) + n_14(elem_14(:,3),1))/3;
    cy = (n_14(elem_14(:,1),2) + n_14(elem_14(:,2),2) + n_14(elem_14(:,3),2))/3;
    th_c = atan2(cy, cx);
    
    % Tìm lớp phần tử nằm sát sát trục X nhất (góc theta rất nhỏ)
    el_y0 = find(th_c < pi/(2*N_user)); 
    [r_fem, idx] = sort(cx(el_y0)); el_y0 = el_y0(idx); % Sắp xếp theo bán kính r
    
    % Xoay tensor ứng suất từ hệ x-y sang hệ r-theta
    sr14 = S_14(:,1).*cos(th_c).^2 + S_14(:,2).*sin(th_c).^2 + 2*S_14(:,3).*sin(th_c).*cos(th_c);
    st14 = S_14(:,1).*sin(th_c).^2 + S_14(:,2).*cos(th_c).^2 - 2*S_14(:,3).*sin(th_c).*cos(th_c);
    sr_fem = sr14(el_y0); st_fem = st14(el_y0);
else
    % Đối với Q4, ứng suất đã được tính ở từng NÚT
    r_fem = r_y0;
    sr14 = S_14(:,1).*cos(th_nodes).^2 + S_14(:,2).*sin(th_nodes).^2 + 2*S_14(:,3).*sin(th_nodes).*cos(th_nodes);
    st14 = S_14(:,1).*sin(th_nodes).^2 + S_14(:,2).*cos(th_nodes).^2 - 2*S_14(:,3).*sin(th_nodes).*cos(th_nodes);
    sr_fem = sr14(n_y0); st_fem = st14(n_y0);
end

% Vẽ 3 đồ thị con so sánh: Chuyển vị, Ứng suất hướng kính (Radial), Ứng suất vòng (Hoop)
figure('Name','So sanh ket qua Giải tich va Số (FEM)');
subplot(1,3,1); plot(r_y0, ur_ex, 'k-', r_y0, ur14(n_y0), 'ro'); title('Chuyen vi u_r'); legend('Giai tich', ['FEM (', elem_name, ')']); grid on;
subplot(1,3,2); plot(r_y0, sr_ex, 'k-', r_fem, sr_fem, 'ro'); title('Ung suat \sigma_r'); grid on;
subplot(1,3,3); plot(r_y0, st_ex, 'k-', r_fem, st_fem, 'ro'); title('Ung suat \sigma_\theta'); grid on;


%% ========================================================================
%  ========================= CÁC HÀM PHỤ TRỢ (FUNCTIONS) ==================
%  ========================================================================

%% Hàm 0: Hàm hiển thị phổ màu (Biểu đồ nhiệt - Contour map)
function draw_6_plots(nodes, elem, U, S, nu, fig_name, is_T3)
    % 1. XỬ LÝ CHUYỂN VỊ: Chuyển vị liên tục từ nút nọ sang nút kia
    th_nodes = atan2(nodes(:,2), nodes(:,1)); % Tính góc theta
    ux = U(1:2:end); uy = U(2:2:end); % Lọc chuyển vị Ux (nút lẻ) và Uy (nút chẵn)
    ur = ux.*cos(th_nodes) + uy.*sin(th_nodes); % Tính chuyển vị hướng kính
    
    % 2. XỬ LÝ ỨNG SUẤT: Gom chung về Tâm phần tử để ép hiển thị phổ màu rỗ hạt (Flat)
    if is_T3
        % T3 có sẵn 1 giá trị ứng suất/phần tử
        cx = (nodes(elem(:,1),1) + nodes(elem(:,2),1) + nodes(elem(:,3),1))/3;
        cy = (nodes(elem(:,1),2) + nodes(elem(:,2),2) + nodes(elem(:,3),2))/3;
        S_plot = S; 
    else
        % Q4 lấy trung bình ứng suất của 4 Nút để quy về Trọng tâm phần tử
        cx = (nodes(elem(:,1),1) + nodes(elem(:,2),1) + nodes(elem(:,3),1) + nodes(elem(:,4),1))/4;
        cy = (nodes(elem(:,1),2) + nodes(elem(:,2),2) + nodes(elem(:,3),2) + nodes(elem(:,4),2))/4;
        S_plot = (S(elem(:,1),:) + S(elem(:,2),:) + S(elem(:,3),:) + S(elem(:,4),:)) / 4;
    end
    
    th_S = atan2(cy, cx); % Góc theta tại tâm phần tử
    
    % Xoay ma trận ứng suất sang tọa độ cực r-theta (Mohr's Circle)
    sr = S_plot(:,1).*cos(th_S).^2 + S_plot(:,2).*sin(th_S).^2 + 2*S_plot(:,3).*sin(th_S).*cos(th_S);
    st = S_plot(:,1).*sin(th_S).^2 + S_plot(:,2).*cos(th_S).^2 - 2*S_plot(:,3).*sin(th_S).*cos(th_S);
    sz = nu*(S_plot(:,1)+S_plot(:,2)); % Ứng suất dọc trục z trong bài toán Plane Strain
    
    % Tính ứng suất tương đương Von Mises (Tiêu chí chảy dẻo)
    svm = sqrt(0.5*((S_plot(:,1)-S_plot(:,2)).^2 + (S_plot(:,2)-sz).^2 + (sz-S_plot(:,1)).^2 + 6*S_plot(:,3).^2));
    
    % Đóng gói dữ liệu để vẽ tự động
    vars = {ur, sr, st, svm, ux, uy};
    tits = {'Chuyen vi u_r', 'Ung suat \sigma_r', 'Ung suat \sigma_\theta', 'Ung suat Von Mises', 'Chuyen vi u_x', 'Chuyen vi u_y'};
    
    figure('Name', fig_name); colormap('jet'); % Bảng màu cầu vồng (jet)
    for i=1:6
        subplot(2,3,i); % Tạo 6 khung ảnh phụ
        
        % Cài đặt cách tô màu:
        % Nếu là Ứng suất (hình 2,3,4) thì ép kiểu 'flat' (hiển thị từng ô rời rạc)
        if i == 2 || i == 3 || i == 4
            fcolor = 'flat'; ecolor = [0.6 0.6 0.6]; % ecolor là màu viền xám
        else
            % Nếu là Chuyển vị (1, 5, 6) thì hiển thị dải màu mượt 'interp'
            fcolor = 'interp'; ecolor = 'none';   
        end
        
        % Lệnh patch vẽ các đa giác, gán dữ liệu màu (FaceVertexCData)
        patch('Faces',elem,'Vertices',nodes,'FaceVertexCData',vars{i},'FaceColor',fcolor,'EdgeColor',ecolor);
        colorbar; title(tits{i}); axis equal off; % Ẩn trục tọa độ
    end
end

%% Hàm 1: Sinh lưới tọa độ và xây dựng ma trận kết nối (Mesh Generator)
function [nodes, eQ, eT] = mesh_gen(N, Ri, Ro, is_full)
    % Tính toán số vòng lặp theo phương góc (theta)
    if is_full, Nth = 4*N; th_max = 2*pi*(1 - 1/Nth); num_th = Nth; 
    else,       Nth = N;   th_max = pi/2;             num_th = Nth + 1; end 
    
    % Hàm ndgrid tạo lưới lưới điểm (r, theta) đều đặn
    [R, T] = ndgrid(linspace(Ri,Ro,N+1), linspace(0, th_max, num_th));
    nodes = [R(:).*cos(T(:)), R(:).*sin(T(:))]; % Chuyển sang tọa độ Descartes (x, y)
    
    % Xây dựng danh sách các Nút (Connectivity) thuộc về 1 phần tử
    [I, J] = ndgrid(1:N, 1:Nth);
    n1 = I(:) + (J(:)-1)*(N+1); % Góc dưới bên trái
    n2 = n1 + 1;                % Góc dưới bên phải
    if is_full, n4 = I(:) + mod(J(:), Nth)*(N+1); % Cấu trúc 360 độ cần nối điểm cuối vòng lại điểm đầu
    else,       n4 = I(:) + J(:)*(N+1); end       % Góc trên bên trái
    n3 = n4 + 1;                % Góc trên bên phải
    
    eQ = [n1, n2, n3, n4]; % Phân tử Tứ giác Q4 được định hình bởi 4 nút này
    
    % Từ phần tử Q4, chẻ đôi theo đường chéo tạo thành 2 phần tử T3
    eT = [eQ(:,1:3); eQ(:,3), eQ(:,4), eQ(:,1)]; 
end

%% Hàm 2: Thuật toán Giải FEM (Lắp ráp & Giải hệ phương trình)
function [U, S] = solveFEM(nodes, elem, type, D, Ri, Ro, pi_p, po_p, is_full)
    nn = size(nodes,1);      % Tổng số nút
    F = zeros(2*nn,1);       % Khởi tạo Vector Tải trọng Toàn cục (Global Force Vector)
    
    % --- Bước 1: Quy đổi Áp suất (Distributed Load) thành Lực tại các Nút (Nodal Force) ---
    for pass=1:2 % Chạy 2 lần: Lần 1 mặt trong, Lần 2 mặt ngoài
        if pass==1, r_t=Ri; p=pi_p; sn=1; else, r_t=Ro; p=po_p; sn=-1; end
        
        % Tìm các nút nằm trên bề mặt chịu áp suất
        bn = find(abs(sqrt(nodes(:,1).^2+nodes(:,2).^2)-r_t)<1e-5);
        [~, id] = sort(atan2(nodes(bn,2), nodes(bn,1))); bn = bn(id); % Xếp theo thứ tự vòng tròn
        n_seg = length(bn) - 1 + is_full; 
        
        % Chiếu áp suất thành phần lực Fx, Fy dựa theo độ dài và phương của các đoạn biên (segment)
        for i=1:n_seg
            n1=bn(i); if i < length(bn), n2 = bn(i+1); else, n2 = bn(1); end
            dx = nodes(n2,1)-nodes(n1,1); dy = nodes(n2,2)-nodes(n1,2);
            F([2*n1-1,2*n1,2*n2-1,2*n2]) = F([2*n1-1,2*n1,2*n2-1,2*n2]) + p*sn*[dy;-dx;dy;-dx]/2;
        end
    end
    
    % --- Bước 2: Tính Ma trận Độ cứng [K] ---
    % Ứng dụng kỹ thuật Triplet (I, J, V) để rải ma trận Sparse cực nhanh trong MATLAB
    ne = size(elem, 1); nd = length(elem(1,:)) * 2; % nd là số bậc tự do 1 phần tử (T3=6, Q4=8)
    I_idx = zeros(ne*nd^2, 1); J_idx = zeros(ne*nd^2, 1); V_val = zeros(ne*nd^2, 1); idx = 1;
    
    for i=1:ne
        en = elem(i,:); % Danh sách nút của phần tử thứ i
        ed = reshape([2*en-1; 2*en], 1, []); % Danh sách các Bậc tự do (Degrees of Freedom)
        
        % Gọi hàm tính Ma trận độ cứng phần tử [Ke]
        if type=='Q', Ke = KeQ(nodes(en,:),D); else, Ke = KeT(nodes(en,:),D); end
        
        % Lắp ráp [Ke] vào hệ tọa độ mảng 1 chiều (Triplet)
        ed_I = repmat(ed', 1, nd); ed_J = repmat(ed, nd, 1); nv = numel(Ke);
        I_idx(idx:idx+nv-1) = ed_I(:); J_idx(idx:idx+nv-1) = ed_J(:); V_val(idx:idx+nv-1) = Ke(:); idx = idx + nv;
    end
    % Tổng hợp thành Ma trận Độ cứng Toàn cục [K]
    K = sparse(I_idx(1:idx-1), J_idx(1:idx-1), V_val(1:idx-1), 2*nn, 2*nn);
    
    % --- Bước 3: Áp dụng Điều kiện biên (Boundary Conditions) và Giải Phương trình ---
    % Mô hình 1/4 bị khóa đối xứng: Nút trên trục y bị khóa phương x, nút trên trục x bị khóa phương y
    fixed = [2*find(abs(nodes(:,1))<1e-5)-1; 2*find(abs(nodes(:,2))<1e-5)];
    free = setdiff(1:2*nn, fixed); % Lấy những bậc tự do (nút) được phép tự do dịch chuyển
    
    U = zeros(2*nn,1); 
    U(free) = K(free,free)\F(free); % GIẢI HỆ PHƯƠNG TRÌNH ĐẠI SỐ TUYẾN TÍNH: [K]{U} = {F}
    
    % --- Bước 4: Trích xuất Ứng suất (Post-Processing) từ Chuyển vị [U] ---
    if type=='Q'
        % Phần tử Q4: Phải tính ứng suất tại các Điểm Gauss, 
        % sau đó ngoại suy ra Nút và lấy trung bình các phần tử dùng chung 1 nút.
        S = zeros(nn,3); count = zeros(nn,1);
        for i=1:ne 
            en = elem(i,:); ed = reshape([2*en-1; 2*en], 1, []);
            % Ma trận đạo hàm của Hàm dạng (Shape functions)
            dN = 0.25*[-1 1 1 -1; -1 -1 1 1]; J = dN*nodes(en,:); dNx = J\dN; B = zeros(3,8);
            B(1,1:2:7)=dNx(1,:); B(2,2:2:8)=dNx(2,:); B(3,1:2:7)=dNx(2,:); B(3,2:2:8)=dNx(1,:);
            S(en,:) = S(en,:) + repmat((D*B*U(ed))',length(en),1); % S = [D]*[B]*{Ue}
            count(en) = count(en) + 1;
        end
        S = S ./ count; % Chia trung bình
    else
        % Phần tử CST (T3): Ứng suất và Biến dạng là một HẰNG SỐ trên toàn bộ diện tích tam giác
        S = zeros(ne,3);
        for i=1:ne
            en = elem(i,:); ed = reshape([2*en-1; 2*en], 1, []);
            A = 0.5*det([1 nodes(en(1),:); 1 nodes(en(2),:); 1 nodes(en(3),:)]); % Diện tích tam giác
            b = [nodes(en(2),2)-nodes(en(3),2), nodes(en(3),2)-nodes(en(1),2), nodes(en(1),2)-nodes(en(2),2)];
            c = [nodes(en(3),1)-nodes(en(2),1), nodes(en(1),1)-nodes(en(3),1), nodes(en(2),1)-nodes(en(1),1)];
            % Xây dựng Ma trận biến dạng [B] (B Matrix)
            B = 1/(2*A)*[b(1) 0 b(2) 0 b(3) 0; 0 c(1) 0 c(2) 0 c(3); c(1) b(1) c(2) b(2) c(3) b(3)];
            S(i,:) = (D*B*U(ed))'; % Ứng suất Element i = [D]*[B]*{Ue}
        end
    end
end

%% Hàm 3: Tính Ma trận Độ cứng [Ke] cho phần tử Tứ giác Q4
function Ke = KeQ(C, D)
    Ke = zeros(8,8); 
    gp = [-0.57735, 0.57735]; % Tọa độ các Điểm tích phân Gauss (Gauss Quadrature 2x2 = 4 điểm)
    for xi = gp
        for eta = gp
            % Đạo hàm hàm dạng bậc 1
            dN = 0.25*[-(1-eta) 1-eta 1+eta -(1+eta); -(1-xi) -(1+xi) 1+xi 1-xi];
            J = dN*C; % Ma trận Jacobian (Ánh xạ tọa độ từ hệ quy chiếu tọa độ tự nhiên sang thực tế)
            dNx = J\dN; B = zeros(3,8); 
            B(1,1:2:7) = dNx(1,:); B(2,2:2:8) = dNx(2,:); B(3,1:2:7) = dNx(2,:); B(3,2:2:8) = dNx(1,:);
            
            % Tích phân biểu thức [Ke] = Tích phân([B]^T * [D] * [B] * det(J)) qua 4 điểm Gauss
            Ke = Ke + B'*D*B*det(J); 
        end
    end
end

%% Hàm 4: Tính Ma trận Độ cứng [Ke] cho phần tử Tam giác CST (T3)
function Ke = KeT(C, D)
    A = 0.5*det([1 C(1,:); 1 C(2,:); 1 C(3,:)]); % Diện tích tam giác
    % Các hệ số hình học b, c
    b = [C(2,2)-C(3,2), C(3,2)-C(1,2), C(1,2)-C(2,2)]; 
    c = [C(3,1)-C(2,1), C(1,1)-C(3,1), C(2,1)-C(1,1)];
    
    % Ma trận B là một hằng số, không phụ thuộc vào x, y
    B = 1/(2*A)*[b(1) 0 b(2) 0 b(3) 0; 0 c(1) 0 c(2) 0 c(3); c(1) b(1) c(2) b(2) c(3) b(3)];
    
    % Vì [B] và [D] là hằng số, tích phân suy biến thành Nhân với diện tích A
    Ke = B'*D*B*A; 
end

%% Hàm 5: Nghiệm Giải tích Chính xác (Exact Lame's Solution)
% Phương trình cổ điển của Lame dành cho bài toán Ống trụ thành dày
function [ur, sr, st] = exact(r, Ri, Ro, pi_p, po_p, E, nu)
    % Các hằng số C1, C2 theo điều kiện áp suất trong (pi) và ngoài (po)
    C1 = (pi_p*Ri^2 - po_p*Ro^2)/(Ro^2-Ri^2); 
    C2 = (pi_p - po_p)*Ri^2*Ro^2/(Ro^2-Ri^2);
    
    sr = C1 - C2./r.^2; % Ứng suất hướng kính (Radial Stress)
    st = C1 + C2./r.^2; % Ứng suất vòng (Hoop Stress)
    ur = (1+nu)/E*((1-2*nu)*C1.*r + C2./r); % Chuyển vị hướng kính (Radial Displacement)
end