function pbb_g = to_global_pbb(pbb_c, R_global, t_global)
% pbb_c:   Nx6  [x y z vx vy vz] in local/sensor frame
% R_global: 3x3 or 3x3xN rotation (local -> global)
% t_global: 3x1 or 1x3 or 3xN translation (local -> global)
% returns: Nx6 in global frame

N = size(pbb_c,1);
pos = pbb_c(:,1:3);    % Nx3
vel = pbb_c(:,4:6);    % Nx3

if ndims(R_global)==2
    % one pose for all detections
    R = R_global;
    t = t_global(:).';                % 1x3
    pos_g = pos*R.' + repmat(t, N, 1);
    vel_g = vel*R.';
else
    % per-detection pose
    pos_g = zeros(N,3);
    vel_g = zeros(N,3);
    for i = 1:N
        R = R_global(:,:,i);
        if size(t_global,2)==N, t = t_global(:,i);
        else,                    t = t_global; % 3x1
        end
        pos_g(i,:) = (R*pos(i,:).').'+ t(:).';
        vel_g(i,:) = (R*vel(i,:).').';
    end
end

pbb_g = [pos_g, vel_g];
end
