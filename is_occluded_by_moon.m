function occ = is_occluded_by_moon(x_tgt, x_obs, mu)
%IS_OCCLUDED_BY_MOON Is 'tgt' occluded by the moon when observed by 'obs'
% vectors should be defined in synodic frame.
% very simple geometric occlusion check using relative angles bw position
% vectors
moon_radius = 1738/384400; % in DU

% decompose vectors
r_tgt = x_tgt(1:3);
r_obs = x_obs(1:3);

%relative position
r_tgt_obs = r_tgt-r_obs;
d_tgt_obs = norm(r_tgt_obs); % dist

% moon position wrt observer
r_lun = [1.0-mu;0;0];
r_lun_obs = r_lun - r_obs;
d_lun_obs = norm(r_lun_obs);

% angles
theta_obs_moon = atan(moon_radius/d_lun_obs); % angle of occluded zone
theta_obs_tgt = acos(dot(r_lun_obs, r_tgt_obs)/(d_lun_obs*d_tgt_obs));

occ = false;
if theta_obs_moon >= theta_obs_tgt && d_tgt_obs > d_lun_obs
    occ = true;
end

end

