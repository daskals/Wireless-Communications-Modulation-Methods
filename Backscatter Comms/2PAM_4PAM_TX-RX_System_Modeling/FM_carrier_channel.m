function [x_FM] = FM_carrier_channel(total_packet_duration, chanParams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis 
% 26/12/2017
% Last Version: 23/1/2018
% Email: daskalakispiros@gmail.com
% Website: www.daskalakispiros.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation code, Rev. 1
% John Kimionis, Sep.13, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_sampling = chanParams.framelength*total_packet_duration;    % Sampling time frame (seconds).
t = 0:chanParams.Ts:t_sampling-chanParams.Ts;

% carrier
   m = cos(2*pi*chanParams.F_sin*t);      % baseband sine
   x_FM = NaN*ones(size(m));              % modulated FM signal

  for tt=1:length(t)
     x_FM(tt) = exp(j*2*pi*(chanParams.f_c+chanParams.CFO)*t(tt) + j*2*pi*chanParams.Df*sum(m(1:tt))*chanParams.Ts);
  end


end




