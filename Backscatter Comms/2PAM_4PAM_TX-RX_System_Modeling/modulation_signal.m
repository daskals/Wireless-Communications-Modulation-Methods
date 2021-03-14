%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis 
% 26/12/2017
% Last Version: 23/1/2018
% Email: daskalakispiros@gmail.com
% Website: www.daskalakispiros.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r] = modulation_signal(x_tag, x_FM,chanParams, t)


 if chanParams.Modtype==1
        %% preset reflection coefficients
        G1 =-0.7245 - j*0.6922;
        G2 =-0.3414 - j*0.2881;
        G3 =+0.0223 + j*0.1779;
        G4 =+0.3079 + j*0.6334;

%         G1 = 0.001*exp(j*unifrnd(0, 2*pi));
%         G2 = 0.33*exp(j*angle(G1));
%         G3 = 0.66*exp(j*angle(G1));
%         G4 = 1*exp(j*angle(G1));
 
        x_tag_complex = x_tag;
        % (-3 || 00)----(-1 || 01)------(1 || 11)-----(3 || 10)
        x_tag_complex(find(x_tag_complex==-3)) = G1 ;
        x_tag_complex(find(x_tag_complex==-1)) = G2 ;
        x_tag_complex(find(x_tag_complex==1)) = G3 ;
        x_tag_complex(find(x_tag_complex==3)) = G4 ;
        
 elseif chanParams.Modtype==0
     
        % random reflection coefficients
        % G0 = unifrnd(0, 1)*exp(j*unifrnd(0, 2*pi));
        % G1 = unifrnd(0, 1)*exp(j*unifrnd(0, 2*pi));

        %% antipodal reflection coefficients
        G0 = 0.5*exp(j*unifrnd(0, 2*pi));
        G1 = 0.5*exp(j*angle(G0) + j*pi);

        %% preset reflection coefficients
        % G0 = 0.95*exp(j*pi/3);
        % G1 = 0.7*exp(j*pi/3 + j*pi + j*pi/10);

        x_tag_complex = x_tag;
        indx = find(x_tag_complex==1);
        x_tag_complex(1:end) = G0;
        x_tag_complex(indx) = G1;
        
 end

  A = sqrt(chanParams.P_TX)*exp(j*unifrnd(0, 2*pi));   % carrier received at SDR
  B = 0.1*exp(j*unifrnd(0, 2*pi));          % modulated signal scaling
 
  %L_padding = 20000;
  %% Random L_padding
  L_padding= randi([0 (length(x_FM)-length(x_tag_complex))], 1,1);
  x_mul=[G1*ones(1, L_padding) x_tag_complex  G1*ones(1, length(x_FM)-length(x_tag_complex)-L_padding)];
  r = A*x_FM + A*B*x_FM.*x_mul;             % full signal
  
  r = awgn(r,chanParams.noiseSNR);
  
  %% Add a small flactuation in our final modulated signal
  F_sin_flactuation = 0.5*rand(1);              % rate of sine
  flactuation = cos(2*pi*F_sin_flactuation*t);     % baseband sine
  r=r.*flactuation;

end

