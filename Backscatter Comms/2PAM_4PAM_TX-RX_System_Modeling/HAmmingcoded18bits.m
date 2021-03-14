%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis 
% 26/12/2017
% Last Version: 23/1/2018
% Email: daskalakispiros@gmail.com
% Website: www.daskalakispiros.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function insert_parity_bits calculates and inserts the parity bits
function [newmes] = HAmmingcoded18bits(message)

        n=length(message);
        
        % mathematical formula that calculates the number of parity bits
        nbp=floor(log2(n+ceil(log2(n))))+1;
        
        %parity bits 1,2,4,8,16
        newmes=zeros(1,n+nbp+1);
        mes_index=1;
        
        for i=1:1:n+nbp
           if (i==1 || i==2 || i==4 || i==8 ||i==16)
              newmes(i)=-1;
           else
              newmes(i)=message(mes_index);
              mes_index=mes_index+1;
           end
        end
     
        %%
        par1=[newmes(17) newmes(15) newmes(13) newmes(11) newmes(9) newmes(7) newmes(5) newmes(3)];
        par1= mod(sum (par1),2);
        
        %%
        par2=[ newmes(18)  newmes(15)  newmes(14)  newmes(11) newmes(10) newmes(7) newmes(6) newmes(3)];
        par2= mod(sum (par2),2);
       
        %%
        par4=[newmes(15) newmes(14)  newmes(13) newmes(12) newmes(7) newmes(6) newmes(5) ];
        par4= mod(sum (par4),2);
        
        %%
        par8=[ newmes(15) newmes(14)  newmes(13) newmes(12) newmes(11) newmes(10) newmes(9)];
        par8= mod(sum (par8),2);
        
        %%
        par16=[newmes(18)  newmes(17) ];
        par16= mod(sum (par16),2);
        
        %%
        newmes (1)=par1;
        newmes (2)=par2;
        newmes (4)=par4;
        newmes (8)=par8;
        newmes (16)=par16;
        SECDED_bit=mod(sum (newmes (1:end-1)),2);
        newmes (n+nbp+1)=SECDED_bit;
end

