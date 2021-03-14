%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spiros Daskalakis 
% 26/12/2017
% Email: daskalakispiros@gmail.com
% Website: www.daskalakispiros.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [packetcor] = HAmmingDecoded18bits(sendbits)
 %function insert_parity_bits calculates and inserts the parity bits

        sendbits_change =sendbits;
        
        %%
        
        res_par1=[sendbits_change(17) sendbits_change(15) sendbits_change(13) sendbits_change(11) sendbits_change(9) sendbits_change(7) sendbits_change(5) sendbits_change(3)];
        res_par1= mod(sum (res_par1),2);
        
        if (res_par1==sendbits_change(1))
             %disp('Correct Par1');
             Par1ind1= 0;
        else 
            %disp('Error in Par1');
            Par1ind1= 1;
        end 
        
        res_par2=[sendbits_change(18) sendbits_change(15) sendbits_change(14) sendbits_change(11) sendbits_change(10) sendbits_change(7) sendbits_change(6) sendbits_change(3)];
        res_par2= mod(sum (res_par2),2);
        
        if (res_par2==sendbits_change(2))
             %disp('Correct Par2');
             Par1ind2= 0;
        else 
            %disp('Error in Par2');
             Par1ind2= 2;
        end 
        
        res_par4=[ sendbits_change(15) sendbits_change(14) sendbits_change(13) sendbits_change(12) sendbits_change(7) sendbits_change(6) sendbits_change(5) ];
        res_par4= mod(sum (res_par4),2);
        
        if (res_par4==sendbits_change(4))
             %disp('Correct Par4');
             Par1ind4= 0;
        else 
            %disp('Error in Par4');
             Par1ind4= 4;
        end 
        
        res_par8=[ sendbits_change(15) sendbits_change(14)  sendbits_change(13) sendbits_change(12) sendbits_change(11) sendbits_change(10) sendbits_change(9)];
        res_par8= mod(sum (res_par8),2);
        
         if (res_par8==sendbits_change(8))
            % disp('Correct Par8');
             Par1ind8= 0;
        else 
           % disp('Error in Par8');
             Par1ind8= 8;
         end 
         
         
        res_par16=[ sendbits_change(18) sendbits_change(17) ];
        res_par16= mod(sum (res_par16),2);
        
         if (res_par16==sendbits_change(16))
            %disp('Correct Par16');
            Par1ind16= 0;
        else 
            %disp('Error in Par16');
           Par1ind16= 16;
         end 
          
         SECDED_bit=mod(sum (sendbits_change (1:end-1)),2);
         c=Par1ind1+Par1ind2+Par1ind4+Par1ind8+Par1ind16;
         
         if(c>length (sendbits_change)) 
             disp('Error');      
         elseif (c~=0 && SECDED_bit~=sendbits_change (end)) 
             disp('One bit Error ');             
             sendbits_change(c)=~sendbits_change(c);
         elseif (c~=0  && SECDED_bit==sendbits_change (end))
            disp('Two bit Error ');
         elseif(c==0 && SECDED_bit==sendbits_change (end)) 
            disp('No Bit  Error');
         elseif(c==0 && SECDED_bit~=sendbits_change (end)) 
            disp('Last Par Bit  Error');
         end 
         packetcor=[sendbits_change(3) sendbits_change(5:7) sendbits_change(9:15) sendbits_change(17:18)];
 
end

