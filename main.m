clc; clearvars; close all;

%% -------------------------------------------------------------------------
% 1. Define TX and RX codebooks used for intial beam selection.
% -------------------------------------------------------------------------
% codebooks span azimuth and elevation
az_deg = [-56:8:56].';
el_deg = [-8:8:8].';

% assume TX and RX codebooks to be the same
txcb_az = az_deg;
txcb_el = el_deg;
rxcb_az = az_deg;
rxcb_el = el_deg;

% full TX codebook in az-el
num_tx_az = length(txcb_az);
num_tx_el = length(txcb_el);
txcb_azel = [repmat(txcb_az,num_tx_el,1) repelem(txcb_el,num_tx_az,1)];

% full RX codebook in az-el
num_rx_az = length(rxcb_az);
num_rx_el = length(rxcb_el);
rxcb_azel = [repmat(rxcb_az,num_rx_el,1) repelem(rxcb_el,num_rx_az,1)];

% size of codebooks (number of beams)
num_tx = length(txcb_azel(:,1));
num_rx = length(rxcb_azel(:,1));
disp(['Number of beams in TX codebook: ' num2str(num_tx)]);
disp(['Number of beams in RX codebook: ' num2str(num_rx)]);
disp(['Total number of TX-RX beam combinations: ' num2str(num_tx*num_rx)]);

%% -------------------------------------------------------------------------
% 2. Set size and resolution of spatial neighborhoods (design parameters).
% -------------------------------------------------------------------------
% neighborhood size
Delta_az = 2;
Delta_el = 2;

% neighborhood resolution
delta_az = 1;
delta_el = 1;

% separate neighborhood sizes for TX and RX
Delta_az_tx = Delta_az;
Delta_el_tx = Delta_el;
Delta_az_rx = Delta_az;
Delta_el_rx = Delta_el;

% separate neighborhood resolutions for TX and RX
delta_az_tx = delta_az;
delta_el_tx = delta_el;
delta_az_rx = delta_az;
delta_el_rx = delta_el;

%% -------------------------------------------------------------------------
% 3. Construct nominal neighborhood.
% -------------------------------------------------------------------------
% TX and RX neighborhoods
nbr_tx = construct_neighborhood(Delta_az_tx,Delta_el_tx,delta_az_tx,delta_el_tx);
nbr_rx = construct_neighborhood(Delta_az_rx,Delta_el_rx,delta_az_rx,delta_el_rx);
disp(['Total number of TX candidates in neighborhood: ' num2str(length(nbr_tx(:,1)))]);
disp(['Total number of RX candidates in neighborhood: ' num2str(length(nbr_rx(:,1)))]);

% full TX-RX neighborhood
nbr = [repmat(nbr_tx,length(nbr_rx(:,1)),1) repelem(nbr_rx,length(nbr_tx(:,1)),1)];
num_nbr = length(nbr(:,1));
disp(['Total number of TX-RX candidates in neighborhood: ' num2str(num_nbr)]);

% sort full neighborhood by distance (can modify distance metric as desired)
[~,idx] = sort(sum(nbr.^2,2));
nbr = nbr(idx,:);

%% -------------------------------------------------------------------------
% 4. For each possible transmit-receive initial beam selection, run STEER.
% -------------------------------------------------------------------------
% set INR target (design parameter)
INR_tgt_dB = -7; % lower = stricter threshold

% reset counter
idx_row = 0;

% reset lookup table
lut = [];

% for each TX beam in codebook
for idx_tx = 1:length(txcb_azel)
    % initial TX steering direction
    tx_azel = txcb_azel(idx_tx,:);
    
    % for each RX beam in codebook
    for idx_rx = 1:length(rxcb_azel)
        % initial RX steering direction
        rx_azel = rxcb_azel(idx_rx,:);
        
        % reset min INR
        INR_min_dB = Inf;
        
        % reset counter
        num_meas = 0;
                
        % search over neighborhood
        for idx_nbr = 1:num_nbr
            % shift TX direction
            tx_azel_shift = nbr(idx_nbr,1:2);
            tx_azel_nbr = tx_azel + tx_azel_shift;
            
            % shift RX direction
            rx_azel_shift = nbr(idx_nbr,3:4);
            rx_azel_nbr = rx_azel + rx_azel_shift;
            
            % measure INR (this depends on either a model or measurements)
            INR_meas_dB = 3*randn(1);
            
            % record number of measurements required
            num_meas = num_meas + 1;
            
            % record nominal INR (with initial selection)
            if idx_nbr == 1
                INR_nom_dB = INR_meas_dB;
            end
            
            % check if new minimum INR found
            if INR_meas_dB < INR_min_dB
                % record new best INR and steering directions
                INR_min_dB = INR_meas_dB;
                tx_azel_opt = tx_azel_nbr;
                rx_azel_opt = rx_azel_nbr;
                
                % check if target met
                if INR_meas_dB <= INR_tgt_dB
                    break;
                end
            end
        end
        
        % add to lookup table
        idx_row = idx_row + 1;
        lut(idx_row,:) = [tx_azel, tx_azel_opt, rx_azel, rx_azel_opt, INR_nom_dB, INR_min_dB, num_meas];
    end
end

%% -------------------------------------------------------------------------
% 5. Print lookup table.
% -------------------------------------------------------------------------
% print lookup table:
% each row is a different initial TX-RX direction pair
% - cols 1-2: initial TX direction
% - cols 3-4: shifted TX direction
% - cols 5-6: initial RX direction
% - cols 7-8: shifted RX direction
% - col 9: INR (in dB) with initial TX-RX directions
% - col 10: INR (in dB) with shifted TX-RX directions
% - col 11: number of measurements taken
lut

%% -------------------------------------------------------------------------
% 6. Plot INR with and without STEER.
% -------------------------------------------------------------------------
figure(1);
[f,x] = ecdf(lut(:,9));
plot(x,f,'k--','DisplayName','Without STEER');
hold on;
[f,x] = ecdf(lut(:,10));
plot(x,f,'k-','DisplayName','With STEER');
hold off;
xlabel('Self-Interference, INR (dB)');
ylabel('Cumulative Probability');
grid on;
grid minor;
legend('Location','Best','FontSize',10);

%% -------------------------------------------------------------------------
% 7. Plot how many measurements are needed.
% -------------------------------------------------------------------------
[f,x] = ecdf(lut(:,11));
figure(2);
plot(x,f,'k-');
xlabel('Number of Measurements');
ylabel('Cumulative Probability');
grid on;
grid minor;