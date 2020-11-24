%% clean

clc; clear all; close all;

M = 2; % antennes receptions
N = 2; % antennes emissions

L = 2; % longueur du mot de code
M_psk = 4;
nb_mdc = 100;
K = L*nb_mdc*M_psk;
phi0 = 0;
var_signal = 1;
%varv = 1; % variance bruit

EbN0dB_min  = -2; % Minimum de EbN0
EbN0dB_max  = 10; % Maximum de EbN0
EbN0dB_step = 1;% Pas de EbN0

nbr_erreur_cap = 200;  % Nombre d'erreurs a observer avant de calculer un BER
nbr_bit_max = 100e6;% Nombre de bits max a simuler
ber_min     = 1e-6; % BER min

EbN0dB = EbN0dB_min:EbN0dB_step:EbN0dB_max;     % Points de EbN0 en dB a simuler
EbN0   = 10.^(EbN0dB/10);% Points de EbN0 a simuler
sigma2 = 1./EbN0;

%% Construction du modulateur
mod_psk = comm.PSKModulator(...
    'ModulationOrder', M_psk, ... % BPSK
    'PhaseOffset'    , phi0, ...
    'SymbolMapping'  , 'Gray',...
    'BitInput'       , true);

%% Construction du demodulateur
demod_psk = comm.PSKDemodulator(...
    'ModulationOrder', M_psk, ...
    'PhaseOffset'    , phi0, ...
    'SymbolMapping'  , 'Gray',...
    'BitOutput'      , true,...
    'DecisionMethod' , 'Log-likelihood ratio');

%% préparation ML
C = zeros(2,2,4^4);
c_seq = 0:255;
c_bin = dec2bin(c_seq);
for c=1:256

    one_code = c_bin(c,:);

    one_code_matrix = reshape(one_code, 2, [])-'0';
    s = zeros(1,4);
    for j=1:length(s)
        s(j) = mod_psk(one_code_matrix(:,j));

    end
    s = reshape(s,2,2);
    C(:,:,c) = s;
end



% vn = sqrt(varv/2)*(randn(M,L) + 1i*randn(M,L)); % autant de bruit thermiques de capteurs recepteurs

% H = sqrt(varv/2)*(randn(M,N) + 1i*randn(M,N)); % h une 2X2 : les coefs des différents "sous-canaux" créés avec un coef par sous-canal
% on les regénère a chaque fois donc on va le foutre dans la boucle

nb_errs = zeros(size(EbN0));

ber = zeros(size(EbN0));
nb_err = 0;
for i_snr = 1:length(EbN0)
    
    i_snr
    varv = sigma2(i_snr);
    nb_bits_sent = 0;
    nb_err = 0;
    while (nb_err < nbr_erreur_cap)

        %% chaine de com
%         for i=1:nb_mdc


            %% emission
%             b = randi([0,1],K,1);    % Generation du message al?atoire
% 
%             x_stream = mod_psk(b);
% 
%             %% attention, voilà le gros encodage spatiotemporel de la mort qui tue V_BLAST
%             meta_X = reshape(x_stream, N, []); % la première 2x2 = le premier mdc=X, il y a nb_mdc matrices 2x2
%             % Voila
% 
%             X = meta_X(1:2, i:i+L-1);
            
            rand_imdc = randi([1,256], 1,1);
            X = C(:,:,rand_imdc);

            %% canal
            H = sqrt(var_signal/2)*(randn(M,N) + 1i*randn(M,N));
            V = sqrt(varv/2)*(randn(M,L) + 1i*randn(M,L)); % autant de bruit thermiques de capteurs recepteurs

            Y = H*X + V;

            %% Reception
            i_min = 1;
            min = norm(Y - H*C(:,:,1), 'fro');
            
            for i_codepossible = 2:4^4 % optimisable

                X_tmp = C(:,:,i_codepossible);
                mynorm = norm(Y - H*X_tmp, 'fro');
                if mynorm < min
                    min = mynorm;
                    i_min = i_codepossible;
                end

            end
            
        X_chap = C(:,:,i_min);
        err = sum(sum(X_chap~=X))>0;
        nb_err = nb_err + err;
        nb_bits_sent = nb_bits_sent + 1;
        nb_errs(i_snr) = nb_err;
        ber(i_snr)= nb_errs(i_snr)/nb_bits_sent;
    end

end

semilogy(EbN0dB, ber);