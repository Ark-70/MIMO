%% clean

clc; clear all; close all;

M = 2; % antennes receptions
N = 2; % antennes emissions

L = 2; % longueur du mot de code
M_psk = 4;
nb_mdc = 100;
K = L*nb_mdc*M_psk;
phi0 = 0;
varv = 1; % variance bruit

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



vn = sqrt(varv/2)*(randn(M,L) + 1i*randn(M,L)); % autant de bruit thermiques de capteurs recepteurs

H = sqrt(varv/2)*(randn(M,N) + 1i*randn(M,N)); % h une 2X2 : les coefs des différents "sous-canaux" créés avec un coef par sous-canal
% on les regénère a chaque fois donc on va le foutre dans la boucle

b = randi([0,1],K,1);    % Generation du message al�atoire

x_stream = mod_psk(b);

%% attention, voilà le gros encodage spatiotemporel de la mort qui tue V_BLAST
meta_X = reshape(x_stream, N, []); % la première 2x2 = le premier mdc=X, il y a nb_mdc matrices 2x2
%% Voila

% for i=1:nb_mdc
%     X = 


% end

