% pulizia workspace
clearvars
clc

%% definizione e aggiunta percorsi

%%%%%%%%% percorso assoluto comune in cui avete le cartelle/sottocartelle con TUTTE
% le funzioni necessarie
% esempio: '/Users/simonetoffoli/funzioniPenna/'
 function_path = '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/utilities';
% function_path = '/Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/TesiTriennaleAnziani/utilities';

% aggiunta percorso con funzioni
addpath(genpath(function_path)) %aggiunge cartella functions al path

%%%%%%%%% percorso assoluto comune delle cartelle/sottocartelle in cui avete tutti i dati
% esempi:
% anziani: '/Users/simonetoffoli/datiAnziani/' -> nella cartella
% 'datiAnziani' ci sono le sotto cartelle 'Spirale', 'FraseCorsivo', 'FraseStampato'

% bambini: '/Users/simonetoffoli/datiBambini/' -> nella cartella
% 'datiBambini' ci sono le sotto cartelle '1_T1, '1_T2' e '1_T3',
% contententi a loro volta 'BVSCO', 'Frase', 'TPV'.
% Ognuna di queste contiene altre sotto cartelle, una per ogni prova del
% protocollo. Per esempio, in 'BVSCO' ci sono 'LE', 'NUMERI_stampatello',
% 'NUMERI_corsivo', 'UNO_stampatello', 'UNO_corsivo'

common_data_path = '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/Dati/';
%common_data_path = '//Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/TesiTriennaleAnziani/Dati/';

%%%%%%%%% percorso assoluto specifico della cartella con i dati che volete
%%%%%%%%% analizzare

% anziani:
task = 'SentenceCursive'; %% una tra 'Spirale', 'FraseCorsivo', 'FraseStampato'
data_folder = [common_data_path  task];


%% creazione cartella in cui salvare gli indicatori

save_folder = fullfile('/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/Indicatori');
% esempio: save_folder = fullfile('/Users/simonetoffoli/datiPenna/Indicatori/');

% if ~isfolder(save_folder) %crea cartella se non esiste già
%     mkdir(save_folder)
% end

%% estrazione elenco dei dati presenti nella cartella specificata (data_folder)

subjects_list= dir(fullfile(data_folder, '*.mat')); %cerco file mat nella cartella data_folder
n_subjects = length(subjects_list); %ottengo numero file mat nella cartella data_folder
subjects = cell(n_subjects, 1); %crea una cella per ogni soggetto nella cartella data_folder
for i = 1:n_subjects
    subjects{i} = strrep(subjects_list(i).name, '.mat', ''); %in ogni cella mantiene solo nome file senza estensione .mat

    %%%%%%%%% INIZIO RIGHE SOLO PER ANZIANI %%%%%%%%%

    subjects_code{i,1} = extractBefore(subjects{i}, '_2023');
    subjects_code{i,1} = extractAfter(subjects_code{i,1}, '_');

    %%%%%%%%% FINE RIGHE SOLO PER ANZIANI %%%%%%%%%
end
%%
unique_codes = unique(subjects_code);
%% flag LeLe
% anziani:
% = 0 se analizzate frasi; = 1 se analizzate spirale

isLeLe = 0;

%% SOLO PER ANZIANI
% percorso assoluto alla cartella che contiene il file
% 'MMSE_usuarios_ESSENCE.xlsx', contenente i punteggi MMSE dei soggetti
% all'inizio e alla fine dello studio
% N.B.: ve l'ho caricato nella cartella condivisa 'TesiTriennaliAnziani'
MMSE_path = ['Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/'];
%MMSE_path = ['Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/TesiTriennaleAnziani/'];

% esempio: '/Users/simonetoffoli/MMSE/'

demo = readtable([MMSE_path 'MMSE_usuarios_ESSENCE.xlsx'], 'Sheet', 'Foglio1'); %leggo file xlsx
demo = table2struct(demo); % converto tabella in struttura
k=1;
%% calcolo indicatori per ogni dato a disposizione


for i = 292:length(subjects) %per ogni file mat

    clear data
    % carico il dato penna (il file mat)
    load(fullfile(data_folder, strcat(subjects{i}, '.mat')));

    %%%%%%%%% INIZIO RIGHE SOLO PER ANZIANI %%%%%%%%%
    clear ind_s
    ind_s = find(contains({demo.User}.',subjects_code{i}));

    if isLeLe ==1
        percentage= 0.15;
        beads_flag = 0; r = 0; fc =0; amp = 0; lam0 = 0 ;cost = 0 ;
    else
        percentage= 0;
        beads_flag = 0; r = 6; fc = 0.006; amp = 0.8; lam0 = 0.5*amp, cost=0.10;
    end

    %%%%%%%%% FINE RIGHE SOLO PER ANZIANI %%%%%%%%%


    % calcolo indicatori per il singolo dato
    Indicators = FunctionIndicatorsComputation(data, 500, 12, percentage, fc, r, amp, lam0, beads_flag, cost, 'NEW',isLeLe);

    % aggiungo nome del soggetto alla lista di indicatori. Questo potrebbe
    % essere da modificare in base al gruppo (anziani o bambini)
    Indicators.subject = subjects_code{i};

    %%%%%%%%% INIZIO RIGHE SOLO PER ANZIANI %%%%%%%%%
    Indicators.age = demo(ind_s).AGE;
    Indicators.MMSE = demo(ind_s).MMSE;
    Indicators.MMSE_age_adj = demo(ind_s).MMSE_AGE_ADJUSTED;
    Indicators.MMSE_post = demo(ind_s).MMSE_post;
    Indicators.MMSE_post_age_adj = demo(ind_s).MMSE_age_adj_post;
    Indicators.years_education = demo(ind_s).YEARS_EDU;
    %%%%%%%%% FINE RIGHE SOLO PER ANZIANI %%%%%%%%%

    % aggiorno struttura che contiene gli indicatori di tutti i soggetti
    S(k,1) = Indicators;
    k=k+1;

    % pulisco lista indicatori del singolo dato e il singolo dato
    clear Indicators data
end


%% salvataggio indicatori come file xlsx

% converto struttura in tabella
T_all=struct2table(S);

% assegno nome alla tabella che voglio salvare
filename = 'Tabella';
% aggiungo estensione al nome della tabella
% N.B.: qui potete includere nel nome del file alcune delle variabili che
% avete assegnato manualmente, come 'task' o 'grade', in modo da salvare
% programmaticamente la tabella
% esempio: save_name_table = [filename '_' task '.xlsx']
save_name_table = [filename 'SentenceCursive.xlsx'];

% salvo file tabella come file xlsx
if exist([save_folder '/' save_name_table], 'file')
    % se file esiste già, modifico il file in modalità append
    writetable(T_all, fullfile([save_folder '/' save_name_table]),"WriteVariableNames",false, 'WriteMode','append');
else
    % se il file non esiste, lo creo per la prima volta
    writetable(T_all, fullfile([save_folder '/' save_name_table]),"WriteVariableNames",true);
end