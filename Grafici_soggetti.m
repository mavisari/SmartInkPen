%% pulizia workspace
clearvars
clc

%% task da analizzare
task= 'SentenceCapitals';

%% definizione e aggiunta percorsi
utente = 'maria vittoria';
utente = upper(utente);

if contains(utente, 'SIMONE')
    root= '/Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/';
    common_data_path = [root '/TesiTriennaleAnziani/Indicatori/'];
    save_folder = [root 'TesiTriennaleAnziani/DirezionePatologia'];
elseif contains(utente,'ALBERTA')
    root = '\Users\alber\OneDrive - Politecnico di Milano\TesiTriennaleAnziani\';
    common_data_path = [root '\Indicatori\'];
    addpath '\Users\alber\OneDrive - Politecnico di Milano\TesiTriennaleAnziani\utilities';
elseif contains(utente,'MARIA VITTORIA')
    root = '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/';
    common_data_path = [root '/TesiTriennaleAnziani/IndicatoriSimone/'];
    addpath '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/utilities'; 
    plot_folder = '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/GraficiSoggetti2/GraficiSentenceCapital';
elseif contains(utente,'MICOL')
    root= '/Users/micol/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Simone Toffoli/TesiTriennali/';
    common_data_path = [root '/TesiTriennaleAnziani/Indicatori/'];
    addpath '/Users/micol/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Simone Toffoli/TesiTriennali/TesiTriennaleAnziani/utilities';
end


addpath(genpath(common_data_path));

%% estrazione elenco dei dati presenti nella cartella specificata (data_folder)
file_list= dir(fullfile(common_data_path, '*.xlsx')); %cerco file excel nella cartella data_folder
n_file = length(file_list); %ottengo numero file mat nella cartella data_folder
file = cell(n_file, 1); %crea una cella per ogni soggetto nella cartella data_folder
for i = 1:n_file
    file{i} = strrep(file_list(i).name, '.xlsx', ''); %in ogni cella mantiene solo nome file senza estensione .mat
end

%% leggere il file excel
% indice = find(contains(file, task));
% fileexcel = readtable([file{indice}]); fileexcel = table2struct(fileexcel);
load('/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/IndicatoriSimone/TabellaSentenceCapitals.mat')
fileexcel = S; 
fileexcelMDC=readtable('TabelleMDCSimone/MDCSentenceCapitals.xlsx');
fileexcelMDC=table2struct(fileexcelMDC); 
fileDirezionePatologia=readtable('DirezionePatologiaSimone/DirezionePatologia_SentenceCapitals.xlsx'); 
fileDirezionePatologia=table2struct(fileDirezionePatologia); 

%% rendo fileexcel con indicatori fileMDC
nomi_campi={fileexcelMDC.Indicatori}; 
nomi_campi=vertcat(nomi_campi', {'subject'}); 

for i=1:numel(fileexcel)
    fileexcel=rmfield(fileexcel, setdiff(fieldnames(fileexcel(i)), nomi_campi')); 
end 

%% estrazione N prove per soggetto

[unique_subjects, IA, IC] = unique({fileexcel.subject}.');
numero_prove_soggetto = groupcounts({fileexcel.subject}.');
numero_righe=size(fileexcel,1); 
indicatori=fieldnames(fileexcel); 
numero_indicatori=size(indicatori,1); 
n_soggetti=size(unique_subjects,1); 

codice_soggetto = cell(n_soggetti, 1); 

for i = 1:n_soggetti
%     codice_soggetto{i,1} = ['soggetto_' num2str(i)];
    codice_soggetto{i, 1} = ['soggetto_' num2str(unique_subjects{i})]; 

end


%%
clear struttura_soggetto
for i = 1:n_soggetti
    
    clear ind_start ind_stop k
    k=1;
    ind_start = IA(i);
     
    if i==n_soggetti
       ind_stop = length(fileexcel);
    else
       ind_stop=IA(i+1)-1;
    end


    for j=ind_start:ind_stop
        struttura_soggetto.(codice_soggetto{i})(k,1) = fileexcel(j);
        k=k+1;
    end

end


%% Grafici per ogni soggetto 

if strcmp(task, 'Spiral')
    indicatori_da_estrarre = {'ExecutionTime', 'S3_median_s', 'RelPOWPC_g_PC_WELCH', 'PeakF_G_PC', 'meanOnsheet', 'meanP', 'ConsPeakDiffPMean'};
elseif strcmp(task, 'SentenceCursive')
    indicatori_da_estrarre = {'S3_median_s', 'PauseNum_Rel', 'ConsPeakDiffGMean', 'OnSheetRatio', 'P_OVS', 'ConsPeakDiffPMean'};
elseif strcmp(task, 'SentenceCapitals')
    indicatori_da_estrarre = {'OnSheetRatio', 'meanOnsheet', 'ConsPeakDiffPMean', 'LDLJ_G_median', 'meanP', 'S2_median_s'};
end

soggetti = fieldnames(struttura_soggetto);

for i = 1 : n_soggetti

    figure1 = figure;
    figure2 = figure;

    soggetto = soggetti{i};
    tabella_soggetto = struttura_soggetto.(soggetto);
    grafico_num_gruppo1 = 1;
    grafico_num_gruppo2 = 1;

    indicatori_gruppo1 = indicatori_da_estrarre(1:4);
    indicatori_gruppo2 = indicatori_da_estrarre(5:end);

    
    for j = 1:length(indicatori_gruppo1)
        indicatore = indicatori_gruppo1{j};
        dati = zeros(numero_prove_soggetto(i), 1); 
        for k = 1:numero_prove_soggetto(i)
            dati(k, 1) = tabella_soggetto(k).(indicatore);
        end
        
        prova1 = dati(1, 1);
        mdc = 0; 
        for z = 1:size(fileexcelMDC, 1)
            if strcmp(indicatore, fileexcelMDC(z).Indicatori)
                if fileexcelMDC(z).trend == 1
                    mdc = fileexcelMDC(z).MDC; 
                elseif fileexcelMDC(z).trend == 0
                    mdc = -fileexcelMDC(z).MDC; 
                end
            end 
        end 

        prova1_mdc = dati(1, 1) + mdc; 
        figure(figure1); % Passa alla prima figura
        subplot(length(indicatori_gruppo1), 1, grafico_num_gruppo1);
        
        plot(dati);
        hold on;
        yline(prova1, 'k-', 'LineWidth', 2);
        yline(prova1_mdc, 'r-', 'LineWidth', 2);

        title(indicatore); 
        grafico_num_gruppo1 = grafico_num_gruppo1 + 1;
    end

    
    for j = 1:length(indicatori_gruppo2)
        indicatore = indicatori_gruppo2{j};
        dati = zeros(numero_prove_soggetto(i), 1); 
        for k = 1:numero_prove_soggetto(i)
            dati(k, 1) = tabella_soggetto(k).(indicatore);
        end
        
        prova1 = dati(1, 1);
        mdc = 0; 
        for z = 1:size(fileexcelMDC, 1)
            if strcmp(indicatore, fileexcelMDC(z).Indicatori)
                if fileexcelMDC(z).trend == 1
                    mdc = fileexcelMDC(z).MDC; 
                elseif fileexcelMDC(z).trend == 0
                    mdc = -fileexcelMDC(z).MDC; 
                end
            end 
        end 

        prova1_mdc = dati(1, 1) + mdc; 
        figure(figure2); 
        subplot(length(indicatori_gruppo2), 1, grafico_num_gruppo2);
        
        plot(dati);
        hold on;
        yline(prova1, 'k-', 'LineWidth', 2);
        yline(prova1_mdc, 'r-', 'LineWidth', 2);

        title(indicatore); 

        grafico_num_gruppo2 = grafico_num_gruppo2 + 1;
    end

  
    nome_file1 = sprintf('%s_%s_gruppo1.fig', soggetto, ''); 
    nome_file2 = sprintf('%s_%s_gruppo2.fig', soggetto, '');
    nome_file1_completo = fullfile(plot_folder, nome_file1);
    nome_file2_completo = fullfile(plot_folder, nome_file2);
    saveas(figure1, nome_file1_completo);
    saveas(figure2, nome_file2_completo);
    close(figure1);
    close(figure2);

end
