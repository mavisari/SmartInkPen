%% pulizia workspace
clearvars
clc

%% task da analizzare
task= 'SentenceCapitals';

%% definizione e aggiunta percorsi
utente = 'simone';
utente = upper(utente);

if contains(utente, 'SIMONE')
    root= '/Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/';
    common_data_path = [root '/TesiTriennaleAnziani/IndicatoriSimone/'];
    MDC_path = [root '/TesiTriennaleAnziani/TabelleMDCSimone/'];
    save_folder = [root 'TesiTriennaleAnziani/DirezionePatologiaSimone'];
elseif contains(utente,'ALBERTA')
    root = '\Users\alber\OneDrive - Politecnico di Milano\TesiTriennaleAnziani\';
    common_data_path = [root '\Indicatori\'];
    addpath '\Users\alber\OneDrive - Politecnico di Milano\TesiTriennaleAnziani\utilities';
elseif contains(utente,'MARIA VITTORIA')
    root = '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/';
    common_data_path = [root '/TesiTriennaleAnziani/Indicatori/'];
    addpath '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/utilities'; 
    
elseif contains(utente,'MICOL')
    root= '/Users/micol/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Simone Toffoli/TesiTriennali/';
    common_data_path = [root '/TesiTriennaleAnziani/Indicatori/'];
    addpath '/Users/micol/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Simone Toffoli/TesiTriennali/TesiTriennaleAnziani/utilities';
end



addpath(genpath(common_data_path));
addpath(MDC_path)

%% estrazione elenco dei dati presenti nella cartella specificata (data_folder)
file_list= dir(fullfile(common_data_path, '*.mat')); %cerco file excel nella cartella data_folder
n_file = length(file_list); %ottengo numero file mat nella cartella data_folder
file = cell(n_file, 1); %crea una cella per ogni soggetto nella cartella data_folder
for i = 1:n_file
    file{i} = strrep(file_list(i).name, '.mat', ''); %in ogni cella mantiene solo nome file senza estensione .mat
end


indice = find(contains(file, task));
% fileexcel = readtable([file{indice}]); fileexcel = table2struct(fileexcel);
load([file{indice}]); fileexcel = S;
%% file MDC
file_list_MDC= dir(fullfile(MDC_path, '*.xlsx')); %cerco file excel nella cartella data_folder
n_file_MDC = length(file_list_MDC); %ottengo numero file mat nella cartella data_folder
file_MDC = cell(n_file_MDC, 1); %crea una cella per ogni soggetto nella cartella data_folder
for i = 1:n_file_MDC
    file_MDC{i} = strrep(file_list_MDC(i).name, '.xlsx', ''); %in ogni cella mantiene solo nome file senza estensione .mat
end
indice_MDC = find(contains(file_MDC, task));
fileexcelMDC=readtable([file_MDC{indice} '.xlsx']);
fileexcelMDC=table2struct(fileexcelMDC); 
%% file MMSE
fileexcelMMSE=readtable('MMSE_usuarios_ESSENCE.xlsx'); 
fileexcelMMSE=table2struct(fileexcelMMSE); 

%% rimuovo indicatori non affidabili da fileexcel
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

for i = 1:n_soggetti
    codice_soggetto{i,1} = ['soggetto_' unique_subjects{i}];
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


%% estrazione prove e differenze 
fields = fieldnames(struttura_soggetto);
array_indicatori={}; 
array_soggetto={}; 
array_differenza=[]; 
array_MDC=[]; 
array_trend=[]; 
n=1; 

for i=1:n_soggetti
    
     
    %scorro tutti i soggetti
    soggetto=fields{i}; %i
    tabella_soggetto=struttura_soggetto.(soggetto); 
    

    differenzaabs = zeros(numero_prove_soggetto(i) - 2, numero_indicatori - 1);
    %scorro tutti gli indicatori e per ognuno faccio prova k - prova 1 
    for j=1:(numero_indicatori-1)
        indicators=indicatori{j};
        valore_prova1=tabella_soggetto(1).(indicators);
        for k=3:numero_prove_soggetto(i)
            valore=tabella_soggetto(k).(indicators);
            differenza(k-2,j)=valore-valore_prova1;
            differenzaabs(k-2,j)=abs(valore-valore_prova1);

            if differenzaabs(k-2,j) >= fileexcelMDC(j).MDC

                if differenza(k-2,j)>0 && fileexcelMDC(j).trend==1
                    array_differenza(n,1)=differenzaabs(k-2,j);
                    array_trend(n,1)=fileexcelMDC(j).trend;
                    array_soggetto{n,1}=soggetto;
                    array_prova(n,1) = k;
                    array_MDC(n,1)=fileexcelMDC(j).MDC;
                    array_indicatori{n,1}=indicatori{j};
                    array_MMSE(n,1)=fileexcelMMSE(i).MMSE;
                    array_MMSE_post(n,1)=fileexcelMMSE(i).MMSE_post; 
                    array_dMMSE(n,1)=fileexcelMMSE(i).MMSE_post-fileexcelMMSE(i).MMSE; 
                    array_MMSE_adj (n,1)=fileexcelMMSE(i).MMSE_AGE_ADJUSTED; 
                    array_MMSE_adj_post(n,1)=fileexcelMMSE(i).MMSE_age_adj_post; 
                    array_MMSE_dadj(n,1)=fileexcelMMSE(i).MMSE_age_adj_post-fileexcelMMSE(i).MMSE_AGE_ADJUSTED; 
                    n=n+1;
                elseif differenza(k-2,j)<0 && fileexcelMDC(j).trend==0
                    array_differenza(n,1)=differenzaabs(k-2,j);
                    array_trend(n,1)=fileexcelMDC(j).trend;
                    array_soggetto{n,1}=soggetto;
                    array_prova(n,1) = k;
                    array_MDC(n,1)=fileexcelMDC(j).MDC;
                    array_indicatori{n,1}=indicatori{j};
                    array_MMSE(n,1)=fileexcelMMSE(i).MMSE;
                    array_MMSE_post(n,1)=fileexcelMMSE(i).MMSE_post; 
                    array_dMMSE(n,1)=fileexcelMMSE(i).MMSE_post-fileexcelMMSE(i).MMSE; 
                    array_MMSE_adj (n,1)=fileexcelMMSE(i).MMSE_AGE_ADJUSTED; 
                    array_MMSE_adj_post(n,1)=fileexcelMMSE(i).MMSE_age_adj_post; 
                    array_MMSE_dadj(n,1)=fileexcelMMSE(i).MMSE_age_adj_post-fileexcelMMSE(i).MMSE_AGE_ADJUSTED; 
                    n=n+1;
                end
            end
        end
    end
   
end 

array_differenza=num2cell(array_differenza); 
array_MDC=num2cell(array_MDC); 
array_trend=num2cell(array_trend); 
array_prova=num2cell(array_prova); 
array_MMSE=num2cell(array_MMSE); 
array_MMSE_post=num2cell(array_MMSE_post); 
array_dMMSE=num2cell(array_dMMSE); 
array_MMSE_adj=num2cell(array_MMSE_adj); 
array_MMSE_adj_post=num2cell(array_MMSE_adj_post); 
array_MMSE_dadj=num2cell(array_MMSE_dadj); 

direzione_patologia=[array_soggetto array_prova array_indicatori array_MDC array_trend array_differenza array_MMSE array_MMSE_post array_dMMSE array_MMSE_adj array_MMSE_adj_post array_MMSE_dadj]; 



%% File excel
if contains(utente, 'MARIA VITTORIA')
        save_folder = fullfile('/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/DirezionePatologia/');
end
data_table=cell2table(direzione_patologia, "VariableNames", ["Soggetti" "Prova" "Indicatori" "MDC" "trend" "differenza" "MMSE" "MMSE_post" "ΔMMSE" "MMSE_age-adj" "MMSE_age_adj_post" "ΔMMSE_age_adj"]);
if ~isfolder(save_folder)
    mkdir(save_folder)
end
filename='DirezionePatologia'; 
save_name_table = [filename '_' task '.xlsx'];
writetable(data_table, fullfile([save_folder '/' save_name_table])); 
