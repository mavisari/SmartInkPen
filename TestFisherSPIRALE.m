%% pulizia workspace
clearvars
clc

%% task da analizzare
task= 'Spiral';

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
    common_data_path = [root '/TesiTriennaleAnziani/Indicatori/'];
    addpath '/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/utilities'; 
    
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
indice = find(contains(file, task));
fileexcel = readtable([file{indice}]); fileexcel = table2struct(fileexcel);
fileexcelMDC=readtable('TabelleMDC/MDCSpiral2.xlsx');
fileexcelMDC=table2struct(fileexcelMDC); 
fileDirezionePatologia=readtable('DirezionePatologia/DirezionePatologia3_Spiral.xlsx'); 
fileDirezionePatologia=table2struct(fileDirezionePatologia); 
fileexcelMMSE=readtable('MMSE_usuarios_ESSENCE.xlsx'); 
fileexcelMMSE=table2struct(fileexcelMMSE); 

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

for i = 1:n_soggetti
    codice_soggetto{i,1} = ['soggetto_' num2str(i)];
end

%% Estraggo deltaMMSE per i soggetti  

for i=1:n_soggetti 
    MMSE(i,1)=fileexcelMMSE(i).MMSE;
    MMSE_post(i,1)=fileexcelMMSE(i).MMSE_post; 
    dMMSE(i,1)=fileexcelMMSE(i).MMSE_post-fileexcelMMSE(i).MMSE;
end 

indici_NaN=isnan(dMMSE);  
dMMSE_noNaN=dMMSE(~indici_NaN); 

%% Indicatori in frequenza SPIRALE
indicatore_freq={'RelPOWVOL_a3D_WELCH';'RelPOWDYS_a3D_WELCH';'RelPOWPDt_a3D_WELCH';'RelPOWPHYt_a3D_WELCH'; 'RelPOWPC_g_PC_WELCH'; 'RelPOWPDt_g_PC_WELCH'; 'PeakF_G_PC'; 'Peak_A3D_WELCH'; 'Peak_G_PC'; 'ampXout_A3D'; 'out_lev_A3D'; 'dfA3D_68'; 'dfA3D_50'; 'fA3D_sort_50'; 'fA3D_sort_CM';'ampXout_G_PC_norm'; 'out_lev_G_PC_norm'; 'dfG_PC_68'; 'fG_PC_sort_68'; 'fG_PC_sort_50'; 'fG_PC_sort_CM'; 'bw_A3D' ; 'median_A3D'; 'bw_G_PC'; 'median_G_PC'}; 
totale_freq=size(indicatore_freq,1);

totale_ind_freq = [];

for i=1:n_soggetti 
    freq_prove=[]; 
    freq_ind={}; 
    z=1; 

    for j=1:size(fileDirezionePatologia,1)
        if strcmp(fileDirezionePatologia(j).Soggetti, codice_soggetto(i))
            freq_prove(z,1)=fileDirezionePatologia(j).Prova; 
            freq_ind{z,1}=fileDirezionePatologia(j).Indicatori; 
            z=z+1; 
        end 
    end 
    
    freq_prove_uniche=unique(freq_prove); 
   
    prova_ind = zeros(numel(freq_prove_uniche), 1);
    tot_ind_freq = zeros(numel(freq_prove_uniche), 1);

    for k=1:numel(freq_prove_uniche)
        freq=0; 
        for z=1:size(freq_prove,1)
            if freq_prove(z)==freq_prove_uniche(k)
                indicatore_corrente=freq_ind{z,1}; 
                

                for p=1:totale_freq
                    if strcmp(indicatore_corrente, indicatore_freq{p,1})
                        freq=freq+1; 
                    end 
                end 
            
                prova_ind(k)=freq_prove_uniche(k); 
                tot_ind_freq(k)=freq; 
                 
            end 
        end 
    end 
   
    totale_ind_freq = [totale_ind_freq; tot_ind_freq];

end 

media=mean(totale_ind_freq); 

 %% Calcolo Prove declinate
 
for i=1:n_soggetti 
    freq_prove=[]; 
    freq_ind={}; 
    z=1; 

    for j=1:size(fileDirezionePatologia,1)
        if strcmp(fileDirezionePatologia(j).Soggetti, codice_soggetto(i))
            freq_prove(z,1)=fileDirezionePatologia(j).Prova; 
            freq_ind{z,1}=fileDirezionePatologia(j).Indicatori; 
            z=z+1; 
        end 
    end 
    
    freq_prove_uniche=unique(freq_prove); 
   
    prova_ind = zeros(numel(freq_prove_uniche), 1);
    tot_ind_freq = zeros(numel(freq_prove_uniche), 1);

    for k=1:numel(freq_prove_uniche)
        freq=0; 
        for z=1:size(freq_prove,1)
            if freq_prove(z)==freq_prove_uniche(k)
                indicatore_corrente=freq_ind{z,1}; 
                

                for p=1:totale_freq
                    if strcmp(indicatore_corrente, indicatore_freq{p,1})
                        freq=freq+1; 
                    end 
                end 
            
                prova_ind(k)=freq_prove_uniche(k); 
                tot_ind_freq(k)=freq;
  
            end 
        end 
    end 

    prova_declinata = [];
    d=1; 
    for s=1:numel(freq_prove_uniche)
        if tot_ind_freq(s)>media 
            prova_declinata(d,1)=prova_ind(s);  
            d=d+1;
        end 
    end 
    
    prove_uniche=unique(prova_declinata); 
    n_prove_declinate_sogg(i,1)=numel(prove_uniche); 
    
    prove_no12=3:numero_prove_soggetto(i);
    n_prove_sogg(i,1)=numel(prove_no12); 

end 

%% Calcolo percentuali prove declinate
n_prove_sogg=n_prove_sogg(~indici_NaN); 
n_prove_declinate_sogg=n_prove_declinate_sogg(~indici_NaN); 
codice_soggetto(indici_NaN) = [] 

k=1; 
for i=1:size(n_prove_sogg)
    if n_prove_sogg(i)>0
        perc_soggetto(k,1)=n_prove_declinate_sogg(i)/n_prove_sogg(i);
        dMMSE_tot(k,1)=dMMSE_noNaN(i); 
        n_sing(k,1)=n_prove_declinate_sogg(i); 
        n_prove_nd(k,1)=n_prove_sogg(i)-n_prove_declinate_sogg(i); 
        soggetto_decl{k,1}=codice_soggetto(i); 
        k=k+1; 
    end 
end 

%% Test di Fisher sui soggetti

k=[0.10 0.25 0.50 0.75]; 

for i=1:4
    condizione_prove=perc_soggetto>=k(i); 
    condizione_MMSE=dMMSE_tot<0; 

   %  tabella = table(condizione_prove,  condizione_MMSE);
   % figure; heatmap(tabella, "condizione_prove", "condizione_MMSE");

    [tabella_contingenza, chi2, p, labels]=crosstab(condizione_prove, condizione_MMSE);
   % % cm= confusionmat(condizione_prove, condizione_MMSE);
   % %  confusionchart(cm)


    [h,p,stats] = fishertest(tabella_contingenza); 
    pvalue(i)=p; 
    stat(i)=stats; 

nomi_righe = { 'Prove declinate < %', 'Prove declinate>= %'};
nomi_colonne = {'dMMSE >= 0','dMMSE < 0'};
tabella_con_titoli{i} = array2table(tabella_contingenza, 'RowNames', nomi_righe, 'VariableNames', nomi_colonne);

end 



%% Test per prove 
prove_decl_neg=0;  
prove_nondec_neg=0; 
prove_nondec_pos=0; 
prove_decl_pos=0; 

for i=1:size(dMMSE_tot);  
    if dMMSE_tot(i)<0 
        prove_nondec_neg=prove_nondec_neg+n_prove_nd(i); 
        prove_decl_neg=prove_decl_neg+n_prove_declinate_sogg(i); 
        
    elseif dMMSE_tot(i)>=0
        prove_nondec_pos=prove_nondec_pos+n_prove_nd(i); 
        prove_decl_pos=prove_decl_pos+n_prove_declinate_sogg(i); 
    
    end 
end 

% tabella_cont=[prove_decl_pos prove_decl_neg ; prove_nondec_pos prove_nondec_neg]; 
% nomi_righe2 = { 'Prove declinate', 'Prove non declinate'};
tabella_cont=[ prove_nondec_pos prove_nondec_neg; prove_decl_pos prove_decl_neg]; 
nomi_righe2 = {'Prove non declinate', 'Prove declinate'};
nomi_colonne2 = {'dMMSE >= 0', 'dMMSE < 0'};
tabella_con_titoli2 = array2table(tabella_cont, 'RowNames', nomi_righe2, 'VariableNames', nomi_colonne2);

[hp,pp,statsp] = fishertest(tabella_cont);
