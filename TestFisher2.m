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
    root= '';
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


%% conto prove declinate rispetto alle prove del soggetto CORSIVO STAMPATELLO
n_prove_sing_soggetto=cell(n_soggetti,1);  
perc_prove_soggetto=zeros(n_soggetti,1); 
 

for i=1:n_soggetti
    prove_singole=[]; 
     

    for j=1:size(fileDirezionePatologia,1)

        if strcmp(fileDirezionePatologia(j).Soggetti, codice_soggetto(i))
            prove_singole = [prove_singole, fileDirezionePatologia(j).Prova];
         
        end 
            
        prove_singole_uniche = unique(prove_singole);
        n_prove_sing_soggetto{i} = numel(prove_singole_uniche);
        perc_prove_soggetto(i)=n_prove_sing_soggetto{i}/numero_prove_soggetto(i); 


end 

    prove_no12=3:numero_prove_soggetto(i);
    n_prove_sogg(i,1)=numel(prove_no12); 
    
end  


n_prove_sogg=n_prove_sogg(~indici_NaN); 
n_prove_sing_soggetto=n_prove_sing_soggetto(~indici_NaN); 
perc_p_sogg_noNaN=perc_prove_soggetto(~indici_NaN);
codice_soggetto(indici_NaN) = [] 

k=1; 
for i=1:size(n_prove_sogg)
    if n_prove_sogg(i)>0
        perc_soggetto(k,1)=n_prove_sing_soggetto{i}/n_prove_sogg(i);
        dMMSE_tot(k,1)=dMMSE_noNaN(i); 
        n_sing(k,1)=n_prove_sing_soggetto{i}; 
        n_prove_nd(k,1)=n_prove_sogg(i)-n_prove_sing_soggetto{i}; 
        soggetto_decl{k,1}=codice_soggetto(i); 
        k=k+1; 
    end 
end 



%% Istogramma delle prove per soggetto

% codice_soggetto_abbreviato = cell(n_soggetti, 1);
% for i = 1:n_soggetti
%     soggetto_numero = extractAfter(codice_soggetto{i}, '_');
%     codice_soggetto_abbreviato{i} = ['sogg_' soggetto_numero];
% end
% 
% figure; 
% bar(perc_prove_soggetto*100); 
% xlabel('Soggetto'); 
% ylabel('Percentuale di prove (%)'); 
% title('Percentuale di prove per soggetto')
% set(gca, 'XTick', 1:n_soggetti, 'XTickLabel', codice_soggetto_abbreviato);


%tabellaRaw = table(perc_p_sogg_noNaN,  dMMSE_noNaN);
%% Test di Fisher su MMSE e perc_prove_soggetti
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
        prove_decl_neg=prove_decl_neg+n_sing(i); 
        
    elseif dMMSE_tot(i)>=0
        prove_nondec_pos=prove_nondec_pos+n_prove_nd(i); 
        prove_decl_pos=prove_decl_pos+n_sing(i); 
    
    end 
end 

% tabella_cont=[prove_decl_pos prove_decl_neg ; prove_nondec_pos prove_nondec_neg]; 
% nomi_righe2 = { 'Prove declinate', 'Prove non declinate'};
tabella_cont=[ prove_nondec_pos prove_nondec_neg; prove_decl_pos prove_decl_neg]; 
nomi_righe2 = {'Prove non declinate', 'Prove declinate'};
nomi_colonne2 = {'dMMSE >= 0', 'dMMSE < 0'};
tabella_con_titoli2 = array2table(tabella_cont, 'RowNames', nomi_righe2, 'VariableNames', nomi_colonne2);

[hp,pp,statsp] = fishertest(tabella_cont);






