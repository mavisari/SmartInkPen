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
    save_folder = fullfile('/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/AnalisiStatisticaTabelle');
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
fileDirezionePatologia=readtable('DirezionePatologia/DirezionePatologia_Spiral.xlsx'); 
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

for i = 1:n_soggetti
    codice_soggetto{i,1} = ['soggetto_' num2str(i)];
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

%% costruzione matrice_declinata_totale e matrice_non_declinata_totale

n_prove_sing_soggetto=zeros(n_soggetti,1);  
fields= fieldnames(struttura_soggetto);


for i=1:n_soggetti
    prove_singole=[]; 
    prove_gen=[]; 
    prove_s=[]; 

    for j=1:size(fileDirezionePatologia,1)

        if strcmp(fileDirezionePatologia(j).Soggetti, codice_soggetto(i))
            prove_singole = [prove_singole, fileDirezionePatologia(j).Prova];
           
        end 

        prove_singole_uniche = unique(prove_singole); 
        n_prove_sing_soggetto(i) = numel(prove_singole_uniche); 
        
        soggetto=fields{i}; 
        tabella_soggetto=struttura_soggetto.(soggetto);

        for k=1:n_prove_sing_soggetto(i)
            indice_da_estrarre = prove_singole_uniche (k); 
            matrice_declinata.(codice_soggetto{i})(k,1)=tabella_soggetto(indice_da_estrarre); 
        end 
    end 
    
    prove_gen=3:numero_prove_soggetto(i);
    comuni=ismember(prove_gen, prove_singole_uniche);
    noncomuni=~comuni;
    prove_s=prove_gen(noncomuni);
    lung_prova_s(i)=numel(prove_s);

    for z=1:lung_prova_s(i)
        indici_nod=prove_s(z);
        matrice_nondeclinata.(codice_soggetto{i})(z,1)=tabella_soggetto(indici_nod);
        
    end
end 

%% Matrici totali 
soggetti_declinati=fieldnames(matrice_declinata); 
soggetti_non_declinati=fieldnames(matrice_nondeclinata); 

valoriCampi=zeros(1, length(nomi_campi)); 
matrice_gen=struct();
matrice_declinata_totale=struct(); 
matrice_gen_nd=struct(); 
matrice_nondeclinata_totale=struct();

for i = 1:length(nomi_campi)
  matrice_gen.(nomi_campi{i}) = valoriCampi(i);
  matrice_declinata_totale.(nomi_campi{i}) = valoriCampi(i);
  matrice_gen_nd.(nomi_campi{i}) = valoriCampi(i);
  matrice_nondeclinata_totale.(nomi_campi{i}) = valoriCampi(i);
end


for i=1:size(soggetti_declinati)
    soggetto_estratto=soggetti_declinati{i}; 
    tabella_soggetto_declinata=matrice_declinata.(soggetto_estratto); 
    matrice_gen=[matrice_gen; tabella_soggetto_declinata]; 

end 

matrice_declinata_totale=[matrice_declinata_totale; matrice_gen];
matrice_declinata_totale = matrice_declinata_totale(3:end, :);


for i=1:size(soggetti_non_declinati)
    soggetto_estratto_nd=soggetti_non_declinati{i}; 
    tabella_soggetto_nondeclinata=matrice_nondeclinata.(soggetto_estratto_nd); 
    matrice_gen_nd=[matrice_gen_nd; tabella_soggetto_nondeclinata]; 

end 

matrice_nondeclinata_totale=[matrice_nondeclinata_totale; matrice_gen_nd];
matrice_nondeclinata_totale = matrice_nondeclinata_totale(3:end, :);


%% calcolo normalità, t-test e Mann-Whiteney test
flag_normale=[]; 

for i=1:(numero_indicatori-1)
    Indicatore=indicatori{i}; 
    colonna_dec=[]; 
    colonna_nd=[]; 

    %ESTRAZIONE COLONNA INDICATORE DA MATRICI 
    for j=1:size(matrice_declinata_totale)
        colonna_dec(j,1)=matrice_declinata_totale(j).(Indicatore); 
    end 

    for k=1:size(matrice_nondeclinata_totale)
        colonna_nd(k,1)=matrice_nondeclinata_totale(k).(Indicatore); 
    end 
    
    indici_NaN=isnan(colonna_dec); 
    colonna_dec(indici_NaN)=0; 

    indici_NaN=isnan(colonna_nd); 
    colonna_nd(indici_NaN)=0; 

    %TEST NORMALITA
    colonna_conc=[colonna_dec; colonna_nd]; 

    [hl,pl]=lillietest(colonna_conc, 'Alpha',0.05); 
    if pl<0.05
        flag_normale(i)=0; %non normale
    else 
        flag_normale(i)=1; %normale
    end

    %SE NORMALE UMPAIRED t-test 
    if flag_normale(i)==1
        [ht, pt, cit, statst] = ttest2(colonna_dec, colonna_nd);
        pvalue_ttest(i,1)=pt;

        media_dec(i,1)=mean(colonna_dec);
        media_nd(i,1)=mean(colonna_nd);

        SD_d(i,1)=std(colonna_dec);
        SD_nd(i,1)=std(colonna_nd);
    else 
         pvalue_ttest(i,1)=0;

        media_dec(i,1)=0;
        media_nd(i,1)=0;

        SD_d(i,1)=0;
        SD_nd(i,1)=0;
    end
  

    %SE NON NORMALE Mann-Whiteney test
    if flag_normale(i)==0
       [pr, ranksum_stat, ranksum_data] = ranksum(colonna_dec, colonna_nd);
       pvalue_ranksum(i,1)=pr; 

       mediana_d(i,1)=median(colonna_dec); 
       mediana_nd(i,1)=median(colonna_nd); 

       IQR_d(i,1)=iqr(colonna_dec); 
       IQR_nd(i,1)=iqr(colonna_nd); 
    
    else 
        pvalue_ranksum(i,1)=0;

        mediana_d(i,1)=0;
        mediana_nd(i,1)=0;

        IQR_d(i,1)=0;
        IQR_nd(i,1)=0;

    end 

end 

pvalue_ttest=num2cell(pvalue_ttest); 
pvalue_ranksum=num2cell(pvalue_ranksum); 

media_dec=num2cell(media_dec); 
media_nd=num2cell(media_nd); 

mediana_nd=num2cell(mediana_nd); 
mediana_d=num2cell(mediana_d); 

SD_d=num2cell(SD_d); 
SD_nd=num2cell(SD_nd); 

IQR_d=num2cell(IQR_d); 
IQR_nd=num2cell(IQR_nd); 

tot_ind=indicatori(1:numero_indicatori-1); 

totale=[tot_ind pvalue_ttest pvalue_ranksum media_dec media_nd SD_d SD_nd mediana_d mediana_nd IQR_d IQR_nd]; 

%% Salvo su excel 
data_table=cell2table(totale, "VariableNames", ["Indicatori" "P-value t-test" "P-value Mann-Whiteney test " "Media Prove Declinate" "Media Prove Non Declinate" "Deviazione Standard Prove declinate" "Deviazione Standard Prove non declinate" "Mediana Prove declinate" "Mediana Prove non declinate" "Range interquartile Prove declinate" "Range interquartile Prove non declinate"]);  
filename='Stat'; 
save_name_table = [filename '_' task '.xlsx'];
writetable(data_table, fullfile([save_folder '/' save_name_table])); 
