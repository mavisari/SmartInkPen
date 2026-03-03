% pulizia workspace
clearvars
clc
%% seleziono prove da considerare per calcolo reliability
n_prove_selezionate = 2;
% task da analizzare
task= 'Spiral';
%% definizione e aggiunta percorsi
utente = 'simone';
utente = upper(utente);

if contains(utente, 'SIMONE')
    root= '/Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/';
    common_data_path = [root '/TesiTriennaleAnziani/IndicatoriSimone/'];
    stats_path = [root '/TesiTriennaleAnziani/utilities_statistica/'];
    reliability_path = [root '/TesiTriennaleAnziani/TabelleSimone/'];
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
addpath(genpath(stats_path));

%% estrazione elenco dei dati presenti nella cartella specificata (data_folder)
file_list= dir(fullfile(common_data_path, '*.mat')); %cerco file excel nella cartella data_folder
n_file = length(file_list); %ottengo numero file mat nella cartella data_folder
file = cell(n_file, 1); %crea una cella per ogni soggetto nella cartella data_folder
for i = 1:n_file
    file{i} = strrep(file_list(i).name, '.mat', ''); %in ogni cella mantiene solo nome file senza estensione .mat
end
%% leggere il file excel
indice = find(contains(file, task));
% fileexcel = readtable([file{indice}]); fileexcel = table2struct(fileexcel);
load(file{indice});
fileexcel=S;
fields_to_remove = {'median_TILT_2'; 'median_TILT_3'; 'median_TILT_4'; 'median_TILT_5'};
for i = 1:size(fields_to_remove,1)
    if isfield(fileexcel, fields_to_remove{i})
        fileexcel=rmfield(fileexcel,fields_to_remove{i});
    end
end
%% estrazione N prove per soggetto
[unique_subjects, IA, IC] = unique({fileexcel.subject}.');
numero_prove_soggetto = groupcounts({fileexcel.subject}.');


for i = 1:max(numero_prove_soggetto)
    numero_prova{i,1} = ['prova_' num2str(i)];
end



if n_prove_selezionate > max(numero_prove_soggetto)
   n_prove_selezionate = max(numero_prove_soggetto);
end


clear struttura_prove
k=1;
for i=1:length(unique_subjects)

    if numero_prove_soggetto(i) >= n_prove_selezionate
       
        for j = 1:n_prove_selezionate
            struttura_prove.(numero_prova{j})(k,1) = fileexcel(IA(i)+j-1);
        end
        
        k=k+1;

    end

end

%% controllo sui soggetti
fields = fieldnames(struttura_prove);
for i = 1:length(struttura_prove.(fields{1,1}))
   
   array_soggetto = {};

   for j = 1:length(fields)
        array_soggetto{j,1} = struttura_prove.(fields{j,1})(i).subject;

   end

   if length(unique(array_soggetto)) > 1
       error('No corrispondenza tra soggetti nelle due prove')
   end

end
%% analisi di affidabilità

%indicators = fieldnames(fileexcel);
indicators = fieldnames(struttura_prove.prova_1); %PROVA MAVI
N = 216; % perchè gli indicatori di interesse sono fino al 216

% ad ogni iterazione calcolo affidabilità per l'indicatore i-esimo
n_soggetti=size(struttura_prove.prova_1, 1);  
SD=zeros(N,1); %deviazione standard per calcolo mdc (round 3)

flag_normale=zeros(N,1); 
flag_ttest=zeros(N,1); 
flag_Wilcoxon=zeros(N,1); 

array_W=zeros(N,1); 
array_r=zeros(N,1); 

pvalue_ttest=zeros(N,1); 
pvalue_Kendall=zeros(N,1);
pvalue_ICC=zeros(N,1);
pvalue_Wilcoxon=zeros(N,1); 

%%
for i = 1:N
    matrice=zeros(n_soggetti,2);

    % ESTRAETE L'INDICATORE DA struttura_prove.prova_1 ecc...
    Indicatore=indicators{i}; 
    for j=1:n_soggetti
         matrice(j,1)=struttura_prove.prova_1(j).(Indicatore);
         matrice(j,2)=struttura_prove.prova_2(j).(Indicatore);   
    end 
    
    indici_NaN=isnan(matrice); 
    matrice(indici_NaN)=0; 
    
    SD(i)=std(matrice(:,1)); 

    matrice_conc=[matrice(:,1); matrice(:,2)]; %creo una colonna per test di normalità


    % TEST di Normalità
    [hl,pl]=lillietest(matrice_conc, 'Alpha',0.05); 
    if pl<0.05 %riufiuto la normalità se il pvalue è sotto il livello di significatività
        flag_normale(i)=0; %rifiuto l'ipotesi nulla di normalità
    else 
        flag_normale(i)=1; %non rifiuto l'ipotesi nulla di normalità
    end
   
    % TEST di ipotesi
    if flag_normale(i)==1
        [ht, pt]=ttest(matrice(:,1), matrice(:,2));
        pvalue_ttest(i,1)=pt; 

        if pt<0.05
            flag_ttest(i)=0; %rifiuto ipotesi nulla 
        else 
            flag_ttest(i)=1; %accetto l'ipotesi nulla 
        end 
    end 


    % TEST di affidabilità [ICC o Kendall W] solò se test di ipotesi non è
    % stato rifiutato
    if flag_ttest(i)==1 
        [r, LB, UB, F, df1, df2, p] = ICC(matrice, 'A-k', 0.05);
        array_r(i,1)=r;
        pvalue_ICC(i,1)=p; 
    end 

%AGGIUNGERE TEST WILCOXON SIGNED RANK TEST se non normale
    if flag_normale(i)==0
        p_wx=signrank(matrice(:,1), matrice(:,2)); 
        pvalue_Wilcoxon(i,1)=p_wx; 

        if p_wx<0.05 
            flag_Wilcoxon(i)=0; %rifiuto l'ipotesi nulla
        else 
            flag_Wilcoxon(i)=1; %accetto ipotesi nulla
        end 
    end 

    if flag_Wilcoxon(i)==1 %se accetto ipotesi nulla allora calcolo Kendall per non normalità
        [W, p_val_q] = KendallCoef_p_val(matrice); 
         array_W(i,1)=W;
         pvalue_Kendall(i,1)=p_val_q; 

    end 
             
end 

indicatori_m=indicators(1:N); 
cell_sd=num2cell(SD); 

cell_W=num2cell(array_W); 
cell_r=num2cell(array_r); 

cell_pvalueW=num2cell(pvalue_Kendall); 
cell_pvaluer=num2cell(pvalue_ICC); 
cell_pvaluettest=num2cell(pvalue_ttest); 
cell_pvalueWilcoxon=num2cell(pvalue_Wilcoxon); 

matrice_affidabilita=[indicatori_m cell_pvaluettest cell_pvalueWilcoxon cell_W cell_pvalueW cell_r cell_pvaluer cell_sd];  


%% Salvo su file

save_folder = fullfile(reliability_path);
if ~isfolder(save_folder)
    mkdir(save_folder)
end
data_table=cell2table(matrice_affidabilita, "VariableNames", ["Indicatori" "P-value t-test" "P-value Wilcoxon" "W (Kendall)" "P-value Kendall" "r (ICC)" "P-value ICC" "Deviazione standard"]);  
filename='MatriceAffidabilità'; 
save_name_table = [filename task '.xlsx'];
writetable(data_table, fullfile([save_folder '/' save_name_table])); 
