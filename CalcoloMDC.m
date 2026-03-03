%% pulizia workspace
clearvars
clc

%% Leggo la tabella affidabilità dei vari task
task = 'Spiral';
root= '/Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/';
common_data_path = [root '/TesiTriennaleAnziani/TabelleSimone/'];


data=readtable([common_data_path 'MatriceAffidabilità' task '.xlsx']);
struttura_excel = table2struct(data);

%% Calcolo MDC

numero_indicatori_task=size(struttura_excel,1);

indicatori={};
mdc=zeros(numero_indicatori_task, 1);
affidabilita=zeros(numero_indicatori_task, 1);
SEM=zeros(numero_indicatori_task,1);
flag_mdc=zeros(numero_indicatori_task,1);
valori_Kendall=[];
valori_ICC=[];
k=1;
for i=1:numero_indicatori_task

    %creo colonne per foglio excel

    %costruisco un vettore con le affidabilità
    if struttura_excel(i).W_Kendall_>0
        affidabilita(i)=struttura_excel(i).W_Kendall_;
    elseif struttura_excel(i).r_ICC_>0
        affidabilita(i)=struttura_excel(i).r_ICC_;
    end

    %estraggo indicatori affidabili (buona e eccellente)
    if struttura_excel(i).W_Kendall_>=0.40
        flag_mdc(i)=1;
    elseif struttura_excel(i).r_ICC_>=0.50
        flag_mdc(i)=1;
    else
        flag_mdc(i)=0;
    end

    %calcolo mdc se affidabili
    if flag_mdc(i)==1
        SEM(i)=struttura_excel(i).DeviazioneStandard*(sqrt(1-affidabilita(i)));
        mdc(i)=SEM(i,1)*1.96*sqrt(2);
        cell_mdc(k,1) = mdc(i);
        valori_Kendall(k,1)=struttura_excel(i).W_Kendall_;
        valori_ICC(k,1)=struttura_excel(i).r_ICC_;
        indicatori{k,1}=struttura_excel(i).Indicatori;
        k=k+1;
    else
        mdc(i)=0;
    end
end

cell_mdc=num2cell(cell_mdc);
cell_W=num2cell(valori_Kendall);
cell_r=num2cell(valori_ICC);

matrice_tot=[indicatori cell_W cell_r cell_mdc];

%% salvo su tabella excel
save_folder = fullfile([root 'TesiTriennaleAnziani/TabelleMDCSimone/']);
if ~isfolder(save_folder)
    mkdir(save_folder)
end
data_table=cell2table(matrice_tot, "VariableNames", ["Indicatori" "AffidabilitàKendall" "AffidabilitàICC" "MDC"]);
filename='MDC';
save_name_table = [filename task '.xlsx'];
writetable(data_table, fullfile([save_folder '/' save_name_table]));

