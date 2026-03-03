function [PEN_test, signals, freq_g,psd_G_d,freq_gm,psd_G_dm,G_Peak_F_max_m,G_Peak_F_m,psd_G_PC1_d,psd_G_PC_d,psd_G_PC_sum_d,G_PC1_peak_F_d,G_PC_peak_F_d,G_PC_sum_peak_F_d]= ...
    FunctionIndicatorsComputation(data, segment_length, stroke_length,percentage, fc, r, amp, lam0, beads_flag, cost, pen_type, isLeLe)
% - baseline removal section compressed - V3
% - code re-organized
% - added asselborn spectral features
%% CONSTANTS
%% COSTANTI PER FAR PARTIRE IL CODICE DI PROVA
isLeLe = 1;
data = json_to_mat_new_pen('data_example.json');
percentage = 0.2;
pen_type = 'NEW';
segment_length = 500;
stroke_length = 12;
%%

%constants for displaying plots
display_Sadikov = 0;display = 0;display_harmonicity = 0;
%constants for ApEn computation
m=2;r_apen=0.2;
%constants for filters cut-off frequency
WnLdir=3; WnL=12; WnH=0.5; lf=0.01; hf=15; hfGyro =20;
%sampling frequency
fs=50;
%threshold for in air moments to be considered as pause, in seconds
threshpause=2;
if isLeLe==1
    [~, pressNcrop, ~, indPressure, indP]= baseline3 (data,percentage, pen_type);
else
    [indPressure, indP, pressNcrop] = extract_handwriting_task(data,fc, r, amp, lam0, cost, beads_flag);
end

%constants for removal of samples at the beginning and at the end of the signal
a=0; b=0; % a = beginning, b=end
P_filt=pressNcrop(indPressure(1)+a:indPressure(length(indP))-b);

if length(P_filt) < segment_length
    segment_length=length(P_filt);
    if mod(segment_length,2)
        segment_length=segment_length-1;
    end
end


%samples overlap for PSD estimation
overlap = ceil(segment_length/2); % 50%

%Frequency bands
bandt=round(segment_length/fs)*2+1; %21; %2Hz

%bande per giro - vecchie
GbandVol =round([(segment_length/fs).*0+1, (segment_length/fs).*1+1]); %[1 11]; %continua - 1Hz
GbandPC = round([(segment_length/fs).*1+2, (segment_length/fs).*3+1]); %[12 31];%1.1 - 3Hz
GbandPD = round([(segment_length/fs).*4+1, (segment_length/fs).*7.9+1]); %[41 71]; %4 - 7Hz
GbandPhy =round([(segment_length/fs).*8+1, (segment_length/fs).*12+1]); %[81 121]; %8 - 12Hz
%Bande Nuove - Van Galen (1990)
band_2=round([(segment_length/fs).*1+1, (segment_length/fs).*4+1]);%1-4
band_5=round([(segment_length/fs).*4+2, (segment_length/fs).*7+1]);%4.1-7
band_8=round([(segment_length/fs).*7+2, (segment_length/fs).*10+1]);%7.1-10
band_11=round([(segment_length/fs).*10+2, (segment_length/fs).*13+1]);%10.1-13
%% BASELINE REMOVAL AND EXTRACTION OF INDEXES WHERE PRESSURE > 0
% the beads parameters, the cost and the flag are passed by the scripts
% which computes the indicators for all files in a folder
%
%Baseline removal depending on the test type

%clearvars -except segment_length data
%% SIGNALS EXTRACTION
Time_adj=data.time/1000;
Time_adj=Time_adj(indP(1)+a : indP(length(indP))-b);
Time_adj=Time_adj-Time_adj(1);
G_adj=data.angVel(indP(1)+a:indP(length(indP))-b,:);
A=Filtering(data.acc(indP(1):indP(length(indP)),:),fs,lf,hf); % acc in mm/s^2 - filtrato con BP  0.1-15Hz
A_Smooth = Filtering(data.acc(indP(1):indP(length(indP)),:),fs,0.5,hf)./1000; %acc in m/s^2 - filttrato con BP 0.5-15Hz
Tremor = Filtering(data.acc(indP(1):indP(end),:),fs,2,hf); %filtrato con BP 2-15 Hz
G_f=Filtering(data.angVel(indP(1):indP(length(indP)),:),50,0.5,hf);
G_s= LPF(data.angVel(indP(1):indP(length(indP)),:),hf);
G_manus=Filtering(data.angVel(indP(1):indP(length(indP)),:),50,2,hf);

meanP=mean(nonzeros(P_filt));
A3D=sqrt(A(:,1).^2 + A(:,2).^2 + A(:,3).^2);
A3Dd=detrend(A3D,0); A3Df=Filtering(A3D,50, 0.5, hf);
%% PC EXTRACTION
% ACCELERAZIONE FILTRATA PASSA ALTO A 2HZ
[~,score,~,~,explained,~]=pca(Tremor);
%PRENDE LE PRIME DUE COMPONENTI PRINCIPALI - S
PC1D=score(:,1); PC2D=score(:,2); PC_sum = PC1D+PC2D;
if explained(1) >= 67
    PC = PC1D;
else
    PC = PC1D+PC2D;
end


if length(Tremor)>= segment_length
    [PSD_tremor,f] = pwelch(Tremor,segment_length, overlap, segment_length, fs);
    [PSD_PC1,f] = pwelch(PC1D,segment_length, overlap, segment_length, fs);
    [PSD_PC,f] = pwelch(PC,segment_length, overlap, segment_length, fs);
    if display
        figure
        subplot 211; plot(f, PSD_tremor); title('Acc HP 2Hz single components - PSD')
        subplot 212; plot(f, PSD_PC1); hold on; plot(f,PSD_PC);title('Acc HP 2Hz PC - PSD'); legend('PC1', 'PC > 67%')
    end
     ind_2hz=find(f==2)+1;
end
   

    % VELOCITA' ANGOLARE FILTRATA PASSA ALTO A 2HZ
    [~,scoreG_tr,~,~,explainedG_tr,~]=pca(G_manus);
    PC1G_tr=scoreG_tr(:,1); PC2G_tr=scoreG_tr(:,2); PCG_tr_sum = PC1G_tr+PC2G_tr;
    if explainedG_tr(1) >= 67
        PCG_tr = PC1G_tr;
    else
        PCG_tr = PC1G_tr+PC2G_tr;
    end
if length(Tremor)>= segment_length
    [PSDG_tr_comp,f] = pwelch(G_manus,segment_length, overlap, segment_length, fs);
    [PSDG_tr_PC,f] = pwelch(PCG_tr,segment_length, overlap, segment_length, fs);
    [PSDG_tr_PC1,f] = pwelch(PC1G_tr,segment_length, overlap, segment_length, fs);
    if display
        figure
        subplot 211; plot(f, PSDG_tr_comp); title('Gyro single components - 2Hz HP')
        subplot 212; hold on; plot(f, PSDG_tr_PC1); plot(f, PSDG_tr_PC); title('Gyro PC - 2Hz HP'); legend('PC1', 'PC > 67%')
    end
end


%% accelerazione e velocità per calcolo jerk
bpx=LPF(A(:,1),WnLdir);bpy=LPF(A(:,2),WnLdir);bpz=LPF(A(:,3),WnLdir);
% bpA3D = sqrt(bpx.^2+bpy.^2+bpz.^2);
% bpA3Df=Filtering(sqrt(bpx.^2+bpy.^2+bpz.^2), fs, WnH, WnLdir);
bpA3D1 = Filtering(A3D, 50, WnH, WnLdir);%uso questa x il jerk
%velox calcola la velocitï¿½ a partire da acc band-passed
V_x=velox(Time_adj,bpx); V_y=velox(Time_adj,bpy);V_z=velox(Time_adj,bpz);
V=sqrt(V_x.^2+V_y.^2+V_z.^2); %USO QUESTA PER IL JERK
% Vd=velox(Time_adj,bpA3Dd);
% Vf=velox(Time_adj,bpA3Df);
%% CREATION OF A MATRIX WITH ALL NEEDED SIGNALS
CELLDATA(:,1)=Time_adj; %TIME in seconds
CELLDATA(:,2)=A(:,1);CELLDATA(:,3)=A(:,2);CELLDATA(:,4)=A(:,3); %BP filtered 0.1-15Hz acc in mm/s^2
CELLDATA(:,5)=G_adj(:,1); CELLDATA(:,6)=G_adj(:,2); CELLDATA(:,7)=G_adj(:,3); %ang vel not filtered
CELLDATA(:,8)=P_filt; %baseline removed pressure
CELLDATA(:,9)=bpA3D1; %BP filtered 3D acc
CELLDATA(:,10)=V; % 3D velocity
CELLDATA(:,11)=A3Dd; % mean detrended 3D acc
CELLDATA(:,12) = HPF(LPF(A3D,3.5),0.1); %other version of BP filtered 3D acc (???)
CELLDATA(:,13) = A_Smooth(:,1);CELLDATA(:,14) = A_Smooth(:,2);CELLDATA(:,15) = A_Smooth(:,3); %acc in m/s^2
CELLDATA(:,16) = G_f(:,1);CELLDATA(:,17) = G_f(:,2);CELLDATA(:,18) = G_f(:,3); % BP filtered 0.5-15Hz ang vel
CELLDATA(:,19) = G_s(:,1);CELLDATA(:,20) = G_s(:,2);CELLDATA(:,21) = G_s(:,3); % LP filtered 15Hz ang vel
CELLDATA(:,22) = G_manus(:,1);CELLDATA(:,23) = G_manus(:,2);CELLDATA(:,24) = G_manus(:,3); % BP filtered 2-15Hz ang vel
CELLDATA(:,25) = Tremor(:,1); CELLDATA(:,26) = Tremor(:,2); CELLDATA(:,27) = Tremor(:,3); %singola componente acc filtrata HP a 2Hz
CELLDATA(:,28) = PC1D; % PC1 acc tr
CELLDATA(:,29) = PC_sum; % PC1+PC2 acc tr
CELLDATA(:,30) = PC; % either PC1 or PC_sum acc tr, according to the amount of explained variance
CELLDATA(:,31) = PC1G_tr; % PC1 gyro tremor
CELLDATA(:,32) = PCG_tr_sum; % % PC1+PC2 gyro tr
CELLDATA(:,33) = PCG_tr; % either PC1 or PC_sum gyro tr, according to the amount of explained variance
%not filtered acc in mm/s^2 -> for TILT
CELLDATA(:,34) = data.acc(indP(1):indP(end),1);CELLDATA(:,35) = data.acc(indP(1):indP(end),2);CELLDATA(:,36) = data.acc(indP(1):indP(end),3);
signals=CELLDATA(:,2:end);
%% check data loss
deltat=diff(data.timestamp(indP(1):indP(length(indP))));%differenza tra timestamp consecutivi
cut=[]; %array in cui salvo gli indici in cui il segnale va tagliato
i=1;
while i < length(deltat)
    if (deltat(i) == 6 && deltat(i+1)~= 4) && (deltat(i+1) + deltat(i+2)~= 9)
        cut=[cut; i+1]; i=i+1;
    elseif deltat(i)>= 7
        cut=[cut; i+1]; i=i+1;
    end
    i=i+1;
end
% divido il segnale in segmenti (DataCut) escludendo i campioni persi
if ~isempty(cut)
    start=1;
    for i=1:length(cut)+1
        if i == length(cut)+1
            DataCut{i,1} = CELLDATA(start:end, :);
        else
            DataCut{i,1} = CELLDATA(start:cut(i)-1, :);
            start=cut(i)+1;
        end
    end %fine for su lunghezza indici di perdita
else
    DataCut{1,1} = CELLDATA(:,:); %se non ci sono perdite considero tutto il segnale
end %fine divisione del segnale in segmenti basata sulle perdite
%fine divisione del segnale in segmenti basata sulle perdite
%Confronto Press - Cut e non Cut
% FirstCut = DataCut{1,1};
% figure
% plot(P_filt, 'Color', 'g'); hold on
% figure
% plot(FirstCut(:,8))
% for i=1:length(cut)
%   xline(cut(i), 'LineWidth', 1.25);
% end
%% STUDIO TREMORE
%500 -> 0.1Hz di risoluzione in frequenza
% [TM_ApEn,TM_RR,TM_DET,RMS_A_TOT,RMS_A_VOL,RMS_A_PD,RMS_A_PHY,RMS_G_TOT,RMS_G_VOL,RMS_G_PD,RMS_G_PHY,RMS_T_TOT,RMS_T_LF,RMS_T_PD,RMS_T_PHY,RPW_A_VOL,RPW_A_PD,RPW_A_PHY,RPW_G_VOL,RPW_G_PD,RPW_G_PHY,RPW_T_LF,RPW_T_PD,RPW_T_PHY,PEAK_POW_A,F_MODAL_A, PEAK_POW_T, F_MODAL_T, F_MODAL_T_IRENE]= HandwritingTREMORstudy(DataCut,duration);
[TM_ApEn,TM_RR,TM_DET,TSI_DB_IMF,TSI_DB_thr_IMF,SNR_IMF, MHP_IMF,TSI_Luft_IMF, ...
    RMS_A_TOT,RMS_A_VOL,RMS_A_DYS,RMS_A_PD,RMS_A_PHY,RMS_G_TOT,~,RMS_G_DYS,RMS_G_PD,RMS_G_PHY,RMS_T_TOT,RMS_T_LF,RMS_T_PD,RMS_T_PHY,...
    RPW_A_VOL, RPW_A_DYS,RPW_A_PD,RPW_A_PHY,~,RPW_G_DYS,RPW_G_PD,RPW_G_PHY,RPW_T_LF,RPW_T_PD,RPW_T_PHY,...
    PEAK_POW_A,F_MODAL_A, PEAK_POW_T, F_MODAL_T, F_MODAL_T_IRENE, PEAK_POW_G, F_MODAL_G]= HandwritingTREMORstudy_2(DataCut,segment_length);
%% FIND STROKES INDEXES
ind_write = find_strokes_indexes(CELLDATA(:,8)); %we consider the pression
% execution time and stroke number
Execution_Time = CELLDATA(ind_write(end),1) - CELLDATA(ind_write(1),1);
Num_Strokes=(length(ind_write)/2);
Relative_Num_Strokes = (length(ind_write)/2)/Execution_Time;

%% inizio analisi su singolo elemento di DataCut
for d = 1 : length(DataCut)
    clear TEMPO G_X G_Y G_Z PRESSURE ACC3D ACC_X ACC_Y ACC_Z
    TEMPO=DataCut{d,1}(:,1);
    ACC_X=DataCut{d,1}(:,2); ACC_Y=DataCut{d,1}(:,3); ACC_Z=DataCut{d,1}(:,4);
    G_X=DataCut{d,1}(:,5); G_Y=DataCut{d,1}(:,6);G_Z=DataCut{d,1}(:,7);
    G_Xm=DataCut{d,1}(:,22); G_Ym=DataCut{d,1}(:,23); G_Zm=DataCut{d,1}(:,24);
    A_Xtr=DataCut{d,1}(:,25); A_Ytr=DataCut{d,1}(:,26); A_Ztr=DataCut{d,1}(:,27);
    PRESSURE=DataCut{d,1}(:,8);
    %     ACC3DJ=DataCut{d,1}(:,9); V=DataCut{d,1}(:,10);
    ACC3D=DataCut{d,1}(:,11);
    %trovo indici strokes del tratto senza perdite
    ind_DC{d}=find(PRESSURE~=0);
    if length(ind_DC{d}) >=2 %continuo solo se ho almeno uno stroke nel tratto senza perdite

        %extract strokes start-stop indexes in datacut
        ind_write_DC{d} = find_strokes_indexes(PRESSURE); %same of ind_write if not data loss
        ind_write_DC{d} = ind_write_DC{d}.'; %trasposto

        % ||||||| SALVATAGGIO HANDWRITING STROKES |||||||
        %Divisione in stroke. Ogni cella di stroke ho TUTTI i dati di uno stroke
        coppie(d)=size(ind_write_DC{d},2)/2;
        on=1; in=1;
        for j=1:2:size(ind_write_DC{d},2)
            %if ind_write_DC{d}(j+1)-ind_write_DC{d}(j)>= 6
            stroke=DataCut{d,1}(ind_write_DC{d}(j):ind_write_DC{d}(j+1),:);
            % Salvo gli in-aire moments e la loro durata.

            if j<(size(ind_write_DC{d},2)-1)
                inair{d,in}=DataCut{d,1}(ind_write_DC{d}(j+1)+1:ind_write_DC{d}(j+2)-1,:); %valori dei segnali corrispondenti all'i-esimo in-air
                inair_duration(d,in)=inair{d,in}(end,1)-inair{d,in}(1,1); %durata in secondi dell'i-esimo in-air
                in=in+1;
            end
            strokes{d,on}=stroke; %in ogni cella ho i segnali del singolo stroke
            onsheet(d,on)=strokes{d,on}(end,1)-strokes{d,on}(1,1); %stroke duration
            on=on+1;
            %end
            %clear stroke
        end
        clear j;
        %in uscita ho stroke, onsheet e inair
        %% INIZIO ANALISI SUL SINGOLO STROKE
        for L=1:coppie(d)
            if size(strokes{d,L},1)>=stroke_length
                % DEFINISCO: Acc XYZ del singolo stroke
                tempo=strokes{d,L}(:,1)-strokes{d,L}(1,1);
                clear P; P=strokes{d,L}(:,8);
                g_x=strokes{d,L}(:,2); g_y=strokes{d,L}(:,3);   g_z=strokes{d,L}(:,4);
                [c_x,c_y,c_z,stimaG(L),errG(L)]=gravitystroke(g_x,g_y,g_z);
                %sostituisco l'acc con quella compensata (per ogni stroke)
                strokes{d,L}(:,2)=c_x; strokes{d,L}(:,3)=c_y; strokes{d,L}(:,4)=c_z;
                V_strokes = strokes{d,L}(:,10);
                %figure
                %subplot 411; plot(g_x); hold on; plot(c_x); yline(0);
                %subplot 412; plot(g_y); hold on; plot(c_y); yline(0);
                %subplot 413; plot(g_z); hold on, plot(c_z); yline(0);
                %subplot 414; plot(sqrt(g_x.^2+g_y.^2+g_z.^2)); hold on; plot(sqrt(c_x.^2+c_y.^2+c_z.^2)); legend('NC', 'C')

                %              AVGerrG(d)=mean(errG); AVGstima(d)=mean(stimaG);
                %              clear errG stimaG

                %%%%%%%%%%%%%%%new smoothness measures%%%%%%%%%%%%%%
                [LDJL_A_s{d,L}, LDJL_G_s{d,L}] = compute_Smoothness_DLJ(strokes{d,L}(:,13),strokes{d,L}(:,14),strokes{d,L}(:,15), ...
                    strokes{d,L}(:,16),strokes{d,L}(:,17),strokes{d,L}(:,18),tempo);

                [~, ~, S1_s{d,L}, S2_s{d,L}, S3_s{d,L}, S4_s{d,L}, S45_s{d,L}, S5_s{d,L}] = ...
                    compute_Smoothness_Stroke_from_kinematic(strokes{d,L}(:,13),strokes{d,L}(:,14),strokes{d,L}(:,15), ...
                    strokes{d,L}(:,16),strokes{d,L}(:,17),strokes{d,L}(:,18),...
                    strokes{d,L}(:,19),strokes{d,L}(:,20),strokes{d,L}(:,21),tempo);

                %%%%%%%%%%%%%%%Signal-to-Noise velocity peaks difference (SNvpd)%%%%%%%%%%%%%%
                bpxs{d,L}=HPF(LPF(strokes{d,L}(:,2),WnL),WnH); %HP per drift
                bpys{d,L}=HPF(LPF(strokes{d,L}(:,3),WnL),WnH);
                bpzs{d,L}=HPF(LPF(strokes{d,L}(:,4),WnL),WnH);
                %velox calcola la velocità a partire da acc band-passed
                V_xs{d,L}=velox(tempo,bpxs{d,L}); V_ys{d,L}=velox(tempo,bpys{d,L});  V_zs{d,L}=velox(tempo,bpzs{d,L});
                Vs{d,L}=sqrt(V_xs{d,L}.^2+V_ys{d,L}.^2+V_zs{d,L}.^2);

                Vs_Filt5{d,L}=LPF(Vs{d,L}, 5); %LP filtered velocity at 5Hz
                Vs_Filt10{d,L}=LPF(Vs{d,L},10);%LP filtered velocity at 10Hz

                [Peaks10{d,L},locs10]=findpeaks(Vs_Filt10{d,L}, 'MinPeakHeight', 0);
                [Peaks5{d,L},locs5]=findpeaks(Vs_Filt5{d,L}, 'MinPeakHeight', 0);

                SNvpd_norm(d,L)=(length(Peaks10{d,L})-length(Peaks5{d,L}))/(tempo(end-1)-tempo(1));

                %%%%%%%%%%%%%%%JERK%%%%%%%%%%%%%%
                GRAD{d,L}=gradient(strokes{d,L}(:,9), tempo);
                %figure
                %plot(GRAD{d,i})
                %title('Gradient')
                integral_a(d,L)=trapz(tempo, GRAD{d,L}.^2);%integrale(squared jerk)
                t_1_2_a(d,L)=tempo(end)-tempo(1);%parametro x normalizzazione numeratore
                DimLessJerk_A(d,L)= (integral_a(d,L))*(t_1_2_a(d,L)^3)/((mean(V_strokes))^2); %DLJ_A

                %%%%%%%%%%%%%%PRESSURE%%%%%%%%%%%%%%
                %             LPP=LPF(P,4);
                [MaximaP{d,L},MinimaP{d,L},maxExtremaP(d,L),avgExtremaP(d,L),NC_press(d,L), forcedNCP_s{d,L}, IndP_s{d,L}]=MinMax(P);
                NCP_REL(d,L) = NC_press(d,L)/onsheet(d,L); % NCP relativo singolo stroke
                P_AVG(d,L) = mean(P);  P_CV(d,L)= std(P)/P_AVG(d,L);  P_OVS(d,L)= max(P)-median(P);
                %Consecutive peak difference in pressure
                [peak_diff_pavg_s(d,L), peak_diff_pCV_s(d,L)] = consPeakDiff(P, IndP_s{d,L});

                %%%%%%%%%%%%%%GYRO%%%%%%%%%%%%%%
                clear G_rms G_S
                G_rms = [rms(strokes{d,L}(:,5)), rms(strokes{d,L}(:,6)), rms(strokes{d,L}(:,7))];
                Max_G_rms=max(G_rms);IND_G=find(G_rms==Max_G_rms);
                G_S=strokes{d,L}(:,IND_G(1)+4);
                clear Max_G_rms IND_G

                LPG = LPF(G_S,4);
                [MaximaG{d,L},MinimaG{d,L},maxExtremaG(d,L),avgExtremaG(d,L),NC_gyro(d,L), forcedNCG_s{d,L}, IndG_s{d,L}]=MinMax(LPG);
                NCG_REL(d,L) = NC_gyro(d,L)/onsheet(d,L); % NCG relativo singolo stroke
                [peak_diff_gavg_s(d,L), peak_diff_gCV_s(d,L)] = consPeakDiff(LPG, IndG_s{d,L});
                %%%%%%%%%%%%%%%%%% VELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [MaximaV{d,L},MinimaV{d,L},maxExtremaV(d,L),avgExtremaV(d,L),NC_vel(d,L), forcedNCV_s{d,L}, IndV_s{d,L}]=MinMax(V_strokes);
                NCV_REL(d,L) = NC_vel(d,L)/onsheet(d,L);
                V_AVG(d,L)=mean(V_strokes); V_CV(d,L)=std(V_strokes)/V_AVG(d,L);
                V_ApEn(d,L)=FastApEn(m,r_apen*std(V_strokes),V_strokes,1);
                [peak_diff_vavg_s(d,L), peak_diff_vCV_s(d,L)] = consPeakDiff(V_strokes, IndV_s{d,L});

                %%%%%%%%%%%%%%%%%% ACCELERAZIONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % SU ACCELERAZIONE TOTALE
                %WnLtot=4;
                % vedi contenuto in freq per decisione soglie si lpf
                %ALP_tot=LPF(A_tot{d,L} ,WnLtot);
                [MaximaALP{d,L},MinimaALP{d,L},maxExtremaALP(d,L),avgExtremaALP(d,L),NC_accLP(d,L)]=MinMax(strokes{d,L}(:,12));
                NCA_REL(d,L) = NC_accLP(d,L)/onsheet(d,L); % NCA relativo singolo stroke

                %%%%%%%%%%%%%%%%%% TILT %%%%%%%%%%%%%%%%%%
                %primo argomento = acc longitudinale a lunghezza penna
                %secondo e terzo argomento = gyro sugli altri due assi
                if contains(pen_type, 'OLD') %penna vecchia
                    theta = computeAttitudeOld(strokes{d,L}(:,35), strokes{d,L}(:,5), strokes{d,L}(:,7));
                else %penna nuova
                    theta = computeAttitudeOld(strokes{d,L}(:,34), strokes{d,L}(:,6), strokes{d,L}(:,7));
                end
                TILT_S(d,L) = mean(theta); %tilt medio dello stroke
                if TILT_S(d,L) == 0
                    TILT_S_CV(d,L) = 0;
                else
                    TILT_S_CV(d,L) = std(theta)/TILT_S(d,L);
                end
                TILT_S_VAR(d,L) = var(theta);

            else %se lo stroke è più corto di stroke length

                LDJL_A_s{d,L} = zeros(1,3); LDJL_G_s{d,L} = zeros(1,3);
                S1_s{d,L} = zeros(1,3); S2_s{d,L} = zeros(1,3);  S3_s{d,L} = zeros(1,3);
                S4_s{d,L} = zeros(1,3);  S45_s{d,L} = zeros(1,3);  S5_s{d,L} = zeros(1,3);

                SNvpd_norm(d,L)=0; integral_a(d,L)=0; DimLessJerk_A(d,L)=0;

                P_AVG(d,L) = 0; P_CV(d,L) = 0; P_OVS(d,L) = 0; peak_diff_pavg_s(d,L)=0; peak_diff_pCV_s(d,L)=0; NCP_REL(d,L)=0;

                NCG_REL(d,L)=0; peak_diff_gavg_s(d,L)=0; peak_diff_gCV_s(d,L)=0;

                NCA_REL(d,L)=0;
                TILT_S(d,L)=0;  TILT_S_CV(d,L)=0; TILT_S_VAR(d,L) = 0;

                NCV_REL(d,L)=0; V_AVG(d,L)=0;  V_CV(d,L)=0; V_ApEn(d,L)=0; peak_diff_vavg_s(d,L)=0; peak_diff_vCV_s(d,L)=0;
            end
        end
        % FINE ANALISI SINGOLO STROKE
        %% ANALISI SU DATACUT
        %figure; plot(sort(nonzeros(DimLessJerk_A(d,:))))
        %yline(median(nonzeros(DimLessJerk_A(d,:)))); yline(median(DimLessJerk_A(d,:)), 'Color', 'g'); yline(mean(nonzeros(DimLessJerk_A(d,:))), 'Color', 'r');
        %title('DLJ tratto senza perdite')

        %CALCOLO MEDIA VARIABILI x STROKE nel tratto senza perdite
        %Smoothness
        LDLJ_A_k(d,:) = computeStrokeMeanAxis(LDJL_A_s(d,:));
        LDLJ_G_k(d,:) = computeStrokeMeanAxis(LDJL_G_s(d,:));

        S1_sm(d,:) = computeStrokeMeanAxis(S1_s(d,:));
        S2_sm(d,:) = computeStrokeMeanAxis(S2_s(d,:));
        S3_sm(d,:) = computeStrokeMeanAxis(S3_s(d,:));
        S4_sm(d,:) = computeStrokeMeanAxis(S4_s(d,:));
        S45_sm(d,:) = computeStrokeMeanAxis(S45_s(d,:));
        S5_sm(d,:) = computeStrokeMeanAxis(S5_s(d,:));
        %se S li calcolo qui non deve più essere calcolata la media
        [S1_k(d,:), S2_k(d,:), S3_k(d,:), S4_k(d,:), S45_k(d,:), S5_k(d,:)] = compute_Smoothness_Sparc(DataCut{d,1}(:,19),DataCut{d,1}(:,20),DataCut{d,1}(:,21));

        clear RPWG RMSG m_manus

        meanP_DC(d) = computeStrokeMean(P_AVG(d,:)); meanP_CV_DC(d) = computeStrokeMean(P_CV(d,:));
        meanNCP_DC(d)=computeStrokeMean(NCP_REL(d,:)); meanPOVS_DC(d)=computeStrokeMean(P_OVS(d,:));
        meanNCG_DC(d)=computeStrokeMean(NCG_REL(d,:)); meanNCA_DC(d)=computeStrokeMean(NCA_REL(d,:)); meanNCV_DC(d)=computeStrokeMean(NCV_REL(d,:));
        meanV_DC(d)=computeStrokeMean(V_AVG(d,:)); V_CV_AVG(d)=computeStrokeMean(V_CV(d,:));  V_ApEn_AVG(d)=computeStrokeMean(V_ApEn(d,:));
        medianDLJ_A(d)=computeStrokeMedian(DimLessJerk_A(d,:)); medianSquaredJerk_A(d)=computeStrokeMedian(integral_a(d,:));

        SNvpd_norm_AVG(d)=computeStrokeMean(SNvpd_norm(d,:));
        meanCPDP_AVG(d)=computeStrokeMean(peak_diff_pavg_s(d,:));
        meanCPDP_CV(d)=computeStrokeMean(peak_diff_pCV_s(d,:));
        meanCPDG_AVG(d)=computeStrokeMean(peak_diff_gavg_s(d,:));
        meanCPDG_CV(d)=computeStrokeMean(peak_diff_gCV_s(d,:));
        meanCPDV_AVG(d)=computeStrokeMean(peak_diff_vavg_s(d,:));
        meanCPDV_CV(d)=computeStrokeMean(peak_diff_vCV_s(d,:));
        meanTILT(d)=computeStrokeMean(TILT_S(d,:));
        meanTILT_CV(d)=computeStrokeMean(TILT_S_CV(d,:));
        meanTILT_VAR(d)=computeStrokeMean(TILT_S_VAR(d,:));

        %%%%%% ANALISI TREMORE SU SEGMENTO SENZA PERDITE %%%%%%%
        L(d)=length(ACC_X);
        if L(d) >= segment_length
            %
            %%%TSI, SNR, MHP di ACC3D %%%
            [TSI_DB(d), TSI_DB_thr(d), SNR(d), MHP(d)] = TSI_computation(ACC3D,TEMPO);
            TSI_Luft(d) = TSI_Luft_computation(ACC3D);
            %feature manus
            ord = 10;
            [RPWG_manus_s(d,:), RMSG_manus_s(d,:), m_manus_s(d,:)] = feature_manus_segment(DataCut{d,1}(:,22), DataCut{d,1}(:,23), ...
                DataCut{d,1}(:,24),segment_length, 50, ord);

            clear A_d G_d A_tr_d G_dm
            %%SPECTRAL FEATURES CON WELCH%%%
            %%%%%%%%%%ACC  -> considero segnale filtrato passa alto a 0.1Hz%%%%%%%%
            A_d(:,1)= ACC_X; A_d(:,2)= ACC_Y; A_d(:,3)= ACC_Z;
            [psd_A_d{d}, freq_a{d}, RPW_A_Vol_d(d,:), RPW_A_PC_d(d,:), RPW_A_PD_d(d,:), RPW_A_Phy_d(d,:), A_peak_d(d,:), A_peak_F_d(d,:)] = ...
                computeSpectralFeaturesAxis(A_d, segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt,1);
            [psd_bw_A_d(d,:) , psd_median_A_d(d,:) ] = asselbornSpectralFeatures(psd_A_d{d}, freq_a{d});
            %%%%% new band %%%%%
            [psd_A_new{d}, freq_A_new{d}, RPW_2_A_ab_d(d,:), RPW_5_A_ab_d(d,:), RPW_8_A_ab_d(d,:), RPW_11_A_ab_d(d,:),RPW_2_A_rel_d(d,:),RPW_5_A_rel_d(d,:),RPW_8_A_rel_d(d,:),RPW_11_A_rel_d(d,:), A_peak_d_new(d,:), A_peak_F_d_new(d,:)] = ...
                computeSpectralFeatures_new(A_d, segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_A_d_new(d,:) , psd_median_A_d_new(d,:) ] = asselbornSpectralFeatures(psd_A_new{d}, freq_A_new{d});



            %%%%%%%%%%ACC_tr  -> considero segnale filtrato passa alto a 2Hz%%%%%%%%
            A_tr_d(:,1)= A_Xtr; A_tr_d(:,2)= A_Ytr; A_tr_d(:,3)= A_Ztr;
            [psd_A_tr_d{d}, freq_A_tr{d}, ~, RPW_A_tr_PC_d(d,:), RPW_A_tr_PD_d(d,:), RPW_A_tr_Phy_d(d,:), A_tr_peak_d(d,:), A_tr_peak_F_d(d,:)] = ...
                computeSpectralFeaturesAxis(A_tr_d, segment_length, overlap, fs, GbandVol, [bandt GbandPD(1)], GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_A_tr_d(d,:) , psd_median_A_tr_d(d,:) ] = asselbornSpectralFeatures(psd_A_tr_d{d}, freq_A_tr{d});
            %%%%% new band %%%%%
            [psd_A_tr_new{d}, freq_A_tr_new{d}, RPW_2_A_tr_ab_d(d,:), RPW_5_A_tr_ab_d(d,:), RPW_8_A_tr_ab_d(d,:), RPW_11_A_tr_ab_d(d,:),RPW_2_A_tr_rel_d(d,:),RPW_5_A_tr_rel_d(d,:),RPW_8_A_tr_rel_d(d,:),RPW_11_A_tr_rel_d(d,:), A_tr_peak_d_new(d,:), A_tr_peak_F_d_new(d,:)] = ...
                computeSpectralFeatures_new(A_tr_d, segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_A_tr_d_new(d,:) , psd_median_A_tr_d_new(d,:) ] = asselbornSpectralFeatures(psd_A_tr_new{d}, freq_A_tr_new{d});

            %%%%%%%%%%GYRO -> considero segnale non filtrato, nelle tre componenti%%%%%%%%
            G_d(:,1) = G_X; G_d(:,2) = G_Y; G_d(:,3) = G_Z;
            [psd_G_d{d}, freq_g{d}, RPW_G_Vol_d(d,:), RPW_G_PC_d(d,:), RPW_G_PD_d(d,:), RPW_G_Phy_d(d,:), G_peak_d(d,:), G_peak_F_d(d,:)] = ...
                computeSpectralFeaturesAxis(G_d, segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, GbandVol(2), GbandVol(1));
            [psd_bw_G_d(d,:) , psd_median_G_d(d,:) ] = asselbornSpectralFeatures(psd_G_d{d}, freq_g{d});
            %%%%% new band %%%%%
            [psd_G_d_new{d}, freq_G_d_new{d}, RPW_2_G_ab_d(d,:), RPW_5_G_ab_d(d,:), RPW_8_G_ab_d(d,:), RPW_11_G_ab_d(d,:),RPW_2_G_rel_d(d,:),RPW_5_G_rel_d(d,:),RPW_8_G_rel_d(d,:),RPW_11_G_rel_d(d,:), G_peak_d_new(d,:), G_peak_F_d_new(d,:)] = ...
                computeSpectralFeatures_new(G_d, segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_G_d_new(d,:) , psd_median_G_d_new(d,:) ] = asselbornSpectralFeatures(psd_G_d_new{d}, freq_G_d_new{d});

            %%%%%%%%%%GYRO_tr -> considero segnale filtrato passa alto a 2Hz, nelle tre componenti%%%%%%%%
            G_dm(:,1) = G_Xm; G_dm(:,2) = G_Ym; G_dm(:,3) = G_Zm;
            [psd_G_dm{d}, freq_gm{d}, ~, RPW_G_PC_dm(d,:), RPW_G_PD_dm(d,:), RPW_G_Phy_dm(d,:), G_peak_dm(d,:), G_peak_F_dm(d,:)] = ...
                computeSpectralFeaturesAxis(G_dm, segment_length, overlap, fs, GbandVol, [bandt GbandPD(1)], GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_G_dm(d,:) , psd_median_G_dm(d,:) ] = asselbornSpectralFeatures(psd_G_dm{d}, freq_gm{d});
            %%%%% new band %%%%%
            [psd_G_dm_new{d}, freq_G_dm_new{d}, RPW_2_G_dm_ab_d(d,:), RPW_5_G_dm_ab_d(d,:), RPW_8_G_dm_ab_d(d,:), RPW_11_G_dm_ab_d(d,:),RPW_2_G_dm_rel_d(d,:),RPW_5_G_dm_rel_d(d,:),RPW_8_G_dm_rel_d(d,:),RPW_11_G_dm_rel_d(d,:), G_dm_peak_d_new(d,:), G_dm_peak_F_d_new(d,:)] = ...
                computeSpectralFeatures_new(G_dm, segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_G_dm_new(d,:) , psd_median_G_dm_new(d,:) ] = asselbornSpectralFeatures(psd_G_dm_new{d}, freq_G_dm_new{d});

            %%%%%%%%%%NORMALIZED Spectra%%%%%%%%%%
            for k=1:3
                psd_A_tr_norm_d{d}(:,k) = pwelch(A_tr_d(:,k)./max(A_tr_d(:,k)), segment_length, overlap, segment_length, fs);
                A_tr_peak_norm_d(d,k)  = max(psd_A_tr_norm_d{d}(:,k));

                psd_G_norm_dm{d}(:,k) = pwelch(G_dm(:,k)./max(G_dm(:,k)), segment_length, overlap, segment_length, fs);
                G_peak_norm_dm(d,k) = max(psd_G_norm_dm{d}(:,k));
            end

            %%%%%%%% TILT %%%%%%%
            theta_data_cut = computeAttitudeOld(DataCut{d}(:,34), DataCut{d}(:,6), DataCut{d}(:,7));
            if length(theta_data_cut) >= segment_length
            [psd_tilt{d},freq_tilt{d}]=TiltSpectrum(theta_data_cut,segment_length,overlap,fs);
            [psd_bw_tilt_d(d,:) , psd_median_tilt_d(d,:) ] = asselbornSpectralFeatures(psd_tilt{d}, freq_tilt{d});
            else
                psd_bw_tilt_d(d,:)=0;
                psd_median_tilt_d(d,:)=0;

            end


            %%%%%%%%%%ACC dominant component%%%%%%%%%% si ordina lo spettro
            [~, Ind_max_A_d(d)] = max(A_tr_peak_d(d,:)); %ricavo indice della componente con potenza maggiore
            RPW_A_tr_dominant_PC_d(d) =  RPW_A_tr_PC_d(d,Ind_max_A_d(d));
            RPW_A_tr_dominant_PD_d(d) =  RPW_A_tr_PD_d(d,Ind_max_A_d(d));
            RPW_A_tr_dominant_Phy_d(d) =  RPW_A_tr_Phy_d(d,Ind_max_A_d(d));
            psd_bw_A_tr_dominant_d(d) = psd_bw_A_tr_d(d,Ind_max_A_d(d));
            psd_median_A_tr_dominant_d(d) = psd_median_A_tr_d(d,Ind_max_A_d(d));

            %%%%ACC dominant component new bands%%%%
            [~, Ind_max_A_d_new(d)] = max(A_tr_peak_d_new(d,:)); %ricavo indice della componente con potenza maggiore
            RPW_A_tr_dominant_2_d(d) =  RPW_2_A_tr_ab_d(d,Ind_max_A_d_new(d));
            RPW_A_tr_dominant_5_d(d) =  RPW_5_A_tr_ab_d(d,Ind_max_A_d_new(d));
            RPW_A_tr_dominant_8_d(d) =  RPW_8_A_tr_ab_d(d,Ind_max_A_d_new(d));
            RPW_A_tr_dominant_11_d(d) =  RPW_11_A_tr_ab_d(d,Ind_max_A_d_new(d));
            psd_bw_A_tr_dominant_d_new(d) = psd_bw_A_tr_d_new(d,Ind_max_A_d_new(d));
            psd_median_A_tr_dominant_d_new(d) = psd_median_A_tr_d_new(d,Ind_max_A_d_new(d));

            %%%%%%Sadikov acc dominant component%%%%%%%%%%
            [~, Ind_max_A_norm_d(d)] = max(A_tr_peak_norm_d(d,:)); %ricavo indice della componente con potenza maggiore
            psd_A_tr_norm_dominant_d(:,d) = psd_A_tr_norm_d{d}(:,Ind_max_A_norm_d(d));
            [ampXout_A_dom_d(d), out_lev_A_dom_d(d)] = peak_ampXout_from_PSD(psd_A_tr_norm_dominant_d(:,d),freq_A_tr{d},display_Sadikov);
            [dfA_dom_68(d), dfA_dom_50(d), fA_dom_sort_68(d), fA_dom_sort_50(d), fA_dom_sort_CM(d)] = harmonicity_from_PSD(psd_A_tr_norm_dominant_d(:,d),freq_A_tr{d},display_harmonicity);

            %%%%%%%%%%Gyro dominant component%%%%%%%%%%
            [max_G_dM(d), Ind_max_G_d(d)] = max(G_peak_dm(d,:)); %ricavo indice della componente con potenza maggiore
            RPW_G_tr_dominant_PC_d(d) =  RPW_G_PC_dm(d,Ind_max_G_d(d));
            RPW_G_tr_dominant_PD_d(d) =  RPW_G_PD_dm(d,Ind_max_G_d(d));
            RPW_G_tr_dominant_Phy_d(d) =  RPW_G_Phy_dm(d,Ind_max_G_d(d));
            psd_bw_G_tr_dominant_d(d) =psd_bw_G_dm(d,Ind_max_G_d(d));
            psd_median_G_tr_dominant_d(d) = psd_median_G_dm(d,Ind_max_G_d(d));

            %%Gyro dominant component new bands%%%%%%%%%%
            [max_G_dM_new(d),Ind_max_G_d_new(d)] = max(G_dm_peak_d_new(d,:)); %ricavo indice della componente con potenza maggiore
            RPW_G_tr_dominant_2_d(d) =  RPW_2_G_dm_ab_d(d,Ind_max_G_d_new(d));
            RPW_G_tr_dominant_5_d(d) =  RPW_5_G_dm_ab_d(d,Ind_max_G_d_new(d));
            RPW_G_tr_dominant_8_d(d) =  RPW_8_G_dm_ab_d(d,Ind_max_G_d_new(d));
            RPW_G_tr_dominant_11_d(d) =  RPW_11_G_dm_ab_d(d,Ind_max_G_d_new(d));
            psd_bw_G_tr_dominant_d_new(d) =psd_bw_G_dm_new(d,Ind_max_G_d_new(d));
            psd_median_G_tr_dominant_d_new(d) = psd_median_G_dm_new(d,Ind_max_G_d_new(d));

            %%%%%Sadikov gyro dominant component%%%%%%%%%%
            [~, Ind_max_G_norm_d(d)] = max(G_peak_norm_dm(d,:)); %ricavo indice della componente con potenza maggiore
            psd_G_tr_norm_dominant_d(:,d) = psd_G_norm_dm{d}(:,Ind_max_G_norm_d(d));
            [ampXout_G_dom_d(d), out_lev_G_dom_d(d)] = peak_ampXout_from_PSD(psd_G_tr_norm_dominant_d(:,d),freq_gm{d},display_Sadikov);
            [dfG_dom_68(d), dfG_dom_50(d), fG_dom_sort_68(d), fG_dom_sort_50(d), fG_dom_sort_CM(d)] = harmonicity_from_PSD(psd_G_tr_norm_dominant_d(:,d),freq_gm{d},display_harmonicity);

            psd_A_sum_d{d} = (psd_A_d{d}(:,1) + psd_A_d{d}(:,2) + psd_A_d{d}(:,3))/3;
            A_sum_peak_d(d) = max(psd_A_sum_d{d});
            A_sum_ind_peak_d(d) = find(psd_A_sum_d{d}== A_sum_peak_d(d));
            A_sum_peak_F_d(d) = freq_a{d}(A_sum_ind_peak_d(d));

            %         figure
            %         subplot 321; plot(freq_g{d}(GbandVol(1):GbandPhy(2)),psd_A_d{d}(GbandVol(1):GbandPhy(2),1)); title('Acc welch - components')
            %         subplot 323; plot(freq_g{d}(GbandVol(1):GbandPhy(2)),psd_A_d{d}(GbandVol(1):GbandPhy(2),2))
            %         subplot 325; plot(freq_g{d}(GbandVol(1):GbandPhy(2)),psd_A_d{d}(GbandVol(1):GbandPhy(2),3))
            %         subplot 322; plot(freq_g{d}(GbandVol(1):GbandPhy(2)),psd_G_d{d}(GbandVol(1):GbandPhy(2),1)); title('Gyro wpsd_bw_G3D_newelch - components')
            %         subplot 324; plot(freq_g{d}(GbandVol(1):GbandPhy(2)),psd_G_d{d}(GbandVol(1):GbandPhy(2),2))
            %         subplot 326; plot(freq_g{d}(GbandVol(1):GbandPhy(2)),psd_G_d{d}(GbandVol(1):GbandPhy(2),3)) s
            %         lfGyro=0.25;
            %         hfGyro=20;

            %velocitÃ  angolare (filtraggio in base ad articolo Heida, Wentink et Marani)
            %gyro_x{d,1}=Filtering(G_X, fs, lfGyro, hfGyro);
            %gyro_y{d,1}=Filtering(G_Y, fs, lfGyro, hfGyro);
            %gyro_z{d,1}=Filtering(G_Z, fs, lfGyro, hfGyro);

            %%%%%%%%%GYRO 3D%%%%%%%%%
            GYRO3D{d,1}=  Filtering(sqrt(G_X.^2  + G_Y.^2 + G_Z.^2), fs, 2, hfGyro);
            % GYRO3D{d,1} =DataCut{d,1}(:,21 + Ind_max_G_d(d)); %componente singola di tremore gyro con picco massimo
            [psd_G3D_d{d}, freq_g3D{d}, RPW_G3D_Vol_d(d), RPW_G3D_PC_d(d), RPW_G3D_PD_d(d), RPW_G3D_Phy_d(d), G3D_peak_d(d), G3D_peak_F_d(d)] = ...
                computeSpectralFeatures(GYRO3D{d,1}, segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, GbandVol(2), GbandVol(1));
            [psd_bw_G3D_d(d,:) , psd_median_G3D_d(d,:) ] = asselbornSpectralFeatures(psd_G3D_d{d}, freq_g3D{d});
            %%%%new band%%%%
            [psd_G3D_new{d}, freq_G3D_new{d}, RPW_2_G3D_ab_d(d), RPW_5_G3D_ab_d(d), RPW_8_G3D_ab_d(d), RPW_11_G3D_ab_d(d),RPW_2_G3D_rel_d(d),RPW_5_G3D_rel_d(d),RPW_8_G3D_rel_d(d),RPW_11_G3D_rel_d(d), G3D_peak_d_new(d), G3D_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(GYRO3D{d,1}, segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_G3D_d_new(d,:) , psd_median_G3D_d_new(d,:) ] = asselbornSpectralFeatures(psd_G3D_new{d}, freq_G3D_new{d});

            %%%%%%%%%ACC3D%%%%%%%%%
            Acc{d,1}= ACC3D;
            [psd_A3D_d{d}, freq_a3D{d}, RPW_A3D_Vol_d(d), RPW_A3D_PC_d(d), RPW_A3D_PD_d(d), RPW_A3D_Phy_d(d), A3D_peak_d(d), A3D_peak_F_d(d)] = ...
                computeSpectralFeatures(Acc{d,1}, segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt, 1);
            [psd_bw_A3D_d(d,:) , psd_median_A3D_d(d,:) ] = asselbornSpectralFeatures(psd_A3D_d{d}, freq_a3D{d});
            %Sadikov
            psd_A3D_norm_d{d} = pwelch(Acc{d,1}./max(Acc{d,1}), segment_length, overlap, segment_length, fs);
            [ampXout_A3D_norm_d(d), out_lev_A3D_norm_d(d)] = peak_ampXout_from_PSD(psd_A3D_norm_d{d},freq_a3D{d},display_Sadikov);
            [dfA3D_68(d), dfA3D_50(d), fA3D_sort_68(d), fA3D_sort_50(d), fA3D_sort_CM(d)] = harmonicity_from_PSD(psd_A3D_norm_d{d},freq_a3D{d},display_harmonicity);
            %%%new band%%%
            [psd_A3D_d_new{d}, freq_A3D_new{d}, RPW_2_A3D_ab_d(d), RPW_5_A3D_ab_d(d), RPW_8_A3D_ab_d(d), RPW_11_A3D_ab_d(d),RPW_2_A3D_rel_d(d),RPW_5_A3D_rel_d(d),RPW_8_A3D_rel_d(d),RPW_11_A3D_rel_d(d), A3D_peak_d_new(d), A3D_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(Acc{d,1}, segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_A3D_d_new(d,:) , psd_median_A3D_d_new(d,:) ] = asselbornSpectralFeatures(psd_A3D_d_new{d}, freq_A3D_new{d});
            %         figure
            %         plot(freq_A3D_new{d},psd_A3D_d_new{d})

            %%%%%%%%%PC%%%%%%%%%%%%%%
            %%%%%%%%%PC1 Gyro%%%%%%%%%
            [psd_G_PC1_d{d}, freq_G_PC1{d}, ~, RPW_G_PC1_PC_d(d), RPW_G_PC1_PD_d(d), RPW_G_PC1_Phy_d(d), G_PC1_peak_d(d), G_PC1_peak_F_d(d)] = ...
                computeSpectralFeatures(DataCut{d,1}(:,31), segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_G_PC1_d(d,:) , psd_median_G_PC1_d(d,:) ] = asselbornSpectralFeatures(psd_G_PC1_d{d}, freq_G_PC1{d});
            %Sadikov
            psd_G_PC1_norm_d{d} = pwelch(DataCut{d,1}(:,31)./max(DataCut{d,1}(:,31)), segment_length, overlap, segment_length, fs);
            [ampXout_G_PC1_norm_d(d), out_lev_G_PC1_norm_d(d)] = peak_ampXout_from_PSD(psd_G_PC1_norm_d{d},freq_G_PC1{d},display_Sadikov);
            [dfG_PC1_68(d), dfG_PC1_50(d), fG_PC1_sort_68(d), fG_PC1_sort_50(d), fG_PC1_sort_CM(d)] = harmonicity_from_PSD(psd_G_PC1_norm_d{d},freq_G_PC1{d},display_harmonicity);
            %%%%new band%%%
            [psd_G_PC1_new{d}, freq_G_PC1_new{d}, RPW_2_G_PC1_ab_d(d), RPW_5_G_PC1_ab_d(d), RPW_8_G_PC1_ab_d(d), RPW_11_G_PC1_ab_d(d),RPW_2_G_PC1_rel_d(d),RPW_5_G_PC1_rel_d(d),RPW_8_G_PC1_rel_d(d),RPW_11_G_PC1_rel_d(d), G_PC1_peak_d_new(d), G_PC1_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(DataCut{d,1}(:,31), segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_G_PC1_d_new(d,:) , psd_median_G_PC1_d_new(d,:) ] = asselbornSpectralFeatures(psd_G_PC1_new{d}, freq_G_PC1_new{d});

            %%%%%%%%%PC1+PC2 Gyro%%%%%%%%%
            [psd_G_PC_sum_d{d}, freq_G_PC_sum{d}, ~, RPW_G_PC_sum_PC_d(d), RPW_G_PC_sum_PD_d(d), RPW_G_PC_sum_Phy_d(d), G_PC_sum_peak_d(d), G_PC_sum_peak_F_d(d)] = ...
                computeSpectralFeatures(DataCut{d,1}(:,32), segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_G_PC_sum_d(d,:) , psd_median_G_PC_sum_d(d,:) ] = asselbornSpectralFeatures(psd_G_PC_sum_d{d}, freq_G_PC_sum{d});
            %Sadikov
            psd_G_PC_sum_norm_d{d} = pwelch(DataCut{d,1}(:,32)./max(DataCut{d,1}(:,32)), segment_length, overlap, segment_length, fs);
            [ampXout_G_PC_sum_norm_d(d), out_lev_G_PC_sum_norm_d(d)] = peak_ampXout_from_PSD(psd_G_PC_sum_norm_d{d},freq_G_PC1{d},display_Sadikov);
            [dfG_PC_sum_68(d), dfG_PC_sum_50(d), fG_PC_sum_sort_68(d), fG_PC_sum_sort_50(d), fG_PC_sum_sort_CM(d)] = harmonicity_from_PSD(psd_G_PC_sum_norm_d{d},freq_G_PC1{d},display_harmonicity);
            %%%new band%%%
            [psd_G_PC_sum_new{d}, freq_G_PC_sum_new{d}, RPW_2_G_PC_sum_ab_d(d), RPW_5_G_PC_sum_ab_d(d), RPW_8_G_PC_sum_ab_d(d), RPW_11_G_PC_sum_ab_d(d),RPW_2_G_PC_sum_rel_d(d),RPW_5_G_PC_sum_rel_d(d),RPW_8_G_PC_sum_rel_d(d),RPW_11_G_PC_sum_rel_d(d), G_PC_sum_peak_d_new(d), G_PC_sum_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(DataCut{d,1}(:,32), segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_G_PC_sum_d_new(d,:) , psd_median_G_PC_sum_d_new(d,:) ] = asselbornSpectralFeatures(psd_G_PC_sum_new{d}, freq_G_PC_sum_new{d});

            %%%%%%%%%PC Gyro%%%%%%%%%
            [psd_G_PC_d{d}, freq_G_PC{d}, ~, RPW_G_PC_PC_d(d), RPW_G_PC_PD_d(d), RPW_G_PC_Phy_d(d), G_PC_peak_d(d), G_PC_peak_F_d(d)] = ...
                computeSpectralFeatures(DataCut{d,1}(:,33), segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_G_PC_d(d,:) , psd_median_G_PC_d(d,:) ] = asselbornSpectralFeatures(psd_G_PC_d{d}, freq_G_PC{d});
            %Sadikov
            psd_G_PC_norm_d{d} = pwelch(DataCut{d,1}(:,33)./max(DataCut{d,1}(:,33)), segment_length, overlap, segment_length, fs);
            [ampXout_G_PC_norm_d(d), out_lev_G_PC_norm_d(d)] = peak_ampXout_from_PSD(psd_G_PC_norm_d{d},freq_G_PC1{d},display_Sadikov);
            [dfG_PC_68(d), dfG_PC_50(d), fG_PC_sort_68(d), fG_PC_sort_50(d), fG_PC_sort_CM(d)] = harmonicity_from_PSD(psd_G_PC_norm_d{d},freq_G_PC1{d},display_harmonicity);
            %%%new band%%%
            [psd_G_PC_new{d}, freq_G_PC_new{d}, RPW_2_G_PC_ab_d(d), RPW_5_G_PC_ab_d(d), RPW_8_G_PC_ab_d(d), RPW_11_G_PC_ab_d(d),RPW_2_G_PC_rel_d(d),RPW_5_G_PC_rel_d(d),RPW_8_G_PC_rel_d(d),RPW_11_G_PC_rel_d(d), G_PC_peak_d_new(d), G_PC_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(DataCut{d,1}(:,33), segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_G_PC_d_new(d,:) , psd_median_G_PC_d_new(d,:) ] = asselbornSpectralFeatures(psd_G_PC_new{d}, freq_G_PC_new{d});

            %%%%%%%%%PC1 Acc%%%%%%%%%
            [psd_A_PC1_d{d}, freq_A_PC1{d}, ~, RPW_A_PC1_PC_d(d), RPW_A_PC1_PD_d(d), RPW_A_PC1_Phy_d(d), A_PC1_peak_d(d), A_PC1_peak_F_d(d)] = ...
                computeSpectralFeatures(DataCut{d,1}(:,28), segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_A_PC1_d(d,:) , psd_median_A_PC1_d(d,:) ] = asselbornSpectralFeatures(psd_A_PC1_d{d}, freq_A_PC1{d});
            %Sadikov
            psd_A_PC1_norm_d{d} = pwelch(DataCut{d,1}(:,28)./max(DataCut{d,1}(:,28)), segment_length, overlap, segment_length, fs);
            [ampXout_A_PC1_norm_d(d), out_lev_A_PC1_norm_d(d)] = peak_ampXout_from_PSD(psd_A_PC1_norm_d{d},freq_A_PC1{d},display_Sadikov);
            [dfA_PC1_68(d), dfA_PC1_50(d), fA_PC1_sort_68(d), fA_PC1_sort_50(d), fA_PC1_sort_CM(d)] = harmonicity_from_PSD(psd_A_PC1_norm_d{d},freq_A_PC1{d},display_harmonicity);
            %%%new band%%%
            [psd_A_PC1_new{d}, freq_A_PC1_new{d}, RPW_2_A_PC1_ab_d(d), RPW_5_A_PC1_ab_d(d), RPW_8_A_PC1_ab_d(d), RPW_11_A_PC1_ab_d(d),RPW_2_A_PC1_rel_d(d),RPW_5_A_PC1_rel_d(d),RPW_8_A_PC1_rel_d(d),RPW_11_A_PC1_rel_d(d), A_PC1_peak_d_new(d), A_PC1_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(DataCut{d,1}(:,28), segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_A_PC1_d_new(d,:) , psd_median_A_PC1_d_new(d,:) ] = asselbornSpectralFeatures(psd_A_PC1_new{d}, freq_A_PC1_new{d});

            %%%%%%%%%PC_sum Acc%%%%%%%%%
            [psd_A_PC_sum_d{d}, freq_A_PC_sum{d}, ~, RPW_A_PC_sum_PC_d(d), RPW_A_PC_sum_PD_d(d), RPW_A_PC_sum_Phy_d(d), A_PC_sum_peak_d(d), A_PC_sum_peak_F_d(d)] = ...
                computeSpectralFeatures(DataCut{d,1}(:,29), segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_A_PC_sum_d(d,:) , psd_median_A_PC_sum_d(d,:) ] = asselbornSpectralFeatures(psd_A_PC_sum_d{d}, freq_A_PC_sum{d});
            %Sadikov
            psd_A_PC_sum_norm_d{d} = pwelch(DataCut{d,1}(:,29)./max(DataCut{d,1}(:,29)), segment_length, overlap, segment_length, fs);
            [ampXout_A_PC_sum_norm_d(d), out_lev_A_PC_sum_norm_d(d)] = peak_ampXout_from_PSD(psd_A_PC_sum_norm_d{d},freq_A_PC1{d},display_Sadikov);
            [dfA_PC_sum_68(d), dfA_PC_sum_50(d), fA_PC_sum_sort_68(d), fA_PC_sum_sort_50(d), fA_PC_sum_sort_CM(d)] = harmonicity_from_PSD(psd_A_PC_sum_norm_d{d},freq_A_PC1{d},display_harmonicity);
            %%%new band%%%
            [psd_A_PC_sum_new{d}, freq_A_PC_sum_new{d}, RPW_2_A_PC_sum_ab_d(d), RPW_5_A_PC_sum_ab_d(d), RPW_8_A_PC_sum_ab_d(d), RPW_11_A_PC_sum_ab_d(d),RPW_2_A_PC_sum_rel_d(d),RPW_5_A_PC_sum_rel_d(d),RPW_8_A_PC_sum_rel_d(d),RPW_11_A_PC_sum_rel_d(d), A_PC_sum_peak_d_new(d), A_PC_sum_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(DataCut{d,1}(:,29), segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_A_PC_sum_d_new(d,:) , psd_median_A_PC_sum_d_new(d,:) ] = asselbornSpectralFeatures(psd_A_PC_sum_new{d}, freq_A_PC_sum_new{d});

            %%%%%%%%%PC Acc%%%%%%%%%
            [psd_A_PC_d{d}, freq_A_PC{d}, ~, RPW_A_PC_PC_d(d), RPW_A_PC_PD_d(d), RPW_A_PC_Phy_d(d), A_PC_peak_d(d), A_PC_peak_F_d(d)] = ...
                computeSpectralFeatures(DataCut{d,1}(:,30), segment_length, overlap, fs, GbandVol, GbandPC, GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_A_PC_d(d,:) , psd_median_A_PC_d(d,:) ] = asselbornSpectralFeatures(psd_A_PC_d{d}, freq_A_PC{d});
            %Sadikov
            psd_A_PC_norm_d{d} = pwelch(DataCut{d,1}(:,30)./max(DataCut{d,1}(:,30)), segment_length, overlap, segment_length, fs);
            [ampXout_A_PC_norm_d(d), out_lev_A_PC_norm_d(d)] = peak_ampXout_from_PSD(psd_A_PC_norm_d{d},freq_A_PC1{d},display_Sadikov);
            [dfA_PC_68(d), dfA_PC_50(d), fA_PC_sort_68(d), fA_PC_sort_50(d), fA_PC_sort_CM(d)] = harmonicity_from_PSD(psd_A_PC_norm_d{d},freq_A_PC1{d},display_harmonicity);
            %%%new band%%%
            [psd_A_PC_new{d}, freq_A_PC_new{d}, RPW_2_A_PC_ab_d(d), RPW_5_A_PC_ab_d(d), RPW_8_A_PC_ab_d(d), RPW_11_A_PC_ab_d(d),RPW_2_A_PC_rel_d(d),RPW_5_A_PC_rel_d(d),RPW_8_A_PC_rel_d(d),RPW_11_A_PC_rel_d(d), A_PC_peak_d_new(d), A_PC_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(DataCut{d,1}(:,30), segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_A_PC_d_new(d,:) , psd_median_A_PC_d_new(d,:) ] = asselbornSpectralFeatures(psd_A_PC_new{d}, freq_A_PC_new{d});
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % manual tremor and voluntary contribution extraction
            vol_tot{d,1}=LPF(Acc{d,1}, 2); %0-2Hz
            tremor_tot{d,1}=HPF(Acc{d,1},2); %2-12Hz

            [psd_tr_d{d}, freq_tr{d}, ~, RPW_tr_PC_d(d), RPW_tr_PD_d(d), RPW_tr_Phy_d(d), tr_peak_d(d), tr_peak_F_d(d)] = ...
                computeSpectralFeatures(tremor_tot{d,1}, segment_length, overlap, fs, GbandVol, [bandt GbandPD(1)], GbandPD, GbandPhy, bandt, bandt);
            [psd_bw_tr_d(d,:) , psd_median_tr_d(d,:) ] = asselbornSpectralFeatures(psd_tr_d{d}, freq_tr{d});
            %%%new bands%%%
            [psd_tr_d_new{d}, freq_tr_new{d}, RPW_2_tr_ab_d(d), RPW_5_tr_ab_d(d), RPW_8_tr_ab_d(d), RPW_11_tr_ab_d(d),RPW_2_tr_rel_d(d),RPW_5_tr_rel_d(d),RPW_8_tr_rel_d(d),RPW_11_tr_rel_d(d), tr_peak_d_new(d), tr_peak_F_d_new(d)] = ...
                computeSpectralFeatures_new(tremor_tot{d,1}, segment_length, overlap, fs, band_2, band_5, band_8, band_11);
            [psd_bw_tr_d_new(d,:) , psd_median_tr_d_new(d,:) ] = asselbornSpectralFeatures(psd_tr_d_new{d}, freq_tr_new{d});

            %%%%%%%%%MISURA FRA%%%%%%%%%
            %assi singoli
            %RPW_Fra_new(d)= compute_RPW_fra(psd_G_dm{d},freq_gm{d},G_peak_F_dm(d,:),G_peak_dm(d,:));
            for k=1:3
                psd_Gmnorm{d}(:,k)= psd_G_dm{d}(:,k)./sum(psd_G_dm{d}(:,k));
            end
            df=mean(diff(freq_gm{d}));
            add2=ceil((1/df)/2);
            [M(d), F(d)] = mode(G_peak_F_dm(d,:));
            flag_moda=0;
            flag_comparison=0;
            if F(d) > 1 %se c'è la moda
                [~,indf_mode(d)]=min(abs(freq_gm{d}-M(d)));
                flag_moda=1;
                for k = 1:3
                    [~,indf(k)]=min(abs(freq_gm{d}-G_peak_F_dm(d,k)));
                    %se freq picco dell'asse è in un intorno della moda di 1Hz,
                    %uso questa frequenza
                    if abs(indf(k)-indf_mode(d))<=add2
                        RPW_Fra(d,k)=sum(psd_Gmnorm{d}(indf(k)-add2:indf(k)+add2,k));
                    else %altrimenti uso la frequenza della moda
                        RPW_Fra(d,k)=sum(psd_Gmnorm{d}(indf_mode(d)-add2:indf_mode(d)+add2,k));
                    end
                end
            else %se non c'è la moda
                %trovo frequenza del picco massimo
                [~, ind_peak_max(d)] = max(G_peak_dm(d,:));
                [~, ind_peak_F_max(d)] = min(abs(freq_gm{d} - G_peak_F_dm(d,ind_peak_max(d))));
                for k = 1:3
                    [~,indf(k)]=min(abs(freq_gm{d}-G_peak_F_dm(d,k)));
                    %se freq picco dell'asse è in un intorno della frequenza di picco massimo di 1Hz,
                    %uso questa frequenza
                    if abs(indf(k)-ind_peak_F_max(d))<=add2
                        RPW_Fra(d,k)=sum(psd_Gmnorm{d}(indf(k)-add2:indf(k)+add2,k));
                    else %altrimenti uso la frequenza del picco massimo
                        RPW_Fra(d,k)=sum(psd_Gmnorm{d}(ind_peak_F_max(d)-add2:ind_peak_F_max(d)+add2,k));
                    end
                end
            end

            %%%new bands
            for k=1:3
                psd_Gmnorm_new{d}(:,k)= psd_G_dm_new{d}(:,k)./sum(psd_G_dm_new{d}(:,k));
            end
            df=mean(diff(freq_G_dm_new{d}));
            add2_new=ceil((1/df)/2);
            [M(d), F(d)] = mode(G_dm_peak_F_d_new(d,:));
            flag_moda=0;
            flag_comparison=0;
            if F(d) > 1 %se c'è la moda
                [~,indf_mode(d)]=min(abs(freq_G_dm_new{d}-M(d)));
                flag_moda=1;
                for k = 1:3
                    [~,indf(k)]=min(abs(freq_G_dm_new{d}-G_dm_peak_F_d_new(d,k)));
                    %se freq picco dell'asse è in un intorno della moda di 1Hz,
                    %uso questa frequenza
                    if abs(indf(k)-indf_mode(d))<=add2_new
                        RPW_Fra_new(d,k)=sum(psd_Gmnorm_new{d}(indf(k)-add2_new:indf(k)+add2_new,k));
                    else %altrimenti uso la frequenza della moda
                        RPW_Fra_new(d,k)=sum(psd_Gmnorm_new{d}(indf_mode(d)-add2_new:indf_mode(d)+add2_new,k));
                    end
                end
            else %se non c'è la moda
                %trovo frequenza del picco massimo
                [~, ind_peak_max(d)] = max(G_dm_peak_F_d_new(d,:));
                [~, ind_peak_F_max(d)] = min(abs(freq_G_dm_new{d} - G_dm_peak_F_d_new(d,ind_peak_max(d))));
                for k = 1:3
                    [~,indf(k)]=min(abs(freq_G_dm_new{d}-G_dm_peak_F_d_new(d,k)));
                    %se freq picco dell'asse è in un intorno della frequenza di picco massimo di 1Hz,
                    %uso questa frequenza
                    if abs(indf(k)-ind_peak_F_max(d))<=add2_new
                        RPW_Fra_new(d,k)=sum(psd_Gmnorm_new{d}(indf(k)-add2_new:indf(k)+add2_new,k));
                    else %altrimenti uso la frequenza del picco massimo
                        RPW_Fra_new(d,k)=sum(psd_Gmnorm_new{d}(ind_peak_F_max(d)-add2_new:ind_peak_F_max(d)+add2_new,k));
                    end
                end
            end

            %Stessa indicatore ma per PC -> non devo calcolare la moda
            psd_G_PC_norm{d}= psd_G_PC_d{d}./sum(psd_G_PC_d{d});
            psd_G_PC1_norm{d}= psd_G_PC1_d{d}./sum(psd_G_PC1_d{d});
            psd_G_PC_sum_norm{d}= psd_G_PC_sum_d{d}./sum(psd_G_PC_sum_d{d});
            %%% new indicators%%%
            psd_G_PC_norm_new{d}= psd_G_PC_new{d}./sum(psd_G_PC_new{d});
            psd_G_PC1_norm_new{d}= psd_G_PC1_new{d}./sum(psd_G_PC1_new{d});
            psd_G_PC_sum_norm_new{d}= psd_G_PC_sum_new{d}./sum(psd_G_PC_sum_new{d});

            [~,indfPC1]=min(abs(freq_gm{d}-G_PC1_peak_F_d(d)));
            RPW_Fra_G_PC1(d)=sum(psd_G_PC1_norm{d}(indfPC1-add2:indfPC1+add2));
            %%% new indicators%%%
            [~,indfPC1_new]=min(abs(freq_G_dm_new{d}-G_PC1_peak_F_d_new(d)));
            RPW_Fra_G_PC1_new(d)=sum(psd_G_PC1_norm_new{d}(indfPC1_new-add2_new:indfPC1_new+add2_new));

            [~,indfPC]=min(abs(freq_gm{d}-G_PC_peak_F_d(d)));
            RPW_Fra_G_PC(d)=sum(psd_G_PC_norm{d}(indfPC-add2:indfPC+add2));
            %%% new indicators%%%
            [~,indfPC_new]=min(abs(freq_G_dm_new{d}-G_PC_peak_F_d(d)));
            RPW_Fra_G_PC_new(d)=sum(psd_G_PC_norm{d}(indfPC_new-add2_new:indfPC_new+add2_new));

            [~,indfPCsum]=min(abs(freq_gm{d}-G_PC_sum_peak_F_d(d)));
            RPW_Fra_G_PC_sum(d)=sum(psd_G_PC_sum_norm{d}(indfPCsum-add2:indfPCsum+add2));
            %%% new indicators%%%
            [~,indfPCsum_new]=min(abs(freq_G_dm_new{d}-G_PC_sum_peak_F_d(d)));
            RPW_Fra_G_PC_sum_new(d)=sum(psd_G_PC_sum_norm{d}(indfPCsum_new-add2_new:indfPCsum_new+add2_new));

        else %COSA FACCIO SE SEGMENTO SENZA PERDITE NON HA ALMENO 500 CAMPIONI
            TSI_DB(d)=0; TSI_DB_thr(d) =0; SNR(d) = 0; MHP(d) = 0;  TSI_Luft(d) =0;

            for k=1:3
                %manus
                RPWG_manus_s(d,k)=0; RMSG_manus_s(d,k)=0;m_manus_s(d,k)=0;

                RPW_A_Vol_d(d,k) = 0; RPW_A_PC_d(d,k) = 0;  RPW_A_PD_d(d,k) = 0; RPW_A_Phy_d(d,k) = 0;
                A_peak_d(d,k) = 0; A_peak_F_d(d,k) = NaN; psd_bw_A_d(d,k) = 0; psd_median_A_d(d,k) = 0;
                RPW_2_A_ab_d(d,k)=0;  RPW_5_A_ab_d(d,k)=0; RPW_8_A_ab_d(d,k)=0; RPW_11_A_ab_d(d,k)=0;RPW_2_A_rel_d(d,k)=0;RPW_5_A_rel_d(d,k)=0;RPW_8_A_rel_d(d,k)=0;RPW_11_A_rel_d(d,k)=0;
                A_peak_d_new(d,k)=0; A_peak_F_d_new(d,k) = NaN; psd_bw_A_d_new(d,k)=0; psd_median_A_d_new(d,k)=0;

                RPW_A_tr_PC_d(d,k) = 0; RPW_A_tr_PD_d(d,k) = 0; RPW_A_tr_Phy_d(d,k) = 0;
                A_tr_peak_d(d,k) = 0; A_tr_peak_F_d(d,k) = NaN; psd_bw_A_tr_d(d,k) = 0; psd_median_A_tr_d(d,k) = 0;
                RPW_2_A_tr_ab_d(d,k)=0; RPW_5_A_tr_ab_d(d,k)=0; RPW_8_A_tr_ab_d(d,k)=0; RPW_11_A_tr_ab_d(d,k)=0;RPW_2_A_tr_rel_d(d,k)=0;RPW_5_A_tr_rel_d(d,k)=0;RPW_8_A_tr_rel_d(d,k)=0;RPW_11_A_tr_rel_d(d,k)=0;
                A_tr_peak_d_new(d,k)=0; A_tr_peak_F_d_new(d,k)= NaN; psd_bw_A_tr_d_new(d,k) =0; psd_median_A_tr_d_new(d,k)=0;

                RPW_G_Vol_d(d,k) = 0; RPW_G_PC_d(d,k) = 0; RPW_G_PD_d(d,k) = 0; RPW_G_Phy_d(d,k) = 0;
                G_peak_d(d,k) = 0; G_peak_F_d(d,k) = NaN; psd_bw_G_d(d,k) = 0; psd_median_G_d(d,k) = 0;
                RPW_2_G_ab_d(d,k)=0; RPW_5_G_ab_d(d,k)=0; RPW_8_G_ab_d(d,k)=0; RPW_11_G_ab_d(d,k)=0;RPW_2_G_rel_d(d,k)=0;RPW_5_G_rel_d(d,k)=0;RPW_8_G_rel_d(d,k)=0;RPW_11_G_rel_d(d,k)=0;
                G_peak_d_new(d,k)=0; G_peak_F_d_new(d,k)=NaN; psd_bw_G_d_new(d,k)=0; psd_median_G_d_new(d,k)=0;

                RPW_G_PC_dm(d,k) = 0;  RPW_G_PD_dm(d,k) = 0; RPW_G_Phy_dm(d,k) = 0;
                G_peak_dm(d,k) = 0;  G_peak_F_dm(d,k) = NaN; psd_bw_G_dm(d,k) = 0; psd_median_G_dm(d,k) = 0;
                RPW_2_G_dm_ab_d(d,k)=0; RPW_5_G_dm_ab_d(d,k)=0; RPW_8_G_dm_ab_d(d,k)=0; RPW_11_G_dm_ab_d(d,k)=0;RPW_2_G_dm_rel_d(d,k)=0;RPW_5_G_dm_rel_d(d,k)=0;RPW_8_G_dm_rel_d(d,k)=0;RPW_11_G_dm_rel_d(d,k)=0;
                G_dm_peak_d_new(d,k)=0; G_dm_peak_F_d_new(d,k)=NaN; psd_bw_G_dm_new(d,k) =0; psd_median_G_dm_new(d,k)=0;

                psd_bw_tilt_d=0; psd_median_tilt_d=0;

                RPW_Fra(d,k) = NaN;
                RPW_Fra_new(d,k)= NaN;

            end

            RPW_A_tr_dominant_PC_d(d)= 0; RPW_A_tr_dominant_PD_d(d)= 0;  RPW_A_tr_dominant_Phy_d(d)= 0;
            psd_bw_A_tr_dominant_d(d) = 0; psd_median_A_tr_dominant_d(d) = 0;
            RPW_A_tr_dominant_2_d(d)=0; RPW_A_tr_dominant_5_d(d)=0; RPW_A_tr_dominant_8_d(d)=0; RPW_A_tr_dominant_11_d(d)=0;
            psd_bw_A_tr_dominant_d_new(d) = 0; psd_median_A_tr_dominant_d_new(d) = 0;


            RPW_G_tr_dominant_PC_d(d)= 0; RPW_G_tr_dominant_PD_d(d)= 0; RPW_G_tr_dominant_Phy_d(d)= 0;
            psd_bw_G_tr_dominant_d(d) = 0; psd_median_G_tr_dominant_d(d) = 0;
            RPW_G_tr_dominant_2_d(d)=0; RPW_G_tr_dominant_5_d(d)=0; RPW_G_tr_dominant_8_d(d)=0; RPW_G_tr_dominant_11_d(d)=0;
            psd_bw_G_tr_dominant_d_new(d) = 0; psd_median_G_tr_dominant_d_new(d) = 0;


            A_sum_peak_d(d)= 0;   A_sum_peak_F_d(d)= 0;

            RPW_A3D_Vol_d(d)= 0; RPW_A3D_PC_d(d)= 0;  RPW_A3D_PD_d(d)= 0; RPW_A3D_Phy_d(d)= 0;
            A3D_peak_d(d)= 0;   A3D_peak_F_d(d)= NaN; psd_bw_A3D_d(d) = 0; psd_median_A3D_d(d) = 0;
            RPW_2_A3D_ab_d(d)=0; RPW_5_A3D_ab_d(d)=0; RPW_8_A3D_ab_d(d)=0; RPW_11_A3D_ab_d(d)=0;RPW_2_A3D_rel_d(d)=0;RPW_5_A3D_rel_d(d)=0;RPW_8_A3D_rel_d(d)=0;RPW_11_A3D_rel_d(d)=0;
            A3D_peak_d_new(d)=0; A3D_peak_F_d_new(d)=NaN; psd_bw_A3D_d_new(d) =0; psd_median_A3D_d_new(d)=0;

            RPW_A_PC1_PC_d(d)= 0; RPW_A_PC1_PD_d(d)= 0; RPW_A_PC1_Phy_d(d)= 0;
            A_PC1_peak_d(d)= 0; A_PC1_peak_F_d(d)= 0; psd_bw_A_PC1_d(d) = 0; psd_median_A_PC1_d(d) = 0;
            RPW_2_A_PC1_ab_d(d)=0; RPW_5_A_PC1_ab_d(d)=0; RPW_8_A_PC1_ab_d(d)=0; RPW_11_A_PC1_ab_d(d)=0;RPW_2_A_PC1_rel_d(d)=0;RPW_5_A_PC1_rel_d(d)=0;RPW_8_A_PC1_rel_d(d)=0;RPW_11_A_PC1_rel_d(d)=0;
            A_PC1_peak_d_new(d)=0; A_PC1_peak_F_d_new(d) = NaN; psd_bw_A_PC1_d_new(d) =0; psd_median_A_PC1_d_new(d) = 0;

            RPW_A_PC_sum_PC_d(d)= 0; RPW_A_PC_sum_PD_d(d)= 0; RPW_A_PC_sum_Phy_d(d)= 0;
            A_PC_sum_peak_d(d)= 0; A_PC_sum_peak_F_d(d)= 0; psd_bw_A_PC_sum_d(d) = 0; psd_median_A_PC_sum_d(d) = 0;
            RPW_2_A_PC_sum_ab_d(d)=0; RPW_5_A_PC_sum_ab_d(d)=0; RPW_8_A_PC_sum_ab_d(d)=0; RPW_11_A_PC_sum_ab_d(d)=0;RPW_2_A_PC_sum_rel_d(d)=0;RPW_5_A_PC_sum_rel_d(d)=0;RPW_8_A_PC_sum_rel_d(d)=0;RPW_11_A_PC_sum_rel_d(d)=0;
            A_PC_sum_peak_d_new(d)=0; A_PC_sum_peak_F_d_new(d) = NaN; psd_bw_A_PC_sum_d_new(d) =0; psd_median_A_PC_sum_d_new(d) = 0;

            RPW_A_PC_PC_d(d)= 0; RPW_A_PC_PD_d(d)= 0;   RPW_A_PC_Phy_d(d)= 0;
            A_PC_peak_d(d)= 0; A_PC_peak_F_d(d)= 0; psd_bw_A_PC_d(d) = 0; psd_median_A_PC_d(d) = 0;
            RPW_2_A_PC_ab_d(d)=0; RPW_5_A_PC_ab_d(d)=0; RPW_8_A_PC_ab_d(d)=0; RPW_11_A_PC_ab_d(d)=0;RPW_2_A_PC_rel_d(d)=0;RPW_5_A_PC_rel_d(d)=0;RPW_8_A_PC_rel_d(d)=0;RPW_11_A_PC_rel_d(d)=0;
            A_PC_peak_d_new(d)=0; A_PC_peak_F_d_new(d) = NaN; psd_bw_A_PC_d_new(d) =0; psd_median_A_PC_d_new(d) = 0;

            RPW_G_PC1_PC_d(d)= 0; RPW_G_PC1_PD_d(d)= 0; RPW_G_PC1_Phy_d(d)= 0;
            G_PC1_peak_d(d)= 0; G_PC1_peak_F_d(d)= 0; psd_bw_G_PC1_d(d) = 0; psd_median_G_PC1_d(d) = 0;
            RPW_2_G_PC1_ab_d(d)=0; RPW_5_G_PC1_ab_d(d)=0; RPW_8_G_PC1_ab_d(d)=0; RPW_11_G_PC1_ab_d(d)=0;RPW_2_G_PC1_rel_d(d)=0;RPW_5_G_PC1_rel_d(d)=0;RPW_8_G_PC1_rel_d(d)=0;RPW_11_G_PC1_rel_d(d)=0;
            G_PC1_peak_d_new(d)=0; G_PC1_peak_F_d_new(d) = NaN; psd_bw_G_PC1_d_new(d) =0; psd_median_G_PC1_d_new(d) = 0;

            RPW_G_PC_sum_PC_d(d)= 0; RPW_G_PC_sum_PD_d(d)= 0;  RPW_G_PC_sum_Phy_d(d)= 0;
            G_PC_sum_peak_d(d)= 0; G_PC_sum_peak_F_d(d)= 0; psd_bw_G_PC_sum_d(d) = 0; psd_median_G_PC_sum_d(d) = 0;
            RPW_2_G_PC_sum_ab_d(d)=0; RPW_5_G_PC_sum_ab_d(d)=0; RPW_8_G_PC_sum_ab_d(d)=0; RPW_11_G_PC_sum_ab_d(d)=0;RPW_2_G_PC_sum_rel_d(d)=0;RPW_5_G_PC_sum_rel_d(d)=0;RPW_8_G_PC_sum_rel_d(d)=0;RPW_11_G_PC_sum_rel_d(d)=0;
            G_PC_sum_peak_d_new(d)=0; G_PC_sum_peak_F_d_new(d) = NaN; psd_bw_G_PC_sum_d_new(d) =0; psd_median_G_PC_sum_d_new(d) = 0;

            RPW_G_PC_PC_d(d)= 0; RPW_G_PC_PD_d(d)= 0; RPW_G_PC_Phy_d(d)= 0;
            G_PC_peak_d(d)= 0;  G_PC_peak_F_d(d)= 0;  psd_bw_G_PC_d(d) = 0; psd_median_G_PC_d(d) = 0;
            RPW_2_G_PC_ab_d(d)=0; RPW_5_G_PC_ab_d(d)=0; RPW_8_G_PC_ab_d(d)=0; RPW_11_G_PC_ab_d(d)=0;RPW_2_G_PC_rel_d(d)=0;RPW_5_G_PC_rel_d(d)=0;RPW_8_G_PC_rel_d(d)=0;RPW_11_G_PC_rel_d(d)=0;
            G_PC_peak_d_new(d)=0; G_PC_peak_F_d_new(d) = NaN; psd_bw_G_PC_d_new(d) =0; psd_median_G_PC_d_new(d) = 0;

            RPW_tr_PC_d(d)= 0; RPW_tr_PD_d(d)= 0;  RPW_tr_Phy_d(d)= 0;
            tr_peak_d(d)= 0; tr_peak_F_d(d)= 0; psd_bw_tr_d(d) = 0; psd_median_tr_d(d) = 0;
            RPW_2_tr_ab_d(d)=0; RPW_5_tr_ab_d(d)=0; RPW_8_tr_ab_d(d)=0; RPW_11_tr_ab_d(d)=0;RPW_2_tr_rel_d(d)=0;RPW_5_tr_rel_d(d)=0;RPW_8_tr_rel_d(d)=0;RPW_11_tr_rel_d(d)=0;
            tr_peak_d_new(d)=0; tr_peak_F_d_new(d) = NaN; psd_bw_tr_d_new(d)=0; psd_median_tr_d_new(d)= 0;

            RPW_G3D_Vol_d(d)= 0; RPW_G3D_PC_d(d)= 0; RPW_G3D_PD_d(d)= 0; RPW_G3D_Phy_d(d)= 0;
            G3D_peak_d(d)= 0; G3D_peak_F_d(d)= NaN; psd_bw_G3D_d(d) = 0; psd_median_G3D_d(d) = 0;
            RPW_2_G3D_ab_d(d)=0; RPW_5_G3D_ab_d(d)=0; RPW_8_G3D_ab_d(d)=0; RPW_11_G3D_ab_d(d)=0;RPW_2_G3D_rel_d(d)=0;RPW_5_G3D_rel_d(d)=0;RPW_8_G3D_rel_d(d)=0;RPW_11_G3D_rel_d(d)=0;
            G3D_peak_d_new(d)=0; G3D_peak_F_d_new(d)= NaN; psd_bw_G3D_d_new(d) =0; psd_median_G3D_d_new(d)= 0;


            %Sadikov
            ampXout_A_dom_d(d)= 0; out_lev_A_dom_d(d)= 0; dfA_dom_68(d)=0;
            dfA_dom_50(d)=NaN;  fA_dom_sort_68(d)=0; fA_dom_sort_50(d)=0;  fA_dom_sort_CM(d)=0;

            ampXout_A3D_norm_d(d)= 0; out_lev_A3D_norm_d(d)= 0;
            dfA3D_68(d)=0; dfA3D_50(d)=NaN;  fA3D_sort_68(d)=0; fA3D_sort_50(d)=0; fA3D_sort_CM(d)=0;

            ampXout_A_PC1_norm_d(d)= 0; out_lev_A_PC1_norm_d(d)= 0;
            dfA_PC1_68(d)=0; dfA_PC1_50(d)=NaN; fA_PC1_sort_68(d)=0; fA_PC1_sort_50(d)=0; fA_PC1_sort_CM(d)=0;

            ampXout_A_PC_sum_norm_d(d)= 0; out_lev_A_PC_sum_norm_d(d)= 0;
            dfA_PC_sum_68(d)=0; dfA_PC_sum_50(d)=NaN; fA_PC_sum_sort_68(d)=0;  fA_PC_sum_sort_50(d)=0;  fA_PC_sum_sort_CM(d)=0;

            ampXout_A_PC_norm_d(d)= 0; out_lev_A_PC_norm_d(d)= 0;
            dfA_PC_68(d)=0; dfA_PC_50(d)=NaN; fA_PC_sort_68(d)=0; fA_PC_sort_50(d)=0; fA_PC_sort_CM(d)=0;

            ampXout_G_dom_d(d)= 0; out_lev_G_dom_d(d)= 0;
            dfG_dom_68(d)=0;  dfG_dom_50(d)=NaN; fG_dom_sort_68(d)=0; fG_dom_sort_50(d)=0; fG_dom_sort_CM(d)=0;

            ampXout_G_PC1_norm_d(d)= 0; out_lev_G_PC1_norm_d(d)= 0;
            dfG_PC1_68(d)=0; dfG_PC1_50(d)=NaN; fG_PC1_sort_68(d)=0; fG_PC1_sort_50(d)=0; fG_PC1_sort_CM(d)=0;

            ampXout_G_PC_sum_norm_d(d)= 0; out_lev_G_PC_sum_norm_d(d)= 0;
            dfG_PC_sum_68(d)=0; dfG_PC_sum_50(d)=NaN; fG_PC_sum_sort_68(d)=0; fG_PC_sum_sort_50(d)=0; fG_PC_sum_sort_CM(d)=0;

            ampXout_G_PC_norm_d(d)= 0; out_lev_G_PC_norm_d(d)= 0;
            dfG_PC_68(d)=0; dfG_PC_50(d)=NaN; fG_PC_sort_68(d)=0; fG_PC_sort_50(d)=0; fG_PC_sort_CM(d)=0;

            RPW_Fra_G_PC1(d)= 0; RPW_Fra_G_PC(d)= 0; RPW_Fra_G_PC_sum(d)= 0;

        end
        %%%%%%%%%%%%%%CALCOLO FEATURE TEMPORALI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if coppie(d) <= 1 %se ho meno di due stroke non ho nessun inair
            inair{d,1} =0;
            inair_duration(d,1)=0;
        end

        %meanOnsheet(d)=mean(onsheet(d,:));
        meanOnsheet(d)=mean(nonzeros(onsheet(d,:)));
        %controllo che dev std sia diversa da 0
        if meanOnsheet(d) > 0
            Onsheet_CV(d)=std(nonzeros(onsheet(d,:)))/meanOnsheet(d);
        else
            Onsheet_CV(d)=0;
        end

        clear i;

        %controllo durata media in air moments
        if mean(nonzeros(inair_duration(d,:))) > 0 %importante visto che potrebbero non esserci in air moments
            %durata complessiva in air
            %           meanInair(d)=mean(inair_duration(d,:));
            meanInair(d)=mean(nonzeros(inair_duration(d,:)));
            %           Inair_CV(d)=std(inair_duration(d,:))/meanInair(d);
            Inair_CV(d)=std(nonzeros(inair_duration(d,:)))/meanInair(d);

            %feature relative agli in air che non sono pause (in air < threshpause)
            ind_NOpause{d}=find(inair_duration(d,:) < threshpause);

            if   isnan(mean(inair_duration(d,ind_NOpause{d})))
                meanInair_nopause(d)=0;
                Inair_nopause_CV(d)=0;
            else
                meanInair_nopause(d)=mean(nonzeros(inair_duration(d,ind_NOpause{d})));
                Inair_nopause_CV(d)=std(nonzeros(inair_duration(d,ind_NOpause{d})))/meanInair_nopause(d);
            end

            %feature relative alle pause (in air >= threshpause)
            ind_pause{d}=find(inair_duration(d,:) >= threshpause);
            if   isnan(mean(inair_duration(d,ind_pause{d})))

                meanpause_duration(d)=0;
                pause_duration_CV(d)=0;
                pause_num(d)=0;
            else
                meanpause_duration(d)=mean(inair_duration(d,ind_pause{d}));
                pause_duration_CV(d)=std(inair_duration(d,ind_pause{d}))/meanpause_duration(d);
                pause_num(d)=length(inair_duration(d,ind_pause{d}));
            end

        else %se non ho in air moments nel tratto senza perdite
            meanInair(d)=0;
            Inair_CV(d)=0;
            meanInair_nopause(d)=0;
            Inair_nopause_CV(d)=0;
            meanpause_duration(d)=0;
            pause_duration_CV(d)=0;
            pause_num(d)=0;
        end
        %end
        %AirSheet_Ratio=NaN;
        %if inair~=0
        AirSheet_Ratio_Nopause(d)=meanInair_nopause(d)/meanOnsheet(d);
        AirSheet_Ratio_pause(d)=meanInair(d)/meanOnsheet(d);
        %end

    else %se non ho stroke nel tratto senza perdite

        meanP_DC(d) = 0; meanP_CV_DC(d) = 0; meanPOVS_DC(d) = 0;
        meanNCP_DC(d)=0; meanNCG_DC(d)=0;  meanNCA_DC(d)=0;  meanNCV_DC(d)=0;
        meanCPDP_AVG(d)=0;   meanCPDP_CV(d)=0; meanCPDG_AVG(d)=0; meanCPDG_CV(d)=0; meanCPDV_AVG(d)=0; meanCPDV_CV(d)=0;
        meanV_DC(d)=0; V_CV_AVG(d)=0; V_ApEn_AVG(d)=0;

        for k=1:3
            LDLJ_A_k(d,k) =0; LDLJ_G_k(d,k) =0;

            S1_sm(d,k) = 0; S2_sm(d,k) = 0; S3_sm(d,k) = 0; S4_sm(d,k) = 0; S45_sm(d,k) = 0; S5_sm(d,k) = 0;
            S1_k(d,k) =0; S2_k(d,k) =0; S3_k(d,k) =0; S4_k(d,k) =0; S45_k(d,k) =0; S5_k(d,k) =0;
        end

        SNvpd_norm_AVG(d)=0; medianDLJ_A(d)=0; medianSquaredJerk_A(d)=0;
        meanTILT(d)=0; meanTILT_CV(d)=0; meanTILT_VAR(d)=0;

        for k=1:3
            %manus
            RPWG_manus_s(d,k)=0; RMSG_manus_s(d,k)=0;m_manus_s(d,k)=0;

            RPW_A_Vol_d(d,k) = 0; RPW_A_PC_d(d,k) = 0; RPW_A_PD_d(d,k) = 0; RPW_A_Phy_d(d,k) = 0;
            A_peak_d(d,k) = 0; A_peak_F_d(d,k) = NaN;  psd_bw_A_d(d,k) = 0; psd_median_A_d(d,k) = 0;
            RPW_2_A_ab_d(d,k)=0;  RPW_5_A_ab_d(d,k)=0; RPW_8_A_ab_d(d,k)=0; RPW_11_A_ab_d(d,k)=0;RPW_2_A_rel_d(d,k)=0;RPW_5_A_rel_d(d,k)=0;RPW_8_A_rel_d(d,k)=0;RPW_11_A_rel_d(d,k)=0;
            A_peak_d_new(d,k)=0; A_peak_F_d_new(d,k) = NaN; psd_bw_A_d_new(d,k)=0; psd_median_A_d_new(d,k)=0;

            RPW_A_tr_PC_d(d,k) = 0; RPW_A_tr_PD_d(d,k) = 0; RPW_A_tr_Phy_d(d,k) = 0;
            A_tr_peak_d(d,k) = 0; A_tr_peak_F_d(d,k) = NaN; psd_bw_A_tr_d(d,k) = 0; psd_median_A_tr_d(d,k) = 0;
            RPW_2_A_tr_ab_d(d,k)=0; RPW_5_A_tr_ab_d(d,k)=0; RPW_8_A_tr_ab_d(d,k)=0; RPW_11_A_tr_ab_d(d,k)=0;RPW_2_A_tr_rel_d(d,k)=0;RPW_5_A_tr_rel_d(d,k)=0;RPW_8_A_tr_rel_d(d,k)=0;RPW_11_A_tr_rel_d(d,k)=0;
            A_tr_peak_d_new(d,k)=0; A_tr_peak_F_d_new(d,k)= NaN; psd_bw_A_tr_d_new(d,k) =0; psd_median_A_tr_d_new(d,k)=0;

            RPW_G_Vol_d(d,k) = 0; RPW_G_PC_d(d,k) = 0; RPW_G_PD_d(d,k) = 0; RPW_G_Phy_d(d,k) = 0;
            G_peak_d(d,k) = 0; G_peak_F_d(d,k) = NaN; psd_bw_G_d(d,k) = 0; psd_median_G_d(d,k) = 0;
            RPW_2_G_ab_d(d,k)=0; RPW_5_G_ab_d(d,k)=0; RPW_8_G_ab_d(d,k)=0; RPW_11_G_ab_d(d,k)=0;RPW_2_G_rel_d(d,k)=0;RPW_5_G_rel_d(d,k)=0;RPW_8_G_rel_d(d,k)=0;RPW_11_G_rel_d(d,k)=0;
            G_peak_d_new(d,k)=0; G_peak_F_d_new(d,k)=NaN; psd_bw_G_d_new(d,k)=0; psd_median_G_d_new(d,k)=0;

            RPW_G_PC_dm(d,k) = 0; RPW_G_PD_dm(d,k) = 0; RPW_G_Phy_dm(d,k) = 0;
            G_peak_dm(d,k) = 0; G_peak_F_dm(d,k) = NaN; psd_bw_G_dm(d,k) = 0; psd_median_G_dm(d,k) = 0;
            RPW_2_G_dm_ab_d(d,k)=0; RPW_5_G_dm_ab_d(d,k)=0; RPW_8_G_dm_ab_d(d,k)=0; RPW_11_G_dm_ab_d(d,k)=0;RPW_2_G_dm_rel_d(d,k)=0;RPW_5_G_dm_rel_d(d,k)=0;RPW_8_G_dm_rel_d(d,k)=0;RPW_11_G_dm_rel_d(d,k)=0;
            G_dm_peak_d_new(d,k)=0; G_dm_peak_F_d_new(d,k)=NaN; psd_bw_G_dm_new(d,k) =0; psd_median_G_dm_new(d,k)=0;

            psd_bw_tilt_d=0; psd_median_tilt_d=0;

            RPW_Fra(d,k) = NaN;
            RPW_Fra_new(d,k)= NaN;
        end

        RPW_A_tr_dominant_PC_d(d)= 0; RPW_A_tr_dominant_PD_d(d)= 0;  RPW_A_tr_dominant_Phy_d(d)= 0;
        psd_bw_A_tr_dominant_d(d) = 0; psd_median_A_tr_dominant_d(d) = 0;
        RPW_A_tr_dominant_2_d(d)=0; RPW_A_tr_dominant_5_d(d)=0; RPW_A_tr_dominant_8_d(d)=0; RPW_A_tr_dominant_11_d(d)=0;
        psd_bw_A_tr_dominant_d_new(d) = 0; psd_median_A_tr_dominant_d_new(d) = 0;


        RPW_G_tr_dominant_PC_d(d)= 0; RPW_G_tr_dominant_PD_d(d)= 0; RPW_G_tr_dominant_Phy_d(d)= 0;
        psd_bw_G_tr_dominant_d(d) = 0; psd_median_G_tr_dominant_d(d) = 0;
        RPW_G_tr_dominant_2_d(d)=0; RPW_G_tr_dominant_5_d(d)=0; RPW_G_tr_dominant_8_d(d)=0; RPW_G_tr_dominant_11_d(d)=0;
        psd_bw_G_tr_dominant_d_new(d) = 0; psd_median_G_tr_dominant_d_new(d) = 0;


        A_sum_peak_d(d)= 0;   A_sum_peak_F_d(d)= 0;

        RPW_A3D_Vol_d(d)= 0; RPW_A3D_PC_d(d)= 0;  RPW_A3D_PD_d(d)= 0; RPW_A3D_Phy_d(d)= 0;
        A3D_peak_d(d)= 0;   A3D_peak_F_d(d)= NaN; psd_bw_A3D_d(d) = 0; psd_median_A3D_d(d) = 0;
        RPW_2_A3D_ab_d(d)=0; RPW_5_A3D_ab_d(d)=0; RPW_8_A3D_ab_d(d)=0; RPW_11_A3D_ab_d(d)=0;RPW_2_A3D_rel_d(d)=0;RPW_5_A3D_rel_d(d)=0;RPW_8_A3D_rel_d(d)=0;RPW_11_A3D_rel_d(d)=0;
        A3D_peak_d_new(d)=0; A3D_peak_F_d_new(d)=NaN; psd_bw_A3D_d_new(d) =0; psd_median_A3D_d_new(d)=0;

        RPW_A_PC1_PC_d(d)= 0; RPW_A_PC1_PD_d(d)= 0; RPW_A_PC1_Phy_d(d)= 0;
        A_PC1_peak_d(d)= 0; A_PC1_peak_F_d(d)= 0; psd_bw_A_PC1_d(d) = 0; psd_median_A_PC1_d(d) = 0;
        RPW_2_A_PC1_ab_d(d)=0; RPW_5_A_PC1_ab_d(d)=0; RPW_8_A_PC1_ab_d(d)=0; RPW_11_A_PC1_ab_d(d)=0;RPW_2_A_PC1_rel_d(d)=0;RPW_5_A_PC1_rel_d(d)=0;RPW_8_A_PC1_rel_d(d)=0;RPW_11_A_PC1_rel_d(d)=0;
        A_PC1_peak_d_new(d)=0; A_PC1_peak_F_d_new(d) = NaN; psd_bw_A_PC1_d_new(d) =0; psd_median_A_PC1_d_new(d) = 0;

        RPW_A_PC_sum_PC_d(d)= 0; RPW_A_PC_sum_PD_d(d)= 0; RPW_A_PC_sum_Phy_d(d)= 0;
        A_PC_sum_peak_d(d)= 0; A_PC_sum_peak_F_d(d)= 0; psd_bw_A_PC_sum_d(d) = 0; psd_median_A_PC_sum_d(d) = 0;
        RPW_2_A_PC_sum_ab_d(d)=0; RPW_5_A_PC_sum_ab_d(d)=0; RPW_8_A_PC_sum_ab_d(d)=0; RPW_11_A_PC_sum_ab_d(d)=0;RPW_2_A_PC_sum_rel_d(d)=0;RPW_5_A_PC_sum_rel_d(d)=0;RPW_8_A_PC_sum_rel_d(d)=0;RPW_11_A_PC_sum_rel_d(d)=0;
        A_PC_sum_peak_d_new(d)=0; A_PC_sum_peak_F_d_new(d) = NaN; psd_bw_A_PC_sum_d_new(d) =0; psd_median_A_PC_sum_d_new(d) = 0;

        RPW_A_PC_PC_d(d)= 0; RPW_A_PC_PD_d(d)= 0;   RPW_A_PC_Phy_d(d)= 0;
        A_PC_peak_d(d)= 0; A_PC_peak_F_d(d)= 0; psd_bw_A_PC_d(d) = 0; psd_median_A_PC_d(d) = 0;
        RPW_2_A_PC_ab_d(d)=0; RPW_5_A_PC_ab_d(d)=0; RPW_8_A_PC_ab_d(d)=0; RPW_11_A_PC_ab_d(d)=0;RPW_2_A_PC_rel_d(d)=0;RPW_5_A_PC_rel_d(d)=0;RPW_8_A_PC_rel_d(d)=0;RPW_11_A_PC_rel_d(d)=0;
        A_PC_peak_d_new(d)=0; A_PC_peak_F_d_new(d) = NaN; psd_bw_A_PC_d_new(d) =0; psd_median_A_PC_d_new(d) = 0;

        RPW_G_PC1_PC_d(d)= 0; RPW_G_PC1_PD_d(d)= 0; RPW_G_PC1_Phy_d(d)= 0;
        G_PC1_peak_d(d)= 0; G_PC1_peak_F_d(d)= 0; psd_bw_G_PC1_d(d) = 0; psd_median_G_PC1_d(d) = 0;
        RPW_2_G_PC1_ab_d(d)=0; RPW_5_G_PC1_ab_d(d)=0; RPW_8_G_PC1_ab_d(d)=0; RPW_11_G_PC1_ab_d(d)=0;RPW_2_G_PC1_rel_d(d)=0;RPW_5_G_PC1_rel_d(d)=0;RPW_8_G_PC1_rel_d(d)=0;RPW_11_G_PC1_rel_d(d)=0;
        G_PC1_peak_d_new(d)=0; G_PC1_peak_F_d_new(d) = NaN; psd_bw_G_PC1_d_new(d) =0; psd_median_G_PC1_d_new(d) = 0;

        RPW_G_PC_sum_PC_d(d)= 0; RPW_G_PC_sum_PD_d(d)= 0;  RPW_G_PC_sum_Phy_d(d)= 0;
        G_PC_sum_peak_d(d)= 0; G_PC_sum_peak_F_d(d)= 0; psd_bw_G_PC_sum_d(d) = 0; psd_median_G_PC_sum_d(d) = 0;
        RPW_2_G_PC_sum_ab_d(d)=0; RPW_5_G_PC_sum_ab_d(d)=0; RPW_8_G_PC_sum_ab_d(d)=0; RPW_11_G_PC_sum_ab_d(d)=0;RPW_2_G_PC_sum_rel_d(d)=0;RPW_5_G_PC_sum_rel_d(d)=0;RPW_8_G_PC_sum_rel_d(d)=0;RPW_11_G_PC_sum_rel_d(d)=0;
        G_PC_sum_peak_d_new(d)=0; G_PC_sum_peak_F_d_new(d) = NaN; psd_bw_G_PC_sum_d_new(d) =0; psd_median_G_PC_sum_d_new(d) = 0;

        RPW_G_PC_PC_d(d)= 0; RPW_G_PC_PD_d(d)= 0; RPW_G_PC_Phy_d(d)= 0;
        G_PC_peak_d(d)= 0;  G_PC_peak_F_d(d)= 0;  psd_bw_G_PC_d(d) = 0; psd_median_G_PC_d(d) = 0;
        RPW_2_G_PC_ab_d(d)=0; RPW_5_G_PC_ab_d(d)=0; RPW_8_G_PC_ab_d(d)=0; RPW_11_G_PC_ab_d(d)=0;RPW_2_G_PC_rel_d(d)=0;RPW_5_G_PC_rel_d(d)=0;RPW_8_G_PC_rel_d(d)=0;RPW_11_G_PC_rel_d(d)=0;
        G_PC_peak_d_new(d)=0; G_PC_peak_F_d_new(d) = NaN; psd_bw_G_PC_d_new(d) =0; psd_median_G_PC_d_new(d) = 0;

        RPW_tr_PC_d(d)= 0; RPW_tr_PD_d(d)= 0;  RPW_tr_Phy_d(d)= 0;
        tr_peak_d(d)= 0; tr_peak_F_d(d)= 0; psd_bw_tr_d(d) = 0; psd_median_tr_d(d) = 0;
        RPW_2_tr_ab_d(d)=0; RPW_5_tr_ab_d(d)=0; RPW_8_tr_ab_d(d)=0; RPW_11_tr_ab_d(d)=0;RPW_2_tr_rel_d(d)=0;RPW_5_tr_rel_d(d)=0;RPW_8_tr_rel_d(d)=0;RPW_11_tr_rel_d(d)=0;
        tr_peak_d_new(d)=0; tr_peak_F_d_new(d) = NaN; psd_bw_tr_d_new(d)=0; psd_median_tr_d_new(d)= 0;

        RPW_G3D_Vol_d(d)= 0; RPW_G3D_PC_d(d)= 0; RPW_G3D_PD_d(d)= 0; RPW_G3D_Phy_d(d)= 0;
        G3D_peak_d(d)= 0; G3D_peak_F_d(d)= NaN; psd_bw_G3D_d(d) = 0; psd_median_G3D_d(d) = 0;
        RPW_2_G3D_ab_d(d)=0; RPW_5_G3D_ab_d(d)=0; RPW_8_G3D_ab_d(d)=0; RPW_11_G3D_ab_d(d)=0;RPW_2_G3D_rel_d(d)=0;RPW_5_GG3D_rel_d(d)=0;RPW_8_G3D_rel_d(d)=0;RPW_11_G3D_rel_d(d)=0;
        G3D_peak_d_new(d)=0; G3D_peak_F_d_new(d)= NaN; psd_bw_G3D_d_new(d) =0; psd_median_G3D_d_new(d)= 0;

        %Sadikov
        ampXout_A_dom_d(d) = 0; out_lev_A_dom_d(d) = 0;
        dfA_dom_68(d)=0; dfA_dom_50(d)=NaN; fA_dom_sort_68(d)=0; fA_dom_sort_50(d)=0; fA_dom_sort_CM(d)=0;

        ampXout_A3D_norm_d(d) = 0; out_lev_A3D_norm_d(d) = 0;
        dfA3D_68(d)=0; dfA3D_50(d)=NaN; fA3D_sort_68(d)=0; fA3D_sort_50(d)=0; fA3D_sort_CM(d)=0;

        ampXout_A_PC1_norm_d(d) = 0; out_lev_A_PC1_norm_d(d) = 0;
        dfA_PC1_68(d)=0; dfA_PC1_50(d)=NaN; fA_PC1_sort_68(d)=0; fA_PC1_sort_50(d)=0; fA_PC1_sort_CM(d)=0;

        ampXout_A_PC_sum_norm_d(d) = 0; out_lev_A_PC_sum_norm_d(d) = 0;
        dfA_PC_sum_68(d)=0; dfA_PC_sum_50(d)=NaN; fA_PC_sum_sort_68(d)=0; fA_PC_sum_sort_50(d)=0; fA_PC_sum_sort_CM(d)=0;

        ampXout_A_PC_norm_d(d) = 0; out_lev_A_PC_norm_d(d) = 0;
        dfA_PC_68(d)=0; dfA_PC_50(d)=NaN; fA_PC_sort_68(d)=0; fA_PC_sort_50(d)=0; fA_PC_sort_CM(d)=0;

        ampXout_G_dom_d(d) = 0; out_lev_G_dom_d(d) = 0;
        dfG_dom_68(d)=0; dfG_dom_50(d)=NaN; fG_dom_sort_68(d)=0; fG_dom_sort_50(d)=0; fG_dom_sort_CM(d)=0;

        ampXout_G_PC1_norm_d(d) = 0; out_lev_G_PC1_norm_d(d) = 0;
        dfG_PC1_68(d)=0; dfG_PC1_50(d)=NaN; fG_PC1_sort_68(d)=0; fG_PC1_sort_50(d)=0; fG_PC1_sort_CM(d)=0;

        ampXout_G_PC_sum_norm_d(d) = 0; out_lev_G_PC_sum_norm_d(d) = 0;
        dfG_PC_sum_68(d)=0;  dfG_PC_sum_50(d)=NaN; fG_PC_sum_sort_68(d)=0; fG_PC_sum_sort_50(d)=0; fG_PC_sum_sort_CM(d)=0;

        ampXout_G_PC_norm_d(d) = 0; out_lev_G_PC_norm_d(d) = 0;
        dfG_PC_68(d)=0; dfG_PC_50(d)=NaN; fG_PC_sort_68(d)=0; fG_PC_sort_50(d)=0; fG_PC_sort_CM(d)=0;

        RPW_Fra_G_PC1(d) = 0; RPW_Fra_G_PC(d) = 0; RPW_Fra_G_PC_sum(d) = 0; out_lev_G_dom_d(d) = 0;
        %%%%%%%%%%%%%%
        meanOnsheet(d)=0; Onsheet_CV(d)=0;
        meanInair(d)=TEMPO(end)-TEMPO(1);
        Inair_CV(d)=0;

        if meanInair(d) < threshpause
            meanInair_nopause(d)=TEMPO(end)-TEMPO(1);
            meanpause_duration(d)=0;
        else
            meanInair_nopause(d)=0;
            meanpause_duration(d)=TEMPO(end)-TEMPO(1);
        end

        Inair_nopause_CV(d)=0; pause_duration_CV(d)=0; pause_num(d)=0;
        AirSheet_Ratio_Nopause(d)=0;  AirSheet_Ratio_pause(d)=0;
    end
end % fine analisi su singolo elemento di DataCut

%figure; plot(data.press); hold on; plot(pressNcrop);
Onsheet_ratio=sum(onsheet, 'all')/Execution_Time;
%% calcolo valori medi su elementi di DataCut

A_sum_peak = mean(nonzeros(A_sum_peak_d));
if isnan(mean(nonzeros(A_sum_peak_F_d)))
    A_sum_peak_F= 0;
else
    A_sum_peak_F= mean(nonzeros(A_sum_peak_F_d));
end

RPW_A3D_Vol= mean(nonzeros(RPW_A3D_Vol_d)); RPW_A3D_PC= mean(nonzeros(RPW_A3D_PC_d)); RPW_A3D_PD= mean(nonzeros(RPW_A3D_PD_d)); RPW_A3D_Phy= mean(nonzeros(RPW_A3D_Phy_d));
A3D_peak=mean(nonzeros(A3D_peak_d)); A3D_peak_F=mean(A3D_peak_F_d, 'omitnan');
psd_bw_A3D = mean(nonzeros(psd_bw_A3D_d)); psd_median_A3D = mean(nonzeros(psd_median_A3D_d));
RPW_2_A3D_ab = mean(nonzeros(RPW_2_A3D_ab_d)); RPW_5_A3D_ab = mean(nonzeros(RPW_5_A3D_ab_d)); RPW_8_A3D_ab = mean(nonzeros(RPW_8_A3D_ab_d)); RPW_11_A3D_ab = mean(nonzeros(RPW_11_A3D_ab_d));
RPW_2_A3D_rel = mean(nonzeros(RPW_2_A3D_rel_d)); RPW_5_A3D_rel = mean(nonzeros(RPW_5_A3D_rel_d)); RPW_8_A3D_rel = mean(nonzeros(RPW_8_A3D_rel_d)); RPW_11_A3D_rel = mean(nonzeros(RPW_11_A3D_rel_d));
A3D_peak_new = mean(nonzeros(A3D_peak_d_new)); A3D_peak_F_new=mean(A3D_peak_F_d_new, 'omitnan');
psd_bw_A3D_new = mean(nonzeros(psd_bw_A3D_d_new)); psd_median_A3D_new = mean(nonzeros(psd_median_A3D_d_new));

RPW_A_PC1_PC= mean(nonzeros(RPW_A_PC1_PC_d)); RPW_A_PC1_PD= mean(nonzeros(RPW_A_PC1_PD_d)); RPW_A_PC1_Phy= mean(nonzeros(RPW_A_PC1_Phy_d));
A_PC1_peak=mean(nonzeros(A_PC1_peak_d)); A_PC1_peak_F=mean(nonzeros(A_PC1_peak_F_d));
psd_bw_A_PC1 = mean(nonzeros(psd_bw_A_PC1_d)); psd_median_A_PC1 = mean(nonzeros(psd_median_A_PC1_d));
RPW_2_A_PC1_ab = mean(nonzeros(RPW_2_A_PC1_ab_d)); RPW_5_A_PC1_ab = mean(nonzeros(RPW_5_A_PC1_ab_d)); RPW_8_A_PC1_ab = mean(nonzeros(RPW_8_A_PC1_ab_d)); RPW_11_A_PC1_ab = mean(nonzeros(RPW_11_A_PC1_ab_d));
RPW_2_A_PC1_rel = mean(nonzeros(RPW_2_A_PC1_rel_d)); RPW_5_A_PC1_rel = mean(nonzeros(RPW_5_A_PC1_rel_d)); RPW_8_A_PC1_rel = mean(nonzeros(RPW_8_A_PC1_rel_d)); RPW_11_A_PC1_rel = mean(nonzeros(RPW_11_A_PC1_rel_d));
A_PC1_peak_new = mean(nonzeros(A_PC1_peak_d_new)); A_PC1_peak_F_new=mean(nonzeros(A_PC1_peak_F_d_new));
psd_bw_A_PC1_new = mean(nonzeros(psd_bw_A_PC1_d_new)); psd_median_A_PC1_new = mean(nonzeros(psd_median_A_PC1_d_new));

RPW_A_PC_sum_PC= mean(nonzeros(RPW_A_PC_sum_PC_d)); RPW_A_PC_sum_PD= mean(nonzeros(RPW_A_PC_sum_PD_d)); RPW_A_PC_sum_Phy= mean(nonzeros(RPW_A_PC_sum_Phy_d));
A_PC_sum_peak=mean(nonzeros(A_PC_sum_peak_d)); A_PC_sum_peak_F=mean(nonzeros(A_PC_sum_peak_F_d));
psd_bw_A_PC_sum = mean(nonzeros(psd_bw_A_PC_sum_d)); psd_median_A_PC_sum = mean(nonzeros(psd_median_A_PC_sum_d));
RPW_2_A_PC_sum_ab = mean(nonzeros(RPW_2_A_PC_sum_ab_d)); RPW_5_A_PC_sum_ab = mean(nonzeros(RPW_5_A_PC_sum_ab_d)); RPW_8_A_PC_sum_ab = mean(nonzeros(RPW_8_A_PC_sum_ab_d)); RPW_11_A_PC_sum_ab = mean(nonzeros(RPW_11_A_PC_sum_ab_d));
RPW_2_A_PC_sum_rel = mean(nonzeros(RPW_2_A_PC_sum_rel_d)); RPW_5_A_PC_sum_rel = mean(nonzeros(RPW_5_A_PC_sum_rel_d)); RPW_8_A_PC_sum_rel = mean(nonzeros(RPW_8_A_PC_sum_rel_d)); RPW_11_A_PC_sum_rel = mean(nonzeros(RPW_11_A_PC_sum_rel_d));
A_PC_sum_peak_new = mean(nonzeros(A_PC_sum_peak_d_new)); A_PC_sum_peak_F_new=mean(nonzeros(A_PC_sum_peak_F_d_new));
psd_bw_A_PC_sum_new = mean(nonzeros(psd_bw_A_PC_sum_d_new)); psd_median_A_PC_sum_new = mean(nonzeros(psd_median_A_PC_sum_d_new));

RPW_A_PC_PC= mean(nonzeros(RPW_A_PC_PC_d)); RPW_A_PC_PD= mean(nonzeros(RPW_A_PC_PD_d)); RPW_A_PC_Phy= mean(nonzeros(RPW_A_PC_Phy_d));
A_PC_peak=mean(nonzeros(A_PC_peak_d)); A_PC_peak_F=mean(nonzeros(A_PC_peak_F_d));
psd_bw_A_PC = mean(nonzeros(psd_bw_A_PC_d)); psd_median_A_PC = mean(nonzeros(psd_median_A_PC_d));
RPW_2_A_PC_ab = mean(nonzeros(RPW_2_A_PC_ab_d)); RPW_5_A_PC_ab = mean(nonzeros(RPW_5_A_PC_ab_d)); RPW_8_A_PC_ab = mean(nonzeros(RPW_8_A_PC_ab_d)); RPW_11_A_PC_ab = mean(nonzeros(RPW_11_A_PC_ab_d));
RPW_2_A_PC_rel = mean(nonzeros(RPW_2_A_PC_rel_d)); RPW_5_A_PC_rel = mean(nonzeros(RPW_5_A_PC_rel_d)); RPW_8_A_PC_rel = mean(nonzeros(RPW_8_A_PC_rel_d)); RPW_11_A_PC_rel = mean(nonzeros(RPW_11_A_PC_rel_d));
A_PC_peak_new = mean(nonzeros(A_PC_peak_d_new)); A_PC_peak_F_new=mean(nonzeros(A_PC_peak_F_d_new));
psd_bw_A_PC_new = mean(nonzeros(psd_bw_A_PC_d_new)); psd_median_A_PC_new = mean(nonzeros(psd_median_A_PC_d_new));

RPW_G_PC1_PC= mean(nonzeros(RPW_G_PC1_PC_d)); RPW_G_PC1_PD= mean(nonzeros(RPW_G_PC1_PD_d)); RPW_G_PC1_Phy= mean(nonzeros(RPW_G_PC1_Phy_d));
G_PC1_peak=mean(nonzeros(G_PC1_peak_d)); G_PC1_peak_F=mean(nonzeros(G_PC1_peak_F_d));
psd_bw_G_PC1 = mean(nonzeros(psd_bw_G_PC1_d)); psd_median_G_PC1 = mean(nonzeros(psd_median_G_PC1_d));
RPW_2_G_PC1_ab = mean(nonzeros(RPW_2_G_PC1_ab_d)); RPW_5_G_PC1_ab = mean(nonzeros(RPW_5_G_PC1_ab_d)); RPW_8_G_PC1_ab = mean(nonzeros(RPW_8_G_PC1_ab_d)); RPW_11_G_PC1_ab = mean(nonzeros(RPW_11_G_PC1_ab_d));
RPW_2_G_PC1_rel = mean(nonzeros(RPW_2_G_PC1_rel_d)); RPW_5_G_PC1_rel = mean(nonzeros(RPW_5_G_PC1_rel_d)); RPW_8_G_PC1_rel = mean(nonzeros(RPW_8_G_PC1_rel_d)); RPW_11_G_PC1_rel = mean(nonzeros(RPW_11_G_PC1_rel_d));
G_PC1_peak_new = mean(nonzeros(G_PC1_peak_d_new)); G_PC1_peak_F_new=mean(nonzeros(G_PC1_peak_F_d_new));
psd_bw_G_PC1_new = mean(nonzeros(psd_bw_G_PC1_d_new)); psd_median_G_PC1_new = mean(nonzeros(psd_median_G_PC1_d_new));

RPW_G_PC_sum_PC= mean(nonzeros(RPW_G_PC_sum_PC_d)); RPW_G_PC_sum_PD= mean(nonzeros(RPW_G_PC_sum_PD_d)); RPW_G_PC_sum_Phy= mean(nonzeros(RPW_G_PC_sum_Phy_d));
G_PC_sum_peak=mean(nonzeros(G_PC_sum_peak_d)); G_PC_sum_peak_F=mean(nonzeros(G_PC_sum_peak_F_d));
psd_bw_G_PC_sum = mean(nonzeros(psd_bw_G_PC_sum_d)); psd_median_G_PC_sum = mean(nonzeros(psd_median_G_PC_sum_d));
RPW_2_G_PC_sum_ab = mean(nonzeros(RPW_2_G_PC_sum_ab_d)); RPW_5_G_PC_sum_ab = mean(nonzeros(RPW_5_G_PC_sum_ab_d)); RPW_8_G_PC_sum_ab = mean(nonzeros(RPW_8_G_PC_sum_ab_d)); RPW_11_G_PC_sum_ab = mean(nonzeros(RPW_11_G_PC_sum_ab_d));
RPW_2_G_PC_sum_rel = mean(nonzeros(RPW_2_G_PC_sum_rel_d)); RPW_5_G_PC_sum_rel = mean(nonzeros(RPW_5_G_PC_sum_rel_d)); RPW_8_G_PC_sum_rel = mean(nonzeros(RPW_8_G_PC_sum_rel_d)); RPW_11_G_PC_sum_rel = mean(nonzeros(RPW_11_G_PC_sum_rel_d));
G_PC_sum_peak_new = mean(nonzeros(G_PC_sum_peak_d_new)); G_PC_sum_peak_F_new=mean(nonzeros(G_PC_sum_peak_F_d_new));
psd_bw_G_PC_sum_new = mean(nonzeros(psd_bw_G_PC_sum_d_new)); psd_median_G_PC_sum_new = mean(nonzeros(psd_median_G_PC_sum_d_new));

RPW_G_PC_PC= mean(nonzeros(RPW_G_PC_PC_d)); RPW_G_PC_PD= mean(nonzeros(RPW_G_PC_PD_d)); RPW_G_PC_Phy= mean(nonzeros(RPW_G_PC_Phy_d));
G_PC_peak=mean(nonzeros(G_PC_peak_d)); G_PC_peak_F=mean(nonzeros(G_PC_peak_F_d));
psd_bw_G_PC = mean(nonzeros(psd_bw_G_PC_d)); psd_median_G_PC = mean(nonzeros(psd_median_G_PC_d));
RPW_2_G_PC_ab = mean(nonzeros(RPW_2_G_PC_ab_d)); RPW_5_G_PC_ab = mean(nonzeros(RPW_5_G_PC_ab_d)); RPW_8_G_PC_ab = mean(nonzeros(RPW_8_G_PC_ab_d)); RPW_11_G_PC_ab = mean(nonzeros(RPW_11_G_PC_ab_d));
RPW_2_G_PC_rel = mean(nonzeros(RPW_2_G_PC_rel_d)); RPW_5_G_PC_rel = mean(nonzeros(RPW_5_G_PC_rel_d)); RPW_8_G_PC_rel = mean(nonzeros(RPW_8_G_PC_rel_d)); RPW_11_G_PC_rel = mean(nonzeros(RPW_11_G_PC_rel_d));
G_PC_peak_new = mean(nonzeros(G_PC_peak_d_new)); G_PC_peak_F_new=mean(nonzeros(G_PC_peak_F_d_new));
psd_bw_G_PC_new = mean(nonzeros(psd_bw_G_PC_d_new)); psd_median_G_PC_new = mean(nonzeros(psd_median_G_PC_d_new));

RPW_tr_PC= mean(nonzeros(RPW_tr_PC_d)); RPW_tr_PD= mean(nonzeros(RPW_tr_PD_d)); RPW_tr_Phy= mean(nonzeros(RPW_tr_Phy_d));
tr_peak=mean(nonzeros(tr_peak_d)); tr_peak_F=mean(nonzeros(tr_peak_F_d));
psd_bw_tr = mean(nonzeros(psd_bw_tr_d)); psd_median_tr = mean(nonzeros(psd_median_tr_d));
RPW_2_tr_ab = mean(nonzeros(RPW_2_tr_ab_d)); RPW_5_tr_ab = mean(nonzeros(RPW_5_tr_ab_d)); RPW_8_tr_ab = mean(nonzeros(RPW_8_tr_ab_d)); RPW_11_tr_ab = mean(nonzeros(RPW_11_tr_ab_d));
RPW_2_tr_rel = mean(nonzeros(RPW_2_tr_rel_d)); RPW_5_tr_rel = mean(nonzeros(RPW_5_tr_rel_d)); RPW_8_tr_rel = mean(nonzeros(RPW_8_tr_rel_d)); RPW_11_tr_rel = mean(nonzeros(RPW_11_tr_rel_d));
tr_peak_new = mean(nonzeros(tr_peak_d_new)); tr_peak_F_new=mean(nonzeros(tr_peak_F_d_new));
psd_bw_tr_new = mean(nonzeros(psd_bw_tr_d_new)); psd_median_tr_new = mean(nonzeros(psd_median_tr_d_new));

RPW_G3D_Vol= mean(nonzeros(RPW_G3D_Vol_d)); RPW_G3D_PC= mean(nonzeros(RPW_G3D_PC_d)); RPW_G3D_PD= mean(nonzeros(RPW_G3D_PD_d)); RPW_G3D_Phy= mean(nonzeros(RPW_G3D_Phy_d));
G3D_peak=mean(nonzeros(G3D_peak_d)); G3D_peak_F=mean(G3D_peak_F_d, 'omitnan');
psd_bw_G3D = mean(nonzeros(psd_bw_G3D_d)); psd_median_G3D = mean(nonzeros(psd_median_G3D_d));
RPW_2_G3D_ab = mean(nonzeros(RPW_2_G3D_ab_d)); RPW_5_G3D_ab = mean(nonzeros(RPW_5_G3D_ab_d)); RPW_8_G3D_ab = mean(nonzeros(RPW_8_G3D_ab_d)); RPW_11_G3D_ab = mean(nonzeros(RPW_11_G3D_ab_d));
RPW_2_G3D_rel = mean(nonzeros(RPW_2_G3D_rel_d)); RPW_5_G3D_rel = mean(nonzeros(RPW_5_G3D_rel_d)); RPW_8_G3D_rel = mean(nonzeros(RPW_8_G3D_rel_d)); RPW_11_G3D_rel = mean(nonzeros(RPW_11_G3D_rel_d));
G3D_peak_new = mean(nonzeros(G3D_peak_d_new)); G3D_peak_F_new=mean(nonzeros(G3D_peak_F_d_new));
psd_bw_G3D_new = mean(nonzeros(psd_bw_G3D_d_new)); psd_median_G3D_new = mean(nonzeros(psd_median_G3D_d_new));

psd_bw_tilt= mean(nonzeros(psd_bw_tilt_d)); psd_median_tilt=mean(nonzeros(psd_median_tilt_d));


%calcolo valori medi su elementi di datacut x feature calcolate su assi
%singoli
for k = 1:3
    LDLJ_A(k) =mean(nonzeros(LDLJ_A_k(:,k)));  LDLJ_G(k) =mean(nonzeros(LDLJ_G_k(:,k)));

    S1s(k) =mean(nonzeros(S1_sm(:,k))); S2s(k) =mean(nonzeros(S2_sm(:,k)));
    S3s(k) =mean(nonzeros(S3_sm(:,k))); S4s(k) =mean(nonzeros(S4_sm(:,k)));
    S45s(k) =mean(nonzeros(S45_sm(:,k))); S5s(k) =mean(nonzeros(S5_sm(:,k)));

    S1(k) =mean(nonzeros(S1_k(:,k))); S2(k) =mean(nonzeros(S2_k(:,k)));
    S3(k) =mean(nonzeros(S3_k(:,k))); S4(k) =mean(nonzeros(S4_k(:,k)));
    S45(k) =mean(nonzeros(S45_k(:,k))); S5(k) =mean(nonzeros(S5_k(:,k)));
    RPWG_M(k) =mean(nonzeros(RPWG_manus_s(:,k)));  RMSG_M(k) =mean(nonzeros(RMSG_manus_s(:,k))); m_M(k) =mean(nonzeros(m_manus_s(:,k)));

    RPW_A_Vol(k)=mean(nonzeros(RPW_A_Vol_d(:,k))); RPW_A_PC(k)=mean(nonzeros(RPW_A_PC_d(:,k)));
    RPW_A_PD_axes(k)=mean(nonzeros(RPW_A_PD_d(:,k))); RPW_A_Phy(k)=mean(nonzeros(RPW_A_Phy_d(:,k)));
    A_Peak(k)=mean(nonzeros(A_peak_d(:,k)));  A_Peak_F(k)=mean(A_peak_F_d(:,k), 'omitnan');
    psd_bw_A(k)=mean(nonzeros(psd_bw_A_d(:,k))); psd_median_A(k)=mean(nonzeros(psd_median_A_d(:,k)));
    RPW_2_A_ab=mean(nonzeros(RPW_2_A_ab_d(:,k))); RPW_5_A_ab=mean(nonzeros(RPW_5_A_ab_d(:,k)));RPW_8_A_ab=mean(nonzeros(RPW_8_A_ab_d(:,k)));RPW_11_A_ab=mean(nonzeros(RPW_11_A_ab_d(:,k)));
    RPW_2_A_rel=mean(nonzeros(RPW_2_A_rel_d(:,k))); RPW_5_A_rel=mean(nonzeros(RPW_5_A_rel_d(:,k)));RPW_8_A_rel=mean(nonzeros(RPW_8_A_rel_d(:,k)));RPW_11_A_rel=mean(nonzeros(RPW_11_A_rel_d(:,k)));
    A_Peak_new(k)=mean(nonzeros(A_peak_d_new(:,k)));  A_Peak_F_new(k)=mean(A_peak_F_d_new(:,k), 'omitnan');
    psd_bw_A_new(k)=mean(nonzeros(psd_bw_A_d_new(:,k))); psd_median_A(k)=mean(nonzeros(psd_median_A_d_new(:,k)));

    RPW_A_tr_PC(k)=mean(nonzeros(RPW_A_tr_PC_d(:,k))); RPW_A_tr_PD(k)=mean(nonzeros(RPW_A_tr_PD_d(:,k))); RPW_A_tr_Phy(k)=mean(nonzeros(RPW_A_tr_Phy_d(:,k)));
    A_tr_Peak(k)=mean(nonzeros(A_tr_peak_d(:,k))); A_tr_Peak_F(k)=mean(A_tr_peak_F_d(:,k), 'omitnan');
    psd_bw_A_tr(k)=mean(nonzeros(psd_bw_A_tr_d(:,k))); psd_median_A_tr(k)=mean(nonzeros(psd_median_A_tr_d(:,k)));
    RPW_2_A_tr_ab=mean(nonzeros(RPW_2_A_tr_ab_d(:,k))); RPW_5_A_tr_ab=mean(nonzeros(RPW_5_A_tr_ab_d(:,k)));RPW_8_A_tr_ab=mean(nonzeros(RPW_8_A_tr_ab_d(:,k)));RPW_11_A_tr_ab=mean(nonzeros(RPW_11_A_tr_ab_d(:,k)));
    RPW_2_A_tr_rel=mean(nonzeros(RPW_2_A_tr_rel_d(:,k))); RPW_5_A_tr_rel=mean(nonzeros(RPW_5_A_tr_rel_d(:,k)));RPW_8_A_tr_rel=mean(nonzeros(RPW_8_A_tr_rel_d(:,k)));RPW_11_A_tr_rel=mean(nonzeros(RPW_11_A_tr_rel_d(:,k)));
    A_tr_Peak_new(k)=mean(nonzeros(A_tr_peak_d_new(:,k)));  A_tr_Peak_F_new(k)=mean(A_tr_peak_F_d_new(:,k), 'omitnan');
    psd_bw_A_tr_new(k)=mean(nonzeros(psd_bw_A_tr_d_new(:,k))); psd_median_A_tr_new(k)=mean(nonzeros(psd_median_A_tr_d_new(:,k)));

    RPW_G_Vol(k)=mean(nonzeros(RPW_G_Vol_d(:,k))); RPW_G_PC(k)=mean(nonzeros(RPW_G_PC_d(:,k)));
    RPW_G_PD_axes(k)=mean(nonzeros(RPW_G_PD_d(:,k))); RPW_G_Phy(k)=mean(nonzeros(RPW_G_Phy_d(:,k)));
    G_Peak(k)=mean(nonzeros(G_peak_d(:,k))); G_Peak_F(k)=mean(G_peak_F_d(:,k),'omitnan');
    psd_bw_G(k)=mean(nonzeros(psd_bw_G_d(:,k))); psd_median_G(k)=mean(nonzeros(psd_median_G_d(:,k)));
    RPW_2_G_ab=mean(nonzeros(RPW_2_G_ab_d(:,k))); RPW_5_G_ab=mean(nonzeros(RPW_5_G_ab_d(:,k)));RPW_8_G_ab=mean(nonzeros(RPW_8_G_ab_d(:,k)));RPW_11_G_ab=mean(nonzeros(RPW_11_G_ab_d(:,k)));
    RPW_2_G_rel=mean(nonzeros(RPW_2_G_rel_d(:,k))); RPW_5_G_rel=mean(nonzeros(RPW_5_G_rel_d(:,k)));RPW_8_G_rel=mean(nonzeros(RPW_8_G_rel_d(:,k)));RPW_11_G_rel=mean(nonzeros(RPW_11_G_rel_d(:,k)));
    G_Peak_new(k)=mean(nonzeros(G_peak_d_new(:,k)));  G_Peak_F_new(k)=mean(G_peak_F_d_new(:,k), 'omitnan');
    psd_bw_G_new(k)=mean(nonzeros(psd_bw_G_d_new(:,k))); psd_median_G(k)=mean(nonzeros(psd_median_G_d_new(:,k)));

    RPW_G_PC_m(k)=mean(nonzeros(RPW_G_PC_dm(:,k))); RPW_G_PD_m(k)=mean(nonzeros(RPW_G_PD_dm(:,k))); RPW_G_Phy_m(k)=mean(nonzeros(RPW_G_Phy_dm(:,k)));
    G_Peak_m(k)=mean(nonzeros(G_peak_dm(:,k))); G_Peak_F_m(k)=mean(G_peak_F_dm(:,k),'omitnan');
    psd_bw_G_m(k)=mean(nonzeros(psd_bw_G_dm(:,k))); psd_median_G_m(k)=mean(nonzeros(psd_median_G_dm(:,k)));
    RPW_2_G_dm_ab=mean(nonzeros(RPW_2_G_dm_ab_d(:,k))); RPW_5_G_dm_ab=mean(nonzeros(RPW_5_G_dm_ab_d(:,k)));RPW_8_G_dm_ab=mean(nonzeros(RPW_8_G_dm_ab_d(:,k)));RPW_11_G_dm_ab=mean(nonzeros(RPW_11_G_dm_ab_d(:,k)));
    RPW_2_G_dm_rel=mean(nonzeros(RPW_2_G_dm_rel_d(:,k))); RPW_5_G_dm_rel=mean(nonzeros(RPW_5_G_dm_rel_d(:,k)));RPW_8_G_dm_rel=mean(nonzeros(RPW_8_G_dm_rel_d(:,k)));RPW_11_G_dm_rel=mean(nonzeros(RPW_11_G_dm_rel_d(:,k)));
    G_dm_Peak_new(k)=mean(nonzeros(G_dm_peak_d_new(:,k)));  G_dm_Peak_F_new(k)=mean(G_dm_peak_F_d_new(:,k), 'omitnan');
    psd_bw_G_m_new(k)=mean(nonzeros(psd_bw_G_dm_new(:,k))); psd_median_G_m_new(k)=mean(nonzeros(psd_median_G_dm_new(:,k)));

    RPW_Fra_m(k) = mean(RPW_Fra(:,k), 'omitnan');
    RPW_Fra_m_new(k) = mean(RPW_Fra_new(:,k), 'omitnan');
end

%indice m manus max in valore assoluto
m_abs = abs(m_M); [~, ind_max_m_abs] = max(m_abs);
%considerando componente con picco massimo
A_Peak_max = max(A_Peak); A_Peak_F_max = A_Peak_F(A_Peak == A_Peak_max);
A_Peak_max_new = max(A_Peak_new); A_Peak_F_max_new = A_Peak_F_new(A_Peak_new == A_Peak_max_new);

%Acc tremore
RPW_A_tr_dominant_PC=mean(nonzeros(RPW_A_tr_dominant_PC_d)); RPW_A_tr_dominant_PD=mean(nonzeros(RPW_A_tr_dominant_PD_d));
RPW_A_tr_dominant_Phy=mean(nonzeros(RPW_A_tr_dominant_Phy_d));
A_tr_Peak_max = max(A_tr_Peak); A_tr_Peak_F_max = A_tr_Peak_F(A_tr_Peak == A_tr_Peak_max);
psd_bw_A_tr_dominant= mean(nonzeros(psd_bw_A_tr_dominant_d)); psd_median_A_tr_dominant= mean(nonzeros(psd_median_A_tr_dominant_d));
RPW_A_tr_dominant_2=mean(nonzeros(RPW_A_tr_dominant_2_d)); RPW_A_tr_dominant_5=mean(nonzeros(RPW_A_tr_dominant_5_d));
RPW_A_tr_dominant_8=mean(nonzeros(RPW_A_tr_dominant_8_d));RPW_A_tr_dominant_11=mean(nonzeros(RPW_A_tr_dominant_11_d));
A_tr_Peak_max_new = max(A_tr_Peak_new); A_tr_Peak_F_max_new = A_tr_Peak_F_new(A_tr_Peak_new == A_tr_Peak_max_new);
psd_bw_A_tr_dominant_new= mean(nonzeros(psd_bw_A_tr_dominant_d_new)); psd_median_A_tr_dominant_new= mean(nonzeros(psd_median_A_tr_dominant_d_new));

%Gyro totale
G_Peak_max = max(G_Peak); G_Peak_F_max = G_Peak_F(G_Peak == G_Peak_max);
G_Peak_max_new = max(G_Peak_new); G_Peak_F_max_new = G_Peak_F_new(G_Peak_new == G_Peak_max_new);

%Gyro tremore
RPW_G_tr_dominant_PC=mean(nonzeros(RPW_G_tr_dominant_PC_d)); RPW_G_tr_dominant_PD=mean(nonzeros(RPW_G_tr_dominant_PD_d));
RPW_G_tr_dominant_Phy=mean(nonzeros(RPW_G_tr_dominant_Phy_d));
G_Peak_max_m = max(G_Peak_m); G_Peak_F_max_m = G_Peak_F_m(G_Peak_m == G_Peak_max_m);
psd_bw_G_tr_dominant= mean(nonzeros(psd_bw_G_tr_dominant_d)); psd_median_G_tr_dominant= mean(nonzeros(psd_median_G_tr_dominant_d));
RPW_G_tr_dominant_2=mean(nonzeros(RPW_G_tr_dominant_2_d)); RPW_G_tr_dominant_5=mean(nonzeros(RPW_G_tr_dominant_5_d));
RPW_G_tr_dominant_8=mean(nonzeros(RPW_G_tr_dominant_8_d));RPW_G_tr_dominant_11=mean(nonzeros(RPW_G_tr_dominant_11_d));
G_Peak_max_m_new = max(G_dm_Peak_new); G_Peak_F_max_m_new = G_dm_Peak_F_new(G_dm_Peak_new == G_Peak_max_m_new);
psd_bw_G_tr_dominant_new= mean(nonzeros(psd_bw_G_tr_dominant_d_new)); psd_median_G_tr_dominant_new= mean(nonzeros(psd_median_G_tr_dominant_d_new));
%% salvataggio indici
% medie sono calcolate sul numero di segmenti senza perdite
%TEMPORALI
PEN_test.ExecutionTime=Execution_Time;
PEN_test.StrokeNum=Num_Strokes;
PEN_test.RelStrokeNum=Relative_Num_Strokes;

PEN_test.meanOnsheet=mean(nonzeros(meanOnsheet));
PEN_test.Onsheet_CV=mean(nonzeros(Onsheet_CV));
PEN_test.OnSheetRatio=Onsheet_ratio;

PEN_test.meanInair=mean(nonzeros(meanInair));
PEN_test.Inair_CV=mean(nonzeros(Inair_CV));
PEN_test.meanInair_nopause=mean(nonzeros(meanInair_nopause));

if mean(Inair_nopause_CV)>0
    PEN_test.Inair_nopause_CV=mean(nonzeros(Inair_nopause_CV));
else
    PEN_test.Inair_nopause_CV=0;
end

if isnan(mean(nonzeros(meanpause_duration)))
    PEN_test.meanPause=0;
else
    PEN_test.meanPause=mean(nonzeros(meanpause_duration));
end

if mean(pause_duration_CV)>0
    PEN_test.Pause_CV=mean(nonzeros(pause_duration_CV));
else
    PEN_test.Pause_CV=0;
end
PEN_test.PauseNum=sum(pause_num);
PEN_test.PauseNum_Rel=sum(pause_num)/Execution_Time;

PEN_test.AirSheet_Ratio_NOpause=mean(nonzeros(AirSheet_Ratio_Nopause));
PEN_test.AirSheet_Ratio_pause=mean(nonzeros(AirSheet_Ratio_pause));
%PRESSIONE
PEN_test.meanPonsheet=meanP;
PEN_test.meanP=mean(nonzeros(meanP_DC));
PEN_test.P_CV=mean(nonzeros(meanP_CV_DC));
PEN_test.P_OVS=mean(nonzeros(meanPOVS_DC));
%FLUIDITA'
PEN_test.NCP_REL=mean(nonzeros(meanNCP_DC));
PEN_test.NCG_REL=mean(nonzeros(meanNCG_DC));
PEN_test.NCA_REL=mean(nonzeros(meanNCA_DC));
%PEN_test.NCV_REL=mean(nonzeros(meanNCV_DC));

PEN_test.ConsPeakDiffPMean=mean(nonzeros(meanCPDP_AVG));
PEN_test.ConsPeakDiffPCV=mean(nonzeros(meanCPDP_CV));
PEN_test.ConsPeakDiffGMean=mean(nonzeros(meanCPDG_AVG));
if isnan(mean(nonzeros(meanCPDG_CV)))
    PEN_test.ConsPeakDiffGCV=0;
else
    PEN_test.ConsPeakDiffGCV=mean(nonzeros(meanCPDG_CV));
end
% PEN_test.ConsPeakDiffVMean=mean(nonzeros(meanCPDV_AVG));
%
% if isnan(mean(nonzeros(meanCPDV_CV)))
%    PEN_test.ConsPeakDiffVCV=0;
% else
%
%    PEN_test.ConsPeakDiffVCV=mean(nonzeros(meanCPDV_CV));
% end
%
% PEN_test.meanV=mean(nonzeros(meanV_DC)); PEN_test.VelCV=mean(nonzeros(V_CV_AVG)); PEN_test.VelApEN=mean(nonzeros(V_ApEn_AVG));

PEN_test.TILT=mean(meanTILT);
PEN_test.TILT_CV=mean(meanTILT_CV);
PEN_test.TILT_VAR=mean(meanTILT_VAR);
% if isnan(mean(nonzeros(SNvpd_norm_AVG)))
%     PEN_test.SNvpd_REL=0;
% else
%     PEN_test.SNvpd_REL=mean(nonzeros(SNvpd_norm_AVG));
% end
%
% PEN_test.medianDLJ=mean(nonzeros(medianDLJ_A));
% PEN_test.medianSJ=mean(nonzeros(medianSquaredJerk_A));

%feature calderon min
% PEN_test.LDLJ_A_min=min(LDLJ_A); PEN_test.LDLJ_G_min=min(LDLJ_G);
% PEN_test.S1_min=min(S1); PEN_test.S2_min=min(S2); PEN_test.S3_min=min(S3);
% PEN_test.S4_min=min(S4); PEN_test.S45_min=min(S45); PEN_test.S5_min=min(S5);

%feature calderon median
PEN_test.LDLJ_A_median=median(LDLJ_A);
PEN_test.LDLJ_G_median=median(LDLJ_G);

PEN_test.S1_median_s=median(S1s);
PEN_test.S2_median_s=median(S2s);
PEN_test.S3_median_s=median(S3s);
PEN_test.S4_median_s=median(S4s);
PEN_test.S45_median_s=median(S45s);
PEN_test.S5_median_s=median(S5s);


PEN_test.S1_median=median(S1);
PEN_test.S2_median=median(S2);
PEN_test.S3_median=median(S3);
PEN_test.S4_median=median(S4);
PEN_test.S45_median=median(S45);
PEN_test.S5_median=median(S5);

%%%%%%%%%%%%%%%% Per dati busto no da qui in avanti
%feature calderon max
% PEN_test.LDLJ_A_max=max(LDLJ_A); PEN_test.LDLJ_G_max=max(LDLJ_G);
% PEN_test.S1_max=max(S1); PEN_test.S2_max=max(S2); PEN_test.S3_max=max(S3);
% PEN_test.S4_max=max(S4); PEN_test.S45_max=max(S45); PEN_test.S5_max=max(S5);

%feature manus median
PEN_test.RPWG_M_median=median(RPWG_M);

if isnan(median(RMSG_M))
    PEN_test.RMSG_M_median=0;
else
    PEN_test.RMSG_M_median=median(RMSG_M);
end

if isnan(median(m_M))
    PEN_test.m_M_median=0;
else
    PEN_test.m_M_median=median(m_M);
end
%feature manus max
PEN_test.RPWG_M_max=max(RPWG_M);
if isnan(max(RMSG_M))
    PEN_test.RMSG_M_max=0;
else
    PEN_test.RMSG_M_max=max(RMSG_M);
end
if isnan(max(m_M))
    PEN_test.m_M_max=0;
else
    PEN_test.m_M_max=m_M(ind_max_m_abs);
end
%PEN_test.V_ApEn=mean(nonzeros(meanV_APEN));
%TREMORE
if isnan(mean(nonzeros(TM_ApEn))) %capita se non ci sono segmenti senza perdite > 500 samples
    PEN_test.TREMOR_ApEn=0;
    PEN_test.meanRR=0;
    PEN_test.meanDET=0;
    %PEN_test.Tremor_RMS=Tremor_Amplitude_RMS;
    PEN_test.TSI_DB_IMF=0;
    PEN_test.TSI_DB_thr_IMF=0;
    PEN_test.TSI_Luft_IMF=0;
    PEN_test.SNR_IMF=0;
    PEN_test.MHP_IMF=0;

    %RMS A
    PEN_test.RMS_A_TOT=0;
    %     PEN_test.RMS_A_VOL=0;  PEN_test.RMS_A_DYS=0; PEN_test.RMS_A_PD=0; PEN_test.RMS_A_PHY=0;
    %RMS G
    PEN_test.RMS_G_TOT=0;
    %     %PEN_test.RMS_G_VOL=0; PEN_test.RMS_G_DYS=0; PEN_test.RMS_G_PD=0; PEN_test.RMS_G_PHY=0;
    %   %RMS T
    %     PEN_test.RMS_T_TOT=0; PEN_test.RMS_T_LF=0; PEN_test.RMS_T_PD=0;  PEN_test.RMS_T_PHY=0;
    %   %RPW A
    %     PEN_test.RPW_A_VOL=0;  PEN_test.RPW_A_DYS=0;  PEN_test.RPW_A_PD=0; PEN_test.RPW_A_PHY=0;
    %   %RPW G
    %     PEN_test.RPW_G_VOL=0; PEN_test.RPW_G_DYS=0; PEN_test.RPW_G_PD=0;  PEN_test.RPW_G_PHY=0;
    %   %RPW T
    %     PEN_test.RPW_T_LF=0; PEN_test.RPW_T_PD=0;  PEN_test.RPW_T_PHY=0;
    %   %PSD A
    %     PEN_test.F_MODAL_A=0; PEN_test.PEAK_POW_A=0;
    %     %PSD T
    %     PEN_test.F_MODAL_T=0; PEN_test.F_MODAL_T_IRENE=0;
    %     PSD G
    %     PEN_test.F_MODAL_G=0; PEN_test.PEAK_POW_G=0;
else
    PEN_test.TREMOR_ApEn=mean(nonzeros(TM_ApEn));
    PEN_test.meanRR=mean(nonzeros(TM_RR));
    PEN_test.meanDET=mean(nonzeros(TM_DET));
    %PEN_test.Tremor_RMS=Tremor_Amplitude_RMS;
    PEN_test.TSI_DB_IMF=mean(nonzeros(TSI_DB_IMF));
    PEN_test.TSI_DB_thr_IMF=mean(nonzeros(TSI_DB_thr_IMF));
    PEN_test.TSI_Luft_IMF=mean(TSI_Luft_IMF, 'all');
    PEN_test.SNR_IMF=mean(nonzeros(SNR_IMF));
    PEN_test.MHP_IMF=mean(nonzeros(MHP_IMF));
    %RMS A
    PEN_test.RMS_A_TOT=mean(nonzeros(RMS_A_TOT));
    %     PEN_test.RMS_A_VOL=mean(nonzeros(RMS_A_VOL)); PEN_test.RMS_A_DYS=mean(nonzeros(RMS_A_DYS));
    %     PEN_test.RMS_A_PD=mean(nonzeros(RMS_A_PD)); PEN_test.RMS_A_PHY=mean(nonzeros(RMS_A_PHY));
    %   %RMS G
    PEN_test.RMS_G_TOT=mean(nonzeros(RMS_G_TOT));
    %     %PEN_test.RMS_G_VOL=mean(nonzeros(RMS_G_VOL)); PEN_test.RMS_G_DYS=mean(nonzeros(RMS_G_DYS));
    %     PEN_test.RMS_G_PD=mean(nonzeros(RMS_G_PD)); PEN_test.RMS_G_PHY=mean(nonzeros(RMS_G_PHY));
    %     %RMS T
    %     PEN_test.RMS_T_TOT=mean(nonzeros(RMS_T_TOT)); PEN_test.RMS_T_LF=mean(nonzeros(RMS_T_LF));
    %     PEN_test.RMS_T_PD=mean(nonzeros(RMS_T_PD)); PEN_test.RMS_T_PHY=mean(nonzeros(RMS_T_PHY));
    %     %RPW A
    %     PEN_test.RPW_A_VOL=mean(nonzeros(RPW_A_VOL)); PEN_test.RPW_A_DYS=mean(nonzeros(RPW_A_DYS));
    %     PEN_test.RPW_A_PD=mean(nonzeros(RPW_A_PD)); PEN_test.RPW_A_PHY=mean(nonzeros(RPW_A_PHY));
    %     %RPW G
    %     PEN_test.RPW_G_VOL=mean(nonzeros(RPW_G_VOL)); PEN_test.RPW_G_DYS=mean(nonzeros(RPW_G_DYS));
    %     PEN_test.RPW_G_PD=mean(nonzeros(RPW_G_PD)); PEN_test.RPW_G_PHY=mean(nonzeros(RPW_G_PHY));
    %     %RPW T
    %     PEN_test.RPW_T_LF=mean(nonzeros(RPW_T_LF)); PEN_test.RPW_T_PD=mean(nonzeros(RPW_T_PD));  PEN_test.RPW_T_PHY=mean(nonzeros(RPW_T_PHY));
    %     %PSD A
    %     PEN_test.F_MODAL_A=mean(nonzeros(F_MODAL_A));PEN_test.PEAK_POW_A=mean(nonzeros(PEAK_POW_A));
    %     %PSD T
    %     PEN_test.F_MODAL_T=mean(nonzeros(F_MODAL_T)); PEN_test.PEAK_POW_T=mean(nonzeros(PEAK_POW_T)); PEN_test.F_MODAL_T_IRENE=mean(nonzeros(F_MODAL_T_IRENE));
    %     %PSD G
    %     PEN_test.F_MODAL_G=mean(nonzeros(F_MODAL_G)); PEN_test.PEAK_POW_G=mean(nonzeros(PEAK_POW_G));
    % end

    PEN_test.RelPOWVOL_a3D_WELCH= RPW_A3D_Vol;
    PEN_test.RelPOWDYS_a3D_WELCH=RPW_A3D_PC;
    PEN_test.RelPOWPDt_a3D_WELCH=RPW_A3D_PD;
    PEN_test.RelPOWPHYt_a3D_WELCH=RPW_A3D_Phy;
    PEN_test.RelPOWband2_a3D_WELCH=RPW_2_A3D_ab;
    PEN_test.RelPOWband5_a3D_WELCH=RPW_5_A3D_ab;
    PEN_test.RelPOWband8_a3D_WELCH=RPW_8_A3D_ab;
    PEN_test.RelPOWband11_a3D_WELCH=RPW_11_A3D_ab;


    %potenza relativa PC acc contributo di tremore
    %PC1
    % PEN_test.RelPOWPC_a_PC1_WELCH=RPW_A_PC1_PC;
    % PEN_test.RelPOWPDt_a_PC1_WELCH=RPW_A_PC1_PD;
    % PEN_test.RelPOWPHYt_a_PC1_WELCH=RPW_A_PC1_Phy;
    % % potenza relativa PC sum
    % PEN_test.RelPOWPC_a_PC_sum_WELCH=RPW_A_PC_sum_PC;
    % PEN_test.RelPOWPDt_a_PC_sum_WELCH=RPW_A_PC_sum_PD;
    % PEN_test.RelPOWPHYt_a_PC_sum_WELCH=RPW_A_PC_sum_Phy;
    % potenza relativa PC
    PEN_test.RelPOWPC_a_PC_WELCH=RPW_A_PC_PC;
    PEN_test.RelPOWPDt_a_PC_WELCH=RPW_A_PC_PD;
    PEN_test.RelPOWPHYt_a_PC_WELCH=RPW_A_PC_Phy;
    PEN_test.RelPOWband2_a_PC_WELCH=RPW_2_A_PC_ab;
    PEN_test.RelPOWband5_a_PC_WELCH=RPW_5_A_PC_ab;
    PEN_test.RelPOWband8_a_PC_WELCH=RPW_8_A_PC_ab;
    PEN_test.RelPOWband11_a_PC_WELCH=RPW_11_A_PC_ab;

    %potenza relativa mediana tra le tre componenti acc
    % PEN_test.RelPOWVOL_a_WELCH_median= median(RPW_A_Vol);
    % PEN_test.RelPOWDYS_a_WELCH_median=median(RPW_A_PC);
    % PEN_test.RelPOWPDt_a_WELCH_median=median(RPW_A_PD_axes);
    % PEN_test.RelPOWPHYt_a_WELCH_median=median(RPW_A_Phy);

    %potenza relativa max tra le tre componenti acc
    % PEN_test.RelPOWVOL_a_WELCH_max= max(RPW_A_Vol);
    % PEN_test.RelPOWDYS_a_WELCH_max=max(RPW_A_PC);
    % PEN_test.RelPOWPDt_a_WELCH_max=max(RPW_A_PD_axes);
    % PEN_test.RelPOWPHYt_a_WELCH_max=max(RPW_A_Phy);

    % %potenza relativa mediana tra le tre componenti acc filtrate passa alto a
    % %2Hz con WELCH
    % PEN_test.RelPOWDYS_a_tr_WELCH_median=median(RPW_A_tr_PC);
    % PEN_test.RelPOWPDt_a_tr_WELCH_median=median(RPW_A_tr_PD);
    % PEN_test.RelPOWPHYt_a_tr_WELCH_median=median(RPW_A_tr_Phy);
    %
    % %potenza relativa max tra le tre componenti acc filtrate passa alto a
    % %2Hz con WELCH
    % PEN_test.RelPOWDYS_a_tr_WELCH_max=max(RPW_A_tr_PC);
    % PEN_test.RelPOWPDt_a_tr_WELCH_max=max(RPW_A_tr_PD);
    % PEN_test.RelPOWPHYt_a_tr_WELCH_max=max(RPW_A_tr_Phy);
    %
    % %potenza relativa della componente di acc filtrata HP a 2Hz con picco
    % %massimo, con WELCH
    PEN_test.RelPOW_a_tr_dominant_PC = RPW_A_tr_dominant_PC;
    PEN_test.RelPOW_a_tr_dominant_PD = RPW_A_tr_dominant_PD;
    PEN_test.RelPOW_a_tr_dominant_Phy = RPW_A_tr_dominant_Phy;
    PEN_test.RelPOW_a_tr_dominant_band2 = RPW_A_tr_dominant_2;
    PEN_test.RelPOW_a_tr_dominant_band5 = RPW_A_tr_dominant_5;
    PEN_test.RelPOW_a_tr_dominant_band8 = RPW_A_tr_dominant_8;
    PEN_test.RelPOW_a_tr_dominant_band11 = RPW_A_tr_dominant_11;

    % PEN_test.RelPOWLF_tr_WELCH=RPW_tr_PC;
    % PEN_test.RelPOWPDt_tr_WELCH=RPW_tr_PD;
    % PEN_test.RelPOWPHYt_tr_WELCH=RPW_tr_Phy;
    %
    % PEN_test.RelPOWVOL_g3D_WELCH= RPW_G3D_Vol;
    % PEN_test.RelPOWDYS_g3D_WELCH=RPW_G3D_PC;
    % PEN_test.RelPOWPDt_g3D_WELCH=RPW_G3D_PD;
    % PEN_test.RelPOWPHYt_g3D_WELCH=RPW_G3D_Phy;

    %potenza relativa PC gyro contributo di tremore
    %PC1
    % PEN_test.RelPOWPC_g_PC1_WELCH=RPW_G_PC1_PC;
    % PEN_test.RelPOWPDt_g_PC1_WELCH=RPW_G_PC1_PD;
    % PEN_test.RelPOWPHYt_g_PC1_WELCH=RPW_G_PC1_Phy;
    % %PC sum
    % PEN_test.RelPOWPC_g_PC_sum_WELCH=RPW_G_PC_sum_PC;
    % PEN_test.RelPOWPDt_g_PC_sum_WELCH=RPW_G_PC_sum_PD;
    % PEN_test.RelPOWPHYt_g_PC_sum_WELCH=RPW_G_PC_sum_Phy;
    % %PC
    PEN_test.RelPOWPC_g_PC_WELCH=RPW_G_PC_PC;
    PEN_test.RelPOWPDt_g_PC_WELCH=RPW_G_PC_PD;
    PEN_test.RelPOWPHYt_g_PC_WELCH=RPW_G_PC_Phy;
    PEN_test.RelPOWband2_g_PC_WELCH = RPW_2_G_PC_ab;
    PEN_test.RelPOWband5_g_PC_WELCH = RPW_5_G_PC_ab;
    PEN_test.RelPOWband8_g_PC_WELCH = RPW_8_G_PC_ab;
    PEN_test.RelPOWband11_g_PC_WELCH = RPW_11_G_PC_ab;


    %potenza relativa mediana tra le tre componenti gyro
    % PEN_test.RelPOWVOL_g_WELCH_median= median(RPW_G_Vol);
    % PEN_test.RelPOWDYS_g_WELCH_median= median(RPW_G_PC);
    % PEN_test.RelPOWPDt_g_WELCH_median= median(RPW_G_PD_axes);
    % PEN_test.RelPOWPHYt_g_WELCH_median= median(RPW_G_Phy);

    %potenza relativa max tra le tre componenti gyro
    % PEN_test.RelPOWVOL_g_WELCH_max= max(RPW_G_Vol);
    % PEN_test.RelPOWDYS_g_WELCH_max= max(RPW_G_PC);
    % PEN_test.RelPOWPDt_g_WELCH_max= max(RPW_G_PD_axes);
    % PEN_test.RelPOWPHYt_g_WELCH_max= max(RPW_G_Phy);
    %
    % % %potenza relativa mediana tra le tre componenti gyro filtrate a 2 Hz
    % PEN_test.RelPOWDYS_g_tr_WELCH_median= median(RPW_G_PC_m);
    % PEN_test.RelPOWPDt_g_tr_WELCH_median= median(RPW_G_PD_m);
    % PEN_test.RelPOWPHYt_g_tr_WELCH_median= median(RPW_G_Phy_m);
    % %
    % % %potenza relativa max tra le tre componenti gyro filtrate a 2Hz
    % PEN_test.RelPOWDYS_g_tr_WELCH_max= max(RPW_G_PC_m);
    % PEN_test.RelPOWPDt_g_tr_WELCH_max= max(RPW_G_PD_m);
    % PEN_test.RelPOWPHYt_g_tr_WELCH_max= max(RPW_G_Phy_m);

    %potenza relativa della componente di gyro filtrata HP a 2Hz con picco
    %massimo
    PEN_test.RelPOW_g_tr_dominant_PC = RPW_G_tr_dominant_PC;
    PEN_test.RelPOW_g_tr_dominant_PD = RPW_G_tr_dominant_PD;
    PEN_test.RelPOW_g_tr_dominant_Phy = RPW_G_tr_dominant_Phy;
    PEN_test.RelPOW_g_tr_dominant_band2 = RPW_G_tr_dominant_2;
    PEN_test.RelPOW_g_tr_dominant_band5 = RPW_G_tr_dominant_5;
    PEN_test.RelPOW_g_tr_dominant_band8 = RPW_G_tr_dominant_8;
    PEN_test.RelPOW_g_tr_dominant_band11 = RPW_G_tr_dominant_11;

    PEN_test.TSI_DB=mean(nonzeros(TSI_DB));
    PEN_test.TSI_DB_thr=mean(nonzeros(TSI_DB_thr));
    PEN_test.TSI_Luft=mean(TSI_Luft);
    PEN_test.SNR=mean(nonzeros(SNR));
    PEN_test.MHP=mean(nonzeros(MHP));

    %FREQUENZE DI PICCO
    PEN_test.PeakF_A3D_WELCH=A3D_peak_F;
    PEN_test.PeakF_A3D_WELCH_new=A3D_peak_F_new;
    % PEN_test.PeakF_A_WELCH=A_Peak_F_max;
    % PEN_test.PeakF_A_sum_WELCH = A_sum_peak_F;
    %Acc tremore Principal Components
    % PEN_test.PeakF_A_PC1 = A_PC1_peak_F;
    % PEN_test.PeakF_A_PC_sum = A_PC_sum_peak_F;
    PEN_test.PeakF_A_PC = A_PC_peak_F;
    PEN_test.PeakF_A_PC_new = A_PC_peak_F_new;%o la prima o la seconda a seconda della varianza spiegata
    %Acc single components
    % PEN_test.PeakF_A_WELCH=A_Peak_F_max;% f modal della componente di acc dominante

    %tremore 3D estratto manualmente, welch
    % PEN_test.PeakF_tr_WELCH=tr_peak_F;
    % %Acc tremor single component
    PEN_test.PeakF_A_tr_WELCH=A_tr_Peak_F_max;
    PEN_test.PeakF_A_tr_WELCH_new=A_tr_Peak_F_max_new; % f modal della componente di acc, filtrata HP a 2Hz, dominante


    % PEN_test.PeakF_G3D_WELCH=G3D_peak_F;
    % PEN_test.PeakF_G_WELCH=G_Peak_F_max; % Singola componente gyro (quella con picco massimo), welch
    % %Gyro tremor principal components
    % PEN_test.PeakF_G_PC1 = G_PC1_peak_F;
    % PEN_test.PeakF_G_PC_sum = G_PC_sum_peak_F;
    PEN_test.PeakF_G_PC = G_PC_peak_F;
    PEN_test.PeakF_G_PC_new = G_PC_peak_F_new;
    % %Gyro tremor single components
    PEN_test.PeakF_G_tr_WELCH=G_Peak_F_max_m; % Singola componente gyro filtrata a 2Hz (quella con picco massimo), welch
    PEN_test.PeakF_G_tr_WELCH_new=G_Peak_F_max_m_new;

    %POTENZE DI PICCO
    PEN_test.Peak_A3D_WELCH=A3D_peak;
    PEN_test.Peak_A3D_WELCH_new=A3D_peak_new;
    % PEN_test.Peak_A_WELCH=A_Peak_max;
    % PEN_test.Peak_A_sum_WELCH = A_sum_peak;
    % %Acc tremor PC
    % PEN_test.Peak_A_PC1 = A_PC1_peak;
    % PEN_test.Peak_A_PC_sum = A_PC_sum_peak;
    PEN_test.Peak_A_PC = A_PC_peak;
    PEN_test.Peak_A_PC_new = A_PC_peak_new;
    %Acc single components
    % PEN_test.Peak_A_WELCH=A_Peak_max; %Singola componente acc (quella con picco massimo), welch

    % PEN_test.Peak_tr_WELCH=tr_peak;
    %Acc tremor Single Components
    PEN_test.Peak_A_tr_WELCH=A_tr_Peak_max;
    PEN_test.Peak_A_tr_WELCH_new=A_tr_Peak_max_new;%Singola componente acc filtrata a 2 Hz (quella con picco massimo), welch

    %%%
    % PEN_test.Peak_G3D_WELCH=G3D_peak;
    % PEN_test.Peak_G_WELCH=G_Peak_max;
    % %Gyro tremor principal components
    % PEN_test.Peak_G_PC1 = G_PC1_peak;
    % PEN_test.Peak_G_PC_sum = G_PC_sum_peak;
    PEN_test.Peak_G_PC = G_PC_peak;
    PEN_test.Peak_G_PC_new = G_PC_peak_new;
    %Gyro tremor single component
    PEN_test.Peak_G_tr_WELCH=G_Peak_max_m;  % Singola componente gyro filtrata a 2Hz (quella con picco massimo), welch
    PEN_test.Peak_G_tr_WELCH_new=G_Peak_max_m_new;

    %Sadikov
    PEN_test.ampXout_A_dom=mean(nonzeros(ampXout_A_dom_d));
    PEN_test.out_lev_A_dom=mean(nonzeros(out_lev_A_dom_d));
    PEN_test.dfA_dom_68=mean(nonzeros(dfA_dom_68));
    PEN_test.dfA_dom_50=mean(dfA_dom_50,'omitnan');
    PEN_test.fA_dom_sort_68=mean(nonzeros(fA_dom_sort_68));
    PEN_test.fA_dom_sort_50=mean(nonzeros(fA_dom_sort_50));
    PEN_test.fA_dom_sort_CM=mean(nonzeros(fA_dom_sort_CM));

    PEN_test.ampXout_A3D=mean(nonzeros(ampXout_A3D_norm_d));
    PEN_test.out_lev_A3D=mean(nonzeros(out_lev_A3D_norm_d));
    PEN_test.dfA3D_68=mean(nonzeros(dfA3D_68));
    PEN_test.dfA3D_50=mean(dfA3D_50,'omitnan');
    PEN_test.fA3D_sort_68=mean(nonzeros(fA3D_sort_68));
    PEN_test.fA3D_sort_50=mean(nonzeros(fA3D_sort_50));
    PEN_test.fA3D_sort_CM=mean(nonzeros(fA3D_sort_CM));

    % PEN_test.ampXout_A_PC1=mean(nonzeros(ampXout_A_PC1_norm_d));
    % PEN_test.out_lev_A_PC1=mean(nonzeros(out_lev_A_PC1_norm_d));
    % PEN_test.dfA_PC1_68=mean(nonzeros(dfA_PC1_68));
    % PEN_test.dfA_PC1_50=mean(dfA_PC1_50,'omitnan');
    % PEN_test.fA_PC1_sort_68=mean(nonzeros(fA_PC1_sort_68));
    % PEN_test.fA_PC1_sort_50=mean(nonzeros(fA_PC1_sort_50));
    % PEN_test.fA_PC1_sort_CM=mean(nonzeros(fA_PC1_sort_CM));
    %
    % PEN_test.ampXout_A_PC=mean(nonzeros(ampXout_A_PC_sum_norm_d));
    % PEN_test.out_lev_A_PC=mean(nonzeros(out_lev_A_PC_sum_norm_d));
    % PEN_test.dfA_PC_sum_68=mean(nonzeros(dfA_PC_sum_68));
    % PEN_test.dfA_PC_sum_50=mean(dfA_PC_sum_50,'omitnan');
    % PEN_test.fA_PC_sum_sort_68=mean(nonzeros(fA_PC_sum_sort_68));
    % PEN_test.fA_PC_sum_sort_50=mean(nonzeros(fA_PC_sum_sort_50));
    % PEN_test.fA_PC_sum_sort_CM=mean(nonzeros(fA_PC_sum_sort_CM));

    PEN_test.ampXout_A_PC_norm=mean(nonzeros(ampXout_A_PC_norm_d));
    PEN_test.out_lev_A_PC_norm=mean(nonzeros(out_lev_A_PC_norm_d));
    PEN_test.dfA_PC_68=mean(nonzeros(dfA_PC_68));
    PEN_test.dfA_PC_50=mean(dfA_PC_50,'omitnan');
    PEN_test.fA_PC_sort_68=mean(nonzeros(fA_PC_sort_68));
    PEN_test.fA_PC_sort_50=mean(nonzeros(fA_PC_sort_50));
    PEN_test.fA_PC_sort_CM=mean(nonzeros(fA_PC_sort_CM));

    PEN_test.ampXout_G_dom=mean(nonzeros(ampXout_G_dom_d));
    PEN_test.out_lev_G_dom=mean(nonzeros(out_lev_G_dom_d));
    PEN_test.dfG_dom_68=mean(nonzeros(dfG_dom_68));
    PEN_test.dfG_dom_50=mean(dfG_dom_50,'omitnan');
    PEN_test.fG_dom_sort_68=mean(nonzeros(fG_dom_sort_68));
    PEN_test.fG_dom_sort_50=mean(nonzeros(fG_dom_sort_50));
    PEN_test.fG_dom_sort_CM=mean(nonzeros(fG_dom_sort_CM));

    % PEN_test.ampXout_G_PC1=mean(nonzeros(ampXout_G_PC1_norm_d));
    % PEN_test.out_lev_G_PC1=mean(nonzeros(out_lev_G_PC1_norm_d));
    % PEN_test.dfG_PC1_68=mean(nonzeros(dfG_PC1_68));
    % PEN_test.dfG_PC1_50=mean(dfG_PC1_50,'omitnan');
    % PEN_test.fG_PC1_sort_68=mean(nonzeros(fG_PC1_sort_68));
    % PEN_test.fG_PC1_sort_50=mean(nonzeros(fG_PC1_sort_50));
    % PEN_test.fG_PC1_sort_CM=mean(nonzeros(fG_PC1_sort_CM));
    %
    % PEN_test.ampXout_G_PC_sum=mean(nonzeros(ampXout_G_PC_sum_norm_d));
    % PEN_test.out_lev_G_PC_sum=mean(nonzeros(out_lev_G_PC_sum_norm_d));
    % PEN_test.dfG_PC_sum_68=mean(nonzeros(dfG_PC_sum_68));
    % PEN_test.dfG_PC_sum_50=mean(dfG_PC_sum_50,'omitnan');
    % PEN_test.fG_PC_sum_sort_68=mean(nonzeros(fG_PC_sum_sort_68));
    % PEN_test.fG_PC_sum_sort_50=mean(nonzeros(fG_PC_sum_sort_50));
    % PEN_test.fG_PC_sum_sort_CM=mean(nonzeros(fG_PC_sum_sort_CM));

    PEN_test.ampXout_G_PC_norm=mean(nonzeros(ampXout_G_PC_norm_d));
    PEN_test.out_lev_G_PC_norm=mean(nonzeros(out_lev_G_PC_norm_d));
    PEN_test.dfG_PC_68=mean(nonzeros(dfG_PC_68));
    PEN_test.dfG_PC_50=mean(dfG_PC_50,'omitnan');
    PEN_test.fG_PC_sort_68=mean(nonzeros(fG_PC_sort_68));
    PEN_test.fG_PC_sort_50=mean(nonzeros(fG_PC_sort_50));
    PEN_test.fG_PC_sort_CM=mean(nonzeros(fG_PC_sort_CM));

    %Asselborn Acc
    PEN_test.bw_A3D = psd_bw_A3D;
    PEN_test.bw_A_PC1 = psd_bw_A_PC1;
    PEN_test.bw_A_PC_sum = psd_bw_A_PC_sum;
    PEN_test.bw_A_PC = psd_bw_A_PC;
    PEN_test.bw_A_tr_dominant = psd_bw_A_tr_dominant;
    PEN_test.bw_A_median = median(psd_bw_A);
    PEN_test.bw_A_tr_median = median(psd_bw_A_tr);
    %%%new
    PEN_test.bw_A3D_new = psd_bw_A3D_new;
    PEN_test.bw_A_PC1_new = psd_bw_A_PC1_new;
    PEN_test.bw_A_PC_sum_new = psd_bw_A_PC_sum_new;
    PEN_test.bw_A_PC_new = psd_bw_A_PC_new;
    PEN_test.bw_A_tr_dominant_new = psd_bw_A_tr_dominant_new;
    PEN_test.bw_A_median_new = median(psd_bw_A_new);
    PEN_test.bw_A_tr_median_new = median(psd_bw_A_tr_new);

    PEN_test.median_A3D = psd_median_A3D;
    PEN_test.median_A_PC1 = psd_median_A_PC1;
    PEN_test.median_A_PC_sum = psd_median_A_PC_sum;
    PEN_test.median_A_PC = psd_median_A_PC;
    PEN_test.median_A_tr_dominant = psd_median_A_tr_dominant;
    PEN_test.median_A_median = median(psd_median_A);
    PEN_test.median_A_tr_median = median(psd_median_A_tr);
    %%%new
    PEN_test.median_A3D_new = psd_median_A3D_new;
    PEN_test.median_A_PC1_new = psd_median_A_PC1_new;
    PEN_test.median_A_PC_sum_new = psd_median_A_PC_sum_new;
    PEN_test.median_A_PC_new = psd_median_A_PC_new;
    PEN_test.median_A_tr_dominant_new = psd_median_A_tr_dominant_new;
    PEN_test.median_A_median_new = median(psd_median_A_tr_new);
    PEN_test.median_A_tr_median_new = median(psd_median_A_tr_new);

    %Asselborn Gyro
    PEN_test.bw_G3D = psd_bw_G3D;
    PEN_test.bw_G_PC1 = psd_bw_G_PC1;
    PEN_test.bw_G_PC_sum = psd_bw_G_PC_sum;
    PEN_test.bw_G_PC = psd_bw_G_PC;
    PEN_test.bw_G_tr_dominant = psd_bw_G_tr_dominant;
    PEN_test.bw_G_median = median(psd_bw_G);
    PEN_test.bw_G_tr_median = median(psd_bw_G_m);
    %%% new
    PEN_test.bw_G3D_new = psd_bw_G3D_new;
    PEN_test.bw_G_PC1_new = psd_bw_G_PC1_new;
    PEN_test.bw_G_PC_sum_new = psd_bw_G_PC_sum_new;
    PEN_test.bw_G_PC_new = psd_bw_G_PC_new;
    PEN_test.bw_G_tr_dominant_new = psd_bw_G_tr_dominant_new;
    PEN_test.bw_G_median_new = median(psd_bw_G_new);
    PEN_test.bw_G_tr_median_new = median(psd_bw_G_m_new);

    PEN_test.median_G3D = psd_median_G3D;
    PEN_test.median_G_PC1 = psd_median_G_PC1;
    PEN_test.median_G_PC_sum = psd_median_G_PC_sum;
    PEN_test.median_G_PC = psd_median_G_PC;
    PEN_test.median_G_tr_dominant = psd_median_G_tr_dominant;
    PEN_test.median_G_median = median(psd_median_G);
    PEN_test.median_G_tr_median = median(psd_median_G_m);
    %%%new
    PEN_test.median_G3D_new = psd_median_G3D_new;
    PEN_test.median_G_PC1_new = psd_median_G_PC1_new;
    PEN_test.median_G_PC_sum_new = psd_median_G_PC_sum_new;
    PEN_test.median_G_PC_new = psd_median_G_PC_new;
    PEN_test.median_G_tr_dominant_new = psd_median_G_tr_dominant_new;
    % PEN_test.median_G_median_new = median(psd_median_G_new);
    % PEN_test.median_G_tr_median_new = median(psd_median_G_m_new);

    %%%TILT
    PEN_test.bw_TILT=psd_bw_tilt;
    PEN_test.median_TILT=psd_median_tilt_d;

    %MISURE FRA
    % PEN_test.RPW_Fra_max = max(RPW_Fra_m);
    PEN_test.RPW_Fra_median = median(RPW_Fra_m);
    PEN_test.RPW_Fra_median_new = median(RPW_Fra_m_new);
    % PEN_test.RPW_Fra_Sara_max = mean(max(RPW_Fra, [], 2),'omitnan');
    PEN_test.RPW_Fra_Sara_mean = mean(mean(RPW_Fra,2, 'omitnan'),1, 'omitnan');
    PEN_test.RPW_Fra_PC = mean(nonzeros(RPW_Fra_G_PC));
    % PEN_test.RPW_Fra_PC1 = mean(nonzeros(RPW_Fra_G_PC1));
    % PEN_test.RPW_Fra_PC_sum = mean(nonzeros(RPW_Fra_G_PC_sum));
end