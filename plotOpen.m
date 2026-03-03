clear all;
clc;
%cd_folder = '/Users/simonetoffoli/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennali/TesiTriennaleAnziani/GraficiSoggetti';
%cd_folder = '/Users/micol/Library/CloudStorage/OneDrive-PolitecnicodiMilano/Simone Toffoli/TesiTriennali/TesiTriennaleAnziani/GraficiSoggetti2';
cd_folder = '\Users\alber\OneDrive - Politecnico di Milano\TesiTriennaleAnziani\GraficiSoggetti2';
%cd_folder='/Users/mariavittoriasari/Library/CloudStorage/OneDrive-PolitecnicodiMilano/TesiTriennaleAnziani/GraficiSoggetti2'; 
cd(cd_folder);

%capitals = [cd_folder '/GraficiSentenceCapitals'];
%cursive = [cd_folder '/GraficiSentenceCursive'];
spiral = [cd_folder '\GraficiSpiral'];

%%
close all
soggetto = '3';

%openfig([ capitals '/soggetto_' soggetto '__gruppo1.fig']); 
%openfig([ capitals '/soggetto_' soggetto '__gruppo2.fig']); 

%openfig([ cursive '\soggetto_' soggetto '__gruppo1.fig']); 
%openfig([ cursive '\soggetto_' soggetto '__gruppo2.fig']); 

openfig([ spiral '/soggetto_' soggetto '__gruppo1.fig']); 
openfig([ spiral '/soggetto_' soggetto '__gruppo2.fig']); 



