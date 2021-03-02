%iDN1004 Script (Created June 12th, 2020) - Daniel Norena Caro
%iDN1004 Script (Edited February 24th, 2021) - Daniel Norena Caro
%Genome scale metabolic model of Anabaena sp. UTEX 2576 
%a.k.a. Nostoc sp. PCC7120

clc
clear

%initCobraToolbox

%Changes default solver to Gurobi
changeCobraSolver ('gurobi', 'all', 1)

%Transforms model from Excel file to matlab model
%Make sure the xls file is stored in the same directory
model = readCbModel('Data_set_1_iDN1004.xls');

%Load model from mat iDN1004.mat file
%load Dataset_2_iDN1004_PAuto.mat or Dataset_3_iDN1004_PDiazo.mat

%Prints the stoiquiometric matrix for the model

%To extract the S matrix from the model
mat = full(model.S); 

%Saves matrix as csv
%csvwrite('stoiqiDN1004.csv',mat) 

%Identify metabolic Dead ends

outputMets = detectDeadEnds(model);
DeadEnds = model.mets(outputMets);

answer = questdlg('What is the nitrogen source?', ...
	'Growth conditions', ...
	'NaNO3 (Photoautotrophic)','N2 (Photodiazotrophic)','Abort');
% Handle response
switch answer
    case 'NaNO3 (Photoautotrophic)'
    %Change objective function (Non-diazotrophic)
    model = changeObjective(model,'BOF_ANA_AUTO',1);


    fprintf ('The Objective reaction is: \n')
    %Print  Objective function
    printRxnFormula(model,'BOF_ANA_AUTO');    
    
    %Change reation limits for cultures in BG11 (Photoautotrophic)

    %Photon flux based on experimental PPFD
    model = changeRxnBounds(model, 'EX_photon_e', -87.09, 'l');
    %Citrate flux from BG11 medium
    model = changeRxnBounds(model, 'EX_cit_e', -0.01, 'b');
    %Pyruvate from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_pyr_e', 0, 'b');
    %L-Methionine flux from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_met__L_e', 0, 'b');
    %L-Glutamate flux from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_gln__L_e', 0, 'b');
    %L-Glutamine flux from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_glu__L_e', 0, 'b');
    %L-Arginine flux from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_arg__L_e', 0, 'b');
    %L-Lysine flux from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_lys__L_e', 0, 'b');
    %L-Histidine flux from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_his__L_e', 0, 'b');
    %Urea flux from medium blocked, but excretion allowed
    model = changeRxnBounds(model, 'EX_urea_e', 0, 'l');
    %Putrescine flux from medium blocked, but excretion allowed
    model = changeRxnBounds(model, 'EX_ptrc_e', 0, 'l');
    %Spermidine flux from medium blocked, but excretion allowed
    model = changeRxnBounds(model, 'EX_spmd_e', 0, 'l');
    %Cyanate flux from medium blocked
    model = changeRxnBounds(model, 'EX_cynt_e', 0, 'b');
    %Ammonium ion flux/uptake from medium
    model = changeRxnBounds(model, 'EX_nh4_e', 0, 'l');
    %model = changeRxnBounds(model, 'EX_nh4_e', 0, 'u');
    %Diazotrophy is blocked in BG11 medium (No N2 uptake)
    model = changeRxnBounds(model, 'EX_n2_e', 0, 'b');
    %Nitrate ion flux/uptake from medium
    model = changeRxnBounds(model, 'EX_no3_e', -0.2438, 'l');
    model = changeRxnBounds(model, 'EX_no3_e', -0.0129, 'u');
    %Biosynthesis of phycocyanobilin pigment
    %model = changeRxnBounds(model, 'PHYCYFX', 2.27E-5, 'l');
    %model = changeRxnBounds(model, 'PHYCYFX', 5.70E-4, 'u');
    %Biosynthesis of c-phycocyanin phycobiliprotein
    model = changeRxnBounds(model, 'CPCPBPS', 5.84E-6, 'l');
    model = changeRxnBounds(model, 'CPCPBPS', 1.47E-4, 'u');
    %Uptake of Vitamin B12/ Calomide from medium blocked
    model = changeRxnBounds(model, 'EX_adocbl_e', 0, 'l');
    %Uptake of Vitamin B2/ Riboflavin from medium blocked
    model = changeRxnBounds(model, 'EX_ribflv_e', 0, 'l');
    %Uptake of Glycerol from medium blocked
    model = changeRxnBounds(model, 'EX_glyc_e', 0, 'b');
    %Constrained HCO3 and CO2 uptake reaction
    %HCO3 uptake is 70% of total Ci required
    %CO2 uptake is 30% of total Ci required
    model = changeRxnBounds(model, 'EX_co2_e', -0.3561, 'b');
    model = changeRxnBounds(model, 'EX_hco3_e', -0.8309, 'b');
    %Constrained secretion of Schizokinen
    model = changeRxnBounds(model, 'SCHIZt', 0.000491, 'l');
    %Fe 2+ avilability constraint (Only Fe 3+ in medium)
    model = changeRxnBounds(model, 'EX_fe2_e', 0, 'b');
     
    
    %Run FBA
    FBAsolution = optimizeCbModel (model, 'max','one');
    obj_f_value = FBAsolution.f;%Objective function

    [{'Reaction ID', 'Objective', 'Flux'};...
    model.rxns, num2cell(model.c), num2cell(FBAsolution.x)]
      
        
    case 'N2 (Photodiazotrophic)'
        
    %Change objective function (Diazotrophic)
    model = changeObjective(model,'BOF_ANA_DIAZO',1);
        


    fprintf ('The Objective reaction is: \n')
    %Print  Objective function
    printRxnFormula(model,'BOF_ANA_DIAZO');
    
     
   
    %Change reation limits for cultures in BG11o (Diazotrophic)

    
    %Photon flux based on experimental PPFD
    model = changeRxnBounds(model, 'EX_photon_e', -102.53, 'l');
    
    %PSI activity constraint for minimum growth of heterocyst
    %model = changeRxnBounds(model, 'PHOt_d', 0.316, 'b');
    %model = changeRxnBounds(model, 'PHOt_e', 19.184643, 'u');
   
        
    %Citrate flux from BG11o medium
    model = changeRxnBounds(model, 'EX_cit_e', -0.01, 'b');
    model = changeRxnBounds(model, 'CITt', 0.0085, 'b');
    model = changeRxnBounds(model, 'CITt_d', 0.0015, 'b');
    %Pyruvate from BG11o medium blocked
    model = changeRxnBounds(model, 'EX_pyr_e', 0, 'b');
    %L-Methionine flux from BG11o medium blocked
    model = changeRxnBounds(model, 'EX_met__L_e', 0, 'b');
    %L-Glutamine flux from BG11 medium blocked
    model = changeRxnBounds(model, 'EX_glu__L_e', 0, 'b');
    %L-Glutamine flux from BG11o medium blocked
    model = changeRxnBounds(model, 'EX_gln__L_e', 0, 'b');
    %L-Arginine flux from BG11o medium blocked
    model = changeRxnBounds(model, 'EX_arg__L_e', 0, 'b');
    %L-Lysine flux from BG11o medium blocked
    model = changeRxnBounds(model, 'EX_lys__L_e', 0, 'b');
    %L-Histidine flux from BG11o medium blocked
    model = changeRxnBounds(model, 'EX_his__L_e', 0, 'b');
    %Urea flux from medium blocked, but excretion allowed
    model = changeRxnBounds(model, 'EX_urea_e', 0, 'l');
    %Putrescine flux from medium blocked, but excretion allowed
    model = changeRxnBounds(model, 'EX_ptrc_e', 0, 'l');
    %Spermidine flux from medium blocked, but excretion allowed
    model = changeRxnBounds(model, 'EX_spmd_e', 0, 'l');
    %Cyanate flux from medium blocked
    model = changeRxnBounds(model, 'EX_cynt_e', 0, 'b');
    %Ammonium ion flux/uptake from medium
    model = changeRxnBounds(model, 'EX_nh4_e', 0, 'l');
    %model = changeRxnBounds(model, 'EX_nh4_e', 0, 'u');
    %Diazotrophy is allowed in BG11o medium (N2 uptake)
    model = changeRxnBounds(model, 'EX_n2_e', -0.1255, 'l');
    model = changeRxnBounds(model, 'EX_n2_e', -0.0104, 'u');
    %BG11o medium lacks Nitrate, no nitrate uptake
    model = changeRxnBounds(model, 'EX_no3_e', 0, 'b');
    %Biosynthesis of phycocyanobilin pigment
    %model = changeRxnBounds(model, 'PHYCYFX', 2.01E-5, 'l');
    %model = changeRxnBounds(model, 'PHYCYFX', 4.05E-4, 'u');
    %Biosynthesis of c-phycocyanin phycobiliprotein
    model = changeRxnBounds(model, 'CPCPBPS', 4.35E-6, 'l');
    model = changeRxnBounds(model, 'CPCPBPS', 8.77E-5, 'u');
    model = changeRxnBounds(model, 'CPCPBPS_d', 7.67E-7, 'l');
    model = changeRxnBounds(model, 'CPCPBPS_d', 1.55E-5, 'u');
    %Uptake of Vitamin B12/ Calomide from medium blocked
    model = changeRxnBounds(model, 'EX_adocbl_e', 0, 'l');
    %Uptake of Vitamin B2/ Riboflavin from medium blocked
    model = changeRxnBounds(model, 'EX_ribflv_e', 0, 'l');
    %Uptake of Glycerol from medium blocked
    model = changeRxnBounds(model, 'EX_glyc_e', 0, 'b');
    %Constrained CO2 uptake reaction
    %HCO3 uptake is 70% of total Ci required
    %CO2 uptake is 30% of total Ci required
    model = changeRxnBounds(model, 'EX_co2_e', -0.3249, 'b');
    model = changeRxnBounds(model, 'EX_hco3_e', -0.7581, 'b');

    %Constrained secretion of Schizokinen
    model = changeRxnBounds(model, 'SCHIZt', 0.000371705, 'l');
    model = changeRxnBounds(model, 'SCHIZt_d', 0.000065595, 'l');
   
    %Fe 2+ avilability constraint (Only Fe 3+ in medium)
    model = changeRxnBounds(model, 'EX_fe2_e', 0, 'b');
        
      
    %Run FBA
    FBAsolution = optimizeCbModel (model, 'max','one');

    %printFluxVector(model,FBAsolution.x,true)
    obj_f_value = FBAsolution.f;%Objective function

    [{'Reaction ID', 'Objective', 'Flux'};...
    model.rxns, num2cell(model.c), num2cell(FBAsolution.x)]
        
        
    case 'Abort'
        disp('Operation aborted')
       
end

%Check Stoichiometry matrix sparsity
S = model.S;
[nMets, nRxns] = size(S);
nElem = numel(S);
nNz = nnz(S);
sparsityRatio = (1-nNz/nElem)*100;
compSparsityRatio = 100 - sparsityRatio;

colDensityAv = 0;

for j=1:nRxns
    colDensityAv = colDensityAv + nnz(S(:,j));
end

colDensityAv = colDensityAv/ nRxns;

colDensityRel = (colDensityAv/nMets)*100;

spyc(S, colormap(advancedColormap('cobratoolbox')));
set(gca, 'fontsize',14);

%Check Mass and charge balance of reactions

[~,imBalancedMass,imBalancedCharge,imBalancedBool,~] = checkMassChargeBalance(model,1);

%writeCbModel(model)





