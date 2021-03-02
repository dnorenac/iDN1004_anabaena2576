%iDN1004 Script FVA (December 13th, 2020) - Daniel N. Caro
%Genome scale metabolic model of Anabaena sp. UTEX 2576 
%a.k.a. Nostoc sp. PCC7120

clc
clear

%initCobraToolbox

%Changes default solver to Gurobi
changeCobraSolver ('gurobi', 'all', 1)

rxnlist = {'ABTA';'ABTA_d';'ADCS';'ADCS_d';'ADCYRS';'ADCYRS_d';'AGTi';'AGTi_d';
    'ALAR';'ALAR_d';'AMPMS3';'AMPMS3_d';'AMPTASECG';'AMPTASECG_d';'ANS';'ANS_d';
    'ARGN';'ARGN_d';'ARGSL';'ARGSL_d';'ARODH';'ARODH_2';'ARODH_2_d';'ARODH_d';
    'AROH';'AROH_d';'ASNN';'ASNN_d';'ASNS1';'ASNS1_d';'ASPCT';'ASPK';'ASPK_d';
    'ASPTA';'ASPTA_d';'ASPTA4';'ASPTA4_d';'BTS6';'BTS6_d';'CBPS';'CBPS_d';'COBNAD';
    'COBNAD_d';'CPPPGO2';'CPPPGO2_d';'CTPS2';'CTPS2_d';'CYPH_BIA';'CYPH_BIA_d';
    'CYSS_2';'CYSS_2_d';'CYSTS';'CYSTS_d';'DAPDC';'DAPDC_d';'DMHDRFS_1';
    'DMHDRFS_1_d';'GDH1';'GDH2_nadp';'GF6PTA';'GF6PTA_d';'GHMT2r';'GHMT2r_d';
    'GLMS_b';'GLNS';'GLNS_d';'GLUN';'GLUPRT';'GLUPRT_d';'GLYTA';'GLYTA_d';
    'GOR1';'GOR1_d';'GTHRDH_syn';'GTHRDH_syn_d';'GTPC';'GTPC_d';'HISTDb';
    'HISTDb_d';'IG3PS_1';'IG3PS_1_d';'ILEDH_nad';'ILEDH_nad_d';'LLEUDr';
    'LLEUDr_d';'METAT';'METAT_d';'METS_1';'METS_1_d';'NACASPAH';'NACASPAH_d';
    'NADS2';'NADS2_d';'OPAH';'OPAH_d';'ORNTA_d';'P5CD';'P5CD_d';'P5CR';'P5CR_d';
    'P5CRx';'P5CRx_d';'PHETA1';'PHETA1_d';'PPHTA';'PPHTA_d';'PRFGS';'PRFGS_d';
    'PSP_L';'PSP_L_d';'R05224';'R05224_d';'SCYSDS';'SCYSDS_d';'SPT_syn';'SPT_syn_d';
    'THRA';'THRA_d';'THRS';'THRS_d';'THZPSN';'TRPS2';'TRPS2_d';'TYRTA';'TYRTA_d';
    'UNK3';'UNK3_d';'VALDHr';'VALDHr_d';'VPAMTr';'VPAMTr_d'};


for i=1:131;

load('Data_set_2_iDN1004_PAuto.mat');
modelfva1= model;
load('Data_set_3_iDN1004_PDiazo.mat');
modelfva2= model;

modelfva1 = changeObjective(modelfva1,'BOF_ANA_AUTO',0);
modelfva1 = changeObjective(modelfva1,rxnlist(i),1);

modelfva2 = changeObjective(modelfva2,'BOF_ANA_DIAZO',0);
modelfva2 = changeObjective(modelfva2,rxnlist(i),1);

gamma = 0.95;
gr_pa = 0.0283;
gr_pd = 0.0258;

gr_pa_2 = gamma*gr_pa;
gr_pd_2 = gamma*gr_pd;

%Change growth constraint to at least gamma=95% of FBA solution
%Change gamma as needed when solution is not feasible
modelfva1 = changeRxnBounds(modelfva1, 'BOF_ANA_AUTO', gr_pa_2, 'l');
modelfva2 = changeRxnBounds(modelfva2, 'BOF_ANA_DIAZO', gr_pd_2, 'l');

%Run FBA (Min/Max reaction)
FVAsolution1 = optimizeCbModel (modelfva1, 'max','one');
FVAsolution2 = optimizeCbModel (modelfva1, 'min','one');
FVAsolution3 = optimizeCbModel (modelfva2, 'max','one');
FVAsolution4 = optimizeCbModel (modelfva2, 'min','one');

PAdataMax(:,i)=FVAsolution1.x;
PAdataMin(:,i)=FVAsolution2.x;
PDdataMax(:,i)=FVAsolution3.x;
PDdataMin(:,i)=FVAsolution4.x;

end
