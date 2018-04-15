%  Multiple Pools of Ser and Gly
%  Disease: Lung Cancer
%  Species: Human
%  Organ & Tissue: Lung, Lung Tissue (Lung Cancer Patient Lung Tissue)
%  Cancer tissue & adjacent non-cancer tissue
%  Hypothesis: Multiple Pools of Ser and Gly  vs.  One pool of Ser and Gly (in Cytoplasm)
%  Other Ser and Gly pools: extracellular (in Medium)  &  Mitochon
%  Data: NMR & MS data
%  Time resolution:  0 h  &  24 h

function MultiPoolsSerGly();

%  Creat Variables 
%  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
%     Compute number of molecules for a Variable (metabolite with all isomers) from its initial concentrations

CellVolume = 0.5E-6;      %  u Liter
MediumVolume = 7.8 * 10^3;       % u Liter		7.8 mL
TissueVolume = 26;    % 0.026 mL = 26 u Liter
CytoplasmVolume = TissueVolume * 0.64;     % 0.64 * 26 u Liter
MitochonVolume = TissueVolume * 0.16;     % 0.16 * 26 u Liter

    CO2 = 1;
    NH3 = 1;
    cCO2 = 1;
    cNH3 = 1;
    mCO2 = 1;
    mNH3 = 1;

    AvogadroConstant = 6.02214129E+23;

% Avogadro Constant * X (pL) = 6.02214129 * 10^11 * X






Kf_GLC_Transport = 0.00073 * 0.3; 				% / min
Kr_GLC_Transport = 0.288 * 0.3; 				% / min

Kf_GLC_SerPool1 = 0.0011 * 0.4;					% / min
Kf_GLC_PYR = 0.1833;					% / min
Kf_GLC_PRPP = 0.000083 * 1.5;   				% / min

Kf_SerPool3_Synthesis = 6.4 * 10^-6 * 1.2;			% M /min 
Kf_SerPool3_Degradation = 1.77 * 10^-2;	 				% /min

Kf_GlyPool3_Synthesis = 6.4 * 10^-6 * 3; 				% M /min 
Kf_GlyPool3_Degradation = 0.53 * 10^-3;    			% /min

Kf_SerPool1_GlyPool1 = 176 * 1.5;     	% /M /min
Kr_SerPool1_GlyPool1 = 3.55 * 0.75; 		% /M /min

Kf_Ser_Transport = 0.000139 * 1.5;			% /min
Kr_Ser_Transport = 0.00667 * 1.5;				% /min

Kf_Gly_Transport = 0.0000603 * 1.5;			% /min
Kr_Gly_Transport = 0.000694 * 1.5;			% /min

Kf_SerPool2_SerPoolMitochon = 0.0013 * 3;				% /min
Kr_SerPool2_SerPoolMitochon = 0.00325 * 3;			% /min

Kf_GlyPool2_GlyPoolMitochon = 0.00195 * 3;			% /min 
Kr_GlyPool2_GlyPoolMitochon = 0.004875 * 3;			% /min 

Kf_SerPoolMitochon_GlyPoolMitochon = 277 * 3;		% /M /min
Kr_SerPoolMitochon_GlyPoolMitochon = 5.5 * 3;		% /M /min

Kf_GlyPool1_CO2 = 17.7 * 3;						% /M /min
Kr_GlyPool1_CO2 = 0.213 * 3;						% /min

Kf_MethyleneTHF_FormylTHF_Cytoplasm = 0.32;			% /min
Kr_MethyleneTHF_FormylTHF_Cytoplasm = 0.8;			% /min 

Kf_FormylTHF_Formate_Cytoplasm = 0.32;						% /min 
Kr_FormylTHF_Formate_Cytoplasm = 1.0 * 10^6; 			% /M /min 

Kf_GlyPoolMitochon_CO2 = 434;						% /M /min
Kr_GlyPoolMitochon_CO2 = 0.42; 					% /min

Kf_MethyleneTHF_FormylTHF_Mitochon = 2.67;					% /min
Kr_MethyleneTHF_FormylTHF_Mitochon = 0.267;					% /min

Kf_FormylTHF_Formate_Mitochon = 2.67;							% /min
Kr_FormylTHF_Formate_Mitochon = 4.44 * 10^4;			% /M /min

Kf_PRPP1_GAR = 44 * 2; 			% /M /min

Kf_GAR_FGAR = 88 * 10^3 * 0.5; 			% /M /min

Kf_FGAR_AMP = 88 * 10^3 * 0.5; 			% /M /min

Kf_IMP_AMP = 0.025 * 0.5;		% /min

Kf_AMP_Degradation = 0.00125 * 0.5; 		% /min

Kr_MethyleneTHF_Transport = 0.013 * 3; 			% /min
Kf_MethyleneTHF_Transport = 0.026 * 0.3; 		% /min

Kr_FormylTHF_Transport = 2.34 * 2; 			% /min
Kf_FormylTHF_Transport = 0.468; 				% /min

Kr_THF_Transport = 0.78; 						% /min
Kf_THF_Transport = 0.156 * 0.5; 		% /min

Kr_Formate_Transport = 4.68; 						% /min
Kf_Formate_Transport = 0.936 * 0.5; 		% /min

Kf_PRPP2_GAR = 22 * 2;  			% /M /min
Kf_PRPP3_GAR = 11 * 2;  			% /M /min


   Kf_SerPool2_GlyPool2 = 176 * 1.5;    		% /M /min
   Kr_SerPool2_GlyPool2 = 3.55 * 0.75;    	% /M /min
   Kf_GlyPool2_CO2 = 17.7 * 3;				% /M /min
   Kr_GlyPool2_CO2 = 0.213 * 3; 			% /min

	Kf_SerPool3_GlyPool3 = Kf_SerPool2_GlyPool2;
	Kr_SerPool3_GlyPool3 = Kr_SerPool2_GlyPool2;
	Kf_GlyPool3_CO2 = Kf_GlyPool2_CO2;
	Kr_GlyPool3_CO2 = Kr_GlyPool2_CO2;
	
	Kf_Ser_Transport3 = Kf_Ser_Transport;
	Kr_Ser_Transport3 = Kr_Ser_Transport;
	Kf_Gly_Transport3 = Kf_Gly_Transport;
	Kr_Gly_Transport3 = Kr_Gly_Transport;





OptimizeParameters = [0.000461657750477753	0.427503086639990	0.000138933974778344	0.912957765988904	0.000302264234311588	0.000264216988034605	0.528433976069209	0.000178629396143585	0.0105686795213842	389.997313257504	3.65988124056536	0.000422747180855367	0.0496192767065513	7.44289150598270e-05	0.0178629396143585	0.0109460124414652	0.0140352064043982	0.00966781987350446	0.0490562874450916	830.998723620212	17.0631986313531	54.1720082557871	1.25536770067575	1.67337425755250	2.03451173487994	0.399124270787164	1000000.00000000	434.000150742915	1.64330339153493	2.22492836752176	0.600557183706194	2.29538774044506	44400.0000000432	89.9889179513641	44000.0000016264	43999.9999979280	0.0266554754467594	0.00151731008995339	0.189003218774087	0.00751136610783774	11.6715973413360	0.495498097191622	2.51038222698186	0.0677700081258078	4.57755449879482	0.579744450285708	42.5463885164964	21.0649543138687	85.0037419701686	120.001190410154	65.0980133729533	1.78629396143585	CellVolume	MediumVolume	TissueVolume	CytoplasmVolume	MitochonVolume];

%noncancer
%OptimizeParameters = [0.000465165222015355  0.242665045280639    0.000036     0.518326766571717        0.0002   0.00003    0.2         0.000075             0.075     130.0               0.6         0.00018               0.075     0.00005               0.12       0.001     0.0008               0.00974237918290052   0.0278492039049790      780.0     0.6         54.2492566486518               0.03       0.95       2.00461377104628          0.398825587055585        999999.999999978               434.005733100667          1.05698640506409          2.24239831194499               0.551424353965021        2.31291778438686          44400.0000092611          80.0               44000.0001051987          44000.0000615634          0.0268575860192817               0.000861408841917836                0.107341967375275        0.00756915225856982               11.6618608584811          0.499324532242705        2.01627621832653               0.0682897697038028      4.57835999241956          0.579450741578679        42.0       21.0               120.0     1.2  65.0  0.72  CellVolume  MediumVolume  TissueVolume     CytoplasmVolume      MitochonVolume];

%  For a 1 u mole substrate, it has 1 u moles = 1 * 10^-6 * 6.02214129*10^23 = 1 * 6.02214129*10^17 molecules

NumberOfIsotopomers = 608;
Xd = zeros(NumberOfIsotopomers, 1);


cGLC_13C0 = 0.3 * 6.02214129 * 10^17;      % 300 n moles unlabeled cGLC
cSerPool1_13C000D000 = 0.01 * 6.02214129 * 10^17;      % 0.01 u moles unlabeled cSerPool1
cPRPP_13C0 = 0.0001 * 6.02214129 * 10^17;      % 0.1 n moles unlabeled cPRPP
cGlyPool1_13C00D00 = 0.1 * 6.02214129 * 10^17;      % 0.1 u moles unlabeled cGlyPool1
cTHF = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled cTHF
cMethyleneTHF_13C0D00 = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled cMethyleneTHF (no carbon label, no deterium label)
cSerPool2_13C000D000 = 0.01 * 6.02214129 * 10^17;      % 0.01 u moles unlabeled cSerPool2 (no carbon label, no deterium label)
cGlyPool2_13C00D00 = 0.1 * 6.02214129 * 10^17;      % 0.1 u moles unlabeled cGlyPool2 (no carbon label, no deterium label)
cSerPool3_13C000D000 = 0.01 * 6.02214129 * 10^17;      % 0.01 u moles unlabeled cSerPool3 (no carbon label, no deterium label)
cGlyPool3_13C00D00 = 0.1 * 6.02214129 * 10^17;      % 0.1 u moles unlabeled cGlyPool3 (no carbon label, no deterium label)
mSerPool_13C000D000 = 0.002 * 6.02214129 * 10^17;      % 0.002 u moles unlabeled mSerPool (no carbon label, no deterium label)
mGlyPool_13C00D00 = 0.02 * 6.02214129 * 10^17;      % 0.02 u moles unlabeled mGlyPool (no carbon label, no deterium label)
mTHF = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled cTHF
mMethyleneTHF_13C0D00 = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled mMethyleneTHF (no carbon label, no deterium label)
cFormylTHF_13C0D0 = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled cFormylTHF (no carbon label, no deterium label)
cFormate_13C0D0 = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled cFormate (no carbon label, no deterium label)
mFormate_13C0D0 = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled mFormate (no carbon label, no deterium label)
mFormylTHF_13C0D0 = 0.00005 * 6.02214129 * 10^17;      % 0.05 n moles unlabeled mFormylTHF (no carbon label, no deterium label)
cGAR_13C000 = 0.0001 * 6.02214129 * 10^17;      % 0.1 n moles unlabeled cGAR (no carbon label, no deterium label)
cFGAR_13C0000D0 = 0.0001 * 6.02214129 * 10^17;      % 0.1 n moles unlabeled cFGAR (no carbon label, no deterium label)
cAMP_13C00000D00 = 0.02 * 6.02214129 * 10^17;      % 20 n moles unlabeled cAMP (no carbon label, no deterium label)


Xd(3) = cGLC_13C0;
Xd(85) = cSerPool1_13C000D000;
Xd(439) = cPRPP_13C0;
Xd(149) = cGlyPool1_13C00D00;
Xd(421) = cTHF;
Xd(405) = cMethyleneTHF_13C0D00;
Xd(165) = cSerPool2_13C000D000;
Xd(229) = cGlyPool2_13C00D00;
Xd(245) = cSerPool3_13C000D000;
Xd(309) = cGlyPool3_13C00D00;
Xd(325) = mSerPool_13C000D000;
Xd(389) = mGlyPool_13C00D00;
Xd(438) = mTHF;
Xd(422) = mMethyleneTHF_13C0D00;
Xd(413) = cFormylTHF_13C0D0;
Xd(417) = cFormate_13C0D0;
Xd(434) = mFormate_13C0D0;
Xd(430) = mFormylTHF_13C0D0;
Xd(441) = cGAR_13C000;
Xd(449) = cFGAR_13C0000D0;
Xd(481) = cAMP_13C00000D00;



RunFlag = 'Cancer';				% 'NonCancer'

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Optimization - No.1:		13C6-GLC only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if strcmp(RunFlag, 'Cancer')

			eGLC_13C1 = 236.310 * 6.02214129 * 10^17;      % 236310 n moles labeled eGLC
			Xd(2) = eGLC_13C1;

			eSer_13C000D000 = 2.987748945 * 6.02214129 * 10^17;      % 2987.748945 n moles unlabeled eSer (no carbon label, no deterium label)
			Xd(5) = eSer_13C000D000;
			
			eGly_13C00D00 = 4.496768784 * 6.02214129 * 10^17;      % 4496.768784 n moles unlabeled eGly (no carbon label, no deterium label)
			Xd(69) = eGly_13C00D00;

			
		elseif strcmp(RunFlag, 'NonCancer')
		
			eGLC_13C1 = 196.292 * 6.02214129 * 10^17;      % 196292 n moles labeled eGLC
			Xd(2) = eGLC_13C1;

			eSer_13C000D000 = 3.310533725 * 6.02214129 * 10^17;      % 3310.533725 n moles unlabeled eSer (no carbon label, no deterium label)
			Xd(5) = eSer_13C000D000;
			
			eGly_13C00D00 = 3.955700675 * 6.02214129 * 10^17;      % 3955.700675 n moles unlabeled eGly (no carbon label, no deterium label)
			Xd(69) = eGly_13C00D00;
			
		end
	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%Optimization - No.1:		13C6-GLC only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%Optimization - No.2:		13C6-GLC + D3-Ser
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 		if strcmp(RunFlag, 'Cancer')
% 
% 			eGLC_13C1 = 218.872 * 6.02214129 * 10^17;      % 218872 n moles labeled eGLC
% 			Xd(2) = eGLC_13C1;
% 
% 			eSer_13C000D111 = 1.988443586 * 6.02214129 * 10^17;      % 1988.443586 n moles labeled eSer (no carbon label, D3 deterium label)
% 			Xd(61) = eSer_13C000D111;
% 			
% 			eGly_13C00D00 = 3.304771955 * 6.02214129 * 10^17;      % 3304.771955 n moles unlabeled eGly (no carbon label, no deterium label)
% 			Xd(69) = eGly_13C00D00;                       %%ODE solver failed.
% 
% 			
% 		elseif strcmp(RunFlag, 'NonCancer')
% 		
% 			eGLC_13C1 = 209.796 * 6.02214129 * 10^17;      % 209796 n moles labeled eGLC
% 			Xd(2) = eGLC_13C1;
% 
% 			eSer_13C000D111 = 2.846508241 * 6.02214129 * 10^17;      % 2846.508241 n moles labeled eSer (no carbon label, D3 deterium label)
% 			Xd(61) = eSer_13C000D111;
% 			
% 			eGly_13C00D00 = 4.070300751 * 6.02214129 * 10^17;      % 4070.300751 n moles unlabeled eGly (no carbon label, no deterium label)
% 			Xd(69) = eGly_13C00D00;
% 			
% 		end                                                %%ODE solver failed.
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%Optimization - No.2:		13C6-GLC + D3-Ser
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%Optimization - No.3:		13C6-GLC + D2-Gly
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 		if strcmp(RunFlag, 'Cancer')
% 
% 			eGLC_13C1 = 235.264 * 6.02214129 * 10^17;      % 235264 n moles labeled eGLC
% 			Xd(2) = eGLC_13C1;
% 
% 			eSer_13C000D000 = 2.99509567 * 6.02214129 * 10^17;      % 2995.09567 n moles unlabeled eSer (no carbon label, no deterium label)
% 			Xd(5) = eSer_13C000D000;
% 			
% 			eGly_13C00D11 = 3.3552663 * 6.02214129 * 10^17;      % 3355.2663 n moles labeled eGly (no carbon label, D2 deterium label)
% 			Xd(81) = eGly_13C00D11;
% 
% 			
% 		elseif strcmp(RunFlag, 'NonCancer')
% 		
% 			eGLC_13C1 = 211.237 * 6.02214129 * 10^17;      % 211237 n moles labeled eGLC
% 			Xd(2) = eGLC_13C1;
% 
% 			eSer_13C000D000 = 4.276838531 * 6.02214129 * 10^17;      % 4276.838531 n moles unlabeled eSer (no carbon label, no deterium label)
% 			Xd(5) = eSer_13C000D000;
% 			
% 			eGly_13C00D11 = 4.693109101 * 6.02214129 * 10^17;      % 4693.109101 n moles labeled eGly (no carbon label, D2 deterium label)
% 			Xd(81) = eGly_13C00D11;
% 			
% 		end
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%Optimization - No.3:		13C6-GLC + D2-Gly
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










t0 = 0 * 60;   % unit: min
tf = 24 * 60;   % unit: min
dt = 60;       % unit: min  
tspan = [t0 : dt : tf];

%options=odeset('RelTol',1e-6,'AbsTol',1e-6);




[t, Yd] = ode15s(@PurineSynthesis, tspan, Xd, [], OptimizeParameters);
    if size(Yd, 1) ~= size(tspan, 2)
%         'ODEsolver failed'
%         size(Yd)
        Yd = zeros(size(tspan, 2), size(Xd, 2));
    end




% Convert number of molecules to real unit
% Then express them with unit of nano mole (10^-9 mole = 10^-9 * 6.02214129*10^23 = 6.02214129*10^14
% ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MoleculeNumberInOneNanoMole = 6.02214129*10^14;

	Yd = Yd ./ MoleculeNumberInOneNanoMole;
% T = array2table(Yd);
% filename = 'data_1.xlsx';
% writetable(T,filename,'Sheet',1,'Range','D1');





eGLC_13C0_24h				  = Yd(end, 1  );
eGLC_13C1_24h                 = Yd(end, 2  );
cGLC_13C0_24h                 = Yd(end, 3  );
cGLC_13C1_24h                 = Yd(end, 4  );
eSer_13C000D000_24h           = Yd(end, 5  );
eSer_13C001D000_24h           = Yd(end, 6  );
eSer_13C010D000_24h           = Yd(end, 7  );
eSer_13C011D000_24h           = Yd(end, 8  );
eSer_13C100D000_24h           = Yd(end, 9  );
eSer_13C101D000_24h           = Yd(end, 10 );
eSer_13C110D000_24h           = Yd(end, 11 );
eSer_13C111D000_24h           = Yd(end, 12 );
eSer_13C000D001_24h           = Yd(end, 13 );
eSer_13C001D001_24h           = Yd(end, 14 );
eSer_13C010D001_24h           = Yd(end, 15 );
eSer_13C011D001_24h           = Yd(end, 16 );
eSer_13C100D001_24h           = Yd(end, 17 );
eSer_13C101D001_24h           = Yd(end, 18 );
eSer_13C110D001_24h           = Yd(end, 19 );
eSer_13C111D001_24h           = Yd(end, 20 );
eSer_13C000D010_24h           = Yd(end, 21 );
eSer_13C001D010_24h           = Yd(end, 22 );
eSer_13C010D010_24h           = Yd(end, 23 );
eSer_13C011D010_24h           = Yd(end, 24 );
eSer_13C100D010_24h           = Yd(end, 25 );
eSer_13C101D010_24h           = Yd(end, 26 );
eSer_13C110D010_24h           = Yd(end, 27 );
eSer_13C111D010_24h           = Yd(end, 28 );
eSer_13C000D011_24h           = Yd(end, 29 );
eSer_13C001D011_24h           = Yd(end, 30 );
eSer_13C010D011_24h           = Yd(end, 31 );
eSer_13C011D011_24h           = Yd(end, 32 );
eSer_13C100D011_24h           = Yd(end, 33 );
eSer_13C101D011_24h           = Yd(end, 34 );
eSer_13C110D011_24h           = Yd(end, 35 );
eSer_13C111D011_24h           = Yd(end, 36 );
eSer_13C000D100_24h           = Yd(end, 37 );
eSer_13C001D100_24h           = Yd(end, 38 );
eSer_13C010D100_24h           = Yd(end, 39 );
eSer_13C011D100_24h           = Yd(end, 40 );
eSer_13C100D100_24h           = Yd(end, 41 );
eSer_13C101D100_24h           = Yd(end, 42 );
eSer_13C110D100_24h           = Yd(end, 43 );
eSer_13C111D100_24h           = Yd(end, 44 );
eSer_13C000D101_24h           = Yd(end, 45 );
eSer_13C001D101_24h           = Yd(end, 46 );
eSer_13C010D101_24h           = Yd(end, 47 );
eSer_13C011D101_24h           = Yd(end, 48 );
eSer_13C100D101_24h           = Yd(end, 49 );
eSer_13C101D101_24h           = Yd(end, 50 );
eSer_13C110D101_24h           = Yd(end, 51 );
eSer_13C111D101_24h           = Yd(end, 52 );
eSer_13C000D110_24h           = Yd(end, 53 );
eSer_13C001D110_24h           = Yd(end, 54 );
eSer_13C010D110_24h           = Yd(end, 55 );
eSer_13C011D110_24h           = Yd(end, 56 );
eSer_13C100D110_24h           = Yd(end, 57 );
eSer_13C101D110_24h           = Yd(end, 58 );
eSer_13C110D110_24h           = Yd(end, 59 );
eSer_13C111D110_24h           = Yd(end, 60 );
eSer_13C000D111_24h           = Yd(end, 61 );
eSer_13C001D111_24h           = Yd(end, 62 );
eSer_13C010D111_24h           = Yd(end, 63 );
eSer_13C011D111_24h           = Yd(end, 64 );
eSer_13C100D111_24h           = Yd(end, 65 );
eSer_13C101D111_24h           = Yd(end, 66 );
eSer_13C110D111_24h           = Yd(end, 67 );
eSer_13C111D111_24h           = Yd(end, 68 );
eGly_13C00D00_24h             = Yd(end, 69 );
eGly_13C01D00_24h             = Yd(end, 70 );
eGly_13C10D00_24h             = Yd(end, 71 );
eGly_13C11D00_24h             = Yd(end, 72 );
eGly_13C00D01_24h             = Yd(end, 73 );
eGly_13C01D01_24h             = Yd(end, 74 );
eGly_13C10D01_24h             = Yd(end, 75 );
eGly_13C11D01_24h             = Yd(end, 76 );
eGly_13C00D10_24h             = Yd(end, 77 );
eGly_13C01D10_24h             = Yd(end, 78 );
eGly_13C10D10_24h             = Yd(end, 79 );
eGly_13C11D10_24h             = Yd(end, 80 );
eGly_13C00D11_24h             = Yd(end, 81 );
eGly_13C01D11_24h             = Yd(end, 82 );
eGly_13C10D11_24h             = Yd(end, 83 );
eGly_13C11D11_24h             = Yd(end, 84 );
cSerPool1_13C000D000_24h      = Yd(end, 85 );
cSerPool1_13C001D000_24h      = Yd(end, 86 );
cSerPool1_13C010D000_24h      = Yd(end, 87 );
cSerPool1_13C011D000_24h      = Yd(end, 88 );
cSerPool1_13C100D000_24h      = Yd(end, 89 );
cSerPool1_13C101D000_24h      = Yd(end, 90 );
cSerPool1_13C110D000_24h      = Yd(end, 91 );
cSerPool1_13C111D000_24h      = Yd(end, 92 );
cSerPool1_13C000D001_24h      = Yd(end, 93 );
cSerPool1_13C001D001_24h      = Yd(end, 94 );
cSerPool1_13C010D001_24h      = Yd(end, 95 );
cSerPool1_13C011D001_24h      = Yd(end, 96 );
cSerPool1_13C100D001_24h      = Yd(end, 97 );
cSerPool1_13C101D001_24h      = Yd(end, 98 );
cSerPool1_13C110D001_24h      = Yd(end, 99 );
cSerPool1_13C111D001_24h      = Yd(end, 100);
cSerPool1_13C000D010_24h      = Yd(end, 101);
cSerPool1_13C001D010_24h      = Yd(end, 102);
cSerPool1_13C010D010_24h      = Yd(end, 103);
cSerPool1_13C011D010_24h      = Yd(end, 104);
cSerPool1_13C100D010_24h      = Yd(end, 105);
cSerPool1_13C101D010_24h      = Yd(end, 106);
cSerPool1_13C110D010_24h      = Yd(end, 107);
cSerPool1_13C111D010_24h      = Yd(end, 108);
cSerPool1_13C000D011_24h      = Yd(end, 109);
cSerPool1_13C001D011_24h      = Yd(end, 110);
cSerPool1_13C010D011_24h      = Yd(end, 111);
cSerPool1_13C011D011_24h      = Yd(end, 112);
cSerPool1_13C100D011_24h      = Yd(end, 113);
cSerPool1_13C101D011_24h      = Yd(end, 114);
cSerPool1_13C110D011_24h      = Yd(end, 115);
cSerPool1_13C111D011_24h      = Yd(end, 116);
cSerPool1_13C000D100_24h      = Yd(end, 117);
cSerPool1_13C001D100_24h      = Yd(end, 118);
cSerPool1_13C010D100_24h      = Yd(end, 119);
cSerPool1_13C011D100_24h      = Yd(end, 120);
cSerPool1_13C100D100_24h      = Yd(end, 121);
cSerPool1_13C101D100_24h      = Yd(end, 122);
cSerPool1_13C110D100_24h      = Yd(end, 123);
cSerPool1_13C111D100_24h      = Yd(end, 124);
cSerPool1_13C000D101_24h      = Yd(end, 125);
cSerPool1_13C001D101_24h      = Yd(end, 126);
cSerPool1_13C010D101_24h      = Yd(end, 127);
cSerPool1_13C011D101_24h      = Yd(end, 128);
cSerPool1_13C100D101_24h      = Yd(end, 129);
cSerPool1_13C101D101_24h      = Yd(end, 130);
cSerPool1_13C110D101_24h      = Yd(end, 131);
cSerPool1_13C111D101_24h      = Yd(end, 132);
cSerPool1_13C000D110_24h      = Yd(end, 133);
cSerPool1_13C001D110_24h      = Yd(end, 134);
cSerPool1_13C010D110_24h      = Yd(end, 135);
cSerPool1_13C011D110_24h      = Yd(end, 136);
cSerPool1_13C100D110_24h      = Yd(end, 137);
cSerPool1_13C101D110_24h      = Yd(end, 138);
cSerPool1_13C110D110_24h      = Yd(end, 139);
cSerPool1_13C111D110_24h      = Yd(end, 140);
cSerPool1_13C000D111_24h      = Yd(end, 141);
cSerPool1_13C001D111_24h      = Yd(end, 142);
cSerPool1_13C010D111_24h      = Yd(end, 143);
cSerPool1_13C011D111_24h      = Yd(end, 144);
cSerPool1_13C100D111_24h      = Yd(end, 145);
cSerPool1_13C101D111_24h      = Yd(end, 146);
cSerPool1_13C110D111_24h      = Yd(end, 147);
cSerPool1_13C111D111_24h      = Yd(end, 148);
cGlyPool1_13C00D00_24h        = Yd(end, 149);
cGlyPool1_13C01D00_24h        = Yd(end, 150);
cGlyPool1_13C10D00_24h        = Yd(end, 151);
cGlyPool1_13C11D00_24h        = Yd(end, 152);
cGlyPool1_13C00D01_24h        = Yd(end, 153);
cGlyPool1_13C01D01_24h        = Yd(end, 154);
cGlyPool1_13C10D01_24h        = Yd(end, 155);
cGlyPool1_13C11D01_24h        = Yd(end, 156);
cGlyPool1_13C00D10_24h        = Yd(end, 157);
cGlyPool1_13C01D10_24h        = Yd(end, 158);
cGlyPool1_13C10D10_24h        = Yd(end, 159);
cGlyPool1_13C11D10_24h        = Yd(end, 160);
cGlyPool1_13C00D11_24h        = Yd(end, 161);
cGlyPool1_13C01D11_24h        = Yd(end, 162);
cGlyPool1_13C10D11_24h        = Yd(end, 163);
cGlyPool1_13C11D11_24h        = Yd(end, 164);
cSerPool2_13C000D000_24h      = Yd(end, 165);
cSerPool2_13C001D000_24h      = Yd(end, 166);
cSerPool2_13C010D000_24h      = Yd(end, 167);
cSerPool2_13C011D000_24h      = Yd(end, 168);
cSerPool2_13C100D000_24h      = Yd(end, 169);
cSerPool2_13C101D000_24h      = Yd(end, 170);
cSerPool2_13C110D000_24h      = Yd(end, 171);
cSerPool2_13C111D000_24h      = Yd(end, 172);
cSerPool2_13C000D001_24h      = Yd(end, 173);
cSerPool2_13C001D001_24h      = Yd(end, 174);
cSerPool2_13C010D001_24h      = Yd(end, 175);
cSerPool2_13C011D001_24h      = Yd(end, 176);
cSerPool2_13C100D001_24h      = Yd(end, 177);
cSerPool2_13C101D001_24h      = Yd(end, 178);
cSerPool2_13C110D001_24h      = Yd(end, 179);
cSerPool2_13C111D001_24h      = Yd(end, 180);
cSerPool2_13C000D010_24h      = Yd(end, 181);
cSerPool2_13C001D010_24h      = Yd(end, 182);
cSerPool2_13C010D010_24h      = Yd(end, 183);
cSerPool2_13C011D010_24h      = Yd(end, 184);
cSerPool2_13C100D010_24h      = Yd(end, 185);
cSerPool2_13C101D010_24h      = Yd(end, 186);
cSerPool2_13C110D010_24h      = Yd(end, 187);
cSerPool2_13C111D010_24h      = Yd(end, 188);
cSerPool2_13C000D011_24h      = Yd(end, 189);
cSerPool2_13C001D011_24h      = Yd(end, 190);
cSerPool2_13C010D011_24h      = Yd(end, 191);
cSerPool2_13C011D011_24h      = Yd(end, 192);
cSerPool2_13C100D011_24h      = Yd(end, 193);
cSerPool2_13C101D011_24h      = Yd(end, 194);
cSerPool2_13C110D011_24h      = Yd(end, 195);
cSerPool2_13C111D011_24h      = Yd(end, 196);
cSerPool2_13C000D100_24h      = Yd(end, 197);
cSerPool2_13C001D100_24h      = Yd(end, 198);
cSerPool2_13C010D100_24h      = Yd(end, 199);
cSerPool2_13C011D100_24h      = Yd(end, 200);
cSerPool2_13C100D100_24h      = Yd(end, 201);
cSerPool2_13C101D100_24h      = Yd(end, 202);
cSerPool2_13C110D100_24h      = Yd(end, 203);
cSerPool2_13C111D100_24h      = Yd(end, 204);
cSerPool2_13C000D101_24h      = Yd(end, 205);
cSerPool2_13C001D101_24h      = Yd(end, 206);
cSerPool2_13C010D101_24h      = Yd(end, 207);
cSerPool2_13C011D101_24h      = Yd(end, 208);
cSerPool2_13C100D101_24h      = Yd(end, 209);
cSerPool2_13C101D101_24h      = Yd(end, 210);
cSerPool2_13C110D101_24h      = Yd(end, 211);
cSerPool2_13C111D101_24h      = Yd(end, 212);
cSerPool2_13C000D110_24h      = Yd(end, 213);
cSerPool2_13C001D110_24h      = Yd(end, 214);
cSerPool2_13C010D110_24h      = Yd(end, 215);
cSerPool2_13C011D110_24h      = Yd(end, 216);
cSerPool2_13C100D110_24h      = Yd(end, 217);
cSerPool2_13C101D110_24h      = Yd(end, 218);
cSerPool2_13C110D110_24h      = Yd(end, 219);
cSerPool2_13C111D110_24h      = Yd(end, 220);
cSerPool2_13C000D111_24h      = Yd(end, 221);
cSerPool2_13C001D111_24h      = Yd(end, 222);
cSerPool2_13C010D111_24h      = Yd(end, 223);
cSerPool2_13C011D111_24h      = Yd(end, 224);
cSerPool2_13C100D111_24h      = Yd(end, 225);
cSerPool2_13C101D111_24h      = Yd(end, 226);
cSerPool2_13C110D111_24h      = Yd(end, 227);
cSerPool2_13C111D111_24h      = Yd(end, 228);
cGlyPool2_13C00D00_24h        = Yd(end, 229);
cGlyPool2_13C01D00_24h        = Yd(end, 230);
cGlyPool2_13C10D00_24h        = Yd(end, 231);
cGlyPool2_13C11D00_24h        = Yd(end, 232);
cGlyPool2_13C00D01_24h        = Yd(end, 233);
cGlyPool2_13C01D01_24h        = Yd(end, 234);
cGlyPool2_13C10D01_24h        = Yd(end, 235);
cGlyPool2_13C11D01_24h        = Yd(end, 236);
cGlyPool2_13C00D10_24h        = Yd(end, 237);
cGlyPool2_13C01D10_24h        = Yd(end, 238);
cGlyPool2_13C10D10_24h        = Yd(end, 239);
cGlyPool2_13C11D10_24h        = Yd(end, 240);
cGlyPool2_13C00D11_24h        = Yd(end, 241);
cGlyPool2_13C01D11_24h        = Yd(end, 242);
cGlyPool2_13C10D11_24h        = Yd(end, 243);
cGlyPool2_13C11D11_24h        = Yd(end, 244);
cSerPool3_13C000D000_24h      = Yd(end, 245);
cSerPool3_13C001D000_24h      = Yd(end, 246);
cSerPool3_13C010D000_24h      = Yd(end, 247);
cSerPool3_13C011D000_24h      = Yd(end, 248);
cSerPool3_13C100D000_24h      = Yd(end, 249);
cSerPool3_13C101D000_24h      = Yd(end, 250);
cSerPool3_13C110D000_24h      = Yd(end, 251);
cSerPool3_13C111D000_24h      = Yd(end, 252);
cSerPool3_13C000D001_24h      = Yd(end, 253);
cSerPool3_13C001D001_24h      = Yd(end, 254);
cSerPool3_13C010D001_24h      = Yd(end, 255);
cSerPool3_13C011D001_24h      = Yd(end, 256);
cSerPool3_13C100D001_24h      = Yd(end, 257);
cSerPool3_13C101D001_24h      = Yd(end, 258);
cSerPool3_13C110D001_24h      = Yd(end, 259);
cSerPool3_13C111D001_24h      = Yd(end, 260);
cSerPool3_13C000D010_24h      = Yd(end, 261);
cSerPool3_13C001D010_24h      = Yd(end, 262);
cSerPool3_13C010D010_24h      = Yd(end, 263);
cSerPool3_13C011D010_24h      = Yd(end, 264);
cSerPool3_13C100D010_24h      = Yd(end, 265);
cSerPool3_13C101D010_24h      = Yd(end, 266);
cSerPool3_13C110D010_24h      = Yd(end, 267);
cSerPool3_13C111D010_24h      = Yd(end, 268);
cSerPool3_13C000D011_24h      = Yd(end, 269);
cSerPool3_13C001D011_24h      = Yd(end, 270);
cSerPool3_13C010D011_24h      = Yd(end, 271);
cSerPool3_13C011D011_24h      = Yd(end, 272);
cSerPool3_13C100D011_24h      = Yd(end, 273);
cSerPool3_13C101D011_24h      = Yd(end, 274);
cSerPool3_13C110D011_24h      = Yd(end, 275);
cSerPool3_13C111D011_24h      = Yd(end, 276);
cSerPool3_13C000D100_24h      = Yd(end, 277);
cSerPool3_13C001D100_24h      = Yd(end, 278);
cSerPool3_13C010D100_24h      = Yd(end, 279);
cSerPool3_13C011D100_24h      = Yd(end, 280);
cSerPool3_13C100D100_24h      = Yd(end, 281);
cSerPool3_13C101D100_24h      = Yd(end, 282);
cSerPool3_13C110D100_24h      = Yd(end, 283);
cSerPool3_13C111D100_24h      = Yd(end, 284);
cSerPool3_13C000D101_24h      = Yd(end, 285);
cSerPool3_13C001D101_24h      = Yd(end, 286);
cSerPool3_13C010D101_24h      = Yd(end, 287);
cSerPool3_13C011D101_24h      = Yd(end, 288);
cSerPool3_13C100D101_24h      = Yd(end, 289);
cSerPool3_13C101D101_24h      = Yd(end, 290);
cSerPool3_13C110D101_24h      = Yd(end, 291);
cSerPool3_13C111D101_24h      = Yd(end, 292);
cSerPool3_13C000D110_24h      = Yd(end, 293);
cSerPool3_13C001D110_24h      = Yd(end, 294);
cSerPool3_13C010D110_24h      = Yd(end, 295);
cSerPool3_13C011D110_24h      = Yd(end, 296);
cSerPool3_13C100D110_24h      = Yd(end, 297);
cSerPool3_13C101D110_24h      = Yd(end, 298);
cSerPool3_13C110D110_24h      = Yd(end, 299);
cSerPool3_13C111D110_24h      = Yd(end, 300);
cSerPool3_13C000D111_24h      = Yd(end, 301);
cSerPool3_13C001D111_24h      = Yd(end, 302);
cSerPool3_13C010D111_24h      = Yd(end, 303);
cSerPool3_13C011D111_24h      = Yd(end, 304);
cSerPool3_13C100D111_24h      = Yd(end, 305);
cSerPool3_13C101D111_24h      = Yd(end, 306);
cSerPool3_13C110D111_24h      = Yd(end, 307);
cSerPool3_13C111D111_24h      = Yd(end, 308);
cGlyPool3_13C00D00_24h        = Yd(end, 309);
cGlyPool3_13C01D00_24h        = Yd(end, 310);
cGlyPool3_13C10D00_24h        = Yd(end, 311);
cGlyPool3_13C11D00_24h        = Yd(end, 312);
cGlyPool3_13C00D01_24h        = Yd(end, 313);
cGlyPool3_13C01D01_24h        = Yd(end, 314);
cGlyPool3_13C10D01_24h        = Yd(end, 315);
cGlyPool3_13C11D01_24h        = Yd(end, 316);
cGlyPool3_13C00D10_24h        = Yd(end, 317);
cGlyPool3_13C01D10_24h        = Yd(end, 318);
cGlyPool3_13C10D10_24h        = Yd(end, 319);
cGlyPool3_13C11D10_24h        = Yd(end, 320);
cGlyPool3_13C00D11_24h        = Yd(end, 321);
cGlyPool3_13C01D11_24h        = Yd(end, 322);
cGlyPool3_13C10D11_24h        = Yd(end, 323);
cGlyPool3_13C11D11_24h        = Yd(end, 324);
mSerPool_13C000D000_24h       = Yd(end, 325);
mSerPool_13C001D000_24h       = Yd(end, 326);
mSerPool_13C010D000_24h       = Yd(end, 327);
mSerPool_13C011D000_24h       = Yd(end, 328);
mSerPool_13C100D000_24h       = Yd(end, 329);
mSerPool_13C101D000_24h       = Yd(end, 330);
mSerPool_13C110D000_24h       = Yd(end, 331);
mSerPool_13C111D000_24h       = Yd(end, 332);
mSerPool_13C000D001_24h       = Yd(end, 333);
mSerPool_13C001D001_24h       = Yd(end, 334);
mSerPool_13C010D001_24h       = Yd(end, 335);
mSerPool_13C011D001_24h       = Yd(end, 336);
mSerPool_13C100D001_24h       = Yd(end, 337);
mSerPool_13C101D001_24h       = Yd(end, 338);
mSerPool_13C110D001_24h       = Yd(end, 339);
mSerPool_13C111D001_24h       = Yd(end, 340);
mSerPool_13C000D010_24h       = Yd(end, 341);
mSerPool_13C001D010_24h       = Yd(end, 342);
mSerPool_13C010D010_24h       = Yd(end, 343);
mSerPool_13C011D010_24h       = Yd(end, 344);
mSerPool_13C100D010_24h       = Yd(end, 345);
mSerPool_13C101D010_24h       = Yd(end, 346);
mSerPool_13C110D010_24h       = Yd(end, 347);
mSerPool_13C111D010_24h       = Yd(end, 348);
mSerPool_13C000D011_24h       = Yd(end, 349);
mSerPool_13C001D011_24h       = Yd(end, 350);
mSerPool_13C010D011_24h       = Yd(end, 351);
mSerPool_13C011D011_24h       = Yd(end, 352);
mSerPool_13C100D011_24h       = Yd(end, 353);
mSerPool_13C101D011_24h       = Yd(end, 354);
mSerPool_13C110D011_24h       = Yd(end, 355);
mSerPool_13C111D011_24h       = Yd(end, 356);
mSerPool_13C000D100_24h       = Yd(end, 357);
mSerPool_13C001D100_24h       = Yd(end, 358);
mSerPool_13C010D100_24h       = Yd(end, 359);
mSerPool_13C011D100_24h       = Yd(end, 360);
mSerPool_13C100D100_24h       = Yd(end, 361);
mSerPool_13C101D100_24h       = Yd(end, 362);
mSerPool_13C110D100_24h       = Yd(end, 363);
mSerPool_13C111D100_24h       = Yd(end, 364);
mSerPool_13C000D101_24h       = Yd(end, 365);
mSerPool_13C001D101_24h       = Yd(end, 366);
mSerPool_13C010D101_24h       = Yd(end, 367);
mSerPool_13C011D101_24h       = Yd(end, 368);
mSerPool_13C100D101_24h       = Yd(end, 369);
mSerPool_13C101D101_24h       = Yd(end, 370);
mSerPool_13C110D101_24h       = Yd(end, 371);
mSerPool_13C111D101_24h       = Yd(end, 372);
mSerPool_13C000D110_24h       = Yd(end, 373);
mSerPool_13C001D110_24h       = Yd(end, 374);
mSerPool_13C010D110_24h       = Yd(end, 375);
mSerPool_13C011D110_24h       = Yd(end, 376);
mSerPool_13C100D110_24h       = Yd(end, 377);
mSerPool_13C101D110_24h       = Yd(end, 378);
mSerPool_13C110D110_24h       = Yd(end, 379);
mSerPool_13C111D110_24h       = Yd(end, 380);
mSerPool_13C000D111_24h       = Yd(end, 381);
mSerPool_13C001D111_24h       = Yd(end, 382);
mSerPool_13C010D111_24h       = Yd(end, 383);
mSerPool_13C011D111_24h       = Yd(end, 384);
mSerPool_13C100D111_24h       = Yd(end, 385);
mSerPool_13C101D111_24h       = Yd(end, 386);
mSerPool_13C110D111_24h       = Yd(end, 387);
mSerPool_13C111D111_24h       = Yd(end, 388);
mGlyPool_13C00D00_24h         = Yd(end, 389);
mGlyPool_13C01D00_24h         = Yd(end, 390);
mGlyPool_13C10D00_24h         = Yd(end, 391);
mGlyPool_13C11D00_24h         = Yd(end, 392);
mGlyPool_13C00D01_24h         = Yd(end, 393);
mGlyPool_13C01D01_24h         = Yd(end, 394);
mGlyPool_13C10D01_24h         = Yd(end, 395);
mGlyPool_13C11D01_24h         = Yd(end, 396);
mGlyPool_13C00D10_24h         = Yd(end, 397);
mGlyPool_13C01D10_24h         = Yd(end, 398);
mGlyPool_13C10D10_24h         = Yd(end, 399);
mGlyPool_13C11D10_24h         = Yd(end, 400);
mGlyPool_13C00D11_24h         = Yd(end, 401);
mGlyPool_13C01D11_24h         = Yd(end, 402);
mGlyPool_13C10D11_24h         = Yd(end, 403);
mGlyPool_13C11D11_24h         = Yd(end, 404);
cMethyleneTHF_13C0D00_24h     = Yd(end, 405);
cMethyleneTHF_13C1D00_24h     = Yd(end, 406);
cMethyleneTHF_13C0D01_24h     = Yd(end, 407);
cMethyleneTHF_13C1D01_24h     = Yd(end, 408);
cMethyleneTHF_13C0D10_24h     = Yd(end, 409);
cMethyleneTHF_13C1D10_24h     = Yd(end, 410);
cMethyleneTHF_13C0D11_24h     = Yd(end, 411);
cMethyleneTHF_13C1D11_24h     = Yd(end, 412);
cFormylTHF_13C0D0_24h         = Yd(end, 413);
cFormylTHF_13C1D0_24h         = Yd(end, 414);
cFormylTHF_13C0D1_24h         = Yd(end, 415);
cFormylTHF_13C1D1_24h         = Yd(end, 416);
cFormate_13C0D0_24h           = Yd(end, 417);
cFormate_13C1D0_24h           = Yd(end, 418);
cFormate_13C0D1_24h           = Yd(end, 419);
cFormate_13C1D1_24h           = Yd(end, 420);
cTHF_24h                      = Yd(end, 421);
mMethyleneTHF_13C0D00_24h     = Yd(end, 422);
mMethyleneTHF_13C1D00_24h     = Yd(end, 423);
mMethyleneTHF_13C0D01_24h     = Yd(end, 424);
mMethyleneTHF_13C1D01_24h     = Yd(end, 425);
mMethyleneTHF_13C0D10_24h     = Yd(end, 426);
mMethyleneTHF_13C1D10_24h     = Yd(end, 427);
mMethyleneTHF_13C0D11_24h     = Yd(end, 428);
mMethyleneTHF_13C1D11_24h     = Yd(end, 429);
mFormylTHF_13C0D0_24h         = Yd(end, 430);
mFormylTHF_13C1D0_24h         = Yd(end, 431);
mFormylTHF_13C0D1_24h         = Yd(end, 432);
mFormylTHF_13C1D1_24h         = Yd(end, 433);
mFormate_13C0D0_24h           = Yd(end, 434);
mFormate_13C1D0_24h           = Yd(end, 435);
mFormate_13C0D1_24h           = Yd(end, 436);
mFormate_13C1D1_24h           = Yd(end, 437);
mTHF_24h                      = Yd(end, 438);
cPRPP_13C0_24h                = Yd(end, 439);
cPRPP_13C1_24h                = Yd(end, 440);
cGAR_13C000_24h               = Yd(end, 441);
cGAR_13C001_24h               = Yd(end, 442);
cGAR_13C010_24h               = Yd(end, 443);
cGAR_13C011_24h               = Yd(end, 444);
cGAR_13C100_24h               = Yd(end, 445);
cGAR_13C101_24h               = Yd(end, 446);
cGAR_13C110_24h               = Yd(end, 447);
cGAR_13C111_24h               = Yd(end, 448);
cFGAR_13C0000D0_24h           = Yd(end, 449);
cFGAR_13C0001D0_24h           = Yd(end, 450);
cFGAR_13C0010D0_24h           = Yd(end, 451);
cFGAR_13C0011D0_24h           = Yd(end, 452);
cFGAR_13C0100D0_24h           = Yd(end, 453);
cFGAR_13C0101D0_24h           = Yd(end, 454);
cFGAR_13C0110D0_24h           = Yd(end, 455);
cFGAR_13C0111D0_24h           = Yd(end, 456);
cFGAR_13C1000D0_24h           = Yd(end, 457);
cFGAR_13C1001D0_24h           = Yd(end, 458);
cFGAR_13C1010D0_24h           = Yd(end, 459);
cFGAR_13C1011D0_24h           = Yd(end, 460);
cFGAR_13C1100D0_24h           = Yd(end, 461);
cFGAR_13C1101D0_24h           = Yd(end, 462);
cFGAR_13C1110D0_24h           = Yd(end, 463);
cFGAR_13C1111D0_24h           = Yd(end, 464);
cFGAR_13C0000D1_24h           = Yd(end, 465);
cFGAR_13C0001D1_24h           = Yd(end, 466);
cFGAR_13C0010D1_24h           = Yd(end, 467);
cFGAR_13C0011D1_24h           = Yd(end, 468);
cFGAR_13C0100D1_24h           = Yd(end, 469);
cFGAR_13C0101D1_24h           = Yd(end, 470);
cFGAR_13C0110D1_24h           = Yd(end, 471);
cFGAR_13C0111D1_24h           = Yd(end, 472);
cFGAR_13C1000D1_24h           = Yd(end, 473);
cFGAR_13C1001D1_24h           = Yd(end, 474);
cFGAR_13C1010D1_24h           = Yd(end, 475);
cFGAR_13C1011D1_24h           = Yd(end, 476);
cFGAR_13C1100D1_24h           = Yd(end, 477);
cFGAR_13C1101D1_24h           = Yd(end, 478);
cFGAR_13C1110D1_24h           = Yd(end, 479);
cFGAR_13C1111D1_24h           = Yd(end, 480);
cAMP_13C00000D00_24h          = Yd(end, 481);
cAMP_13C00001D00_24h          = Yd(end, 482);
cAMP_13C00010D00_24h          = Yd(end, 483);
cAMP_13C00011D00_24h          = Yd(end, 484);
cAMP_13C00100D00_24h          = Yd(end, 485);
cAMP_13C00101D00_24h          = Yd(end, 486);
cAMP_13C00110D00_24h          = Yd(end, 487);
cAMP_13C00111D00_24h          = Yd(end, 488);
cAMP_13C01000D00_24h          = Yd(end, 489);
cAMP_13C01001D00_24h          = Yd(end, 490);
cAMP_13C01010D00_24h          = Yd(end, 491);
cAMP_13C01011D00_24h          = Yd(end, 492);
cAMP_13C01100D00_24h          = Yd(end, 493);
cAMP_13C01101D00_24h          = Yd(end, 494);
cAMP_13C01110D00_24h          = Yd(end, 495);
cAMP_13C01111D00_24h          = Yd(end, 496);
cAMP_13C10000D00_24h          = Yd(end, 497);
cAMP_13C10001D00_24h          = Yd(end, 498);
cAMP_13C10010D00_24h          = Yd(end, 499);
cAMP_13C10011D00_24h          = Yd(end, 500);
cAMP_13C10100D00_24h          = Yd(end, 501);
cAMP_13C10101D00_24h          = Yd(end, 502);
cAMP_13C10110D00_24h          = Yd(end, 503);
cAMP_13C10111D00_24h          = Yd(end, 504);
cAMP_13C11000D00_24h          = Yd(end, 505);
cAMP_13C11001D00_24h          = Yd(end, 506);
cAMP_13C11010D00_24h          = Yd(end, 507);
cAMP_13C11011D00_24h          = Yd(end, 508);
cAMP_13C11100D00_24h          = Yd(end, 509);
cAMP_13C11101D00_24h          = Yd(end, 510);
cAMP_13C11110D00_24h          = Yd(end, 511);
cAMP_13C11111D00_24h          = Yd(end, 512);
cAMP_13C00000D01_24h          = Yd(end, 513);
cAMP_13C00001D01_24h          = Yd(end, 514);
cAMP_13C00010D01_24h          = Yd(end, 515);
cAMP_13C00011D01_24h          = Yd(end, 516);
cAMP_13C00100D01_24h          = Yd(end, 517);
cAMP_13C00101D01_24h          = Yd(end, 518);
cAMP_13C00110D01_24h          = Yd(end, 519);
cAMP_13C00111D01_24h          = Yd(end, 520);
cAMP_13C01000D01_24h          = Yd(end, 521);
cAMP_13C01001D01_24h          = Yd(end, 522);
cAMP_13C01010D01_24h          = Yd(end, 523);
cAMP_13C01011D01_24h          = Yd(end, 524);
cAMP_13C01100D01_24h          = Yd(end, 525);
cAMP_13C01101D01_24h          = Yd(end, 526);
cAMP_13C01110D01_24h          = Yd(end, 527);
cAMP_13C01111D01_24h          = Yd(end, 528);
cAMP_13C10000D01_24h          = Yd(end, 529);
cAMP_13C10001D01_24h          = Yd(end, 530);
cAMP_13C10010D01_24h          = Yd(end, 531);
cAMP_13C10011D01_24h          = Yd(end, 532);
cAMP_13C10100D01_24h          = Yd(end, 533);
cAMP_13C10101D01_24h          = Yd(end, 534);
cAMP_13C10110D01_24h          = Yd(end, 535);
cAMP_13C10111D01_24h          = Yd(end, 536);
cAMP_13C11000D01_24h          = Yd(end, 537);
cAMP_13C11001D01_24h          = Yd(end, 538);
cAMP_13C11010D01_24h          = Yd(end, 539);
cAMP_13C11011D01_24h          = Yd(end, 540);
cAMP_13C11100D01_24h          = Yd(end, 541);
cAMP_13C11101D01_24h          = Yd(end, 542);
cAMP_13C11110D01_24h          = Yd(end, 543);
cAMP_13C11111D01_24h          = Yd(end, 544);
cAMP_13C00000D10_24h          = Yd(end, 545);
cAMP_13C00001D10_24h          = Yd(end, 546);
cAMP_13C00010D10_24h          = Yd(end, 547);
cAMP_13C00011D10_24h          = Yd(end, 548);
cAMP_13C00100D10_24h          = Yd(end, 549);
cAMP_13C00101D10_24h          = Yd(end, 550);
cAMP_13C00110D10_24h          = Yd(end, 551);
cAMP_13C00111D10_24h          = Yd(end, 552);
cAMP_13C01000D10_24h          = Yd(end, 553);
cAMP_13C01001D10_24h          = Yd(end, 554);
cAMP_13C01010D10_24h          = Yd(end, 555);
cAMP_13C01011D10_24h          = Yd(end, 556);
cAMP_13C01100D10_24h          = Yd(end, 557);
cAMP_13C01101D10_24h          = Yd(end, 558);
cAMP_13C01110D10_24h          = Yd(end, 559);
cAMP_13C01111D10_24h          = Yd(end, 560);
cAMP_13C10000D10_24h          = Yd(end, 561);
cAMP_13C10001D10_24h          = Yd(end, 562);
cAMP_13C10010D10_24h          = Yd(end, 563);
cAMP_13C10011D10_24h          = Yd(end, 564);
cAMP_13C10100D10_24h          = Yd(end, 565);
cAMP_13C10101D10_24h          = Yd(end, 566);
cAMP_13C10110D10_24h          = Yd(end, 567);
cAMP_13C10111D10_24h          = Yd(end, 568);
cAMP_13C11000D10_24h          = Yd(end, 569);
cAMP_13C11001D10_24h          = Yd(end, 570);
cAMP_13C11010D10_24h          = Yd(end, 571);
cAMP_13C11011D10_24h          = Yd(end, 572);
cAMP_13C11100D10_24h          = Yd(end, 573);
cAMP_13C11101D10_24h          = Yd(end, 574);
cAMP_13C11110D10_24h          = Yd(end, 575);
cAMP_13C11111D10_24h          = Yd(end, 576);
cAMP_13C00000D11_24h          = Yd(end, 577);
cAMP_13C00001D11_24h          = Yd(end, 578);
cAMP_13C00010D11_24h          = Yd(end, 579);
cAMP_13C00011D11_24h          = Yd(end, 580);
cAMP_13C00100D11_24h          = Yd(end, 581);
cAMP_13C00101D11_24h          = Yd(end, 582);
cAMP_13C00110D11_24h          = Yd(end, 583);
cAMP_13C00111D11_24h          = Yd(end, 584);
cAMP_13C01000D11_24h          = Yd(end, 585);
cAMP_13C01001D11_24h          = Yd(end, 586);
cAMP_13C01010D11_24h          = Yd(end, 587);
cAMP_13C01011D11_24h          = Yd(end, 588);
cAMP_13C01100D11_24h          = Yd(end, 589);
cAMP_13C01101D11_24h          = Yd(end, 590);
cAMP_13C01110D11_24h          = Yd(end, 591);
cAMP_13C01111D11_24h          = Yd(end, 592);
cAMP_13C10000D11_24h          = Yd(end, 593);
cAMP_13C10001D11_24h          = Yd(end, 594);
cAMP_13C10010D11_24h          = Yd(end, 595);
cAMP_13C10011D11_24h          = Yd(end, 596);
cAMP_13C10100D11_24h          = Yd(end, 597);
cAMP_13C10101D11_24h          = Yd(end, 598);
cAMP_13C10110D11_24h          = Yd(end, 599);
cAMP_13C10111D11_24h          = Yd(end, 600);
cAMP_13C11000D11_24h          = Yd(end, 601);
cAMP_13C11001D11_24h          = Yd(end, 602);
cAMP_13C11010D11_24h          = Yd(end, 603);
cAMP_13C11011D11_24h          = Yd(end, 604);
cAMP_13C11100D11_24h          = Yd(end, 605);
cAMP_13C11101D11_24h          = Yd(end, 606);
cAMP_13C11110D11_24h          = Yd(end, 607);
cAMP_13C11111D11_24h          = Yd(end, 608);










cSer13C0D0 = cSerPool1_13C000D000_24h + cSerPool2_13C000D000_24h + cSerPool2_13C000D000_24h + mSerPool_13C000D000_24h
cSer13C1D0 = cSerPool1_13C001D000_24h + cSerPool1_13C010D000_24h + cSerPool1_13C100D000_24h + cSerPool2_13C001D000_24h + cSerPool2_13C010D000_24h + cSerPool2_13C100D000_24h + cSerPool3_13C001D000_24h + cSerPool3_13C010D000_24h + cSerPool3_13C100D000_24h + mSerPool_13C001D000_24h + mSerPool_13C010D000_24h + mSerPool_13C100D000_24h

cSer13C2D0 = cSerPool1_13C011D000_24h + cSerPool1_13C101D000_24h + cSerPool1_13C110D000_24h + cSerPool2_13C011D000_24h + cSerPool2_13C101D000_24h + cSerPool2_13C110D000_24h + cSerPool3_13C011D000_24h + cSerPool3_13C101D000_24h + cSerPool3_13C110D000_24h + mSerPool_13C011D000_24h + mSerPool_13C101D000_24h + mSerPool_13C110D000_24h

cSer13C3D0 = cSerPool1_13C111D000_24h + cSerPool2_13C111D000_24h + cSerPool3_13C111D000_24h + mSerPool_13C111D000_24h


cSer13C0D1 = cSerPool1_13C000D010_24h + cSerPool2_13C000D010_24h + cSerPool3_13C000D010_24h + mSerPool_13C000D010_24h

cSer13C0D1 = cSerPool1_13C000D001_24h + cSerPool1_13C000D100_24h + cSerPool2_13C000D001_24h + cSerPool2_13C000D100_24h + cSerPool3_13C000D001_24h + cSerPool3_13C000D100_24h + mSerPool_13C000D001_24h + mSerPool_13C000D100_24h

cSer13C1D1 = cSerPool1_13C001D001_24h + cSerPool1_13C010D001_24h + cSerPool1_13C100D001_24h + cSerPool1_13C001D010_24h + cSerPool1_13C010D010_24h + cSerPool1_13C100D010_24h + cSerPool1_13C001D100_24h + cSerPool1_13C010D100_24h + cSerPool1_13C100D100_24h + ...
cSerPool2_13C001D001_24h + cSerPool2_13C010D001_24h + cSerPool2_13C100D001_24h + cSerPool2_13C001D010_24h + cSerPool2_13C010D010_24h + cSerPool2_13C100D010_24h + cSerPool2_13C001D100_24h + cSerPool2_13C010D100_24h + cSerPool2_13C100D100_24h + ...
cSerPool3_13C001D001_24h + cSerPool3_13C010D001_24h + cSerPool3_13C100D001_24h + cSerPool3_13C001D010_24h + cSerPool3_13C010D010_24h + cSerPool3_13C100D010_24h + cSerPool3_13C001D100_24h + cSerPool3_13C010D100_24h + cSerPool3_13C100D100_24h + ...
mSerPool_13C001D001_24h + mSerPool_13C010D001_24h + mSerPool_13C100D001_24h + mSerPool_13C001D010_24h + mSerPool_13C010D010_24h + mSerPool_13C100D010_24h + mSerPool_13C001D100_24h + mSerPool_13C010D100_24h + mSerPool_13C100D100_24h

cSer13C2D1 = cSerPool1_13C011D001_24h + cSerPool1_13C101D001_24h + cSerPool1_13C110D001_24h + cSerPool1_13C011D010_24h + cSerPool1_13C101D010_24h + cSerPool1_13C110D010_24h + cSerPool1_13C011D100_24h + cSerPool1_13C101D100_24h + cSerPool1_13C110D100_24h + ...
cSerPool2_13C011D001_24h + cSerPool2_13C101D001_24h + cSerPool2_13C110D001_24h + cSerPool2_13C011D010_24h + cSerPool2_13C101D010_24h + cSerPool2_13C110D010_24h + cSerPool2_13C011D100_24h + cSerPool2_13C101D100_24h + cSerPool2_13C110D100_24h + ...
cSerPool3_13C011D001_24h + cSerPool3_13C101D001_24h + cSerPool3_13C110D001_24h + cSerPool3_13C011D010_24h + cSerPool3_13C101D010_24h + cSerPool3_13C110D010_24h + cSerPool3_13C011D100_24h + cSerPool3_13C101D100_24h + cSerPool3_13C110D100_24h + ...
mSerPool_13C011D001_24h + mSerPool_13C101D001_24h + mSerPool_13C110D001_24h + mSerPool_13C011D010_24h + mSerPool_13C101D010_24h + mSerPool_13C110D010_24h + mSerPool_13C011D100_24h + mSerPool_13C101D100_24h + mSerPool_13C110D100_24h

cSer13C3D1 = cSerPool1_13C111D001_24h + cSerPool1_13C111D010_24h + cSerPool1_13C111D100_24h + cSerPool2_13C111D001_24h + cSerPool2_13C111D010_24h + cSerPool2_13C111D100_24h + cSerPool3_13C111D001_24h + cSerPool3_13C111D010_24h + cSerPool3_13C111D100_24h + mSerPool_13C111D001_24h + mSerPool_13C111D010_24h + mSerPool_13C111D100_24h


cSer13C0D2 = cSerPool1_13C000D011_24h + cSerPool1_13C000D101_24h + cSerPool1_13C000D110_24h + cSerPool2_13C000D011_24h + cSerPool2_13C000D101_24h + cSerPool2_13C000D110_24h + cSerPool3_13C000D011_24h + cSerPool3_13C000D101_24h + cSerPool3_13C000D110_24h + mSerPool_13C000D011_24h + mSerPool_13C000D101_24h + mSerPool_13C000D110_24h
 
cSer13C1D2 = cSerPool1_13C001D011_24h + cSerPool1_13C010D011_24h + cSerPool1_13C100D011_24h + cSerPool1_13C001D101_24h + cSerPool1_13C010D101_24h + cSerPool1_13C100D101_24h + cSerPool1_13C001D110_24h + cSerPool1_13C010D110_24h + cSerPool1_13C100D110_24h + ...
cSerPool2_13C001D011_24h + cSerPool2_13C010D011_24h + cSerPool2_13C100D011_24h + cSerPool2_13C001D101_24h + cSerPool2_13C010D101_24h + cSerPool2_13C100D101_24h + cSerPool2_13C001D110_24h + cSerPool2_13C010D110_24h + cSerPool2_13C100D110_24h + ...
cSerPool3_13C001D011_24h + cSerPool3_13C010D011_24h + cSerPool3_13C100D011_24h + cSerPool3_13C001D101_24h + cSerPool3_13C010D101_24h + cSerPool3_13C100D101_24h + cSerPool3_13C001D110_24h + cSerPool3_13C010D110_24h + cSerPool3_13C100D110_24h + ...
mSerPool_13C001D011_24h + mSerPool_13C010D011_24h + mSerPool_13C100D011_24h + mSerPool_13C001D101_24h + mSerPool_13C010D101_24h + mSerPool_13C100D101_24h + mSerPool_13C001D110_24h + mSerPool_13C010D110_24h + mSerPool_13C100D110_24h

cSer13C2D2 = cSerPool1_13C011D011_24h + cSerPool1_13C101D011_24h + cSerPool1_13C110D011_24h + cSerPool1_13C011D101_24h + cSerPool1_13C101D101_24h + cSerPool1_13C110D101_24h + cSerPool1_13C011D110_24h + cSerPool1_13C101D110_24h + cSerPool1_13C110D110_24h + ...
cSerPool2_13C011D011_24h + cSerPool2_13C101D011_24h + cSerPool2_13C110D011_24h + cSerPool2_13C011D101_24h + cSerPool2_13C101D101_24h + cSerPool2_13C110D101_24h + cSerPool2_13C011D110_24h + cSerPool2_13C101D110_24h + cSerPool2_13C110D110_24h + ...
cSerPool3_13C011D011_24h + cSerPool3_13C101D011_24h + cSerPool3_13C110D011_24h + cSerPool3_13C011D101_24h + cSerPool3_13C101D101_24h + cSerPool3_13C110D101_24h + cSerPool3_13C011D110_24h + cSerPool3_13C101D110_24h + cSerPool3_13C110D110_24h + ...
mSerPool_13C011D011_24h + mSerPool_13C101D011_24h + mSerPool_13C110D011_24h + mSerPool_13C011D101_24h + mSerPool_13C101D101_24h + mSerPool_13C110D101_24h + mSerPool_13C011D110_24h + mSerPool_13C101D110_24h + mSerPool_13C110D110_24h

cSer13C3D2 = cSerPool1_13C111D011_24h + cSerPool1_13C111D101_24h + cSerPool1_13C111D110_24h + cSerPool2_13C111D011_24h + cSerPool2_13C111D101_24h + cSerPool2_13C111D110_24h + cSerPool3_13C111D011_24h + cSerPool3_13C111D101_24h + cSerPool3_13C111D110_24h + mSerPool_13C111D011_24h + mSerPool_13C111D101_24h + mSerPool_13C111D110_24h


cSer13C0D3 = cSerPool1_13C000D111_24h + cSerPool2_13C000D111_24h + cSerPool3_13C000D111_24h + mSerPool_13C000D111_24h

cSer13C1D3 = cSerPool1_13C001D111_24h + cSerPool1_13C010D111_24h + cSerPool1_13C100D111_24h + cSerPool2_13C001D111_24h + cSerPool2_13C010D111_24h + cSerPool2_13C100D111_24h + cSerPool3_13C001D111_24h + cSerPool3_13C010D111_24h + cSerPool3_13C100D111_24h + mSerPool_13C001D111_24h + mSerPool_13C010D111_24h + mSerPool_13C100D111_24h

cSer13C2D3 = cSerPool1_13C011D111_24h + cSerPool1_13C101D111_24h + cSerPool1_13C110D111_24h + cSerPool2_13C011D111_24h + cSerPool2_13C101D111_24h + cSerPool2_13C110D111_24h + cSerPool3_13C011D111_24h + cSerPool3_13C101D111_24h + cSerPool3_13C110D111_24h + mSerPool_13C011D111_24h + mSerPool_13C101D111_24h + mSerPool_13C110D111_24h

cSer13C3D3 = cSerPool1_13C111D111_24h + cSerPool2_13C111D111_24h + cSerPool3_13C111D111_24h + mSerPool_13C111D111_24h





cGly13C0D0 = cGlyPool1_13C00D00_24h + cGlyPool2_13C00D00_24h + cGlyPool3_13C00D00_24h + mGlyPool_13C00D00_24h

cGly13C1D0 = cGlyPool1_13C01D00_24h + cGlyPool1_13C10D00_24h + cGlyPool2_13C01D00_24h + cGlyPool2_13C10D00_24h + cGlyPool3_13C01D00_24h + cGlyPool3_13C10D00_24h + mGlyPool_13C01D00_24h + mGlyPool_13C10D00_24h

cGly13C2D0 = cGlyPool1_13C11D00_24h + cGlyPool2_13C11D00_24h + cGlyPool3_13C11D00_24h + mGlyPool_13C11D00_24h


cGly13C0D1 = cGlyPool1_13C00D01_24h + cGlyPool1_13C00D10_24h + cGlyPool2_13C00D01_24h + cGlyPool2_13C00D10_24h + cGlyPool3_13C00D01_24h + cGlyPool3_13C00D10_24h + mGlyPool_13C00D01_24h + mGlyPool_13C00D10_24h

cGly13C1D1 = cGlyPool1_13C01D01_24h + cGlyPool1_13C10D01_24h + cGlyPool1_13C01D10_24h + cGlyPool1_13C10D10_24h + cGlyPool2_13C01D01_24h + cGlyPool2_13C10D01_24h + cGlyPool2_13C01D10_24h + cGlyPool2_13C10D10_24h + cGlyPool3_13C01D01_24h + cGlyPool3_13C10D01_24h + cGlyPool3_13C01D10_24h + cGlyPool3_13C10D10_24h + mGlyPool_13C01D01_24h + mGlyPool_13C10D01_24h + mGlyPool_13C01D10_24h + mGlyPool_13C10D10_24h

cGly13C2D1 = cGlyPool1_13C11D01_24h + cGlyPool1_13C11D10_24h + cGlyPool2_13C11D01_24h + cGlyPool2_13C11D10_24h + cGlyPool3_13C11D01_24h + cGlyPool3_13C11D10_24h + mGlyPool_13C11D01_24h + mGlyPool_13C11D10_24h

cGly13C0D2 = cGlyPool1_13C00D11_24h + cGlyPool2_13C00D11_24h + cGlyPool3_13C00D11_24h + mGlyPool_13C00D11_24h

cGly13C1D2 = cGlyPool1_13C01D11_24h + cGlyPool1_13C10D11_24h + cGlyPool2_13C01D11_24h + cGlyPool2_13C10D11_24h + cGlyPool3_13C01D11_24h + cGlyPool3_13C10D11_24h + mGlyPool_13C01D11_24h + mGlyPool_13C10D11_24h

cGly13C2D2 = cGlyPool1_13C11D11_24h + cGlyPool2_13C11D11_24h + cGlyPool3_13C11D11_24h + mGlyPool_13C11D11_24h






eSer13C0D0 = eSer_13C000D000_24h

eSer13C1D0 = eSer_13C001D000_24h + eSer_13C010D000_24h + eSer_13C100D000_24h

eSer13C2D0 = eSer_13C011D000_24h + eSer_13C101D000_24h + eSer_13C110D000_24h

eSer13C3D0 = eSer_13C111D000_24h


eSer13C0D1 = eSer_13C000D010_24h

eSer13C0D1 = eSer_13C000D001_24h + eSer_13C000D100_24h

eSer13C1D1 = eSer_13C001D001_24h + eSer_13C010D001_24h + eSer_13C100D001_24h + eSer_13C001D010_24h + eSer_13C010D010_24h + eSer_13C100D010_24h + eSer_13C001D100_24h + eSer_13C010D100_24h + eSer_13C100D100_24h

eSer13C2D1 = eSer_13C011D001_24h + eSer_13C101D001_24h + eSer_13C110D001_24h + eSer_13C011D010_24h + eSer_13C101D010_24h + eSer_13C110D010_24h + eSer_13C011D100_24h + eSer_13C101D100_24h + eSer_13C110D100_24h

eSer13C3D1 = eSer_13C111D001_24h + eSer_13C111D010_24h + eSer_13C111D100_24h


eSer13C0D2 = eSer_13C000D011_24h + eSer_13C000D101_24h + eSer_13C000D110_24h
 
eSer13C1D2 = eSer_13C001D011_24h + eSer_13C010D011_24h + eSer_13C100D011_24h + eSer_13C001D101_24h + eSer_13C010D101_24h + eSer_13C100D101_24h + eSer_13C001D110_24h + eSer_13C010D110_24h + eSer_13C100D110_24h

eSer13C2D2 = eSer_13C011D011_24h + eSer_13C101D011_24h + eSer_13C110D011_24h + eSer_13C011D101_24h + eSer_13C101D101_24h + eSer_13C110D101_24h + eSer_13C011D110_24h + eSer_13C101D110_24h + eSer_13C110D110_24h

eSer13C3D2 = eSer_13C111D011_24h + eSer_13C111D101_24h + eSer_13C111D110_24h


eSer13C0D3 = eSer_13C000D111_24h

eSer13C1D3 = eSer_13C001D111_24h + eSer_13C010D111_24h + eSer_13C100D111_24h

eSer13C2D3 = eSer_13C011D111_24h + eSer_13C101D111_24h + eSer_13C110D111_24h

eSer13C3D3 = eSer_13C111D111_24h






eGly13C0D0 = eGly_13C00D00_24h

eGly13C1D0 = eGly_13C01D00_24h + eGly_13C10D00_24h

eGly13C2D0 = eGly_13C11D00_24h


eGly13C0D1 = eGly_13C00D01_24h + eGly_13C00D10_24h

eGly13C1D1 = eGly_13C01D01_24h + eGly_13C10D01_24h + eGly_13C01D10_24h + eGly_13C10D10_24h

eGly13C2D1 = eGly_13C11D01_24h + eGly_13C11D10_24h


eGly13C0D2 = eGly_13C00D11_24h

eGly13C1D2 = eGly_13C01D11_24h + eGly_13C10D11_24h

eGly13C2D2 = eGly_13C11D11_24h









% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Optimization - No.1:		13C6-GLC only
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 	%---------------------------------------------------------------------
% 	%Cancer tissue and medium
% 	%---------------------------------------------------------------------
% 	%Units: n mole  or  Enrichment (Percentage)
% 	
% 	Target_Cancer_Tissue_GLC_Isotopologues_13C0D0_24h = 2.691;
% 	Target_Cancer_Tissue_GLC_Isotopologues_13C6D0_24h = 97.309;
% 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D0_24h = 99.45;
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D0_24h = 0.10; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D1_24h = 0.03; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D2_24h = 0.01; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D2_24h = 0.15; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D3_24h = 0.07; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D3_24h = 0.19; 
% 	
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D0_24h = 99.90;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C2D0_24h = 0.10; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D2_24h = 0.00;
% 	
% 	Target_Cancer_Tissue_PRPP_Isotopologues_13C0D0_24h = 13.36;
% 	Target_Cancer_Tissue_PRPP_Isotopologues_13C5D0_24h = 86.64;
% 	
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C0D0_24h = 23.19;
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C1D0_24h = 1.01; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C4D0_24h = 2.49; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C5D0_24h = 73.30;
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C6D0_24h = 0.00;
% 	
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D0_24h = 32.64;
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D0_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D1_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D2_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D0_24h = 2.59; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D0_24h = 0.81; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D0_24h = 1.04; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D0_24h = 2.15; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D0_24h = 58.07;
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D0_24h = 2.65; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D0_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D0_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D0_24h = 0.05; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D2_24h = 0.00; 
% 	
% 	
% 	
% 	Target_Cancer_Medium_GLC_Isotopologues_13C0D0_24h = 0.658;
% 	Target_Cancer_Medium_GLC_Isotopologues_13C6D0_24h = 99.342;
% 	
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D0_24h = 99.77;
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D0_24h = 0.16; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C2D1_24h = 0.07; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D3_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D0_24h = 100.00;
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D1_24h = 0.00;  
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D2_24h = 0.00;
% 
% 
% Target_Cancer_Tissue_Ser_24h = 28.12440632;
% Target_Cancer_Tissue_Gly_24h = 228.9665356;
% 
% Target_Cancer_Medium_Ser_24h = 2395.389644;
% Target_Cancer_Medium_Gly_24h = 7019.071506;
%  
% 
% 
% 
% 
%     %---------------------------------------------------------------------
%     %Non-Cancer tissue and medium
%     %---------------------------------------------------------------------
% 	%Units: n mole  or  Enrichment (Percentage)
% 	
% 	Target_NonCancer_Tissue_GLC_Isotopologues_13C0D0_24h = 1.166;	
% 	Target_NonCancer_Tissue_GLC_Isotopologues_13C6D0_24h = 98.834;
% 	
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D0_24h = 98.97;
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D1_24h = 0.15; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D2_24h = 0.86; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D2_24h = 0.02; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D3_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D3_24h = 0.00;
% 	
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D0_24h = 100.00;
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D0_24h = 0.00;  
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C2D0_24h = 0.00;  
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D1_24h = 0.00;  
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D1_24h = 0.00;  
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D2_24h = 0.00;  
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D2_24h = 0.00;
% 	
% %	Target_NonCancer_Tissue_PRPP_Isotopologues_13C0D0_24h = 
% %	Target_NonCancer_Tissue_PRPP_Isotopologues_13C5D0_24h = 
% 	
% 	Target_NonCancer_Tissue_IMP_Isotopologues_13C0D0_24h = 40.71;
% 	Target_NonCancer_Tissue_IMP_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_IMP_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_IMP_Isotopologues_13C4D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_IMP_Isotopologues_13C5D0_24h = 48.27;
% 	Target_NonCancer_Tissue_IMP_Isotopologues_13C6D0_24h = 0.00;
% 	
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D0_24h = 31.99;
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D0_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D1_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D2_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D0_24h = 9.24; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D0_24h = 2.25; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D0_24h = 3.41; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D0_24h = 4.20; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D0_24h = 42.67;
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D0_24h = 6.14; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D0_24h = 0.10; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D2_24h = 0.00; 
% 	
% 	
% 	
% 	Target_NonCancer_Medium_GLC_Isotopologues_13C0D0_24h = 0.672;	
% 	Target_NonCancer_Medium_GLC_Isotopologues_13C6D0_24h = 99.328;
% 	
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D0_24h = 98.58;
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D0_24h = 0.07; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C2D1_24h = 0.04; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D1_24h = 1.31; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D3_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D0_24h = 100.00;
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D1_24h = 0.00;
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D2_24h = 0.00;
% 	
% 	
% Target_NonCancer_Tissue_Ser_24h = 6.568720273;
% Target_NonCancer_Tissue_Gly_24h = 41.46053008;
% 	
% Target_NonCancer_Medium_Ser_24h = 2536.651608;
% Target_NonCancer_Medium_Gly_24h = 5316.458457;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Optimization - No.2:		13C6-GLC + D3-Ser
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 	%---------------------------------------------------------------------
% 	%Cancer tissue and medium
% 	%---------------------------------------------------------------------
% 	%Units: n mole  or  Enrichment (Percentage)
% 	
% 	Target_Cancer_Tissue_GLC_Isotopologues_13C0D0_24h = 6.186;
% 	Target_Cancer_Tissue_GLC_Isotopologues_13C6D0_24h = 93.814;
% 	
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D0_24h = 48.74;
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D0_24h = 0.08; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D1_24h = 9.75; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D1_24h = 0.16; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D1_24h = 0.03; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D2_24h = 19.08;
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D2_24h = 0.44; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D2_24h = 0.03; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D3_24h = 21.68;
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D0_24h = 89.94;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C2D0_24h = 0.02; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D1_24h = 10.03;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D2_24h = 0.01; 
% 	
% 	Target_Cancer_Tissue_PRPP_Isotopologues_13C0D0_24h = 17.27;
% 	Target_Cancer_Tissue_PRPP_Isotopologues_13C5D0_24h = 82.73;
% 	
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C0D0_24h = 24.12;
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C1D0_24h = 0.84; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C4D0_24h = 2.37; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C5D0_24h = 72.67;
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C6D0_24h = 0.00; 
% 	
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D0_24h = 31.38;
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D1_24h = 1.04; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D0_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D1_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D2_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D0_24h = 3.20; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D1_24h = 0.87; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D2_24h = 0.24; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D0_24h = 0.80; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D0_24h = 0.92; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D0_24h = 1.92; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D0_24h = 56.10;
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D1_24h = 0.83; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D2_24h = 0.23; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D0_24h = 2.29; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D0_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D0_24h = 0.10; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D0_24h = 0.10; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D2_24h = 0.00; 
% 	
% 	
% 	
% 	Target_Cancer_Medium_GLC_Isotopologues_13C0D0_24h = 0.567;
% 	Target_Cancer_Medium_GLC_Isotopologues_13C6D0_24h = 99.433;
% 	
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D0_24h = 15.06;
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D0_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D1_24h = 2.07; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C2D1_24h = 0.03; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D1_24h = 0.39; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D2_24h = 7.93; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D2_24h = 0.11; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D3_24h = 74.4;
% 	Target_Cancer_Medium_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D0_24h = 98.86;
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D1_24h = 1.14;
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D2_24h = 0.00;
% 
% 
% Target_Cancer_Tissue_Ser_24h = 11.84526703;
% Target_Cancer_Tissue_Gly_24h = 171.3889876;
% 
% 
% Target_Cancer_Medium_Ser_24h = 1293.090313;
% Target_Cancer_Medium_Gly_24h = 6538.210969;
% 
% 
% 
% 
% 
%     %---------------------------------------------------------------------
%     %Non-Cancer tissue and medium
%     %---------------------------------------------------------------------
% 	%Units: n mole  or  Enrichment (Percentage)
% 	
% 	Target_NonCancer_Tissue_GLC_Isotopologues_13C0D0_24h = 1.227;	
% 	Target_NonCancer_Tissue_GLC_Isotopologues_13C6D0_24h = 98.773;
% 	
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D0_24h = 8.50; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D0_24h = 0.02; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D1_24h = 0.23; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D1_24h = 0.07; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D2_24h = 6.17; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D2_24h = 0.43; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D2_24h = 0.02; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D3_24h = 84.56;
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D0_24h = 85.53;
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D1_24h = 14.47;
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D2_24h = 0.00; 
% 	
% %	Target_NonCancer_Tissue_PRPP_Isotopologues_13C0D0_24h = 
% %	Target_NonCancer_Tissue_PRPP_Isotopologues_13C5D0_24h = 
% 	
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C0D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C1D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C2D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C4D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C5D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C6D0_24h = 
% 	
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D0_24h = 33.44;
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D1_24h = 0.20; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D2_24h = 0.22; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D0_24h =0.15;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D1_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D2_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D0_24h = 2.33; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D1_24h = 2.77; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D2_24h = 1.90; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D0_24h = 0.85; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D0_24h = 1.04; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D0_24h = 1.61; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D0_24h = 48.95;
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D1_24h = 2.65; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D2_24h = 1.82; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D0_24h = 2.09; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D2_24h = 0.00; 
% 	
% 	
% 	
% 	Target_NonCancer_Medium_GLC_Isotopologues_13C0D0_24h = 0.671;
% 	Target_NonCancer_Medium_GLC_Isotopologues_13C6D0_24h = 99.329;
% 	
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D0_24h = 2.45; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D0_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C2D1_24h = 0.32; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D2_24h = 3.35; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D2_24h = 0.13; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D3_24h = 93.75;
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D0_24h = 99.00;
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D1_24h = 1.00;
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D2_24h = 0.00;
% 
% 	
% Target_NonCancer_Tissue_Ser_24h = 6.712839201;
% Target_NonCancer_Tissue_Gly_24h = 31.86986518;
% 
% Target_NonCancer_Medium_Ser_24h = 1756.777659;
% Target_NonCancer_Medium_Gly_24h = 4663.807314;
% 	
% 	
% 
% 
% 
% 
% 
% 	
% 	
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%Optimization - No.3:		13C6-GLC + D2-Gly
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 	%---------------------------------------------------------------------
% 	%Cancer tissue and medium
% 	%---------------------------------------------------------------------
% 	%Units: n mole  or  Enrichment (Percentage)
% 	
% 	Target_Cancer_Tissue_GLC_Isotopologues_13C0D0_24h = 3.903;
% 	Target_Cancer_Tissue_GLC_Isotopologues_13C6D0_24h = 96.097;
% 	
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D0_24h = 95.91;
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D0_24h = 0.06; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D1_24h = 3.24; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D1_24h = 0.11; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D1_24h = 0.05; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D2_24h = 0.01; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C2D2_24h = 0.24; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C0D3_24h = 0.12; 
% 	Target_Cancer_Tissue_Ser_Isotopologues_13C1D3_24h = 0.26; 
% 	
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D0_24h = 94.73;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D0_24h = 0.00;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C2D0_24h = 0.01;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D1_24h = 4.49;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D1_24h = 0.00;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C0D2_24h = 0.76;
% 	Target_Cancer_Tissue_Gly_Isotopologues_13C1D2_24h = 0.00;
% 	
% 	Target_Cancer_Tissue_PRPP_Isotopologues_13C0D0_24h = 13.48;
% 	Target_Cancer_Tissue_PRPP_Isotopologues_13C5D0_24h = 86.52;
% 	
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C0D0_24h = 21.47;
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C1D0_24h = 0.23; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C4D0_24h = 1.23; 
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C5D0_24h = 77.07;
% 	Target_Cancer_Tissue_IMP_Isotopologues_13C6D0_24h = 0.00; 
% 	
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D0_24h = 30.41;
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D0_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D1_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C10D2_24h =0.00;  
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D0_24h = 3.18; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C1D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D0_24h = 0.77; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C2D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D0_24h = 0.93; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D0_24h = 2.15; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C4D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D0_24h = 59.88;
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C5D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D0_24h = 2.69; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C6D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D0_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C7D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D0_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C8D2_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D0_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D1_24h = 0.00; 
% 	Target_Cancer_Tissue_AMP_Isotopologues_13C9D2_24h = 0.00; 
% 	
% 	
% 	
% 	Target_Cancer_Medium_GLC_Isotopologues_13C0D0_24h = 0.465;
% 	Target_Cancer_Medium_GLC_Isotopologues_13C6D0_24h = 99.535;
% 	
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D0_24h = 98.84;
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D0_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D1_24h = 0.91; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C2D1_24h = 0.02; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C0D3_24h = 0.23; 
% 	Target_Cancer_Medium_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D0_24h = 42.01;
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D1_24h = 1.74;	
% 	Target_Cancer_Medium_Gly_Isotopologues_13C0D2_24h = 56.25;
% 
% 
% Target_Cancer_Tissue_Ser_24h = 23.86778678;
% Target_Cancer_Tissue_Gly_24h = 187.4162261;
% 
% 
% Target_Cancer_Medium_Ser_24h = 2584.723899; %3103.757406;
% Target_Cancer_Medium_Gly_24h = 5017.150194; %5630.933977;
% 
% 
% 
% 
% 
%     %---------------------------------------------------------------------
%     %Non-Cancer tissue and medium
%     %---------------------------------------------------------------------
% 	%Units: n mole  or  Enrichment (Percentage)
% 	
% %	Target_NonCancer_Tissue_GLC_Isotopologues = 
% 	
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D0_24h = 98.60;
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D1_24h = 0.73; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D1_24h = 0.11; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C2D2_24h = 0.47; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C0D3_24h = 0.08; 
% 	Target_NonCancer_Tissue_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D0_24h = 67.21;
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C2D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D1_24h = 14.27;
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C0D2_24h = 18.52;
% 	Target_NonCancer_Tissue_Gly_Isotopologues_13C1D2_24h = 0.00; 
% 	
% %	Target_NonCancer_Tissue_PRPP_Isotopologues_13C0D0_24h = 
% %	Target_NonCancer_Tissue_PRPP_Isotopologues_13C5D0_24h = 
% 	
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C0D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C1D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C2D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C4D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C5D0_24h = 
% %	Target_NonCancer_Tissue_IMP_Isotopologues_13C6D0_24h = 
% 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D0_24h = 36.54;
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D0_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D1_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C10D2_24h =0.00;  
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D0_24h = 3.60; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C1D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D0_24h = 0.86; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C2D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D0_24h = 1.21; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D0_24h = 1.82; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C4D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D0_24h = 53.51;
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C5D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D0_24h = 2.46; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C6D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C7D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C8D2_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D0_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D1_24h = 0.00; 
% 	Target_NonCancer_Tissue_AMP_Isotopologues_13C9D2_24h = 0.00; 
% 	
% 	
% 	
% 	Target_NonCancer_Medium_GLC_Isotopologues_13C0D0_24h = 0.634;
% 	Target_NonCancer_Medium_GLC_Isotopologues_13C6D0_24h = 99.366;
% 	
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D0_24h = 98.61;   %99.73;
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D0_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D1_24h = 1.10;    %0.20; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C2D1_24h = 0.02; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D1_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D2_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C3D2_24h = 0.00; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C0D3_24h = 0.27; 
% 	Target_NonCancer_Medium_Ser_Isotopologues_13C1D3_24h = 0.00; 
% 	
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D0_24h = 47.15;   %15.37;
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D1_24h = 1.95;    %2.12;	
% 	Target_NonCancer_Medium_Gly_Isotopologues_13C0D2_24h = 50.90;   %82.50;
% 	
% 
% Target_NonCancer_Tissue_Ser_24h = 8.004947815;
% Target_NonCancer_Tissue_Gly_24h = 24.01158668;
% 
% Target_NonCancer_Medium_Ser_24h = 2606.93612;
% Target_NonCancer_Medium_Gly_24h = 3339.796591;





end
