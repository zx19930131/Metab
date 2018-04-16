
% t0 = 0 * 60;   % unit: min
% tf = 1 * 60;   % unit: min
% dt = 0.5;       % unit: min  
% tspan = [t0 : dt : tf];
% Xd = [7, -10];
% OptimizeParameters = [72, 1, 2, 1, 1];
% 
% [t, Ydd] = ode15s(@Myfunc, tspan, Xd, [], OptimizeParameters);
% 
% Yode = Ydd(end,:);
% Yobs = [6.2333, 20.7001];
% startp = [OptimizeParameters,0.25,0.25]; % epsilon ~ MVN(0, SIGMA)
% muprior = [60,2,2,2,2,0.1,0.1];
% sigmaprior = [300,4,4,4,4,0.01,0.01];
% indexFit = [ones(1,5),0,0];
% 
% [ ps, ps_trial, chi2s, chi2s_trial, acceptance ] = MCMCMetab(Yobs, Yode, startp, muprior,sigmaprior, indexFit, 20, 10, 5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


Optimize_orig = nan(1,57);
Optimize_orig = Optimize_orig(1,:);
Optimize_orig(1) = Kf_GLC_Transport;
Optimize_orig(2) = Kr_GLC_Transport;
Optimize_orig(3) = Kf_GLC_SerPool1;
Optimize_orig(4) = Kf_GLC_PYR;
Optimize_orig(5) = Kf_GLC_PRPP ;
Optimize_orig(6) = Kf_SerPool3_Synthesis                   ;
Optimize_orig(7) = Kf_SerPool3_Degradation                 ;
Optimize_orig(8) = Kf_GlyPool3_Synthesis                   ;
Optimize_orig(9) = Kf_GlyPool3_Degradation                 ;
Optimize_orig(10) = Kf_SerPool1_GlyPool1                      ;
Optimize_orig(11) = Kr_SerPool1_GlyPool1                      ;
Optimize_orig(12) = Kf_Ser_Transport                        ;
Optimize_orig(13) = Kr_Ser_Transport                        ;
Optimize_orig(14) = Kf_Gly_Transport                        ;
Optimize_orig(15) = Kr_Gly_Transport                        ;
Optimize_orig(16) = Kf_SerPool2_SerPoolMitochon             ;
Optimize_orig(17) = Kr_SerPool2_SerPoolMitochon             ;
Optimize_orig(18) = Kf_GlyPool2_GlyPoolMitochon             ;
Optimize_orig(19) = Kr_GlyPool2_GlyPoolMitochon             ;
Optimize_orig(20) = Kf_SerPoolMitochon_GlyPoolMitochon      ;
Optimize_orig(21) = Kr_SerPoolMitochon_GlyPoolMitochon      ;
Optimize_orig(22) = Kf_GlyPool1_CO2                         ;
Optimize_orig(23) = Kr_GlyPool1_CO2                          ;
Optimize_orig(24) = Kf_MethyleneTHF_FormylTHF_Cytoplasm     ;
Optimize_orig(25) = Kr_MethyleneTHF_FormylTHF_Cytoplasm     ;
Optimize_orig(26) = Kf_FormylTHF_Formate_Cytoplasm          ;
Optimize_orig(27) = Kr_FormylTHF_Formate_Cytoplasm          ;
Optimize_orig(28) = Kf_GlyPoolMitochon_CO2                  ;
Optimize_orig(29) = Kr_GlyPoolMitochon_CO2                  ;
Optimize_orig(30) = Kf_MethyleneTHF_FormylTHF_Mitochon      ;
Optimize_orig(31) = Kr_MethyleneTHF_FormylTHF_Mitochon      ;
Optimize_orig(32) = Kf_FormylTHF_Formate_Mitochon           ;
Optimize_orig(33) = Kr_FormylTHF_Formate_Mitochon           ;
Optimize_orig(34) = Kf_PRPP1_GAR                            ;
Optimize_orig(35) = Kf_GAR_FGAR                             ;
Optimize_orig(36) = Kf_FGAR_AMP                             ;
Optimize_orig(37) = Kf_IMP_AMP                              ;
Optimize_orig(38) = Kf_AMP_Degradation                      ;
Optimize_orig(39) = Kr_MethyleneTHF_Transport               ;
Optimize_orig(40) = Kf_MethyleneTHF_Transport               ;
Optimize_orig(41) = Kr_FormylTHF_Transport                  ;
Optimize_orig(42) = Kf_FormylTHF_Transport                  ;
Optimize_orig(43) = Kr_THF_Transport                        ;
Optimize_orig(44) = Kf_THF_Transport                        ;
Optimize_orig(45) = Kr_Formate_Transport                    ;
Optimize_orig(46) = Kf_Formate_Transport                    ;
Optimize_orig(47) = Kf_PRPP2_GAR                            ;
Optimize_orig(48) = Kf_PRPP3_GAR                            ;

Optimize_orig(49) = Kf_SerPool2_GlyPool2                    ;
Optimize_orig(50) = Kr_SerPool2_GlyPool2                    ;
Optimize_orig(51) = Kf_GlyPool2_CO2                         ;
Optimize_orig(52) = Kr_GlyPool2_CO2                         ;


Optimize_orig(53) = CellVolume                                ;
Optimize_orig(54) = MediumVolume                            ;
Optimize_orig(55) = TissueVolume                            ;
Optimize_orig(56) = CytoplasmVolume                         ;
Optimize_orig(57) = MitochonVolume                          ;


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

t0 = 0 * 60;   % unit: min
tf = 24 * 60;   % unit: min
dt = 60;       % unit: min  
tspan = [t0 : dt : tf];

%options=odeset('RelTol',1e-6,'AbsTol',1e-6);




%[t, Yd_orig] = ode15s(@PurineSynthesis, tspan, Xd, [], Optimize_orig);
[t, Yd] = ode15s(@PurineSynthesis, tspan, Xd, [], OptimizeParameters); % 75s
%    if size(Yd, 1) ~= size(tspan, 2)
%         'ODEsolver failed'
%         size(Yd)
%        Yd = zeros(size(tspan, 2), size(Xd, 2));
%   end
%[t, Yd_45] = ode45(@PurineSynthesis, tspan, Xd, [], OptimizeParameters); % faster

MoleculeNumberInOneNanoMole = 6.02214129*10^14;

%Yd_orig = Yd_orig ./ MoleculeNumberInOneNanoMole;
Yd = Yd ./ MoleculeNumberInOneNanoMole;% unit is nano mole
%Yd_45 = Yd_45 ./ MoleculeNumberInOneNanoMole;

format short e;

%observed = (~Yd(end,:))==0; %exclude unobserved data with 0
observed = (~Yd(end,:))==0;
Yode = Yd(end,observed); % number of actually observed Y is 95.
%Yode = Yd_orig(end,observed); %now without data, just manipulate some
Y_var = (0.05*Yode).^2;
%startp = [OptimizeParameters(1:(end-5)),Y_var]; % epsilon ~ MVN(0, SIGMA)
startp = [OptimizeParameters(1:(end-5)),Y_var];
% muprior = [0.005,0.5,0.05,0.05,0.005,0.5e-4,0.05,0.5e-4,2.5e-3,200,5,0.005,0.05,0.0005,0.005,0.005,0.015,0.005,0.025,200,5,25,1,1.5,5,...
%     1.5,5e6,500,2,5,0.5,5,5e4,50,2.5e5,2.5e5,0.15,0.005,0.05,0.15,2.5,2.5,4,0.5,5,5,50,25,200,5,25,0.1,Y_var];
muprior = [OptimizeParameters(1:(end-5)),Y_var];
sigma_par = ([0.005,0.5,0.05,0.5,0.005,0.5e-4,0.05,0.5e-4,2.5e-3,225,5,0.005,0.05,0.0005,0.005,0.005,0.015,0.005,0.025,225,5,20,1,1.5,5,...
    1.5,5e6,500,2,5,0.5,5,5e4,50,2.5e5,2.5e5,0.15,0.005,0.05,0.15,2.5,2.5,4,0.5,5,5,50,25,225,5,25,0.1]/3).^2;
sigmaprior = [sigma_par,0.01*Y_var]; % sigma^2
%indexFit = [ones(1,52),repmat(0,1,length(Yode))];

% indexFit = [ones(1,5),repmat(0,1,47+length(Yode))];
% startp(1:5) = Optimize_orig(1:5);

indexFit = [0,0,1,1,0,repmat(0,1,47+length(Yode))];
% startp(1) = Optimize_orig(1);
% startp(2) = Optimize_orig(2);
startp(3) = Optimize_orig(3);
startp(4) = Optimize_orig(4);

ode_fun = @PurineSynthesis;
sigma = diag(Y_var);
R = chol(sigma);

matlabpool open 16

parfor seed = 1:30
rng(seed);
Yobs = Yode + randn(1,length(Yode))*R;
%[ps, ps_trial, chi2s, chi2s_trial, acceptance ] = MCMCMetab(Yobs, Yode, startp, muprior,sigmaprior, indexFit,1000, 500, 5,Xd,observed,MoleculeNumberInOneNanoMole,ode_fun);
[ps, ps_trial, chi2s, chi2s_trial, acceptance ] = Metab_GibbsWithMH(Yobs, Yode, startp, muprior,sigmaprior, indexFit,2000, 600, 5,Xd,observed,MoleculeNumberInOneNanoMole,ode_fun);

sample_mm(seed,:) = mean(ps(:,1:5));
variation_m(seed,:) = (sample_mm(seed,1:5)-OptimizeParameters(1:5))./sqrt(sigmaprior(1:5));
se_sample(seed,:) = std(ps(:,3:4)); % se of the sample of par3 and par4

conf_interv(:,:,seed) = prctile(ps(:,3:4),[5 95],1);

psTosave34 = transpose(ps(:,3:4));

name0 = 'Par3to4_sample';
name1 = strcat(name0,num2str(seed),'.txt');

fileIDs = fopen(name1,'w');
fprintf(fileIDs,'%12s %12s\n','par3','par4');
fprintf(fileIDs,'%12.6e %12.6e\n',psTosave34);
fclose(fileIDs);

end

fileID = fopen('Par3to4_stat01to30.txt','w');
fprintf(fileID,'%12s\n','samplemean');
for ii = 1:size(sample_mm,1)
    fprintf(fileID,'%12.6e\t',sample_mm(ii,:));
    fprintf(fileID,'\n');
end
fprintf(fileID,'%12s\n','variation');
for jj = 1:size(variation_m,1)
    fprintf(fileID,'%12.6e\t',variation_m(jj,:));
    fprintf(fileID,'\n');
end
fprintf(fileID,'%12s\n','sampleSD');
for kk = 1:size(se_sample,1)
    fprintf(fileID,'%12.6e\t',se_sample(kk,:));
    fprintf(fileID,'\n');
end
fclose(fileID);

fileID = fopen('Par3to4_CI01to30.txt','w');
fprintf(fileIDs,'%12s %12s %12s %12s\n','par3_low','par3_upp','par4_low','par4_upp');
for ii = 1:3
    
   fprintf(fileIDs,'%12.6e %12.6e\t',conf_interv(:,:,ii));
   fprintf(fileID,'\n');
    
end
fclose(fileID);

%type Par3to4_CI01to30.txt

%dlmwrite('test.txt',sample_mm);
%save('MyMatrix.txt', 'sample_mm', '-ascii', '-double', '-tabs');



% filename = 'Par3to4_CI01to30.txt';
% delimiterIn = ' ';
% headerlinesIn = 1;
% A = importdata(filename,delimiterIn,headerlinesIn);


% fid = fopen('Par3to4_stat01to30.txt');
% numlines = 3; % number of seeds or samples
% thissectionCell = textscan(fid, '%f %f %f %f %f', numlines, 'HeaderLines', 5, 'CollectOutput', 1);
% fclose(fid);
% celldisp(thissectionCell)
% thedataMatrix = thissectionCell{1};









% YYY = ps(:,4);
% XXX = (1:1000);
% plot(XXX,YYY)
% line([1,1000],[OptimizeParameters(4),OptimizeParameters(4)])
% line([1,1000],[sample_m(4),sample_m(4)])
% 
% XX4 = (0.2:0.04:1.12);
% post_ddd = nan(1,length(XX4));
% for i=1:length(XX4)
%     ppp = XX4(i);
%     OptimizeParameters_test = OptimizeParameters;
%     OptimizeParameters_test(4) = ppp;
%     [t, Ydd] = ode15s(@PurineSynthesis, tspan, Xd, [], OptimizeParameters_test);
%     Yode_ddd = Ydd(end,observed)./ MoleculeNumberInOneNanoMole;
%     post_ddd(1,i) = mvnpdf(Yobs, Yode_ddd, diag(Y_var))*normpdf(XX4(i),muprior(4),sigmaprior(4));
% end
% 
% plot(XX4,post_ddd)
% 
% 
% XXX3 = (0.5E-04:0.02E-04:2E-04);
% XXX4 = (0.4:0.01:1.15);
% llikli_h = nan(length(XXX3),length(XXX4));
% lpost_ddd = nan(length(XXX3),length(XXX4));
% for i=1:length(XXX4)
%     for j=1:length(XXX3)
%     ppp4 = XXX4(i);
%     ppp3 = XXX3(j);
%     OptimizeParameters_test = OptimizeParameters;
%     OptimizeParameters_test(4) = ppp4;
%     OptimizeParameters_test(3) = ppp3;
%     [t, Ydd] = ode15s(@PurineSynthesis, tspan, Xd, [], OptimizeParameters_test);
%     Yode_ddd = Ydd(end,observed)./ MoleculeNumberInOneNanoMole;
%     llikli_h(j,i) = log(mvnpdf(Yobs, Yode_ddd, diag(Y_var)));
%     lpost_ddd(j,i) = llikli_h(j,i)+log(normpdf(XXX3(j),muprior(3),sigmaprior(3)))...
%        +log(normpdf(XXX4(i),muprior(4),sigmaprior(4)));
%     [i,j]
%     end
% end
% 
% surf(XXX3,XXX4,lpost_ddd);
% 
% [max_pc,I_c] = max(lpost_ddd);
% [max_pr,I_r] = max(transpose(lpost_ddd));
% [max_t,I_t] = max(max_pc); %52 col
% [max_tt,I_tt] = max(max_pr); %44 row
% 
% surf(XXX3,XXX4,llikli_h); % results are cool
% 



