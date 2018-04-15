function [dxdt] = PurineSynthesis(t, Xd,	OptimizeParameters)
% -----------------------------------------------------------------------------------------------------------------------
    CO2 = 1;
    NH3 = 1;
    cCO2 = 1;
    cNH3 = 1;
    mCO2 = 1;
    mNH3 = 1;

    AvogadroConstant = 6.02214129E+23;


Kf_GLC_Transport						= OptimizeParameters(1);
Kr_GLC_Transport                        = OptimizeParameters(2);           
Kf_GLC_SerPool1                         = OptimizeParameters(3);           
Kf_GLC_PYR                              = OptimizeParameters(4);           
Kf_GLC_PRPP                             = OptimizeParameters(5);           
Kf_SerPool3_Synthesis                   = OptimizeParameters(6);           
Kf_SerPool3_Degradation                 = OptimizeParameters(7);           
Kf_GlyPool3_Synthesis                   = OptimizeParameters(8);           
Kf_GlyPool3_Degradation                 = OptimizeParameters(9);           
Kf_SerPool1_GlyPool1                      = OptimizeParameters(10);        
Kr_SerPool1_GlyPool1                      = OptimizeParameters(11);        
Kf_Ser_Transport                        = OptimizeParameters(12);          
Kr_Ser_Transport                        = OptimizeParameters(13);          
Kf_Gly_Transport                        = OptimizeParameters(14);          
Kr_Gly_Transport                        = OptimizeParameters(15);          
Kf_SerPool2_SerPoolMitochon             = OptimizeParameters(16);          
Kr_SerPool2_SerPoolMitochon             = OptimizeParameters(17);          
Kf_GlyPool2_GlyPoolMitochon             = OptimizeParameters(18);          
Kr_GlyPool2_GlyPoolMitochon             = OptimizeParameters(19);          
Kf_SerPoolMitochon_GlyPoolMitochon      = OptimizeParameters(20);          
Kr_SerPoolMitochon_GlyPoolMitochon      = OptimizeParameters(21);          
Kf_GlyPool1_CO2                          = OptimizeParameters(22);         
Kr_GlyPool1_CO2                          = OptimizeParameters(23);         
Kf_MethyleneTHF_FormylTHF_Cytoplasm     = OptimizeParameters(24);          
Kr_MethyleneTHF_FormylTHF_Cytoplasm     = OptimizeParameters(25);          
Kf_FormylTHF_Formate_Cytoplasm          = OptimizeParameters(26);          
Kr_FormylTHF_Formate_Cytoplasm          = OptimizeParameters(27);          
Kf_GlyPoolMitochon_CO2                  = OptimizeParameters(28);          
Kr_GlyPoolMitochon_CO2                  = OptimizeParameters(29);          
Kf_MethyleneTHF_FormylTHF_Mitochon      = OptimizeParameters(30);          
Kr_MethyleneTHF_FormylTHF_Mitochon      = OptimizeParameters(31);          
Kf_FormylTHF_Formate_Mitochon           = OptimizeParameters(32);          
Kr_FormylTHF_Formate_Mitochon           = OptimizeParameters(33);          
Kf_PRPP1_GAR                            = OptimizeParameters(34);          
Kf_GAR_FGAR                             = OptimizeParameters(35);          
Kf_FGAR_AMP                             = OptimizeParameters(36);          
Kf_IMP_AMP                              = OptimizeParameters(37);          
Kf_AMP_Degradation                      = OptimizeParameters(38);          
Kr_MethyleneTHF_Transport               = OptimizeParameters(39);          
Kf_MethyleneTHF_Transport               = OptimizeParameters(40);          
Kr_FormylTHF_Transport                  = OptimizeParameters(41);          
Kf_FormylTHF_Transport                  = OptimizeParameters(42);          
Kr_THF_Transport                        = OptimizeParameters(43);          
Kf_THF_Transport                        = OptimizeParameters(44);          
Kr_Formate_Transport                    = OptimizeParameters(45);          
Kf_Formate_Transport                    = OptimizeParameters(46);          
Kf_PRPP2_GAR                            = OptimizeParameters(47);          
Kf_PRPP3_GAR                            = OptimizeParameters(48);          

Kf_SerPool2_GlyPool2										= OptimizeParameters(49); 
Kr_SerPool2_GlyPool2                    = OptimizeParameters(50); 
Kf_GlyPool2_CO2                         = OptimizeParameters(51); 
Kr_GlyPool2_CO2                         = OptimizeParameters(52); 


CellVolume															= OptimizeParameters(53);   
MediumVolume                            = OptimizeParameters(54);   
TissueVolume                            = OptimizeParameters(55);   
CytoplasmVolume                         = OptimizeParameters(56);   
MitochonVolume                          = OptimizeParameters(57);   


%%%%%%%%%%%%%%%%%%%%%%%%
	Kf_SerPool3_GlyPool3 = Kf_SerPool2_GlyPool2;
	Kr_SerPool3_GlyPool3 = Kr_SerPool2_GlyPool2;
	Kf_GlyPool3_CO2 = Kf_GlyPool2_CO2;
	Kr_GlyPool3_CO2 = Kr_GlyPool2_CO2;
	
	Kf_Ser_Transport3 = Kf_Ser_Transport;
	Kr_Ser_Transport3 = Kr_Ser_Transport;
	Kf_Gly_Transport3 = Kf_Gly_Transport;
	Kr_Gly_Transport3 = Kr_Gly_Transport;
%%%%%%%%%%%%%%%%%%%%%%%%%





eGLC_13C0					   = Xd(1  );
eGLC_13C1                      = Xd(2  );
cGLC_13C0                      = Xd(3  );
cGLC_13C1                      = Xd(4  );
eSer_13C000D000                = Xd(5  );
eSer_13C001D000                = Xd(6  );
eSer_13C010D000                = Xd(7  );
eSer_13C011D000                = Xd(8  );
eSer_13C100D000                = Xd(9  );
eSer_13C101D000                = Xd(10 );
eSer_13C110D000                = Xd(11 );
eSer_13C111D000                = Xd(12 );
eSer_13C000D001                = Xd(13 );
eSer_13C001D001                = Xd(14 );
eSer_13C010D001                = Xd(15 );
eSer_13C011D001                = Xd(16 );
eSer_13C100D001                = Xd(17 );
eSer_13C101D001                = Xd(18 );
eSer_13C110D001                = Xd(19 );
eSer_13C111D001                = Xd(20 );
eSer_13C000D010                = Xd(21 );
eSer_13C001D010                = Xd(22 );
eSer_13C010D010                = Xd(23 );
eSer_13C011D010                = Xd(24 );
eSer_13C100D010                = Xd(25 );
eSer_13C101D010                = Xd(26 );
eSer_13C110D010                = Xd(27 );
eSer_13C111D010                = Xd(28 );
eSer_13C000D011                = Xd(29 );
eSer_13C001D011                = Xd(30 );
eSer_13C010D011                = Xd(31 );
eSer_13C011D011                = Xd(32 );
eSer_13C100D011                = Xd(33 );
eSer_13C101D011                = Xd(34 );
eSer_13C110D011                = Xd(35 );
eSer_13C111D011                = Xd(36 );
eSer_13C000D100                = Xd(37 );
eSer_13C001D100                = Xd(38 );
eSer_13C010D100                = Xd(39 );
eSer_13C011D100                = Xd(40 );
eSer_13C100D100                = Xd(41 );
eSer_13C101D100                = Xd(42 );
eSer_13C110D100                = Xd(43 );
eSer_13C111D100                = Xd(44 );
eSer_13C000D101                = Xd(45 );
eSer_13C001D101                = Xd(46 );
eSer_13C010D101                = Xd(47 );
eSer_13C011D101                = Xd(48 );
eSer_13C100D101                = Xd(49 );
eSer_13C101D101                = Xd(50 );
eSer_13C110D101                = Xd(51 );
eSer_13C111D101                = Xd(52 );
eSer_13C000D110                = Xd(53 );
eSer_13C001D110                = Xd(54 );
eSer_13C010D110                = Xd(55 );
eSer_13C011D110                = Xd(56 );
eSer_13C100D110                = Xd(57 );
eSer_13C101D110                = Xd(58 );
eSer_13C110D110                = Xd(59 );
eSer_13C111D110                = Xd(60 );
eSer_13C000D111                = Xd(61 );
eSer_13C001D111                = Xd(62 );
eSer_13C010D111                = Xd(63 );
eSer_13C011D111                = Xd(64 );
eSer_13C100D111                = Xd(65 );
eSer_13C101D111                = Xd(66 );
eSer_13C110D111                = Xd(67 );
eSer_13C111D111                = Xd(68 );
eGly_13C00D00                  = Xd(69 );
eGly_13C01D00                  = Xd(70 );
eGly_13C10D00                  = Xd(71 );
eGly_13C11D00                  = Xd(72 );
eGly_13C00D01                  = Xd(73 );
eGly_13C01D01                  = Xd(74 );
eGly_13C10D01                  = Xd(75 );
eGly_13C11D01                  = Xd(76 );
eGly_13C00D10                  = Xd(77 );
eGly_13C01D10                  = Xd(78 );
eGly_13C10D10                  = Xd(79 );
eGly_13C11D10                  = Xd(80 );
eGly_13C00D11                  = Xd(81 );
eGly_13C01D11                  = Xd(82 );
eGly_13C10D11                  = Xd(83 );
eGly_13C11D11                  = Xd(84 );
cSerPool1_13C000D000           = Xd(85 );
cSerPool1_13C001D000           = Xd(86 );
cSerPool1_13C010D000           = Xd(87 );
cSerPool1_13C011D000           = Xd(88 );
cSerPool1_13C100D000           = Xd(89 );
cSerPool1_13C101D000           = Xd(90 );
cSerPool1_13C110D000           = Xd(91 );
cSerPool1_13C111D000           = Xd(92 );
cSerPool1_13C000D001           = Xd(93 );
cSerPool1_13C001D001           = Xd(94 );
cSerPool1_13C010D001           = Xd(95 );
cSerPool1_13C011D001           = Xd(96 );
cSerPool1_13C100D001           = Xd(97 );
cSerPool1_13C101D001           = Xd(98 );
cSerPool1_13C110D001           = Xd(99 );
cSerPool1_13C111D001           = Xd(100);
cSerPool1_13C000D010           = Xd(101);
cSerPool1_13C001D010           = Xd(102);
cSerPool1_13C010D010           = Xd(103);
cSerPool1_13C011D010           = Xd(104);
cSerPool1_13C100D010           = Xd(105);
cSerPool1_13C101D010           = Xd(106);
cSerPool1_13C110D010           = Xd(107);
cSerPool1_13C111D010           = Xd(108);
cSerPool1_13C000D011           = Xd(109);
cSerPool1_13C001D011           = Xd(110);
cSerPool1_13C010D011           = Xd(111);
cSerPool1_13C011D011           = Xd(112);
cSerPool1_13C100D011           = Xd(113);
cSerPool1_13C101D011           = Xd(114);
cSerPool1_13C110D011           = Xd(115);
cSerPool1_13C111D011           = Xd(116);
cSerPool1_13C000D100           = Xd(117);
cSerPool1_13C001D100           = Xd(118);
cSerPool1_13C010D100           = Xd(119);
cSerPool1_13C011D100           = Xd(120);
cSerPool1_13C100D100           = Xd(121);
cSerPool1_13C101D100           = Xd(122);
cSerPool1_13C110D100           = Xd(123);
cSerPool1_13C111D100           = Xd(124);
cSerPool1_13C000D101           = Xd(125);
cSerPool1_13C001D101           = Xd(126);
cSerPool1_13C010D101           = Xd(127);
cSerPool1_13C011D101           = Xd(128);
cSerPool1_13C100D101           = Xd(129);
cSerPool1_13C101D101           = Xd(130);
cSerPool1_13C110D101           = Xd(131);
cSerPool1_13C111D101           = Xd(132);
cSerPool1_13C000D110           = Xd(133);
cSerPool1_13C001D110           = Xd(134);
cSerPool1_13C010D110           = Xd(135);
cSerPool1_13C011D110           = Xd(136);
cSerPool1_13C100D110           = Xd(137);
cSerPool1_13C101D110           = Xd(138);
cSerPool1_13C110D110           = Xd(139);
cSerPool1_13C111D110           = Xd(140);
cSerPool1_13C000D111           = Xd(141);
cSerPool1_13C001D111           = Xd(142);
cSerPool1_13C010D111           = Xd(143);
cSerPool1_13C011D111           = Xd(144);
cSerPool1_13C100D111           = Xd(145);
cSerPool1_13C101D111           = Xd(146);
cSerPool1_13C110D111           = Xd(147);
cSerPool1_13C111D111           = Xd(148);
cGlyPool1_13C00D00             = Xd(149);
cGlyPool1_13C01D00             = Xd(150);
cGlyPool1_13C10D00             = Xd(151);
cGlyPool1_13C11D00             = Xd(152);
cGlyPool1_13C00D01             = Xd(153);
cGlyPool1_13C01D01             = Xd(154);
cGlyPool1_13C10D01             = Xd(155);
cGlyPool1_13C11D01             = Xd(156);
cGlyPool1_13C00D10             = Xd(157);
cGlyPool1_13C01D10             = Xd(158);
cGlyPool1_13C10D10             = Xd(159);
cGlyPool1_13C11D10             = Xd(160);
cGlyPool1_13C00D11             = Xd(161);
cGlyPool1_13C01D11             = Xd(162);
cGlyPool1_13C10D11             = Xd(163);
cGlyPool1_13C11D11             = Xd(164);
cSerPool2_13C000D000           = Xd(165);
cSerPool2_13C001D000           = Xd(166);
cSerPool2_13C010D000           = Xd(167);
cSerPool2_13C011D000           = Xd(168);
cSerPool2_13C100D000           = Xd(169);
cSerPool2_13C101D000           = Xd(170);
cSerPool2_13C110D000           = Xd(171);
cSerPool2_13C111D000           = Xd(172);
cSerPool2_13C000D001           = Xd(173);
cSerPool2_13C001D001           = Xd(174);
cSerPool2_13C010D001           = Xd(175);
cSerPool2_13C011D001           = Xd(176);
cSerPool2_13C100D001           = Xd(177);
cSerPool2_13C101D001           = Xd(178);
cSerPool2_13C110D001           = Xd(179);
cSerPool2_13C111D001           = Xd(180);
cSerPool2_13C000D010           = Xd(181);
cSerPool2_13C001D010           = Xd(182);
cSerPool2_13C010D010           = Xd(183);
cSerPool2_13C011D010           = Xd(184);
cSerPool2_13C100D010           = Xd(185);
cSerPool2_13C101D010           = Xd(186);
cSerPool2_13C110D010           = Xd(187);
cSerPool2_13C111D010           = Xd(188);
cSerPool2_13C000D011           = Xd(189);
cSerPool2_13C001D011           = Xd(190);
cSerPool2_13C010D011           = Xd(191);
cSerPool2_13C011D011           = Xd(192);
cSerPool2_13C100D011           = Xd(193);
cSerPool2_13C101D011           = Xd(194);
cSerPool2_13C110D011           = Xd(195);
cSerPool2_13C111D011           = Xd(196);
cSerPool2_13C000D100           = Xd(197);
cSerPool2_13C001D100           = Xd(198);
cSerPool2_13C010D100           = Xd(199);
cSerPool2_13C011D100           = Xd(200);
cSerPool2_13C100D100           = Xd(201);
cSerPool2_13C101D100           = Xd(202);
cSerPool2_13C110D100           = Xd(203);
cSerPool2_13C111D100           = Xd(204);
cSerPool2_13C000D101           = Xd(205);
cSerPool2_13C001D101           = Xd(206);
cSerPool2_13C010D101           = Xd(207);
cSerPool2_13C011D101           = Xd(208);
cSerPool2_13C100D101           = Xd(209);
cSerPool2_13C101D101           = Xd(210);
cSerPool2_13C110D101           = Xd(211);
cSerPool2_13C111D101           = Xd(212);
cSerPool2_13C000D110           = Xd(213);
cSerPool2_13C001D110           = Xd(214);
cSerPool2_13C010D110           = Xd(215);
cSerPool2_13C011D110           = Xd(216);
cSerPool2_13C100D110           = Xd(217);
cSerPool2_13C101D110           = Xd(218);
cSerPool2_13C110D110           = Xd(219);
cSerPool2_13C111D110           = Xd(220);
cSerPool2_13C000D111           = Xd(221);
cSerPool2_13C001D111           = Xd(222);
cSerPool2_13C010D111           = Xd(223);
cSerPool2_13C011D111           = Xd(224);
cSerPool2_13C100D111           = Xd(225);
cSerPool2_13C101D111           = Xd(226);
cSerPool2_13C110D111           = Xd(227);
cSerPool2_13C111D111           = Xd(228);
cGlyPool2_13C00D00             = Xd(229);
cGlyPool2_13C01D00             = Xd(230);
cGlyPool2_13C10D00             = Xd(231);
cGlyPool2_13C11D00             = Xd(232);
cGlyPool2_13C00D01             = Xd(233);
cGlyPool2_13C01D01             = Xd(234);
cGlyPool2_13C10D01             = Xd(235);
cGlyPool2_13C11D01             = Xd(236);
cGlyPool2_13C00D10             = Xd(237);
cGlyPool2_13C01D10             = Xd(238);
cGlyPool2_13C10D10             = Xd(239);
cGlyPool2_13C11D10             = Xd(240);
cGlyPool2_13C00D11             = Xd(241);
cGlyPool2_13C01D11             = Xd(242);
cGlyPool2_13C10D11             = Xd(243);
cGlyPool2_13C11D11             = Xd(244);
cSerPool3_13C000D000           = Xd(245);
cSerPool3_13C001D000           = Xd(246);
cSerPool3_13C010D000           = Xd(247);
cSerPool3_13C011D000           = Xd(248);
cSerPool3_13C100D000           = Xd(249);
cSerPool3_13C101D000           = Xd(250);
cSerPool3_13C110D000           = Xd(251);
cSerPool3_13C111D000           = Xd(252);
cSerPool3_13C000D001           = Xd(253);
cSerPool3_13C001D001           = Xd(254);
cSerPool3_13C010D001           = Xd(255);
cSerPool3_13C011D001           = Xd(256);
cSerPool3_13C100D001           = Xd(257);
cSerPool3_13C101D001           = Xd(258);
cSerPool3_13C110D001           = Xd(259);
cSerPool3_13C111D001           = Xd(260);
cSerPool3_13C000D010           = Xd(261);
cSerPool3_13C001D010           = Xd(262);
cSerPool3_13C010D010           = Xd(263);
cSerPool3_13C011D010           = Xd(264);
cSerPool3_13C100D010           = Xd(265);
cSerPool3_13C101D010           = Xd(266);
cSerPool3_13C110D010           = Xd(267);
cSerPool3_13C111D010           = Xd(268);
cSerPool3_13C000D011           = Xd(269);
cSerPool3_13C001D011           = Xd(270);
cSerPool3_13C010D011           = Xd(271);
cSerPool3_13C011D011           = Xd(272);
cSerPool3_13C100D011           = Xd(273);
cSerPool3_13C101D011           = Xd(274);
cSerPool3_13C110D011           = Xd(275);
cSerPool3_13C111D011           = Xd(276);
cSerPool3_13C000D100           = Xd(277);
cSerPool3_13C001D100           = Xd(278);
cSerPool3_13C010D100           = Xd(279);
cSerPool3_13C011D100           = Xd(280);
cSerPool3_13C100D100           = Xd(281);
cSerPool3_13C101D100           = Xd(282);
cSerPool3_13C110D100           = Xd(283);
cSerPool3_13C111D100           = Xd(284);
cSerPool3_13C000D101           = Xd(285);
cSerPool3_13C001D101           = Xd(286);
cSerPool3_13C010D101           = Xd(287);
cSerPool3_13C011D101           = Xd(288);
cSerPool3_13C100D101           = Xd(289);
cSerPool3_13C101D101           = Xd(290);
cSerPool3_13C110D101           = Xd(291);
cSerPool3_13C111D101           = Xd(292);
cSerPool3_13C000D110           = Xd(293);
cSerPool3_13C001D110           = Xd(294);
cSerPool3_13C010D110           = Xd(295);
cSerPool3_13C011D110           = Xd(296);
cSerPool3_13C100D110           = Xd(297);
cSerPool3_13C101D110           = Xd(298);
cSerPool3_13C110D110           = Xd(299);
cSerPool3_13C111D110           = Xd(300);
cSerPool3_13C000D111           = Xd(301);
cSerPool3_13C001D111           = Xd(302);
cSerPool3_13C010D111           = Xd(303);
cSerPool3_13C011D111           = Xd(304);
cSerPool3_13C100D111           = Xd(305);
cSerPool3_13C101D111           = Xd(306);
cSerPool3_13C110D111           = Xd(307);
cSerPool3_13C111D111           = Xd(308);
cGlyPool3_13C00D00             = Xd(309);
cGlyPool3_13C01D00             = Xd(310);
cGlyPool3_13C10D00             = Xd(311);
cGlyPool3_13C11D00             = Xd(312);
cGlyPool3_13C00D01             = Xd(313);
cGlyPool3_13C01D01             = Xd(314);
cGlyPool3_13C10D01             = Xd(315);
cGlyPool3_13C11D01             = Xd(316);
cGlyPool3_13C00D10             = Xd(317);
cGlyPool3_13C01D10             = Xd(318);
cGlyPool3_13C10D10             = Xd(319);
cGlyPool3_13C11D10             = Xd(320);
cGlyPool3_13C00D11             = Xd(321);
cGlyPool3_13C01D11             = Xd(322);
cGlyPool3_13C10D11             = Xd(323);
cGlyPool3_13C11D11             = Xd(324);
mSerPool_13C000D000            = Xd(325);
mSerPool_13C001D000            = Xd(326);
mSerPool_13C010D000            = Xd(327);
mSerPool_13C011D000            = Xd(328);
mSerPool_13C100D000            = Xd(329);
mSerPool_13C101D000            = Xd(330);
mSerPool_13C110D000            = Xd(331);
mSerPool_13C111D000            = Xd(332);
mSerPool_13C000D001            = Xd(333);
mSerPool_13C001D001            = Xd(334);
mSerPool_13C010D001            = Xd(335);
mSerPool_13C011D001            = Xd(336);
mSerPool_13C100D001            = Xd(337);
mSerPool_13C101D001            = Xd(338);
mSerPool_13C110D001            = Xd(339);
mSerPool_13C111D001            = Xd(340);
mSerPool_13C000D010            = Xd(341);
mSerPool_13C001D010            = Xd(342);
mSerPool_13C010D010            = Xd(343);
mSerPool_13C011D010            = Xd(344);
mSerPool_13C100D010            = Xd(345);
mSerPool_13C101D010            = Xd(346);
mSerPool_13C110D010            = Xd(347);
mSerPool_13C111D010            = Xd(348);
mSerPool_13C000D011            = Xd(349);
mSerPool_13C001D011            = Xd(350);
mSerPool_13C010D011            = Xd(351);
mSerPool_13C011D011            = Xd(352);
mSerPool_13C100D011            = Xd(353);
mSerPool_13C101D011            = Xd(354);
mSerPool_13C110D011            = Xd(355);
mSerPool_13C111D011            = Xd(356);
mSerPool_13C000D100            = Xd(357);
mSerPool_13C001D100            = Xd(358);
mSerPool_13C010D100            = Xd(359);
mSerPool_13C011D100            = Xd(360);
mSerPool_13C100D100            = Xd(361);
mSerPool_13C101D100            = Xd(362);
mSerPool_13C110D100            = Xd(363);
mSerPool_13C111D100            = Xd(364);
mSerPool_13C000D101            = Xd(365);
mSerPool_13C001D101            = Xd(366);
mSerPool_13C010D101            = Xd(367);
mSerPool_13C011D101            = Xd(368);
mSerPool_13C100D101            = Xd(369);
mSerPool_13C101D101            = Xd(370);
mSerPool_13C110D101            = Xd(371);
mSerPool_13C111D101            = Xd(372);
mSerPool_13C000D110            = Xd(373);
mSerPool_13C001D110            = Xd(374);
mSerPool_13C010D110            = Xd(375);
mSerPool_13C011D110            = Xd(376);
mSerPool_13C100D110            = Xd(377);
mSerPool_13C101D110            = Xd(378);
mSerPool_13C110D110            = Xd(379);
mSerPool_13C111D110            = Xd(380);
mSerPool_13C000D111            = Xd(381);
mSerPool_13C001D111            = Xd(382);
mSerPool_13C010D111            = Xd(383);
mSerPool_13C011D111            = Xd(384);
mSerPool_13C100D111            = Xd(385);
mSerPool_13C101D111            = Xd(386);
mSerPool_13C110D111            = Xd(387);
mSerPool_13C111D111            = Xd(388);
mGlyPool_13C00D00              = Xd(389);
mGlyPool_13C01D00              = Xd(390);
mGlyPool_13C10D00              = Xd(391);
mGlyPool_13C11D00              = Xd(392);
mGlyPool_13C00D01              = Xd(393);
mGlyPool_13C01D01              = Xd(394);
mGlyPool_13C10D01              = Xd(395);
mGlyPool_13C11D01              = Xd(396);
mGlyPool_13C00D10              = Xd(397);
mGlyPool_13C01D10              = Xd(398);
mGlyPool_13C10D10              = Xd(399);
mGlyPool_13C11D10              = Xd(400);
mGlyPool_13C00D11              = Xd(401);
mGlyPool_13C01D11              = Xd(402);
mGlyPool_13C10D11              = Xd(403);
mGlyPool_13C11D11              = Xd(404);
cMethyleneTHF_13C0D00          = Xd(405);
cMethyleneTHF_13C1D00          = Xd(406);
cMethyleneTHF_13C0D01          = Xd(407);
cMethyleneTHF_13C1D01          = Xd(408);
cMethyleneTHF_13C0D10          = Xd(409);
cMethyleneTHF_13C1D10          = Xd(410);
cMethyleneTHF_13C0D11          = Xd(411);
cMethyleneTHF_13C1D11          = Xd(412);
cFormylTHF_13C0D0              = Xd(413);
cFormylTHF_13C1D0              = Xd(414);
cFormylTHF_13C0D1              = Xd(415);
cFormylTHF_13C1D1              = Xd(416);
cFormate_13C0D0                = Xd(417);
cFormate_13C1D0                = Xd(418);
cFormate_13C0D1                = Xd(419);
cFormate_13C1D1                = Xd(420);
cTHF                           = Xd(421);
mMethyleneTHF_13C0D00          = Xd(422);
mMethyleneTHF_13C1D00          = Xd(423);
mMethyleneTHF_13C0D01          = Xd(424);
mMethyleneTHF_13C1D01          = Xd(425);
mMethyleneTHF_13C0D10          = Xd(426);
mMethyleneTHF_13C1D10          = Xd(427);
mMethyleneTHF_13C0D11          = Xd(428);
mMethyleneTHF_13C1D11          = Xd(429);
mFormylTHF_13C0D0              = Xd(430);
mFormylTHF_13C1D0              = Xd(431);
mFormylTHF_13C0D1              = Xd(432);
mFormylTHF_13C1D1              = Xd(433);
mFormate_13C0D0                = Xd(434);
mFormate_13C1D0                = Xd(435);
mFormate_13C0D1                = Xd(436);
mFormate_13C1D1                = Xd(437);
mTHF                           = Xd(438);
cPRPP_13C0                     = Xd(439);
cPRPP_13C1                     = Xd(440);
cGAR_13C000                    = Xd(441);
cGAR_13C001                    = Xd(442);
cGAR_13C010                    = Xd(443);
cGAR_13C011                    = Xd(444);
cGAR_13C100                    = Xd(445);
cGAR_13C101                    = Xd(446);
cGAR_13C110                    = Xd(447);
cGAR_13C111                    = Xd(448);
cFGAR_13C0000D0                = Xd(449);
cFGAR_13C0001D0                = Xd(450);
cFGAR_13C0010D0                = Xd(451);
cFGAR_13C0011D0                = Xd(452);
cFGAR_13C0100D0                = Xd(453);
cFGAR_13C0101D0                = Xd(454);
cFGAR_13C0110D0                = Xd(455);
cFGAR_13C0111D0                = Xd(456);
cFGAR_13C1000D0                = Xd(457);
cFGAR_13C1001D0                = Xd(458);
cFGAR_13C1010D0                = Xd(459);
cFGAR_13C1011D0                = Xd(460);
cFGAR_13C1100D0                = Xd(461);
cFGAR_13C1101D0                = Xd(462);
cFGAR_13C1110D0                = Xd(463);
cFGAR_13C1111D0                = Xd(464);
cFGAR_13C0000D1                = Xd(465);
cFGAR_13C0001D1                = Xd(466);
cFGAR_13C0010D1                = Xd(467);
cFGAR_13C0011D1                = Xd(468);
cFGAR_13C0100D1                = Xd(469);
cFGAR_13C0101D1                = Xd(470);
cFGAR_13C0110D1                = Xd(471);
cFGAR_13C0111D1                = Xd(472);
cFGAR_13C1000D1                = Xd(473);
cFGAR_13C1001D1                = Xd(474);
cFGAR_13C1010D1                = Xd(475);
cFGAR_13C1011D1                = Xd(476);
cFGAR_13C1100D1                = Xd(477);
cFGAR_13C1101D1                = Xd(478);
cFGAR_13C1110D1                = Xd(479);
cFGAR_13C1111D1                = Xd(480);
cAMP_13C00000D00               = Xd(481);
cAMP_13C00001D00               = Xd(482);
cAMP_13C00010D00               = Xd(483);
cAMP_13C00011D00               = Xd(484);
cAMP_13C00100D00               = Xd(485);
cAMP_13C00101D00               = Xd(486);
cAMP_13C00110D00               = Xd(487);
cAMP_13C00111D00               = Xd(488);
cAMP_13C01000D00               = Xd(489);
cAMP_13C01001D00               = Xd(490);
cAMP_13C01010D00               = Xd(491);
cAMP_13C01011D00               = Xd(492);
cAMP_13C01100D00               = Xd(493);
cAMP_13C01101D00               = Xd(494);
cAMP_13C01110D00               = Xd(495);
cAMP_13C01111D00               = Xd(496);
cAMP_13C10000D00               = Xd(497);
cAMP_13C10001D00               = Xd(498);
cAMP_13C10010D00               = Xd(499);
cAMP_13C10011D00               = Xd(500);
cAMP_13C10100D00               = Xd(501);
cAMP_13C10101D00               = Xd(502);
cAMP_13C10110D00               = Xd(503);
cAMP_13C10111D00               = Xd(504);
cAMP_13C11000D00               = Xd(505);
cAMP_13C11001D00               = Xd(506);
cAMP_13C11010D00               = Xd(507);
cAMP_13C11011D00               = Xd(508);
cAMP_13C11100D00               = Xd(509);
cAMP_13C11101D00               = Xd(510);
cAMP_13C11110D00               = Xd(511);
cAMP_13C11111D00               = Xd(512);
cAMP_13C00000D01               = Xd(513);
cAMP_13C00001D01               = Xd(514);
cAMP_13C00010D01               = Xd(515);
cAMP_13C00011D01               = Xd(516);
cAMP_13C00100D01               = Xd(517);
cAMP_13C00101D01               = Xd(518);
cAMP_13C00110D01               = Xd(519);
cAMP_13C00111D01               = Xd(520);
cAMP_13C01000D01               = Xd(521);
cAMP_13C01001D01               = Xd(522);
cAMP_13C01010D01               = Xd(523);
cAMP_13C01011D01               = Xd(524);
cAMP_13C01100D01               = Xd(525);
cAMP_13C01101D01               = Xd(526);
cAMP_13C01110D01               = Xd(527);
cAMP_13C01111D01               = Xd(528);
cAMP_13C10000D01               = Xd(529);
cAMP_13C10001D01               = Xd(530);
cAMP_13C10010D01               = Xd(531);
cAMP_13C10011D01               = Xd(532);
cAMP_13C10100D01               = Xd(533);
cAMP_13C10101D01               = Xd(534);
cAMP_13C10110D01               = Xd(535);
cAMP_13C10111D01               = Xd(536);
cAMP_13C11000D01               = Xd(537);
cAMP_13C11001D01               = Xd(538);
cAMP_13C11010D01               = Xd(539);
cAMP_13C11011D01               = Xd(540);
cAMP_13C11100D01               = Xd(541);
cAMP_13C11101D01               = Xd(542);
cAMP_13C11110D01               = Xd(543);
cAMP_13C11111D01               = Xd(544);
cAMP_13C00000D10               = Xd(545);
cAMP_13C00001D10               = Xd(546);
cAMP_13C00010D10               = Xd(547);
cAMP_13C00011D10               = Xd(548);
cAMP_13C00100D10               = Xd(549);
cAMP_13C00101D10               = Xd(550);
cAMP_13C00110D10               = Xd(551);
cAMP_13C00111D10               = Xd(552);
cAMP_13C01000D10               = Xd(553);
cAMP_13C01001D10               = Xd(554);
cAMP_13C01010D10               = Xd(555);
cAMP_13C01011D10               = Xd(556);
cAMP_13C01100D10               = Xd(557);
cAMP_13C01101D10               = Xd(558);
cAMP_13C01110D10               = Xd(559);
cAMP_13C01111D10               = Xd(560);
cAMP_13C10000D10               = Xd(561);
cAMP_13C10001D10               = Xd(562);
cAMP_13C10010D10               = Xd(563);
cAMP_13C10011D10               = Xd(564);
cAMP_13C10100D10               = Xd(565);
cAMP_13C10101D10               = Xd(566);
cAMP_13C10110D10               = Xd(567);
cAMP_13C10111D10               = Xd(568);
cAMP_13C11000D10               = Xd(569);
cAMP_13C11001D10               = Xd(570);
cAMP_13C11010D10               = Xd(571);
cAMP_13C11011D10               = Xd(572);
cAMP_13C11100D10               = Xd(573);
cAMP_13C11101D10               = Xd(574);
cAMP_13C11110D10               = Xd(575);
cAMP_13C11111D10               = Xd(576);
cAMP_13C00000D11               = Xd(577);
cAMP_13C00001D11               = Xd(578);
cAMP_13C00010D11               = Xd(579);
cAMP_13C00011D11               = Xd(580);
cAMP_13C00100D11               = Xd(581);
cAMP_13C00101D11               = Xd(582);
cAMP_13C00110D11               = Xd(583);
cAMP_13C00111D11               = Xd(584);
cAMP_13C01000D11               = Xd(585);
cAMP_13C01001D11               = Xd(586);
cAMP_13C01010D11               = Xd(587);
cAMP_13C01011D11               = Xd(588);
cAMP_13C01100D11               = Xd(589);
cAMP_13C01101D11               = Xd(590);
cAMP_13C01110D11               = Xd(591);
cAMP_13C01111D11               = Xd(592);
cAMP_13C10000D11               = Xd(593);
cAMP_13C10001D11               = Xd(594);
cAMP_13C10010D11               = Xd(595);
cAMP_13C10011D11               = Xd(596);
cAMP_13C10100D11               = Xd(597);
cAMP_13C10101D11               = Xd(598);
cAMP_13C10110D11               = Xd(599);
cAMP_13C10111D11               = Xd(600);
cAMP_13C11000D11               = Xd(601);
cAMP_13C11001D11               = Xd(602);
cAMP_13C11010D11               = Xd(603);
cAMP_13C11011D11               = Xd(604);
cAMP_13C11100D11               = Xd(605);
cAMP_13C11101D11               = Xd(606);
cAMP_13C11110D11               = Xd(607);
cAMP_13C11111D11               = Xd(608);











    Vf_Vin_GLC_13C0_13C0 = Kf_GLC_Transport * eGLC_13C0;
    Vf_Vin_GLC_13C1_13C1 = Kf_GLC_Transport * eGLC_13C1;



    Vr_Vin_GLC_13C0_13C0 = Kr_GLC_Transport * cGLC_13C0;
    Vr_Vin_GLC_13C1_13C1 = Kr_GLC_Transport * cGLC_13C1;



    Vf_GLC_SerPool1_13C0_13C0D000 = Kf_GLC_SerPool1 * cGLC_13C0;
    Vf_GLC_SerPool1_13C1_13C3D000 = Kf_GLC_SerPool1 * cGLC_13C1;



    Vf_GLC_PYR_13C0 = Kf_GLC_PYR * cGLC_13C0;
    Vf_GLC_PYR_13C1 = Kf_GLC_PYR * cGLC_13C1;



    Vf_GLC_PRPP_13C0_13C0 = Kf_GLC_PRPP * cGLC_13C0;
    Vf_GLC_PRPP_13C1_13C1 = Kf_GLC_PRPP * cGLC_13C1;



    Vf_Vin_Ser_13C000D000_13C000D000 = Kf_Ser_Transport * eSer_13C000D000;
    Vf_Vin_Ser_13C001D000_13C001D000 = Kf_Ser_Transport * eSer_13C001D000;
    Vf_Vin_Ser_13C010D000_13C010D000 = Kf_Ser_Transport * eSer_13C010D000;
    Vf_Vin_Ser_13C011D000_13C011D000 = Kf_Ser_Transport * eSer_13C011D000;
    Vf_Vin_Ser_13C100D000_13C100D000 = Kf_Ser_Transport * eSer_13C100D000;
    Vf_Vin_Ser_13C101D000_13C101D000 = Kf_Ser_Transport * eSer_13C101D000;
    Vf_Vin_Ser_13C110D000_13C110D000 = Kf_Ser_Transport * eSer_13C110D000;
    Vf_Vin_Ser_13C111D000_13C111D000 = Kf_Ser_Transport * eSer_13C111D000;
    Vf_Vin_Ser_13C000D001_13C000D001 = Kf_Ser_Transport * eSer_13C000D001;
    Vf_Vin_Ser_13C001D001_13C001D001 = Kf_Ser_Transport * eSer_13C001D001;
    Vf_Vin_Ser_13C010D001_13C010D001 = Kf_Ser_Transport * eSer_13C010D001;
    Vf_Vin_Ser_13C011D001_13C011D001 = Kf_Ser_Transport * eSer_13C011D001;
    Vf_Vin_Ser_13C100D001_13C100D001 = Kf_Ser_Transport * eSer_13C100D001;
    Vf_Vin_Ser_13C101D001_13C101D001 = Kf_Ser_Transport * eSer_13C101D001;
    Vf_Vin_Ser_13C110D001_13C110D001 = Kf_Ser_Transport * eSer_13C110D001;
    Vf_Vin_Ser_13C111D001_13C111D001 = Kf_Ser_Transport * eSer_13C111D001;
    Vf_Vin_Ser_13C000D010_13C000D010 = Kf_Ser_Transport * eSer_13C000D010; 
    Vf_Vin_Ser_13C001D010_13C001D010 = Kf_Ser_Transport * eSer_13C001D010; 
    Vf_Vin_Ser_13C010D010_13C010D010 = Kf_Ser_Transport * eSer_13C010D010; 
    Vf_Vin_Ser_13C011D010_13C011D010 = Kf_Ser_Transport * eSer_13C011D010; 
    Vf_Vin_Ser_13C100D010_13C100D010 = Kf_Ser_Transport * eSer_13C100D010; 
    Vf_Vin_Ser_13C101D010_13C101D010 = Kf_Ser_Transport * eSer_13C101D010; 
    Vf_Vin_Ser_13C110D010_13C110D010 = Kf_Ser_Transport * eSer_13C110D010; 
    Vf_Vin_Ser_13C111D010_13C111D010 = Kf_Ser_Transport * eSer_13C111D010; 
    Vf_Vin_Ser_13C000D011_13C000D011 = Kf_Ser_Transport * eSer_13C000D011; 
    Vf_Vin_Ser_13C001D011_13C001D011 = Kf_Ser_Transport * eSer_13C001D011; 
    Vf_Vin_Ser_13C010D011_13C010D011 = Kf_Ser_Transport * eSer_13C010D011; 
    Vf_Vin_Ser_13C011D011_13C011D011 = Kf_Ser_Transport * eSer_13C011D011; 
    Vf_Vin_Ser_13C100D011_13C100D011 = Kf_Ser_Transport * eSer_13C100D011; 
    Vf_Vin_Ser_13C101D011_13C101D011 = Kf_Ser_Transport * eSer_13C101D011; 
    Vf_Vin_Ser_13C110D011_13C110D011 = Kf_Ser_Transport * eSer_13C110D011; 
    Vf_Vin_Ser_13C111D011_13C111D011 = Kf_Ser_Transport * eSer_13C111D011; 
    Vf_Vin_Ser_13C000D100_13C000D100 = Kf_Ser_Transport * eSer_13C000D100; 
    Vf_Vin_Ser_13C001D100_13C001D100 = Kf_Ser_Transport * eSer_13C001D100; 
    Vf_Vin_Ser_13C010D100_13C010D100 = Kf_Ser_Transport * eSer_13C010D100; 
    Vf_Vin_Ser_13C011D100_13C011D100 = Kf_Ser_Transport * eSer_13C011D100; 
    Vf_Vin_Ser_13C100D100_13C100D100 = Kf_Ser_Transport * eSer_13C100D100; 
    Vf_Vin_Ser_13C101D100_13C101D100 = Kf_Ser_Transport * eSer_13C101D100; 
    Vf_Vin_Ser_13C110D100_13C110D100 = Kf_Ser_Transport * eSer_13C110D100; 
    Vf_Vin_Ser_13C111D100_13C111D100 = Kf_Ser_Transport * eSer_13C111D100; 
    Vf_Vin_Ser_13C000D101_13C000D101 = Kf_Ser_Transport * eSer_13C000D101; 
    Vf_Vin_Ser_13C001D101_13C001D101 = Kf_Ser_Transport * eSer_13C001D101; 
    Vf_Vin_Ser_13C010D101_13C010D101 = Kf_Ser_Transport * eSer_13C010D101; 
    Vf_Vin_Ser_13C011D101_13C011D101 = Kf_Ser_Transport * eSer_13C011D101; 
    Vf_Vin_Ser_13C100D101_13C100D101 = Kf_Ser_Transport * eSer_13C100D101; 
    Vf_Vin_Ser_13C101D101_13C101D101 = Kf_Ser_Transport * eSer_13C101D101; 
    Vf_Vin_Ser_13C110D101_13C110D101 = Kf_Ser_Transport * eSer_13C110D101; 
    Vf_Vin_Ser_13C111D101_13C111D101 = Kf_Ser_Transport * eSer_13C111D101; 
    Vf_Vin_Ser_13C000D110_13C000D110 = Kf_Ser_Transport * eSer_13C000D110; 
    Vf_Vin_Ser_13C001D110_13C001D110 = Kf_Ser_Transport * eSer_13C001D110; 
    Vf_Vin_Ser_13C010D110_13C010D110 = Kf_Ser_Transport * eSer_13C010D110; 
    Vf_Vin_Ser_13C011D110_13C011D110 = Kf_Ser_Transport * eSer_13C011D110; 
    Vf_Vin_Ser_13C100D110_13C100D110 = Kf_Ser_Transport * eSer_13C100D110; 
    Vf_Vin_Ser_13C101D110_13C101D110 = Kf_Ser_Transport * eSer_13C101D110; 
    Vf_Vin_Ser_13C110D110_13C110D110 = Kf_Ser_Transport * eSer_13C110D110; 
    Vf_Vin_Ser_13C111D110_13C111D110 = Kf_Ser_Transport * eSer_13C111D110; 
    Vf_Vin_Ser_13C000D111_13C000D111 = Kf_Ser_Transport * eSer_13C000D111; 
    Vf_Vin_Ser_13C001D111_13C001D111 = Kf_Ser_Transport * eSer_13C001D111; 
    Vf_Vin_Ser_13C010D111_13C010D111 = Kf_Ser_Transport * eSer_13C010D111; 
    Vf_Vin_Ser_13C011D111_13C011D111 = Kf_Ser_Transport * eSer_13C011D111; 
    Vf_Vin_Ser_13C100D111_13C100D111 = Kf_Ser_Transport * eSer_13C100D111; 
    Vf_Vin_Ser_13C101D111_13C101D111 = Kf_Ser_Transport * eSer_13C101D111; 
    Vf_Vin_Ser_13C110D111_13C110D111 = Kf_Ser_Transport * eSer_13C110D111; 
    Vf_Vin_Ser_13C111D111_13C111D111 = Kf_Ser_Transport * eSer_13C111D111; 



    Vr_Vin_Ser_13C000D000_13C000D000 = Kr_Ser_Transport * cSerPool2_13C000D000;
    Vr_Vin_Ser_13C001D000_13C001D000 = Kr_Ser_Transport * cSerPool2_13C001D000;      
    Vr_Vin_Ser_13C010D000_13C010D000 = Kr_Ser_Transport * cSerPool2_13C010D000;  
    Vr_Vin_Ser_13C011D000_13C011D000 = Kr_Ser_Transport * cSerPool2_13C011D000;  
    Vr_Vin_Ser_13C100D000_13C100D000 = Kr_Ser_Transport * cSerPool2_13C100D000;   
    Vr_Vin_Ser_13C101D000_13C101D000 = Kr_Ser_Transport * cSerPool2_13C101D000;   
    Vr_Vin_Ser_13C110D000_13C110D000 = Kr_Ser_Transport * cSerPool2_13C110D000;   
    Vr_Vin_Ser_13C111D000_13C111D000 = Kr_Ser_Transport * cSerPool2_13C111D000;   
    Vr_Vin_Ser_13C000D001_13C000D001 = Kr_Ser_Transport * cSerPool2_13C000D001;   
    Vr_Vin_Ser_13C001D001_13C001D001 = Kr_Ser_Transport * cSerPool2_13C001D001;   
    Vr_Vin_Ser_13C010D001_13C010D001 = Kr_Ser_Transport * cSerPool2_13C010D001;   
    Vr_Vin_Ser_13C011D001_13C011D001 = Kr_Ser_Transport * cSerPool2_13C011D001;   
    Vr_Vin_Ser_13C100D001_13C100D001 = Kr_Ser_Transport * cSerPool2_13C100D001;   
    Vr_Vin_Ser_13C101D001_13C101D001 = Kr_Ser_Transport * cSerPool2_13C101D001;   
    Vr_Vin_Ser_13C110D001_13C110D001 = Kr_Ser_Transport * cSerPool2_13C110D001;   
    Vr_Vin_Ser_13C111D001_13C111D001 = Kr_Ser_Transport * cSerPool2_13C111D001;   
    Vr_Vin_Ser_13C000D010_13C000D010 = Kr_Ser_Transport * cSerPool2_13C000D010;
    Vr_Vin_Ser_13C001D010_13C001D010 = Kr_Ser_Transport * cSerPool2_13C001D010;      
    Vr_Vin_Ser_13C010D010_13C010D010 = Kr_Ser_Transport * cSerPool2_13C010D010;  
    Vr_Vin_Ser_13C011D010_13C011D010 = Kr_Ser_Transport * cSerPool2_13C011D010;  
    Vr_Vin_Ser_13C100D010_13C100D010 = Kr_Ser_Transport * cSerPool2_13C100D010;   
    Vr_Vin_Ser_13C101D010_13C101D010 = Kr_Ser_Transport * cSerPool2_13C101D010;   
    Vr_Vin_Ser_13C110D010_13C110D010 = Kr_Ser_Transport * cSerPool2_13C110D010;   
    Vr_Vin_Ser_13C111D010_13C111D010 = Kr_Ser_Transport * cSerPool2_13C111D010;   
    Vr_Vin_Ser_13C000D011_13C000D011 = Kr_Ser_Transport * cSerPool2_13C000D011;   
    Vr_Vin_Ser_13C001D011_13C001D011 = Kr_Ser_Transport * cSerPool2_13C001D011;   
    Vr_Vin_Ser_13C010D011_13C010D011 = Kr_Ser_Transport * cSerPool2_13C010D011;   
    Vr_Vin_Ser_13C011D011_13C011D011 = Kr_Ser_Transport * cSerPool2_13C011D011;   
    Vr_Vin_Ser_13C100D011_13C100D011 = Kr_Ser_Transport * cSerPool2_13C100D011;   
    Vr_Vin_Ser_13C101D011_13C101D011 = Kr_Ser_Transport * cSerPool2_13C101D011;   
    Vr_Vin_Ser_13C110D011_13C110D011 = Kr_Ser_Transport * cSerPool2_13C110D011;   
    Vr_Vin_Ser_13C111D011_13C111D011 = Kr_Ser_Transport * cSerPool2_13C111D011;   
    Vr_Vin_Ser_13C000D100_13C000D100 = Kr_Ser_Transport * cSerPool2_13C000D100;
    Vr_Vin_Ser_13C001D100_13C001D100 = Kr_Ser_Transport * cSerPool2_13C001D100;      
    Vr_Vin_Ser_13C010D100_13C010D100 = Kr_Ser_Transport * cSerPool2_13C010D100;  
    Vr_Vin_Ser_13C011D100_13C011D100 = Kr_Ser_Transport * cSerPool2_13C011D100;  
    Vr_Vin_Ser_13C100D100_13C100D100 = Kr_Ser_Transport * cSerPool2_13C100D100;   
    Vr_Vin_Ser_13C101D100_13C101D100 = Kr_Ser_Transport * cSerPool2_13C101D100;   
    Vr_Vin_Ser_13C110D100_13C110D100 = Kr_Ser_Transport * cSerPool2_13C110D100;   
    Vr_Vin_Ser_13C111D100_13C111D100 = Kr_Ser_Transport * cSerPool2_13C111D100;   
    Vr_Vin_Ser_13C000D101_13C000D101 = Kr_Ser_Transport * cSerPool2_13C000D101;   
    Vr_Vin_Ser_13C001D101_13C001D101 = Kr_Ser_Transport * cSerPool2_13C001D101;   
    Vr_Vin_Ser_13C010D101_13C010D101 = Kr_Ser_Transport * cSerPool2_13C010D101;   
    Vr_Vin_Ser_13C011D101_13C011D101 = Kr_Ser_Transport * cSerPool2_13C011D101;   
    Vr_Vin_Ser_13C100D101_13C100D101 = Kr_Ser_Transport * cSerPool2_13C100D101;   
    Vr_Vin_Ser_13C101D101_13C101D101 = Kr_Ser_Transport * cSerPool2_13C101D101;   
    Vr_Vin_Ser_13C110D101_13C110D101 = Kr_Ser_Transport * cSerPool2_13C110D101;   
    Vr_Vin_Ser_13C111D101_13C111D101 = Kr_Ser_Transport * cSerPool2_13C111D101;   
    Vr_Vin_Ser_13C000D110_13C000D110 = Kr_Ser_Transport * cSerPool2_13C000D110;
    Vr_Vin_Ser_13C001D110_13C001D110 = Kr_Ser_Transport * cSerPool2_13C001D110;      
    Vr_Vin_Ser_13C010D110_13C010D110 = Kr_Ser_Transport * cSerPool2_13C010D110;  
    Vr_Vin_Ser_13C011D110_13C011D110 = Kr_Ser_Transport * cSerPool2_13C011D110;  
    Vr_Vin_Ser_13C100D110_13C100D110 = Kr_Ser_Transport * cSerPool2_13C100D110;   
    Vr_Vin_Ser_13C101D110_13C101D110 = Kr_Ser_Transport * cSerPool2_13C101D110;   
    Vr_Vin_Ser_13C110D110_13C110D110 = Kr_Ser_Transport * cSerPool2_13C110D110;   
    Vr_Vin_Ser_13C111D110_13C111D110 = Kr_Ser_Transport * cSerPool2_13C111D110;   
    Vr_Vin_Ser_13C000D111_13C000D111 = Kr_Ser_Transport * cSerPool2_13C000D111;   
    Vr_Vin_Ser_13C001D111_13C001D111 = Kr_Ser_Transport * cSerPool2_13C001D111;   
    Vr_Vin_Ser_13C010D111_13C010D111 = Kr_Ser_Transport * cSerPool2_13C010D111;   
    Vr_Vin_Ser_13C011D111_13C011D111 = Kr_Ser_Transport * cSerPool2_13C011D111;   
    Vr_Vin_Ser_13C100D111_13C100D111 = Kr_Ser_Transport * cSerPool2_13C100D111;   
    Vr_Vin_Ser_13C101D111_13C101D111 = Kr_Ser_Transport * cSerPool2_13C101D111;   
    Vr_Vin_Ser_13C110D111_13C110D111 = Kr_Ser_Transport * cSerPool2_13C110D111;   
    Vr_Vin_Ser_13C111D111_13C111D111 = Kr_Ser_Transport * cSerPool2_13C111D111;   



    Vf_Vin_Ser3_13C000D000_13C000D000 = Kf_Ser_Transport3 * eSer_13C000D000;
    Vf_Vin_Ser3_13C001D000_13C001D000 = Kf_Ser_Transport3 * eSer_13C001D000;
    Vf_Vin_Ser3_13C010D000_13C010D000 = Kf_Ser_Transport3 * eSer_13C010D000;
    Vf_Vin_Ser3_13C011D000_13C011D000 = Kf_Ser_Transport3 * eSer_13C011D000;
    Vf_Vin_Ser3_13C100D000_13C100D000 = Kf_Ser_Transport3 * eSer_13C100D000;
    Vf_Vin_Ser3_13C101D000_13C101D000 = Kf_Ser_Transport3 * eSer_13C101D000;
    Vf_Vin_Ser3_13C110D000_13C110D000 = Kf_Ser_Transport3 * eSer_13C110D000;
    Vf_Vin_Ser3_13C111D000_13C111D000 = Kf_Ser_Transport3 * eSer_13C111D000;
    Vf_Vin_Ser3_13C000D001_13C000D001 = Kf_Ser_Transport3 * eSer_13C000D001;
    Vf_Vin_Ser3_13C001D001_13C001D001 = Kf_Ser_Transport3 * eSer_13C001D001;
    Vf_Vin_Ser3_13C010D001_13C010D001 = Kf_Ser_Transport3 * eSer_13C010D001;
    Vf_Vin_Ser3_13C011D001_13C011D001 = Kf_Ser_Transport3 * eSer_13C011D001;
    Vf_Vin_Ser3_13C100D001_13C100D001 = Kf_Ser_Transport3 * eSer_13C100D001;
    Vf_Vin_Ser3_13C101D001_13C101D001 = Kf_Ser_Transport3 * eSer_13C101D001;
    Vf_Vin_Ser3_13C110D001_13C110D001 = Kf_Ser_Transport3 * eSer_13C110D001;
    Vf_Vin_Ser3_13C111D001_13C111D001 = Kf_Ser_Transport3 * eSer_13C111D001;
    Vf_Vin_Ser3_13C000D010_13C000D010 = Kf_Ser_Transport3 * eSer_13C000D010;
    Vf_Vin_Ser3_13C001D010_13C001D010 = Kf_Ser_Transport3 * eSer_13C001D010;
    Vf_Vin_Ser3_13C010D010_13C010D010 = Kf_Ser_Transport3 * eSer_13C010D010;
    Vf_Vin_Ser3_13C011D010_13C011D010 = Kf_Ser_Transport3 * eSer_13C011D010;
    Vf_Vin_Ser3_13C100D010_13C100D010 = Kf_Ser_Transport3 * eSer_13C100D010;
    Vf_Vin_Ser3_13C101D010_13C101D010 = Kf_Ser_Transport3 * eSer_13C101D010;
    Vf_Vin_Ser3_13C110D010_13C110D010 = Kf_Ser_Transport3 * eSer_13C110D010;
    Vf_Vin_Ser3_13C111D010_13C111D010 = Kf_Ser_Transport3 * eSer_13C111D010;
    Vf_Vin_Ser3_13C000D011_13C000D011 = Kf_Ser_Transport3 * eSer_13C000D011;
    Vf_Vin_Ser3_13C001D011_13C001D011 = Kf_Ser_Transport3 * eSer_13C001D011;
    Vf_Vin_Ser3_13C010D011_13C010D011 = Kf_Ser_Transport3 * eSer_13C010D011;
    Vf_Vin_Ser3_13C011D011_13C011D011 = Kf_Ser_Transport3 * eSer_13C011D011;
    Vf_Vin_Ser3_13C100D011_13C100D011 = Kf_Ser_Transport3 * eSer_13C100D011;
    Vf_Vin_Ser3_13C101D011_13C101D011 = Kf_Ser_Transport3 * eSer_13C101D011;
    Vf_Vin_Ser3_13C110D011_13C110D011 = Kf_Ser_Transport3 * eSer_13C110D011;
    Vf_Vin_Ser3_13C111D011_13C111D011 = Kf_Ser_Transport3 * eSer_13C111D011;
    Vf_Vin_Ser3_13C000D100_13C000D100 = Kf_Ser_Transport3 * eSer_13C000D100;
    Vf_Vin_Ser3_13C001D100_13C001D100 = Kf_Ser_Transport3 * eSer_13C001D100;
    Vf_Vin_Ser3_13C010D100_13C010D100 = Kf_Ser_Transport3 * eSer_13C010D100;
    Vf_Vin_Ser3_13C011D100_13C011D100 = Kf_Ser_Transport3 * eSer_13C011D100;
    Vf_Vin_Ser3_13C100D100_13C100D100 = Kf_Ser_Transport3 * eSer_13C100D100;
    Vf_Vin_Ser3_13C101D100_13C101D100 = Kf_Ser_Transport3 * eSer_13C101D100;
    Vf_Vin_Ser3_13C110D100_13C110D100 = Kf_Ser_Transport3 * eSer_13C110D100;
    Vf_Vin_Ser3_13C111D100_13C111D100 = Kf_Ser_Transport3 * eSer_13C111D100;
    Vf_Vin_Ser3_13C000D101_13C000D101 = Kf_Ser_Transport3 * eSer_13C000D101;
    Vf_Vin_Ser3_13C001D101_13C001D101 = Kf_Ser_Transport3 * eSer_13C001D101;
    Vf_Vin_Ser3_13C010D101_13C010D101 = Kf_Ser_Transport3 * eSer_13C010D101;
    Vf_Vin_Ser3_13C011D101_13C011D101 = Kf_Ser_Transport3 * eSer_13C011D101;
    Vf_Vin_Ser3_13C100D101_13C100D101 = Kf_Ser_Transport3 * eSer_13C100D101;
    Vf_Vin_Ser3_13C101D101_13C101D101 = Kf_Ser_Transport3 * eSer_13C101D101;
    Vf_Vin_Ser3_13C110D101_13C110D101 = Kf_Ser_Transport3 * eSer_13C110D101;
    Vf_Vin_Ser3_13C111D101_13C111D101 = Kf_Ser_Transport3 * eSer_13C111D101;
    Vf_Vin_Ser3_13C000D110_13C000D110 = Kf_Ser_Transport3 * eSer_13C000D110;
    Vf_Vin_Ser3_13C001D110_13C001D110 = Kf_Ser_Transport3 * eSer_13C001D110;
    Vf_Vin_Ser3_13C010D110_13C010D110 = Kf_Ser_Transport3 * eSer_13C010D110;
    Vf_Vin_Ser3_13C011D110_13C011D110 = Kf_Ser_Transport3 * eSer_13C011D110;
    Vf_Vin_Ser3_13C100D110_13C100D110 = Kf_Ser_Transport3 * eSer_13C100D110;
    Vf_Vin_Ser3_13C101D110_13C101D110 = Kf_Ser_Transport3 * eSer_13C101D110;
    Vf_Vin_Ser3_13C110D110_13C110D110 = Kf_Ser_Transport3 * eSer_13C110D110;
    Vf_Vin_Ser3_13C111D110_13C111D110 = Kf_Ser_Transport3 * eSer_13C111D110;
    Vf_Vin_Ser3_13C000D111_13C000D111 = Kf_Ser_Transport3 * eSer_13C000D111;
    Vf_Vin_Ser3_13C001D111_13C001D111 = Kf_Ser_Transport3 * eSer_13C001D111;
    Vf_Vin_Ser3_13C010D111_13C010D111 = Kf_Ser_Transport3 * eSer_13C010D111;
    Vf_Vin_Ser3_13C011D111_13C011D111 = Kf_Ser_Transport3 * eSer_13C011D111;
    Vf_Vin_Ser3_13C100D111_13C100D111 = Kf_Ser_Transport3 * eSer_13C100D111;
    Vf_Vin_Ser3_13C101D111_13C101D111 = Kf_Ser_Transport3 * eSer_13C101D111;
    Vf_Vin_Ser3_13C110D111_13C110D111 = Kf_Ser_Transport3 * eSer_13C110D111;
    Vf_Vin_Ser3_13C111D111_13C111D111 = Kf_Ser_Transport3 * eSer_13C111D111;



    Vr_Vin_Ser3_13C000D000_13C000D000 = Kr_Ser_Transport3 * cSerPool3_13C000D000;
    Vr_Vin_Ser3_13C001D000_13C001D000 = Kr_Ser_Transport3 * cSerPool3_13C001D000;      
    Vr_Vin_Ser3_13C010D000_13C010D000 = Kr_Ser_Transport3 * cSerPool3_13C010D000;  
    Vr_Vin_Ser3_13C011D000_13C011D000 = Kr_Ser_Transport3 * cSerPool3_13C011D000;  
    Vr_Vin_Ser3_13C100D000_13C100D000 = Kr_Ser_Transport3 * cSerPool3_13C100D000;   
    Vr_Vin_Ser3_13C101D000_13C101D000 = Kr_Ser_Transport3 * cSerPool3_13C101D000;   
    Vr_Vin_Ser3_13C110D000_13C110D000 = Kr_Ser_Transport3 * cSerPool3_13C110D000;   
    Vr_Vin_Ser3_13C111D000_13C111D000 = Kr_Ser_Transport3 * cSerPool3_13C111D000;   
    Vr_Vin_Ser3_13C000D001_13C000D001 = Kr_Ser_Transport3 * cSerPool3_13C000D001;   
    Vr_Vin_Ser3_13C001D001_13C001D001 = Kr_Ser_Transport3 * cSerPool3_13C001D001;   
    Vr_Vin_Ser3_13C010D001_13C010D001 = Kr_Ser_Transport3 * cSerPool3_13C010D001;   
    Vr_Vin_Ser3_13C011D001_13C011D001 = Kr_Ser_Transport3 * cSerPool3_13C011D001;   
    Vr_Vin_Ser3_13C100D001_13C100D001 = Kr_Ser_Transport3 * cSerPool3_13C100D001;   
    Vr_Vin_Ser3_13C101D001_13C101D001 = Kr_Ser_Transport3 * cSerPool3_13C101D001;   
    Vr_Vin_Ser3_13C110D001_13C110D001 = Kr_Ser_Transport3 * cSerPool3_13C110D001;   
    Vr_Vin_Ser3_13C111D001_13C111D001 = Kr_Ser_Transport3 * cSerPool3_13C111D001;   
    Vr_Vin_Ser3_13C000D010_13C000D010 = Kr_Ser_Transport3 * cSerPool3_13C000D010;
    Vr_Vin_Ser3_13C001D010_13C001D010 = Kr_Ser_Transport3 * cSerPool3_13C001D010;      
    Vr_Vin_Ser3_13C010D010_13C010D010 = Kr_Ser_Transport3 * cSerPool3_13C010D010;  
    Vr_Vin_Ser3_13C011D010_13C011D010 = Kr_Ser_Transport3 * cSerPool3_13C011D010;  
    Vr_Vin_Ser3_13C100D010_13C100D010 = Kr_Ser_Transport3 * cSerPool3_13C100D010;   
    Vr_Vin_Ser3_13C101D010_13C101D010 = Kr_Ser_Transport3 * cSerPool3_13C101D010;   
    Vr_Vin_Ser3_13C110D010_13C110D010 = Kr_Ser_Transport3 * cSerPool3_13C110D010;   
    Vr_Vin_Ser3_13C111D010_13C111D010 = Kr_Ser_Transport3 * cSerPool3_13C111D010;   
    Vr_Vin_Ser3_13C000D011_13C000D011 = Kr_Ser_Transport3 * cSerPool3_13C000D011;   
    Vr_Vin_Ser3_13C001D011_13C001D011 = Kr_Ser_Transport3 * cSerPool3_13C001D011;   
    Vr_Vin_Ser3_13C010D011_13C010D011 = Kr_Ser_Transport3 * cSerPool3_13C010D011;   
    Vr_Vin_Ser3_13C011D011_13C011D011 = Kr_Ser_Transport3 * cSerPool3_13C011D011;   
    Vr_Vin_Ser3_13C100D011_13C100D011 = Kr_Ser_Transport3 * cSerPool3_13C100D011;   
    Vr_Vin_Ser3_13C101D011_13C101D011 = Kr_Ser_Transport3 * cSerPool3_13C101D011;   
    Vr_Vin_Ser3_13C110D011_13C110D011 = Kr_Ser_Transport3 * cSerPool3_13C110D011;   
    Vr_Vin_Ser3_13C111D011_13C111D011 = Kr_Ser_Transport3 * cSerPool3_13C111D011;   
    Vr_Vin_Ser3_13C000D100_13C000D100 = Kr_Ser_Transport3 * cSerPool3_13C000D100;
    Vr_Vin_Ser3_13C001D100_13C001D100 = Kr_Ser_Transport3 * cSerPool3_13C001D100;      
    Vr_Vin_Ser3_13C010D100_13C010D100 = Kr_Ser_Transport3 * cSerPool3_13C010D100;  
    Vr_Vin_Ser3_13C011D100_13C011D100 = Kr_Ser_Transport3 * cSerPool3_13C011D100;  
    Vr_Vin_Ser3_13C100D100_13C100D100 = Kr_Ser_Transport3 * cSerPool3_13C100D100;   
    Vr_Vin_Ser3_13C101D100_13C101D100 = Kr_Ser_Transport3 * cSerPool3_13C101D100;   
    Vr_Vin_Ser3_13C110D100_13C110D100 = Kr_Ser_Transport3 * cSerPool3_13C110D100;   
    Vr_Vin_Ser3_13C111D100_13C111D100 = Kr_Ser_Transport3 * cSerPool3_13C111D100;   
    Vr_Vin_Ser3_13C000D101_13C000D101 = Kr_Ser_Transport3 * cSerPool3_13C000D101;   
    Vr_Vin_Ser3_13C001D101_13C001D101 = Kr_Ser_Transport3 * cSerPool3_13C001D101;   
    Vr_Vin_Ser3_13C010D101_13C010D101 = Kr_Ser_Transport3 * cSerPool3_13C010D101;   
    Vr_Vin_Ser3_13C011D101_13C011D101 = Kr_Ser_Transport3 * cSerPool3_13C011D101;   
    Vr_Vin_Ser3_13C100D101_13C100D101 = Kr_Ser_Transport3 * cSerPool3_13C100D101;   
    Vr_Vin_Ser3_13C101D101_13C101D101 = Kr_Ser_Transport3 * cSerPool3_13C101D101;   
    Vr_Vin_Ser3_13C110D101_13C110D101 = Kr_Ser_Transport3 * cSerPool3_13C110D101;   
    Vr_Vin_Ser3_13C111D101_13C111D101 = Kr_Ser_Transport3 * cSerPool3_13C111D101;   
    Vr_Vin_Ser3_13C000D110_13C000D110 = Kr_Ser_Transport3 * cSerPool3_13C000D110;
    Vr_Vin_Ser3_13C001D110_13C001D110 = Kr_Ser_Transport3 * cSerPool3_13C001D110;      
    Vr_Vin_Ser3_13C010D110_13C010D110 = Kr_Ser_Transport3 * cSerPool3_13C010D110;  
    Vr_Vin_Ser3_13C011D110_13C011D110 = Kr_Ser_Transport3 * cSerPool3_13C011D110;  
    Vr_Vin_Ser3_13C100D110_13C100D110 = Kr_Ser_Transport3 * cSerPool3_13C100D110;   
    Vr_Vin_Ser3_13C101D110_13C101D110 = Kr_Ser_Transport3 * cSerPool3_13C101D110;   
    Vr_Vin_Ser3_13C110D110_13C110D110 = Kr_Ser_Transport3 * cSerPool3_13C110D110;   
    Vr_Vin_Ser3_13C111D110_13C111D110 = Kr_Ser_Transport3 * cSerPool3_13C111D110;   
    Vr_Vin_Ser3_13C000D111_13C000D111 = Kr_Ser_Transport3 * cSerPool3_13C000D111;   
    Vr_Vin_Ser3_13C001D111_13C001D111 = Kr_Ser_Transport3 * cSerPool3_13C001D111;   
    Vr_Vin_Ser3_13C010D111_13C010D111 = Kr_Ser_Transport3 * cSerPool3_13C010D111;   
    Vr_Vin_Ser3_13C011D111_13C011D111 = Kr_Ser_Transport3 * cSerPool3_13C011D111;   
    Vr_Vin_Ser3_13C100D111_13C100D111 = Kr_Ser_Transport3 * cSerPool3_13C100D111;   
    Vr_Vin_Ser3_13C101D111_13C101D111 = Kr_Ser_Transport3 * cSerPool3_13C101D111;   
    Vr_Vin_Ser3_13C110D111_13C110D111 = Kr_Ser_Transport3 * cSerPool3_13C110D111;   
    Vr_Vin_Ser3_13C111D111_13C111D111 = Kr_Ser_Transport3 * cSerPool3_13C111D111;   



   Vf_Vin_Gly_13C00D00_13C00D00 = Kf_Gly_Transport * eGly_13C00D00;
   Vf_Vin_Gly_13C01D00_13C01D00 = Kf_Gly_Transport * eGly_13C01D00;
   Vf_Vin_Gly_13C10D00_13C10D00 = Kf_Gly_Transport * eGly_13C10D00;
   Vf_Vin_Gly_13C11D00_13C11D00 = Kf_Gly_Transport * eGly_13C11D00;
   Vf_Vin_Gly_13C00D01_13C00D01 = Kf_Gly_Transport * eGly_13C00D01;
   Vf_Vin_Gly_13C01D01_13C01D01 = Kf_Gly_Transport * eGly_13C01D01;
   Vf_Vin_Gly_13C10D01_13C10D01 = Kf_Gly_Transport * eGly_13C10D01;
   Vf_Vin_Gly_13C11D01_13C11D01 = Kf_Gly_Transport * eGly_13C11D01;
   Vf_Vin_Gly_13C00D10_13C00D10 = Kf_Gly_Transport * eGly_13C00D10;
   Vf_Vin_Gly_13C01D10_13C01D10 = Kf_Gly_Transport * eGly_13C01D10;
   Vf_Vin_Gly_13C10D10_13C10D10 = Kf_Gly_Transport * eGly_13C10D10;
   Vf_Vin_Gly_13C11D10_13C11D10 = Kf_Gly_Transport * eGly_13C11D10;
   Vf_Vin_Gly_13C00D11_13C00D11 = Kf_Gly_Transport * eGly_13C00D11;
   Vf_Vin_Gly_13C01D11_13C01D11 = Kf_Gly_Transport * eGly_13C01D11;
   Vf_Vin_Gly_13C10D11_13C10D11 = Kf_Gly_Transport * eGly_13C10D11;
   Vf_Vin_Gly_13C11D11_13C11D11 = Kf_Gly_Transport * eGly_13C11D11;



    Vr_Vin_Gly_13C00D00_13C00D00 = Kr_Gly_Transport * cGlyPool2_13C00D00;
    Vr_Vin_Gly_13C01D00_13C01D00 = Kr_Gly_Transport * cGlyPool2_13C01D00;
    Vr_Vin_Gly_13C10D00_13C10D00 = Kr_Gly_Transport * cGlyPool2_13C10D00;
    Vr_Vin_Gly_13C11D00_13C11D00 = Kr_Gly_Transport * cGlyPool2_13C11D00;
    Vr_Vin_Gly_13C00D01_13C00D01 = Kr_Gly_Transport * cGlyPool2_13C00D01;
    Vr_Vin_Gly_13C01D01_13C01D01 = Kr_Gly_Transport * cGlyPool2_13C01D01;
    Vr_Vin_Gly_13C10D01_13C10D01 = Kr_Gly_Transport * cGlyPool2_13C10D01;
    Vr_Vin_Gly_13C11D01_13C11D01 = Kr_Gly_Transport * cGlyPool2_13C11D01;
    Vr_Vin_Gly_13C00D10_13C00D10 = Kr_Gly_Transport * cGlyPool2_13C00D10;
    Vr_Vin_Gly_13C01D10_13C01D10 = Kr_Gly_Transport * cGlyPool2_13C01D10;
    Vr_Vin_Gly_13C10D10_13C10D10 = Kr_Gly_Transport * cGlyPool2_13C10D10;
    Vr_Vin_Gly_13C11D10_13C11D10 = Kr_Gly_Transport * cGlyPool2_13C11D10;
    Vr_Vin_Gly_13C00D11_13C00D11 = Kr_Gly_Transport * cGlyPool2_13C00D11;
    Vr_Vin_Gly_13C01D11_13C01D11 = Kr_Gly_Transport * cGlyPool2_13C01D11;
    Vr_Vin_Gly_13C10D11_13C10D11 = Kr_Gly_Transport * cGlyPool2_13C10D11;
    Vr_Vin_Gly_13C11D11_13C11D11 = Kr_Gly_Transport * cGlyPool2_13C11D11;



    Vf_Vin_Gly3_13C00D00_13C00D00 = Kf_Gly_Transport3 * eGly_13C00D00;
    Vf_Vin_Gly3_13C01D00_13C01D00 = Kf_Gly_Transport3 * eGly_13C01D00;
    Vf_Vin_Gly3_13C10D00_13C10D00 = Kf_Gly_Transport3 * eGly_13C10D00;
    Vf_Vin_Gly3_13C11D00_13C11D00 = Kf_Gly_Transport3 * eGly_13C11D00;
    Vf_Vin_Gly3_13C00D01_13C00D01 = Kf_Gly_Transport3 * eGly_13C00D01;
    Vf_Vin_Gly3_13C01D01_13C01D01 = Kf_Gly_Transport3 * eGly_13C01D01;
    Vf_Vin_Gly3_13C10D01_13C10D01 = Kf_Gly_Transport3 * eGly_13C10D01;
    Vf_Vin_Gly3_13C11D01_13C11D01 = Kf_Gly_Transport3 * eGly_13C11D01;
    Vf_Vin_Gly3_13C00D10_13C00D10 = Kf_Gly_Transport3 * eGly_13C00D10;
    Vf_Vin_Gly3_13C01D10_13C01D10 = Kf_Gly_Transport3 * eGly_13C01D10;
    Vf_Vin_Gly3_13C10D10_13C10D10 = Kf_Gly_Transport3 * eGly_13C10D10;
    Vf_Vin_Gly3_13C11D10_13C11D10 = Kf_Gly_Transport3 * eGly_13C11D10;
    Vf_Vin_Gly3_13C00D11_13C00D11 = Kf_Gly_Transport3 * eGly_13C00D11;
    Vf_Vin_Gly3_13C01D11_13C01D11 = Kf_Gly_Transport3 * eGly_13C01D11;
    Vf_Vin_Gly3_13C10D11_13C10D11 = Kf_Gly_Transport3 * eGly_13C10D11;
    Vf_Vin_Gly3_13C11D11_13C11D11 = Kf_Gly_Transport3 * eGly_13C11D11;



    Vr_Vin_Gly3_13C00D00_13C00D00 = Kr_Gly_Transport3 * cGlyPool3_13C00D00;
    Vr_Vin_Gly3_13C01D00_13C01D00 = Kr_Gly_Transport3 * cGlyPool3_13C01D00;
    Vr_Vin_Gly3_13C10D00_13C10D00 = Kr_Gly_Transport3 * cGlyPool3_13C10D00;
    Vr_Vin_Gly3_13C11D00_13C11D00 = Kr_Gly_Transport3 * cGlyPool3_13C11D00;
    Vr_Vin_Gly3_13C00D01_13C00D01 = Kr_Gly_Transport3 * cGlyPool3_13C00D01;
    Vr_Vin_Gly3_13C01D01_13C01D01 = Kr_Gly_Transport3 * cGlyPool3_13C01D01;
    Vr_Vin_Gly3_13C10D01_13C10D01 = Kr_Gly_Transport3 * cGlyPool3_13C10D01;
    Vr_Vin_Gly3_13C11D01_13C11D01 = Kr_Gly_Transport3 * cGlyPool3_13C11D01;
    Vr_Vin_Gly3_13C00D10_13C00D10 = Kr_Gly_Transport3 * cGlyPool3_13C00D10;
    Vr_Vin_Gly3_13C01D10_13C01D10 = Kr_Gly_Transport3 * cGlyPool3_13C01D10;
    Vr_Vin_Gly3_13C10D10_13C10D10 = Kr_Gly_Transport3 * cGlyPool3_13C10D10;
    Vr_Vin_Gly3_13C11D10_13C11D10 = Kr_Gly_Transport3 * cGlyPool3_13C11D10;
    Vr_Vin_Gly3_13C00D11_13C00D11 = Kr_Gly_Transport3 * cGlyPool3_13C00D11;
    Vr_Vin_Gly3_13C01D11_13C01D11 = Kr_Gly_Transport3 * cGlyPool3_13C01D11;
    Vr_Vin_Gly3_13C10D11_13C10D11 = Kr_Gly_Transport3 * cGlyPool3_13C10D11;
    Vr_Vin_Gly3_13C11D11_13C11D11 = Kr_Gly_Transport3 * cGlyPool3_13C11D11;



		Vf_SerPool1_GlyPool1_13C000D000_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D000_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D000_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D000_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D000_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D000_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D000_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D000_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C000D001_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D001_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D001_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D001_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D001_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D001_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D001_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D001_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C000D010_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D010_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D010_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D010_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D010_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D010_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D010_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D010_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C000D011_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D011_13C00D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D011_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D011_13C01D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D011_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D011_13C10D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D011_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D011_13C11D00 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C000D100_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D100_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D100_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D100_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D100_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D100_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D100_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D100_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C000D101_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D101_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D101_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D101_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D101_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D101_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D101_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D101_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C000D110_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D110_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D110_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D110_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D110_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D110_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D110_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D110_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C000D111_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C000D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C001D111_13C00D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C001D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C010D111_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C010D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C011D111_13C01D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C011D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C100D111_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C100D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C101D111_13C10D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C101D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C110D111_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C110D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool1_GlyPool1_13C111D111_13C11D10 = Kf_SerPool1_GlyPool1 * cSerPool1_13C111D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



		Vr_SerPool1_GlyPool1_13C00D00_13C0D00__13C000D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C0D00__13C010D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C0D00__13C100D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C0D00__13C110D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C0D00__13C000D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C0D00__13C010D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C0D00__13C100D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C0D00__13C110D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C0D00__13C000D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C0D00__13C010D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C0D00__13C100D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C0D00__13C110D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C0D00__13C000D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C0D00__13C010D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C0D00__13C100D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C0D00__13C110D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D00_13C1D00__13C001D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C1D00__13C011D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C1D00__13C101D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C1D00__13C111D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C1D00__13C001D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C1D00__13C011D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C1D00__13C101D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C1D00__13C111D000 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C1D00__13C001D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C1D00__13C011D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C1D00__13C101D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C1D00__13C111D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C1D00__13C001D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C1D00__13C011D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C1D00__13C101D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C1D00__13C111D100 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D00_13C0D01__13C000D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C0D01__13C010D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C0D01__13C100D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C0D01__13C110D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C0D01__13C000D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C0D01__13C010D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C0D01__13C100D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C0D01__13C110D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C0D01__13C000D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C0D01__13C010D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C0D01__13C100D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C0D01__13C110D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C0D01__13C000D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C0D01__13C010D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C0D01__13C100D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C0D01__13C110D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D00_13C1D01__13C001D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C1D01__13C011D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C1D01__13C101D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C1D01__13C111D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C1D01__13C001D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C1D01__13C011D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C1D01__13C101D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C1D01__13C111D001 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C1D01__13C001D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C1D01__13C011D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C1D01__13C101D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C1D01__13C111D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C1D01__13C001D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C1D01__13C011D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C1D01__13C101D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C1D01__13C111D101 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D00_13C0D10__13C000D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C0D10__13C010D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C0D10__13C100D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C0D10__13C110D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C0D10__13C000D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C0D10__13C010D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C0D10__13C100D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C0D10__13C110D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C0D10__13C000D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C0D10__13C010D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C0D10__13C100D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C0D10__13C110D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C0D10__13C000D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C0D10__13C010D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C0D10__13C100D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C0D10__13C110D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D00_13C1D10__13C001D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C1D10__13C011D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C1D10__13C101D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C1D10__13C111D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C1D10__13C001D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C1D10__13C011D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C1D10__13C101D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C1D10__13C111D010 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C1D10__13C001D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C1D10__13C011D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C1D10__13C101D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C1D10__13C111D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C1D10__13C001D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C1D10__13C011D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C1D10__13C101D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C1D10__13C111D110 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D00_13C0D11__13C000D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C0D11__13C010D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C0D11__13C100D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C0D11__13C110D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C0D11__13C000D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C0D11__13C010D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C0D11__13C100D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C0D11__13C110D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C0D11__13C000D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C0D11__13C010D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C0D11__13C100D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C0D11__13C110D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C0D11__13C000D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C0D11__13C010D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C0D11__13C100D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C0D11__13C110D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D00_13C1D11__13C001D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D00_13C1D11__13C011D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D00_13C1D11__13C101D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D00_13C1D11__13C111D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D01_13C1D11__13C001D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D01_13C1D11__13C011D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D01_13C1D11__13C101D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D01_13C1D11__13C111D011 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D10_13C1D11__13C001D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D10_13C1D11__13C011D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D10_13C1D11__13C101D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D10_13C1D11__13C111D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C00D11_13C1D11__13C001D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C01D11_13C1D11__13C011D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C10D11_13C1D11__13C101D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool1_GlyPool1_13C11D11_13C1D11__13C111D111 = Kr_SerPool1_GlyPool1 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



		Vf_SerPool2_GlyPool2_13C000D000_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D000_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D000_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D000_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D000_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D000_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D000_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D000_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C000D001_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D001_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D001_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D001_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D001_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D001_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D001_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D001_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C000D010_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D010_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D010_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D010_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D010_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D010_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D010_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D010_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C000D011_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D011_13C00D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D011_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D011_13C01D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D011_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D011_13C10D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D011_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D011_13C11D00 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C000D100_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D100_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D100_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D100_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D100_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D100_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D100_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D100_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C000D101_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D101_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D101_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D101_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D101_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D101_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D101_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D101_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C000D110_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D110_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D110_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D110_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D110_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D110_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D110_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D110_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C000D111_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C000D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C001D111_13C00D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C001D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C010D111_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C010D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C011D111_13C01D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C011D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C100D111_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C100D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C101D111_13C10D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C101D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C110D111_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C110D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool2_GlyPool2_13C111D111_13C11D10 = Kf_SerPool2_GlyPool2 * cSerPool2_13C111D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



		Vr_SerPool2_GlyPool2_13C00D00_13C0D00__13C000D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C0D00__13C010D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C0D00__13C100D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C0D00__13C110D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C0D00__13C000D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C0D00__13C010D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C0D00__13C100D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C0D00__13C110D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C0D00__13C000D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C0D00__13C010D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C0D00__13C100D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C0D00__13C110D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C0D00__13C000D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C0D00__13C010D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C0D00__13C100D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C0D00__13C110D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D00_13C1D00__13C001D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C1D00__13C011D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C1D00__13C101D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C1D00__13C111D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C1D00__13C001D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C1D00__13C011D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C1D00__13C101D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C1D00__13C111D000 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C1D00__13C001D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C1D00__13C011D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C1D00__13C101D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C1D00__13C111D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C1D00__13C001D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C1D00__13C011D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C1D00__13C101D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C1D00__13C111D100 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D00_13C0D01__13C000D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C0D01__13C010D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C0D01__13C100D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C0D01__13C110D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C0D01__13C000D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C0D01__13C010D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C0D01__13C100D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C0D01__13C110D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C0D01__13C000D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C0D01__13C010D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C0D01__13C100D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C0D01__13C110D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C0D01__13C000D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C0D01__13C010D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C0D01__13C100D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C0D01__13C110D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D00_13C1D01__13C001D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C1D01__13C011D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C1D01__13C101D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C1D01__13C111D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C1D01__13C001D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C1D01__13C011D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C1D01__13C101D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C1D01__13C111D001 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C1D01__13C001D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C1D01__13C011D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C1D01__13C101D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C1D01__13C111D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C1D01__13C001D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C1D01__13C011D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C1D01__13C101D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C1D01__13C111D101 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D00_13C0D10__13C000D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C0D10__13C010D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C0D10__13C100D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C0D10__13C110D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C0D10__13C000D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C0D10__13C010D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C0D10__13C100D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C0D10__13C110D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C0D10__13C000D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C0D10__13C010D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C0D10__13C100D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C0D10__13C110D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C0D10__13C000D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C0D10__13C010D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C0D10__13C100D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C0D10__13C110D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D00_13C1D10__13C001D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C1D10__13C011D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C1D10__13C101D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C1D10__13C111D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C1D10__13C001D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C1D10__13C011D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C1D10__13C101D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C1D10__13C111D010 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C1D10__13C001D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C1D10__13C011D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C1D10__13C101D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C1D10__13C111D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C1D10__13C001D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C1D10__13C011D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C1D10__13C101D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C1D10__13C111D110 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D00_13C0D11__13C000D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C0D11__13C010D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C0D11__13C100D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C0D11__13C110D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C0D11__13C000D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C0D11__13C010D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C0D11__13C100D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C0D11__13C110D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C0D11__13C000D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C0D11__13C010D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C0D11__13C100D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C0D11__13C110D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C0D11__13C000D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C0D11__13C010D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C0D11__13C100D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C0D11__13C110D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D00_13C1D11__13C001D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D00_13C1D11__13C011D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D00_13C1D11__13C101D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D00_13C1D11__13C111D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D01_13C1D11__13C001D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D01_13C1D11__13C011D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D01_13C1D11__13C101D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D01_13C1D11__13C111D011 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D10_13C1D11__13C001D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D10_13C1D11__13C011D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D10_13C1D11__13C101D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D10_13C1D11__13C111D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C00D11_13C1D11__13C001D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C01D11_13C1D11__13C011D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C10D11_13C1D11__13C101D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool2_GlyPool2_13C11D11_13C1D11__13C111D111 = Kr_SerPool2_GlyPool2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D000_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D000_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D000_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D000_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D000_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D000_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D000_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D000_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D000 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D001_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D001_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D001_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D001_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D001_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D001_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D001_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D001_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D001 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D010_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D010_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D010_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D010_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D010_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D010_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D010_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D010_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D010 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D011_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D011_13C00D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D011_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D011_13C01D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D011_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D011_13C10D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D011_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D011_13C11D00 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D011 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D100_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D100_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D100_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D100_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D100_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D100_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D100_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D100_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D100 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D101_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D101_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D101_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D101_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D101_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D101_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D101_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D101_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D101 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D110_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D110_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D110_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D110_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D110_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D110_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D110_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D110_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D110 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C000D111_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C000D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C001D111_13C00D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C001D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C010D111_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C010D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C011D111_13C01D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C011D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C100D111_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C100D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C101D111_13C10D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C101D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C110D111_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C110D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPool3_GlyPool3_13C111D111_13C11D10 = Kf_SerPool3_GlyPool3 * cSerPool3_13C111D111 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



		Vr_SerPool3_GlyPool3_13C00D00_13C0D00__13C000D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C0D00__13C010D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C0D00__13C100D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C0D00__13C110D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C0D00__13C000D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C0D00__13C010D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C0D00__13C100D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C0D00__13C110D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C0D00__13C000D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C0D00__13C010D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C0D00__13C100D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C0D00__13C110D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C0D00__13C000D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C0D00__13C010D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C0D00__13C100D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C0D00__13C110D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D00_13C1D00__13C001D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C1D00__13C011D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C1D00__13C101D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C1D00__13C111D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C1D00__13C001D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C1D00__13C011D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C1D00__13C101D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C1D00__13C111D000 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C1D00__13C001D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C1D00__13C011D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C1D00__13C101D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C1D00__13C111D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C1D00__13C001D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C1D00__13C011D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C1D00__13C101D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C1D00__13C111D100 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D00_13C0D01__13C000D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C0D01__13C010D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C0D01__13C100D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C0D01__13C110D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C0D01__13C000D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C0D01__13C010D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C0D01__13C100D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C0D01__13C110D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C0D01__13C000D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C0D01__13C010D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C0D01__13C100D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C0D01__13C110D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C0D01__13C000D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C0D01__13C010D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C0D01__13C100D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C0D01__13C110D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D00_13C1D01__13C001D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C1D01__13C011D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C1D01__13C101D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C1D01__13C111D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C1D01__13C001D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C1D01__13C011D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C1D01__13C101D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C1D01__13C111D001 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C1D01__13C001D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C1D01__13C011D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C1D01__13C101D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C1D01__13C111D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C1D01__13C001D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C1D01__13C011D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C1D01__13C101D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C1D01__13C111D101 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D00_13C0D10__13C000D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C0D10__13C010D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C0D10__13C100D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C0D10__13C110D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C0D10__13C000D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C0D10__13C010D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C0D10__13C100D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C0D10__13C110D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C0D10__13C000D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C0D10__13C010D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C0D10__13C100D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C0D10__13C110D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C0D10__13C000D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C0D10__13C010D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C0D10__13C100D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C0D10__13C110D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D00_13C1D10__13C001D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C1D10__13C011D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C1D10__13C101D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C1D10__13C111D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C1D10__13C001D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C1D10__13C011D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C1D10__13C101D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C1D10__13C111D010 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C1D10__13C001D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C1D10__13C011D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C1D10__13C101D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C1D10__13C111D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C1D10__13C001D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C1D10__13C011D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C1D10__13C101D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C1D10__13C111D110 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D00_13C0D11__13C000D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C0D11__13C010D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C0D11__13C100D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C0D11__13C110D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C0D11__13C000D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C0D11__13C010D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C0D11__13C100D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C0D11__13C110D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C0D11__13C000D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C0D11__13C010D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C0D11__13C100D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C0D11__13C110D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C0D11__13C000D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C0D11__13C010D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C0D11__13C100D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C0D11__13C110D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D00_13C1D11__13C001D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D00_13C1D11__13C011D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D00_13C1D11__13C101D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D00_13C1D11__13C111D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D01_13C1D11__13C001D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D01_13C1D11__13C011D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D01_13C1D11__13C101D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D01_13C1D11__13C111D011 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D10_13C1D11__13C001D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D10_13C1D11__13C011D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D10_13C1D11__13C101D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D10_13C1D11__13C111D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C00D11_13C1D11__13C001D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C01D11_13C1D11__13C011D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C10D11_13C1D11__13C101D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPool3_GlyPool3_13C11D11_13C1D11__13C111D111 = Kr_SerPool3_GlyPool3 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);




        Vf_SerPoolMitochon_GlyPoolMitochon_13C000D000_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D000_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D000_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D000_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D000_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D000_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D000_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D000_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D000 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C000D001_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D001_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D001_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D001_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D001_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D001_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D001_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D001_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D001 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C000D010_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D010_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D010_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D010_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D010_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D010_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D010_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D010_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D010 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C000D011_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D011_13C00D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D011_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D011_13C01D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D011_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D011_13C10D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D011_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D011_13C11D00 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D011 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C000D100_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D100_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D100_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D100_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D100_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D100_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D100_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D100_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D100 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C000D101_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D101_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D101_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D101_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D101_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D101_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D101_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D101_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D101 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C000D110_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D110_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D110_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D110_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D110_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D110_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D110_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D110_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D110 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C000D111_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C000D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C001D111_13C00D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C001D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C010D111_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C010D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C011D111_13C01D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C011D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C100D111_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C100D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C101D111_13C10D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C101D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C110D111_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C110D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vf_SerPoolMitochon_GlyPoolMitochon_13C111D111_13C11D10 = Kf_SerPoolMitochon_GlyPoolMitochon * mSerPool_13C111D111 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D00__13C000D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D00__13C010D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D00__13C100D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D00__13C110D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D00__13C000D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D00__13C010D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D00__13C100D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D00__13C110D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D00__13C000D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D00__13C010D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D00__13C100D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D00__13C110D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D00__13C000D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D00__13C010D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D00__13C100D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D00__13C110D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D00__13C001D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D00__13C011D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D00__13C101D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D00__13C111D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D00__13C001D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D00__13C011D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D00__13C101D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D00__13C111D000 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D00__13C001D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D00__13C011D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D00__13C101D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D00__13C111D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D00__13C001D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D00__13C011D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D00__13C101D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D00__13C111D100 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D01__13C000D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D01__13C010D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D01__13C100D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D01__13C110D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D01__13C000D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D01__13C010D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D01__13C100D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D01__13C110D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D01__13C000D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D01__13C010D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D01__13C100D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D01__13C110D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D01__13C000D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D01__13C010D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D01__13C100D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D01__13C110D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D01__13C001D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D01__13C011D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D01__13C101D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D01__13C111D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D01__13C001D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D01__13C011D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D01__13C101D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D01__13C111D001 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D01__13C001D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D01__13C011D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D01__13C101D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D01__13C111D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D01__13C001D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D01__13C011D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D01__13C101D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D01__13C111D101 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D10__13C000D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D10__13C010D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D10__13C100D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D10__13C110D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D10__13C000D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D10__13C010D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D10__13C100D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D10__13C110D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D10__13C000D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D10__13C010D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D10__13C100D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D10__13C110D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D10__13C000D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D10__13C010D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D10__13C100D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D10__13C110D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D10__13C001D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D10__13C011D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D10__13C101D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D10__13C111D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D10__13C001D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D10__13C011D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D10__13C101D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D10__13C111D010 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D10__13C001D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D10__13C011D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D10__13C101D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D10__13C111D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D10__13C001D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D10__13C011D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D10__13C101D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D10__13C111D110 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D11__13C000D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D11__13C010D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D11__13C100D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D11__13C110D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D11__13C000D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D11__13C010D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D11__13C100D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D11__13C110D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D11__13C000D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D11__13C010D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D11__13C100D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D11__13C110D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D11__13C000D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D11__13C010D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D11__13C100D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D11__13C110D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C0D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D11__13C001D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D11__13C011D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D11__13C101D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D11__13C111D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D11__13C001D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D11__13C011D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D11__13C101D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D11__13C111D011 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D11__13C001D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D11__13C011D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D11__13C101D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D11__13C111D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D11__13C001D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D11__13C011D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D11__13C101D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
		Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D11__13C111D111 = Kr_SerPoolMitochon_GlyPoolMitochon * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mMethyleneTHF_13C1D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vf_GlyPool1_CO2_13C00D00_13C0D00 = Kf_GlyPool1_CO2 * cGlyPool1_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C01D00_13C1D00 = Kf_GlyPool1_CO2 * cGlyPool1_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C10D00_13C0D00 = Kf_GlyPool1_CO2 * cGlyPool1_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C11D00_13C1D00 = Kf_GlyPool1_CO2 * cGlyPool1_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C00D01_13C0D01 = Kf_GlyPool1_CO2 * cGlyPool1_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C01D01_13C1D01 = Kf_GlyPool1_CO2 * cGlyPool1_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C10D01_13C0D01 = Kf_GlyPool1_CO2 * cGlyPool1_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C11D01_13C1D01 = Kf_GlyPool1_CO2 * cGlyPool1_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C00D10_13C0D10 = Kf_GlyPool1_CO2 * cGlyPool1_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C01D10_13C1D10 = Kf_GlyPool1_CO2 * cGlyPool1_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C10D10_13C0D10 = Kf_GlyPool1_CO2 * cGlyPool1_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C11D10_13C1D10 = Kf_GlyPool1_CO2 * cGlyPool1_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C00D11_13C0D11 = Kf_GlyPool1_CO2 * cGlyPool1_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C01D11_13C1D11 = Kf_GlyPool1_CO2 * cGlyPool1_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C10D11_13C0D11 = Kf_GlyPool1_CO2 * cGlyPool1_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool1_CO2_13C11D11_13C1D11 = Kf_GlyPool1_CO2 * cGlyPool1_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vr_GlyPool1_CO2_13C0D00_13C00D00 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D00;
    Vr_GlyPool1_CO2_13C1D00_13C01D00 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D00;
    Vr_GlyPool1_CO2_13C0D01_13C00D01 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D01;
    Vr_GlyPool1_CO2_13C1D01_13C01D01 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D01;
    Vr_GlyPool1_CO2_13C0D10_13C00D10 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D10;
    Vr_GlyPool1_CO2_13C1D10_13C01D10 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D10;
    Vr_GlyPool1_CO2_13C0D11_13C00D11 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D11;
    Vr_GlyPool1_CO2_13C1D11_13C01D11 = Kr_GlyPool1_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D11;



    Vf_GlyPool2_CO2_13C00D00_13C0D00 = Kf_GlyPool2_CO2 * cGlyPool2_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C01D00_13C1D00 = Kf_GlyPool2_CO2 * cGlyPool2_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C10D00_13C0D00 = Kf_GlyPool2_CO2 * cGlyPool2_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C11D00_13C1D00 = Kf_GlyPool2_CO2 * cGlyPool2_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C00D01_13C0D01 = Kf_GlyPool2_CO2 * cGlyPool2_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C01D01_13C1D01 = Kf_GlyPool2_CO2 * cGlyPool2_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C10D01_13C0D01 = Kf_GlyPool2_CO2 * cGlyPool2_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C11D01_13C1D01 = Kf_GlyPool2_CO2 * cGlyPool2_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C00D10_13C0D10 = Kf_GlyPool2_CO2 * cGlyPool2_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C01D10_13C1D10 = Kf_GlyPool2_CO2 * cGlyPool2_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C10D10_13C0D10 = Kf_GlyPool2_CO2 * cGlyPool2_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C11D10_13C1D10 = Kf_GlyPool2_CO2 * cGlyPool2_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C00D11_13C0D11 = Kf_GlyPool2_CO2 * cGlyPool2_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C01D11_13C1D11 = Kf_GlyPool2_CO2 * cGlyPool2_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C10D11_13C0D11 = Kf_GlyPool2_CO2 * cGlyPool2_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool2_CO2_13C11D11_13C1D11 = Kf_GlyPool2_CO2 * cGlyPool2_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vr_GlyPool2_CO2_13C0D00_13C00D00 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D00;
    Vr_GlyPool2_CO2_13C1D00_13C01D00 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D00;
    Vr_GlyPool2_CO2_13C0D01_13C00D01 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D01;
    Vr_GlyPool2_CO2_13C1D01_13C01D01 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D01;
    Vr_GlyPool2_CO2_13C0D10_13C00D10 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D10;
    Vr_GlyPool2_CO2_13C1D10_13C01D10 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D10;
    Vr_GlyPool2_CO2_13C0D11_13C00D11 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D11;
    Vr_GlyPool2_CO2_13C1D11_13C01D11 = Kr_GlyPool2_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D11;



    Vf_GlyPool3_CO2_13C00D00_13C0D00 = Kf_GlyPool3_CO2 * cGlyPool3_13C00D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C01D00_13C1D00 = Kf_GlyPool3_CO2 * cGlyPool3_13C01D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C10D00_13C0D00 = Kf_GlyPool3_CO2 * cGlyPool3_13C10D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C11D00_13C1D00 = Kf_GlyPool3_CO2 * cGlyPool3_13C11D00 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C00D01_13C0D01 = Kf_GlyPool3_CO2 * cGlyPool3_13C00D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C01D01_13C1D01 = Kf_GlyPool3_CO2 * cGlyPool3_13C01D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C10D01_13C0D01 = Kf_GlyPool3_CO2 * cGlyPool3_13C10D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C11D01_13C1D01 = Kf_GlyPool3_CO2 * cGlyPool3_13C11D01 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C00D10_13C0D10 = Kf_GlyPool3_CO2 * cGlyPool3_13C00D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C01D10_13C1D10 = Kf_GlyPool3_CO2 * cGlyPool3_13C01D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C10D10_13C0D10 = Kf_GlyPool3_CO2 * cGlyPool3_13C10D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C11D10_13C1D10 = Kf_GlyPool3_CO2 * cGlyPool3_13C11D10 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C00D11_13C0D11 = Kf_GlyPool3_CO2 * cGlyPool3_13C00D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C01D11_13C1D11 = Kf_GlyPool3_CO2 * cGlyPool3_13C01D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C10D11_13C0D11 = Kf_GlyPool3_CO2 * cGlyPool3_13C10D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPool3_CO2_13C11D11_13C1D11 = Kf_GlyPool3_CO2 * cGlyPool3_13C11D11 / (6.02214129 * 10^11) / CytoplasmVolume * cTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vr_GlyPool3_CO2_13C0D00_13C00D00 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D00;
    Vr_GlyPool3_CO2_13C1D00_13C01D00 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D00;
    Vr_GlyPool3_CO2_13C0D01_13C00D01 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D01;
    Vr_GlyPool3_CO2_13C1D01_13C01D01 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D01;
    Vr_GlyPool3_CO2_13C0D10_13C00D10 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D10;
    Vr_GlyPool3_CO2_13C1D10_13C01D10 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D10;
    Vr_GlyPool3_CO2_13C0D11_13C00D11 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C0D11;
    Vr_GlyPool3_CO2_13C1D11_13C01D11 = Kr_GlyPool3_CO2 * cCO2 * cNH3 * cMethyleneTHF_13C1D11;



    Vf_GlyPoolMitochon_CO2_13C00D00_13C0D00 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C00D00 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C01D00_13C1D00 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C01D00 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C10D00_13C0D00 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C10D00 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C11D00_13C1D00 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C11D00 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C00D01_13C0D01 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C00D01 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C01D01_13C1D01 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C01D01 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C10D01_13C0D01 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C10D01 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C11D01_13C1D01 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C11D01 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C00D10_13C0D10 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C00D10 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C01D10_13C1D10 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C01D10 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C10D10_13C0D10 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C10D10 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C11D10_13C1D10 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C11D10 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C00D11_13C0D11 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C00D11 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C01D11_13C1D11 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C01D11 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C10D11_13C0D11 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C10D11 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GlyPoolMitochon_CO2_13C11D11_13C1D11 = Kf_GlyPoolMitochon_CO2 * mGlyPool_13C11D11 / (6.02214129 * 10^11) / MitochonVolume * mTHF / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vr_GlyPoolMitochon_CO2_13C0D00_13C00D00 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C0D00;
    Vr_GlyPoolMitochon_CO2_13C1D00_13C01D00 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C1D00;
    Vr_GlyPoolMitochon_CO2_13C0D01_13C00D01 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C0D01;
    Vr_GlyPoolMitochon_CO2_13C1D01_13C01D01 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C1D01;
    Vr_GlyPoolMitochon_CO2_13C0D10_13C00D10 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C0D10;
    Vr_GlyPoolMitochon_CO2_13C1D10_13C01D10 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C1D10;
    Vr_GlyPoolMitochon_CO2_13C0D11_13C00D11 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C0D11;
    Vr_GlyPoolMitochon_CO2_13C1D11_13C01D11 = Kr_GlyPoolMitochon_CO2 * mCO2 * mNH3 * mMethyleneTHF_13C1D11;



		Vf_SerPool3_Synthesis_13C000D000_13C000D000 = Kf_SerPool3_Synthesis * (6.02214129 * 10^17) * CytoplasmVolume;



    Vf_SerPool3_Degradation_13C000D000_13C000D000 = Kf_SerPool3_Degradation * cSerPool3_13C000D000;
    Vf_SerPool3_Degradation_13C001D000_13C001D000 = Kf_SerPool3_Degradation * cSerPool3_13C001D000;
    Vf_SerPool3_Degradation_13C010D000_13C010D000 = Kf_SerPool3_Degradation * cSerPool3_13C010D000;
    Vf_SerPool3_Degradation_13C011D000_13C011D000 = Kf_SerPool3_Degradation * cSerPool3_13C011D000;
    Vf_SerPool3_Degradation_13C100D000_13C100D000 = Kf_SerPool3_Degradation * cSerPool3_13C100D000;
    Vf_SerPool3_Degradation_13C101D000_13C101D000 = Kf_SerPool3_Degradation * cSerPool3_13C101D000;
    Vf_SerPool3_Degradation_13C110D000_13C110D000 = Kf_SerPool3_Degradation * cSerPool3_13C110D000;
    Vf_SerPool3_Degradation_13C111D000_13C111D000 = Kf_SerPool3_Degradation * cSerPool3_13C111D000;
    Vf_SerPool3_Degradation_13C000D001_13C000D001 = Kf_SerPool3_Degradation * cSerPool3_13C000D001;
    Vf_SerPool3_Degradation_13C001D001_13C001D001 = Kf_SerPool3_Degradation * cSerPool3_13C001D001;
    Vf_SerPool3_Degradation_13C010D001_13C010D001 = Kf_SerPool3_Degradation * cSerPool3_13C010D001;
    Vf_SerPool3_Degradation_13C011D001_13C011D001 = Kf_SerPool3_Degradation * cSerPool3_13C011D001;
    Vf_SerPool3_Degradation_13C100D001_13C100D001 = Kf_SerPool3_Degradation * cSerPool3_13C100D001;
    Vf_SerPool3_Degradation_13C101D001_13C101D001 = Kf_SerPool3_Degradation * cSerPool3_13C101D001;
    Vf_SerPool3_Degradation_13C110D001_13C110D001 = Kf_SerPool3_Degradation * cSerPool3_13C110D001;
    Vf_SerPool3_Degradation_13C111D001_13C111D001 = Kf_SerPool3_Degradation * cSerPool3_13C111D001;
    Vf_SerPool3_Degradation_13C000D010_13C000D010 = Kf_SerPool3_Degradation * cSerPool3_13C000D010;
    Vf_SerPool3_Degradation_13C001D010_13C001D010 = Kf_SerPool3_Degradation * cSerPool3_13C001D010;
    Vf_SerPool3_Degradation_13C010D010_13C010D010 = Kf_SerPool3_Degradation * cSerPool3_13C010D010;
    Vf_SerPool3_Degradation_13C011D010_13C011D010 = Kf_SerPool3_Degradation * cSerPool3_13C011D010;
    Vf_SerPool3_Degradation_13C100D010_13C100D010 = Kf_SerPool3_Degradation * cSerPool3_13C100D010;
    Vf_SerPool3_Degradation_13C101D010_13C101D010 = Kf_SerPool3_Degradation * cSerPool3_13C101D010;
    Vf_SerPool3_Degradation_13C110D010_13C110D010 = Kf_SerPool3_Degradation * cSerPool3_13C110D010;
    Vf_SerPool3_Degradation_13C111D010_13C111D010 = Kf_SerPool3_Degradation * cSerPool3_13C111D010;
    Vf_SerPool3_Degradation_13C000D011_13C000D011 = Kf_SerPool3_Degradation * cSerPool3_13C000D011;
    Vf_SerPool3_Degradation_13C001D011_13C001D011 = Kf_SerPool3_Degradation * cSerPool3_13C001D011;
    Vf_SerPool3_Degradation_13C010D011_13C010D011 = Kf_SerPool3_Degradation * cSerPool3_13C010D011;
    Vf_SerPool3_Degradation_13C011D011_13C011D011 = Kf_SerPool3_Degradation * cSerPool3_13C011D011;
    Vf_SerPool3_Degradation_13C100D011_13C100D011 = Kf_SerPool3_Degradation * cSerPool3_13C100D011;
    Vf_SerPool3_Degradation_13C101D011_13C101D011 = Kf_SerPool3_Degradation * cSerPool3_13C101D011;
    Vf_SerPool3_Degradation_13C110D011_13C110D011 = Kf_SerPool3_Degradation * cSerPool3_13C110D011;
    Vf_SerPool3_Degradation_13C111D011_13C111D011 = Kf_SerPool3_Degradation * cSerPool3_13C111D011;
    Vf_SerPool3_Degradation_13C000D100_13C000D100 = Kf_SerPool3_Degradation * cSerPool3_13C000D100;
    Vf_SerPool3_Degradation_13C001D100_13C001D100 = Kf_SerPool3_Degradation * cSerPool3_13C001D100;
    Vf_SerPool3_Degradation_13C010D100_13C010D100 = Kf_SerPool3_Degradation * cSerPool3_13C010D100;
    Vf_SerPool3_Degradation_13C011D100_13C011D100 = Kf_SerPool3_Degradation * cSerPool3_13C011D100;
    Vf_SerPool3_Degradation_13C100D100_13C100D100 = Kf_SerPool3_Degradation * cSerPool3_13C100D100;
    Vf_SerPool3_Degradation_13C101D100_13C101D100 = Kf_SerPool3_Degradation * cSerPool3_13C101D100;
    Vf_SerPool3_Degradation_13C110D100_13C110D100 = Kf_SerPool3_Degradation * cSerPool3_13C110D100;
    Vf_SerPool3_Degradation_13C111D100_13C111D100 = Kf_SerPool3_Degradation * cSerPool3_13C111D100;
    Vf_SerPool3_Degradation_13C000D101_13C000D101 = Kf_SerPool3_Degradation * cSerPool3_13C000D101;
    Vf_SerPool3_Degradation_13C001D101_13C001D101 = Kf_SerPool3_Degradation * cSerPool3_13C001D101;
    Vf_SerPool3_Degradation_13C010D101_13C010D101 = Kf_SerPool3_Degradation * cSerPool3_13C010D101;
    Vf_SerPool3_Degradation_13C011D101_13C011D101 = Kf_SerPool3_Degradation * cSerPool3_13C011D101;
    Vf_SerPool3_Degradation_13C100D101_13C100D101 = Kf_SerPool3_Degradation * cSerPool3_13C100D101;
    Vf_SerPool3_Degradation_13C101D101_13C101D101 = Kf_SerPool3_Degradation * cSerPool3_13C101D101;
    Vf_SerPool3_Degradation_13C110D101_13C110D101 = Kf_SerPool3_Degradation * cSerPool3_13C110D101;
    Vf_SerPool3_Degradation_13C111D101_13C111D101 = Kf_SerPool3_Degradation * cSerPool3_13C111D101;
    Vf_SerPool3_Degradation_13C000D110_13C000D110 = Kf_SerPool3_Degradation * cSerPool3_13C000D110;
    Vf_SerPool3_Degradation_13C001D110_13C001D110 = Kf_SerPool3_Degradation * cSerPool3_13C001D110;
    Vf_SerPool3_Degradation_13C010D110_13C010D110 = Kf_SerPool3_Degradation * cSerPool3_13C010D110;
    Vf_SerPool3_Degradation_13C011D110_13C011D110 = Kf_SerPool3_Degradation * cSerPool3_13C011D110;
    Vf_SerPool3_Degradation_13C100D110_13C100D110 = Kf_SerPool3_Degradation * cSerPool3_13C100D110;
    Vf_SerPool3_Degradation_13C101D110_13C101D110 = Kf_SerPool3_Degradation * cSerPool3_13C101D110;
    Vf_SerPool3_Degradation_13C110D110_13C110D110 = Kf_SerPool3_Degradation * cSerPool3_13C110D110;
    Vf_SerPool3_Degradation_13C111D110_13C111D110 = Kf_SerPool3_Degradation * cSerPool3_13C111D110;
    Vf_SerPool3_Degradation_13C000D111_13C000D111 = Kf_SerPool3_Degradation * cSerPool3_13C000D111;
    Vf_SerPool3_Degradation_13C001D111_13C001D111 = Kf_SerPool3_Degradation * cSerPool3_13C001D111;
    Vf_SerPool3_Degradation_13C010D111_13C010D111 = Kf_SerPool3_Degradation * cSerPool3_13C010D111;
    Vf_SerPool3_Degradation_13C011D111_13C011D111 = Kf_SerPool3_Degradation * cSerPool3_13C011D111;
    Vf_SerPool3_Degradation_13C100D111_13C100D111 = Kf_SerPool3_Degradation * cSerPool3_13C100D111;
    Vf_SerPool3_Degradation_13C101D111_13C101D111 = Kf_SerPool3_Degradation * cSerPool3_13C101D111;
    Vf_SerPool3_Degradation_13C110D111_13C110D111 = Kf_SerPool3_Degradation * cSerPool3_13C110D111;
    Vf_SerPool3_Degradation_13C111D111_13C111D111 = Kf_SerPool3_Degradation * cSerPool3_13C111D111;



    Vf_GlyPool3_Synthesis_13C00D00_13C00D00 = Kf_GlyPool3_Synthesis * (6.02214129 * 10^17) * CytoplasmVolume;



    Vf_GlyPool3_Degradation_13C00D00_13C00D00 = Kf_GlyPool3_Degradation * cGlyPool3_13C00D00;
    Vf_GlyPool3_Degradation_13C01D00_13C01D00 = Kf_GlyPool3_Degradation * cGlyPool3_13C01D00;
    Vf_GlyPool3_Degradation_13C10D00_13C10D00 = Kf_GlyPool3_Degradation * cGlyPool3_13C10D00;
    Vf_GlyPool3_Degradation_13C11D00_13C11D00 = Kf_GlyPool3_Degradation * cGlyPool3_13C11D00;
    Vf_GlyPool3_Degradation_13C00D01_13C00D01 = Kf_GlyPool3_Degradation * cGlyPool3_13C00D01;
    Vf_GlyPool3_Degradation_13C01D01_13C01D01 = Kf_GlyPool3_Degradation * cGlyPool3_13C01D01;
    Vf_GlyPool3_Degradation_13C10D01_13C10D01 = Kf_GlyPool3_Degradation * cGlyPool3_13C10D01;
    Vf_GlyPool3_Degradation_13C11D01_13C11D01 = Kf_GlyPool3_Degradation * cGlyPool3_13C11D01;
    Vf_GlyPool3_Degradation_13C00D10_13C00D10 = Kf_GlyPool3_Degradation * cGlyPool3_13C00D10;
    Vf_GlyPool3_Degradation_13C01D10_13C01D10 = Kf_GlyPool3_Degradation * cGlyPool3_13C01D10;
    Vf_GlyPool3_Degradation_13C10D10_13C10D10 = Kf_GlyPool3_Degradation * cGlyPool3_13C10D10;
    Vf_GlyPool3_Degradation_13C11D10_13C11D10 = Kf_GlyPool3_Degradation * cGlyPool3_13C11D10;
    Vf_GlyPool3_Degradation_13C00D11_13C00D11 = Kf_GlyPool3_Degradation * cGlyPool3_13C00D11;
    Vf_GlyPool3_Degradation_13C01D11_13C01D11 = Kf_GlyPool3_Degradation * cGlyPool3_13C01D11;
    Vf_GlyPool3_Degradation_13C10D11_13C10D11 = Kf_GlyPool3_Degradation * cGlyPool3_13C10D11;
    Vf_GlyPool3_Degradation_13C11D11_13C11D11 = Kf_GlyPool3_Degradation * cGlyPool3_13C11D11;



    Vf_SerPool2_SerPoolMitochon_13C000D000_13C000D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D000;
    Vf_SerPool2_SerPoolMitochon_13C001D000_13C001D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D000;
    Vf_SerPool2_SerPoolMitochon_13C010D000_13C010D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D000;
    Vf_SerPool2_SerPoolMitochon_13C011D000_13C011D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D000;
    Vf_SerPool2_SerPoolMitochon_13C100D000_13C100D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D000;
    Vf_SerPool2_SerPoolMitochon_13C101D000_13C101D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D000;
    Vf_SerPool2_SerPoolMitochon_13C110D000_13C110D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D000;
    Vf_SerPool2_SerPoolMitochon_13C111D000_13C111D000 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D000;
    Vf_SerPool2_SerPoolMitochon_13C000D001_13C000D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D001;
    Vf_SerPool2_SerPoolMitochon_13C001D001_13C001D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D001;
    Vf_SerPool2_SerPoolMitochon_13C010D001_13C010D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D001;
    Vf_SerPool2_SerPoolMitochon_13C011D001_13C011D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D001;
    Vf_SerPool2_SerPoolMitochon_13C100D001_13C100D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D001;
    Vf_SerPool2_SerPoolMitochon_13C101D001_13C101D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D001;
    Vf_SerPool2_SerPoolMitochon_13C110D001_13C110D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D001;
    Vf_SerPool2_SerPoolMitochon_13C111D001_13C111D001 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D001;
    Vf_SerPool2_SerPoolMitochon_13C000D010_13C000D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D010;
    Vf_SerPool2_SerPoolMitochon_13C001D010_13C001D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D010;
    Vf_SerPool2_SerPoolMitochon_13C010D010_13C010D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D010;
    Vf_SerPool2_SerPoolMitochon_13C011D010_13C011D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D010;
    Vf_SerPool2_SerPoolMitochon_13C100D010_13C100D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D010;
    Vf_SerPool2_SerPoolMitochon_13C101D010_13C101D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D010;
    Vf_SerPool2_SerPoolMitochon_13C110D010_13C110D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D010;
    Vf_SerPool2_SerPoolMitochon_13C111D010_13C111D010 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D010;
    Vf_SerPool2_SerPoolMitochon_13C000D011_13C000D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D011;
    Vf_SerPool2_SerPoolMitochon_13C001D011_13C001D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D011;
    Vf_SerPool2_SerPoolMitochon_13C010D011_13C010D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D011;
    Vf_SerPool2_SerPoolMitochon_13C011D011_13C011D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D011;
    Vf_SerPool2_SerPoolMitochon_13C100D011_13C100D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D011;
    Vf_SerPool2_SerPoolMitochon_13C101D011_13C101D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D011;
    Vf_SerPool2_SerPoolMitochon_13C110D011_13C110D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D011;
    Vf_SerPool2_SerPoolMitochon_13C111D011_13C111D011 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D011;
    Vf_SerPool2_SerPoolMitochon_13C000D100_13C000D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D100;
    Vf_SerPool2_SerPoolMitochon_13C001D100_13C001D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D100;
    Vf_SerPool2_SerPoolMitochon_13C010D100_13C010D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D100;
    Vf_SerPool2_SerPoolMitochon_13C011D100_13C011D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D100;
    Vf_SerPool2_SerPoolMitochon_13C100D100_13C100D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D100;
    Vf_SerPool2_SerPoolMitochon_13C101D100_13C101D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D100;
    Vf_SerPool2_SerPoolMitochon_13C110D100_13C110D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D100;
    Vf_SerPool2_SerPoolMitochon_13C111D100_13C111D100 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D100;
    Vf_SerPool2_SerPoolMitochon_13C000D101_13C000D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D101;
    Vf_SerPool2_SerPoolMitochon_13C001D101_13C001D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D101;
    Vf_SerPool2_SerPoolMitochon_13C010D101_13C010D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D101;
    Vf_SerPool2_SerPoolMitochon_13C011D101_13C011D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D101;
    Vf_SerPool2_SerPoolMitochon_13C100D101_13C100D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D101;
    Vf_SerPool2_SerPoolMitochon_13C101D101_13C101D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D101;
    Vf_SerPool2_SerPoolMitochon_13C110D101_13C110D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D101;
    Vf_SerPool2_SerPoolMitochon_13C111D101_13C111D101 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D101;
    Vf_SerPool2_SerPoolMitochon_13C000D110_13C000D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D110;
    Vf_SerPool2_SerPoolMitochon_13C001D110_13C001D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D110;
    Vf_SerPool2_SerPoolMitochon_13C010D110_13C010D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D110;
    Vf_SerPool2_SerPoolMitochon_13C011D110_13C011D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D110;
    Vf_SerPool2_SerPoolMitochon_13C100D110_13C100D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D110;
    Vf_SerPool2_SerPoolMitochon_13C101D110_13C101D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D110;
    Vf_SerPool2_SerPoolMitochon_13C110D110_13C110D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D110;
    Vf_SerPool2_SerPoolMitochon_13C111D110_13C111D110 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D110;
    Vf_SerPool2_SerPoolMitochon_13C000D111_13C000D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C000D111;
    Vf_SerPool2_SerPoolMitochon_13C001D111_13C001D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C001D111;
    Vf_SerPool2_SerPoolMitochon_13C010D111_13C010D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C010D111;
    Vf_SerPool2_SerPoolMitochon_13C011D111_13C011D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C011D111;
    Vf_SerPool2_SerPoolMitochon_13C100D111_13C100D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C100D111;
    Vf_SerPool2_SerPoolMitochon_13C101D111_13C101D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C101D111;
    Vf_SerPool2_SerPoolMitochon_13C110D111_13C110D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C110D111;
    Vf_SerPool2_SerPoolMitochon_13C111D111_13C111D111 = Kf_SerPool2_SerPoolMitochon * cSerPool2_13C111D111;



    Vr_SerPool2_SerPoolMitochon_13C000D000_13C000D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D000;
    Vr_SerPool2_SerPoolMitochon_13C001D000_13C001D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D000;
    Vr_SerPool2_SerPoolMitochon_13C010D000_13C010D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D000;
    Vr_SerPool2_SerPoolMitochon_13C011D000_13C011D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D000;
    Vr_SerPool2_SerPoolMitochon_13C100D000_13C100D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D000;
    Vr_SerPool2_SerPoolMitochon_13C101D000_13C101D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D000;
    Vr_SerPool2_SerPoolMitochon_13C110D000_13C110D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D000;
    Vr_SerPool2_SerPoolMitochon_13C111D000_13C111D000 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D000;
    Vr_SerPool2_SerPoolMitochon_13C000D001_13C000D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D001;
    Vr_SerPool2_SerPoolMitochon_13C001D001_13C001D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D001;
    Vr_SerPool2_SerPoolMitochon_13C010D001_13C010D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D001;
    Vr_SerPool2_SerPoolMitochon_13C011D001_13C011D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D001;
    Vr_SerPool2_SerPoolMitochon_13C100D001_13C100D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D001;
    Vr_SerPool2_SerPoolMitochon_13C101D001_13C101D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D001;
    Vr_SerPool2_SerPoolMitochon_13C110D001_13C110D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D001;
    Vr_SerPool2_SerPoolMitochon_13C111D001_13C111D001 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D001;
    Vr_SerPool2_SerPoolMitochon_13C000D010_13C000D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D010;
    Vr_SerPool2_SerPoolMitochon_13C001D010_13C001D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D010;
    Vr_SerPool2_SerPoolMitochon_13C010D010_13C010D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D010;
    Vr_SerPool2_SerPoolMitochon_13C011D010_13C011D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D010;
    Vr_SerPool2_SerPoolMitochon_13C100D010_13C100D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D010;
    Vr_SerPool2_SerPoolMitochon_13C101D010_13C101D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D010;
    Vr_SerPool2_SerPoolMitochon_13C110D010_13C110D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D010;
    Vr_SerPool2_SerPoolMitochon_13C111D010_13C111D010 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D010;
    Vr_SerPool2_SerPoolMitochon_13C000D011_13C000D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D011;
    Vr_SerPool2_SerPoolMitochon_13C001D011_13C001D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D011;
    Vr_SerPool2_SerPoolMitochon_13C010D011_13C010D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D011;
    Vr_SerPool2_SerPoolMitochon_13C011D011_13C011D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D011;
    Vr_SerPool2_SerPoolMitochon_13C100D011_13C100D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D011;
    Vr_SerPool2_SerPoolMitochon_13C101D011_13C101D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D011;
    Vr_SerPool2_SerPoolMitochon_13C110D011_13C110D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D011;
    Vr_SerPool2_SerPoolMitochon_13C111D011_13C111D011 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D011;
    Vr_SerPool2_SerPoolMitochon_13C000D100_13C000D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D100;
    Vr_SerPool2_SerPoolMitochon_13C001D100_13C001D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D100;
    Vr_SerPool2_SerPoolMitochon_13C010D100_13C010D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D100;
    Vr_SerPool2_SerPoolMitochon_13C011D100_13C011D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D100;
    Vr_SerPool2_SerPoolMitochon_13C100D100_13C100D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D100;
    Vr_SerPool2_SerPoolMitochon_13C101D100_13C101D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D100;
    Vr_SerPool2_SerPoolMitochon_13C110D100_13C110D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D100;
    Vr_SerPool2_SerPoolMitochon_13C111D100_13C111D100 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D100;
    Vr_SerPool2_SerPoolMitochon_13C000D101_13C000D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D101;
    Vr_SerPool2_SerPoolMitochon_13C001D101_13C001D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D101;
    Vr_SerPool2_SerPoolMitochon_13C010D101_13C010D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D101;
    Vr_SerPool2_SerPoolMitochon_13C011D101_13C011D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D101;
    Vr_SerPool2_SerPoolMitochon_13C100D101_13C100D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D101;
    Vr_SerPool2_SerPoolMitochon_13C101D101_13C101D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D101;
    Vr_SerPool2_SerPoolMitochon_13C110D101_13C110D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D101;
    Vr_SerPool2_SerPoolMitochon_13C111D101_13C111D101 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D101;
    Vr_SerPool2_SerPoolMitochon_13C000D110_13C000D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D110;
    Vr_SerPool2_SerPoolMitochon_13C001D110_13C001D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D110;
    Vr_SerPool2_SerPoolMitochon_13C010D110_13C010D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D110;
    Vr_SerPool2_SerPoolMitochon_13C011D110_13C011D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D110;
    Vr_SerPool2_SerPoolMitochon_13C100D110_13C100D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D110;
    Vr_SerPool2_SerPoolMitochon_13C101D110_13C101D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D110;
    Vr_SerPool2_SerPoolMitochon_13C110D110_13C110D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D110;
    Vr_SerPool2_SerPoolMitochon_13C111D110_13C111D110 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D110;
    Vr_SerPool2_SerPoolMitochon_13C000D111_13C000D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C000D111;
    Vr_SerPool2_SerPoolMitochon_13C001D111_13C001D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C001D111;
    Vr_SerPool2_SerPoolMitochon_13C010D111_13C010D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C010D111;
    Vr_SerPool2_SerPoolMitochon_13C011D111_13C011D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C011D111;
    Vr_SerPool2_SerPoolMitochon_13C100D111_13C100D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C100D111;
    Vr_SerPool2_SerPoolMitochon_13C101D111_13C101D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C101D111;
    Vr_SerPool2_SerPoolMitochon_13C110D111_13C110D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C110D111;
    Vr_SerPool2_SerPoolMitochon_13C111D111_13C111D111 = Kr_SerPool2_SerPoolMitochon * mSerPool_13C111D111;



    Vf_GlyPool2_GlyPoolMitochon_13C00D00_13C00D00 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C00D00;
    Vf_GlyPool2_GlyPoolMitochon_13C01D00_13C01D00 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C01D00;
    Vf_GlyPool2_GlyPoolMitochon_13C10D00_13C10D00 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C10D00;
    Vf_GlyPool2_GlyPoolMitochon_13C11D00_13C11D00 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C11D00;
    Vf_GlyPool2_GlyPoolMitochon_13C00D01_13C00D01 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C00D01;
    Vf_GlyPool2_GlyPoolMitochon_13C01D01_13C01D01 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C01D01;
    Vf_GlyPool2_GlyPoolMitochon_13C10D01_13C10D01 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C10D01;
    Vf_GlyPool2_GlyPoolMitochon_13C11D01_13C11D01 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C11D01;
    Vf_GlyPool2_GlyPoolMitochon_13C00D10_13C00D10 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C00D10;
    Vf_GlyPool2_GlyPoolMitochon_13C01D10_13C01D10 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C01D10;
    Vf_GlyPool2_GlyPoolMitochon_13C10D10_13C10D10 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C10D10;
    Vf_GlyPool2_GlyPoolMitochon_13C11D10_13C11D10 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C11D10;
    Vf_GlyPool2_GlyPoolMitochon_13C00D11_13C00D11 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C00D11;
    Vf_GlyPool2_GlyPoolMitochon_13C01D11_13C01D11 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C01D11;
    Vf_GlyPool2_GlyPoolMitochon_13C10D11_13C10D11 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C10D11;
    Vf_GlyPool2_GlyPoolMitochon_13C11D11_13C11D11 = Kf_GlyPool2_GlyPoolMitochon * cGlyPool2_13C11D11;



    Vr_GlyPool2_GlyPoolMitochon_13C00D00_13C00D00 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C00D00;
    Vr_GlyPool2_GlyPoolMitochon_13C01D00_13C01D00 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C01D00;
    Vr_GlyPool2_GlyPoolMitochon_13C10D00_13C10D00 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C10D00;
    Vr_GlyPool2_GlyPoolMitochon_13C11D00_13C11D00 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C11D00;
    Vr_GlyPool2_GlyPoolMitochon_13C00D01_13C00D01 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C00D01;
    Vr_GlyPool2_GlyPoolMitochon_13C01D01_13C01D01 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C01D01;
    Vr_GlyPool2_GlyPoolMitochon_13C10D01_13C10D01 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C10D01;
    Vr_GlyPool2_GlyPoolMitochon_13C11D01_13C11D01 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C11D01;
    Vr_GlyPool2_GlyPoolMitochon_13C00D10_13C00D10 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C00D10;
    Vr_GlyPool2_GlyPoolMitochon_13C01D10_13C01D10 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C01D10;
    Vr_GlyPool2_GlyPoolMitochon_13C10D10_13C10D10 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C10D10;
    Vr_GlyPool2_GlyPoolMitochon_13C11D10_13C11D10 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C11D10;
    Vr_GlyPool2_GlyPoolMitochon_13C00D11_13C00D11 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C00D11;
    Vr_GlyPool2_GlyPoolMitochon_13C01D11_13C01D11 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C01D11;
    Vr_GlyPool2_GlyPoolMitochon_13C10D11_13C10D11 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C10D11;
    Vr_GlyPool2_GlyPoolMitochon_13C11D11_13C11D11 = Kr_GlyPool2_GlyPoolMitochon * mGlyPool_13C11D11;



    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D00_13C0D0 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C0D00;
    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D00_13C1D0 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C1D00;
    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D01_13C0D0 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C0D01;
    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D01_13C1D0 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C1D01;
    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D10_13C0D1 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C0D10;
    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D10_13C1D1 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C1D10;
    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D11_13C0D1 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C0D11;
    Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D11_13C1D1 = Kf_MethyleneTHF_FormylTHF_Cytoplasm * cMethyleneTHF_13C1D11;



    Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C0D0_13C0D00 = Kr_MethyleneTHF_FormylTHF_Cytoplasm * cFormylTHF_13C0D0;
    Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C1D0_13C1D00 = Kr_MethyleneTHF_FormylTHF_Cytoplasm * cFormylTHF_13C1D0;
    Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C0D1_13C0D10 = Kr_MethyleneTHF_FormylTHF_Cytoplasm * cFormylTHF_13C0D1;
    Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C1D1_13C1D10 = Kr_MethyleneTHF_FormylTHF_Cytoplasm * cFormylTHF_13C1D1;



    Vf_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0 = Kf_FormylTHF_Formate_Cytoplasm * cFormylTHF_13C0D0;
    Vf_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0 = Kf_FormylTHF_Formate_Cytoplasm * cFormylTHF_13C1D0;
    Vf_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1 = Kf_FormylTHF_Formate_Cytoplasm * cFormylTHF_13C0D1;
    Vf_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1 = Kf_FormylTHF_Formate_Cytoplasm * cFormylTHF_13C1D1;



    Vr_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0 = Kr_FormylTHF_Formate_Cytoplasm * cTHF / (6.02214129 * 10^11) / CytoplasmVolume * cFormate_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vr_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0 = Kr_FormylTHF_Formate_Cytoplasm * cTHF / (6.02214129 * 10^11) / CytoplasmVolume * cFormate_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vr_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1 = Kr_FormylTHF_Formate_Cytoplasm * cTHF / (6.02214129 * 10^11) / CytoplasmVolume * cFormate_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vr_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1 = Kr_FormylTHF_Formate_Cytoplasm * cTHF / (6.02214129 * 10^11) / CytoplasmVolume * cFormate_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vf_Formate_Transport_13C0D0 = Kf_Formate_Transport * cFormate_13C0D0;
    Vf_Formate_Transport_13C1D0 = Kf_Formate_Transport * cFormate_13C1D0;
    Vf_Formate_Transport_13C0D1 = Kf_Formate_Transport * cFormate_13C0D1;
    Vf_Formate_Transport_13C1D1 = Kf_Formate_Transport * cFormate_13C1D1;



    Vr_Formate_Transport_13C0D0 = Kr_Formate_Transport * mFormate_13C0D0;
    Vr_Formate_Transport_13C1D0 = Kr_Formate_Transport * mFormate_13C1D0;
    Vr_Formate_Transport_13C0D1 = Kr_Formate_Transport * mFormate_13C0D1;
    Vr_Formate_Transport_13C1D1 = Kr_Formate_Transport * mFormate_13C1D1;



    Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D00_13C0D0 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C0D00;
    Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D00_13C1D0 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C1D00;
    Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D10_13C0D1 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C0D10;
    Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D10_13C1D1 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C1D10;
    Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D01_13C0D0 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C0D01;
    Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D01_13C1D0 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C1D01;
    Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D11_13C0D1 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C0D11;
    Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D11_13C1D1 = Kf_MethyleneTHF_FormylTHF_Mitochon * mMethyleneTHF_13C1D11;



    Vr_MethyleneTHF_FormylTHF_Mitochon_13C0D0_13C0D00 = Kr_MethyleneTHF_FormylTHF_Mitochon * mFormylTHF_13C0D0;
    Vr_MethyleneTHF_FormylTHF_Mitochon_13C1D0_13C1D00 = Kr_MethyleneTHF_FormylTHF_Mitochon * mFormylTHF_13C1D0;
    Vr_MethyleneTHF_FormylTHF_Mitochon_13C0D1_13C0D10 = Kr_MethyleneTHF_FormylTHF_Mitochon * mFormylTHF_13C0D1;
    Vr_MethyleneTHF_FormylTHF_Mitochon_13C1D1_13C1D10 = Kr_MethyleneTHF_FormylTHF_Mitochon * mFormylTHF_13C1D1;



    Vf_FormylTHF_Formate_Mitochon_13C0D0_13C0D0 = Kf_FormylTHF_Formate_Mitochon * mFormylTHF_13C0D0;
    Vf_FormylTHF_Formate_Mitochon_13C1D0_13C1D0 = Kf_FormylTHF_Formate_Mitochon * mFormylTHF_13C1D0;
    Vf_FormylTHF_Formate_Mitochon_13C0D1_13C0D1 = Kf_FormylTHF_Formate_Mitochon * mFormylTHF_13C0D1;
    Vf_FormylTHF_Formate_Mitochon_13C1D1_13C1D1 = Kf_FormylTHF_Formate_Mitochon * mFormylTHF_13C1D1;



    Vr_FormylTHF_Formate_Mitochon_13C0D0_13C0D0 = Kr_FormylTHF_Formate_Mitochon * mTHF / (6.02214129 * 10^11) / MitochonVolume * mFormate_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vr_FormylTHF_Formate_Mitochon_13C1D0_13C1D0 = Kr_FormylTHF_Formate_Mitochon * mTHF / (6.02214129 * 10^11) / MitochonVolume * mFormate_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vr_FormylTHF_Formate_Mitochon_13C0D1_13C0D1 = Kr_FormylTHF_Formate_Mitochon * mTHF / (6.02214129 * 10^11) / MitochonVolume * mFormate_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vr_FormylTHF_Formate_Mitochon_13C1D1_13C1D1 = Kr_FormylTHF_Formate_Mitochon * mTHF / (6.02214129 * 10^11) / MitochonVolume * mFormate_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vf_MethyleneTHF_Transport_13C0D00_13C0D00 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C0D00;
    Vf_MethyleneTHF_Transport_13C1D00_13C1D00 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C1D00;
    Vf_MethyleneTHF_Transport_13C0D01_13C0D01 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C0D01;
    Vf_MethyleneTHF_Transport_13C1D01_13C1D01 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C1D01;
    Vf_MethyleneTHF_Transport_13C0D10_13C0D10 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C0D10;
    Vf_MethyleneTHF_Transport_13C1D10_13C1D10 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C1D10;
    Vf_MethyleneTHF_Transport_13C0D11_13C0D11 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C0D11;
    Vf_MethyleneTHF_Transport_13C1D11_13C1D11 = Kf_MethyleneTHF_Transport * cMethyleneTHF_13C1D11;



    Vr_MethyleneTHF_Transport_13C0D00_13C0D00 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C0D00;
    Vr_MethyleneTHF_Transport_13C1D00_13C1D00 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C1D00;
    Vr_MethyleneTHF_Transport_13C0D01_13C0D01 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C0D01;
    Vr_MethyleneTHF_Transport_13C1D01_13C1D01 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C1D01;
    Vr_MethyleneTHF_Transport_13C0D10_13C0D10 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C0D10;
    Vr_MethyleneTHF_Transport_13C1D10_13C1D10 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C1D10;
    Vr_MethyleneTHF_Transport_13C0D11_13C0D11 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C0D11;
    Vr_MethyleneTHF_Transport_13C1D11_13C1D11 = Kr_MethyleneTHF_Transport * mMethyleneTHF_13C1D11;



    Vf_FormylTHF_Transport_13C0D0_13C0D0 = Kf_FormylTHF_Transport * cFormylTHF_13C0D0;
    Vf_FormylTHF_Transport_13C1D0_13C1D0 = Kf_FormylTHF_Transport * cFormylTHF_13C1D0;
    Vf_FormylTHF_Transport_13C0D1_13C0D1 = Kf_FormylTHF_Transport * cFormylTHF_13C0D1;
    Vf_FormylTHF_Transport_13C1D1_13C1D1 = Kf_FormylTHF_Transport * cFormylTHF_13C1D1;



    Vr_FormylTHF_Transport_13C0D0_13C0D0 = Kr_FormylTHF_Transport * mFormylTHF_13C0D0;
    Vr_FormylTHF_Transport_13C1D0_13C1D0 = Kr_FormylTHF_Transport * mFormylTHF_13C1D0;
    Vr_FormylTHF_Transport_13C0D1_13C0D1 = Kr_FormylTHF_Transport * mFormylTHF_13C0D1;
    Vr_FormylTHF_Transport_13C1D1_13C1D1 = Kr_FormylTHF_Transport * mFormylTHF_13C1D1;



    Vf_THF_Transport = Kf_THF_Transport * cTHF;



    Vr_THF_Transport = Kr_THF_Transport * mTHF;



    Vf_PRPP_GAR_GlyPool1_13C0_13C00D00__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C00D00__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C01D00__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C01D00__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C10D00__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C10D00__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C11D00__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C11D00__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C00D01__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C00D01__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C01D01__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C01D01__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C10D01__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C10D01__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C11D01__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C11D01__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C00D10__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C00D10__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C01D10__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C01D10__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C10D10__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C10D10__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C11D10__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C11D10__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C00D11__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C00D11__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C00D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C01D11__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C01D11__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C01D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C10D11__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C10D11__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C10D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C0_13C11D11__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool1_13C1_13C11D11__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool1_13C11D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vf_PRPP_GAR_GlyPool2_13C0_13C00D00__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C00D00__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C01D00__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C01D00__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C10D00__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C10D00__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C11D00__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C11D00__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C00D01__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C00D01__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C01D01__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C01D01__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C10D01__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C10D01__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C11D01__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C11D01__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C00D10__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C00D10__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C01D10__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C01D10__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C10D10__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C10D10__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C11D10__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C11D10__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C00D11__13C000 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C00D11__13C100 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C00D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C01D11__13C001 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C01D11__13C101 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C01D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C10D11__13C010 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C10D11__13C110 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C10D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C0_13C11D11__13C011 = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool2_13C1_13C11D11__13C111 = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool2_13C11D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vf_PRPP_GAR_GlyPool3_13C0_13C00D00__13C000  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C00D00__13C100  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C01D00__13C001  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C01D00__13C101  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C10D00__13C010  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C10D00__13C110  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C11D00__13C011  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C11D00__13C111  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D00 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C00D01__13C000  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C00D01__13C100  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C01D01__13C001  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C01D01__13C101  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C10D01__13C010  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C10D01__13C110  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C11D01__13C011  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C11D01__13C111  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D01 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C00D10__13C000  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C00D10__13C100  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C01D10__13C001  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C01D10__13C101  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C10D10__13C010  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C10D10__13C110  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C11D10__13C011  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C11D10__13C111  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D10 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C00D11__13C000  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C00D11__13C100  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C00D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C01D11__13C001  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C01D11__13C101  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C01D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C10D11__13C010  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C10D11__13C110  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C10D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C0_13C11D11__13C011  = Kf_PRPP1_GAR * cPRPP_13C0 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_PRPP_GAR_GlyPool3_13C1_13C11D11__13C111  = Kf_PRPP1_GAR * cPRPP_13C1 / (6.02214129 * 10^11) / CytoplasmVolume * cGlyPool3_13C11D11 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vf_GAR_FGAR_13C000_13C0D0__13C0000D0 = Kf_GAR_FGAR * cGAR_13C000 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C001_13C0D0__13C0010D0 = Kf_GAR_FGAR * cGAR_13C001 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C010_13C0D0__13C0100D0 = Kf_GAR_FGAR * cGAR_13C010 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C011_13C0D0__13C0110D0 = Kf_GAR_FGAR * cGAR_13C011 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C100_13C0D0__13C1000D0 = Kf_GAR_FGAR * cGAR_13C100 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C101_13C0D0__13C1010D0 = Kf_GAR_FGAR * cGAR_13C101 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C110_13C0D0__13C1100D0 = Kf_GAR_FGAR * cGAR_13C110 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C111_13C0D0__13C1110D0 = Kf_GAR_FGAR * cGAR_13C111 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C000_13C1D0__13C0001D0 = Kf_GAR_FGAR * cGAR_13C000 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C001_13C1D0__13C0011D0 = Kf_GAR_FGAR * cGAR_13C001 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C010_13C1D0__13C0101D0 = Kf_GAR_FGAR * cGAR_13C010 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C011_13C1D0__13C0111D0 = Kf_GAR_FGAR * cGAR_13C011 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C100_13C1D0__13C1001D0 = Kf_GAR_FGAR * cGAR_13C100 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C101_13C1D0__13C1011D0 = Kf_GAR_FGAR * cGAR_13C101 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C110_13C1D0__13C1101D0 = Kf_GAR_FGAR * cGAR_13C110 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C111_13C1D0__13C1111D0 = Kf_GAR_FGAR * cGAR_13C111 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C000_13C0D1__13C0000D1 = Kf_GAR_FGAR * cGAR_13C000 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C001_13C0D1__13C0010D1 = Kf_GAR_FGAR * cGAR_13C001 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C010_13C0D1__13C0100D1 = Kf_GAR_FGAR * cGAR_13C010 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C011_13C0D1__13C0110D1 = Kf_GAR_FGAR * cGAR_13C011 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C100_13C0D1__13C1000D1 = Kf_GAR_FGAR * cGAR_13C100 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C101_13C0D1__13C1010D1 = Kf_GAR_FGAR * cGAR_13C101 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C110_13C0D1__13C1100D1 = Kf_GAR_FGAR * cGAR_13C110 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C111_13C0D1__13C1110D1 = Kf_GAR_FGAR * cGAR_13C111 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C000_13C1D1__13C0001D1 = Kf_GAR_FGAR * cGAR_13C000 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C001_13C1D1__13C0011D1 = Kf_GAR_FGAR * cGAR_13C001 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C010_13C1D1__13C0101D1 = Kf_GAR_FGAR * cGAR_13C010 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C011_13C1D1__13C0111D1 = Kf_GAR_FGAR * cGAR_13C011 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C100_13C1D1__13C1001D1 = Kf_GAR_FGAR * cGAR_13C100 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C101_13C1D1__13C1011D1 = Kf_GAR_FGAR * cGAR_13C101 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C110_13C1D1__13C1101D1 = Kf_GAR_FGAR * cGAR_13C110 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_GAR_FGAR_13C111_13C1D1__13C1111D1 = Kf_GAR_FGAR * cGAR_13C111 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);




    Vf_FGAR_AMP_13C0000D0_13C0D0__13C00000D00 = Kf_FGAR_AMP * cFGAR_13C0000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D0_13C0D0__13C00010D00 = Kf_FGAR_AMP * cFGAR_13C0001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D0_13C0D0__13C00100D00 = Kf_FGAR_AMP * cFGAR_13C0010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D0_13C0D0__13C00110D00 = Kf_FGAR_AMP * cFGAR_13C0011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D0_13C0D0__13C01000D00 = Kf_FGAR_AMP * cFGAR_13C0100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D0_13C0D0__13C01010D00 = Kf_FGAR_AMP * cFGAR_13C0101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D0_13C0D0__13C01100D00 = Kf_FGAR_AMP * cFGAR_13C0110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D0_13C0D0__13C01110D00 = Kf_FGAR_AMP * cFGAR_13C0111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D0_13C0D0__13C10000D00 = Kf_FGAR_AMP * cFGAR_13C1000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D0_13C0D0__13C10010D00 = Kf_FGAR_AMP * cFGAR_13C1001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D0_13C0D0__13C10100D00 = Kf_FGAR_AMP * cFGAR_13C1010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D0_13C0D0__13C10110D00 = Kf_FGAR_AMP * cFGAR_13C1011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D0_13C0D0__13C11000D00 = Kf_FGAR_AMP * cFGAR_13C1100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D0_13C0D0__13C11010D00 = Kf_FGAR_AMP * cFGAR_13C1101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D0_13C0D0__13C11100D00 = Kf_FGAR_AMP * cFGAR_13C1110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D0_13C0D0__13C11110D00 = Kf_FGAR_AMP * cFGAR_13C1111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0000D1_13C0D0__13C00000D10 = Kf_FGAR_AMP * cFGAR_13C0000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D1_13C0D0__13C00010D10 = Kf_FGAR_AMP * cFGAR_13C0001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D1_13C0D0__13C00100D10 = Kf_FGAR_AMP * cFGAR_13C0010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D1_13C0D0__13C00110D10 = Kf_FGAR_AMP * cFGAR_13C0011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D1_13C0D0__13C01000D10 = Kf_FGAR_AMP * cFGAR_13C0100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D1_13C0D0__13C01010D10 = Kf_FGAR_AMP * cFGAR_13C0101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D1_13C0D0__13C01100D10 = Kf_FGAR_AMP * cFGAR_13C0110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D1_13C0D0__13C01110D10 = Kf_FGAR_AMP * cFGAR_13C0111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D1_13C0D0__13C10000D10 = Kf_FGAR_AMP * cFGAR_13C1000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D1_13C0D0__13C10010D10 = Kf_FGAR_AMP * cFGAR_13C1001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D1_13C0D0__13C10100D10 = Kf_FGAR_AMP * cFGAR_13C1010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D1_13C0D0__13C10110D10 = Kf_FGAR_AMP * cFGAR_13C1011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D1_13C0D0__13C11000D10 = Kf_FGAR_AMP * cFGAR_13C1100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D1_13C0D0__13C11010D10 = Kf_FGAR_AMP * cFGAR_13C1101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D1_13C0D0__13C11100D10 = Kf_FGAR_AMP * cFGAR_13C1110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D1_13C0D0__13C11110D10 = Kf_FGAR_AMP * cFGAR_13C1111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0000D0_13C1D0__13C00001D00 = Kf_FGAR_AMP * cFGAR_13C0000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D0_13C1D0__13C00011D00 = Kf_FGAR_AMP * cFGAR_13C0001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D0_13C1D0__13C00101D00 = Kf_FGAR_AMP * cFGAR_13C0010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D0_13C1D0__13C00111D00 = Kf_FGAR_AMP * cFGAR_13C0011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D0_13C1D0__13C01001D00 = Kf_FGAR_AMP * cFGAR_13C0100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D0_13C1D0__13C01011D00 = Kf_FGAR_AMP * cFGAR_13C0101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D0_13C1D0__13C01101D00 = Kf_FGAR_AMP * cFGAR_13C0110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D0_13C1D0__13C01111D00 = Kf_FGAR_AMP * cFGAR_13C0111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D0_13C1D0__13C10001D00 = Kf_FGAR_AMP * cFGAR_13C1000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D0_13C1D0__13C10011D00 = Kf_FGAR_AMP * cFGAR_13C1001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D0_13C1D0__13C10101D00 = Kf_FGAR_AMP * cFGAR_13C1010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D0_13C1D0__13C10111D00 = Kf_FGAR_AMP * cFGAR_13C1011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D0_13C1D0__13C11001D00 = Kf_FGAR_AMP * cFGAR_13C1100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D0_13C1D0__13C11011D00 = Kf_FGAR_AMP * cFGAR_13C1101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D0_13C1D0__13C11101D00 = Kf_FGAR_AMP * cFGAR_13C1110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D0_13C1D0__13C11111D00 = Kf_FGAR_AMP * cFGAR_13C1111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0000D1_13C1D0__13C00001D10 = Kf_FGAR_AMP * cFGAR_13C0000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D1_13C1D0__13C00011D10 = Kf_FGAR_AMP * cFGAR_13C0001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D1_13C1D0__13C00101D10 = Kf_FGAR_AMP * cFGAR_13C0010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D1_13C1D0__13C00111D10 = Kf_FGAR_AMP * cFGAR_13C0011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D1_13C1D0__13C01001D10 = Kf_FGAR_AMP * cFGAR_13C0100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D1_13C1D0__13C01011D10 = Kf_FGAR_AMP * cFGAR_13C0101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D1_13C1D0__13C01101D10 = Kf_FGAR_AMP * cFGAR_13C0110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D1_13C1D0__13C01111D10 = Kf_FGAR_AMP * cFGAR_13C0111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D1_13C1D0__13C10001D10 = Kf_FGAR_AMP * cFGAR_13C1000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D1_13C1D0__13C10011D10 = Kf_FGAR_AMP * cFGAR_13C1001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D1_13C1D0__13C10101D10 = Kf_FGAR_AMP * cFGAR_13C1010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D1_13C1D0__13C10111D10 = Kf_FGAR_AMP * cFGAR_13C1011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D1_13C1D0__13C11001D10 = Kf_FGAR_AMP * cFGAR_13C1100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D1_13C1D0__13C11011D10 = Kf_FGAR_AMP * cFGAR_13C1101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D1_13C1D0__13C11101D10 = Kf_FGAR_AMP * cFGAR_13C1110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D1_13C1D0__13C11111D10 = Kf_FGAR_AMP * cFGAR_13C1111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D0 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0000D0_13C0D1__13C00000D01 = Kf_FGAR_AMP * cFGAR_13C0000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D0_13C0D1__13C00010D01 = Kf_FGAR_AMP * cFGAR_13C0001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D0_13C0D1__13C00100D01 = Kf_FGAR_AMP * cFGAR_13C0010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D0_13C0D1__13C00110D01 = Kf_FGAR_AMP * cFGAR_13C0011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D0_13C0D1__13C01000D01 = Kf_FGAR_AMP * cFGAR_13C0100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D0_13C0D1__13C01010D01 = Kf_FGAR_AMP * cFGAR_13C0101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D0_13C0D1__13C01100D01 = Kf_FGAR_AMP * cFGAR_13C0110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D0_13C0D1__13C01110D01 = Kf_FGAR_AMP * cFGAR_13C0111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D0_13C0D1__13C10000D01 = Kf_FGAR_AMP * cFGAR_13C1000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D0_13C0D1__13C10010D01 = Kf_FGAR_AMP * cFGAR_13C1001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D0_13C0D1__13C10100D01 = Kf_FGAR_AMP * cFGAR_13C1010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D0_13C0D1__13C10110D01 = Kf_FGAR_AMP * cFGAR_13C1011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D0_13C0D1__13C11000D01 = Kf_FGAR_AMP * cFGAR_13C1100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D0_13C0D1__13C11010D01 = Kf_FGAR_AMP * cFGAR_13C1101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D0_13C0D1__13C11100D01 = Kf_FGAR_AMP * cFGAR_13C1110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D0_13C0D1__13C11110D01 = Kf_FGAR_AMP * cFGAR_13C1111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0000D1_13C0D1__13C00000D11 = Kf_FGAR_AMP * cFGAR_13C0000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D1_13C0D1__13C00010D11 = Kf_FGAR_AMP * cFGAR_13C0001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D1_13C0D1__13C00100D11 = Kf_FGAR_AMP * cFGAR_13C0010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D1_13C0D1__13C00110D11 = Kf_FGAR_AMP * cFGAR_13C0011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D1_13C0D1__13C01000D11 = Kf_FGAR_AMP * cFGAR_13C0100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D1_13C0D1__13C01010D11 = Kf_FGAR_AMP * cFGAR_13C0101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D1_13C0D1__13C01100D11 = Kf_FGAR_AMP * cFGAR_13C0110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D1_13C0D1__13C01110D11 = Kf_FGAR_AMP * cFGAR_13C0111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D1_13C0D1__13C10000D11 = Kf_FGAR_AMP * cFGAR_13C1000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D1_13C0D1__13C10010D11 = Kf_FGAR_AMP * cFGAR_13C1001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D1_13C0D1__13C10100D11 = Kf_FGAR_AMP * cFGAR_13C1010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D1_13C0D1__13C10110D11 = Kf_FGAR_AMP * cFGAR_13C1011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D1_13C0D1__13C11000D11 = Kf_FGAR_AMP * cFGAR_13C1100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D1_13C0D1__13C11010D11 = Kf_FGAR_AMP * cFGAR_13C1101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D1_13C0D1__13C11100D11 = Kf_FGAR_AMP * cFGAR_13C1110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D1_13C0D1__13C11110D11 = Kf_FGAR_AMP * cFGAR_13C1111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C0D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0000D0_13C1D1__13C00001D01 = Kf_FGAR_AMP * cFGAR_13C0000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D0_13C1D1__13C00011D01 = Kf_FGAR_AMP * cFGAR_13C0001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D0_13C1D1__13C00101D01 = Kf_FGAR_AMP * cFGAR_13C0010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D0_13C1D1__13C00111D01 = Kf_FGAR_AMP * cFGAR_13C0011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D0_13C1D1__13C01001D01 = Kf_FGAR_AMP * cFGAR_13C0100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D0_13C1D1__13C01011D01 = Kf_FGAR_AMP * cFGAR_13C0101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D0_13C1D1__13C01101D01 = Kf_FGAR_AMP * cFGAR_13C0110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D0_13C1D1__13C01111D01 = Kf_FGAR_AMP * cFGAR_13C0111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D0_13C1D1__13C10001D01 = Kf_FGAR_AMP * cFGAR_13C1000D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D0_13C1D1__13C10011D01 = Kf_FGAR_AMP * cFGAR_13C1001D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D0_13C1D1__13C10101D01 = Kf_FGAR_AMP * cFGAR_13C1010D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D0_13C1D1__13C10111D01 = Kf_FGAR_AMP * cFGAR_13C1011D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D0_13C1D1__13C11001D01 = Kf_FGAR_AMP * cFGAR_13C1100D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D0_13C1D1__13C11011D01 = Kf_FGAR_AMP * cFGAR_13C1101D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D0_13C1D1__13C11101D01 = Kf_FGAR_AMP * cFGAR_13C1110D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D0_13C1D1__13C11111D01 = Kf_FGAR_AMP * cFGAR_13C1111D0 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0000D1_13C1D1__13C00001D11 = Kf_FGAR_AMP * cFGAR_13C0000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0001D1_13C1D1__13C00011D11 = Kf_FGAR_AMP * cFGAR_13C0001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0010D1_13C1D1__13C00101D11 = Kf_FGAR_AMP * cFGAR_13C0010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0011D1_13C1D1__13C00111D11 = Kf_FGAR_AMP * cFGAR_13C0011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0100D1_13C1D1__13C01001D11 = Kf_FGAR_AMP * cFGAR_13C0100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0101D1_13C1D1__13C01011D11 = Kf_FGAR_AMP * cFGAR_13C0101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0110D1_13C1D1__13C01101D11 = Kf_FGAR_AMP * cFGAR_13C0110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C0111D1_13C1D1__13C01111D11 = Kf_FGAR_AMP * cFGAR_13C0111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1000D1_13C1D1__13C10001D11 = Kf_FGAR_AMP * cFGAR_13C1000D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1001D1_13C1D1__13C10011D11 = Kf_FGAR_AMP * cFGAR_13C1001D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1010D1_13C1D1__13C10101D11 = Kf_FGAR_AMP * cFGAR_13C1010D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1011D1_13C1D1__13C10111D11 = Kf_FGAR_AMP * cFGAR_13C1011D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1100D1_13C1D1__13C11001D11 = Kf_FGAR_AMP * cFGAR_13C1100D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1101D1_13C1D1__13C11011D11 = Kf_FGAR_AMP * cFGAR_13C1101D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1110D1_13C1D1__13C11101D11 = Kf_FGAR_AMP * cFGAR_13C1110D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);
    Vf_FGAR_AMP_13C1111D1_13C1D1__13C11111D11 = Kf_FGAR_AMP * cFGAR_13C1111D1 / (6.02214129 * 10^11) / CytoplasmVolume * cFormylTHF_13C1D1 / (6.02214129 * 10^11) * (6.02214129 * 10^5);



    Vf_AMP_Degradation_13C00000D00_13C00000D00 = Kf_AMP_Degradation * cAMP_13C00000D00;
    Vf_AMP_Degradation_13C00001D00_13C00001D00 = Kf_AMP_Degradation * cAMP_13C00001D00;
    Vf_AMP_Degradation_13C00010D00_13C00010D00 = Kf_AMP_Degradation * cAMP_13C00010D00;
    Vf_AMP_Degradation_13C00011D00_13C00011D00 = Kf_AMP_Degradation * cAMP_13C00011D00;
    Vf_AMP_Degradation_13C00100D00_13C00100D00 = Kf_AMP_Degradation * cAMP_13C00100D00;
    Vf_AMP_Degradation_13C00101D00_13C00101D00 = Kf_AMP_Degradation * cAMP_13C00101D00;
    Vf_AMP_Degradation_13C00110D00_13C00110D00 = Kf_AMP_Degradation * cAMP_13C00110D00;
    Vf_AMP_Degradation_13C00111D00_13C00111D00 = Kf_AMP_Degradation * cAMP_13C00111D00;
    Vf_AMP_Degradation_13C01000D00_13C01000D00 = Kf_AMP_Degradation * cAMP_13C01000D00;
    Vf_AMP_Degradation_13C01001D00_13C01001D00 = Kf_AMP_Degradation * cAMP_13C01001D00;
    Vf_AMP_Degradation_13C01010D00_13C01010D00 = Kf_AMP_Degradation * cAMP_13C01010D00;
    Vf_AMP_Degradation_13C01011D00_13C01011D00 = Kf_AMP_Degradation * cAMP_13C01011D00;
    Vf_AMP_Degradation_13C01100D00_13C01100D00 = Kf_AMP_Degradation * cAMP_13C01100D00;
    Vf_AMP_Degradation_13C01101D00_13C01101D00 = Kf_AMP_Degradation * cAMP_13C01101D00;
    Vf_AMP_Degradation_13C01110D00_13C01110D00 = Kf_AMP_Degradation * cAMP_13C01110D00;
    Vf_AMP_Degradation_13C01111D00_13C01111D00 = Kf_AMP_Degradation * cAMP_13C01111D00;
    Vf_AMP_Degradation_13C10000D00_13C10000D00 = Kf_AMP_Degradation * cAMP_13C10000D00;
    Vf_AMP_Degradation_13C10001D00_13C10001D00 = Kf_AMP_Degradation * cAMP_13C10001D00;
    Vf_AMP_Degradation_13C10010D00_13C10010D00 = Kf_AMP_Degradation * cAMP_13C10010D00;
    Vf_AMP_Degradation_13C10011D00_13C10011D00 = Kf_AMP_Degradation * cAMP_13C10011D00;
    Vf_AMP_Degradation_13C10100D00_13C10100D00 = Kf_AMP_Degradation * cAMP_13C10100D00;
    Vf_AMP_Degradation_13C10101D00_13C10101D00 = Kf_AMP_Degradation * cAMP_13C10101D00;
    Vf_AMP_Degradation_13C10110D00_13C10110D00 = Kf_AMP_Degradation * cAMP_13C10110D00;
    Vf_AMP_Degradation_13C10111D00_13C10111D00 = Kf_AMP_Degradation * cAMP_13C10111D00;
    Vf_AMP_Degradation_13C11000D00_13C11000D00 = Kf_AMP_Degradation * cAMP_13C11000D00;
    Vf_AMP_Degradation_13C11001D00_13C11001D00 = Kf_AMP_Degradation * cAMP_13C11001D00;
    Vf_AMP_Degradation_13C11010D00_13C11010D00 = Kf_AMP_Degradation * cAMP_13C11010D00;
    Vf_AMP_Degradation_13C11011D00_13C11011D00 = Kf_AMP_Degradation * cAMP_13C11011D00;
    Vf_AMP_Degradation_13C11100D00_13C11100D00 = Kf_AMP_Degradation * cAMP_13C11100D00;
    Vf_AMP_Degradation_13C11101D00_13C11101D00 = Kf_AMP_Degradation * cAMP_13C11101D00;
    Vf_AMP_Degradation_13C11110D00_13C11110D00 = Kf_AMP_Degradation * cAMP_13C11110D00;
    Vf_AMP_Degradation_13C11111D00_13C11111D00 = Kf_AMP_Degradation * cAMP_13C11111D00;
    Vf_AMP_Degradation_13C00000D01_13C00000D01 = Kf_AMP_Degradation * cAMP_13C00000D01;
    Vf_AMP_Degradation_13C00001D01_13C00001D01 = Kf_AMP_Degradation * cAMP_13C00001D01;
    Vf_AMP_Degradation_13C00010D01_13C00010D01 = Kf_AMP_Degradation * cAMP_13C00010D01;
    Vf_AMP_Degradation_13C00011D01_13C00011D01 = Kf_AMP_Degradation * cAMP_13C00011D01;
    Vf_AMP_Degradation_13C00100D01_13C00100D01 = Kf_AMP_Degradation * cAMP_13C00100D01;
    Vf_AMP_Degradation_13C00101D01_13C00101D01 = Kf_AMP_Degradation * cAMP_13C00101D01;
    Vf_AMP_Degradation_13C00110D01_13C00110D01 = Kf_AMP_Degradation * cAMP_13C00110D01;
    Vf_AMP_Degradation_13C00111D01_13C00111D01 = Kf_AMP_Degradation * cAMP_13C00111D01;
    Vf_AMP_Degradation_13C01000D01_13C01000D01 = Kf_AMP_Degradation * cAMP_13C01000D01;
    Vf_AMP_Degradation_13C01001D01_13C01001D01 = Kf_AMP_Degradation * cAMP_13C01001D01;
    Vf_AMP_Degradation_13C01010D01_13C01010D01 = Kf_AMP_Degradation * cAMP_13C01010D01;
    Vf_AMP_Degradation_13C01011D01_13C01011D01 = Kf_AMP_Degradation * cAMP_13C01011D01;
    Vf_AMP_Degradation_13C01100D01_13C01100D01 = Kf_AMP_Degradation * cAMP_13C01100D01;
    Vf_AMP_Degradation_13C01101D01_13C01101D01 = Kf_AMP_Degradation * cAMP_13C01101D01;
    Vf_AMP_Degradation_13C01110D01_13C01110D01 = Kf_AMP_Degradation * cAMP_13C01110D01;
    Vf_AMP_Degradation_13C01111D01_13C01111D01 = Kf_AMP_Degradation * cAMP_13C01111D01;
    Vf_AMP_Degradation_13C10000D01_13C10000D01 = Kf_AMP_Degradation * cAMP_13C10000D01;
    Vf_AMP_Degradation_13C10001D01_13C10001D01 = Kf_AMP_Degradation * cAMP_13C10001D01;
    Vf_AMP_Degradation_13C10010D01_13C10010D01 = Kf_AMP_Degradation * cAMP_13C10010D01;
    Vf_AMP_Degradation_13C10011D01_13C10011D01 = Kf_AMP_Degradation * cAMP_13C10011D01;
    Vf_AMP_Degradation_13C10100D01_13C10100D01 = Kf_AMP_Degradation * cAMP_13C10100D01;
    Vf_AMP_Degradation_13C10101D01_13C10101D01 = Kf_AMP_Degradation * cAMP_13C10101D01;
    Vf_AMP_Degradation_13C10110D01_13C10110D01 = Kf_AMP_Degradation * cAMP_13C10110D01;
    Vf_AMP_Degradation_13C10111D01_13C10111D01 = Kf_AMP_Degradation * cAMP_13C10111D01;
    Vf_AMP_Degradation_13C11000D01_13C11000D01 = Kf_AMP_Degradation * cAMP_13C11000D01;
    Vf_AMP_Degradation_13C11001D01_13C11001D01 = Kf_AMP_Degradation * cAMP_13C11001D01;
    Vf_AMP_Degradation_13C11010D01_13C11010D01 = Kf_AMP_Degradation * cAMP_13C11010D01;
    Vf_AMP_Degradation_13C11011D01_13C11011D01 = Kf_AMP_Degradation * cAMP_13C11011D01;
    Vf_AMP_Degradation_13C11100D01_13C11100D01 = Kf_AMP_Degradation * cAMP_13C11100D01;
    Vf_AMP_Degradation_13C11101D01_13C11101D01 = Kf_AMP_Degradation * cAMP_13C11101D01;
    Vf_AMP_Degradation_13C11110D01_13C11110D01 = Kf_AMP_Degradation * cAMP_13C11110D01;
    Vf_AMP_Degradation_13C11111D01_13C11111D01 = Kf_AMP_Degradation * cAMP_13C11111D01;
    Vf_AMP_Degradation_13C00000D10_13C00000D10 = Kf_AMP_Degradation * cAMP_13C00000D10;
    Vf_AMP_Degradation_13C00001D10_13C00001D10 = Kf_AMP_Degradation * cAMP_13C00001D10;
    Vf_AMP_Degradation_13C00010D10_13C00010D10 = Kf_AMP_Degradation * cAMP_13C00010D10;
    Vf_AMP_Degradation_13C00011D10_13C00011D10 = Kf_AMP_Degradation * cAMP_13C00011D10;
    Vf_AMP_Degradation_13C00100D10_13C00100D10 = Kf_AMP_Degradation * cAMP_13C00100D10;
    Vf_AMP_Degradation_13C00101D10_13C00101D10 = Kf_AMP_Degradation * cAMP_13C00101D10;
    Vf_AMP_Degradation_13C00110D10_13C00110D10 = Kf_AMP_Degradation * cAMP_13C00110D10;
    Vf_AMP_Degradation_13C00111D10_13C00111D10 = Kf_AMP_Degradation * cAMP_13C00111D10;
    Vf_AMP_Degradation_13C01000D10_13C01000D10 = Kf_AMP_Degradation * cAMP_13C01000D10;
    Vf_AMP_Degradation_13C01001D10_13C01001D10 = Kf_AMP_Degradation * cAMP_13C01001D10;
    Vf_AMP_Degradation_13C01010D10_13C01010D10 = Kf_AMP_Degradation * cAMP_13C01010D10;
    Vf_AMP_Degradation_13C01011D10_13C01011D10 = Kf_AMP_Degradation * cAMP_13C01011D10;
    Vf_AMP_Degradation_13C01100D10_13C01100D10 = Kf_AMP_Degradation * cAMP_13C01100D10;
    Vf_AMP_Degradation_13C01101D10_13C01101D10 = Kf_AMP_Degradation * cAMP_13C01101D10;
    Vf_AMP_Degradation_13C01110D10_13C01110D10 = Kf_AMP_Degradation * cAMP_13C01110D10;
    Vf_AMP_Degradation_13C01111D10_13C01111D10 = Kf_AMP_Degradation * cAMP_13C01111D10;
    Vf_AMP_Degradation_13C10000D10_13C10000D10 = Kf_AMP_Degradation * cAMP_13C10000D10;
    Vf_AMP_Degradation_13C10001D10_13C10001D10 = Kf_AMP_Degradation * cAMP_13C10001D10;
    Vf_AMP_Degradation_13C10010D10_13C10010D10 = Kf_AMP_Degradation * cAMP_13C10010D10;
    Vf_AMP_Degradation_13C10011D10_13C10011D10 = Kf_AMP_Degradation * cAMP_13C10011D10;
    Vf_AMP_Degradation_13C10100D10_13C10100D10 = Kf_AMP_Degradation * cAMP_13C10100D10;
    Vf_AMP_Degradation_13C10101D10_13C10101D10 = Kf_AMP_Degradation * cAMP_13C10101D10;
    Vf_AMP_Degradation_13C10110D10_13C10110D10 = Kf_AMP_Degradation * cAMP_13C10110D10;
    Vf_AMP_Degradation_13C10111D10_13C10111D10 = Kf_AMP_Degradation * cAMP_13C10111D10;
    Vf_AMP_Degradation_13C11000D10_13C11000D10 = Kf_AMP_Degradation * cAMP_13C11000D10;
    Vf_AMP_Degradation_13C11001D10_13C11001D10 = Kf_AMP_Degradation * cAMP_13C11001D10;
    Vf_AMP_Degradation_13C11010D10_13C11010D10 = Kf_AMP_Degradation * cAMP_13C11010D10;
    Vf_AMP_Degradation_13C11011D10_13C11011D10 = Kf_AMP_Degradation * cAMP_13C11011D10;
    Vf_AMP_Degradation_13C11100D10_13C11100D10 = Kf_AMP_Degradation * cAMP_13C11100D10;
    Vf_AMP_Degradation_13C11101D10_13C11101D10 = Kf_AMP_Degradation * cAMP_13C11101D10;
    Vf_AMP_Degradation_13C11110D10_13C11110D10 = Kf_AMP_Degradation * cAMP_13C11110D10;
    Vf_AMP_Degradation_13C11111D10_13C11111D10 = Kf_AMP_Degradation * cAMP_13C11111D10;
    Vf_AMP_Degradation_13C00000D11_13C00000D11 = Kf_AMP_Degradation * cAMP_13C00000D11;
    Vf_AMP_Degradation_13C00001D11_13C00001D11 = Kf_AMP_Degradation * cAMP_13C00001D11;
    Vf_AMP_Degradation_13C00010D11_13C00010D11 = Kf_AMP_Degradation * cAMP_13C00010D11;
    Vf_AMP_Degradation_13C00011D11_13C00011D11 = Kf_AMP_Degradation * cAMP_13C00011D11;
    Vf_AMP_Degradation_13C00100D11_13C00100D11 = Kf_AMP_Degradation * cAMP_13C00100D11;
    Vf_AMP_Degradation_13C00101D11_13C00101D11 = Kf_AMP_Degradation * cAMP_13C00101D11;
    Vf_AMP_Degradation_13C00110D11_13C00110D11 = Kf_AMP_Degradation * cAMP_13C00110D11;
    Vf_AMP_Degradation_13C00111D11_13C00111D11 = Kf_AMP_Degradation * cAMP_13C00111D11;
    Vf_AMP_Degradation_13C01000D11_13C01000D11 = Kf_AMP_Degradation * cAMP_13C01000D11;
    Vf_AMP_Degradation_13C01001D11_13C01001D11 = Kf_AMP_Degradation * cAMP_13C01001D11;
    Vf_AMP_Degradation_13C01010D11_13C01010D11 = Kf_AMP_Degradation * cAMP_13C01010D11;
    Vf_AMP_Degradation_13C01011D11_13C01011D11 = Kf_AMP_Degradation * cAMP_13C01011D11;
    Vf_AMP_Degradation_13C01100D11_13C01100D11 = Kf_AMP_Degradation * cAMP_13C01100D11;
    Vf_AMP_Degradation_13C01101D11_13C01101D11 = Kf_AMP_Degradation * cAMP_13C01101D11;
    Vf_AMP_Degradation_13C01110D11_13C01110D11 = Kf_AMP_Degradation * cAMP_13C01110D11;
    Vf_AMP_Degradation_13C01111D11_13C01111D11 = Kf_AMP_Degradation * cAMP_13C01111D11;
    Vf_AMP_Degradation_13C10000D11_13C10000D11 = Kf_AMP_Degradation * cAMP_13C10000D11;
    Vf_AMP_Degradation_13C10001D11_13C10001D11 = Kf_AMP_Degradation * cAMP_13C10001D11;
    Vf_AMP_Degradation_13C10010D11_13C10010D11 = Kf_AMP_Degradation * cAMP_13C10010D11;
    Vf_AMP_Degradation_13C10011D11_13C10011D11 = Kf_AMP_Degradation * cAMP_13C10011D11;
    Vf_AMP_Degradation_13C10100D11_13C10100D11 = Kf_AMP_Degradation * cAMP_13C10100D11;
    Vf_AMP_Degradation_13C10101D11_13C10101D11 = Kf_AMP_Degradation * cAMP_13C10101D11;
    Vf_AMP_Degradation_13C10110D11_13C10110D11 = Kf_AMP_Degradation * cAMP_13C10110D11;
    Vf_AMP_Degradation_13C10111D11_13C10111D11 = Kf_AMP_Degradation * cAMP_13C10111D11;
    Vf_AMP_Degradation_13C11000D11_13C11000D11 = Kf_AMP_Degradation * cAMP_13C11000D11;
    Vf_AMP_Degradation_13C11001D11_13C11001D11 = Kf_AMP_Degradation * cAMP_13C11001D11;
    Vf_AMP_Degradation_13C11010D11_13C11010D11 = Kf_AMP_Degradation * cAMP_13C11010D11;
    Vf_AMP_Degradation_13C11011D11_13C11011D11 = Kf_AMP_Degradation * cAMP_13C11011D11;
    Vf_AMP_Degradation_13C11100D11_13C11100D11 = Kf_AMP_Degradation * cAMP_13C11100D11;
    Vf_AMP_Degradation_13C11101D11_13C11101D11 = Kf_AMP_Degradation * cAMP_13C11101D11;
    Vf_AMP_Degradation_13C11110D11_13C11110D11 = Kf_AMP_Degradation * cAMP_13C11110D11;
    Vf_AMP_Degradation_13C11111D11_13C11111D11 = Kf_AMP_Degradation * cAMP_13C11111D11;









    % eGLC	
    %  dxdt(1, 1) = (eGLC) = 
    % - Vf_Vin_GLC
    % + Vr_Vin_GLC
    
    dxdt(1, 1) =  - Vf_Vin_GLC_13C0_13C0	 +  Vr_Vin_GLC_13C0_13C0;					%  eGLC_13C0 
    dxdt(2, 1) =  - Vf_Vin_GLC_13C1_13C1	 +  Vr_Vin_GLC_13C1_13C1;					%  eGLC_13C1 

    
    
    % cGLC
    %  dxdt(1, 1) = (cGLC) = 
    % - Vf_GLC_PYR
    % - Vf_GLC_PRPP
    % - Vf_GLC_SerPool1
    % - Vr_Vin_GLC 		
    % + Vf_Vin_GLC

		dxdt(3, 1) =  - Vf_GLC_PYR_13C0  -  Vf_GLC_PRPP_13C0_13C0  -  Vf_GLC_SerPool1_13C0_13C0D000  -  Vr_Vin_GLC_13C0_13C0  +  Vf_Vin_GLC_13C0_13C0;							% cGLC_13C0 
		dxdt(4, 1) =  - Vf_GLC_PYR_13C1  -  Vf_GLC_PRPP_13C1_13C1  -  Vf_GLC_SerPool1_13C1_13C3D000  -  Vr_Vin_GLC_13C1_13C1  +  Vf_Vin_GLC_13C1_13C1;							% cGLC_13C1 
    


    % eSer
    %  dxdt(1, 1) = (eSer) = 
    % - Vf_Vin_Ser 	
    % + Vr_Vin_Ser	
    % - Vf_Vin_Ser3	
    % + Vr_Vin_Ser3	

		dxdt(5, 1)  =  -	Vf_Vin_Ser_13C000D000_13C000D000  +  Vr_Vin_Ser_13C000D000_13C000D000	 -  Vf_Vin_Ser3_13C000D000_13C000D000  +  Vr_Vin_Ser3_13C000D000_13C000D000;														%  eSer_13C000D000
		dxdt(6, 1)  =  -	Vf_Vin_Ser_13C001D000_13C001D000  +  Vr_Vin_Ser_13C001D000_13C001D000	 -  Vf_Vin_Ser3_13C001D000_13C001D000  +  Vr_Vin_Ser3_13C001D000_13C001D000;                            %  eSer_13C001D000
		dxdt(7, 1)  =  -	Vf_Vin_Ser_13C010D000_13C010D000  +  Vr_Vin_Ser_13C010D000_13C010D000	 -  Vf_Vin_Ser3_13C010D000_13C010D000  +  Vr_Vin_Ser3_13C010D000_13C010D000;                            %  eSer_13C010D000
		dxdt(8, 1)  =  -	Vf_Vin_Ser_13C011D000_13C011D000  +  Vr_Vin_Ser_13C011D000_13C011D000	 -  Vf_Vin_Ser3_13C011D000_13C011D000  +  Vr_Vin_Ser3_13C011D000_13C011D000;                            %  eSer_13C011D000
		dxdt(9, 1)  =  -	Vf_Vin_Ser_13C100D000_13C100D000  +  Vr_Vin_Ser_13C100D000_13C100D000	 -  Vf_Vin_Ser3_13C100D000_13C100D000  +  Vr_Vin_Ser3_13C100D000_13C100D000;                            %  eSer_13C100D000
		dxdt(10, 1) =  -	Vf_Vin_Ser_13C101D000_13C101D000  +  Vr_Vin_Ser_13C101D000_13C101D000	 -  Vf_Vin_Ser3_13C101D000_13C101D000  +  Vr_Vin_Ser3_13C101D000_13C101D000;														%  eSer_13C101D000 
		dxdt(11, 1) =  -	Vf_Vin_Ser_13C110D000_13C110D000  +  Vr_Vin_Ser_13C110D000_13C110D000	 -  Vf_Vin_Ser3_13C110D000_13C110D000  +  Vr_Vin_Ser3_13C110D000_13C110D000;                            %  eSer_13C110D000 
		dxdt(12, 1) =  -	Vf_Vin_Ser_13C111D000_13C111D000  +  Vr_Vin_Ser_13C111D000_13C111D000	 -  Vf_Vin_Ser3_13C111D000_13C111D000  +  Vr_Vin_Ser3_13C111D000_13C111D000;                            %  eSer_13C111D000 
		dxdt(13, 1) =  -	Vf_Vin_Ser_13C000D001_13C000D001  +  Vr_Vin_Ser_13C000D001_13C000D001	 -  Vf_Vin_Ser3_13C000D001_13C000D001  +  Vr_Vin_Ser3_13C000D001_13C000D001;                            %  eSer_13C000D001 
		dxdt(14, 1) =  -	Vf_Vin_Ser_13C001D001_13C001D001  +  Vr_Vin_Ser_13C001D001_13C001D001	 -  Vf_Vin_Ser3_13C001D001_13C001D001  +  Vr_Vin_Ser3_13C001D001_13C001D001;                            %  eSer_13C001D001 
		dxdt(15, 1) =  -	Vf_Vin_Ser_13C010D001_13C010D001  +  Vr_Vin_Ser_13C010D001_13C010D001	 -  Vf_Vin_Ser3_13C010D001_13C010D001  +  Vr_Vin_Ser3_13C010D001_13C010D001;                            %  eSer_13C010D001 
		dxdt(16, 1) =  -	Vf_Vin_Ser_13C011D001_13C011D001  +  Vr_Vin_Ser_13C011D001_13C011D001	 -  Vf_Vin_Ser3_13C011D001_13C011D001  +  Vr_Vin_Ser3_13C011D001_13C011D001;                            %  eSer_13C011D001 
		dxdt(17, 1) =  -	Vf_Vin_Ser_13C100D001_13C100D001  +  Vr_Vin_Ser_13C100D001_13C100D001	 -  Vf_Vin_Ser3_13C100D001_13C100D001  +  Vr_Vin_Ser3_13C100D001_13C100D001;                            %  eSer_13C100D001 
		dxdt(18, 1) =  -	Vf_Vin_Ser_13C101D001_13C101D001  +  Vr_Vin_Ser_13C101D001_13C101D001	 -  Vf_Vin_Ser3_13C101D001_13C101D001  +  Vr_Vin_Ser3_13C101D001_13C101D001;                            %  eSer_13C101D001 
		dxdt(19, 1) =  -	Vf_Vin_Ser_13C110D001_13C110D001  +  Vr_Vin_Ser_13C110D001_13C110D001	 -  Vf_Vin_Ser3_13C110D001_13C110D001  +  Vr_Vin_Ser3_13C110D001_13C110D001;                            %  eSer_13C110D001 
		dxdt(20, 1) =  -	Vf_Vin_Ser_13C111D001_13C111D001  +  Vr_Vin_Ser_13C111D001_13C111D001	 -  Vf_Vin_Ser3_13C111D001_13C111D001  +  Vr_Vin_Ser3_13C111D001_13C111D001;                            %  eSer_13C111D001 
		dxdt(21, 1) =  -	Vf_Vin_Ser_13C000D010_13C000D010  +  Vr_Vin_Ser_13C000D010_13C000D010	 -  Vf_Vin_Ser3_13C000D010_13C000D010  +  Vr_Vin_Ser3_13C000D010_13C000D010;                            %  eSer_13C000D010 
		dxdt(22, 1) =  -	Vf_Vin_Ser_13C001D010_13C001D010  +  Vr_Vin_Ser_13C001D010_13C001D010	 -  Vf_Vin_Ser3_13C001D010_13C001D010  +  Vr_Vin_Ser3_13C001D010_13C001D010;                            %  eSer_13C001D010 
		dxdt(23, 1) =  -	Vf_Vin_Ser_13C010D010_13C010D010  +  Vr_Vin_Ser_13C010D010_13C010D010	 -  Vf_Vin_Ser3_13C010D010_13C010D010  +  Vr_Vin_Ser3_13C010D010_13C010D010;                            %  eSer_13C010D010 
		dxdt(24, 1) =  -	Vf_Vin_Ser_13C011D010_13C011D010  +  Vr_Vin_Ser_13C011D010_13C011D010	 -  Vf_Vin_Ser3_13C011D010_13C011D010  +  Vr_Vin_Ser3_13C011D010_13C011D010;                            %  eSer_13C011D010 
		dxdt(25, 1) =  -	Vf_Vin_Ser_13C100D010_13C100D010  +  Vr_Vin_Ser_13C100D010_13C100D010	 -  Vf_Vin_Ser3_13C100D010_13C100D010  +  Vr_Vin_Ser3_13C100D010_13C100D010;                            %  eSer_13C100D010 
		dxdt(26, 1) =  -	Vf_Vin_Ser_13C101D010_13C101D010  +  Vr_Vin_Ser_13C101D010_13C101D010	 -  Vf_Vin_Ser3_13C101D010_13C101D010  +  Vr_Vin_Ser3_13C101D010_13C101D010;                            %  eSer_13C101D010 
		dxdt(27, 1) =  -	Vf_Vin_Ser_13C110D010_13C110D010  +  Vr_Vin_Ser_13C110D010_13C110D010	 -  Vf_Vin_Ser3_13C110D010_13C110D010  +  Vr_Vin_Ser3_13C110D010_13C110D010;                            %  eSer_13C110D010 
		dxdt(28, 1) =  -	Vf_Vin_Ser_13C111D010_13C111D010  +  Vr_Vin_Ser_13C111D010_13C111D010	 -  Vf_Vin_Ser3_13C111D010_13C111D010  +  Vr_Vin_Ser3_13C111D010_13C111D010;                            %  eSer_13C111D010 
		dxdt(29, 1) =  -	Vf_Vin_Ser_13C000D011_13C000D011  +  Vr_Vin_Ser_13C000D011_13C000D011	 -  Vf_Vin_Ser3_13C000D011_13C000D011  +  Vr_Vin_Ser3_13C000D011_13C000D011;                            %  eSer_13C000D011 
		dxdt(30, 1) =  -	Vf_Vin_Ser_13C001D011_13C001D011  +  Vr_Vin_Ser_13C001D011_13C001D011	 -  Vf_Vin_Ser3_13C001D011_13C001D011  +  Vr_Vin_Ser3_13C001D011_13C001D011;                            %  eSer_13C001D011 
		dxdt(31, 1) =  -	Vf_Vin_Ser_13C010D011_13C010D011  +  Vr_Vin_Ser_13C010D011_13C010D011	 -  Vf_Vin_Ser3_13C010D011_13C010D011  +  Vr_Vin_Ser3_13C010D011_13C010D011;                            %  eSer_13C010D011 
		dxdt(32, 1) =  -	Vf_Vin_Ser_13C011D011_13C011D011  +  Vr_Vin_Ser_13C011D011_13C011D011	 -  Vf_Vin_Ser3_13C011D011_13C011D011  +  Vr_Vin_Ser3_13C011D011_13C011D011;                            %  eSer_13C011D011 
		dxdt(33, 1) =  -	Vf_Vin_Ser_13C100D011_13C100D011  +  Vr_Vin_Ser_13C100D011_13C100D011	 -  Vf_Vin_Ser3_13C100D011_13C100D011  +  Vr_Vin_Ser3_13C100D011_13C100D011;                            %  eSer_13C100D011 
		dxdt(34, 1) =  -	Vf_Vin_Ser_13C101D011_13C101D011  +  Vr_Vin_Ser_13C101D011_13C101D011	 -  Vf_Vin_Ser3_13C101D011_13C101D011  +  Vr_Vin_Ser3_13C101D011_13C101D011;                            %  eSer_13C101D011 
		dxdt(35, 1) =  -	Vf_Vin_Ser_13C110D011_13C110D011  +  Vr_Vin_Ser_13C110D011_13C110D011	 -  Vf_Vin_Ser3_13C110D011_13C110D011  +  Vr_Vin_Ser3_13C110D011_13C110D011;                            %  eSer_13C110D011 
		dxdt(36, 1) =  -	Vf_Vin_Ser_13C111D011_13C111D011  +  Vr_Vin_Ser_13C111D011_13C111D011	 -  Vf_Vin_Ser3_13C111D011_13C111D011  +  Vr_Vin_Ser3_13C111D011_13C111D011;                            %  eSer_13C111D011 
		dxdt(37, 1) =  -	Vf_Vin_Ser_13C000D100_13C000D100  +  Vr_Vin_Ser_13C000D100_13C000D100	 -  Vf_Vin_Ser3_13C000D100_13C000D100  +  Vr_Vin_Ser3_13C000D100_13C000D100;                            %  eSer_13C000D100 
		dxdt(38, 1) =  -	Vf_Vin_Ser_13C001D100_13C001D100  +  Vr_Vin_Ser_13C001D100_13C001D100	 -  Vf_Vin_Ser3_13C001D100_13C001D100  +  Vr_Vin_Ser3_13C001D100_13C001D100;                            %  eSer_13C001D100 
		dxdt(39, 1) =  -	Vf_Vin_Ser_13C010D100_13C010D100  +  Vr_Vin_Ser_13C010D100_13C010D100	 -  Vf_Vin_Ser3_13C010D100_13C010D100  +  Vr_Vin_Ser3_13C010D100_13C010D100;                            %  eSer_13C010D100 
		dxdt(40, 1) =  -	Vf_Vin_Ser_13C011D100_13C011D100  +  Vr_Vin_Ser_13C011D100_13C011D100	 -  Vf_Vin_Ser3_13C011D100_13C011D100  +  Vr_Vin_Ser3_13C011D100_13C011D100;                            %  eSer_13C011D100 
		dxdt(41, 1) =  -	Vf_Vin_Ser_13C100D100_13C100D100  +  Vr_Vin_Ser_13C100D100_13C100D100	 -  Vf_Vin_Ser3_13C100D100_13C100D100  +  Vr_Vin_Ser3_13C100D100_13C100D100;                            %  eSer_13C100D100 
		dxdt(42, 1) =  -	Vf_Vin_Ser_13C101D100_13C101D100  +  Vr_Vin_Ser_13C101D100_13C101D100	 -  Vf_Vin_Ser3_13C101D100_13C101D100  +  Vr_Vin_Ser3_13C101D100_13C101D100;                            %  eSer_13C101D100 
		dxdt(43, 1) =  -	Vf_Vin_Ser_13C110D100_13C110D100  +  Vr_Vin_Ser_13C110D100_13C110D100	 -  Vf_Vin_Ser3_13C110D100_13C110D100  +  Vr_Vin_Ser3_13C110D100_13C110D100;                            %  eSer_13C110D100 
		dxdt(44, 1) =  -	Vf_Vin_Ser_13C111D100_13C111D100  +  Vr_Vin_Ser_13C111D100_13C111D100	 -  Vf_Vin_Ser3_13C111D100_13C111D100  +  Vr_Vin_Ser3_13C111D100_13C111D100;                            %  eSer_13C111D100 
		dxdt(45, 1) =  -	Vf_Vin_Ser_13C000D101_13C000D101  +  Vr_Vin_Ser_13C000D101_13C000D101	 -  Vf_Vin_Ser3_13C000D101_13C000D101  +  Vr_Vin_Ser3_13C000D101_13C000D101;                            %  eSer_13C000D101 
		dxdt(46, 1) =  -	Vf_Vin_Ser_13C001D101_13C001D101  +  Vr_Vin_Ser_13C001D101_13C001D101	 -  Vf_Vin_Ser3_13C001D101_13C001D101  +  Vr_Vin_Ser3_13C001D101_13C001D101;                            %  eSer_13C001D101 
		dxdt(47, 1) =  -	Vf_Vin_Ser_13C010D101_13C010D101  +  Vr_Vin_Ser_13C010D101_13C010D101	 -  Vf_Vin_Ser3_13C010D101_13C010D101  +  Vr_Vin_Ser3_13C010D101_13C010D101;                            %  eSer_13C010D101 
		dxdt(48, 1) =  -	Vf_Vin_Ser_13C011D101_13C011D101  +  Vr_Vin_Ser_13C011D101_13C011D101	 -  Vf_Vin_Ser3_13C011D101_13C011D101  +  Vr_Vin_Ser3_13C011D101_13C011D101;                            %  eSer_13C011D101 
		dxdt(49, 1) =  -	Vf_Vin_Ser_13C100D101_13C100D101  +  Vr_Vin_Ser_13C100D101_13C100D101	 -  Vf_Vin_Ser3_13C100D101_13C100D101  +  Vr_Vin_Ser3_13C100D101_13C100D101;                            %  eSer_13C100D101 
		dxdt(50, 1) =  -	Vf_Vin_Ser_13C101D101_13C101D101  +  Vr_Vin_Ser_13C101D101_13C101D101	 -  Vf_Vin_Ser3_13C101D101_13C101D101  +  Vr_Vin_Ser3_13C101D101_13C101D101;                            %  eSer_13C101D101 
		dxdt(51, 1) =  -	Vf_Vin_Ser_13C110D101_13C110D101  +  Vr_Vin_Ser_13C110D101_13C110D101	 -  Vf_Vin_Ser3_13C110D101_13C110D101  +  Vr_Vin_Ser3_13C110D101_13C110D101;                            %  eSer_13C110D101 
		dxdt(52, 1) =  -	Vf_Vin_Ser_13C111D101_13C111D101  +  Vr_Vin_Ser_13C111D101_13C111D101	 -  Vf_Vin_Ser3_13C111D101_13C111D101  +  Vr_Vin_Ser3_13C111D101_13C111D101;                            %  eSer_13C111D101 
		dxdt(53, 1) =  -	Vf_Vin_Ser_13C000D110_13C000D110  +  Vr_Vin_Ser_13C000D110_13C000D110	 -  Vf_Vin_Ser3_13C000D110_13C000D110  +  Vr_Vin_Ser3_13C000D110_13C000D110;                            %  eSer_13C000D110 
		dxdt(54, 1) =  -	Vf_Vin_Ser_13C001D110_13C001D110  +  Vr_Vin_Ser_13C001D110_13C001D110	 -  Vf_Vin_Ser3_13C001D110_13C001D110  +  Vr_Vin_Ser3_13C001D110_13C001D110;                            %  eSer_13C001D110 
		dxdt(55, 1) =  -	Vf_Vin_Ser_13C010D110_13C010D110  +  Vr_Vin_Ser_13C010D110_13C010D110	 -  Vf_Vin_Ser3_13C010D110_13C010D110  +  Vr_Vin_Ser3_13C010D110_13C010D110;                            %  eSer_13C010D110 
		dxdt(56, 1) =  -	Vf_Vin_Ser_13C011D110_13C011D110  +  Vr_Vin_Ser_13C011D110_13C011D110	 -  Vf_Vin_Ser3_13C011D110_13C011D110  +  Vr_Vin_Ser3_13C011D110_13C011D110;                            %  eSer_13C011D110 
		dxdt(57, 1) =  -	Vf_Vin_Ser_13C100D110_13C100D110  +  Vr_Vin_Ser_13C100D110_13C100D110	 -  Vf_Vin_Ser3_13C100D110_13C100D110  +  Vr_Vin_Ser3_13C100D110_13C100D110;                            %  eSer_13C100D110 
		dxdt(58, 1) =  -	Vf_Vin_Ser_13C101D110_13C101D110  +  Vr_Vin_Ser_13C101D110_13C101D110	 -  Vf_Vin_Ser3_13C101D110_13C101D110  +  Vr_Vin_Ser3_13C101D110_13C101D110;                            %  eSer_13C101D110 
		dxdt(59, 1) =  -	Vf_Vin_Ser_13C110D110_13C110D110  +  Vr_Vin_Ser_13C110D110_13C110D110	 -  Vf_Vin_Ser3_13C110D110_13C110D110  +  Vr_Vin_Ser3_13C110D110_13C110D110;                            %  eSer_13C110D110 
		dxdt(60, 1) =  -	Vf_Vin_Ser_13C111D110_13C111D110  +  Vr_Vin_Ser_13C111D110_13C111D110	 -  Vf_Vin_Ser3_13C111D110_13C111D110  +  Vr_Vin_Ser3_13C111D110_13C111D110;                            %  eSer_13C111D110 
		dxdt(61, 1) =  -	Vf_Vin_Ser_13C000D111_13C000D111  +  Vr_Vin_Ser_13C000D111_13C000D111	 -  Vf_Vin_Ser3_13C000D111_13C000D111  +  Vr_Vin_Ser3_13C000D111_13C000D111;                            %  eSer_13C000D111 
		dxdt(62, 1) =  -	Vf_Vin_Ser_13C001D111_13C001D111  +  Vr_Vin_Ser_13C001D111_13C001D111	 -  Vf_Vin_Ser3_13C001D111_13C001D111  +  Vr_Vin_Ser3_13C001D111_13C001D111;                            %  eSer_13C001D111 
		dxdt(63, 1) =  -	Vf_Vin_Ser_13C010D111_13C010D111  +  Vr_Vin_Ser_13C010D111_13C010D111	 -  Vf_Vin_Ser3_13C010D111_13C010D111  +  Vr_Vin_Ser3_13C010D111_13C010D111;                            %  eSer_13C010D111 
		dxdt(64, 1) =  -	Vf_Vin_Ser_13C011D111_13C011D111  +  Vr_Vin_Ser_13C011D111_13C011D111	 -  Vf_Vin_Ser3_13C011D111_13C011D111  +  Vr_Vin_Ser3_13C011D111_13C011D111;                            %  eSer_13C011D111 
		dxdt(65, 1) =  -	Vf_Vin_Ser_13C100D111_13C100D111  +  Vr_Vin_Ser_13C100D111_13C100D111	 -  Vf_Vin_Ser3_13C100D111_13C100D111  +  Vr_Vin_Ser3_13C100D111_13C100D111;                            %  eSer_13C100D111 
		dxdt(66, 1) =  -	Vf_Vin_Ser_13C101D111_13C101D111  +  Vr_Vin_Ser_13C101D111_13C101D111	 -  Vf_Vin_Ser3_13C101D111_13C101D111  +  Vr_Vin_Ser3_13C101D111_13C101D111;                            %  eSer_13C101D111 
		dxdt(67, 1) =  -	Vf_Vin_Ser_13C110D111_13C110D111  +  Vr_Vin_Ser_13C110D111_13C110D111	 -  Vf_Vin_Ser3_13C110D111_13C110D111  +  Vr_Vin_Ser3_13C110D111_13C110D111;                            %  eSer_13C110D111 
		dxdt(68, 1) =  -	Vf_Vin_Ser_13C111D111_13C111D111  +  Vr_Vin_Ser_13C111D111_13C111D111	 -  Vf_Vin_Ser3_13C111D111_13C111D111  +  Vr_Vin_Ser3_13C111D111_13C111D111;                            %  eSer_13C111D111 
          


    % eGly
    %  dxdt(1, 1) = (eGly) = 
    % - Vf_Vin_Gly 
    % + Vr_Vin_Gly	
    % - Vf_Vin_Gly3 
    % + Vr_Vin_Gly3	

		dxdt(69, 1)  =  - Vf_Vin_Gly_13C00D00_13C00D00  +  Vr_Vin_Gly_13C00D00_13C00D00  -  Vf_Vin_Gly3_13C00D00_13C00D00  +  Vr_Vin_Gly3_13C00D00_13C00D00;							%  eGly_13C00D00
		dxdt(70, 1)  =  - Vf_Vin_Gly_13C01D00_13C01D00  +  Vr_Vin_Gly_13C01D00_13C01D00  -  Vf_Vin_Gly3_13C01D00_13C01D00  +  Vr_Vin_Gly3_13C01D00_13C01D00;              %  eGly_13C01D00
		dxdt(71, 1)  =  - Vf_Vin_Gly_13C10D00_13C10D00  +  Vr_Vin_Gly_13C10D00_13C10D00  -  Vf_Vin_Gly3_13C10D00_13C10D00  +  Vr_Vin_Gly3_13C10D00_13C10D00;              %  eGly_13C10D00
		dxdt(72, 1)  =  - Vf_Vin_Gly_13C11D00_13C11D00  +  Vr_Vin_Gly_13C11D00_13C11D00  -  Vf_Vin_Gly3_13C11D00_13C11D00  +  Vr_Vin_Gly3_13C11D00_13C11D00;              %  eGly_13C11D00
		dxdt(73, 1)  =  - Vf_Vin_Gly_13C00D01_13C00D01  +  Vr_Vin_Gly_13C00D01_13C00D01  -  Vf_Vin_Gly3_13C00D01_13C00D01  +  Vr_Vin_Gly3_13C00D01_13C00D01;              %  eGly_13C00D01
		dxdt(74, 1)  =  - Vf_Vin_Gly_13C01D01_13C01D01  +  Vr_Vin_Gly_13C01D01_13C01D01  -  Vf_Vin_Gly3_13C01D01_13C01D01  +  Vr_Vin_Gly3_13C01D01_13C01D01;              %  eGly_13C01D01
		dxdt(75, 1)  =  - Vf_Vin_Gly_13C10D01_13C10D01  +  Vr_Vin_Gly_13C10D01_13C10D01  -  Vf_Vin_Gly3_13C10D01_13C10D01  +  Vr_Vin_Gly3_13C10D01_13C10D01;              %  eGly_13C10D01
		dxdt(76, 1)  =  - Vf_Vin_Gly_13C11D01_13C11D01  +  Vr_Vin_Gly_13C11D01_13C11D01  -  Vf_Vin_Gly3_13C11D01_13C11D01  +  Vr_Vin_Gly3_13C11D01_13C11D01;              %  eGly_13C11D01
		dxdt(77, 1)  =  - Vf_Vin_Gly_13C00D10_13C00D10  +  Vr_Vin_Gly_13C00D10_13C00D10  -  Vf_Vin_Gly3_13C00D10_13C00D10  +  Vr_Vin_Gly3_13C00D10_13C00D10;              %  eGly_13C00D10
		dxdt(78, 1)  =  - Vf_Vin_Gly_13C01D10_13C01D10  +  Vr_Vin_Gly_13C01D10_13C01D10  -  Vf_Vin_Gly3_13C01D10_13C01D10  +  Vr_Vin_Gly3_13C01D10_13C01D10;              %  eGly_13C01D10
		dxdt(79, 1)  =  - Vf_Vin_Gly_13C10D10_13C10D10  +  Vr_Vin_Gly_13C10D10_13C10D10  -  Vf_Vin_Gly3_13C10D10_13C10D10  +  Vr_Vin_Gly3_13C10D10_13C10D10;              %  eGly_13C10D10
		dxdt(80, 1)  =  - Vf_Vin_Gly_13C11D10_13C11D10  +  Vr_Vin_Gly_13C11D10_13C11D10  -  Vf_Vin_Gly3_13C11D10_13C11D10  +  Vr_Vin_Gly3_13C11D10_13C11D10;              %  eGly_13C11D10
		dxdt(81, 1)  =  - Vf_Vin_Gly_13C00D11_13C00D11  +  Vr_Vin_Gly_13C00D11_13C00D11  -  Vf_Vin_Gly3_13C00D11_13C00D11  +  Vr_Vin_Gly3_13C00D11_13C00D11;              %  eGly_13C00D11
		dxdt(82, 1)  =  - Vf_Vin_Gly_13C01D11_13C01D11  +  Vr_Vin_Gly_13C01D11_13C01D11  -  Vf_Vin_Gly3_13C01D11_13C01D11  +  Vr_Vin_Gly3_13C01D11_13C01D11;              %  eGly_13C01D11
		dxdt(83, 1)  =  - Vf_Vin_Gly_13C10D11_13C10D11  +  Vr_Vin_Gly_13C10D11_13C10D11  -  Vf_Vin_Gly3_13C10D11_13C10D11  +  Vr_Vin_Gly3_13C10D11_13C10D11;              %  eGly_13C10D11
		dxdt(84, 1)  =  - Vf_Vin_Gly_13C11D11_13C11D11  +  Vr_Vin_Gly_13C11D11_13C11D11  -  Vf_Vin_Gly3_13C11D11_13C11D11  +  Vr_Vin_Gly3_13C11D11_13C11D11;              %  eGly_13C11D11



    % cSerPool1
    %  dxdt(1, 1) = (cSerPool1) = 
    % - Vf_SerPool1_GlyPool1 
    % + Vr_SerPool1_GlyPool1 
    % + Vf_GLC_SerPool1 

		dxdt(85, 1) =  - Vf_SerPool1_GlyPool1_13C000D000_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C0D00__13C000D000  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D00__13C000D000  +  Vf_GLC_SerPool1_13C0_13C0D000;								%  cSerPool1_13C000D000 
		dxdt(86, 1) =  - Vf_SerPool1_GlyPool1_13C001D000_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D00__13C001D000  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D00__13C001D000;                                                  %  cSerPool1_13C001D000 
		dxdt(87, 1) =  - Vf_SerPool1_GlyPool1_13C010D000_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D00__13C010D000  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D00__13C010D000;                                                  %  cSerPool1_13C010D000 
		dxdt(88, 1) =  - Vf_SerPool1_GlyPool1_13C011D000_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D00__13C011D000  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D00__13C011D000;                                                  %  cSerPool1_13C011D000 
		dxdt(89, 1) =  - Vf_SerPool1_GlyPool1_13C100D000_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D00__13C100D000  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D00__13C100D000;                                                  %  cSerPool1_13C100D000 
		dxdt(90, 1) =  - Vf_SerPool1_GlyPool1_13C101D000_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D00__13C101D000  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D00__13C101D000;                                                  %  cSerPool1_13C101D000 
		dxdt(91, 1) =  - Vf_SerPool1_GlyPool1_13C110D000_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D00__13C110D000  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D00__13C110D000;                                                  %  cSerPool1_13C110D000 
		dxdt(92, 1) =  - Vf_SerPool1_GlyPool1_13C111D000_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D00__13C111D000  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D00__13C111D000  +  Vf_GLC_SerPool1_13C1_13C3D000;                %  cSerPool1_13C111D000 
		dxdt(93, 1) =  - Vf_SerPool1_GlyPool1_13C000D001_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C0D01__13C000D001  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D01__13C000D001;                                                  %  cSerPool1_13C000D001 
		dxdt(94, 1) =  - Vf_SerPool1_GlyPool1_13C001D001_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D01__13C001D001  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D01__13C001D001;                                                  %  cSerPool1_13C001D001 
		dxdt(95, 1) =  - Vf_SerPool1_GlyPool1_13C010D001_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D01__13C010D001  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D01__13C010D001;                                                  %  cSerPool1_13C010D001 
		dxdt(96, 1) =  - Vf_SerPool1_GlyPool1_13C011D001_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D01__13C011D001  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D01__13C011D001;                                                  %  cSerPool1_13C011D001 
		dxdt(97, 1) =  - Vf_SerPool1_GlyPool1_13C100D001_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D01__13C100D001  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D01__13C100D001;                                                  %  cSerPool1_13C100D001 
		dxdt(98, 1) =  - Vf_SerPool1_GlyPool1_13C101D001_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D01__13C101D001  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D01__13C101D001;                                                  %  cSerPool1_13C101D001 
		dxdt(99, 1) =  - Vf_SerPool1_GlyPool1_13C110D001_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D01__13C110D001  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D01__13C110D001;                                                  %  cSerPool1_13C110D001 
		dxdt(100, 1) =  - Vf_SerPool1_GlyPool1_13C111D001_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D01__13C111D001  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D01__13C111D001;																									%  cSerPool1_13C111D001 
		dxdt(101, 1) =  - Vf_SerPool1_GlyPool1_13C000D010_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C0D10__13C000D010  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D10__13C000D010;                                                 %  cSerPool1_13C000D010 
		dxdt(102, 1) =  - Vf_SerPool1_GlyPool1_13C001D010_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D10__13C001D010  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D10__13C001D010;                                                 %  cSerPool1_13C001D010 
		dxdt(103, 1) =  - Vf_SerPool1_GlyPool1_13C010D010_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D10__13C010D010  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D10__13C010D010;                                                 %  cSerPool1_13C010D010 
		dxdt(104, 1) =  - Vf_SerPool1_GlyPool1_13C011D010_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D10__13C011D010  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D10__13C011D010;                                                 %  cSerPool1_13C011D010 
		dxdt(105, 1) =  - Vf_SerPool1_GlyPool1_13C100D010_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D10__13C100D010  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D10__13C100D010;                                                 %  cSerPool1_13C100D010 
		dxdt(106, 1) =  - Vf_SerPool1_GlyPool1_13C101D010_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D10__13C101D010  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D10__13C101D010;                                                 %  cSerPool1_13C101D010 
		dxdt(107, 1) =  - Vf_SerPool1_GlyPool1_13C110D010_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D10__13C110D010  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D10__13C110D010;                                                 %  cSerPool1_13C110D010 
		dxdt(108, 1) =  - Vf_SerPool1_GlyPool1_13C111D010_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D10__13C111D010  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D10__13C111D010;                                                 %  cSerPool1_13C111D010 
		dxdt(109, 1) =  - Vf_SerPool1_GlyPool1_13C000D011_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C0D11__13C000D011  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D11__13C000D011;                                                 %  cSerPool1_13C000D011 
		dxdt(110, 1) =  - Vf_SerPool1_GlyPool1_13C001D011_13C00D00  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D11__13C001D011  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D11__13C001D011;                                                 %  cSerPool1_13C001D011 
		dxdt(111, 1) =  - Vf_SerPool1_GlyPool1_13C010D011_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D11__13C010D011  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D11__13C010D011;                                                 %  cSerPool1_13C010D011 
		dxdt(112, 1) =  - Vf_SerPool1_GlyPool1_13C011D011_13C01D00  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D11__13C011D011  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D11__13C011D011;                                                 %  cSerPool1_13C011D011 
		dxdt(113, 1) =  - Vf_SerPool1_GlyPool1_13C100D011_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D11__13C100D011  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D11__13C100D011;                                                 %  cSerPool1_13C100D011 
		dxdt(114, 1) =  - Vf_SerPool1_GlyPool1_13C101D011_13C10D00  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D11__13C101D011  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D11__13C101D011;                                                 %  cSerPool1_13C101D011 
		dxdt(115, 1) =  - Vf_SerPool1_GlyPool1_13C110D011_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D11__13C110D011  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D11__13C110D011;                                                 %  cSerPool1_13C110D011 
		dxdt(116, 1) =  - Vf_SerPool1_GlyPool1_13C111D011_13C11D00  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D11__13C111D011  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D11__13C111D011;                                                 %  cSerPool1_13C111D011 
		dxdt(117, 1) =  - Vf_SerPool1_GlyPool1_13C000D100_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D00__13C000D100  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D00__13C000D100;                                                 %  cSerPool1_13C000D100 
		dxdt(118, 1) =  - Vf_SerPool1_GlyPool1_13C001D100_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D00__13C001D100  +  Vr_SerPool1_GlyPool1_13C00D11_13C1D00__13C001D100;                                                 %  cSerPool1_13C001D100 
		dxdt(119, 1) =  - Vf_SerPool1_GlyPool1_13C010D100_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D00__13C010D100  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D00__13C010D100;                                                 %  cSerPool1_13C010D100 
		dxdt(120, 1) =  - Vf_SerPool1_GlyPool1_13C011D100_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D00__13C011D100  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D00__13C011D100;                                                 %  cSerPool1_13C011D100 
		dxdt(121, 1) =  - Vf_SerPool1_GlyPool1_13C100D100_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D00__13C100D100  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D00__13C100D100;                                                 %  cSerPool1_13C100D100 
		dxdt(122, 1) =  - Vf_SerPool1_GlyPool1_13C101D100_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D00__13C101D100  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D00__13C101D100;                                                 %  cSerPool1_13C101D100 
		dxdt(123, 1) =  - Vf_SerPool1_GlyPool1_13C110D100_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D00__13C110D100  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D00__13C110D100;                                                 %  cSerPool1_13C110D100 
		dxdt(124, 1) =  - Vf_SerPool1_GlyPool1_13C111D100_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D00__13C111D100  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D00__13C111D100;                                                 %  cSerPool1_13C111D100 
		dxdt(125, 1) =  - Vf_SerPool1_GlyPool1_13C000D101_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D01__13C000D101  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D01__13C000D101;                                                 %  cSerPool1_13C000D101 
		dxdt(126, 1) =  - Vf_SerPool1_GlyPool1_13C001D101_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D01__13C001D101  +  Vr_SerPool1_GlyPool1_13C00D11_13C1D01__13C001D101;                                                 %  cSerPool1_13C001D101 
		dxdt(127, 1) =  - Vf_SerPool1_GlyPool1_13C010D101_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D01__13C010D101  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D01__13C010D101;                                                 %  cSerPool1_13C010D101 
		dxdt(128, 1) =  - Vf_SerPool1_GlyPool1_13C011D101_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D01__13C011D101  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D01__13C011D101;                                                 %  cSerPool1_13C011D101 
		dxdt(129, 1) =  - Vf_SerPool1_GlyPool1_13C100D101_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D01__13C100D101  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D01__13C100D101;                                                 %  cSerPool1_13C100D101 
		dxdt(130, 1) =  - Vf_SerPool1_GlyPool1_13C101D101_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D01__13C101D101  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D01__13C101D101;                                                 %  cSerPool1_13C101D101 
		dxdt(131, 1) =  - Vf_SerPool1_GlyPool1_13C110D101_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D01__13C110D101  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D01__13C110D101;                                                 %  cSerPool1_13C110D101 
		dxdt(132, 1) =  - Vf_SerPool1_GlyPool1_13C111D101_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D01__13C111D101  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D01__13C111D101;                                                 %  cSerPool1_13C111D101 
		dxdt(133, 1) =  - Vf_SerPool1_GlyPool1_13C000D110_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D10__13C000D110  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D10__13C000D110;                                                 %  cSerPool1_13C000D110 
		dxdt(134, 1) =  - Vf_SerPool1_GlyPool1_13C001D110_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D10__13C001D110  +  Vr_SerPool1_GlyPool1_13C00D11_13C1D10__13C001D110;                                                 %  cSerPool1_13C001D110 
		dxdt(135, 1) =  - Vf_SerPool1_GlyPool1_13C010D110_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D10__13C010D110  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D10__13C010D110;                                                 %  cSerPool1_13C010D110 
		dxdt(136, 1) =  - Vf_SerPool1_GlyPool1_13C011D110_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D10__13C011D110  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D10__13C011D110;                                                 %  cSerPool1_13C011D110 
		dxdt(137, 1) =  - Vf_SerPool1_GlyPool1_13C100D110_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D10__13C100D110  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D10__13C100D110;                                                 %  cSerPool1_13C100D110 
		dxdt(138, 1) =  - Vf_SerPool1_GlyPool1_13C101D110_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D10__13C101D110  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D10__13C101D110;                                                 %  cSerPool1_13C101D110 
		dxdt(139, 1) =  - Vf_SerPool1_GlyPool1_13C110D110_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D10__13C110D110  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D10__13C110D110;                                                 %  cSerPool1_13C110D110 
		dxdt(140, 1) =  - Vf_SerPool1_GlyPool1_13C111D110_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D10__13C111D110  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D10__13C111D110;                                                 %  cSerPool1_13C111D110 
		dxdt(141, 1) =  - Vf_SerPool1_GlyPool1_13C000D111_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D11__13C000D111  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D11__13C000D111;                                                 %  cSerPool1_13C000D111 
		dxdt(142, 1) =  - Vf_SerPool1_GlyPool1_13C001D111_13C00D10  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D11__13C001D111  +  Vr_SerPool1_GlyPool1_13C00D11_13C1D11__13C001D111;                                                 %  cSerPool1_13C001D111 
		dxdt(143, 1) =  - Vf_SerPool1_GlyPool1_13C010D111_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D11__13C010D111  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D11__13C010D111;                                                 %  cSerPool1_13C010D111 
		dxdt(144, 1) =  - Vf_SerPool1_GlyPool1_13C011D111_13C01D10  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D11__13C011D111  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D11__13C011D111;                                                 %  cSerPool1_13C011D111 
		dxdt(145, 1) =  - Vf_SerPool1_GlyPool1_13C100D111_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D11__13C100D111  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D11__13C100D111;                                                 %  cSerPool1_13C100D111 
		dxdt(146, 1) =  - Vf_SerPool1_GlyPool1_13C101D111_13C10D10  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D11__13C101D111  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D11__13C101D111;                                                 %  cSerPool1_13C101D111 
		dxdt(147, 1) =  - Vf_SerPool1_GlyPool1_13C110D111_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D11__13C110D111  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D11__13C110D111;                                                 %  cSerPool1_13C110D111 
		dxdt(148, 1) =  - Vf_SerPool1_GlyPool1_13C111D111_13C11D10  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D11__13C111D111  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D11__13C111D111;                                                 %  cSerPool1_13C111D111 
           


    % cGlyPool1
    %  dxdt(1, 1) = (cGlyPool1) = 
    % - Vf_GlyPool1_CO2 
    % - Vf_PRPP_GAR_GlyPool1 
    % - Vr_SerPool1_GlyPool1 
    % + Vr_GlyPool1_CO2 
    % + Vf_SerPool1_GlyPool1 

		dxdt(149, 1) =  - Vf_GlyPool1_CO2_13C00D00_13C0D00  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D00__13C000  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D00__13C100  -	Vr_SerPool1_GlyPool1_13C00D00_13C0D00__13C000D000  -	Vr_SerPool1_GlyPool1_13C00D00_13C1D00__13C001D000  - 	Vr_SerPool1_GlyPool1_13C00D00_13C0D01__13C000D001  -	Vr_SerPool1_GlyPool1_13C00D00_13C1D01__13C001D001  -	Vr_SerPool1_GlyPool1_13C00D00_13C0D10__13C000D010  -	Vr_SerPool1_GlyPool1_13C00D00_13C1D10__13C001D010  -	Vr_SerPool1_GlyPool1_13C00D00_13C0D11__13C000D011  -	Vr_SerPool1_GlyPool1_13C00D00_13C1D11__13C001D011  +  Vr_GlyPool1_CO2_13C0D00_13C00D00  +  Vf_SerPool1_GlyPool1_13C000D000_13C00D00  +  Vf_SerPool1_GlyPool1_13C001D000_13C00D00  +  Vf_SerPool1_GlyPool1_13C000D001_13C00D00  +  Vf_SerPool1_GlyPool1_13C001D001_13C00D00  +  Vf_SerPool1_GlyPool1_13C000D010_13C00D00  +  Vf_SerPool1_GlyPool1_13C001D010_13C00D00  +  Vf_SerPool1_GlyPool1_13C000D011_13C00D00  +  Vf_SerPool1_GlyPool1_13C001D011_13C00D00;								%  cGlyPool1_13C00D00 
		dxdt(150, 1) =  - Vf_GlyPool1_CO2_13C01D00_13C1D00  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D00__13C001  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D00__13C101  -	Vr_SerPool1_GlyPool1_13C01D00_13C0D00__13C010D000  -	Vr_SerPool1_GlyPool1_13C01D00_13C1D00__13C011D000  - 	Vr_SerPool1_GlyPool1_13C01D00_13C0D01__13C010D001  -	Vr_SerPool1_GlyPool1_13C01D00_13C1D01__13C011D001  -	Vr_SerPool1_GlyPool1_13C01D00_13C0D10__13C010D010  -	Vr_SerPool1_GlyPool1_13C01D00_13C1D10__13C011D010  -	Vr_SerPool1_GlyPool1_13C01D00_13C0D11__13C010D011  -	Vr_SerPool1_GlyPool1_13C01D00_13C1D11__13C011D011  +  Vr_GlyPool1_CO2_13C1D00_13C01D00  +  Vf_SerPool1_GlyPool1_13C010D000_13C01D00  +  Vf_SerPool1_GlyPool1_13C011D000_13C01D00  +  Vf_SerPool1_GlyPool1_13C010D001_13C01D00  +  Vf_SerPool1_GlyPool1_13C011D001_13C01D00  +  Vf_SerPool1_GlyPool1_13C010D010_13C01D00  +  Vf_SerPool1_GlyPool1_13C011D010_13C01D00  +  Vf_SerPool1_GlyPool1_13C010D011_13C01D00  +  Vf_SerPool1_GlyPool1_13C011D011_13C01D00;               %  cGlyPool1_13C01D00 
		dxdt(151, 1) =  - Vf_GlyPool1_CO2_13C10D00_13C0D00  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D00__13C110  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D00__13C010  -	Vr_SerPool1_GlyPool1_13C10D00_13C0D00__13C100D000  -	Vr_SerPool1_GlyPool1_13C10D00_13C1D00__13C101D000  - 	Vr_SerPool1_GlyPool1_13C10D00_13C0D01__13C100D001  -	Vr_SerPool1_GlyPool1_13C10D00_13C1D01__13C101D001  -	Vr_SerPool1_GlyPool1_13C10D00_13C0D10__13C100D010  -	Vr_SerPool1_GlyPool1_13C10D00_13C1D10__13C101D010  -	Vr_SerPool1_GlyPool1_13C10D00_13C0D11__13C100D011  -	Vr_SerPool1_GlyPool1_13C10D00_13C1D11__13C101D011  +                                       Vf_SerPool1_GlyPool1_13C100D000_13C10D00  +  Vf_SerPool1_GlyPool1_13C101D000_13C10D00  +  Vf_SerPool1_GlyPool1_13C100D001_13C10D00  +  Vf_SerPool1_GlyPool1_13C101D001_13C10D00  +  Vf_SerPool1_GlyPool1_13C100D010_13C10D00  +  Vf_SerPool1_GlyPool1_13C101D010_13C10D00  +  Vf_SerPool1_GlyPool1_13C100D011_13C10D00  +  Vf_SerPool1_GlyPool1_13C101D011_13C10D00;               %  cGlyPool1_13C10D00 
		dxdt(152, 1) =  - Vf_GlyPool1_CO2_13C11D00_13C1D00  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D00__13C011  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D00__13C111  -	Vr_SerPool1_GlyPool1_13C11D00_13C0D00__13C110D000  -	Vr_SerPool1_GlyPool1_13C11D00_13C1D00__13C111D000  - 	Vr_SerPool1_GlyPool1_13C11D00_13C0D01__13C110D001  -	Vr_SerPool1_GlyPool1_13C11D00_13C1D01__13C111D001  -	Vr_SerPool1_GlyPool1_13C11D00_13C0D10__13C110D010  -	Vr_SerPool1_GlyPool1_13C11D00_13C1D10__13C111D010  -	Vr_SerPool1_GlyPool1_13C11D00_13C0D11__13C110D011  -	Vr_SerPool1_GlyPool1_13C11D00_13C1D11__13C111D011  +                                       Vf_SerPool1_GlyPool1_13C110D000_13C11D00  +  Vf_SerPool1_GlyPool1_13C111D000_13C11D00  +  Vf_SerPool1_GlyPool1_13C110D001_13C11D00  +  Vf_SerPool1_GlyPool1_13C111D001_13C11D00  +  Vf_SerPool1_GlyPool1_13C110D010_13C11D00  +  Vf_SerPool1_GlyPool1_13C111D010_13C11D00  +  Vf_SerPool1_GlyPool1_13C110D011_13C11D00  +  Vf_SerPool1_GlyPool1_13C111D011_13C11D00;               %  cGlyPool1_13C11D00 
		dxdt(153, 1) =  - Vf_GlyPool1_CO2_13C00D01_13C0D01  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D01__13C000  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D01__13C100  -	Vr_SerPool1_GlyPool1_13C00D01_13C0D00__13C000D000  -	Vr_SerPool1_GlyPool1_13C00D01_13C1D00__13C001D000  - 	Vr_SerPool1_GlyPool1_13C00D01_13C0D01__13C000D001  -	Vr_SerPool1_GlyPool1_13C00D01_13C1D01__13C001D001  -	Vr_SerPool1_GlyPool1_13C00D01_13C0D10__13C000D010  -	Vr_SerPool1_GlyPool1_13C00D01_13C1D10__13C001D010  -	Vr_SerPool1_GlyPool1_13C00D01_13C0D11__13C000D011  -	Vr_SerPool1_GlyPool1_13C00D01_13C1D11__13C001D011  +  Vr_GlyPool1_CO2_13C0D01_13C00D01;                                                                                                                                                                                                                                                                                                                                                                                       %  cGlyPool1_13C00D01 
		dxdt(154, 1) =  - Vf_GlyPool1_CO2_13C01D01_13C1D01  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D01__13C001  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D01__13C101  -	Vr_SerPool1_GlyPool1_13C01D01_13C0D00__13C010D000  -	Vr_SerPool1_GlyPool1_13C01D01_13C1D00__13C011D000  - 	Vr_SerPool1_GlyPool1_13C01D01_13C0D01__13C010D001  -	Vr_SerPool1_GlyPool1_13C01D01_13C1D01__13C011D001  -	Vr_SerPool1_GlyPool1_13C01D01_13C0D10__13C010D010  -	Vr_SerPool1_GlyPool1_13C01D01_13C1D10__13C011D010  -	Vr_SerPool1_GlyPool1_13C01D01_13C0D11__13C010D011  -	Vr_SerPool1_GlyPool1_13C01D01_13C1D11__13C011D011  +  Vr_GlyPool1_CO2_13C1D01_13C01D01;                                                                                                                                                                                                                                                                                                                                                                                       %  cGlyPool1_13C01D01 
		dxdt(155, 1) =  - Vf_GlyPool1_CO2_13C10D01_13C0D01  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D01__13C010  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D01__13C110  -	Vr_SerPool1_GlyPool1_13C10D01_13C0D00__13C100D000  -	Vr_SerPool1_GlyPool1_13C10D01_13C1D00__13C101D000  - 	Vr_SerPool1_GlyPool1_13C10D01_13C0D01__13C100D001  -	Vr_SerPool1_GlyPool1_13C10D01_13C1D01__13C101D001  -	Vr_SerPool1_GlyPool1_13C10D01_13C0D10__13C100D010  -	Vr_SerPool1_GlyPool1_13C10D01_13C1D10__13C101D010  -	Vr_SerPool1_GlyPool1_13C10D01_13C0D11__13C100D011  -	Vr_SerPool1_GlyPool1_13C10D01_13C1D11__13C101D011  ;                                                                                                                                                                                                                                                                                                                                                                                                                          %  cGlyPool1_13C10D01 
		dxdt(156, 1) =  - Vf_GlyPool1_CO2_13C11D01_13C1D01  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D01__13C011  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D01__13C111  -	Vr_SerPool1_GlyPool1_13C11D01_13C0D00__13C110D000  -	Vr_SerPool1_GlyPool1_13C11D01_13C1D00__13C111D000  - 	Vr_SerPool1_GlyPool1_13C11D01_13C0D01__13C110D001  -	Vr_SerPool1_GlyPool1_13C11D01_13C1D01__13C111D001  -	Vr_SerPool1_GlyPool1_13C11D01_13C0D10__13C110D010  -	Vr_SerPool1_GlyPool1_13C11D01_13C1D10__13C111D010  -	Vr_SerPool1_GlyPool1_13C11D01_13C0D11__13C110D011  -	Vr_SerPool1_GlyPool1_13C11D01_13C1D11__13C111D011  ;                                                                                                                                                                                                                                                                                                                                                                                                                          %  cGlyPool1_13C11D01 
		dxdt(157, 1) =  - Vf_GlyPool1_CO2_13C00D10_13C0D10  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D10__13C000  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D10__13C100  -	Vr_SerPool1_GlyPool1_13C00D10_13C0D00__13C000D100  -	Vr_SerPool1_GlyPool1_13C00D10_13C1D00__13C001D100  - 	Vr_SerPool1_GlyPool1_13C00D10_13C0D01__13C000D101  -	Vr_SerPool1_GlyPool1_13C00D10_13C1D01__13C001D101  -	Vr_SerPool1_GlyPool1_13C00D10_13C0D10__13C000D110  -	Vr_SerPool1_GlyPool1_13C00D10_13C1D10__13C001D110  -	Vr_SerPool1_GlyPool1_13C00D10_13C0D11__13C000D111  -	Vr_SerPool1_GlyPool1_13C00D10_13C1D11__13C001D111  +  Vr_GlyPool1_CO2_13C0D10_13C00D10  +  Vf_SerPool1_GlyPool1_13C000D100_13C00D10  +  Vf_SerPool1_GlyPool1_13C001D100_13C00D10  +  Vf_SerPool1_GlyPool1_13C000D101_13C00D10  +  Vf_SerPool1_GlyPool1_13C001D101_13C00D10  +  Vf_SerPool1_GlyPool1_13C000D110_13C00D10  +  Vf_SerPool1_GlyPool1_13C001D110_13C00D10  +  Vf_SerPool1_GlyPool1_13C000D111_13C00D10  +  Vf_SerPool1_GlyPool1_13C001D111_13C00D10;               %  cGlyPool1_13C00D10 
		dxdt(158, 1) =  - Vf_GlyPool1_CO2_13C01D10_13C1D10  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D10__13C001  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D10__13C101  -	Vr_SerPool1_GlyPool1_13C01D10_13C0D00__13C010D100  -	Vr_SerPool1_GlyPool1_13C01D10_13C1D00__13C011D100  - 	Vr_SerPool1_GlyPool1_13C01D10_13C0D01__13C010D101  -	Vr_SerPool1_GlyPool1_13C01D10_13C1D01__13C011D101  -	Vr_SerPool1_GlyPool1_13C01D10_13C0D10__13C010D110  -	Vr_SerPool1_GlyPool1_13C01D10_13C1D10__13C011D110  -	Vr_SerPool1_GlyPool1_13C01D10_13C0D11__13C010D111  -	Vr_SerPool1_GlyPool1_13C01D10_13C1D11__13C011D111  +  Vr_GlyPool1_CO2_13C1D10_13C01D10  +  Vf_SerPool1_GlyPool1_13C010D100_13C01D10  +  Vf_SerPool1_GlyPool1_13C011D100_13C01D10  +  Vf_SerPool1_GlyPool1_13C010D101_13C01D10  +  Vf_SerPool1_GlyPool1_13C011D101_13C01D10  +  Vf_SerPool1_GlyPool1_13C010D110_13C01D10  +  Vf_SerPool1_GlyPool1_13C011D110_13C01D10  +  Vf_SerPool1_GlyPool1_13C010D111_13C01D10  +  Vf_SerPool1_GlyPool1_13C011D111_13C01D10;               %  cGlyPool1_13C01D10 
		dxdt(159, 1) =  - Vf_GlyPool1_CO2_13C10D10_13C0D10  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D10__13C010  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D10__13C110  -	Vr_SerPool1_GlyPool1_13C10D10_13C0D00__13C100D100  -	Vr_SerPool1_GlyPool1_13C10D10_13C1D00__13C101D100  - 	Vr_SerPool1_GlyPool1_13C10D10_13C0D01__13C100D101  -	Vr_SerPool1_GlyPool1_13C10D10_13C1D01__13C101D101  -	Vr_SerPool1_GlyPool1_13C10D10_13C0D10__13C100D110  -	Vr_SerPool1_GlyPool1_13C10D10_13C1D10__13C101D110  -	Vr_SerPool1_GlyPool1_13C10D10_13C0D11__13C100D111  -	Vr_SerPool1_GlyPool1_13C10D10_13C1D11__13C101D111  +                                       Vf_SerPool1_GlyPool1_13C100D100_13C10D10  +  Vf_SerPool1_GlyPool1_13C101D100_13C10D10  +  Vf_SerPool1_GlyPool1_13C100D101_13C10D10  +  Vf_SerPool1_GlyPool1_13C101D101_13C10D10  +  Vf_SerPool1_GlyPool1_13C100D110_13C10D10  +  Vf_SerPool1_GlyPool1_13C101D110_13C10D10  +  Vf_SerPool1_GlyPool1_13C100D111_13C10D10  +  Vf_SerPool1_GlyPool1_13C101D111_13C10D10;               %  cGlyPool1_13C10D10 
		dxdt(160, 1) =  - Vf_GlyPool1_CO2_13C11D10_13C1D10  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D10__13C011  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D10__13C111  -	Vr_SerPool1_GlyPool1_13C11D10_13C0D00__13C110D100  -	Vr_SerPool1_GlyPool1_13C11D10_13C1D00__13C111D100  - 	Vr_SerPool1_GlyPool1_13C11D10_13C0D01__13C110D101  -	Vr_SerPool1_GlyPool1_13C11D10_13C1D01__13C111D101  -	Vr_SerPool1_GlyPool1_13C11D10_13C0D10__13C110D110  -	Vr_SerPool1_GlyPool1_13C11D10_13C1D10__13C111D110  -	Vr_SerPool1_GlyPool1_13C11D10_13C0D11__13C110D111  -	Vr_SerPool1_GlyPool1_13C11D10_13C1D11__13C111D111  +                                       Vf_SerPool1_GlyPool1_13C110D100_13C11D10  +  Vf_SerPool1_GlyPool1_13C111D100_13C11D10  +  Vf_SerPool1_GlyPool1_13C110D101_13C11D10  +  Vf_SerPool1_GlyPool1_13C111D101_13C11D10  +  Vf_SerPool1_GlyPool1_13C110D110_13C11D10  +  Vf_SerPool1_GlyPool1_13C111D110_13C11D10  +  Vf_SerPool1_GlyPool1_13C110D111_13C11D10  +  Vf_SerPool1_GlyPool1_13C111D111_13C11D10;               %  cGlyPool1_13C11D10 
		dxdt(161, 1) =  - Vf_GlyPool1_CO2_13C00D11_13C0D11  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D11__13C000  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D11__13C100  -	Vr_SerPool1_GlyPool1_13C00D11_13C0D00__13C000D100  -	Vr_SerPool1_GlyPool1_13C00D11_13C1D00__13C001D100  - 	Vr_SerPool1_GlyPool1_13C00D11_13C0D01__13C000D101  -	Vr_SerPool1_GlyPool1_13C00D11_13C1D01__13C001D101  -	Vr_SerPool1_GlyPool1_13C00D11_13C0D10__13C000D110  -	Vr_SerPool1_GlyPool1_13C00D11_13C1D10__13C001D110  -	Vr_SerPool1_GlyPool1_13C00D11_13C0D11__13C000D111  -	Vr_SerPool1_GlyPool1_13C00D11_13C1D11__13C001D111  +  Vr_GlyPool1_CO2_13C0D11_13C00D11;                                                                                                                                                                                                                                                                                                                                                                                       %  cGlyPool1_13C00D11 
		dxdt(162, 1) =  - Vf_GlyPool1_CO2_13C01D11_13C1D11  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D11__13C001  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D11__13C101  -	Vr_SerPool1_GlyPool1_13C01D11_13C0D00__13C010D100  -	Vr_SerPool1_GlyPool1_13C01D11_13C1D00__13C011D100  - 	Vr_SerPool1_GlyPool1_13C01D11_13C0D01__13C010D101  -	Vr_SerPool1_GlyPool1_13C01D11_13C1D01__13C011D101  -	Vr_SerPool1_GlyPool1_13C01D11_13C0D10__13C010D110  -	Vr_SerPool1_GlyPool1_13C01D11_13C1D10__13C011D110  -	Vr_SerPool1_GlyPool1_13C01D11_13C0D11__13C010D111  -	Vr_SerPool1_GlyPool1_13C01D11_13C1D11__13C011D111  +  Vr_GlyPool1_CO2_13C1D11_13C01D11;                                                                                                                                                                                                                                                                                                                                                                                       %  cGlyPool1_13C01D11 
		dxdt(163, 1) =  - Vf_GlyPool1_CO2_13C10D11_13C0D11  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D11__13C010  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D11__13C110  -	Vr_SerPool1_GlyPool1_13C10D11_13C0D00__13C100D100  -	Vr_SerPool1_GlyPool1_13C10D11_13C1D00__13C101D100  - 	Vr_SerPool1_GlyPool1_13C10D11_13C0D01__13C100D101  -	Vr_SerPool1_GlyPool1_13C10D11_13C1D01__13C101D101  -	Vr_SerPool1_GlyPool1_13C10D11_13C0D10__13C100D110  -	Vr_SerPool1_GlyPool1_13C10D11_13C1D10__13C101D110  -	Vr_SerPool1_GlyPool1_13C10D11_13C0D11__13C100D111  -	Vr_SerPool1_GlyPool1_13C10D11_13C1D11__13C101D111  ;                                                                                                                                                                                                                                                                                                                                                                                                                          %  cGlyPool1_13C10D11 
		dxdt(164, 1) =  - Vf_GlyPool1_CO2_13C11D11_13C1D11  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D11__13C011  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D11__13C111  -	Vr_SerPool1_GlyPool1_13C11D11_13C0D00__13C110D100  -	Vr_SerPool1_GlyPool1_13C11D11_13C1D00__13C111D100  - 	Vr_SerPool1_GlyPool1_13C11D11_13C0D01__13C110D101  -	Vr_SerPool1_GlyPool1_13C11D11_13C1D01__13C111D101  -	Vr_SerPool1_GlyPool1_13C11D11_13C0D10__13C110D110  -	Vr_SerPool1_GlyPool1_13C11D11_13C1D10__13C111D110  -	Vr_SerPool1_GlyPool1_13C11D11_13C0D11__13C110D111  -	Vr_SerPool1_GlyPool1_13C11D11_13C1D11__13C111D111  ;                                                                                                                                                                                                                                                                                                                                                                                                                          %  cGlyPool1_13C11D11 



    % cSerPool2
    %  dxdt(1, 1) = (cSerPool2) = 
    % - Vf_SerPool2_GlyPool2 
    % - Vf_SerPool2_SerPoolMitochon 
    % - Vr_Vin_Ser 
    % + Vr_SerPool2_GlyPool2 
    % + Vr_SerPool2_SerPoolMitochon 
    % + Vf_Vin_Ser 
    
		dxdt(165, 1) =  - Vf_SerPool2_GlyPool2_13C000D000_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C000D000_13C000D000  -  Vr_Vin_Ser_13C000D000_13C000D000  +  Vr_SerPool2_GlyPool2_13C00D00_13C0D00__13C000D000  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D00__13C000D000  +  Vr_SerPool2_SerPoolMitochon_13C000D000_13C000D000  +  Vf_Vin_Ser_13C000D000_13C000D000;						%  cSerPool2_13C000D000 
		dxdt(166, 1) =  - Vf_SerPool2_GlyPool2_13C001D000_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C001D000_13C001D000  -  Vr_Vin_Ser_13C001D000_13C001D000  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D00__13C001D000  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D00__13C001D000  +  Vr_SerPool2_SerPoolMitochon_13C001D000_13C001D000  +  Vf_Vin_Ser_13C001D000_13C001D000;           %  cSerPool2_13C001D000 
		dxdt(167, 1) =  - Vf_SerPool2_GlyPool2_13C010D000_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C010D000_13C010D000  -  Vr_Vin_Ser_13C010D000_13C010D000  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D00__13C010D000  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D00__13C010D000  +  Vr_SerPool2_SerPoolMitochon_13C010D000_13C010D000  +  Vf_Vin_Ser_13C010D000_13C010D000;           %  cSerPool2_13C010D000 
		dxdt(168, 1) =  - Vf_SerPool2_GlyPool2_13C011D000_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C011D000_13C011D000  -  Vr_Vin_Ser_13C011D000_13C011D000  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D00__13C011D000  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D00__13C011D000  +  Vr_SerPool2_SerPoolMitochon_13C011D000_13C011D000  +  Vf_Vin_Ser_13C011D000_13C011D000;           %  cSerPool2_13C011D000 
		dxdt(169, 1) =  - Vf_SerPool2_GlyPool2_13C100D000_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C100D000_13C100D000  -  Vr_Vin_Ser_13C100D000_13C100D000  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D00__13C100D000  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D00__13C100D000  +  Vr_SerPool2_SerPoolMitochon_13C100D000_13C100D000  +  Vf_Vin_Ser_13C100D000_13C100D000;           %  cSerPool2_13C100D000 
		dxdt(170, 1) =  - Vf_SerPool2_GlyPool2_13C101D000_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C101D000_13C101D000  -  Vr_Vin_Ser_13C101D000_13C101D000  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D00__13C101D000  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D00__13C101D000  +  Vr_SerPool2_SerPoolMitochon_13C101D000_13C101D000  +  Vf_Vin_Ser_13C101D000_13C101D000;           %  cSerPool2_13C101D000 
		dxdt(171, 1) =  - Vf_SerPool2_GlyPool2_13C110D000_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C110D000_13C110D000  -  Vr_Vin_Ser_13C110D000_13C110D000  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D00__13C110D000  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D00__13C110D000  +  Vr_SerPool2_SerPoolMitochon_13C110D000_13C110D000  +  Vf_Vin_Ser_13C110D000_13C110D000;           %  cSerPool2_13C110D000 
		dxdt(172, 1) =  - Vf_SerPool2_GlyPool2_13C111D000_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C111D000_13C111D000  -  Vr_Vin_Ser_13C111D000_13C111D000  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D00__13C111D000  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D00__13C111D000  +  Vr_SerPool2_SerPoolMitochon_13C111D000_13C111D000  +  Vf_Vin_Ser_13C111D000_13C111D000;           %  cSerPool2_13C111D000 
		dxdt(173, 1) =  - Vf_SerPool2_GlyPool2_13C000D001_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C000D001_13C000D001  -  Vr_Vin_Ser_13C000D001_13C000D001  +  Vr_SerPool2_GlyPool2_13C00D00_13C0D01__13C000D001  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D01__13C000D001  +  Vr_SerPool2_SerPoolMitochon_13C000D001_13C000D001  +  Vf_Vin_Ser_13C000D001_13C000D001;           %  cSerPool2_13C000D001 
		dxdt(174, 1) =  - Vf_SerPool2_GlyPool2_13C001D001_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C001D001_13C001D001  -  Vr_Vin_Ser_13C001D001_13C001D001  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D01__13C001D001  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D01__13C001D001  +  Vr_SerPool2_SerPoolMitochon_13C001D001_13C001D001  +  Vf_Vin_Ser_13C001D001_13C001D001;           %  cSerPool2_13C001D001 
		dxdt(175, 1) =  - Vf_SerPool2_GlyPool2_13C010D001_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C010D001_13C010D001  -  Vr_Vin_Ser_13C010D001_13C010D001  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D01__13C010D001  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D01__13C010D001  +  Vr_SerPool2_SerPoolMitochon_13C010D001_13C010D001  +  Vf_Vin_Ser_13C010D001_13C010D001;           %  cSerPool2_13C010D001 
		dxdt(176, 1) =  - Vf_SerPool2_GlyPool2_13C011D001_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C011D001_13C011D001  -  Vr_Vin_Ser_13C011D001_13C011D001  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D01__13C011D001  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D01__13C011D001  +  Vr_SerPool2_SerPoolMitochon_13C011D001_13C011D001  +  Vf_Vin_Ser_13C011D001_13C011D001;           %  cSerPool2_13C011D001 
		dxdt(177, 1) =  - Vf_SerPool2_GlyPool2_13C100D001_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C100D001_13C100D001  -  Vr_Vin_Ser_13C100D001_13C100D001  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D01__13C100D001  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D01__13C100D001  +  Vr_SerPool2_SerPoolMitochon_13C100D001_13C100D001  +  Vf_Vin_Ser_13C100D001_13C100D001;           %  cSerPool2_13C100D001 
		dxdt(178, 1) =  - Vf_SerPool2_GlyPool2_13C101D001_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C101D001_13C101D001  -  Vr_Vin_Ser_13C101D001_13C101D001  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D01__13C101D001  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D01__13C101D001  +  Vr_SerPool2_SerPoolMitochon_13C101D001_13C101D001  +  Vf_Vin_Ser_13C101D001_13C101D001;           %  cSerPool2_13C101D001 
		dxdt(179, 1) =  - Vf_SerPool2_GlyPool2_13C110D001_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C110D001_13C110D001  -  Vr_Vin_Ser_13C110D001_13C110D001  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D01__13C110D001  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D01__13C110D001  +  Vr_SerPool2_SerPoolMitochon_13C110D001_13C110D001  +  Vf_Vin_Ser_13C110D001_13C110D001;           %  cSerPool2_13C110D001 
		dxdt(180, 1) =  - Vf_SerPool2_GlyPool2_13C111D001_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C111D001_13C111D001  -  Vr_Vin_Ser_13C111D001_13C111D001  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D01__13C111D001  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D01__13C111D001  +  Vr_SerPool2_SerPoolMitochon_13C111D001_13C111D001  +  Vf_Vin_Ser_13C111D001_13C111D001;           %  cSerPool2_13C111D001 
		dxdt(181, 1) =  - Vf_SerPool2_GlyPool2_13C000D010_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C000D010_13C000D010  -  Vr_Vin_Ser_13C000D010_13C000D010  +  Vr_SerPool2_GlyPool2_13C00D00_13C0D10__13C000D010  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D10__13C000D010  +  Vr_SerPool2_SerPoolMitochon_13C000D010_13C000D010  +  Vf_Vin_Ser_13C000D010_13C000D010;           %  cSerPool2_13C000D010 
		dxdt(182, 1) =  - Vf_SerPool2_GlyPool2_13C001D010_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C001D010_13C001D010  -  Vr_Vin_Ser_13C001D010_13C001D010  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D10__13C001D010  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D10__13C001D010  +  Vr_SerPool2_SerPoolMitochon_13C001D010_13C001D010  +  Vf_Vin_Ser_13C001D010_13C001D010;           %  cSerPool2_13C001D010 
		dxdt(183, 1) =  - Vf_SerPool2_GlyPool2_13C010D010_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C010D010_13C010D010  -  Vr_Vin_Ser_13C010D010_13C010D010  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D10__13C010D010  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D10__13C010D010  +  Vr_SerPool2_SerPoolMitochon_13C010D010_13C010D010  +  Vf_Vin_Ser_13C010D010_13C010D010;           %  cSerPool2_13C010D010 
		dxdt(184, 1) =  - Vf_SerPool2_GlyPool2_13C011D010_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C011D010_13C011D010  -  Vr_Vin_Ser_13C011D010_13C011D010  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D10__13C011D010  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D10__13C011D010  +  Vr_SerPool2_SerPoolMitochon_13C011D010_13C011D010  +  Vf_Vin_Ser_13C011D010_13C011D010;           %  cSerPool2_13C011D010 
		dxdt(185, 1) =  - Vf_SerPool2_GlyPool2_13C100D010_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C100D010_13C100D010  -  Vr_Vin_Ser_13C100D010_13C100D010  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D10__13C100D010  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D10__13C100D010  +  Vr_SerPool2_SerPoolMitochon_13C100D010_13C100D010  +  Vf_Vin_Ser_13C100D010_13C100D010;           %  cSerPool2_13C100D010 
		dxdt(186, 1) =  - Vf_SerPool2_GlyPool2_13C101D010_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C101D010_13C101D010  -  Vr_Vin_Ser_13C101D010_13C101D010  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D10__13C101D010  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D10__13C101D010  +  Vr_SerPool2_SerPoolMitochon_13C101D010_13C101D010  +  Vf_Vin_Ser_13C101D010_13C101D010;           %  cSerPool2_13C101D010 
		dxdt(187, 1) =  - Vf_SerPool2_GlyPool2_13C110D010_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C110D010_13C110D010  -  Vr_Vin_Ser_13C110D010_13C110D010  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D10__13C110D010  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D10__13C110D010  +  Vr_SerPool2_SerPoolMitochon_13C110D010_13C110D010  +  Vf_Vin_Ser_13C110D010_13C110D010;           %  cSerPool2_13C110D010 
		dxdt(188, 1) =  - Vf_SerPool2_GlyPool2_13C111D010_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C111D010_13C111D010  -  Vr_Vin_Ser_13C111D010_13C111D010  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D10__13C111D010  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D10__13C111D010  +  Vr_SerPool2_SerPoolMitochon_13C111D010_13C111D010  +  Vf_Vin_Ser_13C111D010_13C111D010;           %  cSerPool2_13C111D010 
		dxdt(189, 1) =  - Vf_SerPool2_GlyPool2_13C000D011_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C000D011_13C000D011  -  Vr_Vin_Ser_13C000D011_13C000D011  +  Vr_SerPool2_GlyPool2_13C00D00_13C0D11__13C000D011  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D11__13C000D011  +  Vr_SerPool2_SerPoolMitochon_13C000D011_13C000D011  +  Vf_Vin_Ser_13C000D011_13C000D011;           %  cSerPool2_13C000D011 
		dxdt(190, 1) =  - Vf_SerPool2_GlyPool2_13C001D011_13C00D00  -  Vf_SerPool2_SerPoolMitochon_13C001D011_13C001D011  -  Vr_Vin_Ser_13C001D011_13C001D011  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D11__13C001D011  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D11__13C001D011  +  Vr_SerPool2_SerPoolMitochon_13C001D011_13C001D011  +  Vf_Vin_Ser_13C001D011_13C001D011;           %  cSerPool2_13C001D011 
		dxdt(191, 1) =  - Vf_SerPool2_GlyPool2_13C010D011_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C010D011_13C010D011  -  Vr_Vin_Ser_13C010D011_13C010D011  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D11__13C010D011  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D11__13C010D011  +  Vr_SerPool2_SerPoolMitochon_13C010D011_13C010D011  +  Vf_Vin_Ser_13C010D011_13C010D011;           %  cSerPool2_13C010D011 
		dxdt(192, 1) =  - Vf_SerPool2_GlyPool2_13C011D011_13C01D00  -  Vf_SerPool2_SerPoolMitochon_13C011D011_13C011D011  -  Vr_Vin_Ser_13C011D011_13C011D011  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D11__13C011D011  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D11__13C011D011  +  Vr_SerPool2_SerPoolMitochon_13C011D011_13C011D011  +  Vf_Vin_Ser_13C011D011_13C011D011;           %  cSerPool2_13C011D011 
		dxdt(193, 1) =  - Vf_SerPool2_GlyPool2_13C100D011_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C100D011_13C100D011  -  Vr_Vin_Ser_13C100D011_13C100D011  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D11__13C100D011  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D11__13C100D011  +  Vr_SerPool2_SerPoolMitochon_13C100D011_13C100D011  +  Vf_Vin_Ser_13C100D011_13C100D011;           %  cSerPool2_13C100D011 
		dxdt(194, 1) =  - Vf_SerPool2_GlyPool2_13C101D011_13C10D00  -  Vf_SerPool2_SerPoolMitochon_13C101D011_13C101D011  -  Vr_Vin_Ser_13C101D011_13C101D011  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D11__13C101D011  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D11__13C101D011  +  Vr_SerPool2_SerPoolMitochon_13C101D011_13C101D011  +  Vf_Vin_Ser_13C101D011_13C101D011;           %  cSerPool2_13C101D011 
		dxdt(195, 1) =  - Vf_SerPool2_GlyPool2_13C110D011_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C110D011_13C110D011  -  Vr_Vin_Ser_13C110D011_13C110D011  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D11__13C110D011  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D11__13C110D011  +  Vr_SerPool2_SerPoolMitochon_13C110D011_13C110D011  +  Vf_Vin_Ser_13C110D011_13C110D011;           %  cSerPool2_13C110D011 
		dxdt(196, 1) =  - Vf_SerPool2_GlyPool2_13C111D011_13C11D00  -  Vf_SerPool2_SerPoolMitochon_13C111D011_13C111D011  -  Vr_Vin_Ser_13C111D011_13C111D011  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D11__13C111D011  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D11__13C111D011  +  Vr_SerPool2_SerPoolMitochon_13C111D011_13C111D011  +  Vf_Vin_Ser_13C111D011_13C111D011;           %  cSerPool2_13C111D011 
		dxdt(197, 1) =  - Vf_SerPool2_GlyPool2_13C000D100_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C000D100_13C000D100  -  Vr_Vin_Ser_13C000D100_13C000D100  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D00__13C000D100  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D00__13C000D100  +  Vr_SerPool2_SerPoolMitochon_13C000D100_13C000D100  +  Vf_Vin_Ser_13C000D100_13C000D100;           %  cSerPool2_13C000D100 
		dxdt(198, 1) =  - Vf_SerPool2_GlyPool2_13C001D100_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C001D100_13C001D100  -  Vr_Vin_Ser_13C001D100_13C001D100  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D00__13C001D100  +  Vr_SerPool2_GlyPool2_13C00D11_13C1D00__13C001D100  +  Vr_SerPool2_SerPoolMitochon_13C001D100_13C001D100  +  Vf_Vin_Ser_13C001D100_13C001D100;           %  cSerPool2_13C001D100 
		dxdt(199, 1) =  - Vf_SerPool2_GlyPool2_13C010D100_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C010D100_13C010D100  -  Vr_Vin_Ser_13C010D100_13C010D100  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D00__13C010D100  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D00__13C010D100  +  Vr_SerPool2_SerPoolMitochon_13C010D100_13C010D100  +  Vf_Vin_Ser_13C010D100_13C010D100;           %  cSerPool2_13C010D100 
		dxdt(200, 1) =  - Vf_SerPool2_GlyPool2_13C011D100_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C011D100_13C011D100  -  Vr_Vin_Ser_13C011D100_13C011D100  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D00__13C011D100  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D00__13C011D100  +  Vr_SerPool2_SerPoolMitochon_13C011D100_13C011D100  +  Vf_Vin_Ser_13C011D100_13C011D100;           %  cSerPool2_13C011D100 
		dxdt(201, 1) =  - Vf_SerPool2_GlyPool2_13C100D100_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C100D100_13C100D100  -  Vr_Vin_Ser_13C100D100_13C100D100  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D00__13C100D100  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D00__13C100D100  +  Vr_SerPool2_SerPoolMitochon_13C100D100_13C100D100  +  Vf_Vin_Ser_13C100D100_13C100D100;           %  cSerPool2_13C100D100 
		dxdt(202, 1) =  - Vf_SerPool2_GlyPool2_13C101D100_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C101D100_13C101D100  -  Vr_Vin_Ser_13C101D100_13C101D100  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D00__13C101D100  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D00__13C101D100  +  Vr_SerPool2_SerPoolMitochon_13C101D100_13C101D100  +  Vf_Vin_Ser_13C101D100_13C101D100;           %  cSerPool2_13C101D100 
		dxdt(203, 1) =  - Vf_SerPool2_GlyPool2_13C110D100_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C110D100_13C110D100  -  Vr_Vin_Ser_13C110D100_13C110D100  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D00__13C110D100  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D00__13C110D100  +  Vr_SerPool2_SerPoolMitochon_13C110D100_13C110D100  +  Vf_Vin_Ser_13C110D100_13C110D100;           %  cSerPool2_13C110D100 
		dxdt(204, 1) =  - Vf_SerPool2_GlyPool2_13C111D100_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C111D100_13C111D100  -  Vr_Vin_Ser_13C111D100_13C111D100  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D00__13C111D100  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D00__13C111D100  +  Vr_SerPool2_SerPoolMitochon_13C111D100_13C111D100  +  Vf_Vin_Ser_13C111D100_13C111D100;           %  cSerPool2_13C111D100 
		dxdt(205, 1) =  - Vf_SerPool2_GlyPool2_13C000D101_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C000D101_13C000D101  -  Vr_Vin_Ser_13C000D101_13C000D101  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D01__13C000D101  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D01__13C000D101  +  Vr_SerPool2_SerPoolMitochon_13C000D101_13C000D101  +  Vf_Vin_Ser_13C000D101_13C000D101;           %  cSerPool2_13C000D101 
		dxdt(206, 1) =  - Vf_SerPool2_GlyPool2_13C001D101_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C001D101_13C001D101  -  Vr_Vin_Ser_13C001D101_13C001D101  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D01__13C001D101  +  Vr_SerPool2_GlyPool2_13C00D11_13C1D01__13C001D101  +  Vr_SerPool2_SerPoolMitochon_13C001D101_13C001D101  +  Vf_Vin_Ser_13C001D101_13C001D101;           %  cSerPool2_13C001D101 
		dxdt(207, 1) =  - Vf_SerPool2_GlyPool2_13C010D101_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C010D101_13C010D101  -  Vr_Vin_Ser_13C010D101_13C010D101  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D01__13C010D101  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D01__13C010D101  +  Vr_SerPool2_SerPoolMitochon_13C010D101_13C010D101  +  Vf_Vin_Ser_13C010D101_13C010D101;           %  cSerPool2_13C010D101 
		dxdt(208, 1) =  - Vf_SerPool2_GlyPool2_13C011D101_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C011D101_13C011D101  -  Vr_Vin_Ser_13C011D101_13C011D101  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D01__13C011D101  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D01__13C011D101  +  Vr_SerPool2_SerPoolMitochon_13C011D101_13C011D101  +  Vf_Vin_Ser_13C011D101_13C011D101;           %  cSerPool2_13C011D101 
		dxdt(209, 1) =  - Vf_SerPool2_GlyPool2_13C100D101_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C100D101_13C100D101  -  Vr_Vin_Ser_13C100D101_13C100D101  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D01__13C100D101  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D01__13C100D101  +  Vr_SerPool2_SerPoolMitochon_13C100D101_13C100D101  +  Vf_Vin_Ser_13C100D101_13C100D101;           %  cSerPool2_13C100D101 
		dxdt(210, 1) =  - Vf_SerPool2_GlyPool2_13C101D101_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C101D101_13C101D101  -  Vr_Vin_Ser_13C101D101_13C101D101  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D01__13C101D101  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D01__13C101D101  +  Vr_SerPool2_SerPoolMitochon_13C101D101_13C101D101  +  Vf_Vin_Ser_13C101D101_13C101D101;           %  cSerPool2_13C101D101 
		dxdt(211, 1) =  - Vf_SerPool2_GlyPool2_13C110D101_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C110D101_13C110D101  -  Vr_Vin_Ser_13C110D101_13C110D101  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D01__13C110D101  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D01__13C110D101  +  Vr_SerPool2_SerPoolMitochon_13C110D101_13C110D101  +  Vf_Vin_Ser_13C110D101_13C110D101;           %  cSerPool2_13C110D101 
		dxdt(212, 1) =  - Vf_SerPool2_GlyPool2_13C111D101_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C111D101_13C111D101  -  Vr_Vin_Ser_13C111D101_13C111D101  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D01__13C111D101  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D01__13C111D101  +  Vr_SerPool2_SerPoolMitochon_13C111D101_13C111D101  +  Vf_Vin_Ser_13C111D101_13C111D101;           %  cSerPool2_13C111D101 
		dxdt(213, 1) =  - Vf_SerPool2_GlyPool2_13C000D110_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C000D110_13C000D110  -  Vr_Vin_Ser_13C000D110_13C000D110  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D10__13C000D110  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D10__13C000D110  +  Vr_SerPool2_SerPoolMitochon_13C000D110_13C000D110  +  Vf_Vin_Ser_13C000D110_13C000D110;           %  cSerPool2_13C000D110 
		dxdt(214, 1) =  - Vf_SerPool2_GlyPool2_13C001D110_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C001D110_13C001D110  -  Vr_Vin_Ser_13C001D110_13C001D110  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D10__13C001D110  +  Vr_SerPool2_GlyPool2_13C00D11_13C1D10__13C001D110  +  Vr_SerPool2_SerPoolMitochon_13C001D110_13C001D110  +  Vf_Vin_Ser_13C001D110_13C001D110;           %  cSerPool2_13C001D110 
		dxdt(215, 1) =  - Vf_SerPool2_GlyPool2_13C010D110_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C010D110_13C010D110  -  Vr_Vin_Ser_13C010D110_13C010D110  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D10__13C010D110  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D10__13C010D110  +  Vr_SerPool2_SerPoolMitochon_13C010D110_13C010D110  +  Vf_Vin_Ser_13C010D110_13C010D110;           %  cSerPool2_13C010D110 
		dxdt(216, 1) =  - Vf_SerPool2_GlyPool2_13C011D110_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C011D110_13C011D110  -  Vr_Vin_Ser_13C011D110_13C011D110  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D10__13C011D110  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D10__13C011D110  +  Vr_SerPool2_SerPoolMitochon_13C011D110_13C011D110  +  Vf_Vin_Ser_13C011D110_13C011D110;           %  cSerPool2_13C011D110 
		dxdt(217, 1) =  - Vf_SerPool2_GlyPool2_13C100D110_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C100D110_13C100D110  -  Vr_Vin_Ser_13C100D110_13C100D110  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D10__13C100D110  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D10__13C100D110  +  Vr_SerPool2_SerPoolMitochon_13C100D110_13C100D110  +  Vf_Vin_Ser_13C100D110_13C100D110;           %  cSerPool2_13C100D110 
		dxdt(218, 1) =  - Vf_SerPool2_GlyPool2_13C101D110_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C101D110_13C101D110  -  Vr_Vin_Ser_13C101D110_13C101D110  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D10__13C101D110  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D10__13C101D110  +  Vr_SerPool2_SerPoolMitochon_13C101D110_13C101D110  +  Vf_Vin_Ser_13C101D110_13C101D110;           %  cSerPool2_13C101D110 
		dxdt(219, 1) =  - Vf_SerPool2_GlyPool2_13C110D110_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C110D110_13C110D110  -  Vr_Vin_Ser_13C110D110_13C110D110  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D10__13C110D110  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D10__13C110D110  +  Vr_SerPool2_SerPoolMitochon_13C110D110_13C110D110  +  Vf_Vin_Ser_13C110D110_13C110D110;           %  cSerPool2_13C110D110 
		dxdt(220, 1) =  - Vf_SerPool2_GlyPool2_13C111D110_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C111D110_13C111D110  -  Vr_Vin_Ser_13C111D110_13C111D110  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D10__13C111D110  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D10__13C111D110  +  Vr_SerPool2_SerPoolMitochon_13C111D110_13C111D110  +  Vf_Vin_Ser_13C111D110_13C111D110;           %  cSerPool2_13C111D110 
		dxdt(221, 1) =  - Vf_SerPool2_GlyPool2_13C000D111_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C000D111_13C000D111  -  Vr_Vin_Ser_13C000D111_13C000D111  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D11__13C000D111  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D11__13C000D111  +  Vr_SerPool2_SerPoolMitochon_13C000D111_13C000D111  +  Vf_Vin_Ser_13C000D111_13C000D111;           %  cSerPool2_13C000D111 
		dxdt(222, 1) =  - Vf_SerPool2_GlyPool2_13C001D111_13C00D10  -  Vf_SerPool2_SerPoolMitochon_13C001D111_13C001D111  -  Vr_Vin_Ser_13C001D111_13C001D111  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D11__13C001D111  +  Vr_SerPool2_GlyPool2_13C00D11_13C1D11__13C001D111  +  Vr_SerPool2_SerPoolMitochon_13C001D111_13C001D111  +  Vf_Vin_Ser_13C001D111_13C001D111;           %  cSerPool2_13C001D111 
		dxdt(223, 1) =  - Vf_SerPool2_GlyPool2_13C010D111_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C010D111_13C010D111  -  Vr_Vin_Ser_13C010D111_13C010D111  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D11__13C010D111  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D11__13C010D111  +  Vr_SerPool2_SerPoolMitochon_13C010D111_13C010D111  +  Vf_Vin_Ser_13C010D111_13C010D111;           %  cSerPool2_13C010D111 
		dxdt(224, 1) =  - Vf_SerPool2_GlyPool2_13C011D111_13C01D10  -  Vf_SerPool2_SerPoolMitochon_13C011D111_13C011D111  -  Vr_Vin_Ser_13C011D111_13C011D111  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D11__13C011D111  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D11__13C011D111  +  Vr_SerPool2_SerPoolMitochon_13C011D111_13C011D111  +  Vf_Vin_Ser_13C011D111_13C011D111;           %  cSerPool2_13C011D111 
		dxdt(225, 1) =  - Vf_SerPool2_GlyPool2_13C100D111_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C100D111_13C100D111  -  Vr_Vin_Ser_13C100D111_13C100D111  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D11__13C100D111  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D11__13C100D111  +  Vr_SerPool2_SerPoolMitochon_13C100D111_13C100D111  +  Vf_Vin_Ser_13C100D111_13C100D111;           %  cSerPool2_13C100D111 
		dxdt(226, 1) =  - Vf_SerPool2_GlyPool2_13C101D111_13C10D10  -  Vf_SerPool2_SerPoolMitochon_13C101D111_13C101D111  -  Vr_Vin_Ser_13C101D111_13C101D111  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D11__13C101D111  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D11__13C101D111  +  Vr_SerPool2_SerPoolMitochon_13C101D111_13C101D111  +  Vf_Vin_Ser_13C101D111_13C101D111;           %  cSerPool2_13C101D111 
		dxdt(227, 1) =  - Vf_SerPool2_GlyPool2_13C110D111_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C110D111_13C110D111  -  Vr_Vin_Ser_13C110D111_13C110D111  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D11__13C110D111  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D11__13C110D111  +  Vr_SerPool2_SerPoolMitochon_13C110D111_13C110D111  +  Vf_Vin_Ser_13C110D111_13C110D111;           %  cSerPool2_13C110D111 
		dxdt(228, 1) =  - Vf_SerPool2_GlyPool2_13C111D111_13C11D10  -  Vf_SerPool2_SerPoolMitochon_13C111D111_13C111D111  -  Vr_Vin_Ser_13C111D111_13C111D111  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D11__13C111D111  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D11__13C111D111  +  Vr_SerPool2_SerPoolMitochon_13C111D111_13C111D111  +  Vf_Vin_Ser_13C111D111_13C111D111;           %  cSerPool2_13C111D111 
           

  
    % cGlyPool2
    %  dxdt(1, 1) = (cGlyPool2) = 
    % - Vf_GlyPool2_GlyPoolMitochon     
    % - Vf_GlyPool2_CO2     
    % - Vf_PRPP_GAR_GlyPool2     
    % - Vr_Vin_Gly     
    % - Vr_SerPool2_GlyPool2     
    % + Vr_GlyPool2_GlyPoolMitochon 
    % + Vr_GlyPool2_CO2 
    % + Vf_Vin_Gly 
    % + Vf_SerPool2_GlyPool2 

		dxdt(229, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C00D00_13C00D00  -  Vf_GlyPool2_CO2_13C00D00_13C0D00  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D00__13C000  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D00__13C100  -  Vr_Vin_Gly_13C00D00_13C00D00  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D00__13C000D000  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D00__13C001D000  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D01__13C000D001  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D01__13C001D001  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D10__13C000D010  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D10__13C001D010  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D11__13C000D011  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D11__13C001D011  +  Vr_GlyPool2_GlyPoolMitochon_13C00D00_13C00D00  +  Vr_GlyPool2_CO2_13C0D00_13C00D00  +  Vf_Vin_Gly_13C00D00_13C00D00  +  Vf_SerPool2_GlyPool2_13C000D000_13C00D00  +  Vf_SerPool2_GlyPool2_13C001D000_13C00D00  +  Vf_SerPool2_GlyPool2_13C000D001_13C00D00  +  Vf_SerPool2_GlyPool2_13C001D001_13C00D00  +  Vf_SerPool2_GlyPool2_13C000D010_13C00D00  +  Vf_SerPool2_GlyPool2_13C001D010_13C00D00  +  Vf_SerPool2_GlyPool2_13C000D011_13C00D00  +  Vf_SerPool2_GlyPool2_13C001D011_13C00D00;					%  cGlyPool2_13C00D00 
		dxdt(230, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C01D00_13C01D00  -  Vf_GlyPool2_CO2_13C01D00_13C1D00  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D00__13C001  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D00__13C101  -  Vr_Vin_Gly_13C01D00_13C01D00  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D00__13C010D000  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D00__13C011D000  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D01__13C010D001  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D01__13C011D001  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D10__13C010D010  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D10__13C011D010  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D11__13C010D011  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D11__13C011D011  +  Vr_GlyPool2_GlyPoolMitochon_13C01D00_13C01D00  +  Vr_GlyPool2_CO2_13C1D00_13C01D00  +  Vf_Vin_Gly_13C01D00_13C01D00  +  Vf_SerPool2_GlyPool2_13C010D000_13C01D00  +  Vf_SerPool2_GlyPool2_13C011D000_13C01D00  +  Vf_SerPool2_GlyPool2_13C010D001_13C01D00  +  Vf_SerPool2_GlyPool2_13C011D001_13C01D00  +  Vf_SerPool2_GlyPool2_13C010D010_13C01D00  +  Vf_SerPool2_GlyPool2_13C011D010_13C01D00  +  Vf_SerPool2_GlyPool2_13C010D011_13C01D00  +  Vf_SerPool2_GlyPool2_13C011D011_13C01D00;          %  cGlyPool2_13C01D00 
		dxdt(231, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C10D00_13C10D00  -  Vf_GlyPool2_CO2_13C10D00_13C0D00  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D00__13C010  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D00__13C110  -  Vr_Vin_Gly_13C10D00_13C10D00  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D00__13C100D000  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D00__13C101D000  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D01__13C100D001  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D01__13C101D001  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D10__13C100D010  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D10__13C101D010  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D11__13C100D011  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D11__13C101D011  +  Vr_GlyPool2_GlyPoolMitochon_13C10D00_13C10D00                                       +  Vf_Vin_Gly_13C10D00_13C10D00  +  Vf_SerPool2_GlyPool2_13C100D000_13C10D00  +  Vf_SerPool2_GlyPool2_13C101D000_13C10D00  +  Vf_SerPool2_GlyPool2_13C100D001_13C10D00  +  Vf_SerPool2_GlyPool2_13C101D001_13C10D00  +  Vf_SerPool2_GlyPool2_13C100D010_13C10D00  +  Vf_SerPool2_GlyPool2_13C101D010_13C10D00  +  Vf_SerPool2_GlyPool2_13C100D011_13C10D00  +  Vf_SerPool2_GlyPool2_13C101D011_13C10D00;          %  cGlyPool2_13C10D00 
		dxdt(232, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C11D00_13C11D00  -  Vf_GlyPool2_CO2_13C11D00_13C1D00  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D00__13C011  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D00__13C111  -  Vr_Vin_Gly_13C11D00_13C11D00  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D00__13C110D000  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D00__13C111D000  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D01__13C110D001  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D01__13C111D001  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D10__13C110D010  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D10__13C111D010  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D11__13C110D011  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D11__13C111D011  +  Vr_GlyPool2_GlyPoolMitochon_13C11D00_13C11D00                                       +  Vf_Vin_Gly_13C11D00_13C11D00  +  Vf_SerPool2_GlyPool2_13C110D000_13C11D00  +  Vf_SerPool2_GlyPool2_13C111D000_13C11D00  +  Vf_SerPool2_GlyPool2_13C110D001_13C11D00  +  Vf_SerPool2_GlyPool2_13C111D001_13C11D00  +  Vf_SerPool2_GlyPool2_13C110D010_13C11D00  +  Vf_SerPool2_GlyPool2_13C111D010_13C11D00  +  Vf_SerPool2_GlyPool2_13C110D011_13C11D00  +  Vf_SerPool2_GlyPool2_13C111D011_13C11D00;          %  cGlyPool2_13C11D00 
		dxdt(233, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C00D01_13C00D01  -  Vf_GlyPool2_CO2_13C00D01_13C0D01  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D01__13C000  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D01__13C100  -  Vr_Vin_Gly_13C00D01_13C00D01  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D00__13C000D000  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D00__13C001D000  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D01__13C000D001  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D01__13C001D001  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D10__13C000D010  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D10__13C001D010  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D11__13C000D011  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D11__13C001D011  +  Vr_GlyPool2_GlyPoolMitochon_13C00D01_13C00D01  +  Vr_GlyPool2_CO2_13C0D01_13C00D01  +  Vf_Vin_Gly_13C00D01_13C00D01;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C00D01 
		dxdt(234, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C01D01_13C01D01  -  Vf_GlyPool2_CO2_13C01D01_13C1D01  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D01__13C001  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D01__13C101  -  Vr_Vin_Gly_13C01D01_13C01D01  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D00__13C010D000  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D00__13C011D000  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D01__13C010D001  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D01__13C011D001  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D10__13C010D010  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D10__13C011D010  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D11__13C010D011  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D11__13C011D011  +  Vr_GlyPool2_GlyPoolMitochon_13C01D01_13C01D01  +  Vr_GlyPool2_CO2_13C1D01_13C01D01  +  Vf_Vin_Gly_13C01D01_13C01D01;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C01D01 
		dxdt(235, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C10D01_13C10D01  -  Vf_GlyPool2_CO2_13C10D01_13C0D01  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D01__13C010  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D01__13C110  -  Vr_Vin_Gly_13C10D01_13C10D01  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D00__13C100D000  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D00__13C101D000  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D01__13C100D001  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D01__13C101D001  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D10__13C100D010  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D10__13C101D010  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D11__13C100D011  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D11__13C101D011  +  Vr_GlyPool2_GlyPoolMitochon_13C10D01_13C10D01                                       +  Vf_Vin_Gly_13C10D01_13C10D01;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C10D01 
		dxdt(236, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C11D01_13C11D01  -  Vf_GlyPool2_CO2_13C11D01_13C1D01  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D01__13C011  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D01__13C111  -  Vr_Vin_Gly_13C11D01_13C11D01  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D00__13C110D000  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D00__13C111D000  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D01__13C110D001  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D01__13C111D001  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D10__13C110D010  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D10__13C111D010  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D11__13C110D011  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D11__13C111D011  +  Vr_GlyPool2_GlyPoolMitochon_13C11D01_13C11D01                                       +  Vf_Vin_Gly_13C11D01_13C11D01;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C11D01 
		dxdt(237, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C00D10_13C00D10  -  Vf_GlyPool2_CO2_13C00D10_13C0D10  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D10__13C000  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D10__13C100  -  Vr_Vin_Gly_13C00D10_13C00D10  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D00__13C000D100  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D00__13C001D100  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D01__13C000D101  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D01__13C001D101  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D10__13C000D110  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D10__13C001D110  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D11__13C000D111  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D11__13C001D111  +  Vr_GlyPool2_GlyPoolMitochon_13C00D10_13C00D10  +  Vr_GlyPool2_CO2_13C0D10_13C00D10  +  Vf_Vin_Gly_13C00D10_13C00D10  +  Vf_SerPool2_GlyPool2_13C000D100_13C00D10  +  Vf_SerPool2_GlyPool2_13C001D100_13C00D10  +  Vf_SerPool2_GlyPool2_13C000D101_13C00D10  +  Vf_SerPool2_GlyPool2_13C001D101_13C00D10  +  Vf_SerPool2_GlyPool2_13C000D110_13C00D10  +  Vf_SerPool2_GlyPool2_13C001D110_13C00D10  +  Vf_SerPool2_GlyPool2_13C000D111_13C00D10  +  Vf_SerPool2_GlyPool2_13C001D111_13C00D10;          %  cGlyPool2_13C00D10 
		dxdt(238, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C01D10_13C01D10  -  Vf_GlyPool2_CO2_13C01D10_13C1D10  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D10__13C001  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D10__13C101  -  Vr_Vin_Gly_13C01D10_13C01D10  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D00__13C010D100  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D00__13C011D100  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D01__13C010D101  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D01__13C011D101  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D10__13C010D110  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D10__13C011D110  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D11__13C010D111  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D11__13C011D111  +  Vr_GlyPool2_GlyPoolMitochon_13C01D10_13C01D10  +  Vr_GlyPool2_CO2_13C1D10_13C01D10  +  Vf_Vin_Gly_13C01D10_13C01D10  +  Vf_SerPool2_GlyPool2_13C010D100_13C01D10  +  Vf_SerPool2_GlyPool2_13C011D100_13C01D10  +  Vf_SerPool2_GlyPool2_13C010D101_13C01D10  +  Vf_SerPool2_GlyPool2_13C011D101_13C01D10  +  Vf_SerPool2_GlyPool2_13C010D110_13C01D10  +  Vf_SerPool2_GlyPool2_13C011D110_13C01D10  +  Vf_SerPool2_GlyPool2_13C010D111_13C01D10  +  Vf_SerPool2_GlyPool2_13C011D111_13C01D10;          %  cGlyPool2_13C01D10 
		dxdt(239, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C10D10_13C10D10  -  Vf_GlyPool2_CO2_13C10D10_13C0D10  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D10__13C010  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D10__13C110  -  Vr_Vin_Gly_13C10D10_13C10D10  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D00__13C100D100  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D00__13C101D100  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D01__13C100D101  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D01__13C101D101  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D10__13C100D110  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D10__13C101D110  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D11__13C100D111  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D11__13C101D111  +  Vr_GlyPool2_GlyPoolMitochon_13C10D10_13C10D10                                       +  Vf_Vin_Gly_13C10D10_13C10D10  +  Vf_SerPool2_GlyPool2_13C100D100_13C10D10  +  Vf_SerPool2_GlyPool2_13C101D100_13C10D10  +  Vf_SerPool2_GlyPool2_13C100D101_13C10D10  +  Vf_SerPool2_GlyPool2_13C101D101_13C10D10  +  Vf_SerPool2_GlyPool2_13C100D110_13C10D10  +  Vf_SerPool2_GlyPool2_13C101D110_13C10D10  +  Vf_SerPool2_GlyPool2_13C100D111_13C10D10  +  Vf_SerPool2_GlyPool2_13C101D111_13C10D10;          %  cGlyPool2_13C10D10 
		dxdt(240, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C11D10_13C11D10  -  Vf_GlyPool2_CO2_13C11D10_13C1D10  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D10__13C011  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D10__13C111  -  Vr_Vin_Gly_13C11D10_13C11D10  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D00__13C110D100  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D00__13C111D100  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D01__13C110D101  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D01__13C111D101  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D10__13C110D110  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D10__13C111D110  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D11__13C110D111  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D11__13C111D111  +  Vr_GlyPool2_GlyPoolMitochon_13C11D10_13C11D10                                       +  Vf_Vin_Gly_13C11D10_13C11D10  +  Vf_SerPool2_GlyPool2_13C110D100_13C11D10  +  Vf_SerPool2_GlyPool2_13C111D100_13C11D10  +  Vf_SerPool2_GlyPool2_13C110D101_13C11D10  +  Vf_SerPool2_GlyPool2_13C111D101_13C11D10  +  Vf_SerPool2_GlyPool2_13C110D110_13C11D10  +  Vf_SerPool2_GlyPool2_13C111D110_13C11D10  +  Vf_SerPool2_GlyPool2_13C110D111_13C11D10  +  Vf_SerPool2_GlyPool2_13C111D111_13C11D10;          %  cGlyPool2_13C11D10 
		dxdt(241, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C00D11_13C00D11  -  Vf_GlyPool2_CO2_13C00D11_13C0D11  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D11__13C000  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D11__13C100  -  Vr_Vin_Gly_13C00D11_13C00D11  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D00__13C000D100  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D00__13C001D100  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D01__13C000D101  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D01__13C001D101  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D10__13C000D110  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D10__13C001D110  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D11__13C000D111  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D11__13C001D111  +  Vr_GlyPool2_GlyPoolMitochon_13C00D11_13C00D11  +  Vr_GlyPool2_CO2_13C0D11_13C00D11  +  Vf_Vin_Gly_13C00D11_13C00D11;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C00D11 
		dxdt(242, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C01D11_13C01D11  -  Vf_GlyPool2_CO2_13C01D11_13C1D11  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D11__13C001  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D11__13C101  -  Vr_Vin_Gly_13C01D11_13C01D11  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D00__13C010D100  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D00__13C011D100  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D01__13C010D101  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D01__13C011D101  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D10__13C010D110  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D10__13C011D110  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D11__13C010D111  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D11__13C011D111  +  Vr_GlyPool2_GlyPoolMitochon_13C01D11_13C01D11  +  Vr_GlyPool2_CO2_13C1D11_13C01D11  +  Vf_Vin_Gly_13C01D11_13C01D11;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C01D11 
		dxdt(243, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C10D11_13C10D11  -  Vf_GlyPool2_CO2_13C10D11_13C0D11  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D11__13C010  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D11__13C110  -  Vr_Vin_Gly_13C10D11_13C10D11  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D00__13C100D100  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D00__13C101D100  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D01__13C100D101  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D01__13C101D101  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D10__13C100D110  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D10__13C101D110  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D11__13C100D111  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D11__13C101D111  +  Vr_GlyPool2_GlyPoolMitochon_13C10D11_13C10D11                                       +  Vf_Vin_Gly_13C10D11_13C10D11;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C10D11 
		dxdt(244, 1) =  - Vf_GlyPool2_GlyPoolMitochon_13C11D11_13C11D11  -  Vf_GlyPool2_CO2_13C11D11_13C1D11  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D11__13C011  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D11__13C111  -  Vr_Vin_Gly_13C11D11_13C11D11  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D00__13C110D100  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D00__13C111D100  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D01__13C110D101  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D01__13C111D101  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D10__13C110D110  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D10__13C111D110  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D11__13C110D111  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D11__13C111D111  +  Vr_GlyPool2_GlyPoolMitochon_13C11D11_13C11D11                                       +  Vf_Vin_Gly_13C11D11_13C11D11;                                                                                                                                                                                                                                                                                                                                                                                  %  cGlyPool2_13C11D11 

 
    
    % cSerPool3
    %  dxdt(1, 1) = (cSerPool3) = 
    % - Vf_SerPool3_Degradation
    % - Vf_SerPool3_GlyPool3
    % + Vr_SerPool3_GlyPool3
    % + Vf_SerPool3_Synthesis;
    % - Vr_Vin_Ser3 
    % + Vf_Vin_Ser3 
    
		dxdt(245, 1) =  -  Vf_SerPool3_Degradation_13C000D000_13C000D000  -  Vf_SerPool3_GlyPool3_13C000D000_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C0D00__13C000D000  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D00__13C000D000  +  Vf_SerPool3_Synthesis_13C000D000_13C000D000  -  Vr_Vin_Ser3_13C000D000_13C000D000  +  Vf_Vin_Ser3_13C000D000_13C000D000;  							%  cSerPool3_13C000D000 
		dxdt(246, 1) =  -  Vf_SerPool3_Degradation_13C001D000_13C001D000  -  Vf_SerPool3_GlyPool3_13C001D000_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D00__13C001D000  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D00__13C001D000                                                  -  Vr_Vin_Ser3_13C001D000_13C001D000  +  Vf_Vin_Ser3_13C001D000_13C001D000;                %  cSerPool3_13C001D000 
		dxdt(247, 1) =  -  Vf_SerPool3_Degradation_13C010D000_13C010D000  -  Vf_SerPool3_GlyPool3_13C010D000_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D00__13C010D000  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D00__13C010D000                                                  -  Vr_Vin_Ser3_13C010D000_13C010D000  +  Vf_Vin_Ser3_13C010D000_13C010D000;                %  cSerPool3_13C010D000 
		dxdt(248, 1) =  -  Vf_SerPool3_Degradation_13C011D000_13C011D000  -  Vf_SerPool3_GlyPool3_13C011D000_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D00__13C011D000  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D00__13C011D000                                                  -  Vr_Vin_Ser3_13C011D000_13C011D000  +  Vf_Vin_Ser3_13C011D000_13C011D000;                %  cSerPool3_13C011D000 
		dxdt(249, 1) =  -  Vf_SerPool3_Degradation_13C100D000_13C100D000  -  Vf_SerPool3_GlyPool3_13C100D000_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D00__13C100D000  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D00__13C100D000                                                  -  Vr_Vin_Ser3_13C100D000_13C100D000  +  Vf_Vin_Ser3_13C100D000_13C100D000;                %  cSerPool3_13C100D000 
		dxdt(250, 1) =  -  Vf_SerPool3_Degradation_13C101D000_13C101D000  -  Vf_SerPool3_GlyPool3_13C101D000_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D00__13C101D000  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D00__13C101D000                                                  -  Vr_Vin_Ser3_13C101D000_13C101D000  +  Vf_Vin_Ser3_13C101D000_13C101D000;                %  cSerPool3_13C101D000 
		dxdt(251, 1) =  -  Vf_SerPool3_Degradation_13C110D000_13C110D000  -  Vf_SerPool3_GlyPool3_13C110D000_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D00__13C110D000  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D00__13C110D000                                                  -  Vr_Vin_Ser3_13C110D000_13C110D000  +  Vf_Vin_Ser3_13C110D000_13C110D000;                %  cSerPool3_13C110D000 
		dxdt(252, 1) =  -  Vf_SerPool3_Degradation_13C111D000_13C111D000  -  Vf_SerPool3_GlyPool3_13C111D000_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D00__13C111D000  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D00__13C111D000                                                  -  Vr_Vin_Ser3_13C111D000_13C111D000  +  Vf_Vin_Ser3_13C111D000_13C111D000;                %  cSerPool3_13C111D000 
		dxdt(253, 1) =  -  Vf_SerPool3_Degradation_13C000D001_13C000D001  -  Vf_SerPool3_GlyPool3_13C000D001_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C0D01__13C000D001  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D01__13C000D001                                                  -  Vr_Vin_Ser3_13C000D001_13C000D001  +  Vf_Vin_Ser3_13C000D001_13C000D001;                %  cSerPool3_13C000D001 
		dxdt(254, 1) =  -  Vf_SerPool3_Degradation_13C001D001_13C001D001  -  Vf_SerPool3_GlyPool3_13C001D001_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D01__13C001D001  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D01__13C001D001                                                  -  Vr_Vin_Ser3_13C001D001_13C001D001  +  Vf_Vin_Ser3_13C001D001_13C001D001;                %  cSerPool3_13C001D001 
		dxdt(255, 1) =  -  Vf_SerPool3_Degradation_13C010D001_13C010D001  -  Vf_SerPool3_GlyPool3_13C010D001_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D01__13C010D001  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D01__13C010D001                                                  -  Vr_Vin_Ser3_13C010D001_13C010D001  +  Vf_Vin_Ser3_13C010D001_13C010D001;                %  cSerPool3_13C010D001 
		dxdt(256, 1) =  -  Vf_SerPool3_Degradation_13C011D001_13C011D001  -  Vf_SerPool3_GlyPool3_13C011D001_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D01__13C011D001  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D01__13C011D001                                                  -  Vr_Vin_Ser3_13C011D001_13C011D001  +  Vf_Vin_Ser3_13C011D001_13C011D001;                %  cSerPool3_13C011D001 
		dxdt(257, 1) =  -  Vf_SerPool3_Degradation_13C100D001_13C100D001  -  Vf_SerPool3_GlyPool3_13C100D001_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D01__13C100D001  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D01__13C100D001                                                  -  Vr_Vin_Ser3_13C100D001_13C100D001  +  Vf_Vin_Ser3_13C100D001_13C100D001;                %  cSerPool3_13C100D001 
		dxdt(258, 1) =  -  Vf_SerPool3_Degradation_13C101D001_13C101D001  -  Vf_SerPool3_GlyPool3_13C101D001_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D01__13C101D001  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D01__13C101D001                                                  -  Vr_Vin_Ser3_13C101D001_13C101D001  +  Vf_Vin_Ser3_13C101D001_13C101D001;                %  cSerPool3_13C101D001 
		dxdt(259, 1) =  -  Vf_SerPool3_Degradation_13C110D001_13C110D001  -  Vf_SerPool3_GlyPool3_13C110D001_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D01__13C110D001  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D01__13C110D001                                                  -  Vr_Vin_Ser3_13C110D001_13C110D001  +  Vf_Vin_Ser3_13C110D001_13C110D001;                %  cSerPool3_13C110D001 
		dxdt(260, 1) =  -  Vf_SerPool3_Degradation_13C111D001_13C111D001  -  Vf_SerPool3_GlyPool3_13C111D001_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D01__13C111D001  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D01__13C111D001                                                  -  Vr_Vin_Ser3_13C111D001_13C111D001  +  Vf_Vin_Ser3_13C111D001_13C111D001;                %  cSerPool3_13C111D001 
		dxdt(261, 1) =  -  Vf_SerPool3_Degradation_13C000D010_13C000D010  -  Vf_SerPool3_GlyPool3_13C000D010_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C0D10__13C000D010  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D10__13C000D010                                                  -  Vr_Vin_Ser3_13C000D010_13C000D010  +  Vf_Vin_Ser3_13C000D010_13C000D010;                %  cSerPool3_13C000D010 
		dxdt(262, 1) =  -  Vf_SerPool3_Degradation_13C001D010_13C001D010  -  Vf_SerPool3_GlyPool3_13C001D010_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D10__13C001D010  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D10__13C001D010                                                  -  Vr_Vin_Ser3_13C001D010_13C001D010  +  Vf_Vin_Ser3_13C001D010_13C001D010;                %  cSerPool3_13C001D010 
		dxdt(263, 1) =  -  Vf_SerPool3_Degradation_13C010D010_13C010D010  -  Vf_SerPool3_GlyPool3_13C010D010_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D10__13C010D010  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D10__13C010D010                                                  -  Vr_Vin_Ser3_13C010D010_13C010D010  +  Vf_Vin_Ser3_13C010D010_13C010D010;                %  cSerPool3_13C010D010 
		dxdt(264, 1) =  -  Vf_SerPool3_Degradation_13C011D010_13C011D010  -  Vf_SerPool3_GlyPool3_13C011D010_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D10__13C011D010  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D10__13C011D010                                                  -  Vr_Vin_Ser3_13C011D010_13C011D010  +  Vf_Vin_Ser3_13C011D010_13C011D010;                %  cSerPool3_13C011D010 
		dxdt(265, 1) =  -  Vf_SerPool3_Degradation_13C100D010_13C100D010  -  Vf_SerPool3_GlyPool3_13C100D010_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D10__13C100D010  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D10__13C100D010                                                  -  Vr_Vin_Ser3_13C100D010_13C100D010  +  Vf_Vin_Ser3_13C100D010_13C100D010;                %  cSerPool3_13C100D010 
		dxdt(266, 1) =  -  Vf_SerPool3_Degradation_13C101D010_13C101D010  -  Vf_SerPool3_GlyPool3_13C101D010_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D10__13C101D010  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D10__13C101D010                                                  -  Vr_Vin_Ser3_13C101D010_13C101D010  +  Vf_Vin_Ser3_13C101D010_13C101D010;                %  cSerPool3_13C101D010 
		dxdt(267, 1) =  -  Vf_SerPool3_Degradation_13C110D010_13C110D010  -  Vf_SerPool3_GlyPool3_13C110D010_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D10__13C110D010  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D10__13C110D010                                                  -  Vr_Vin_Ser3_13C110D010_13C110D010  +  Vf_Vin_Ser3_13C110D010_13C110D010;                %  cSerPool3_13C110D010 
		dxdt(268, 1) =  -  Vf_SerPool3_Degradation_13C111D010_13C111D010  -  Vf_SerPool3_GlyPool3_13C111D010_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D10__13C111D010  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D10__13C111D010                                                  -  Vr_Vin_Ser3_13C111D010_13C111D010  +  Vf_Vin_Ser3_13C111D010_13C111D010;                %  cSerPool3_13C111D010 
		dxdt(269, 1) =  -  Vf_SerPool3_Degradation_13C000D011_13C000D011  -  Vf_SerPool3_GlyPool3_13C000D011_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C0D11__13C000D011  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D11__13C000D011                                                  -  Vr_Vin_Ser3_13C000D011_13C000D011  +  Vf_Vin_Ser3_13C000D011_13C000D011;                %  cSerPool3_13C000D011 
		dxdt(270, 1) =  -  Vf_SerPool3_Degradation_13C001D011_13C001D011  -  Vf_SerPool3_GlyPool3_13C001D011_13C00D00  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D11__13C001D011  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D11__13C001D011                                                  -  Vr_Vin_Ser3_13C001D011_13C001D011  +  Vf_Vin_Ser3_13C001D011_13C001D011;                %  cSerPool3_13C001D011 
		dxdt(271, 1) =  -  Vf_SerPool3_Degradation_13C010D011_13C010D011  -  Vf_SerPool3_GlyPool3_13C010D011_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D11__13C010D011  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D11__13C010D011                                                  -  Vr_Vin_Ser3_13C010D011_13C010D011  +  Vf_Vin_Ser3_13C010D011_13C010D011;                %  cSerPool3_13C010D011 
		dxdt(272, 1) =  -  Vf_SerPool3_Degradation_13C011D011_13C011D011  -  Vf_SerPool3_GlyPool3_13C011D011_13C01D00  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D11__13C011D011  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D11__13C011D011                                                  -  Vr_Vin_Ser3_13C011D011_13C011D011  +  Vf_Vin_Ser3_13C011D011_13C011D011;                %  cSerPool3_13C011D011 
		dxdt(273, 1) =  -  Vf_SerPool3_Degradation_13C100D011_13C100D011  -  Vf_SerPool3_GlyPool3_13C100D011_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D11__13C100D011  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D11__13C100D011                                                  -  Vr_Vin_Ser3_13C100D011_13C100D011  +  Vf_Vin_Ser3_13C100D011_13C100D011;                %  cSerPool3_13C100D011 
		dxdt(274, 1) =  -  Vf_SerPool3_Degradation_13C101D011_13C101D011  -  Vf_SerPool3_GlyPool3_13C101D011_13C10D00  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D11__13C101D011  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D11__13C101D011                                                  -  Vr_Vin_Ser3_13C101D011_13C101D011  +  Vf_Vin_Ser3_13C101D011_13C101D011;                %  cSerPool3_13C101D011 
		dxdt(275, 1) =  -  Vf_SerPool3_Degradation_13C110D011_13C110D011  -  Vf_SerPool3_GlyPool3_13C110D011_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D11__13C110D011  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D11__13C110D011                                                  -  Vr_Vin_Ser3_13C110D011_13C110D011  +  Vf_Vin_Ser3_13C110D011_13C110D011;                %  cSerPool3_13C110D011 
		dxdt(276, 1) =  -  Vf_SerPool3_Degradation_13C111D011_13C111D011  -  Vf_SerPool3_GlyPool3_13C111D011_13C11D00  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D11__13C111D011  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D11__13C111D011                                                  -  Vr_Vin_Ser3_13C111D011_13C111D011  +  Vf_Vin_Ser3_13C111D011_13C111D011;                %  cSerPool3_13C111D011 
		dxdt(277, 1) =  -  Vf_SerPool3_Degradation_13C000D100_13C000D100  -  Vf_SerPool3_GlyPool3_13C000D100_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D00__13C000D100  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D00__13C000D100                                                  -  Vr_Vin_Ser3_13C000D100_13C000D100  +  Vf_Vin_Ser3_13C000D100_13C000D100;                %  cSerPool3_13C000D100 
		dxdt(278, 1) =  -  Vf_SerPool3_Degradation_13C001D100_13C001D100  -  Vf_SerPool3_GlyPool3_13C001D100_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D00__13C001D100  +  Vr_SerPool3_GlyPool3_13C00D11_13C1D00__13C001D100                                                  -  Vr_Vin_Ser3_13C001D100_13C001D100  +  Vf_Vin_Ser3_13C001D100_13C001D100;                %  cSerPool3_13C001D100 
		dxdt(279, 1) =  -  Vf_SerPool3_Degradation_13C010D100_13C010D100  -  Vf_SerPool3_GlyPool3_13C010D100_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D00__13C010D100  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D00__13C010D100                                                  -  Vr_Vin_Ser3_13C010D100_13C010D100  +  Vf_Vin_Ser3_13C010D100_13C010D100;                %  cSerPool3_13C010D100 
		dxdt(280, 1) =  -  Vf_SerPool3_Degradation_13C011D100_13C011D100  -  Vf_SerPool3_GlyPool3_13C011D100_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D00__13C011D100  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D00__13C011D100                                                  -  Vr_Vin_Ser3_13C011D100_13C011D100  +  Vf_Vin_Ser3_13C011D100_13C011D100;                %  cSerPool3_13C011D100 
		dxdt(281, 1) =  -  Vf_SerPool3_Degradation_13C100D100_13C100D100  -  Vf_SerPool3_GlyPool3_13C100D100_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D00__13C100D100  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D00__13C100D100                                                  -  Vr_Vin_Ser3_13C100D100_13C100D100  +  Vf_Vin_Ser3_13C100D100_13C100D100;                %  cSerPool3_13C100D100 
		dxdt(282, 1) =  -  Vf_SerPool3_Degradation_13C101D100_13C101D100  -  Vf_SerPool3_GlyPool3_13C101D100_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D00__13C101D100  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D00__13C101D100                                                  -  Vr_Vin_Ser3_13C101D100_13C101D100  +  Vf_Vin_Ser3_13C101D100_13C101D100;                %  cSerPool3_13C101D100 
		dxdt(283, 1) =  -  Vf_SerPool3_Degradation_13C110D100_13C110D100  -  Vf_SerPool3_GlyPool3_13C110D100_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D00__13C110D100  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D00__13C110D100                                                  -  Vr_Vin_Ser3_13C110D100_13C110D100  +  Vf_Vin_Ser3_13C110D100_13C110D100;                %  cSerPool3_13C110D100 
		dxdt(284, 1) =  -  Vf_SerPool3_Degradation_13C111D100_13C111D100  -  Vf_SerPool3_GlyPool3_13C111D100_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D00__13C111D100  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D00__13C111D100                                                  -  Vr_Vin_Ser3_13C111D100_13C111D100  +  Vf_Vin_Ser3_13C111D100_13C111D100;                %  cSerPool3_13C111D100 
		dxdt(285, 1) =  -  Vf_SerPool3_Degradation_13C000D101_13C000D101  -  Vf_SerPool3_GlyPool3_13C000D101_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D01__13C000D101  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D01__13C000D101                                                  -  Vr_Vin_Ser3_13C000D101_13C000D101  +  Vf_Vin_Ser3_13C000D101_13C000D101;                %  cSerPool3_13C000D101 
		dxdt(286, 1) =  -  Vf_SerPool3_Degradation_13C001D101_13C001D101  -  Vf_SerPool3_GlyPool3_13C001D101_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D01__13C001D101  +  Vr_SerPool3_GlyPool3_13C00D11_13C1D01__13C001D101                                                  -  Vr_Vin_Ser3_13C001D101_13C001D101  +  Vf_Vin_Ser3_13C001D101_13C001D101;                %  cSerPool3_13C001D101 
		dxdt(287, 1) =  -  Vf_SerPool3_Degradation_13C010D101_13C010D101  -  Vf_SerPool3_GlyPool3_13C010D101_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D01__13C010D101  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D01__13C010D101                                                  -  Vr_Vin_Ser3_13C010D101_13C010D101  +  Vf_Vin_Ser3_13C010D101_13C010D101;                %  cSerPool3_13C010D101 
		dxdt(288, 1) =  -  Vf_SerPool3_Degradation_13C011D101_13C011D101  -  Vf_SerPool3_GlyPool3_13C011D101_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D01__13C011D101  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D01__13C011D101                                                  -  Vr_Vin_Ser3_13C011D101_13C011D101  +  Vf_Vin_Ser3_13C011D101_13C011D101;                %  cSerPool3_13C011D101 
		dxdt(289, 1) =  -  Vf_SerPool3_Degradation_13C100D101_13C100D101  -  Vf_SerPool3_GlyPool3_13C100D101_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D01__13C100D101  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D01__13C100D101                                                  -  Vr_Vin_Ser3_13C100D101_13C100D101  +  Vf_Vin_Ser3_13C100D101_13C100D101;                %  cSerPool3_13C100D101 
		dxdt(290, 1) =  -  Vf_SerPool3_Degradation_13C101D101_13C101D101  -  Vf_SerPool3_GlyPool3_13C101D101_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D01__13C101D101  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D01__13C101D101                                                  -  Vr_Vin_Ser3_13C101D101_13C101D101  +  Vf_Vin_Ser3_13C101D101_13C101D101;                %  cSerPool3_13C101D101 
		dxdt(291, 1) =  -  Vf_SerPool3_Degradation_13C110D101_13C110D101  -  Vf_SerPool3_GlyPool3_13C110D101_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D01__13C110D101  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D01__13C110D101                                                  -  Vr_Vin_Ser3_13C110D101_13C110D101  +  Vf_Vin_Ser3_13C110D101_13C110D101;                %  cSerPool3_13C110D101 
		dxdt(292, 1) =  -  Vf_SerPool3_Degradation_13C111D101_13C111D101  -  Vf_SerPool3_GlyPool3_13C111D101_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D01__13C111D101  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D01__13C111D101                                                  -  Vr_Vin_Ser3_13C111D101_13C111D101  +  Vf_Vin_Ser3_13C111D101_13C111D101;                %  cSerPool3_13C111D101 
		dxdt(293, 1) =  -  Vf_SerPool3_Degradation_13C000D110_13C000D110  -  Vf_SerPool3_GlyPool3_13C000D110_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D10__13C000D110  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D10__13C000D110                                                  -  Vr_Vin_Ser3_13C000D110_13C000D110  +  Vf_Vin_Ser3_13C000D110_13C000D110;                %  cSerPool3_13C000D110 
		dxdt(294, 1) =  -  Vf_SerPool3_Degradation_13C001D110_13C001D110  -  Vf_SerPool3_GlyPool3_13C001D110_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D10__13C001D110  +  Vr_SerPool3_GlyPool3_13C00D11_13C1D10__13C001D110                                                  -  Vr_Vin_Ser3_13C001D110_13C001D110  +  Vf_Vin_Ser3_13C001D110_13C001D110;                %  cSerPool3_13C001D110 
		dxdt(295, 1) =  -  Vf_SerPool3_Degradation_13C010D110_13C010D110  -  Vf_SerPool3_GlyPool3_13C010D110_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D10__13C010D110  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D10__13C010D110                                                  -  Vr_Vin_Ser3_13C010D110_13C010D110  +  Vf_Vin_Ser3_13C010D110_13C010D110;                %  cSerPool3_13C010D110 
		dxdt(296, 1) =  -  Vf_SerPool3_Degradation_13C011D110_13C011D110  -  Vf_SerPool3_GlyPool3_13C011D110_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D10__13C011D110  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D10__13C011D110                                                  -  Vr_Vin_Ser3_13C011D110_13C011D110  +  Vf_Vin_Ser3_13C011D110_13C011D110;                %  cSerPool3_13C011D110 
		dxdt(297, 1) =  -  Vf_SerPool3_Degradation_13C100D110_13C100D110  -  Vf_SerPool3_GlyPool3_13C100D110_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D10__13C100D110  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D10__13C100D110                                                  -  Vr_Vin_Ser3_13C100D110_13C100D110  +  Vf_Vin_Ser3_13C100D110_13C100D110;                %  cSerPool3_13C100D110 
		dxdt(298, 1) =  -  Vf_SerPool3_Degradation_13C101D110_13C101D110  -  Vf_SerPool3_GlyPool3_13C101D110_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D10__13C101D110  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D10__13C101D110                                                  -  Vr_Vin_Ser3_13C101D110_13C101D110  +  Vf_Vin_Ser3_13C101D110_13C101D110;                %  cSerPool3_13C101D110 
		dxdt(299, 1) =  -  Vf_SerPool3_Degradation_13C110D110_13C110D110  -  Vf_SerPool3_GlyPool3_13C110D110_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D10__13C110D110  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D10__13C110D110                                                  -  Vr_Vin_Ser3_13C110D110_13C110D110  +  Vf_Vin_Ser3_13C110D110_13C110D110;                %  cSerPool3_13C110D110 
		dxdt(300, 1) =  -  Vf_SerPool3_Degradation_13C111D110_13C111D110  -  Vf_SerPool3_GlyPool3_13C111D110_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D10__13C111D110  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D10__13C111D110                                                  -  Vr_Vin_Ser3_13C111D110_13C111D110  +  Vf_Vin_Ser3_13C111D110_13C111D110;                %  cSerPool3_13C111D110 
		dxdt(301, 1) =  -  Vf_SerPool3_Degradation_13C000D111_13C000D111  -  Vf_SerPool3_GlyPool3_13C000D111_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D11__13C000D111  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D11__13C000D111                                                  -  Vr_Vin_Ser3_13C000D111_13C000D111  +  Vf_Vin_Ser3_13C000D111_13C000D111;                %  cSerPool3_13C000D111 
		dxdt(302, 1) =  -  Vf_SerPool3_Degradation_13C001D111_13C001D111  -  Vf_SerPool3_GlyPool3_13C001D111_13C00D10  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D11__13C001D111  +  Vr_SerPool3_GlyPool3_13C00D11_13C1D11__13C001D111                                                  -  Vr_Vin_Ser3_13C001D111_13C001D111  +  Vf_Vin_Ser3_13C001D111_13C001D111;                %  cSerPool3_13C001D111 
		dxdt(303, 1) =  -  Vf_SerPool3_Degradation_13C010D111_13C010D111  -  Vf_SerPool3_GlyPool3_13C010D111_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D11__13C010D111  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D11__13C010D111                                                  -  Vr_Vin_Ser3_13C010D111_13C010D111  +  Vf_Vin_Ser3_13C010D111_13C010D111;                %  cSerPool3_13C010D111 
		dxdt(304, 1) =  -  Vf_SerPool3_Degradation_13C011D111_13C011D111  -  Vf_SerPool3_GlyPool3_13C011D111_13C01D10  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D11__13C011D111  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D11__13C011D111                                                  -  Vr_Vin_Ser3_13C011D111_13C011D111  +  Vf_Vin_Ser3_13C011D111_13C011D111;                %  cSerPool3_13C011D111 
		dxdt(305, 1) =  -  Vf_SerPool3_Degradation_13C100D111_13C100D111  -  Vf_SerPool3_GlyPool3_13C100D111_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D11__13C100D111  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D11__13C100D111                                                  -  Vr_Vin_Ser3_13C100D111_13C100D111  +  Vf_Vin_Ser3_13C100D111_13C100D111;                %  cSerPool3_13C100D111 
		dxdt(306, 1) =  -  Vf_SerPool3_Degradation_13C101D111_13C101D111  -  Vf_SerPool3_GlyPool3_13C101D111_13C10D10  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D11__13C101D111  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D11__13C101D111                                                  -  Vr_Vin_Ser3_13C101D111_13C101D111  +  Vf_Vin_Ser3_13C101D111_13C101D111;                %  cSerPool3_13C101D111 
		dxdt(307, 1) =  -  Vf_SerPool3_Degradation_13C110D111_13C110D111  -  Vf_SerPool3_GlyPool3_13C110D111_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D11__13C110D111  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D11__13C110D111                                                  -  Vr_Vin_Ser3_13C110D111_13C110D111  +  Vf_Vin_Ser3_13C110D111_13C110D111;                %  cSerPool3_13C110D111 
		dxdt(308, 1) =  -  Vf_SerPool3_Degradation_13C111D111_13C111D111  -  Vf_SerPool3_GlyPool3_13C111D111_13C11D10  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D11__13C111D111  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D11__13C111D111                                                  -  Vr_Vin_Ser3_13C111D111_13C111D111  +  Vf_Vin_Ser3_13C111D111_13C111D111;                %  cSerPool3_13C111D111 
           


    % cGlyPool3
    %  dxdt(1, 1) = (cGlyPool3) = 
    % - Vf_GlyPool3_Degradation     
    % - Vf_GlyPool3_CO2     
    % - Vf_PRPP_GAR_GlyPool3     
    % - Vr_SerPool3_GlyPool3     
    % + Vr_GlyPool3_CO2 
    % + Vf_SerPool3_GlyPool3 
    % + Vf_GlyPool3_Synthesis;
    % - Vr_Vin_Gly3     
    % + Vf_Vin_Gly3 
    

		dxdt(309, 1) =  - Vf_GlyPool3_Degradation_13C00D00_13C00D00  -  Vf_GlyPool3_CO2_13C00D00_13C0D00  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D00__13C000  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D00__13C100  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D00__13C000D000  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D00__13C001D000  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D01__13C000D001  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D01__13C001D001  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D10__13C000D010  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D10__13C001D010  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D11__13C000D011  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D11__13C001D011  +  Vr_GlyPool3_CO2_13C0D00_13C00D00  +  Vf_SerPool3_GlyPool3_13C000D000_13C00D00  +  Vf_SerPool3_GlyPool3_13C001D000_13C00D00  +  Vf_SerPool3_GlyPool3_13C000D001_13C00D00  +  Vf_SerPool3_GlyPool3_13C001D001_13C00D00  +  Vf_SerPool3_GlyPool3_13C000D010_13C00D00  +  Vf_SerPool3_GlyPool3_13C001D010_13C00D00  +  Vf_SerPool3_GlyPool3_13C000D011_13C00D00  +  Vf_SerPool3_GlyPool3_13C001D011_13C00D00  +  Vf_GlyPool3_Synthesis_13C00D00_13C00D00  -  Vr_Vin_Gly3_13C00D00_13C00D00  +  Vf_Vin_Gly3_13C00D00_13C00D00;							%  cGlyPool3_13C00D00 
		dxdt(310, 1) =  - Vf_GlyPool3_Degradation_13C01D00_13C01D00  -  Vf_GlyPool3_CO2_13C01D00_13C1D00  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D00__13C001  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D00__13C101  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D00__13C010D000  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D00__13C011D000  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D01__13C010D001  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D01__13C011D001  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D10__13C010D010  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D10__13C011D010  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D11__13C010D011  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D11__13C011D011  +  Vr_GlyPool3_CO2_13C1D00_13C01D00  +  Vf_SerPool3_GlyPool3_13C010D000_13C01D00  +  Vf_SerPool3_GlyPool3_13C011D000_13C01D00  +  Vf_SerPool3_GlyPool3_13C010D001_13C01D00  +  Vf_SerPool3_GlyPool3_13C011D001_13C01D00  +  Vf_SerPool3_GlyPool3_13C010D010_13C01D00  +  Vf_SerPool3_GlyPool3_13C011D010_13C01D00  +  Vf_SerPool3_GlyPool3_13C010D011_13C01D00  +  Vf_SerPool3_GlyPool3_13C011D011_13C01D00                                              -  Vr_Vin_Gly3_13C01D00_13C01D00  +  Vf_Vin_Gly3_13C01D00_13C01D00;              %  cGlyPool3_13C01D00 
		dxdt(311, 1) =  - Vf_GlyPool3_Degradation_13C10D00_13C10D00  -  Vf_GlyPool3_CO2_13C10D00_13C0D00  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D00__13C010  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D00__13C110  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D00__13C100D000  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D00__13C101D000  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D01__13C100D001  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D01__13C101D001  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D10__13C100D010  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D10__13C101D010  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D11__13C100D011  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D11__13C101D011                                       +  Vf_SerPool3_GlyPool3_13C100D000_13C10D00  +  Vf_SerPool3_GlyPool3_13C101D000_13C10D00  +  Vf_SerPool3_GlyPool3_13C100D001_13C10D00  +  Vf_SerPool3_GlyPool3_13C101D001_13C10D00  +  Vf_SerPool3_GlyPool3_13C100D010_13C10D00  +  Vf_SerPool3_GlyPool3_13C101D010_13C10D00  +  Vf_SerPool3_GlyPool3_13C100D011_13C10D00  +  Vf_SerPool3_GlyPool3_13C101D011_13C10D00                                              -  Vr_Vin_Gly3_13C10D00_13C10D00  +  Vf_Vin_Gly3_13C10D00_13C10D00;              %  cGlyPool3_13C10D00 
		dxdt(312, 1) =  - Vf_GlyPool3_Degradation_13C11D00_13C11D00  -  Vf_GlyPool3_CO2_13C11D00_13C1D00  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D00__13C011  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D00__13C111  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D00__13C110D000  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D00__13C111D000  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D01__13C110D001  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D01__13C111D001  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D10__13C110D010  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D10__13C111D010  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D11__13C110D011  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D11__13C111D011                                       +  Vf_SerPool3_GlyPool3_13C110D000_13C11D00  +  Vf_SerPool3_GlyPool3_13C111D000_13C11D00  +  Vf_SerPool3_GlyPool3_13C110D001_13C11D00  +  Vf_SerPool3_GlyPool3_13C111D001_13C11D00  +  Vf_SerPool3_GlyPool3_13C110D010_13C11D00  +  Vf_SerPool3_GlyPool3_13C111D010_13C11D00  +  Vf_SerPool3_GlyPool3_13C110D011_13C11D00  +  Vf_SerPool3_GlyPool3_13C111D011_13C11D00                                              -  Vr_Vin_Gly3_13C11D00_13C11D00  +  Vf_Vin_Gly3_13C11D00_13C11D00;              %  cGlyPool3_13C11D00 
		dxdt(313, 1) =  - Vf_GlyPool3_Degradation_13C00D01_13C00D01  -  Vf_GlyPool3_CO2_13C00D01_13C0D01  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D01__13C000  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D01__13C100  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D00__13C000D000  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D00__13C001D000  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D01__13C000D001  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D01__13C001D001  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D10__13C000D010  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D10__13C001D010  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D11__13C000D011  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D11__13C001D011  +  Vr_GlyPool3_CO2_13C0D01_13C00D01                                                                                                                                                                                                                                                                                                                                                                                                                      -  Vr_Vin_Gly3_13C00D01_13C00D01  +  Vf_Vin_Gly3_13C00D01_13C00D01;              %  cGlyPool3_13C00D01 
		dxdt(314, 1) =  - Vf_GlyPool3_Degradation_13C01D01_13C01D01  -  Vf_GlyPool3_CO2_13C01D01_13C1D01  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D01__13C001  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D01__13C101  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D00__13C010D000  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D00__13C011D000  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D01__13C010D001  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D01__13C011D001  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D10__13C010D010  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D10__13C011D010  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D11__13C010D011  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D11__13C011D011  +  Vr_GlyPool3_CO2_13C1D01_13C01D01                                                                                                                                                                                                                                                                                                                                                                                                                      -  Vr_Vin_Gly3_13C01D01_13C01D01  +  Vf_Vin_Gly3_13C01D01_13C01D01;              %  cGlyPool3_13C01D01 
		dxdt(315, 1) =  - Vf_GlyPool3_Degradation_13C10D01_13C10D01  -  Vf_GlyPool3_CO2_13C10D01_13C0D01  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D01__13C010  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D01__13C110  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D00__13C100D000  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D00__13C101D000  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D01__13C100D001  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D01__13C101D001  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D10__13C100D010  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D10__13C101D010  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D11__13C100D011  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D11__13C101D011                                                                                                                                                                                                                                                                                                                                                                                                                                                           -  Vr_Vin_Gly3_13C10D01_13C10D01  +  Vf_Vin_Gly3_13C10D01_13C10D01;              %  cGlyPool3_13C10D01 
		dxdt(316, 1) =  - Vf_GlyPool3_Degradation_13C11D01_13C11D01  -  Vf_GlyPool3_CO2_13C11D01_13C1D01  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D01__13C011  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D01__13C111  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D00__13C110D000  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D00__13C111D000  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D01__13C110D001  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D01__13C111D001  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D10__13C110D010  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D10__13C111D010  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D11__13C110D011  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D11__13C111D011                                                                                                                                                                                                                                                                                                                                                                                                                                                           -  Vr_Vin_Gly3_13C11D01_13C11D01  +  Vf_Vin_Gly3_13C11D01_13C11D01;              %  cGlyPool3_13C11D01 
		dxdt(317, 1) =  - Vf_GlyPool3_Degradation_13C00D10_13C00D10  -  Vf_GlyPool3_CO2_13C00D10_13C0D10  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D10__13C000  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D10__13C100  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D00__13C000D100  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D00__13C001D100  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D01__13C000D101  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D01__13C001D101  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D10__13C000D110  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D10__13C001D110  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D11__13C000D111  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D11__13C001D111  +  Vr_GlyPool3_CO2_13C0D10_13C00D10  +  Vf_SerPool3_GlyPool3_13C000D100_13C00D10  +  Vf_SerPool3_GlyPool3_13C001D100_13C00D10  +  Vf_SerPool3_GlyPool3_13C000D101_13C00D10  +  Vf_SerPool3_GlyPool3_13C001D101_13C00D10  +  Vf_SerPool3_GlyPool3_13C000D110_13C00D10  +  Vf_SerPool3_GlyPool3_13C001D110_13C00D10  +  Vf_SerPool3_GlyPool3_13C000D111_13C00D10  +  Vf_SerPool3_GlyPool3_13C001D111_13C00D10                                              -  Vr_Vin_Gly3_13C00D10_13C00D10  +  Vf_Vin_Gly3_13C00D10_13C00D10;              %  cGlyPool3_13C00D10 
		dxdt(318, 1) =  - Vf_GlyPool3_Degradation_13C01D10_13C01D10  -  Vf_GlyPool3_CO2_13C01D10_13C1D10  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D10__13C001  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D10__13C101  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D00__13C010D100  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D00__13C011D100  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D01__13C010D101  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D01__13C011D101  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D10__13C010D110  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D10__13C011D110  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D11__13C010D111  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D11__13C011D111  +  Vr_GlyPool3_CO2_13C1D10_13C01D10  +  Vf_SerPool3_GlyPool3_13C010D100_13C01D10  +  Vf_SerPool3_GlyPool3_13C011D100_13C01D10  +  Vf_SerPool3_GlyPool3_13C010D101_13C01D10  +  Vf_SerPool3_GlyPool3_13C011D101_13C01D10  +  Vf_SerPool3_GlyPool3_13C010D110_13C01D10  +  Vf_SerPool3_GlyPool3_13C011D110_13C01D10  +  Vf_SerPool3_GlyPool3_13C010D111_13C01D10  +  Vf_SerPool3_GlyPool3_13C011D111_13C01D10                                              -  Vr_Vin_Gly3_13C01D10_13C01D10  +  Vf_Vin_Gly3_13C01D10_13C01D10;              %  cGlyPool3_13C01D10 
		dxdt(319, 1) =  - Vf_GlyPool3_Degradation_13C10D10_13C10D10  -  Vf_GlyPool3_CO2_13C10D10_13C0D10  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D10__13C010  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D10__13C110  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D00__13C100D100  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D00__13C101D100  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D01__13C100D101  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D01__13C101D101  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D10__13C100D110  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D10__13C101D110  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D11__13C100D111  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D11__13C101D111                                       +  Vf_SerPool3_GlyPool3_13C100D100_13C10D10  +  Vf_SerPool3_GlyPool3_13C101D100_13C10D10  +  Vf_SerPool3_GlyPool3_13C100D101_13C10D10  +  Vf_SerPool3_GlyPool3_13C101D101_13C10D10  +  Vf_SerPool3_GlyPool3_13C100D110_13C10D10  +  Vf_SerPool3_GlyPool3_13C101D110_13C10D10  +  Vf_SerPool3_GlyPool3_13C100D111_13C10D10  +  Vf_SerPool3_GlyPool3_13C101D111_13C10D10                                              -  Vr_Vin_Gly3_13C10D10_13C10D10  +  Vf_Vin_Gly3_13C10D10_13C10D10;              %  cGlyPool3_13C10D10 
		dxdt(320, 1) =  - Vf_GlyPool3_Degradation_13C11D10_13C11D10  -  Vf_GlyPool3_CO2_13C11D10_13C1D10  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D10__13C011  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D10__13C111  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D00__13C110D100  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D00__13C111D100  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D01__13C110D101  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D01__13C111D101  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D10__13C110D110  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D10__13C111D110  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D11__13C110D111  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D11__13C111D111                                       +  Vf_SerPool3_GlyPool3_13C110D100_13C11D10  +  Vf_SerPool3_GlyPool3_13C111D100_13C11D10  +  Vf_SerPool3_GlyPool3_13C110D101_13C11D10  +  Vf_SerPool3_GlyPool3_13C111D101_13C11D10  +  Vf_SerPool3_GlyPool3_13C110D110_13C11D10  +  Vf_SerPool3_GlyPool3_13C111D110_13C11D10  +  Vf_SerPool3_GlyPool3_13C110D111_13C11D10  +  Vf_SerPool3_GlyPool3_13C111D111_13C11D10                                              -  Vr_Vin_Gly3_13C11D10_13C11D10  +  Vf_Vin_Gly3_13C11D10_13C11D10;              %  cGlyPool3_13C11D10 
		dxdt(321, 1) =  - Vf_GlyPool3_Degradation_13C00D11_13C00D11  -  Vf_GlyPool3_CO2_13C00D11_13C0D11  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D11__13C000  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D11__13C100  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D00__13C000D100  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D00__13C001D100  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D01__13C000D101  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D01__13C001D101  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D10__13C000D110  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D10__13C001D110  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D11__13C000D111  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D11__13C001D111  +  Vr_GlyPool3_CO2_13C0D11_13C00D11                                                                                                                                                                                                                                                                                                                                                                                                                      -  Vr_Vin_Gly3_13C00D11_13C00D11  +  Vf_Vin_Gly3_13C00D11_13C00D11;              %  cGlyPool3_13C00D11 
		dxdt(322, 1) =  - Vf_GlyPool3_Degradation_13C01D11_13C01D11  -  Vf_GlyPool3_CO2_13C01D11_13C1D11  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D11__13C001  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D11__13C101  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D00__13C010D100  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D00__13C011D100  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D01__13C010D101  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D01__13C011D101  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D10__13C010D110  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D10__13C011D110  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D11__13C010D111  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D11__13C011D111  +  Vr_GlyPool3_CO2_13C1D11_13C01D11                                                                                                                                                                                                                                                                                                                                                                                                                      -  Vr_Vin_Gly3_13C01D11_13C01D11  +  Vf_Vin_Gly3_13C01D11_13C01D11;              %  cGlyPool3_13C01D11 
		dxdt(323, 1) =  - Vf_GlyPool3_Degradation_13C10D11_13C10D11  -  Vf_GlyPool3_CO2_13C10D11_13C0D11  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D11__13C010  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D11__13C110  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D00__13C100D100  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D00__13C101D100  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D01__13C100D101  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D01__13C101D101  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D10__13C100D110  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D10__13C101D110  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D11__13C100D111  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D11__13C101D111                                                                                                                                                                                                                                                                                                                                                                                                                                                           -  Vr_Vin_Gly3_13C10D11_13C10D11  +  Vf_Vin_Gly3_13C10D11_13C10D11;              %  cGlyPool3_13C10D11 
		dxdt(324, 1) =  - Vf_GlyPool3_Degradation_13C11D11_13C11D11  -  Vf_GlyPool3_CO2_13C11D11_13C1D11  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D11__13C011  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D11__13C111  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D00__13C110D100  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D00__13C111D100  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D01__13C110D101  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D01__13C111D101  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D10__13C110D110  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D10__13C111D110  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D11__13C110D111  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D11__13C111D111                                                                                                                                                                                                                                                                                                                                                                                                                                                           -  Vr_Vin_Gly3_13C11D11_13C11D11  +  Vf_Vin_Gly3_13C11D11_13C11D11;              %  cGlyPool3_13C11D11 



    % mSerPool
    %  dxdt(1, 1) = (mSerPool) = 
    % - Vf_SerPoolMitochon_GlyPoolMitochon 
    % - Vr_SerPool2_SerPoolMitochon 
    % + Vr_SerPoolMitochon_GlyPoolMitochon 
    % + Vf_SerPool2_SerPoolMitochon;
    
		dxdt(325, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D000_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C000D000_13C000D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D00__13C000D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D00__13C000D000  +  Vf_SerPool2_SerPoolMitochon_13C000D000_13C000D000; 							%  mSerPool_13C000D000 
		dxdt(326, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D000_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C001D000_13C001D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D00__13C001D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D00__13C001D000  +  Vf_SerPool2_SerPoolMitochon_13C001D000_13C001D000;               %  mSerPool_13C001D000 
		dxdt(327, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D000_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C010D000_13C010D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D00__13C010D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D00__13C010D000  +  Vf_SerPool2_SerPoolMitochon_13C010D000_13C010D000;               %  mSerPool_13C010D000 
		dxdt(328, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D000_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C011D000_13C011D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D00__13C011D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D00__13C011D000  +  Vf_SerPool2_SerPoolMitochon_13C011D000_13C011D000;               %  mSerPool_13C011D000 
		dxdt(329, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D000_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C100D000_13C100D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D00__13C100D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D00__13C100D000  +  Vf_SerPool2_SerPoolMitochon_13C100D000_13C100D000;               %  mSerPool_13C100D000 
		dxdt(330, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D000_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C101D000_13C101D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D00__13C101D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D00__13C101D000  +  Vf_SerPool2_SerPoolMitochon_13C101D000_13C101D000;               %  mSerPool_13C101D000 
		dxdt(331, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D000_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C110D000_13C110D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D00__13C110D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D00__13C110D000  +  Vf_SerPool2_SerPoolMitochon_13C110D000_13C110D000;               %  mSerPool_13C110D000 
		dxdt(332, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D000_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C111D000_13C111D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D00__13C111D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D00__13C111D000  +  Vf_SerPool2_SerPoolMitochon_13C111D000_13C111D000;               %  mSerPool_13C111D000 
		dxdt(333, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D001_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C000D001_13C000D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D01__13C000D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D01__13C000D001  +  Vf_SerPool2_SerPoolMitochon_13C000D001_13C000D001;               %  mSerPool_13C000D001 
		dxdt(334, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D001_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C001D001_13C001D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D01__13C001D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D01__13C001D001  +  Vf_SerPool2_SerPoolMitochon_13C001D001_13C001D001;               %  mSerPool_13C001D001 
		dxdt(335, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D001_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C010D001_13C010D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D01__13C010D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D01__13C010D001  +  Vf_SerPool2_SerPoolMitochon_13C010D001_13C010D001;               %  mSerPool_13C010D001 
		dxdt(336, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D001_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C011D001_13C011D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D01__13C011D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D01__13C011D001  +  Vf_SerPool2_SerPoolMitochon_13C011D001_13C011D001;               %  mSerPool_13C011D001 
		dxdt(337, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D001_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C100D001_13C100D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D01__13C100D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D01__13C100D001  +  Vf_SerPool2_SerPoolMitochon_13C100D001_13C100D001;               %  mSerPool_13C100D001 
		dxdt(338, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D001_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C101D001_13C101D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D01__13C101D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D01__13C101D001  +  Vf_SerPool2_SerPoolMitochon_13C101D001_13C101D001;               %  mSerPool_13C101D001 
		dxdt(339, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D001_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C110D001_13C110D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D01__13C110D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D01__13C110D001  +  Vf_SerPool2_SerPoolMitochon_13C110D001_13C110D001;               %  mSerPool_13C110D001 
		dxdt(340, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D001_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C111D001_13C111D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D01__13C111D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D01__13C111D001  +  Vf_SerPool2_SerPoolMitochon_13C111D001_13C111D001;               %  mSerPool_13C111D001 
		dxdt(341, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D010_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C000D010_13C000D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D10__13C000D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D10__13C000D010  +  Vf_SerPool2_SerPoolMitochon_13C000D010_13C000D010;               %  mSerPool_13C000D010 
		dxdt(342, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D010_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C001D010_13C001D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D10__13C001D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D10__13C001D010  +  Vf_SerPool2_SerPoolMitochon_13C001D010_13C001D010;               %  mSerPool_13C001D010 
		dxdt(343, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D010_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C010D010_13C010D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D10__13C010D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D10__13C010D010  +  Vf_SerPool2_SerPoolMitochon_13C010D010_13C010D010;               %  mSerPool_13C010D010 
		dxdt(344, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D010_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C011D010_13C011D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D10__13C011D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D10__13C011D010  +  Vf_SerPool2_SerPoolMitochon_13C011D010_13C011D010;               %  mSerPool_13C011D010 
		dxdt(345, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D010_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C100D010_13C100D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D10__13C100D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D10__13C100D010  +  Vf_SerPool2_SerPoolMitochon_13C100D010_13C100D010;               %  mSerPool_13C100D010 
		dxdt(346, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D010_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C101D010_13C101D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D10__13C101D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D10__13C101D010  +  Vf_SerPool2_SerPoolMitochon_13C101D010_13C101D010;               %  mSerPool_13C101D010 
		dxdt(347, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D010_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C110D010_13C110D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D10__13C110D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D10__13C110D010  +  Vf_SerPool2_SerPoolMitochon_13C110D010_13C110D010;               %  mSerPool_13C110D010 
		dxdt(348, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D010_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C111D010_13C111D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D10__13C111D010  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D10__13C111D010  +  Vf_SerPool2_SerPoolMitochon_13C111D010_13C111D010;               %  mSerPool_13C111D010 
		dxdt(349, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D011_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C000D011_13C000D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D11__13C000D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D11__13C000D011  +  Vf_SerPool2_SerPoolMitochon_13C000D011_13C000D011;               %  mSerPool_13C000D011 
		dxdt(350, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D011_13C00D00  -  Vr_SerPool2_SerPoolMitochon_13C001D011_13C001D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D11__13C001D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D11__13C001D011  +  Vf_SerPool2_SerPoolMitochon_13C001D011_13C001D011;               %  mSerPool_13C001D011 
		dxdt(351, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D011_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C010D011_13C010D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D11__13C010D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D11__13C010D011  +  Vf_SerPool2_SerPoolMitochon_13C010D011_13C010D011;               %  mSerPool_13C010D011 
		dxdt(352, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D011_13C01D00  -  Vr_SerPool2_SerPoolMitochon_13C011D011_13C011D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D11__13C011D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D11__13C011D011  +  Vf_SerPool2_SerPoolMitochon_13C011D011_13C011D011;               %  mSerPool_13C011D011 
		dxdt(353, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D011_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C100D011_13C100D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D11__13C100D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D11__13C100D011  +  Vf_SerPool2_SerPoolMitochon_13C100D011_13C100D011;               %  mSerPool_13C100D011 
		dxdt(354, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D011_13C10D00  -  Vr_SerPool2_SerPoolMitochon_13C101D011_13C101D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D11__13C101D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D11__13C101D011  +  Vf_SerPool2_SerPoolMitochon_13C101D011_13C101D011;               %  mSerPool_13C101D011 
		dxdt(355, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D011_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C110D011_13C110D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D11__13C110D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D11__13C110D011  +  Vf_SerPool2_SerPoolMitochon_13C110D011_13C110D011;               %  mSerPool_13C110D011 
		dxdt(356, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D011_13C11D00  -  Vr_SerPool2_SerPoolMitochon_13C111D011_13C111D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D11__13C111D011  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D11__13C111D011  +  Vf_SerPool2_SerPoolMitochon_13C111D011_13C111D011;               %  mSerPool_13C111D011 
		dxdt(357, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D100_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C000D100_13C000D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D00__13C000D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D00__13C000D100  +  Vf_SerPool2_SerPoolMitochon_13C000D100_13C000D100;               %  mSerPool_13C000D100 
		dxdt(358, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D100_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C001D100_13C001D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D00__13C001D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D00__13C001D100  +  Vf_SerPool2_SerPoolMitochon_13C001D100_13C001D100;               %  mSerPool_13C001D100 
		dxdt(359, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D100_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C010D100_13C010D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D00__13C010D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D00__13C010D100  +  Vf_SerPool2_SerPoolMitochon_13C010D100_13C010D100;               %  mSerPool_13C010D100 
		dxdt(360, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D100_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C011D100_13C011D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D00__13C011D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D00__13C011D100  +  Vf_SerPool2_SerPoolMitochon_13C011D100_13C011D100;               %  mSerPool_13C011D100 
		dxdt(361, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D100_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C100D100_13C100D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D00__13C100D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D00__13C100D100  +  Vf_SerPool2_SerPoolMitochon_13C100D100_13C100D100;               %  mSerPool_13C100D100 
		dxdt(362, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D100_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C101D100_13C101D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D00__13C101D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D00__13C101D100  +  Vf_SerPool2_SerPoolMitochon_13C101D100_13C101D100;               %  mSerPool_13C101D100 
		dxdt(363, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D100_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C110D100_13C110D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D00__13C110D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D00__13C110D100  +  Vf_SerPool2_SerPoolMitochon_13C110D100_13C110D100;               %  mSerPool_13C110D100 
		dxdt(364, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D100_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C111D100_13C111D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D00__13C111D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D00__13C111D100  +  Vf_SerPool2_SerPoolMitochon_13C111D100_13C111D100;               %  mSerPool_13C111D100 
		dxdt(365, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D101_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C000D101_13C000D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D01__13C000D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D01__13C000D101  +  Vf_SerPool2_SerPoolMitochon_13C000D101_13C000D101;               %  mSerPool_13C000D101 
		dxdt(366, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D101_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C001D101_13C001D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D01__13C001D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D01__13C001D101  +  Vf_SerPool2_SerPoolMitochon_13C001D101_13C001D101;               %  mSerPool_13C001D101 
		dxdt(367, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D101_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C010D101_13C010D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D01__13C010D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D01__13C010D101  +  Vf_SerPool2_SerPoolMitochon_13C010D101_13C010D101;               %  mSerPool_13C010D101 
		dxdt(368, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D101_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C011D101_13C011D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D01__13C011D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D01__13C011D101  +  Vf_SerPool2_SerPoolMitochon_13C011D101_13C011D101;               %  mSerPool_13C011D101 
		dxdt(369, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D101_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C100D101_13C100D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D01__13C100D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D01__13C100D101  +  Vf_SerPool2_SerPoolMitochon_13C100D101_13C100D101;               %  mSerPool_13C100D101 
		dxdt(370, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D101_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C101D101_13C101D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D01__13C101D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D01__13C101D101  +  Vf_SerPool2_SerPoolMitochon_13C101D101_13C101D101;               %  mSerPool_13C101D101 
		dxdt(371, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D101_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C110D101_13C110D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D01__13C110D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D01__13C110D101  +  Vf_SerPool2_SerPoolMitochon_13C110D101_13C110D101;               %  mSerPool_13C110D101 
		dxdt(372, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D101_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C111D101_13C111D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D01__13C111D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D01__13C111D101  +  Vf_SerPool2_SerPoolMitochon_13C111D101_13C111D101;               %  mSerPool_13C111D101 
		dxdt(373, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D110_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C000D110_13C000D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D10__13C000D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D10__13C000D110  +  Vf_SerPool2_SerPoolMitochon_13C000D110_13C000D110;               %  mSerPool_13C000D110 
		dxdt(374, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D110_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C001D110_13C001D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D10__13C001D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D10__13C001D110  +  Vf_SerPool2_SerPoolMitochon_13C001D110_13C001D110;               %  mSerPool_13C001D110 
		dxdt(375, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D110_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C010D110_13C010D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D10__13C010D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D10__13C010D110  +  Vf_SerPool2_SerPoolMitochon_13C010D110_13C010D110;               %  mSerPool_13C010D110 
		dxdt(376, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D110_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C011D110_13C011D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D10__13C011D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D10__13C011D110  +  Vf_SerPool2_SerPoolMitochon_13C011D110_13C011D110;               %  mSerPool_13C011D110 
		dxdt(377, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D110_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C100D110_13C100D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D10__13C100D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D10__13C100D110  +  Vf_SerPool2_SerPoolMitochon_13C100D110_13C100D110;               %  mSerPool_13C100D110 
		dxdt(378, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D110_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C101D110_13C101D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D10__13C101D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D10__13C101D110  +  Vf_SerPool2_SerPoolMitochon_13C101D110_13C101D110;               %  mSerPool_13C101D110 
		dxdt(379, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D110_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C110D110_13C110D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D10__13C110D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D10__13C110D110  +  Vf_SerPool2_SerPoolMitochon_13C110D110_13C110D110;               %  mSerPool_13C110D110 
		dxdt(380, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D110_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C111D110_13C111D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D10__13C111D110  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D10__13C111D110  +  Vf_SerPool2_SerPoolMitochon_13C111D110_13C111D110;               %  mSerPool_13C111D110 
		dxdt(381, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D111_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C000D111_13C000D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D11__13C000D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D11__13C000D111  +  Vf_SerPool2_SerPoolMitochon_13C000D111_13C000D111;               %  mSerPool_13C000D111 
		dxdt(382, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C001D111_13C00D10  -  Vr_SerPool2_SerPoolMitochon_13C001D111_13C001D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D11__13C001D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D11__13C001D111  +  Vf_SerPool2_SerPoolMitochon_13C001D111_13C001D111;               %  mSerPool_13C001D111 
		dxdt(383, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C010D111_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C010D111_13C010D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D11__13C010D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D11__13C010D111  +  Vf_SerPool2_SerPoolMitochon_13C010D111_13C010D111;               %  mSerPool_13C010D111 
		dxdt(384, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C011D111_13C01D10  -  Vr_SerPool2_SerPoolMitochon_13C011D111_13C011D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D11__13C011D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D11__13C011D111  +  Vf_SerPool2_SerPoolMitochon_13C011D111_13C011D111;               %  mSerPool_13C011D111 
		dxdt(385, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C100D111_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C100D111_13C100D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D11__13C100D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D11__13C100D111  +  Vf_SerPool2_SerPoolMitochon_13C100D111_13C100D111;               %  mSerPool_13C100D111 
		dxdt(386, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C101D111_13C10D10  -  Vr_SerPool2_SerPoolMitochon_13C101D111_13C101D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D11__13C101D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D11__13C101D111  +  Vf_SerPool2_SerPoolMitochon_13C101D111_13C101D111;               %  mSerPool_13C101D111 
		dxdt(387, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C110D111_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C110D111_13C110D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D11__13C110D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D11__13C110D111  +  Vf_SerPool2_SerPoolMitochon_13C110D111_13C110D111;               %  mSerPool_13C110D111 
		dxdt(388, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C111D111_13C11D10  -  Vr_SerPool2_SerPoolMitochon_13C111D111_13C111D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D11__13C111D111  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D11__13C111D111  +  Vf_SerPool2_SerPoolMitochon_13C111D111_13C111D111;               %  mSerPool_13C111D111 
           

    
    % mGlyPool
    %  dxdt(1, 1) = (mGlyPool) = 
    % - Vf_GlyPoolMitochon_CO2 
    % - Vr_SerPoolMitochon_GlyPoolMitochon 
    % - Vr_GlyPool2_GlyPoolMitochon 
    % + Vr_GlyPoolMitochon_CO2 
    % + Vf_SerPoolMitochon_GlyPoolMitochon 
    % + Vf_GlyPool2_GlyPoolMitochon;
    
		dxdt(389, 1) =  - Vf_GlyPoolMitochon_CO2_13C00D00_13C0D00  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D00__13C000D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D00__13C001D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D01__13C000D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D01__13C001D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D10__13C000D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D10__13C001D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D11__13C000D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D11__13C001D011  -  Vr_GlyPool2_GlyPoolMitochon_13C00D00_13C00D00  +  Vr_GlyPoolMitochon_CO2_13C0D00_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D000_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D000_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D001_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D001_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D010_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D010_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D011_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D011_13C00D00  +  Vf_GlyPool2_GlyPoolMitochon_13C00D00_13C00D00;    				%  mGlyPool_13C00D00 
		dxdt(390, 1) =  - Vf_GlyPoolMitochon_CO2_13C01D00_13C1D00  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D00__13C010D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D00__13C011D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D01__13C010D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D01__13C011D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D10__13C010D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D10__13C011D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D11__13C010D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D11__13C011D011  -  Vr_GlyPool2_GlyPoolMitochon_13C01D00_13C01D00  +  Vr_GlyPoolMitochon_CO2_13C1D00_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D000_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D000_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D001_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D001_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D010_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D010_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D011_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D011_13C01D00  +  Vf_GlyPool2_GlyPoolMitochon_13C01D00_13C01D00;            %  mGlyPool_13C01D00 
		dxdt(391, 1) =  - Vf_GlyPoolMitochon_CO2_13C10D00_13C0D00  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D00__13C100D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D00__13C101D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D01__13C100D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D01__13C101D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D10__13C100D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D10__13C101D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D11__13C100D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D11__13C101D011  -  Vr_GlyPool2_GlyPoolMitochon_13C10D00_13C10D00                                              +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D000_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D000_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D001_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D001_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D010_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D010_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D011_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D011_13C10D00  +  Vf_GlyPool2_GlyPoolMitochon_13C10D00_13C10D00;            %  mGlyPool_13C10D00 
		dxdt(392, 1) =  - Vf_GlyPoolMitochon_CO2_13C11D00_13C1D00  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D00__13C110D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D00__13C111D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D01__13C110D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D01__13C111D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D10__13C110D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D10__13C111D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D11__13C110D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D11__13C111D011  -  Vr_GlyPool2_GlyPoolMitochon_13C11D00_13C11D00                                              +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D000_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D000_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D001_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D001_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D010_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D010_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D011_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D011_13C11D00  +  Vf_GlyPool2_GlyPoolMitochon_13C11D00_13C11D00;            %  mGlyPool_13C11D00 
		dxdt(393, 1) =  - Vf_GlyPoolMitochon_CO2_13C00D01_13C0D01  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D00__13C000D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D00__13C001D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D01__13C000D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D01__13C001D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D10__13C000D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D10__13C001D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D11__13C000D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D11__13C001D011  -  Vr_GlyPool2_GlyPoolMitochon_13C00D01_13C00D01  +  Vr_GlyPoolMitochon_CO2_13C0D01_13C00D01                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          +  Vf_GlyPool2_GlyPoolMitochon_13C00D01_13C00D01;            %  mGlyPool_13C00D01 
		dxdt(394, 1) =  - Vf_GlyPoolMitochon_CO2_13C01D01_13C1D01  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D00__13C010D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D00__13C011D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D01__13C010D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D01__13C011D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D10__13C010D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D10__13C011D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D11__13C010D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D11__13C011D011  -  Vr_GlyPool2_GlyPoolMitochon_13C01D01_13C01D01  +  Vr_GlyPoolMitochon_CO2_13C1D01_13C01D01                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          +  Vf_GlyPool2_GlyPoolMitochon_13C01D01_13C01D01;            %  mGlyPool_13C01D01 
		dxdt(395, 1) =  - Vf_GlyPoolMitochon_CO2_13C10D01_13C0D01  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D00__13C100D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D00__13C101D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D01__13C100D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D01__13C101D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D10__13C100D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D10__13C101D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D11__13C100D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D11__13C101D011  -  Vr_GlyPool2_GlyPoolMitochon_13C10D01_13C10D01                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      +  Vf_GlyPool2_GlyPoolMitochon_13C10D01_13C10D01;            %  mGlyPool_13C10D01 
		dxdt(396, 1) =  - Vf_GlyPoolMitochon_CO2_13C11D01_13C1D01  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D00__13C110D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D00__13C111D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D01__13C110D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D01__13C111D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D10__13C110D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D10__13C111D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D11__13C110D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D11__13C111D011  -  Vr_GlyPool2_GlyPoolMitochon_13C11D01_13C11D01                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      +  Vf_GlyPool2_GlyPoolMitochon_13C11D01_13C11D01;            %  mGlyPool_13C11D01 
		dxdt(397, 1) =  - Vf_GlyPoolMitochon_CO2_13C00D10_13C0D10  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D00__13C000D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D00__13C001D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D01__13C000D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D01__13C001D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D10__13C000D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D10__13C001D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D11__13C000D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D11__13C001D111  -  Vr_GlyPool2_GlyPoolMitochon_13C00D10_13C00D10  +  Vr_GlyPoolMitochon_CO2_13C0D10_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D100_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D100_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D101_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D101_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D110_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D110_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D111_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D111_13C00D10  +  Vf_GlyPool2_GlyPoolMitochon_13C00D10_13C00D10;            %  mGlyPool_13C00D10 
		dxdt(398, 1) =  - Vf_GlyPoolMitochon_CO2_13C01D10_13C1D10  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D00__13C010D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D00__13C011D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D01__13C010D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D01__13C011D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D10__13C010D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D10__13C011D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D11__13C010D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D11__13C011D111  -  Vr_GlyPool2_GlyPoolMitochon_13C01D10_13C01D10  +  Vr_GlyPoolMitochon_CO2_13C1D10_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D100_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D100_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D101_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D101_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D110_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D110_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D111_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D111_13C01D10  +  Vf_GlyPool2_GlyPoolMitochon_13C01D10_13C01D10;            %  mGlyPool_13C01D10 
		dxdt(399, 1) =  - Vf_GlyPoolMitochon_CO2_13C10D10_13C0D10  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D00__13C100D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D00__13C101D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D01__13C100D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D01__13C101D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D10__13C100D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D10__13C101D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D11__13C100D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D11__13C101D111  -  Vr_GlyPool2_GlyPoolMitochon_13C10D10_13C10D10                                              +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D100_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D100_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D101_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D101_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D110_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D110_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D111_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D111_13C10D10  +  Vf_GlyPool2_GlyPoolMitochon_13C10D10_13C10D10;            %  mGlyPool_13C10D10 
		dxdt(400, 1) =  - Vf_GlyPoolMitochon_CO2_13C11D10_13C1D10  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D00__13C110D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D00__13C111D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D01__13C110D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D01__13C111D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D10__13C110D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D10__13C111D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D11__13C110D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D11__13C111D111  -  Vr_GlyPool2_GlyPoolMitochon_13C11D10_13C11D10                                              +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D100_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D100_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D101_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D101_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D110_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D110_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D111_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D111_13C11D10  +  Vf_GlyPool2_GlyPoolMitochon_13C11D10_13C11D10;            %  mGlyPool_13C11D10 
		dxdt(401, 1) =  - Vf_GlyPoolMitochon_CO2_13C00D11_13C0D11  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D00__13C000D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D00__13C001D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D01__13C000D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D01__13C001D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D10__13C000D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D10__13C001D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D11__13C000D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D11__13C001D111  -  Vr_GlyPool2_GlyPoolMitochon_13C00D11_13C00D11  +  Vr_GlyPoolMitochon_CO2_13C0D11_13C00D11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          +  Vf_GlyPool2_GlyPoolMitochon_13C00D11_13C00D11;            %  mGlyPool_13C00D11 
		dxdt(402, 1) =  - Vf_GlyPoolMitochon_CO2_13C01D11_13C1D11  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D00__13C010D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D00__13C011D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D01__13C010D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D01__13C011D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D10__13C010D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D10__13C011D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D11__13C010D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D11__13C011D111  -  Vr_GlyPool2_GlyPoolMitochon_13C01D11_13C01D11  +  Vr_GlyPoolMitochon_CO2_13C1D11_13C01D11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          +  Vf_GlyPool2_GlyPoolMitochon_13C01D11_13C01D11;            %  mGlyPool_13C01D11 
		dxdt(403, 1) =  - Vf_GlyPoolMitochon_CO2_13C10D11_13C0D11  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D00__13C100D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D00__13C101D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D01__13C100D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D01__13C101D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D10__13C100D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D10__13C101D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D11__13C100D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D11__13C101D111  -  Vr_GlyPool2_GlyPoolMitochon_13C10D11_13C10D11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      +  Vf_GlyPool2_GlyPoolMitochon_13C10D11_13C10D11;            %  mGlyPool_13C10D11 
		dxdt(404, 1) =  - Vf_GlyPoolMitochon_CO2_13C11D11_13C1D11  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D00__13C110D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D00__13C111D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D01__13C110D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D01__13C111D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D10__13C110D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D10__13C111D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D11__13C110D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D11__13C111D111  -  Vr_GlyPool2_GlyPoolMitochon_13C11D11_13C11D11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      +  Vf_GlyPool2_GlyPoolMitochon_13C11D11_13C11D11;            %  mGlyPool_13C11D11 



    % cMethyleneTHF
    %  dxdt(1, 1) = (cMethyleneTHF) = 
    % - Vf_MethyleneTHF_FormylTHF_Cytoplasm 
    % - Vf_MethyleneTHF_Transport 
    % - Vr_SerPool1_GlyPool1    
    % - Vr_SerPool2_GlyPool2     
    % - Vr_SerPool3_GlyPool3     
    % - Vr_GlyPool1_CO2     
    % - Vr_GlyPool2_CO2     
    % - Vr_GlyPool3_CO2     
    % + Vr_MethyleneTHF_FormylTHF_Cytoplasm 
    % + Vr_MethyleneTHF_Transport 
    % + Vf_SerPool1_GlyPool1 
    % + Vf_SerPool2_GlyPool2 
    % + Vf_SerPool3_GlyPool3
    % + Vf_GlyPool1_CO2 
    % + Vf_GlyPool2_CO2;
    % + Vf_GlyPool3_CO2;

		dxdt(405, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D00_13C0D0  -  Vf_MethyleneTHF_Transport_13C0D00_13C0D00  -  Vr_SerPool1_GlyPool1_13C00D00_13C0D00__13C000D000  -  Vr_SerPool1_GlyPool1_13C01D00_13C0D00__13C010D000  -  Vr_SerPool1_GlyPool1_13C10D00_13C0D00__13C100D000  -  Vr_SerPool1_GlyPool1_13C11D00_13C0D00__13C110D000  -  Vr_SerPool1_GlyPool1_13C00D01_13C0D00__13C000D000  -  Vr_SerPool1_GlyPool1_13C01D01_13C0D00__13C010D000  -  Vr_SerPool1_GlyPool1_13C10D01_13C0D00__13C100D000  -  Vr_SerPool1_GlyPool1_13C11D01_13C0D00__13C110D000  -  Vr_SerPool1_GlyPool1_13C00D10_13C0D00__13C000D100  -  Vr_SerPool1_GlyPool1_13C01D10_13C0D00__13C010D100  -  Vr_SerPool1_GlyPool1_13C10D10_13C0D00__13C100D100  -  Vr_SerPool1_GlyPool1_13C11D10_13C0D00__13C110D100  -  Vr_SerPool1_GlyPool1_13C00D11_13C0D00__13C000D100  -  Vr_SerPool1_GlyPool1_13C01D11_13C0D00__13C010D100  -  Vr_SerPool1_GlyPool1_13C10D11_13C0D00__13C100D100  -  Vr_SerPool1_GlyPool1_13C11D11_13C0D00__13C110D100  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D00__13C000D000  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D00__13C010D000  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D00__13C100D000  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D00__13C110D000  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D00__13C000D000  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D00__13C010D000  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D00__13C100D000  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D00__13C110D000  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D00__13C000D100  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D00__13C010D100  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D00__13C100D100  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D00__13C110D100  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D00__13C000D100  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D00__13C010D100  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D00__13C100D100  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D00__13C110D100  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D00__13C000D000  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D00__13C010D000  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D00__13C100D000  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D00__13C110D000  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D00__13C000D000  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D00__13C010D000  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D00__13C100D000  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D00__13C110D000  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D00__13C000D100  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D00__13C010D100  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D00__13C100D100  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D00__13C110D100  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D00__13C000D100  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D00__13C010D100  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D00__13C100D100  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D00__13C110D100  -  Vr_GlyPool1_CO2_13C0D00_13C00D00  -  Vr_GlyPool2_CO2_13C0D00_13C00D00  -  Vr_GlyPool3_CO2_13C0D00_13C00D00  +  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C0D0_13C0D00  +  Vr_MethyleneTHF_Transport_13C0D00_13C0D00  +  Vf_SerPool1_GlyPool1_13C000D000_13C00D00  +  Vf_SerPool1_GlyPool1_13C010D000_13C01D00  +  Vf_SerPool1_GlyPool1_13C100D000_13C10D00  +  Vf_SerPool1_GlyPool1_13C110D000_13C11D00  +  Vf_SerPool1_GlyPool1_13C000D100_13C00D10  +  Vf_SerPool1_GlyPool1_13C010D100_13C01D10  +  Vf_SerPool1_GlyPool1_13C100D100_13C10D10  +  Vf_SerPool1_GlyPool1_13C110D100_13C11D10  +  Vf_SerPool2_GlyPool2_13C000D000_13C00D00  +  Vf_SerPool2_GlyPool2_13C010D000_13C01D00  +  Vf_SerPool2_GlyPool2_13C100D000_13C10D00  +  Vf_SerPool2_GlyPool2_13C110D000_13C11D00  +  Vf_SerPool2_GlyPool2_13C000D100_13C00D10  +  Vf_SerPool2_GlyPool2_13C010D100_13C01D10  +  Vf_SerPool2_GlyPool2_13C100D100_13C10D10  +  Vf_SerPool2_GlyPool2_13C110D100_13C11D10  +  Vf_SerPool3_GlyPool3_13C000D000_13C00D00  +  Vf_SerPool3_GlyPool3_13C010D000_13C01D00  +  Vf_SerPool3_GlyPool3_13C100D000_13C10D00  +  Vf_SerPool3_GlyPool3_13C110D000_13C11D00  +  Vf_SerPool3_GlyPool3_13C000D100_13C00D10  +  Vf_SerPool3_GlyPool3_13C010D100_13C01D10  +  Vf_SerPool3_GlyPool3_13C100D100_13C10D10  +  Vf_SerPool3_GlyPool3_13C110D100_13C11D10  + ... 
Vf_GlyPool1_CO2_13C00D00_13C0D00  +  Vf_GlyPool1_CO2_13C10D00_13C0D00  +  Vf_GlyPool2_CO2_13C00D00_13C0D00  +  Vf_GlyPool2_CO2_13C10D00_13C0D00  +  Vf_GlyPool3_CO2_13C00D00_13C0D00  +  Vf_GlyPool3_CO2_13C10D00_13C0D00;
		dxdt(406, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D00_13C1D0  -  Vf_MethyleneTHF_Transport_13C1D00_13C1D00  -  Vr_SerPool1_GlyPool1_13C00D00_13C1D00__13C001D000  -  Vr_SerPool1_GlyPool1_13C01D00_13C1D00__13C011D000  -  Vr_SerPool1_GlyPool1_13C10D00_13C1D00__13C101D000  -  Vr_SerPool1_GlyPool1_13C11D00_13C1D00__13C111D000  -  Vr_SerPool1_GlyPool1_13C00D01_13C1D00__13C001D000  -  Vr_SerPool1_GlyPool1_13C01D01_13C1D00__13C011D000  -  Vr_SerPool1_GlyPool1_13C10D01_13C1D00__13C101D000  -  Vr_SerPool1_GlyPool1_13C11D01_13C1D00__13C111D000  -  Vr_SerPool1_GlyPool1_13C00D10_13C1D00__13C001D100  -  Vr_SerPool1_GlyPool1_13C01D10_13C1D00__13C011D100  -  Vr_SerPool1_GlyPool1_13C10D10_13C1D00__13C101D100  -  Vr_SerPool1_GlyPool1_13C11D10_13C1D00__13C111D100  -  Vr_SerPool1_GlyPool1_13C00D11_13C1D00__13C001D100  -  Vr_SerPool1_GlyPool1_13C01D11_13C1D00__13C011D100  -  Vr_SerPool1_GlyPool1_13C10D11_13C1D00__13C101D100  -  Vr_SerPool1_GlyPool1_13C11D11_13C1D00__13C111D100  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D00__13C001D000  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D00__13C011D000  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D00__13C101D000  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D00__13C111D000  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D00__13C001D000  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D00__13C011D000  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D00__13C101D000  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D00__13C111D000  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D00__13C001D100  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D00__13C011D100  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D00__13C101D100  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D00__13C111D100  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D00__13C001D100  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D00__13C011D100  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D00__13C101D100  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D00__13C111D100  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D00__13C001D000  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D00__13C011D000  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D00__13C101D000  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D00__13C111D000  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D00__13C001D000  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D00__13C011D000  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D00__13C101D000  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D00__13C111D000  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D00__13C001D100  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D00__13C011D100  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D00__13C101D100  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D00__13C111D100  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D00__13C001D100  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D00__13C011D100  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D00__13C101D100  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D00__13C111D100  -  Vr_GlyPool1_CO2_13C1D00_13C01D00  -  Vr_GlyPool2_CO2_13C1D00_13C01D00  -  Vr_GlyPool3_CO2_13C1D00_13C01D00  +  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C1D0_13C1D00  +  Vr_MethyleneTHF_Transport_13C1D00_13C1D00  +  Vf_SerPool1_GlyPool1_13C001D100_13C00D10  +  Vf_SerPool1_GlyPool1_13C011D100_13C01D10  +  Vf_SerPool1_GlyPool1_13C101D100_13C10D10  +  Vf_SerPool1_GlyPool1_13C111D100_13C11D10  +  Vf_SerPool1_GlyPool1_13C001D000_13C00D00  +  Vf_SerPool1_GlyPool1_13C011D000_13C01D00  +  Vf_SerPool1_GlyPool1_13C101D000_13C10D00  +  Vf_SerPool1_GlyPool1_13C111D000_13C11D00  +  Vf_SerPool2_GlyPool2_13C001D100_13C00D10  +  Vf_SerPool2_GlyPool2_13C011D100_13C01D10  +  Vf_SerPool2_GlyPool2_13C101D100_13C10D10  +  Vf_SerPool2_GlyPool2_13C111D100_13C11D10  +  Vf_SerPool2_GlyPool2_13C001D000_13C00D00  +  Vf_SerPool2_GlyPool2_13C011D000_13C01D00  +  Vf_SerPool2_GlyPool2_13C101D000_13C10D00  +  Vf_SerPool2_GlyPool2_13C111D000_13C11D00  +  Vf_SerPool3_GlyPool3_13C001D100_13C00D10  +  Vf_SerPool3_GlyPool3_13C011D100_13C01D10  +  Vf_SerPool3_GlyPool3_13C101D100_13C10D10  +  Vf_SerPool3_GlyPool3_13C111D100_13C11D10  +  Vf_SerPool3_GlyPool3_13C001D000_13C00D00  +  Vf_SerPool3_GlyPool3_13C011D000_13C01D00  +  Vf_SerPool3_GlyPool3_13C101D000_13C10D00  +  Vf_SerPool3_GlyPool3_13C111D000_13C11D00  + ... 
Vf_GlyPool1_CO2_13C01D00_13C1D00  +  Vf_GlyPool1_CO2_13C11D00_13C1D00  +  Vf_GlyPool2_CO2_13C01D00_13C1D00  +  Vf_GlyPool2_CO2_13C11D00_13C1D00  +  Vf_GlyPool3_CO2_13C01D00_13C1D00  +  Vf_GlyPool3_CO2_13C11D00_13C1D00;
		dxdt(407, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D01_13C0D0  -  Vf_MethyleneTHF_Transport_13C0D01_13C0D01  -  Vr_SerPool1_GlyPool1_13C00D00_13C0D01__13C000D001  -  Vr_SerPool1_GlyPool1_13C01D00_13C0D01__13C010D001  -  Vr_SerPool1_GlyPool1_13C10D00_13C0D01__13C100D001  -  Vr_SerPool1_GlyPool1_13C11D00_13C0D01__13C110D001  -  Vr_SerPool1_GlyPool1_13C00D01_13C0D01__13C000D001  -  Vr_SerPool1_GlyPool1_13C01D01_13C0D01__13C010D001  -  Vr_SerPool1_GlyPool1_13C10D01_13C0D01__13C100D001  -  Vr_SerPool1_GlyPool1_13C11D01_13C0D01__13C110D001  -  Vr_SerPool1_GlyPool1_13C00D10_13C0D01__13C000D101  -  Vr_SerPool1_GlyPool1_13C01D10_13C0D01__13C010D101  -  Vr_SerPool1_GlyPool1_13C10D10_13C0D01__13C100D101  -  Vr_SerPool1_GlyPool1_13C11D10_13C0D01__13C110D101  -  Vr_SerPool1_GlyPool1_13C00D11_13C0D01__13C000D101  -  Vr_SerPool1_GlyPool1_13C01D11_13C0D01__13C010D101  -  Vr_SerPool1_GlyPool1_13C10D11_13C0D01__13C100D101  -  Vr_SerPool1_GlyPool1_13C11D11_13C0D01__13C110D101  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D01__13C000D001  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D01__13C010D001  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D01__13C100D001  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D01__13C110D001  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D01__13C000D001  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D01__13C010D001  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D01__13C100D001  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D01__13C110D001  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D01__13C000D101  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D01__13C010D101  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D01__13C100D101  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D01__13C110D101  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D01__13C000D101  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D01__13C010D101  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D01__13C100D101  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D01__13C110D101  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D01__13C000D001  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D01__13C010D001  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D01__13C100D001  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D01__13C110D001  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D01__13C000D001  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D01__13C010D001  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D01__13C100D001  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D01__13C110D001  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D01__13C000D101  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D01__13C010D101  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D01__13C100D101  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D01__13C110D101  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D01__13C000D101  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D01__13C010D101  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D01__13C100D101  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D01__13C110D101  -  Vr_GlyPool1_CO2_13C0D01_13C00D01  -  Vr_GlyPool2_CO2_13C0D01_13C00D01  -  Vr_GlyPool3_CO2_13C0D01_13C00D01                                                         +  Vr_MethyleneTHF_Transport_13C0D01_13C0D01  +  Vf_SerPool1_GlyPool1_13C000D001_13C00D00  +  Vf_SerPool1_GlyPool1_13C010D001_13C01D00  +  Vf_SerPool1_GlyPool1_13C100D001_13C10D00  +  Vf_SerPool1_GlyPool1_13C110D001_13C11D00  +  Vf_SerPool1_GlyPool1_13C000D101_13C00D10  +  Vf_SerPool1_GlyPool1_13C010D101_13C01D10  +  Vf_SerPool1_GlyPool1_13C100D101_13C10D10  +  Vf_SerPool1_GlyPool1_13C110D101_13C11D10  +  Vf_SerPool2_GlyPool2_13C000D001_13C00D00  +  Vf_SerPool2_GlyPool2_13C010D001_13C01D00  +  Vf_SerPool2_GlyPool2_13C100D001_13C10D00  +  Vf_SerPool2_GlyPool2_13C110D001_13C11D00  +  Vf_SerPool2_GlyPool2_13C000D101_13C00D10  +  Vf_SerPool2_GlyPool2_13C010D101_13C01D10  +  Vf_SerPool2_GlyPool2_13C100D101_13C10D10  +  Vf_SerPool2_GlyPool2_13C110D101_13C11D10  +  Vf_SerPool3_GlyPool3_13C000D001_13C00D00  +  Vf_SerPool3_GlyPool3_13C010D001_13C01D00  +  Vf_SerPool3_GlyPool3_13C100D001_13C10D00  +  Vf_SerPool3_GlyPool3_13C110D001_13C11D00  +  Vf_SerPool3_GlyPool3_13C000D101_13C00D10  +  Vf_SerPool3_GlyPool3_13C010D101_13C01D10  +  Vf_SerPool3_GlyPool3_13C100D101_13C10D10  +  Vf_SerPool3_GlyPool3_13C110D101_13C11D10  + ... 
Vf_GlyPool1_CO2_13C00D01_13C0D01  +  Vf_GlyPool1_CO2_13C10D01_13C0D01  +  Vf_GlyPool2_CO2_13C00D01_13C0D01  +  Vf_GlyPool2_CO2_13C10D01_13C0D01  +  Vf_GlyPool3_CO2_13C00D01_13C0D01  +  Vf_GlyPool3_CO2_13C10D01_13C0D01;
		dxdt(408, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D01_13C1D0  -  Vf_MethyleneTHF_Transport_13C1D01_13C1D01  -  Vr_SerPool1_GlyPool1_13C00D00_13C1D01__13C001D001  -  Vr_SerPool1_GlyPool1_13C01D00_13C1D01__13C011D001  -  Vr_SerPool1_GlyPool1_13C10D00_13C1D01__13C101D001  -  Vr_SerPool1_GlyPool1_13C11D00_13C1D01__13C111D001  -  Vr_SerPool1_GlyPool1_13C00D01_13C1D01__13C001D001  -  Vr_SerPool1_GlyPool1_13C01D01_13C1D01__13C011D001  -  Vr_SerPool1_GlyPool1_13C10D01_13C1D01__13C101D001  -  Vr_SerPool1_GlyPool1_13C11D01_13C1D01__13C111D001  -  Vr_SerPool1_GlyPool1_13C00D10_13C1D01__13C001D101  -  Vr_SerPool1_GlyPool1_13C01D10_13C1D01__13C011D101  -  Vr_SerPool1_GlyPool1_13C10D10_13C1D01__13C101D101  -  Vr_SerPool1_GlyPool1_13C11D10_13C1D01__13C111D101  -  Vr_SerPool1_GlyPool1_13C00D11_13C1D01__13C001D101  -  Vr_SerPool1_GlyPool1_13C01D11_13C1D01__13C011D101  -  Vr_SerPool1_GlyPool1_13C10D11_13C1D01__13C101D101  -  Vr_SerPool1_GlyPool1_13C11D11_13C1D01__13C111D101  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D01__13C001D001  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D01__13C011D001  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D01__13C101D001  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D01__13C111D001  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D01__13C001D001  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D01__13C011D001  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D01__13C101D001  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D01__13C111D001  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D01__13C001D101  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D01__13C011D101  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D01__13C101D101  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D01__13C111D101  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D01__13C001D101  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D01__13C011D101  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D01__13C101D101  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D01__13C111D101  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D01__13C001D001  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D01__13C011D001  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D01__13C101D001  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D01__13C111D001  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D01__13C001D001  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D01__13C011D001  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D01__13C101D001  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D01__13C111D001  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D01__13C001D101  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D01__13C011D101  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D01__13C101D101  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D01__13C111D101  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D01__13C001D101  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D01__13C011D101  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D01__13C101D101  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D01__13C111D101  -  Vr_GlyPool1_CO2_13C1D01_13C01D01  -  Vr_GlyPool2_CO2_13C1D01_13C01D01  -  Vr_GlyPool3_CO2_13C1D01_13C01D01                                                         +  Vr_MethyleneTHF_Transport_13C1D01_13C1D01  +  Vf_SerPool1_GlyPool1_13C001D101_13C00D10  +  Vf_SerPool1_GlyPool1_13C011D101_13C01D10  +  Vf_SerPool1_GlyPool1_13C101D101_13C10D10  +  Vf_SerPool1_GlyPool1_13C111D101_13C11D10  +  Vf_SerPool1_GlyPool1_13C001D001_13C00D00  +  Vf_SerPool1_GlyPool1_13C011D001_13C01D00  +  Vf_SerPool1_GlyPool1_13C101D001_13C10D00  +  Vf_SerPool1_GlyPool1_13C111D001_13C11D00  +  Vf_SerPool2_GlyPool2_13C001D101_13C00D10  +  Vf_SerPool2_GlyPool2_13C011D101_13C01D10  +  Vf_SerPool2_GlyPool2_13C101D101_13C10D10  +  Vf_SerPool2_GlyPool2_13C111D101_13C11D10  +  Vf_SerPool2_GlyPool2_13C001D001_13C00D00  +  Vf_SerPool2_GlyPool2_13C011D001_13C01D00  +  Vf_SerPool2_GlyPool2_13C101D001_13C10D00  +  Vf_SerPool2_GlyPool2_13C111D001_13C11D00  +  Vf_SerPool3_GlyPool3_13C001D101_13C00D10  +  Vf_SerPool3_GlyPool3_13C011D101_13C01D10  +  Vf_SerPool3_GlyPool3_13C101D101_13C10D10  +  Vf_SerPool3_GlyPool3_13C111D101_13C11D10  +  Vf_SerPool3_GlyPool3_13C001D001_13C00D00  +  Vf_SerPool3_GlyPool3_13C011D001_13C01D00  +  Vf_SerPool3_GlyPool3_13C101D001_13C10D00  +  Vf_SerPool3_GlyPool3_13C111D001_13C11D00  + ... 
Vf_GlyPool1_CO2_13C01D01_13C1D01  +  Vf_GlyPool1_CO2_13C11D01_13C1D01  +  Vf_GlyPool2_CO2_13C01D01_13C1D01  +  Vf_GlyPool2_CO2_13C11D01_13C1D01  +  Vf_GlyPool3_CO2_13C01D01_13C1D01  +  Vf_GlyPool3_CO2_13C11D01_13C1D01;
		dxdt(409, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D10_13C0D1  -  Vf_MethyleneTHF_Transport_13C0D10_13C0D10  -  Vr_SerPool1_GlyPool1_13C00D00_13C0D10__13C000D010  -  Vr_SerPool1_GlyPool1_13C01D00_13C0D10__13C010D010  -  Vr_SerPool1_GlyPool1_13C10D00_13C0D10__13C100D010  -  Vr_SerPool1_GlyPool1_13C11D00_13C0D10__13C110D010  -  Vr_SerPool1_GlyPool1_13C00D01_13C0D10__13C000D010  -  Vr_SerPool1_GlyPool1_13C01D01_13C0D10__13C010D010  -  Vr_SerPool1_GlyPool1_13C10D01_13C0D10__13C100D010  -  Vr_SerPool1_GlyPool1_13C11D01_13C0D10__13C110D010  -  Vr_SerPool1_GlyPool1_13C00D10_13C0D10__13C000D110  -  Vr_SerPool1_GlyPool1_13C01D10_13C0D10__13C010D110  -  Vr_SerPool1_GlyPool1_13C10D10_13C0D10__13C100D110  -  Vr_SerPool1_GlyPool1_13C11D10_13C0D10__13C110D110  -  Vr_SerPool1_GlyPool1_13C00D11_13C0D10__13C000D110  -  Vr_SerPool1_GlyPool1_13C01D11_13C0D10__13C010D110  -  Vr_SerPool1_GlyPool1_13C10D11_13C0D10__13C100D110  -  Vr_SerPool1_GlyPool1_13C11D11_13C0D10__13C110D110  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D10__13C000D010  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D10__13C010D010  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D10__13C100D010  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D10__13C110D010  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D10__13C000D010  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D10__13C010D010  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D10__13C100D010  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D10__13C110D010  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D10__13C000D110  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D10__13C010D110  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D10__13C100D110  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D10__13C110D110  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D10__13C000D110  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D10__13C010D110  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D10__13C100D110  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D10__13C110D110  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D10__13C000D010  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D10__13C010D010  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D10__13C100D010  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D10__13C110D010  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D10__13C000D010  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D10__13C010D010  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D10__13C100D010  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D10__13C110D010  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D10__13C000D110  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D10__13C010D110  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D10__13C100D110  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D10__13C110D110  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D10__13C000D110  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D10__13C010D110  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D10__13C100D110  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D10__13C110D110  -  Vr_GlyPool1_CO2_13C0D10_13C00D10  -  Vr_GlyPool2_CO2_13C0D10_13C00D10  -  Vr_GlyPool3_CO2_13C0D10_13C00D10  +  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C0D1_13C0D10  +  Vr_MethyleneTHF_Transport_13C0D10_13C0D10  +  Vf_SerPool1_GlyPool1_13C000D010_13C00D00  +  Vf_SerPool1_GlyPool1_13C010D010_13C01D00  +  Vf_SerPool1_GlyPool1_13C100D010_13C10D00  +  Vf_SerPool1_GlyPool1_13C110D010_13C11D00  +  Vf_SerPool1_GlyPool1_13C000D110_13C00D10  +  Vf_SerPool1_GlyPool1_13C010D110_13C01D10  +  Vf_SerPool1_GlyPool1_13C100D110_13C10D10  +  Vf_SerPool1_GlyPool1_13C110D110_13C11D10  +  Vf_SerPool2_GlyPool2_13C000D010_13C00D00  +  Vf_SerPool2_GlyPool2_13C010D010_13C01D00  +  Vf_SerPool2_GlyPool2_13C100D010_13C10D00  +  Vf_SerPool2_GlyPool2_13C110D010_13C11D00  +  Vf_SerPool2_GlyPool2_13C000D110_13C00D10  +  Vf_SerPool2_GlyPool2_13C010D110_13C01D10  +  Vf_SerPool2_GlyPool2_13C100D110_13C10D10  +  Vf_SerPool2_GlyPool2_13C110D110_13C11D10  +  Vf_SerPool3_GlyPool3_13C000D010_13C00D00  +  Vf_SerPool3_GlyPool3_13C010D010_13C01D00  +  Vf_SerPool3_GlyPool3_13C100D010_13C10D00  +  Vf_SerPool3_GlyPool3_13C110D010_13C11D00  +  Vf_SerPool3_GlyPool3_13C000D110_13C00D10  +  Vf_SerPool3_GlyPool3_13C010D110_13C01D10  +  Vf_SerPool3_GlyPool3_13C100D110_13C10D10  +  Vf_SerPool3_GlyPool3_13C110D110_13C11D10  + ... 
Vf_GlyPool1_CO2_13C00D10_13C0D10  +  Vf_GlyPool1_CO2_13C10D10_13C0D10  +  Vf_GlyPool2_CO2_13C00D10_13C0D10  +  Vf_GlyPool2_CO2_13C10D10_13C0D10  +  Vf_GlyPool3_CO2_13C00D10_13C0D10  +  Vf_GlyPool3_CO2_13C10D10_13C0D10;
		dxdt(410, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D10_13C1D1  -  Vf_MethyleneTHF_Transport_13C1D10_13C1D10  -  Vr_SerPool1_GlyPool1_13C00D00_13C1D10__13C001D010  -  Vr_SerPool1_GlyPool1_13C01D00_13C1D10__13C011D010  -  Vr_SerPool1_GlyPool1_13C10D00_13C1D10__13C101D010  -  Vr_SerPool1_GlyPool1_13C11D00_13C1D10__13C111D010  -  Vr_SerPool1_GlyPool1_13C00D01_13C1D10__13C001D010  -  Vr_SerPool1_GlyPool1_13C01D01_13C1D10__13C011D010  -  Vr_SerPool1_GlyPool1_13C10D01_13C1D10__13C101D010  -  Vr_SerPool1_GlyPool1_13C11D01_13C1D10__13C111D010  -  Vr_SerPool1_GlyPool1_13C00D10_13C1D10__13C001D110  -  Vr_SerPool1_GlyPool1_13C01D10_13C1D10__13C011D110  -  Vr_SerPool1_GlyPool1_13C10D10_13C1D10__13C101D110  -  Vr_SerPool1_GlyPool1_13C11D10_13C1D10__13C111D110  -  Vr_SerPool1_GlyPool1_13C00D11_13C1D10__13C001D110  -  Vr_SerPool1_GlyPool1_13C01D11_13C1D10__13C011D110  -  Vr_SerPool1_GlyPool1_13C10D11_13C1D10__13C101D110  -  Vr_SerPool1_GlyPool1_13C11D11_13C1D10__13C111D110  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D10__13C001D010  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D10__13C011D010  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D10__13C101D010  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D10__13C111D010  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D10__13C001D010  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D10__13C011D010  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D10__13C101D010  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D10__13C111D010  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D10__13C001D110  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D10__13C011D110  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D10__13C101D110  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D10__13C111D110  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D10__13C001D110  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D10__13C011D110  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D10__13C101D110  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D10__13C111D110  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D10__13C001D010  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D10__13C011D010  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D10__13C101D010  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D10__13C111D010  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D10__13C001D010  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D10__13C011D010  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D10__13C101D010  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D10__13C111D010  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D10__13C001D110  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D10__13C011D110  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D10__13C101D110  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D10__13C111D110  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D10__13C001D110  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D10__13C011D110  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D10__13C101D110  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D10__13C111D110  -  Vr_GlyPool1_CO2_13C1D10_13C01D10  -  Vr_GlyPool2_CO2_13C1D10_13C01D10  -  Vr_GlyPool3_CO2_13C1D10_13C01D10  +  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C1D1_13C1D10  +  Vr_MethyleneTHF_Transport_13C1D10_13C1D10  +  Vf_SerPool1_GlyPool1_13C001D110_13C00D10  +  Vf_SerPool1_GlyPool1_13C011D110_13C01D10  +  Vf_SerPool1_GlyPool1_13C101D110_13C10D10  +  Vf_SerPool1_GlyPool1_13C111D110_13C11D10  +  Vf_SerPool1_GlyPool1_13C001D010_13C00D00  +  Vf_SerPool1_GlyPool1_13C011D010_13C01D00  +  Vf_SerPool1_GlyPool1_13C101D010_13C10D00  +  Vf_SerPool1_GlyPool1_13C111D010_13C11D00  +  Vf_SerPool2_GlyPool2_13C001D110_13C00D10  +  Vf_SerPool2_GlyPool2_13C011D110_13C01D10  +  Vf_SerPool2_GlyPool2_13C101D110_13C10D10  +  Vf_SerPool2_GlyPool2_13C111D110_13C11D10  +  Vf_SerPool2_GlyPool2_13C001D010_13C00D00  +  Vf_SerPool2_GlyPool2_13C011D010_13C01D00  +  Vf_SerPool2_GlyPool2_13C101D010_13C10D00  +  Vf_SerPool2_GlyPool2_13C111D010_13C11D00  +  Vf_SerPool3_GlyPool3_13C001D110_13C00D10  +  Vf_SerPool3_GlyPool3_13C011D110_13C01D10  +  Vf_SerPool3_GlyPool3_13C101D110_13C10D10  +  Vf_SerPool3_GlyPool3_13C111D110_13C11D10  +  Vf_SerPool3_GlyPool3_13C001D010_13C00D00  +  Vf_SerPool3_GlyPool3_13C011D010_13C01D00  +  Vf_SerPool3_GlyPool3_13C101D010_13C10D00  +  Vf_SerPool3_GlyPool3_13C111D010_13C11D00  + ... 
Vf_GlyPool1_CO2_13C01D10_13C1D10  +  Vf_GlyPool1_CO2_13C11D10_13C1D10  +  Vf_GlyPool2_CO2_13C01D10_13C1D10  +  Vf_GlyPool2_CO2_13C11D10_13C1D10  +  Vf_GlyPool3_CO2_13C01D10_13C1D10  +  Vf_GlyPool3_CO2_13C11D10_13C1D10;
		dxdt(411, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D11_13C0D1  -  Vf_MethyleneTHF_Transport_13C0D11_13C0D11  -  Vr_SerPool1_GlyPool1_13C00D00_13C0D11__13C000D011  -  Vr_SerPool1_GlyPool1_13C01D00_13C0D11__13C010D011  -  Vr_SerPool1_GlyPool1_13C10D00_13C0D11__13C100D011  -  Vr_SerPool1_GlyPool1_13C11D00_13C0D11__13C110D011  -  Vr_SerPool1_GlyPool1_13C00D01_13C0D11__13C000D011  -  Vr_SerPool1_GlyPool1_13C01D01_13C0D11__13C010D011  -  Vr_SerPool1_GlyPool1_13C10D01_13C0D11__13C100D011  -  Vr_SerPool1_GlyPool1_13C11D01_13C0D11__13C110D011  -  Vr_SerPool1_GlyPool1_13C00D10_13C0D11__13C000D111  -  Vr_SerPool1_GlyPool1_13C01D10_13C0D11__13C010D111  -  Vr_SerPool1_GlyPool1_13C10D10_13C0D11__13C100D111  -  Vr_SerPool1_GlyPool1_13C11D10_13C0D11__13C110D111  -  Vr_SerPool1_GlyPool1_13C00D11_13C0D11__13C000D111  -  Vr_SerPool1_GlyPool1_13C01D11_13C0D11__13C010D111  -  Vr_SerPool1_GlyPool1_13C10D11_13C0D11__13C100D111  -  Vr_SerPool1_GlyPool1_13C11D11_13C0D11__13C110D111  -  Vr_SerPool2_GlyPool2_13C00D00_13C0D11__13C000D011  -  Vr_SerPool2_GlyPool2_13C01D00_13C0D11__13C010D011  -  Vr_SerPool2_GlyPool2_13C10D00_13C0D11__13C100D011  -  Vr_SerPool2_GlyPool2_13C11D00_13C0D11__13C110D011  -  Vr_SerPool2_GlyPool2_13C00D01_13C0D11__13C000D011  -  Vr_SerPool2_GlyPool2_13C01D01_13C0D11__13C010D011  -  Vr_SerPool2_GlyPool2_13C10D01_13C0D11__13C100D011  -  Vr_SerPool2_GlyPool2_13C11D01_13C0D11__13C110D011  -  Vr_SerPool2_GlyPool2_13C00D10_13C0D11__13C000D111  -  Vr_SerPool2_GlyPool2_13C01D10_13C0D11__13C010D111  -  Vr_SerPool2_GlyPool2_13C10D10_13C0D11__13C100D111  -  Vr_SerPool2_GlyPool2_13C11D10_13C0D11__13C110D111  -  Vr_SerPool2_GlyPool2_13C00D11_13C0D11__13C000D111  -  Vr_SerPool2_GlyPool2_13C01D11_13C0D11__13C010D111  -  Vr_SerPool2_GlyPool2_13C10D11_13C0D11__13C100D111  -  Vr_SerPool2_GlyPool2_13C11D11_13C0D11__13C110D111  -  Vr_SerPool3_GlyPool3_13C00D00_13C0D11__13C000D011  -  Vr_SerPool3_GlyPool3_13C01D00_13C0D11__13C010D011  -  Vr_SerPool3_GlyPool3_13C10D00_13C0D11__13C100D011  -  Vr_SerPool3_GlyPool3_13C11D00_13C0D11__13C110D011  -  Vr_SerPool3_GlyPool3_13C00D01_13C0D11__13C000D011  -  Vr_SerPool3_GlyPool3_13C01D01_13C0D11__13C010D011  -  Vr_SerPool3_GlyPool3_13C10D01_13C0D11__13C100D011  -  Vr_SerPool3_GlyPool3_13C11D01_13C0D11__13C110D011  -  Vr_SerPool3_GlyPool3_13C00D10_13C0D11__13C000D111  -  Vr_SerPool3_GlyPool3_13C01D10_13C0D11__13C010D111  -  Vr_SerPool3_GlyPool3_13C10D10_13C0D11__13C100D111  -  Vr_SerPool3_GlyPool3_13C11D10_13C0D11__13C110D111  -  Vr_SerPool3_GlyPool3_13C00D11_13C0D11__13C000D111  -  Vr_SerPool3_GlyPool3_13C01D11_13C0D11__13C010D111  -  Vr_SerPool3_GlyPool3_13C10D11_13C0D11__13C100D111  -  Vr_SerPool3_GlyPool3_13C11D11_13C0D11__13C110D111  -  Vr_GlyPool1_CO2_13C0D11_13C00D11  -  Vr_GlyPool2_CO2_13C0D11_13C00D11  -  Vr_GlyPool3_CO2_13C0D11_13C00D11                                                         +  Vr_MethyleneTHF_Transport_13C0D11_13C0D11  +  Vf_SerPool1_GlyPool1_13C000D011_13C00D00  +  Vf_SerPool1_GlyPool1_13C010D011_13C01D00  +  Vf_SerPool1_GlyPool1_13C100D011_13C10D00  +  Vf_SerPool1_GlyPool1_13C110D011_13C11D00  +  Vf_SerPool1_GlyPool1_13C000D111_13C00D10  +  Vf_SerPool1_GlyPool1_13C010D111_13C01D10  +  Vf_SerPool1_GlyPool1_13C100D111_13C10D10  +  Vf_SerPool1_GlyPool1_13C110D111_13C11D10  +  Vf_SerPool2_GlyPool2_13C000D011_13C00D00  +  Vf_SerPool2_GlyPool2_13C010D011_13C01D00  +  Vf_SerPool2_GlyPool2_13C100D011_13C10D00  +  Vf_SerPool2_GlyPool2_13C110D011_13C11D00  +  Vf_SerPool2_GlyPool2_13C000D111_13C00D10  +  Vf_SerPool2_GlyPool2_13C010D111_13C01D10  +  Vf_SerPool2_GlyPool2_13C100D111_13C10D10  +  Vf_SerPool2_GlyPool2_13C110D111_13C11D10  +  Vf_SerPool3_GlyPool3_13C000D011_13C00D00  +  Vf_SerPool3_GlyPool3_13C010D011_13C01D00  +  Vf_SerPool3_GlyPool3_13C100D011_13C10D00  +  Vf_SerPool3_GlyPool3_13C110D011_13C11D00  +  Vf_SerPool3_GlyPool3_13C000D111_13C00D10  +  Vf_SerPool3_GlyPool3_13C010D111_13C01D10  +  Vf_SerPool3_GlyPool3_13C100D111_13C10D10  +  Vf_SerPool3_GlyPool3_13C110D111_13C11D10  + ... 
Vf_GlyPool1_CO2_13C00D11_13C0D11  +  Vf_GlyPool1_CO2_13C10D11_13C0D11  +  Vf_GlyPool2_CO2_13C00D11_13C0D11  +  Vf_GlyPool2_CO2_13C10D11_13C0D11  +  Vf_GlyPool3_CO2_13C00D11_13C0D11  +  Vf_GlyPool3_CO2_13C10D11_13C0D11;
		dxdt(412, 1) =  - Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D11_13C1D1  -  Vf_MethyleneTHF_Transport_13C1D11_13C1D11  -  Vr_SerPool1_GlyPool1_13C00D00_13C1D11__13C001D011  -  Vr_SerPool1_GlyPool1_13C01D00_13C1D11__13C011D011  -  Vr_SerPool1_GlyPool1_13C10D00_13C1D11__13C101D011  -  Vr_SerPool1_GlyPool1_13C11D00_13C1D11__13C111D011  -  Vr_SerPool1_GlyPool1_13C00D01_13C1D11__13C001D011  -  Vr_SerPool1_GlyPool1_13C01D01_13C1D11__13C011D011  -  Vr_SerPool1_GlyPool1_13C10D01_13C1D11__13C101D011  -  Vr_SerPool1_GlyPool1_13C11D01_13C1D11__13C111D011  -  Vr_SerPool1_GlyPool1_13C00D10_13C1D11__13C001D111  -  Vr_SerPool1_GlyPool1_13C01D10_13C1D11__13C011D111  -  Vr_SerPool1_GlyPool1_13C10D10_13C1D11__13C101D111  -  Vr_SerPool1_GlyPool1_13C11D10_13C1D11__13C111D111  -  Vr_SerPool1_GlyPool1_13C00D11_13C1D11__13C001D111  -  Vr_SerPool1_GlyPool1_13C01D11_13C1D11__13C011D111  -  Vr_SerPool1_GlyPool1_13C10D11_13C1D11__13C101D111  -  Vr_SerPool1_GlyPool1_13C11D11_13C1D11__13C111D111  -  Vr_SerPool2_GlyPool2_13C00D00_13C1D11__13C001D011  -  Vr_SerPool2_GlyPool2_13C01D00_13C1D11__13C011D011  -  Vr_SerPool2_GlyPool2_13C10D00_13C1D11__13C101D011  -  Vr_SerPool2_GlyPool2_13C11D00_13C1D11__13C111D011  -  Vr_SerPool2_GlyPool2_13C00D01_13C1D11__13C001D011  -  Vr_SerPool2_GlyPool2_13C01D01_13C1D11__13C011D011  -  Vr_SerPool2_GlyPool2_13C10D01_13C1D11__13C101D011  -  Vr_SerPool2_GlyPool2_13C11D01_13C1D11__13C111D011  -  Vr_SerPool2_GlyPool2_13C00D10_13C1D11__13C001D111  -  Vr_SerPool2_GlyPool2_13C01D10_13C1D11__13C011D111  -  Vr_SerPool2_GlyPool2_13C10D10_13C1D11__13C101D111  -  Vr_SerPool2_GlyPool2_13C11D10_13C1D11__13C111D111  -  Vr_SerPool2_GlyPool2_13C00D11_13C1D11__13C001D111  -  Vr_SerPool2_GlyPool2_13C01D11_13C1D11__13C011D111  -  Vr_SerPool2_GlyPool2_13C10D11_13C1D11__13C101D111  -  Vr_SerPool2_GlyPool2_13C11D11_13C1D11__13C111D111  -  Vr_SerPool3_GlyPool3_13C00D00_13C1D11__13C001D011  -  Vr_SerPool3_GlyPool3_13C01D00_13C1D11__13C011D011  -  Vr_SerPool3_GlyPool3_13C10D00_13C1D11__13C101D011  -  Vr_SerPool3_GlyPool3_13C11D00_13C1D11__13C111D011  -  Vr_SerPool3_GlyPool3_13C00D01_13C1D11__13C001D011  -  Vr_SerPool3_GlyPool3_13C01D01_13C1D11__13C011D011  -  Vr_SerPool3_GlyPool3_13C10D01_13C1D11__13C101D011  -  Vr_SerPool3_GlyPool3_13C11D01_13C1D11__13C111D011  -  Vr_SerPool3_GlyPool3_13C00D10_13C1D11__13C001D111  -  Vr_SerPool3_GlyPool3_13C01D10_13C1D11__13C011D111  -  Vr_SerPool3_GlyPool3_13C10D10_13C1D11__13C101D111  -  Vr_SerPool3_GlyPool3_13C11D10_13C1D11__13C111D111  -  Vr_SerPool3_GlyPool3_13C00D11_13C1D11__13C001D111  -  Vr_SerPool3_GlyPool3_13C01D11_13C1D11__13C011D111  -  Vr_SerPool3_GlyPool3_13C10D11_13C1D11__13C101D111  -  Vr_SerPool3_GlyPool3_13C11D11_13C1D11__13C111D111  -  Vr_GlyPool1_CO2_13C1D11_13C01D11  -  Vr_GlyPool2_CO2_13C1D11_13C01D11  -  Vr_GlyPool3_CO2_13C1D11_13C01D11                                                         +  Vr_MethyleneTHF_Transport_13C1D11_13C1D11  +  Vf_SerPool1_GlyPool1_13C001D011_13C00D00  +  Vf_SerPool1_GlyPool1_13C011D011_13C01D00  +  Vf_SerPool1_GlyPool1_13C101D011_13C10D00  +  Vf_SerPool1_GlyPool1_13C111D011_13C11D00  +  Vf_SerPool1_GlyPool1_13C001D111_13C00D10  +  Vf_SerPool1_GlyPool1_13C011D111_13C01D10  +  Vf_SerPool1_GlyPool1_13C101D111_13C10D10  +  Vf_SerPool1_GlyPool1_13C111D111_13C11D10  +  Vf_SerPool2_GlyPool2_13C001D011_13C00D00  +  Vf_SerPool2_GlyPool2_13C011D011_13C01D00  +  Vf_SerPool2_GlyPool2_13C101D011_13C10D00  +  Vf_SerPool2_GlyPool2_13C111D011_13C11D00  +  Vf_SerPool2_GlyPool2_13C001D111_13C00D10  +  Vf_SerPool2_GlyPool2_13C011D111_13C01D10  +  Vf_SerPool2_GlyPool2_13C101D111_13C10D10  +  Vf_SerPool2_GlyPool2_13C111D111_13C11D10  +  Vf_SerPool3_GlyPool3_13C001D011_13C00D00  +  Vf_SerPool3_GlyPool3_13C011D011_13C01D00  +  Vf_SerPool3_GlyPool3_13C101D011_13C10D00  +  Vf_SerPool3_GlyPool3_13C111D011_13C11D00  +  Vf_SerPool3_GlyPool3_13C001D111_13C00D10  +  Vf_SerPool3_GlyPool3_13C011D111_13C01D10  +  Vf_SerPool3_GlyPool3_13C101D111_13C10D10  +  Vf_SerPool3_GlyPool3_13C111D111_13C11D10  + ... 
Vf_GlyPool1_CO2_13C01D11_13C1D11  +  Vf_GlyPool1_CO2_13C11D11_13C1D11  +  Vf_GlyPool2_CO2_13C01D11_13C1D11  +  Vf_GlyPool2_CO2_13C11D11_13C1D11  +  Vf_GlyPool3_CO2_13C01D11_13C1D11  +  Vf_GlyPool3_CO2_13C11D11_13C1D11;

% cMethyleneTHF_13C0D00 
% cMethyleneTHF_13C1D00 
% cMethyleneTHF_13C0D01 
% cMethyleneTHF_13C1D01 
% cMethyleneTHF_13C0D10 
% cMethyleneTHF_13C1D10 
% cMethyleneTHF_13C0D11 
% cMethyleneTHF_13C1D11 



    % cFormylTHF
    %  dxdt(1, 1) = (cFormylTHF) = 
    % - Vf_FormylTHF_Formate_Cytoplasm 
    % - Vf_GAR_FGAR    
    % - Vf_FGAR_AMP
    % - Vf_FormylTHF_Transport 
    % - Vr_MethyleneTHF_FormylTHF_Cytoplasm 
    % + Vr_FormylTHF_Formate_Cytoplasm 
    % + Vr_FormylTHF_Transport 
    % + Vf_MethyleneTHF_FormylTHF_Cytoplasm;

		dxdt(413, 1) =  - Vf_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0  -  Vf_GAR_FGAR_13C000_13C0D0__13C0000D0  -  Vf_GAR_FGAR_13C001_13C0D0__13C0010D0  -  Vf_GAR_FGAR_13C010_13C0D0__13C0100D0  -  Vf_GAR_FGAR_13C011_13C0D0__13C0110D0  -  Vf_GAR_FGAR_13C100_13C0D0__13C1000D0  -  Vf_GAR_FGAR_13C101_13C0D0__13C1010D0  -  Vf_GAR_FGAR_13C110_13C0D0__13C1100D0  -  Vf_GAR_FGAR_13C111_13C0D0__13C1110D0  -  Vf_FGAR_AMP_13C0000D0_13C0D0__13C00000D00  -  Vf_FGAR_AMP_13C0001D0_13C0D0__13C00010D00  -  Vf_FGAR_AMP_13C0010D0_13C0D0__13C00100D00  -  Vf_FGAR_AMP_13C0011D0_13C0D0__13C00110D00  -  Vf_FGAR_AMP_13C0100D0_13C0D0__13C01000D00  -  Vf_FGAR_AMP_13C0101D0_13C0D0__13C01010D00  -  Vf_FGAR_AMP_13C0110D0_13C0D0__13C01100D00  -  Vf_FGAR_AMP_13C0111D0_13C0D0__13C01110D00  -  Vf_FGAR_AMP_13C1000D0_13C0D0__13C10000D00  -  Vf_FGAR_AMP_13C1001D0_13C0D0__13C10010D00  -  Vf_FGAR_AMP_13C1010D0_13C0D0__13C10100D00  -  Vf_FGAR_AMP_13C1011D0_13C0D0__13C10110D00  -  Vf_FGAR_AMP_13C1100D0_13C0D0__13C11000D00  -  Vf_FGAR_AMP_13C1101D0_13C0D0__13C11010D00  -  Vf_FGAR_AMP_13C1110D0_13C0D0__13C11100D00  -  Vf_FGAR_AMP_13C1111D0_13C0D0__13C11110D00  -  Vf_FGAR_AMP_13C0000D1_13C0D0__13C00000D10  -  Vf_FGAR_AMP_13C0001D1_13C0D0__13C00010D10  -  Vf_FGAR_AMP_13C0010D1_13C0D0__13C00100D10  -  Vf_FGAR_AMP_13C0011D1_13C0D0__13C00110D10  -  Vf_FGAR_AMP_13C0100D1_13C0D0__13C01000D10  -  Vf_FGAR_AMP_13C0101D1_13C0D0__13C01010D10  -  Vf_FGAR_AMP_13C0110D1_13C0D0__13C01100D10  -  Vf_FGAR_AMP_13C0111D1_13C0D0__13C01110D10  -  Vf_FGAR_AMP_13C1000D1_13C0D0__13C10000D10  -  Vf_FGAR_AMP_13C1001D1_13C0D0__13C10010D10  -  Vf_FGAR_AMP_13C1010D1_13C0D0__13C10100D10  -  Vf_FGAR_AMP_13C1011D1_13C0D0__13C10110D10  -  Vf_FGAR_AMP_13C1100D1_13C0D0__13C11000D10  -  Vf_FGAR_AMP_13C1101D1_13C0D0__13C11010D10  -  Vf_FGAR_AMP_13C1110D1_13C0D0__13C11100D10  -  Vf_FGAR_AMP_13C1111D1_13C0D0__13C11110D10  -  Vf_FormylTHF_Transport_13C0D0_13C0D0  -  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C0D0_13C0D00  +  Vr_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0  +  Vr_FormylTHF_Transport_13C0D0_13C0D0  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D00_13C0D0  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D01_13C0D0;					%  cFormylTHF_13C0D0 
		dxdt(414, 1) =  - Vf_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0  -  Vf_GAR_FGAR_13C000_13C1D0__13C0001D0  -  Vf_GAR_FGAR_13C001_13C1D0__13C0011D0  -  Vf_GAR_FGAR_13C010_13C1D0__13C0101D0  -  Vf_GAR_FGAR_13C011_13C1D0__13C0111D0  -  Vf_GAR_FGAR_13C100_13C1D0__13C1001D0  -  Vf_GAR_FGAR_13C101_13C1D0__13C1011D0  -  Vf_GAR_FGAR_13C110_13C1D0__13C1101D0  -  Vf_GAR_FGAR_13C111_13C1D0__13C1111D0  -  Vf_FGAR_AMP_13C0000D0_13C1D0__13C00001D00  -  Vf_FGAR_AMP_13C0001D0_13C1D0__13C00011D00  -  Vf_FGAR_AMP_13C0010D0_13C1D0__13C00101D00  -  Vf_FGAR_AMP_13C0011D0_13C1D0__13C00111D00  -  Vf_FGAR_AMP_13C0100D0_13C1D0__13C01001D00  -  Vf_FGAR_AMP_13C0101D0_13C1D0__13C01011D00  -  Vf_FGAR_AMP_13C0110D0_13C1D0__13C01101D00  -  Vf_FGAR_AMP_13C0111D0_13C1D0__13C01111D00  -  Vf_FGAR_AMP_13C1000D0_13C1D0__13C10001D00  -  Vf_FGAR_AMP_13C1001D0_13C1D0__13C10011D00  -  Vf_FGAR_AMP_13C1010D0_13C1D0__13C10101D00  -  Vf_FGAR_AMP_13C1011D0_13C1D0__13C10111D00  -  Vf_FGAR_AMP_13C1100D0_13C1D0__13C11001D00  -  Vf_FGAR_AMP_13C1101D0_13C1D0__13C11011D00  -  Vf_FGAR_AMP_13C1110D0_13C1D0__13C11101D00  -  Vf_FGAR_AMP_13C1111D0_13C1D0__13C11111D00  -  Vf_FGAR_AMP_13C0000D1_13C1D0__13C00001D10  -  Vf_FGAR_AMP_13C0001D1_13C1D0__13C00011D10  -  Vf_FGAR_AMP_13C0010D1_13C1D0__13C00101D10  -  Vf_FGAR_AMP_13C0011D1_13C1D0__13C00111D10  -  Vf_FGAR_AMP_13C0100D1_13C1D0__13C01001D10  -  Vf_FGAR_AMP_13C0101D1_13C1D0__13C01011D10  -  Vf_FGAR_AMP_13C0110D1_13C1D0__13C01101D10  -  Vf_FGAR_AMP_13C0111D1_13C1D0__13C01111D10  -  Vf_FGAR_AMP_13C1000D1_13C1D0__13C10001D10  -  Vf_FGAR_AMP_13C1001D1_13C1D0__13C10011D10  -  Vf_FGAR_AMP_13C1010D1_13C1D0__13C10101D10  -  Vf_FGAR_AMP_13C1011D1_13C1D0__13C10111D10  -  Vf_FGAR_AMP_13C1100D1_13C1D0__13C11001D10  -  Vf_FGAR_AMP_13C1101D1_13C1D0__13C11011D10  -  Vf_FGAR_AMP_13C1110D1_13C1D0__13C11101D10  -  Vf_FGAR_AMP_13C1111D1_13C1D0__13C11111D10  -  Vf_FormylTHF_Transport_13C1D0_13C1D0  -  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C1D0_13C1D00  +  Vr_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0  +  Vr_FormylTHF_Transport_13C1D0_13C1D0  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D00_13C1D0  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D01_13C1D0;         %  cFormylTHF_13C1D0 
		dxdt(415, 1) =  - Vf_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1  -  Vf_GAR_FGAR_13C000_13C0D1__13C0000D1  -  Vf_GAR_FGAR_13C001_13C0D1__13C0010D1  -  Vf_GAR_FGAR_13C010_13C0D1__13C0100D1  -  Vf_GAR_FGAR_13C011_13C0D1__13C0110D1  -  Vf_GAR_FGAR_13C100_13C0D1__13C1000D1  -  Vf_GAR_FGAR_13C101_13C0D1__13C1010D1  -  Vf_GAR_FGAR_13C110_13C0D1__13C1100D1  -  Vf_GAR_FGAR_13C111_13C0D1__13C1110D1  -  Vf_FGAR_AMP_13C0000D0_13C0D1__13C00000D01  -  Vf_FGAR_AMP_13C0001D0_13C0D1__13C00010D01  -  Vf_FGAR_AMP_13C0010D0_13C0D1__13C00100D01  -  Vf_FGAR_AMP_13C0011D0_13C0D1__13C00110D01  -  Vf_FGAR_AMP_13C0100D0_13C0D1__13C01000D01  -  Vf_FGAR_AMP_13C0101D0_13C0D1__13C01010D01  -  Vf_FGAR_AMP_13C0110D0_13C0D1__13C01100D01  -  Vf_FGAR_AMP_13C0111D0_13C0D1__13C01110D01  -  Vf_FGAR_AMP_13C1000D0_13C0D1__13C10000D01  -  Vf_FGAR_AMP_13C1001D0_13C0D1__13C10010D01  -  Vf_FGAR_AMP_13C1010D0_13C0D1__13C10100D01  -  Vf_FGAR_AMP_13C1011D0_13C0D1__13C10110D01  -  Vf_FGAR_AMP_13C1100D0_13C0D1__13C11000D01  -  Vf_FGAR_AMP_13C1101D0_13C0D1__13C11010D01  -  Vf_FGAR_AMP_13C1110D0_13C0D1__13C11100D01  -  Vf_FGAR_AMP_13C1111D0_13C0D1__13C11110D01  -  Vf_FGAR_AMP_13C0000D1_13C0D1__13C00000D11  -  Vf_FGAR_AMP_13C0001D1_13C0D1__13C00010D11  -  Vf_FGAR_AMP_13C0010D1_13C0D1__13C00100D11  -  Vf_FGAR_AMP_13C0011D1_13C0D1__13C00110D11  -  Vf_FGAR_AMP_13C0100D1_13C0D1__13C01000D11  -  Vf_FGAR_AMP_13C0101D1_13C0D1__13C01010D11  -  Vf_FGAR_AMP_13C0110D1_13C0D1__13C01100D11  -  Vf_FGAR_AMP_13C0111D1_13C0D1__13C01110D11  -  Vf_FGAR_AMP_13C1000D1_13C0D1__13C10000D11  -  Vf_FGAR_AMP_13C1001D1_13C0D1__13C10010D11  -  Vf_FGAR_AMP_13C1010D1_13C0D1__13C10100D11  -  Vf_FGAR_AMP_13C1011D1_13C0D1__13C10110D11  -  Vf_FGAR_AMP_13C1100D1_13C0D1__13C11000D11  -  Vf_FGAR_AMP_13C1101D1_13C0D1__13C11010D11  -  Vf_FGAR_AMP_13C1110D1_13C0D1__13C11100D11  -  Vf_FGAR_AMP_13C1111D1_13C0D1__13C11110D11  -  Vf_FormylTHF_Transport_13C0D1_13C0D1  -  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C0D1_13C0D10  +  Vr_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1  +  Vr_FormylTHF_Transport_13C0D1_13C0D1  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D10_13C0D1  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C0D11_13C0D1;         %  cFormylTHF_13C0D1 
		dxdt(416, 1) =  - Vf_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1  -  Vf_GAR_FGAR_13C000_13C1D1__13C0001D1  -  Vf_GAR_FGAR_13C001_13C1D1__13C0011D1  -  Vf_GAR_FGAR_13C010_13C1D1__13C0101D1  -  Vf_GAR_FGAR_13C011_13C1D1__13C0111D1  -  Vf_GAR_FGAR_13C100_13C1D1__13C1001D1  -  Vf_GAR_FGAR_13C101_13C1D1__13C1011D1  -  Vf_GAR_FGAR_13C110_13C1D1__13C1101D1  -  Vf_GAR_FGAR_13C111_13C1D1__13C1111D1  -  Vf_FGAR_AMP_13C0000D0_13C1D1__13C00001D01  -  Vf_FGAR_AMP_13C0001D0_13C1D1__13C00011D01  -  Vf_FGAR_AMP_13C0010D0_13C1D1__13C00101D01  -  Vf_FGAR_AMP_13C0011D0_13C1D1__13C00111D01  -  Vf_FGAR_AMP_13C0100D0_13C1D1__13C01001D01  -  Vf_FGAR_AMP_13C0101D0_13C1D1__13C01011D01  -  Vf_FGAR_AMP_13C0110D0_13C1D1__13C01101D01  -  Vf_FGAR_AMP_13C0111D0_13C1D1__13C01111D01  -  Vf_FGAR_AMP_13C1000D0_13C1D1__13C10001D01  -  Vf_FGAR_AMP_13C1001D0_13C1D1__13C10011D01  -  Vf_FGAR_AMP_13C1010D0_13C1D1__13C10101D01  -  Vf_FGAR_AMP_13C1011D0_13C1D1__13C10111D01  -  Vf_FGAR_AMP_13C1100D0_13C1D1__13C11001D01  -  Vf_FGAR_AMP_13C1101D0_13C1D1__13C11011D01  -  Vf_FGAR_AMP_13C1110D0_13C1D1__13C11101D01  -  Vf_FGAR_AMP_13C1111D0_13C1D1__13C11111D01  -  Vf_FGAR_AMP_13C0000D1_13C1D1__13C00001D11  -  Vf_FGAR_AMP_13C0001D1_13C1D1__13C00011D11  -  Vf_FGAR_AMP_13C0010D1_13C1D1__13C00101D11  -  Vf_FGAR_AMP_13C0011D1_13C1D1__13C00111D11  -  Vf_FGAR_AMP_13C0100D1_13C1D1__13C01001D11  -  Vf_FGAR_AMP_13C0101D1_13C1D1__13C01011D11  -  Vf_FGAR_AMP_13C0110D1_13C1D1__13C01101D11  -  Vf_FGAR_AMP_13C0111D1_13C1D1__13C01111D11  -  Vf_FGAR_AMP_13C1000D1_13C1D1__13C10001D11  -  Vf_FGAR_AMP_13C1001D1_13C1D1__13C10011D11  -  Vf_FGAR_AMP_13C1010D1_13C1D1__13C10101D11  -  Vf_FGAR_AMP_13C1011D1_13C1D1__13C10111D11  -  Vf_FGAR_AMP_13C1100D1_13C1D1__13C11001D11  -  Vf_FGAR_AMP_13C1101D1_13C1D1__13C11011D11  -  Vf_FGAR_AMP_13C1110D1_13C1D1__13C11101D11  -  Vf_FGAR_AMP_13C1111D1_13C1D1__13C11111D11  -  Vf_FormylTHF_Transport_13C1D1_13C1D1  -  Vr_MethyleneTHF_FormylTHF_Cytoplasm_13C1D1_13C1D10  +  Vr_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1  +  Vr_FormylTHF_Transport_13C1D1_13C1D1  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D10_13C1D1  +  Vf_MethyleneTHF_FormylTHF_Cytoplasm_13C1D11_13C1D1;         %  cFormylTHF_13C1D1 


    
    % cFormate
    %  dxdt(1, 1) = (cFormate) = 
    % - Vf_Formate_Transport 
    % - Vr_FormylTHF_Formate_Cytoplasm 
    % + Vr_Formate_Transport 
    % + Vf_FormylTHF_Formate_Cytoplasm;
    
		dxdt(417, 1) =  - Vf_Formate_Transport_13C0D0  -  Vr_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0  +  Vr_Formate_Transport_13C0D0  +  Vf_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0;					%  cFormate_13C0D0 
		dxdt(418, 1) =  - Vf_Formate_Transport_13C1D0  -  Vr_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0  +  Vr_Formate_Transport_13C1D0  +  Vf_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0;          %  cFormate_13C1D0 
		dxdt(419, 1) =  - Vf_Formate_Transport_13C0D1  -  Vr_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1  +  Vr_Formate_Transport_13C0D1  +  Vf_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1;          %  cFormate_13C0D1 
		dxdt(420, 1) =  - Vf_Formate_Transport_13C1D1  -  Vr_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1  +  Vr_Formate_Transport_13C1D1  +  Vf_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1;          %  cFormate_13C1D1 



    % cTHF
    %  dxdt(1, 1) = (cTHF) =  
    % - Vf_SerPool1_GlyPool1 
    % - Vf_SerPool2_GlyPool2 
    % - Vf_SerPool3_GlyPool3 
    % - Vf_GlyPool1_CO2 
    % - Vf_THF_Transport 
    % - Vf_GlyPool2_CO2 
    % - Vf_GlyPool3_CO2 
    % - Vr_FormylTHF_Formate_Cytoplasm 
    % + Vr_SerPool1_GlyPool1 
    % + Vr_SerPool2_GlyPool2 
    % + Vr_SerPool3_GlyPool3
    % + Vr_GlyPool1_CO2 
    % + Vr_THF_Transport 
    % + Vr_GlyPool2_CO2 
    % + Vr_GlyPool3_CO2 
    % + Vf_FormylTHF_Formate_Cytoplasm 
    % + Vf_GAR_FGAR 
    % + Vf_FGAR_AMP;

		dxdt(421, 1) =  - Vf_SerPool1_GlyPool1_13C000D000_13C00D00  -  Vf_SerPool1_GlyPool1_13C001D000_13C00D00  -  Vf_SerPool1_GlyPool1_13C010D000_13C01D00  -  Vf_SerPool1_GlyPool1_13C011D000_13C01D00  -  Vf_SerPool1_GlyPool1_13C100D000_13C10D00  -  Vf_SerPool1_GlyPool1_13C101D000_13C10D00  -  Vf_SerPool1_GlyPool1_13C110D000_13C11D00  -  Vf_SerPool1_GlyPool1_13C111D000_13C11D00  -  Vf_SerPool1_GlyPool1_13C000D001_13C00D00  -  Vf_SerPool1_GlyPool1_13C001D001_13C00D00  -  Vf_SerPool1_GlyPool1_13C010D001_13C01D00  -  Vf_SerPool1_GlyPool1_13C011D001_13C01D00  -  Vf_SerPool1_GlyPool1_13C100D001_13C10D00  -  Vf_SerPool1_GlyPool1_13C101D001_13C10D00  -  Vf_SerPool1_GlyPool1_13C110D001_13C11D00  -  Vf_SerPool1_GlyPool1_13C111D001_13C11D00  -  Vf_SerPool1_GlyPool1_13C000D010_13C00D00  -  Vf_SerPool1_GlyPool1_13C001D010_13C00D00  -  Vf_SerPool1_GlyPool1_13C010D010_13C01D00  -  Vf_SerPool1_GlyPool1_13C011D010_13C01D00  -  Vf_SerPool1_GlyPool1_13C100D010_13C10D00  -  Vf_SerPool1_GlyPool1_13C101D010_13C10D00  -  Vf_SerPool1_GlyPool1_13C110D010_13C11D00  -  Vf_SerPool1_GlyPool1_13C111D010_13C11D00  -  Vf_SerPool1_GlyPool1_13C000D011_13C00D00  -  Vf_SerPool1_GlyPool1_13C001D011_13C00D00  -  Vf_SerPool1_GlyPool1_13C010D011_13C01D00  -  Vf_SerPool1_GlyPool1_13C011D011_13C01D00  -  Vf_SerPool1_GlyPool1_13C100D011_13C10D00  -  Vf_SerPool1_GlyPool1_13C101D011_13C10D00  -  Vf_SerPool1_GlyPool1_13C110D011_13C11D00  -  Vf_SerPool1_GlyPool1_13C111D011_13C11D00  -  Vf_SerPool1_GlyPool1_13C000D100_13C00D10  -  Vf_SerPool1_GlyPool1_13C001D100_13C00D10  -  Vf_SerPool1_GlyPool1_13C010D100_13C01D10  -  Vf_SerPool1_GlyPool1_13C011D100_13C01D10  -  Vf_SerPool1_GlyPool1_13C100D100_13C10D10  -  Vf_SerPool1_GlyPool1_13C101D100_13C10D10  -  Vf_SerPool1_GlyPool1_13C110D100_13C11D10  -  Vf_SerPool1_GlyPool1_13C111D100_13C11D10  -  Vf_SerPool1_GlyPool1_13C000D101_13C00D10  -  Vf_SerPool1_GlyPool1_13C001D101_13C00D10  -  Vf_SerPool1_GlyPool1_13C010D101_13C01D10  -  Vf_SerPool1_GlyPool1_13C011D101_13C01D10  -  Vf_SerPool1_GlyPool1_13C100D101_13C10D10  -  Vf_SerPool1_GlyPool1_13C101D101_13C10D10  -  Vf_SerPool1_GlyPool1_13C110D101_13C11D10  -  Vf_SerPool1_GlyPool1_13C111D101_13C11D10  -  Vf_SerPool1_GlyPool1_13C000D110_13C00D10  -  Vf_SerPool1_GlyPool1_13C001D110_13C00D10  -  Vf_SerPool1_GlyPool1_13C010D110_13C01D10  -  Vf_SerPool1_GlyPool1_13C011D110_13C01D10  -  Vf_SerPool1_GlyPool1_13C100D110_13C10D10  -  Vf_SerPool1_GlyPool1_13C101D110_13C10D10  -  Vf_SerPool1_GlyPool1_13C110D110_13C11D10  -  Vf_SerPool1_GlyPool1_13C111D110_13C11D10  -  Vf_SerPool1_GlyPool1_13C000D111_13C00D10  -  Vf_SerPool1_GlyPool1_13C001D111_13C00D10  -  Vf_SerPool1_GlyPool1_13C010D111_13C01D10  -  Vf_SerPool1_GlyPool1_13C011D111_13C01D10  -  Vf_SerPool1_GlyPool1_13C100D111_13C10D10  -  Vf_SerPool1_GlyPool1_13C101D111_13C10D10  -  Vf_SerPool1_GlyPool1_13C110D111_13C11D10  -  Vf_SerPool1_GlyPool1_13C111D111_13C11D10  ...
-  Vf_SerPool2_GlyPool2_13C000D000_13C00D00  -  Vf_SerPool2_GlyPool2_13C001D000_13C00D00  -  Vf_SerPool2_GlyPool2_13C010D000_13C01D00  -  Vf_SerPool2_GlyPool2_13C011D000_13C01D00  -  Vf_SerPool2_GlyPool2_13C100D000_13C10D00  -  Vf_SerPool2_GlyPool2_13C101D000_13C10D00  -  Vf_SerPool2_GlyPool2_13C110D000_13C11D00  -  Vf_SerPool2_GlyPool2_13C111D000_13C11D00  -  Vf_SerPool2_GlyPool2_13C000D001_13C00D00  -  Vf_SerPool2_GlyPool2_13C001D001_13C00D00  -  Vf_SerPool2_GlyPool2_13C010D001_13C01D00  -  Vf_SerPool2_GlyPool2_13C011D001_13C01D00  -  Vf_SerPool2_GlyPool2_13C100D001_13C10D00  -  Vf_SerPool2_GlyPool2_13C101D001_13C10D00  -  Vf_SerPool2_GlyPool2_13C110D001_13C11D00  -  Vf_SerPool2_GlyPool2_13C111D001_13C11D00  -  Vf_SerPool2_GlyPool2_13C000D010_13C00D00  -  Vf_SerPool2_GlyPool2_13C001D010_13C00D00  -  Vf_SerPool2_GlyPool2_13C010D010_13C01D00  -  Vf_SerPool2_GlyPool2_13C011D010_13C01D00  -  Vf_SerPool2_GlyPool2_13C100D010_13C10D00  -  Vf_SerPool2_GlyPool2_13C101D010_13C10D00  -  Vf_SerPool2_GlyPool2_13C110D010_13C11D00  -  Vf_SerPool2_GlyPool2_13C111D010_13C11D00  -  Vf_SerPool2_GlyPool2_13C000D011_13C00D00  -  Vf_SerPool2_GlyPool2_13C001D011_13C00D00  -  Vf_SerPool2_GlyPool2_13C010D011_13C01D00  -  Vf_SerPool2_GlyPool2_13C011D011_13C01D00  -  Vf_SerPool2_GlyPool2_13C100D011_13C10D00  -  Vf_SerPool2_GlyPool2_13C101D011_13C10D00  -  Vf_SerPool2_GlyPool2_13C110D011_13C11D00  -  Vf_SerPool2_GlyPool2_13C111D011_13C11D00  -  Vf_SerPool2_GlyPool2_13C000D100_13C00D10  -  Vf_SerPool2_GlyPool2_13C001D100_13C00D10  -  Vf_SerPool2_GlyPool2_13C010D100_13C01D10  -  Vf_SerPool2_GlyPool2_13C011D100_13C01D10  -  Vf_SerPool2_GlyPool2_13C100D100_13C10D10  -  Vf_SerPool2_GlyPool2_13C101D100_13C10D10  -  Vf_SerPool2_GlyPool2_13C110D100_13C11D10  -  Vf_SerPool2_GlyPool2_13C111D100_13C11D10  -  Vf_SerPool2_GlyPool2_13C000D101_13C00D10  -  Vf_SerPool2_GlyPool2_13C001D101_13C00D10  -  Vf_SerPool2_GlyPool2_13C010D101_13C01D10  -  Vf_SerPool2_GlyPool2_13C011D101_13C01D10  -  Vf_SerPool2_GlyPool2_13C100D101_13C10D10  -  Vf_SerPool2_GlyPool2_13C101D101_13C10D10  -  Vf_SerPool2_GlyPool2_13C110D101_13C11D10  -  Vf_SerPool2_GlyPool2_13C111D101_13C11D10  -  Vf_SerPool2_GlyPool2_13C000D110_13C00D10  -  Vf_SerPool2_GlyPool2_13C001D110_13C00D10  -  Vf_SerPool2_GlyPool2_13C010D110_13C01D10  -  Vf_SerPool2_GlyPool2_13C011D110_13C01D10  -  Vf_SerPool2_GlyPool2_13C100D110_13C10D10  -  Vf_SerPool2_GlyPool2_13C101D110_13C10D10  -  Vf_SerPool2_GlyPool2_13C110D110_13C11D10  -  Vf_SerPool2_GlyPool2_13C111D110_13C11D10  -  Vf_SerPool2_GlyPool2_13C000D111_13C00D10  -  Vf_SerPool2_GlyPool2_13C001D111_13C00D10  -  Vf_SerPool2_GlyPool2_13C010D111_13C01D10  -  Vf_SerPool2_GlyPool2_13C011D111_13C01D10  -  Vf_SerPool2_GlyPool2_13C100D111_13C10D10  -  Vf_SerPool2_GlyPool2_13C101D111_13C10D10  -  Vf_SerPool2_GlyPool2_13C110D111_13C11D10  -  Vf_SerPool2_GlyPool2_13C111D111_13C11D10  ...
-  Vf_SerPool3_GlyPool3_13C000D000_13C00D00  -  Vf_SerPool3_GlyPool3_13C001D000_13C00D00  -  Vf_SerPool3_GlyPool3_13C010D000_13C01D00  -  Vf_SerPool3_GlyPool3_13C011D000_13C01D00  -  Vf_SerPool3_GlyPool3_13C100D000_13C10D00  -  Vf_SerPool3_GlyPool3_13C101D000_13C10D00  -  Vf_SerPool3_GlyPool3_13C110D000_13C11D00  -  Vf_SerPool3_GlyPool3_13C111D000_13C11D00  -  Vf_SerPool3_GlyPool3_13C000D001_13C00D00  -  Vf_SerPool3_GlyPool3_13C001D001_13C00D00  -  Vf_SerPool3_GlyPool3_13C010D001_13C01D00  -  Vf_SerPool3_GlyPool3_13C011D001_13C01D00  -  Vf_SerPool3_GlyPool3_13C100D001_13C10D00  -  Vf_SerPool3_GlyPool3_13C101D001_13C10D00  -  Vf_SerPool3_GlyPool3_13C110D001_13C11D00  -  Vf_SerPool3_GlyPool3_13C111D001_13C11D00  -  Vf_SerPool3_GlyPool3_13C000D010_13C00D00  -  Vf_SerPool3_GlyPool3_13C001D010_13C00D00  -  Vf_SerPool3_GlyPool3_13C010D010_13C01D00  -  Vf_SerPool3_GlyPool3_13C011D010_13C01D00  -  Vf_SerPool3_GlyPool3_13C100D010_13C10D00  -  Vf_SerPool3_GlyPool3_13C101D010_13C10D00  -  Vf_SerPool3_GlyPool3_13C110D010_13C11D00  -  Vf_SerPool3_GlyPool3_13C111D010_13C11D00  -  Vf_SerPool3_GlyPool3_13C000D011_13C00D00  -  Vf_SerPool3_GlyPool3_13C001D011_13C00D00  -  Vf_SerPool3_GlyPool3_13C010D011_13C01D00  -  Vf_SerPool3_GlyPool3_13C011D011_13C01D00  -  Vf_SerPool3_GlyPool3_13C100D011_13C10D00  -  Vf_SerPool3_GlyPool3_13C101D011_13C10D00  -  Vf_SerPool3_GlyPool3_13C110D011_13C11D00  -  Vf_SerPool3_GlyPool3_13C111D011_13C11D00  -  Vf_SerPool3_GlyPool3_13C000D100_13C00D10  -  Vf_SerPool3_GlyPool3_13C001D100_13C00D10  -  Vf_SerPool3_GlyPool3_13C010D100_13C01D10  -  Vf_SerPool3_GlyPool3_13C011D100_13C01D10  -  Vf_SerPool3_GlyPool3_13C100D100_13C10D10  -  Vf_SerPool3_GlyPool3_13C101D100_13C10D10  -  Vf_SerPool3_GlyPool3_13C110D100_13C11D10  -  Vf_SerPool3_GlyPool3_13C111D100_13C11D10  -  Vf_SerPool3_GlyPool3_13C000D101_13C00D10  -  Vf_SerPool3_GlyPool3_13C001D101_13C00D10  -  Vf_SerPool3_GlyPool3_13C010D101_13C01D10  -  Vf_SerPool3_GlyPool3_13C011D101_13C01D10  -  Vf_SerPool3_GlyPool3_13C100D101_13C10D10  -  Vf_SerPool3_GlyPool3_13C101D101_13C10D10  -  Vf_SerPool3_GlyPool3_13C110D101_13C11D10  -  Vf_SerPool3_GlyPool3_13C111D101_13C11D10  -  Vf_SerPool3_GlyPool3_13C000D110_13C00D10  -  Vf_SerPool3_GlyPool3_13C001D110_13C00D10  -  Vf_SerPool3_GlyPool3_13C010D110_13C01D10  -  Vf_SerPool3_GlyPool3_13C011D110_13C01D10  -  Vf_SerPool3_GlyPool3_13C100D110_13C10D10  -  Vf_SerPool3_GlyPool3_13C101D110_13C10D10  -  Vf_SerPool3_GlyPool3_13C110D110_13C11D10  -  Vf_SerPool3_GlyPool3_13C111D110_13C11D10  -  Vf_SerPool3_GlyPool3_13C000D111_13C00D10  -  Vf_SerPool3_GlyPool3_13C001D111_13C00D10  -  Vf_SerPool3_GlyPool3_13C010D111_13C01D10  -  Vf_SerPool3_GlyPool3_13C011D111_13C01D10  -  Vf_SerPool3_GlyPool3_13C100D111_13C10D10  -  Vf_SerPool3_GlyPool3_13C101D111_13C10D10  -  Vf_SerPool3_GlyPool3_13C110D111_13C11D10  -  Vf_SerPool3_GlyPool3_13C111D111_13C11D10  ...
-  Vf_GlyPool1_CO2_13C00D00_13C0D00  -  Vf_GlyPool1_CO2_13C01D00_13C1D00  -  Vf_GlyPool1_CO2_13C10D00_13C0D00  -  Vf_GlyPool1_CO2_13C11D00_13C1D00  -  Vf_GlyPool1_CO2_13C00D01_13C0D01  -  Vf_GlyPool1_CO2_13C01D01_13C1D01  -  Vf_GlyPool1_CO2_13C10D01_13C0D01  -  Vf_GlyPool1_CO2_13C11D01_13C1D01  -  Vf_GlyPool1_CO2_13C00D10_13C0D10  -  Vf_GlyPool1_CO2_13C01D10_13C1D10  -  Vf_GlyPool1_CO2_13C10D10_13C0D10  -  Vf_GlyPool1_CO2_13C11D10_13C1D10  -  Vf_GlyPool1_CO2_13C00D11_13C0D11  -  Vf_GlyPool1_CO2_13C01D11_13C1D11  -  Vf_GlyPool1_CO2_13C10D11_13C0D11  -  Vf_GlyPool1_CO2_13C11D11_13C1D11  ...
-  Vf_THF_Transport  ...
-  Vf_GlyPool2_CO2_13C00D00_13C0D00  -  Vf_GlyPool2_CO2_13C01D00_13C1D00  -  Vf_GlyPool2_CO2_13C10D00_13C0D00  -  Vf_GlyPool2_CO2_13C11D00_13C1D00  -  Vf_GlyPool2_CO2_13C00D01_13C0D01  -  Vf_GlyPool2_CO2_13C01D01_13C1D01  -  Vf_GlyPool2_CO2_13C10D01_13C0D01  -  Vf_GlyPool2_CO2_13C11D01_13C1D01  -  Vf_GlyPool2_CO2_13C00D10_13C0D10  -  Vf_GlyPool2_CO2_13C01D10_13C1D10  -  Vf_GlyPool2_CO2_13C10D10_13C0D10  -  Vf_GlyPool2_CO2_13C11D10_13C1D10  -  Vf_GlyPool2_CO2_13C00D11_13C0D11  -  Vf_GlyPool2_CO2_13C01D11_13C1D11  -  Vf_GlyPool2_CO2_13C10D11_13C0D11  -  Vf_GlyPool2_CO2_13C11D11_13C1D11  ...
-  Vf_GlyPool3_CO2_13C00D00_13C0D00  -  Vf_GlyPool3_CO2_13C01D00_13C1D00  -  Vf_GlyPool3_CO2_13C10D00_13C0D00  -  Vf_GlyPool3_CO2_13C11D00_13C1D00  -  Vf_GlyPool3_CO2_13C00D01_13C0D01  -  Vf_GlyPool3_CO2_13C01D01_13C1D01  -  Vf_GlyPool3_CO2_13C10D01_13C0D01  -  Vf_GlyPool3_CO2_13C11D01_13C1D01  -  Vf_GlyPool3_CO2_13C00D10_13C0D10  -  Vf_GlyPool3_CO2_13C01D10_13C1D10  -  Vf_GlyPool3_CO2_13C10D10_13C0D10  -  Vf_GlyPool3_CO2_13C11D10_13C1D10  -  Vf_GlyPool3_CO2_13C00D11_13C0D11  -  Vf_GlyPool3_CO2_13C01D11_13C1D11  -  Vf_GlyPool3_CO2_13C10D11_13C0D11  -  Vf_GlyPool3_CO2_13C11D11_13C1D11  ...
-  Vr_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0  -  Vr_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0  -  Vr_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1  -  Vr_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1  ...
+  Vr_SerPool1_GlyPool1_13C00D00_13C0D00__13C000D000  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D00__13C010D000  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D00__13C100D000  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D00__13C110D000  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D00__13C000D000  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D00__13C010D000  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D00__13C100D000  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D00__13C110D000  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D00__13C000D100  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D00__13C010D100  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D00__13C100D100  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D00__13C110D100  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D00__13C000D100  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D00__13C010D100  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D00__13C100D100  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D00__13C110D100  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D00__13C001D000  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D00__13C011D000  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D00__13C101D000  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D00__13C111D000  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D00__13C001D000  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D00__13C011D000  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D00__13C101D000  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D00__13C111D000  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D00__13C001D100  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D00__13C011D100  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D00__13C101D100  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D00__13C111D100  +  Vr_SerPool1_GlyPool1_13C00D11_13C1D00__13C001D100  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D00__13C011D100  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D00__13C101D100  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D00__13C111D100  +  Vr_SerPool1_GlyPool1_13C00D00_13C0D01__13C000D001  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D01__13C010D001  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D01__13C100D001  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D01__13C110D001  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D01__13C000D001  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D01__13C010D001  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D01__13C100D001  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D01__13C110D001  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D01__13C000D101  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D01__13C010D101  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D01__13C100D101  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D01__13C110D101  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D01__13C000D101  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D01__13C010D101  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D01__13C100D101  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D01__13C110D101  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D01__13C001D001  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D01__13C011D001  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D01__13C101D001  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D01__13C111D001  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D01__13C001D001  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D01__13C011D001  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D01__13C101D001  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D01__13C111D001  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D01__13C001D101  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D01__13C011D101  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D01__13C101D101  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D01__13C111D101  ...
+  Vr_SerPool1_GlyPool1_13C00D11_13C1D01__13C001D101  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D01__13C011D101  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D01__13C101D101  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D01__13C111D101  +  Vr_SerPool1_GlyPool1_13C00D00_13C0D10__13C000D010  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D10__13C010D010  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D10__13C100D010  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D10__13C110D010  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D10__13C000D010  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D10__13C010D010  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D10__13C100D010  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D10__13C110D010  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D10__13C000D110  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D10__13C010D110  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D10__13C100D110  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D10__13C110D110  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D10__13C000D110  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D10__13C010D110  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D10__13C100D110  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D10__13C110D110  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D10__13C001D010  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D10__13C011D010  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D10__13C101D010  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D10__13C111D010  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D10__13C001D010  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D10__13C011D010  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D10__13C101D010  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D10__13C111D010  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D10__13C001D110  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D10__13C011D110  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D10__13C101D110  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D10__13C111D110  +  Vr_SerPool1_GlyPool1_13C00D11_13C1D10__13C001D110  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D10__13C011D110  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D10__13C101D110  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D10__13C111D110  +  Vr_SerPool1_GlyPool1_13C00D00_13C0D11__13C000D011  +  Vr_SerPool1_GlyPool1_13C01D00_13C0D11__13C010D011  +  Vr_SerPool1_GlyPool1_13C10D00_13C0D11__13C100D011  +  Vr_SerPool1_GlyPool1_13C11D00_13C0D11__13C110D011  +  Vr_SerPool1_GlyPool1_13C00D01_13C0D11__13C000D011  +  Vr_SerPool1_GlyPool1_13C01D01_13C0D11__13C010D011  +  Vr_SerPool1_GlyPool1_13C10D01_13C0D11__13C100D011  +  Vr_SerPool1_GlyPool1_13C11D01_13C0D11__13C110D011  +  Vr_SerPool1_GlyPool1_13C00D10_13C0D11__13C000D111  +  Vr_SerPool1_GlyPool1_13C01D10_13C0D11__13C010D111  +  Vr_SerPool1_GlyPool1_13C10D10_13C0D11__13C100D111  +  Vr_SerPool1_GlyPool1_13C11D10_13C0D11__13C110D111  +  Vr_SerPool1_GlyPool1_13C00D11_13C0D11__13C000D111  +  Vr_SerPool1_GlyPool1_13C01D11_13C0D11__13C010D111  +  Vr_SerPool1_GlyPool1_13C10D11_13C0D11__13C100D111  +  Vr_SerPool1_GlyPool1_13C11D11_13C0D11__13C110D111  +  Vr_SerPool1_GlyPool1_13C00D00_13C1D11__13C001D011  +  Vr_SerPool1_GlyPool1_13C01D00_13C1D11__13C011D011  +  Vr_SerPool1_GlyPool1_13C10D00_13C1D11__13C101D011  +  Vr_SerPool1_GlyPool1_13C11D00_13C1D11__13C111D011  +  Vr_SerPool1_GlyPool1_13C00D01_13C1D11__13C001D011  +  Vr_SerPool1_GlyPool1_13C01D01_13C1D11__13C011D011  +  Vr_SerPool1_GlyPool1_13C10D01_13C1D11__13C101D011  +  Vr_SerPool1_GlyPool1_13C11D01_13C1D11__13C111D011  +  Vr_SerPool1_GlyPool1_13C00D10_13C1D11__13C001D111  +  Vr_SerPool1_GlyPool1_13C01D10_13C1D11__13C011D111  +  Vr_SerPool1_GlyPool1_13C10D10_13C1D11__13C101D111  +  Vr_SerPool1_GlyPool1_13C11D10_13C1D11__13C111D111  +  Vr_SerPool1_GlyPool1_13C00D11_13C1D11__13C001D111  +  Vr_SerPool1_GlyPool1_13C01D11_13C1D11__13C011D111  +  Vr_SerPool1_GlyPool1_13C10D11_13C1D11__13C101D111  +  Vr_SerPool1_GlyPool1_13C11D11_13C1D11__13C111D111  ...
+  Vr_SerPool2_GlyPool2_13C00D00_13C0D00__13C000D000  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D00__13C010D000  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D00__13C100D000  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D00__13C110D000  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D00__13C000D000  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D00__13C010D000  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D00__13C100D000  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D00__13C110D000  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D00__13C000D100  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D00__13C010D100  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D00__13C100D100  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D00__13C110D100  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D00__13C000D100  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D00__13C010D100  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D00__13C100D100  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D00__13C110D100  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D00__13C001D000  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D00__13C011D000  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D00__13C101D000  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D00__13C111D000  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D00__13C001D000  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D00__13C011D000  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D00__13C101D000  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D00__13C111D000  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D00__13C001D100  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D00__13C011D100  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D00__13C101D100  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D00__13C111D100  +  Vr_SerPool2_GlyPool2_13C00D11_13C1D00__13C001D100  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D00__13C011D100  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D00__13C101D100  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D00__13C111D100  +  Vr_SerPool2_GlyPool2_13C00D00_13C0D01__13C000D001  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D01__13C010D001  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D01__13C100D001  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D01__13C110D001  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D01__13C000D001  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D01__13C010D001  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D01__13C100D001  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D01__13C110D001  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D01__13C000D101  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D01__13C010D101  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D01__13C100D101  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D01__13C110D101  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D01__13C000D101  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D01__13C010D101  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D01__13C100D101  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D01__13C110D101  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D01__13C001D001  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D01__13C011D001  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D01__13C101D001  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D01__13C111D001  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D01__13C001D001  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D01__13C011D001  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D01__13C101D001  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D01__13C111D001  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D01__13C001D101  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D01__13C011D101  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D01__13C101D101  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D01__13C111D101  ...
+  Vr_SerPool2_GlyPool2_13C00D11_13C1D01__13C001D101  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D01__13C011D101  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D01__13C101D101  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D01__13C111D101  +  Vr_SerPool2_GlyPool2_13C00D00_13C0D10__13C000D010  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D10__13C010D010  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D10__13C100D010  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D10__13C110D010  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D10__13C000D010  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D10__13C010D010  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D10__13C100D010  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D10__13C110D010  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D10__13C000D110  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D10__13C010D110  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D10__13C100D110  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D10__13C110D110  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D10__13C000D110  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D10__13C010D110  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D10__13C100D110  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D10__13C110D110  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D10__13C001D010  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D10__13C011D010  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D10__13C101D010  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D10__13C111D010  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D10__13C001D010  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D10__13C011D010  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D10__13C101D010  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D10__13C111D010  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D10__13C001D110  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D10__13C011D110  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D10__13C101D110  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D10__13C111D110  +  Vr_SerPool2_GlyPool2_13C00D11_13C1D10__13C001D110  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D10__13C011D110  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D10__13C101D110  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D10__13C111D110  +  Vr_SerPool2_GlyPool2_13C00D00_13C0D11__13C000D011  +  Vr_SerPool2_GlyPool2_13C01D00_13C0D11__13C010D011  +  Vr_SerPool2_GlyPool2_13C10D00_13C0D11__13C100D011  +  Vr_SerPool2_GlyPool2_13C11D00_13C0D11__13C110D011  +  Vr_SerPool2_GlyPool2_13C00D01_13C0D11__13C000D011  +  Vr_SerPool2_GlyPool2_13C01D01_13C0D11__13C010D011  +  Vr_SerPool2_GlyPool2_13C10D01_13C0D11__13C100D011  +  Vr_SerPool2_GlyPool2_13C11D01_13C0D11__13C110D011  +  Vr_SerPool2_GlyPool2_13C00D10_13C0D11__13C000D111  +  Vr_SerPool2_GlyPool2_13C01D10_13C0D11__13C010D111  +  Vr_SerPool2_GlyPool2_13C10D10_13C0D11__13C100D111  +  Vr_SerPool2_GlyPool2_13C11D10_13C0D11__13C110D111  +  Vr_SerPool2_GlyPool2_13C00D11_13C0D11__13C000D111  +  Vr_SerPool2_GlyPool2_13C01D11_13C0D11__13C010D111  +  Vr_SerPool2_GlyPool2_13C10D11_13C0D11__13C100D111  +  Vr_SerPool2_GlyPool2_13C11D11_13C0D11__13C110D111  +  Vr_SerPool2_GlyPool2_13C00D00_13C1D11__13C001D011  +  Vr_SerPool2_GlyPool2_13C01D00_13C1D11__13C011D011  +  Vr_SerPool2_GlyPool2_13C10D00_13C1D11__13C101D011  +  Vr_SerPool2_GlyPool2_13C11D00_13C1D11__13C111D011  +  Vr_SerPool2_GlyPool2_13C00D01_13C1D11__13C001D011  +  Vr_SerPool2_GlyPool2_13C01D01_13C1D11__13C011D011  +  Vr_SerPool2_GlyPool2_13C10D01_13C1D11__13C101D011  +  Vr_SerPool2_GlyPool2_13C11D01_13C1D11__13C111D011  +  Vr_SerPool2_GlyPool2_13C00D10_13C1D11__13C001D111  +  Vr_SerPool2_GlyPool2_13C01D10_13C1D11__13C011D111  +  Vr_SerPool2_GlyPool2_13C10D10_13C1D11__13C101D111  +  Vr_SerPool2_GlyPool2_13C11D10_13C1D11__13C111D111  +  Vr_SerPool2_GlyPool2_13C00D11_13C1D11__13C001D111  +  Vr_SerPool2_GlyPool2_13C01D11_13C1D11__13C011D111  +  Vr_SerPool2_GlyPool2_13C10D11_13C1D11__13C101D111  +  Vr_SerPool2_GlyPool2_13C11D11_13C1D11__13C111D111  ...
+  Vr_SerPool3_GlyPool3_13C00D00_13C0D00__13C000D000  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D00__13C010D000  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D00__13C100D000  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D00__13C110D000  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D00__13C000D000  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D00__13C010D000  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D00__13C100D000  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D00__13C110D000  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D00__13C000D100  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D00__13C010D100  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D00__13C100D100  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D00__13C110D100  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D00__13C000D100  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D00__13C010D100  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D00__13C100D100  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D00__13C110D100  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D00__13C001D000  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D00__13C011D000  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D00__13C101D000  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D00__13C111D000  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D00__13C001D000  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D00__13C011D000  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D00__13C101D000  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D00__13C111D000  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D00__13C001D100  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D00__13C011D100  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D00__13C101D100  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D00__13C111D100  +  Vr_SerPool3_GlyPool3_13C00D11_13C1D00__13C001D100  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D00__13C011D100  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D00__13C101D100  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D00__13C111D100  +  Vr_SerPool3_GlyPool3_13C00D00_13C0D01__13C000D001  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D01__13C010D001  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D01__13C100D001  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D01__13C110D001  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D01__13C000D001  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D01__13C010D001  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D01__13C100D001  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D01__13C110D001  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D01__13C000D101  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D01__13C010D101  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D01__13C100D101  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D01__13C110D101  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D01__13C000D101  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D01__13C010D101  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D01__13C100D101  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D01__13C110D101  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D01__13C001D001  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D01__13C011D001  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D01__13C101D001  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D01__13C111D001  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D01__13C001D001  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D01__13C011D001  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D01__13C101D001  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D01__13C111D001  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D01__13C001D101  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D01__13C011D101  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D01__13C101D101  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D01__13C111D101  ...
+  Vr_SerPool3_GlyPool3_13C00D11_13C1D01__13C001D101  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D01__13C011D101  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D01__13C101D101  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D01__13C111D101  +  Vr_SerPool3_GlyPool3_13C00D00_13C0D10__13C000D010  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D10__13C010D010  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D10__13C100D010  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D10__13C110D010  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D10__13C000D010  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D10__13C010D010  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D10__13C100D010  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D10__13C110D010  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D10__13C000D110  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D10__13C010D110  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D10__13C100D110  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D10__13C110D110  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D10__13C000D110  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D10__13C010D110  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D10__13C100D110  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D10__13C110D110  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D10__13C001D010  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D10__13C011D010  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D10__13C101D010  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D10__13C111D010  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D10__13C001D010  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D10__13C011D010  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D10__13C101D010  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D10__13C111D010  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D10__13C001D110  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D10__13C011D110  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D10__13C101D110  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D10__13C111D110  +  Vr_SerPool3_GlyPool3_13C00D11_13C1D10__13C001D110  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D10__13C011D110  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D10__13C101D110  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D10__13C111D110  +  Vr_SerPool3_GlyPool3_13C00D00_13C0D11__13C000D011  +  Vr_SerPool3_GlyPool3_13C01D00_13C0D11__13C010D011  +  Vr_SerPool3_GlyPool3_13C10D00_13C0D11__13C100D011  +  Vr_SerPool3_GlyPool3_13C11D00_13C0D11__13C110D011  +  Vr_SerPool3_GlyPool3_13C00D01_13C0D11__13C000D011  +  Vr_SerPool3_GlyPool3_13C01D01_13C0D11__13C010D011  +  Vr_SerPool3_GlyPool3_13C10D01_13C0D11__13C100D011  +  Vr_SerPool3_GlyPool3_13C11D01_13C0D11__13C110D011  +  Vr_SerPool3_GlyPool3_13C00D10_13C0D11__13C000D111  +  Vr_SerPool3_GlyPool3_13C01D10_13C0D11__13C010D111  +  Vr_SerPool3_GlyPool3_13C10D10_13C0D11__13C100D111  +  Vr_SerPool3_GlyPool3_13C11D10_13C0D11__13C110D111  +  Vr_SerPool3_GlyPool3_13C00D11_13C0D11__13C000D111  +  Vr_SerPool3_GlyPool3_13C01D11_13C0D11__13C010D111  +  Vr_SerPool3_GlyPool3_13C10D11_13C0D11__13C100D111  +  Vr_SerPool3_GlyPool3_13C11D11_13C0D11__13C110D111  +  Vr_SerPool3_GlyPool3_13C00D00_13C1D11__13C001D011  +  Vr_SerPool3_GlyPool3_13C01D00_13C1D11__13C011D011  +  Vr_SerPool3_GlyPool3_13C10D00_13C1D11__13C101D011  +  Vr_SerPool3_GlyPool3_13C11D00_13C1D11__13C111D011  +  Vr_SerPool3_GlyPool3_13C00D01_13C1D11__13C001D011  +  Vr_SerPool3_GlyPool3_13C01D01_13C1D11__13C011D011  +  Vr_SerPool3_GlyPool3_13C10D01_13C1D11__13C101D011  +  Vr_SerPool3_GlyPool3_13C11D01_13C1D11__13C111D011  +  Vr_SerPool3_GlyPool3_13C00D10_13C1D11__13C001D111  +  Vr_SerPool3_GlyPool3_13C01D10_13C1D11__13C011D111  +  Vr_SerPool3_GlyPool3_13C10D10_13C1D11__13C101D111  +  Vr_SerPool3_GlyPool3_13C11D10_13C1D11__13C111D111  +  Vr_SerPool3_GlyPool3_13C00D11_13C1D11__13C001D111  +  Vr_SerPool3_GlyPool3_13C01D11_13C1D11__13C011D111  +  Vr_SerPool3_GlyPool3_13C10D11_13C1D11__13C101D111  +  Vr_SerPool3_GlyPool3_13C11D11_13C1D11__13C111D111  ...
+  Vr_GlyPool1_CO2_13C0D00_13C00D00  +  Vr_GlyPool1_CO2_13C1D00_13C01D00  +  Vr_GlyPool1_CO2_13C0D01_13C00D01  +  Vr_GlyPool1_CO2_13C1D01_13C01D01  +  Vr_GlyPool1_CO2_13C0D10_13C00D10  +  Vr_GlyPool1_CO2_13C1D10_13C01D10  +  Vr_GlyPool1_CO2_13C0D11_13C00D11  +  Vr_GlyPool1_CO2_13C1D11_13C01D11  ...
+  Vr_THF_Transport   ...
+  Vr_GlyPool2_CO2_13C0D00_13C00D00  +  Vr_GlyPool2_CO2_13C1D00_13C01D00  +  Vr_GlyPool2_CO2_13C0D01_13C00D01  +  Vr_GlyPool2_CO2_13C1D01_13C01D01  +  Vr_GlyPool2_CO2_13C0D10_13C00D10  +  Vr_GlyPool2_CO2_13C1D10_13C01D10  +  Vr_GlyPool2_CO2_13C0D11_13C00D11  +  Vr_GlyPool2_CO2_13C1D11_13C01D11  ...
+  Vr_GlyPool3_CO2_13C0D00_13C00D00  +  Vr_GlyPool3_CO2_13C1D00_13C01D00  +  Vr_GlyPool3_CO2_13C0D01_13C00D01  +  Vr_GlyPool3_CO2_13C1D01_13C01D01  +  Vr_GlyPool3_CO2_13C0D10_13C00D10  +  Vr_GlyPool3_CO2_13C1D10_13C01D10  +  Vr_GlyPool3_CO2_13C0D11_13C00D11  +  Vr_GlyPool3_CO2_13C1D11_13C01D11  ...
+  Vf_FormylTHF_Formate_Cytoplasm_13C0D0_13C0D0  +  Vf_FormylTHF_Formate_Cytoplasm_13C1D0_13C1D0  +  Vf_FormylTHF_Formate_Cytoplasm_13C0D1_13C0D1  +  Vf_FormylTHF_Formate_Cytoplasm_13C1D1_13C1D1  ...
+  Vf_GAR_FGAR_13C000_13C0D0__13C0000D0  +  Vf_GAR_FGAR_13C001_13C0D0__13C0010D0  +  Vf_GAR_FGAR_13C010_13C0D0__13C0100D0  +  Vf_GAR_FGAR_13C011_13C0D0__13C0110D0  +  Vf_GAR_FGAR_13C100_13C0D0__13C1000D0  +  Vf_GAR_FGAR_13C101_13C0D0__13C1010D0  +  Vf_GAR_FGAR_13C110_13C0D0__13C1100D0  +  Vf_GAR_FGAR_13C111_13C0D0__13C1110D0  +  Vf_GAR_FGAR_13C000_13C1D0__13C0001D0  +  Vf_GAR_FGAR_13C001_13C1D0__13C0011D0  +  Vf_GAR_FGAR_13C010_13C1D0__13C0101D0  +  Vf_GAR_FGAR_13C011_13C1D0__13C0111D0  +  Vf_GAR_FGAR_13C100_13C1D0__13C1001D0  +  Vf_GAR_FGAR_13C101_13C1D0__13C1011D0  +  Vf_GAR_FGAR_13C110_13C1D0__13C1101D0  +  Vf_GAR_FGAR_13C111_13C1D0__13C1111D0  +  Vf_GAR_FGAR_13C000_13C0D1__13C0000D1  +  Vf_GAR_FGAR_13C001_13C0D1__13C0010D1  +  Vf_GAR_FGAR_13C010_13C0D1__13C0100D1  +  Vf_GAR_FGAR_13C011_13C0D1__13C0110D1  +  Vf_GAR_FGAR_13C100_13C0D1__13C1000D1  +  Vf_GAR_FGAR_13C101_13C0D1__13C1010D1  +  Vf_GAR_FGAR_13C110_13C0D1__13C1100D1  +  Vf_GAR_FGAR_13C111_13C0D1__13C1110D1  +  Vf_GAR_FGAR_13C000_13C1D1__13C0001D1  +  Vf_GAR_FGAR_13C001_13C1D1__13C0011D1  +  Vf_GAR_FGAR_13C010_13C1D1__13C0101D1  +  Vf_GAR_FGAR_13C011_13C1D1__13C0111D1  +  Vf_GAR_FGAR_13C100_13C1D1__13C1001D1  +  Vf_GAR_FGAR_13C101_13C1D1__13C1011D1  +  Vf_GAR_FGAR_13C110_13C1D1__13C1101D1  +  Vf_GAR_FGAR_13C111_13C1D1__13C1111D1  ...
+  Vf_FGAR_AMP_13C0000D0_13C0D0__13C00000D00  +  Vf_FGAR_AMP_13C0001D0_13C0D0__13C00010D00  +  Vf_FGAR_AMP_13C0010D0_13C0D0__13C00100D00  +  Vf_FGAR_AMP_13C0011D0_13C0D0__13C00110D00  +  Vf_FGAR_AMP_13C0100D0_13C0D0__13C01000D00  +  Vf_FGAR_AMP_13C0101D0_13C0D0__13C01010D00  +  Vf_FGAR_AMP_13C0110D0_13C0D0__13C01100D00  +  Vf_FGAR_AMP_13C0111D0_13C0D0__13C01110D00  +  Vf_FGAR_AMP_13C1000D0_13C0D0__13C10000D00  +  Vf_FGAR_AMP_13C1001D0_13C0D0__13C10010D00  +  Vf_FGAR_AMP_13C1010D0_13C0D0__13C10100D00  +  Vf_FGAR_AMP_13C1011D0_13C0D0__13C10110D00  +  Vf_FGAR_AMP_13C1100D0_13C0D0__13C11000D00  +  Vf_FGAR_AMP_13C1101D0_13C0D0__13C11010D00  +  Vf_FGAR_AMP_13C1110D0_13C0D0__13C11100D00  +  Vf_FGAR_AMP_13C1111D0_13C0D0__13C11110D00  +  Vf_FGAR_AMP_13C0000D1_13C0D0__13C00000D10  +  Vf_FGAR_AMP_13C0001D1_13C0D0__13C00010D10  +  Vf_FGAR_AMP_13C0010D1_13C0D0__13C00100D10  +  Vf_FGAR_AMP_13C0011D1_13C0D0__13C00110D10  +  Vf_FGAR_AMP_13C0100D1_13C0D0__13C01000D10  +  Vf_FGAR_AMP_13C0101D1_13C0D0__13C01010D10  +  Vf_FGAR_AMP_13C0110D1_13C0D0__13C01100D10  +  Vf_FGAR_AMP_13C0111D1_13C0D0__13C01110D10  +  Vf_FGAR_AMP_13C1000D1_13C0D0__13C10000D10  +  Vf_FGAR_AMP_13C1001D1_13C0D0__13C10010D10  +  Vf_FGAR_AMP_13C1010D1_13C0D0__13C10100D10  +  Vf_FGAR_AMP_13C1011D1_13C0D0__13C10110D10  +  Vf_FGAR_AMP_13C1100D1_13C0D0__13C11000D10  +  Vf_FGAR_AMP_13C1101D1_13C0D0__13C11010D10  +  Vf_FGAR_AMP_13C1110D1_13C0D0__13C11100D10  +  Vf_FGAR_AMP_13C1111D1_13C0D0__13C11110D10  +  Vf_FGAR_AMP_13C0000D0_13C1D0__13C00001D00  +  Vf_FGAR_AMP_13C0001D0_13C1D0__13C00011D00  +  Vf_FGAR_AMP_13C0010D0_13C1D0__13C00101D00  +  Vf_FGAR_AMP_13C0011D0_13C1D0__13C00111D00  +  Vf_FGAR_AMP_13C0100D0_13C1D0__13C01001D00  +  Vf_FGAR_AMP_13C0101D0_13C1D0__13C01011D00  +  Vf_FGAR_AMP_13C0110D0_13C1D0__13C01101D00  +  Vf_FGAR_AMP_13C0111D0_13C1D0__13C01111D00  +  Vf_FGAR_AMP_13C1000D0_13C1D0__13C10001D00  +  Vf_FGAR_AMP_13C1001D0_13C1D0__13C10011D00  +  Vf_FGAR_AMP_13C1010D0_13C1D0__13C10101D00  +  Vf_FGAR_AMP_13C1011D0_13C1D0__13C10111D00  +  Vf_FGAR_AMP_13C1100D0_13C1D0__13C11001D00  +  Vf_FGAR_AMP_13C1101D0_13C1D0__13C11011D00  +  Vf_FGAR_AMP_13C1110D0_13C1D0__13C11101D00  +  Vf_FGAR_AMP_13C1111D0_13C1D0__13C11111D00  +  Vf_FGAR_AMP_13C0000D1_13C1D0__13C00001D10  +  Vf_FGAR_AMP_13C0001D1_13C1D0__13C00011D10  +  Vf_FGAR_AMP_13C0010D1_13C1D0__13C00101D10  +  Vf_FGAR_AMP_13C0011D1_13C1D0__13C00111D10  +  Vf_FGAR_AMP_13C0100D1_13C1D0__13C01001D10  +  Vf_FGAR_AMP_13C0101D1_13C1D0__13C01011D10  +  Vf_FGAR_AMP_13C0110D1_13C1D0__13C01101D10  +  Vf_FGAR_AMP_13C0111D1_13C1D0__13C01111D10  +  Vf_FGAR_AMP_13C1000D1_13C1D0__13C10001D10  +  Vf_FGAR_AMP_13C1001D1_13C1D0__13C10011D10  +  Vf_FGAR_AMP_13C1010D1_13C1D0__13C10101D10  +  Vf_FGAR_AMP_13C1011D1_13C1D0__13C10111D10  +  Vf_FGAR_AMP_13C1100D1_13C1D0__13C11001D10  +  Vf_FGAR_AMP_13C1101D1_13C1D0__13C11011D10  +  Vf_FGAR_AMP_13C1110D1_13C1D0__13C11101D10  +  Vf_FGAR_AMP_13C1111D1_13C1D0__13C11111D10  +  Vf_FGAR_AMP_13C0000D0_13C0D1__13C00000D01  +  Vf_FGAR_AMP_13C0001D0_13C0D1__13C00010D01  +  Vf_FGAR_AMP_13C0010D0_13C0D1__13C00100D01  +  Vf_FGAR_AMP_13C0011D0_13C0D1__13C00110D01  +  Vf_FGAR_AMP_13C0100D0_13C0D1__13C01000D01  +  Vf_FGAR_AMP_13C0101D0_13C0D1__13C01010D01  +  Vf_FGAR_AMP_13C0110D0_13C0D1__13C01100D01  +  Vf_FGAR_AMP_13C0111D0_13C0D1__13C01110D01  ...
+  Vf_FGAR_AMP_13C1000D0_13C0D1__13C10000D01  +  Vf_FGAR_AMP_13C1001D0_13C0D1__13C10010D01  +  Vf_FGAR_AMP_13C1010D0_13C0D1__13C10100D01  +  Vf_FGAR_AMP_13C1011D0_13C0D1__13C10110D01  +  Vf_FGAR_AMP_13C1100D0_13C0D1__13C11000D01  +  Vf_FGAR_AMP_13C1101D0_13C0D1__13C11010D01  +  Vf_FGAR_AMP_13C1110D0_13C0D1__13C11100D01  +  Vf_FGAR_AMP_13C1111D0_13C0D1__13C11110D01  +  Vf_FGAR_AMP_13C0000D1_13C0D1__13C00000D11  +  Vf_FGAR_AMP_13C0001D1_13C0D1__13C00010D11  +  Vf_FGAR_AMP_13C0010D1_13C0D1__13C00100D11  +  Vf_FGAR_AMP_13C0011D1_13C0D1__13C00110D11  +  Vf_FGAR_AMP_13C0100D1_13C0D1__13C01000D11  +  Vf_FGAR_AMP_13C0101D1_13C0D1__13C01010D11  +  Vf_FGAR_AMP_13C0110D1_13C0D1__13C01100D11  +  Vf_FGAR_AMP_13C0111D1_13C0D1__13C01110D11  +  Vf_FGAR_AMP_13C1000D1_13C0D1__13C10000D11  +  Vf_FGAR_AMP_13C1001D1_13C0D1__13C10010D11  +  Vf_FGAR_AMP_13C1010D1_13C0D1__13C10100D11  +  Vf_FGAR_AMP_13C1011D1_13C0D1__13C10110D11  +  Vf_FGAR_AMP_13C1100D1_13C0D1__13C11000D11  +  Vf_FGAR_AMP_13C1101D1_13C0D1__13C11010D11  +  Vf_FGAR_AMP_13C1110D1_13C0D1__13C11100D11  +  Vf_FGAR_AMP_13C1111D1_13C0D1__13C11110D11  +  Vf_FGAR_AMP_13C0000D0_13C1D1__13C00001D01  +  Vf_FGAR_AMP_13C0001D0_13C1D1__13C00011D01  +  Vf_FGAR_AMP_13C0010D0_13C1D1__13C00101D01  +  Vf_FGAR_AMP_13C0011D0_13C1D1__13C00111D01  +  Vf_FGAR_AMP_13C0100D0_13C1D1__13C01001D01  +  Vf_FGAR_AMP_13C0101D0_13C1D1__13C01011D01  +  Vf_FGAR_AMP_13C0110D0_13C1D1__13C01101D01  +  Vf_FGAR_AMP_13C0111D0_13C1D1__13C01111D01  +  Vf_FGAR_AMP_13C1000D0_13C1D1__13C10001D01  +  Vf_FGAR_AMP_13C1001D0_13C1D1__13C10011D01  +  Vf_FGAR_AMP_13C1010D0_13C1D1__13C10101D01  +  Vf_FGAR_AMP_13C1011D0_13C1D1__13C10111D01  +  Vf_FGAR_AMP_13C1100D0_13C1D1__13C11001D01  +  Vf_FGAR_AMP_13C1101D0_13C1D1__13C11011D01  +  Vf_FGAR_AMP_13C1110D0_13C1D1__13C11101D01  +  Vf_FGAR_AMP_13C1111D0_13C1D1__13C11111D01  +  Vf_FGAR_AMP_13C0000D1_13C1D1__13C00001D11  +  Vf_FGAR_AMP_13C0001D1_13C1D1__13C00011D11  +  Vf_FGAR_AMP_13C0010D1_13C1D1__13C00101D11  +  Vf_FGAR_AMP_13C0011D1_13C1D1__13C00111D11  +  Vf_FGAR_AMP_13C0100D1_13C1D1__13C01001D11  +  Vf_FGAR_AMP_13C0101D1_13C1D1__13C01011D11  +  Vf_FGAR_AMP_13C0110D1_13C1D1__13C01101D11  +  Vf_FGAR_AMP_13C0111D1_13C1D1__13C01111D11  +  Vf_FGAR_AMP_13C1000D1_13C1D1__13C10001D11  +  Vf_FGAR_AMP_13C1001D1_13C1D1__13C10011D11  +  Vf_FGAR_AMP_13C1010D1_13C1D1__13C10101D11  +  Vf_FGAR_AMP_13C1011D1_13C1D1__13C10111D11  +  Vf_FGAR_AMP_13C1100D1_13C1D1__13C11001D11  +  Vf_FGAR_AMP_13C1101D1_13C1D1__13C11011D11  +  Vf_FGAR_AMP_13C1110D1_13C1D1__13C11101D11  +  Vf_FGAR_AMP_13C1111D1_13C1D1__13C11111D11;

%  cTHF 


    
    % mMethyleneTHF
    %  dxdt(1, 1) = (mMethyleneTHF) = 
    % - Vf_MethyleneTHF_FormylTHF_Mitochon 
    % - Vr_MethyleneTHF_Transport 
    % - Vr_SerPoolMitochon_GlyPoolMitochon 
    % - Vr_GlyPoolMitochon_CO2 
    % + Vr_MethyleneTHF_FormylTHF_Mitochon 
    % + Vf_MethyleneTHF_Transport 
    % + Vf_SerPoolMitochon_GlyPoolMitochon 
    % + Vf_GlyPoolMitochon_CO2;
    
		dxdt(422, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D00_13C0D0  -  Vr_MethyleneTHF_Transport_13C0D00_13C0D00  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D00__13C000D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D00__13C010D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D00__13C100D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D00__13C110D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D00__13C000D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D00__13C010D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D00__13C100D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D00__13C110D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D00__13C000D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D00__13C010D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D00__13C100D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D00__13C110D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D00__13C000D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D00__13C010D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D00__13C100D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D00__13C110D100  -  Vr_GlyPoolMitochon_CO2_13C0D00_13C00D00  +  Vr_MethyleneTHF_FormylTHF_Mitochon_13C0D0_13C0D00  +  Vf_MethyleneTHF_Transport_13C0D00_13C0D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D000_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D000_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D000_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D000_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D100_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D100_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D100_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D100_13C11D10  +  Vf_GlyPoolMitochon_CO2_13C00D00_13C0D00  +  Vf_GlyPoolMitochon_CO2_13C10D00_13C0D00;								%  mMethyleneTHF_13C0D00
		dxdt(423, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D00_13C1D0  -  Vr_MethyleneTHF_Transport_13C1D00_13C1D00  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D00__13C001D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D00__13C011D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D00__13C101D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D00__13C111D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D00__13C001D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D00__13C011D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D00__13C101D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D00__13C111D000  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D00__13C001D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D00__13C011D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D00__13C101D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D00__13C111D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D00__13C001D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D00__13C011D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D00__13C101D100  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D00__13C111D100  -  Vr_GlyPoolMitochon_CO2_13C1D00_13C01D00  +  Vr_MethyleneTHF_FormylTHF_Mitochon_13C1D0_13C1D00  +  Vf_MethyleneTHF_Transport_13C1D00_13C1D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D100_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D100_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D100_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D100_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D000_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D000_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D000_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D000_13C11D00  +  Vf_GlyPoolMitochon_CO2_13C01D00_13C1D00  +  Vf_GlyPoolMitochon_CO2_13C11D00_13C1D00;         			%  mMethyleneTHF_13C1D00
		dxdt(424, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D10_13C0D1  -  Vr_MethyleneTHF_Transport_13C0D01_13C0D01  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D01__13C000D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D01__13C010D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D01__13C100D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D01__13C110D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D01__13C000D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D01__13C010D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D01__13C100D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D01__13C110D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D01__13C000D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D01__13C010D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D01__13C100D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D01__13C110D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D01__13C000D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D01__13C010D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D01__13C100D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D01__13C110D101  -  Vr_GlyPoolMitochon_CO2_13C0D01_13C00D01                                                        +  Vf_MethyleneTHF_Transport_13C0D01_13C0D01  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D001_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D001_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D001_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D001_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D101_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D101_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D101_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D101_13C11D10  +  Vf_GlyPoolMitochon_CO2_13C00D01_13C0D01  +  Vf_GlyPoolMitochon_CO2_13C10D01_13C0D01;         			%  mMethyleneTHF_13C0D01
		dxdt(425, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D10_13C1D1  -  Vr_MethyleneTHF_Transport_13C1D01_13C1D01  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D01__13C001D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D01__13C011D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D01__13C101D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D01__13C111D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D01__13C001D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D01__13C011D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D01__13C101D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D01__13C111D001  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D01__13C001D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D01__13C011D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D01__13C101D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D01__13C111D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D01__13C001D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D01__13C011D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D01__13C101D101  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D01__13C111D101  -  Vr_GlyPoolMitochon_CO2_13C1D01_13C01D01                                                        +  Vf_MethyleneTHF_Transport_13C1D01_13C1D01  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D101_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D101_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D101_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D101_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D001_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D001_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D001_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D001_13C11D00  +  Vf_GlyPoolMitochon_CO2_13C01D01_13C1D01  +  Vf_GlyPoolMitochon_CO2_13C11D01_13C1D01;         			%  mMethyleneTHF_13C1D01
		dxdt(426, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D01_13C0D0  -  Vr_MethyleneTHF_Transport_13C0D10_13C0D10  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D10__13C000D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D10__13C010D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D10__13C100D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D10__13C110D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D10__13C000D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D10__13C010D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D10__13C100D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D10__13C110D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D10__13C000D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D10__13C010D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D10__13C100D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D10__13C110D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D10__13C000D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D10__13C010D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D10__13C100D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D10__13C110D110  -  Vr_GlyPoolMitochon_CO2_13C0D10_13C00D10  +  Vr_MethyleneTHF_FormylTHF_Mitochon_13C0D1_13C0D10  +  Vf_MethyleneTHF_Transport_13C0D10_13C0D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D010_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D010_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D010_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D010_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D110_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D110_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D110_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D110_13C11D10  +  Vf_GlyPoolMitochon_CO2_13C00D10_13C0D10  +  Vf_GlyPoolMitochon_CO2_13C10D10_13C0D10;         			%  mMethyleneTHF_13C0D10
		dxdt(427, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D01_13C1D0  -  Vr_MethyleneTHF_Transport_13C1D10_13C1D10  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D10__13C001D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D10__13C011D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D10__13C101D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D10__13C111D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D10__13C001D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D10__13C011D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D10__13C101D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D10__13C111D010  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D10__13C001D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D10__13C011D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D10__13C101D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D10__13C111D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D10__13C001D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D10__13C011D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D10__13C101D110  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D10__13C111D110  -  Vr_GlyPoolMitochon_CO2_13C1D10_13C01D10  +  Vr_MethyleneTHF_FormylTHF_Mitochon_13C1D1_13C1D10  +  Vf_MethyleneTHF_Transport_13C1D10_13C1D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D110_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D110_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D110_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D110_13C11D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D010_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D010_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D010_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D010_13C11D00  +  Vf_GlyPoolMitochon_CO2_13C01D10_13C1D10  +  Vf_GlyPoolMitochon_CO2_13C11D10_13C1D10;         			%  mMethyleneTHF_13C1D10
		dxdt(428, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D11_13C0D1  -  Vr_MethyleneTHF_Transport_13C0D11_13C0D11  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D11__13C000D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D11__13C010D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D11__13C100D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D11__13C110D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D11__13C000D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D11__13C010D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D11__13C100D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D11__13C110D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D11__13C000D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D11__13C010D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D11__13C100D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D11__13C110D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D11__13C000D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D11__13C010D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D11__13C100D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D11__13C110D111  -  Vr_GlyPoolMitochon_CO2_13C0D11_13C00D11                                                        +  Vf_MethyleneTHF_Transport_13C0D11_13C0D11  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D011_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D011_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D011_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D011_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D111_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D111_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D111_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D111_13C11D10  +  Vf_GlyPoolMitochon_CO2_13C00D11_13C0D11  +  Vf_GlyPoolMitochon_CO2_13C10D11_13C0D11;         			%  mMethyleneTHF_13C0D11
		dxdt(429, 1) = - Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D11_13C1D1  -  Vr_MethyleneTHF_Transport_13C1D11_13C1D11  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D11__13C001D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D11__13C011D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D11__13C101D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D11__13C111D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D11__13C001D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D11__13C011D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D11__13C101D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D11__13C111D011  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D11__13C001D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D11__13C011D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D11__13C101D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D11__13C111D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D11__13C001D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D11__13C011D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D11__13C101D111  -  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D11__13C111D111  -  Vr_GlyPoolMitochon_CO2_13C1D11_13C01D11                                                        +  Vf_MethyleneTHF_Transport_13C1D11_13C1D11  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D011_13C00D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D011_13C01D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D011_13C10D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D011_13C11D00  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D111_13C00D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D111_13C01D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D111_13C10D10  +  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D111_13C11D10  +  Vf_GlyPoolMitochon_CO2_13C01D11_13C1D11  +  Vf_GlyPoolMitochon_CO2_13C11D11_13C1D11;         			%  mMethyleneTHF_13C1D11


    
    % mFormylTHF
    %  dxdt(1, 1) = (mFormylTHF) = 
    % - Vf_FormylTHF_Formate_Mitochon 
    % - Vr_FormylTHF_Transport 
    % - Vr_MethyleneTHF_FormylTHF_Mitochon 
    % + Vr_FormylTHF_Formate_Mitochon 
    % + Vf_FormylTHF_Transport 
    % + Vf_MethyleneTHF_FormylTHF_Mitochon;

		dxdt(430, 1) = - Vf_FormylTHF_Formate_Mitochon_13C0D0_13C0D0  -  Vr_FormylTHF_Transport_13C0D0_13C0D0  -  Vr_MethyleneTHF_FormylTHF_Mitochon_13C0D0_13C0D00  +  Vr_FormylTHF_Formate_Mitochon_13C0D0_13C0D0  +  Vf_FormylTHF_Transport_13C0D0_13C0D0  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D00_13C0D0  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D01_13C0D0;							%  mFormylTHF_13C0D0 
		dxdt(431, 1) = - Vf_FormylTHF_Formate_Mitochon_13C1D0_13C1D0  -  Vr_FormylTHF_Transport_13C1D0_13C1D0  -  Vr_MethyleneTHF_FormylTHF_Mitochon_13C1D0_13C1D00  +  Vr_FormylTHF_Formate_Mitochon_13C1D0_13C1D0  +  Vf_FormylTHF_Transport_13C1D0_13C1D0  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D00_13C1D0  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D01_13C1D0;             %  mFormylTHF_13C1D0 
		dxdt(432, 1) = - Vf_FormylTHF_Formate_Mitochon_13C0D1_13C0D1  -  Vr_FormylTHF_Transport_13C0D1_13C0D1  -  Vr_MethyleneTHF_FormylTHF_Mitochon_13C0D1_13C0D10  +  Vr_FormylTHF_Formate_Mitochon_13C0D1_13C0D1  +  Vf_FormylTHF_Transport_13C0D1_13C0D1  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D10_13C0D1  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C0D11_13C0D1;             %  mFormylTHF_13C0D1 
		dxdt(433, 1) = - Vf_FormylTHF_Formate_Mitochon_13C1D1_13C1D1  -  Vr_FormylTHF_Transport_13C1D1_13C1D1  -  Vr_MethyleneTHF_FormylTHF_Mitochon_13C1D1_13C1D10  +  Vr_FormylTHF_Formate_Mitochon_13C1D1_13C1D1  +  Vf_FormylTHF_Transport_13C1D1_13C1D1  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D10_13C1D1  +  Vf_MethyleneTHF_FormylTHF_Mitochon_13C1D11_13C1D1;             %  mFormylTHF_13C1D1 

    
    
    % mFormate
    %  dxdt(1, 1) = (mFormate) = 
    % - Vr_Formate_Transport 
    % - Vr_FormylTHF_Formate_Mitochon 
    % + Vf_Formate_Transport 
    % + Vf_FormylTHF_Formate_Mitochon;
    
		dxdt(434, 1) =  - Vr_Formate_Transport_13C0D0  -  Vr_FormylTHF_Formate_Mitochon_13C0D0_13C0D0  +  Vf_Formate_Transport_13C0D0  +  Vf_FormylTHF_Formate_Mitochon_13C0D0_13C0D0;				%  mFormate_13C0D0 
		dxdt(435, 1) =  - Vr_Formate_Transport_13C1D0  -  Vr_FormylTHF_Formate_Mitochon_13C1D0_13C1D0  +  Vf_Formate_Transport_13C1D0  +  Vf_FormylTHF_Formate_Mitochon_13C1D0_13C1D0;        %  mFormate_13C1D0 
		dxdt(436, 1) =  - Vr_Formate_Transport_13C0D1  -  Vr_FormylTHF_Formate_Mitochon_13C0D1_13C0D1  +  Vf_Formate_Transport_13C0D1  +  Vf_FormylTHF_Formate_Mitochon_13C0D1_13C0D1;        %  mFormate_13C0D1 
		dxdt(437, 1) =  - Vr_Formate_Transport_13C1D1  -  Vr_FormylTHF_Formate_Mitochon_13C1D1_13C1D1  +  Vf_Formate_Transport_13C1D1  +  Vf_FormylTHF_Formate_Mitochon_13C1D1_13C1D1;        %  mFormate_13C1D1 


    

    % mTHF
    %  dxdt(1, 1) = (mTHF) = 
    % - Vf_SerPoolMitochon_GlyPoolMitochon     
    % - Vf_GlyPoolMitochon_CO2     
    % - Vr_FormylTHF_Formate_Mitochon     
    % - Vr_THF_Transport 
    % + Vr_SerPoolMitochon_GlyPoolMitochon 
    % + Vr_GlyPoolMitochon_CO2 
    % + Vf_FormylTHF_Formate_Mitochon 
    % + Vf_THF_Transport;
    
		dxdt(438, 1) =  - Vf_SerPoolMitochon_GlyPoolMitochon_13C000D000_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D000_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D000_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D000_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D000_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D000_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D000_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D000_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D001_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D001_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D001_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D001_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D001_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D001_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D001_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D001_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D010_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D010_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D010_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D010_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D010_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D010_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D010_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D010_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D011_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D011_13C00D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D011_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D011_13C01D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D011_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D011_13C10D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D011_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D011_13C11D00  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D100_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D100_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D100_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D100_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D100_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D100_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D100_13C11D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D100_13C11D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D101_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D101_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D101_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D101_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D101_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D101_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D101_13C11D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D101_13C11D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D110_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D110_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D110_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D110_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D110_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D110_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D110_13C11D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D110_13C11D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C000D111_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C001D111_13C00D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C010D111_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C011D111_13C01D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C100D111_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C101D111_13C10D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C110D111_13C11D10  -  Vf_SerPoolMitochon_GlyPoolMitochon_13C111D111_13C11D10  ...
-  Vf_GlyPoolMitochon_CO2_13C00D00_13C0D00  -  Vf_GlyPoolMitochon_CO2_13C01D00_13C1D00  -  Vf_GlyPoolMitochon_CO2_13C10D00_13C0D00  -  Vf_GlyPoolMitochon_CO2_13C11D00_13C1D00  -  Vf_GlyPoolMitochon_CO2_13C00D01_13C0D01  -  Vf_GlyPoolMitochon_CO2_13C01D01_13C1D01  -  Vf_GlyPoolMitochon_CO2_13C10D01_13C0D01  -  Vf_GlyPoolMitochon_CO2_13C11D01_13C1D01  -  Vf_GlyPoolMitochon_CO2_13C00D10_13C0D10  -  Vf_GlyPoolMitochon_CO2_13C01D10_13C1D10  -  Vf_GlyPoolMitochon_CO2_13C10D10_13C0D10  -  Vf_GlyPoolMitochon_CO2_13C11D10_13C1D10  -  Vf_GlyPoolMitochon_CO2_13C00D11_13C0D11  -  Vf_GlyPoolMitochon_CO2_13C01D11_13C1D11  -  Vf_GlyPoolMitochon_CO2_13C10D11_13C0D11  -  Vf_GlyPoolMitochon_CO2_13C11D11_13C1D11  ...
-  Vr_FormylTHF_Formate_Mitochon_13C0D0_13C0D0  -  Vr_FormylTHF_Formate_Mitochon_13C1D0_13C1D0  -  Vr_FormylTHF_Formate_Mitochon_13C0D1_13C0D1  -  Vr_FormylTHF_Formate_Mitochon_13C1D1_13C1D1  ...
-  Vr_THF_Transport  ...
+  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D00__13C000D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D00__13C010D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D00__13C100D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D00__13C110D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D00__13C000D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D00__13C010D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D00__13C100D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D00__13C110D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D00__13C000D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D00__13C010D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D00__13C100D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D00__13C110D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D00__13C000D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D00__13C010D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D00__13C100D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D00__13C110D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D00__13C001D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D00__13C011D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D00__13C101D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D00__13C111D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D00__13C001D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D00__13C011D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D00__13C101D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D00__13C111D000  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D00__13C001D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D00__13C011D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D00__13C101D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D00__13C111D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C1D00__13C001D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C1D00__13C011D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C1D00__13C101D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C1D00__13C111D100  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C0D01__13C000D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C0D01__13C010D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C0D01__13C100D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C0D01__13C110D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C0D01__13C000D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C0D01__13C010D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C0D01__13C100D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C0D01__13C110D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C0D01__13C000D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C0D01__13C010D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C0D01__13C100D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C0D01__13C110D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D11_13C0D01__13C000D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D11_13C0D01__13C010D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D11_13C0D01__13C100D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D11_13C0D01__13C110D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D00_13C1D01__13C001D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D00_13C1D01__13C011D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D00_13C1D01__13C101D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D00_13C1D01__13C111D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D01_13C1D01__13C001D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D01_13C1D01__13C011D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D01_13C1D01__13C101D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D01_13C1D01__13C111D001  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C00D10_13C1D01__13C001D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C01D10_13C1D01__13C011D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C10D10_13C1D01__13C101D101  +  Vr_SerPoolMitochon_GlyPoolMitochon_13C11D10_13C1D01__13C111D101  ...
+  Vr_GlyPoolMitochon_CO2_13C0D00_13C00D00  +  Vr_GlyPoolMitochon_CO2_13C1D00_13C01D00  +  Vr_GlyPoolMitochon_CO2_13C0D01_13C00D01  +  Vr_GlyPoolMitochon_CO2_13C1D01_13C01D01  +  Vr_GlyPoolMitochon_CO2_13C0D10_13C00D10  +  Vr_GlyPoolMitochon_CO2_13C1D10_13C01D10  +  Vr_GlyPoolMitochon_CO2_13C0D11_13C00D11  +  Vr_GlyPoolMitochon_CO2_13C1D11_13C01D11  ...
+  Vf_FormylTHF_Formate_Mitochon_13C0D0_13C0D0  +  Vf_FormylTHF_Formate_Mitochon_13C1D0_13C1D0  +  Vf_FormylTHF_Formate_Mitochon_13C0D1_13C0D1  +  Vf_FormylTHF_Formate_Mitochon_13C1D1_13C1D1  ...
+  Vf_THF_Transport;

%  mTHF


  
    % cPRPP
    %  dxdt(1, 1) = (cPRPP) = 
    % - Vf_PRPP_GAR_GlyPool2 
    % - Vf_PRPP_GAR_GlyPool1 
    % - Vf_PRPP_GAR_GlyPool3 
    % + Vf_GLC_PRPP;

		dxdt(439, 1) =  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D00__13C000  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D00__13C001  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D00__13C010  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D00__13C011  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D01__13C000  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D01__13C001  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D01__13C010  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D01__13C011  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D10__13C000  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D10__13C001  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D10__13C010  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D10__13C011  -  Vf_PRPP_GAR_GlyPool2_13C0_13C00D11__13C000  -  Vf_PRPP_GAR_GlyPool2_13C0_13C01D11__13C001  -  Vf_PRPP_GAR_GlyPool2_13C0_13C10D11__13C010  -  Vf_PRPP_GAR_GlyPool2_13C0_13C11D11__13C011  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D00__13C000  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D00__13C001  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D00__13C010  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D00__13C011  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D01__13C000  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D01__13C001  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D01__13C010  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D01__13C011  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D10__13C000  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D10__13C001  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D10__13C010  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D10__13C011  -  Vf_PRPP_GAR_GlyPool1_13C0_13C00D11__13C000  -  Vf_PRPP_GAR_GlyPool1_13C0_13C01D11__13C001  -  Vf_PRPP_GAR_GlyPool1_13C0_13C10D11__13C010  -  Vf_PRPP_GAR_GlyPool1_13C0_13C11D11__13C011  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D00__13C000  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D00__13C001  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D00__13C010  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D00__13C011  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D01__13C000  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D01__13C001  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D01__13C010  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D01__13C011  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D10__13C000  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D10__13C001  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D10__13C010  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D10__13C011  -  Vf_PRPP_GAR_GlyPool3_13C0_13C00D11__13C000  -  Vf_PRPP_GAR_GlyPool3_13C0_13C01D11__13C001  -  Vf_PRPP_GAR_GlyPool3_13C0_13C10D11__13C010  -  Vf_PRPP_GAR_GlyPool3_13C0_13C11D11__13C011  +  Vf_GLC_PRPP_13C0_13C0;					 %  cPRPP_13C0
		dxdt(440, 1) =  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D00__13C100  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D00__13C101  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D00__13C110  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D00__13C111  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D01__13C100  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D01__13C101  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D01__13C110  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D01__13C111  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D10__13C100  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D10__13C101  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D10__13C110  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D10__13C111  -  Vf_PRPP_GAR_GlyPool2_13C1_13C00D11__13C100  -  Vf_PRPP_GAR_GlyPool2_13C1_13C01D11__13C101  -  Vf_PRPP_GAR_GlyPool2_13C1_13C10D11__13C110  -  Vf_PRPP_GAR_GlyPool2_13C1_13C11D11__13C111  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D00__13C100  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D00__13C101  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D00__13C110  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D00__13C111  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D01__13C100  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D01__13C101  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D01__13C110  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D01__13C111  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D10__13C100  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D10__13C101  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D10__13C110  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D10__13C111  -  Vf_PRPP_GAR_GlyPool1_13C1_13C00D11__13C100  -  Vf_PRPP_GAR_GlyPool1_13C1_13C01D11__13C101  -  Vf_PRPP_GAR_GlyPool1_13C1_13C10D11__13C110  -  Vf_PRPP_GAR_GlyPool1_13C1_13C11D11__13C111  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D00__13C100  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D00__13C101  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D00__13C110  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D00__13C111  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D01__13C100  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D01__13C101  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D01__13C110  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D01__13C111  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D10__13C100  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D10__13C101  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D10__13C110  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D10__13C111  -  Vf_PRPP_GAR_GlyPool3_13C1_13C00D11__13C100  -  Vf_PRPP_GAR_GlyPool3_13C1_13C01D11__13C101  -  Vf_PRPP_GAR_GlyPool3_13C1_13C10D11__13C110  -  Vf_PRPP_GAR_GlyPool3_13C1_13C11D11__13C111  +  Vf_GLC_PRPP_13C1_13C1;          %  cPRPP_13C1



    % cGAR
    %  dxdt(1, 1) = (cGAR) = 
    % - Vf_GAR_FGAR 
    % + Vf_PRPP_GAR_GlyPool1 
    % + Vf_PRPP_GAR_GlyPool2;
    % + Vf_PRPP_GAR_GlyPool3;
    
		dxdt(441, 1) = -  Vf_GAR_FGAR_13C000_13C0D0__13C0000D0  -  Vf_GAR_FGAR_13C000_13C1D0__13C0001D0  -  Vf_GAR_FGAR_13C000_13C0D1__13C0000D1  -  Vf_GAR_FGAR_13C000_13C1D1__13C0001D1  +  Vf_PRPP_GAR_GlyPool1_13C0_13C00D00__13C000  +  Vf_PRPP_GAR_GlyPool1_13C0_13C00D01__13C000  +  Vf_PRPP_GAR_GlyPool1_13C0_13C00D10__13C000  +  Vf_PRPP_GAR_GlyPool1_13C0_13C00D11__13C000  +  Vf_PRPP_GAR_GlyPool2_13C0_13C00D00__13C000  +  Vf_PRPP_GAR_GlyPool2_13C0_13C00D01__13C000  +  Vf_PRPP_GAR_GlyPool2_13C0_13C00D10__13C000  +  Vf_PRPP_GAR_GlyPool2_13C0_13C00D11__13C000  +  Vf_PRPP_GAR_GlyPool3_13C0_13C00D00__13C000  +  Vf_PRPP_GAR_GlyPool3_13C0_13C00D01__13C000  +  Vf_PRPP_GAR_GlyPool3_13C0_13C00D10__13C000  +  Vf_PRPP_GAR_GlyPool3_13C0_13C00D11__13C000;					%  cGAR_13C000 
		dxdt(442, 1) = -  Vf_GAR_FGAR_13C001_13C0D0__13C0010D0  -  Vf_GAR_FGAR_13C001_13C1D0__13C0011D0  -  Vf_GAR_FGAR_13C001_13C0D1__13C0010D1  -  Vf_GAR_FGAR_13C001_13C1D1__13C0011D1  +  Vf_PRPP_GAR_GlyPool1_13C0_13C01D00__13C001  +  Vf_PRPP_GAR_GlyPool1_13C0_13C01D01__13C001  +  Vf_PRPP_GAR_GlyPool1_13C0_13C01D10__13C001  +  Vf_PRPP_GAR_GlyPool1_13C0_13C01D11__13C001  +  Vf_PRPP_GAR_GlyPool2_13C0_13C01D00__13C001  +  Vf_PRPP_GAR_GlyPool2_13C0_13C01D01__13C001  +  Vf_PRPP_GAR_GlyPool2_13C0_13C01D10__13C001  +  Vf_PRPP_GAR_GlyPool2_13C0_13C01D11__13C001  +  Vf_PRPP_GAR_GlyPool3_13C0_13C01D00__13C001  +  Vf_PRPP_GAR_GlyPool3_13C0_13C01D01__13C001  +  Vf_PRPP_GAR_GlyPool3_13C0_13C01D10__13C001  +  Vf_PRPP_GAR_GlyPool3_13C0_13C01D11__13C001;          %  cGAR_13C001 
		dxdt(443, 1) = -  Vf_GAR_FGAR_13C010_13C0D0__13C0100D0  -  Vf_GAR_FGAR_13C010_13C1D0__13C0101D0  -  Vf_GAR_FGAR_13C010_13C0D1__13C0100D1  -  Vf_GAR_FGAR_13C010_13C1D1__13C0101D1  +  Vf_PRPP_GAR_GlyPool1_13C0_13C10D00__13C010  +  Vf_PRPP_GAR_GlyPool1_13C0_13C10D01__13C010  +  Vf_PRPP_GAR_GlyPool1_13C0_13C10D10__13C010  +  Vf_PRPP_GAR_GlyPool1_13C0_13C10D11__13C010  +  Vf_PRPP_GAR_GlyPool2_13C0_13C10D00__13C010  +  Vf_PRPP_GAR_GlyPool2_13C0_13C10D01__13C010  +  Vf_PRPP_GAR_GlyPool2_13C0_13C10D10__13C010  +  Vf_PRPP_GAR_GlyPool2_13C0_13C10D11__13C010  +  Vf_PRPP_GAR_GlyPool3_13C0_13C10D00__13C010  +  Vf_PRPP_GAR_GlyPool3_13C0_13C10D01__13C010  +  Vf_PRPP_GAR_GlyPool3_13C0_13C10D10__13C010  +  Vf_PRPP_GAR_GlyPool3_13C0_13C10D11__13C010;          %  cGAR_13C010 
		dxdt(444, 1) = -  Vf_GAR_FGAR_13C011_13C0D0__13C0110D0  -  Vf_GAR_FGAR_13C011_13C1D0__13C0111D0  -  Vf_GAR_FGAR_13C011_13C0D1__13C0110D1  -  Vf_GAR_FGAR_13C011_13C1D1__13C0111D1  +  Vf_PRPP_GAR_GlyPool1_13C0_13C11D00__13C011  +  Vf_PRPP_GAR_GlyPool1_13C0_13C11D01__13C011  +  Vf_PRPP_GAR_GlyPool1_13C0_13C11D10__13C011  +  Vf_PRPP_GAR_GlyPool1_13C0_13C11D11__13C011  +  Vf_PRPP_GAR_GlyPool2_13C0_13C11D00__13C011  +  Vf_PRPP_GAR_GlyPool2_13C0_13C11D01__13C011  +  Vf_PRPP_GAR_GlyPool2_13C0_13C11D10__13C011  +  Vf_PRPP_GAR_GlyPool2_13C0_13C11D11__13C011  +  Vf_PRPP_GAR_GlyPool3_13C0_13C11D00__13C011  +  Vf_PRPP_GAR_GlyPool3_13C0_13C11D01__13C011  +  Vf_PRPP_GAR_GlyPool3_13C0_13C11D10__13C011  +  Vf_PRPP_GAR_GlyPool3_13C0_13C11D11__13C011;          %  cGAR_13C011 
		dxdt(445, 1) = -  Vf_GAR_FGAR_13C100_13C0D0__13C1000D0  -  Vf_GAR_FGAR_13C100_13C1D0__13C1001D0  -  Vf_GAR_FGAR_13C100_13C0D1__13C1000D1  -  Vf_GAR_FGAR_13C100_13C1D1__13C1001D1  +  Vf_PRPP_GAR_GlyPool1_13C1_13C00D00__13C100  +  Vf_PRPP_GAR_GlyPool1_13C1_13C00D01__13C100  +  Vf_PRPP_GAR_GlyPool1_13C1_13C00D10__13C100  +  Vf_PRPP_GAR_GlyPool1_13C1_13C00D11__13C100  +  Vf_PRPP_GAR_GlyPool2_13C1_13C00D00__13C100  +  Vf_PRPP_GAR_GlyPool2_13C1_13C00D01__13C100  +  Vf_PRPP_GAR_GlyPool2_13C1_13C00D10__13C100  +  Vf_PRPP_GAR_GlyPool2_13C1_13C00D11__13C100  +  Vf_PRPP_GAR_GlyPool3_13C1_13C00D00__13C100  +  Vf_PRPP_GAR_GlyPool3_13C1_13C00D01__13C100  +  Vf_PRPP_GAR_GlyPool3_13C1_13C00D10__13C100  +  Vf_PRPP_GAR_GlyPool3_13C1_13C00D11__13C100;          %  cGAR_13C100 
		dxdt(446, 1) = -  Vf_GAR_FGAR_13C101_13C0D0__13C1010D0  -  Vf_GAR_FGAR_13C101_13C1D0__13C1011D0  -  Vf_GAR_FGAR_13C101_13C0D1__13C1010D1  -  Vf_GAR_FGAR_13C101_13C1D1__13C1011D1  +  Vf_PRPP_GAR_GlyPool1_13C1_13C01D00__13C101  +  Vf_PRPP_GAR_GlyPool1_13C1_13C01D01__13C101  +  Vf_PRPP_GAR_GlyPool1_13C1_13C01D10__13C101  +  Vf_PRPP_GAR_GlyPool1_13C1_13C01D11__13C101  +  Vf_PRPP_GAR_GlyPool2_13C1_13C01D00__13C101  +  Vf_PRPP_GAR_GlyPool2_13C1_13C01D01__13C101  +  Vf_PRPP_GAR_GlyPool2_13C1_13C01D10__13C101  +  Vf_PRPP_GAR_GlyPool2_13C1_13C01D11__13C101  +  Vf_PRPP_GAR_GlyPool3_13C1_13C01D00__13C101  +  Vf_PRPP_GAR_GlyPool3_13C1_13C01D01__13C101  +  Vf_PRPP_GAR_GlyPool3_13C1_13C01D10__13C101  +  Vf_PRPP_GAR_GlyPool3_13C1_13C01D11__13C101;          %  cGAR_13C101 
		dxdt(447, 1) = -  Vf_GAR_FGAR_13C110_13C0D0__13C1100D0  -  Vf_GAR_FGAR_13C110_13C1D0__13C1101D0  -  Vf_GAR_FGAR_13C110_13C0D1__13C1100D1  -  Vf_GAR_FGAR_13C110_13C1D1__13C1101D1  +  Vf_PRPP_GAR_GlyPool1_13C1_13C10D00__13C110  +  Vf_PRPP_GAR_GlyPool1_13C1_13C10D01__13C110  +  Vf_PRPP_GAR_GlyPool1_13C1_13C10D10__13C110  +  Vf_PRPP_GAR_GlyPool1_13C1_13C10D11__13C110  +  Vf_PRPP_GAR_GlyPool2_13C1_13C10D00__13C110  +  Vf_PRPP_GAR_GlyPool2_13C1_13C10D01__13C110  +  Vf_PRPP_GAR_GlyPool2_13C1_13C10D10__13C110  +  Vf_PRPP_GAR_GlyPool2_13C1_13C10D11__13C110  +  Vf_PRPP_GAR_GlyPool3_13C1_13C10D00__13C110  +  Vf_PRPP_GAR_GlyPool3_13C1_13C10D01__13C110  +  Vf_PRPP_GAR_GlyPool3_13C1_13C10D10__13C110  +  Vf_PRPP_GAR_GlyPool3_13C1_13C10D11__13C110;          %  cGAR_13C110 
		dxdt(448, 1) = -  Vf_GAR_FGAR_13C111_13C0D0__13C1110D0  -  Vf_GAR_FGAR_13C111_13C1D0__13C1111D0  -  Vf_GAR_FGAR_13C111_13C0D1__13C1110D1  -  Vf_GAR_FGAR_13C111_13C1D1__13C1111D1  +  Vf_PRPP_GAR_GlyPool1_13C1_13C11D00__13C111  +  Vf_PRPP_GAR_GlyPool1_13C1_13C11D01__13C111  +  Vf_PRPP_GAR_GlyPool1_13C1_13C11D10__13C111  +  Vf_PRPP_GAR_GlyPool1_13C1_13C11D11__13C111  +  Vf_PRPP_GAR_GlyPool2_13C1_13C11D00__13C111  +  Vf_PRPP_GAR_GlyPool2_13C1_13C11D01__13C111  +  Vf_PRPP_GAR_GlyPool2_13C1_13C11D10__13C111  +  Vf_PRPP_GAR_GlyPool2_13C1_13C11D11__13C111  +  Vf_PRPP_GAR_GlyPool3_13C1_13C11D00__13C111  +  Vf_PRPP_GAR_GlyPool3_13C1_13C11D01__13C111  +  Vf_PRPP_GAR_GlyPool3_13C1_13C11D10__13C111  +  Vf_PRPP_GAR_GlyPool3_13C1_13C11D11__13C111;          %  cGAR_13C111 


    
    % cFGAR
    %  dxdt(1, 1) = (cFGAR) = 
    % - Vf_FGAR_AMP 
    % + Vf_GAR_FGAR;
    
		dxdt(449, 1) =  -  Vf_FGAR_AMP_13C0000D0_13C0D0__13C00000D00  -  Vf_FGAR_AMP_13C0000D0_13C1D0__13C00001D00  -  Vf_FGAR_AMP_13C0000D0_13C0D1__13C00000D01  -  Vf_FGAR_AMP_13C0000D0_13C1D1__13C00001D01  +  Vf_GAR_FGAR_13C000_13C0D0__13C0000D0;							%  cFGAR_13C0000D0 
		dxdt(450, 1) =  -  Vf_FGAR_AMP_13C0001D0_13C0D0__13C00010D00  -  Vf_FGAR_AMP_13C0001D0_13C1D0__13C00011D00  -  Vf_FGAR_AMP_13C0001D0_13C0D1__13C00010D01  -  Vf_FGAR_AMP_13C0001D0_13C1D1__13C00011D01  +  Vf_GAR_FGAR_13C000_13C1D0__13C0001D0;              %  cFGAR_13C0001D0 
		dxdt(451, 1) =  -  Vf_FGAR_AMP_13C0010D0_13C0D0__13C00100D00  -  Vf_FGAR_AMP_13C0010D0_13C1D0__13C00101D00  -  Vf_FGAR_AMP_13C0010D0_13C0D1__13C00100D01  -  Vf_FGAR_AMP_13C0010D0_13C1D1__13C00101D01  +  Vf_GAR_FGAR_13C001_13C0D0__13C0010D0;              %  cFGAR_13C0010D0 
		dxdt(452, 1) =  -  Vf_FGAR_AMP_13C0011D0_13C0D0__13C00110D00  -  Vf_FGAR_AMP_13C0011D0_13C1D0__13C00111D00  -  Vf_FGAR_AMP_13C0011D0_13C0D1__13C00110D01  -  Vf_FGAR_AMP_13C0011D0_13C1D1__13C00111D01  +  Vf_GAR_FGAR_13C001_13C1D0__13C0011D0;              %  cFGAR_13C0011D0 
		dxdt(453, 1) =  -  Vf_FGAR_AMP_13C0100D0_13C0D0__13C01000D00  -  Vf_FGAR_AMP_13C0100D0_13C1D0__13C01001D00  -  Vf_FGAR_AMP_13C0100D0_13C0D1__13C01000D01  -  Vf_FGAR_AMP_13C0100D0_13C1D1__13C01001D01  +  Vf_GAR_FGAR_13C010_13C0D0__13C0100D0;              %  cFGAR_13C0100D0 
		dxdt(454, 1) =  -  Vf_FGAR_AMP_13C0101D0_13C0D0__13C01010D00  -  Vf_FGAR_AMP_13C0101D0_13C1D0__13C01011D00  -  Vf_FGAR_AMP_13C0101D0_13C0D1__13C01010D01  -  Vf_FGAR_AMP_13C0101D0_13C1D1__13C01011D01  +  Vf_GAR_FGAR_13C010_13C1D0__13C0101D0;              %  cFGAR_13C0101D0 
		dxdt(455, 1) =  -  Vf_FGAR_AMP_13C0110D0_13C0D0__13C01100D00  -  Vf_FGAR_AMP_13C0110D0_13C1D0__13C01101D00  -  Vf_FGAR_AMP_13C0110D0_13C0D1__13C01100D01  -  Vf_FGAR_AMP_13C0110D0_13C1D1__13C01101D01  +  Vf_GAR_FGAR_13C011_13C0D0__13C0110D0;              %  cFGAR_13C0110D0 
		dxdt(456, 1) =  -  Vf_FGAR_AMP_13C0111D0_13C0D0__13C01110D00  -  Vf_FGAR_AMP_13C0111D0_13C1D0__13C01111D00  -  Vf_FGAR_AMP_13C0111D0_13C0D1__13C01110D01  -  Vf_FGAR_AMP_13C0111D0_13C1D1__13C01111D01  +  Vf_GAR_FGAR_13C011_13C1D0__13C0111D0;              %  cFGAR_13C0111D0 
		dxdt(457, 1) =  -  Vf_FGAR_AMP_13C1000D0_13C0D0__13C10000D00  -  Vf_FGAR_AMP_13C1000D0_13C1D0__13C10001D00  -  Vf_FGAR_AMP_13C1000D0_13C0D1__13C10000D01  -  Vf_FGAR_AMP_13C1000D0_13C1D1__13C10001D01  +  Vf_GAR_FGAR_13C100_13C0D0__13C1000D0;              %  cFGAR_13C1000D0 
		dxdt(458, 1) =  -  Vf_FGAR_AMP_13C1001D0_13C0D0__13C10010D00  -  Vf_FGAR_AMP_13C1001D0_13C1D0__13C10011D00  -  Vf_FGAR_AMP_13C1001D0_13C0D1__13C10010D01  -  Vf_FGAR_AMP_13C1001D0_13C1D1__13C10011D01  +  Vf_GAR_FGAR_13C100_13C1D0__13C1001D0;              %  cFGAR_13C1001D0 
		dxdt(459, 1) =  -  Vf_FGAR_AMP_13C1010D0_13C0D0__13C10100D00  -  Vf_FGAR_AMP_13C1010D0_13C1D0__13C10101D00  -  Vf_FGAR_AMP_13C1010D0_13C0D1__13C10100D01  -  Vf_FGAR_AMP_13C1010D0_13C1D1__13C10101D01  +  Vf_GAR_FGAR_13C101_13C0D0__13C1010D0;              %  cFGAR_13C1010D0 
		dxdt(460, 1) =  -  Vf_FGAR_AMP_13C1011D0_13C0D0__13C10110D00  -  Vf_FGAR_AMP_13C1011D0_13C1D0__13C10111D00  -  Vf_FGAR_AMP_13C1011D0_13C0D1__13C10110D01  -  Vf_FGAR_AMP_13C1011D0_13C1D1__13C10111D01  +  Vf_GAR_FGAR_13C101_13C1D0__13C1011D0;              %  cFGAR_13C1011D0 
		dxdt(461, 1) =  -  Vf_FGAR_AMP_13C1100D0_13C0D0__13C11000D00  -  Vf_FGAR_AMP_13C1100D0_13C1D0__13C11001D00  -  Vf_FGAR_AMP_13C1100D0_13C0D1__13C11000D01  -  Vf_FGAR_AMP_13C1100D0_13C1D1__13C11001D01  +  Vf_GAR_FGAR_13C110_13C0D0__13C1100D0;              %  cFGAR_13C1100D0 
		dxdt(462, 1) =  -  Vf_FGAR_AMP_13C1101D0_13C0D0__13C11010D00  -  Vf_FGAR_AMP_13C1101D0_13C1D0__13C11011D00  -  Vf_FGAR_AMP_13C1101D0_13C0D1__13C11010D01  -  Vf_FGAR_AMP_13C1101D0_13C1D1__13C11011D01  +  Vf_GAR_FGAR_13C110_13C1D0__13C1101D0;              %  cFGAR_13C1101D0 
		dxdt(463, 1) =  -  Vf_FGAR_AMP_13C1110D0_13C0D0__13C11100D00  -  Vf_FGAR_AMP_13C1110D0_13C1D0__13C11101D00  -  Vf_FGAR_AMP_13C1110D0_13C0D1__13C11100D01  -  Vf_FGAR_AMP_13C1110D0_13C1D1__13C11101D01  +  Vf_GAR_FGAR_13C111_13C0D0__13C1110D0;              %  cFGAR_13C1110D0 
		dxdt(464, 1) =  -  Vf_FGAR_AMP_13C1111D0_13C0D0__13C11110D00  -  Vf_FGAR_AMP_13C1111D0_13C1D0__13C11111D00  -  Vf_FGAR_AMP_13C1111D0_13C0D1__13C11110D01  -  Vf_FGAR_AMP_13C1111D0_13C1D1__13C11111D01  +  Vf_GAR_FGAR_13C111_13C1D0__13C1111D0;              %  cFGAR_13C1111D0 
		dxdt(465, 1) =  -  Vf_FGAR_AMP_13C0000D1_13C0D0__13C00000D10  -  Vf_FGAR_AMP_13C0000D1_13C1D0__13C00001D10  -  Vf_FGAR_AMP_13C0000D1_13C0D1__13C00000D11  -  Vf_FGAR_AMP_13C0000D1_13C1D1__13C00001D11  +  Vf_GAR_FGAR_13C000_13C0D1__13C0000D1;              %  cFGAR_13C0000D1 
		dxdt(466, 1) =  -  Vf_FGAR_AMP_13C0001D1_13C0D0__13C00010D10  -  Vf_FGAR_AMP_13C0001D1_13C1D0__13C00011D10  -  Vf_FGAR_AMP_13C0001D1_13C0D1__13C00010D11  -  Vf_FGAR_AMP_13C0001D1_13C1D1__13C00011D11  +  Vf_GAR_FGAR_13C000_13C1D1__13C0001D1;              %  cFGAR_13C0001D1 
		dxdt(467, 1) =  -  Vf_FGAR_AMP_13C0010D1_13C0D0__13C00100D10  -  Vf_FGAR_AMP_13C0010D1_13C1D0__13C00101D10  -  Vf_FGAR_AMP_13C0010D1_13C0D1__13C00100D11  -  Vf_FGAR_AMP_13C0010D1_13C1D1__13C00101D11  +  Vf_GAR_FGAR_13C001_13C0D1__13C0010D1;              %  cFGAR_13C0010D1 
		dxdt(468, 1) =  -  Vf_FGAR_AMP_13C0011D1_13C0D0__13C00110D10  -  Vf_FGAR_AMP_13C0011D1_13C1D0__13C00111D10  -  Vf_FGAR_AMP_13C0011D1_13C0D1__13C00110D11  -  Vf_FGAR_AMP_13C0011D1_13C1D1__13C00111D11  +  Vf_GAR_FGAR_13C001_13C1D1__13C0011D1;              %  cFGAR_13C0011D1 
		dxdt(469, 1) =  -  Vf_FGAR_AMP_13C0100D1_13C0D0__13C01000D10  -  Vf_FGAR_AMP_13C0100D1_13C1D0__13C01001D10  -  Vf_FGAR_AMP_13C0100D1_13C0D1__13C01000D11  -  Vf_FGAR_AMP_13C0100D1_13C1D1__13C01001D11  +  Vf_GAR_FGAR_13C010_13C0D1__13C0100D1;              %  cFGAR_13C0100D1 
		dxdt(470, 1) =  -  Vf_FGAR_AMP_13C0101D1_13C0D0__13C01010D10  -  Vf_FGAR_AMP_13C0101D1_13C1D0__13C01011D10  -  Vf_FGAR_AMP_13C0101D1_13C0D1__13C01010D11  -  Vf_FGAR_AMP_13C0101D1_13C1D1__13C01011D11  +  Vf_GAR_FGAR_13C010_13C1D1__13C0101D1;              %  cFGAR_13C0101D1 
		dxdt(471, 1) =  -  Vf_FGAR_AMP_13C0110D1_13C0D0__13C01100D10  -  Vf_FGAR_AMP_13C0110D1_13C1D0__13C01101D10  -  Vf_FGAR_AMP_13C0110D1_13C0D1__13C01100D11  -  Vf_FGAR_AMP_13C0110D1_13C1D1__13C01101D11  +  Vf_GAR_FGAR_13C011_13C0D1__13C0110D1;              %  cFGAR_13C0110D1 
		dxdt(472, 1) =  -  Vf_FGAR_AMP_13C0111D1_13C0D0__13C01110D10  -  Vf_FGAR_AMP_13C0111D1_13C1D0__13C01111D10  -  Vf_FGAR_AMP_13C0111D1_13C0D1__13C01110D11  -  Vf_FGAR_AMP_13C0111D1_13C1D1__13C01111D11  +  Vf_GAR_FGAR_13C011_13C1D1__13C0111D1;              %  cFGAR_13C0111D1 
		dxdt(473, 1) =  -  Vf_FGAR_AMP_13C1000D1_13C0D0__13C10000D10  -  Vf_FGAR_AMP_13C1000D1_13C1D0__13C10001D10  -  Vf_FGAR_AMP_13C1000D1_13C0D1__13C10000D11  -  Vf_FGAR_AMP_13C1000D1_13C1D1__13C10001D11  +  Vf_GAR_FGAR_13C100_13C0D1__13C1000D1;              %  cFGAR_13C1000D1 
		dxdt(474, 1) =  -  Vf_FGAR_AMP_13C1001D1_13C0D0__13C10010D10  -  Vf_FGAR_AMP_13C1001D1_13C1D0__13C10011D10  -  Vf_FGAR_AMP_13C1001D1_13C0D1__13C10010D11  -  Vf_FGAR_AMP_13C1001D1_13C1D1__13C10011D11  +  Vf_GAR_FGAR_13C100_13C1D1__13C1001D1;              %  cFGAR_13C1001D1 
		dxdt(475, 1) =  -  Vf_FGAR_AMP_13C1010D1_13C0D0__13C10100D10  -  Vf_FGAR_AMP_13C1010D1_13C1D0__13C10101D10  -  Vf_FGAR_AMP_13C1010D1_13C0D1__13C10100D11  -  Vf_FGAR_AMP_13C1010D1_13C1D1__13C10101D11  +  Vf_GAR_FGAR_13C101_13C0D1__13C1010D1;              %  cFGAR_13C1010D1 
		dxdt(476, 1) =  -  Vf_FGAR_AMP_13C1011D1_13C0D0__13C10110D10  -  Vf_FGAR_AMP_13C1011D1_13C1D0__13C10111D10  -  Vf_FGAR_AMP_13C1011D1_13C0D1__13C10110D11  -  Vf_FGAR_AMP_13C1011D1_13C1D1__13C10111D11  +  Vf_GAR_FGAR_13C101_13C1D1__13C1011D1;              %  cFGAR_13C1011D1 
		dxdt(477, 1) =  -  Vf_FGAR_AMP_13C1100D1_13C0D0__13C11000D10  -  Vf_FGAR_AMP_13C1100D1_13C1D0__13C11001D10  -  Vf_FGAR_AMP_13C1100D1_13C0D1__13C11000D11  -  Vf_FGAR_AMP_13C1100D1_13C1D1__13C11001D11  +  Vf_GAR_FGAR_13C110_13C0D1__13C1100D1;              %  cFGAR_13C1100D1 
		dxdt(478, 1) =  -  Vf_FGAR_AMP_13C1101D1_13C0D0__13C11010D10  -  Vf_FGAR_AMP_13C1101D1_13C1D0__13C11011D10  -  Vf_FGAR_AMP_13C1101D1_13C0D1__13C11010D11  -  Vf_FGAR_AMP_13C1101D1_13C1D1__13C11011D11  +  Vf_GAR_FGAR_13C110_13C1D1__13C1101D1;              %  cFGAR_13C1101D1 
		dxdt(479, 1) =  -  Vf_FGAR_AMP_13C1110D1_13C0D0__13C11100D10  -  Vf_FGAR_AMP_13C1110D1_13C1D0__13C11101D10  -  Vf_FGAR_AMP_13C1110D1_13C0D1__13C11100D11  -  Vf_FGAR_AMP_13C1110D1_13C1D1__13C11101D11  +  Vf_GAR_FGAR_13C111_13C0D1__13C1110D1;              %  cFGAR_13C1110D1 
		dxdt(480, 1) =  -  Vf_FGAR_AMP_13C1111D1_13C0D0__13C11110D10  -  Vf_FGAR_AMP_13C1111D1_13C1D0__13C11111D10  -  Vf_FGAR_AMP_13C1111D1_13C0D1__13C11110D11  -  Vf_FGAR_AMP_13C1111D1_13C1D1__13C11111D11  +  Vf_GAR_FGAR_13C111_13C1D1__13C1111D1;              %  cFGAR_13C1111D1 
		

    
    % cAMP
    %  dxdt(1, 1) = (cAMP) = 
    % - Vf_AMP_Degradation 
    % + Vf_FGAR_AMP;
    
		dxdt(481, 1)  =  - Vf_AMP_Degradation_13C00000D00_13C00000D00  +  Vf_FGAR_AMP_13C0000D0_13C0D0__13C00000D00;								%  cAMP_13C00000D00
		dxdt(482, 1)  =  - Vf_AMP_Degradation_13C00001D00_13C00001D00  +  Vf_FGAR_AMP_13C0000D0_13C1D0__13C00001D00;               %  cAMP_13C00001D00
		dxdt(483, 1)  =  - Vf_AMP_Degradation_13C00010D00_13C00010D00  +  Vf_FGAR_AMP_13C0001D0_13C0D0__13C00010D00;               %  cAMP_13C00010D00
		dxdt(484, 1)  =  - Vf_AMP_Degradation_13C00011D00_13C00011D00  +  Vf_FGAR_AMP_13C0001D0_13C1D0__13C00011D00;               %  cAMP_13C00011D00
		dxdt(485, 1)  =  - Vf_AMP_Degradation_13C00100D00_13C00100D00  +  Vf_FGAR_AMP_13C0010D0_13C0D0__13C00100D00;               %  cAMP_13C00100D00
		dxdt(486, 1)  =  - Vf_AMP_Degradation_13C00101D00_13C00101D00  +  Vf_FGAR_AMP_13C0010D0_13C1D0__13C00101D00;               %  cAMP_13C00101D00
		dxdt(487, 1)  =  - Vf_AMP_Degradation_13C00110D00_13C00110D00  +  Vf_FGAR_AMP_13C0011D0_13C0D0__13C00110D00;               %  cAMP_13C00110D00
		dxdt(488, 1)  =  - Vf_AMP_Degradation_13C00111D00_13C00111D00  +  Vf_FGAR_AMP_13C0011D0_13C1D0__13C00111D00;               %  cAMP_13C00111D00
		dxdt(489, 1)  =  - Vf_AMP_Degradation_13C01000D00_13C01000D00  +  Vf_FGAR_AMP_13C0100D0_13C0D0__13C01000D00;               %  cAMP_13C01000D00
		dxdt(490, 1)  =  - Vf_AMP_Degradation_13C01001D00_13C01001D00  +  Vf_FGAR_AMP_13C0100D0_13C1D0__13C01001D00;               %  cAMP_13C01001D00
		dxdt(491, 1)  =  - Vf_AMP_Degradation_13C01010D00_13C01010D00  +  Vf_FGAR_AMP_13C0101D0_13C0D0__13C01010D00;               %  cAMP_13C01010D00
		dxdt(492, 1)  =  - Vf_AMP_Degradation_13C01011D00_13C01011D00  +  Vf_FGAR_AMP_13C0101D0_13C1D0__13C01011D00;               %  cAMP_13C01011D00
		dxdt(493, 1)  =  - Vf_AMP_Degradation_13C01100D00_13C01100D00  +  Vf_FGAR_AMP_13C0110D0_13C0D0__13C01100D00;               %  cAMP_13C01100D00
		dxdt(494, 1)  =  - Vf_AMP_Degradation_13C01101D00_13C01101D00  +  Vf_FGAR_AMP_13C0110D0_13C1D0__13C01101D00;               %  cAMP_13C01101D00
		dxdt(495, 1)  =  - Vf_AMP_Degradation_13C01110D00_13C01110D00  +  Vf_FGAR_AMP_13C0111D0_13C0D0__13C01110D00;               %  cAMP_13C01110D00
		dxdt(496, 1)  =  - Vf_AMP_Degradation_13C01111D00_13C01111D00  +  Vf_FGAR_AMP_13C0111D0_13C1D0__13C01111D00;               %  cAMP_13C01111D00
		dxdt(497, 1)  =  - Vf_AMP_Degradation_13C10000D00_13C10000D00  +  Vf_FGAR_AMP_13C1000D0_13C0D0__13C10000D00;               %  cAMP_13C10000D00
		dxdt(498, 1)  =  - Vf_AMP_Degradation_13C10001D00_13C10001D00  +  Vf_FGAR_AMP_13C1000D0_13C1D0__13C10001D00;               %  cAMP_13C10001D00
		dxdt(499, 1)  =  - Vf_AMP_Degradation_13C10010D00_13C10010D00  +  Vf_FGAR_AMP_13C1001D0_13C0D0__13C10010D00;               %  cAMP_13C10010D00
		dxdt(500, 1)  =  - Vf_AMP_Degradation_13C10011D00_13C10011D00  +  Vf_FGAR_AMP_13C1001D0_13C1D0__13C10011D00;               %  cAMP_13C10011D00
		dxdt(501, 1)  =  - Vf_AMP_Degradation_13C10100D00_13C10100D00  +  Vf_FGAR_AMP_13C1010D0_13C0D0__13C10100D00;               %  cAMP_13C10100D00
		dxdt(502, 1)  =  - Vf_AMP_Degradation_13C10101D00_13C10101D00  +  Vf_FGAR_AMP_13C1010D0_13C1D0__13C10101D00;               %  cAMP_13C10101D00
		dxdt(503, 1)  =  - Vf_AMP_Degradation_13C10110D00_13C10110D00  +  Vf_FGAR_AMP_13C1011D0_13C0D0__13C10110D00;               %  cAMP_13C10110D00
		dxdt(504, 1)  =  - Vf_AMP_Degradation_13C10111D00_13C10111D00  +  Vf_FGAR_AMP_13C1011D0_13C1D0__13C10111D00;               %  cAMP_13C10111D00
		dxdt(505, 1)  =  - Vf_AMP_Degradation_13C11000D00_13C11000D00  +  Vf_FGAR_AMP_13C1100D0_13C0D0__13C11000D00;               %  cAMP_13C11000D00
		dxdt(506, 1)  =  - Vf_AMP_Degradation_13C11001D00_13C11001D00  +  Vf_FGAR_AMP_13C1100D0_13C1D0__13C11001D00;               %  cAMP_13C11001D00
		dxdt(507, 1)  =  - Vf_AMP_Degradation_13C11010D00_13C11010D00  +  Vf_FGAR_AMP_13C1101D0_13C0D0__13C11010D00;               %  cAMP_13C11010D00
		dxdt(508, 1)  =  - Vf_AMP_Degradation_13C11011D00_13C11011D00  +  Vf_FGAR_AMP_13C1101D0_13C1D0__13C11011D00;               %  cAMP_13C11011D00
		dxdt(509, 1)  =  - Vf_AMP_Degradation_13C11100D00_13C11100D00  +  Vf_FGAR_AMP_13C1110D0_13C0D0__13C11100D00;               %  cAMP_13C11100D00
		dxdt(510, 1)  =  - Vf_AMP_Degradation_13C11101D00_13C11101D00  +  Vf_FGAR_AMP_13C1110D0_13C1D0__13C11101D00;               %  cAMP_13C11101D00
		dxdt(511, 1)  =  - Vf_AMP_Degradation_13C11110D00_13C11110D00  +  Vf_FGAR_AMP_13C1111D0_13C0D0__13C11110D00;               %  cAMP_13C11110D00
		dxdt(512, 1)  =  - Vf_AMP_Degradation_13C11111D00_13C11111D00  +  Vf_FGAR_AMP_13C1111D0_13C1D0__13C11111D00;               %  cAMP_13C11111D00
		dxdt(513, 1)  =  - Vf_AMP_Degradation_13C00000D01_13C00000D01  +  Vf_FGAR_AMP_13C0000D0_13C0D1__13C00000D01;               %  cAMP_13C00000D01
		dxdt(514, 1)  =  - Vf_AMP_Degradation_13C00001D01_13C00001D01  +  Vf_FGAR_AMP_13C0000D0_13C1D1__13C00001D01;               %  cAMP_13C00001D01
		dxdt(515, 1)  =  - Vf_AMP_Degradation_13C00010D01_13C00010D01  +  Vf_FGAR_AMP_13C0001D0_13C0D1__13C00010D01;               %  cAMP_13C00010D01
		dxdt(516, 1)  =  - Vf_AMP_Degradation_13C00011D01_13C00011D01  +  Vf_FGAR_AMP_13C0001D0_13C1D1__13C00011D01;               %  cAMP_13C00011D01
		dxdt(517, 1)  =  - Vf_AMP_Degradation_13C00100D01_13C00100D01  +  Vf_FGAR_AMP_13C0010D0_13C0D1__13C00100D01;               %  cAMP_13C00100D01
		dxdt(518, 1)  =  - Vf_AMP_Degradation_13C00101D01_13C00101D01  +  Vf_FGAR_AMP_13C0010D0_13C1D1__13C00101D01;               %  cAMP_13C00101D01
		dxdt(519, 1)  =  - Vf_AMP_Degradation_13C00110D01_13C00110D01  +  Vf_FGAR_AMP_13C0011D0_13C0D1__13C00110D01;               %  cAMP_13C00110D01
		dxdt(520, 1)  =  - Vf_AMP_Degradation_13C00111D01_13C00111D01  +  Vf_FGAR_AMP_13C0011D0_13C1D1__13C00111D01;               %  cAMP_13C00111D01
		dxdt(521, 1)  =  - Vf_AMP_Degradation_13C01000D01_13C01000D01  +  Vf_FGAR_AMP_13C0100D0_13C0D1__13C01000D01;               %  cAMP_13C01000D01
		dxdt(522, 1)  =  - Vf_AMP_Degradation_13C01001D01_13C01001D01  +  Vf_FGAR_AMP_13C0100D0_13C1D1__13C01001D01;               %  cAMP_13C01001D01
		dxdt(523, 1)  =  - Vf_AMP_Degradation_13C01010D01_13C01010D01  +  Vf_FGAR_AMP_13C0101D0_13C0D1__13C01010D01;               %  cAMP_13C01010D01
		dxdt(524, 1)  =  - Vf_AMP_Degradation_13C01011D01_13C01011D01  +  Vf_FGAR_AMP_13C0101D0_13C1D1__13C01011D01;               %  cAMP_13C01011D01
		dxdt(525, 1)  =  - Vf_AMP_Degradation_13C01100D01_13C01100D01  +  Vf_FGAR_AMP_13C0110D0_13C0D1__13C01100D01;               %  cAMP_13C01100D01
		dxdt(526, 1)  =  - Vf_AMP_Degradation_13C01101D01_13C01101D01  +  Vf_FGAR_AMP_13C0110D0_13C1D1__13C01101D01;               %  cAMP_13C01101D01
		dxdt(527, 1)  =  - Vf_AMP_Degradation_13C01110D01_13C01110D01  +  Vf_FGAR_AMP_13C0111D0_13C0D1__13C01110D01;               %  cAMP_13C01110D01
		dxdt(528, 1)  =  - Vf_AMP_Degradation_13C01111D01_13C01111D01  +  Vf_FGAR_AMP_13C0111D0_13C1D1__13C01111D01;               %  cAMP_13C01111D01
		dxdt(529, 1)  =  - Vf_AMP_Degradation_13C10000D01_13C10000D01  +  Vf_FGAR_AMP_13C1000D0_13C0D1__13C10000D01;               %  cAMP_13C10000D01
		dxdt(530, 1)  =  - Vf_AMP_Degradation_13C10001D01_13C10001D01  +  Vf_FGAR_AMP_13C1000D0_13C1D1__13C10001D01;               %  cAMP_13C10001D01
		dxdt(531, 1)  =  - Vf_AMP_Degradation_13C10010D01_13C10010D01  +  Vf_FGAR_AMP_13C1001D0_13C0D1__13C10010D01;               %  cAMP_13C10010D01
		dxdt(532, 1)  =  - Vf_AMP_Degradation_13C10011D01_13C10011D01  +  Vf_FGAR_AMP_13C1001D0_13C1D1__13C10011D01;               %  cAMP_13C10011D01
		dxdt(533, 1)  =  - Vf_AMP_Degradation_13C10100D01_13C10100D01  +  Vf_FGAR_AMP_13C1010D0_13C0D1__13C10100D01;               %  cAMP_13C10100D01
		dxdt(534, 1)  =  - Vf_AMP_Degradation_13C10101D01_13C10101D01  +  Vf_FGAR_AMP_13C1010D0_13C1D1__13C10101D01;               %  cAMP_13C10101D01
		dxdt(535, 1)  =  - Vf_AMP_Degradation_13C10110D01_13C10110D01  +  Vf_FGAR_AMP_13C1011D0_13C0D1__13C10110D01;               %  cAMP_13C10110D01
		dxdt(536, 1)  =  - Vf_AMP_Degradation_13C10111D01_13C10111D01  +  Vf_FGAR_AMP_13C1011D0_13C1D1__13C10111D01;               %  cAMP_13C10111D01
		dxdt(537, 1)  =  - Vf_AMP_Degradation_13C11000D01_13C11000D01  +  Vf_FGAR_AMP_13C1100D0_13C0D1__13C11000D01;               %  cAMP_13C11000D01
		dxdt(538, 1)  =  - Vf_AMP_Degradation_13C11001D01_13C11001D01  +  Vf_FGAR_AMP_13C1100D0_13C1D1__13C11001D01;               %  cAMP_13C11001D01
		dxdt(539, 1)  =  - Vf_AMP_Degradation_13C11010D01_13C11010D01  +  Vf_FGAR_AMP_13C1101D0_13C0D1__13C11010D01;               %  cAMP_13C11010D01
		dxdt(540, 1)  =  - Vf_AMP_Degradation_13C11011D01_13C11011D01  +  Vf_FGAR_AMP_13C1101D0_13C1D1__13C11011D01;               %  cAMP_13C11011D01
		dxdt(541, 1)  =  - Vf_AMP_Degradation_13C11100D01_13C11100D01  +  Vf_FGAR_AMP_13C1110D0_13C0D1__13C11100D01;               %  cAMP_13C11100D01
		dxdt(542, 1)  =  - Vf_AMP_Degradation_13C11101D01_13C11101D01  +  Vf_FGAR_AMP_13C1110D0_13C1D1__13C11101D01;               %  cAMP_13C11101D01
		dxdt(543, 1)  =  - Vf_AMP_Degradation_13C11110D01_13C11110D01  +  Vf_FGAR_AMP_13C1111D0_13C0D1__13C11110D01;               %  cAMP_13C11110D01
		dxdt(544, 1)  =  - Vf_AMP_Degradation_13C11111D01_13C11111D01  +  Vf_FGAR_AMP_13C1111D0_13C1D1__13C11111D01;               %  cAMP_13C11111D01
		dxdt(545, 1)  =  - Vf_AMP_Degradation_13C00000D10_13C00000D10  +  Vf_FGAR_AMP_13C0000D1_13C0D0__13C00000D10;               %  cAMP_13C00000D10
		dxdt(546, 1)  =  - Vf_AMP_Degradation_13C00001D10_13C00001D10  +  Vf_FGAR_AMP_13C0000D1_13C1D0__13C00001D10;               %  cAMP_13C00001D10
		dxdt(547, 1)  =  - Vf_AMP_Degradation_13C00010D10_13C00010D10  +  Vf_FGAR_AMP_13C0001D1_13C0D0__13C00010D10;               %  cAMP_13C00010D10
		dxdt(548, 1)  =  - Vf_AMP_Degradation_13C00011D10_13C00011D10  +  Vf_FGAR_AMP_13C0001D1_13C1D0__13C00011D10;               %  cAMP_13C00011D10
		dxdt(549, 1)  =  - Vf_AMP_Degradation_13C00100D10_13C00100D10  +  Vf_FGAR_AMP_13C0010D1_13C0D0__13C00100D10;               %  cAMP_13C00100D10
		dxdt(550, 1)  =  - Vf_AMP_Degradation_13C00101D10_13C00101D10  +  Vf_FGAR_AMP_13C0010D1_13C1D0__13C00101D10;               %  cAMP_13C00101D10
		dxdt(551, 1)  =  - Vf_AMP_Degradation_13C00110D10_13C00110D10  +  Vf_FGAR_AMP_13C0011D1_13C0D0__13C00110D10;               %  cAMP_13C00110D10
		dxdt(552, 1)  =  - Vf_AMP_Degradation_13C00111D10_13C00111D10  +  Vf_FGAR_AMP_13C0011D1_13C1D0__13C00111D10;               %  cAMP_13C00111D10
		dxdt(553, 1)  =  - Vf_AMP_Degradation_13C01000D10_13C01000D10  +  Vf_FGAR_AMP_13C0100D1_13C0D0__13C01000D10;               %  cAMP_13C01000D10
		dxdt(554, 1)  =  - Vf_AMP_Degradation_13C01001D10_13C01001D10  +  Vf_FGAR_AMP_13C0100D1_13C1D0__13C01001D10;               %  cAMP_13C01001D10
		dxdt(555, 1)  =  - Vf_AMP_Degradation_13C01010D10_13C01010D10  +  Vf_FGAR_AMP_13C0101D1_13C0D0__13C01010D10;               %  cAMP_13C01010D10
		dxdt(556, 1)  =  - Vf_AMP_Degradation_13C01011D10_13C01011D10  +  Vf_FGAR_AMP_13C0101D1_13C1D0__13C01011D10;               %  cAMP_13C01011D10
		dxdt(557, 1)  =  - Vf_AMP_Degradation_13C01100D10_13C01100D10  +  Vf_FGAR_AMP_13C0110D1_13C0D0__13C01100D10;               %  cAMP_13C01100D10
		dxdt(558, 1)  =  - Vf_AMP_Degradation_13C01101D10_13C01101D10  +  Vf_FGAR_AMP_13C0110D1_13C1D0__13C01101D10;               %  cAMP_13C01101D10
		dxdt(559, 1)  =  - Vf_AMP_Degradation_13C01110D10_13C01110D10  +  Vf_FGAR_AMP_13C0111D1_13C0D0__13C01110D10;               %  cAMP_13C01110D10
		dxdt(560, 1)  =  - Vf_AMP_Degradation_13C01111D10_13C01111D10  +  Vf_FGAR_AMP_13C0111D1_13C1D0__13C01111D10;               %  cAMP_13C01111D10
		dxdt(561, 1)  =  - Vf_AMP_Degradation_13C10000D10_13C10000D10  +  Vf_FGAR_AMP_13C1000D1_13C0D0__13C10000D10;               %  cAMP_13C10000D10
		dxdt(562, 1)  =  - Vf_AMP_Degradation_13C10001D10_13C10001D10  +  Vf_FGAR_AMP_13C1000D1_13C1D0__13C10001D10;               %  cAMP_13C10001D10
		dxdt(563, 1)  =  - Vf_AMP_Degradation_13C10010D10_13C10010D10  +  Vf_FGAR_AMP_13C1001D1_13C0D0__13C10010D10;               %  cAMP_13C10010D10
		dxdt(564, 1)  =  - Vf_AMP_Degradation_13C10011D10_13C10011D10  +  Vf_FGAR_AMP_13C1001D1_13C1D0__13C10011D10;               %  cAMP_13C10011D10
		dxdt(565, 1)  =  - Vf_AMP_Degradation_13C10100D10_13C10100D10  +  Vf_FGAR_AMP_13C1010D1_13C0D0__13C10100D10;               %  cAMP_13C10100D10
		dxdt(566, 1)  =  - Vf_AMP_Degradation_13C10101D10_13C10101D10  +  Vf_FGAR_AMP_13C1010D1_13C1D0__13C10101D10;               %  cAMP_13C10101D10
		dxdt(567, 1)  =  - Vf_AMP_Degradation_13C10110D10_13C10110D10  +  Vf_FGAR_AMP_13C1011D1_13C0D0__13C10110D10;               %  cAMP_13C10110D10
		dxdt(568, 1)  =  - Vf_AMP_Degradation_13C10111D10_13C10111D10  +  Vf_FGAR_AMP_13C1011D1_13C1D0__13C10111D10;               %  cAMP_13C10111D10
		dxdt(569, 1)  =  - Vf_AMP_Degradation_13C11000D10_13C11000D10  +  Vf_FGAR_AMP_13C1100D1_13C0D0__13C11000D10;               %  cAMP_13C11000D10
		dxdt(570, 1)  =  - Vf_AMP_Degradation_13C11001D10_13C11001D10  +  Vf_FGAR_AMP_13C1100D1_13C1D0__13C11001D10;               %  cAMP_13C11001D10
		dxdt(571, 1)  =  - Vf_AMP_Degradation_13C11010D10_13C11010D10  +  Vf_FGAR_AMP_13C1101D1_13C0D0__13C11010D10;               %  cAMP_13C11010D10
		dxdt(572, 1)  =  - Vf_AMP_Degradation_13C11011D10_13C11011D10  +  Vf_FGAR_AMP_13C1101D1_13C1D0__13C11011D10;               %  cAMP_13C11011D10
		dxdt(573, 1)  =  - Vf_AMP_Degradation_13C11100D10_13C11100D10  +  Vf_FGAR_AMP_13C1110D1_13C0D0__13C11100D10;               %  cAMP_13C11100D10
		dxdt(574, 1)  =  - Vf_AMP_Degradation_13C11101D10_13C11101D10  +  Vf_FGAR_AMP_13C1110D1_13C1D0__13C11101D10;               %  cAMP_13C11101D10
		dxdt(575, 1)  =  - Vf_AMP_Degradation_13C11110D10_13C11110D10  +  Vf_FGAR_AMP_13C1111D1_13C0D0__13C11110D10;               %  cAMP_13C11110D10
		dxdt(576, 1)  =  - Vf_AMP_Degradation_13C11111D10_13C11111D10  +  Vf_FGAR_AMP_13C1111D1_13C1D0__13C11111D10;               %  cAMP_13C11111D10
		dxdt(577, 1)  =  - Vf_AMP_Degradation_13C00000D11_13C00000D11  +  Vf_FGAR_AMP_13C0000D1_13C0D1__13C00000D11;               %  cAMP_13C00000D11
		dxdt(578, 1)  =  - Vf_AMP_Degradation_13C00001D11_13C00001D11  +  Vf_FGAR_AMP_13C0000D1_13C1D1__13C00001D11;               %  cAMP_13C00001D11
		dxdt(579, 1)  =  - Vf_AMP_Degradation_13C00010D11_13C00010D11  +  Vf_FGAR_AMP_13C0001D1_13C0D1__13C00010D11;               %  cAMP_13C00010D11
		dxdt(580, 1)  =  - Vf_AMP_Degradation_13C00011D11_13C00011D11  +  Vf_FGAR_AMP_13C0001D1_13C1D1__13C00011D11;               %  cAMP_13C00011D11
		dxdt(581, 1)  =  - Vf_AMP_Degradation_13C00100D11_13C00100D11  +  Vf_FGAR_AMP_13C0010D1_13C0D1__13C00100D11;               %  cAMP_13C00100D11
		dxdt(582, 1)  =  - Vf_AMP_Degradation_13C00101D11_13C00101D11  +  Vf_FGAR_AMP_13C0010D1_13C1D1__13C00101D11;               %  cAMP_13C00101D11
		dxdt(583, 1)  =  - Vf_AMP_Degradation_13C00110D11_13C00110D11  +  Vf_FGAR_AMP_13C0011D1_13C0D1__13C00110D11;               %  cAMP_13C00110D11
		dxdt(584, 1)  =  - Vf_AMP_Degradation_13C00111D11_13C00111D11  +  Vf_FGAR_AMP_13C0011D1_13C1D1__13C00111D11;               %  cAMP_13C00111D11
		dxdt(585, 1)  =  - Vf_AMP_Degradation_13C01000D11_13C01000D11  +  Vf_FGAR_AMP_13C0100D1_13C0D1__13C01000D11;               %  cAMP_13C01000D11
		dxdt(586, 1)  =  - Vf_AMP_Degradation_13C01001D11_13C01001D11  +  Vf_FGAR_AMP_13C0100D1_13C1D1__13C01001D11;               %  cAMP_13C01001D11
		dxdt(587, 1)  =  - Vf_AMP_Degradation_13C01010D11_13C01010D11  +  Vf_FGAR_AMP_13C0101D1_13C0D1__13C01010D11;               %  cAMP_13C01010D11
		dxdt(588, 1)  =  - Vf_AMP_Degradation_13C01011D11_13C01011D11  +  Vf_FGAR_AMP_13C0101D1_13C1D1__13C01011D11;               %  cAMP_13C01011D11
		dxdt(589, 1)  =  - Vf_AMP_Degradation_13C01100D11_13C01100D11  +  Vf_FGAR_AMP_13C0110D1_13C0D1__13C01100D11;               %  cAMP_13C01100D11
		dxdt(590, 1)  =  - Vf_AMP_Degradation_13C01101D11_13C01101D11  +  Vf_FGAR_AMP_13C0110D1_13C1D1__13C01101D11;               %  cAMP_13C01101D11
		dxdt(591, 1)  =  - Vf_AMP_Degradation_13C01110D11_13C01110D11  +  Vf_FGAR_AMP_13C0111D1_13C0D1__13C01110D11;               %  cAMP_13C01110D11
		dxdt(592, 1)  =  - Vf_AMP_Degradation_13C01111D11_13C01111D11  +  Vf_FGAR_AMP_13C0111D1_13C1D1__13C01111D11;               %  cAMP_13C01111D11
		dxdt(593, 1)  =  - Vf_AMP_Degradation_13C10000D11_13C10000D11  +  Vf_FGAR_AMP_13C1000D1_13C0D1__13C10000D11;               %  cAMP_13C10000D11
		dxdt(594, 1)  =  - Vf_AMP_Degradation_13C10001D11_13C10001D11  +  Vf_FGAR_AMP_13C1000D1_13C1D1__13C10001D11;               %  cAMP_13C10001D11
		dxdt(595, 1)  =  - Vf_AMP_Degradation_13C10010D11_13C10010D11  +  Vf_FGAR_AMP_13C1001D1_13C0D1__13C10010D11;               %  cAMP_13C10010D11
		dxdt(596, 1)  =  - Vf_AMP_Degradation_13C10011D11_13C10011D11  +  Vf_FGAR_AMP_13C1001D1_13C1D1__13C10011D11;               %  cAMP_13C10011D11
		dxdt(597, 1)  =  - Vf_AMP_Degradation_13C10100D11_13C10100D11  +  Vf_FGAR_AMP_13C1010D1_13C0D1__13C10100D11;               %  cAMP_13C10100D11
		dxdt(598, 1)  =  - Vf_AMP_Degradation_13C10101D11_13C10101D11  +  Vf_FGAR_AMP_13C1010D1_13C1D1__13C10101D11;               %  cAMP_13C10101D11
		dxdt(599, 1)  =  - Vf_AMP_Degradation_13C10110D11_13C10110D11  +  Vf_FGAR_AMP_13C1011D1_13C0D1__13C10110D11;               %  cAMP_13C10110D11
		dxdt(600, 1)  =  - Vf_AMP_Degradation_13C10111D11_13C10111D11  +  Vf_FGAR_AMP_13C1011D1_13C1D1__13C10111D11;               %  cAMP_13C10111D11
		dxdt(601, 1)  =  - Vf_AMP_Degradation_13C11000D11_13C11000D11  +  Vf_FGAR_AMP_13C1100D1_13C0D1__13C11000D11;               %  cAMP_13C11000D11
		dxdt(602, 1)  =  - Vf_AMP_Degradation_13C11001D11_13C11001D11  +  Vf_FGAR_AMP_13C1100D1_13C1D1__13C11001D11;               %  cAMP_13C11001D11
		dxdt(603, 1)  =  - Vf_AMP_Degradation_13C11010D11_13C11010D11  +  Vf_FGAR_AMP_13C1101D1_13C0D1__13C11010D11;               %  cAMP_13C11010D11
		dxdt(604, 1)  =  - Vf_AMP_Degradation_13C11011D11_13C11011D11  +  Vf_FGAR_AMP_13C1101D1_13C1D1__13C11011D11;               %  cAMP_13C11011D11
		dxdt(605, 1)  =  - Vf_AMP_Degradation_13C11100D11_13C11100D11  +  Vf_FGAR_AMP_13C1110D1_13C0D1__13C11100D11;               %  cAMP_13C11100D11
		dxdt(606, 1)  =  - Vf_AMP_Degradation_13C11101D11_13C11101D11  +  Vf_FGAR_AMP_13C1110D1_13C1D1__13C11101D11;               %  cAMP_13C11101D11
		dxdt(607, 1)  =  - Vf_AMP_Degradation_13C11110D11_13C11110D11  +  Vf_FGAR_AMP_13C1111D1_13C0D1__13C11110D11;               %  cAMP_13C11110D11
		dxdt(608, 1)  =  - Vf_AMP_Degradation_13C11111D11_13C11111D11  +  Vf_FGAR_AMP_13C1111D1_13C1D1__13C11111D11;               %  cAMP_13C11111D11
           

end

