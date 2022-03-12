from numba import njit
import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter.ttk import Checkbutton
from tkinter import messagebox
import PyPDF3 as py3
import re
import os
import unidecode
import nltk
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
from numpy.core.numeric import full
import pubchempy as pcp
from pubchempy import Compound, get_compounds
import pandas as pd
from numbers import Real
from time import sleep
import sys


#executable windows
exe=tk.Tk()
exe.title('HPSTC')
exe.geometry('1280x1024')
exe.configure(background='#FFFFFF')

##icon
file_icone=os.path.dirname(__file__)
exe.iconbitmap(file_icone+"\\HPSTC.ico")

#logo
file_image_GDEM=os.path.dirname(__file__)
Logo_GDEM=PhotoImage(file=file_image_GDEM+"\\GDEM.png")
logo1=Label(exe,image=Logo_GDEM,background="#FFFFFF",foreground="#FFFFFF")
logo1.place(x=25,y=5)

file_image_UEPG=os.path.dirname(__file__)
Logo_UEPG=PhotoImage(file=file_image_GDEM+"\\UEPG.png")
logo2=Label(exe,image=Logo_UEPG,background="#FFFFFF",foreground="#FFFFFF")
logo2.place(x=1000,y=20)
                 
#display
display_box=Text(exe,foreground="#000000", font='arial,14',bd=2,relief='solid',cursor='arrow')
display_box.place(x=250,y=120,width=1000,height=550)
scrollbar_1=Scrollbar(display_box,orient=VERTICAL)
scrollbar_1.pack(side=RIGHT,fill=Y)
scrollbar_1.config(command=display_box.yview)
display_box.config(yscrollcommand=scrollbar_1.set)

# Box buttons
Select_file_box1=Button(exe,text="Select Article",background='#003333',foreground="#FFFFFF", cursor='hand2', font="arial",command= lambda:open_file())
Select_file_box1.place(x=20,y=150,width=200,height=60)


#Reset button
Select_file_box6=Button(exe,text="Reset",background='#003333',foreground="#FFFFFF", cursor='hand2', font="arial",command= lambda:reset_systen())
Select_file_box6.place(x=20,y=120,width=200,height=20)

def reset_systen():
    sys.stdout.flush()
    os.execl(sys.executable, sys.executable, *sys.argv[1:])
    

#open_file_path
root=tk.Tk()
root.withdraw()

chemical_elements=['','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']
chemical_elements_name=['', "hydrogen", "helium", "lithium", "beryllium", "boron", "carbon", "nitrogen", "oxygen", "fluorine", "neon", "sodium", "magnesium", "aluminium", "silicon", "phosphorus", "sulfur", "chlorine", "argon", "potassium", "calcium", "scandium", "titanium", "vanadium", "chromium", "manganese", "iron", "cobalt", "nickel", "copper", "zinc", "gallium", "germanium", "arsenic", "selenium",  "bromine", "krypton", "rubidium", "strontium", "yttrium", "zirconium", "niobium", "molybdenum", "technetium", "ruthenium", "rhodium", "palladium", "silver", "cadmium", "Indium", "tin", "antimony", "tellurium", "iodine", "xenon", "caesium", "barium",  "lanthanum", "cerium", "praseodymium", "neodymium", "Promethium", "samarium", "europium", "gadolinium", "terbium", "dysprosium", "holmium", "erbium", "thulium", "ytterbium", "lutetium", "hafnium", "tantalum", "tungsten", "rhenium", "osmium", "iridium",  "platinum", "gold", "mercury", "thallium", "lead", "bismuth", "polonium", "astatine", "radon", "francium", "radium", "actinium", "thorium", "protactinium", "uranium", "neptunium",  "plutonium", "americium",  "curium", "berkelium", "californium", "einsteinium", "fermium", "mendelevium", "nobelium", "lawrencium", "rutherfordium", "dubnium", "seaborgium", "bohrium", "hassium",  "mitnerium", "darmstadtium", "roentgenium", "copernicium", "nihonium", "flerovium", "moscovium", "livermorium", "tennessine", "oganesson"]

#-------------------------------------------------------Library Lr ------------------------------------------------------------
library_Lr=['Lrtext']

#-------------------------------------------------------Library Md ------------------------------------------------------------
library_Md=['Mdacid_ino','Mdbase_ino','Mdnoble_gases','Mdgases','Mdinert_gases','Mdwater','Mdperoxide','Mdsuperoxide','Mdalotropic_carbon', 'Mdmettalic_nanoparticle','Mdmaterials_inorganic','Mdacid_mine_ino','Mdsilane','Mdoxidant','Mdreducer','Mdbuffer','Mdelectrolyte']
acid_inorganic=['H2S','H2SO3','HNO2','HClO','HClO2','H3AsO4','H2AsO4','HClO3','H2CO3','H3BO3','HIO3','HN3','HBrO','HCN','HIO','HSCN','H2CrO4,','H3O+',  'H3PO3','HSO3F',  'H5IO6', 'H2SeO3', 'H3PO2', 'HPF6', 'H2SiF6', 'HPO3',  'H2PtCl6', 'H2PtCl6.6H2O', 'H3PMo12O40.12H2O', 'H4P2O7', 'H2MoO4', 'ClSO3H', 'HBF4', 'H2ZrF6', 'HAuCl4', 'H4SiO4', 'HBF4', 'H2ZrF6', 'HAuCl4', 'H4SiO4', 'H2WO4', 'HReO4', 'H2TiO3', 'HSbF6', 'regia']
acid_mine_ino=['HCl','HNO3','H2SO4','H3BO3','HF','HBr','HI','HClO4','H3PO4','hydrochloric','nitric','phosphoric','shulphuric','boric','hydrofluoric','hydrobromic','perchloric','hydroiodic']
base_inorganic=['OH-','LiOH','NaOH','KOH','RbOH','CsOH','FrOH','BeOH2','MgOH2','CaOH2','SrOH2','BaOH2','RaOH2','NH4OH','BOH3','Ag(OH)2','Fe(OH)2','Ni(OH)2','Al(OH)3','Fe(OH)3','Sn(OH)4','Pb(OH)4','Mn(OH)4','NH3','NaH','LiH','NaNH2','KNH2','Zn(OH)2','LiAlH4',' KAlH4','NaAlH4 ','Mg(AlH4)2 ', 'alkali']
noble_gases=['He','Ne','Ar','Kr','Xe','Rn','Og']
gases=['O2','He','Ne','Ar','Kr','Xe','Rn','H2','N2','CO','CO2','NH3','H2S','AsH3','Br2','F2','Cl2','ClO2','B2H6','HCN','HCl','HF','CH3SH','CH4','NO','NO2','O3','COCl2','PH3','SiH4','SO2','SO3','SF6','NF3','C2F6','CF4','CHClF2', 'N2O5']
inert_gases=['He','Ne','Ar','Kr','Xe','Rn','N2','CO2']
water=['water','H2O']
peroxide=['H2O2','Na2O2','K2O2','BaO2','ZnO2','H3PO5','H2SO5','peroxyacetic','PAA','K2S2O8']
superoxide=['KO2','NaO2','CsO2','RbO2','LiO2']
alotropic_carbon=['graphene','nanotubes','nanobuds','diamond','graphite','buckminsterfullerene','Fullerite','Fullerene','glassy']
mettalic_nanoparticle=['CuNps', 'AgNps', 'AuNps', 'PtNps', 'FeNps', 'FeONps', 'F2O3Nps', 'Fe3O4Nps', 'PdNps', 'NiNps', 'ZnONps', 'GdNps', 'TiO2Nps', 'RuNps', 'ReNps', 'CoNps']
materials_inorganic=['zeolite', 'nanoparticle', 'ceramic', 'ILs', 'boehmite','silica', 'alumina', 'titania', 'aluminosilicate', 'QD', 'QDs' 'organometallic', 'coordination', 'superconductor', 'quartz', 'steel', 'spine', 'semiconductor', 'fenton', 'radioactivity', 'conductor', 'alloy', 'intermetallic', 'bronze', 'colloid', 'aerosols', 'gels','MOF']
silane=['TEOS','MTEOS','MTMOS','MTMS','MPTMS','APTMOS']
agent_oxidant=['H2O2','NaClO','Na2S2O8','HClO4','KMnO4','KHSO5','oxone','(NH4)2S2O8','K2Cr2O7','NaIO4','Na2Cr2O7','H2SO4','O3','HNO3','H2S2O8','H2SO5','PCC','KNO3','NaBiO3','PbO2','PCC','NaH2BO04','(NH4)2Ce(NO3)6','Ce(SO4)2','LiClO4','I2','I2O5','MnO2','K2S2O8','KIO4','SeO2','NaClO3','NaBrO3','NaICl2']
agent_reducer=['NaBH4','NaH','LiAlH4','NaBH3CN','Na2S2O6','Na2S2O3','KI','N2H4','DIBAL-H','LiBH4','SmI2','CaH2','LiCl','Mg','K2BH4','Mg(BH4)2','KH']
buffer=['PBS','KH2PO4','K2HPO4','TAPS','Bicine','Tris','TAPSO','HEPES','TES','MOPS','PIPES','MEES','Cacodylate','citric','acetic','NaHPO4','Na2H2PO4','CHES','BARBITAL','Bis-Tris','MES','ADA','ACES','PIPES','MOPSO','BES','MOPS','TES','Trizma','Na2CO3','NaHCO3','Britton–Robinson', 'Sorensen', 'NH4Cl', 'HCN', 'KCN', 'Tris-HCl']
electrolyte=['NaCl','KCl','K2SO4','Na2SO4','LiCl','NaOH','H2SO4','tetrabutylammonium','tetraethylammonium','NaClO4','LiClO4']
reagent=acid_inorganic+acid_mine_ino+base_inorganic+noble_gases+gases+inert_gases+water+peroxide+superoxide+alotropic_carbon+mettalic_nanoparticle+materials_inorganic+silane+agent_oxidant+agent_reducer+buffer+electrolyte
#-------------------------------------------------------Library Bk ------------------------------------------------------------
library_Bk=['Bksolvent_organic','Bksolvent_aqueous']
solvent_organic=['acetone', 'methanol', 'ethanol', 'tetrachloroethylene', 'toluene', 'hexane', 'benzene', 'acetonitrile', '1-butanol','2-butanol', '2-butanone', 'chlorobenzene', 'chloroform', 'cyclohexane', 'dichloroethane', 'diglyme', 'DMSO', 'DMF', 'DME', 'glycerin', 'HMPA', 'heptane', 'HMPT', 'MTBE', 'NMP', 'nitromethane', 'pentane', 'ligroine', '1-propanol', '2-propanol', 'THF', 'm-xylene', 'o-xylene', 'p-xylene', 'acetic', 't-BuOH', 'CCl4', 'DEG', 'ether', 'glycol', 'CH2Cl2', 'ethyl-acetate', 'anisole', 'tetralin', 'pyridine', 'nitrobenzene', 'CS2', 'decane', 'octanol', 'aniline', 'benzonitrile', 'i-butanol', 'thinner', 'dioxane', 'NH3', 'n-Butanol', 'IPA', 'isooctane', 'Propylene-carbonate','CH3OH', 'CH3CH2OH', 'EtOH', 'MeOH']
solvent_aqueous=['water','H2O']

solvent=solvent_organic+solvent_aqueous
#-------------------------------------------------------Library Fm ------------------------------------------------------------
library_Fm=['Fmelectrode_modifield']
mettalic_nanoparticle=['CuNps', 'AgNps', 'AuNps', 'PtNps', 'FeNps', 'FeONps', 'F2O3Nps', 'Fe3O4Nps', 'PdNps', 'NiNps', 'ZnONps', 'GdNps', 'TiO2Nps', 'RuNps', 'ReNps', 'CoNps']
graphene=['GO','RGO']
nanotube=['MWCNT','SWCNT']
fullerene=['C60','C20','Buckyball','C70','Buckinsterfullerene','B80','C20','C72','C74','C76','C78','C80','C84','C86','C88','C90']
oxide=['Re2O7','V2O3','V2O5','ZnO','SnO2','Nb2O5','Al4O3','La2Al4O3','CeO2','LaGaO3','Ba2In2O5','ZrO2','Y2O3','TiO2','WO3','CuO','Cu2O','F2O3','UO2','Bi2O3','UO3','La2CuO4','VO2','MnO2','La2O3','Yb2O3','NiO','EuO','FeO']
titanate=['BaTiO3','CaTiO3','SrTiO3','PbTiO3']
niobate=['LiNbO3']
perovskite=titanate+niobate+['MgSiO3','LaMnO3','LaFeo3','CH3NH3PbI3','BiFeO3','LaAlO3','LaYbO3','MgFeSiO3','CaSiO3','LaMnO3','LaSrMno3','LSAT','PST','PbZrO3','MALHs','YAP','LuAP']
QDs=['CdS','CdSe','CdTe','ZnS','ZnS','ZnSe','GaAs','InAs','InP','PbS','PbSe','PbTe','InGaAs','ZnTe','HgS','HgSe','HgTe','BN','BP','BAs','BSb','AlN','AlP','AlAs','AlSb','GaN','GaP','GaSb','InN','InSb','TlN','TlP','TlAs','TlSb','GeS','GeSe','GeTe','SnS','SnSe','SnTe','InCuS','B12As2','Cu2S','SnS2','Tl2SnTe5','Tl2GeTe5','Bi2Te3','Cd3P2','Cd3As2','Cd3Sb2','Zn3P2','ZnP2','Zn2As3','Zn3Sb2','PbI2','MoS2','Bi2S','GaMnAs','InMnAs','CdMnTe','PbMnTe','CuInSe2','AgGaS2','ZnSiP2','As2S3','PtSi','BiI3', 'HgI2','TlBr','FeS','Cu2ZnSnS4','Cu2SnS3','EuS']
zeolite=['zeolite']
polymer=['PANi','PEDOT','PPy','PPS','polythiophene','PAC','PPV','PTh','Poly(dialkylfluorene)','Poly(thiophene)','Poly(pyrrole)','polyacetylene','SiPy','Nafion', 'polyphenylene','PPO','PPP']
precussor_Sc=['Sc(ClO4)3','ScCl3','Sc(SO3CF3)3','ScF3','Sc(OCH(CH3)2)3','ScCl3.6H2O','Sc(NO3)3.xH2O','[LSc(Me)Cl','[ScF6]','Sc(acac)3']
precussor_Y=['YF3','YBO3','(CF3SO3)3Y','YCl3','Y(CF3CO2)3','Y(acac)3','Y(CH3CO2)3','Y(hf-acac)3','Y(2-eha)3','Y(NO3)3','Y(naph)3']
precussor_La=['LaF3','LaB6','[La(acac)3(H2O)2]','La2(SO4)3.9H20','[La(EDTA)(H2O)4]','La(NO3)3.6H2O ','La(NO3)3','La2(SO4)3','LaBO3','LaCl3','LaI3','Ti(acac)2OiPr','TiOPc','TiPcCl2']
precussor_Ti=['[Et3NH]2[Ti(O2C6H4)3]','K4[TiO(O2C6H4)2].9H2O','TiBr4','TiCl4','TiCl3','Ti(OCH3)4','TiF4','TiO(acac)2','TEOT','TiCl4(THF)2','TDEAT','Cp2TiCl2','C10H15Cl3Ti']
precussor_Zr=['ZrCl4(THF)2','ZrCl4','TEMAZ','Zr(SO4)2','Zr[OC(CH3)3]4','ZrO(NO3)2','ZrF4','Bis(indenyl)zirconium','Zr(CH2CMe3)3','Zr(acac)4','[ZrF8]','']
precussor_Hf=['[HfCl4(THT)2]','HfCl4','HfCl2',' HfOCl2','HfI4','HfF4','TDMAH','TEMAH','Hf(acac)4','(Hf(tfac)4)','[HfF8]','Cyclopentadienylhafnium(IV)']  
precussor_Ce=['CeCl3','CeF4','CeCl3.7H2O','Ce(SO4)2','Ce2(SO4)3','Ce(NO3)3.6H2O','CeBr3','Ce2(SO4)3.8H2O','Ce(NH4)4(SO4)4.2H2O','Ce(ClO4)3','Ce(OTf)3']
precussor_Th=['ThCl4','Th(BH4)4','Th(NO3)4.4H2O']
precussor_V=['[VO(acac)2]','VCl3','OV(OC2H5)3','VCl2','VTIP','VOSO4','Bis(cyclopentadienyl)vanadium(II)','[VO(NCS)4]','VOPc']
precussor_Nb=['NbCl5','Nb(OCH2CH3)5','NbF5','((CH3)6N)3NbN(CH3)3','KNbO3','Nb(OEt)5','[NEt4][NbCl6]']
precussor_Ta=['Nb(OEt)5','TaCl2','TaCl5','Ta(OC2H5)5','TaF5','Ta(OCH3)5','Ta(OCH2CH2CH2CH3)5','Ta(N(CH3)2)5','Tris(diethylamino)(tert-butylimino)tantalum(V)']
precussor_Pr=['Pr(NO3)3.6H2O','PrCl3','Pr(acac)3','PrCl3.6H2O','Pr(NO3)3']
precussor_Pa=['Pa4O(O2)6F1','PaVO2(C2O4)']
precussor_Cr=['CrPic','Cr(CO)6','Cr(NO3)3.9H2O','CrCl2','CrCl3','Cr(acac)3','CrCl3.6H2O','Cr(NO3)3.9H2O','Cr(C5H5)2','CrF3.4H2O','Dibenzenechromium','Benchrotrene','[Cr(CN)6]','K2Cr2O7']
precussor_Mo=['MoCl5','MoOCl4','bis(cyclopentadienyl)molybdenum(IV)','[MoO2(acac)2]','Molybdenumhexacarbonyl','Na2MoO4.2H2O','(NH4)2MoO4','H2MoO4','[MoO2(C5H7O2)2]']
precussor_W=['W(CO)6','WCl4','WOCl4','WCl6','Bis(cyclopentadienyl)tungsten(IV)','Tetracarbonyl(1,5-cyclooctadiene)tungsten(0)','TpW(NO)Br2','[W(CN)8]']
precussor_Nd=['NdCl3','NdCl3','Nd(CF3SO3)3','NdF3','NdCl3.6H2O','Nd(NO3)3.6H2O','Nd(acac)3','Nd(Oi-Pr)3']
precussor_U=['UO2(NO3)2.6H2O','[UO2(NO3)2.6H2O]','[UI3(THF)4]','UCl4','UI3(MeCN)4','U(N(SiMe3)3']
precussor_Mn=['MnSO4 ','MnCl2','MnCl2','KMnO4','MnI2','Mn(acac)3','MnF2','Mn2(CO)10','Mn(NO3)2.4H2O',' Mn(acac)2','Mn(TMHD)3','Mn(OAc)2·4H2O','MMT','Decamethylmanganocene','Bis(isopropylcyclopentadienyl)manganesee','MnIIIClPc','MnIIClPc','MnPcCl','Bis(isopropylcyclopentadienyl)manganese','MnPc','MnPcCl']
precussor_Re=['ReCl3','Trioxo(triphenylsilyloxy)rhenium(VII)','ReCl5','[(C6H5)3P]2ReO2I','C20H21Cl3O2PReS','[(C6H5)3P]2ReOCl3','[(C6H5)3P]2Re(CH3CN)Cl3','[ReBr(CO)3(NCCH3)2]','Re2(CO)10','Re(CO)5Br','Re(CO)5Cl','Re(CO)5OSO2CF3','[Re(CO)3(dmso)3](CF3SO3)']
precussor_Fe=['FeCl3','FeCl2','FeBr2','Fe(acac)3','Fe(OTf)3','FePO4.2H2O','FeBr3','Decamethylferrocene','ferrocene','FeCl3.6H2O','FeSO4.7H2O','Fe(NO3)3.9H2O','Pentacarbonyliron(0)','Tetracarbonylbis(cyclopentadienyl)diiron','Tricarbonyl(cyclooctatetraene)iron(II)','Na3[FeF6]','Fe(dibm)3','Ferrocenium','Ferroceneboronic','Diironnonacarbonyl','Ferrocenecarboxaldehyde','FePc','FePcCl']
precussor_Ru=['Diiodo(p-cymene)ruthenium(II) ','Ru-PNN','RuI3','Ru(acac)3','RuCl3','Dihydridotetrakis(triphenylphosphine)ruthenium(II)','RuCl3.3H2O','[(C6H5)3P]3RuCl2','Bis(cyclopentadienyl)ruthenium(II)','Diethylruthenocene','RuCl(OAc)(PPh3)3','[RuCl2(mesitylene)]2','Dichloro(mesitylene)ruthenium(II)','[Ru(p-cymene)Cl2]2',]
precussor_Os=['OsCl3','OsCl4','K2OsCl6','Os3(CO)12','(NH4)2OsCl6','[OsBr2(PPh3)3]','OsH6(PiPr3)2 ']
precussor_Co=['(CoCl2⋅6H2O)','(Co(NO3)2⋅6H2O','Co(CH3COO)2·4H2O','CoPc','CoTsPc']
precussor_Rh=['RhCl3','Rh(acac)3','Rh(OAc)3','Rh(NO3)3','Rh2(SO4)3','Dicarbonyl(pentamethylcyclopentadienyl)rhodium(I)','Rh2(OOCCH3)4','[(CF3COO)2Rh]2','RhCl(PPh3)3','Rh(acac)(COD)','Rh2(TPA)4']   
precussor_Ir=['IrCl3','Ir(acac)3','Chlorobis(cyclooctene)iridium(I)dimer','[Ir(cod)(PCy3)(py)]PF6','Ir(acac)(COD)','Ir(Fppy)3','Ir(pq)2acac','F2Irpic','Ir(dFppy)3',' Ir(acac)(COD)','[Ir{dF(CF3)ppy}2(dtbpy)]PF6','(NH4)2IrCl6','Ir(acac)(CO2)','K2IrCl6','(PPZ)2Ir(acac)']
precussor_Sm=['SmCl3·6H2O','SmI3','Sm(SO3CF3)3','SmI2','Sm[OCH(CH3)2]3','Sm(NO3)3.6H2O','SmCl3','Sm(OAc)3','Sm(acac)2','Sm2(SO4)3.8H2O']
precussor_Ni=['NiI2','NiBr2','NiCl2','Ni(OAc)2.4H2O','Ni(NO3)2.6H2O','[((CH3)3P)]2NICl2','[(C6H5)3P]4Ni','Ni(acac)2','NiSO4.6H2O','Nickelocene','NiBr2.3H2O','Dibromobis(triphenylphosphine)nickel(II)','NiPc','NiTsPc']
precussor_Pd=['PdBr2','Pd(NO3)2','PdSO4','PdI2','Pd(acac)2','Pd(CN)2','Pd(NO3)2.2H2O','PdCl2','Pd(TFA)2','Pd(PPh3)4','[Pd(1-phenylallyl)Cl]2','Pd(OAc)2','PdCl2(PPh3)2','Pd(OAc)2(PPh3)2','Pd(dba)2','PdCl2[P(cy)3]2','PdCl2(cod)','Pd(OAc)2','Pd(dba)2','CX22','CX31','CX32','Pd[(o-tol)3P]2']
precussor_Pt=['PtI2','PtBr2','cis-Dichlorobis(pyridine)platinum(II)','Pt(acac)2','Dichloro(ethylenediamine)platinum(II)','Platinum(0)-tetrakis(triphenylphosphine)','PtCl2','Trimethyl(methylcyclopentadienyl)platinum(IV)','Pd(NH3)4(NO3)2','Pt(DMSO)2Cl2','Dichloro(1,10-phenanthroline)platinum(II)','HS161','Pt(NH3)2(NO2)2','H2PtCl6.6H2O','HS425','HS432','K2PtCl6','H2PtCl6','Pt(COD)Cl2','Pt(NH2)Cl2','PtOEP']
precussor_Eu=['EuI2','EuCl3','Eu(OAc)3','Eu(NO3)3.5H2O','EuCl2', 'Eu(hfc)3','Eu(NO3)3','EuCl3.6H2O','Eu(acac)3', 'Eu(hfc)3','Eu(facam)3','Eu(tfc)3','Eu(dbm)3(phen)']    
precussor_Cu=['CuCl2','CuCl','CuCO2CH3','Bromotris(triphenylphosphine)copper(I)','Cu(acac)2','Phenylthiocopper(I)','CuI','CuBr2','CuSO4','CuSO4.5H2O','CuCl2.2H2O','Cu(OTs)2','(Ethylcyclopentadienyl)(triphenylphosphine)copper(I)','[(CH3CN)4Cu]PF6','Bis(ethylenediamine)copper(II)','[Cu(bpy)2]','CuPc','CuTsPc','Cu-TMEDA','(Ethylcyclopentadienyl)(triphenylphosphine)copper(I)',]
precussor_Ag=['AgNO3','Ag(NH2)2','Ag2(SO4)','Ag3PO4','AgNO2','AgBF4','AgCN','AgPF6','AgF2','AgSCN','Ag2CO3','AgOCN','AgClO4','Ag2CrO4','AgI','AgClO3','AgReO4','K[Ag(CN)2]']
precussor_Au=['AuCl','Chloro(triethylphosphine)gold(I)','AuI','AuCl3','Chloro(trimethylphosphine)gold(I)','Chloro(dimethylsulfide)gold(I)','AuBr3','Chloro(methyldiphenylphosphine)gold(I)','Methyl(triphenylphosphine)gold(I)','Chloro(triphenylphosphine)gold(I)','Dichloro(2-pyridinecarboxylato)gold','[Tris(para-trifluoromethylphenyl)phosphine]gold(I)','(tBu3P)AuCl','NaAuCl4.2H2O','HAuCl4','AuF','AuF3','AuF5','[Au(SIPr)-(alkoxide)]','[AuCl(PPh3)]','CAuClO','Dichloro[(+-)−BINAP]digold(I)']
precussor_Zn=['ZnCl2',' Ziram','ZrBr2','ZnI2','ZnF2','Zn3(PO4)2','ZnCl2-TMEDA','ZnCl2-TMEDA','Zn(NTf2)2','Dichlorobis(tetrahydrofuran)zinc','ZnSO4.7H2O','Zn(OAc)2','ZnPc','Zn(NTf2)2']
precussor_Er=['ErCl3']
precussor_Gd=['GdCl3']
mettalic_complex_precussor=precussor_Sc+precussor_Y+precussor_La+precussor_Ti+precussor_Zr+precussor_Hf+precussor_Ce+precussor_Th+precussor_V+precussor_Nb+precussor_Ta+precussor_Pa+precussor_Cr+precussor_Mo+precussor_W+precussor_U+precussor_Mn+precussor_Re+precussor_Fe+precussor_Ru+precussor_Os+precussor_Co+precussor_Rh+precussor_Ir+precussor_Sm+precussor_Ni+precussor_Pd+precussor_Pt+precussor_Eu+precussor_Cu+precussor_Ag+precussor_Au+precussor_Zn+precussor_Er+precussor_Gd
ILs=['ILs']
Enzyme=['Enzyme']
MOF=['MOF']
physical=['irradiation','colloidal']
mettalic=['alloy','powder','Ti','Fe','Ru','Ni','Co','Pt','Cu','Au','Pd','cluster']
electrode_modifield=mettalic_nanoparticle+graphene+nanotube+fullerene+oxide+perovskite+QDs+zeolite+polymer+mettalic_complex_precussor+ILs+Enzyme
#-------------------------------------------------------Library Am ------------------------------------------------------------
library_Am=['Amtoken_chemical','Amanalyt_all', 'Amanalyt_phenol_enol','Amanalyt_hetcyclicaro','Amanalyt_thiophenol/tiol','Amanalyt_amine','Amanalyt_thiocarbamate','Amanalyt_sulfonamide','Amanalyt_imine','Amanalyt_aldehyde']
het_cycle=['thiazin', 'thiadiazol', 'imidazol', 'indol', 'acridin', 'purin', 'benzo', 'furan', 'benzofuran', 'pyrido', 'triazine', 'quinolin', 'pyrimidin', 'pteridin', 'imidazol', 'pyrazol', 'pyridine', 'pyrazol', 'piperidin', 'pyrimidin', 'purine', 'benzotriazin', 'triazol', 'thiadiazol', 'thiophene', 'benzotriazol', 'pyrrole', 'thiophen', 'aniline', 'carbazole', 'quinoxal', 'benzoate', 'acridine', 'phenarsazinine', 'tetrazole', 'quinolin', 'triazin', 'thioxanthen', 'cinnoline', 'tetrazine', 'thiazol', 'pyrazol', 'pyran', 'thiopyran', 'pterid', 'quino', 'thianthrene', 'oxazol', 'porphyrin', 'chromen', 'furo', 'oxazol', 'triazine', 'thioxanthene', 'quino', 'thione', 'thianthrene', 'thiazole', 'thiadiazole', 'dioxin', 'pyrazine', 'phenanthroline', 'oxadiazole']

#-------------------------------------------------------Library Cf ------------------------------------------------------------
library_Cf=['Cfanalyt_solubility','Cfarticle_name']
#-------------------------------------------------------Library No ------------------------------------------------------------
library_No=['Noscan_rate','Noconcentration']


#open file pdf

def open_file():
    file_pdf=filedialog.askopenfiles(mode='rb',title="Select file",filetypes=(("Text files","*.pdf"),("All files","*.*")))

#----------------------------------------------------------Attributes----------------------------------------------------------------------------------------------
    #button file attribute
    Select_file_box2=Button(exe,text="Select Attributes ",background='#003333',foreground="#FFFFFF", cursor='hand2', font="arial",command= lambda:attribute_file())
    Select_file_box2.place(x=20,y=230,width=200,height=60)

    def attribute_file():
        attributes_open = filedialog.askopenfile(mode='rb',title="Select File Attributes",filetypes=(("Text files","*.pdf"),("All files","*.*")))
      
        attributes_read =py3.PdfFileReader(attributes_open,strict=False)
        attributes_name=os.path.basename(attributes_open.name)
        display_box.insert(tk.END,'\nSuccessfully imported the  ' + attributes_name +'\n')
         
        page_attributes=attributes_read.getPage(0)                                                                                                                                                                                                                                                                                                                       
        attributes_extract=page_attributes.extractText()
        attributes_extract=re.sub('\n','',attributes_extract) 
  
        word_token_attribute=nltk.tokenize.word_tokenize(str(attributes_extract))
        list(word_token_attribute)
      

        Save_button=Button(exe,text="Save Extract",background='#003333',foreground="#FFFFFF", cursor='hand2', font="arial",command= lambda:save_path_extract())
        Save_button.place(x=20,y=320,width=200,height=60)
                       
        def save_path_extract():
            save_arff=filedialog.asksaveasfile(title="Save extract",defaultextension='.arff',filetypes=[("Weka File","*.*")])
            o=len(word_token_attribute)
            save_arff.writelines('%'+'Text fractions obtained by HPSTC'+"\n"+'\n')
            save_arff.writelines('@relation '+'HPSTC'+"\n"+'\n')
        
            for p in range(0,o): 
                
                if word_token_attribute[p] in library_Lr[0]:
                    result_attribute=' string'
                    save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                tokenchemical_ini=0
                if word_token_attribute[p] in library_Am:
                    if word_token_attribute[p] in library_Am[0]:
                        tokenchemical_ini=1
                        
                    if word_token_attribute[p] in library_Am[1]:
                        result_attribute=' { phenol/enol, het_cyclic_aro, thiophenol/thiol, amine, thiocarbamate, sulfonamide, imine, aldehyde, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n") 
                        analyt_soluble=0
                        space_analyt=1
                    if word_token_attribute[p] in library_Am[2]:
                        result_attribute=' { phenol/enol, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n") 
                        space_analyt=0
                    if word_token_attribute[p] in library_Am[3]:
                        result_attribute=' { het_cyclic_aro, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")     
                        space_analyt=0
                    if word_token_attribute[p] in library_Am[4]:
                        result_attribute=' { thiophenol/thiol, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n") 
                        space_analyt=0
                    if word_token_attribute[p] in library_Am[5]:
                        result_attribute=' { amine, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        space_analyt=0
                    if word_token_attribute[p] in library_Am[6]:
                        result_attribute=' { thiocarbamate, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")     
                        space_analyt=0
                    if word_token_attribute[p] in library_Am[7]:
                        result_attribute=' { sulfonamide, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n") 
                        space_analyt=0

                    if word_token_attribute[p] in library_Am[8]:
                        result_attribute=' { imine, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        space_analyt=0
                    if word_token_attribute[p] in library_Am[9]:
                        result_attribute=' { aldehyde, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        space_analyt=0   
                if word_token_attribute[p] in library_Cf:
                    name_article=0                 
                    if word_token_attribute[p] in library_Cf[0]:
                        analyt_soluble=1
                        result_attribute=' { Water , THF, No}'
                        save_arff.writelines('@attribute '+'Cf1'+word_token_attribute[p]+result_attribute+"\n")  
                        result_attribute=' { KCl, LiClO4, No}'
                        save_arff.writelines('@attribute '+'Cf2'+word_token_attribute[p]+result_attribute+"\n")                                   
                    if word_token_attribute[p] in library_Cf[1]: 
                        result_attribute=' string';name_article=1
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                research_number=0
                if word_token_attribute[p] in library_No:
                    if word_token_attribute[p] in library_No[0]:
                        result_attribute=' Real';research_number=1   
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n") 
                    
                    if word_token_attribute[p] in library_No[1]:
                        result_attribute=' Real';research_number=1
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n") 
                       
                if word_token_attribute[p] in library_Fm[0]:
                    result_attribute=' { mettalic_nanoparticle, graphene, oxide, nanotubes, fullerene, ILs, QDs,perovskite, zeolite, polymer, physical, enzyme, MOF, mettalic,mettalic_complex_precussor, No_or_other}'
                    save_arff.writelines('@attribute '+'At1'+word_token_attribute[p]+result_attribute+"\n")
                    save_arff.writelines('@attribute '+'As1'+word_token_attribute[p]+' string '+"\n")
                    save_arff.writelines('@attribute '+'At2'+word_token_attribute[p]+result_attribute+"\n")
                    save_arff.writelines('@attribute '+'As2'+word_token_attribute[p]+' string '+"\n")
                    save_arff.writelines('@attribute '+'At3'+word_token_attribute[p]+result_attribute+"\n")
                    save_arff.writelines('@attribute '+'As3'+word_token_attribute[p]+' string '+"\n")
          
                if word_token_attribute[p] in library_Bk:
                    if word_token_attribute[p] in library_Bk[0]:
                        result_attribute=' { acetone, methanol, ethanol, tetrachloroethylene, toluene, hexane, benzene, acetonitrile, 1-butanol,2-butanol, 2-butanone, chlorobenzene, chloroform, cyclohexane, dichloroethane, diglyme, DMSO, DMF, DME, glycerin, HMPA, heptane, HMPT, MTBE, NMP, nitromethane, pentane, ligroine, 1-propanol, 2-propanol, THF, m-xylene, o-xylene, p-xylene, acetic, t-BuOH, CCl4, DEG, ether, glycol, CH2Cl2, ethyl-acetate, anisole, tetralin, pyridine, nitrobenzene, CS2, decane, octanol, aniline, benzonitrile, i-butanol, thinner, dioxane, NH3, n-Butanol, IPA, isooctane, Propylene-carbonate, CH3OH, CH3CH2OH, EtOH, MeOH, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                
                    if word_token_attribute[p] in library_Bk[1]:
                        result_attribute=' { water, H2O, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        
                if word_token_attribute[p] in library_Md:
                    if word_token_attribute[p] in library_Md[0]:
                        result_attribute=' { HF, HBr, HI, H2S, H2SO3, HNO2, HClO, HClO2, H3AsO4, H2AsO4, HClO3, HClO4, H2CO3, H3BO3, HIO3, HN3, HBrO, HCN, HIO, HSCN, H2CrO4, H3O+, H3PO3, HSO3F,  H5IO6, H2SeO3, H3PO2, HPF6, H2SiF6, HPO3,  H2PtCl6, H2PtCl6.6H2O, H3PMo12O40.12H2O, H4P2O7, H2MoO4, ClSO3H, HBF4, H2ZrF6, HAuCl4, H4SiO4, H2WO4, HReO4, H2TiO3, HSbF6, regia, acid No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                    
                    if word_token_attribute[p] in library_Md[1]:
                        result_attribute=' { OH-, LiOH, NaOH, KOH, RbOH, CsOH, FrOH, Be(OH)2, Mg(OH)2, Ca(OH)2, Sr(OH2), Ba(OH2), Ra(OH2), NH4OH, B(OH)3, Ag(OH)2, Fe(OH)2, Ni(OH)2, Al(OH)3, Fe(OH3), Sn(OH)4, Pb(OH)4, Mn(OH)4, NH3, NaH, LiH, NaNH2, KNH2, Zn(OH)2, LiAlH4, KAlH4, NaAlH4, Mg(AlH4)2, alkali, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        
                    if word_token_attribute[p] in library_Md[2]:
                        result_attribute='{ He, Ne,Ar, Kr, Xe, Rn, Og, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        
                    if word_token_attribute[p] in library_Md[3]:
                        result_attribute=' { O2, He, Ne, Ar, Kr, Xe, Rn, H2, N2, CO, CO2, NH3, H2S, AsH3, Br2, F2, Cl2, ClO2, B2H6, HCN, HCl, HF, CH3SH, CH4, NO, NO2, O3, COCl2, PH3, SiH4, SO2, SO3, SF6, NF3, C2F6, CF4, CHClF2, N2O5, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                    
                    if word_token_attribute[p] in library_Md[4]:
                        result_attribute=' { He, Ne, Ar, Kr, Xe, Rn, N2, CO2, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                    
                    if word_token_attribute[p] in library_Md[5]:
                        result_attribute=' { water, H2O, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        
                    if word_token_attribute[p] in library_Md[6]:
                        result_attribute=' { H2O2, Na2O2, K2O2, BaO2, ZnO2, H3PO5, H2SO5, peroxyacetic, PAA, K2S2O8, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                        
                    if word_token_attribute[p] in library_Md[7]:
                        result_attribute=' { KO2, NaO2, CsO2, RbO2, LiO2, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                    
                    if word_token_attribute[p] in library_Md[8]:
                        result_attribute=' { graphene, nanotubes, nanobuds, diamond, graphite, buckminsterfullerene, Fullerite, Fullerene, glassy, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                    
                    if word_token_attribute[p] in library_Md[9]:
                        result_attribute=' { CuNps, AgNps, AuNps, PtNps, FeNps, FeONps, F2O3Nps, Fe3O4, PdNps, NiNps, ZnONps, GdNps, TiO2Nps, RuNps, ReNps, CoNps, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                   
                    if word_token_attribute[p] in library_Md[10]:
                        result_attribute=' { zeolite,nanoparticle, ceramic, ILs, boehmite, silica, alumina, titania,  aluminosilicate, QD, QDs organometallic, coordination,  superconductor, quartz, steel, spinel, semiconductor, fenton, radioactivity, conductor, alloy, intermetallic, bronze, colloid, aerosols, gels, MOF, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                     
                    if word_token_attribute[p] in library_Md[11]:
                        result_attribute=' { HCl, HNO3, H2SO4, H3BO3, HF, HBr, HI, HClO4, H3PO4, hydrochloric, nitric, phosphoric, shulphuric, boric, hydrofluoric, hydrobromic, perchloric, hydroiodic, No or other}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")   
                    
                    if word_token_attribute[p] in library_Md[12]:
                        result_attribute=' { TEOS, MTEOS, MTMOS, MTMS, MPTMS, APTMOS, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                    
                    if word_token_attribute[p] in library_Md[13]:
                        result_attribute=' { H2O2, NaClO, Na2S2O8, HClO4, KMnO4, KHSO5, oxone, (NH4)2S2O8, K2Cr2O7, NaIO4, Na2Cr2O7, H2SO4, O3, HNO3, H2S2O8, H2SO5, PCC, KNO3, NaBiO3, PbO2, PCC, NaH2BO04, (NH4)2Ce(NO3)6, Ce(SO4)2, LiClO4, I2, I2O5, MnO2, K2S2O8, KIO4, SeO2, NaClO3, NaBrO3, NaICl2, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                                             
                    if word_token_attribute[p] in library_Md[14]:
                        result_attribute=' { NaBH4, NaH, LiAlH4, NaBH3CN, Na2S2O6, Na2S2O3, KI, N2H4, DIBAL-H, LiBH4, SmI2, CaH2, LiCl , Mg, K2BH4, Mg(BH4)2, KH, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")                         
                    
                    if word_token_attribute[p] in library_Md[15]:
                        result_attribute=' { PBS, KH2PO4, K2HPO4, TAPS, Bicine, Tris, TAPSO, HEPES, TES, MOPS, PIPES, MEES, Cacodylate, citric, acetic, NaHPO4, Na2H2PO4,C HES, BARBITAL, Bis-Tris, MES, ADA, ACES, PIPES, MOPSO, BES, MOPS, TES, Trizma, Na2CO3, NaHCO3, BrittonRobinson, Sorensen, NH4Cl, HCN, KCN, Tris-HCl, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                   
                    if word_token_attribute[p] in library_Md[16]:
                        result_attribute=' { NaCl, KCl, K2SO4, Na2SO4, LiCl, NaOH, H2SO4, HCl , tetrabutylammonium, tetraethylammonium, NaClO4, LiClO4, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")
                
                else:
                    library_all=library_Bk+library_Md+library_Lr+library_Fm+library_Am+library_Cf+library_No
                    if word_token_attribute[p] not in library_all:
                        result_attribute=' { Yes, No}'
                        save_arff.writelines('@attribute '+word_token_attribute[p]+result_attribute+"\n")    
                        tokenchemical_ini=0
                   
            save_arff.writelines("\n\n"+'@data'+"\n\n")
#------------------------------------------------------Start-------------------------------------------------------------------------------------------------------------
    
            Start_button=Button(exe,text="Start",background='#003333',foreground="#FFFFFF", cursor='hand2', font="arial",command= lambda:start_read())
            Start_button.place(x=20,y=390,width=200,height=60)
            
            def  start_read():   
                word_token_all=[]
                list(word_token_all)
           
                #number pdf files 
                k=len(file_pdf) 
                
                for x in range(0,k):
                    pdf_read=py3.PdfFileReader(file_pdf[x],strict=False)
                    pdf_n=os.path.basename(file_pdf[x].name)
                    display_box.insert(tk.END, '\nAnalyzing '+ str(x+1)+ ' article of '+ str(k))
                    display_box.insert(tk.END,'\nInitializing reading of: ' + pdf_n)
                    
                    #number pages
                    n=pdf_read.numPages 
                      
                    for i in range(0,n):
                        page=pdf_read.getPage(i)                                                                                                                                                                                                                                                                                                                       
                        text_page=page.extractText()
                        text_page=re.sub('\n','',text_page)
                        
                        #extract name pdf file
                        text_page=unidecode.unidecode(text_page)
                        text_page=text_page.encode('utf8')               
                    
                    #-------------------------------------Tokenize--------------------------------------#
                        word_token=nltk.tokenize.word_tokenize(str(text_page))
                        word_token_all=(word_token_all+word_token)
                         
                        if i==(n-1):
                            
                            display_box.insert(tk.END,'\nEnd of reading, initializing text tokenize')
                            
                            #add keyword - chemical informations
                            processed_text_6=[];
                            processed_text_27_1=[];unity_concen=['mol','mmmol','pmol','nmol','M','mM','pM','nM','dm']
                            processed_text_27_2=[];processed_text_35=[];processed_text_36=[];processed_text_37=[]
                            pontuaction=[':','/','=',',',';','!','*','&','(',')','[',']','{','}','@','$','<','>','|','?','%','and','of','p.A']
                            
                            #concentration
                            for kw1 in range(0,50) :
                                processed_text_27_1.append(word_token_all[kw1])
                            
                            for kw2 in range(0,400):
                                if word_token_all[kw2] not in pontuaction and research_number==1:
                                    processed_text_35.append(word_token_all[kw2])
                                    
                            for kw22 in range(400,len(word_token_all)):
                                if word_token_all[kw22] in unity_concen and research_number==1:
                                    processed_text_35.append(word_token_all[kw22-3]);
                                    processed_text_35.append(word_token_all[kw22-2])
                                    processed_text_35.append(word_token_all[kw22-1])
                                    processed_text_35.append(word_token_all[kw22])
                          
                            for kw3 in range(0,len(processed_text_35)):
                                if kw3 < (len(processed_text_35)) and research_number==1:
                                    if (bool(re.search(r'\d', processed_text_35[kw3])))==True and str(processed_text_35[kw3])!=None and (processed_text_35[kw3])!='': 
                                        if '-' in  processed_text_35[kw3]:  
                                            if (bool(re.search(r'\d', processed_text_35[kw3]))):
                                                    processed_text_37=processed_text_35[kw3].split('-')
                                                    if len(processed_text_37)==2:
                                                        processed_text_35[kw3-1]=''
                                                        processed_text_35[kw3-1],processed_text_35[kw3]=processed_text_37[0],processed_text_37[1]
                                    if 'M.' in processed_text_35[kw3] or 'mol.' in processed_text_35[kw3]:
                                        processed_text_38=(str(processed_text_35[kw3]).split('.'))
                                        processed_text_35[kw3]=processed_text_38[0]
                                      
                                    if 'x' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split('x'))
                                        processed_text_35[kw3]=processed_text_38[0]
                                    if 'X' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split('X'))
                                        processed_text_35[kw3]=processed_text_38[0]
                                    if 'limit' in processed_text_35[kw3]:
                                        del(processed_text_35[kw3])
                                    if ',' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split(','))
                                        processed_text_35[kw3]=processed_text_38[0]  
                                    if '=' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split('='))
                                        processed_text_35[kw3]=processed_text_38[0]
                                    if '/' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split('/'))
                                        processed_text_35[kw3]=processed_text_38[0]
                                    if ':' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split(':'))
                                        processed_text_35[kw3]=processed_text_38[0]
                                    if '^' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split('^'))
                                        processed_text_35[kw3]=processed_text_38[0]
                                    if '*' in processed_text_35[kw3] :
                                        processed_text_38=(str(processed_text_35[kw3]).split('*'))
                                        processed_text_35[kw3]=processed_text_38[0]
                            
                            for kw4 in range(0,len(processed_text_35)):
                                if processed_text_35[kw4] in unity_concen and research_number==1:
                                    if processed_text_35[kw4]=='M' or processed_text_35[kw4]=='mol' :
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-2]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-2]))==False and (bool(re.search(r'\d', processed_text_35[kw4-2]))) :
                                            processed_text_36.append(float(processed_text_35[kw4-2]))
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-1]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-1]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-2]))==False and (bool(re.search(r'\d', processed_text_35[kw4-1])))  :
                                            processed_text_36.append(float(processed_text_35[kw4-1]))   
                                        if processed_text_35[kw4-2]=='10' or processed_text_35[kw4-2]=='x10' or processed_text_35[kw4-2]=='*10' or processed_text_35[kw4-2]=='X10':
                                            if bool(re.search(r'[a-y]',processed_text_35[kw4-3]))==False and bool(re.search(r'[a-y]',processed_text_35[kw4-1]))==False and (bool(re.search(r'\d', processed_text_35[kw4-3]))) and  processed_text_35[kw4-3] not in pontuaction:
                                                processed_text_36.append(float(processed_text_35[kw4-3])*10**-(int(processed_text_35[kw4-1])))
                                                      
                                    if processed_text_35[kw4]=='mM' or processed_text_35[kw4]=='mmol' or processed_text_35[kw4]=='dm':
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-2]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-2]))==False and (bool(re.search(r'\d', processed_text_35[kw4-2])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-2]))==False :
                                            processed_text_36.append(round(((float(processed_text_35[kw4-2]))*10**-3),5))  
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-1]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-1]))==False and (bool(re.search(r'\d', processed_text_35[kw4-1])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-1]))==False:
                                            processed_text_36.append(round(((float(processed_text_35[kw4-1]))*10**-3),5)) 
                                                              
                                    if processed_text_35[kw4]=='uM' or processed_text_35[kw4]=='umol':
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-2]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-2]))==False and (bool(re.search(r'\d', processed_text_35[kw4-2])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-2]))==False :
                                            processed_text_36.append(round(((float(processed_text_35[kw4-2]))*10**-6),8))   
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-1]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-1]))==False and (bool(re.search(r'\d', processed_text_35[kw4-1])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-1]))==False:
                                            processed_text_36.append(round(((float(processed_text_35[kw4-1]))*10**-6),8))
                                                            
                                    if processed_text_35[kw4]=='nM' or processed_text_35[kw4]=='nmol':
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-2]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-2]))==False and (bool(re.search(r'\d', processed_text_35[kw4-2])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-2]))==False:
                                            processed_text_36.append(round(((float(processed_text_35[kw4-2]))*10**-9),11))   
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-1]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-1]))==False and (bool(re.search(r'\d', processed_text_35[kw4-1])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-1]))==False:
                                            processed_text_36.append(round(((float(processed_text_35[kw4-1]))*10**-9),11))
                                    
                                    if processed_text_35[kw4]=='pM' or processed_text_35[kw4]=='pmol':
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-2]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-2]))==False and (bool(re.search(r'\d', processed_text_35[kw4-2])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-2]))==False :
                                            processed_text_36.append(round(((float(processed_text_35[kw4-2]))*10**-12),14))   
                                        if bool(re.search(r'[a-z]',processed_text_35[kw4-1]))==False and bool(re.search(r'[A-Z]',processed_text_35[kw4-1]))==False and (bool(re.search(r'\d', processed_text_35[kw4-1])))==True and bool(re.search(r'[:/=;,!*^&()[]{}@$<>|%-+}#]',processed_text_35[kw4-1]))==False:
                                            processed_text_36.append(round(((float(processed_text_35[kw4-1]))*10**-12),14))
                                if processed_text_36==[] and research_number==1:
                                    processed_text_36.append(0)

                            if processed_text_36==[]:
                                processed_text_36.append(0)
                            
                    
                            if 'Keywords' in word_token_all:
                                kw=word_token_all.index('Keywords')
                                processed_text_27_2=(word_token_all[kw:kw+15]);
                            processed_text_6=processed_text_27_1+processed_text_27_2
                            
                            #removed Abstract
                            if 'Abstract' in word_token_all:
                                a=word_token_all.index('Abstract')
                                del(word_token_all[0:(a+1)])
                            if 'ABSTRACT' in word_token_all:
                                a=word_token_all.index('ABSTRACT')
                                del(word_token_all[0:(a+1)])
                            
                            if 'Introduction' in word_token_all:
                                key_1=0
                                #removed Introduction
                                if 'Introduction' and 'Experimental' in word_token_all:
                                    b=word_token_all.index('Introduction')
                                    c=word_token_all.index('Experimental')
                                    del(word_token_all[b:(c)])
                                    key_1=1
                                if 'Introduction' and 'Experiment' in word_token_all:
                                    b_0=word_token_all.index('Introduction')
                                    c_0=word_token_all.index('Experiment')
                                    del(word_token_all[b_0:(c_0)])
                                    key_1=1
                                
                                if key_1==0 :
                                    if 'Materials'in word_token_all:
                                        #removed Introduction
                                        if 'Introduction' and 'Materials' and 'methods' in word_token_all:
                                            b_1=word_token_all.index('Introduction')
                                            c_1=word_token_all.index('Materials')
                                            c_2=word_token_all.index('methods')
                                            if c_2-c_1==2 :
                                                del(word_token_all[b_1:(c_2+1)])
                                                 
                                        if 'Introduction' and 'methods' in word_token_all:
                                            b_2=word_token_all.index('methods')
                                            del(word_token_all[0:(b_2)])

                                        else:
                                            if 'Introduction' in word_token_all:
                                                b_3=word_token_all.index('Introduction')
                                                del(word_token_all[0:b_3])
                                              
                            #removal of camouflaged introduction               
                            if 'EXPERIMENTAL' or 'Experimental' in word_token_all:
                                if 'EXPERIMENTAL'and 'SECTION' in  word_token_all:
                                    c_3=word_token_all.index('EXPERIMENTAL')
                                    del(word_token_all[0:(c_3+2)])
                                if 'Experimental' in word_token_all:
                                    c_4=word_token_all.index('Experimental')
                                    del(word_token_all[0:(c_4)])
                                    
                            #removed References
                            if 'References' in word_token_all:
                                d=word_token_all.index('References')
                                y=len(word_token_all)
                                del(word_token_all[d:y])
                                
                            #removed Acknowledgments
                            if 'Acknowledgment' in word_token_all:
                                e=word_token_all.index('Acknowledgment')
                                y=len(word_token_all)
                                del(word_token_all[e:y])
                                
                            #removed Acknowledgments
                            if 'Acknowledgments' in word_token_all:
                                e=word_token_all.index('Acknowledgments')
                                y=len(word_token_all)
                                del(word_token_all[e:y])
                                
                            if 'Conflicts' in word_token_all:
                            #removed Conflicts of interest
                                if 'Conflicts'  and 'of' and 'interest' in word_token_all:
                                    f=word_token_all.index('Conflicts')
                                    h=word_token_all.index('interest')
                                    y=len(word_token_all)
                                    if h-f==2 :
                                        del(word_token_all[(f-1):y])

                            if 'Author' in word_token_all:
                                if 'Author' and 'INFORMATION' in word_token_all:
                                    h_1=word_token_all.index('Author')
                                    y=len(word_token_all)
                                    del(word_token_all[h_1:y])
                            
                            if 'ASSOCIATED' in word_token_all:
                                if 'ASSOCIATED' and 'CONTENT' in word_token_all:
                                    h_2=word_token_all.index('ASSOCIATED')
                                    y=len(word_token_all)
                                    del(word_token_all[h_2:y])
                            
                            if 'Figure' in word_token_all:
                                h_3=word_token_all.index('Figure')
                                del(word_token_all[h_3:h_3+1])
                                    
                            #removing stop words
                            stop_words=set(stopwords.words('english'))
                            processed_text=[]
                            for j in  word_token_all:
                                    if j not in stop_words:
                                        processed_text.append(j)  
                
                            #removing pontuaction
                            pontuaction=[':','/','=',',',';','!','*','&','.','(',')','[',']','{','}','@','$','<','>','|','?','%']
                            processed_text_2=[]
                            for l in  processed_text:
                                    if l not in pontuaction:
                                        processed_text_2.append(l) 
                             
                            #removing  manufacturer and others
                            manufacture=['mixture','Metrohm-Autolab','Autolab','Germany','g','Millipore','discussion','Sigma', 'Aldrich','Sigma-Aldrich','Merck','Dinamica','Fluka,','Synth','et','al','Autolab','Materials','Figure','ACCEPTED','MANUSCRIPT','in','In','UK','USA','BR','buffer','solutions','solution','All','Henrifarma']
                            processed_text_3=[]
                            for m in  processed_text_2:
                                    if m not in manufacture:
                                        processed_text_3.append(m)
                                       
                            UNITY=['mg','mm','ml','ul','SS','mL','uL','cm','um','Km','L','h','min','s','mg/mL','ug/mL','mg/L','a','A','V','mV','mA','uA','Hz']
                            processed_text_3_1=[]
                            for m01 in  processed_text_3:
                                    if m01 not in UNITY:
                                        processed_text_3_1.append(m01)
                            
                            if tokenchemical_ini==1:
                                chem00=reagent+solvent+electrode_modifield
                                processed_text_4=[]
                                processed_text_5=[]
                                for m0 in processed_text_3_1:
                                    if m0 in chem00:
                                        processed_text_4.append(m0)
                                    if m0 not in chem00:
                                        processed_text_5.append(m0)
                                if 'Results' in  processed_text_5:
                                    chem01=processed_text_5.index('Results')
                                    chem02=len(processed_text_5)
                                    del(processed_text_5[chem01:chem02])  
                                if 'RESULTS' in  processed_text_5:
                                    chem03=processed_text_5.index('RESULTS')
                                    chem04=len(processed_text_5)
                                    del(processed_text_5[chem03:chem04])     
                                for par in nltk.pos_tag(processed_text_5):
                                    if par[1] in ['NN','NNP','NNS','NNPS']:
                                        processed_text_6.append(par[0])
                                
                                #removing repeat
                                processed_text_6_1=[]
                                processed_text_6_1=list(set(processed_text_6))
                                processed_text_6_2=[]
                                for m001 in processed_text_6_1:
                                    if m001 not in chemical_elements_name:
                                            processed_text_6_2.append(m001)
                                
                                processed_text_6_2_1=[]
                                for m002 in processed_text_6_2:
                                    if m002 not in chemical_elements:
                                            processed_text_6_2_1.append(m002)
#--------------------------------------------------------analyt-------------------------------------------------------------------------------------------------                                            
                                
                                processed_text_7=[]
                                acid_org1=['gallic','ascorbic','acetic','glutamic','aspartic','uric','latic','malic,','citric','oxalic','tartaric','salicylic']    
                                acid_org2=['Gallic','Ascorbic','Acetic','Glutamic','Aspartic','Uric','Latic','Malic,','Citric','Oxalic','Tartaric','Salicylic']    
                                for m002_2 in processed_text_6_2_1:
                                    if m002_2 in acid_org1 or m002_2 in acid_org2:
                                        acid_pos=processed_text_6_2_1.index(m002_2)
                                        processed_text_6_2_1[acid_pos]=m002_2+" "+'acid'

                                #Seacrh to pesticide
                                file_FDAA=os.path.dirname(__file__)
                                list_drug=pd.read_excel(file_FDAA+"\\Pesticide.xlsx",engine='openpyxl')
                                for m3 in  processed_text_6_2_1: 
                                    search_drug_5=list_drug[list_drug['Pesticide'].isin([m3])].any().any()
                                    if search_drug_5:
                                        processed_text_7.append(m3)
                                
                                #remove siglas
                                processed_text_6_3=[]   
                                m003=len(processed_text_6_2_1)
                                for m004 in range(0,m003):
                                    if (len(processed_text_6_2_1[m004])) > 3:
                                        processed_text_6_3.append(processed_text_6_2_1[m004]) 
                                
                                if 'bisphenol' in processed_text_6_3:
                                    processed_text_7.append('bisphenol')
                                
                                #Search for peroxide
                                H2O2=['H2O2','peroxide']
                                if H2O2 in processed_text_6_3:
                                    processed_text_7.append('H2O2')
                                
                                #Search for drug
                                file_FDAA=os.path.dirname(__file__)
                                list_drug=pd.read_excel(file_FDAA+"\\Drug.xlsx",engine='openpyxl')
                                for m2 in  processed_text_6_3: 
                                    search_drug_1=list_drug[list_drug['DrugName'].isin([m2])].any().any()
                                    search_drug_2=list_drug[list_drug['ActiveIngredient'].isin([m2])].any().any()
                                    search_drug_3=list_drug[list_drug['DRUG_NAME'].isin([m2])].any().any()
                                    search_drug_4=list_drug[list_drug['TARGET_NAME'].isin([m2])].any().any()
                                    if search_drug_1 or search_drug_2 or search_drug_3 or search_drug_4 :
                                        processed_text_7.append(m2)                                   
                                
                                #Search of Biomolecule
                                list_biomol=pd.read_excel(file_FDAA+"\\COCONUT.xlsx",engine='openpyxl')
                                for m3 in processed_text_6_3:
                                    chem3=pcp.get_compounds(m3,'name')
                                    chem4=len(chem3)
                                    chem5=False
                                    if chem5==False:
                                        for m4 in range(0,chem4):
                                            SMILE_BIOMOL=chem3[m4].canonical_smiles
                                            search_biomol_1=list_biomol[list_biomol['SMILE1'].isin([SMILE_BIOMOL])].any().any()
                                            search_biomol_2=list_biomol[list_biomol['SMILE2'].isin([SMILE_BIOMOL])].any().any()                                        
                                            if search_biomol_1 or search_biomol_2:
                                                processed_text_7.append(m3)
                                                chem5=True
                                
                                processed_text_29=[]
                                processed_text_30=[]
                                processed_text_31=[]
                                processed_text_34=[];processed_text_34_1=[];processed_text_34_3=[];

                                for m32 in range (0,len(processed_text_7)):
                                    processed_text_34.append((str(processed_text_7[m32])).lower())
                                if 'phosphate' in  processed_text_34:
                                    del(processed_text_34[processed_text_34.index('phosphate')])
                                if 'tetrabutylammonium' in  processed_text_34:
                                    del(processed_text_34[processed_text_34.index('tetrabutylammonium')])
                                for m34 in processed_text_34:
                                    if m34 not in chemical_elements_name:
                                        processed_text_34_1.append(m34)
                                processed_text_34_2=list(set(processed_text_34_1))     
                                
                                for m33 in processed_text_34_2 :
                                    chem1=pcp.get_compounds(m33, 'name')
                                    if len(chem1)!=0:
                                        processed_text_34_3.append(chem1[0].cid)
                                processed_text_34_4=list(set(processed_text_34_3)) 
                               
                                for m30 in processed_text_34_4 :
                                    chem2=Compound.from_cid(m30)
                                    processed_text_29.append(chem2.iupac_name)
                                    processed_text_31.append(chem2.canonical_smiles)
                                    processed_text_30.append(chem2.xlogp) 
                                
                                if None in processed_text_30:
                                    for m34 in range(0,len(processed_text_30)):
                                        if (str(processed_text_30[m34]))=='None' :
                                            processed_text_30[m34]=-1
                            display_box.insert(tk.END,'\nEnd of Tokenize'+'\nExtracting information')
                            display_box.insert(tk.END,'\n')      
#-------------------------------------------------------------------------------------------------------------------------------------------                   
                         
                            for q in range(0,o):
                                key_start=0

                                if word_token_attribute[q] in library_Lr:
                                    write_text=" ".join(processed_text_3)
                                    write_text=unidecode.unidecode(write_text)
                                    save_arff.writelines(write_text+',')
                                    key_start=1
                                                  
                                if word_token_attribute[q] in library_Am:  
                                    if library_Am[0] in word_token_attribute[q]:
                                        key_start=1;scan_rate=0
                                        AM1=0;phenol_stop=0;hca_stop=0;thiol_stop=0;amine_stop=0;thiocarbamate_stop=0;sulfonamide_stop=0;imine_stop=0;aldehyde_stop=0 
                                    if library_Am[1] in word_token_attribute[q]:
                                        key_start=1;scan_rate=0
                                        AM1=1;
                                        phenol_stop=0;hca_stop=0;thiol_stop=0;amine_stop=0;thiocarbamate_stop=0;sulfonamide_stop=0;imine_stop=0;aldehyde_stop=0 
                                    
                                    for ph001 in processed_text_31:
                                        chem_token=[];
                                        chem_token=nltk.tokenize.word_tokenize(ph001) 
                                        ph002=len(chem_token)
                                        ph004=0
                                        phenol_stop=0;scan_rate=0
                                        for ph003 in range(0,ph002):
                                            if '=C' in str(chem_token[ph003]):
                                                ph004=1; 
                                            if ph004==1 and phenol_stop==0:
                                                ph005=len(chem_token)
                                                
                                                if chem_token[ph005-2]==')' and chem_token[ph005-1]=='O' :
                                                    if (library_Am[2] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q]):                           
                                                        save_arff.writelines( 'phenol/enol'+',');
                                                        key_start=1;scan_rate=1
                                                        processed_text_32=(processed_text_31.index(ph001))
                                                        if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1) and space_analyt==1:
                                                                save_arff.writelines 
                                                                key_start=1
                                                        if analyt_soluble==1:
                                                            if processed_text_30[processed_text_32] < 0:
                                                                save_arff.writelines( 'Water'+',');
                                                                save_arff.writelines( 'KCl'+',')
                                                                key_start=1
                                                            if processed_text_30[processed_text_32] >= 0:
                                                                save_arff.writelines( 'THF'+',')   
                                                                save_arff.writelines( 'LiClO4'+',')
                                                                key_start=1
                                                                    
                                                        if  library_No[0] in word_token_attribute:
                                                            position_scan=[];
                                                            for s_r in range(0,len( processed_text_3)):
                                                                if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                    for s_r2 in range(0,7):
                                                                        if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                            position_scan.append(s_r);key_start=1;posis_r=s_r2
                                                                            s_r2=8
                                                            if position_scan==[]:
                                                                position_scan=[0];
                                                            if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                Volt='E0,';
                                                            elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                Volt='E-3,';
                                                            else:
                                                                Volt='E0,';
                                                            
                                                            save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                        if  library_No[1] in word_token_attribute:
                                                            save_arff.writelines(str(min(processed_text_36))+',')
                                                        if name_article==1:
                                                            save_arff.writelines(str(pdf_n+','))
                                                        if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                            save_arff.writelines('\n') 
                                                        key_start=1
                                                        phenol_stop=1
                                            if phenol_stop==0 and ph003==(ph002-1) and AM1==0 and (library_Am[2] in word_token_attribute[q]):    
                                                save_arff.writelines( 'No'+','+'No'+','+'No'+',');key_start=1
                                                if  library_No[0] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0,');
                                                if  library_No[1] in word_token_attribute:
                                                    if scan_rate==0:
                                                        save_arff.writelines('0,');
                                                if name_article==1:
                                                    save_arff.writelines(str(pdf_n+','))
                                                if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                    save_arff.writelines('\n')
                                                phenol_stop=1;key_start=1
                                                
                                    for hca001 in processed_text_29:               
                                        if phenol_stop==0 :
                                            hca003=len(het_cycle)
                                            hca004=processed_text_29[(len(processed_text_29))-1]
                                            hca_stop=0;scan_rate=0
                                            for hc002 in range(0,hca003):
                                                if hca_stop==0 :
                                                    if str(het_cycle[hc002]) in str(hca001):
                                                        if (library_Am[3] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q])  :                           
                                                            save_arff.writelines( 'het_cyclic_aro'+',')
                                                            key_start=1;scan_rate=1
                                                            if  len(processed_text_29) > 1 and (processed_text_29.index(hca001)) < ((len(processed_text_29))-1) and space_analyt==1:
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            processed_text_32=(processed_text_29.index(hca001))
                                                            if analyt_soluble==1:
                                                                if processed_text_30[processed_text_32] < 0:
                                                                    save_arff.writelines( 'Water'+',')
                                                                    save_arff.writelines( 'KCl'+',')
                                                                if processed_text_30[processed_text_32] >= 0:
                                                                    save_arff.writelines( 'THF'+',')   
                                                                    save_arff.writelines( 'LiClO4'+',')
                                                                key_start=1
                                                                
                                                            if  library_No[0] in word_token_attribute:
                                                                position_scan=[];key_start=1;
                                                                for s_r in range(0,len( processed_text_3)):
                                                                    if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                        for s_r2 in range(0,7):
                                                                            if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                                position_scan.append(s_r);key_start=1;posis_r=s_r2
                                                                                s_r2=8
                                                                if position_scan==[]:
                                                                    position_scan=[0];
                                                                if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                    Volt='E0,';
                                                                elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                    Volt='E-3,';
                                                                else:
                                                                    Volt='E0,';
                                                                
                                                                save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                            if  library_No[1] in word_token_attribute:
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                            if name_article==1:
                                                                save_arff.writelines(str(pdf_n+','))
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                                save_arff.writelines('\n') 
                                                            key_start=1 
                                                            
                                                            hca_stop=1                  
                                                    if hca_stop==0 and hca001==(hca004) and AM1==0 and (library_Am[3] in word_token_attribute[q]):    
                                                        save_arff.writelines( 'No'+','+'No'+','+'No'+',')  
                                                        if  library_No[0] in word_token_attribute:
                                                            if  scan_rate==0:
                                                                save_arff.writelines('0');
                                                        if  library_No[1] in word_token_attribute:
                                                            if  scan_rate==0:
                                                                save_arff.writelines('0,');
                                                        if name_article==1:
                                                            save_arff.writelines(str(pdf_n+','))
                                                        if  len(processed_text_29) > 1 and (processed_text_29.index(hca001))  < ((len(processed_text_29))-1) :
                                                            save_arff.writelines('\n') 
                                                            key_start=1
                                                        hca_stop=1  
                                    
                                    for th001 in processed_text_31:
                                        if phenol_stop==0 and  hca_stop==0:
                                            chem_token=[];key_start=1
                                            chem_token=nltk.tokenize.word_tokenize(th001) 
                                            th002=len(chem_token)
                                            th004=0;th005=0;
                                            thiol_stop=0;scan_rate=0
                                            for th003 in range(0,th002):
                                                if '=S' in str(chem_token[th003]):
                                                    th004=1   
                                                if 'NS' in str(chem_token[th003]):
                                                    th005=1     
                                                if th004==0 and th005==0  and thiol_stop==0:
                                                    if 'S' in str(chem_token[th003]):
                                                        if (library_Am[4] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q]):                           
                                                            save_arff.writelines( 'thiophenol/thiol'+',')
                                                            key_start=1;scan_rate=1
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(th001)) != ((len(processed_text_31))-1) and space_analyt==1:
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            processed_text_32=(processed_text_31.index(th001))
                                                            if analyt_soluble==1:
                                                                if processed_text_30[processed_text_32] < 0:
                                                                    save_arff.writelines( 'Water'+',')
                                                                    save_arff.writelines( 'KCl'+',')
                                                                    key_start=1
                                                                if processed_text_30[processed_text_32] >= 0:
                                                                    save_arff.writelines( 'THF'+',')   
                                                                    save_arff.writelines( 'LiClO4'+',')
                                                                    key_start=1
                                                                
                                                            
                                                            if  library_No[0] in word_token_attribute:
                                                                position_scan=[];
                                                                for s_r in range(0,len( processed_text_3)):
                                                                    if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                        for s_r2 in range(0,7):
                                                                            if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                                position_scan.append(s_r);key_start=1;posis_r=s_r2
                                                                                s_r2=8
                                                                if position_scan==[]:
                                                                    position_scan=[0];
                                                                if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                    Volt='E0,';
                                                                elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                    Volt='E-3,';
                                                                else:
                                                                    Volt='E0,';
                                                                
                                                                save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                            if  library_No[1] in word_token_attribute:
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                            if name_article==1:
                                                                save_arff.writelines(str(pdf_n+','))
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            thiol_stop=1
                                            if thiol_stop==0 and th003==(th002-1) and AM1==0 and (library_Am[4] in word_token_attribute[q]):    
                                                save_arff.writelines( 'No'+','+'No'+','+'No'+',')  
                                                key_start=1
                                                if  library_No[0] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0');
                                                if  library_No[1] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0,');
                                                if name_article==1:
                                                    save_arff.writelines(str(pdf_n+','))
                                                if  len(processed_text_31) > 1 and (processed_text_31.index(th001)) != ((len(processed_text_31))-1):
                                                    save_arff.writelines('\n')
                                                thiol_stop=1           
                               
                                    for am001 in processed_text_31  :
                                        if phenol_stop==0 and  hca_stop==0 and thiol_stop==0:
                                            chem_token=[]
                                            chem_token=nltk.tokenize.word_tokenize(am001) 
                                            am002=len(chem_token)
                                            am004=0;am005=0;am006=0;
                                            amine_stop=0;scan_rate=0
                                            for am003 in range(0,am002):
                                                if '=N' in str(chem_token[am003]):
                                                    am004=1   
                                                if 'NS' in str(chem_token[am003]):
                                                    am005=1  
                                                if '#N' in str(chem_token[am003]):
                                                    am006=1        
                                                if am004==0 and am005==0  and am006==0 and amine_stop==0:
                                                    if 'N' in str(chem_token[am003]):
                                                        if (library_Am[5] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q]):                           
                                                            save_arff.writelines( 'amine'+',')
                                                            key_start=1;scan_rate=1
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(am001)) != ((len(processed_text_31))-1) and space_analyt==1:
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            processed_text_32=(processed_text_31.index(am001))
                                                            if analyt_soluble==1:
                                                                if processed_text_30[processed_text_32] < 0:
                                                                    save_arff.writelines( 'Water'+',')
                                                                    save_arff.writelines( 'KCl'+',')
                                                                    key_start=1
                                                                if processed_text_30[processed_text_32] >= 0:
                                                                    save_arff.writelines( 'THF'+',')   
                                                                    save_arff.writelines( 'LiClO4'+',')
                                                                    key_start=1
                                                                if  len(processed_text_31) > 1 and (processed_text_31.index(am001)) != ((len(processed_text_31))-1):
                                                                    save_arff.writelines('\n') 
                                                                    key_start=1
                                                            
                                                            if  library_No[0] in word_token_attribute:
                                                                position_scan=[];
                                                                for s_r in range(0,len( processed_text_3)):
                                                                    if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                        for s_r2 in range(0,7):
                                                                            if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                                position_scan.append(s_r);key_start=1;posis_r=s_r2
                                                                                s_r2=8
                                                                if position_scan==[]:
                                                                    position_scan=[0];
                                                                if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                    Volt='E0,';
                                                                elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                    Volt='E-3,';
                                                                else:
                                                                    Volt='E0,';
                                                                
                                                                save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                            if  library_No[1] in word_token_attribute:
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                            if name_article==1:
                                                                save_arff.writelines(str(pdf_n+','))
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            amine_stop=1
                                            if  amine_stop==0 and am003==(am002-1) and AM1==0 and (library_Am[5] in word_token_attribute[q]):    
                                                save_arff.writelines( 'No'+','+'No'+','+'No'+',')  
                                                key_start=1
                                                if  library_No[0] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0');
                                                if  library_No[1] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0,');
                                                if name_article==1:
                                                    save_arff.writelines(str(pdf_n+','))
                                                if len(processed_text_31) > 1 and (processed_text_31.index(am001)) != ((len(processed_text_31))-1):
                                                    save_arff.writelines('\n')
                                                amine_stop=1             
                                    
                                    for tc001 in processed_text_31  :
                                        if phenol_stop==0 and  hca_stop==0 and thiol_stop==0 and amine_stop==0:
                                            chem_token=[]
                                            chem_token=nltk.tokenize.word_tokenize(tc001) 
                                            tc002=len(chem_token)
                                            tc004=0;tc005=0;tc006=0;
                                            thiocarbamate_stop=0;scan_rate=0
                                            for tc003 in range(0,tc002):
                                                if '=S' in str(chem_token[tc003]):
                                                    tc004=1   
                                                if 'N' in str(chem_token[tc003]):
                                                    tc005=1  
                                                if 'O' in str(chem_token[tc003]):
                                                    tc006=1        
                                                if tc004==1 and tc005==1  and tc006==1 and thiocarbamate_stop==0:
                                                    if (library_Am[6] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q]):                           
                                                            save_arff.writelines( 'thiocarbamate'+',')
                                                            key_start=1;scan_rate=1
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(tc001)) != ((len(processed_text_31))-1) and space_analyt==1:
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            processed_text_32=(processed_text_31.index(tc001))
                                                            if analyt_soluble==1:
                                                                if processed_text_30[processed_text_32] < 0:
                                                                    save_arff.writelines( 'Water'+',')
                                                                    save_arff.writelines( 'KCl'+',')
                                                                    key_start=1
                                                                if processed_text_30[processed_text_32] >= 0:
                                                                    save_arff.writelines( 'THF'+',')   
                                                                    save_arff.writelines( 'LiClO4'+',')
                                                                    key_start=1
                                                                if  len(processed_text_31) > 1 and (processed_text_31.index(tc001)) != ((len(processed_text_31))-1):
                                                                    save_arff.writelines('\n') 
                                                                    key_start=1
                                                            
                                                            if  library_No[0] in word_token_attribute:
                                                                position_scan=[];
                                                                for s_r in range(0,len( processed_text_3)):
                                                                    if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                        for s_r2 in range(0,7):
                                                                            if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                                position_scan.append(s_r);key_start=1;posis_r=s_r2
                                                                                s_r2=8
                                                                if position_scan==[]:
                                                                    position_scan=[0];
                                                                if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                    Volt='E0,';
                                                                elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                    Volt='E-3,';
                                                                else:
                                                                    Volt='E0,';
                                                                
                                                                save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                            if  library_No[1] in word_token_attribute:
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                            if name_article==1:
                                                                save_arff.writelines(str(pdf_n+','))
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            thiocarbamate_stop=1       
                                            if  thiocarbamate_stop==0 and tc003==(tc002-1) and AM1==0 and (library_Am[6] in word_token_attribute[q]):    
                                                save_arff.writelines( 'No'+','+'No'+','+'No'+',')  
                                                key_start=1
                                                if  library_No[0] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0');
                                                if  library_No[1] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0,');
                                                if name_article==1:
                                                    save_arff.writelines(str(pdf_n+','))
                                                if  len(processed_text_31) > 1 and (processed_text_31.index(tc001)) != ((len(processed_text_31))-1):
                                                    save_arff.writelines('\n')
                                                thiocarbamate_stop=1                   
                                   
                                    for sfa001 in processed_text_31  :
                                        if phenol_stop==0 and  hca_stop==0 and thiol_stop==0 and amine_stop==0 and thiocarbamate_stop==0   :
                                            chem_token=[]
                                            chem_token=nltk.tokenize.word_tokenize(sfa001) 
                                            sfa002=len(chem_token)
                                            sulfonamide_stop=0;scan_rate=0
                                            for sfa003 in range(0,sfa002):        
                                                if 'NS' in str(chem_token[sfa003]) and  thiocarbamate_stop==0:
                                                    if (library_Am[7] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q]):                           
                                                            save_arff.writelines( 'sulfonamide'+',')
                                                            key_start=1;scan_rate=1
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(sfa001)) != ((len(processed_text_31))-1) and space_analyt==1:
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            processed_text_32=(processed_text_31.index(sfa001))
                                                            if analyt_soluble==1:
                                                                if processed_text_30[processed_text_32] < 0:
                                                                    save_arff.writelines( 'Water'+',')
                                                                    save_arff.writelines( 'KCl'+',')
                                                                    key_start=1
                                                                if processed_text_30[processed_text_32] >= 0:
                                                                    save_arff.writelines( 'THF'+',')   
                                                                    save_arff.writelines( 'LiClO4'+',')
                                                                    key_start=1
                                                                if  len(processed_text_31) > 1 and (processed_text_31.index(sfa001)) != ((len(processed_text_31))-1):
                                                                    save_arff.writelines('\n') 
                                                                    key_start=1
                                                            
                                                            if  library_No[0] in word_token_attribute:
                                                                position_scan=[];
                                                                for s_r in range(0,len( processed_text_3)):
                                                                    if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                        for s_r2 in range(0,7):
                                                                            if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                                position_scan.append(s_r);key_start=1;posis_r=s_r2
                                                                                s_r2=8
                                                                if position_scan==[]:
                                                                    position_scan=[0];
                                                                if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                    Volt='E0,';
                                                                elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                    Volt='E-3,';
                                                                else:
                                                                    Volt='E0,';
                                                                
                                                                save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                            if  library_No[1] in word_token_attribute:
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                            if name_article==1:
                                                                save_arff.writelines(str(pdf_n+','))
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                                save_arff.writelines('\n') 
                                                            key_start=1
                                                            sulfonamide_stop=1       
                                            if   sulfonamide_stop==0 and sfa003==(sfa002-1) and AM1==0 and (library_Am[7] in word_token_attribute[q]):    
                                                save_arff.writelines( 'No'+','+'No'+','+'No'+',')  
                                                key_start=1
                                                if  library_No[0] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0');
                                                if  library_No[1] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0,');
                                                if name_article==1:
                                                    save_arff.writelines(str(pdf_n+','))
                                                if  len(processed_text_31) > 1 and (processed_text_31.index(sfa001)) != ((len(processed_text_31))-1):
                                                    save_arff.writelines('\n')
                                                sulfonamide_stop=1     
                                    
                                    for imi001 in processed_text_31  :
                                        if phenol_stop==0 and  hca_stop==0 and thiol_stop==0 and amine_stop==0 and thiocarbamate_stop==0 and sulfonamide_stop==0  :
                                            chem_token=[];key_start=1
                                            chem_token=nltk.tokenize.word_tokenize(imi001) 
                                            imi002=len(chem_token)
                                            imine_stop=0;scan_rate=0
                                            for imi003 in range(0,imi002):        
                                                if 'N=C' in str(chem_token[imi003]) and  imine_stop==0:
                                                    if (library_Am[8] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q]):                           
                                                            save_arff.writelines( 'imine'+',')
                                                            key_start=1;scan_rate=1
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(imi001)) != ((len(processed_text_31))-1) and space_analyt==1:
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            processed_text_32=(processed_text_31.index(imi001))
                                                            if analyt_soluble==1:
                                                                if processed_text_30[processed_text_32] < 0:
                                                                    save_arff.writelines( 'Water'+',')
                                                                    save_arff.writelines( 'KCl'+',')
                                                                    key_start=1
                                                                if processed_text_30[processed_text_32] >= 0:
                                                                    save_arff.writelines( 'THF'+',')   
                                                                    save_arff.writelines( 'LiClO4'+',')
                                                                    key_start=1
                                                                if  len(processed_text_31) > 1 and (processed_text_31.index(imi001)) != ((len(processed_text_31))-1):
                                                                    save_arff.writelines('\n') 
                                                                    key_start=1
                                                            
                                                            if  library_No[0] in word_token_attribute:
                                                                position_scan=[];
                                                                for s_r in range(0,len( processed_text_3)):
                                                                    if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                        for s_r2 in range(0,7):
                                                                            if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                                position_scan.append(s_r);key_start=1;posis_r=s_r2
                                                                                s_r2=8
                                                                if position_scan==[]:
                                                                    position_scan=[0];
                                                                if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                    Volt='E0,';
                                                                elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                    Volt='E-3,';
                                                                else:
                                                                    Volt='E0,';
                                                                
                                                                save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                            if  library_No[1] in word_token_attribute:
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                            if name_article==1:
                                                                save_arff.writelines(str(pdf_n+','))
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                                save_arff.writelines('\n') 
                                                                key_start=1
                                                            imine_stop=1       
                                            if   imine_stop==0 and imi003==(imi002-1) and AM1==0 and (library_Am[8] in word_token_attribute[q]):    
                                                save_arff.writelines( 'No'+','+'No'+','+'No'+',')  
                                                key_start=1
                                                if  library_No[0] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0');
                                                if  library_No[1] in word_token_attribute:
                                                    if  scan_rate==0:
                                                        save_arff.writelines('0,');
                                                if name_article==1:
                                                    save_arff.writelines(str(pdf_n+','))
                                                if  len(processed_text_31) > 1 and (processed_text_31.index(imi001)) != ((len(processed_text_31))-1):
                                                    save_arff.writelines('\n')
                                                imine_stop=1;key_start=1  
                                    
                                    for ald001 in processed_text_31  :
                                        if phenol_stop==0 and  hca_stop==0 and thiol_stop==0 and amine_stop==0 and thiocarbamate_stop==0 and sulfonamide_stop==0 and imine_stop==0:
                                            chem_token=[];key_start=1
                                            chem_token=nltk.tokenize.word_tokenize(ald001) 
                                            ald002=len(chem_token)
                                            ald004=0;ald005=0;ald006=0
                                            aldehyde_stop=0;scan_rate=0
                                            for ald003 in range(0,ald002):
                                                if 'C=O' in str(chem_token[ald003]):
                                                    ald004=1   
                                                if  str(chem_token[ald003])=='(' and str(chem_token[ald003+1])=='O1' and str(chem_token[ald003+2])==')' :
                                                    ald005=1       
                                                if  str(chem_token[ald003])=='(' and str(chem_token[ald003+1])=='O2' and str(chem_token[ald003+2])==')' :
                                                    ald006=1   
                                                if ald004==1 or ald005==1  and ald006==0 and aldehyde_stop==0 and ald003==(ald002-1):
                                                    if (library_Am[9] in word_token_attribute[q]) or (library_Am[1] in word_token_attribute[q]):                           
                                                            save_arff.writelines( 'aldehyde'+',')
                                                            key_start=1;scan_rate=1
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ald001)) != ((len(processed_text_31))-1) and space_analyt==1:
                                                                save_arff.writelines('\n') 
                                                                key_start=1;
                                                            processed_text_32=(processed_text_31.index(ald001))
                                                            if analyt_soluble==1:
                                                                if processed_text_30[processed_text_32] < 0:
                                                                    save_arff.writelines( 'Water'+',')
                                                                    save_arff.writelines( 'KCl'+',')
                                                                    key_start=1
                                                                if processed_text_30[processed_text_32] >= 0:
                                                                    save_arff.writelines( 'THF'+',')   
                                                                    save_arff.writelines( 'LiClO4'+',')
                                                                    key_start=1
                                                                if  len(processed_text_31) > 1 and (processed_text_31.index(ald001)) != ((len(processed_text_31))-1):
                                                                    save_arff.writelines('\n') 
                                                                    key_start=1
                                                            
                                                            if  library_No[0] in word_token_attribute:
                                                                position_scan=[];
                                                                for s_r in range(0,len( processed_text_3)):
                                                                    if processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rate' or processed_text_3[s_r]=='scan' and processed_text_3[s_r+1]=='rates':
                                                                        for s_r2 in range(0,7):
                                                                            if (bool(re.search(r'\d', processed_text_3[s_r+s_r2])))==True and bool(re.search(r'[a-z]',processed_text_3[s_r+s_r2]))==False:
                                                                                if '-' in  str(processed_text_3[s_r+s_r2]) and (bool(re.search(r'\d', processed_text_3[s_r+s_r2]))):
                                                                                    processed_text_38_1=processed_text_3[s_r+s_r2].split('-')
                                                                                    if len(processed_text_38_1)==2 :
                                                                                        processed_text_3[s_r+s_r2]='';processed_text_3[s_r+s_r2]=processed_text_38_1[0]
                                                                                position_scan.append(s_r);key_start=1;posis_r=s_r2;s_r2=8
                                                                if position_scan==[]:
                                                                    position_scan=[0];
                                                                if 'V' in processed_text_3  or 'Vs' in processed_text_3 or 'Vs1'  in processed_text_3 or 'V/s' in processed_text_3  or 'v'  in processed_text_3:
                                                                    Volt='E0,';
                                                                elif 'mV' in processed_text_3 or 'mVs' in processed_text_3 or 'mVs1' in processed_text_3 or 'mV/s' in processed_text_3 :       
                                                                    Volt='E-3,';
                                                                else:
                                                                    Volt='E0,';
                                                                
                                                                save_arff.writelines(processed_text_3[(position_scan[-1])+posis_r]+Volt)
                                                            if  library_No[1] in word_token_attribute:
                                                                save_arff.writelines(str(min(processed_text_36))+',')
                                                            if name_article==1:
                                                                save_arff.writelines(str(pdf_n+','))
                                                            if  len(processed_text_31) > 1 and (processed_text_31.index(ph001)) != ((len(processed_text_31))-1):
                                                                save_arff.writelines('\n') 
                                                            key_start=1
                                                            aldehyde_stop=1       
                                                if  aldehyde_stop==0 and ald003==(ald002-1) and AM1==0 and (library_Am[9] in word_token_attribute[q]):    
                                                    save_arff.writelines( 'No'+','+'No'+','+'No'+',')  
                                                    key_start=1
                                                    if  library_No[0] in word_token_attribute:
                                                        if  scan_rate==0:
                                                            save_arff.writelines('0');
                                                    if  library_No[1] in word_token_attribute:
                                                        if  scan_rate==0:
                                                            save_arff.writelines('0,');
                                                    if name_article==1:
                                                        save_arff.writelines(str(pdf_n+','))
                                                    if  len(processed_text_31) > 1 and (processed_text_31.index(ald001)) != ((len(processed_text_31))-1):
                                                        save_arff.writelines('\n')
                                                    save_arff.writelines('\n');aldehyde_stop=1;                 
                                         
                                    if  phenol_stop==0 and  hca_stop==0 and thiol_stop==0 and amine_stop==0 and thiocarbamate_stop==0 and sulfonamide_stop==0 and imine_stop==0 and  aldehyde_stop==0 :
                                        save_arff.writelines( 'No'+','+'No'+','+'No'+',') 
                                        key_start=1
                                        if  library_No[0] in word_token_attribute:
                                            if  scan_rate==0:
                                                save_arff.writelines('0,');
                                        if  library_No[1] in word_token_attribute:
                                            if  scan_rate==0:
                                                save_arff.writelines('0,');
                                        if name_article==1:
                                            save_arff.writelines(str(pdf_n+','))
                                        save_arff.writelines('\n')
                                        
                                if word_token_attribute[q] in library_Cf:  
                                    if library_Cf[0] in word_token_attribute[q]:
                                        key_start=1 
                                    if library_Cf[0] in word_token_attribute[q]:
                                        key_start=1 
                                if word_token_attribute[q] in library_No:  
                                    if library_No[0] in word_token_attribute[q]:
                                        key_start=1 
                                    if library_No[1] in word_token_attribute[q]:
                                        key_start=1
                                         
                                if word_token_attribute[q] in library_Bk:
                                    if library_Bk[0] in word_token_attribute[q]:
                                        key_bk=0
                                        for bk_0 in processed_text_3:
                                            if key_bk==0:
                                                if bk_0 in solvent_organic:
                                                    solvorg=[]
                                                    solvorg.append(bk_0)
                                                    write_text="".join(solvorg)
                                                    save_arff.writelines(write_text+',')
                                                    key_bk=1 
                                                    key_start=1
                                            
                                        for bk_01 in processed_text_3:
                                            if key_bk==0:
                                                if bk_01 not in solvent_organic:    
                                                        save_arff.writelines(' No'+',')
                                                        key_bk=2  
                                                        key_start=1
                                    
                                    if library_Bk[1] in word_token_attribute[q]:
                                        key_bk1=0
                                        for bk1_0 in processed_text_3:
                                            if key_bk1==0:
                                                if bk1_0 in solvent_aqueous:
                                                    solaqu=[]
                                                    solaqu.append(bk1_0)
                                                    write_text="".join(solaqu)
                                                    save_arff.writelines(write_text+',')
                                                    key_bk1=1 
                                                    key_start=1
                                        for bk1_01 in processed_text_3:
                                            if key_bk1==0:
                                                if bk1_01 not in solvent_aqueous:    
                                                        save_arff.writelines(' No'+',')
                                                        key_bk1=2  
                                                        key_start=1  
                                
                                
                                if word_token_attribute[q] in library_Fm:
                                    add_Fm=0
                                    if library_Fm[0] in word_token_attribute[q]:
                                        key_Fm0=0
                                        for Fm_0 in processed_text_3:
                                            if key_Fm0==0:
                                                elecm0=[]
                                                if Fm_0 in electrode_modifield:
                                                    elecm0.append(Fm_0)
                                                    write_text="".join(elecm0)
                                                    if Fm_0 in mettalic_nanoparticle:
                                                        save_arff.writelines('mettalic_nanoparticle'+',')
                                                    if Fm_0 in graphene:
                                                        save_arff.writelines('graphene'+',')
                                                    if Fm_0 in polymer:
                                                        save_arff.writelines('polymer'+',')
                                                    if Fm_0 in fullerene:
                                                        save_arff.writelines('fullerene'+',')
                                                    if Fm_0 in nanotube:
                                                        save_arff.writelines('nanotubes'+',')
                                                    if Fm_0 in perovskite:
                                                        save_arff.writelines('perovskite'+',')    
                                                    if Fm_0 in oxide:
                                                        save_arff.writelines('oxide'+',')
                                                    if Fm_0 in QDs:
                                                        save_arff.writelines('QDs'+',')
                                                    if Fm_0 in zeolite:
                                                        save_arff.writelines('zeolite'+',')           
                                                    if Fm_0 in mettalic_complex_precussor:
                                                        save_arff.writelines('mettalic_complex_precussor'+',')           
                                                    if Fm_0 in ILs:
                                                        save_arff.writelines('ILs'+',')  
                                                    if Fm_0 in Enzyme:
                                                        save_arff.writelines('enzyme'+',')  
                                                    if Fm_0 in mettalic:
                                                        save_arff.writelines('mettalic'+',') 
                                                    if Fm_0 in MOF:
                                                        save_arff.writelines('MOF'+',') 
                                                    if Fm_0 in physical:
                                                        save_arff.writelines('phisical'+',') 
                                                    save_arff.writelines(write_text+',')
                                                    key_Fm0=1
                                                    key_start=1
                                                    add_Fm=add_Fm+1              
                                        for Fm_01 in processed_text_3:
                                            if key_Fm0==0:
                                                if Fm_01 not in electrode_modifield:    
                                                        save_arff.writelines(' No_or_other'+',')
                                                        save_arff.writelines(' No'+',')
                                                        key_Fm0=2  
                                                        key_start=1
                                                        add_Fm=add_Fm+1
                                        key_Fm1=0
                                        for Fm_1 in processed_text_3:
                                            if key_Fm0==1 and key_Fm1==0 and add_Fm !=3:
                                                elecm1=[]
                                                if Fm_1 in electrode_modifield and Fm_1 not in elecm0:
                                                    elecm1.append(Fm_1)
                                                    write_text="".join(elecm1)
                                                    if Fm_1 in mettalic_nanoparticle:
                                                        save_arff.writelines('mettalic_nanoparticle'+',')
                                                    if Fm_1 in graphene:
                                                        save_arff.writelines('graphene'+',')
                                                    if Fm_1 in polymer:
                                                        save_arff.writelines('polymer'+',')
                                                    if Fm_1 in fullerene:
                                                        save_arff.writelines('fullerene'+',')
                                                    if Fm_1 in nanotube:
                                                        save_arff.writelines('nanotubes'+',')
                                                    if Fm_1 in perovskite:
                                                        save_arff.writelines('perovskite'+',')    
                                                    if Fm_1 in oxide:
                                                        save_arff.writelines('oxide'+',')
                                                    if Fm_1 in QDs:
                                                        save_arff.writelines('QDs'+',')
                                                    if Fm_1 in zeolite:
                                                        save_arff.writelines('zeolite'+',') 
                                                    if Fm_1 in ILs:
                                                        save_arff.writelines('ILs'+',')  
                                                    if Fm_1 in Enzyme:
                                                        save_arff.writelines('enzyme'+',')  
                                                    if Fm_1 in mettalic:
                                                        save_arff.writelines('mettalic'+',') 
                                                    if Fm_1 in MOF:
                                                        save_arff.writelines('MOF'+',') 
                                                    if Fm_1 in physical:
                                                        save_arff.writelines('phisical'+',')     
                                                    if Fm_1 in mettalic_complex_precussor:
                                                        save_arff.writelines('mettalic_complex_precussor'+',')     
                                                    save_arff.writelines(write_text+',')
                                                    key_Fm1=1 
                                                    key_start=1
                                                    add_Fm=add_Fm+1
                                        for Fm_11 in processed_text_3:
                                            if key_Fm0==0 or key_Fm1==0 and add_Fm !=3:
                                                if Fm_11 not in electrode_modifield:
                                                        save_arff.writelines(' No_or_other'+',')
                                                        save_arff.writelines(' No'+',')    
                                                        key_Fm1=2  
                                                        key_start=1
                                                        add_Fm=add_Fm+1
                                        
                                        key_Fm2=0
                                        for Fm_2 in processed_text_3:
                                            if key_Fm1==1 and add_Fm !=3:
                                                elecm2=[]
                                                if Fm_2 in electrode_modifield and Fm_2 not in elecm0 and Fm_2 not in elecm1:
                                                    elecm2.append(Fm_2)
                                                    write_text="".join(elecm2)
                                                    if Fm_2 in mettalic_nanoparticle:
                                                        save_arff.writelines('mettalic_nanoparticle'+',')
                                                    if Fm_2 in graphene:
                                                        save_arff.writelines('graphene'+',')
                                                    if Fm_2 in polymer:
                                                        save_arff.writelines('polymer'+',')
                                                    if Fm_2 in fullerene:
                                                        save_arff.writelines('fullerene'+',')
                                                    if Fm_2 in nanotube:
                                                        save_arff.writelines('nanotubes'+',')
                                                    if Fm_2 in perovskite:
                                                        save_arff.writelines('perovskite'+',')    
                                                    if Fm_2 in oxide:
                                                        save_arff.writelines('oxide'+',')
                                                    if Fm_2 in QDs:
                                                        save_arff.writelines('QDs'+',')
                                                    if Fm_2 in zeolite:
                                                        save_arff.writelines('zeolite'+',') 
                                                    if Fm_2 in mettalic_complex_precussor:
                                                        save_arff.writelines('mettalic_complex_precussor'+',')   
                                                    if Fm_2 in ILs:
                                                        save_arff.writelines('ILs'+',')  
                                                    if Fm_2 in Enzyme:
                                                        save_arff.writelines('enzyme'+',')  
                                                    if Fm_2 in mettalic:
                                                        save_arff.writelines('mettalic'+',') 
                                                    if Fm_2 in MOF:
                                                        save_arff.writelines('MOF'+',') 
                                                    if Fm_2 in physical:
                                                        save_arff.writelines('phisical'+',')   
                                                    save_arff.writelines(write_text+',')
                                                    key_Fm2=1 
                                                    key_start=1
                                                    add_Fm=add_Fm+1
                                                    
                                        for Fm_21 in processed_text_3:
                                            if  key_Fm0==0 or key_Fm1==0 or key_Fm2==0 and add_Fm !=4:
                                                if Fm_21 not in electrode_modifield:
                                                        save_arff.writelines(' No_or_other'+',')
                                                        save_arff.writelines(' No'+',')    
                                                        key_Fm2=2  
                                                        key_start=1
                                                        add_Fm=add_Fm+1
                 
                                if word_token_attribute[q] in library_Md:
                                    if library_Md[0] in word_token_attribute[q]:
                                        key_r0=0
                                        for r_0 in processed_text_3:
                                            if key_r0==0:
                                                if r_0 in acid_inorganic:
                                                    aci=[]
                                                    aci.append(r_0)
                                                    write_text="".join(aci)
                                                    save_arff.writelines(write_text+',')
                                                    key_r0=1 
                                                    key_start=1
                                        for r_01 in processed_text_3:
                                            if key_r0==0:
                                                if r_01 not in acid_inorganic:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r0=2  
                                                        key_start=1

                                    if library_Md[1] in word_token_attribute[q]:
                                        key_r1=0
                                        for r_1 in processed_text_3:
                                            if key_r1==0:
                                                if r_1 in base_inorganic:
                                                    bai=[]
                                                    bai.append(r_1)
                                                    write_text="".join(bai)
                                                    save_arff.writelines(write_text+',')
                                                    key_r1=1 
                                                    key_start=1
                                        for r_11 in processed_text_3:
                                            if key_r1==0:
                                                if r_11 not in base_inorganic:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r1=2 
                                                        key_start=1
                                        
                                    if library_Md[2] in word_token_attribute[q]:
                                        key_r2=0
                                        for r_2 in processed_text_3:
                                            if key_r2==0:
                                                if r_2 in noble_gases:
                                                    ng=[]
                                                    ng.append(r_2)
                                                    write_text="".join(ng)
                                                    save_arff.writelines(write_text+',')
                                                    key_r2=1 
                                                    key_start=1
                                        for r_21 in processed_text_3:
                                            if key_r2==0:
                                                if r_21 not in noble_gases:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r2=2  
                                                        key_start=1  
                                                                      
                                    if library_Md[3] in word_token_attribute[q]:
                                        key_r3=0
                                        for r_3 in processed_text_3:
                                            if key_r3==0:
                                                if r_3 in gases:
                                                    gas=[]
                                                    gas.append(r_3)
                                                    write_text="".join(gas)
                                                    save_arff.writelines(write_text+',')
                                                    key_r3=1 
                                                    key_start=1
                                        for r_31 in processed_text_3:
                                            if key_r3==0:
                                                if r_31 not in gases:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r3=2 
                                                        key_start=1 
                                           
                                    if library_Md[4] in word_token_attribute[q]:
                                        key_r4=0
                                        for r_4 in processed_text_3:
                                            if key_r4==0:
                                                if r_4 in inert_gases:
                                                    ing=[]
                                                    ing.append(r_4)
                                                    write_text="".join(ing)
                                                    save_arff.writelines(write_text+',')
                                                    key_r4=1 
                                                    key_start=1
                                        for r_41 in processed_text_3:
                                            if key_r4==0:
                                                if r_41 not in inert_gases:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r4=2 
                                                        key_start=1   
                                        
                                    if library_Md[5] in word_token_attribute[q]:
                                        key_r5=0
                                        for r_5 in processed_text_3:
                                            if key_r5==0:
                                                if r_5 in water:
                                                    wat=[]
                                                    wat.append(r_5)
                                                    write_text="".join(wat)
                                                    save_arff.writelines(write_text+',')
                                                    key_r5=1 
                                                    key_start=1
                                        for r_51 in processed_text_3:
                                            if key_r5==0:
                                                if r_51 not in water:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r5=2 
                                                        key_start=1
                                                                    
                                    if library_Md[6] in word_token_attribute[q]:
                                        key_r6=0
                                        for r_6 in processed_text_3:
                                            if key_r6==0:
                                                if r_6 in peroxide:
                                                    per=[]
                                                    per.append(r_6)
                                                    write_text="".join(per)
                                                    save_arff.writelines(write_text+',')
                                                    key_r6=1
                                                    key_start=1
                                        for r_61 in processed_text_3:
                                            if key_r6==0:
                                                if r_61 not in peroxide:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r6=2 
                                                        key_start=1
                                                        
                                    if library_Md[7] in word_token_attribute[q]:
                                        key_r7=0
                                        for r_7 in processed_text_3:
                                            if key_r7==0:
                                                if r_7 in superoxide:
                                                    sper=[]
                                                    sper.append(r_7)
                                                    write_text="".join(sper)
                                                    save_arff.writelines(write_text+',')
                                                    key_r7=1 
                                                    key_start=1
                                        for r_71 in processed_text_3:
                                            if key_r7==0:
                                                if r_71 not in superoxide:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r7=2 
                                                        key_start=1                    
                                                        
                                    if library_Md[8] in word_token_attribute[q]:
                                        key_r8=0
                                        for r_8 in processed_text_3:
                                            if key_r8==0:
                                                if r_8 in alotropic_carbon:
                                                    alc=[]
                                                    alc.append(r_8)
                                                    write_text="".join(alc)
                                                    save_arff.writelines(write_text+',')
                                                    key_r8=1  
                                                    key_start=1
                                        for r_81 in processed_text_3:
                                            if key_r8==0:
                                                if r_81 not in alotropic_carbon:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r8=2 
                                                        key_start=1

                                    if library_Md[9] in word_token_attribute[q]:
                                        key_r9=0
                                        for r_9 in processed_text_3:
                                            if key_r9==0:
                                                if r_9 in mettalic_nanoparticle:
                                                    MNps=[]
                                                    MNps.append(r_9)
                                                    write_text="".join(MNps)
                                                    save_arff.writelines(write_text+',')
                                                    key_r9=1 
                                                    key_start=1
                                        for r_91 in processed_text_3:
                                            if key_r9==0:
                                                if r_91 not in mettalic_nanoparticle:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r9=2   
                                                        key_start=1
                                                          
                                    if library_Md[10] in word_token_attribute[q]:
                                        key_r10=0
                                        for r_10 in processed_text_3:
                                            if key_r10==0:
                                                if r_10 in materials_inorganic:
                                                    mati=[]
                                                    mati.append(r_10)
                                                    write_text="".join(mati)
                                                    save_arff.writelines(write_text+',')
                                                    key_r10=1 
                                                    key_start=1
                                        for r_101 in processed_text_3:
                                            if key_r10==0:
                                                if r_101 not in materials_inorganic:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r10=2  
                                                        key_start=1            
                                
 
                                    if library_Md[11] in word_token_attribute[q]:
                                        key_r11=0
                                        for r_11 in processed_text_3:
                                            if key_r11==0:
                                                if r_11 in acid_mine_ino:
                                                    ami=[]
                                                    ami.append(r_11)
                                                    write_text="".join(ami)
                                                    save_arff.writelines(write_text+',')
                                                    key_r11=1 
                                                    key_start=1
                                        for r_111 in processed_text_3:
                                            if key_r11==0:
                                                if r_111 not in materials_inorganic:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r11=2       
                                                        key_start=1     
                                                            
                                                            
                                    if library_Md[12] in word_token_attribute[q]:
                                        key_r12=0
                                        for r_12 in processed_text_3:
                                            if key_r12==0:
                                                if r_12 in silane:
                                                    sil=[]
                                                    sil.append(r_12)
                                                    write_text="".join(sil)
                                                    save_arff.writelines(write_text+',')
                                                    key_r12=1 
                                                    key_start=1
                                        for r_121 in processed_text_3:
                                            if key_r12==0:
                                                if r_121 not in silane:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r12=2  
                                                        key_start=1 
                                    
                                    if library_Md[13] in word_token_attribute[q]:
                                        key_r13=0
                                        for r_13 in processed_text_3:
                                            if key_r13==0:
                                                if r_13 in agent_oxidant:
                                                    oxid=[]
                                                    oxid.append(r_13)
                                                    write_text="".join(oxid)
                                                    save_arff.writelines(write_text+',')
                                                    key_r13=1 
                                                    key_start=1
                                        for r_131 in processed_text_3:
                                            if key_r13==0:
                                                if r_131 not in agent_oxidant:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r13=2  
                                                        key_start=1     
                                                        
                                    if library_Md[14] in word_token_attribute[q]:
                                        key_r14=0
                                        for r_14 in processed_text_3:
                                            if key_r14==0:
                                                if r_14 in agent_reducer:
                                                    redu=[]
                                                    redu.append(r_14)
                                                    write_text="".join(redu)
                                                    save_arff.writelines(write_text+',')
                                                    key_r14=1 
                                                    key_start=1
                                        for r_141 in processed_text_3:
                                            if key_r14==0:
                                                if r_141 not in agent_reducer:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r14=2  
                                                        key_start=1 
                                    
                                    if library_Md[15] in word_token_attribute[q]:
                                        key_r15=0
                                        for r_15 in processed_text_3:
                                            if key_r15==0:
                                                if r_15 in buffer:
                                                    buff=[]
                                                    buff.append(r_15)
                                                    write_text="".join(buff)
                                                    save_arff.writelines(write_text+',')
                                                    key_r15=1 
                                                    key_start=1
                                        for r_151 in processed_text_3:
                                            if key_r15==0:
                                                if r_151 not in buffer:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r15=2  
                                                        key_start=1
                                    
                                    if library_Md[16] in word_token_attribute[q]:
                                        key_r16=0
                                        for r_16 in processed_text_3:
                                            if key_r16==0:
                                                if r_16 in electrolyte:
                                                    elect=[]
                                                    elect.append(r_16)
                                                    write_text="".join(elect)
                                                    save_arff.writelines(write_text+',')
                                                    key_r16=1 
                                                    key_start=1
                                        for r_161 in processed_text_3:
                                            if key_r16==0:
                                                if r_161 not in electrolyte:    
                                                        save_arff.writelines(' No'+',')
                                                        key_r6=2  
                                                        key_start=1                                       
                                
                                else:
                                    key_2=0
                                    if key_2==0  and key_start==0:
                                        if word_token_attribute[q] in processed_text_3:
                                            write_text=' Yes,'
                                            save_arff.writelines(write_text)
                                            key_2=1

                                        if word_token_attribute[q] not in  processed_text_3:
                                            write_text=' No,'
                                            save_arff.writelines(write_text)
                                            key_2=2
                             
                                if word_token_attribute[-1]==word_token_attribute[q]:
                                    save_arff.writelines('\n')  
                                    word_token_all.clear()
                            save_arff.writelines('\n')
                                              
                                
#-------------------------------------------------------------------------------------------------------------------------------                   
                       
                              
                word_token_all.clear()
                attributes_open.close()
                file_pdf[x].close()  
                Menssager=messagebox.showinfo('HPSTC','Extraction text finished')
                print(Menssager)

#root.mainloop()
exe.mainloop()