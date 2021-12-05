
#>>> =====================================================================================#
#>>>                            IMAGING SCRIPT                                            #
#>>> =====================================================================================#
#>>>
#>>> Helpful tip: Use the commands %cpaste or %paste to copy
#>>> and paste indented sections of code into the casa command line.

#>>> The commands below serve as a guide to best practices for imaging
#>>> ALMA data. It does not replace careful thought on your part while
#>>> imaging the data. Before imaging, you should use the commands in
#>>> the first section of this script to prep the data for imaging.
#>>> The commands in both sections should be able to be run as
#>>> standard Python script. However, the cleaning in this script is
#>>> done interactively making the final product somewhat dependent on
#>>> the individual doing the clean. The final data products
#>>> are the primary beam corrected images (*.pbcor), and the primary
#>>> beams (*.pb). These images are converted to fits at the end
#>>> of the script. 

########################################
# Check CASA version

import re

try:
    import casalith
except:
    print("Script requires CASA 6.0 or greater")


if casalith.compare_version("<",[6,1,1,15]):
    print("Please use CASA version greater than or equal to 6.1.1-15 with this script")


##################################################
# Create an Averaged Continuum MS

#>>> Continuum images can be sped up considerably by averaging the data
#>>> together to reduce overall volume. Since the sensitivity of a
#>>> continuum image depends on its bandwidth, continuum images are
#>>> typically made by including as much bandwidth as possible in the
#>>> data while excluding any line emission. The following plotms command
#>>> pages through the spectral windows in a project allowing you to
#>>> identify channel ranges within spectral windows that do not include
#>>> *strong* line emission. You will form a continuum image by averaging
#>>> the line-free spws and/or channel ranges within spws. In most cases,
#>>> you will not need to create an image to select line channels,
#>>> although this is a possible path for exploration for cases where there is
#>>> wide-spread line emission.
#>>>
#>>> In general, it is not necessary to include narrow spectral windows 
#>>> (<250MHz) in the continuum image.

finalvis='calibrated_final.ms' # This is your output ms from the data
                               # preparation script.

# Use plotms to identify line and continuum spectral windows.
#>>> If you have a project with multiple fields, you will want to run
#>>> the following plotms command separately for each source. If the
#>>> spectra for each field are significantly different from each other,
#>>> it may be necessary to make separate average continuum  and
#>>> continuum-subtracted measurement sets for each field.
plotms(vis=finalvis, xaxis='channel', yaxis='amplitude',
       ydatacolumn='data',
       avgtime='1e8', avgscan=True, avgchannel='1', 
       iteraxis='spw' )


#>>> Note that when you average channels in plotms, it displays
#>>> the "bin" number rather than the average channel number of each
#>>> bin. 

#>>> If you don't see any obvious lines in the above plot, you may to try
#>>> to set avgbaseline=True with uvrange (e.g., <100m). Limiting the
#>>> uvrange to the short baselines greatly improves the visibility of
#>>> lines with extended emission.

#>>> If you have multiple sources, the line channel ranges may be
#>>> different for different sources. Thus you would need to repeat the
#>>> process below for each source.

# Set spws to be used to form continuum
contspws = '25,27,29,31'

# If you have complex line emission and no dedicated continuum
# windows, you will need to flag the line channels prior to averaging.
flagmanager(vis=finalvis,mode='save',
            versionname='before_cont_flags')

initweights(vis=finalvis,wtmode='weight',dowtsp=True)

# Flag the "line channels"
flagchannels='25:2~5;25~36;61~62;66~67;84~85;88~88;90~92;94~95;98~101;106~106;136~140;142~143;150~154;160~160;164~164;169~174;180~185;189~194;213~217;231~238;283~286;292~295;306~308;311~311;334~340;380~380;391~395;452~460;493~497;535~537;556~559;569~572;580~582;590~592;630~639;702~706;713~716;763~769;779~781;834~834;838~839;848~850;858~861;872~878;891~891;915~915;919~922;934~934;939~948;955~959;962~967;978~981;986~989;1016~1018;1123~1127;1139~1147;1155~1165;1173~1180;1185~1187;1197~1198;1218~1222;1227~1229;1240~1247;1254~1260;1285~1291;1307~1307;1313~1318;1336~1340;1346~1347;1363~1363;1368~1370;1380~1380;1449~1454;1457~1461;1496~1499;1504~1504;1513~1522;1539~1546;1553~1555;1562~1562;1573~1576;1600~1606;1619~1631;1636~1644;1662~1662;1666~1666;1672~1672;1683~1683;1689~1692;1702~1703;1706~1706;1735~1737;1746~1747;1776~1781;1795~1798;1802~1804;1839~1839;1849~1855;1863~1863;1919~1919,'\
        '27:26~31;92~92;102~105;152~154;158~162;172~176;183~188;194~195;199~200;206~206;226~229;234~234;240~240;257~257;279~280;284~289;294~295;312~315;328~328;334~337;369~371;395~395;408~413;416~416;424~424;428~434;445~446;469~472;477~480;487~490;509~511;519~522;540~542;555~561;611~614;626~629;639~643;711~713;721~726;735~739;742~744;754~758;768~770;786~788;800~804;856~859;868~878;912~914;945~948;1032~1035;1052~1057;1151~1151;1167~1168;1194~1199;1217~1217;1243~1250;1266~1273;1345~1348;1392~1396;1408~1416;1446~1450;1453~1456;1466~1467;1504~1508;1516~1518;1538~1539;1593~1596;1667~1689;1697~1701;1783~1788;1813~1817;1870~1874;1880~1883;1919~1919,'\
        '29:34~38;52~52;60~61;68~69;84~94;131~134;152~152;158~160;167~172;189~192;205~205;243~249;256~256;268~299;304~312;343~348;352~361;368~370;376~382;392~392;406~406;416~416;423~431;434~434;445~447;452~456;468~474;477~477;491~492;512~514;519~519;525~537;547~547;561~561;565~565;570~574;591~594;598~602;607~607;610~613;625~631;633~633;647~649;652~653;664~666;678~678;703~703;709~709;712~714;727~729;735~735;743~743;747~750;811~811;860~862;871~878;885~885;892~892;899~899;917~918;930~933;935~935;959~961;1019~1019;1057~1057;1070~1073;1105~1105;1121~1127;1141~1142;1158~1158;1172~1173;1179~1182;1185~1185;1211~1211;1220~1224;1229~1230;1247~1247;1249~1250;1259~1259;1273~1276;1301~1306;1319~1324;1336~1345;1372~1375;1416~1416;1442~1443;1461~1464;1471~1479;1503~1503;1513~1516;1537~1541;1582~1589;1606~1607;1624~1624;1628~1628;1637~1642;1669~1689;1709~1712;1721~1722;1727~1727;1736~1738;1758~1761;1777~1777;1783~1783;1809~1809;1814~1823;1825~1826;1831~1832;1834~1842;1844~1846;1849~1849;1858~1858;1861~1861;1888~1888;1905~1905;1913~1913;1916~1916,'\
        '31:7~10;16~20;35~35;38~38;44~44;46~47;49~49;53~55;60~62;73~77;80~84;102~106;130~132;148~148;150~150;158~165;172~172;178~181;184~186;191~193;228~228;245~248;253~253;290~290;306~310;347~350;401~403;410~413;440~445;452~454;467~467;471~471;473~473;508~508;524~524;540~542;551~551;561~562;565~569;628~629;661~661;664~664;682~682;693~694;707~719;724~731;745~753;760~761;766~766;775~783;794~795;798~799;816~823;827~828;840~840;857~870;872~881;888~888;922~925;927~927;966~968;986~989;999~1001;1019~1019;1054~1082;1097~1106;1112~1113;1121~1123;1175~1177;1211~1213;1218~1221;1255~1255;1265~1268;1272~1283;1291~1297;1303~1308;1316~1319;1343~1348;1392~1399;1405~1406;1433~1440;1444~1446;1479~1482;1491~1493;1518~1538;1567~1572;1583~1585;1662~1662;1669~1671;1674~1686;1715~1718;1741~1751;1761~1763;1775~1775;1777~1778;1808~1808;1825~1828;1879~1884;1903~1906'

flagdata(vis=finalvis,mode='manual',
          spw=flagchannels,flagbackup=False)

# check that flags are as expected, NOTE must check reload on plotms
# gui if its still open.
plotms(vis=finalvis,yaxis='amp',xaxis='channel',
       avgchannel='1',avgtime='1e8',avgscan=True,iteraxis='spw') 

# Average the channels within spws
contvis='calibrated_final_cont.ms'
rmtables(contvis)
os.system('rm -rf ' + contvis + '.flagversions')

#>>> Note that to mitigate bandwidth smearing, please keep the width
#>>> of averaged channels less than 125MHz in Bands 3, 4, and 6, and 250MHz
#>>> in Band 7 for both TDM and FDM modes. For example, for a 2GHz TDM window
#>>> with 15.625 MHz channels, this means that the maximum width parameter
#>>> should be 8 channels for Bands 3, 4, and 6 and 16 channels for Band 7.
#>>> This is especially important for any long baseline data. These limits
#>>> have been designed to have minimize the reduction of the peak flux to
#>>> 95%. 
split(vis=finalvis,
     spw=contspws,      
     outputvis=contvis,
      width=[128,128,128,128], # number of channels to average together. The final channel width should be less than 125MHz in Bands 3, 4, and 6 and 250MHz in Band 7.
     datacolumn='data')


# Check the weights. You will need to change antenna and field to
# appropriate values
plotms(vis=contvis, yaxis='wtsp',xaxis='freq',spw='',antenna='DA42',field='3')

# If you flagged any line channels, restore the previous flags
flagmanager(vis=finalvis,mode='restore',
            versionname='before_cont_flags')

# Inspect continuum for any problems
plotms(vis=contvis,xaxis='uvwave',yaxis='amp',coloraxis='spw')

# #############################################
# Image Parameters

#>>> You're now ready to image. 

# source parameters
# ------------------

field='3' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL POTENTIALLY TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.
gridder='standard' # uncomment if single field 


# image parameters.
# ----------------

#>>> Generally, you want 5-8 cells (i.e., pixels) across the narrowest
#>>> part of the beam. You can estimate the beam size using the following
#>>> equation: 206265.0/(longest baseline in wavelengths).  To determine
#>>> the longest baseline, use plotms with xaxis='uvwave' and
#>>> yaxis='amp'. Divide the estimated beam size by five to eight to get
#>>> your cell size. It's better to error on the side of slightly too
#>>> many cells per beam than too few. Once you have made an image,
#>>> please re-assess the cell size based on the beam of the image. 

#>>> To determine the image size (i.e., the imsize parameter),  
#>>> an imsize equal to the size of the primary beam is
#>>> usually sufficient. The ALMA 12m primary beam in arcsec scales as
#>>> 6300 / nu[GHz] and the ALMA 7m primary beam in arcsec scales as
#>>> 10608 / nu[GHz], where nu[GHz] is the sky frequency. However, if
#>>> there is significant point source and/or extended emission beyond
#>>> the edges of your initial images, you should increase the imsize to
#>>> incorporate more emission.  Note that the imsize parameter is
#>>> in PIXELS, not arcsec, so you will need to divide the image size
#>>> in arcsec by the pixel size to determine a value for imsize.

cell='0.0042arcsec' # cell size for imaging.
imsize = [2000,2000] # size of image in pixels.

# velocity parameters
# -------------------

outframe='lsrk' # velocity reference frame. 
veltype='radio' # velocity type. 

# imaging control
# ----------------

# The cleaning below is done interactively, so niter and threshold can
# be controlled within clean. 

weighting = 'briggs'
robust=0.0
niter=1000
threshold = '0.0mJy'

#>>> Guidelines for setting robust:

#>>> Robust < 0.0 is not recommended for mosaics with poor-uv
#>>> coverage. Using values of robust less than or equal to 0.0 will
#>>> lead to major artifacts in the images including uneven noise
#>>> across the image.

#>>> If you are uv-tapering the data, you should set robust=2 (natural
#>>> weighting) to avoid upweighting points that are going to be
#>>> downweighted by uv-taper.

#############################################
# Imaging the Continuuum

# Set the ms and continuum image name.
contvis = 'calibrated_final_cont.ms'         

contimagename ='S255IR-SMA1_sci.spw25_27_29_31.mfs.I.manual' 

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=contvis)
#delmod(vis=contvis)

for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
    rmtables(contimagename+ext)

#>>> If you're going be be imaging with nterms>1, then you also need
#>>> to removed the *.tt0, and *.tt1 images in additional to those
#>>> listed above.

#>>> If the fractional bandwidth for the aggregate continuum is
#>>> greater than 10%, set deconvolver='mtmfs' to use multi-term,
#>>> multi-frequency synthesis. This algorithm takes into account the
#>>> spatial spectral index variations in an image.  Note that only
#>>> ALMA Band 3 and the lower end of Band 4 can have fractional
#>>> bandwidths of greater than 10% and only when both sidebands are
#>>> employed.

tclean(vis=contvis,
       imagename=contimagename,
       field=field,
       specmode='mfs',
       deconvolver='hogbom', 
       imsize = imsize, 
       cell= cell, 
       weighting = weighting,
       robust = robust,
       niter = niter, 
       threshold = threshold,
       interactive = True,
       gridder = gridder,
       pbcor = True,
       usepointing=False)
       
#>>> If interactively cleaning (interactive=True), then note number of
#>>> iterations at which you stop. This number will help
#>>> replicate images. Do not clean empty
#>>> images. Just click the red X to stop the interactive and note the
#>>> RMS.

#>>> Note RMS.
# ~2400 iterations
# RMS ~ 34 uJy
# Beam = 0.025 x 0.019 arcsec

# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an optional step.
#contmaskname = 'cont.mask'
##rmtables(contmaskname) # if you want to delete the old mask
#os.system('cp -ir ' + contimagename + '.mask ' + contmaskname)



########################################
# Continuum Subtraction for Line Imaging

#>>> If you have observations that include both line and strong (>3 sigma
#>>> per final line image channel) continuum emission, you need to
#>>> subtract the continuum from the line data. You should not continuum
#>>> subtract if the line of interest is in absorption.

#>>> You can use au.invertChannelRanges(flagchannels,vis=finalvis) to
#>>> get the fitspw below. You will need to insert any continuum spws
#>>> that weren't included in flagchannels. For example, if your continuum
#>>> spws are '0,1,2' and flagchannels='1:260~500', au.invertChannelRanges will return
#>>> '1:0~259,1:501~3839'. The fitspw parameter should be '0,1:0~259,1:501~3839,2'
#>>> Make sure to cut and paste the output in fitspw below since PIs don't have
#>>> analysisUtilities by default.

fitspw = '25:0~1;6~24;37~60;63~65;68~83;86~87;89~89;93~93;96~97;102~105;107~135;141~141;144~149;155~159;161~163;165~168;175~179;186~188;195~212;218~230;239~282;287~291;296~305;309~310;312~333;341~379;381~390;396~451;461~492;498~534;538~555;560~568;573~579;583~589;593~629;640~701;707~712;717~762;770~778;782~833;835~837;840~847;851~857;862~871;879~890;892~914;916~918;923~933;935~938;949~954;960~961;968~977;982~985;990~1015;1019~1122;1128~1138;1148~1154;1166~1172;1181~1184;1188~1196;1199~1217;1223~1226;1230~1239;1248~1253;1261~1284;1292~1306;1308~1312;1319~1335;1341~1345;1348~1362;1364~1367;1371~1379;1381~1448;1455~1456;1462~1495;1500~1503;1505~1512;1523~1538;1547~1552;1556~1561;1563~1572;1577~1599;1607~1618;1632~1635;1645~1661;1663~1665;1667~1671;1673~1682;1684~1688;1693~1701;1704~1705;1707~1734;1738~1745;1748~1775;1782~1794;1799~1801;1805~1838;1840~1848;1856~1862;1864~1918,'\
        '27:0~25;32~91;93~101;106~151;155~157;163~171;177~182;189~193;196~198;201~205;207~225;230~233;235~239;241~256;258~278;281~283;290~293;296~311;316~327;329~333;338~368;372~394;396~407;414~415;417~423;425~427;435~444;447~468;473~476;481~486;491~508;512~518;523~539;543~554;562~610;615~625;630~638;644~710;714~720;727~734;740~741;745~753;759~767;771~785;789~799;805~855;860~867;879~911;915~944;949~1031;1036~1051;1058~1150;1152~1166;1169~1193;1200~1216;1218~1242;1251~1265;1274~1344;1349~1391;1397~1407;1417~1445;1451~1452;1457~1465;1468~1503;1509~1515;1519~1537;1540~1592;1597~1666;1690~1696;1702~1782;1789~1812;1818~1869;1875~1879;1884~1918,'\
        '29:0~33;39~51;53~59;62~67;70~83;95~130;135~151;153~157;161~166;173~188;193~204;206~242;250~255;257~267;300~303;313~342;349~351;362~367;371~375;383~391;393~405;407~415;417~422;432~433;435~444;448~451;457~467;475~476;478~490;493~511;515~518;520~524;538~546;548~560;562~564;566~569;575~590;595~597;603~606;608~609;614~624;632~632;634~646;650~651;654~663;667~677;679~702;704~708;710~711;715~726;730~734;736~742;744~746;751~810;812~859;863~870;879~884;886~891;893~898;900~916;919~929;934~934;936~958;962~1018;1020~1056;1058~1069;1074~1104;1106~1120;1128~1140;1143~1157;1159~1171;1174~1178;1183~1184;1186~1210;1212~1219;1225~1228;1231~1246;1248~1248;1251~1258;1260~1272;1277~1300;1307~1318;1325~1335;1346~1371;1376~1415;1417~1441;1444~1460;1465~1470;1480~1502;1504~1512;1517~1536;1542~1581;1590~1605;1608~1623;1625~1627;1629~1636;1643~1668;1690~1708;1713~1720;1723~1726;1728~1735;1739~1757;1762~1776;1778~1782;1784~1808;1810~1813;1824~1824;1827~1830;1833~1833;1843~1843;1847~1848;1850~1857;1859~1860;1862~1887;1889~1904;1906~1912;1914~1915;1917~1919,'\
        '31:0~6;11~15;21~34;36~37;39~43;45~45;48~48;50~52;56~59;63~72;78~79;85~101;107~129;133~147;149~149;151~157;166~171;173~177;182~183;187~190;194~227;229~244;249~252;254~289;291~305;311~346;351~400;404~409;414~439;446~451;455~466;468~470;472~472;474~507;509~523;525~539;543~550;552~560;563~564;570~627;630~660;662~663;665~681;683~692;695~706;720~723;732~744;754~759;762~765;767~774;784~793;796~797;800~815;824~826;829~839;841~856;871~871;882~887;889~921;926~926;928~965;969~985;990~998;1002~1018;1020~1053;1083~1096;1107~1111;1114~1120;1124~1174;1178~1210;1214~1217;1222~1254;1256~1264;1269~1271;1284~1290;1298~1302;1309~1315;1320~1342;1349~1391;1400~1404;1407~1432;1441~1443;1447~1478;1483~1490;1494~1517;1539~1566;1573~1582;1586~1661;1663~1668;1672~1673;1687~1714;1719~1740;1752~1760;1764~1774;1776~1776;1779~1807;1809~1824;1829~1878;1885~1902;1907~1919' # *line-free* channels for fitting continuum
linespw = '25' # line spectral windows. You can subtract the continuum from multiple spectral line windows at once. The linespw must be in the fitspw unless combine='spw'

finalvis='calibrated_final.ms'

uvcontsub(vis=finalvis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # regions without lines.
          excludechans=False, # fit the regions in fitspw
          #combine='spw', #uncomment if there are no line-free channels in the line spectral window.  
          solint='int',
          fitorder=1,
          want_cont=False) # This value should not be changed.

#>>> Note that the continuum subtraction is done for each field in 
#>>> turn. However, if the fields have different line-free channels, you
#>>> will need to do the continuum subtraction separately for each field.

# NOTE: Imaging the continuum produced by uvcontsub with
# want_cont=True will lead to extremely poor continuum images because
# of bandwidth smearing effects. For imaging the continuum, you should
# always create a line-free continuum data set using the process
# outlined above.


##############################################
# Image line emission [REPEAT AS NECESSARY]

#>>> If you did an mstransform/cvel, use the same velocity parameters in
#>>> the clean that you did for the regridding. If you did not do an
#>>> mstransform and have multiple executions of a scheduling block,
#>>> select the spws with the same rest frequency using the spw parameter
#>>> (currently commented out below). DO NOT INCLUDE SPWS WITH DIFFERENT
#>>> REST FREQUENCIES IN THE SAME RUN OF CLEAN: THEY WILL SLOW DOWN
#>>> IMAGING CONSIDERABLY.

finalvis = 'calibrated_final.ms'
linevis = finalvis + '.contsub' # uncomment if continuum subtracted


restfreq='234.50GHz' # Typically the rest frequency of the line of
                        # interest. If the source has a significant
                        # redshift (z>0.2), use the observed sky
                        # frequency (nu_rest/(1+z)) instead of the
                        # rest frequency of the
                        # line.

spw='0' # uncomment and replace with appropriate spw 

lineimagename =  'S255IR-SMA1_sci.spw25.cube.I.manual'

start='' # start velocity. See science goals for appropriate value.
width='5km/s' # velocity width. See science goals.
nchan = -1  # number of channels. See science goals for appropriate value.

# If necessary, run the following commands to get rid of older clean
# data.

#clearcal(vis=linevis)
#delmod(vis=linevis)

for ext in ['.image','.mask','.model','.image.pbcor','.psf','.residual','.pb','.sumwt']:
    rmtables(lineimagename + ext)

tclean(vis=linevis,
       imagename=lineimagename, 
       field=field,
       spw=spw,
       specmode='cube', # comment this if observing an ephemeris source
       perchanweightdensity=False, # uncomment if you are running in CASA >= 5.5.0. 
       start=start,
       width=width,
       nchan=nchan, 
       outframe=outframe,
       veltype=veltype, 
       restfreq=restfreq, 
       niter=niter,  
       threshold=threshold, 
       interactive=True,
       cell=cell,
       imsize=imsize, 
       weighting=weighting,
       robust=robust,
       gridder=gridder,
       pbcor=True,
       restoringbeam='common',
       usepointing=False) 

#>>> If interactively cleaning (interactive=True), then note number of
#>>> iterations at which you stop for the PI. This number will help the
#>>> PI replicate the delivered images. Do not clean empty
#>>> images. Just click the red X to stop the interactive and note the
#>>> RMS.

# No iterations performed
# RMS ~ 0.7 mJy/beam
# Beam = 0.025 x 0.019 arcsec

# If you'd like to redo your clean, but don't want to make a new mask
# use the following commands to save your original mask. This is an
# optional step.
# linemaskname = 'line.mask'
## rmtables(linemaskname) # uncomment if you want to overwrite the mask.
# os.system('cp -ir ' + lineimagename + '.mask ' + linemaskname)

##############################################
# Export the images

import glob

myimages = glob.glob("*manual*.pbcor")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True)

myimages = glob.glob("*manual*.pb")
for image in myimages:
    exportfits(imagename=image, fitsimage=image+'.fits',overwrite=True) 

##############################################
# Create Diagnostic PNGs

os.system("rm -rf *.png")
mycontimages = glob.glob("*mfs*manual.image")
for cimage in mycontimages:
    mymax=imstat(cimage)['max'][0]
    mymin=-0.1*mymax
    outimage = cimage+'.png'
    os.system('rm -rf '+outimage)
    imview(raster={'file':cimage,'range':[mymin,mymax]},out=outimage)

mylineimages = glob.glob("*cube*manual.image")
for limage in mylineimages:
    mom8=limage+'.mom8'
    os.system("rm -rf "+mom8)
    immoments(limage,moments=[8],outfile=mom8)
    mymax=imstat(mom8)['max'][0]
    mymin=-0.1*mymax
    os.system("rm -rf "+mom8+".png")
    imview(raster={'file':mom8,'range':[mymin,mymax]},out=mom8+'.png')


##############################################
# Analysis

# For examples of how to get started analyzing your data, see
#     https://casaguides.nrao.edu/index.php/TWHydraBand7_Imaging_4.3
#     
