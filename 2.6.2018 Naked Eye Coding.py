import ephem
import math
import time
print(" New Moon Visibility Authentication")
lokasi = ephem.Observer()

print('Enter Longitude')
lokasi.lon = input()

print('Enter Latitude')
lokasi.lat = input()

print("Masukkan Jam")
t = float(input())

print("Masukkan Hari")
d = int(input())

print("Masukkan Bulan")
m = int(input())

print("Masukkan Tahun")
y = int(input ())

print("Masukkan Time Zone")
tz = int(input())
tz1 = t - tz

print("Altitude lokasi cerapan (km)")
lokasi.elevation = float(input())

print("Populasi Manusia di Lokasi Cerapan")
popu = int(input())

print("Jarak Lokasi Cerapan dari Pusat Bandar")
dis = int(input())

print("Anak bulan menghadap ke Bandar (yes/no)")
ab = input()

print("Snellen Ration (x/20)")
slnrto = int(input())

#Pengiraan Time Zone#
if tz1 < 0:
    dtz = d - 1
    ttz = (24 + t) - tz

elif tz1 <= 24:
    ttz = t-tz
    dtz=d
else :
    dtz = d+1
    ttz = ((t-tz)-24)
#print('')
#print('Waktu UTC mengikut timezone %i adalah %i dan hari %i'% (tz,ttz,dtz))

#Pengiraan Tahun Lompat
x = y%4
if x == 0:
    k = 1
    dt = 366
else:
    k = 2
    dt = 365

#Pengiraan Hari dalam Setahun dan Nisbah Hari dalam setahun#
n =  round(((275*m)/9)-0.5) - k * round(((m+9)/12)-0.5) + dtz -30 - 1 + ttz/24
y1 = str(y+ (n/dt))

#print('Nilai y adalah %s'%(y1))

d = ephem.Date(y1)
#print(ephem.Date(d))

lokasi.date=d

bulan =  ephem.Moon()
mata = ephem.Sun()

bulan.compute(lokasi)
mata.compute(lokasi)
bulanalt = float(math.degrees(bulan.alt))
mataalt = float(math.degrees(mata.alt))
d = str(lokasi.date)


table_data = [['Masa Cerapan','Altitud Matahari','Altitud Bulan','Beza Azimuth','Elongaasi', 'Kecerahan Langit Senja','Magnitud Bulan','Kontras Cahaya', 'Sempadan Kontras Cahaya', 'Kenampakkan']]
for row in table_data:
    print("{: >10} {: >22} {: >10} {: >10} {: >10} {: >20} {: >15} {: >17} {: >24} {: >15}  ".format(*row))
while mataalt >= -20:


#Main Geometrical Component
    bulan.compute(lokasi)
    mata.compute(lokasi)
    bulanalt = float(math.degrees(bulan.alt))
    mataaz = float(math.degrees(mata.az))
    bulanaz = float(math.degrees(bulan.az))
    bulanelon = float(math.degrees(bulan.elong))
    difaz = abs(bulanaz - mataaz)



#semi diameter geocentric
    sd = float((bulan.size/2)/60)
    belong = float(math.degrees(bulan.elong))
    mataalt= float(math.degrees(mata.alt))
    w= float((sd *(1-math.cos(math.radians(abs(mataalt-bulanalt)))))*60)

#Kecerahan Permukaan Bulan
    bulanmag = float(bulan.mag)#kene uji ketepatan equation ni.#
    bulanluas = float(math.pi*(pow(sd/60,2)))
    bulanphase = float(bulan.moon_phase)
    bulanluassinar = bulanluas*bulanphase
    bulansuf = (pow(2.51,(10-bulanmag)))/bulanluassinar
    bulansufmag = -((math.log10(bulansuf))*2.5-26.33)

#Kecerahan Permukaan Bulan
    bulanmag = float(bulan.mag)#kene uji ketepatan equation ni.#
    bulanluas = float(math.pi*(pow(sd/60,2)))
    bulanphase = float(bulan.moon_phase)
    bulanluassinar = bulanluas*bulanphase
    bulansuf = (pow(2.51,(10-bulanmag)))/bulanluassinar
    bulansufmag = -((math.log10(bulansuf))*2.5-26.33)
    #Rayleigh
    bulanzen = float(90-(abs(bulanalt)))
    jure1 = float(math.cos(math.radians(bulanzen)))
    jur1 = float(0.01*math.sqrt(8.2))
    jur2 = (-30*math.cos(math.radians(bulanzen)))/(math.sqrt(8.2))
    jur = pow(jure1+jur1*math.exp(jur2),-1)
    #Aerosol
    jue1 = float(0.01*math.sqrt(1.5))
    jue2 = (-30*math.cos(math.radians(bulanzen)))/(math.sqrt(1.5))
    jue = pow(jure1+jue1*math.exp(jue2),-1)
    #Ozone
    juo1 = (math.sin(math.radians(bulanzen)))/(1+(20/6378.14))
    juo2 = 1-pow(juo1,2)
    juo = pow(juo2,-0.5)
    #Rayleigh
    ele = float(lokasi.elevation)
    kr = 0.1066 * math.exp(-ele/8.2)
    #Ozone
    lat = float(math.degrees(lokasi.lat))
    matara = float(math.degrees(mata.ra))
    ko1 = (lat*math.cos(math.radians(matara))-math.cos(math.radians(3*lat)))
    ko2 = 3.0 + 0.4*ko1
    ko = 0.031 *(ko2/3.0)
    #Aerosol
    ke1 = 0.12*pow(1,1.3)
    ke2 = math.exp(-ele/15)
    ke3 = pow((1-(0.32/math.log(0.5))),4/3)
    ke4 = 1+0.33 * math.sin(math.radians(matara))
    ke = ke1 * ke2 *ke3 * ke4
    #Kecerahan Bulan dalam Atmosfera
    bulanmagfc = pow(10,((-16.8-bulansufmag)/2.5))
    difmag = float(kr*jur + ke*jue + ko*juo)
    maghlg = bulanmagfc*(pow(10,(-difmag/2.5)))
    bulanmagket = float(-16.8-2.5*(math.log10(maghlg)))

#Pencemaran Cahaya
    B = 50
    U = 2.59
    V = 0.008
    h = 2.4
    k = 0.026
    lp1 = B*math.sqrt(popu)*(U*pow(pow(dis,2)+pow(h,2),-1))
    lp2 = V*pow((pow(dis,2)+pow(h,2)),-0.5)
    lp3 = math.exp(-k*pow((pow(dis,2)+pow(h,2)),0.5))
    lp = lp1 + lp2 + lp3
    maglp = -(math.log10(lp)*2.5-27.78)

#Kecerahan Langit Senja yang Dipengaruhi Pencemaran Cahaya
    #Kecerahan Langit Senja Kastner
    tw1 = -(7.5*(pow(10,-5))*bulanzen+5.05*pow(10,-3))*difaz
    tw2 = (3.67*pow(10,-4)*bulanzen-0.458)*-1*(mataalt)
    tw3 = (9.17*pow(10,-3))*bulanzen+3.525
    tw = tw1+tw2+tw3
    twmag1 = pow(10,tw)
    twmag2 = 9.17*pow(10,4)*twmag1
    twmag = -(math.log10(twmag2)*2.5-27.78)
    #lPencemaran Cahaya Stopage Point
    if(ab == 'yes'):
      ab = 1
    else:
     ab = -1
    bulanzenlp = bulanzen*ab
    alptw = 319.186511755673
    blptw = -46.9419919726453
    clptw = 0.434369283318305
    dlptw = 2.39138122541257
    elptw = -0.00121311966265414
    flptw = -0.0488735671014705
    glptw = -0.0396262989070584
    hlptw = -9.28758741258742E-07
    ilptw = 0.0000548215089277455
    jlptw = 0.00134667559232511
    lptw = alptw + blptw*maglp+ clptw*bulanzenlp + dlptw*pow(maglp,2)+elptw*pow(bulanzenlp,2)+flptw*maglp*bulanzenlp+glptw*pow(maglp,3)+hlptw*pow(bulanzenlp,3)+ilptw*maglp*pow(bulanzenlp,2)+jlptw*pow(maglp,2)*bulanzenlp
    #Kecerahan Langit Senja yang dicemari Pencemaran Cahaya
    if (lptw<twmag):
        twlp = lptw
    else:
     twlp = twmag

#Sempadan Kontrast Cahaya mata
    dia = (sd/60*2)/3600
    twlpno = pow(10,(-0.4*(twlp-26.33)))
    eyesns1 = (40/slnrto)
    eyesns2 = 8.28*pow(twlpno,-0.29)
    eyesns3 = pow(10,eyesns2)
    eyesns4 = (eyesns1 * eyesns3)
    eyesns = eyesns4/3600
    bulanmagketno = pow(10,(-0.4*(bulanmagket-26.33)))
    cont = (abs(bulanmagketno-twlpno))/twlpno
    contsch1 = pow((eyesns/dia),2)
    contsch2 = 2.4*pow(twlpno,-0.1)
    contsch = 0.0028+contsch2*contsch1
    if cont>contsch:
        kenam = 'Ya'
    else:
        kenam = 'Tidak'



#untuk kegunaan table
    bulanalt = '%0.2f' %bulanalt
    mataalt = '%0.2f' % mataalt
    difaz = '%0.2f' % difaz
    bulanelon = '%0.2f' % bulanelon
    w = '%0.2f' % w
    bulanphase = ('%0.2f') % bulan.phase
    bulanmag = '%0.2f' % bulanmag
    difmag = '%0.2f'% difmag
    bulanmagket = '%0.2f'% bulanmagket
    maglp = '%0.2f'% maglp
    twlp ='%0.2f'% twlp
    lokasi.date = lokasi.date + ephem.hour * tz
    cont = '%0.2f'%cont
    contsch = '%0.2f' %contsch
    d = str(lokasi.date)

    table_data = [[d,mataalt,bulanalt,difaz,bulanelon,twlp,bulanmagket,cont,contsch,kenam]]
    for row in table_data:
        print("{: >10} {: >10} {: >15} {: >13} {: >12} {: >16} {: >17} {: >22} {: >19} {: >20}".format(*row))
    lokasi.date = lokasi.date - ephem.hour * tz
    lokasi.date = lokasi.date + ephem.minute * 1
    bulanalt = float(math.degrees(bulan.alt))
    mataalt = float(math.degrees(mata.alt))

else:
    print("Bulan Sudah Terbenam")

time.sleep(6000)