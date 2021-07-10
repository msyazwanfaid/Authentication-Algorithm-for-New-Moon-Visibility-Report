#!/usr/bin/env python
# coding: utf-8

# In[4]:


import os,glob, pandas as pd
import datetime as dt
from pytz  import timezone
from skyfield import almanac
from skyfield.api import Topos,load
from skyfield.api import EarthSatellite
from skyfield import api
import json
from skyfield.api import GREGORIAN_START
import numpy as np
import math

planets = load('de422.bsp')
earth = planets['earth']
sun = planets['sun']
moon = planets['moon']
d=1


#skyfield variable
planets = load('de422.bsp')
earth = planets['earth']
sun = planets['sun']
moon = planets['moon']
d=1

ts = load.timescale()
eph = api.load('de422.bsp')
path = 'F:/OneDrive - Universiti Malaya/PhD/References/Source of Moon Sighting Report/Raw Data/Calculated Data/calculatedtest.csv'

df = pd.read_csv(path)

#Schaefer Contrast Threshold

from skyfield.api import load, Topos
from skyfield.trigonometry import position_angle_of

def Atmospheric_MoonBrightness(x):
    
    try:

        ts.julian_calendar_cutoff = GREGORIAN_START
        location = Topos(latitude_degrees=x['Lat'], longitude_degrees=x['Long'])
        t0 = ts.utc(x['Year'], x['Month'], x['Day']-x['TZ']/24)
        t1 = ts.utc(x['Year'], x['Month'], x['Day']+1-x['TZ']/24)


        f = almanac.sunrise_sunset(eph,location)
        t, y = almanac.find_discrete(t0, t1, f)

        if x['O'] == "M":
            z = True
        else:
            z = False
        for ti, yi in zip(t, y):
            if yi == z:
                #print(ti.utc_iso())
                #return('%s' % (ti.utc_strftime('%H:%M:%S')))
                ysunset = ti.utc.year
                mosunset = ti.utc.month
                dsunset = ti.utc.day
                hsunset = ti.utc.hour
                msunset = ti.utc.minute
                ssunset = ti.utc.second
            else:
                None

        boston = earth + Topos(latitude_degrees=x['Lat'],longitude_degrees=x['Long'], elevation_m=0)
        sun_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(sun)
        sun_app = sun_astro.apparent()
        sun_alt, sun_az, sun_distance = sun_app.altaz()
        sun_ra, sun_dec, sun_distance = sun_app.radec()


        moon_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(moon)
        moon_app = moon_astro.apparent()
        moon_alt, moon_az, sun_distance = moon_app.altaz()
        daz= abs(moon_az.degrees-sun_az.degrees)

        #print(moon_alt.degrees)


        while moon_alt.degrees > 0:
            boston = earth + Topos(latitude_degrees=x['Lat'],longitude_degrees=x['Long'], elevation_m=0)
            sun_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(sun)
            sun_app = sun_astro.apparent()
            sun_alt, az, distance = sun_app.altaz()

            moon_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(moon)
            moon_app = moon_astro.apparent()
            moon_alt, az, distance = moon_app.altaz()


            #Calculating Width
            moon_earth_distance = (moon.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))-earth.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))).distance().km
            sun_earth_distance = (sun.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))-earth.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))).distance().km
            horizontal_parallax = math.degrees(math.asin(6378.14/moon_earth_distance))
            moon_sun_distance_ratio=0.272481
            semidiameter_geocentric = math.degrees( math.asin( moon_sun_distance_ratio * math.sin(math.radians(horizontal_parallax)) ) )
            semidiameter_topocentric = semidiameter_geocentric * ( 1 + math.sin( math.radians(( moon_alt.degrees ) )) * math.sin( math.radians(( horizontal_parallax ))) )
            moonsun_relative_altitude = abs(moon_alt.degrees) + abs(sun_alt.degrees)
            width = semidiameter_topocentric * (1- math.cos( math.radians ( sun_app.separation_from(moon_app).degrees)))



            #Calculating Phase Angle
            arcl=sun_app.separation_from(moon_app).degrees
            v = planets['moon'].at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))) - planets['Sun'].at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))
            moon_sun_distance = v.distance().km
            sun_earth_distance = sun_distance.km
            phase_angle=180+math.degrees(math.atan(((moon_sun_distance*math.sin(math.radians(arcl)))/(sun_earth_distance-moon_sun_distance*math.cos(math.radians(arcl))))))


            #Calculating Illuminated Fraction
            k = (1+math.cos(math.radians(phase_angle)))/2

            #Calculating Moon Area
            moon_area = math.pi*semidiameter_topocentric**2

            #Calculating Illuminated Area
            ilum_area = k*moon_area

            #Extra Atmospheric Moon Brightness
            m= -12.73+0.026*phase_angle+4*10**-9*phase_angle**4
            out_moonbrightness_s10 = (1/ilum_area)*2.51**(10-m)
            out_moonbrightness_mag = -((math.log10(out_moonbrightness_s10))*2.5-26.33)

            #Light Loss from Atmosphere
            ##Air Mass Effect
            ###Raylieh
            moon_zenith = float(90-((moon_alt.degrees)))
            jure1 = float(math.cos(math.radians(moon_zenith)))
            jur1 = float(0.01*math.sqrt(8.2))
            jur2 = (-30*math.cos(math.radians(moon_zenith)))/(math.sqrt(8.2))
            jur3 = jure1+jur1*math.exp(jur2)
            jur = pow(jur3,-1)

            ###Aerosol
            jue1 = float(0.01*math.sqrt(1.5))
            jue2 = (-30*math.cos(math.radians(moon_zenith)))/(math.sqrt(1.5))
            jue3 = jure1+jue1*math.exp(jue2)
            jue = pow(jue3,-1)

            ###Ozone
            juo1 = (math.sin(math.radians(moon_zenith)))/(1+(20/6378.14))
            juo2 = 1-pow(juo1,2)
            juo = pow(juo2,-0.5)

            ##Atmosphere Effect
            ###Rayleigh
            ele = (float(10))
            kr = 0.1066 * math.exp(-(ele/1000)/8.2)

            ###Ozone
            lat = float(x['Lat'])
            matara = float((sun_ra.hours*15))
            ko1 = (lat*math.cos(math.radians(matara))-math.cos(math.radians(3*lat)))
            ko2 = 3.0 + 0.4*ko1
            ko = 0.031 *(ko2/3.0)

            #Aerosol
            ke1 = 0.12*pow(1,1.3)
            ke2 = math.exp(-ele/1.5)
            ke3 = math.exp(-ele/1.5)*(1-(0.32/math.log(0.5)))**(4/3)
            ke4 = 1+0.33 * math.sin(math.radians(sun_ra.hours*15))
            ke5= ke2*ke3
            ke = ke1 *ke3 * ke4

            #Inner Atmospheric Moon Brightness
            bulanmagfc = pow(10,((-16.8-out_moonbrightness_mag)/2.5))
            difmag = float(kr*jur + ke*jue + ko*juo)
            maghlg = bulanmagfc*(pow(10,(-difmag/2.5)))
            bulanmagket = float(-16.8-2.5*(math.log10(maghlg)))
            sunra = sun_ra.hours*15

            #Calculate Twilight Brighntess
            moon_zenith_distance = float(90-moon_alt.degrees)
            tw1 = -(7.5*(pow(10,-5))*moon_zenith_distance+5.05*pow(10,-3))*daz 
            tw2 = (3.67*pow(10,-4)*moon_zenith_distance-0.458)*-1*(sun_alt.degrees) 
            tw3 = (9.17*pow(10,-3))*moon_zenith_distance+3.525 
            tw = tw1+tw2+tw3 
            twmag1 = pow(10,tw) 
            twmag2 = 9.17*pow(10,4)*twmag1 
            tw_no_lightpollution = -(math.log10(twmag2)*2.5-27.78)

            a=319.186511755673
            b=-46.9419919726453
            c=0.434369283318305
            d=2.39138122541257
            e=-0.00121311966265414
            f=-0.0488735671014705
            g=-0.0396262989070584
            h=-0.000000928758741258742
            i=0.0000548215089277455
            j=0.001346675592325110
            θ=moon_zenith_distance
            Zlp=x['LP']
            tw_lightpollution= a+b*Zlp+c*θ+d*pow(Zlp,2)+e*pow(θ,2)+f*Zlp*θ+g*pow(Zlp,3)+h*pow(θ,3)+i*Zlp*pow(θ,2)+j*pow(Zlp,2)*θ

            if(tw_no_lightpollution<=tw_lightpollution):
                twmag=tw_no_lightpollution
            else:
                twmag=tw_lightpollution


            #Calculate Contrast
            moon_brightness_nL = pow(10,(-0.4*(bulanmagket-26.33)))
            twilight_brightness_nL = pow(10,(-0.4*(twmag-26.33)))

            c= ((moon_brightness_nL- twilight_brightness_nL)/twilight_brightness_nL)
            msunset=msunset+1

            dia = (semidiameter_topocentric/60*2)/60
            #twlpno = pow(10,(-0.4*(twlp-26.33)))
            eyesns1 = (40/(20))
            eyesns2 = 8.28*pow(twilight_brightness_nL,-0.29)
            eyesns3 = pow(10,eyesns2)
            eyesns4 = (eyesns1 * eyesns3)
            eyesns = eyesns4/3600
            #bulanmagketno = pow(10,(-0.4*(bulanmagket-26.33)))
            cont = (abs(moon_brightness_nL-twilight_brightness_nL))/twilight_brightness_nL
            contsch1 = pow((eyesns/dia),2)
            contsch2 = 2.4*pow(twilight_brightness_nL,-0.1)
            contsch = 0.0028+contsch2*contsch1

            #print(twmag,bulanmagket)
            if cont > contsch:
                #print(twmag,bulanmagket)
                return "Moon Sighted"
            else :
                #print(twmag,bulanmagket)
                None


    except:
        return ''

df['Modified Schaefer Prediction Model'] = df.apply(Atmospheric_MoonBrightness, axis=1)

#Crumey Contrast Threshold

from skyfield.api import load, Topos
from skyfield.trigonometry import position_angle_of

def Atmospheric_MoonBrightness(x):
    
    try:

        ts.julian_calendar_cutoff = GREGORIAN_START
        location = Topos(latitude_degrees=x['Lat'], longitude_degrees=x['Long'])
        t0 = ts.utc(x['Year'], x['Month'], x['Day']-x['TZ']/24)
        t1 = ts.utc(x['Year'], x['Month'], x['Day']+1-x['TZ']/24)


        f = almanac.sunrise_sunset(eph,location)
        t, y = almanac.find_discrete(t0, t1, f)

        if x['O'] == "M":
            z = True
        else:
            z = False
        for ti, yi in zip(t, y):
            if yi == z:
                #print(ti.utc_iso())
                #return('%s' % (ti.utc_strftime('%H:%M:%S')))
                ysunset = ti.utc.year
                mosunset = ti.utc.month
                dsunset = ti.utc.day
                hsunset = ti.utc.hour
                msunset = ti.utc.minute
                ssunset = ti.utc.second
            else:
                None

        boston = earth + Topos(latitude_degrees=x['Lat'],longitude_degrees=x['Long'], elevation_m=0)
        sun_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(sun)
        sun_app = sun_astro.apparent()
        sun_alt, sun_az, sun_distance = sun_app.altaz()
        sun_ra, sun_dec, sun_distance = sun_app.radec()


        moon_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(moon)
        moon_app = moon_astro.apparent()
        moon_alt, moon_az, sun_distance = moon_app.altaz()
        daz= abs(moon_az.degrees-sun_az.degrees)

        #print(moon_alt.degrees)


        while moon_alt.degrees > 0:
            boston = earth + Topos(latitude_degrees=x['Lat'],longitude_degrees=x['Long'], elevation_m=0)
            sun_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(sun)
            sun_app = sun_astro.apparent()
            sun_alt, az, distance = sun_app.altaz()

            moon_astro = boston.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))).observe(moon)
            moon_app = moon_astro.apparent()
            moon_alt, az, distance = moon_app.altaz()


            #Calculating Width
            moon_earth_distance = (moon.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))-earth.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))).distance().km
            sun_earth_distance = (sun.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))-earth.at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))).distance().km
            horizontal_parallax = math.degrees(math.asin(6378.14/moon_earth_distance))
            moon_sun_distance_ratio=0.272481
            semidiameter_geocentric = math.degrees( math.asin( moon_sun_distance_ratio * math.sin(math.radians(horizontal_parallax)) ) )
            semidiameter_topocentric = semidiameter_geocentric * ( 1 + math.sin( math.radians(( moon_alt.degrees ) )) * math.sin( math.radians(( horizontal_parallax ))) )
            moonsun_relative_altitude = abs(moon_alt.degrees) + abs(sun_alt.degrees)
            width = semidiameter_topocentric * (1- math.cos( math.radians ( sun_app.separation_from(moon_app).degrees)))



            #Calculating Phase Angle
            arcl=sun_app.separation_from(moon_app).degrees
            v = planets['moon'].at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset))) - planets['Sun'].at(ts.utc((ysunset), (mosunset), (dsunset),(hsunset),(msunset)))
            moon_sun_distance = v.distance().km
            sun_earth_distance = sun_distance.km
            phase_angle=180+math.degrees(math.atan(((moon_sun_distance*math.sin(math.radians(arcl)))/(sun_earth_distance-moon_sun_distance*math.cos(math.radians(arcl))))))


            #Calculating Illuminated Fraction
            k = (1+math.cos(math.radians(phase_angle)))/2

            #Calculating Moon Area
            moon_area = math.pi*semidiameter_topocentric**2

            #Calculating Illuminated Area
            ilum_area = k*moon_area

            #Extra Atmospheric Moon Brightness
            m= -12.73+0.026*phase_angle+4*10**-9*phase_angle**4
            out_moonbrightness_s10 = (1/ilum_area)*2.51**(10-m)
            out_moonbrightness_mag = -((math.log10(out_moonbrightness_s10))*2.5-26.33)

            #Light Loss from Atmosphere
            ##Air Mass Effect
            ###Raylieh
            moon_zenith = float(90-((moon_alt.degrees)))
            jure1 = float(math.cos(math.radians(moon_zenith)))
            jur1 = float(0.01*math.sqrt(8.2))
            jur2 = (-30*math.cos(math.radians(moon_zenith)))/(math.sqrt(8.2))
            jur3 = jure1+jur1*math.exp(jur2)
            jur = pow(jur3,-1)

            ###Aerosol
            jue1 = float(0.01*math.sqrt(1.5))
            jue2 = (-30*math.cos(math.radians(moon_zenith)))/(math.sqrt(1.5))
            jue3 = jure1+jue1*math.exp(jue2)
            jue = pow(jue3,-1)

            ###Ozone
            juo1 = (math.sin(math.radians(moon_zenith)))/(1+(20/6378.14))
            juo2 = 1-pow(juo1,2)
            juo = pow(juo2,-0.5)

            ##Atmosphere Effect
            ###Rayleigh
            ele = (float(10))
            kr = 0.1066 * math.exp(-(ele/1000)/8.2)

            ###Ozone
            lat = float(x['Lat'])
            matara = float((sun_ra.hours*15))
            ko1 = (lat*math.cos(math.radians(matara))-math.cos(math.radians(3*lat)))
            ko2 = 3.0 + 0.4*ko1
            ko = 0.031 *(ko2/3.0)

            #Aerosol
            ke1 = 0.12*pow(1,1.3)
            ke2 = math.exp(-ele/1.5)
            ke3 = math.exp(-ele/1.5)*(1-(0.32/math.log(0.5)))**(4/3)
            ke4 = 1+0.33 * math.sin(math.radians(sun_ra.hours*15))
            ke5= ke2*ke3
            ke = ke1 *ke3 * ke4

            #Inner Atmospheric Moon Brightness
            bulanmagfc = pow(10,((-16.8-out_moonbrightness_mag)/2.5))
            difmag = float(kr*jur + ke*jue + ko*juo)
            maghlg = bulanmagfc*(pow(10,(-difmag/2.5)))
            bulanmagket = float(-16.8-2.5*(math.log10(maghlg)))
            sunra = sun_ra.hours*15

            #Calculate Twilight Brighntess
            moon_zenith_distance = float(90-moon_alt.degrees)
            tw1 = -(7.5*(pow(10,-5))*moon_zenith_distance+5.05*pow(10,-3))*daz 
            tw2 = (3.67*pow(10,-4)*moon_zenith_distance-0.458)*-1*(sun_alt.degrees) 
            tw3 = (9.17*pow(10,-3))*moon_zenith_distance+3.525 
            tw = tw1+tw2+tw3 
            twmag1 = pow(10,tw) 
            twmag2 = 9.17*pow(10,4)*twmag1
            
            tw_no_lightpollution = -(math.log10(twmag2)*2.5-27.78)

            a=319.186511755673
            b=-46.9419919726453
            c=0.434369283318305
            d=2.39138122541257
            e=-0.00121311966265414
            f=-0.0488735671014705
            g=-0.0396262989070584
            h=-0.000000928758741258742
            i=0.0000548215089277455
            j=0.001346675592325110
            θ=moon_zenith_distance
            Zlp=x['LP']
            tw_lightpollution= a+b*Zlp+c*θ+d*pow(Zlp,2)+e*pow(θ,2)+f*Zlp*θ+g*pow(Zlp,3)+h*pow(θ,3)+i*Zlp*pow(θ,2)+j*pow(Zlp,2)*θ

            if(tw_no_lightpollution<=tw_lightpollution):
                twmag=tw_no_lightpollution
            else:
                twmag=tw_lightpollution

            #Calculate Contrast
            moon_brightness_nL = pow(10,(-0.4*(bulanmagket-26.33)))
            twilight_brightness_nL = pow(10,(-0.4*(twmag-26.33)))

            c= ((moon_brightness_nL- twilight_brightness_nL)/twilight_brightness_nL)

            twilight_brightness_cd=pow(10,((12.58-twmag)/2.5))
            r1=6.505*pow(10,-4)
            r2=-8.461*pow(10,-4)
            k1=7.633*pow(10,-3)
            k2=-7.174*pow(10,-3)


            contrast_crumey = pow(pow((pow((r1*pow(twilight_brightness_cd,(-1/4))+r2),2)/semidiameter_topocentric),3/5)+
                                       pow(k1*pow(twilight_brightness_cd,(-1/4))+k2,3/5),5/3)
            msunset=msunset+1



            #print(semidiameter_topocentric*60,c,contrast_crumey,twmag,twilight_brightness_cd)
            try:
                if c > contrast_crumey:
                    return "Moon Sighted"
                else :
                    None
            except:
                pass

    except:
        return ''

df['Crumey Prediction Model'] = df.apply(Atmospheric_MoonBrightness, axis=1)

df=df.round(decimals=2)


df.to_csv( "F:/OneDrive - Universiti Malaya/PhD/References/Source of Moon Sighting Report/Raw Data/Calculated Data/calculatedcontrast.csv", index=False, encoding='utf-8-sig')


# In[6]:


df.to_csv( "F:/OneDrive - Universiti Malaya/PhD/References/Source of Moon Sighting Report/Raw Data/Calculated Data/calculatedcontrast.csv", index=False, encoding='utf-8-sig')


# In[ ]:





# In[7]:


try:
    import winsound
    winsound.Beep(400, 1000)
except RuntimeError:
    print("The system is not able to beep the speaker")
except ImportError:
    print("Can't import winsound module")


# In[ ]:





# In[ ]:





# In[ ]:




