#!/usr/bin/env python
# coding: utf-8

# In[75]:


import numpy as np
import math


# In[76]:


#telescope parameters
mirrorsize = 2.4 #m maybe scale up to check number
readOut = 2.97 #e-/pixel
darkCurrent= 5 #e-/hr //estimate
quantumEfficiency = 0.8 #estimate
pixelScale= 0.05
gain = 1.57 #e-/ADU


# In[77]:



#user input
targetType = input("Enter target type: 1 for point source, 2 for extended source: ")
if (targetType == '2'):
    targetType = input("Enter galaxy type: Elliptical or Spiral: ")

System = input("Enter magnitude system: (Vega or AB): ")

targetMag = float(input("Enter target magnitude: "))
       
if(targetType == '1'):
    seeing = float(input("Enter seeing: "))
    FWHM = seeing
else:
    FWHM = float(input("Enter FWHM: "))

filt = input("Enter filter of choice: options: UBVRI: ")

moon = input("Enter moon phase: ") #for sky brightness

airmass = float(input("Enter air mass: "))

if(targetType!='1'):
    seeing = float(input("Enter star condition (seeing): "))

if(targetType =='1'):
    print("Aperture: (Suggested aperture : 1.6*seeing/2.35): ")
else:
    print("Aperture: (Suggested aperture = 1.5 ~ 2 * FWHM): ")
   
aperture = float(input("Enter aperture: "))

#way = input("find (1) S/N ratio or (2) find time")


# #Ns
# # sky brightness/ 0 days
# if(moon=='0'):
#     SkyR = 0
#
# if(moon=='3'):
#     if(filt == 'U'):
#         SkyR = 21.5-22
#     if (filt =='B'):
#         SkyR = 22.4-22.7
#     if (filt =='V'):
#         SkyR = 21.7-21.8
#     if(filt == 'R'):
#         20.8-20.9
#     if(filt == 'I'):
#         SkyR = 0
#        
# if(moon =='7'):
#     if(filt == 'U'):
#         SkyR = 19.9-22
#     if (filt =='B'):
#         SkyR = 21.6-22.7
#     if (filt =='V'):
#         SkyR = 21.4-21.8
#     if(filt == 'R'):
#         SkyR = 20.6-20.9
#     if(filt == 'I'):
#         SkyR = 19.7-19.9
#
# if(moon =='10'):
#     if(filt == 'U'):
#         SkyR = 18.5-22
#     if (filt =='B'):
#         SkyR = 20.7-22.7
#     if (filt =='V'):
#         SkyR = 20.7-21.8
#     if(filt == 'R'):
#         SkyR = 20.3-20.9
#     if(filt == 'I'):
#         SkyR = 19.5-19.9
#
# else:#14
#     if(filt == 'U'):
#         SkyR = 17-22
#     if (filt =='B'):
#         SkyR = 19.5-22.7
#     if (filt =='V'):
#         SkyR = 20-21.8
#     if(filt == 'R'):
#         SkyR = 19.9-20.9
#     if(filt == 'I'):
#         SkyR = 19.2-19.9
#
# SkyR = pow(10,SkyR/(-2.5))
#
#
# #different for 4m Mosaic
# #0.5,0.06,0.1,0.05
# if(filt == 'U'):
#     Sky = SkyR * 0.5
#     correction = 0.5
# if (filt =='B'):
#     Sky = SkyR*1.5
#     correction = 0.25
# if (filt =='V'):
#     Sky = SkyR * 4.2
#     correction = 0.15
# if(filt == 'R'):
#     Sky = SkyR *9.9
#     correction = 0.1
# if(filt == 'I'):
#     Sky = SkyR * 13
#     correction = 0.07
#

# In[78]:



'''step 1'''
#point source = Vega
#spiral = AB
#elliptical = AB
print("works")

#filters file
if(filt == 'U'):
    data2= np.genfromtxt("./U.txt",delimiter ="",names = ["A","t"])
if (filt =='B'):
    data2= np.genfromtxt("./B.txt",delimiter ="",names = ["A","t"])
if (filt =='V'):
    data2= np.genfromtxt("./V.txt",delimiter ="",names = ["A","t"])
if(filt == 'R'):
    data2= np.genfromtxt("./R.txt",delimiter ="",names = ["A","t"])
if(filt == 'I'):
    data2= np.genfromtxt("./I.txt",delimiter ="",names = ["A","t"])

wave2 = data2["A"]
t = data2["t"]


#magnitude system file (fzp)
if(System=="Vega"):
    #data
    data1 = np.genfromtxt("./Vega_SED.dat",delimiter ="",names = ["A","flux"])
    wave1 = data1["A"]
    f = data1["flux"]

    #calculation
    i = 0
    j = 0
    s = 0
    while(i<len(wave1)-1):
        while(j<len(wave2)-1):
            #go thru each wavelength in transmission list
            if(wave2[j]<=wave1[i]<wave2[j+1]):
                #wavelength = wave1[i]
                transmission = (t[j] + t[j+1])/2 #find average
                transmission = transmission/100 #gives percentage in decimal
                s = s+transmission*f[i]*(wave1[i+1]-wave1[i]) #delta wavelength
            j=j+1
        i=i+1
        j=0

   
else:
    #AB
    fAB = []
    fv= 3.63 * pow(10,-20) #erg/cm s Hz
    c = 2.997 * pow(10,18)#A/s
   
    i =0
    while(i<len(wave2)):
        fAB.append(fv*c/wave2[i]/wave2[i])
        i=i+1
    #set up fAB
       
    i = 0
    s = 0
    while(i<len(wave2)-1):
        #wavelength = wave3[i]
        transmission = t[i]/100 #gives percentage in decimal
        s = s+transmission*fAB[i]*(wave2[i+1]-wave2[i]) #delta wavelength
        i=i+1
    #s = f temp
   
fzp = s

#for f temp


if(targetType=="1"):
    #point source
    data1 = np.genfromtxt("./Vega_SED.dat",delimiter ="",names = ["A","flux"])
    wave4 = data1["A"]
    f = data1["flux"]

if (targetType =="Spiral"):
    data4 = np.genfromtxt("./Sc_SED.dat",delimiter ="",names = ["A","flux"])
    wave4 = data4["A"]
    f = data4["flux"]
    print(wave4[0])
    for i in range(len(wave4)):
        wave4[i] = wave4[i]*10
        f[i]=f[i]/10
       
    print(wave4[0])


   
if (targetType == "Elliptical"):
    data4 = np.genfromtxt("./E_SED.dat",delimiter ="",names = ["A","flux"])
    wave4 = data4["A"]
    f = data4["flux"]
    for i in range(len(wave4)):
        wave4[i] = wave4[i]*10
        f[i]=f[i]/10

#calculation for galaxies
i = 0
j = 0
s = 0
while(i<len(wave4)-1):
    while(j<len(wave2)-1):
        #go thru each wavelength in transmission list
        #print(wave2[j],wave4[i],wave2[j+1])
        if wave2[j]<=wave4[i] and wave4[i]<wave2[j+1]:
            #wavelength = wave1[i]
            #print("hello")
            transmission = (t[j] + t[j+1])/2 #find average
            transmission = transmission/100 #gives percentage in decimal
            s = s+transmission*f[i]*(wave4[i+1]-wave4[i]) #delta wavelength
        j=j+1
    i=i+1
    j=0

ftemp = s

print("ftemp",ftemp)
print("fzp",fzp)
#done with ftemp

ratio = ftemp/fzp
mtemp = -2.5 * math.log(ratio,10)
norm = (targetMag - mtemp )/-2.5
norm = pow(10,norm)
print("norm",norm)


# In[ ]:





# In[79]:



'''step 2'''
s = 0
h = 6.62607004*pow(10,-27)#erg * s
c = 2.997 * pow(10,18)#A/s
i = 0
j = 0

#wave1,f for given target type (ftemp)
if (targetType == "1"):
    data1 = np.genfromtxt("./Vega_SED.dat",delimiter ="",names = ["A","flux"])
    wave1 = data1["A"]
    f = data1["flux"]
   
if(targetType == "Elliptical"):
    data1 = np.genfromtxt("./E_SED.dat",delimiter ="",names = ["A","flux"])
    wave1 = data1["A"]
    f = data1["flux"]

if(targetType =="Spiral"):
    data1 = np.genfromtxt("./Sc_SED.dat",delimiter ="",names = ["A","flux"])
    wave1 = data1["A"]
    f = data1["flux"]

   
photonNum = 0

while(i<len(wave1)-1):
    while(j<len(wave2)-1):
        #go thru each wavelength in transmission list
        if(wave2[j]<=wave1[i]<wave2[j+1]):
            #wavelength = wave1[i]
            transmission = (t[j] + t[j+1])/2 #find average
            transmission = transmission/100 #gives percentage in decimal
            photonNum = photonNum+norm * transmission*f[i]*(wave1[i+1]-wave1[i])/(h*c/wave1[i]) #delta wavelength
        j=j+1
    i=i+1
    j=0
#photon_num is result
result = 0
result = photonNum
Ntemp = photonNum
print(result)


# In[80]:


result = 0
result = photonNum
Ntemp = photonNum
print(result)


# In[81]:


print("quantum efficiency",quantumEfficiency)


# In[82]:



'''step 3'''
#x = airmass
#k depends on filter

aExtinction = correction* airmass / (-2.5)
aExtinction = pow(10,aExtinction)

result = result * aExtinction
print("extinction",aExtinction)
print("result",result)

'''step 4'''
if(targetType ==1): #for point source
    sigma = seeing/2.36

else:
    sigma =  FWHM/2.36
   
e = 2.71828
exponent = - aperture*aperture/2/sigma/sigma
faperture = 1-pow(e,exponent)
result = result * faperture
print("result",result)

print("aperture fraction",faperture)
'''step 5'''
result = result * quantumEfficiency
print("result",result)

'''step 6'''
mirrorsize = 4
mirrorsize = mirrorsize*100/2 #in cm
area = mirrorsize*mirrorsize*3.14
result= result * area #now total photon # and not photon # per unit area


# In[83]:


Nstar= result
Nd = darkCurrent/3600
Nr = readOut
Ns = Sky

#n = n pix, where r = aperture, convert
n = 3.14 * (aperture/pixelScale)*(aperture/pixelScale)


# In[84]:




way = input("find (1) S/N ratio or (2) find time: ")
#consider t to N/s

if (way == "1"):
    time = float(input("Enter time for observation: "))
    NS = time * Nstar / (math.sqrt(time*Nstar + n*(Nd*time + Ns*time + pow(Nr,2))))
    print("Signal to noise ratio is:",NS)
   

if (way == "2"):
    NS = float(input("Enter S/N for observation: "))
    a = pow(Nstar,2) #star is per second
    b = -pow(NS,2)*(Nstar + n*Ns + n*Nd)
    c = -pow(NS,2)*n*pow(Nr,2)
    time = (math.sqrt(b*b - 4*a*c) - b)/2/a
   
    print("time is:",time)

