from visual import *
import numpy as np
import Conversions as cv

def crossDot(a,b,c):
    return np.dot(np.cross(a,b),c)

def OD():
    inFile = open("JPL Test Data.txt","r")
    arrTime = np.zeros(3)
    arrRA = np.zeros(3)
    arrDec = np.zeros(3)
    arrayR = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
    arrR = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])

    ## NOTE: parallax correction done by using SBO-Sun vector (SBO: Sommers–Bausch Observatory)
    
    for a in range(3):
        data = inFile.readline().split()
        arrTime[a] = data[0]
        arrRA[a] = data[1]
        arrDec[a] = data[2]
        arrR[a,0] = float(data[3])
        arrR[a,1] = float(data[4])
        arrR[a,2] = float(data[5])

    t = np.zeros(3)
    ra = np.zeros(3)
    dec = np.zeros(3)
    R = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
    
    for i in range(3):
        t[i] = arrTime[i]
        ra[i] = cv.toRad(arrRA[i]*15)
        dec[i] = cv.toRad(arrDec[i])
        R[i] = arrR[i]

    # convert time to Gaussian days
    k = 0.017202099
    tau = np.zeros(3)

    tau[0] = k*(t[0]-t[1])
    tau[1] = k*(t[2]-t[0])
    tau[2] = k*(t[2]-t[1])

    print "tau:" 
    print tau[0],tau[1],tau[2]

    # first guess for aNeg and aPos
    a1 = abs(t[2]-t[1])/(t[2]-t[0])
    a3 = abs(t[0]-t[1])/(t[2]-t[0])

    print "a1,a3:"
    print a1,a3

    # calculate rhoHat from RA and Dec
    rhoHat = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
    for i in range(3):
        rhoHat[i] = np.array([cos(ra[i])*cos(dec[i]),sin(ra[i])*cos(dec[i]),sin(dec[i])])

    print "rhoHat:"
    print rhoHat[0],rhoHat[1],rhoHat[2]
    
    # initial guess for r and rDot
    rhoMag = np.zeros(3)
    rhoMag[0] = (a1*crossDot(R[0],rhoHat[1],rhoHat[2])-crossDot(R[1],rhoHat[1],rhoHat[2])+a3*crossDot(R[2],rhoHat[1],rhoHat[2]))/(a1*crossDot(rhoHat[0],rhoHat[1],rhoHat[2]))
    rhoMag[1] = (a1*crossDot(rhoHat[0],R[0],rhoHat[2])-crossDot(rhoHat[0],R[1],rhoHat[2])+a3*crossDot(rhoHat[0],R[2],rhoHat[2]))/(-crossDot(rhoHat[0],rhoHat[1],rhoHat[2]))
    rhoMag[2] = (a1*crossDot(rhoHat[1],R[0],rhoHat[0])-crossDot(rhoHat[1],R[1],rhoHat[0])+a3*crossDot(rhoHat[1],R[2],rhoHat[0]))/(a3*crossDot(rhoHat[1],rhoHat[2],rhoHat[0]))

    print "rhoMag:"
    print rhoMag[0],rhoMag[1],rhoMag[2]
    
    rhoVec = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
    rVec = np.array([np.zeros(3),np.zeros(3),np.zeros(3)])
    for i in range(3):
        rhoVec[i] = rhoMag[i]*rhoHat[i]
        rVec[i] = rhoVec[i]-R[i]

    # light time correction
    c = 173.1446
    tCorr = np.zeros(3)
    tCorr[0] = t[0] - mag(rhoVec[0])/c
    tCorr[1] = t[1] - mag(rhoVec[1])/c
    tCorr[2] = t[2] - mag(rhoVec[2])/c

    tau[0] = k*(tCorr[0]-tCorr[1])
    tau[1] = k*(tCorr[2]-tCorr[0])
    tau[2] = k*(tCorr[2]-tCorr[1])

    # calculate r0 and rDot0
    r0 = mag(rVec[1])
    rDotVec0 = (rVec[2]-rVec[0])/tau[1]

    print "rhoVec:"
    print rhoVec[0],rhoVec[1],rhoVec[2]
    print "rVec:"
    print rVec[0],rVec[1],rVec[2]
    print "rMag:"
    print r0
    print "rDot:"
    print rDotVec0
    print "rDotMag:"
    print mag(rDotVec0)

    # recalculate more precise values
    f1 = 1-1./(2*r0**3)*tau[0]**2
    f3 = 1-1./(2*r0**3)*tau[2]**2
    g1 = tau[0]-1./(6*r0**3)*tau[0]**3
    g3 = tau[2]-1./(6*r0**3)*tau[2]**3
    a1 = g3/(f1*g3-f3*g1)
    a3 = -g1/(f1*g3-f3*g1)
    print "f1,f3:"
    print f1,f3
    print "g1,g3:"
    print g1,g3
    print "a1,a3:"
    print a1,a3

    rhoMag[0] = (a1*crossDot(R[0],rhoHat[1],rhoHat[2])-crossDot(R[1],rhoHat[1],rhoHat[2])+a3*crossDot(R[2],rhoHat[1],rhoHat[2]))/(a1*crossDot(rhoHat[0],rhoHat[1],rhoHat[2]))
    rhoMag[1] = (a1*crossDot(rhoHat[0],R[0],rhoHat[2])-crossDot(rhoHat[0],R[1],rhoHat[2])+a3*crossDot(rhoHat[0],R[2],rhoHat[2]))/(-crossDot(rhoHat[0],rhoHat[1],rhoHat[2]))
    rhoMag[2] = (a1*crossDot(rhoHat[1],R[0],rhoHat[0])-crossDot(rhoHat[1],R[1],rhoHat[0])+a3*crossDot(rhoHat[1],R[2],rhoHat[0]))/(a3*crossDot(rhoHat[1],rhoHat[2],rhoHat[0]))
    for i in range(3):
        rhoVec[i] = rhoMag[i]*rhoHat[i]
        rVec[i] = rhoVec[i]-R[i]

    # light time correction
    c = 173.1446
    tCorr = np.zeros(3)
    tCorr[0] = t[0] - mag(rhoVec[0])/c
    tCorr[1] = t[1] - mag(rhoVec[1])/c
    tCorr[2] = t[2] - mag(rhoVec[2])/c

    tau[0] = k*(tCorr[0]-tCorr[1])
    tau[1] = k*(tCorr[2]-tCorr[0])
    tau[2] = k*(tCorr[2]-tCorr[1])

    r0 = mag(rVec[1])
    rDotVec0 = (rVec[2]-rVec[0])/tau[1]

    print "rhoMag:"
    print rhoMag

    # iterations
    while true:
        last_r0 = r0
        
        f1 = 1-1./(2*r0**3)*tau[0]**2+np.dot(rVec[1],rDotVec0)/(2*r0**5)*tau[0]**3+1./24*(3./r0**3*(np.dot(rDotVec0,rDotVec0)/r0**2-1./r0**3)-15*np.dot(rVec[1],rDotVec0)**2/r0**7+1./r0**6)*tau[0]**4
        f3 = 1-1./(2*r0**3)*tau[2]**2+np.dot(rVec[1],rDotVec0)/(2*r0**5)*tau[2]**3+1./24*(3./r0**3*(np.dot(rDotVec0,rDotVec0)/r0**2-1./r0**3)-15*np.dot(rVec[1],rDotVec0)**2/r0**7+1./r0**6)*tau[2]**4
        g1 = tau[0]-1./(6*r0**3)*tau[0]**3+np.dot(rVec[1],rDotVec0)/(4*r0**5)*tau[0]**4
        g3 = tau[2]-1./(6*r0**3)*tau[2]**3+np.dot(rVec[1],rDotVec0)/(4*r0**5)*tau[2]**4
        a1 = g3/(f1*g3-f3*g1)
        a3 = -g1/(f1*g3-f3*g1)

        # FIND r0 USING METHOD #1:
        rhoMag[0] = (a1*crossDot(R[0],rhoHat[1],rhoHat[2])-crossDot(R[1],rhoHat[1],rhoHat[2])+a3*crossDot(R[2],rhoHat[1],rhoHat[2]))/(a1*crossDot(rhoHat[0],rhoHat[1],rhoHat[2]))
        rhoMag[1] = (a1*crossDot(rhoHat[0],R[0],rhoHat[2])-crossDot(rhoHat[0],R[1],rhoHat[2])+a3*crossDot(rhoHat[0],R[2],rhoHat[2]))/(-crossDot(rhoHat[0],rhoHat[1],rhoHat[2]))
        rhoMag[2] = (a1*crossDot(rhoHat[1],R[0],rhoHat[0])-crossDot(rhoHat[1],R[1],rhoHat[0])+a3*crossDot(rhoHat[1],R[2],rhoHat[0]))/(a3*crossDot(rhoHat[1],rhoHat[2],rhoHat[0]))
        for i in range(3):
            rhoVec[i] = rhoMag[i]*rhoHat[i]
            rVec[i] = rhoVec[i]-R[i]
        r0 = mag(rVec[1])

        # light time correction
        tCorr[0] = t[0] - mag(rhoVec[0])/c
        tCorr[1] = t[1] - mag(rhoVec[1])/c
        tCorr[2] = t[2] - mag(rhoVec[2])/c

        tau[0] = k*(tCorr[0]-tCorr[1])
        tau[1] = k*(tCorr[2]-tCorr[0])
        tau[2] = k*(tCorr[2]-tCorr[1])
        
        rDotVec0 = f3/(g1*f3-g3*f1)*rVec[0]-f1/(g1*f3-g3*f1)*rVec[2]

        if abs(r0-last_r0) < 1e-9:
            break

    print "r (AU, equatorial) = " + str(rVec[1])
    print "rDot (AU, equatorial) = " + str(rDotVec0*k)
    print "light-corrected time = " + str(tCorr[1])

    # equatorial to ecliptic
    ep = 23.45027755*pi/180.
    rVec0 = np.array([rVec[1][0],rVec[1][1]*cos(ep)+rVec[1][2]*sin(ep),rVec[1][1]*-sin(ep)+rVec[1][2]*cos(ep)])
    rDVec0 = np.array([rDotVec0[0],rDotVec0[1]*cos(ep)+rDotVec0[2]*sin(ep),rDotVec0[1]*-sin(ep)+rDotVec0[2]*cos(ep)])
    print "r (AU, ecliptic) = " + str(rVec0)
    print "rDot (AU, ecliptic) = " + str(rDVec0*k*365.2505337327707)
    
    # final r and rDot vectors
    r = rVec0
    rDot = rDVec0

    # JPL orbital elements
    JPL = np.array([1.82032556817991,0.326561881175707,23.2570743997991,164.851008520451,154.367006465303,351.363189589307])
    
    # classical orbital elements
    h = np.cross(r,rDot)
    pos = mag(r)
    vel = mag(rDot)*2*pi

    # semimajor axis
    a = 1./(2./pos-vel**2/(4*pi**2))
    aErr = abs(a-JPL[0])/(0.5*(a+JPL[0]))*100
    print "a = " + str(a) + ", % diff = " + str(aErr)

    # eccentricity
    lPerM = mag(h)*2*pi
    e = sqrt(1-(lPerM**2)/(a*4*pi**2))
    eErr = abs(e-JPL[1])/(0.5*(e+JPL[1]))*100
    print "e = " + str(e) + ", % diff = " + str(eErr)

    # inclination angle
    i = acos(h[2]/mag(h))
    iErr = abs(cv.toDeg(i)-JPL[2])/(0.5*(cv.toDeg(i)+JPL[2]))*100
    print "i = " + str(cv.toDeg(i)) + ", % diff = " + str(iErr)

    # longitude of ascending node
    cosO = -h[1]/(mag(h)*sin(i))
    sinO = h[0]/(mag(h)*sin(i))
    o = cv.getAngle(cosO,sinO)
    oErr = abs(cv.toDeg(o)-JPL[3])/(0.5*(cv.toDeg(o)+JPL[3]))*100
    print "o = " + str(cv.toDeg(o)) + ", % diff = " + str(oErr)

    # argument of periapsis
    cosU = (r[0]*cosO+r[1]*sinO)/mag(r)
    sinU = r[2]/(mag(r)*sin(i))
    u = cv.getAngle(cosU,sinU)
    sinV = a*(1-e**2)/(e*mag(h))*(np.dot(r,rDot)/mag(r))
    cosV = 1/e*(a*(1-e**2)/mag(r)-1)
    v = cv.getAngle(cosV,sinV)
    w = (u-v)%(2*pi)
    wErr = abs(cv.toDeg(w)-JPL[4])/(0.5*(cv.toDeg(w)+JPL[4]))*100
    print "w = " + str(cv.toDeg(w)) + ", % diff = " + str(wErr)

    # eccentric anomaly
    cosE = (e+cos(v))/(1+e*cos(v))
    sinE = (sqrt(1-e**2)*sin(v))/(1+e*cos(v))
    E = cv.getAngle(cosE,sinE)
    
    # mean anomaly
    M = E-e*sin(E)
    mErr = abs(cv.toDeg(M)-JPL[5])/(0.5*(cv.toDeg(M)+JPL[5]))*100
    print "M = " + str(cv.toDeg(M)) + ", % diff = " + str(mErr)

OD()
