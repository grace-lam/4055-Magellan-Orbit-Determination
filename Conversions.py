from visual import *
import numpy as np

def eulerRodriguez(phi,k,x):
    k = k*1./sqrt(k[0]**2+k[1]**2+k[2]**2)
    
    a = cos(phi/2.)
    b = k[0]*sin(phi/2.)
    c = k[1]*sin(phi/2.)
    d = k[2]*sin(phi/2.)

    rot = np.array([[a**2+b**2-c**2-d**2,2*(b*c-a*d),2*(b*d+a*c)],
                    [2*(b*c+a*d),a**2+c**2-b**2-d**2,2*(c*d-a*b)],
                    [2*(b*d-a*c),2*(c*d+a*b),a**2+d**2-b**2-c**2]])

    return np.dot(rot,x)

def toRad(degrees):
    return degrees*pi/180.

def toDeg(radians):
    return radians*180./pi

def toGaussian(jd):
    return jd*2*pi/365.2568983

def toAU(km):
    return km*6.68458712e-9

def toM(AU):
    return AU/6.68458712e-9*1000

def toKM(AU):
    return AU/6.68458712e-9

def eccentricAnomaly(e,M):
    E1 = 0
    while True:
        E0 = E1
        E1 += (M+e*sin(E1)-E1)/(1-e*cos(E1))
        if abs(E1-E0)<1e-10:
            break
    return E1

def toJD(year,month,day,hour,minute):
    if month==1 or month==2:
        year -= 1
        month += 12
    a = int(year/100.)
    b = int(a/4.)
    c = int(2-a+b)
    e = int(365.25*(year+4716))
    f = int(30.6001*(month+1))
    time = hour+minute/60.
    jd = c+day+e+f-1524.5+time/24.
    return jd

def getAngle(cosVal,sinVal):
    if cosVal >= 0:
        if sinVal >= 0:
            return acos(cosVal)
        else:
            return 2*pi - acos(cosVal)
    else:
        if sinVal >= 0:
            return acos(cosVal)
        else:
            return 2*pi - acos(cosVal)

def timeFromJD2000(year,month,day,hour,minute):
    if month==1 or month==2:
        year -= 1
        month += 12
    a = int(year/100.)
    b = int(a/4.)
    c = int(2-a+b)
    e = int(365.25*(year+4716))
    f = int(30.6001*(month+1))
    time = hour+minute/60.
    jd = c+day+e+f-1524.5+time/24.
    return jd-2451544.5

def earthSun(timeJD):
    G = 6.67e-11
    scalefactor = 50

    JD = timeFromJD2000(2015,6,21,21,58) #5650.69305555569 at 6/21/15 16:38 UTC=CT

    timeLabel = label(pos=(0,0,0),text="JD: " + str(JD+2451544.5),yoffset=-150)

    # variables for sun
    sun = sphere(pos=(0,0,0),color=color.yellow)
    sun.mass = 1.989e30 #mass of sun
    sun.p = vector(0,0,0)*sun.mass
    sun.radius = 695800.5e3*scalefactor #radius of sun
    sunLabel = label(pos=sun.pos,text="Sun")

    # variables for mercury
    mercury = sphere(pos=(4.278505337872649e10,-4.514486125017954e10,-7.614075244270340e9),color=color.red)
    mercury.mass = 4.8676e24 #mass of mercury
    mercury.p = vector(2.569053176558061e4,3.582894594784849e4,5.705914421209162e2)*mercury.mass
    mercury.radius = 6051.8e3*scalefactor #radius of mercury
    mercuryLabel = label(pos=mercury.pos,text="Mercury",xoffset=20,yoffset=40)

    # variables for venus
    venus = sphere(pos=(-6.196934157171567E+10,-8.888229875533229E+10,2.357927279994033E+09),color=color.orange)
    venus.mass = 0.3301e24 #mass of venus
    venus.p = vector(2.848320677199765E+04,-2.018984403111513E+04,-1.920458188769198E+03)*venus.mass
    venus.radius = 2439.7e3*scalefactor #radius of venus
    venusLabel = label(pos=venus.pos,text="Venus",xoffset=20,yoffset=40)

    # variables for earth
    semimajora = 149.60e9
    ecc = 0.016771123
    theta = (timeFromJD2000(2015,7,6,19,41)-timeFromJD2000(2015,6,21,21,58))/365.25*2*pi
    earth = sphere(pos=(0,-semimajora*(1-ecc**2)/(1-ecc*cos(theta)),0), color=color.green)
    earth.mass = 5.97e24 #mass of earth
    earth.p = vector(sqrt(G*sun.mass*(2./mag(earth.pos)-1./semimajora)),0,0)*earth.mass
    earth.radius = 6378.5e3*scalefactor #radius of earth
    earthLabel = label(pos=earth.pos,text="Earth",xoffset=20,yoffset=-30)

    # variables for moon
    moon = sphere(pos=(-9.078264635857833E+08,-1.518263560147420E+11,-1.551609456056356E+07),color=color.white)
    moon.mass = 0.07342e24
    moon.p = vector(2.879129441152029E+04,-1.036225983934763E+03,6.697590931749625E+01)*moon.mass
    moon.radius = 1738.1e3*scalefactor
    moonLabel = label(pos=moon.pos,text="Moon",yoffset=40)

    # variables for mars
    mars = sphere(pos=(1.387813253417010e10,2.331461012957870e11,4.544579412856832e9),color=color.red)
    mars.mass = 0.64174e24 #mass of mars
    mars.p = vector(-2.326621410943059e4,3.497163697777045e3,6.443278308533560e2)*mars.mass
    mars.radius = 3396.2e3*scalefactor #radius of mars
    marsLabel = label(pos=mars.pos,text="Mars",xoffset=20,yoffset=40)

    # variables for jupiter
    jupiter = sphere(pos=(-6.813325578777800E+11,4.257433806049283E+11,1.347768662908298E+10),color=color.orange)
    jupiter.mass = 1898.3e24 #mass of jupiter
    jupiter.p = vector(-7.086706064901995E+03,-1.047332330332310E+04,2.020503554452096E+02)*jupiter.mass
    jupiter.radius = 71492e3*scalefactor #radius of jupiter
    jupiterLabel = label(pos=jupiter.pos,text="Jupiter",xoffset=20,yoffset=40)

    # creating the trails for the minor bodies
    for minorbody in [mercury,venus,earth,moon,mars,jupiter]:
        minorbody.orbit=curve(pos=[minorbody.pos],color=minorbody.color,radius = 1e9)

    planets = [mercury,venus,earth,moon,mars,jupiter]
    planetLabels = [mercuryLabel,venusLabel,earthLabel,moonLabel,marsLabel,jupiterLabel]

    # creating the arrows for vectors
    vE = arrow(pos=earth.pos,axis=earth.p/earth.mass*5e5,shaftwidth=100,color=color.red)
    rES = arrow(pos=earth.pos,axis=sun.pos-earth.pos,shaftwidth=0.01,color=color.yellow)
    vELabel = label(pos=vE.pos+vE.axis/2,text="Earth Velocity",xoffset=20,yoffset=20,height=7,border=0.5)
    rESLabel = label(pos=rES.pos+rES.axis/2,text="Earth Sun Vector",xoffset=20,yoffset=20,height=7,border=0.5)

    changeTime = true
    userJD = timeJD-2451544.5

    dt = 100
    counter = 0

    # loop for orbits
    while true:
        JD += dt/86400.
        timeLabel.text = "JD: " + str(JD+2451544.5)

        if changeTime == true:
            rate(1e100)
        else:
            rate(1e20)
        
        if scene.mouse.clicked:
            scene.center = scene.mouse.getclick().pos

        for minorbody in [mercury,venus,earth,moon,mars,jupiter]:
            force = vector(0,0,0)
            for actor in [sun,mercury,venus,earth,moon,mars,jupiter]:
                if minorbody != actor:
                    dist1 = minorbody.pos - actor.pos
                    force += 6.67e-11 * minorbody.mass * actor.mass * dist1 / mag(dist1)**3
            if counter == 0:
                minorbody.new = minorbody.pos + dt*(minorbody.p/minorbody.mass) + 0.5/minorbody.mass*force*dt**2
                minorbody.old = minorbody.pos*1.
                minorbody.pos = minorbody.new*1.
                if minorbody == jupiter and actor == jupiter:
                    counter+=1
            else:
                minorbody.new = 2*minorbody.pos - minorbody.old - force/minorbody.mass*dt**2
                minorbody.old = minorbody.pos*1.
                minorbody.pos = minorbody.new*1.
                vE.pos = earth.pos
                vE.axis = 10000000*(earth.pos-earth.old)/(2.*dt)
                vELabel.pos = vE.pos+vE.axis/2

        for i in range(0,6):
            planetLabels[i].pos = planets[i].pos
            
        for a in [mercury,venus,earth,moon,mars,jupiter]:
            a.orbit.append(pos=a.pos)

        rES.pos=earth.pos
        rES.axis=sun.pos-earth.pos
        rESLabel.pos = rES.pos+rES.axis/2

        if changeTime == true and userJD <= JD:
            changeTime = false
            return rES.axis

def det(matrix,dim):
    if dim==2:
        a = matrix[0,0]
        b = matrix[0,1]
        c = matrix[1,0]
        d = matrix[1,1]
        return a*d-b*c
    else:
        deter = 0
        for a in range(len(matrix)):
            if a%2==0:
                deter += matrix[0,a]*det(np.delete(np.delete(matrix.copy(),a,1),0,0),len(matrix)-1)
            else:
                deter += -matrix[0,a]*det(np.delete(np.delete(matrix.copy(),a,1),0,0),len(matrix)-1)
        return deter

def cramers(mVar,mCon):
    d = det(mVar,len(mVar))*1.
    result = np.zeros(len(mVar))
    for a in range(len(mVar[0,:])):
        mCopy = mVar.copy()
        for i in range(len(mCopy)):
            mCopy[i,a] = mCon[i]
        result[a] = det(mCopy,len(mVar))*1./d
    return result

def toAltAndAz(haDeg,decDeg):
    lat = degToRad(40.003611)
    ha = degToRad(haDeg)
    dec = degToRad(decDeg)
    
    altRad = asin(sin(dec)*sin(lat)+cos(dec)*cos(lat)*cos(ha))
    cosA = (sin(dec)-sin(altRad)*sin(lat))/(cos(altRad)*cos(lat))
    aRad = acos(cosA)
    sinA = -1*cos(dec)*sin(ha)/cos(altRad)

    if sinA<0:
        aRad = 2*pi-aRad

    alt = str(radToDeg(altRad))
    az = str(radToDeg(aRad))
    
    return "altitude: " + alt + ", azimuth: " + az

def equatorialToEcliptic(raDeg,decDeg):
    e = degToRad(23.4)
    ra = degToRad(raDeg)
    dec = degToRad(decDeg)
    
    latRad = asin(-1*sin(ra)*cos(dec)*sin(e)+sin(dec)*cos(e))
    cosLong = cos(ra)*cos(dec)/cos(latRad)
    longRad = acos(cosLong)
    sinLong = (sin(ra)*cos(dec)*cos(e)+sin(dec)*sin(e))/cos(latRad)

    if sinLong<0:
        longRad = 2*pi-longRad

    longitude = radToDeg(longRad)
    latitude = radToDeg(latRad)
    while (longitude > 360):
        longitude -= 360
    while (longitude < 0):
        longitude += 360
    while (latitude > 90):
        latitude -= 360
    while (latitude < -90):
        latitude += 360
    
    return "longitude: " + str(longitude) + ", latitude: " + str(latitude)

def sphericalRad(a,c,B):
    cosb = cos(a)*cos(c)+sin(a)*sin(c)*cos(B)
    b = acos(cosb)
    sinA = sin(a)*sin(B)/sin(b)
    cosA = (cos(a)-cos(b)*cos(c))/(sin(b)*sin(c))
    A = getAngle(cosA,sinA)
    sinC = sin(c)*sin(B)/sin(b)
    cosC = (cos(c)-cos(a)*cos(b))/(sin(a)*sin(b))
    C = getAngle(cosC,sinC)
    print np.array([A,C,b])
    return np.array([A,C,b])

def sphericalDeg():
    aDeg = input("Enter the dec degree val of one side: ")
    cDeg = input("Enter the dec degree val of another side: ")
    BDeg = input("Enter the included angle of the spherical triangle: ")
    ans = sphericalRad(degToRad(aDeg),degToRad(cDeg),degToRad(BDeg))
    print "A: " + str(radToDeg(ans[0])) + "\nC: " + str(radToDeg(ans[1])) + "\nb: " + str(radToDeg(ans[2])) 
