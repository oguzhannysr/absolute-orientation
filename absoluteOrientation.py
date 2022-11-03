import numpy as np

def AbsoluteOrientation(coorTxtPath):
    
    veriler=np.genfromtxt(coorTxtPath,dtype=str).astype(float)
    model_x1=veriler[:,0]
    model_y1=veriler[:,1]
    model_z1=veriler[:,2]
    cisim_x1=veriler[:,3]
    cisim_y1=veriler[:,4]
    cisim_z1=veriler[:,5]
    model_koord=np.stack((model_x1,model_y1,model_z1),axis=1)
    cisim_koord=np.stack((cisim_x1,cisim_y1,cisim_z1),axis=1)
    
    l=l=np.vstack(veriler[:,3:6]).reshape((veriler.shape[0]*3,1))
    birim_matris=np.identity(3)
    katsayılar_m=np.zeros((3*veriler.shape[0],12))
    
    for i in range(0,12,3):
        katsayılar_m[i:i+3,0:3]=birim_matris
    
    
    katsayılar_m[0:1,3:6]=veriler[0:1,0:3]
    katsayılar_m[3:4,3:6]=veriler[1:2,0:3]
    katsayılar_m[6:7,3:6]=veriler[2:3,0:3]
    katsayılar_m[9:10,3:6]=veriler[3:4,0:3]
    
    katsayılar_m[1:2,6:9]=veriler[0:1,0:3]
    katsayılar_m[4:5,6:9]=veriler[1:2,0:3]
    katsayılar_m[7:8,6:9]=veriler[2:3,0:3]
    katsayılar_m[10:11,6:9]=veriler[3:4,0:3]
    
    katsayılar_m[2:3,9:12]=veriler[0:1,0:3]
    katsayılar_m[5:6,9:12]=veriler[1:2,0:3]
    katsayılar_m[8:9,9:12]=veriler[2:3,0:3]
    katsayılar_m[11:12,9:12]=veriler[3:4,0:3]
    
    if  veriler.shape[0] == 4:
        x=np.dot(np.linalg.inv(katsayılar_m),l)
    
    elif veriler.shape[0] > 4: 
        N=np.dot(np.transpose(katsayılar_m),katsayılar_m)
        n=np.dot(np.transpose(katsayılar_m),l)
        x=np.dot(np.linalg.inv(N),n)
    
    x0=x[0,0]
    y0=x[1,0]
    z0=x[2,0]
    a11=x[3,0]
    a12=x[4,0]
    a13=x[5,0]
    a21=x[6,0]
    a22=x[7,0]
    a23=x[8,0]
    a31=x[9,0]
    a32=x[10,0]
    a33=x[11,0]
    
    #transformasyon_matrisi
    t=np.stack((a11,a12,a13,a21,a22,a23,a31,a32,a33)).reshape(3,3)
    lamda=np.sqrt((a11**2+a21**2+a31**2))
    #y. dönüklük matris
    d=((1/lamda)*t)
    
    #yenimodel
    
    cisim_yeni = []
    for i in range(0,veriler.shape[0]):
        sonuc=np.vstack((x0,y0,z0))+((np.dot(t,np.vstack((veriler[i:i+1,0],veriler[i:i+1,1],veriler[i:i+1,2])))))
        cisim_yeni.append(sonuc)
    
    v = []
    for i,j in zip(range(0,len(cisim_yeni)),range(0,len(cisim_yeni))):
        sonuc2=(np.vstack((cisim_x1[i],
                           cisim_y1[i],
                           cisim_z1[i]))-cisim_yeni[j])
        v.append(sonuc2)
    v=np.vstack((v))
    
    v = []
    for i in (range(0,veriler.shape[0])):
        bir = np.vstack((cisim_x1[i],
                           cisim_y1[i],
                           cisim_z1[i]))
        sonuc2 = bir - cisim_yeni[i]
        v.append(sonuc2)
    v=np.vstack((v))
    
    
    if  veriler.shape[0] == 4:
        qxx=np.linalg.inv(katsayılar_m)
        m0=np.sqrt(np.dot((np.transpose(v)),v)/(veriler.shape[0]+1))
    
    elif veriler.shape[0] > 4: 
        qxx=np.linalg.inv(N)
        m0=np.sqrt(np.dot((np.transpose(v)),v)/(3*veriler.shape[0]-7))
    
    
    #son cisim korrdinatı hesabı istenen nokta
    def cisimler(x,y,z):
        yeni_nokta_cisim=np.vstack((x0,y0,z0)) + np.dot(t,np.vstack((x,y,z)))
        return yeni_nokta_cisim

