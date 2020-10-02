import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import glob
import sys
import os.path

def HueCalculation(r,g,b,delta,Cmax):
    if delta==0:
        return 0
    elif Cmax==r:
        return int(360+ 60*(((g-b)/delta)))%360
    elif Cmax==g:
        return int(120+60*(b-r)/delta)
    else :
        return int(240+60*(r-g)/delta)

def SaturationCalculation(delta,Cmax):
    return 0 if Cmax==0 else delta/Cmax 

def RgbToHSVPixel(array,x,y):
    r=array[x][y][0]/255
    g=array[x][y][1]/255
    b=array[x][y][2]/255
    Cmax=max(r,g,b)
    Cmin=min(r,g,b)
    delta=Cmax-Cmin
    hue=HueCalculation(r,g,b,delta,Cmax)
    saturation=SaturationCalculation(delta,Cmax)
    return hue,saturation,Cmax

def NormalizeHSV(H,S,V):
    return int(H/2),int(255*S),int(255*V) 

def RGBtoHSV(RGBarray):
    HSVarray=np.zeros(RGBarray.shape,dtype=int)
    for i in range(len(RGBarray)):
        for j in range(len(RGBarray[i])):
            H,S,V=RgbToHSVPixel(RGBarray,i,j)
            H,S,V=NormalizeHSV(H,S,V)
            HSVarray[i][j][0]=H
            HSVarray[i][j][1]=S
            HSVarray[i][j][2]=V
    return HSVarray

def RGBtoHSVchangeVal(RGBarray,optionalFlag,value,boundary):
    HSVarray=np.zeros(RGBarray.shape,dtype=int)
    for i in range(len(RGBarray)):
        for j in range(len(RGBarray[i])):
            H,S,V=RgbToHSVPixel(RGBarray,i,j)
            if value>=0:
                if optionalFlag==1:
                    if H<=boundary:
                        H+=value
                elif optionalFlag==2:
                    if S*100<=boundary:
                        S+=float(value)/100
                elif optionalFlag==3:
                    if V*100<=boundary:
                       V+=float(value)/100
            else :
                if optionalFlag==1:
                    if H>=boundary:
                        H+=value
                elif optionalFlag==2:
                    if S*100>=boundary:
                        S+=float(value)/100
                elif optionalFlag==3:
                    if V*100>=boundary:
                        V+=float(value)/100
            H,S,V=NormalizeHSV(H,S,V)
            HSVarray[i][j][0]=H
            HSVarray[i][j][1]=S
            HSVarray[i][j][2]=V
    return HSVarray
    
def HSVToRGBPixel(array,x,y):
    H=array[x][y][0]*2
    S=array[x][y][1]/255
    V=array[x][y][2]/255
    C=V*S
    X=C*(1-abs((H/60)%2-1))
    m=V-C
    if 0<=H and H<60 :
        return int((C+m)*255),int( (X+m)*255),int(m*255)
    elif 60<=H and H<120 :
        return int((m+X)*255),int((m+C)*255),int(m*255)
    elif 120<=H and H<180 :
        return int(m*255),int((m+C)*255),int((m+X)*255)
    elif 180<=H and H<240 :
        return int(m*255),int((m+X)*255),int((m+C)*255)
    elif 240<=H and H<300 :
        return int((m+X)*255),int(m*255),int((m+C)*255)
    elif 300<=H and H<360 :
        return int((m+C)*255),int(m*255),int((m+X)*255)
    return 

def HSVtoRGB(HSVarray):
    RGBarray=np.zeros(HSVarray.shape,dtype=int)
    for i in range(len(HSVarray)):
        for j in range(len(HSVarray[i])):
            R,G,B=HSVToRGBPixel(HSVarray,i,j)
            RGBarray[i][j][0]=R
            RGBarray[i][j][1]=G
            RGBarray[i][j][2]=B
    return RGBarray

def RGBtoXYZPixel(array,x,y):
    R=array[x][y][0]/255
    G=array[x][y][1]/255
    B=array[x][y][2]/255

    X=0.4124564*R + 0.3575761*G + 0.1804375*B 
    Y=0.2126729*R + 0.7151522*G + 0.0721750*B 
    Z=0.0193339*R + 0.1191920*G + 0.9503041*B 
    return X,Y,Z

def RGBtoXYZ(RGBarray):
    XYZarray=np.zeros(RGBarray.shape)
    for i in range(len(RGBarray)):
        for j in range(len(RGBarray[i])):
            X,Y,Z=RGBtoXYZPixel(RGBarray,i,j)
            XYZarray[i][j][0]=X
            XYZarray[i][j][1]=Y
            XYZarray[i][j][2]=Z
    return XYZarray

def XYZtoLABpixel(array,x,y):
    Xn=0.950456
    Zn=1.088754

    X=array[x][y][0]
    Y=array[x][y][1]
    Z=array[x][y][2]
    L=116*FunctionHelperXYZtoLAB(Y) -16
    a=500*(FunctionHelperXYZtoLAB(X/Xn)-FunctionHelperXYZtoLAB(Y))
    b=200*(FunctionHelperXYZtoLAB(Y)-FunctionHelperXYZtoLAB(Z/Zn))
    return L,a,b


def FunctionHelperXYZtoLAB(value):
    if value > pow(6/29,3):
        return pow(value,1/3)
    else:
        return 1/3* pow(29/6,2)*value +4/29

def XYZtoLAB(XYZarray):
    LABarray=np.zeros(XYZarray.shape,dtype=int)
    for i in range(len(XYZarray)):
        for j in range(len(XYZarray[i])):
            L,A,B=XYZtoLABpixel(XYZarray,i,j)
            LABarray[i][j][0]=int(L*255/100)
            LABarray[i][j][1]=int(A+128)
            LABarray[i][j][2]=int(B+128)
    return LABarray

def CalculateArrays(RGBarray):
    HSVarray= RGBtoHSV(RGBarray)
    XYZarray=RGBtoXYZ(RGBarray)
    LABarray= XYZtoLAB(XYZarray)
    return HSVarray,XYZarray,LABarray

def PlotGraph(array):
    plt.imshow(array)
    plt.show()

def SetChangeHSVvalue(min,max):
    option=int(input("choose\n1:increase\n2:decrease\n"))
    if option==1:
        boundary=int(input("enter boundary(action will use for values smaller than boundary number): ("
            +str(min)+"-"+str(max)+"):"))
        value=int(input("enter value for increasing("+str(min)+"-"+ str(max-boundary)+") :"))
        if value+boundary>max:
            print("incorrect range, it must be "+str(min)+"-"+str(max))
            return np.nan
        else :
            return value,boundary
    elif option==2:
        boundary=int(input("enter boundary(action will use for values more than boundary number): ("
            +str(min)+"-"+str(max)+"):"))
        value=int(input("enter value for decreasing("+str(min)+ "-"+str(boundary)+") :"))
        if boundary - value <min:
            print("incorrect range, it must be "+min+"-"+max)
            return
        else :
            return -value,boundary        
    else :
        print("incorrect option")
        return


def ChangeHSVvalue(RGBarray,optionFlag,value,boundary):
    if value is not None:
        HSVarray=RGBtoHSVchangeVal(RGBarray,optionFlag,value,boundary)
        RGBarray=HSVtoRGB(HSVarray)
        XYZarray=RGBtoXYZ(RGBarray)
        LABarray=XYZtoLAB(XYZarray)
        return RGBarray,HSVarray,XYZarray,LABarray
    else :
        print("incorrect parameters")
        return

def SetHSVvalueOptions():
    changedParam=int(input("1:Change hue\n"+
        "2:Change saturation\n"+
        "3:Change value\n"+
        "4:Cancel"))
    if changedParam == 4:
        return
    else:
        plt.close('all')
        if changedParam==1:
            value,boundary=SetChangeHSVvalue(0,360)
            return value,boundary,changedParam
        elif changedParam==2:
            value,boundary=SetChangeHSVvalue(0,100)
            return value,boundary,changedParam
        elif changedParam==3:
            value,boundary=SetChangeHSVvalue(0,100)
            return value,boundary,changedParam
        print("incorrect option")
        return

def ConvertToPureL(Larray):
    for i in range(len(Larray)):
        Larray[i]*=100/255 
    return Larray

def GetLComponentArray(LabArray):
    return ConvertToPureL(LabArray[:,:,0].flatten())

def PlotHistogram(Larray):
    plt.hist(Larray)
    plt.show()

def OptionalAction(RGBarray,HSVarray,XYZarray,LABarray):
    while 1:
        optionNumber=input("Choose next action: \n"+
        "1:Show RGB\n"+
        "2:Show HSV\n"+
        "3:Show LAB\n"+
        "4:Show LAB and HSV for concrete pixel\n"+
        "5:Change HSV value\n"+
        "6:Plot L component histogram\n"+
        "7:Exit\n")
        if int(optionNumber)==1:
            PlotGraph(RGBarray)
        elif int(optionNumber)==2:
            PlotGraph(HSVarray)
        elif int(optionNumber)==3:
            PlotGraph(LABarray)
        elif int(optionNumber)==4:
            x,y=map(int, 
                input("Enter concrete pixel(x,y) where \n"+
                    "0<=x<="+str(len(RGBarray))+
                    "\n0<=y<="+str(len(RGBarray[0]))+"\nformat: x y :").split())
            print("HSV: ",RgbToHSVPixel(RGBarray,x,y))
            print("LAB: ",XYZtoLABpixel(XYZarray,x,y))
        elif int(optionNumber)==5:
            value,boundary,changedParam = SetHSVvalueOptions()
            tmpRgb,tmpHSV,tmpXYZ,tmpLAB=ChangeHSVvalue(RGBarray,changedParam,value,boundary)
            if tmpRgb is not None:
                RGBarray,HSVarray,XYZarray,LABarray=tmpRgb,tmpHSV,tmpXYZ,tmpLAB
            else :
                print("incorrect change result")
        elif int(optionNumber)==6:
            PlotHistogram(GetLComponentArray(LABarray))
        elif int(optionNumber)==7:
            return
        else :
            print("incorrect action")

def CheckIsFileExist(filepath):
    return os.path.isfile(filepath)

while 1:
    option=int(input("1.Run\n2.Exit\n"))
    if option==1:
        imagePath=input("Please enter full image path\n")
        if CheckIsFileExist(imagePath) or (len(glob.glob(imagePath))>0 and 
                CheckIsFileExist(glob.glob(imagePath)[0])):
            RGBarray =mpimg.imread(imagePath)
            HSVarray,XYZarray,LABarray=CalculateArrays(RGBarray)
            OptionalAction(RGBarray,HSVarray,XYZarray,LABarray)
        else:
            print("incorrect filepath")
    elif option==2:
        sys.exit(0)