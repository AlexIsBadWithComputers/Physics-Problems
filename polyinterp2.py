import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import chisquare as tealatte
import math

def ListRipper(list_,*args): #Deletes list pieces by index.
	indexes = sorted(list(args), reverse = True)
	for index in indexes:
		del list_[index]
	return list_

def fileyoink(filename): #grab the data I want
	with open(str(filename),"r") as f: 
		search = f.readlines()
		f.close()
	data = [] 
	ReadBegin = True  #Modified as I've filtered data manually.
	for i, line in enumerate(search):
		if ReadBegin ==True:
			info = line.split()
			if "-99.99" in line:
				continue #No data available, do not record.
			else:
				ListRipper(info,0,1,4,5,6)#I onlt want two columns of data
				info = list(map(float,info))
				data.append(info)
			if "#             date" in line:
				ReadBegin = True
	return data
	f.close()


co2 = np.array(fileyoink("co2dat.txt")) #Data to fit to
co2points = np.array(fileyoink("co2dat2.txt")) #Extra data to plot... but not to fit to
co2points[:,0] = co2points[:,0]-1998
co2[:,0] = co2[:,0]-1998 #Normalize data but avoid negative values below datas


xp = np.linspace(0,10,1000)

#Set up data transformations

trans = np.array([np.log(co2[:,1]),1.0/co2[:,1],np.arcsin(co2[:,1]/max(co2points[:,1])),np.sqrt(co2[:,1])])
#print(trans[3])
trans2 = np.array([np.log(co2points[:,1]), 1.0/co2points[:,1],np.arcsin(co2points[:,1]/max(co2points[:,1])),np.sqrt(co2points[:,1])])

labs =["Ln(Co2)","1/Co2)","Arcsin(Co2/MaxVal(Co2))","Sqrt(Co2)"]

colormap = plt.cm.gist_ncar # I demand a pretty plot. No uggos here. 
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, 8)])

labels = []  #Just set some junk up. 

for j in range(0,len(trans)):
    plt.subplot2grid((5,3),(j,0),colspan=1)
    a,b =[min(trans2[j]),max(trans2[j])]  #set axes range to be reasonable
    print(a,b)
    plt.ylim((a,b))
    plt.scatter(co2points[:,0],trans2[j])
    if j == len(trans)-1:
        plt.xlabel("Year - 1998")
    
   
    plt.ylabel(labs[j])
    residlist = []
    orderlist = []
    for i in range(1,50,5):
        if i == 0:
            continue
        plt.xlim((0,10))
	    #plt.subplot2grid((5,3),(j,0),colspan=2)
        fit = np.polyfit(co2[:,0],trans[j],i)
        iterfit = np.poly1d(fit)
        plt.plot(xp,iterfit(xp))
        labels.append(r'$n = %i$' % i) #Make pretty ledged 
        residlist.append(iterfit)
        orderlist.append(i)
    plt.subplot2grid((5,3),(j,1))
    chi2 = []
    lines = []
    for i in range(0,len(residlist)):  #Plot the residual and chi2
        plt.xlim((2,8))
        plt.ylabel("Residual")
	    #calculate residuals
        resids = np.subtract(trans[j],residlist[i](co2[:,0])) #resjdhals
        lines.append(plt.plot(co2[:,0],resids))
        chi2.append(tealatte(residlist[i](co2[:,0]).tolist(),trans[j])[0])  #chi2
        if j == len(trans)-1:
            plt.xlabel("Year - 1998")
    
    residlist = []
    plt.subplot2grid((5,3),(j,2))
    plt.ylabel("Chi Squared")
    plt.scatter(orderlist,chi2)
    if j == len(trans)-1:
        plt.xlabel("Polynomial Order")
    chi2 = []
    orderlist = []


plt.show()



