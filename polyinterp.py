import numpy as np 
import matplotlib.pyplot as plt

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
def chi(p,x,y,dof):
	return np.sum((np.polyval(p,x)-y)**2.0)/dof

co2 = np.array(fileyoink("co2dat.txt")) #Data to fit to
co2points = np.array(fileyoink("co2dat2.txt")) #Extra data to plot... but not to fit to
co2points[:,0] = co2points[:,0]-1998
co2[:,0] = co2[:,0]-1998 #Normalize data but avoid negative values below datas


xp = np.linspace(0,10,1000)


colormap = plt.cm.gist_ncar # I demand a pretty plot. No uggos here. 
plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, 8)])

labels = []  #Just set some junk up. 
plt.subplot(2,2,1)
plt.scatter(co2points[:,0],co2points[:,1])
plt.ylim((360,390))
plt.xlim((0,10))
plt.xlabel(" Year - 1998")
plt.ylabel("Carbon Dioxide (ppm)")
location  = 1
for i in range(0,73,2):
	if i == 0:
		continue
	
	fit = np.polyfit(co2[:,0],co2[:,1],i)
	print(chi(fit,co2[:,0],co2[:,1],1), i,"Hellllloooooo")
	iterfit = np.poly1d(fit)
	plt.plot(xp,iterfit(xp))
	labels.append(r'$n = %i$' % i) #Make pretty ledged 

	
	if i == 16 or i == 36 or i == 54: #Switch which subplot to put the interpolations
		location = location + 1       #would have been prettier with modulus division
		plt.subplot(2,2,location)	   #but I didn't feel like thinking that hard
		plt.scatter(co2points[:,0],co2points[:,1])
		plt.ylim((360,390))
		plt.xlim((0,10))
		plt.xlabel(" Year - 1998")
		plt.ylabel("Carbon Dioxide (ppm)")
	

#Below is just setting up the legend

plt.ylim((360,390))
plt.xlim((0,10))
plt.subplot(2,2,1)
plt.legend(labels[0:8], loc= "upper center", ncol = 4)
plt.subplot(2,2,2)
plt.legend(labels[9:17], loc= "upper center", ncol = 4)
plt.subplot(2,2,3)

plt.legend(labels[18:26],loc = "upper center", ncol =4)
plt.subplot(2,2,4)
plt.legend(labels[27:35], loc= "upper center", ncol = 4)

plt.show()



