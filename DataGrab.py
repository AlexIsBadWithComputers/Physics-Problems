#python thing to grab C02 data and performa couple cute little calculations

import math

def fileyoink(filename): #grab the data I want
	with open(str(filename),"r") as f: 
		search = f.readlines()
    	f.close

	data = []
	ReadBegin = False
	for i, line in enumerate(search):
		if ReadBegin ==True:
			info = line.split()
			if "-99.99" in line:
				continue #vad data, go away
			else:
				ListRipper(info,0,1,6)
				info = map(float,info)
				data.append(info)
		if "#             date" in line:
			ReadBegin = True
	return data
	f.close()

def ListRipper(list_,*args): #Deletes list pieces by index.
	indexes = sorted(list(args), reverse = True)
	for index in indexes:
		del list_[index]
	return list_


co2 = fileyoink("co2dat.txt")
for i in range(0,10):
	print co2[i][0], co2[i][1],co2[i][2],co2[i][3]
f = open("Co2Plots.txt","w")
for i in range(0,len(co2)):
	f.write("{} {}\n".format(co2[i][0],co2[i][1]))

f.write("\n")
f.write("\n")

for i in range(0,len(co2)):
	f.write("{} {}\n".format(co2[i][0],math.log10(co2[i][1])))

f.write("\n")
f.write("\n")

for i in range(0,len(co2)):
	f.write("{} {}\n".format(co2[i][0],math.sqrt(co2[i][1])))

f.write("\n")
f.write("\n")

for i in range(0,len(co2)):
	f.write("{} {}\n".format(co2[i][0],math.log10(co2[i][3])))

f.close()