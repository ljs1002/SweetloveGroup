#################################################################################
#These functions are used to process biological data
#Author: Sanu Shameer 					PI: Lee J. Sweetlove	#
#email: sanushameer@gmail.com				Last updated: 05/09/2017#
#################################################################################

#This function processes raw biomass data and returns upper and lower bounds to 
#be used in constraint scanning as done by Colombie et al 2015
#	!!!Anaconda users need to uncomment %matplotlib inline!!!		                  #
#			!!!Script is still in development!!!			                                #
def generateBoundsFromBiomass(datafile="/home/sanu/ColobieDataRaw.csv",met="sucrose",Rsampling = 1000):
  import pandas as pd
  import math
  import numpy as np
  
  
  
  ################################Main function####################################
  
  #import data
  df = pd.read_csv(datafile,sep="\t")
  Rsampling = 1000
  num2rem = int(0.025*Rsampling)
  
  x_values = list()
  y_values = list()
  
  
  for i in range(0,len(df)):
    if not df[met][i]== 0:
      x_values.append(df["DPA"][i])
      y_values.append(df[met][i])
  
  
  log_y_values = list()
  for i in y_values:
    log_y_values.append(math.log(i))
  
  
  
  x2 = range(min(x_values),max(x_values)+1)
  y2 = np.poly1d(np.polyfit(x_values,log_y_values,3))
  
  ys = dict()
  
  for i in range(0,Rsampling):
    ind = np.random.choice(range(0,len(x_values)),size=50,replace=False)
    #print ind
    x=list()
    y=list()
    for j in ind:
      x.append(x_values[j])
      y.append(log_y_values[j])
    ys[i]=np.poly1d(np.polyfit(x,y,3))
  
  
  maxys = list()
  max95s = list()
  minys = list()
  min95s = list()
  temp=dict()
  for i in ys.keys():
    temp[i]=ys[i](x2)
  
  yi = dict()
  for j in range(0,len(x2)):
    templist=list()
    for i in temp.keys():
      templist.append(temp[i][j])
    temp2list=sorted(templist)[num2rem:Rsampling-num2rem]
    yi[j]=templist
    maxys.append(max(templist))
    max95s.append(max(temp2list))
    minys.append(min(templist))
    min95s.append(min(temp2list))
  
  
  
  
  y3 = list()
  for i in y2(x2):
    y3.append(math.exp(i))
  
  maxy1 = list()
  for i in maxys:
    maxy1.append(math.exp(i))
  
  
  miny1 = list()
  for i in minys:
    miny1.append(math.exp(i))
  
  maxy95_1 = list()
  for i in max95s:
    maxy95_1.append(math.exp(i))
  
  
  miny95_1 = list()
  for i in min95s:
    miny95_1.append(math.exp(i))
  
  
  
  
  import matplotlib.pyplot as plt
  #%matplotlib inline								#uncomment this line when using matplotlib
  plt.rcParams.update({'font.size': 20}) #sets a global fontsize
  plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
  plt.rcParams['xtick.major.width'] = 1
  plt.rcParams['ytick.major.size'] = 5
  plt.rcParams['ytick.major.width'] = 1
  plt.rcParams['axes.linewidth']=3 # makes axes line thicker
  
  
  ax = plt.subplot()
  ax.plot(x_values,y_values,".",label="raw data")
  ax.plot(x2,y3,"g-",label="fitted curve")
  ax.plot(x2,maxy1,"r--",label="fitted curve with max. y")
  ax.plot(x2,miny1,"b--",label="fitted curve with min. y")
  plt.xlabel("Time(days)")
  plt.ylabel("Concentration (micromol/fruit)")
  plt.legend(loc="best",fontsize=15)
  plt.show()
  
  ##############################derivatives#####################
  
  y2_deriv = y2.deriv()
  
  
  ys_deriv = dict()
  for i in ys.keys():
    ys_deriv[i] = ys[i].deriv()
  
  maxys_deriv = list()
  maxy95_deriv = list()
  minys_deriv = list()
  miny95_deriv = list()
  temp=dict()
  for i in ys.keys():
    temp[i]=ys_deriv[i](x2)
  
  
  yi_deriv = dict()
  for j in range(0,len(x2)):
    templist=list()
    for i in temp.keys():
      templist.append(temp[i][j])
    temp2list=sorted(templist)[num2rem:Rsampling-num2rem]
    yi_deriv[j]=templist
    maxys_deriv.append(max(templist))
    maxy95_deriv.append(max(temp2list))
    minys_deriv.append(min(templist))
    miny95_deriv.append(min(temp2list))
    
  
  
  
  y4 = list()
  for i in range(0,len(x2)):
    y4.append(math.exp(y2(x2[i]))*y2_deriv(x2)[i])
  
  maxy2 = list()
  maxy95_2 = list()
  for i in range(0,len(x2)):
    maxy2.append(maxy1[i]*maxys_deriv[i])
    maxy95_2.append(maxy95_1[i]*maxy95_deriv[i])
  
  
  miny2 = list()
  miny95_2 = list()
  for i in range(0,len(x2)):
    miny2.append(miny1[i]*minys_deriv[i])
    miny95_2.append(miny95_1[i]*miny95_deriv[i])
  
  
  
  import matplotlib.pyplot as plt
  #%matplotlib inline								#uncomment this line when using matplotlib
  plt.rcParams.update({'font.size': 20}) #sets a global fontsize
  plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
  plt.rcParams['xtick.major.width'] = 1
  plt.rcParams['ytick.major.size'] = 5
  plt.rcParams['ytick.major.width'] = 1
  plt.rcParams['axes.linewidth']=3 # makes axes line thicker
  
  
  ax = plt.subplot()
  ax.plot(x2,y4,"g-",label="derivitative of fitted curve")
  ax.plot(x2,maxy2,"r-",label="derivative of fitted curve with max.y")
  ax.plot(x2,miny2,"b-",label="derivative of fitted curve with min.y")
  ax.plot(x2,maxy95_2,"r--",label="derivative of fitted curve with max. y 95%")
  ax.plot(x2,miny95_2,"b--",label="derivative of fitted curve with min. y 95%")
  plt.xlabel("Time(days)")
  plt.ylabel("Flux (micromol/fruit/day)")
  plt.legend(loc="best",fontsize=15)
  plt.show()
  
  
