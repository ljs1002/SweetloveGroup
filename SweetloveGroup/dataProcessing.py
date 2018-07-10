#################################################################################
#This function processes raw biomass data and returns upper and lower bounds to #
#be used in constraint scanning as done by Colombie et al 2015                  #
#                                                                               #
#################################################################################
def generateBoundsFromBiomass(datafile="/home/sanu/ColobieDataRaw.csv",met="sucrose",Nsampling = 1000,DPA=(11.4,14.8,18.2,21.6,25,28.4,31.8,25.2,38.6,41.3),show_plots=True,start=4.0,stop=57.0,degree=3,Ssampling=0.75):
  import pandas as pd
  import math
  import numpy as np
  
  conc=list()
  rate_max=list()
  rate_min=list()
    
  #import data
  df = pd.read_csv(datafile,sep="\t")
  num2rem = int(0.025*Nsampling)
  
  x_values = list()
  y_values = list()
  
  
  for i in range(0,len(df)):
    if not df[met][i]== 0:
      x_values.append(df["DPA"][i])
      y_values.append(df[met][i])
  
  
  log_y_values = list()
  for i in y_values:
    log_y_values.append(math.log(i))
  
  
  
  x2 = np.arange(start,stop+0.1,0.1)
  y2 = np.poly1d(np.polyfit(x_values,log_y_values,degree))
  
  ys = dict()
  
  for i in range(0,Nsampling):
    ind = np.random.choice(range(0,len(x_values)),size=int(len(x_values)*Ssampling),replace=False)
    #print ind
    x=list()
    y=list()
    for j in ind:
      x.append(x_values[j])
      y.append(log_y_values[j])
    ys[i]=np.poly1d(np.polyfit(x,y,degree))
  
  
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
    temp2list=sorted(templist)[num2rem:Nsampling-num2rem]
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
  
  #print y3
  #print "----"
  #print x2
  
  for i in DPA:
    #print i
    conc.append(y3[int((i-min(x_values))*10)])
  
  if show_plots:
    import matplotlib.pyplot as plt
    %matplotlib inline
    plt.rcParams.update({'font.size': 20}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=3 # makes axes line thicker
    plt.figure(figsize=(5,5))
    
    
    ax = plt.subplot()
    ax.plot(x_values,y_values,".",label="raw data")
    ax.plot(x2,y3,"g-",label="fitted curve")
    ax.plot(x2,maxy1,"r--",label="fitted curve with max. y")
    ax.plot(x2,miny1,"b--",label="fitted curve with min. y")
    plt.xlabel("Time(days)")
    plt.ylabel("Concentration (micromol/fruit)")
    plt.legend(bbox_to_anchor=(2, 1),fontsize=15)
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
    templist=list()
    for j in range(0,len(x2)):
      templist.append(y3[j]*ys_deriv[i](x2[j]))
    temp[i]=templist
  
  
  yi_deriv = dict()
  for j in range(0,len(x2)):
    templist=list()
    for i in temp.keys():
      templist.append(temp[i][j])
    temp2list=sorted(templist)[num2rem:Nsampling-num2rem]
    yi_deriv[j]=templist
    maxys_deriv.append(max(templist))
    maxy95_deriv.append(max(temp2list))
    minys_deriv.append(min(templist))
    miny95_deriv.append(min(temp2list))
    
  
  
  
  y4 = list()
  for i in range(0,len(x2)):
    y4.append(y3[i]*y2_deriv(x2)[i])
  
  rate_max = list()
  rate_min = list()
  if not len(DPA)==0:
    j=0
    for i in DPA:
      rate_max.append(maxy95_deriv[int((i-min(x_values))*10)])
      rate_min.append(miny95_deriv[int((i-min(x_values))*10)])
      j=j+1
  
  
  if show_plots:
    import matplotlib.pyplot as plt
    
    %matplotlib inline
    plt.rcParams.update({'font.size': 20}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=3 # makes axes line thicker
    plt.figure(figsize=(5,5))
    
    
    ax = plt.subplot()
    ax.plot(x2,y4,"g-",label="derivitative of fitted curve")
    ax.plot(x2,maxys_deriv,"r-",label="derivative of fitted curve with max.y")
    ax.plot(x2,minys_deriv,"b-",label="derivative of fitted curve with min.y")
    ax.plot(x2,maxy95_deriv,"r--",label="derivative of fitted curve with max. y 95%")
    ax.plot(x2,miny95_deriv,"b--",label="derivative of fitted curve with min. y 95%")
    ax.set_title(met+" flux")
    plt.xlabel("Time(days)")
    plt.ylabel("Flux (micromol/fruit/day)")
    plt.legend(bbox_to_anchor=(2.55, 1),fontsize=15)
    plt.show()
  
  return (conc,rate_max,rate_min)


  
###################################################################
#This function removes gene and protein associations from a given #
#sbml file.                                                       #
###################################################################
def removeGeneProteinAssociations(orignal_sbml_file,final_sbml_file):
  fin =open(orignal_sbml_file)
  fout = open(final_sbml_file,"w")
  ignore=False
  Type=""
  reset_ignore_gene=False
  reset_ignore_protein=False
  for line in fin:
    if line.__contains__("<html:p>GENE_ASSOCIATION:") and not line.__contains__("<html:p>GENE_ASSOCIATION: </html:p>"):
      ignore=True
      Type="gene"
    if line.__contains__("<html:p>PROTEIN_ASSOCIATION:") and not line.__contains__("<html:p>PROTEIN_ASSOCIATION: </html:p>"):
      ignore=True
      Type="protein"
    if ignore and line.__contains__("</html:p>"):
      if Type=="gene":
        reset_ignore_gene = True
      elif Type == "protein":
        reset_ignore_protein=True
    if not ignore: 
      #print line
      fout.write(line)
    if reset_ignore_gene:
      fout.write("      <html:p>GENE_ASSOCIATION: </html:p>\n")
      ignore=False
      Type=""
      reset_ignore_gene = False
    if reset_ignore_protein:
      fout.write("      <html:p>PROTEIN_ASSOCIATION: </html:p>\n")
      ignore=False
      Type=""
      reset_ignore_protein = False  
  fout.close()
