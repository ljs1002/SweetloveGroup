

#Function to identify net synthesis/consumption of metabolites in a given set of
#reactions 
#args: 1) a solved cobra model 2) a list of reactions from which net metabolite s-
#-toichiometry should be calculated
#output: a dictionay file with metabolite ids as key and net soichiometry as value
def netMetaboliteStoich(cobra_model,rxnlist):
  netMet = dict()
  for rxn in rxnlist:
    rxn = cobra_model.reactions.get_by_id(rxn)
    if round(rxn.x,6)==0:
      print rxn.id+" flux is 0."
      netMet = dict()
      break
    for met in rxn.metabolites:
      if netMet.keys().__contains__(met.id):
        netMet[met.id]=netMet[met.id]+((rxn.x/abs(rxn.x))*rxn.metabolites.get(met))
      else:
        netMet[met.id]=((rxn.x/abs(rxn.x))*rxn.metabolites.get(met))
  return netMet


#Function to print out all reactions generating/consuming a metabolite of inter-
#est 
#args: 1) a solved cobra model 2) metabolite ID 3) ID of alternate charged state
#(use "" if none) 4) output file (use "" if no output file is required)
#output: none
def writeMetabSummary(cobra_model, met, Amet, outfile):
  met=cobra_model.metabolites.get_by_id(met)
  if not Amet == "":
    Amet=cobra_model.metabolites.get_by_id(Amet)
  if not outfile=="":
    fout = open(outfile,"w")
    fout.write("rxn ID\treaction\tmetabolite flux\n")
  for rxn in met.reactions:
    sto=rxn.metabolites.get(met)
    if Amet=="":
      sto1=0
    else:
      sto1=rxn.metabolites.get(Amet)
    if outfile=="":
      print  rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*(sto+sto1))
    else:
      fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*(sto+sto1))+"\n")
  if not outfile=="":
    fout.close()
      
      
#Function to calculate night time carbon conversion efficiency in diel model
#args: 1) a solved cobra model 2) day-night accumulation tag (default = "diel-
#-Transfer") 3) reaction ID for night time output (default = phloem_output_tx-
#-2), 4) reaction ID representing CO2 respired
#output: carbon conversion efficiency
def predictCCE(C3_model,accumulation_tag="dielTransfer",output="Phloem_output_tx2",CO2rxn = "CO2_tx2"):
  import re
  
  for met in C3_model.metabolites:
    if not met.formula:
      met.formula=""
  
  Cin = 0
  Cout = 0
  
  for rxn in C3_model.reactions.query(accumulation_tag):
    if round(rxn.x,5)>0:
      for met in rxn.products:
        #print met
        if met.formula.__contains__("C"):
          #print str(Cin)+"---"+met.id
          #print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cin = Cin + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  for rxn in C3_model.reactions.query(accumulation_tag):
    if round(rxn.x,5)<0:
      for met in rxn.reactants:
        if met.formula.__contains__("C"):
          #print str(Cout)+"---"+met.id
          #print str(rxn.x)+"\t"+str(rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
          Cout = Cout + (rxn.x * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  rxn = C3_model.reactions.get_by_id(output)
  for met in rxn.reactants:
    if met.formula.__contains__("C"):
      #print str(Cout)+"---"+met.id
      #print str(rxn.x)+"\t"+str(-1 * rxn.metabolites.get(met))+"\t"+str(int(re.split(r"[C,H]",met.formula)[1]))
      Cout = Cout + (rxn.x * -1 * rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1]))
  
  Cout = Cout + (-1*C3_model.reactions.get_by_id(CO2rxn).x)
  
  if(not round(Cin,5) == round(Cout,5)):
    print "Error, Cin = "+str(Cin)+" and Cout = "+str(Cout)
    return 0
  else:
    print "Cin = "+str(Cin)+"\tCO2 ="+str(C3_model.reactions.get_by_id(CO2rxn).x)+"\t"+str(1 + ((C3_model.reactions.get_by_id(CO2rxn).x)/Cin))
    return 1 + ((C3_model.reactions.get_by_id("CO2_tx2").x)/Cin)



####################################################################
#This function generates a tab seperated file that can be used with#
#Cytoscape to visualize metabolic flux                             #
#                                                                  #
#inputs: 1) a cobra model with feasible solution 2) the name of the#
#output file  3)  the number of cells in the model (eg: 2 for diel #
#C3 and 4 for diel C4)                                             #
#                                                                  #
####################################################################

def generateFluxMap(cobra_model, outfile,phases = 2):
    import cobra
    solution = cobra.flux_analysis.parsimonious.optimize_minimal_flux(cobra_model)
    #solution = cobra.flux_analysis.parsimonious.pfba(cobra_model)          #If the previous line returns error comment it out and uncomment this line instead

    #open output file for writing
    f = open(outfile,"w");

    #use rxnSet to identify reaction that have already been processed
    rxnSet = set()

    mult=set()
    #Looping through all reactions in the model
    for rxn in cobra_model.reactions:
        #Get the ID
        RXN=rxn.id
        #declare a boolean variable multFlag to keep track of whether this reaction is present in multiple models
        multFlag=False

        #check if the reaction has already been processed before and if yes skip this run in the loop
        if(rxnSet.__contains__(RXN)):
            continue
        if rxn.id.__contains__("EX") or rxn.id.__contains__("Transfer"):
            multFlag=False
        #check if the reaction ends with one or two i.e it is present more than once in the model
        elif(["1","2","3","4","5","6","7","8","9"].__contains__(rxn.id[len(rxn.id)-1])):
            #change the id to without the suffix 1-9 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-1]
            multFlag=True
        elif rxn.id[len(rxn.id)-2:] == "10":
            #change the id to without the suffix 10 and declare it as a reaction which has multiple instances
            RXN = rxn.id[0:len(rxn.id)-2]
            multFlag=True

        #if metabolite has multiple instances
        values = dict()
        status1 = dict()
        status2 = dict()
        if(multFlag):
            tempvalue = list()
            temp1 = list()
            temp2 = list()
            mult.add(RXN)
            #add the reaction we are about to process to the reactions processed list
            for i in range(1,phases+1):
                rxnSet.add(RXN+str(i))
                if(round(float(cobra_model.reactions.get_by_id(RXN+str(i)).x)*10000) == 0):
                    tempvalue.append(0)
                    temp1.append("none")
                    temp2.append("none")
                elif(float(cobra_model.reactions.get_by_id(RXN+str(i)).x)*10000 > 0):
                    tempvalue.append(cobra_model.reactions.get_by_id(RXN+str(i)).x*1000)
                    temp1.append("produced")
                    temp2.append("consumed")
                elif(float(cobra_model.reactions.get_by_id(RXN+str(i)).x)*10000 < 0):
                    tempvalue.append(cobra_model.reactions.get_by_id(RXN+str(i)).x*1000)
                    temp1.append("consumed")
                    temp2.append("produced")
            values[RXN] = tempvalue
            status1[RXN] = temp1
            status2[RXN] = temp2

            #select 1 reaction so that we can identify the reactants and products which can be then used to generate the edge shared_name
            rxn=cobra_model.reactions.get_by_id(RXN+"1")

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):
                        REAC=REAC[0:len(REAC)-1]
                    f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                    f.write("\n")
                if(RXN.__contains__("biomass")):
                    f.write("R_"+RXN+" (reaction-product)) M_"+REAC)
                    for i in range(1,phases+1):
                        f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                    f.write("\n")
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_"+RXN+" (reaction-product) M_"+PROD)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status1[RXN][i-1]))
                f.write("\n")
            if(RXN.__contains__("biomass")):
                f.write("R_"+RXN+" (reaction-reactant) M_"+REAC)
                for i in range(1,phases+1):
                    f.write("\t"+str(values[RXN][i-1])+"\t"+str(status2[RXN][i-1]))
                f.write("\n")
        else:
            #add the reaction we are about to process to the reactions processed list
            rxnSet.add(RXN)
            if(round(float(solution.x_dict.get(rxn.id))*10000) == 0):
                value = 0;
                status1= "none";
                status0= "none";
            elif(solution.x_dict.get(rxn.id)*10000 > 0):
                value = solution.x_dict.get(rxn.id)*1000;
                status1= "produced";
                status0= "consumed";
            elif(solution.x_dict.get(rxn.id)*10000 < 0):
                value = solution.x_dict.get(rxn.id)*1000;
                status1= "consumed";
                status0= "produced";

            for reac in rxn.reactants:
                REAC=reac.id
                if(REAC.__contains__("1")):
                    if(REAC.rindex("1")==len(REAC)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        REAC=REAC[0:len(REAC)-1]
                f.write("R_%s (reaction-reactant) M_%s\t%s\t%s\t0\tnone\n" % (RXN,REAC,value,status0));
            for prod in rxn.products:
                PROD=prod.id
                if(PROD.__contains__("1")):
                    if(PROD.rindex("1")==len(PROD)-1):# or (met.id.rindex("2")==len(rxn.id)-1):
                        PROD=PROD[0:len(PROD)-1]
                f.write("R_%s (reaction-product) M_%s\t%s\t%s\t0\tnone\n" % (RXN,PROD,value,status1));

    f.close();



####################################################
# This function estimates Rubisco carboxylase flux #
# at which the net CO2 uptake rate is equal to the #
# user defined value                               #
####################################################
def estimateVcFromNetCO2(model,netCO2uptake,Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1",CO2in_ID="CO2_tx1",verbose=False):
    
    from cobra import flux_analysis
    # Initally constraint Vc flux to net CO2 uptake rate
    model.reactions.get_by_id(Vc_ID).lower_bound = netCO2uptake
    model.reactions.get_by_id(Vc_ID).upper_bound = netCO2uptake
    
    #perform pFBA
    flux_analysis.parsimonious.optimize_minimal_flux(model)
    
    #set loop counter
    i=0
    
    #Use a while loop to increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
    while((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake > 0.001 and i<10):
        i=i+1
        prev = model.reactions.get_by_id(Vc_ID).x
        # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
        now = prev + ((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x))
        model.reactions.get_by_id(Vc_ID).lower_bound = now
        model.reactions.get_by_id(Vc_ID).upper_bound = now
        flux_analysis.parsimonious.optimize_minimal_flux(model)
        if verbose:
            print("----"+str(i)+"----")
            print("Vc flux ="+str(now))
            print("net CO2 uptake ="+str(model.reactions.get_by_id(CO2in_ID).x))
            print("Target CO2 uptake ="+str(netCO2uptake))
    return prev
  

######################################################
# This function estimates biomass/phloem output flux #
# at which the net CO2 uptake rate is equal to the   #
# user defined value                                 #
####################################################
def estimateOutputFromNetCO2(model,netCO2uptake,Output_ID="diel_biomass",Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1",CO2in_ID="CO2_tx1",verbose=False):
    
    from cobra import flux_analysis
    # Initally constraint Vc flux to net CO2 uptake rate
    model.reactions.get_by_id(Vc_ID).lower_bound = netCO2uptake
    model.reactions.get_by_id(Vc_ID).upper_bound = netCO2uptake
    
    #perform pFBA
    flux_analysis.parsimonious.pfba(model)
    
    #unconstrain Vc
    model.reactions.get_by_id(Vc_ID).lower_bound = 0
    model.reactions.get_by_id(Vc_ID).upper_bound = 1000
    
    #set loop counter
    i=0
    
    #Use a while loop to increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
    while((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake > 0.001):# and i<10):
        i=i+1
        prev = model.reactions.get_by_id(Output_ID).x
        # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
        now = prev + (prev*((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake))
        model.reactions.get_by_id(Output_ID).lower_bound = now
        model.reactions.get_by_id(Output_ID).upper_bound = now
        
        flux_analysis.parsimonious.pfba(model)
        if verbose:
            print("----"+str(i)+"----")
            print("Vc flux ="+str(model.reactions.get_by_id(Vc_ID).x))
            print("net CO2 uptake ="+str(model.reactions.get_by_id(CO2in_ID).x))
            print("Target CO2 uptake ="+str(netCO2uptake))
            print("Before:"+str(prev))
            print("After:"+str(now))
    return prev


  
#####################################################################
# This function prints information on CO2 generating fluxes and can #
# be used to study dark respiration                                 #
# inputs: 1) an FBA solved model, 2) a list of fluxs to ignore, 3)  #
# the tag used to mark night time metabolites ex: "_night", 4)      #
# a list of solution objects                                        #
#####################################################################
def printDarkRespirationFluxes(model,rxn2avoid=["CO2_tx1","CO2_ec1","CO2_mc1","CO2_pc1","GCVMULTI_RXN_m1"],night_tag = "2",custom_solutions=[]):
  for met in model.metabolites.query("CARBON_DIOXIDE"):
    if met.id.__contains__(night_tag):
      continue
    for rxn in met.reactions:
      if rxn2avoid.__contains__(rxn.id):
        continue
      if len(custom_solutions) == 0:
        if (rxn.metabolites.get(met) > 0 and rxn.x > 0) or (rxn.metabolites.get(met) < 0 and rxn.x < 0):
          print rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*rxn.metabolites.get(met))
        else:
          if (rxn.metabolites.get(met) > 0 and rxn.upper_bound > 0) or (rxn.metabolites.get(met) < 0 and rxn.lower_bound < 0):
            tot=0
            for sol in custom_solutions:
              tot=tot+abs(sol.x_dict.get(rxn.id))
              if tot==0:
                continue
            
            print rxn.id+"\t"+rxn.reaction,
            for sol in custom_solutions:
              print "\t"+str(rxn.metabolites.get(met)*sol.x_dict.get(rxn.id)),
            print ""

#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose wether to use percentage or absolute values in the plot.   #
#####################################################################
def generateATPbudget(model,solution,outfile="",show_plot=True,percentage=False):
  if outfile!="":
    fout = open(outfile,"w")
  ATPdict = dict()
  total = 0
  for p in ("c","p","m","x"):
    met=model.metabolites.get_by_id("ATP_"+p+"1")
    met1=model.metabolites.get_by_id("aATP_"+p+"1")
    for rxn in met.reactions:
      if rxn.id.__contains__("ATP_AMP_mc") or rxn.id.__contains__("ATP_ADP_mc") or rxn.id.__contains__("ATP_pc") or rxn.id.__contains__("AMP_ATP_xc") or rxn.id.__contains__("ATP_ADP_Pi_pc"): 
        continue
      sto=rxn.metabolites.get(met)
      sto1=rxn.metabolites.get(met1)
      if outfile!="":
        fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
      ATPdict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
      if solution.get(rxn.id)*(sto+sto1) > 0:
        total = total + (solution.get(rxn.id)*(sto+sto1))
  if outfile!="":
    fout.close()
  
  ATPdict2 = dict()
  ATPdict2["Others-pos"]=0
  ATPdict2["Others-neg"]=0
  baseline = dict()
  pos_base=0
  neg_base=0
  i=0
  for rxn in ATPdict.keys():
    if ATPdict[rxn]>0:
      if ATPdict[rxn] < total*0.05:
        if percentage:
          ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+float(ATPdict[rxn]*100)/total
        else:
          ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+ATPdict[rxn]
        continue
      base = pos_base
      if percentage:
        ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
        pos_base = pos_base + float(ATPdict[rxn]*100)/total
      else:
        pos_base = pos_base + ATPdict[rxn]
        ATPdict2[rxn]=ATPdict[rxn]
    else:
      if abs(ATPdict[rxn]) < total*0.05:
        if percentage:
          ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+float(ATPdict[rxn]*100)/total
        else:
          ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+ATPdict[rxn]
        continue
      base = neg_base
      if percentage:
        ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
        neg_base = neg_base + float(ATPdict[rxn]*100)/total
      else:
        neg_base = neg_base + ATPdict[rxn]
        ATPdict2[rxn]=ATPdict[rxn]
    i=i+1
    baseline[rxn]=base
  baseline["Others-pos"]=pos_base
  baseline["Others-neg"]=neg_base
  
  if show_plot:
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.figure(figsize=(3,4))
    for rxn in ATPdict2.keys():
      plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
    plt.axhline(0,linestyle="--",color="black")
    plt.xlim(0.8,1.2)
    if percentage:
      plt.ylabel("ATP produced/consumed (%)")
    else:
      plt.ylabel("ATP produced/consumed (in moles)")
    lgd=plt.legend(bbox_to_anchor=(1,1))
    plt.tight_layout
    plt.savefig('temp.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
 
