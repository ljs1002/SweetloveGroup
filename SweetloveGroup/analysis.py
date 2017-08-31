

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







  
