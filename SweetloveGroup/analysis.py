

#Function to identify net synthesis/consumption of metabolites in a given set of
#reactions 
#args: 1) a solved cobra model 2) a list of reactions from which net metabolite s-
#-toichiometry should be calculated
#output: a dictionay file with metabolite ids as key and net soichiometry as value
def netMetaboliteStoich(cobra_model,rxnlist):
  netMet = dict()
  for rxn in a:
    rxn = cobra_mode.reactions.get_by_id(rxn)
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
    fout.write("rxn ID\treaction\tmetabolite flux")
  for rxn in met.reactions:
    sto=rxn.metabolites.get(met)
    if Amet=="":
      sto1=0
    else:
      sto1=rxn.metabolites.get(met1)
    if outfile=="":
      print  rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*(sto+sto1))
    else:
      fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(rxn.x*(sto+sto1))+"\n")
  if not outfile=="":
    fout.close()
      
      


  
