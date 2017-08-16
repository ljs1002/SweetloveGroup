

#Function to identify net synthesis/consumption of metabolites in a given set of
#reactions 
#args: 1) a solved cobra model 2) a list of reactions from which net metabolite s-
#-toichiometry should be calculated
#output: a dictionay file with metabolite ids as key and net soichiometry as value
def netMetaboliteStoich(cobra_model,rxnlist):
  netMet = dict()
  for rxn in a:
    rxn = cobra_model.reactions.get_by_id(rxn)
    for met in rxn.metabolites:
      if netMet.keys().__contains__(met.id):
        netMet[met.id]=netMet[met.id]+((rxn.x/abs(rxn.x))*rxn.metabolites.get(met))
      else:
        netMet[met.id]=((rxn.x/abs(rxn.x))*rxn.metabolites.get(met))
  return netMet
