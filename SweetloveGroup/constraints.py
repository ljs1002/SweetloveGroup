# This function can be used to set night time carbon conversion efficiency (CCE) in diel leaf model
# inputs: a model, CCE (default is 0.5), transfer reaction tag (default is dielTransfer)
# output: a model
# WARNING!! Please check IDs of CO2 exchange referenced in the function to prevent errors
def setCCE(C3_model,CCE = 0.5,tag = "dielTransfer"):
  import re
  MET = Metabolite("Cin")
  for met in C3_model.metabolites:
    if not met.formula:
      met.formula=""
  Cin = 0
  Cout = 0
  C3_model = rev2irrev(C3_model)
  for rxn in C3_model.reactions.query(tag):
    if not rxn.id.__contains__("_reverse"):
      for met in rxn.products:
        #print met
        if met.formula.__contains__("C"):
          rxn.add_metabolites({MET:rxn.metabolites.get(met) * int(re.split(r"[C,H]",met.formula)[1])})
      print rxn.reaction
  
  rxn = C3_model.reactions.get_by_id("CO2_tx2_reverse")
  met = C3_model.metabolites.get_by_id("CARBON_DIOXIDE_e2")
  rxn.add_metabolites({MET:1 * rxn.metabolites.get(met) * (float(1)/(1-CCE))})
  print rxn.reaction
  
  return C3_model


#Function to constraint sum of fluxes when performing FBA
#args: 1) a cobra model, 2) a python list of reactions to leave out from constrai-
#-nt, 3) the float value that sum of fluxes must be constrained to & 4) value obj-
#-ective function needs to be constraint to (provide "" to avoid constraining obj-
#ective function)
#output: a cobra model with sum of fluxes constrained to 
def constrainSumOfFluxes(cobra_model, rxn2avoid,SFvalue,objvalue):
  from cobra.core import Metabolite, Reaction
  
  temp=cobra_model.copy()
  SFMet = Metabolite("SFMet",name="Sum of fluxes pseudometabolite",compartment="c2")
  for rxn in cobra_model.reactions:
    if not rxn2avoid.__contains__(rxn.id):
      if rxn.id.__contains__("reverse"):
	temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:-1})
      else:
	temp.reactions.get_by_id(rxn.id).add_metabolites({SFMet:1})
  SFRxn = Reaction("SFRxn",name="Sum of fluxes pseudoreaction")
  SFRxn.add_metabolites({SFMet:-1})
  SFRxn.lower_bound=SFvalue
  SFRxn.upper_bound=SFvalue
  temp.add_reaction(SFRxn)
  if (not objvalue=="") and (len(temp.objective) == 1):
    for rxn in temp.objective.keys():
      rxn.lower_bound=objvalue
      rxn.upper_bound=objvalue
  return temp
