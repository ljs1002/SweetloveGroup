#Function to perform FVA analysis which maintains sum of fluxes at a minimal val-
#ue
#args: 1) a cobra model 2) Objective 3) reaction to avoid when constraining sum 
#of fluxes 4) reaction list for FVA
#output: a cobra model with FVA as an attribute called fva
def FBA_FVA_run(cobra_model,obj,rxn2avoid = [],rxnlist=[]):
  from SweetloveGroup.transform import rev2irrev
  from SweetloveGroup.constraints import constrainSumOfFluxes
  from cobra import flux_analysis
  
  
  if len(rxnlist)==0:
    rxnlist = cobra_model.reactions
  #print("Rxn list ="+str(rxnlist))
  print("Runing pFBA")
  solution = flux_analysis.parsimonious.pfba(cobra_model)
  objvalue = solution.x_dict.get(obj.id)
  a = 0
  for i in cobra_model.reactions:
    a = a + abs(solution.x_dict.get(i.id))
  
  sumOfFluxes = a
  
  cobra_model2 = cobra_model.copy()
  irr_model = rev2irrev(cobra_model2)
  print("Setting SOF model")
  sfmodel = constrainSumOfFluxes(irr_model,rxn2avoid,sumOfFluxes,objvalue)
  rxnlist2 = list()
  if rxnlist == cobra_model.reactions:
    rxnlist2 = sfmodel.reactions
  else:
    for rxn in rxnlist:
      if rxn.lower_bound<0 and rxn.upper_bound>0:
        rxnlist2.append(sfmodel.reactions.get_by_id(rxn.id+"_reverse"))
      rxnlist2.append(sfmodel.reactions.get_by_id(rxn.id))
  #print("Rxn list ="+str(rxnlist2))
  print("Running FVA")
  sfmodel.solver="cplex"
  fva = flux_analysis.flux_variability_analysis(sfmodel,reaction_list = rxnlist2)
  print("Processing results")
  #print("FVA ="+str(fva))
  FVArxnSet = set()
  tempdict=dict()
  for rxn in fva.keys():
    if FVArxnSet.__contains__(rxn):
      continue
    if rxn.__contains__("reverse"):
      rxn = rxn.replace("_reverse","")
    FVArxnSet.add(rxn)
    if not fva.keys().__contains__(rxn+"_reverse"):
      tempdict[rxn]=fva.get(rxn)
      continue
    FVArxnSet.add(rxn+"_reverse")
    maxi = fva.get(rxn).get("maximum")# + fva.get(rxn+"_reverse").get("minimum")
    mini = fva.get(rxn+"_reverse").get("minimum")# + fva.get(rxn).get("maximum")
    if mini<maxi:
      tempdict[rxn]={"minimum":mini,"maximum":maxi}
    else:
      tempdict[rxn]={"minimum":maxi,"maximum":mini}
  
  sfmodel.fva = fva
  cobra_model.fva = tempdict
  return cobra_model


