
#Function to perform FVA analysis which maintains sum of fluxes at a minimal val-
#ue
#args: 1) a cobra model
#output: a cobra model with FVA as an attribute called fva
#Function to perform FVA analysis which maintains sum of fluxes at a minimal val-
#ue
#args: 1) a cobra model
#output: a cobra model with FVA as an attribute called fva
def FBA_FVA_run(cobra_model):
  objvalue = cobra_model.solution.f
  a = 0
  for i in solution.x_dict.keys():
    a = a + abs(solution.x_dict.get(i))
  
  sumOfFluxes = a
  
  cobra_model2 = cobra_model.copy()
  irr_model = rev2irrev(cobra_model2)
  sfmodel = constrainSumOfFluxes(irr_model,[],sumOfFluxes,objvalue)
  sfmodel.optimize()
  
  fva = flux_analysis.flux_variability_analysis(sfmodel)
  
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
    maxi = fva.get(rxn).get("maximum") + fva.get(rxn+"_reverse").get("minimum")
    mini = fva.get(rxn).get("minimum") + fva.get(rxn+"_reverse").get("maximum")
    if mini<maxi:
      tempdict[rxn]={"minimum":mini,"maximum":maxi}
    else:
      tempdict[rxn]={"minimum":maxi,"maximum":mini}
  
  sfmodel.fva = fva
  cobra_model.fva = tempdict
  return cobra_model


