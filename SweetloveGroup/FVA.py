#Function to perform FVA analysis which maintains sum of fluxes at a minimal val-
#ue
#args: 1) a cobra model 2) Objective 3) reaction to avoid when constraining sum 
#of fluxes 4) reaction list for FVA 5) solver used to perform FVA
#output: a cobra model with FVA as an attribute called fva
def FBA_FVA_run(cobra_model,obj,rxn2avoid = [],rxnlist=[],solver="",weightings={}):
  from SweetloveGroup.transform import rev2irrev
  from SweetloveGroup.constraints import constrainSumOfFluxes
  from SweetloveGroup.FBA import pfba_Weighted
  from cobra import flux_analysis
  
  
  if len(rxnlist)==0:
    rxnlist = cobra_model.reactions
  for rxn in cobra_model.reactions:
    if not rxn.id in weightings.keys():
      weightings[rxn.id]=1
      print("Warning")
      print rxn.id
      return
  print("Runing pFBA")
  solution = pfba_Weighted(cobra_model,weightings)
  objvalue = solution.x_dict.get(obj.id)
  a = 0
  for i in cobra_model.reactions:
    a = a + abs(solution.x_dict.get(i.id)*weightings[i.id])
  
  sumOfFluxes = a
  
  cobra_model2 = cobra_model.copy()
  irr_model = rev2irrev(cobra_model2)
  print("Setting SOF model")
  sfmodel = constrainSumOfFluxes(irr_model,rxn2avoid,sumOfFluxes,objvalue,weightings)
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
  
  if solver != "":
    import optlang
    if optlang.available_solvers.keys().__contains__(solver) and optlang.available_solvers[solver]:
      sfmodel.solver=solver
    else:
      print("Requested solver "+solver+" not available, using current model solver...")
  fva = flux_analysis.flux_variability_analysis(sfmodel,reaction_list = rxnlist2)
  print("Processing results")
  
  fva2=dict()
  for mode in fva.keys():
    if mode == "maximum":
      tempdict = dict()
      FVArxnSet = set()
      for rxn in fva[mode].keys():
        if rxn.__contains__("_reverse"):
          rxn = rxn.replace("_reverse","")
        if FVArxnSet.__contains__(rxn):
          continue
        FVArxnSet.add(rxn)
        if not fva[mode].keys().__contains__(rxn+"_reverse"):
          maxi = fva[mode][rxn]
        else:
          maxi = fva[mode][rxn]+fva[mode][rxn+"_reverse"]
        tempdict[rxn]=maxi
    else:
      tempdict=dict()
      FVArxnSet = set()
      for rxn in fva[mode].keys():
        if rxn.__contains__("_reverse"):
          rxn = rxn.replace("_reverse","")
        if FVArxnSet.__contains__(rxn):
          continue
        FVArxnSet.add(rxn)
        if not fva[mode].keys().__contains__(rxn+"_reverse"):
          mini = fva[mode][rxn]
        else:
          mini = fva[mode][rxn]+fva[mode][rxn+"_reverse"]
        tempdict[rxn]=mini
    fva2[mode]=tempdict
  
  sfmodel.fva = fva
  cobra_model.fva = fva2
  return cobra_model
