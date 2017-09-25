def rev2irrev(cobra_model):
  exp_model=cobra_model.copy()
  
  for RXN in cobra_model.reactions:
    rxn=exp_model.reactions.get_by_id(RXN.id)
    if (rxn.lower_bound < 0):
      rxn_reverse = rxn.copy()
      rxn_reverse.id = "%s_reverse" %(rxn.id)
      rxn.lower_bound = 0
      rxn_reverse.upper_bound = 0
      exp_model.add_reaction(rxn_reverse)
  
  return exp_model
