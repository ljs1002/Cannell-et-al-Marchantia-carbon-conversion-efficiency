### Cannell et al 2020.
# Functions included:
#   SetupDielModel is used in conjunction with list of metabolites provided on individual lines in a text file to produce a diel_model from a single phase model
#   fixModelCompatibilityIssueCobra015 is used to fix a compatibility issue with our models and cobrapy v0.15
#   estimateOutputFromNetCO2 is used to run pFBA on a stoichiometric model while disallowing the night-time fixation of CO2
###

def SetupDielModel(core_model,transferMets):
    from cobra.core import Metabolite, Reaction
    import re

    #create two copies of model elements for day and night
    cobra_model2 = core_model.copy()
    for met in cobra_model2.metabolites:
        met.id = met.id+"1"
        met.compartment = met.compartment+"1"
    for rxn in cobra_model2.reactions:
        rxn.id = rxn.id+"1"

    cobra_model3 = core_model.copy()
    for met in cobra_model3.metabolites:
        met.id = met.id+"2"
        met.compartment = met.compartment+"2"
    for rxn in cobra_model3.reactions:
        rxn.id = rxn.id+"2"

    #merge the day and night model
    cobra_model = cobra_model2+cobra_model3
    for met in cobra_model3.metabolites:
        if not cobra_model.metabolites.__contains__(met.id):
            cobra_model.add_metabolites(met.copy())

    met1 = Metabolite("Biomass_t1",name="Biomass during the day")
    cobra_model.reactions.get_by_id("AraCore_Biomass_tx1").add_metabolites({met1:1})
    met2 = Metabolite("Biomass_t2",name="Biomass during at night")
    cobra_model.reactions.get_by_id("AraCore_Biomass_tx2").add_metabolites({met2:1})
    
    rxn = Reaction("diel_biomass")
    rxn.add_metabolites({met1:-3,met2:-1})
    #rxn.add_metabolites({met1:-1,met2:-1})
    rxn.lower_bound = 0
    rxn.upper_bound = 1000
    cobra_model.add_reaction(rxn)

    #Adding reactions to allow for day-night metabolite accumulations
    tmfile = open(transferMets,"r")
    tmset=set()
    for line in tmfile:
        tmset.add(line.replace("\n",""))

    for met in tmset:
        if met == "AMMONIUM_v" or met=="FRUCTAN_v":
            continue
        tempRxn = Reaction(met+"_dielTransfer")
        tempRxn.add_metabolites({cobra_model.metabolites.get_by_id(met+"1"):-1,cobra_model.metabolites.get_by_id(met+"2"):1})
        tempRxn.lower_bound=-1000
        if not ((met == "STARCH_p") or (met == "SUCROSE_v") or (met == "MAL_v") or (met == "aMAL_v") or (met == "NITRATE_v") or (met == "CIT_v") or (met == "aCIT_v") or (met == "PROTON_v")):
            tempRxn.lower_bound=0
        tempRxn.upper_bound=1000
        cobra_model.add_reaction(tempRxn)

    fractionMets=dict()
    for rxn in cobra_model.reactions:
        for met in rxn.metabolites.keys():
            prefix=""
            a=re.search("^a{1,3}",met.id)
            anion=""
            if a:
                anion=a.group(0)
                prefix=anion
            b=re.search("^b{1,3}",met.id)
            basic=""
            if b:
                basic=b.group(0)
                prefix=basic
            if ((not prefix == "") and met.compartment == "v1"):
                fractionMets[met]=prefix

    temp=cobra_model.copy()
    for met in fractionMets.keys():
        for rxn in met.reactions:
            if rxn.id.__contains__("_dielTransfer"):
                continue
            else:
                mainMet = met.id[len(fractionMets[met]):]
                coeff1 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(mainMet))
                coeff2 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(met.id))
                if not coeff1:
                    coeff1=0
                if not coeff2:
                    coeff2=0
                total = coeff1 + coeff2
                coeff1 = float(coeff1)/total
                coeff2 = float(coeff2)/total
                if cobra_model.reactions.has_id(met.id[0:len(met.id)-1]+"_dielTransfer"):
                    ub = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").upper_bound
                    lb = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").lower_bound
                    temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").remove_from_model()
                    temp.reactions.get_by_id(mainMet[0:len(mainMet)-1]+"_dielTransfer").remove_from_model()
                    Reac = Reaction(mainMet[0:len(mainMet)-1]+"_dielTransfer",name=mainMet+"_dielTransfer")
                    Reac.add_metabolites({temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"1"):-coeff2,temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"2"):coeff2,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"1"):-coeff1,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"2"):coeff1})
                    Reac.lower_bound=lb
                    Reac.upper_bound=ub
                    temp.add_reaction(Reac)
                    print Reac.reaction
                break
    ####ADD CONSTRAINTS TO MODEL####
    cobra_model = temp.copy()

    #Leaves - light
    cobra_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("CO2_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").upper_bound=0
    #Leaves - dark
    cobra_model.reactions.get_by_id("Sucrose_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("CO2_tx2").upper_bound=0

    #Set pG6P transporter to 0
    cobra_model.reactions.get_by_id("G6P_Pi_pc1").lower_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc1").upper_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc2").lower_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc2").upper_bound=0

    #Turn off PTOX
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").upper_bound=0

    #nitrate uptake constrain
    Nitrate_balance = Metabolite("Nitrate_bal_c", name = "Weights to balance nitrate uptake", compartment = "c1")
    cobra_model.reactions.get_by_id("Nitrate_ec1").add_metabolites({Nitrate_balance:-2})
    cobra_model.reactions.get_by_id("Nitrate_ec2").add_metabolites({Nitrate_balance:3})

    
    #Rubisco balance
    Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
    cobra_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:3})
    cobra_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})

    #generic ATPase and NADPH oxidase
    Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
    Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
    Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
    cobra_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
    cobra_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
    cobra_model.reactions.get_by_id("NADPHoxc_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxc_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx2").add_metabolites({Maintenance_constraint2:-3})

    #Plastid enolase was not detected in Arabidopsis mesophyll tissue
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").upper_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").upper_bound=0

    #Setting chloroplastic NADPH dehydrogenase to 0  ((Yamamoto et al., 2011)
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").upper_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").upper_bound=0

    #Set biomass to zero
    cobra_model.reactions.get_by_id("Biomass_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("Biomass_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("Biomass_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Biomass_tx2").upper_bound=0

    #ATP_ADP_Pi constrained to 0 because while there is evidence for its existance, it does not carry high flux
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").upper_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").upper_bound = 0

    #turn off chlorophyll a/b cycling for energy dissipation
    cobra_model.reactions.get_by_id("RXN_7674_p1").lower_bound = 0
    cobra_model.reactions.get_by_id("RXN_7674_p1").upper_bound = 0


    #turn off cytosolic ferric chelate reductase cycle for NADH dissipation
    cobra_model.reactions.get_by_id("FERRIC_CHELATE_REDUCTASE_RXN_c1").lower_bound = 0
    cobra_model.reactions.get_by_id("FERRIC_CHELATE_REDUCTASE_RXN_c1").upper_bound = 0

    #Adding a H_mc reaction to allow protons into mitochondria
    for i in range(1,3):
        rxn = Reaction("H_mc"+str(i))
        rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_c"+str(i)):-1,cobra_model.metabolites.get_by_id("PROTON_m"+str(i)):1})
        rxn.lower_bound=0
        rxn.upper_bound=1000
        cobra_model.add_reactions({rxn})

    return cobra_model


def fixModelCompatibilityIssueCobra015(model,sbml_file):
    '''
    This function can correct model imported in cobra v0.15.0
    '''

    import libsbml

    #remove EX reactions
    rxn2remove = list()
    for rxn in model.reactions:
        if rxn.id.startswith("EX_"):
            rxn2remove.append(rxn.id)
    for rxn in rxn2remove:
        rxn = model.reactions.get_by_id(rxn)
        rxn.remove_from_model()

    #remove boundary metabolites (metabolite with boundaryCondition=True)
    met2remove = list()
    for met in model.metabolites:
        if met.compartment == "t":
            met2remove.append(met.id)
    for met in met2remove:
        met = model.metabolites.get_by_id(met)
        met.remove_from_model()

    #add metbolite InChI and formula from an sbml file
    sbmlDoc = libsbml.readSBMLFromFile(sbml_file)
    InChiDict = dict()
    for x in sbmlDoc.getListOfAllElements():
        if type(x) == libsbml.Species:
            notes = x.getNotes().toXMLString()
            if "INCHI" in notes:
                inchi = notes.split("INCHI:")[1].split("<")[0].replace(" ","")
                if inchi == "NA" or inchi == "":
                    continue
                else:
                    InChiDict[x.id[2:]] = inchi

    for met in model.metabolites:
        if met.id in InChiDict.keys():
            met.notes["INCHI"] = InChiDict[met.id]
        else:
            met.notes["INCHI"] = ""

    ChargeDict = dict()
    for x in sbmlDoc.getListOfAllElements():
        if type(x) == libsbml.Species:
            notes = x.getNotes().toXMLString()
            if "CHARGE" in notes:
                charge = notes.split("CHARGE:")[1].split("<")[0].replace(" ","")
                if charge == "":
                    continue
                else:
                    ChargeDict[x.id[2:]] = float(charge)

    for met in model.metabolites:
        if met.id in ChargeDict.keys():
            met.charge = ChargeDict[met.id]

    ECDict = dict()
    for x in sbmlDoc.getListOfAllElements():
        if type(x) == libsbml.Reaction:
            EC = ""
            notes = x.getNotes().toXMLString()
            if "PROTEIN_CLASS" in notes:
                EC = notes.split("PROTEIN_CLASS:")[1].split("<")[0].replace(" ","")
                if EC.count(".")==2:
                    EC=EC+".-"
                elif EC.count(".")==1:
                    EC=EC+".-.-"
                ECDict[x.id[2:]] = EC

    PathDict = dict()
    for x in sbmlDoc.getListOfAllElements():
        if type(x) == libsbml.Reaction:
            temp = list()
            notes = x.getNotes().toXMLString()
            if "SUBSYSTEM" in notes:
                for i in range(1,len(notes.split("SUBSYSTEM:"))):
                    subsystem = notes.split("SUBSYSTEM:")[i].split("<")[0].replace(" ","")
                    if subsystem == "":
                        continue
                    else:
                        temp.append(subsystem)
                    PathDict[x.id[2:]] = temp
    for rxn in model.reactions:
        if rxn.id in PathDict.keys():
            rxn.notes["SUBSYSTEM"] = str(PathDict[rxn.id]).replace("[","").replace("]","").replace("\'","")
            rxn.notes["PROTEIN_CLASS"] = str(ECDict[rxn.id])
    return model

def estimateOutputFromNetCO2(model,netCO2uptake,Output_ID="diel_biomass",Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1",CO2in_ID="CO2_tx1",verbose=False):

    from cobra import flux_analysis
    # Initally constraint Vc flux to net CO2 uptake rate
    print(model.reactions.get_by_id("Photon_tx1").upper_bound)
    print(model.reactions.get_by_id("ATPase_tx1").upper_bound)
    model.reactions.get_by_id(Vc_ID).lower_bound = netCO2uptake
    model.reactions.get_by_id(Vc_ID).upper_bound = netCO2uptake
    print(model.reactions.get_by_id(Vc_ID).upper_bound)
    #perform pFBA
    flux_analysis.parsimonious.pfba(model)

    #unconstrain Vc
    model.reactions.get_by_id(Vc_ID).lower_bound = 0
    model.reactions.get_by_id(Vc_ID).upper_bound = 1000

    #set loop counter
    i=0

    #Use a while loop to increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
    while((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake > 0.00001 and i<100):
        i=i+1
        prev = model.reactions.get_by_id(Output_ID).x
        # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
        now = prev + (prev*((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake))
        model.reactions.get_by_id(Output_ID).lower_bound = now
        model.reactions.get_by_id(Output_ID).upper_bound = now
        print(now)
        flux_analysis.parsimonious.pfba(model)
        if verbose:
            print("----"+str(i)+"----")
            print("Vc flux ="+str(model.reactions.get_by_id(Vc_ID).x))
            print("net CO2 uptake ="+str(model.reactions.get_by_id(CO2in_ID).x))
            print("Target CO2 uptake ="+str(netCO2uptake))
            print("Before:"+str(prev))
            print("After:"+str(now))
            print("photon flux = "+str(model.reactions.Photon_tx1.x))
    return prev