
def is_liver_disease(pathdf, corpus='en_core_sci_lg', term_set='en_clinical', update=True, only_liv_biopsy=True):

    '''
    corpus: en_ner_bc5cdr_md, en_core_sci_md, en_core_sci_lg
    termset: en, en_clinical, en_clinical_sensitive

    '''
    
    fil_subset = pathdf.MRN_Type.isin(['MGH','BWH','NWH','FH','NSM'])
    df_path = pathdf[fil_subset].copy()
    
    
    if only_liv_biopsy:
        # check the column exists first:
        if 'is_liver_biopsy' in df_path.columns.tolist():
            fil_finaldx_trunc = df_path.is_liver_biopsy == True
            df_path = df_path[fil_finaldx_trunc]
        else:
            print('The flag *only_finaldx=True* was passed, however truncate_finaldx() has not been called. Aborting...')
            return None
        
#     filter_keywords = df_path['Report_Text'].str.contains('steato|balloon|baloon|ballon|inflam|hepatitis|hepatic|fibrosis|bridging|cirrhosis|aih', case=False, na=False)
#     df_path = df_path[filter_keywords]

    import spacy
    from negspacy.negation import Negex
    from negspacy.termsets import termset
    import numpy as np
    import pandas as pd
    import re
    #from spacy.pipeline import EntityRuler

    ts = termset(term_set)

    config={
        "neg_termset":{
            "pseudo_negations": ts.terms['pseudo_negations'] + ['not limited to', 'not excluded', 'needs to be ruled out', 'although not apparent'],
            "preceding_negations": ts.terms['preceding_negations'] + ['negative', 'insufficient', 'without evidence of', 'rather than', 'precludes'], #'grade 0'
            "following_negations": ts.terms['following_negations'] + ['negative', 'ruled out', 'less likely', 'is not', 'are not', 'does not', 'have not', 'was not', 'were not', 'absent', 'not present'], #unremarkable
            "termination": ts.terms['termination'] + ['note:', ';', ', negative', ',negative'] #'negative for', 'with'
        }
    }


    nlp_2 = spacy.load(corpus) 

    # ruler = EntityRuler(nlp_2, overwrite_ents=True)
    # patterns = [
    #     {"label": "ENTITY", "pattern": [{"LOWER": "chronic inflammation"}]}
    #         ]
    # ruler.add_patterns(patterns)

    nlp_2.add_pipe(
        "negex",
        config = config
    )

    num_reports = df_path.shape[0]
    steatosis_col = []
    ballooning_col = []
    inflammation_col = []
    lobular_inflammation_col = []
#     portal_inflammation_col = []
    zone3_inflammation_col = []
    lobular_hepatitis_col= []
    zone3_hepatitis_col = []
    
    fibrosis_col = []
    bridging_fibrosis_col = []
    sinusoidal_fibrosis_col = []
    portal_fibrosis_col = []
    periportal_fibrosis_col = []
    pericellular_fibrosis_col = []
    perivenular_fibrosis_col = []
    septal_fibrosis_col = []
    central_fibrosis_col = []
    zone3_fibrosis_col = []
    zone1_fibrosis_col = []
    centrilob_fibrosis_col = []
    hepatitis_col = []
    autoimmune_hepatitis_col = []
    mallory_col = []
    
    pbc_col = []
    cirrhosis_col = []
    steatohepatitis_col = []
    hepatitisa_col = []
    hepatitisb_col = []
    hepatitisc_col = []
    drug_hepatitis_col = []
    interface_hepatitis_col = []
    viral_hepatitis_col = []
    granulomatous_hepatitis_col = []
    hepatic_parenchyma_col = []
    
    hemochromatosis_col = []
    antitrypsin_col = []
    cholangitis_col = []
    wilsons_col = []
    drug_ind_liv_inj_col = []
    budd_chiari_col = []
    alcoholic_col = []
    carcinoma_col = []
    methotrexate_col = []
    ext_hematopoiesis_col = []
    
    nafld_col = []
    nash_col = []
    
    fibrosis_stage_4_col = []
    fibrosis_stage_6_col = []
    
    disease_list_col = []

    for i in range(0,num_reports):

        # extract path report for this entry
        disease_list = []
        report_text = df_path.iloc[i,:].Report_Text
        result_text = entities(report_text, nlp=nlp_2)

        steatosis = np.nan
        ballooning = np.nan
        inflammation = np.nan
        lobular_inflammation = np.nan
#         portal_inflammation = False
        zone3_inflammation = np.nan
        lobular_hepatitis = np.nan
        zone3_hepatitis = np.nan
        
        fibrosis = np.nan
        bridging_fibrosis = np.nan
        sinusoidal_fibrosis = np.nan
        portal_fibrosis = np.nan
        periportal_fibrosis = np.nan
        pericellular_fibrosis = np.nan
        perivenular_fibrosis = np.nan
        septal_fibrosis = np.nan
        central_fibrosis = np.nan
        zone3_fibrosis = np.nan
        zone1_fibrosis = np.nan
        centrilob_fibrosis = np.nan
        
        hepatitis = np.nan
        autoimmune_hepatitis = np.nan
        mallory = np.nan
        
        pbc = np.nan
        cirrhosis = np.nan
        steatohepatitis = np.nan
        hepatitisa = np.nan
        hepatitisb = np.nan
        hepatitisc = np.nan
        drug_hepatitis = np.nan
        interface_hepatitis = np.nan
        viral_hepatitis = np.nan
        granulomatous_hepatitis = np.nan

        hepatic_parenchyma = np.nan
        hemochromatosis = np.nan
        antitrypsin = np.nan
        cholangitis = np.nan
        wilsons = np.nan
        drug_ind_liv_inj = np.nan
        budd_chiari = np.nan
        alcoholic = np.nan
        carcinoma = np.nan
        methotrexate = np.nan
        ext_hematopoiesis = np.nan
        
        nafld = np.nan
        nash = np.nan
        
        fibrosis_stage_4 = np.nan
        fibrosis_stage_6 = np.nan
        
        other_liv_diseases = np.nan
        
        steatosis_lt5 = sum([1 for ent_text in result_text.split('\n') if '<5%' in ent_text and 'steatosis' in ent_text])==0
        
        fib_stage = -1
        fib_ref = -1
        
        for x in result_text.split('\n'):
            
            is_disease = False
            
            if 'steatosis' in x and (np.isnan(steatosis) or 'True' in x): #and steatosis_lt5
                if 'True' in x: steatosis, is_disease = True, True
                else: steatosis = False
            if ('balloon' in x or 'baloon' in x or 'ballon' in x) and (np.isnan(ballooning) or 'True' in x):
                if 'True' in x: ballooning, is_disease = True, True
                else: ballooning = False
            if 'inflammation' in x and (np.isnan(inflammation) or 'True' in x):
                if 'True' in x: inflammation, is_disease = True, True
                else: inflammation = False
            if 'lobular' in x and ('inflammation' in x or 'activity' in x or 'infiltrate' in x) and (np.isnan(lobular_inflammation) or 'True' in x):
                if 'True' in x: lobular_inflammation, is_disease = True, True
                else: lobular_inflammation = False
            if 'zone-3' in x and 'inflammation' in x and (np.isnan(zone3_inflammation) or 'True' in x):
                if 'True' in x: zone3_inflammation, is_disease = True, True
                else: zone3_inflammation = False
            if ('lobular' in x and 'hepatitis' in x and not 'steato' in x) and (np.isnan(lobular_hepatitis) or 'True' in x):
                if 'True' in x: lobular_hepatitis, is_disease = True, True
                else: lobular_hepatitis = False
            if ('zone-3' in x and 'hepatitis' in x and not 'steato' in x) and (np.isnan(zone3_hepatitis) or 'True' in x):
                if 'True' in x: zone3_hepatitis, is_disease = True, True
                else: zone3_hepatitis = False
            
            if 'fibrosis' in x and (np.isnan(fibrosis) or 'True' in x):
                if 'True' in x: fibrosis, is_disease = True, True
                else: fibrosis = False
            if ('bridging' in x and not 'bridging-necrosis' in x) and (np.isnan(bridging_fibrosis) or 'True' in x):
                if 'True' in x: bridging_fibrosis, is_disease = True, True
                else: bridging_fibrosis = False
            if ('fibrosis' in x and 'sinusoidal' in x) and (np.isnan(sinusoidal_fibrosis) or 'True' in x):
                if 'True' in x: sinusoidal_fibrosis, is_disease = True, True
                else: sinusoidal_fibrosis = False
            if ('fibrosis' in x and 'portal' in x and not 'periportal' in x) and (np.isnan(portal_fibrosis) or 'True' in x):
                if 'True' in x: portal_fibrosis, is_disease = True, True
                else: portal_fibrosis = False
            if ('fibrosis' in x and 'periportal' in x) and (np.isnan(periportal_fibrosis) or 'True' in x):
                if 'True' in x: periportal_fibrosis, is_disease = True, True
                else: periportal_fibrosis = False
            if ('fibrosis' in x and 'pericellular' in x) and (np.isnan(pericellular_fibrosis) or 'True' in x):
                if 'True' in x: pericellular_fibrosis, is_disease = True, True
                else: pericellular_fibrosis = False
            if ('fibrosis' in x and 'perivenular' in x) and (np.isnan(perivenular_fibrosis) or 'True' in x):
                if 'True' in x: perivenular_fibrosis, is_disease = True, True
                else: perivenular_fibrosis = False
            if ('fibrosis' in x and 'septal' in x) and (np.isnan(septal_fibrosis) or 'True' in x):
                if 'True' in x: septal_fibrosis, is_disease = True, True
                else: septal_fibrosis = False
            if ('fibrosis' in x and 'central' in x) and (np.isnan(central_fibrosis) or 'True' in x):
                if 'True' in x: central_fibrosis, is_disease = True, True
                else: central_fibrosis = False
            if ('fibrosis' in x and 'zone-3' in x) and (np.isnan(zone3_fibrosis) or 'True' in x):
                if 'True' in x: zone3_fibrosis, is_disease = True, True
                else: zone3_fibrosis = False
            if ('fibrosis' in x and 'zone-1' in x) and (np.isnan(zone1_fibrosis) or 'True' in x):
                if 'True' in x: zone1_fibrosis, is_disease = True, True
                else: zone1_fibrosis = False
            if ('fibrosis' in x and 'centrilobular' in x) and (np.isnan(centrilob_fibrosis) or 'True' in x):
                if 'True' in x: centrilob_fibrosis, is_disease = True, True
                else: centrilob_fibrosis = False
                

            if ('hepatitis' in x and not 'steato' in x) and (np.isnan(hepatitis) or 'True' in x):
                if 'True' in x: hepatitis, is_disease = True, True
                else: hepatitis = False
            if ('autoimmune hepatitis' in x or bool(re.search(r'\baih\b', x)) or 'auto-immune hepatitis' in x) and (np.isnan(autoimmune_hepatitis) or 'True' in x):
                if 'True' in x: autoimmune_hepatitis, is_disease = True, True
                else: autoimmune_hepatitis = False
            if 'mallory' in x and (np.isnan(mallory) or 'True' in x):
                if 'True' in x: mallory, is_disease = True, True
                else: mallory = False
            
            if ('biliary' in x and 'cirrhosis' in x) and (np.isnan(pbc) or 'True' in x):
                if 'True' in x: pbc, is_disease = True, True
                else: pbc = False
            if ('cirrhosis' in x and 'biliary' not in x) and (np.isnan(cirrhosis) or 'True' in x):
                if 'True' in x: cirrhosis, is_disease = True, True
                else: cirrhosis = False
            if 'steatohepatitis' in x and (np.isnan(steatohepatitis) or 'True' in x):
                if 'True' in x: steatohepatitis, is_disease = True, True
                else: steatohepatitis = False
                    
#             if 'portal' in x and 'inflammation' in x and 'True' in x:
#                 portal_inflammation = True
#             if ('chronic' in x and 'hepatitis' in x) and (np.isnan(chronic_hepatitis) or 'True' in x):
#                 if 'True' in x: chronic_hepatitis, is_disease = True, True
#                 else: chronic_hepatitis = False

            if (bool(re.search(r'\bhepatitis a\b', x)) or bool(re.search(r'\bhepatitis-a\b', x)) or bool(re.search(r'\bhep a\b', x))) and (np.isnan(hepatitisa) or 'True' in x):
                if 'True' in x: hepatitisa, is_disease = True, True
                else: hepatitisa = False
            if (bool(re.search(r'\bhepatitis b\b', x)) or bool(re.search(r'\bhepatitis-b\b', x)) 
                    or bool(re.search(r'\bhep b\b', x)) or bool(re.search(r'\bhbv\b', x)) 
                    or '_hbv' in x) and (np.isnan(hepatitisb) or 'True' in x):
                if 'True' in x: hepatitisb, is_disease = True, True
                else: hepatitisb = False
            # Need to add bool(re.search(r'\bchronichepatitis c\b', x)) bool(re.search(r'\bhepatitis c\b', x))
            # bool(re.search(r'\bhepatitis c\b', x))
            if (bool(re.search(r'\bhepatitis c\b', x)) or bool(re.search(r'\bhepatitis-c\b', x)) 
                    or bool(re.search(r'\bhcv\b', x)) or bool(re.search(r'\bhep c\b', x))
                    or 'ishak' in x or '_hcv' in x) and (np.isnan(hepatitisc) or 'True' in x):
                if 'True' in x: hepatitisc, is_disease = True, True
                else: hepatitisc = False
            if ('drug' in x and 'hepatitis' in x) and (np.isnan(drug_hepatitis) or 'True' in x):
                if 'True' in x: drug_hepatitis, is_disease = True, True
                else: drug_hepatitis = False
            if ('interface' in x and 'hepatitis' in x) and (np.isnan(interface_hepatitis) or 'True' in x):
                if 'True' in x: interface_hepatitis, is_disease = True, True
                else: inflammation = False
            if ('viral' in x and 'hepatitis' in x) and (np.isnan(viral_hepatitis) or 'True' in x):
                if 'True' in x: viral_hepatitis, is_disease = True, True
                else: viral_hepatitis = False
            if ('granulomatous' in x and 'hepatitis' in x) and (np.isnan(granulomatous_hepatitis) or 'True' in x):
                if 'True' in x: granulomatous_hepatitis, is_disease = True, True
                else: granulomatous_hepatitis = False
            
            if ('hepatic' in x and 'parenchyma' in x) and (np.isnan(hepatic_parenchyma) or 'True' in x):
                if 'True' in x: hepatic_parenchyma, is_disease = True, True
                else: hepatic_parenchyma = False
            if 'hemochromatosis' in x and (np.isnan(hemochromatosis) or 'True' in x):
                if 'True' in x: hemochromatosis, is_disease = True, True
                else: hemochromatosis = False
            if 'antitrypsin' in x and (np.isnan(antitrypsin) or 'True' in x):
                if 'True' in x: antitrypsin, is_disease = True, True
                else: antitrypsin = False
            if 'cholangitis' in x and (np.isnan(cholangitis) or 'True' in x):
                if 'True' in x: cholangitis, is_disease = True, True
                else: cholangitis = False
            if "wilson's" in x and (np.isnan(wilsons) or 'True' in x):
                if 'True' in x: wilsons, is_disease = True, True
                else: wilsons = False
            if ('drug-induced' in x or bool(re.search(r'\bdili\b', x)) or 'drug-related' in x) and (np.isnan(drug_ind_liv_inj) or 'True' in x):
                if 'True' in x: drug_ind_liv_inj, is_disease = True, True
                else: drug_ind_liv_inj = False
            if 'budd-chiari' in x and (np.isnan(budd_chiari) or 'True' in x):
                if 'True' in x: budd_chiari, is_disease = True, True
                else: budd_chiari = False
            if (bool(re.search(r'\balcoholic\b', x)) or 'ethanol' in x) and (np.isnan(alcoholic) or 'True' in x):
                if 'True' in x: alcoholic, is_disease = True, True
                else: alcoholic = False
            if ('metastatic' in x or 'metastases' in x or 'metastasis' in x or 'carcinoma' in x or 'lymphoma' in x
                       or bool(re.search(r'\bhcc\b', x)) or 'malign' in x or 'cancer' in x or 'carcinoid' in x
                       or 'angiosarcoma' in x) and (np.isnan(carcinoma) or 'True' in x):
                if 'True' in x: carcinoma, is_disease = True, True
                else: carcinoma = False
                    
            'metastatic', 'metastases', 'metastasis', 'carcinoma', 'lymphoma', 'hcc', 'malign', 'cancer', 'carcinoid', 'angiosarcoma'
                    
            if 'methotrexate' in x and (np.isnan(methotrexate) or 'True' in x):
                if 'True' in x: methotrexate, is_disease = True, True
                else: methotrexate = False
            if 'extramedullary-hematopoiesis' in x and (np.isnan(ext_hematopoiesis) or 'True' in x):
                if 'True' in x: ext_hematopoiesis, is_disease = True, True
                else: ext_hematopoiesis = False
            
            if (bool(re.search(r'\bnafld\b', x)) or 'nonalcoholic fatty liver disease' in x) and (np.isnan(nafld) or 'True' in x):
                if 'True' in x: nafld, is_disease = True, True
                else: nafld = False
            if (bool(re.search(r'\bnash\b', x)) or 'nonalcoholic steatohepatitis' in x) and (np.isnan(nash) or 'True' in x):
                if 'True' in x: nash, is_disease = True, True
                else: nash = False
                    
            if ('fibrosis' in x or 'bridging' in x or 'cirrhosis' in x) and ' stage:' in x:
                try:
                    if 'True' in x:
                        fib_stage = float(x[-12:-9])
                        fib_ref = float(x[-8:-5])
                    elif 'False' in x:
                        fib_stage = float(x[-13:-10])
                        fib_ref = float(x[-9:-6])

                    if fib_ref<4.0:
                        fib_ref = 4.0
                        
                    if fib_ref==4 and np.isnan(fibrosis_stage_4):
                        fibrosis_stage_4 = fib_stage
                    elif fib_ref==6 and np.isnan(fibrosis_stage_6):
                        fibrosis_stage_6 = fib_stage
                    
                    is_disease = True
                    
                except:
                    pass
                
            if is_disease:
                disease_list.append(x)
                
        
        steatosis_col.append(steatosis)
        ballooning_col.append(ballooning)
        inflammation_col.append(inflammation)
        lobular_inflammation_col.append(lobular_inflammation)
#         portal_inflammation_col.append(portal_inflammation)
        zone3_inflammation_col.append(zone3_inflammation)
        lobular_hepatitis_col.append(lobular_hepatitis)
        zone3_hepatitis_col.append(zone3_hepatitis)
        fibrosis_col.append(fibrosis)
        
        
        bridging_fibrosis_col.append(bridging_fibrosis)
        sinusoidal_fibrosis_col.append(sinusoidal_fibrosis)
        portal_fibrosis_col.append(portal_fibrosis)
        periportal_fibrosis_col.append(periportal_fibrosis)
        pericellular_fibrosis_col.append(pericellular_fibrosis)
        perivenular_fibrosis_col.append(perivenular_fibrosis)
        septal_fibrosis_col.append(septal_fibrosis)
        central_fibrosis_col.append(central_fibrosis)
        zone3_fibrosis_col.append(zone3_fibrosis)
        zone1_fibrosis_col.append(zone1_fibrosis)
        centrilob_fibrosis_col.append(centrilob_fibrosis)
        hepatitis_col.append(hepatitis)
        autoimmune_hepatitis_col.append(autoimmune_hepatitis)
        mallory_col.append(mallory)
        
        pbc_col.append(pbc)
        cirrhosis_col.append(cirrhosis)
        steatohepatitis_col.append(steatohepatitis)
        hepatitisa_col.append(hepatitisa)
        hepatitisb_col.append(hepatitisb)
        hepatitisc_col.append(hepatitisc)
        drug_hepatitis_col.append(drug_hepatitis)
        interface_hepatitis_col.append(interface_hepatitis)
        viral_hepatitis_col.append(viral_hepatitis)
        granulomatous_hepatitis_col.append(granulomatous_hepatitis)
        hepatic_parenchyma_col.append(hepatic_parenchyma)
        
        hemochromatosis_col.append(hemochromatosis)
        antitrypsin_col.append(antitrypsin)
        cholangitis_col.append(cholangitis)
        wilsons_col.append(wilsons)
        drug_ind_liv_inj_col.append(drug_ind_liv_inj)
        budd_chiari_col.append(budd_chiari)
        alcoholic_col.append(alcoholic)
        carcinoma_col.append(carcinoma)
        methotrexate_col.append(methotrexate)
        ext_hematopoiesis_col.append(ext_hematopoiesis)
        
        nafld_col.append(nafld)
        nash_col.append(nash)
        
        fibrosis_stage_4_col.append(fibrosis_stage_4)
        fibrosis_stage_6_col.append(fibrosis_stage_6)
        
        disease_list_col.append(disease_list)
        
        
    df_path['steatosis'] = steatosis_col
    df_path['ballooning'] = ballooning_col
    df_path['inflammation'] = inflammation_col
    df_path['lobular_inflammation'] = lobular_inflammation_col
#     df_path['portal_inflammation'] = portal_inflammation_col
    df_path['zone3_inflammation'] = zone3_inflammation_col
    df_path['lobular_hepatitis'] = lobular_hepatitis_col
    df_path['zone3_hepatitis'] = zone3_hepatitis_col
    
    df_path['fibrosis'] = fibrosis_col
    df_path['bridging_fibrosis'] = bridging_fibrosis_col
    df_path['sinusoidal_fibrosis'] = sinusoidal_fibrosis_col
    df_path['portal_fibrosis'] = portal_fibrosis_col
    df_path['periportal_fibrosis'] = periportal_fibrosis_col
    df_path['pericellular_fibrosis'] = pericellular_fibrosis_col
    df_path['perivenular_fibrosis'] = perivenular_fibrosis_col
    df_path['septal_fibrosis'] = septal_fibrosis_col
    df_path['central_fibrosis'] = central_fibrosis_col
    df_path['zone3_fibrosis'] = zone3_fibrosis_col
    df_path['zone1_fibrosis'] = zone1_fibrosis_col
    df_path['centrilob_fibrosis'] = centrilob_fibrosis_col
    df_path['hepatitis'] = hepatitis_col
    df_path['autoimmune_hepatitis'] = autoimmune_hepatitis_col
    df_path['mallory'] = mallory_col
    
    df_path['pbc'] = pbc_col
    df_path['cirrhosis'] = cirrhosis_col
    df_path['steatohepatitis'] = steatohepatitis_col
    df_path['hepatitisa'] = hepatitisa_col
    df_path['hepatitisb'] = hepatitisb_col
    df_path['hepatitisc'] = hepatitisc_col
    df_path['drug_hepatitis'] = drug_hepatitis_col
    df_path['interface_hepatitis'] = interface_hepatitis_col
    df_path['viral_hepatitis'] = viral_hepatitis_col
    df_path['granulomatous_hepatitis'] = granulomatous_hepatitis_col
    df_path['hepatic_parenchyma'] = hepatic_parenchyma_col
    
    df_path['hemochromatosis'] = hemochromatosis_col
    df_path['antitrypsin'] = antitrypsin_col
    df_path['cholangitis'] = cholangitis_col
    df_path['wilsons'] = wilsons_col
    df_path['drug_ind_liv_inj'] = drug_ind_liv_inj_col
    df_path['budd_chiari'] = budd_chiari_col
    df_path['alcoholic'] = alcoholic_col
    df_path['carcinoma'] = carcinoma_col
    df_path['methotrexate'] = methotrexate_col
    df_path['ext_hematopoiesis'] = ext_hematopoiesis_col
    
    df_path['nafld'] = nafld_col
    df_path['nash'] = nash_col
    
    df_path['fibrosis_stage_4'] = fibrosis_stage_4_col
    df_path['fibrosis_stage_6'] = fibrosis_stage_6_col
        
    df_path['disease_list'] = disease_list_col
   
    if update:
        # re-merge with original data
        print('Updating input path dataframe')
        pathdf['steatosis'] = np.nan
        pathdf['ballooning'] = np.nan
        pathdf['inflammation'] = np.nan
        pathdf['lobular_inflammation'] = np.nan
#         pathdf['portal_inflammation'] = np.nan
        pathdf['zone3_inflammation'] = np.nan
        pathdf['lobular_hepatitis'] = np.nan
        pathdf['zone3_hepatitis'] = np.nan
        
        pathdf['fibrosis'] = np.nan
        pathdf['bridging_fibrosis'] = np.nan
        pathdf['sinusoidal_fibrosis'] = np.nan
        pathdf['portal_fibrosis'] = np.nan
        pathdf['periportal_fibrosis'] = np.nan
        pathdf['pericellular_fibrosis'] = np.nan
        pathdf['perivenular_fibrosis'] = np.nan
        pathdf['septal_fibrosis'] = np.nan
        pathdf['central_fibrosis'] = np.nan
        pathdf['zone3_fibrosis'] = np.nan
        pathdf['zone1_fibrosis'] = np.nan
        pathdf['centrilob_fibrosis'] = np.nan
        pathdf['hepatitis'] = np.nan
        pathdf['autoimmune_hepatitis'] = np.nan
        pathdf['mallory'] = np.nan
        
        pathdf['pbc'] = np.nan
        pathdf['cirrhosis'] = np.nan
        pathdf['steatohepatitis'] = np.nan
        pathdf['hepatitisa'] = np.nan
        pathdf['hepatitisb'] = np.nan
        pathdf['hepatitisc'] = np.nan
        pathdf['drug_hepatitis'] = np.nan
        pathdf['interface_hepatitis'] = np.nan
        pathdf['viral_hepatitis'] = np.nan
        pathdf['granulomatous_hepatitis'] = np.nan
        pathdf['hepatic_parenchyma'] = np.nan
        
        pathdf['hemochromatosis'] = np.nan
        pathdf['antitrypsin'] = np.nan
        pathdf['cholangitis'] = np.nan
        pathdf['wilsons'] = np.nan
        pathdf['drug_ind_liv_inj'] = np.nan
        pathdf['budd_chiari'] = np.nan
        pathdf['alcoholic'] = np.nan
        pathdf['carcinoma'] = np.nan
        pathdf['methotrexate'] = np.nan
        pathdf['ext_hematopoiesis'] = np.nan
        
        pathdf['nafld'] = np.nan
        pathdf['nash'] = np.nan
        
        pathdf['fibrosis_stage_4'] = np.nan
        pathdf['fibrosis_stage_6'] = np.nan
            
        pathdf['disease_list'] = np.nan
        pathdf.update(df_path)
        return_df = pathdf.copy()
    else:
        # return this mgh path only file
        #print('Returning MGH, BWH only entries with truncated path reports')
        return_df = df_path
        

    return return_df


def entities(text, nlp):
    
    import re
    
    text = (text
            .replace('Hep. A', 'Hep A').replace('Hep. B', 'Hep B').replace('Hep. C', 'Hep C')
            .replace('/IV', '/4').replace('/VI', '/6')
            .replace('PBC', 'PBC (primary biliary cirrhosis)')
            .replace('PSC', 'PSC (primary sclerosing cholangitis)')
           )
    
    text = text.lower().replace(' bridging.', ' active-bridging.') #.replace(';',' ')
    
    entity_result = ''
    
    for line in text.split('.'):
        
        line = " ".join(line.split())
        line = line.strip()
        line = (line
                .replace('+/-', ',')
                .replace(' no ', ' , no ')
                .replace('no present', 'not present')
                .replace('is not present', 'not present').replace('not present', ' is not present')
                .replace(' minimal ', ' ,minimal ')
                .replace('noted in the', ' in ')
                .replace('is noted in', ' in ')
                .replace(' as well as', ', ')
                .replace('may not', 'will not')
                .replace('neither', 'no').replace('nor', 'no')
                .replace('very', '')
                .replace('mildly', 'mild').replace('mildl', 'mild')
                .replace('non alcoholic', 'nonalcoholic')
                .replace('non-alcoholic', 'nonalcoholic')
                .replace('steato-hepatitis', 'steatohepatitis')
                .replace('steato hepatitis', 'steatohepatitis')
                .replace('nonalcoholic steatohepatitis', 'nonalcoholic-steatohepatitis')
                .replace('inflammatory infiltrate', 'inflammatory-infiltrate')
                .replace('necroinflammatory', 'necroinflammatory (inflammation)')
                .replace('centric inflammation', 'centric-inflammation')
                .replace('inflammatory', 'inflammation')
                .replace('inflamed', 'inflammation')
                .replace('severely', 'severe')
                .replace('moderately', 'moderate')
                .replace('moderate ', 'moderate_').replace('(moderate) ', 'moderate_').replace('(moderate)', 'moderate_')
                .replace('mild to moderate', 'mild&moderate')
                .replace('moderate to severe', 'moderate&severe')
                .replace('mild to severe', 'mild&severe')
                .replace('mild active', 'mild-active')
                .replace('mild chronic', 'mild-chronic')
                .replace('mild ', 'mild_').replace('(mild) ', 'mild_').replace('(mild)', 'mild_')
                .replace('severe ', 'severe_').replace('(severe) ', 'severe_').replace('(severe)', 'severe_')
                .replace('minimal ', 'minimal_')
                .replace('chronic ', '')
                .replace('focal ', 'focal_')
                .replace(' areas', ' area')
                .replace('-area ', ' area ')
#                 .replace('portal area', 'portalarea')
                .replace(' tracts', ' tract')
                .replace('-tracts', ' tract')
#                 .replace('portal tract', 'portaltract')
#                 .replace('centrilobular', 'centri-lobular')
                .replace('peri-portal', 'periportal')
                .replace('periportal and lobular', 'portal&lobular')
                .replace('portal and lobular', 'portal&lobular')
                .replace('portal or lobular', 'portal&lobular')
                .replace('lobular and portal', 'lobular&portal')
                .replace('lobular or portal', 'lobular&portal')
                .replace('portal&lobular', 'lobular&portal')
                .replace('lobular inflammat', 'lobular-inflammat')
                .replace('mixed ', 'mixed-')
                .replace('kupffer cell', 'kupffer-cell')
                .replace('zone 3', 'zone-3').replace('zone3', 'zone-3')
                .replace('zone 1', 'zone-1').replace('zone1', 'zone-1')
                .replace('hepatic plate', 'hepatic-plate')
                .replace('hepatic parenchyma', 'hepatic-parenchyma')
#                 .replace('hepatic ', 'hepatitis ')
                .replace('steatotic', 'steatosis')
#                 .replace('microvesicular ', 'microvesicular-')
#                 .replace('macrovesicular ', 'macrovesicular-')
                .replace('< 5%', '<5%')
                .replace('non classical', 'nonclassical')
                .replace('non-classical', 'nonclassical')
                .replace('ballooning degeneration', 'ballooning_degeneration')
                .replace('hepatocyte ballooning', 'hepatocyte-ballooning')
                .replace('hepatocytic ballooning', 'hepatocytic-ballooning')
                .replace('heptocellular ballooning', 'heptocellular-ballooning')
                .replace('ballooned', 'ballooning')
                .replace('ballooning', 'baloning')
                .replace('portal to portal', 'portal-portal')
                .replace('central to central', 'central-central')
                .replace('portal to central', 'portal-central')
                .replace('bridging necrosis', 'bridging-necrosis')
                .replace('fibrosis bridging', 'fibrosis-bridging')
                .replace('fibrous bridging', 'fibrosis-bridging')
                .replace('bridging fibrosis', 'bridging-fibrosis')
                .replace('sinusoidal fibrosis', 'sinusoidal-fibrosis')
                .replace('portal fibrosis', 'portal-fibrosis')
                .replace('pericellular fibrosis', 'pericellular-fibrosis')
                .replace('perivenular fibrosis', 'perivenular-fibrosis')
                .replace('centrilobular fibrosis', 'centrilobular-fibrosis')
#                 .replace('septal fibrosis', 'septal-fibrosis')
#                 .replace('fibrous septa', 'septal-fibrosis')
#                 .replace('portal septa', 'portal-septa')
                .replace('central fibrosis', 'central-fibrosis')
                .replace('ductal fibrosis', 'ductal-fibrosis')
                .replace('portal bridging', 'portal-bridging')
                .replace('portal expansion', 'portal-expansion')
                .replace('central bridging', 'central-bridging')
                .replace(' bridging ', ' active-bridging ')
                .replace(' bridging;', ' active-bridging;')
                .replace(' bridging,', ' active-bridging,')
                .replace('(p-p)','')
                .replace(' hep ', ' hepatitis ')
                .replace('a1at', 'alpha-1-antitrypsin')
#                 .replace('alpha-1 antitrypsin', 'alpha-1-antitrypsin')
#                 .replace('alpha 1 antitrypsin', 'alpha-1-antitrypsin')
#                 .replace('alpha - 1 antitrypsin', 'alpha-1-antitrypsin')
#                 .replace('alpha 1-antitrypsin', 'alpha-1-antitrypsin')
                .replace('wilsons disease', "wilson's disease")
                .replace('wilson disease', "wilson's disease")
#                 .replace('drug induced liver injury', 'drug-induced-liver-injury')
#                 .replace('drug-induced liver injury', 'drug-induced-liver-injury')
#                 .replace('drug induced cholestatic liver injury', 'drug-induced-liver-injury')
#                 .replace('drug induced injury', 'drug-induced-liver-injury')
#                 .replace('drug induced hepatocellular injury', 'drug-induced-liver-injury')
#                 .replace('drug induced lesion', 'drug-induced-liver-injury')
                .replace('drug induced', 'drug-induced').replace('drug related', 'drug-related')
                .replace('budd chiari', 'budd-chiari')
                .replace('hepatocellular carcinoma', 'hepatocellular-carcinoma')
                .replace('cirrhoses', 'cirrhosis')
                .replace('hepatitis-a', 'hepatitis a')
                .replace('hepatitis-b', 'hepatitis b')
                .replace('hepatitis-c', 'hepatitis c')
                .replace('steato-hepatitis', 'steatohepatitis').replace('steato hepatitis', 'steatohepatitis')
                .replace('centri-lobular', 'centrilobular')
                .replace('centrilobular necrosis', 'centrilobular-necrosis')
                .replace('centrilobular hepatic necrosis', 'centrilobular-hepatic-necrosis')
                .replace('droplet steatosis', 'droplet-steatosis')
#                 .replace('_inflammation', ' iniflammation')
                .replace('_lobular', ' lobular').replace('_portal', ' portal')
                .replace('_fibrosis', ' fibrosis').replace('_inflammation', ' inflammation')
                .replace('early cirrhotic', 'early cirrhotic (cirrhosis)')
                .replace('primary biliary cirrhosis', 'primary_biliary_cirrhosis')
                .replace('biliary cirrhosis', 'biliary_cirrhosis')
                .replace('billiary cirrhosis', 'biliary_cirrhosis')
                .replace('AMA negative', 'AMA-negative')
                .replace('nash- ', 'nash ')
                .replace('extramedullary hematopoiesis', 'extramedullary-hematopoiesis')
                
               )
        
#         if 'fibrosis' in line:
#             line = (line
#                     .replace('pericellular ', 'pericellular-fibrosis ').replace('pericellular,', 'pericellular-fibrosis ')
#                     .replace('sinusoidal ', 'sinusoidal-fibrosis ').replace('sinusoidal,', 'sinusoidal-fibrosis ')
# #                     .replace(' perisinusoidal ', ' perisinusoidal-fibrosis ').replace(' perisinusoidal,', ' perisinusoidal-fibrosis ')
#                     .replace('periportal ', 'periportal-fibrosis ').replace('periportal, ', 'periportal-fibrosis, ')
#                     .replace('central ', 'central-fibrosis ').replace('central,', 'central-fibrosis ')
#                     .replace('septal ', 'septal-fibrosis ').replace('septal,', 'septal-fibrosis ')
#                     .replace('ductal ', 'ductal-fibrosis ').replace('ductal,', 'ductal-fibrosis ')
#                     .replace('perivenular ', 'perivenular-fibrosis ').replace('perivenular,', 'perivenular-fibrosis ')
#                    )
        if 'bridging' in line:
            line = (line
                    .replace('portal-portal-fibrosis', 'portal-portal-fibrosis-bridging')
                    .replace('portal-portal ', 'portal-portal-bridging')
                    .replace('portal-central-fibrosis', 'portal-central-fibrosis-bridging')
                    .replace('portal-central ', 'portal-central-bridging')
                    .replace('central-central-fibrosis', 'central-central-fibrosis-bridging')
                    .replace('central-central ', 'central-central-bridging')
                   )
            
        
        #global doc, e
   
        doc = nlp(line)
    
#         print(line)
    
        for e in doc.ents:
            
            e_text = e.text
            e_text = re.sub(' +', ' ', e_text)
            e_bool = e._.negex
            
            # Replace negation words in the entity and adjust sentiment
            if e_text.startswith(('no ', 'non-', 'non ')):
                to_match = ['^no ', '^non-', '^non ']
                e_text = re.sub('|'.join(to_match), '', e_text)
                e_bool = not e_bool
            
#             if 'lobular' in e_text or 'portal' in e_text:
#                 e_text = re.sub('(mild|moderate|mixed)', '', e_text).strip()

            inf_count = line.count('inflammation')
            
            lobu_inf1 = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?inflammation)\b', line))
            lobu_inf2 = bool(re.search(r'\b(?:lobules\W+(?:\w+\W+){0,6}?inflammation|inflammation\W+(?:\w+\W+){0,3}?lobules)\b', line))
            lobu_inf3 = bool(re.search(r'\b(?:lobule\W+(?:\w+\W+){0,6}?inflammation|inflammation\W+(?:\w+\W+){0,3}?lobule)\b', line))
            lobu_inf4 = bool(re.search(r'\b(?:centrilobular\W+(?:\w+\W+){0,3}?inflammation|inflammation\W+(?:\w+\W+){0,6}?centrilobular)\b', line))
            
            lobu_inf = (lobu_inf1 or lobu_inf2 or lobu_inf3 or lobu_inf4)
            
#             lobu_inf_2 = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?necroinflammatory)\b', line))
#             lobu_inf = lobu_inf_1 or lobu_inf_2
            
            lobu_infil1 = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,4}?infiltrate)\b', line))
            lobu_infil2 = bool(re.search(r'\b(?:infiltrate\W+(?:\w+\W+){0,6}?lobules)\b', line))
            lobu_infil3 = bool(re.search(r'\b(?:infiltrate\W+(?:\w+\W+){0,6}?lobular)\b', line))
            lobu_infil4 = bool(re.search(r'\b(?:infiltrate\W+(?:\w+\W+){0,3}?centrilobular)\b', line))
            
            lobu_infil = (lobu_infil1 or lobu_infil2 or lobu_infil3 or lobu_infil4)

            
            zone3_inf = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,2}?inflammation|inflammation\W+(?:\w+\W+){0,1}?zone-3)\b', line))
            
            # port_inf = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?inflammation|inflammation\W+(?:\w+\W+){0,1}?portal)\b', line))
            port_inf = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,4}?inflammation)\b', line))
            
#             port_tract_inf = bool(re.search(r'\b(?:portaltract\W+(?:\w+\W+){0,5}?inflammation)\b', line))
            
            lob_port_inf = bool(re.search(r'\b(?:lobular&portal\W+(?:\w+\W+){0,3}?inflammation)\b', line))

            lobu_act = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?activity)\b', line))
            port_act = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?activity)\b', line))

            lobu_hep = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,2}?hepatitis)\b', line))
            
            zone3_hep = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,3}?hepatitis)\b', line))
            
            lobu_dis = bool(re.search(r'\b(?:lobular\W+(?:\w+\W+){0,3}?disarray)\b', line)) 
            
            zone3_fib = bool(re.search(r'\b(?:zone-3\W+(?:\w+\W+){0,3}?fibrosis|fibrosis\W+(?:\w+\W+){0,3}?zone-3)\b', line))
            
            zone1_fib = bool(re.search(r'\b(?:zone-1\W+(?:\w+\W+){0,3}?fibrosis|fibrosis\W+(?:\w+\W+){0,3}?zone-1)\b', line))
            
            cent_scar_1 = bool(re.search(r'\b(?:central\W+(?:\w+\W+){0,2}?scarring|scarring\W+(?:\w+\W+){0,4}?central)\b', line))
            cent_scar_2 = bool(re.search(r'\b(?:centrilobular\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,4}?centrilobular)\b', line))
            cent_scar_3 = bool(re.search(r'\b(?:pericentral\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,4}?pericentral)\b', line))
            cent_scar = cent_scar_1 or cent_scar_3 or cent_scar_3
            
            port_scar_1 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?portal)\b', line))
            port_scar_2 = bool(re.search(r'\b(?:portal tract\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?portal tract)\b', line))
            port_scar_3 = bool(re.search(r'\b(?:portal area\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?portal area)\b', line))
            port_scar_4 = bool(re.search(r'\b(?:periportal\W+(?:\w+\W+){0,3}?scarring|scarring\W+(?:\w+\W+){0,1}?periportal)\b', line))
            port_scar = port_scar_1 or port_scar_2 or port_scar_3 or port_scar_4
            
#             stg_fib = bool(re.search(r'\b(?:fibrosis\W+(?:\w+\W+){0,8}?stage|stage\W+(?:\w+\W+){0,6}?fibrosis)\b', line))
#             stg_bri = bool(re.search(r'\b(?:bridging\W+(?:\w+\W+){0,8}?stage|stage\W+(?:\w+\W+){0,6}?bridging)\b', line))
#             stg_cir = bool(re.search(r'\b(?:cirrhosis\W+(?:\w+\W+){0,8}?stage|stage\W+(?:\w+\W+){0,6}?cirrhosis)\b', line))
            
#             sin_fib = bool(re.search(r'\b(?:sinusoidal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?sinusoidal)\b', line))
#             perisin_fib = bool(re.search(r'\b(?:perisinusoidal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?perisinusoidal)\b', line))
#             periport_fib = bool(re.search(r'\b(?:periportal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?periportal)\b', line))
#             bridg_fib = bool(re.search(r'\b(?:bridging\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?bridging)\b', line))
#             cent_fib = bool(re.search(r'\b(?:central\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?central)\b', line))
#             sept_fib = bool(re.search(r'\b(?:septal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?septal)\b', line))
#             periven_fib = bool(re.search(r'\b(?:perivenular\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?perivenular)\b', line))
#             pericel_fib = bool(re.search(r'\b(?:pericellular\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?pericellular)\b', line))
            centrilob_fib = bool(re.search(r'\b(?:centrilobular\W+(?:\w+\W+){0,2}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?centrilobular)\b', line))
            porttract_fib = bool(re.search(r'\b(?:portal tract\W+(?:\w+\W+){0,2}?fibrosis|fibrosis\W+(?:\w+\W+){0,3}?portaltract)\b', line))

            port_fib1 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,4}?fibrosis|fibrosis\W+(?:\w+\W+){0,2}?portal)\b', line))
            port_fib2 = bool(re.search(r'\b(?:fibrous\W+(?:\w+\W+){0,3}?portal)\b', line))
            
            port_exp_1 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,3}?expanded|expanded\W+(?:\w+\W+){0,1}?portal)\b', line))
            port_exp_2 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,2}?expansion|expansion\W+(?:\w+\W+){0,2}?portal)\b', line))
            port_exp = port_exp_1 or port_exp_2
            
            periport_exp_1 = bool(re.search(r'\b(?:periportal\W+(?:\w+\W+){0,2}?expanded|expanded\W+(?:\w+\W+){0,1}?periportal)\b', line))
            periport_exp_2 = bool(re.search(r'\b(?:periportal\W+(?:\w+\W+){0,1}?expansion|expansion\W+(?:\w+\W+){0,2}?periportal)\b', line))
            periport_exp = periport_exp_1 or periport_exp_2
            
            port_sep_1 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,4}?septa)\b', line))
            port_sep_2 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,4}?septal formation)\b', line))
            port_sep_3 = bool(re.search(r'\b(?:portal\W+(?:\w+\W+){0,4}?septae)\b', line))
            port_sep = (port_sep_1 or port_sep_2 or port_sep_3)
        
            
            
            if 'steatosis' in line and ('<5%' in line or 'less than 5%' in line) and 'steatosis' in e_text: #('5-33%' not in line)
                e_text = '<5% ' + e_text
                
            if 'lobular&portal' in line and 'inflammation' in e_text and not 'lobular&portal' in e_text and lob_port_inf:
                e_text = 'lobular&portal ' + e_text
            
            if 'inflammation' in e_text and not 'lobular' in e_text and lobu_inf and inf_count<=1:
                e_text = 'lobular ' + e_text
            
            if 'inflammation' in e_text and not 'lobular' in e_text and not 'portal' in e_text and not 'lobular-inflammation' in line and lobu_inf and inf_count>1:
                e_text = 'lobular ' + e_text
            
            if 'inflammation' in e_text and not 'portal' in e_text and port_inf and inf_count<=1:
                e_text = 'portal ' + e_text
                
            if 'inflammation' in e_text and not 'portal' in e_text and not 'lobular' in e_text and not 'portal inflammation' in line and port_inf and inf_count>1:
                e_text = 'portal ' + e_text
            
            if ('portal tract' in line or 'portal area' in line) and 'inflammation' in e_text and not 'portal' in e_text  and port_inf: #(port_inf or port_tract_inf)
                e_text = 'portal ' + e_text
                
            if 'lobular' in e_text and not 'activity' in e_text and lobu_act:
                e_text = e_text + ' (lobular activity)'
            if 'portal' in e_text and not 'activity' in e_text and port_act:
                e_text = e_text + ' (portal activity)'
            
            if 'lobular' in line and 'disarray' in e_text and not 'lobular' in e_text and lobu_dis:
                e_text = 'lobular ' + e_text
            
            if 'lobular&portal' in line and 'hepatitis' in e_text and not 'lobular&portal' in e_text and lobu_hep:
                e_text = 'lobular&portal ' + e_text
            
            if 'lobular' in line and ('hepatitis' in e_text and not 'lobular' in e_text) and lobu_hep:
                e_text = 'lobular ' + e_text
                
            if 'zone-3' in line and not 'zone-3 injury' in line and 'hepatitis' in e_text and not 'zone-3' in e_text and zone3_hep:
                e_text = 'zone-3 ' + e_text
                
            if 'zone-3' in line and not 'zone-3 injury' in line and 'inflammation' in e_text and not 'zone-3' in e_text and zone3_inf:
                e_text = 'zone-3 ' + e_text
                
            if 'fibrosis' in e_text and zone3_fib and not 'zone-3' in e_text:
                e_text = 'zone-3 ' + e_text
            
            if 'fibrosis' in e_text and zone1_fib and not 'zone-1' in e_text:
                e_text = 'zone-1 ' + e_text
                
            if 'fibrosis' in e_text and centrilob_fib and not 'centrilobular' in e_text:
                e_text = 'centrilobular ' + e_text
                
            if 'scarring' in e_text and cent_scar and not 'central' in e_text:
                e_text = 'central ' + e_text
                
            if 'scarring' in e_text and port_scar and not 'portal' in e_text:
                e_text = 'portal ' + e_text
                
            if 'infiltrate' in e_text and lobu_infil and not 'lobul' in e_text:
                e_text = 'lobular ' + e_text
                
            if 'fibrosis' in e_text and (port_fib1 or porttract_fib) and not 'portal' in e_text:
                e_text = e_text + ' (portal-fibrosis)'
                
            if 'fibrous' in e_text and port_fib2:
                e_text = e_text + ' (portal-fibrosis)'
                
            if 'expan' in e_text and (port_exp or periport_exp):
                e_text = e_text + ' (portal-fibrosis)'
            
            if 'brunt' in e_text and not 'fibrosis' in line:
                e_text = e_text + ' (fibrosis)'
                
            if 'septa' in e_text and port_sep and not 'periportal' in e_text:
                e_text = e_text + ' (periportal-fibrosis)'
                
            if 'fibrosis' in e_text:
                if 'pericellular' in line and not 'pericellular' in e_text:
                    e_text = e_text + ' (pericellular-fibrosis)'
                if 'sinusoidal' in line and not 'sinusoidal' in e_text:
                    e_text = e_text + ' (sinusoidal-fibrosis)'
                if 'periportal' in line and not 'periportal' in e_text:
                    e_text = e_text + ' (periportal-fibrosis)'
                if 'central' in line and not 'central' in e_text:
                    e_text = e_text + ' (central-fibrosis)'
                if 'septal' in line and not 'septal' in e_text:
                    e_text = e_text + ' (septal-fibrosis)'
                if 'ductal' in line and not 'ductal' in e_text:
                    e_text = e_text + ' (ductal-fibrosis)'
                if 'perivenular' in line and not 'perivenular' in e_text:
                    e_text = e_text + ' (perivenular-fibrosis)'
                
           
            if ('fibrosis' in e_text or 'bridging' in e_text or 'cirrhosis' in e_text) and 'stage' in line:
                
                line = (line
                        .replace('stage i-ii', 'stage 1-2').replace('stage ii-iii', 'stage 2-3')
                        .replace('stage iii-iv', 'stage 3-4').replace('stage iv-v', 'stage 4-5')
                        .replace('stage v-vi', 'stage 5-6')
                        
                        .replace('stage iii', 'stage 3').replace('stage ii', 'stage 2').replace('stage iv', 'stage 4')
                        .replace('stage vi', 'stage 6').replace('stage v', 'stage 5').replace('stage i)', 'stage 1)').replace('stage i ', 'stage 1 ')
                        
                        .replace('1/2', '1-2').replace('2/3', '2-3').replace('3/4 of 4', '3-4 of 4')
                        .replace('0 to 1', '0-1').replace('1 to 2', '1-2').replace('2 to 3', '2-3')
                        .replace('3 to 4', '3-4').replace('4 to 5', '4-5').replace('5 to 6', '5-6')
                        .replace('0 to focally 1', '0-1').replace('1 to focally 2', '1-2')
                        .replace('2 to focally 3', '2-3').replace('3 to focally 4', '3-4')
                        .replace('4 to focally 5', '4-5').replace('5 to focally 6', '5-6')
                        .replace('0-1', '0.5').replace('1-2', '1.5').replace('2-3', '2.5')
                        .replace('3-4', '3.5').replace('4-5', '4.5').replace('5-6', '5.5')
                        .replace('1b-2', '1.5')
                        .replace(' 1a ', ' 1 ').replace(' 2a ', ' 2 ').replace(' 3a ', ' 3 ')
                        .replace(' 4a ', ' 4 ').replace(' 5a ', ' 5 ').replace(' 6a ', ' 6 ')
                   )
                        
                stagelist_1 = re.findall(r'stage.*?(\d+(?:\.\d+)?).*?(\d+(?:\.\d+)?)', line)
                stagelist_2 = re.findall(r'stage.*?(\d+(?:\.\d+)?)', line)
                
                stagelist_ishak = re.findall(r'ishak.*?(\d+(?:\.\d+)?).*?(\d+(?:\.\d+)?)', line)
                
                dist_fib = abs(line.find('fibrosis')-line.find('stage'))
                dist_pbc = abs(line.find('pbc')-line.find('stage'))
                
                
                if not (len(stagelist_2)==0 or ('pbc' in line and dist_pbc<dist_fib)):
                    
                    stage_val = float(stagelist_2[0])
                    ref_val = 4.0

                    if len(stagelist_1)==1:

                        bool_1 = '.' in stagelist_1[0][0] and '.' in stagelist_2[0]
                        bool_2 = '.' not in stagelist_2[0]

                        if bool_1 or bool_2:
                            stage_val = float(stagelist_1[0][0])
                            ref_val = float(stagelist_1[0][1])
                    
                    if len(stagelist_ishak)==1:
                        stage_val = float(stagelist_ishak[0][0])
                        ref_val = float(stagelist_ishak[0][1])

                    if stage_val>=5:
                        ref_val = 6.0
                    
                    if ref_val<=6 and stage_val<=ref_val and (ref_val in [4,6]):
                        
                        if 'ishak' in line:
                            e_text = e_text + ' ishak'
                            ref_val = 6.0

                        e_text = e_text + ' stage: ' + str(stage_val) + '/' + str(ref_val)

                    if stage_val==0:
                        e_bool = True
                    else:
                        e_bool = False

            e_text = " ".join(e_text.split())
            
            
            entity_result = entity_result + e_text + ' ' + str(not e_bool) + '\n'
            
            entity_result = entity_result.replace('baloon', 'balloon').replace('ballon', 'balloon').replace('balon', 'balloon')
            
#         print(line)
        
    return entity_result
