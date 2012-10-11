from paste.deploy.loadwsgi import appconfig
from ptmscout.database import DBSession  # base declarative object
import os, sys
from ptmscout.database import gene_expression as ge
import traceback
from DB_init import DatabaseInitialization

def getExpressionTissues(table_name):
    result = DBSession.execute("SHOW COLUMNS FROM `%s`" % table_name)
    
    cols = [ row[0] for row in result ]
    cols.remove('probeset_id')
    cols.remove('id')
    
    return cols

if __name__ == '__main__':
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'development.ini'))
            
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        
        arg_dict = {}
        for i, item in enumerate(sys.argv):
            arg_dict[i] = item
            
        run_type = arg_dict.get(1, "test")
        
        if(run_type == "build_tables"):
            print "Getting tissue names..."
            
            human_tissues = getExpressionTissues('expression_human')
            mouse_tissues = getExpressionTissues('expression_mouse')
            NCI60_tissues = getExpressionTissues('expression_NCI60')
            
            tissues = set()
            for t in human_tissues + mouse_tissues + NCI60_tissues:
                tissues.add(t)
            
            tissues.remove("CD8+T-cells")
            tissues.remove("CD4+T-cells")
            
            
            print "Inserting tissue and collection types"
            
            for tissue in tissues:
                DBSession.execute("INSERT INTO `expression_tissue` (`name`) VALUES ('%s')" % (tissue))
            
            DBSession.execute("INSERT INTO `expression_collection` (`name`) VALUES ('Human')")
            DBSession.execute("INSERT INTO `expression_collection` (`name`) VALUES ('Mouse')")
            DBSession.execute("INSERT INTO `expression_collection` (`name`) VALUES ('NCI60')")
            
            
            print "Getting DB objects for inserted tissue and collection types"
            
            exp_collections = DBSession.query(ge.ExpressionCollection).all()
            
            collection_map = {}
            for collection in exp_collections:
                collection_map[collection.name] = collection
            
            tissue_refs = {}
            
            exp_tissues = DBSession.query(ge.ExpressionTissue).all()
            
            for t in exp_tissues:
                tissue_refs[t.name.lower()] = t
            
            
            print "Remapping probeset accessions"
            
            probesets = DBSession.query(ge.ExpressionProbeset).all()
            
            def addAllSyms(probe, t, field, delimiter, replace = None):
                symbols = getattr(probe, field)
                
                if replace != None:
                    symbols = symbols.replace(replace[0], replace[1])
                
                symbols = symbols.split(delimiter)
                
                for sym in symbols:
                    sym = sym.strip()
                    if len(sym) > 0:
                        acc = ge.ExpressionAccession(t, sym, probe.id)
                        DBSession.add(acc)
    
            probeset_id_map = {}
    
            c = 0        
            for probe in probesets:
                addAllSyms(probe, 'gene_symbol', 'symbol', ';', ('///',';'))
                addAllSyms(probe, 'refseq', 'refseq', ';', ('///',';'))
                addAllSyms(probe, 'uniprot', 'uniprot', ';', ('///',';'))
                addAllSyms(probe, 'alias', 'aliases', ';', ('///',';'))
                DBSession.flush()
                
                probeset_id_map[probe.probeset_id] = probe
                
                c+=1
                if c % 1000 == 0:
                    print c/1000,
                    
            print ""
            
            
            def remapTable(table_name, collection, tissue_set):
                c = 0
                for row in (DBSession.execute("SELECT * FROM %s" % table_name)):
                    probeset_id = row[0]
                    #id = row[-1]
                    
                    for i, t in enumerate(tissue_set):
                        value = row[i+1]
                        
                        sample = ge.ExpressionSample()
                        sample.probeset_id = probeset_id_map[probeset_id].id
                        sample.collection_id = collection_map[collection].id 
                        sample.tissue_id = tissue_refs[t.lower()].id
                        sample.value = value
                        
                        DBSession.add(sample)
                        
                    DBSession.flush()
                        
                    c+=1
                    if c % 1000 == 0:
                        print c/1000,
                print ""
            
            print "Remapping human expression table..."
            remapTable('expression_human', 'Human', human_tissues)
            print "Remapping mouse expression table..."
            remapTable('expression_mouse', 'Mouse', mouse_tissues)
            print "Remapping NCI60 expression table..."
            remapTable('expression_NCI60', 'NCI60', NCI60_tissues)

        if run_type == "build_tables" or run_type == "test":
            print "Testing results..."
            
            human_sample = DBSession.query(ge.ExpressionProbeset).filter_by(probeset_id='1007_s_at').first()
            accs = [ acc.value for acc in human_sample.accessions ]
            expected_accs = ['DDR1', 'NP_001945.3', 'NP_054699.2', 'NP_054700.2', 'A2ABK8', 'A2ABK9', 'A2ABL0', 'A2ABL1', 'A2ABL2', 'A2ABL3', 'A2ABL4', 'A2ABM8', 'Q08345', 'Q2L6H3', 'Q4LE50', 'Q5ST11', 'Q6NSK4', 'Q6ZNR9', 'Q96T61', 'Q96T62', 'Q9UD36', 'Q9UD86']
            
            sort_order = '721 B lymphoblasts    Adipocyte    Adrenal Cortex    Amygdala    Appendix    BM-CD105+Endothelial    BM-CD33+Myeloid    BM-CD34+    BM-CD71+EarlyErythroid    CD4+T-Cells    CD8+T-Cells    Cardiac Myocytes    Cerebellum    Cerebellum Peduncles    Ciliary Ganglion    Cingulate Cortex    Colorectal Adenocarcinoma    Heart    Hypothalamus    Kidney    Liver    Lung    Medulla Oblongata    Occipital Lobe    Olfactory Bulb    Ovary    PB-BDCA4+Dentritic_Cells    PB-CD14+Monocytes    PB-CD19+Bcells    PB-CD56+NKCells    Placenta    Pancreas    Pancreatic Islets    Parietal Lobe    Pituitary    Pons    Prefrontal Cortex    Prostate    Skeletal Muscle    Smooth Muscle    Superior Cervical Ganglion    Tongue    Temporal Lobe    Testis    Testis Germ Cell    Testis Interstitial    Testis Leydig Cell    Testis Seminiferous Tubule    Thalamus    Thymus    Thyroid    Tonsil    Trachea    Trigeminal Ganglion    Uterus    Uterus Corpus    Whole Blood    Whole Brain    Adrenal Gland    Atrioventricular Node    Bone Marrow    Bronchial Epithelial Cells    Caudate Nucleus    Dorsal Root Ganglia    Fetal Thyroid    Fetal Brain    Fetal Liver    Fetal Lung    Globus Pallidus    Leukemia Chronicmyelogenous(k562)    Leukemia Lymphoblastic(molt4)    Leukemia Myelocytic(hl60)    Lymph Node    Lymphoma Burkitts Daudi    Lymphoma Burkitts Raji    Salivary Gland    Skin    Spinal Cord    subthalamicnucleus'.split('    ')
            exp_values = [ float(v) for v in '3170.95    904.3    1094.65    11161.8    636.55    303.2    364    712.4    201.65    683.9    494.1    1616.7    3902.7    3222.6    423.1    2965.05    10165.6    1666.65    18497.3    2815.05    1594.05    6480.1    3547.3    5168.8    5634.95    672.6    540.85    221.5    1065.05    269.9    12693.5    2951.1    10562.3    5229.6    5178.5    2128.5    12819.3    11181.2    386.05    1894.65    728.5    2981.2    1251.6    2259.7    1297.85    732.8    1316.15    1233.3    6476.35    7216.7    20480.2    1171    6831.7    576.35    2207.55    910.2    1197.85    7684.4    2112.9    368.85    813.95    18162.3    5532.8    1069.9    10331.2    11092.5    769.35    3914.5    1227.4    1314    760.75    844.1    1129    818.8    1391.4    2479.6    1582.25    26996.8    1766.1'.split() ]
            
            exp_val_dict = {}
            for i,s in enumerate(sort_order):
                exp_val_dict[s] = exp_values[i]
            
            exp_values = sorted( exp_val_dict.items(), key=lambda item: item[0] )
            values = sorted( [ (s.tissue.name, s.value) for s in human_sample.samples if s.collection.name == 'Human' ], key=lambda item: item[0] )
            
            print "Testing sample 1..."
            assert human_sample.species.name == "homo sapiens"
            assert sorted(expected_accs) == sorted(accs), "Lists do not match:\n%s\n%s\n" % (str(sorted(expected_accs)), str(sorted(accs)))
            assert exp_values == values, "Lists do not match:\n%s\n%s\n" % (exp_values, values)
            
            mouse_sample = DBSession.query(ge.ExpressionProbeset).filter_by(probeset_id='gnf1m00001_at').first()
            accs = [ acc.value for acc in mouse_sample.accessions ]
            expected_accs = ['Serpinb3a']
            
            sort_order = 'B220+B-cells    CD4+T-Cells    CD8+T-Cells    Adipose Tissue    Adrenal Gland    Amygdala    Bladder    Blastocysts    Bone    Bone Marrow    Brown Fat    Cerebellum    Cerebral Cortex    Digits    Dorsal Root Ganglia    Dorsal Striatum    embryo day 10.5    embryo day 6.5    embryo day 7.5    embryo day 8.5    embryo day 9.5    Epidermis    Fertilized Egg    Frontal Cortex    Heart    Hippocampus    Hypothalamus    Kidney    Large Intestine    Liver    Lung    Lymph Node    Mammary Gland (lact)    Medial Olfactory Epithelium    Olfactory Bulb    Oocyte    Ovary    Pancreas    Pituitary    Placenta    Preoptic    Prostate    Retina    Salivary Gland    Skeletal Muscle    Small Intestine    Snout Epidermis    Spinal Cord Lower    Spinal Cord Upper    Spleen    Stomach    Substantial Nigra    Testis    Thymus    Thyroid    Tongue    Trachea    Trigeminal    Umbilical Cord    Uterus    Vomeralnasal Organ'.split('    ')            
            exp_values = [ float(v) for v in '107.7    133.1    155.35    101.95    102.5    103.95    146.35    113.15    140.25    261.15    121    113.25    105.45    982.05    102.3    198.2    139.15    159.95    136.4    127.05    136.1    181.5    146.8    94.6    128.7    134.95    101.45    115.5    99.65    157    104.75    136.4    126.05    86.15    105.2    166.1    104.25    224.8    96.95    147.55    90.75    116.3    352.6    171.45    141.4    116.3    1739.75    93.8    94.95    288.2    120.45    97    119.85    128.15    166.6    106108    19875    97.35    471.95    110.45    100.6'.split() ]
            
            exp_val_dict = {}
            for i,s in enumerate(sort_order):
                exp_val_dict[s] = exp_values[i]
            
            exp_values = sorted( exp_val_dict.items(), key=lambda item: item[0] )
            
            values = sorted( [ (s.tissue.name, s.value) for s in mouse_sample.samples if s.collection.name == 'Mouse' ], key=lambda item: item[0] )
            
            print "Testing sample 2..."
            assert mouse_sample.species.name == "mus musculus"
            assert sorted(expected_accs) == sorted(accs), "Lists do not match:\n%s\n%s\n" % (str(sorted(expected_accs)), str(sorted(accs)))
            assert exp_values == values, "Lists do not match:\n%s\n%s\n" % (exp_values, values)
            
            NCI60_sample = DBSession.query(ge.ExpressionProbeset).filter_by(probeset_id='1255_g_at').first()
            accs = [ acc.value for acc in NCI60_sample.accessions ]
            expected_accs = ['GUCA1A', 'NP_000400.2', 'A6NIT3', 'P43080']
            
            sort_order = '786-0    A172    A204    A361    A498    A549    ACC3    ACHN    ALVA31    CAKI1    CCRT_CEM    COLO205    DLD1    EKVX    HCC2998    HCT116    HCT15    HEK 293    HEK 293T    HELA    HEPG2    HOP62    HOP92    HS578T    HSG    HT1080    IGROV1    IOSE80    JURKAT    K562    KM12    LN18    LNCAP    LOX_IMVI    M14    MALME_3M    MCF7    MDA_MB231    MDA_MB435    MOLT4    NCI-460    NCI-ADR-RES    NCI-H226    NCI-H23    NCI-H322    OVCAR3    OVCAR4    OVCAR5    OVCAR8    PANC1    RL7    RPMI_8226    RS11846    RXF393    SAOS2    SF268    SF295    SF539    SHSYSY+RA    SHSYSY-RA    SKMEL2    SKMEL28    SKMEL5    SKOV3    SN12C    SNB19    SNB75    SR    SW620    T24    T3M4    TK10    U118    U138    U20S    U251    U87    UACC257    UACC62    UO31    ZR75    astrocytes    huh-7'.split('    ')
            exp_values = [ float(v) for v in '9.28498    8.48497    10.2919    9.70591    9.32519    9.26441    9.96471    9.23604    10.7122    9.71672    13.762    9.39302    9.29439    8.77076    8.48583    9.28798    8.82024    9.27564    9.28584    8.92726    9.47159    9.35291    9.07691    9.34741    9.25657    9.08692    9.47966    9.33211    9.37901    10.4807    9.68824    8.52402    9.56472    9.6189    8.94914    10.8523    8.71031    10.3445    9.94654    10.0382    10.6916    8.49493    9.66228    9.48966    9.17454    9.62523    9.26272    9.22584    8.47142    9.51672    9.822    9.65119    10.0936    9.3166    9.32509    9.20707    9.18708    8.75416    53.6318    74.5769    9.24819    9.27232    8.84921    9.38438    9.0581    8.51357    9.2233    10.199    9.8988    9.36692    9.70918    9.34133    8.89915    8.64985    9.61672    9.1199    8.5732    9.17579    9.32824    8.57601    9.86936    10.0235    9.52266'.split() ]
            
            exp_val_dict = {}
            for i,s in enumerate(sort_order):
                exp_val_dict[s] = exp_values[i]
            
            exp_values = sorted( exp_val_dict.items(), key=lambda item: item[0] )
            
            values = sorted( [ (s.tissue.name, s.value) for s in NCI60_sample.samples if s.collection.name == 'NCI60' ], key=lambda item: item[0] )
            
            print "Testing sample 3..."
            assert NCI60_sample.species.name == "homo sapiens"
            assert sorted(expected_accs) == sorted(accs), "Lists do not match:\n%s\n%s\n" % (str(sorted(expected_accs)), str(sorted(accs)))
            assert exp_values == values, "Lists do not match:\n%s\n%s\n" % (exp_values, values)
            
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()