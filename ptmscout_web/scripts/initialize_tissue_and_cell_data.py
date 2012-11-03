from paste.deploy.loadwsgi import appconfig
from ptmscout.database import DBSession  # base declarative object
import os
from ptmscout.database import experiment
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
        
        tissue_types = []
        cell_types = []
        
        print "Loading cell and tissue types..."
        
        items = '721 B lymphoblasts    Adipocyte    Adrenal Cortex    Amygdala    Appendix    BM-CD105+Endothelial    BM-CD33+Myeloid    BM-CD34+    BM-CD71+EarlyErythroid    CD4+T-Cells    CD8+T-Cells    Cardiac Myocytes    Cerebellum    Cerebellum Peduncles    Ciliary Ganglion    Cingulate Cortex    Colorectal Adenocarcinoma    Heart    Hypothalamus    Kidney    Liver    Lung    Medulla Oblongata    Occipital Lobe    Olfactory Bulb    Ovary    PB-BDCA4+Dentritic_Cells    PB-CD14+Monocytes    PB-CD19+Bcells    PB-CD56+NKCells    Placenta    Pancreas    Pancreatic Islets    Parietal Lobe    Pituitary    Pons    Prefrontal Cortex    Prostate    Skeletal Muscle    Smooth Muscle    Superior Cervical Ganglion    Tongue    Temporal Lobe    Testis    Testis Germ Cell    Testis Interstitial    Testis Leydig Cell    Testis Seminiferous Tubule    Thalamus    Thymus    Thyroid    Tonsil    Trachea    Trigeminal Ganglion    Uterus    Uterus Corpus    Whole Blood    Whole Brain    Adrenal Gland    Atrioventricular Node    Bone Marrow    Bronchial Epithelial Cells    Caudate Nucleus    Dorsal Root Ganglia    Fetal Thyroid    Fetal Brain    Fetal Liver    Fetal Lung    Globus Pallidus    Leukemia Chronicmyelogenous(k562)    Leukemia Lymphoblastic(molt4)    Leukemia Myelocytic(hl60)    Lymph Node    Lymphoma Burkitts Daudi    Lymphoma Burkitts Raji    Salivary Gland    Skin    Spinal Cord    subthalamicnucleus'.split('    ')
        tissue_types.extend(items)
        items = 'B220+B-cells    CD4+T-Cells    CD8+T-Cells    Adipose Tissue    Adrenal Gland    Amygdala    Bladder    Blastocysts    Bone    Bone Marrow    Brown Fat    Cerebellum    Cerebral Cortex    Digits    Dorsal Root Ganglia    Dorsal Striatum    embryo day 10.5    embryo day 6.5    embryo day 7.5    embryo day 8.5    embryo day 9.5    Epidermis    Fertilized Egg    Frontal Cortex    Heart    Hippocampus    Hypothalamus    Kidney    Large Intestine    Liver    Lung    Lymph Node    Mammary Gland (lact)    Medial Olfactory Epithelium    Olfactory Bulb    Oocyte    Ovary    Pancreas    Pituitary    Placenta    Preoptic    Prostate    Retina    Salivary Gland    Skeletal Muscle    Small Intestine    Snout Epidermis    Spinal Cord Lower    Spinal Cord Upper    Spleen    Stomach    Substantial Nigra    Testis    Thymus    Thyroid    Tongue    Trachea    Trigeminal    Umbilical Cord    Uterus    Vomeralnasal Organ'.split('    ')            
        tissue_types.extend(items)
        items = '786-0    A172    A204    A361    A498    A549    ACC3    ACHN    ALVA31    CAKI1    CCRT_CEM    COLO205    DLD1    EKVX    HCC2998    HCT116    HCT15    HEK 293    HEK 293T    HELA    HEPG2    HOP62    HOP92    HS578T    HSG    HT1080    IGROV1    IOSE80    JURKAT    K562    KM12    LN18    LNCAP    LOX_IMVI    M14    MALME_3M    MCF7    MDA_MB231    MDA_MB435    MOLT4    NCI-460    NCI-ADR-RES    NCI-H226    NCI-H23    NCI-H322    OVCAR3    OVCAR4    OVCAR5    OVCAR8    PANC1    RL7    RPMI_8226    RS11846    RXF393    SAOS2    SF268    SF295    SF539    SHSYSY+RA    SHSYSY-RA    SKMEL2    SKMEL28    SKMEL5    SKOV3    SN12C    SNB19    SNB75    SR    SW620    T24    T3M4    TK10    U118    U138    U20S    U251    U87    UACC257    UACC62    UO31    ZR75    astrocytes    huh-7'.split('    ')
        cell_types.extend(items)
        
        for t in tissue_types:
            expc = experiment.ExperimentCondition()
            expc.type = 'tissue'
            expc.value = t
            DBSession.add(expc)
            
        for c in cell_types:
            expc = experiment.ExperimentCondition()
            expc.type = 'cell'
            expc.value = c
            DBSession.add(expc)
            
#        print "Loading drugbank..."
#        with open('data/drugbank.txt', 'r') as drugbank:
#            for line in drugbank:
#                line = line.strip()
#                if line == '':
#                    continue
#                
#                if line[0] == "#":
#                    capture = False
#                    
#                if capture:
#                    expc = experiment.ExperimentCondition()
#                    expc.type = 'drug'
#                    expc.value = line
#                    DBSession.add(expc)
#                    
#                if line == "# Brand_Names:":
#                    capture = True
#                if line == "# Generic_Name:":
#                    capture = True
        
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()