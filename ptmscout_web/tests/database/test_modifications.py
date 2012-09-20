from tests.DBTestCase import DBTestCase
from ptmscout.database.protein import Protein, Species, getSpeciesByName
from ptmscout.database.experiment import Experiment, getExperimentById
from ptmscout.database.modifications import Phosphopep, Modification,\
    getModificationsByProtein
from ptmscout.database.user import User, getUserById


class TestModifications(DBTestCase):
    def test_get_modifications_should_return_mods_allowed_by_user_credentials(self):
        u = User("username", "name", "email", "institution")
        u.createUser("password")
        u.saveUser()
        
        p = Protein()
        p.sequence = "ABCDEFG"
        p.species = getSpeciesByName("homo sapiens")
        p.acc_gene = "ZOMG"
        p.name = "zomg this is a protein"
        p.date = "12-1986"
        
        p.saveProtein()
        
        exp = Experiment()
        exp.name = "A test experiment"
        exp.experiment_id=0
        exp.public = 0
        exp.author = ""
        exp.date = ""
        exp.published = ""
        exp.ambiguity = ""
        exp.export = 0
        exp.dataset = ""
        exp.submitter = ""
        exp.primaryModification='Y,T'
        exp.saveExperiment()
        
        exp2 = Experiment()
        exp2.name = "Another test experiment"
        exp2.experiment_id=0
        exp2.public = 0
        exp2.author = ""
        exp2.date = ""
        exp2.published = ""
        exp2.ambiguity = ""
        exp2.export = 0
        exp2.dataset = ""
        exp2.submitter = ""
        exp2.primaryModification='S'
        exp2.saveExperiment()
        
        self.session.flush()
        
        exp.grantPermission(u, 'view')
        exp.saveExperiment()
        
        self.session.flush()
        
        exp = getExperimentById(exp.id, u)
        u = getUserById(u.id)
        
        p1 = Phosphopep()
        p1.pep_tryps="blag"
        p1.pep_aligned="blag"
        p1.site_pos=1
        p1.site_type='Y'
        p1.pfam_site="SH3_1"
        p1.protein_id = p.id
        
        p2 = Phosphopep()
        p2.pep_tryps="blag2"
        p2.pep_aligned="blag2"
        p2.site_pos=1
        p2.site_type='T'
        p2.pfam_site="SH3_1"
        p2.protein_id = p.id
        
        p3 = Phosphopep()
        p3.pep_tryps="blag3"
        p3.pep_aligned="blag3"
        p3.site_pos=1
        p3.site_type='S'
        p3.pfam_site="SH3_1"
        p3.protein_id = p.id
        
        self.session.add(p1)
        self.session.add(p2)
        self.session.add(p3)
        
        mod1 = Modification()
        mod1.experiment_id = exp.id
        mod1.protein_id = p.id
        mod1.phosphopep = "blag"
        self.session.add(mod1)
        
        mod2 = Modification()
        mod2.experiment_id = exp.id
        mod2.protein_id = p.id
        mod2.phosphopep = "blag2"        
        self.session.add(mod2)
        
        mod3 = Modification()
        mod3.experiment_id = exp2.id
        mod3.protein_id = p.id
        mod3.phosphopep = "blag3"        
        self.session.add(mod3)
        
        self.session.flush()
        
        mod1.phosphopeps.append(p1)
        mod2.phosphopeps.append(p2)
        mod3.phosphopeps.append(p3)
        self.session.add(mod1)
        self.session.add(mod2)
        self.session.add(mod3)
        self.session.flush()
        
        modifications = getModificationsByProtein(p.id, u)
        
        phosphopep_ids = [ p.id for m in modifications for p in m.phosphopeps ]
        self.assertEqual([p1.id, p2.id], phosphopep_ids)