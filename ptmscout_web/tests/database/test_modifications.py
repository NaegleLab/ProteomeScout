from tests.DBTestCase import DBTestCase
from ptmscout.database.protein import Protein
from ptmscout.database.taxonomies import getSpeciesByName
from ptmscout.database.experiment import Experiment, getExperimentById,\
    ExperimentData
from ptmscout.database.modifications import Peptide, MeasuredPeptide,\
    getMeasuredPeptidesByProtein, findMatchingPTM, \
    countMeasuredPeptidesForExperiment, getModificationById, \
    getMeasuredPeptide, countProteinsForExperiment
from ptmscout.database.user import User, getUserById


class TestModifications(DBTestCase):

    def test_findPTM_should_get_ptms_with_strain_specific(self):
        taxonomy = ['Plasmodium falciparum (isolate 3D7)']
        modname = 'Sulfothreonine'
        target = 't'

        ptms, found_match, found_match_residue = findMatchingPTM(modname, residue=target, taxons=taxonomy)

        self.assertTrue(found_match)
        self.assertTrue(found_match_residue)
        self.assertEqual(1, len(ptms))
        self.assertEqual(modname, ptms[0].name)
        self.assertEqual('T', ptms[0].target)

    def test_countProteinsForExperiment(self):
        exp_id = 28
        cnt = countProteinsForExperiment(exp_id)

        self.assertEqual(52, cnt)

    def test_countMeasuredPeptidesForExperiment(self):
        exp_id = self.test_get_modifications_should_return_mods_allowed_by_user_credentials()
        cnt = countMeasuredPeptidesForExperiment(exp_id)

        self.assertEqual(2, cnt)

    def test_getMeasuredPeptide_should_allow_filter_by_modifications(self):
        ms1 = getMeasuredPeptide(1, 'AAAACLsRQAssDS', 2193)

        ptm = getModificationById(4575)
        ptm2 = getModificationById(4599)

        pep1 = 'AAAAACLsRQASSDS'
        pep2 = 'ACLSRQAsSDSDSIL'
        pep3 = 'CLSRQASsDSDSILS'
        ms2 = getMeasuredPeptide(1, 'AAAACLsRQAssDS', 2193, [(pep1,ptm),(pep2,ptm),(pep3,ptm)])
        ms3 = getMeasuredPeptide(1, 'AAAACLsRQAssDS', 2193, [(pep1,ptm),(pep2,ptm),(pep3,ptm2)])


        self.assertEqual(ms1,ms2)
        self.assertEqual(None, ms3)

    def test_get_modifications_should_return_mods_allowed_by_user_credentials(self):
        u = User("username", "name", "email", "institution")
        u.createUser("password")
        u.saveUser()
        
        p = Protein()
        p.sequence = "ABCDEFG"
        p.species = getSpeciesByName("homo sapiens")
        p.acc_gene = "ZOMG"
        p.name = "zomg this is a protein"
        p.locus = "ZOMG_HUMAN"
        p.date = "12-1986"
        
        p.saveProtein()
        
        exp = Experiment()
        exp.name = "A test experiment"
        exp.experiment_id=None
        exp.public = 0
        exp.author = ""
        exp.status = 'loaded'
        exp.published = 0
        exp.ambiguity = 0
        exp.export = 0
        exp.dataset = ""
        exp.submitter = ""
        exp.primaryModification='Y,T'
        exp.submitter_id = None
        exp.saveExperiment()
        
        exp2 = Experiment()
        exp2.name = "Another test experiment"
        exp2.experiment_id=None
        exp2.public = 0
        exp2.author = ""
        exp2.status = 'loaded'
        exp2.published = 0
        exp2.ambiguity = 0
        exp2.export = 0
        exp2.dataset = ""
        exp2.submitter = ""
        exp2.primaryModification='S'
        exp2.submitter_id = None
        exp2.saveExperiment()
        

        exp3 = Experiment()
        exp3.name = "A third test experiment"
        exp3.experiment_id=None
        exp3.public = 0
        exp3.author = ""
        exp3.published = 0
        exp3.status = 'loading'
        exp3.ambiguity = 0
        exp3.export = 0
        exp3.dataset = ""
        exp3.submitter = ""
        exp3.primaryModification='S'
        exp3.submitter_id = None
        exp3.saveExperiment()
        
        self.session.flush()
        
        exp.grantPermission(u, 'view')
        exp.saveExperiment()
        
        exp3.grantPermission(u, 'owner')
        exp.saveExperiment()
        
        self.session.flush()
        
        exp = getExperimentById(exp.id, u)
        u = getUserById(u.id)
        
        phos_mods, _exist, _exist_residue = findMatchingPTM('Phosphorylation')
        
        tyr_mod = None
        the_mod = None
        ser_mod = None
        
        for mod in phos_mods:
            if mod.target=='Y':
                tyr_mod = mod
            if mod.target=='T':
                the_mod = mod
            if mod.target=='S':
                ser_mod = mod
        
        p1 = Peptide()
        p1.pep_tryps="blag"
        p1.pep_aligned="blag"
        p1.site_pos=1
        p1.site_type='Y'
        p1.protein_id = p.id
        
        p2 = Peptide()
        p2.pep_tryps="blag2"
        p2.pep_aligned="blag2"
        p2.site_pos=1
        p2.site_type='T'
        p2.protein_id = p.id
        
        p3 = Peptide()
        p3.pep_tryps="blag3"
        p3.pep_aligned="blag3"
        p3.site_pos=1
        p3.site_type='S'
        p3.protein_id = p.id
        
        p4 = Peptide()
        p4.pep_tryps="blag4"
        p4.pep_aligned="blag4"
        p4.site_pos=1
        p4.site_type='S'
        p4.protein_id = p.id
        
        self.session.add(p1)
        self.session.add(p2)
        self.session.add(p3)
        
        mod1 = MeasuredPeptide()
        mod1.query_accession='blah1'
        mod1.experiment_id = exp.id
        mod1.protein_id = p.id
        mod1.peptide = "blag"
        self.session.add(mod1)
        
        mod2 = MeasuredPeptide()
        mod2.query_accession='blah2'
        mod2.experiment_id = exp.id
        mod2.protein_id = p.id
        mod2.peptide = "blag2"        
        self.session.add(mod2)
        
        mod3 = MeasuredPeptide()
        mod3.query_accession='blah3'
        mod3.experiment_id = exp2.id
        mod3.protein_id = p.id
        mod3.peptide = "blag3"        
        self.session.add(mod3)
        
        mod4 = MeasuredPeptide()
        mod4.query_accession='blah4'
        mod4.experiment_id = exp3.id
        mod4.protein_id = p.id
        mod4.peptide = "blag4"        
        self.session.add(mod4)
        
        self.session.flush()
        
        d1 = ExperimentData()
        d1.type = "time(min)"
        d1.MS_id = mod1.id
        d1.priority = 1
        d1.label = 0
        d1.value = 10
        
        d2 = ExperimentData()
        d2.type = "time(min)"
        d2.MS_id = mod1.id
        d2.priority = 2
        d2.label = 1
        d2.value = 20
        
        d3 = ExperimentData()
        d3.type = "time(min)"
        d3.MS_id = mod1.id
        d3.priority = 3
        d3.label = 10
        d3.value = 30
        
        mod1.data.append(d1)
        mod1.data.append(d2)
        mod1.data.append(d3)
        
        mod1.addPeptideModification(p1, tyr_mod)
        mod2.addPeptideModification(p2, the_mod)
        mod3.addPeptideModification(p3, ser_mod)
        mod4.addPeptideModification(p4, ser_mod)
        
        self.session.add(mod1)
        self.session.add(mod2)
        self.session.add(mod3)
        self.session.flush()
        
        modifications = getMeasuredPeptidesByProtein(p.id, u)
        
        Peptide_ids = [ p.peptide.id for m in modifications for p in m.peptides ]
        self.assertEqual([p1.id, p2.id], Peptide_ids)
        
        self.assertEqual(['time(min)']*3, [d.type for d in modifications[0].data])
        self.assertEqual([0,1,10], sorted([d.label for d in modifications[0].data]))

        return exp.id

if __name__=='__main__':
    import unittest
    unittest.main()
