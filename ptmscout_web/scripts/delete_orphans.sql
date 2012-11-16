DELETE FROM `peptide` WHERE id NOT IN (SELECT peptide_id FROM `MS_modifications`);
DELETE FROM `protein` WHERE id NOT IN (SELECT protein_id FROM MS);
DELETE FROM `GO` WHERE id NOT IN (SELECT GO_id FROM protein_GO);