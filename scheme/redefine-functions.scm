
;; scm aliases.  
;; 
(define enhanced-ligand-coot? enhanced-ligand-coot-p)
(define drag-intermediate-atom drag-intermediate-atom-scm)
(define mark-atom-as-fixed mark-atom-as-fixed-scm)
(define ncs-chain-ids ncs-chain-ids-scm)
(define ncs-chain-differences ncs-chain-differences-scm)
(define refmac-parameters refmac-parameters-scm)
(define water-chain water-chain-scm)
(define water-chain-from-shelx-ins water-chain-from-shelx-ins-scm)
(define generic-object-name generic-object-name-scm)
(define generic-object-is-displayed? generic-object-is-displayed-p)
(define is-closed-generic-object? is-closed-generic-object-p)
(define additional-representation-info additional-representation-info-scm)
(define missing-atom-info missing-atom-info-scm)
(define key-sym-code key-sym-code-scm)
(define ncs-ghosts ncs-ghosts-scm)
(define inverse-rtop inverse-rtop-scm)
(define coot-has-python? coot-has-python-p)
(define map-sigma map-sigma-scm)
(define pucker-info pucker-info-scm)
(define map-parameters map-parameters-scm)
(define cell cell-scm)
(define map-cell cell)
(define ccp4i-projects ccp4i-projects-scm)
(define map-sigma map-sigma-scm)
(define get-rotamer-name get-rotamer-name-scm)
(define test-internal test-internal-scm)
(define atom-info-string atom-info-string-scm)
(define get-refmac-sad-atom-info get-refmac-sad-atom-info-scm)
;; fix typo of set-find-hydrogen-torsions (backward compatibility in case anyone was using that)
(define set-find-hydrogen-torsion set-find-hydrogen-torsions)
(define residues-near-position residues-near-position-scm)
(define residues-near-residues residues-near-residues-scm)
(define non-standard-residue-names non-standard-residue-names-scm)
(define refine-residues refine-residues-scm)
(define refine-residues-with-alt-conf refine-residues-with-alt-conf-scm)
(define regularize-residues regularize-residues-scm)
(define refine-zone-with-score refine-zone-with-score-scm)
(define regularize-zone-with-score regularize-zone-with-score-scm)
(define set-merge-molecules-ligand-spec set-merge-molecules-ligand-spec-scm)
(define map-peaks map-peaks-scm)
(define map-peaks-near-point map-peaks-near-point-scm)
(define add-dipole add-dipole-scm)
(define add-dipole-for-residues add-dipole-for-residues-scm)
(define get-torsion get-torsion-scm)
(define test-internal-single test-internal-single-scm)
(define user-defined-click user-defined-click-scm)
(define add-lsq-atom-pair add-lsq-atom-pair-scm)
(define coot-sys-build-type coot-sys-build-type-scm)
(define add-alt-conf add-alt-conf-scm)
(define origin-pre-shift origin-pre-shift-scm)
(define alignment-mismatches alignment-mismatches-scm)
(define average-map average-map-scm)
(define symmetry-operators symmetry-operators-scm)
(define symmetry-operators->xHM symmetry-operators-to-xHM-scm)
(define user-mods user-mods-scm)
(define refine-zone-with-full-residue-spec refine-zone-with-full-residue-spec-scm)
(define pkgdatadir get-pkgdatadir-scm)
(define matching-compound-names-from-sbase matching-compound-names-from-sbase-scm)
(define chain-id chain-id-scm)
(define highly-coordinated-waters highly-coordinated-waters-scm)
(define handle-pisa-interfaces handle-pisa-interfaces-scm)
(define space-group space-group-scm)
(define nearest-residue-by-sequence nearest-residue-by-sequence-scm)
(define set-torsion set-torsion-scm)
(define get-map-colour get-map-colour-scm)
(define list-extra-restraints list-extra-restraints-scm)
(define delete-extra-restraint delete-extra-restraint-scm)
(define delete-extra-restraints-for-residue-spec delete-extra-restraints-for-residue-spec-scm)
(define do-clipped-surface do-clipped-surface-scm)
(define copy-residue-range-from-ncs-master-to-chains copy-residue-range-from-ncs-master-to-chains-scm)
(define copy-from-ncs-master-to-chains copy-from-ncs-master-to-chains-scm)
(define remarks remarks-scm)
(define protein-db-loops protein-db-loops-scm)
(define residue-centre residue-centre-scm)
(define make-link make-link-scm)
(define ncs-master-chains ncs-master-chains-scm)
(define get-lsq-matrix get-lsq-matrix-scm)
(define list-nomenclature-errors list-nomenclature-errors-scm)
(define matching-compound-names-from-dictionary matching-compound-names-from-dictionary-scm)
(define comp-id->name comp-id-to-name-scm)
(define multi-residue-torsion multi-residue-torsion-scm)
(define multi-residue-torsion-fit multi-residue-torsion-fit-scm)
(define add-linked-residue add-linked-residue-scm)
(define set-go-to-atom-from-res-spec set-go-to-atom-from-res-spec-scm)
(define test-function test-function-scm)
(define make-variance-map make-variance-map-scm)
;; (define residue->sdf-file residue-to-sdf-file)
(define cif-file-for-comp-id cif-file-for-comp-id-scm)
(define display-maps display-maps-scm)
(define ligand-search-make-conformers ligand-search-make-conformers-scm)
(define coot-can-do-lidia? coot-can-do-lidia-p) ;; a c++ bool, which SWIG translates
(define dictionary-entries dictionary-entries-scm)
(define SMILES-for-comp-id SMILES-for-comp-id-scm)
(define (molecule-is-drawn-as-surface? imol) (= (molecule-is-drawn-as-surface-int imol) 1))
(define molecule-name-stub molecule-name-stub-scm)
(define map-to-model-correlation map-to-model-correlation-scm)
(define map-to-model-correlation-per-residue map-to-model-correlation-per-residue-scm)
(define het-group-residues het-group-residues-scm)
(define score-rotamers score-rotamers-scm)
(define link-info link-info-scm)
(define kolmogorov-smirnov     kolmogorov-smirnov-scm)
(define kullback-liebler       kullback-liebler-scm)
(define residue-name           residue-name-scm)
(define chain-fragments        chain-fragments-scm)
(define find-blobs             find-blobs-scm)
(define morph-fit-residues     morph-fit-residues-scm)
(define density-score-residue  density-score-residue-scm)
(define align-to-closest-chain align-to-closest-chain-scm)
(define map-statistics         map-statistics-scm)
(define goto-next-atom-maybe   goto-next-atom-maybe-scm)
(define goto-prev-atom-maybe   goto-prev-atom-maybe-scm)
(define glyco-tree-residues   glyco-tree-residues-scm)
(define glyco-tree-residue-id glyco-tree-residue-id-scm)
(define add-extra-bond-restraints add-extra-bond-restraints-scm)
(define new-molecule-by-residue-specs new-molecule-by-residue-specs-scm)
(define delete-residues delete-residues-scm)
(define qq-plot-map-and-model qq-plot-map-and-model-scm)
(define refinement-already-ongoing? (lambda() (= (refinement-already-ongoing-p) 1)))
(define linked-residues linked-residues-scm)
(define molecule-atom-overlaps molecule-atom-overlaps-scm)
(define CG-spin-search CG-spin-search-scm)
(define molecule-atom-overlaps molecule-atom-overlaps-scm)
(define c-beta-deviations c-beta-deviations-scm)
(define set-go-to-atom-from-atom-spec set-go-to-atom-from-atom-spec-scm)
(define chiral-volume-errors chiral-volume-errors-scm)

;; (define residues-distortions residues-distortions-scm)

;; I changed the function name - save those (just a few) with scripts that I've handed out
(define toggle-idle-ligand-interactions toggle-flev-idle-ligand-interactions)
