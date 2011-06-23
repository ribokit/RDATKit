import rna
import settings
import secondary_structure
import view
import mapping

s = 'UGCGCUUUUUUUUGCGCUUUUUUUUGCGCU'
molecule = rna.RNA(s)
data = mapping.MappingData(data=[0,0], seqpos=[2,26])
structures = secondary_structure.fold(molecule.sequence, mapping_data=data,algorithm='rnastructure')
ba = molecule.bootstrap(data, 5, nsamples=1, replacement=True)
for b in ba:
    ba[b] = str(ba[b]) +'%'
applet = view.VARNA(sequences=[molecule.sequence], structures=structures)
print applet.render(base_annotations=ba, annotation_by_helix=True, \
                    helix_function=(lambda x, y: max(x,y)))

