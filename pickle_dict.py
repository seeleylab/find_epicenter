import pickle
import nibabel as nib

parcel_dict = {parcel_index: nib.load('/data/mridata/jbrown/brains/brainnetome_suit_comb/vol_%d.nii' % parcel_index).get_data() for parcel_index in range(1, 274)}

pickle.dump(parcel_dict, open('brainnetome_suit_comb_274_parcel_dict.p', 'w'))