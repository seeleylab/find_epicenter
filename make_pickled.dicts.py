import pickle
import nibabel as nib

epicenter_parcel_dict = {parcel_index: nib.load('/data/mridata/jbrown/brains/brainnetome_suit_comb/vol_%d.nii' % parcel_index).get_data() for parcel_index in range(1, 274)}
pickle.dump(epicenter_parcel_dict, open('epicenter_parcel_dict.p', 'w'))

epicenter_seedmap_dict = {parcel_index: nib.load('/data/mridata/jdeng/sd_bvftd/100_controls/1ST_vol_%s/spmT_0001.nii' % parcel_index).get_data() for parcel_index in range(1, 274)}
pickle.dump(epicenter_seedmap_dict, open('epicenter_seedmap_dict.p', 'w'))
