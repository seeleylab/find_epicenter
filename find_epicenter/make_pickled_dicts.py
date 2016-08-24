import pickle
import nibabel as nib

epicenter_parcel_dict = {parcel_index: nib.load('/data/mridata/jbrown/brains/brainnetome_suit_comb/vol_%d.nii' % parcel_index).get_data() for parcel_index in range(1, 274)}
with open('epicenter_parcel_dict.p', 'w') as f:
  pickle.dump(epicenter_parcel_dict, f)

epicenter_seedmap_dict_all = {parcel_index: nib.load('/data/mridata/jdeng/sd_bvftd/100_controls/1ST_vol_%s/spmT_0001.nii' % parcel_index).get_data() for parcel_index in range(1, 274)}
with open('epicenter_seedmap_dict_all.p', 'w') as f:
  pickle.dump(epicenter_seedmap_dict_all, f)

epicenter_seedmap_dict_all_v2 = {parcel_index: nib.load('/data/mridata/jdeng/sd_bvftd/100_controls_v2/1ST_vol_%s/spmT_0001.nii' % parcel_index).get_data() for parcel_index in range(1, 274)}
with open('epicenter_seedmap_dict_all_v2.p', 'w') as f:
  pickle.dump(epicenter_seedmap_dict_all_v2, f)
