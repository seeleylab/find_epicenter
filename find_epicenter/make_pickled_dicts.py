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

epicenter_seedmap_dict_all_75hc = {parcel_index: nib.load('/data/mridata/jdeng/sd_bvftd/75_controls/1ST_vol_%s/spmT_0001.nii' % parcel_index).get_data() for parcel_index in range(1, 274)}
with open('epicenter_seedmap_dict_all_75hc.p', 'w') as f:
  pickle.dump(epicenter_seedmap_dict_all_75hc, f)

########
# After trimming epicenters from 273 to 252...
parcels_273 = range(1,274)
excluded_parcels = '''93 101 102 118 254 256 257 258 259 260 261 262 263 264 265 266 267 268 270 271 273'''.split()
excluded_parcels = [int(item) for item in excluded_parcels]
parcels_252_orig_indices = sorted(list(set(parcels_273) - set(excluded_parcels)))
parcels_252_new_indices = range(1,253)
lut = dict(zip(parcels_252_new_indices, parcels_252_orig_indices)) 

epicenter_parcel_dict_252_75hc = {parcel_index: nib.load('/data/mridata/jbrown/brains/brainnetome_suit_comb/vol_%s.nii'
                                    % lut[parcel_index]).get_data() for parcel_index in lut.keys()}

epicenter_parcelindices_dict_252_75hc = {k: np.where(v.flatten() > 0)[0] for
                                         k, v in epicenter_parcel_dict_252_75hc.items()}

with open('epicenter_parcelindices_dict_252_75hc.p', 'w') as f:
  pickle.dump(epicenter_parcelindices_dict_252_75hc, f)



epicenter_seedmap_dict_252_75hc = {parcel_index: nib.load('/data/mridata/jdeng/sd_bvftd/75_controls/1ST_vol_%s/spmT_0001.nii'
                                    % lut[parcel_index]).get_data() for parcel_index in lut.keys()}

with open('epicenter_seedmap_dict_252_75hc.p', 'w') as f:
  pickle.dump(epicenter_seedmap_dict_252_75hc, f)