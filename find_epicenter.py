import nibabel as nib
import numpy as np
import pickle

def threshold_image_using_percentile_threshold(img_data, percentile_threshold_level):
	"""Threshold image using a percentile threshold.
	Parameters
	----------
	img_data : ndarray
		Image in the form of a 1D array of voxel values.
	threshold_level : int or float
		Defines cutoff for which any voxels greater than or equal to this value will be kept; otherwise the voxels will be zeroed out.
	Returns
	-------
	thresholded_image_indices: ndarray
		Indices of thresholded image as a 1D array of voxel values.
	"""
	nonzero_img_data = img_data[img_data != 0]
	img_threshold = np.percentile(nonzero_img_data, percentile_threshold_level)
	thresholded_image = img_data >= img_threshold
	thresholded_image_indices = np.where(thresholded_image.ravel())[0]
	return thresholded_image_indices

def return_indices_of_masked_thr_voxels(img_path, threshold_level, mask_path='/data/mridata/jbrown/brains/gm_mask/merged_ho_cereb_stn_comb.nii'):
	"""
	"""
	raw_image = nib.load(img_path).get_data().ravel()
	mask = nib.load(mask_path).get_data().ravel()
	in_mask_voxels_boolean = np.ma.make_mask(mask)
	above_threshold_voxels_boolean = raw_image > threshold_level
	in_mask_and_above_threshold_voxels_boolean = in_mask_voxels_boolean & above_threshold_voxels_boolean
	in_mask_and_above_threshold_voxels_indices = np.where(in_mask_and_above_threshold_voxels_boolean)[0]
	return in_mask_and_above_threshold_voxels_indices

def filter_epicenters_by_parcel_overlap(img_indices, min_overlap):
	"""Keep the parcel as an epicenter candidate if it shares at least the specified number of voxels of overlap with the image.
	
	Parameters
	----------
	img_indices : ndarray
		Indices in the image.
	min_overlap : int or float
		Defines cutoff for which a parcel that overlaps the image by any number of voxels greater than or equal to this cutoff will be kept.

	Returns
	-------
	epicenter_candidates : list
		List of candidate epicenter indices.
	"""
	epicenter_candidates = []
	epicenter_parcel_dict = pickle.load(open('/data/mridata/jdeng/tools/find_epicenter/epicenter_parcel_dict.p', 'r'))
	
	for i in epicenter_parcel_dict.keys():
		parcel = epicenter_parcel_dict[i].ravel()
		parcel_indices = np.nonzero(parcel)[0]
		if len(set(parcel_indices) & set(img_indices)) >= min_overlap:
			epicenter_candidates.append(i)
	
	return epicenter_candidates

def create_epicenter_thr_seedmap_dict(candidates_list, percentile_threshold_level):
	"""
	"""
	epicenter_seedmap_dict_all = pickle.load(open('/data/mridata/jdeng/tools/find_epicenter/epicenter_seedmap_dict.p', 'r'))
	epicenter_seedmap_dict = {i: epicenter_seedmap_dict_all[i] for i in candidates_list}
	
	epicenter_thr_seedmap_dict = {i: threshold_image_using_percentile_threshold(epicenter_seedmap_dict[i], percentile_threshold_level) for i in epicenter_seedmap_dict.keys()}
	
	return epicenter_thr_seedmap_dict
	
def dicecoef(a, b):
	"""Return the Dice coefficient of a and b.

	Parameters
	----------
	a, b : lists or arrays
	
	Returns
	-------
	out : float
		Dice coefficient of a and b.
	"""
	if type(a) == str or type(b) == str:
		raise(TypeError, 'This implementation of the Dice coefficient does not handle strings.')
	else:
		a_voxels = set(a)
		b_voxels = set(b)
		overlap = len(a_voxels & b_voxels)
		total = len(a_voxels) + len(b_voxels)
		assert total > 0
		return round((overlap * 2.0) / total, 3)

def find_epicenter(img_indices, epicenter_thr_seedmap_dict):
	"""
	"""
	epicenter_dicecoef_dict = {i: dicecoef(img_indices, epicenter_thr_seedmap_dict[i]) for i in epicenter_thr_seedmap_dict.keys()}

	epicenters = epicenter_dicecoef_dict.keys()
	dice_coefs = epicenter_dicecoef_dict.values()
	sorted_indices = np.argsort(dice_coefs)

	print '%s %s' % (np.array(epicenters)[sorted_indices][-1:-11:-1], np.array(dice_coefs)[sorted_indices][-1:-11:-1])

if __name__ == '__main__':
	mask_thr_indices = return_indices_of_masked_thr_voxels('tests/wmap.nii', 2.0)
	epicenter_candidates = filter_epicenters_by_parcel_overlap(mask_thr_indices, 10)
	epicenter_thr_seedmap_dict = create_epicenter_thr_seedmap_dict(epicenter_candidates, 90)
	find_epicenter(mask_thr_indices, epicenter_thr_seedmap_dict)