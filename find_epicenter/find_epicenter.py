import nibabel as nib
import numpy as np
import pickle
from multiprocessing import Lock, Pool
import sys

def mask_and_absolute_threshold(img_path, absolute_threshold_level, \
	mask_path='find_epicenter/merged_ho_cereb_stn_comb.nii'):
	"""Return indices of voxels in an image that are within a mask and greater
	than or equal to an absolute threshold.
	
	Parameters
	----------
	img_path : str
		Absolute path to an image.
	threshold_level : int or float
		Absolute cutoff for which any voxels with values greater than or equal
		to this value will be kept; otherwise the voxels will be zeroed out.
	mask_path : str
		Absolute path to a mask.
		
	Returns
	-------
	in_mask_and_above_threshold_voxels_indices : 1darray
		Indices of voxels in image that are within the mask and greater than or
		equal to the threshold.
	"""
	raw_image = nib.load(img_path).get_data().ravel()
	mask = nib.load(mask_path).get_data().ravel()
	in_mask_voxels_boolean = np.ma.make_mask(mask)
	above_threshold_voxels_boolean = (raw_image >= absolute_threshold_level)
	in_mask_and_above_threshold_voxels_boolean = in_mask_voxels_boolean & \
			above_threshold_voxels_boolean
	in_mask_and_above_threshold_voxels_indices = \
			np.where(in_mask_and_above_threshold_voxels_boolean)[0]
	return in_mask_and_above_threshold_voxels_indices

def filter_parcels(img_indices, min_overlap):
	"""Keep the parcel as an epicenter candidate if it shares at least the
	specified number of voxels of overlap with the image.
	
	Parameters
	----------
	img_indices : 1darray
		Indices of voxels in the image.
	min_overlap : int or float
		Defines cutoff for which a parcel that overlaps the image by any number
		of voxels greater than or equal to this cutoff will be kept.
	epicenter_parcel_dict : pickled dict

	Returns
	-------
	epicenter_candidates : list
		List of candidate epicenters.
	"""
	epicenter_candidates = []
	with open('find_epicenter/epicenter_parcelindices_dict.p', 'r') as f:
		epicenter_parcel_dict = pickle.load(f)
	
	for i in epicenter_parcel_dict.keys():
		parcel = epicenter_parcel_dict[i]
		if len(set(parcel) & set(img_indices)) >= min_overlap:
			epicenter_candidates.append(i)
	
	return epicenter_candidates

def percentile_threshold(img_data, percentile_threshold_level):
	"""Return indices of voxels in an image that are greater than or equal
	to a percentile threshold.
	Parameters
	----------
	img_data : ndarray
		Values of every voxel in an image.
	threshold_level : int or float
		Percentile cutoff for which any voxels with values greater than or
		equal to this value will be kept; otherwise the voxels will be
		zeroed out.
		
	Returns
	-------
	thresholded_image_indices: 1darray
		Indices of voxels in image that are greater than or equal to the
		threshold.
	"""
	nonzero_img_data = img_data[img_data != 0]
	img_threshold = np.percentile(nonzero_img_data, percentile_threshold_level)
	thresholded_image = (img_data >= img_threshold)
	thresholded_image_indices = np.where(thresholded_image.ravel())[0]
	return thresholded_image_indices

def create_epicenter_thr_seedmap_dict(candidates_list,
									  percentile_threshold_level):
	"""Return a dictionary of key-value pairs where keys are epicenters and
	values are the indices of voxels in the epicenter-seeded functional
	connectivity map that are greater than or equal to the threshold.
	
	Parameters
	----------
	candidates_list : list
		Epicenters.
	percentile_threshold_level : int or float
		Percentile cutoff for which any voxels greater than or equal to this
		value will be kept; otherwise the voxels will be zeroed out.
		
	Returns
	-------
	epicenter_thr_seedmap_dict : dict
		Epicenters and the indices of voxels in the epicenter-seeded functional
		connectivity map that are greater than or equal to the threshold.
	"""
	with open('find_epicenter/epicenter_seedmap_dict_all_v2.p', 'r') as f:
		epicenter_seedmap_dict_all_v2 = pickle.load(f)
	epicenter_seedmap_dict = {i: epicenter_seedmap_dict_all_v2[i] for i in\
							  candidates_list}
	
	epicenter_thr_seedmap_dict = \
	{i: percentile_threshold(epicenter_seedmap_dict[i],
	percentile_threshold_level) for i in epicenter_seedmap_dict.keys()}
	
	return epicenter_thr_seedmap_dict
	
def dicecoef(a, b):
	"""Return the Dice coefficient of a and b.

	Parameters
	----------
	a, b : lists or arrays.
	
	Returns
	-------
	out : float
		Dice coefficient of a and b.
	"""
	if type(a) == str or type(b) == str:
		raise(TypeError, 'This implementation of the Dice coefficient does not '
			  'handle strings.')
	else:
		a_voxels = set(a)
		b_voxels = set(b)
		overlap = len(a_voxels & b_voxels)
		total = len(a_voxels) + len(b_voxels)
		assert total > 0
		dc = (overlap * 2.0) / total
		return round(dc, 3)

lock = Lock()

def write_top_epicenters(img_path, img_indices, epicenter_thr_seedmap_dict,\
						 outfile):
	"""Find the Dice coefficient between the image and each epicenter's
	masked+thresholded functional connectivty map and print the top
	epicenters with the highest Dice coefficients.
	
	Parameters
	----------
	img_path : str
		Absolute path to an image.
	img_indices : 1darray
		Indices of voxels in the image.
	epicenter_thr_seedmap_dict : dict
		Epicenters and the indices of voxels in the epicenter-seeded functional
		connectivity map that are greater than or equal to the threshold.
		
	Returns
	-------
	out : str
		Sorted epicenters and the Dice coefficients of their functional
		connectivity maps to the image.
	"""
	epicenter_dicecoef_dict = {i: dicecoef(img_indices, \
	epicenter_thr_seedmap_dict[i]) for i in epicenter_thr_seedmap_dict.keys()}

	epicenters = epicenter_dicecoef_dict.keys()
	dice_coefs = epicenter_dicecoef_dict.values()
	sorted_indices = np.argsort(dice_coefs)
	lock.acquire()
	with open(outfile, 'a') as f:
		f.write('%s %s %s' % (img_path,
							  np.array(epicenters)[sorted_indices][-1:-11:-1],
							  np.array(dice_coefs)[sorted_indices][-1:-11:-1]))
		f.write('\n')
	lock.release()
	
def find_epicenter(subj):
	"""Return the top epicenters.
	
	Parameters
	----------
	img_path : str
		Absolute path to an image.
		
	Returns
	-------
	out : str
		Top epicenters and Dice coefficients.
	"""
	mask_thr_indices = mask_and_absolute_threshold(img_path=subj,
					absolute_threshold_level=2.0)
	epicenter_candidates = filter_parcels(img_indices=mask_thr_indices,
			min_overlap=10)
	epicenter_thr_seedmap_dict = \
	create_epicenter_thr_seedmap_dict(candidates_list=epicenter_candidates,
				percentile_threshold_level=90)
	write_top_epicenters(img_path=subj, img_indices=mask_thr_indices,\
					epicenter_thr_seedmap_dict=epicenter_thr_seedmap_dict,\
					outfile='/data/mridata/jdeng/sd_bvftd/v3_new_registration/bvsd_epicenters_top10.txt')

if __name__ == '__main__':
	with open(sys.argv[1], 'r') as f:
		subjs = [line.strip() for line in f]
	pool = Pool(processes=10)
	pool.map(find_epicenter, subjs)
	pool.close()
	pool.join()