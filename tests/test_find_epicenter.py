from find_epicenter import percentile_threshold, mask_and_absolute_threshold, filter_parcels, create_epicenter_thr_seedmap_dict, find_epicenter
import nibabel as nib
import numpy as np
from numpy.testing.utils import assert_allclose, assert_equal

mask = '/data/mridata/jbrown/brains/gm_mask/merged_ho_cereb_stn_comb.nii'

def test_percentile_threshold():
	img = nib.load('100_avg_vol_1_seedmap.nii').get_data()
	result = percentile_threshold(img, 90)
	expected = nib.load('100_avg_vol_1_seedmap_pct_thr_90_using_fsl.nii').get_data().ravel()
	expected_indices = np.where(expected != 0)[0]
	assert_equal(result, expected_indices)

def test_mask_and_absolute_threshold():
	result = mask_and_absolute_threshold('wmap.nii', 2.0)
	expected = nib.load('wmap_gm_masked_abs_thr_at_2_using_fsl.nii').get_data().ravel()
	expected_indices = np.where(expected != 0)[0]
	assert_equal(result, expected_indices)

def test_filter_parcels():
	pass

def test_create_epicenter_thr_seedmap_dict():
	pass

def test_find_epicenter():
	pass