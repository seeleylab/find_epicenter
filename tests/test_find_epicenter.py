from find_epicenter import mask_image, threshold_image, keep_if_overlap_by, find_epicenter
import nibabel as nib
from numpy.testing.utils import assert_allclose, assert_equal

mask = '/data/mridata/jbrown/brains/gm_mask/merged_ho_cereb_stn_comb.nii'

def test_mask_image():
	img2mask = 'wmap.nii'
	result = mask_image(img2mask, mask)
	expected = nib.load('wmap_gm_masked_using_fsl.nii').get_data()
	expected = expected[expected != 0].ravel()
	assert_equal(result, expected)

def test_threshold_image_with_float_threshold_level():
	img2threshold = mask_image('wmap.nii', mask)
	result = threshold_image(img2threshold, 2.0)
	expected = nib.load('wmap_gm_masked_thr_using_fsl.nii').get_data()
	expected = expected[expected != 0].ravel()
	assert_equal(result, expected)

def test_threshold_image_with_int_threshold_level():
	img2threshold = mask_image('wmap.nii', mask)
	result = threshold_image(img2threshold, 2)
	expected = nib.load('wmap_gm_masked_thr_using_fsl.nii').get_data()
	expected = expected[expected != 0].ravel()
	assert_equal(result, expected)


