from find_epicenter import return_indices_of_masked_thr_voxels
import nibabel as nib
import numpy as np
from numpy.testing.utils import assert_allclose, assert_equal

mask = '/data/mridata/jbrown/brains/gm_mask/merged_ho_cereb_stn_comb.nii'

def test_return_indices_of_masked_thr_voxels():
	img2process = 'wmap.nii'
	result = return_indices_of_masked_thr_voxels(img2process, 2.0)
	expected = nib.load('wmap_gm_masked_thr_using_fsl.nii').get_data().ravel()
	expected_indices = np.where(expected != 0)[0]
	assert_equal(result, expected_indices)

