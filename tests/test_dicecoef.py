from find_epicenter import dicecoef
import numpy as np
import nibabel as nib
import random
from numpy.testing.utils import assert_allclose, assert_equal, assert_raises

def test_identical_arrays():
	a = np.array([-5.0, 0.0, 2.5, 10.0])
	b = np.array([-5.0, 0.0, 2.5, 10.0])
	result = dicecoef(a, b)
	expected = 1.0
	assert_equal(result, expected)

def test_half_same_arrays():
	a = np.array([234, 74, 980, 28])
	b = np.array([483, 74, 980, 4281])
	result = dicecoef(a, b)
	expected = 0.5
	assert_equal(result, expected)

def test_different_arrays():
	np.random.seed(42)
	a = np.random.rand(100)
	b = np.negative(a)  
	result = dicecoef(a, b)
	expected = 0
	assert_equal(result, expected)

def test_one_input_is_empty_array():
	a = np.array([])
	np.random.seed(42)
	b = np.random.rand(100)
	result = dicecoef(a, b)
	expected = 0
	assert_equal(result, expected)

def test_both_inputs_are_empty_arrays():
	a = np.array([])
	b = np.array([])	
	assert_raises(AssertionError, dicecoef, a, b) 

def test_identical_lists():
	a = [-5.0, 0.0, 2.5, 10.0]
	b = [-5.0, 0.0, 2.5, 10.0]
	result = dicecoef(a, b)
	expected = 1.0
	assert_equal(result, expected)

def test_half_same_lists():
	a = [234, 74, 980, 28]
	b = [483, 74, 980, 4281]
	result = dicecoef(a, b)
	expected = 0.5
	assert_equal(result, expected)

def test_different_lists():
	random.seed(42)
	a = random.sample(range(1, 101), 100)
	b = [-x for x in a]
	result = dicecoef(a, b)
	expected = 0
	assert_equal(result, expected)

def test_one_list_is_empty_array():
	a = []
	b = [1, 2, 3, 4]
	result = dicecoef(a, b)
	expected = 0
	assert_equal(result, expected)

def test_both_inputs_are_empty_arrays():
	a = []
	b = []
	assert_raises(AssertionError, dicecoef, a, b)

def test_identical_wmaps():
#	a = nib.load('wmap_for_testing.nii').get_data()
#	a_indices = np.where(a >= 2) 
#	b = nib.load('wmap_for_testing_plus_1.nii').get_data()
#	b_indices = np.where(b >= 2)
#	result = dicecoef(a, b)
#	expected = 1.0
#	assert_equal(result, expected)
	pass

def test_different_wmaps():
	pass

def test_str_and_str():
	a = 'meow'
	b = 'mrow'
	assert_raises(TypeError, dicecoef, a, b)

def test_array_and_str():
	a = np.array([2, 4, 6, 8])
	b = 'who do we appreciate?'
	assert_raises(TypeError, dicecoef, a, b)

def test_list_and_str():
	a = [2, 4, 6, 8]
	b = 'who do we appreciate?'
	assert_raises(TypeError, dicecoef, a, b)	
