import nibabel as nib
import numpy as np

# There are 482 parcels in v14, which will be stored in a dictionary
L = {x: nib.load('/data/mridata/jbrown/brains/brainnetome_suit_comb/vol_%d.nii' % x).get_data() for x in range(1, 273)}

def dicecoef(a, b):
    """Returns the Dice coefficient of a and b.
    
    Parameters
    ----------
    a, b : lists or arrays
    
    Returns
    -------
    out : float
        Dice coefficient of a and b
    """
    a_voxels = set(a)
    b_voxels = set(b)
    overlap = len(a_voxels & b_voxels)
    return round(overlap * 2.0/(len(a_voxels) + len(b_voxels)), 3)

def find_epicenter(wmap_path, w_thr, n_overlap, FC_thr, mask_path='/data/mridata/jbrown/brains/gm_mask/merged_ho_cereb_stn_comb.nii'):
    """Returns the parcel that is the epicenter for a subject.
    
    Parameters
    ----------
    wmap_path : string
        The absolute path to the wmap.
    w_thr : int
        The lower w-score threshold of the wmap, inclusive.
    n_overlap : int
        The minimum number of overlapping voxels between the wmap and a parcel for that parcel to be considered as an epicenter candidate.
    FC_thr : int
        The lower percentile threshold of the parcel-seeded mean functional connectivity seed maps, inclusive.
    mask_path : string, optional
        The absolute path to the image which will be used to mask the wmap.
        
    Returns
    -------
    out : string
        The indices of the parcels with the best-fitting seed maps to the wmap and their respective Dice coefficients.
    """
    
    # Mask and threshold wmap
    wmap = nib.load(wmap_path).get_data()
    thr = wmap.ravel() >= w_thr
    mask = nib.load(mask_path).get_data()
    mask = np.ma.make_mask(mask.ravel())
    mask_thr = np.logical_and(mask, thr)
    mask_thr_wmap = np.where(mask_thr == True)[0]
    
    # Find parcels that contain at least "n_overlap" voxels of overlap with the masked and thresholded wmap and keep their FC seed maps
    for i in L.keys():
        parcel = L[i].ravel()
        parcel = np.nonzero(parcel)[0]
        if len(set(parcel) & set(mask_thr_wmap)) >= n_overlap:
            # Replace parcel with seed map data
            L[i] = nib.load('/data/mridata/jdeng/sd_bvftd/100_controls/1ST_vol_%s/spmT_0001.nii' % i).get_data().ravel()
        else:
            # otherwise remove parcel from consideration as epicenter
            del L[i]
        
    # Threshold each seed map
    for i in L.keys():
        # FC_thr is percentile cutoff
        FC_thr_val = np.percentile(L[i][L[i] != 0], FC_thr)
        L[i] = np.nonzero(L[i] >= FC_thr_val)[0]

    # Get the best-fitting thresholded seed maps to the masked & thresholded wmap by taking the thresholded seed maps
    # with the greatest Dice coefficients to the wmap
    for i in L.keys():
        L[i] = dicecoef(mask_thr_wmap, L[i])
        
    # Separate L into a list of parcel IDs and a list of Dice coefficients, then sort them
    parcel_IDs = L.keys()
    dice_coefficients = L.values()
    sorted_indices = np.argsort(dice_coefficients)
    
    # top 10 epicenters and Dice coefficients
    pid = wmap_path.split("/")[4].split("_")[0]
    date = wmap_path.split("/")[4].split("_")[1]
    print '%s %s' % (pid, date), np.array(parcel_IDs)[sorted_indices][-1:-10:-1], np.array(dice_coefficients)[sorted_indices][-1:-10:-1]