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
    """
    Returns the parcel that is the epicenter for a subject.
    
    Parameters
    ----------
    wmap_path : string
        The absolute path to the wmap.
    w_thr : int
        The lower threshold of the wmap, inclusive.
    n_overlap : int
        The minimum number of overlapping voxels between the wmap and a parcel for that parcel to be considered as an epicenter candidate.
    FC_thr : int
        The lower threshold of the parcel-seeded mean functional connectivity seed maps, inclusive. Appropriate lines should be commented out of the
        source code depending on whether the threshold is a percentile or an absolute z-score.
    mask_path : string, optional
        The absolute path to the image which will be used to mask the wmap.
        
    Returns
    -------
    out : string
        A statement about which parcel has the best-fitting seed map to the wmap and the Dice coefficient of the seed map and the wmap OR
        A statement about the top three parcels with the best-fitting seed maps to the wmap and their Dice coefficients.
        Appropriate lines should be commented out of the source code depending on which output is desired.
        
    Examples
    --------
    >>> fe.find_epicenter('/data5/patientNIC/3521/3521_20100114/struc/SPM12_SEG_Full/wmap_123SD_vs_288HC/wmap_2mm.nii', 1.5, 1, 0.5)
    3521 20100114 vol_158.nii 0.256000
    """
    # Find parcels that contain at least "n_overlap" voxels of overlap with the "mask"-masked "w_thr"-thresholded wmap and keep their FC seed maps
    wmap = nib.load(wmap_path).get_data()
    thr = wmap.ravel() >= w_thr
    mask = nib.load(mask_path).get_data()
    mask = np.ma.make_mask(mask.ravel())
    mask_thr = np.logical_and(mask, thr)
    mask_thr_wmap = np.where(mask_thr == True)[0]
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
        # uncomment if FC_thr is an absolute cutoff
        #L[i] = np.nonzero(L[i] >= FC_thr)[0]
        
        # uncomment if FC_thr is a percentile cutoff
        FC_thr_val = np.percentile(L[i][L[i] != 0], FC_thr)
        L[i] = np.nonzero(L[i] >= FC_thr_val)[0]

    # Get the best-fitting thresholded seed map(s) to the masked & thresholded wmap by taking the thresholded seed map with the greatest Dice coefficient to the wmap
    for i in L.keys():
        L[i] = dicecoef(mask_thr_wmap, L[i])
    
    L_z = L.copy()
    for i in L_z:
        L_z[i] = (L[i] - np.mean(L.values())) / np.std(L.values())
        
    pid = wmap_path.split("/")[4].split("_")[0]
    date = wmap_path.split("/")[4].split("_")[1]
    # uncomment if only the top epicenter and its Dice coefficient is wanted
    print('%s %s vol_%d.nii %f' % (pid, date, max(L, key=L.get), max(L.values())))    
    
    # # uncomment if top 3 epicenters and their Dice coefficients are wanted
    # first = sorted(L, key=L.get)[-1]; second = sorted(L, key=L.get)[-2]; third = sorted(L, key=L.get)[-3]
    # print('%s %s vol_%d.nii vol_%d.nii vol_%d.nii %f %f %f %f %f %f' % (pid, date, first, second, third, L[first], L[second], L[third], L_z[first], L_z[second], L_z[third]))
    
    # # uncomment if a file containing all Dice coefficients is wanted
    # np.savetxt("/".join(wmap_path.split("/")[:5])+'/struc/epicenters/dice_coefficients.txt', sorted(L_z.values()), fmt='%1.3f')
    
    # # addendum code to return top 4 epicenters for Jesse's group map
    # first = sorted(L, key=L.get)[-1]; second = sorted(L, key=L.get)[-2]; third = sorted(L, key=L.get)[-3]; fourth = sorted(L, key=L.get)[-4]
    # print('vol_%d.nii vol_%d.nii vol_%d.nii vol_%d.nii %f %f %f %f %f %f %f %f' % (first, second, third, fourth, L[first], L[second], L[third], L[fourth], L_z[first], L_z[second], L_z[third], L_z[fourth]))
