import nibabel as nib
import numpy as np

# There are 482 parcels in v14
parcel_list = np.zeros((91, 109, 91, 482))
for i in range(0,482):
    parcel_list[:,:,:,i] = nib.load('/data/mridata/jbrown/hcp/v14/vol_%d.nii' %(i+1)).get_data()

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

def find_epicenter(wmap_path, w_thr, n_overlap, FC_thr, mask_path='/data/mridata/jbrown/brains/merged_ho_cereb_stn_max_bin.nii'):
    """Returns the parcel that is the epicenter for a subject.
    
    Parameters
    ----------
    wmap_path : string
        The absolute path to the wmap.
    w_thr : int
        The lower threshold of the wmap.
    n_overlap : int
        The minimum number of overlapping voxels between the wmap and a parcel for that parcel to be considered as an epicenter candidate.
    FC_thr : int
        The lower threshold of the parcel-seeded mean functional connectivity seed maps. Appropriate lines should be commented out of the source
        code depending on whether the threshold is a percentile or an absolute z-score.
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
    # Find parcels that contain at least "n_overlap" voxels of overlap with the "mask"-masked "w_thr"-thresholded wmap
    # and keep their FC seed maps
    wmap = nib.load(wmap_path).get_data()
    thr = wmap.ravel() >= w_thr
    mask = nib.load(mask_path).get_data()
    mask = np.ma.make_mask(mask.ravel())
    mask_thr = np.logical_and(mask, thr)
    mask_thr_wmap = np.where(mask_thr == True)[0]
    seed_maps = []
    epicenter_candidates = []
    for i in range(0,482):
        parcel = parcel_list[:,:,:,i].ravel()
        parcel = np.nonzero(parcel)[0]
        if len(set(parcel) & set(mask_thr_wmap)) >= n_overlap:
            epicenter_candidates.append(i)
            seed_map = nib.load('/data/mridata/jbrown/hcp/v14/vol_%d_mean_fc_v2.nii' %(i+1)).get_data()
            seed_map = seed_map.ravel()
            seed_maps.append(seed_map)
    
    # Threshold each seed map
    thr_seed_maps = []
    for i in range(0, len(seed_maps)):
        seed_map = seed_maps[i]
        #thr_seed_maps.append(np.nonzero(seed_map >= FC_thr)[0])     # uncomment if FC_thr is an absolute cutoff
        FC_thr_val = np.percentile(seed_map[seed_map != 0], FC_thr)     # uncomment if FC_thr is a percentile cutoff
        thr_seed_maps.append(np.nonzero(seed_map >= FC_thr_val)[0])     # uncomment if FC_thr is a percentile cutoff

    # Get the best-fitting thresholded seed map to the masked & thresholded wmap by taking the seed map with the greatest Dice coefficient to the wmap
    numbers = np.zeros(len(thr_seed_maps))
    
    for i in range(0,len(numbers)):
        numbers[i] = dicecoef(mask_thr_wmap, thr_seed_maps[i])
    
    pid = wmap_path.split("/")[4].split("_")[0]
    date = wmap_path.split("/")[4].split("_")[1]
    #print('%s %s vol_%d.nii %f' % (pid, date, epicenter_candidates[np.argmax(numbers)]+1, numbers[np.argmax(numbers)]))    # uncomment if only the
                                                                                                                            # top epicenter and its
                                                                                                                            # Dice coefficient is
                                                                                                                            # wanted
                                                                                                                            
    # numbers_z = (numbers - np.mean(numbers)) / np.std(numbers)                                                             # uncomment if top 3
    # first = np.argsort(numbers)[-1]; second = np.argsort(numbers)[-2]; third = np.argsort(numbers)[-3]                     # epicenters and their 
    # print('%s %s vol_%d.nii vol_%d.nii vol_%d.nii %f %f %f %f %f %f' % (pid, date, epicenter_candidates[first]+1,          # their Dice coefficients
    # epicenter_candidates[second]+1, epicenter_candidates[third]+1, numbers[first], numbers[second], numbers[third],        # are
    # numbers_z[first], numbers_z[second], numbers_z[third]))                                                                # wanted
    
    numbers_z = (numbers - np.mean(numbers)) / np.std(numbers)                                                           # uncomment if a file
    numbers_z_sorted = numbers_z[np.argsort(numbers_z)]                                                                  # containing all Dice
    np.savetxt("/".join(wmap_path.split("/")[:5])+'/struc/epicenters/dice_coefficients.txt', numbers_z_sorted, fmt='%1.3f') # coefficients is wanted