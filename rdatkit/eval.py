from numpy import *
from random import *

if __package__ is None or not __package__:
    from . import secstr
else:
    from . import secstr


def _bppm_thresh_indices(bppm_pred, bppm_true, thresh):
    return logical_or(bppm_pred >= thresh, bppm_true >= thresh)


def get_helices_from_structures(structures):
    helices = []
    for struct in structures:
        obj = secstr.SecondaryStructure(dbn=struct) if isinstance(struct, str) else struct

        struct_helices = obj.helices()
        for struct_helix in struct_helices:
            add_helix = True
            for i, helix in enumerate(helices):
                count = 0
                for bp in struct_helix:
                    if bp in helix:
                        count += 1

                if count >= len(helix) / 2 or count >= len(struct_helix) / 2:
                    if len(struct_helix) > len(helix):
                        helices[i] = struct_helix
                    else:
                        add_helix = False
            if add_helix:
                helices.append(struct_helix)

    return helices


def _find_helices_from_bppm(bppm):
    visited = []
    helices = []
    for i in range(bppm.shape[0]):
        for j in range(i + 1, bppm.shape[1]):
            if (i, j) not in visited and bppm[i, j] != 0:
                helix = [(i, j)]
                ii = i + 1
                jj = j - 1
                while ii < bppm.shape[0] and jj >= i and (ii, jj) not in visited and bppm[ii, jj] != 0:
                    helix.append((ii, jj))
                    visited.append((ii, jj))
                    ii += 1
                    jj -= 1

                ii = i - 1
                jj = j + 1
                while ii <= 0 and jj < bppm.shape[1] and (ii, jj) not in visited and bppm[ii, jj] != 0:
                    helix.append((ii, jj))
                    visited.append((ii, jj))
                    ii -= 1
                    jj += 1
                helices.append(helix)
    return helices


def get_indv_bppm_tp_fp_tn_fn(bppm_pred, bppm_true, thresh, diff_thresh=None, thresh2=None, helices=False):
    if thresh2 is None:
        thresh1 = thresh

    if helices:
        (tp, fp, tn, fn) = (0., 0., 0., 0.)
        helices = _find_helices_from_bppm(bppm_true >= thresh2)
        helices_pred = _find_helices_from_bppm(bppm_pred >= thresh2)
        for h in helices:
            if len(h) < 3: continue
            (bps_tp, bps_tn, bps_fn) = (0., 0., 0.)
            for n1, n2 in h:
                if bppm_pred[n1, n2] >= thresh or bppm_true[n1, n2] >= thresh2:
                    if bppm_pred[n1, n2] >= thresh:
                        if bppm_true[n1, n2] >= thresh2:
                            bps_tp += 1
                    else:
                        if bppm_true[n1, n2] >= thresh2:
                            bps_fn += 1
                        else:
                            bps_tn += 1
            if bps_tp >= len(h) / 2:
                tp += 1
            if bps_tn >= len(h) / 2:
                tn += 1
            if bps_fn > len(h) / 2:
                fn += 1

        for h in helices_pred:
            bps_fp = 0.
            if len(h) < 3:
                continue
            for n1, n2 in h:
                if bppm_pred[n1, n2] >= thresh:
                    if bppm_pred[n1, n2] >= thresh and bppm_true[n1, n2] < thresh2:
                        bps_fp += 1
            if bps_fp > len(h)/2:
                fp += 1
        return tp, fp, tn, fn

    else:
        true_indices = bppm_true >= thresh2
        pred_indices = bppm_pred >= thresh
        not_pred_indices = logical_not(pred_indices)
        not_true_indices = logical_not(true_indices)
        #tp = (abs(bppm_pred[pred_indices] - bppm_true[pred_indices]) <= diff_thresh).sum()
        #fp = (abs(bppm_pred[pred_indices] - bppm_true[pred_indices]) > diff_thresh).sum()

        #tn = logical_and(not_pred_indices, not_true_indices).sum()
        #fn = (abs(bppm_pred[not_pred_indices] - bppm_true[not_pred_indices]) > diff_thresh).sum()

        tp = logical_and(pred_indices, true_indices).sum() / 2.
        fp = logical_and(pred_indices, not_true_indices).sum() / 2.
        tn = logical_and(not_pred_indices, not_true_indices).sum() / 2.
        fn = logical_and(not_pred_indices, true_indices).sum() / 2.
    return (tp, fp, tn, fn)


def get_bppm_tp_fp_tn_fn(bppm_pred, bppm_true, thresh, diff_thresh=None, thresh2=None):
    if thresh2 is None:
        thresh2 = thresh
    true_indices = _bppm_thresh_indices(bppm_true, bppm_true, thresh)
    pred_indices = _bppm_thresh_indices(bppm_pred, bppm_pred, thresh2)
    not_pred_indices = logical_not(pred_indices)
    not_true_indices = logical_not(true_indices)
    #tp = (abs(bppm_pred[pred_indices] - bppm_true[pred_indices]) <= diff_thresh).sum()
    #fp = (abs(bppm_pred[pred_indices] - bppm_true[pred_indices]) > diff_thresh).sum()

    #tn = logical_and(not_pred_indices, not_true_indices).sum()
    #fn = (abs(bppm_pred[not_pred_indices] - bppm_true[not_pred_indices]) > diff_thresh).sum()

    tp = logical_and(pred_indices, true_indices).sum()
    fp = logical_and(pred_indices, not_true_indices).sum()
    tn = logical_and(not_pred_indices, not_true_indices).sum()
    fn = logical_and(not_pred_indices, true_indices).sum()
    return (tp, fp, tn, fn)


def bpp_rmsd(bppm_pred, bppm_true, thresh=0.01, indices=None):
    if indices is None:
        indices = _bppm_thresh_indices(bppm_pred, bppm_true, thresh)
    return sqrt(((bppm_pred[indices] - bppm_true[indices]) ** 2).mean())


def _find_mean_helix_value(helixm, i, j, indices, return_visited=False):
    ii = i - 1
    jj = j + 1
    visited = [(i, j)]
    while True:
        if ii < 0 or jj >= helixm.shape[1] or not indices[ii, jj] or helixm[ii, jj] == 0:
            break
        visited.append((ii, jj))
        ii -= 1
        jj += 1

    ii = i + 1
    jj = j - 1
    while True:
        if  ii >= helixm.shape[0] or jj < 0 or not indices[ii, jj] or helixm[ii, jj] == 0 :
            break
        visited.append((ii, jj))
        ii += 1
        jj -= 1

    mean_val = array([helixm[i, j] for i, j in visited]).min()
    if return_visited:
        return mean_val, visited

    return mean_val


def helix_rmsd(bppm_pred, bppm_true, helices, thresh=0.01, indices=None):
    rmsd = 0.
    if indices is None:
        indices = _bppm_thresh_indices(bppm_pred, bppm_true, thresh)

    def _get_helix_vals(matrix, helix):
        return array([matrix[n1, n2] for n1, n2 in helix])

    for helix in helices:
        hp_true = _get_helix_vals(bppm_true, helix).mean()
        hp_pred = _get_helix_vals(bppm_pred, helix).mean()
        rmsd += (hp_true - hp_pred) ** 2

    rmsd = sqrt(rmsd / len(helices))
    return rmsd
"""
    helix_errm = zeros(bppm_true.shape)
    helix_errm[indices] = (bppm_true[indices] - bppm_pred[indices])**2
    visited = []
    err = 0.
    count = 0.
    for i in xrange(helix_errm.shape[0]):
        for j in xrange(i+1, helix_errm.shape[1]):
            if helix_errm[i,j] != 0 and (i,j) not in visited:
                val, helix_visited = _find_mean_helix_value(helix_errm, i, j, indices, return_visited=True)
                err += val
                count += 1.
                visited += helix_visited
    if count == 0:
        return 0.
    else:
        return sqrt(err/count)
"""


def bpp_ppv(bppm_pred, bppm_true, thresh=0.01, diff_thresh=0.1):
    (tp, fp, tn, fn) = get_bppm_tp_fp_tn_fn(bppm_pred, bppm_true, thresh, diff_thresh)
    return tp / float(tp + fp)


def bpp_sensitivity(bppm_pred, bppm_true, thresh=0.01, diff_thresh=0.1):
    (tp, fp, tn, fn) = get_bppm_tp_fp_tn_fn(bppm_pred, bppm_true, thresh,  diff_thresh)
    return tp/float(tp + fn)


def bpp_roc(bppm_pred, bppm_true, interval=None):
    if interval is None:
        max_thresh = min(bppm_pred.max(), bppm_true.max())
        interval = arange(0.01, max_thresh, 0.1)
    ppv = []
    sens = []
    for thresh in interval:
        ppv.append(bpp_ppv(bppm_pred, bppm_true, diff_thresh=thresh))
        sens.append(bpp_sensitivity(bppm_pred, bppm_true, diff_thresh=thresh))
    sorted_indices = [i for i, p in sorted(enumerate(ppv), key=lambda x: x[1])]
    return [ppv[i] for i in sorted_indices], [sens[i] for i in sorted_indices]


def auc(ppv, sens):
    res = 0
    sorted_indices = [i for i, p in sorted(enumerate(ppv), key=lambda x: x[1])]
    for idx, i in enumerate(sorted_indices[:-1]):
        j = sorted_indices[idx + 1]
        base = (ppv[j] - ppv[i])
        h = max(sens[j], sens[i])
        res += base * h - (base * (h - min(sens[j], sens[i])) / 2.)
    return res
