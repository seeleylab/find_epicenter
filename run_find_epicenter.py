#!/usr/bin/env python
import sys
import find_epicenter as fe

def main():
    wmap_path=sys.argv[1]
    w_thr=float(sys.argv[2])
    n_overlap=int(sys.argv[3])
    FC_thr=float(sys.argv[4])
    fe.find_epicenter(wmap_path, w_thr, n_overlap, FC_thr)
if __name__ == "__main__":
    main()