import numpy as np
import sys
from scipy.sparse import bsr_matrix, csr_matrix, coo_matrix

from enum import Enum

class mat_type(Enum):
    ParCSR = 0
    ParBSR = 1

class mat_header:
    def __init__(self, arr):
        assert(arr.shape[0] == 4 or arr.shape[0] == 6)
        self.mat_type = mat_type(arr[0])
        self.shape = (arr[1], arr[2])
        self.nprocs = arr[3]
        if (self.mat_type is mat_type.ParBSR):
            self.block_shape = (arr[4], arr[5])


def parse(fname, b_size):
        offset = 0
        nrows = np.fromfile(fname, dtype='int32', count=1)[0]
        offset += 4
        rowptr = np.fromfile(fname, dtype='int32', offset=offset, count=nrows+1)
        offset += 4 * rowptr.shape[0]
        nnz = rowptr[-1]
        colind = np.fromfile(fname, dtype='int32', offset=offset, count=nnz)
        offset += 4 * colind.shape[0]
        values = np.fromfile(fname, dtype=np.double, offset=offset, count=nnz*b_size)
        return [rowptr, colind, values]

def read(base):
    info = mat_header(np.fromfile(f'{base}.hdr', dtype='int32'))
    print(f'{info.mat_type} on {info.nprocs} ranks')

    nprocs = info.nprocs
    b_size = np.prod(info.block_shape) if info.mat_type is mat_type.ParBSR else 1

    proc_data = [parse(f'{base}.{r}', b_size) for r in range(nprocs)]

    nrows = sum([p[0].shape[0] - 1 for p in proc_data])
    nnz = sum([p[1].shape[0] for p in proc_data])

    for r in range(nprocs - 1):
        proc_data[r+1][0] = proc_data[r+1][0] + proc_data[r][0][-1]

    rowptr = np.concat([proc_data[i][0] if i == 0 else proc_data[i][0][1:] for i in range(nprocs)])
    colind = np.concat([p[1] for p in proc_data])

    values = np.concat([p[2] for p in proc_data])
    if info.mat_type is mat_type.ParBSR:
        values = values.reshape(nnz, info.block_shape[0], info.block_shape[1])
        return bsr_matrix((values, colind, rowptr))
    elif info.mat_type is mat_type.ParCSR:
        return csr_matrix((values, colind, rowptr))

def read_hypre(base):
    import glob
    nprocs = len(glob.glob(f'{base}.*'))
    def parse_hypre(fname):
        data = np.loadtxt(fname, skiprows=1)
        return [np.array(data[:, 0], dtype=int),
                np.array(data[:, 1], dtype=int),
                np.array(data[:, 2], dtype=np.float64)]

    fnames = [f'{base}.{"0" * (5 - len(str(r))) + str(r)}' for r in range(nprocs)]
    proc_data = [parse_hypre(fname) for fname in fnames]
    I = np.concat([p[0] for p in proc_data])
    J = np.concat([p[1] for p in proc_data])
    data = np.concat([p[2] for p in proc_data])
    return coo_matrix((data, (I, J)))
