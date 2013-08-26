import numpy as np
import sys

def load_binary(fname):

    hdr = np.dtype(
        {'names':
            ['magic',
             'major_version',
             'minor_version',
             'sizeof_header',
             'creation_time',
             'sizeof_pos_t',
             'sizeof_vel_t',
             'sizeof_energy_t',
             'sizeof_id',
             'simulation_step',
             'simulation_time',
             'output_index',
             'Nparticles',
             'sizeof_options'],
         'formats':
            ['S8', '<u2', '<u2', '<u4', '<u8', '<u1',
             '<u1', '<u1','<u1','<u8', 'f8', '<u8', '<u8', '<u8']},
        align = True)

    H = np.memmap(fname, hdr, mode='c', offset=0, order='C', shape=1)[0]

    particle = np.dtype( 
        {'names':
            ['x', 'v', 'E', 'id'],
         'formats':
            ['3f%i' % H['sizeof_pos_t'],
             '3f%i' % H['sizeof_vel_t'],
             '1f%i' % H['sizeof_energy_t'], '1i%i' % H['sizeof_id']]}, align = False)
    D = np.memmap(fname, particle, mode='c', offset=H['sizeof_header'], order='C')
    return {'header': H, 'particles': D}

if __name__ == '__main__':
    import sys
    D = load_binary(sys.argv[1])
    print D['header']
    print D['particles']
