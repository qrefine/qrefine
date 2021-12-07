import numpy as np
import zmq
import struct
from ase.calculators import calculator
from ase import units

default_socket = '/tmp/aniserver.sock'
in_header_fmt = 'I'
out_header_fmt = 'I'

import logging
logging.basicConfig(level=logging.INFO)


class ANIFrontEnd(object):
    default_pbc = np.array([False, False, False])
    default_cell = np.eye(3, dtype=np.float32)

    def __init__(self, socket):
        self.socket = socket
        self.msgpack = MSGPack(in_header_fmt, out_header_fmt)

    def __call__(self, flag, arrays):
        arrays = self._prepare_in(arrays)
        in_msg_parts = self.msgpack.pack_input((flag,), arrays)
        self.socket.send_multipart(in_msg_parts, copy=False)
        out_msg_parts = self.socket.recv_multipart(copy=False)
        out_flag, arrays = self.msgpack.unpack_output(out_msg_parts)
        assert out_flag[0] == 0, 'Calculation failed.'
        return arrays

    def _prepare_in(self, arrays):
        R, Z = arrays[:2]
        R = R.astype(np.float32, copy=False)
        Z = Z.astype(np.uint8, copy=False)
        natom = Z.shape[-1]
        R = R.reshape(-1, natom, 3)
        Z = Z.reshape(-1, natom)
        nbatch = R.shape[0]
        assert Z.shape[0] == nbatch

        if len(arrays) > 2:
            assert len(arrays) == 4
            pbc = arrays[2]
            pbc = pbc.astype(np.bool, copy=False)
            pbc = pbc.reshape(3)
            cell = arrays[3]
            cell = cell.astype(np.float32, copy=False)
            cell = cell.reshape(3, 3)
        else:
            pbc = self.default_pbc
            cell = self.default_cell

        return R, Z, pbc, cell


class ANIRPCCalculator(calculator.Calculator):
    implemented_properties = ['energy', 'forces', 'stress', 'free_energy']

    def __init__(self, address='ipc://{}'.format(default_socket)):
        super(ANIRPCCalculator, self).__init__()
        ctx = zmq.Context()
        logging.debug('Opening connection with %s..', address)
        socket = ctx.socket(zmq.REQ)
        rc = socket.connect(address)
        if rc:
            raise zmq.error.ZMQError('Connection failed.')
        logging.debug('Done.')
        self.interface = ANIFrontEnd(socket)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=calculator.all_changes):
        super(ANIRPCCalculator, self).calculate(
            atoms, properties, system_changes)

        calc_type = 0
        if 'forces' in properties:
            calc_type = 1
        if 'stress' in properties:
            calc_type = 2

        R = self.atoms.get_positions()[None, ...].astype(np.float32)
        Z = self.atoms.get_atomic_numbers()[None, ...].astype(np.uint8)
        pbc = self.atoms.get_pbc().astype(np.bool)
        cell = self.atoms.get_cell(complete=True).astype(np.float32)

        results = self.interface(calc_type, (R, Z, pbc, cell))
        self.results['energy'] = self.results['free_energy'] = results[0][0] * units.Hartree
        if calc_type == 1:
            self.results['forces'] = results[1][0] * units.Hartree
        if calc_type == 2:
            self.results['stress'] = results[2] * units.Hartree


class ArrPack(object):
    """ Protocol to pack/unpack numpy arrays.
    Only simple dtype allowed, array.dtype.num
    Up to 6 dimensions
    Number of elements (total and in any dimension) is p to 2^32 (uint32)
    Change `header_fmt` to overwrite these limits
    """

    header_fmt = 'I' * 8  # dtype, numel, 6*shape

    def __init__(self):
        self.packer = struct.Struct(self.header_fmt)

    def pack_header(self, array):
        dtype = array.dtype.num
        numel = array.size
        shape = [0] * 6
        shape[-array.ndim:] = array.shape
        return self.packer.pack(dtype, numel, *shape)

    def unpack_header(self, data):
        values = self.packer.unpack(data)
        dtype, numel = values[:2]
        shape = np.trim_zeros(values[2:], 'f')
        return dtype, numel, shape

    def pack_array(self, array):
        return self.pack_header(array), array.data

    def unpack_array(self, header_data, array_data):
        dtype, numel, shape = self.unpack_header(header_data)
        array = np.frombuffer(
            array_data, dtype=np.typeDict[dtype]).reshape(shape)
        return array

        
class MSGPack(object):
    """ Protocol to pack/unpack messages containing header + arrays.
    Use multipart messages for 0mq.
    Arrays are packed/unpacked with arrpack.ArrPack
    """

    def __init__(self, in_header_fmt, out_header_fmt):
        self.inh_packer = struct.Struct(in_header_fmt)
        self.ouh_packer = struct.Struct(out_header_fmt)
        self.arrpack = ArrPack()

    def _pack(self, packer, header, arrays):
        msg_parts = [packer.pack(*header)]
        for arr in arrays:
            msg_parts.extend(self.arrpack.pack_array(arr))
        return msg_parts

    def _unpack(self, packer, msg_parts):
        header = packer.unpack(msg_parts[0])
        arrays = []
        i = 1
        while i < len(msg_parts):
            arrays.append(self.arrpack.unpack_array(
                msg_parts[i], msg_parts[i+1]))
            i += 2
        return header, arrays

    def pack_input(self, header, arrays):
        return self._pack(self.inh_packer, header, arrays)

    def unpack_input(self, msg_parts):
        return self._unpack(self.inh_packer, msg_parts)

    def pack_output(self, header, arrays):
        return self._pack(self.ouh_packer, header, arrays)

    def unpack_output(self, msg_parts):
        return self._unpack(self.ouh_packer, msg_parts)
    
