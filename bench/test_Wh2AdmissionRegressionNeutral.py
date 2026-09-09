"""Neutral probe tests. No timing calls, scientific launches or DSO mutation."""
import ctypes as T
import os
import struct
import unittest
from unittest.mock import patch

import Wh2AdmissionRegressionNeutral as C


class ProbeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw=[C.A.read_regular(p,4*1024**2) for p,_ in C.LIBRARIES]

    def test_exact_existing_elf_rosters(self):
        for raw,(path,digest) in zip(self.raw,C.LIBRARIES):
            self.assertEqual(C.A.sha(raw),digest)
            elf=C.Elf(raw)
            self.assertEqual(len(elf.exports),53)
            self.assertEqual(len(elf.slots),6)
            for name in elf.exports:
                offset,size=elf.symbol(name,2)
                self.assertTrue(elf.allocated(offset,size,6))
            self.assertEqual(elf.symbol('GF256Ctx',1)[1],0x22810)

    def test_header_and_section_bounds_fail_closed(self):
        for data in (b'',self.raw[0][:64],self.raw[0][:-1],b'X'+self.raw[0][1:]):
            with self.assertRaises(ValueError):C.Elf(data)
        for offset,format,value in ((16,'H',1),(18,'H',183),(40,'Q',2**63),(60,'H',0)):
            data=bytearray(self.raw[0]);struct.pack_into('<'+format,data,offset,value)
            with self.assertRaises(ValueError):C.Elf(bytes(data))

    def test_private_symbol_type_and_bounds(self):
        elf=C.Elf(self.raw[0])
        with self.assertRaises(ValueError):elf.symbol('GF256Ctx',2)
        with self.assertRaises(ValueError):elf.symbol('not_a_symbol',1)
        name='gf256_get_active_x86_cpu_features'
        old=elf.defined[name][0]
        elf.defined[name]=[old[:-2]+(2**63,old[-1])]
        with self.assertRaises(ValueError):elf.symbol(name,2)
        self.assertFalse(elf.allocated(2**63,8,3))

    def test_authentication_precedes_loading(self):
        altered=((C.LIBRARIES[0][0],'0'*64),C.LIBRARIES[1])
        with patch.object(C,'LIBRARIES',altered),patch.object(C.T,'CDLL') as load:
            with self.assertRaisesRegex(ValueError,'exact qualified DSO'):C.Library(0)
            load.assert_not_called()

    def test_bad_environment_and_order_do_not_load(self):
        with patch.object(C,'Library') as load:
            for key in C.ENV_KEYS:
                with patch.dict(os.environ,{key:''},clear=True):
                    with self.assertRaisesRegex(ValueError,'clean loader/allocator'):C.probe('old-new')
            with self.assertRaisesRegex(ValueError,'explicit load order'):C.probe('invalid')
            load.assert_not_called()

    def test_whole_case_roster(self):
        rows=list(C.cases())
        self.assertEqual(len(rows),216)
        self.assertEqual(len(set(rows)),216)
        self.assertEqual([sum(r[0]==f for r in rows) for f in ('certified','wh1','small','k6')],
                         [96,96,12,12])

    def test_failed_encode_frees_its_encoder(self):
        class Fake:
            def __init__(self):self.freed=[]
            def call(self,name,result,*arguments):
                if name.endswith('encoder_create'):
                    def create(source,m,b,policy,profile,capacity):
                        T.memmove(profile,struct.pack('<4sHHQQII',b'WHK3',1,32,0x5748324b33544d31,m,b,0),32)
                        return C.Created(0,4096)
                    return create
                if name.endswith('encode'):return lambda *args:C.Result(2,0,0)
                if name.endswith('free'):return lambda handle:self.freed.append(handle.value)
                raise AssertionError(name)
        fake=Fake()
        with self.assertRaisesRegex(ValueError,'small encode'):
            C.exercise(fake,('small',3,2,2,2))
        self.assertEqual(fake.freed,[4096])

    def test_failed_recover_frees_both_handles(self):
        class Fake:
            def __init__(self):self.freed=[];self.feeds=0
            def call(self,name,result,*arguments):
                if name.endswith('encoder_create'):
                    def create(source,m,b,policy,profile,capacity):
                        T.memmove(profile,struct.pack('<4sHHQQII',b'WHK3',1,32,0x5748324b33544d31,m,b,0),32)
                        return C.Created(0,4096)
                    return create
                if name.endswith('decoder_create'):return lambda *args:C.Created(0,8192)
                if name.endswith('encode'):
                    def encode(handle,packet_id,output,capacity):
                        T.memset(output,0,capacity);return C.Result(0,capacity,capacity)
                    return encode
                if name.endswith('decode'):
                    def feed(*args):
                        self.feeds+=1;return 0 if self.feeds==3 else 1
                    return feed
                if name.endswith('recover'):return lambda *args:C.Result(2,0,0)
                if name.endswith('free'):return lambda handle:self.freed.append(handle.value)
                raise AssertionError(name)
        fake=Fake()
        with self.assertRaisesRegex(ValueError,'small recover'):
            C.exercise(fake,('small',3,2,2,2))
        self.assertEqual(fake.freed,[4096,8192])


if __name__ == '__main__':unittest.main()
