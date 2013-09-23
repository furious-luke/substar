import os, tempfile, shutil, unittest
from ..Package import *

class TestSearch(unittest.TestCase):

    def test_locate_files_all(self,):
        res = self.ver._locate_files(self.all_hdrs, self.all_dirs)
        self.assertTrue(res[0])
        self.assertSetEqual(res[1], set(self.all_dirs))

    def test_locate_files_some(self,):
        res = self.ver._locate_files(self.all_hdrs[2:], self.all_dirs)
        self.assertTrue(res[0])
        self.assertSetEqual(res[1], set(self.all_dirs[1:]))

    def test_locate_files_fail(self,):
        res = self.ver._locate_files(self.all_hdrs, self.all_dirs[1:])
        self.assertFalse(res[0])
        self.assertEquals(res[1], 'hdr1.h')

    def setUp(self):
        self.pkg = Package()
        self.ver = Version(self.pkg, None)

        # Create the directory structure.
        self.path = tempfile.mkdtemp()
        os.makedirs(os.path.join(self.path, 'a', 'c'))
        os.makedirs(os.path.join(self.path, 'b'))
        file(os.path.join(self.path, 'a', 'c', 'hdr1.h'), 'a');
        file(os.path.join(self.path, 'a', 'c', 'hdr2.h'), 'a');
        file(os.path.join(self.path, 'b', 'hdr3.h'), 'a');
        file(os.path.join(self.path, 'hdr4.h'), 'a');
        self.all_hdrs = ['hdr1.h', 'hdr2.h', 'hdr3.h', 'hdr4.h']
        self.all_dirs = [
            os.path.join(self.path, 'a', 'c'),
            os.path.join(self.path, 'b'),
            os.path.join(self.path),
        ]

    def tearDown(self):
        shutil.rmtree(self.path)

if __name__ == '__main__':
    unittest.main()
