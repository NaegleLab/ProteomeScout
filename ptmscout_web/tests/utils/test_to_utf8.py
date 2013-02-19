import unittest
import os
import codecs
from ptmscout.utils import to_utf8
from ptmscout.config import settings

class UTF8TestCase(unittest.TestCase):
    def test_convert_encoding(self):
        source = os.path.join(settings.ptmscout_path, settings.experiment_data_file_path, 'test', 'Rikova4.txt')
        dest = source + ".tmp"

        to_utf8.convert_encoding_to_utf8(source, dest)

        with codecs.open(dest, 'rb', encoding='utf8') as f:
            for line in f:
                pass

        os.remove(dest)

if __name__ == '__main__':
    unittest.main()
