"""Download/extraction checks with tiny local archives, without MITECO traffic."""
import hashlib
import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/download_inventories.py'
spec = importlib.util.spec_from_file_location('download_inventories', SCRIPT)
downloader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(downloader)


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.payload = b'study database contents'
        self.archive = self.root / 'source.zip'
        with zipfile.ZipFile(self.archive, 'w') as z:
            z.writestr('nested/Ifn4_Cataluña.accdb', self.payload)
            z.writestr('../../unwanted.txt', b'never extract this')
        self.entry = {'id': 'fixture', 'destination': 'ifn4/Ifn4_Cataluna.accdb',
                      'bytes': len(self.payload), 'sha256': hashlib.sha256(self.payload).hexdigest(),
                      'url': self.archive.as_uri()}
        self.data_root = self.root / 'data'

    def test_download_renames_unicode_member_and_skips_verified_file(self):
        result = downloader.acquire(self.entry, self.data_root)
        self.assertEqual(result['status'], 'downloaded')
        self.assertEqual((self.data_root / self.entry['destination']).read_bytes(), self.payload)
        self.assertFalse((self.root / 'unwanted.txt').exists())
        self.archive.unlink()  # a second run must work without reaching the source
        self.assertEqual(downloader.acquire(self.entry, self.data_root)['status'], 'verified-existing')

    def test_hash_mismatch_does_not_install_database(self):
        entry = dict(self.entry, sha256='0' * 64)
        with self.assertRaisesRegex(ValueError, 'SHA-256'):
            downloader.acquire(entry, self.data_root)
        self.assertFalse((self.data_root / entry['destination']).exists())
        self.assertFalse(list(self.data_root.rglob('*.part')))

    def test_existing_different_file_is_preserved(self):
        target = self.data_root / self.entry['destination']
        target.parent.mkdir(parents=True)
        target.write_bytes(b'local changes')
        with self.assertRaisesRegex(ValueError, 'will not be overwritten'):
            downloader.acquire(self.entry, self.data_root)
        self.assertEqual(target.read_bytes(), b'local changes')

    def test_verify_only_does_not_download_or_create_directories(self):
        with self.assertRaises(FileNotFoundError):
            downloader.acquire(self.entry, self.data_root, verify_only=True)
        self.assertFalse(self.data_root.exists())

    def test_ambiguous_archive_is_rejected(self):
        with zipfile.ZipFile(self.archive, 'a') as z:
            z.writestr('other/Ifn4_Cataluna.accdb', self.payload)
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            downloader.acquire(self.entry, self.data_root)
        self.assertFalse((self.data_root / self.entry['destination']).exists())


if __name__ == '__main__':
    unittest.main()
