import hashlib
import json
import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path

import numpy as np

from gritic import timingio


class TimingArchiveTest(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)

    def tearDown(self):
        self.temporary_directory.cleanup()

    @staticmethod
    def simple_hierarchy():
        return {'route': np.array([[1, 2], [3, 4]], dtype=np.int16)}

    def write_simple_archive(self, segment_id='segment'):
        return timingio.write_timing_archive(
            self.simple_hierarchy(),
            self.directory,
            segment_id,
        )

    @staticmethod
    def read_manifest(path):
        with Path(path).open('r', encoding='utf-8') as manifest_file:
            return json.load(manifest_file)

    @staticmethod
    def write_manifest(path, manifest):
        with Path(path).open('w', encoding='utf-8') as manifest_file:
            json.dump(manifest, manifest_file)
            manifest_file.write('\n')

    @staticmethod
    def update_archive_hash(manifest, archive_path):
        manifest['archive_sha256'] = hashlib.sha256(
            Path(archive_path).read_bytes()
        ).hexdigest()

    def rewrite_archive(self, archive_path, manifest_path, **arrays):
        with Path(archive_path).open('wb') as archive_file:
            np.savez_compressed(archive_file, **arrays)
        manifest = self.read_manifest(manifest_path)
        self.update_archive_hash(manifest, archive_path)
        self.write_manifest(manifest_path, manifest)

    def assert_manifest_rejected(self, mutator, message):
        archive_path, manifest_path = self.write_simple_archive()
        manifest = self.read_manifest(manifest_path)
        mutator(manifest)
        self.write_manifest(manifest_path, manifest)
        with self.assertRaisesRegex(ValueError, message):
            timingio.load_timing_archive(archive_path, manifest_path)

    def test_canonical_paths_and_manifest_inference(self):
        archive_path, manifest_path = timingio.get_timing_archive_paths(
            self.directory,
            17,
        )
        self.assertEqual(archive_path.name, '17_timing_dict.npz')
        self.assertEqual(
            manifest_path.name,
            '17_timing_dict.manifest.json',
        )
        self.assertEqual(
            timingio.get_timing_manifest_path(archive_path),
            manifest_path,
        )

        with self.assertRaisesRegex(ValueError, 'segment_id must not be empty'):
            timingio.get_timing_archive_paths(self.directory, '')
        with self.assertRaisesRegex(ValueError, 'must end in .npz'):
            timingio.get_timing_manifest_path(self.directory / 'archive.zip')

    def test_round_trip_preserves_nested_key_types_order_shapes_and_numeric_dtypes(self):
        hierarchy = OrderedDict([
            ('boolean', np.array([True, False], dtype=np.bool_)),
            ('signed', np.array([-2, 3], dtype=np.int32)),
            ('unsigned', np.array([2, 3], dtype=np.uint16)),
            ('float', np.array([[0.25, np.nan]], dtype=np.float64)),
            ('complex', np.array(1 + 2j, dtype=np.complex64)),
            (7, OrderedDict([
                ('empty', np.empty((0, 3), dtype=np.float32)),
                (-2, np.array([5], dtype=np.int8)),
            ])),
        ])

        archive_path, manifest_path = timingio.write_timing_archive(
            hierarchy,
            self.directory / 'nested' / 'output',
            'sample:segment',
        )
        restored = timingio.load_timing_archive(archive_path)

        self.assertTrue(archive_path.is_file())
        self.assertTrue(manifest_path.is_file())
        self.assertEqual(list(restored), list(hierarchy))
        self.assertEqual(list(restored[7]), list(hierarchy[7]))
        for key in ('boolean', 'signed', 'unsigned', 'float', 'complex'):
            self.assertEqual(restored[key].dtype, hierarchy[key].dtype)
            self.assertEqual(restored[key].shape, hierarchy[key].shape)
            np.testing.assert_equal(restored[key], hierarchy[key])
        for key in ('empty', -2):
            self.assertEqual(restored[7][key].dtype, hierarchy[7][key].dtype)
            self.assertEqual(restored[7][key].shape, hierarchy[7][key].shape)
            np.testing.assert_equal(restored[7][key], hierarchy[7][key])

        manifest = self.read_manifest(manifest_path)
        self.assertEqual(manifest['format'], timingio.TIMING_ARCHIVE_FORMAT)
        self.assertEqual(manifest['version'], timingio.TIMING_ARCHIVE_VERSION)
        self.assertEqual(manifest['archive'], archive_path.name)
        self.assertEqual(manifest['table_count'], 7)
        self.assertEqual(
            manifest['archive_sha256'],
            hashlib.sha256(archive_path.read_bytes()).hexdigest(),
        )

    def test_empty_mapping_round_trips_as_a_valid_zero_table_archive(self):
        archive_path, manifest_path = timingio.write_timing_archive(
            {},
            self.directory,
            'empty',
        )
        self.assertEqual(timingio.load_timing_archive(archive_path), {})
        manifest = self.read_manifest(manifest_path)
        self.assertEqual(manifest['table_count'], 0)
        self.assertEqual(manifest['tables'], [])
        with np.load(archive_path, allow_pickle=False) as archive:
            self.assertEqual(archive.files, [])

    def test_overwrite_replaces_both_members_and_leaves_no_temporary_files(self):
        archive_path, manifest_path = timingio.write_timing_archive(
            {'old': np.array([1], dtype=np.int8)},
            self.directory,
            'replace',
        )
        returned_paths = timingio.write_timing_archive(
            {'new': np.array([1000], dtype=np.int64)},
            self.directory,
            'replace',
        )

        self.assertEqual(returned_paths, (archive_path, manifest_path))
        restored = timingio.load_timing_archive(archive_path, manifest_path)
        self.assertEqual(list(restored), ['new'])
        np.testing.assert_array_equal(restored['new'], [1000])
        self.assertEqual(
            [path for path in self.directory.iterdir() if path.name.startswith('.')],
            [],
        )

    def test_writer_rejects_invalid_roots_keys_and_leaves(self):
        invalid_cases = (
            (np.array([1]), 'hierarchy must be a mapping'),
            ({True: np.array([1])}, 'non-boolean integers'),
            ({1.5: np.array([1])}, 'observed float'),
            ({'leaf': [1, 2]}, 'leaves must be numpy arrays'),
            ({'leaf': np.array(['text'])}, 'plain numeric arrays'),
            ({'leaf': np.array([object()], dtype=object)}, 'plain numeric arrays'),
            ({'leaf': np.array([(1,)], dtype=[('value', 'i4')])}, 'plain numeric arrays'),
        )
        for hierarchy, message in invalid_cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex((TypeError, ValueError), message):
                    timingio.write_timing_archive(
                        hierarchy,
                        self.directory,
                        'invalid',
                    )

    def test_missing_archive_or_manifest_propagates_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            timingio.load_timing_archive(self.directory / 'missing.npz')

        archive_path, manifest_path = self.write_simple_archive()
        manifest_path.unlink()
        with self.assertRaises(FileNotFoundError):
            timingio.load_timing_archive(archive_path)

    def test_manifest_rejects_nonobject_missing_and_unexpected_root_fields(self):
        archive_path, manifest_path = self.write_simple_archive()
        self.write_manifest(manifest_path, [])
        with self.assertRaisesRegex(ValueError, 'manifest must be a JSON object'):
            timingio.load_timing_archive(archive_path, manifest_path)

        self.assert_manifest_rejected(
            lambda manifest: manifest.pop('format'),
            'invalid fields.*missing: format',
        )
        self.assert_manifest_rejected(
            lambda manifest: manifest.__setitem__('extra', 1),
            'invalid fields.*unexpected: extra',
        )

    def test_manifest_rejects_duplicate_json_keys_at_any_level(self):
        archive_path, manifest_path = self.write_simple_archive()
        manifest_text = manifest_path.read_text(encoding='utf-8')
        manifest_path.write_text(
            manifest_text.replace(
                '"format": "gritic_timing_archive",',
                '"format": "first",\n  "format": "gritic_timing_archive",',
                1,
            ),
            encoding='utf-8',
        )
        with self.assertRaisesRegex(ValueError, "duplicate JSON key 'format'"):
            timingio.load_timing_archive(archive_path, manifest_path)

    def test_manifest_rejects_format_version_and_archive_name_mismatches(self):
        cases = (
            (
                lambda manifest: manifest.__setitem__('format', 'other'),
                'Unsupported GRITIC timing archive format',
            ),
            (
                lambda manifest: manifest.__setitem__(
                    'version',
                    timingio.TIMING_ARCHIVE_VERSION + 1,
                ),
                'Unsupported GRITIC timing archive version',
            ),
            (
                lambda manifest: manifest.__setitem__('archive', 'other.npz'),
                'archive name does not match',
            ),
        )
        for mutator, message in cases:
            with self.subTest(message=message):
                self.assert_manifest_rejected(mutator, message)

    def test_manifest_rejects_malformed_hashes_and_archive_tampering(self):
        invalid_hashes = (None, '', '0' * 63, 'A' * 64, 'g' * 64)
        for invalid_hash in invalid_hashes:
            with self.subTest(invalid_hash=invalid_hash):
                self.assert_manifest_rejected(
                    lambda manifest, value=invalid_hash: manifest.__setitem__(
                        'archive_sha256',
                        value,
                    ),
                    'archive_sha256 must be a lowercase SHA-256',
                )

        archive_path, manifest_path = self.write_simple_archive()
        archive_path.write_bytes(archive_path.read_bytes() + b'tampered')
        with self.assertRaisesRegex(ValueError, 'SHA-256 does not match'):
            timingio.load_timing_archive(archive_path, manifest_path)

    def test_table_count_and_tables_container_are_strict(self):
        cases = (
            (
                lambda manifest: manifest.__setitem__('table_count', True),
                'table_count must be a non-negative integer',
            ),
            (
                lambda manifest: manifest.__setitem__('table_count', -1),
                'table_count must be a non-negative integer',
            ),
            (
                lambda manifest: manifest.__setitem__('tables', {}),
                'tables must be a JSON array',
            ),
            (
                lambda manifest: manifest.__setitem__('table_count', 2),
                'table_count does not match tables length',
            ),
        )
        for mutator, message in cases:
            with self.subTest(message=message):
                self.assert_manifest_rejected(mutator, message)

    def test_table_metadata_requires_exact_fields_and_linear_index_and_key(self):
        cases = (
            (
                lambda manifest: manifest['tables'].__setitem__(0, []),
                'table metadata 0 must be a JSON object',
            ),
            (
                lambda manifest: manifest['tables'][0].pop('dtype'),
                'invalid fields.*missing: dtype',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__('extra', 1),
                'invalid fields.*unexpected: extra',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__('index', True),
                'index must be a non-negative integer',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__('index', 1),
                'indexes must be consecutive',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__(
                    'archive_key',
                    'table_000001',
                ),
                'archive keys must match linear table order',
            ),
        )
        for mutator, message in cases:
            with self.subTest(message=message):
                self.assert_manifest_rejected(mutator, message)

    def test_table_metadata_rejects_invalid_unsupported_dtype_and_shape(self):
        cases = (
            (
                lambda manifest: manifest['tables'][0].__setitem__('dtype', 'not-a-dtype'),
                'invalid dtype',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__('dtype', '|O'),
                'unsupported dtype',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__(
                    'dtype',
                    np.dtype([('value', 'i4')]).str,
                ),
                'unsupported dtype',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__('shape', {}),
                'shape must be a JSON array',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__('shape', [True]),
                'shape dimension must be a non-negative integer',
            ),
            (
                lambda manifest: manifest['tables'][0].__setitem__('shape', [-1]),
                'shape dimension must be a non-negative integer',
            ),
        )
        for mutator, message in cases:
            with self.subTest(message=message):
                self.assert_manifest_rejected(mutator, message)

    def test_npz_entries_must_exactly_match_manifest_order(self):
        archive_path, manifest_path = self.write_simple_archive()
        original = self.simple_hierarchy()['route']
        self.rewrite_archive(
            archive_path,
            manifest_path,
            unexpected=original,
        )
        with self.assertRaisesRegex(ValueError, 'entries do not match'):
            timingio.load_timing_archive(archive_path, manifest_path)

        hierarchy = {
            'first': np.array([1]),
            'second': np.array([2]),
        }
        archive_path, manifest_path = timingio.write_timing_archive(
            hierarchy,
            self.directory,
            'ordered',
        )
        self.rewrite_archive(
            archive_path,
            manifest_path,
            table_000001=hierarchy['second'],
            table_000000=hierarchy['first'],
        )
        with self.assertRaisesRegex(ValueError, 'entries do not match'):
            timingio.load_timing_archive(archive_path, manifest_path)

    def test_npz_table_must_match_manifest_dtype_and_shape(self):
        archive_path, manifest_path = self.write_simple_archive('dtype')
        self.rewrite_archive(
            archive_path,
            manifest_path,
            table_000000=np.array([[1, 2], [3, 4]], dtype=np.float32),
        )
        with self.assertRaisesRegex(ValueError, 'dtype does not match'):
            timingio.load_timing_archive(archive_path, manifest_path)

        archive_path, manifest_path = self.write_simple_archive('shape')
        self.rewrite_archive(
            archive_path,
            manifest_path,
            table_000000=np.array([1, 2, 3, 4], dtype=np.int16),
        )
        with self.assertRaisesRegex(ValueError, 'shape does not match'):
            timingio.load_timing_archive(archive_path, manifest_path)

    def test_object_npz_is_never_loaded_with_pickle(self):
        archive_path, manifest_path = self.write_simple_archive()
        self.rewrite_archive(
            archive_path,
            manifest_path,
            table_000000=np.array([[object(), object()], [object(), object()]], dtype=object),
        )
        with self.assertRaisesRegex(ValueError, 'cannot be loaded without pickle'):
            timingio.load_timing_archive(archive_path, manifest_path)

    def test_hierarchy_nodes_require_supported_type_and_exact_fields(self):
        cases = (
            (
                lambda manifest: manifest.__setitem__('hierarchy', []),
                'hierarchy must be an object',
            ),
            (
                lambda manifest: manifest.__setitem__('hierarchy', {'type': 'other'}),
                'unsupported node type',
            ),
            (
                lambda manifest: manifest['hierarchy'].pop('entries'),
                'invalid fields.*missing: entries',
            ),
            (
                lambda manifest: manifest['hierarchy'].__setitem__('extra', 1),
                'invalid fields.*unexpected: extra',
            ),
            (
                lambda manifest: manifest['hierarchy'].__setitem__('entries', {}),
                'entries must be an array',
            ),
        )
        for mutator, message in cases:
            with self.subTest(message=message):
                self.assert_manifest_rejected(mutator, message)

    def test_hierarchy_entries_and_encoded_keys_are_strict(self):
        def first_entry(manifest):
            return manifest['hierarchy']['entries'][0]

        cases = (
            (
                lambda manifest: manifest['hierarchy']['entries'].__setitem__(0, []),
                r'entries\[0\] must be a JSON object',
            ),
            (
                lambda manifest: first_entry(manifest).pop('value'),
                'invalid fields.*missing: value',
            ),
            (
                lambda manifest: first_entry(manifest)['key'].pop('value'),
                'invalid fields.*missing: value',
            ),
            (
                lambda manifest: first_entry(manifest)['key'].__setitem__('type', 'other'),
                'unsupported key type',
            ),
            (
                lambda manifest: first_entry(manifest)['key'].__setitem__('value', 1),
                'string value must be text',
            ),
            (
                lambda manifest: first_entry(manifest)['key'].update({
                    'type': 'integer',
                    'value': True,
                }),
                'integer value must be an integer',
            ),
        )
        for mutator, message in cases:
            with self.subTest(message=message):
                self.assert_manifest_rejected(mutator, message)

    def test_hierarchy_rejects_duplicate_decoded_mapping_keys(self):
        def duplicate_entry(manifest):
            entry = manifest['hierarchy']['entries'][0]
            manifest['hierarchy']['entries'].append(json.loads(json.dumps(entry)))

        self.assert_manifest_rejected(
            duplicate_entry,
            "has duplicate key 'route'",
        )

    def test_hierarchy_table_references_are_bounded_unique_and_exhaustive(self):
        cases = (
            (
                lambda manifest: manifest['hierarchy']['entries'][0]['value'].__setitem__(
                    'table_index',
                    True,
                ),
                'table_index must be a non-negative integer',
            ),
            (
                lambda manifest: manifest['hierarchy']['entries'][0]['value'].__setitem__(
                    'table_index',
                    2,
                ),
                'references missing table 2',
            ),
            (
                lambda manifest: manifest['hierarchy']['entries'][0]['value'].__setitem__(
                    'extra',
                    1,
                ),
                'invalid fields.*unexpected: extra',
            ),
        )
        for mutator, message in cases:
            with self.subTest(message=message):
                self.assert_manifest_rejected(mutator, message)

        archive_path, manifest_path = timingio.write_timing_archive(
            {
                'first': np.array([1]),
                'second': np.array([2]),
            },
            self.directory,
            'references',
        )
        manifest = self.read_manifest(manifest_path)
        entries = manifest['hierarchy']['entries']
        entries[1]['value']['table_index'] = 0
        self.write_manifest(manifest_path, manifest)
        with self.assertRaisesRegex(ValueError, 'references table 0 more than once'):
            timingio.load_timing_archive(archive_path, manifest_path)

        archive_path, manifest_path = timingio.write_timing_archive(
            {
                'first': np.array([1]),
                'second': np.array([2]),
            },
            self.directory,
            'unused',
        )
        manifest = self.read_manifest(manifest_path)
        manifest['hierarchy']['entries'].pop()
        self.write_manifest(manifest_path, manifest)
        with self.assertRaisesRegex(ValueError, r'does not reference every table.*unused: \[1\]'):
            timingio.load_timing_archive(archive_path, manifest_path)

    def test_manifest_root_must_restore_to_a_mapping(self):
        archive_path, manifest_path = self.write_simple_archive()
        manifest = self.read_manifest(manifest_path)
        manifest['hierarchy'] = {'type': 'table', 'table_index': 0}
        self.write_manifest(manifest_path, manifest)
        with self.assertRaisesRegex(ValueError, 'root must be a mapping'):
            timingio.load_timing_archive(archive_path, manifest_path)


if __name__ == '__main__':
    unittest.main()
