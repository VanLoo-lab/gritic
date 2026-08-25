"""Safe, linear storage for GRITIC timing-array hierarchies."""

import hashlib
import hmac
import json
import os
import tempfile
from collections.abc import Mapping
from numbers import Integral
from pathlib import Path

import numpy as np


TIMING_ARCHIVE_FORMAT = 'gritic_timing_archive'
TIMING_ARCHIVE_VERSION = 1
TIMING_ARCHIVE_SUFFIX = '_timing_dict.npz'
TIMING_MANIFEST_SUFFIX = '_timing_dict.manifest.json'

_MANIFEST_KEYS = {
    'format',
    'version',
    'archive',
    'archive_sha256',
    'table_count',
    'hierarchy',
    'tables',
}
_TABLE_KEYS = {'index', 'archive_key', 'dtype', 'shape'}
_MAPPING_NODE_KEYS = {'type', 'entries'}
_TABLE_NODE_KEYS = {'type', 'table_index'}
_ENTRY_KEYS = {'key', 'value'}
_KEY_KEYS = {'type', 'value'}


def get_timing_archive_paths(directory, segment_id):
    """Return the canonical NPZ and JSON paths for one logical timing object."""
    segment_id = str(segment_id)
    if not segment_id:
        raise ValueError('segment_id must not be empty')
    directory = Path(directory)
    return (
        directory / f'{segment_id}{TIMING_ARCHIVE_SUFFIX}',
        directory / f'{segment_id}{TIMING_MANIFEST_SUFFIX}',
    )


def get_timing_manifest_path(archive_path):
    """Return the same-stem manifest path paired with an NPZ archive path."""
    archive_path = Path(archive_path)
    if archive_path.suffix != '.npz':
        raise ValueError('GRITIC timing archive path must end in .npz')
    return archive_path.with_suffix('.manifest.json')


def _encode_mapping_key(key):
    if isinstance(key, str):
        return {'type': 'string', 'value': key}
    if isinstance(key, Integral) and not isinstance(key, (bool, np.bool_)):
        return {'type': 'integer', 'value': int(key)}
    raise TypeError(
        'GRITIC timing hierarchy keys must be strings or non-boolean integers; '
        f'observed {type(key).__name__}'
    )


def _linearize_hierarchy(value, tables):
    if isinstance(value, Mapping):
        entries = []
        for key, child in value.items():
            entries.append({
                'key': _encode_mapping_key(key),
                'value': _linearize_hierarchy(child, tables),
            })
        return {'type': 'mapping', 'entries': entries}

    if not isinstance(value, np.ndarray):
        raise TypeError(
            'GRITIC timing hierarchy leaves must be numpy arrays; observed '
            f'{type(value).__name__}'
        )
    if value.dtype.kind not in 'biufc' or value.dtype.fields is not None:
        raise TypeError(
            'GRITIC timing archives support only plain numeric arrays'
        )

    table_index = len(tables)
    tables.append(value)
    return {'type': 'table', 'table_index': table_index}


def _temporary_path(target_path):
    target_path = Path(target_path)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=target_path.parent,
        prefix=f'.{target_path.name}.',
        suffix='.tmp',
    )
    os.close(descriptor)
    return Path(temporary_name)


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def write_timing_archive(timing_hierarchy, directory, segment_id):
    """Write one compressed linear NPZ and its hierarchy manifest."""
    if not isinstance(timing_hierarchy, Mapping):
        raise TypeError('GRITIC timing hierarchy must be a mapping')

    archive_path, manifest_path = get_timing_archive_paths(
        directory,
        segment_id,
    )
    archive_path.parent.mkdir(parents=True, exist_ok=True)

    tables = []
    hierarchy = _linearize_hierarchy(timing_hierarchy, tables)
    table_metadata = []
    archive_tables = {}
    for table_index, table in enumerate(tables):
        archive_key = f'table_{table_index:06d}'
        archive_tables[archive_key] = table
        table_metadata.append({
            'index': table_index,
            'archive_key': archive_key,
            'dtype': table.dtype.str,
            'shape': [int(dimension) for dimension in table.shape],
        })

    temporary_archive = _temporary_path(archive_path)
    temporary_manifest = _temporary_path(manifest_path)
    try:
        with temporary_archive.open('wb') as archive_file:
            np.savez_compressed(archive_file, **archive_tables)

        manifest = {
            'format': TIMING_ARCHIVE_FORMAT,
            'version': TIMING_ARCHIVE_VERSION,
            'archive': archive_path.name,
            'archive_sha256': _file_sha256(temporary_archive),
            'table_count': len(tables),
            'hierarchy': hierarchy,
            'tables': table_metadata,
        }
        with temporary_manifest.open('w', encoding='utf-8') as manifest_file:
            json.dump(
                manifest,
                manifest_file,
                ensure_ascii=False,
                allow_nan=False,
                indent=2,
            )
            manifest_file.write('\n')

        os.replace(temporary_archive, archive_path)
        os.replace(temporary_manifest, manifest_path)
    finally:
        temporary_archive.unlink(missing_ok=True)
        temporary_manifest.unlink(missing_ok=True)

    return archive_path, manifest_path


def _require_exact_keys(value, expected_keys, description):
    if not isinstance(value, dict):
        raise ValueError(f'{description} must be a JSON object')
    observed_keys = set(value)
    if observed_keys != expected_keys:
        missing = sorted(expected_keys - observed_keys)
        unexpected = sorted(observed_keys - expected_keys)
        details = []
        if missing:
            details.append('missing: ' + ', '.join(missing))
        if unexpected:
            details.append('unexpected: ' + ', '.join(unexpected))
        raise ValueError(
            f'{description} has invalid fields ({"; ".join(details)})'
        )


def _require_nonnegative_integer(value, description):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f'{description} must be a non-negative integer')
    return value


def _reject_duplicate_json_keys(pairs):
    value = {}
    for key, child in pairs:
        if key in value:
            raise ValueError(
                f'GRITIC timing manifest contains duplicate JSON key {key!r}'
            )
        value[key] = child
    return value


def _validate_table_metadata(manifest):
    table_count = _require_nonnegative_integer(
        manifest['table_count'],
        'GRITIC timing manifest table_count',
    )
    table_metadata = manifest['tables']
    if not isinstance(table_metadata, list):
        raise ValueError('GRITIC timing manifest tables must be a JSON array')
    if len(table_metadata) != table_count:
        raise ValueError(
            'GRITIC timing manifest table_count does not match tables length'
        )

    expected_archive_keys = []
    for expected_index, metadata in enumerate(table_metadata):
        _require_exact_keys(
            metadata,
            _TABLE_KEYS,
            f'GRITIC timing table metadata {expected_index}',
        )
        table_index = _require_nonnegative_integer(
            metadata['index'],
            f'GRITIC timing table metadata {expected_index} index',
        )
        if table_index != expected_index:
            raise ValueError(
                'GRITIC timing manifest table indexes must be consecutive and '
                'match their linear order'
            )

        expected_archive_key = f'table_{expected_index:06d}'
        if metadata['archive_key'] != expected_archive_key:
            raise ValueError(
                'GRITIC timing manifest archive keys must match linear table '
                f'order; expected {expected_archive_key!r}'
            )
        expected_archive_keys.append(expected_archive_key)

        try:
            dtype = np.dtype(metadata['dtype'])
        except (TypeError, ValueError) as error:
            raise ValueError(
                f'GRITIC timing table {expected_index} has an invalid dtype'
            ) from error
        if dtype.kind not in 'biufc' or dtype.fields is not None:
            raise ValueError(
                f'GRITIC timing table {expected_index} has an unsupported dtype'
            )

        shape = metadata['shape']
        if not isinstance(shape, list):
            raise ValueError(
                f'GRITIC timing table {expected_index} shape must be a JSON array'
            )
        for dimension in shape:
            _require_nonnegative_integer(
                dimension,
                f'GRITIC timing table {expected_index} shape dimension',
            )

    return table_metadata, expected_archive_keys


def _decode_mapping_key(encoded_key, description):
    _require_exact_keys(encoded_key, _KEY_KEYS, description)
    key_type = encoded_key['type']
    key_value = encoded_key['value']
    if key_type == 'string':
        if not isinstance(key_value, str):
            raise ValueError(f'{description} string value must be text')
        return key_value
    if key_type == 'integer':
        if isinstance(key_value, bool) or not isinstance(key_value, int):
            raise ValueError(f'{description} integer value must be an integer')
        return key_value
    raise ValueError(f'{description} has unsupported key type {key_type!r}')


def _restore_hierarchy(node, tables, used_indexes, description='hierarchy'):
    if not isinstance(node, dict):
        raise ValueError(f'GRITIC timing manifest {description} must be an object')
    node_type = node.get('type')
    if node_type == 'mapping':
        _require_exact_keys(
            node,
            _MAPPING_NODE_KEYS,
            f'GRITIC timing manifest {description}',
        )
        entries = node['entries']
        if not isinstance(entries, list):
            raise ValueError(
                f'GRITIC timing manifest {description} entries must be an array'
            )
        restored = {}
        for entry_index, entry in enumerate(entries):
            entry_description = f'{description}.entries[{entry_index}]'
            _require_exact_keys(
                entry,
                _ENTRY_KEYS,
                f'GRITIC timing manifest {entry_description}',
            )
            key = _decode_mapping_key(
                entry['key'],
                f'GRITIC timing manifest {entry_description}.key',
            )
            if key in restored:
                raise ValueError(
                    f'GRITIC timing manifest {description} has duplicate key '
                    f'{key!r}'
                )
            restored[key] = _restore_hierarchy(
                entry['value'],
                tables,
                used_indexes,
                f'{entry_description}.value',
            )
        return restored

    if node_type == 'table':
        _require_exact_keys(
            node,
            _TABLE_NODE_KEYS,
            f'GRITIC timing manifest {description}',
        )
        table_index = _require_nonnegative_integer(
            node['table_index'],
            f'GRITIC timing manifest {description} table_index',
        )
        if table_index >= len(tables):
            raise ValueError(
                f'GRITIC timing manifest {description} references missing '
                f'table {table_index}'
            )
        if table_index in used_indexes:
            raise ValueError(
                f'GRITIC timing manifest references table {table_index} more '
                'than once'
            )
        used_indexes.add(table_index)
        return tables[table_index]

    raise ValueError(
        f'GRITIC timing manifest {description} has unsupported node type '
        f'{node_type!r}'
    )


def load_timing_archive(archive_path, manifest_path=None):
    """Load and validate a GRITIC NPZ/JSON pair without enabling pickle."""
    archive_path = Path(archive_path)
    if manifest_path is None:
        manifest_path = get_timing_manifest_path(archive_path)
    else:
        manifest_path = Path(manifest_path)

    with manifest_path.open('r', encoding='utf-8') as manifest_file:
        manifest = json.load(
            manifest_file,
            object_pairs_hook=_reject_duplicate_json_keys,
        )
    _require_exact_keys(
        manifest,
        _MANIFEST_KEYS,
        'GRITIC timing manifest',
    )
    if manifest['format'] != TIMING_ARCHIVE_FORMAT:
        raise ValueError(
            f'Unsupported GRITIC timing archive format {manifest["format"]!r}'
        )
    if manifest['version'] != TIMING_ARCHIVE_VERSION:
        raise ValueError(
            f'Unsupported GRITIC timing archive version {manifest["version"]!r}'
        )
    if manifest['archive'] != archive_path.name:
        raise ValueError(
            'GRITIC timing manifest archive name does not match its NPZ pair'
        )
    archive_sha256 = manifest['archive_sha256']
    if (
        not isinstance(archive_sha256, str)
        or len(archive_sha256) != 64
        or any(character not in '0123456789abcdef' for character in archive_sha256)
    ):
        raise ValueError(
            'GRITIC timing manifest archive_sha256 must be a lowercase SHA-256'
        )
    if not hmac.compare_digest(_file_sha256(archive_path), archive_sha256):
        raise ValueError(
            'GRITIC timing NPZ SHA-256 does not match its JSON manifest'
        )

    table_metadata, expected_archive_keys = _validate_table_metadata(manifest)
    tables = []
    with np.load(archive_path, allow_pickle=False) as archive:
        if archive.files != expected_archive_keys:
            raise ValueError(
                'GRITIC timing NPZ entries do not match manifest linear order'
            )
        for table_index, metadata in enumerate(table_metadata):
            archive_key = metadata['archive_key']
            try:
                table = archive[archive_key]
            except ValueError as error:
                raise ValueError(
                    f'GRITIC timing table {table_index} cannot be loaded '
                    'without pickle'
                ) from error
            if table.dtype.str != metadata['dtype']:
                raise ValueError(
                    f'GRITIC timing table {table_index} dtype does not match '
                    'the manifest'
                )
            if list(table.shape) != metadata['shape']:
                raise ValueError(
                    f'GRITIC timing table {table_index} shape does not match '
                    'the manifest'
                )
            if table.dtype.kind not in 'biufc' or table.dtype.fields is not None:
                raise ValueError(
                    f'GRITIC timing table {table_index} has an unsupported dtype'
                )
            tables.append(table)

    used_indexes = set()
    timing_hierarchy = _restore_hierarchy(
        manifest['hierarchy'],
        tables,
        used_indexes,
    )
    if used_indexes != set(range(len(tables))):
        unused_indexes = sorted(set(range(len(tables))) - used_indexes)
        raise ValueError(
            'GRITIC timing manifest hierarchy does not reference every table; '
            f'unused: {unused_indexes}'
        )
    if not isinstance(timing_hierarchy, dict):
        raise ValueError('GRITIC timing manifest root must be a mapping')
    return timing_hierarchy
