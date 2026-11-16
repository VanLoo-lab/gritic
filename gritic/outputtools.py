"""Shared output collision checks and staged file installation."""

from pathlib import Path

from gritic import validation


def check_output_paths(*, directories=(), files=(), overwrite=False):
    """Reject existing outputs unless the caller explicitly permits reuse."""
    overwrite = validation.validate_boolean(overwrite, 'overwrite')
    for paths, kind in ((directories, 'directory'), (files, 'file')):
        for path in map(Path, paths):
            if not (path.exists() or path.is_symlink()):
                continue
            if kind == 'directory' and not path.is_dir():
                raise FileExistsError(
                    f'Output directory path is not a directory: {path}'
                )
            if kind == 'file' and path.is_dir():
                raise FileExistsError(
                    f'Output file path is a directory: {path}'
                )
            if not overwrite:
                raise FileExistsError(
                    f'Output {kind} already exists: {path}. '
                    'Use --overwrite to reuse existing output locations and '
                    'overwrite files with matching names.'
                )


def replace_files_with_rollback(replacements, backup_directory, *, removals=()):
    """Install files and retire named outputs together, restoring on failure."""
    replacements = tuple(
        (Path(staged_path), Path(target_path))
        for staged_path, target_path in replacements
    )
    removals = tuple(map(Path, removals))
    targets = tuple(target for _, target in replacements) + removals
    if not targets or len(set(targets)) != len(targets):
        raise ValueError('Replacement and removal targets must be nonempty and unique')
    if any(not staged_path.is_file() for staged_path, _ in replacements):
        raise ValueError('Every staged replacement must be an existing file')
    check_output_paths(
        files=targets, overwrite=True,
    )
    backup_directory = Path(backup_directory)
    backup_directory.mkdir()
    backed_up = []
    installed = []
    try:
        for index, target_path in enumerate(targets):
            if target_path.exists() or target_path.is_symlink():
                backup_path = backup_directory / f'{index}-{target_path.name}'
                target_path.replace(backup_path)
                backed_up.append((target_path, backup_path))
        for staged_path, target_path in replacements:
            staged_path.replace(target_path)
            installed.append(target_path)
    except BaseException:
        for target_path in reversed(installed):
            target_path.unlink()
        for target_path, backup_path in reversed(backed_up):
            backup_path.replace(target_path)
        raise
