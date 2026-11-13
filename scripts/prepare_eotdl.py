"""Build a local EOTDL upload folder with a checksum manifest; never upload."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
FILES = ['README.md', 'LICENSE', 'LICENSE-DATA', 'DATA.md', 'REPRODUCING.md',
         'EOTDL.md', 'requirements.txt', 'requirements-lock.txt']
DIRECTORIES = ['analysis', 'models', 'scripts', 'provenance', 'expected', 'tests']
SUFFIXES = {'.py', '.md', '.csv', '.json', '.geojson', '.png', '.pdf', '.txt'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'dist/iberian-pines-ifn-validation')
    args = parser.parse_args()
    readme = (ROOT / 'README.md').read_text()
    if not readme.startswith('---\n'):
        raise ValueError('README must start with EOTDL YAML metadata')
    metadata = readme.split('---', 2)[1]
    match = re.search(r'^name: ([a-zA-Z0-9-]{3,45})$', metadata, re.M)
    if not match or 'license: CC-BY-4.0' not in metadata:
        raise ValueError('Invalid dataset name or unexpected dataset licence')
    sources = [ROOT / name for name in FILES]
    for directory in DIRECTORIES:
        sources.extend(p for p in (ROOT / directory).rglob('*')
                       if p.is_file() and p.suffix in SUFFIXES
                       and not any(part.startswith('.') or part == '__pycache__'
                                   for part in p.relative_to(ROOT).parts))
    for source in sources:
        if source.is_symlink() or not source.is_file():
            raise ValueError(f'Expected a regular release file: {source}')
    # Refuse to overwrite a previous upload snapshot.
    args.output.mkdir(parents=True, exist_ok=False)
    entries = []
    for source in sorted(sources):
        relative = source.relative_to(ROOT)
        target = args.output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        entries.append({'path': relative.as_posix(), 'bytes': target.stat().st_size,
                        'sha256': hashlib.sha256(target.read_bytes()).hexdigest()})
    manifest = {'name': match[1], 'files': entries,
                'note': 'This manifest lists every payload file except itself.'}
    (args.output / 'release-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Prepared {len(entries)} files ({sum(e["bytes"] for e in entries):,} bytes) in {args.output}')
    print('Nothing uploaded. See EOTDL.md for authentication and ingestion.')


if __name__ == '__main__':
    main()
