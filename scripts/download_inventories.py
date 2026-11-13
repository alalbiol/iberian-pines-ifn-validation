"""Download and verify the exact MITECO IFN databases used by this study.

Uses only Python's standard library. Existing verified databases are skipped;
changed publisher files are rejected rather than silently changing the study.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import time
import unicodedata
from urllib.request import Request, urlopen
import zipfile

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'provenance/inventory_sources.json'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verified(path, entry):
    return path.is_file() and path.stat().st_size == entry['bytes'] and digest(path) == entry['sha256']


def normal_name(name):
    text = Path(name.replace('\\', '/')).name
    return ''.join(c for c in unicodedata.normalize('NFKD', text) if not unicodedata.combining(c)).casefold()


def extract_database(archive, destination, entry):
    """Stream only the selected database; ZIP paths never become output paths."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        matches = [n for n in z.namelist() if normal_name(n) == normal_name(entry['destination'])]
        if len(matches) != 1:
            raise ValueError(f'Expected exactly one {Path(entry["destination"]).name} in {archive}; found {matches}')
        member = z.getinfo(matches[0])
        if member.file_size != entry['bytes']:
            raise ValueError('Publisher database size differs from study input; no file installed')
        with tempfile.NamedTemporaryFile(dir=destination.parent, suffix='.part', delete=False) as out:
            temporary = Path(out.name)
            try:
                with z.open(member) as incoming:
                    shutil.copyfileobj(incoming, out)
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
    try:
        if not verified(temporary, entry):
            raise ValueError('Publisher database SHA-256 differs from study input; no file installed')
        if destination.exists():
            raise FileExistsError(f'Destination appeared during download: {destination}')
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return matches[0]


def download(url, target):
    last_error = None
    for attempt in range(3):
        temporary = target.with_suffix('.part')
        try:
            request = Request(url, headers={'User-Agent': 'iberian-pines-ifn-validation/1.0'})
            with urlopen(request, timeout=120) as response, temporary.open('wb') as out:
                final_url = response.url
                shutil.copyfileobj(response, out)
            if not zipfile.is_zipfile(temporary):
                raise ValueError('MITECO response is not a ZIP archive')
            temporary.replace(target)
            return final_url
        except Exception as error:
            temporary.unlink(missing_ok=True)
            last_error = error
            if attempt < 2:
                time.sleep(attempt + 1)
    raise RuntimeError(f'Download failed for {url}: {last_error}') from last_error


def acquire(entry, data_root, verify_only=False):
    destination = data_root / entry['destination']
    if verified(destination, entry):
        return {'id': entry['id'], 'status': 'verified-existing', 'destination': entry['destination'],
                'sha256': entry['sha256'], 'source_url': entry['url']}
    if destination.exists():
        raise ValueError(f'{destination} differs from the study input. Move it aside before downloading; it will not be overwritten.')
    if verify_only:
        raise FileNotFoundError(destination)
    cache = data_root / 'downloads'
    cache.mkdir(parents=True, exist_ok=True)
    archive = cache / (entry['id'] + '.zip')
    final_url = download(entry['url'], archive) if not archive.exists() else None
    member = extract_database(archive, destination, entry)
    return {'id': entry['id'], 'status': 'downloaded', 'destination': entry['destination'],
            'sha256': entry['sha256'], 'source_url': entry['url'], 'final_url': final_url,
            'archive_sha256': digest(archive), 'archive_member': member}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=Path(os.environ.get('IFN_DATA_ROOT', '~/ifn-data')).expanduser())
    parser.add_argument('--ids', nargs='+', help='Optional manifest IDs; default: the full study coverage')
    parser.add_argument('--list', action='store_true', help='List selected sources without downloading')
    parser.add_argument('--verify-only', action='store_true', help='Verify local files without network access or writes')
    parser.add_argument('--workers', type=int, default=2, choices=range(1, 5), metavar='1-4')
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    entries = manifest['files']
    if args.ids:
        unknown = set(args.ids) - {e['id'] for e in entries}
        if unknown:
            parser.error(f'Unknown IDs: {sorted(unknown)}')
        entries = [e for e in entries if e['id'] in args.ids]
    if args.list:
        for e in entries:
            print(f'{e["id"]}: {e["destination"]} ({e["bytes"]:,} bytes)\n  {e["url"]}')
        print(f'{len(entries)} databases; {sum(e["bytes"] for e in entries):,} bytes uncompressed; ZIP cache additional.')
        return
    results, errors = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending = {pool.submit(acquire, e, args.data_root, args.verify_only): e for e in entries}
        for future in as_completed(pending):
            entry = pending[future]
            try:
                result = future.result()
                results.append(result)
                print(f'{entry["id"]}: {result["status"]}', flush=True)
            except Exception as error:
                errors.append({'id': entry['id'], 'error': str(error)})
                print(f'{entry["id"]}: FAILED: {error}', flush=True)
    if not args.verify_only:
        args.data_root.mkdir(parents=True, exist_ok=True)
        receipt = {'date_utc': datetime.now(timezone.utc).isoformat(),
                   'manifest_sha256': digest(MANIFEST), 'files': results, 'errors': errors}
        (args.data_root / 'inventory_download_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        (args.data_root / 'ifn_validation').mkdir(exist_ok=True)
    print(f'{len(results)}/{len(entries)} verified; {len(errors)} failures.')
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
