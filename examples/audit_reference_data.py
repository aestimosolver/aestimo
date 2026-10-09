"""Inventory reference CSV provenance claims without certifying those claims."""

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from aeslibs.validation_policy import SYNTHETIC_FILES, reference_status


def inventory(directory):
    records = []
    for path in sorted(Path(directory).glob('*.csv')):
        content = path.read_bytes()
        comments = [line[1:].strip() for line in content.decode('utf-8').splitlines()
                    if line.startswith('#')]
        kind = ('SYNTHETIC_OR_MODEL_REFERENCE' if path.name in SYNTHETIC_FILES else
                'LITERATURE_BASED_UNVERIFIED' if 'Literature-Based' in '\n'.join(comments) else
                'CLAIMED_MEASUREMENT_OR_DIGITIZATION_UNVERIFIED')
        records.append({
            'filename': path.name,
            'sha256': hashlib.sha256(content).hexdigest(),
            'reference_kind': kind,
            'comparison_status': reference_status([path]),
            'original_source_review': 'NOT PERFORMED',
            'claimed_sources': [line for line in comments if any(
                token in line.lower() for token in ('reference:', 'source:', 'doi:', 'digitization'))],
            'review_required': ('Document the model/generator and its parameters.'
                                if path.name in SYNTHETIC_FILES else
                                'Match these numeric points to the original figure/table/page; '
                                'retain extraction project, method, uncertainty and reviewer evidence.'),
        })
    return records


if __name__ == '__main__':
    destination = ROOT / 'docs' / 'reference-data-audit.json'
    destination.write_text(json.dumps(inventory(ROOT / 'examples' / 'experimental_data'),
                                      indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'Reference inventory saved: {destination}')
