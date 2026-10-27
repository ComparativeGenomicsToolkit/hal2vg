#!/usr/bin/env python3
"""Regression tests for node cuts/orientation invariance and damaged-graph rejection."""
import argparse
import pathlib
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fingerprint', default='./hal2vg-fingerprint')
    parser.add_argument('--fixture', default='./fingerprint-fixture')
    args = parser.parse_args()
    fingerprint = str(pathlib.Path(args.fingerprint).resolve())
    fixture = str(pathlib.Path(args.fixture).resolve())
    compare = str(pathlib.Path(__file__).with_name('compare_graph_fingerprints.py').resolve())
    with tempfile.TemporaryDirectory(prefix='hal2vg-verify-') as directory:
        work = pathlib.Path(directory)

        def generate(mode, valid=True):
            graph = work / f'{mode}.vg'
            manifest = work / f'{mode}.jsonl'
            subprocess.run([fixture, mode, str(graph)], check=True)
            with manifest.open('w') as out:
                result = subprocess.run([fingerprint, str(graph)], stdout=out, stderr=subprocess.PIPE, text=True)
            assert (result.returncode == 0) == valid, (mode, result.stderr)
            return manifest

        def check(old, new, equal=True):
            result = subprocess.run(['python3', compare, str(old), str(new)], capture_output=True, text=True)
            assert (result.returncode == 0) == equal, result.stdout + result.stderr

        baseline = generate('baseline')
        for mode in ['split', 'renumber_flip']:
            check(baseline, generate(mode))
            print(f'PASS invariant: {mode}')
        for mode in ['mutation', 'sharing_loss']:
            check(baseline, generate(mode), False)
            print(f'PASS rejection: {mode}')
        for mode in ['orphan_node', 'unused_edge', 'missing_edge']:
            check(baseline, generate(mode, False), False)
            print(f'PASS rejection: {mode}')
        check(generate('orientation'), generate('orientation_bad'), False)
        print('PASS rejection: same spelling, reversed homology orientation')
        for old, new in [('loop', 'loop_split'), ('circular', 'circular_split')]:
            check(generate(old), generate(new))
            print(f'PASS invariant: {old}')
        hairpin = generate('hairpin')
        check(hairpin, hairpin)
        print('PASS hairpin edge')
        incomplete = work / 'truncated.jsonl'
        incomplete.write_text('\n'.join(baseline.read_text().splitlines()[:-1]) + '\n')
        check(baseline, incomplete, False)
        print('PASS rejection: truncated manifest')
        print('PASS 12 verifier regression cases')


if __name__ == '__main__':
    main()
