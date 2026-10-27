#!/usr/bin/env python3
"""Compare complete graph manifests; node IDs and linear subdivisions may differ."""
import argparse
import json
import sys


def read(path):
    paths, summary = {}, None
    with open(path, encoding='utf-8') as stream:
        for line in stream:
            row = json.loads(line)
            if summary is not None:
                raise ValueError('records after the completion summary')
            if row['type'] == 'path':
                if row['name'] in paths:
                    raise ValueError('duplicate path record')
                paths[row['name']] = row
            elif row['type'] == 'summary':
                summary = row
            else:
                raise ValueError('unknown manifest record')
    if not summary or not summary.get('complete') or summary.get('format_version') != '1':
        raise ValueError('manifest is incomplete or has an unsupported version')
    if summary.get('unused_nodes') != '0' or summary.get('unused_edges') != '0':
        raise ValueError('manifest did not establish full path coverage')
    if int(summary['path_count']) != len(paths) or int(summary['total_bases']) != sum(int(x['length']) for x in paths.values()):
        raise ValueError('manifest summary does not match its path records')
    return paths, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('old')
    parser.add_argument('new')
    args = parser.parse_args()
    old, old_summary = read(args.old)
    new, new_summary = read(args.new)
    if old.keys() != new.keys():
        raise ValueError(f"path sets differ: missing={list(old.keys() - new.keys())[:5]}, extra={list(new.keys() - old.keys())[:5]}")
    fields = ['metadata', 'length', 'sequence_sha256', 'uppercase_sha256', 'homology_sha256', 'homology_intervals']
    failures = []
    for name in sorted(old):
        different = [field for field in fields if old[name][field] != new[name][field]]
        if different:
            failures.append((name, different))
    if failures:
        for name, fields in failures[:20]:
            print(f'FAIL {name}: {", ".join(fields)}', file=sys.stderr)
        raise ValueError(f'{len(failures)} paths differ')
    print(f"PASS: {len(old)} paths; {old_summary['total_bases']} path bases")
    print('Identical metadata, case-sensitive sequences and oriented base-level homology maps.')
    print('All nodes/edges are used by paths, and every path transition has an edge.')
    print(f"Node counts (diagnostic only): {old_summary['node_count']} / {new_summary['node_count']}")


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        print(f'FAIL: {error}', file=sys.stderr)
        sys.exit(1)
