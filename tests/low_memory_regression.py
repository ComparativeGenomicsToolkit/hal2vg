#!/usr/bin/env python3
"""Generate small HALs and compare base-level homology, edges and paths, independent of node IDs."""
import argparse
import collections
import json
import pathlib
import random
import shutil
import subprocess

COMPLEMENT = str.maketrans('ACGTNRYKMSWBDHV', 'TGCANYRMKSWVHDB')


def complement(base):
    return base.upper().translate(COMPLEMENT)


def read_graph(path):
    text = path.read_text()
    decoder = json.JSONDecoder()
    nodes, edges, paths = {}, [], {}
    while text.strip():
        text = text.lstrip()
        chunk, end = decoder.raw_decode(text)
        text = text[end:]
        for node in chunk.get('node', []):
            nodes[int(node['id'])] = node['sequence'].upper()
        edges.extend(chunk.get('edge', []))
        for path in chunk.get('path', []):
            if path['name'] not in paths:
                paths[path['name']] = {**path, 'mapping': []}
            paths[path['name']]['mapping'].extend(path.get('mapping', []))
    return nodes, edges, paths


def signature(path):
    nodes, edges, paths = read_graph(path)
    occurrences = collections.defaultdict(list)
    spellings, metadata = {}, {}
    for name, path in paths.items():
        offset = 0
        spelling = []
        for step in path.get('mapping', []):
            node = int(step['position']['node_id'])
            reverse = step['position'].get('is_reverse', False)
            sequence = nodes[node]
            if reverse:
                sequence = complement(sequence)[::-1]
            spelling.append(sequence)
            for j in range(len(sequence)):
                forward = len(sequence) - 1 - j if reverse else j
                occurrences[node, forward].append((name, offset + j, -1 if reverse else 1))
            offset += len(sequence)
        spellings[name] = ''.join(spelling)
        metadata[name] = {k: v for k, v in path.items() if k != 'mapping'}
    keys, flips = {}, {}
    for node, sequence in nodes.items():
        for j, base in enumerate(sequence):
            members = tuple(sorted(occurrences[node, j]))
            assert members, 'unexpected pathless node'
            original = (base, members)
            reversed_key = (complement(base), tuple(sorted((n, p, -s) for n, p, s in members)))
            keys[node, j] = min(original, reversed_key)
            flips[node, j] = 1 if original <= reversed_key else -1
    canonical_edges = set()

    def add_edge(a, sa, b, sb):
        edge = (keys[a], sa * flips[a], keys[b], sb * flips[b])
        opposite = (edge[2], -edge[3], edge[0], -edge[1])
        canonical_edges.add(min(edge, opposite))

    for node, sequence in nodes.items():
        for j in range(len(sequence) - 1):
            add_edge((node, j), 1, (node, j + 1), 1)
    for edge in edges:
        a, b = int(edge['from']), int(edge['to'])
        reverse_a, reverse_b = edge.get('from_start', False), edge.get('to_end', False)
        add_edge((a, 0 if reverse_a else len(nodes[a]) - 1), -1 if reverse_a else 1,
                 (b, len(nodes[b]) - 1 if reverse_b else 0), -1 if reverse_b else 1)
    return spellings, metadata, collections.Counter(keys.values()), canonical_edges


def make_maf(path, seed, identical=False):
    rng = random.Random(seed)
    blocks, slot = 5, 25
    total = blocks * slot
    orders = {}
    for name in ['h0', 'h1', 'h2', 'h3', 'dup']:
        order = list(range(blocks))
        rng.shuffle(order)
        orders[name] = order
    lines = ['##maf version=1', '']
    for block in range(blocks):
        ref = ''.join(rng.choice('ACGTN') for _ in range(15))
        if seed % 3 == 0:
            ref = ''.join(base.lower() if j % 3 == 0 else base for j, base in enumerate(ref))
        if not identical:
            ref = ref[:7] + '--' + ref[7:]
        lines += ['a score=0', f's ref.chr {block * slot} 15 + {total} {ref}']
        for name in orders:
            for copy in range(2 if name == 'dup' else 1):
                aligned = list(ref)
                if not identical:
                    for j, base in enumerate(aligned):
                        if base == '-':
                            aligned[j] = rng.choice('ACGT') if rng.random() < .7 else '-'
                        elif rng.random() < .2:
                            aligned[j] = rng.choice('ACGTN-')
                aligned = ''.join(aligned)
                size = len(aligned.replace('-', ''))
                forward_start = orders[name][block] * slot + copy * total
                source_size = total * (2 if name == 'dup' else 1)
                reverse = not identical and rng.random() < .45
                start = source_size - forward_start - size if reverse else forward_start
                strand = '-' if reverse else '+'
                lines.append(f's {name}.chr {start} {size} {strand} {source_size} {aligned}')
        lines.append('')
    path.write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', default='./hal2vg')
    parser.add_argument('--maf2hal', default='./deps/hal/bin/maf2hal')
    parser.add_argument('--vg', default='vg')
    parser.add_argument('--work-dir', default='./low_memory_test')
    parser.add_argument('--hal', help='optional additional small HAL; never pass a full chromosome here')
    args = parser.parse_args()
    binary = str(pathlib.Path(args.binary).resolve())
    importer = str(pathlib.Path(args.maf2hal).resolve())
    work = pathlib.Path(args.work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    cases = []
    for seed in range(12):
        maf, hal = work / f'fixture_{seed}.maf', work / f'fixture_{seed}.hal'
        make_maf(maf, seed, identical=(seed == 0))
        subprocess.run([importer, str(maf), str(hal), '--refGenome', 'ref'], check=True)
        cases.append((f'fixture_{seed}', hal, []))
    # Ancestor exclusion and selected targets also change which virtual nodes are used.
    cases += [('no_ancestors', cases[1][1], ['--noAncestors']),
              ('ignore_root', cases[1][1], ['--ignoreGenomes', 'ref']),
              ('targets', cases[2][1], ['--targetGenomes', 'h0,h1']),
              ('references', cases[3][1], ['--refGenomes', 'ref,h0']),
              ('one_target', cases[2][1], ['--targetGenomes', 'h0']),
              ('chopped', cases[3][1], ['--chop', '3']),
              ('hash_graph', cases[4][1], ['--outputFormat', 'hg'])]
    renamed = work / 'phased_names.hal'
    shutil.copyfile(cases[5][1], renamed)
    rename_map = work / 'rename.tsv'
    rename_map.write_text('h0\tsample.1\nh1\tsample.2\n')
    rename_tool = str(pathlib.Path(importer).parent / 'halRenameGenomes')
    subprocess.run([rename_tool, str(renamed), str(rename_map)], check=True)
    cases.append(('phased_names', renamed, ['--refGenomes', 'ref']))
    pansn = work / 'pansn_names.hal'
    shutil.copyfile(cases[6][1], pansn)
    rename_map.write_text('h0\tsample#1\nh1\tsample#2\nh2\tassembly.1#0\nh3\tname.1J\n')
    subprocess.run([rename_tool, str(pansn), str(rename_map)], check=True)
    cases.append(('pansn_names', pansn, ['--refGenomes', 'ref']))
    if args.hal:
        cases.append(('provided_hal', pathlib.Path(args.hal).resolve(), []))
    for name, hal, flags in cases:
        signatures = []
        for mode in ['legacy', 'planned']:
            graph = work / f'{name}.{mode}.vg'
            log = work / f'{name}.{mode}.log'
            command = [binary, str(hal), '--inMemory', '--progress', *flags]
            if mode == 'legacy':
                command += ['--legacyConversion']
            else:
                # a tiny DNA window, so that blocks are read in many pieces
                command += ['--baseBuffer', '17']
            with graph.open('wb') as out, log.open('w') as err:
                subprocess.run(command, stdout=out, stderr=err, check=True)
            subprocess.run([args.vg, 'validate', str(graph)], check=True)
            exported = graph.with_suffix('.json')
            with exported.open('w') as out:
                subprocess.run([args.vg, 'view', '-j', str(graph)], stdout=out, check=True)
            signatures.append(signature(exported))
        assert signatures[0] == signatures[1], f'graph homology/path mismatch: {name}'
        print(f'PASS {name}: paths, metadata, base-level homology and oriented edges', flush=True)
    print(f'PASS {len(cases)} HAL comparisons; fixtures retained in {work}')


if __name__ == '__main__':
    main()
