# ----------------------------------------------------------------------------
# Copyright (c) 2024, BIOTA, inc..
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import tempfile
from pathlib import Path
import pandas as pd


from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)

def run_commands(cmds, verbose=True):
    for cmd in cmds:
        print('\nCommand:', end=' ')
        print(' '.join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


# BBDuk の代表的パラメータ。必要に応じて拡張してください。
_trim_defaults = {
    'threads': 1,          # bbduk.sh の threads
    'k': 27,
    'ktrim': 'f',          # f/l/r
    'mink': None,
    'qtrim': 'f',          # f/l/r/rl/w
    'trimq': 6.0,
    'minlength': 10,
    'tpe': False,
    'tbo': False,
    'minoverlap': 14,
    'ref': None,           # "adapters,phix" や "path1.fa,path2.fa"
    'literal': None,       # "ACGT,NNNN..."
    'forcetrimleft': 0,
    'forcetrimright': 0,
    'kmask': None,         # e.g., 'lc' or 'X'
    'overwrite': True,
}

def _build_bbduk_cmd(
    in1,
    out1,
    in2=None,
    out2=None,
    threads=_trim_defaults['threads'],
    k=_trim_defaults['k'],
    ktrim=_trim_defaults['ktrim'],
    mink=_trim_defaults['mink'],
    qtrim=_trim_defaults['qtrim'],
    trimq=_trim_defaults['trimq'],
    minlength=_trim_defaults['minlength'],
    tpe=_trim_defaults['tpe'],
    tbo=_trim_defaults['tbo'],
    minoverlap=_trim_defaults['minoverlap'],
    ref=_trim_defaults['ref'],
    literal=_trim_defaults['literal'],
    forcetrimleft=_trim_defaults['forcetrimleft'],
    forcetrimright=_trim_defaults['forcetrimright'],
    kmask=_trim_defaults['kmask'],
    overwrite=_trim_defaults['overwrite'],
    stats_fp=None,
):
    cmd = [
        'bbduk.sh',
        f'in={in1}',
        f'out={out1}',
        f'k={k}',
        f'ktrim={ktrim}',
        f'qtrim={qtrim}',
        f'trimq={trimq}',
        f'minlength={minlength}',
        f'threads={threads}',
        f'overwrite={"t" if overwrite else "f"}',
        'ziplevel=2',
        'ordered=f',
    ]

    if mink is not None and mink > 0:
        cmd += [f'mink={mink}']

    if in2 is not None and out2 is not None:
        cmd += [f'in2={in2}', f'out2={out2}']
        if tpe:
            cmd += ['tpe=t']
        if tbo:
            # tbo を使うなら minoverlap も合わせて渡す
            cmd += ['tbo=t', f'minoverlap={minoverlap}']

    if ref:
        cmd += [f'ref={ref}']
    if literal:
        cmd += [f'literal={literal}']
    if forcetrimleft and forcetrimleft > 0:
        cmd += [f'forcetrimleft={forcetrimleft}']
    if forcetrimright and forcetrimright > 0:
        cmd += [f'forcetrimright={forcetrimright}']
    if kmask:
        cmd += [f'kmask={kmask}']
    if stats_fp:
        cmd += [f'stats={stats_fp}']

    return cmd


def trim_single(
    demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    k: int = _trim_defaults['k'],
    ktrim: str = _trim_defaults['ktrim'],
    mink: int = _trim_defaults['mink'],
    qtrim: str = _trim_defaults['qtrim'],
    trimq: float = _trim_defaults['trimq'],
    minlength: int = _trim_defaults['minlength'],
    ref: str = _trim_defaults['ref'],
    literal: str = _trim_defaults['literal'],
    forcetrimleft: int = _trim_defaults['forcetrimleft'],
    forcetrimright: int = _trim_defaults['forcetrimright'],
    kmask: str = _trim_defaults['kmask'],
    threads: int = _trim_defaults['threads'],
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    cmds = []
    with tempfile.TemporaryDirectory() as td:
        for _, fwd in df.itertuples():
            # 出力ファイル名は入力の basename を維持
            out1 = trimmed_sequences.path / os.path.basename(fwd)
            stats_fp = Path(td) / (os.path.basename(fwd) + '.bbduk.stats.txt')
            cmds.append(_build_bbduk_cmd(
                in1=fwd, out1=str(out1),
                threads=threads, k=k, ktrim=ktrim, mink=mink,
                qtrim=qtrim, trimq=trimq, minlength=minlength,
                ref=ref, literal=literal,
                forcetrimleft=forcetrimleft, forcetrimright=forcetrimright,
                kmask=kmask, stats_fp=str(stats_fp),
            ))

        run_commands(cmds)

    return trimmed_sequences


def trim_paired(
    demultiplexed_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
    k: int = _trim_defaults['k'],
    ktrim: str = _trim_defaults['ktrim'],
    mink: int = _trim_defaults['mink'],
    qtrim: str = _trim_defaults['qtrim'],
    trimq: float = _trim_defaults['trimq'],
    minlength: int = _trim_defaults['minlength'],
    tpe: bool = _trim_defaults['tpe'],
    tbo: bool = _trim_defaults['tbo'],
    minoverlap: int = _trim_defaults['minoverlap'],
    ref: str = _trim_defaults['ref'],
    literal: str = _trim_defaults['literal'],
    forcetrimleft: int = _trim_defaults['forcetrimleft'],
    forcetrimright: int = _trim_defaults['forcetrimright'],
    kmask: str = _trim_defaults['kmask'],
    threads: int = _trim_defaults['threads'],
) -> CasavaOneEightSingleLanePerSampleDirFmt:

    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    cmds = []
    with tempfile.TemporaryDirectory() as td:
        for _, fwd, rev in df.itertuples():
            out1 = trimmed_sequences.path / os.path.basename(fwd)
            out2 = trimmed_sequences.path / os.path.basename(rev)
            stats_fp = Path(td) / (os.path.basename(fwd) + '.bbduk.stats.txt')
            cmds.append(_build_bbduk_cmd(
                in1=fwd, out1=str(out1),
                in2=rev, out2=str(out2),
                threads=threads, k=k, ktrim=ktrim, mink=mink,
                qtrim=qtrim, trimq=trimq, minlength=minlength,
                tpe=tpe, tbo=tbo, minoverlap=minoverlap,
                ref=ref, literal=literal,
                forcetrimleft=forcetrimleft, forcetrimright=forcetrimright,
                kmask=kmask, stats_fp=str(stats_fp),
            ))

        run_commands(cmds)

    return trimmed_sequences
