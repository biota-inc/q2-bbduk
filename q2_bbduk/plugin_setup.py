# ----------------------------------------------------------------------------
# Copyright (c) 2024, BIOTA, inc..
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (
    Plugin, Choices, Citations, Float, Range, Int, Str, Bool, Threads
)
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality
)

import q2_bbduk
import q2_bbduk._trim

plugin = Plugin(
    name='bbduk',
    version=q2_bbduk.__version__,
    website='https://github.com/yourname/q2-bbduk',
    package='q2_bbduk',
    description='QIIME 2 plugin wrapping BBMap/BBDuk for k-mer based '
                'adapter/contaminant trimming and filtering.',
    short_description='BBDuk-based trimming/filtering for FASTQ.',
    citations=Citations.load('citations.bib', package='q2_bbduk')
)

# single-end
plugin.methods.register_function(
    function=q2_bbduk._trim.trim_single,
    inputs={'demultiplexed_sequences': SampleData[SequencesWithQuality]},
    parameters={
        'k': Int % Range(1, None),
        'ktrim': Str % Choices(['f', 'r', 'l']),
        'mink': Int % Range(0, None),
        'qtrim': Str % Choices(['rl', 'f', 'r', 'l', 'w']),
        'trimq': Float % Range(0, None, inclusive_start=True),
        'minlength': Int % Range(1, None),
        'ref': Str,              # "adapters,phix" or CSV of paths
        'literal': Str,
        'forcetrimleft': Int % Range(0, None),
        'forcetrimright': Int % Range(0, None),
        'kmask': Str,
        'threads': Threads,
    },
    outputs=[('trimmed_sequences', SampleData[SequencesWithQuality])],
    input_descriptions={
        'demultiplexed_sequences': 'Demultiplexed single-end reads.'
    },
    parameter_descriptions={
        'k': 'K-mer length.',
        'ktrim': 'Trim direction on k-mer matches: f, l, or r.',
        'mink': 'Minimum k at read tips.',
        'qtrim': 'Quality trimming mode: f/l/r/rl/w.',
        'trimq': 'Quality threshold for qtrim.',
        'minlength': 'Discard reads shorter than this length (post-trim).',
        'ref': 'Reference files or keywords passed to BBDuk (comma-delimited).',
        'literal': 'Literal adapter sequences (comma-delimited).',
        'forcetrimleft': 'Hard-trim bases left of this position.',
        'forcetrimright': 'Hard-trim bases right of this position.',
        'kmask': 'Mask matched k-mers with given symbol or lc.',
        'threads': 'Number of threads for BBDuk.',
    },
    output_descriptions={
        'trimmed_sequences': 'Trimmed single-end reads.'
    },
    name='Trim/filter demultiplexed single-end sequences using BBDuk.',
    description='Run bbduk.sh per-sample over demultiplexed single-end FASTQ.'
)

# paired-end
plugin.methods.register_function(
    function=q2_bbduk._trim.trim_paired,
    inputs={'demultiplexed_sequences': SampleData[PairedEndSequencesWithQuality]},
    parameters={
        'k': Int % Range(1, None),
        'ktrim': Str % Choices(['f', 'r', 'l']),
        'mink': Int % Range(0, None),
        'qtrim': Str % Choices(['rl', 'f', 'r', 'l', 'w']),
        'trimq': Float % Range(0, None, inclusive_start=True),
        'minlength': Int % Range(1, None),
        'tpe': Bool,
        'tbo': Bool,
        'minoverlap': Int % Range(1, None),
        'ref': Str,
        'literal': Str,
        'forcetrimleft': Int % Range(0, None),
        'forcetrimright': Int % Range(0, None),
        'kmask': Str,
        'threads': Threads,
    },
    outputs=[('trimmed_sequences', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions={
        'demultiplexed_sequences': 'Demultiplexed paired-end reads.'
    },
    parameter_descriptions={
        'k': 'K-mer length.',
        'ktrim': 'Trim direction on k-mer matches: f, l, or r.',
        'mink': 'Minimum k at read tips.',
        'qtrim': 'Quality trimming mode: f/l/r/rl/w.',
        'trimq': 'Quality threshold for qtrim.',
        'minlength': 'Discard reads shorter than this length (post-trim).',
        'tpe': 'Trim pairs evenly when right-trimming by k-mers.',
        'tbo': 'Trim adapters based on read overlap.',
        'minoverlap': 'Minimum overlap when using tbo.',
        'ref': 'Reference files or keywords passed to BBDuk (comma-delimited).',
        'literal': 'Literal adapter sequences (comma-delimited).',
        'forcetrimleft': 'Hard-trim bases left of this position.',
        'forcetrimright': 'Hard-trim bases right of this position.',
        'kmask': 'Mask matched k-mers with given symbol or lc.',
        'threads': 'Number of threads for BBDuk.',
    },
    output_descriptions={
        'trimmed_sequences': 'Trimmed paired-end reads.'
    },
    name='Trim/filter demultiplexed paired-end sequences using BBDuk.',
    description='Run bbduk.sh per-sample over demultiplexed paired-end FASTQ.'
)
