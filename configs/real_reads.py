import pandas as pd
import json, os, glob, copy
from collections import defaultdict

wd = os.path.dirname(os.getcwd())

# Ecoli

fname = os.path.join(wd, 'dataset/isolates/ecoli_mtbird/accessions.tsv')
df = pd.read_csv(fname, sep='\t')
df = df[(df['Chromsome or plasmid'].str.strip() == 'Chromosome') | (df['Chromsome or plasmid'] == 'Chromsome')]
gb = df.set_index('Isolate')['Genbank accession'].str.strip().to_dict()
gb = {str(isolate): gbid for isolate, gbid in gb.items()}

nanopore = os.path.join(wd, 'dataset/isolates/ecoli_mtbird/nanopore/*/*.fastq.gz')
nanopore = sorted(glob.glob(nanopore))
nanopore = {os.path.basename(os.path.dirname(fname)): fname for fname in nanopore}

fnames = os.path.join(wd, 'dataset/isolates/ecoli_mtbird/illumina/*/*.cleaned.fastq.gz')
fnames = glob.glob(fnames)
illumina = defaultdict(list)

for fname in fnames:
    sample = os.path.basename(os.path.dirname(fname))
    illumina[sample].append(fname)

for sample, fnames in illumina.items():
    illumina[sample] = sorted(fnames)
    assert len(fnames) == 2

references = {
    isolate: os.path.join(wd, f'dataset/isolates/ecoli_mtbird/genbank/{isolate}/{gbi}.fasta.gz')
    for isolate, gbi in gb.items()
}

# Klebsiella and B. licheniformis

fname = os.path.join(wd, 'dataset/isolates/kpbl/config.json')
with open(fname) as f:
    jdata = json.load(f)

for species_name, samples in jdata.items():
    for sampleid, rsa_id in samples.items():
        nanopore[sampleid] = os.path.join(wd, f'dataset/isolates/kpbl/reads/{species_name}/{rsa_id}.fastq.gz')
        references[sampleid] = os.path.join(wd, f'dataset/isolates/kpbl/maincontigs/{species_name}/{sampleid}.fasta')


# Merge everything together in the configuration file

subgroups = {
    'ecoli_mtbird': ['786605', '899091', '542093', '623214'],
    'lklebsiella': ['N319_barcode04', 'N326_barcode18', 'N325_barcode01', 'N319_barcode12'],
    'blicheniformis': ['N413_barcode11', 'N419_barcode16', 'N419_barcode20', 'N413_barcode03']
}

refnames = {
    'ecoli_mtbird': '786605',
    'lklebsiella': 'N319_barcode04',
    'blicheniformis': 'N413_barcode11'
}

# Make the files

coverages = [25, 20, 15, 30]

jdata = {'samples': {}, 'outputs': {}}

for species, isolates in subgroups.items():
    for i in range(2, 5):
        name = f'{species}_{i}'
        config = {}
        for j in range(i):
            isolate = isolates[j]

            config[isolate] = {
                'refpath': references[isolate],
                'mapping_ref': "True" if isolate in refnames[species] else "False",
                'nanopore': nanopore[isolate],
                'seed': 50,
                'quantity': coverages[j]
            }

            if isolate in illumina:
                config[isolate]['illumina_r1'] = illumina[isolate][0]
                config[isolate]['illumina_r2'] = illumina[isolate][1]

            jdata['samples'][name] = config

# ----------------------------------------------------------------------
# Regular outputs

jdata['outputs'] = [
    {
        'group': 'all',
        'assembly': {
            'name': 'inpref'
        },
        'vcalling': 'longshot',
        'phasing': {
            'name': 'floria',
            'readtype': 'nanopore',
            'fmode': 'none',
            'post_assembler': 'flye',
            'assembler_rtype': 'long_reads',
            'assembler_preset': 'none'
        },
        'ref_filtering':{
            'length': 'none',
            'similarity': 'none'
        }
    },

    {
        'group': 'all',
        'assembly': {
            'name': 'inpref'
        },
        'vcalling': 'longshot',
        'phasing': {
            'name': 'floria',
            'readtype': 'nanopore',
            'fmode': 'none',
            'post_assembler': 'wtdbg2',
            'assembler_rtype': 'long_reads',
            'assembler_preset': 'nano'
        },
        'ref_filtering':{
            'length': 'none',
            'similarity': 'none'
        }
    },

    {
        'group': 'all',
        'assembly': {
            'name': 'inpref'
        },
        'vcalling': 'freebayes_ploidy_auto',
        'phasing': {
            'name': 'floria',
            'readtype': 'nanopore',
            'fmode': 'none',
            'post_assembler': 'flye',
            'assembler_rtype': 'long_reads',
            'assembler_preset': 'none'
        },
        'ref_filtering':{
            'length': 'none',
            'similarity': 'none'
        }
    },

    {
        'group': 'all',
        'assembly': {
            'name': 'inpref'
        },
        'vcalling': 'longshot',
        'phasing': {
            'name': 'strainberry',
            'readtype': 'nanopore'
        },
        'ref_filtering':{
            'length': 'none',
            'similarity': 'none'
        }
    },

    {
        'group': 'all',
        'assembly': {
            'name': 'flye',
            'read': 'nanopore',
            'mode': 'raw'
        },
        'ref_filtering':{
            'length': 'none',
            'similarity': 'none'
        }
    },

    {
        'group': 'all',
        'assembly': {
            'name': 'flye_single',
            'read': 'nanopore',
            'mode': 'raw'
        },
        'ref_filtering':{
            'length': 'none',
            'similarity': 'none'
        }
    }
]

jdata['miscs'] = {
    'GENERIC_THREADS': 8,
    'GENERIC_THREADS_SPLIT': 8
}

# ----------------------------------------------------------------------
# Floria multiple configurations

presets = {
    ** {f'fp_e{idx + 1}': f'-e {epsilon / 100}' for idx, epsilon in enumerate(range(1, 11, 1))},
    ** {f'fp_ps{idx + 1}': f'-s {ploidy_sensitivity}' for idx, ploidy_sensitivity in enumerate((1, 2, 3))}
}

jdata['miscs']['floria_presets'] = presets

template = {
    'group': 'all',
    'assembly': {
        'name': 'inpref'
    },
    'vcalling': 'longshot',
    'phasing': {
        'name': 'floria',
        'readtype': 'nanopore',
        'fmode': 'none',
        'post_assembler': 'flye',
        'assembler_rtype': 'long_reads',
        'assembler_preset': 'none'
    },
    'ref_filtering':{
        'length': 'none',
        'similarity': 'none'
    }
}

for preset in presets:
    preset_config = copy.deepcopy(template)
    preset_config['phasing']['fmode'] = preset
    jdata['outputs'].append(preset_config)

# ----------------------------------------------------------------------
# WhatsHap

wh_config = {
    'group': 'all',
    'assembly': {
        'name': 'inpref'
    },
    'vcalling': 'freebayes_ploidy_auto',
    'phasing': {
        'name': 'whatshap',
        'readtype': 'nanopore',
        'post_assembler': 'flye',
        'assembler_rtype': 'long_reads',
        'assembler_preset': 'none'
    },
    'ref_filtering':{
        'length': 'none',
        'similarity': 'none'
    }
}

jdata['outputs'].append(wh_config)
# jdata['outputs'] = [wh_config]

# ----------------------------------------------------------------------
# Saving the file on disk

outfile = './real_reads.json'
with open(outfile, 'w') as f:
    json.dump(jdata, f, indent=4, sort_keys=True)