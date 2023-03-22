import json

# ---------------------------------------------------------
# Basic

ref = 'NC_000913.3'
rname = 'NC_000913'
fracs = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
outfile = './ecoli_sim_basic.json'

jdata = {}

default = {
	"ncbinuc": ref,
	"mapping_ref": "True",
	"quantity": 20,
	"circular": "true",
	"nanopore": "average",
	"seed": 1
}

for idx, frac in enumerate(fracs, start=1):
	prc = float(frac * 100)
	name = f'simulated_{rname}_{idx}'

	mut = {key: value for key, value in default.items()
		if key != 'mapping_ref'}
	mut["simulate"] = {
		"snp_prc": prc,
		"indel_count": 0,
		"seed": 1
		}

	jdata[name] = {
		'ref': default,
		'mut': mut
	}

with open(outfile, 'w') as f:
	json.dump(jdata, f, indent=4)


# ---------------------------------------------------------
# Cov diff

ref = 'NC_000913.3'
rname = 'NC_000913'
fracs = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
outfile = './ecoli_sim_covdiff.json'

jdata = {}

default = {
	"ncbinuc": ref,
	"mapping_ref": "True",
	"quantity": 20,
	"circular": "true",
	"nanopore": "average",
	"seed": 1
}

for idx, frac in enumerate(fracs, start=1):
	prc = float(frac * 100)
	name = f'simulated_{rname}_{idx}'

	mut = {key: value for key, value in default.items()
		if key != 'mapping_ref'}

	mut['quantity'] = 30

	mut["simulate"] = {
		"snp_prc": prc,
		"indel_count": 0,
		"seed": 1
		}

	jdata[name] = {
		'ref': default,
		'mut': mut
	}

with open(outfile, 'w') as f:
	json.dump(jdata, f, indent=4)

# ---------------------------------------------------------
# Cov gradient

ref = 'NC_000913.3'
rname = 'NC_000913'
gradient = [5, 10, 15, 20, 25, 30]
outfile = './ecoli_sim_covgrad.json'

jdata = {}

default = {
	"ncbinuc": ref,
	"mapping_ref": "True",
	"quantity": 20,
	"circular": "true",
	"nanopore": "average",
	"seed": 1
}

for idx, cov in enumerate(gradient, start=1):
	name = f'simulated_{rname}_{idx}'

	mut = {key: value for key, value in default.items()
		if key != 'mapping_ref'}

	mut['quantity'] = cov

	mut["simulate"] = {
		"snp_prc": 2,
		"indel_count": 0,
		"seed": 1
		}

	jdata[name] = {
		'ref': default,
		'mut': mut
	}

with open(outfile, 'w') as f:
	json.dump(jdata, f, indent=4)

# ---------------------------------------------------------
# Cov gradient 2

ref = 'NC_000913.3'
rname = 'NC_000913'
gradient = [5, 10, 15, 20, 25, 30]
outfile = './ecoli_sim_covgrad_3strains.json'

jdata = {}

default = {
	"ncbinuc": ref,
	"mapping_ref": "True",
	"quantity": 20,
	"circular": "true",
	"nanopore": "average",
	"seed": 1
}

for idx, cov in enumerate(gradient, start=1):
	
	mut = {key: value for key, value in default.items()
		if key != 'mapping_ref'}

	mut['quantity'] = cov

	mut["simulate"] = {
		"snp_prc": 1,
		"indel_count": 0,
		"seed": 1
		}

	for idx2, cov2 in enumerate(gradient, start=1):
		name = f'simulated_{rname}_{idx}_{idx2}'

		mut2 = {key: value for key, value in default.items()
			if key != 'mapping_ref'}

		mut2.pop('ncbinuc')
		mut2['refpath'] = f"references/simulated/{name}/mut.simseq.genome.fa"
		mut2['quantity'] = cov2

		mut2["simulate"] = {
			"snp_prc": 1,
			"indel_count": 0,
			"seed": 1
			}

		jdata[name] = {
			'ref': default,
			'mut': mut,
			'mut2': mut2
	}

with open(outfile, 'w') as f:
	json.dump(jdata, f, indent=4)


# ---------------------------------------------------------
# Cov gradient 3 with klebsiella

ref = 'ASM1990034v1'
rname = 'kleb_cg'
gradient = [5, 10, 15, 20, 25, 30]
outfile = './kleb_covgrad_3strains.json'

jdata = {}

default = {
	"ncbinuc": ref,
	"mapping_ref": "True",
	"quantity": 15,
	"circular": "true",
	"nanopore": "average",
	"seed": 1
}

for idx, cov in enumerate(gradient, start=1):
	
	mut = {key: value for key, value in default.items()
		if key != 'mapping_ref'}

	mut['quantity'] = cov
	mut['ncbinuc'] = 'ASM1990030v1'

	for idx2, cov2 in enumerate(gradient, start=1):
		name = f'simulated_{rname}_{idx}_{idx2}'

		mut2 = {key: value for key, value in default.items()
			if key != 'mapping_ref'}

		mut2['ncbinuc'] = 'ASM1990032v1'
		mut2['quantity'] = cov2

		jdata[name] = {
			'kleb_MN': default,
			'kleb_NC': mut,
			'kleb_PA': mut2
	}

with open(outfile, 'w') as f:
	json.dump(jdata, f, indent=4)