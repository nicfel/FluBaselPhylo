# JanusTree

Input files etc. to get a time tree (using treetime and raxml) from janus https://github.com/nextstrain/janus.

## install janus

download janus and associated code from github:
`git clone --recursive https://github.com/nextstrain/janus.git`

`cd janus`

`git submodule update --init --recursive`

Check out the code and branch.



`git submodule update --init --recursive`

Install Python environment.

`./install.sh`



## run augur

Once Janus is installed, use the following python commands to build a timetree:

* Install requirements with `pip install -r requirements.txt`. Also requires that `fasttree`, `raxml` and `mafft` are all callable from the command line. This means relabeling your RAxML binary to `raxml`. You probably want to use `raxmlHPC-PTHREADS-AVX`.
* Go to `flu/`.
* Move `h3n2_ha.fasta` into `flu/` directory.
* Prepare analysis with `python flu.prepare.py --sequences h3n2_ha.fasta -v 1000000000 --time_interval 2010-01-01 2018-01-01 --file_prefix h3n2_ha`
* Run with `python flu.process.py --no_mut_freqs --json prepared/h3n2_ha.json`
* This will produce JSON files in `flu/auspice/` that are prefixed with `flu_basel_`. (edited)
To visualize in auspice:
* Install requirements with `npm install`.
* Move JSONs into `data/` directory.
* Run with `npm run start:local`.
* Go to http://localhost:4000/flu/h3n2/ha/3y
Can convert JSON to NEWICK with `scripts/json_tree_to_nexus.py`.
* to get tree:
* to then get the tree, use: `python json_tree_to_nexus.py -t ../builds/flu/auspice/h3n2_ha_tree.json --temporal -o ../../time_tree.tree`
** to get the tree in genetic distance (i.e the raxml tree), use: `python json_tree_to_nexus.py -t ../builds/flu/auspice/h3n2_ha_tree.json -o ../../time_tree_raxml.tree`
* note that gettingt the tree requires for some reason a changed utils.py class
* to get a time resolved tree: `python json_tree_to_nexus.py -t ../builds/flu/auspice/h3n2_ha_tree.json --temporal --resolve -o ../../time_tree_resolved.tree`

