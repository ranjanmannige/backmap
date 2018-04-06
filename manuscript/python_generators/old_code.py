"""
KEYS FOR DSSP:
G = 3-turn helix (3_10-helix). Min length 3 residues.
H = 4-turn helix (alpha-helix). Min length 4 residues.
I = 5-turn helix (pi-helix). Min length 5 residues.
E = extended strand in parallel and/or anti-parallel beta-sheet conformation. Min length 2 residues.
----- Less "exciting" motifs
T = hydrogen bonded turn (3, 4 or 5 turn)
B = residue in isolated beta-bridge (single pair beta-sheet hydrogen bond formation)
S = bend (the only non-hydrogen-bond based assignment).
C = coil (residues which are not in any of the above conformations).
"""
def get_resid_to_dssp(fn):
	# Create the dssp report
	dssp_output = fn+".dssp"
	os.system("./dssp-2.0.4-linux-amd64 -i "+str(fn)+" -o "+dssp_output)
	
	f = open(dssp_output)
	dssp_block = f.read()
	f.close()
	
	# populating the DSSP secondary structure records for this domain
	resno_to_sstype = {}
	dssp_block = dssp_block.split("#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA")[-1].rstrip()
	for line in dssp_block.split("\n"):
		line = line.rstrip()
		#['  304        !              0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0'] 
		if len(line):
			residue_str = line[5:10].lstrip().rstrip()
			if len(residue_str):
				#print [line],dsspfn
				#print line[5:10].lstrip().rstrip()
				residue = int(residue_str)
				ss_type = line[16]
				resno_to_sstype[residue]=ss_type
	
	return resno_to_sstype
#

