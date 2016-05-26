import Geometry, PeptideBuilder, plotmap, os

A = 75.0
B = 135.0
positionOne  = (-A,B)
positionTwo  = (-B,A)
positionOneP = (B,-A)
positionTwoP = (A,-B)

N = 20

def build_structure(angles):
	geo = Geometry.geometry ("G")
	geo.phi     = angles[0][0]
	geo.psi_im1 = angles[0][1] # Ignored in the first instance
	geo.omega   = angles[0][2] # Ignored in the first instance
	
	structure = PeptideBuilder.initialize_res(geo)
	for i in range(1,len(angles)):
		geo = Geometry.geometry ("G")
		geo.phi     = angles[i][0] # current phi
		geo.psi_im1 = angles[i-1][1] # previous residue's psi
		geo.omega   = angles[i-1][2] # Ignored in the first instance
		structure = PeptideBuilder.add_residue(structure, geo)
	return structure

#structure_start = PeptideBuilder.initialize_res("G", -60, -40)

# FIRST, we compare structures created by just one value repeated N times
PHIPSIs = [positionOne,positionTwo,positionOneP,positionTwoP]
names = ["i ","j ","ip","jp"]

f = open("study_binaryness.output","w")

f.write("\n")
f.write("A:"+str(A)+"\n")
f.write("B:"+str(B)+"\n")
f.write("N:"+str(N)+"\n")

f.write("\n")
f.write("----- \t ---    \t-\t --- \n")
f.write("index:\t(phi    \t,\t psi)\n")
f.write("----- \t ---    \t-\t --- \n")
for i in range(len(PHIPSIs)):
	f.write(str(names[i])+": \t ("+str(PHIPSIs[i][0])+"  \t,\t"+str(PHIPSIs[i][1])+"\t)\n")

f.write("\n\t|-----------------------------")
f.write("\n\t|")
for i in range(len(PHIPSIs)):
	f.write(str(names[i])+"\t")
f.write("\n\t|-----------------------------\n")
atoms_to_consider = range(0,N+1)
for i in range(len(PHIPSIs)):
	phi1,psi1 = PHIPSIs[i]
	angles = []
	for n in range(N):
		angles.append((phi1,psi1,180.0))
	st1 = build_structure(angles)
	rg1 = plotmap.calculate_rg(st1,atoms_to_consider)
	f.write(str(names[i])+"\t|")
	for j in range(len(PHIPSIs)):
		phi2,psi2 = PHIPSIs[j]
		angles = []
		for n in range(N):
			angles.append((phi2,psi2,180.0))
		st2 = build_structure(angles)
		
		rmsd, st2_aligned = plotmap.calculate_rmsd(st1,st2,atoms_to_consider) # this ignores the first and last residues
		rg2   = plotmap.calculate_rg(st2,atoms_to_consider)
		
		
		alignment = ""
		
		basedir = "alignments"
		if not os.path.isdir(basedir):
			os.makedirs(basedir)
		basedir = os.path.abspath(basedir)
		
		"""
		import Bio.PDB # import Biopython's PDB module .
		out = Bio.PDB.PDBIO()
		out.set_structure(st1[0])
		out.save("aligned1.pdb")
		out.set_structure(st2[0])
		out.save("aligned2.pdb")
		os.system("sed -i 's/ A / B /g' aligned2.pdb")
		os.system('head -n -2 aligned1.pdb > aligned1b.pdb; mv aligned1b.pdb aligned1.pdb')
		os.system('head -n -2 aligned2.pdb > aligned2b.pdb; mv aligned2b.pdb aligned2.pdb')
		os.system('cat aligned1.pdb aligned2.pdb > aligned.pdb')

		pdbbase = basedir+"/st_"+names[i].rstrip().lstrip()+"_vs_"+names[j].rstrip().lstrip()
		
		os.system("cp aligned.pdb "+pdbbase+".pdb")
		command = "vmd -e compare_template.vmd -args "+pdbbase
		os.system(command)
		os.system("convert -trim "+pdbbase+".dat.png "+pdbbase+".dat.png")
		#os.system("eog "+pdbbase+".dat.png")
		#raw_input()
		"""
		
		rmsd = round(rmsd,2)
		
		#if i == j:
		#	rmsd = "-"
		
		f.write(str(rmsd)+"\t")
	f.write("\n")
f.write("\n\n")

PHIPSIs = [(positionOne,positionTwo),(positionOneP,positionTwoP),(positionTwo,positionOne),(positionTwoP,positionOneP)]
names = ["i-j","ip-jp","j-i","jp-ip"]

f.write("\n")
f.write("----- \t   ------         \t \t ------ \n")
f.write("index:\t(  state1         \t,\t state2)\n")
f.write("----- \t   ------         \t \t ------ \n")
for i in range(len(PHIPSIs)):
	f.write(str(names[i])+": \t ("+str(PHIPSIs[i][0])+"  \t,\t"+str(PHIPSIs[i][1])+")\n")

f.write("\n\t|--------------------------------------------------------------")
f.write("\n\t|")
for i in range(len(PHIPSIs)):
	f.write(str(names[i])+"\t \t")
f.write("\n\t|")
for i in range(len(PHIPSIs)):
	f.write("full\t skew\t")
f.write("\n\t|--------------------------------------------------------------\n")
atoms_to_consider = range(0,N+1)
for i in range(len(PHIPSIs)):
	phipsis1 = PHIPSIs[i]
	angles = []
	for n in range(N):
		phi1 = phipsis1[1][0]
		psi1 = phipsis1[1][1]
		if n%2 == 0:
			phi1 = phipsis1[0][0]
			psi1 = phipsis1[0][1]
		angles.append((phi1,psi1,180.0))
	st1 = build_structure(angles)
	rg1 = plotmap.calculate_rg(st1,atoms_to_consider)
	f.write(str(names[i])+"\t|")
	
	for j in range(len(PHIPSIs)):
		phipsis2 = PHIPSIs[j]
		angles = []
		for n in range(N):
			phi2 = phipsis2[1][0]
			psi2 = phipsis2[1][1]
			if n%2 == 0:
				phi2 = phipsis2[0][0]
				psi2 = phipsis2[0][1]
			angles.append((phi2,psi2,180.0))
		
		st2 = build_structure(angles)
		
		for range1,range2,comparisontype in [(range(1,N+1),range(1,N+1),"all"),(range(1,N),range(2,N+1),"translated")]:
			rmsd, st2_aligned = plotmap.calculate_rmsd(st1,st2,range1,range2)
			rg2   = plotmap.calculate_rg(st2,atoms_to_consider)
			
			alignment = ""
			basedir = "alignments"
			if not os.path.isdir(basedir):
				os.makedirs(basedir)
			basedir = os.path.abspath(basedir)
			
			
			import Bio.PDB # import Biopython's PDB module .
			out = Bio.PDB.PDBIO()
			out.set_structure(st1[0])
			out.save("aligned1.pdb")
			out.set_structure(st2[0])
			out.save("aligned2.pdb")
			os.system("sed -i 's/ A / B /g' aligned2.pdb")
			os.system('head -n -2 aligned1.pdb > aligned1b.pdb; mv aligned1b.pdb aligned1.pdb')
			os.system('head -n -2 aligned2.pdb > aligned2b.pdb; mv aligned2b.pdb aligned2.pdb')
			os.system('cat aligned1.pdb aligned2.pdb > aligned.pdb')
			
			pdbbase = basedir+"/binary_"+names[i].rstrip().lstrip()+"_vs_"+names[j].rstrip().lstrip()+"_"+comparisontype
			
			os.system("cp aligned.pdb "+pdbbase+".pdb")
			command = "vmd -e compare_template.vmd -args "+pdbbase
			os.system(command)
			os.system("convert -trim "+pdbbase+".dat.png "+pdbbase+".dat.png")
			os.system("eog "+pdbbase+".dat.png")
			raw_input()
			
			rmsd = round(rmsd,2)
			f.write(str(rmsd)+"\t")
	f.write("\n")
f.write("\t|--------------------------------------------------------------\n")
f.write("\n\n")

