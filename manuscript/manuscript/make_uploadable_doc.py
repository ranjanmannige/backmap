import os

current_filename = "ramachandran_v7.tex"
target_filename_base = "ramachandran_number_"
target_dir = "for_upload"


old_to_new_fig_names = {
"figures/fig1":  "fig1",
"figures/fig2":  "fig2",
"figures/pathogenicity2":  "fig3",
"figures/fig_ss_line":  "fig4",
"figures/collage2":  "fig5",
"figures/fig_classes4":  "fig6",
"figures/binary_birth":  "fig7",
"figures/fig_coords2":  "fig8",
"figures/sigma_vs_rmsa_and_rmsd":  "fig9"
}

f = open(current_filename,"r")
block = f.read()
f.close()

if not os.path.isdir(target_dir):
	os.makedirs(target_dir)

for oldfn,newfn in old_to_new_fig_names.items():
	command = "cp "+oldfn+".pdf "+target_dir+"/"+target_filename_base+newfn+".pdf"
	print command
	os.system(command)
	
	block = block.replace(oldfn,target_filename_base+newfn)
fn = target_dir+"/"+target_filename_base.rstrip("_")+".tex"
f=open(fn,"w")
f.write(block)
f.close()

FILENAME = target_filename_base.rstrip("_")

runmetemplate = """import os
fn = "FILENAME"
os.system("pdflatex "+fn)
os.system("bibtex "+fn)
os.system("pdflatex "+fn)
os.system("pdflatex "+fn)
os.system("evince "+fn+".pdf")
"""

fn = target_dir+"/runme.py"
f = open(fn,"w")
f.write(runmetemplate.replace("FILENAME",FILENAME))
f.close()
print fn