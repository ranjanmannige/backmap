import os

current_filename = "backmap.tex"
target_filename_base = "backmap_"
target_dir = "for_upload"

old_to_new_fig_names = {
'figures/peptide_small':'fig1',
'automated_figures/fig_rama_intro':'fig2',
'automated_figures/fig_ss_2d_1d':'fig3',
'automated_figures/fig_r_intro':'fig4',
'figures/fig2_metrics':'fig5',
'figures/nano_r_demo':'fig6',
'automated_figures/fig_ramachandran_plots_vs_numbers':'fig7',
'automated_figures/fig_ramachandran_numbers_are_useful1':'fig8',
'automated_figures/example1':'figA',
'automated_figures/example2':'figB',
'figures/cmaps':'fig9',
'figures/1xqq_spread':'fig10',
'figures/2fft_spread':'fig11',
'figures/signed4':'fig12',
'automated_figures/fig_R_vs_R2':'fig13'
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
	
os.system('cp wlpeerj.* '+target_dir+'/')
os.system('cp all.bib   '+target_dir+'/')

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