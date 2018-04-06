import os
fn = "plotmap"
os.system("pdflatex "+fn)
os.system("bibtex "+fn)
os.system("pdflatex "+fn)
os.system("pdflatex "+fn)

#fn+".bbl "+
os.system("rm "+fn+".blg "+fn+".aux "+fn+".log "+fn+".out ")

"""
manuscriptendingpage=5

os.system("pdftk "+fn+".pdf cat 1-"+str(manuscriptendingpage)+" output "+fn+"_MS"+".pdf")
os.system("pdftk "+fn+".pdf cat "+str(manuscriptendingpage+1)+"-end output "+fn+"_SI"+".pdf")
"""

os.system("evince "+fn+".pdf")
