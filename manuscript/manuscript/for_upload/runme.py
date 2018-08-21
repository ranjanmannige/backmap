import os
fn = "backmap"
os.system("pdflatex "+fn)
os.system("bibtex "+fn)
os.system("pdflatex "+fn)
os.system("pdflatex "+fn)
os.system("evince "+fn+".pdf")
