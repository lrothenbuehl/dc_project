# This was a major pain!!!!!
# copied from https://stackoverflow.com/questions/75632445/anaconda-jupyter-notebook-and-missingidfieldwarning
import nbformat

with open("./02 - Data Curation/DB creation.ipynb", "r") as file:
    nb_corrupted = nbformat.reader.read(file)


nb_fixed = nbformat.validator.normalize(nb_corrupted)
nbformat.validator.validate(nb_fixed[1])


with open("./02 - Data Curation/DB creation normalized.ipynb", "w") as file:
    nbformat.write(nb_fixed[1], file)

