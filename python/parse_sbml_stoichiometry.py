import libsbml
import argparse
from numpy import empty
import ipdb

"""
Load an SBML model from file and print the stoichiometric matrix to stdout.
"""

def _parser():

    parser = argparse.ArgumentParser(description="Parse stoichiometry matrix of SBML file")
    parser.add_argument('file', metavar="filename", type=argparse.FileType('r'),
                        help="Filename of SBML file to parse")

    return parser

def _main():
    parser = _parser()
    args = parser.parse_args()

    file_ = args.file
    species, reactions, stoichiometry_matrix = parse_file(file_)
    _print_stoichiometry_matrix(species, reactions, stoichiometry_matrix)


def parse_file(open_file_):
    reader = libsbml.SBMLReader()
    document = reader.readSBML(open_file_.name)

    sbml_model = document.getModel()

    species = [s.getName() for s in sbml_model.getListOfSpecies()]
    reactions = [r.getId() for r in sbml_model.getListOfReactions()]
    parameters = [p.getName() for p in sbml_model.getListOfParameters()]

    # stoichiometry_matrix = {}
    Nr = len(sbml_model.getListOfReactions())
    Ns = len(sbml_model.getListOfSpecies())

    stoichiometry_matrix = empty((Ns,Nr))

    for reaction_index, reaction in enumerate(sbml_model.getListOfReactions()):
        reactants = {r.getSpecies(): r.getStoichiometry() for r in reaction.getListOfReactants()}
        products = {p.getSpecies(): p.getStoichiometry() for p in reaction.getListOfProducts()}

        for species_index, species_node in enumerate(sbml_model.getListOfSpecies()):
            species_id = species_node.getId()
            net_stoichiometry = int(products.get(species_id, 0)) - int(reactants.get(species_id, 0))
            stoichiometry_matrix[species_index, reaction_index] = net_stoichiometry
            
    # get rate laws as a lambda function
    reactionObjs = sbml_model.getListOfReactions()
    reactionFormulas = [r.getKineticLaw().getFormula() for r in reactionObjs]

    return species, reactions, parameters, reactionFormulas, stoichiometry_matrix
        

def _print_stoichiometry_matrix(species, reactions, stoichiometry_matrix):
    print('\t'.join(['---'] + reactions))
    for species_ix, species_label in enumerate(species):
        to_print = [species_label]
        for reaction_ix in range(len(reactions)):
            stoichiometry = stoichiometry_matrix[species_ix, reaction_ix]
            to_print.append(stoichiometry)

        print('\t'.join(map(str, to_print)))

if __name__ == '__main__':
    _main()

