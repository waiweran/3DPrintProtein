from pdb_model import PDBFile
import pdb_to_stl


def make_stl(pdb_input, output_file, model_type, sidechains, hydrogens, hets, resolution=10):
    """Creates an STL file from the given PDB file data.
       STL file saved at the location specified by output_file
       The model_type can be 'space', 'ballstick' or 'stick' to produce a
       space-filling model, a ball-and-stick model, or a model that
       shows bonds but does not include atom balls.
       Boolean values for sidechains, hydrogens, and hets specify whether
       each of those items is included in the model.
       Note that currently, sidechain hydrogens are always left out of
       ball-and-stick and stick only models."""
    
    # Initialize
    pdb_file = PDBFile()
    pdb_file.load_pdb(pdb_input)

    # Select components for output
    exclude_tags = list()
    if not sidechains:
        exclude_tags.append('sidechain')
    if not hydrogens:
        exclude_tags.append('hydrogens')
    if not hets:
        exclude_tags.append('hets')

    # Make and Save Model
    if model_type == 'space':
        model = pdb_to_stl.make_space_fill(pdb_file, resolution=resolution)
    elif model_type == 'ballstick':
        model = pdb_to_stl.make_ball_stick(pdb_file, resolution=resolution)
    else:
        model = pdb_to_stl.make_line_model(pdb_file, resolution=resolution)
    model.to_stl(exclude_tags=exclude_tags)
    model.save(output_file)


if __name__ == '__main__':
    with open('Files/6j5i.pdb') as opened_file:
        pdb_input_lines = opened_file.readlines()
        make_stl(pdb_input_lines, 'Files/output.stl', 'space', True, False, False, resolution=10)
