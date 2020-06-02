import mol_props
from stl_gen import Model3D

DEFAULT_ATOM_SIZE = 0.4
DEFAULT_H_SIZE = 0.25
DEFAULT_BOND_WIDTH = 0.2
DEFAULT_BOND_WIDTH_H = 0.1


def make_space_fill(pdb_file):
    """Creates STL file data from the given PDB file data.
       Bonds are filled in as cylinders and atoms are not included.
       resolution specifies how many faces are included in spheres and cylinders.
       Note that currently, adding sidechain hydrogens are not supported."""
    model3d = Model3D()

    for model in pdb_file.models:
        for atom in model.mainchain_atoms:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'mainchain', 'atoms']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, tags)
        for atom in model.sidechain_atoms:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'sidechain', 'atoms']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, tags)
        for atom in model.hets_atoms:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'het_{}'.format(atom.residue_code), 'hets', 'atoms']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, tags)
        for atom in model.mainchain_hydrogens:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'mainchain', 'hydrogens']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, tags)
        for atom in model.sidechain_hydrogens:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'sidechain', 'hydrogens']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, tags)

    return model3d


def make_ball_stick(pdb_file, atom_size=DEFAULT_ATOM_SIZE, hydrogen_size=DEFAULT_H_SIZE, bond_width=DEFAULT_BOND_WIDTH,
                    bond_width_h=DEFAULT_BOND_WIDTH_H):
    """Creates STL file data from the given PDB file data.
       Bonds are filled in as cylinders and atoms are not included.
       ball_size specifies the diameter of atoms in the model.
       ball_size_h specifies diameter of hydrogen atoms.
       bond_width specifies the diameter of bond cylinders.
       bond_width_h specifies diameter of bonds with hydrogen.
       resolution specifies how many faces are included in spheres and cylinders.
       Note that currently, adding sidechain hydrogens are not supported."""
    model3d = Model3D()

    for model in pdb_file.models:
        for atom in model.mainchain_atoms:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'mainchain', 'atoms']
            model3d.add_sphere(atom_size, atom.coords, tags)
        for connection in model.mainchain_connections:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'mainchain', 'atoms']
            model3d.add_cylinder(bond_width, connection[0].coords, connection[1].coords, tags)
        for atom in model.sidechain_atoms:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'sidechain', 'atoms']
            model3d.add_sphere(atom_size, atom.coords, tags)
        for connection in model.sidechain_connections:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'sidechain', 'atoms']
            model3d.add_cylinder(bond_width, connection[0].coords, connection[1].coords, tags)
        for atom in model.hets_atoms:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'het_{}'.format(atom.residue_code), 'hets']
            if atom.element == 'H':
                tags.append('hydrogens')
                model3d.add_sphere(hydrogen_size, atom.coords, tags)
            else:
                tags.append('atoms')
                model3d.add_sphere(atom_size, atom.coords, tags)
        for connection in model.hets_connections:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id)]
            if connection[0] in model.mainchain_atoms or connection[1] in model.mainchain_atoms:
                tags.append('mainchain')
            if connection[0] in model.sidechain_atoms or connection[1] in model.sidechain_atoms:
                tags.append('sidechain')
            if connection[0] in model.hets_atoms or connection[1] in model.hets_atoms:
                tags.append('hets')
                if connection[0] in model.hets_atoms:
                    tags.append('het_{}'.format(connection[0].residue_code))
                if connection[1] in model.hets_atoms:
                    tags.append('het_{}'.format(connection[1].residue_code))
            if connection[0].element == 'H' or connection[1].element == 'H':
                tags.append('hydrogens')
                model3d.add_cylinder(bond_width_h, connection[0].coords, connection[1].coords, tags)
            else:
                tags.append('atoms')
                model3d.add_cylinder(bond_width, connection[0].coords, connection[1].coords, tags)
        for atom in model.mainchain_hydrogens:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'mainchain', 'hydrogens']
            model3d.add_sphere(hydrogen_size, atom.coords, tags)
        for connection in model.mainchain_h_connections:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'mainchain', 'hydrogens']
            model3d.add_cylinder(bond_width_h, connection[0].coords, connection[1].coords, tags)
        for atom in model.sidechain_hydrogens:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, atom.chain_id),
                    'element_{}'.format(atom.element), 'sidechain', 'hydrogens']
            model3d.add_sphere(hydrogen_size, atom.coords, tags)
        for connection in model.sidechain_h_connections:
            tags = ['model_{}'.format(model.serial_number),
                    'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'sidechain', 'hydrogens']
            model3d.add_cylinder(bond_width_h, connection[0].coords, connection[1].coords, tags)

    return model3d


def make_line_model(pdb_file, bond_width=DEFAULT_BOND_WIDTH, bond_width_h=DEFAULT_BOND_WIDTH_H):
    """Creates STL file data from the given PDB file data.
       Bonds are filled in as cylinders and atoms are not included.
       bond_width specifies the diameter of bond cylinders.
       bond_width_h specifies diameter of bonds with hydrogen.
       resolution specifies how many faces are included in spheres and cylinders.
       Note that currently, adding sidechain hydrogens are not supported."""
    return make_ball_stick(pdb_file, atom_size=bond_width, hydrogen_size=bond_width_h, bond_width=bond_width,
                           bond_width_h=bond_width_h)
