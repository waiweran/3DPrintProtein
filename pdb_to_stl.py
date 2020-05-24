import mol_props
from stl_gen import Model3D

def make_space_fill(pdb_file, resolution=10):
    '''Creates STL file data from the given PDB file data.
       Bonds are filled in as cylinders and atoms are not included.
       resolution specifies how many faces are included in spheres and cylinders.
       Note that currently, adding sidechain hydrogens are not supported.'''
    model3d = Model3D()

    for model in pdb_file.models:
        for atom in model.mainchain_atoms:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'mainchain', 'atoms']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, resolution, tags)
        for atom in model.sidechain_atoms:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'sidechain', 'atoms']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, resolution, tags)
        for atom in model.hets_atoms:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'hets', 'atoms']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, resolution, tags)
        for atom in model.mainchain_hydrogens:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'mainchain', 'hydrogens']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, resolution, tags)
        for atom in model.sidechain_hydrogens:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'sidechain', 'hydrogens']
            model3d.add_sphere(mol_props.get_space_size(atom.element), atom.coords, resolution, tags)

    return model3d


def make_ball_stick(pdb_file, ball_size=0.4, ball_size_h=0.25, bond_width=0.2, bond_width_h=0.1, resolution=10):
    '''Creates STL file data from the given PDB file data.
       Bonds are filled in as cylinders and atoms are not included.
       ball_size specifies the diameter of atoms in the model.
       ball_size_h specifies diameter of hydrogen atoms.
       bond_width specifies the diameter of bond cylinders.
       bond_width_h specifies diameter of bonds with hydrogen.
       resolution specifies how many faces are included in spheres and cylinders.
       Note that currently, adding sidechain hydrogens are not supported.'''
    model3d = Model3D()

    for model in pdb_file.models:
        for atom in model.mainchain_atoms:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'mainchain', 'atoms']
            model3d.add_sphere(ball_size, atom.coords, resolution, tags)
        for connection in model.mainchain_connections:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'mainchain', 'atoms']
            model3d.add_cylinder(bond_width, connection[0].coords, connection[1].coords, resolution, tags)
        for atom in model.sidechain_atoms:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'sidechain', 'atoms']
            model3d.add_sphere(ball_size, atom.coords, resolution, tags)
        for connection in model.sidechain_connections:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'sidechain', 'atoms']
            model3d.add_cylinder(bond_width, connection[0].coords, connection[1].coords, resolution, tags)
        for atom in model.hets_atoms:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'hets', 'atoms']
            if(atom.element == 'H'):
                model3d.add_sphere(ball_size_h, atom.coords, resolution, tags)
            else:
                model3d.add_sphere(ball_size, atom.coords, resolution, tags)
        for connection in model.hets_connections:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'hets', 'atoms']
            if connection[0].element == 'H' or connection[1].element == 'H':
                model3d.add_cylinder(bond_width_h, connection[0].coords, connection[1].coords, resolution, tags)
            else:
                model3d.add_cylinder(bond_width, connection[0].coords, connection[1].coords, resolution, tags)
        for atom in model.mainchain_hydrogens:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'mainchain', 'hydrogens']
            model3d.add_sphere(ball_size_h, atom.coords, resolution, tags)
        for connection in model.mainchain_h_connections:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'mainchain', 'hydrogens']
            model3d.add_cylinder(bond_width_h, connection[0].coords, connection[1].coords, resolution, tags)
        for atom in model.sidechain_hydrogens:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, atom.chain_id), 'element_{}'.format(atom.element), 'sidechain', 'hydrogens']
            model3d.add_sphere(ball_size_h, atom.coords, resolution, tags)
        for connection in model.sidechain_h_connections:
            tags = ['model_{}'.format(model.serial_number), 'model_{}_chain_{}'.format(model.serial_number, connection[0].chain_id), 'element_{}'.format(atom.element), 'sidechain', 'hydrogens']
            model3d.add_cylinder(bond_width_h, connection[0].coords, connection[1].coords, resolution, tags)

    return model3d

def make_line_model(pdb_file, bond_width=0.2, bond_width_h=0.1, resolution=10):
    '''Creates STL file data from the given PDB file data.
       Bonds are filled in as cylinders and atoms are not included.
       bond_width specifies the diameter of bond cylinders.
       bond_width_h specifies diameter of bonds with hydrogen.
       resolution specifies how many faces are included in spheres and cylinders.
       Note that currently, adding sidechain hydrogens are not supported.'''
    return make_ball_stick(pdb_file, ball_size=bond_width, ball_size_h=bond_width_h, bond_width=bond_width, bond_width_h=bond_width_h, resolution=resolution)
