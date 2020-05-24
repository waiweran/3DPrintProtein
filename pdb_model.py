class PDBFile:
    '''Represents a PDB file'''

    def __init__(self):
        self.models = [PDBModel()]
        self.molecule_name = ''
        self.molecule_source = ''
        self.alpha_helices = list()
        self.beta_sheets = list()

    def load_pdb(self, pdb_data):
        '''Loads data from a PDB file, passed in as a list of lines'''
        for line in pdb_data:
            self._process_line(line)
        for model in self.models:
            model._load_residue_connections()
            model._load_peptide_bonds()

    def _process_line(self, line):
        '''Processes one line of a PDB file'''
        key = line[0:6]
        # Title Section
        if key == 'COMPND':
            pass # Contains molecule name
        elif key =='SOURCE':
            pass # Contains organism
        # Secondary Structure Section
        elif key == 'HELIX ':
            self.alpha_helices.append(Helix(line))
        elif key == 'SHEET ':
            self.beta_sheets.append(Sheet(line))
        # Connectivity Annotation Section
        elif key == 'SSBOND':
            pass # Recategorize CONECT from hets_connections to sidechain_connections
        elif key == 'LINK  ':
            pass # Recategorize CONECT from hets_connections to mainchain/sidechain_connections
        # Coordinates Section
        elif key == 'ATOM  ':
            self.models[-1]._add_atom(line)
        elif key == 'HETATM':
            self.models[-1]._add_het(line)
        elif key == 'MODEL ':
            self._start_model(line)
        # Connectivity Section
        elif key == 'CONECT':
            self.models[-1]._add_extra_connection(line)

    def _start_model(self, line):
        '''Starts recording a model from the PDB file'''
        serial_number = int(line[10:14])
        if not self.models[-1].serial_number:
            self.models = list()
        self.models.append(PDBModel(serial_number))


class Helix:
    '''Represents an alpha helix from a PDB file'''

    def __init__(self, line):
        '''Reads data from a HELIX entry of a PDB file'''
        if line[0:6] == 'HELIX ':
            self.record_type = line[0:6].strip()
            self.serial_number = int(line[7:10])
            self.helix_id = line[11:14].strip()
            self.initial_residue_name = line[15:18].strip()
            self.initial_residue_chain = line[19].strip()
            self.initial_residue_number = int(line[21:25])
            self.initial_insertion_code = line[25].strip()
            self.end_residue_name = line[27:30].strip()
            self.end_residue_chain = line[31].strip()
            self.end_residue_number = int(line[33:37])
            self.end_insertion_code = line[37].strip()
            self.helix_class = int(line[38:40])
            self.comment = line[40:70].strip()
            self.length = int(line[71:76])
        else:
            raise ValueError('Given line is not a helix')


class Sheet:
    '''Represents a beta sheet from a PDB file'''

    def __init__(self, line):
        '''Reads data from a SHEET entry of a PDB file'''
        if line[0:6] == 'SHEET ':
            self.record_type = line[0:6].strip()
            self.strand_number = int(line[7:10])
            self.sheet_id = line[11:14].strip()
            self.num_strands = int(line[14:16])
            self.initial_residue_name = line[17:20].strip()
            self.initial_residue_chain = line[21].strip()
            self.initial_residue_number = int(line[22:26])
            self.initial_insertion_code = line[26].strip()
            self.end_residue_name = line[28:31].strip()
            self.end_residue_chain = line[32].strip()
            self.end_residue_number = int(line[34:37])
            self.end_insertion_code = line[37].strip()
        else:
            raise ValueError('Given line is not a sheet')


class PDBModel:
    '''Represents a model from a PDB file'''

    def __init__(self, serial=0):
        self.serial_number = serial
        self.chains = list()
        self.all_atoms = list()
        self.mainchain_atoms = list()
        self.mainchain_hydrogens = list()
        self.sidechain_atoms = list()
        self.sidechain_hydrogens = list()
        self.hets_atoms = list()
        self.mainchain_connections = list()
        self.mainchain_h_connections = list()
        self.sidechain_connections = list()
        self.sidechain_h_connections = list() # Not fully filled in
        self.hets_connections = list()

    def _add_atom(self, line):
        '''Adds an ATOM entry to the model representation'''
        atom = Atom()._from_line(line)
        if atom.alt_loc_id in ('', 'A'): # Currently forces alt A only
            self.all_atoms.append(atom)
            if atom.chain_id not in self.chains:
                self.chains.append(atom.chain_id)
            if atom.atom_name in ('N', 'C', 'CA', 'O'):
                self.mainchain_atoms.append(atom)
            elif atom.atom_name in ('H', 'HA'):
                self.mainchain_hydrogens.append(atom)
            elif atom.element == 'H':
                self.sidechain_hydrogens.append(atom)
            else:
                self.sidechain_atoms.append(atom)

    def _add_het(self, line):
        '''Adds a HETATM entry to the model representation'''
        atom = Atom()._from_line(line)
        if atom.alt_loc_id in ('', 'A'): # Currently forces alt A only
            self.all_atoms.append(atom);
            self.hets_atoms.append(atom)

    def _add_extra_connection(self, line):
        '''Adds a CONECT record to the model representation'''
        connections = list()
        serial_number = int(line[6:11])
        for atom in self.all_atoms:
            if atom.serial_number == serial_number:
                break
        if line[11:16].strip() != '':
            serial_number_2 = int(line[11:16])
            for atom2 in self.all_atoms:
                if atom2.serial_number == serial_number_2:
                    connections.append((atom, atom2))
                    break
        if line[16:21].strip() != '':
            serial_number_2 = int(line[16:21])
            for atom2 in self.all_atoms:
                if atom2.serial_number == serial_number_2:
                    connections.append((atom, atom2))
                    break
        if line[21:26].strip() != '':
            serial_number_2 = int(line[21:26])
            for atom2 in self.all_atoms:
                if atom2.serial_number == serial_number_2:
                    connections.append((atom, atom2))
                    break
        if line[27:31].strip() != '':
            serial_number_2 = int(line[26:31])
            for atom2 in self.all_atoms:
                if atom2.serial_number == serial_number_2:
                    connections.append((atom, atom2))
                    break
        for atom, atom2 in connections:
            if atom in self.mainchain_atoms:
                if atom2 in self.mainchain_atoms and not (atom2, atom) in self.mainchian_connections:
                    self.mainchain_connections.append((atom, atom2))
                elif atom2 in self.sidechain_atoms and not (atom2, atom) in self.sidechain_connections:
                    self.sidechain_connections.append((atom, atom2))
            if atom in self.sidechain_atoms and (atom2 in self.mainchain_atoms or atom2 in self.sidechain_atoms):
                if not (atom2, atom) in self.sidechain_connections:
                    self.sidechain_connections.append((atom, atom2))
            else:
                if not (atom2, atom) in self.sidechain_connections:
                    self.hets_connections.append((atom, atom2))

    def _load_residue_connections(self):
        '''Adds default bonds within residues'''
        index = 0
        while index < len(self.all_atoms):
            residue_num = self.all_atoms[index].residue_number
            i_code = self.all_atoms[index].insertion_code
            res_atoms = dict()
            residue_type = self.all_atoms[index].residue_code
            while index < len(self.all_atoms) and self.all_atoms[index].residue_number == residue_num and self.all_atoms[index].insertion_code == i_code:
                res_atoms[self.all_atoms[index].atom_name] = self.all_atoms[index]
                index += 1

            self._try_connection_add(self.mainchain_connections, res_atoms, 'CA', 'N')
            self._try_connection_add(self.mainchain_connections, res_atoms, 'CA', 'C')
            self._try_connection_add(self.mainchain_connections, res_atoms, 'C', 'O')
            self._try_connection_add(self.mainchain_h_connections, res_atoms, 'N', 'H')
            self._try_connection_add(self.mainchain_h_connections, res_atoms, 'CA', 'HA')
            if residue_type == "GLY":
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CA', 'HA2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CA', 'HA3')
            if residue_type == "ALA":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
            if residue_type == "VAL":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HG11')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HB12')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HB13')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HB21')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HB22')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HB23')
            if residue_type == "LEU":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD11')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD12')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD13')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD21')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD22')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD23')
            if residue_type == "ILE":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG1', 'CD1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HG12')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HG13')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG21')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG22')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG23')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD11')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD12')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD13')
            if residue_type == "MET":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'SD')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'SD', 'CE')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE3')
            if residue_type == "PHE":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD1', 'CE1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE1', 'CZ')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE2', 'CZ')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE1', 'HE1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE2', 'HE2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CZ', 'HZ')
            if residue_type == "TYR":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD1', 'CE1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE1', 'CZ')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE2', 'CZ')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CZ', 'OH')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE1', 'HE1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE2', 'HE2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'OH', 'HH')
            if residue_type == "SER":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'OG')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'OG', 'HG')
            if residue_type == "THR":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'OG1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'OG1', 'HG1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG21')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG22')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG23')
            if residue_type == "CYS":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'SG')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'SG', 'HG')
            if residue_type == "ASP":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'OD1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'OD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'OD2', 'HD2')
            if residue_type == "ASN":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'OD1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'ND2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'ND2', 'HD21')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'ND2', 'HD22')
            if residue_type == "GLU":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'OE1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'OE2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'OE2', 'HE2')
            if residue_type == "GLN":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'OE1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'NE2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NE2', 'HE21')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NE2', 'HE22')
            if residue_type == "LYS":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'CE')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE', 'NZ')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NZ', 'HZ1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NZ', 'HZ2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NZ', 'HZ3')
            if residue_type == "ARG":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'NE')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'NE', 'CZ')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CZ', 'NH1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CZ', 'NH2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NH1', 'HH11')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NH1', 'HH12')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NH2', 'HH21')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NH2', 'HH22')
            if residue_type == "TRP":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD1', 'NE1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE3')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'NE1', 'CE2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE2', 'CZ2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE3', 'CZ3')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CZ2', 'CH2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CZ3', 'CH2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'NE1', 'HE1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE3', 'HE3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CZ2', 'HZ2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CZ3', 'HZ3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CH2', 'HH2')
            if residue_type == "HIS":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'ND1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'ND1', 'CE1')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'NE2')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CE1', 'NE2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'ND1', 'HD1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE1', 'HE1')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CE2', 'HE2')
            if residue_type == "PRO":
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                self._try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'N')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD2')
                self._try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD3')

    def _try_connection_add(self, list, dict, key1, key2):
        '''Helper function to add a connection only if both atoms are present'''
        if key1 in dict and key2 in dict:
            list.append((dict[key1], dict[key2]))

    def _load_peptide_bonds(self):
        '''Adds peptide bonds between residues'''
        c_atom = None
        for atom in self.all_atoms:
            if atom.atom_name == 'C':
                c_atom = atom
            elif atom.atom_name == 'N' and c_atom and atom.chain_id == c_atom.chain_id:
                self.mainchain_connections.append((c_atom, atom))

class Atom:
    '''Represents an atom from a PDB file'''

    def __init__(self):
        pass
    
    def _from_line(self, line):
        '''Reads data from an ATOM or HETATM entry in a PDB file'''
        if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
            self.record_type = line[0:6].strip()
            self.serial_number = int(line[6:11])
            self.atom_name = line[12:16].strip()
            self.alt_loc_id = line[16].strip()
            self.residue_code = line[17:20].strip()
            self.chain_id = line[21].strip()
            self.residue_number = int(line[22:26])
            self.insertion_code = line[26]
            self.coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            self.occupancy = float(line[54:60])
            self.temp_factor = float(line[60:66])
            self.element = line[76:78].strip()
            self.charge = line[78:80].strip()
            return self
        else:
            raise ValueError('Given line is not an atom')
