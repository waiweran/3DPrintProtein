class PDBFile:
    """Represents a PDB file"""

    def __init__(self):
        self.models = [PDBModel()]
        self.molecule_name = ''
        self.molecule_source = ''
        self.alpha_helices = list()
        self.beta_sheets = list()

    def load_pdb(self, pdb_data):
        """Loads data from a PDB file, passed in as a list of lines"""
        for line in pdb_data:
            self._process_line(line)
        for model in self.models:
            model.load_residue_connections()
            model.load_peptide_bonds()

    def _process_line(self, line):
        """Processes one line of a PDB file"""
        key = line[0:6]
        # Title Section
        if key == 'COMPND':
            pass  # Contains molecule name
        elif key == 'SOURCE':
            pass  # Contains organism
        # Secondary Structure Section
        elif key == 'HELIX ':
            self.alpha_helices.append(Helix(line))
        elif key == 'SHEET ':
            self.beta_sheets.append(Sheet(line))
        # Connectivity Annotation Section
        elif key == 'SSBOND':
            pass  # Disulfides re-categorized in processing of hets connections
        elif key == 'LINK  ':
            pass  # Bonds indicated here re-categorized in processing of hets connections
        # Coordinates Section
        elif key == 'ATOM  ':
            self.models[-1].add_atom(line)
        elif key == 'HETATM':
            self.models[-1].add_het(line)
        elif key == 'MODEL ':
            self._start_model(line)
        # Connectivity Section
        elif key == 'CONECT':
            self.models[-1].add_extra_connection(line)

    def _start_model(self, line):
        """Starts recording a model from the PDB file"""
        serial_number = int(line[10:14])
        if not self.models[-1].serial_number:
            self.models = list()
        self.models.append(PDBModel(serial_number))


class Helix:
    """Represents an alpha helix from a PDB file"""

    def __init__(self, line):
        """Reads data from a HELIX entry of a PDB file"""
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
    """Represents a beta sheet from a PDB file"""

    def __init__(self, line):
        """Reads data from a SHEET entry of a PDB file"""
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
    """Represents a model from a PDB file"""

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
        self.sidechain_h_connections = list()  # Not fully filled in
        self.hets_connections = list()

    def add_atom(self, line):
        """Adds an ATOM entry to the model representation"""
        atom = Atom(line)
        if atom.alt_loc_id in ('', 'A'):  # Currently forces alt A only
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

    def add_het(self, line):
        """Adds a HETATM entry to the model representation"""
        atom = Atom(line)
        if atom.alt_loc_id in ('', 'A'):  # Currently forces alt A only
            self.all_atoms.append(atom)
            self.hets_atoms.append(atom)

    def add_extra_connection(self, line):
        """Adds a CONECT record to the model representation"""
        connections = list()
        serial_number = int(line[6:11])
        atom = None
        for atm in self.all_atoms:
            if atm.serial_number == serial_number:
                atom = atm
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
                if atom2 in self.mainchain_atoms and not (atom2, atom) in self.mainchain_connections:
                    self.mainchain_connections.append((atom, atom2))
                elif atom2 in self.sidechain_atoms and not (atom2, atom) in self.sidechain_connections:
                    self.sidechain_connections.append((atom, atom2))
            if atom in self.sidechain_atoms and (atom2 in self.mainchain_atoms or atom2 in self.sidechain_atoms):
                if not (atom2, atom) in self.sidechain_connections:
                    self.sidechain_connections.append((atom, atom2))
            else:
                if not (atom2, atom) in self.sidechain_connections:
                    self.hets_connections.append((atom, atom2))

    def load_residue_connections(self):
        """Adds default bonds within residues"""
        index = 0
        while index < len(self.all_atoms):
            residue_num = self.all_atoms[index].residue_number
            i_code = self.all_atoms[index].insertion_code
            res_atoms = dict()
            residue_type = self.all_atoms[index].residue_code
            while index < len(self.all_atoms) and self.all_atoms[index].residue_number == residue_num \
                    and self.all_atoms[index].insertion_code == i_code:
                res_atoms[self.all_atoms[index].atom_name] = self.all_atoms[index]
                index += 1

            _try_connection_add(self.mainchain_connections, res_atoms, 'CA', 'N')
            _try_connection_add(self.mainchain_connections, res_atoms, 'CA', 'C')
            _try_connection_add(self.mainchain_connections, res_atoms, 'C', 'O')
            _try_connection_add(self.mainchain_h_connections, res_atoms, 'N', 'H')
            _try_connection_add(self.mainchain_h_connections, res_atoms, 'CA', 'HA')
            if residue_type == "GLY":
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CA', 'HA2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CA', 'HA3')
            if residue_type == "ALA":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
            if residue_type == "VAL":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HG11')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HB12')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HB13')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HB21')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HB22')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HB23')
            if residue_type == "LEU":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD11')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD12')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD13')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD21')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD22')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD23')
            if residue_type == "ILE":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG1', 'CD1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HG12')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG1', 'HG13')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG21')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG22')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG23')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD11')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD12')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD13')
            if residue_type == "MET":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'SD')
                _try_connection_add(self.sidechain_connections, res_atoms, 'SD', 'CE')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE3')
            if residue_type == "PHE":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD1', 'CE1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE1', 'CZ')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE2', 'CZ')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE1', 'HE1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE2', 'HE2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CZ', 'HZ')
            if residue_type == "TYR":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD1', 'CE1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE1', 'CZ')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE2', 'CZ')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CZ', 'OH')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE1', 'HE1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE2', 'HE2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'OH', 'HH')
            if residue_type == "SER":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'OG')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'OG', 'HG')
            if residue_type == "THR":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'OG1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'OG1', 'HG1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG21')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG22')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG2', 'HG23')
            if residue_type == "CYS":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'SG')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'SG', 'HG')
            if residue_type == "ASP":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'OD1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'OD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'OD2', 'HD2')
            if residue_type == "ASN":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'OD1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'ND2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'ND2', 'HD21')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'ND2', 'HD22')
            if residue_type == "GLU":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'OE1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'OE2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'OE2', 'HE2')
            if residue_type == "GLN":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'OE1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'NE2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NE2', 'HE21')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NE2', 'HE22')
            if residue_type == "LYS":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'CE')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE', 'NZ')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE', 'HE3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NZ', 'HZ1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NZ', 'HZ2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NZ', 'HZ3')
            if residue_type == "ARG":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'NE')
                _try_connection_add(self.sidechain_connections, res_atoms, 'NE', 'CZ')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CZ', 'NH1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CZ', 'NH2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NH1', 'HH11')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NH1', 'HH12')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NH2', 'HH21')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NH2', 'HH22')
            if residue_type == "TRP":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD1', 'NE1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'CE3')
                _try_connection_add(self.sidechain_connections, res_atoms, 'NE1', 'CE2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE2', 'CZ2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE3', 'CZ3')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CZ2', 'CH2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CZ3', 'CH2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD1', 'HD1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'NE1', 'HE1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE3', 'HE3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CZ2', 'HZ2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CZ3', 'HZ3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CH2', 'HH2')
            if residue_type == "HIS":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'ND1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'ND1', 'CE1')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD2', 'NE2')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CE1', 'NE2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'ND1', 'HD1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD2', 'HD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE1', 'HE1')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CE2', 'HE2')
            if residue_type == "PRO":
                _try_connection_add(self.sidechain_connections, res_atoms, 'CA', 'CB')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CB', 'CG')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CG', 'CD')
                _try_connection_add(self.sidechain_connections, res_atoms, 'CD', 'N')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CB', 'HB3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CG', 'HG3')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD2')
                _try_connection_add(self.sidechain_h_connections, res_atoms, 'CD', 'HD3')

    def load_peptide_bonds(self):
        """Adds peptide bonds between residues"""
        c_atom = None
        for atom in self.all_atoms:
            if atom.atom_name == 'C':
                c_atom = atom
            elif atom.atom_name == 'N' and c_atom and atom.chain_id == c_atom.chain_id:
                self.mainchain_connections.append((c_atom, atom))


def _try_connection_add(bond_list, atom_dict, key1, key2):
    """Helper function to add a connection only if both atoms are present"""
    if key1 in atom_dict and key2 in atom_dict:
        bond_list.append((atom_dict[key1], atom_dict[key2]))


class Atom:
    """Represents an atom from a PDB file"""

    def __init__(self, line):
        """Reads data from an ATOM or HETATM entry in a PDB file"""
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
        else:
            raise ValueError('Given line is not an atom')
