import torch, os, warnings, io, subprocess, tempfile
from Bio.PDB import PDBIO, Chain, Residue, Polypeptide, Atom, PDBParser
from Bio import pairwise2
import numpy as np
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqUtils import seq1, seq3
from .logging import get_logger
from .protein_residues import normal as RESIDUES

logger = get_logger(__name__)
parser = PDBParser()

def PROCESS_RESIDUES(d):
    d['HIS'] = d['HIP']
    d = {key: val for key, val in d.items() if seq1(key) != 'X'}
    for key in d:
        atoms = d[key]['atoms']
        d[key] = {'CA': 'C'} | {key: val['symbol'] for key, val in atoms.items() if val['symbol'] != 'H' and key != 'CA'}
    return d

RESIDUES = PROCESS_RESIDUES(RESIDUES)
my_dir = tempfile.mkdtemp()

def pdb_to_npy(pdb_path, model_num=0, chain_id=None, seqres=None):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        try:
            model = parser.get_structure('', pdb_path)[model_num]
            chain = model.child_dict[chain_id] if chain_id else model.child_list[0]
        except:
            logger.warning(f'Error opening PDB file {pdb_path}')
            return

    coords, seq = [], ''
    try:
        for resi in list(chain):
            if (resi.id[0] != ' ') or ('CA' not in resi.child_dict) or (resi.resname not in RESIDUES):
                continue
            co = np.zeros((14, 3)) * np.nan
            atoms = RESIDUES[resi.resname]
            seq += resi.resname
            for i, at in enumerate(atoms):
                try: co[i] = resi.child_dict[at].coord
                except: pass
            coords.append(co)
        coords = np.stack(coords)
        seq = seq1(seq)
    except:
        logger.warning(f'Error parsing chain {pdb_path}')
        return

    if not seqres:
        return coords, seq

    coords_ = np.zeros((len(seqres), 14, 3)) * np.nan
    alignment = pairwise2.align.globalxx(seq, seqres)[0]
    if '-' in alignment.seqB:
        logger.warning(f'Alignment gaps {pdb_path}')
        return

    mask = np.array([a == b for a, b in zip(alignment.seqA, alignment.seqB)])
    coords_[mask] = coords
    return coords_, mask

def tmscore(X_path, Y_path, molseq=None, lddt=False, lddt_start=1):
    """Windows-compatible TMscore only. LDDT is skipped."""
    if not isinstance(X_path, str):
        PDBFile(molseq).add(X_path).write(os.path.join(my_dir, 'X.pdb'))
        X_path = os.path.join(my_dir, 'X.pdb')
    if not isinstance(Y_path, str):
        PDBFile(molseq).add(Y_path).write(os.path.join(my_dir, 'Y.pdb'))
        Y_path = os.path.join(my_dir, 'Y.pdb')

    X_path = os.path.abspath(X_path)
    Y_path = os.path.abspath(Y_path)

    try:
        out = subprocess.check_output(['TMscore', '-seq', Y_path, X_path],
                                      stderr=subprocess.DEVNULL,
                                      cwd=my_dir)
    except FileNotFoundError as e:
        logger.error(f"TMscore binary not found: {e}")
        return {'rmsd': -1, 'tm': -1, 'gdt_ts': -1, 'gdt_ha': -1, 'lddt': -1}

    start = out.find(b'RMSD')
    end = out.find(b'rotation')
    out = out[start:end]

    try:
        rmsd, _, tm, _, gdt_ts, gdt_ha, *_ = out.split(b'\n')
        return {
            'rmsd': float(rmsd.split(b'=')[-1]),
            'tm': float(tm.split(b'=')[1].split()[0]),
            'gdt_ts': float(gdt_ts.split(b'=')[1].split()[0]),
            'gdt_ha': float(gdt_ha.split(b'=')[1].split()[0]),
            'lddt': -1  # skipped on Windows
        }
    except Exception as e:
        logger.warning(f"TMscore parse failed: {e}")
        return {'rmsd': -1, 'tm': -1, 'gdt_ts': -1, 'gdt_ha': -1, 'lddt': -1}

class PDBFile:
    def __init__(self, molseq):
        self.molseq = molseq
        self.blocks = []
        self.chain = chain = Chain.Chain('A')
        j = 1
        for i, aa in enumerate(molseq):
            aa = Residue.Residue(id=(' ', i+1, ' '), resname=seq3(aa).upper(), segid='    ')
            for atom in RESIDUES[aa.resname.upper()]:
                at = Atom.Atom(
                    name=atom, coord=None, bfactor=0, occupancy=1, altloc=' ',
                    fullname=f' {atom} ', serial_number=j, element=RESIDUES[aa.resname][atom]
                )
                j += 1
                aa.add(at)
            chain.add(aa)

    def add(self, coords):
        coords = coords.cpu().double().numpy() if isinstance(coords, torch.Tensor) else coords.astype(np.float64)
        for i, resi in enumerate(self.chain):
            resi['CA'].coord = coords[i]

        k = len(self.chain)
        for i, resi in enumerate(self.chain):
            for j, at in enumerate(RESIDUES[resi.resname]):
                if at != 'CA':
                    try:
                        resi[at].coord = coords[k] + resi['CA'].coord
                    except:
                        resi[at].coord = (np.nan, np.nan, np.nan)
                    k += 1

        pdbio = PDBIO()
        pdbio.set_structure(self.chain)
        stringio = io.StringIO()
        pdbio.save(stringio)
        block = stringio.getvalue().split('\n')[:-2]
        self.blocks.append(block)
        return self

    def clear(self):
        self.blocks = []
        return self

    def write(self, path=None, idx=None, reverse=False):
        is_first, str_ = True, ''
        blocks = self.blocks[::-1] if reverse else self.blocks
        for block in ([blocks[idx]] if idx is not None else blocks):
            block = [line for line in block if 'nan' not in line]
            if not is_first:
                block = [line for line in block if 'CONECT' not in line]
            else:
                is_first = False
            str_ += 'MODEL\n' + '\n'.join(block) + '\nENDMDL\n'
        if not path:
            return str_
        with open(path, 'w') as f:
            f.write(str_)

def renumber_pdb(molseq, X_path, X_renum, start=1):
    assert isinstance(molseq, str)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        model = parser.get_structure('', X_path)[0]
    chain = model.child_list[0]
    polypeptide = Polypeptide.Polypeptide(chain.get_residues())
    seq = polypeptide.get_sequence()
    alignment = pairwise2.align.globalxx(seq, molseq)[0]
    assert '-' not in alignment.seqB
    numbering = [i + start for i, c in enumerate(alignment.seqA) if c != '-']
    for n, resi in enumerate(chain):
        resi.id = (' ', 10000 + n, ' ')
    for n, resi in zip(numbering, chain):
        resi.id = (' ', n, ' ')
    pdbio = PDBIO()
    pdbio.set_structure(chain)
    pdbio.save(X_renum)
