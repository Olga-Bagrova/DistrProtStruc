import numpy as np
import pandas as pd

from dataclasses import dataclass

import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb


class Domain():
    """
    A class for whole domain data and processing
    
    Attributes
    ----------
    PDBid : str
        identification number of PDB-structure
    length_dom : int
        length of the domain in amino acids
    chain_dom : str
        chain of the domain in PDB structure
    start_dom : int
        first amino acid in domain
    stop_dom : int
        last amino acid in domain
    struc_distr_df : pd.DataFrame
        data about specific type of secondary structure
    pdbx_file : biotite.structure.io.pdbx.PDBxFile
        pdbx-file of the PDB-structure

    Methods
    -------
    get_suitable_structures():
        filter structures for the domain
    calculate_distr(column_name: str):
        calculate the distribution of a specific type of secondary structure in the domain
    """

    def __init__(self, PDBid: str, length_dom: int, chain_dom: str, start_dom: int, stop_dom: int):
        """Constructs all the necessary attributes for the object."""
        
        self.PDBid = PDBid
        self.length_dom = length_dom
        self.chain_dom = chain_dom
        self.start_dom = start_dom
        self.stop_dom = stop_dom
        self.struc_distr_df = pd.DataFrame({'PDBid': [self.PDBid], 'chain': [self.chain_dom], 
                                            'start': [start_dom], 'end': [self.stop_dom], 
                                            'length': [self.length_dom]})
        self.pdbx_file = pdbx.PDBxFile.read(rcsb.fetch(self.PDBid, 'pdbx'))

    def get_suitable_structures(self):
        """Filter structures for the domain from PDB-structure."""
        
        self.structures_df = self.structures_df[self.structures_df.chain_str == self.chain_dom]
        self.structures_df = self.structures_df[self.structures_df.start_str >= self.start_dom]
        self.structures_df = self.structures_df[self.structures_df.stop_str <= self.stop_dom]
        self.structures_df.index = np.arange(len(self.structures_df))

    def calculate_distr(self, column_name: str):
        """Calculate the distribution of a specific type of secondary structure in the domain."""
        
        if self.structures_df.empty:
            distr = [0] * 100
            distr_str = (';'.join(map(str, distr)))
            self.struc_distr_df.loc[len(self.struc_distr_df) - 1, column_name] = distr_str
        else:
            distr = []
            distr = calculate_distribution(self.structures_df, self.length_dom, self.start_dom, self)
            distr_str = (';'.join(map(str, distr)))
            self.struc_distr_df.loc[len(self.struc_distr_df) - 1, column_name] = distr_str
        
        return self.struc_distr_df


@dataclass
class SecondaryStructure():
    """
    Dataclass for secondary structure data
    
    Attributes
    ----------
    id_str : list
        identification numbers of specific secondary structure for the domain
    start_str : list
        first amino acids of the secondary structure for the domain
    stop_str : list
        last amino acids of the secondary structure for the domain
    class_str : list
        classes of the secondary structure for the domain
    length_str : list
        lengths of the secondary structure for the domain
    chain_str : list
        chain identificators of the secondary structure for the domain
    """
    
    id_str: list
    start_str: list
    stop_str: list
    class_str: list
    length_str: list
    chain_str: list


class Helix(Domain):
    """
    Child class of Domain class for all types of helices
    
    Attributes
    ----------
    secondary_structure_obj: SecondaryStructure
        object of SecondaryStructure
    structures_df: pd.DataFrame
        dataframe for all helices data

    Methods
    -------
    get_suitable_helices(class_str: str = '1'):
        filter helices structures (class_str the same as in pdbx annotation)
    """
    
    def __init__(self, PDBid: str, length_dom: int, chain_dom: str, start_dom: int, stop_dom: int):
        """Constructs all the necessary attributes for the object."""
        
        super().__init__(PDBid=PDBid, length_dom=length_dom, chain_dom=chain_dom, start_dom=start_dom, stop_dom=stop_dom)
        self.secondary_structure_obj = SecondaryStructure(id_str=self.pdbx_file [self.PDBid, 'struct_conf']['pdbx_PDB_helix_id'],
                                                          start_str = self.pdbx_file[self.PDBid, 'struct_conf']['beg_label_seq_id'],
                                                          stop_str = self.pdbx_file [self.PDBid, 'struct_conf']['end_label_seq_id'],
                                                          class_str = self.pdbx_file [self.PDBid, 'struct_conf']['pdbx_PDB_helix_class'],
                                                          length_str = self.pdbx_file [self.PDBid, 'struct_conf']['pdbx_PDB_helix_length'],
                                                          chain_str = self.pdbx_file [self.PDBid, 'struct_conf']['beg_auth_asym_id'])
        self.secondary_structure_obj.start_str = np.array([int(self.secondary_structure_obj.start_str)] if isinstance(self.secondary_structure_obj.start_str, str) 
                                                          else [int(t) for t in self.secondary_structure_obj.start_str])
        self.secondary_structure_obj.stop_str = np.array([int(self.secondary_structure_obj.stop_str)] if isinstance(self.secondary_structure_obj.stop_str, str) 
                                                         else [int(t) for t in self.secondary_structure_obj.stop_str])   
        self.structures_df = pd.DataFrame({'id_str' : self.secondary_structure_obj.id_str,
                                           'start_str' : self.secondary_structure_obj.start_str, 
                                           'stop_str': self.secondary_structure_obj.stop_str, 
                                           'class_str': self.secondary_structure_obj.class_str,
                                           'length_str': self.secondary_structure_obj.length_str, 
                                           'chain_str': self.secondary_structure_obj.chain_str})

    def get_suitable_helices(self, class_str: str = '1'):
        """Filter helices structures."""
        
        self.get_suitable_structures()
        self.structures_df = self.structures_df[self.structures_df.class_str == str(class_str)]


class Sheet(Domain):
    """
    Child class of Domain class for all sheets
    
    Attributes
    ----------
    secondary_structure_obj: SecondaryStructure
        object of SecondaryStructure
    structures_df: pd.DataFrame
        dataframe for all sheets data

    Methods
    -------
    get_suitable_sheets(class_str: str = '1'):
        filter sheets structures
    """
    
    def __init__(self,PDBid: str, length_dom: int, chain_dom: str, start_dom: int, stop_dom: int):
        """Constructs all the necessary attributes for the object."""
        
        super().__init__(PDBid=PDBid, length_dom=length_dom, chain_dom=chain_dom, start_dom=start_dom, stop_dom=stop_dom)

        self.secondary_structure_obj = SecondaryStructure(id_str=self.pdbx_file[self.PDBid, 'struct_sheet_range']['sheet_id'],
                                                          start_str = self.pdbx_file[self.PDBid, 'struct_sheet_range']['beg_label_seq_id'],
                                                          stop_str = self.pdbx_file [self.PDBid, 'struct_sheet_range']['end_label_seq_id'],
                                                          class_str = None,
                                                          length_str = None,
                                                          chain_str = self.pdbx_file [self.PDBid, 'struct_sheet_range']['beg_auth_asym_id'])
        self.secondary_structure_obj.start_str = np.array([int(self.secondary_structure_obj.start_str)] if isinstance(self.secondary_structure_obj.start_str, str) 
                                                          else [int(t) for t in self.secondary_structure_obj.start_str])
        self.secondary_structure_obj.stop_str = np.array([int(self.secondary_structure_obj.stop_str)] if isinstance(self.secondary_structure_obj.stop_str, str) 
                                                         else [int(t) for t in self.secondary_structure_obj.stop_str])   
        self.structures_df = pd.DataFrame({'id_str': self.secondary_structure_obj.id_str,
                                           'start_str' : self.secondary_structure_obj.start_str, 
                                           'stop_str': self.secondary_structure_obj.stop_str, 
                                           'chain_str': self.secondary_structure_obj.chain_str})

    def get_suitable_sheets(self):
        """Filter suitable for domain sheets."""
        
        self.get_suitable_structures()


class AlphaHelix(Helix):
    """
    Child class of Helix class for alpha-helices

    Methods
    -------
    get_alpha_helices():
        filter suitable for domain alpha-helices
    """
    
    def __init__(self, PDBid: str, length_dom: int, chain_dom: str, start_dom: int, stop_dom: int):
        """Constructs all the necessary attributes for the object."""
        
        super().__init__(PDBid=PDBid, length_dom=length_dom, chain_dom=chain_dom, start_dom=start_dom, stop_dom=stop_dom)
        
    def get_alpha_helices(self):
        """Filter suitable for domain alpha-helices."""
        
        self.get_suitable_helices('1')        


class ThreeTenHelix(Helix):
    """
    Child class of Helix class for 3_10-helices

    Methods
    -------
    get_threeten_helices():
        filter suitable for domain 3_10-helices
    """
    
    def __init__(self, PDBid: str, length_dom: int, chain_dom: str, start_dom: int, stop_dom: int):
        """Constructs all the necessary attributes for the object."""
        
        super().__init__(PDBid=PDBid, length_dom=length_dom, chain_dom=chain_dom, start_dom=start_dom, stop_dom=stop_dom)
        
    def get_threeten_helices(self):
        """Filter suitable for domain 3_10-helices."""
        
        self.get_suitable_helices('5')


class BetaSheet(Sheet):
    """
    Child class of Helix class for alpha-helices
    """
    
    def __init__(self, PDBid: str, length_dom: int, chain_dom: str, start_dom: int, stop_dom: int):
        """Constructs all the necessary attributes for the object."""
        
        super().__init__(PDBid=PDBid, length_dom=length_dom, chain_dom=chain_dom, start_dom=start_dom, stop_dom=stop_dom)

 
class IrregularStructure(Domain):
    """
    Child class of Domain class for irregular structure
    
    Attributes
    ----------
    ah_df_res: pd.DataFrame
        alpha-helix distribution in the domain dataframe
    tth_df_res: pd.DataFrame
        3_10-helix distribution in the domain dataframe
    bs_df_res: pd.DataFrame
        bta-sheet distribution in the domain dataframe
    total_df_res: pd.DataFrame
        all typesof secondary structure distribution in the domain dataframe

    Methods
    -------
    calculate_irregular_distr():
        calculate the distribution of irregular structure in the domain
    """

    def __init__(self, ah_df_res, tth_df_res, bs_df_res):
        """Constructs all the necessary attributes for the object."""
        
        self.ah_df_res = ah_df_res
        self.tth_df_res = tth_df_res
        self.bs_df_res = bs_df_res
        self.total_df_res = ah_df_res[['PDBid', 'chain', 'start', 'end', 'length']].copy()

    def calculate_irregular_distr(self):
        """Calculate the distribution of irregular structure in the domain."""
        
        for i in range(len(self.ah_df_res)):
            list_alpha = self.ah_df_res.iloc[i, -1].split(';')
            list_beta = self.bs_df_res.iloc[i, -1].split(';')
            list_threeten = self.tth_df_res.iloc[i, -1].split(';')
            
            distr = [0] * 100 
            for index in range(100):
                distr[index] = 1 - float(list_alpha[index]) - float(list_beta[index]) - float(list_threeten[index])
                distr[index] = round(distr[index], 2)
                distr_str = (';'.join(map(str,distr)))
                self.total_df_res.loc[i, 'Irregular_structure_distr'] = distr_str 
                
        return self.total_df_res


def read_prot_set(path: str, mode: str = 'xlsx'):
    """
    Read protein set data from file through. 
    The following data-file structure is required: PDBid, chain, start, end. Header is needed.
    arguments:
        - path (str): path to the file for reading
        - mode (str): file format to read. Only xlsx, csv and tsv formats are allowed
    return:
        - prot_set_df (pandas.DataFrame): dataframe of the protein set
    """
    
    match mode:
        case 'xlsx':
            prot_set_df = pd.read_excel(path)
        case 'csv':
            prot_set_df = pd.read_csv(path)
        case 'tsv':
            prot_set_df = pd.read_csv(path, sep='\t')

    prot_set_df['length'] = prot_set_df['end'] - prot_set_df['start'] + 1

    return prot_set_df


def calculate_distribution(structures_df, length_dom: int, start_dom: int, structure_obj) -> list:
    """
    Calculates the distribution of a specific type of secondary structure.
    arguments:
        - structures_df (pd.DataFrame): dataframe with suitable secondary structure
        - length_dom (int): length of the domain chain
        - start_dom (int): first amino acid in polypeptide chain
        - structure_obj (Domain): object of the Domain class
    return:
        - struc_distr_scaled (list): secondary structure distribution (scaled)
    """
     
    struc_distr = [0] * length_dom
    
    for i in range(0, len(structures_df)):
        start = int(structures_df.iloc[i]['start_str']) - structure_obj.start_dom
        stop = int(structures_df.iloc[i]['stop_str']) - structure_obj.start_dom
        
        for aa in range (start, stop + 1):
            struc_distr[aa] = 1
    
    struc_distr_scaled = scale(struc_distr, length_dom)
    
    return struc_distr_scaled


def scale(struc_distr: list, length_dom: int) -> list:
    """
    Scale the distribution of a specific type of secondary structure to 100 units.
    arguments:
        - struc_distr (list): domain-long distribution of the secondary structure
        - length_dom (int): length of the domain
    return:
        - scaled_distr: scaled distribution of a specific type of secondary structure
    """
    
    scaled_distr = [0] * 100
    weight = length_dom / 100
    
    if weight != 1:
        k, t = 1, weight
        i,j  = 1, 1
        while i < length_dom and j < 100:
            a = min(k, t)
            scaled_distr[j] = round(scaled_distr[j] + a * struc_distr[i], 2)
            k, t = round(k - a, 2), round(t - a, 2)
            if k == 0:
                i += 1
                k = 1
            if t == 0:
                j = j + 1
                t = weight
        scaled_distr = [round(x / weight, 2) for x in scaled_distr]
    else:
        scaled_distr = struc_distr
    
    return scaled_distr


def process_alpha_helix(df_set, res_file_name: str = 'Result_ah.xlsx', _column_name: str = 'Alpha_helix_distr'):
    """
    Get alpha-helix distribution for all proteins in the set. 
    arguments:
        - df_set (pd.DataFrame): dataframe of the protein set
        - res_file_name (str): name for file with results 
        - _column_name (str): name of the column in results
    return:
        - ah_df_res (pandas.DataFrame): dataframe of the alpha-helix distribution
    """
    
    for i, row in df_set.iterrows():
        ah = AlphaHelix(PDBid=row['PDBid'], length_dom=row['length'], chain_dom=row['chain'], start_dom=row['start'], stop_dom=row['end'])
        ah.get_alpha_helices()
        ah.calculate_distr(column_name=_column_name)
        ah_df_res = ah.struc_distr_df if i == 0 else pd.concat([ah_df_res, ah.struc_distr_df])
        
    ah_df_res = ah_df_res.reset_index(drop=True)
    ah_df_res.to_excel(res_file_name)
    
    return ah_df_res


def process_threeten_helix(df_set, res_file_name: str = 'Result_tth.xlsx', _column_name: str = 'Three_ten_helix_distr'):
    """
    Get 3_10-helix distribution for all proteins in the set. 
    arguments:
        - df_set (pd.DataFrame): dataframe of the protein set
        - res_file_name (str): name for file with results 
        - _column_name (str): name of the column in results
    return:
        - tth_df_res (pandas.DataFrame): dataframe of the 3_10-helix distribution
    """
    
    for i, row in df_set.iterrows():
        tth = ThreeTenHelix(PDBid=row['PDBid'], length_dom=row['length'], chain_dom=row['chain'], start_dom=row['start'], stop_dom=row['end'])
        tth.get_threeten_helices()
        tth.calculate_distr(column_name=_column_name)
        tth_df_res = tth.struc_distr_df if i == 0 else pd.concat([tth_df_res, tth.struc_distr_df])
        
    tth_df_res = tth_df_res.reset_index(drop=True)
    tth_df_res.to_excel(res_file_name)
    
    return tth_df_res


def process_beta_sheets(df_set, res_file_name: str = 'Result_bs.xlsx', _column_name: str = 'Beta_sheet_distr'):
    """
    Get beta-sheet distribution for all proteins in the set. 
    arguments:
        - df_set (pd.DataFrame): dataframe of the protein set
        - res_file_name (str): name for file with results 
        - _column_name (str): name of the column in results
    return:
        - bs_df_res (pandas.DataFrame): dataframe of the beta-sheet distribution
    """
    
    for i, row in df_set.iterrows():
        bs = BetaSheet(PDBid=row['PDBid'], length_dom=row['length'], chain_dom=row['chain'], start_dom=row['start'], stop_dom=row['end'])
        bs.get_suitable_sheets()
        bs.calculate_distr(column_name=_column_name)
        bs_df_res = bs.struc_distr_df if i == 0 else pd.concat([bs_df_res, bs.struc_distr_df])
        
    bs_df_res = bs_df_res.reset_index(drop=True)
    
    bs_df_res.to_excel(res_file_name)

    return bs_df_res


def process_irregular_struc(ah_df_res, tth_df_res, bs_df_res, res_file_name: str = 'Result_iss.xlsx'):
    """
    Get alpha-helix distribution for all proteins in the set. 
    arguments:
        - ah_df_res (pandas.DataFrame): dataframe of the alpha-helix distribution
        - tth_df_res (pandas.DataFrame): dataframe of the 3_10-helix distribution
        - bs_df_res (pandas.DataFrame): dataframe of the beta-sheet distribution
        - res_file_name (str): name for file with results 
    return:
        - iss_df_res (pandas.DataFrame): dataframe of the irregular structure distribution
    """
    
    iss = IrregularStructure(ah_df_res, tth_df_res, bs_df_res)
    iss_df_res = iss.calculate_irregular_distr()
    iss_df_res.to_excel(res_file_name)
    
    return iss_df_res
