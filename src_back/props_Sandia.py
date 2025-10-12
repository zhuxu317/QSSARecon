import cantera as ct
import os

props_list = []

props = {
    'file_name': "Sandia_F.csv",
    'T': 1000.,
    'P': 0.993 * ct.one_atm,
    'fuel': 'CH4',
    'oxidizer': 'O2:1.0, N2:3.762',
    'strain_factor': 1.2,
    'width': 0.05,
    'dt': 1e-2,
    'time_end': 1e-1,
    'perturb': 0,
    'mixing': None,
    'mech_file': "/data/ZhuXu/Cantera/CEMA/mechanism/gri12/gri12.cti",
    'mech_path': '/data/ZhuXu/Cantera/CEMA/mechanism/gr12',
    'type': "Sandia_flame",
    'element_num': 5,
    'inert_specie': 'AR'
}

props_list.append(props)
