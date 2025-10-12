
import sys
import os
from utils import *

# = = = = = = = = = =
props = {
    'T': 1000.,                # Temperature [K]
    'P': 1.013e5,              # Pressure [Pa]
    'phi': 1.0,                # equivalence ratio
    'fuel': 'CH4:1.0',         # fuel components [molar ratio]
    'oxyd': 'O2:1.0, N2:3.76', # oxydizer components [molar ratio]
    'width': 0.05,             # domain width for calculating flame [m]
    'mixing': None
}
# for CH4 non-premixed (case Sandia-D)
# mech_file = "src/Skeletal29_N/Skeletal29_N.cti"  # mech file, in .cti or .yaml format
mech_file = "mechanism/gri12.yaml"

props["type"] = "nonpremixed_counterflow"  #
props["path"] = (
    "./data/case_CH4_counterflow_PP/"  # data path to save 1D flame data and figures
)
props["fuel"] = "CH4:0.25, N2:0.75"  # fuel composition [molar ratio]
props["oxid"] = "O2:1.0, N2:3.762"  # oxidizer composition [molar ratio]
props["P"] = 1.0 * ct.one_atm  # pressure [Pa]
props["Tfuel"] = 294  # fuel temperature [K]
props["Toxid"] = 291  # oxidizer temperature [K]
props["Text"] = 300  # extinction temperature [K]
props["width"] = 0.05  # width [m]
props["element_num"] = 5
props['inert_specie'] = 'AR'

# TODO: 
# 1. 把所有的需要的参数一次储存, 增加一个储存的csv文件(cantera_input and cantera_output)
# 2. 画图，画出不同拉伸率下的温度分布，CEM值分布以及其他的参数

def get_partially_premixd_counterflow(gas, props, left_phi, right_phi, left_U, right_U_list):
    f = ct.CounterflowDiffusionFlame(gas, width=props["width"])
    f.P = props["P"]  # 1 bar
    
    gas.set_equivalence_ratio(left_phi, props["fuel"], props["oxid"])
    f.fuel_inlet.mdot = gas.density * left_U  # kg/m^2/s
    f.fuel_inlet.X= gas.X
    f.fuel_inlet.T = 300  # K
    
    gas.set_equivalence_ratio(right_phi, props["fuel"], props["oxid"])
    f.oxidizer_inlet.mdot = gas.density * right_U_list[0]  # kg/m^2/s
    f.oxidizer_inlet.T = 300  # K
    f.oxidizer_inlet.X = gas.X
    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
    f.set_initial_guess()
    temperature_limit_extinction = 310 # [K]
    def interrupt_extinction(t):
        if np.max(f.T) < temperature_limit_extinction:
            raise FlameExtinguished('Flame extinguished')
        return 0.
    f.set_interrupt(interrupt_extinction)
    print(f'Running simulation for left_phi={left_phi}, right_phi={right_phi}, left_U={left_U}, right_U={right_U_list[0]}')
    try:
        f.solve(loglevel=1, auto=True)
        output_dir = os.path.join(props["path"], "run_partially_premixed")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        file_name = os.path.join(
            output_dir,
            f"left_phi_{left_phi}_right_phi_{right_phi}_left_U_{left_U:.4f}_right_U_{right_U_list[0]:.4f}.csv"
        )
        f.write_csv(file_name, quiet=False)
        mass_fractions, MF = get_mass_fractions_from_1D_flame(f, gas, props["fuel"], props["oxid"])
        CEM_values = get_cem_from_1D_flame(f, gas, props["P"], element_num=props["element_num"])
        df = pd.read_csv(file_name)
        mole_fraction_columns = [col for col in df.columns if col.startswith('X_')]
        df.drop(columns=mole_fraction_columns, inplace=True)
        species_names = gas.species_names
        for j, species in enumerate(species_names):
            df[species] = [mf[j] for mf in mass_fractions]
        df['CEM'] = CEM_values
        df['MF'] = MF
        normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
        df['normalized_grid'] = normalized_grid
        df.to_csv(file_name, index=False)
    except FlameExtinguished:
        print(f'Flame extinguished for left_phi={left_phi}, right_phi={right_phi}, left_U={left_U:.4f}')
    except ct.CanteraError as e:
        print(f'Error during simulation for left_phi={left_phi}, right_phi={right_phi}, left_U={left_U:.4f}', e)

    try:
        exp_d_a = -1. / 2.
        for i, right_U in enumerate(right_U_list[1:]):
            strain_factor = (left_U + right_U)/(left_U + right_U_list[i])
            while np.max(f.T) > temperature_limit_extinction:
                f.flame.grid *= strain_factor ** exp_d_a
                f.oxidizer_inlet.mdot = gas.density * right_U  # kg/m^2/ss
                f.solve(loglevel=1, auto=True)
                output_dir = os.path.join(props["path"], "run_partially_premixed")
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                file_name = os.path.join(
                    output_dir,
                    f"left_phi_{left_phi}_right_phi_{right_phi}_left_U_{left_U:.4f}_right_U_{right_U:.4f}.csv"
                )
                f.write_csv(file_name, quiet=False)
                mass_fractions, MF = get_mass_fractions_from_1D_flame(f, gas, props["fuel"], props["oxid"])
                CEM_values = get_cem_from_1D_flame(f, gas, props["premixed_P"], element_num=props["element_num"])
                df = pd.read_csv(file_name)
                mole_fraction_columns = [col for col in df.columns if col.startswith('X_')]
                df.drop(columns=mole_fraction_columns, inplace=True)
                species_names = gas.species_names
                for j, species in enumerate(species_names):
                    df[species] = [mf[j] for mf in mass_fractions]
                df['CEM'] = CEM_values
                df['MF'] = MF
                normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
                df['normalized_grid'] = normalized_grid
                df.to_csv(file_name, index=False)
    except FlameExtinguished:
        print(f'Flame extinguished for left_phi={left_phi}, right_phi={right_phi}, left_U={left_U:.4f}')
    except ct.CanteraError as e:
        print(f'Error during simulation for left_phi={left_phi}, right_phi={right_phi}, left_U={left_U:.4f}', e)



# = = = = = = = = = =
if __name__ == "__main__":
    # prepare paths
    figs_path = os.path.join(props["path"], "figs/")
    check_path(props["path"])
    check_path(figs_path)

    # prepare cases, only one case is needed for counterflow
    gas = ct.Solution(mech_file)
    last_spec = props['inert_specie']
    specs = gas.species()[:]
    N2_ind = gas.species_index(last_spec)
    gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
            species=specs[:N2_ind] + specs[N2_ind + 1:] + [specs[N2_ind]],
            reactions=gas.reactions())

    gas.set_equivalence_ratio(1.0, props["fuel"], props["oxid"])
    Zst = gas.mixture_fraction(props["fuel"], props["oxid"])
    print("Zst = ", Zst)
    # if you want simulate the 1D flames, turn on this comment
    left_phi=1.4
    right_phi=0.7
    left_U=0.39
    right_U_list=np.linspace(0.25,0.35,5) 
    print("right_U_list=", right_U_list)
    get_partially_premixd_counterflow(gas, props, left_phi, right_phi, left_U, right_U_list)

    # terminate useless threads
    if rank != 0:
        sys.exit(0)