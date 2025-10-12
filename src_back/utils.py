import os, sys, time, json
import argparse
import threading
import logging
import numpy as np
import pandas as pd
import cantera as ct
from copy import deepcopy
from scipy import interpolate
from scipy.stats import beta
from scipy.special import erfcinv, inv_boxcox
from itertools import accumulate
# import torch
import warnings
import re
import shutil
import csv
import math
# import torch.nn as nn
# import torch.nn.functional as F
from scipy import stats 
# = = = = = = = = = =
# Matplotlib settings
try:
    import matplotlib
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt
    matplotlib.rc('font', **{'size':13})
    matplotlib.rc('text', usetex=False)
    matplotlib.rc('legend', **{'frameon':False} )
except:
    print("NOTE: matplotlib failed, are you on a server with no GUI?")

# = = = = = = = = = =
# MPI settings
try:
    # if MPI is supported
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
except:
    # if MPI is not supported
    print("NOTE: mpi4py failed, you may need to check the installation of mpi4py")
    comm = None
    rank = 0
    size = 1

# = = = = = = = = = =
# global settings
t0 = time.time()
threadLock = threading.Lock()
np.random.seed(0x7777777>>7)
line_arr = ('-','--','-.',':')
color_arr = ('k','r','b','m','y','g', 'c')
symbol_arr = ('s','o','v','^','*')
os.environ["MKL_NUM_THREADS"] = '4'
os.environ["NUMEXPR_NUM_THREADS"] = '4'
os.environ["OMP_NUM_THREADS"] = '4'

# for graphical / saving figures
def c2i(*tupl):
    """ Convert centermeter to inch
    """
    return tuple(i/2.54 for i in tupl)

def check_path(directory_path):
    """
    Checks if a directory exists, and if not, creates it.

    Parameters:
    - directory_path: str, the path of the directory to check/create.
    """
    if not os.path.exists(directory_path):
        try:
            os.makedirs(directory_path)
            print(f"Directory {directory_path} created.")
        except OSError as e:
            print(f"Error creating directory {directory_path}: {e}")
    else:
        print(f"Directory {directory_path} already exists.")


def checkexists(filename, delete=False):
    """ Check if file exists, if yes, delete it by default input
    """
    if os.path.exists(filename):
        if delete:
            try:
                os.remove(filename)
            except:
                pass
        return True
    return False

# = = = = = = = = = =
# find conditional rank
def first(iterable, condition = lambda x: True):
    """ Find first value satisfying the condition
    """
    for i,li in enumerate(iterable):
        if condition(li):
            return i, li
    return None, None

def accum(iterable, condition = lambda x: x**2, sum_limit = 0.99):
    """ Find fist position satisfying the cumsum condition
    """
    adder = 0.
    for i,li in enumerate(iterable):
        adder += condition(li)
        if adder > sum_limit:
            return i

# = = = = = = = = = =
# algebra
def normalize(data):
    """ Normalize a vector
    """
    return data/np.linalg.norm(data)

def cosine(a, b):
    """ Calculate consine of two vector
    """
    return np.dot(normalize(a), normalize(b))

def cdiff(data):
    """ Central differential of data
    """
    if len(data)==1: return np.array([0])
    d = np.diff(data)
    return np.array([d[0]] + list((d[1:]+d[:-1])/2) + [d[-1]])

def diffMax(x, y):
    """ Find max position of dy/dx
    """
    dydx = np.abs(cdiff(y)/cdiff(x))
    pos = np.argmax(dydx)
    return pos, dydx[pos]

def curvMax(x, y, n=5, N=100):
    """ Find max position of dy/dx with high resolution
    """
    l = len(x)-1
    pos, x_val = diffMax(x, y)
    lpos = 0 if pos-n<0 else pos-n
    rpos = l if pos+n>l else pos+n
    x, y = x[lpos:rpos], y[lpos:rpos]
    x_new = np.linspace(x[0], x[-1], N)
    f = interpolate.interp1d(x, y, 'quadratic')
    return diffMax(x_new,f(x_new))

def check_cuda_and_print_info():
    if torch.cuda.is_available():
        print("CUDA is available. Using GPU.")
        print(f"Device name: {torch.cuda.get_device_name(0)}")
        print(f"Memory Allocated: {torch.cuda.memory_allocated(0)} bytes")
        print(f"Memory Cached: {torch.cuda.memory_reserved(0)} bytes")
    else:
        warnings.warn("CUDA is not available. Please check your CUDA installation and GPU settings.")
        print("Using CPU.")


def save_checkpoint(state, filename='checkpoint.pth.tar'):
    """
    保存模型和训练状态。
    """
    torch.save(state, filename)

def find_checkpoint(model_type, output_dir):
    """
    Finds the latest checkpoint based on the model type.
    """
    if model_type == 'MLP':
        checkpoint_path = find_latest_checkpoint(output_dir)
        return checkpoint_path, None
    elif model_type == 'KANs':
        return KANs_find_latest_checkpoint(output_dir)
    
    
def load_checkpoint(model, model_type, checkpoint_path, model_dir, checkpoint_epoch, optimizer, device):
    """
    Loads the checkpoint based on the model type.
    """
    # Correct logging format using logging's built-in string formatting
    logging.info("=> Loading checkpoint '%s'", checkpoint_path)
    
    if model_type == 'MLP':
        checkpoint = torch.load(checkpoint_path)
        model.load_state_dict(checkpoint['state_dict'])
        model.to(device)  # Move model to the correct device before loading optimizer state
        optimizer.load_state_dict(checkpoint['optimizer'])
        
        # Move optimizer state to GPU if needed
        for state in optimizer.state.values():  
            for k, v in state.items():
                if isinstance(v, torch.Tensor):
                    state[k] = v.to(device)
        
        start_epoch = checkpoint['epoch']
        train_losses = checkpoint['train_losses']
        val_losses = checkpoint['val_losses']
        
    elif model_type == 'KANs':
        logging.info("checkpoint_path=%s", checkpoint_path)
        folder_dir, name = os.path.split(checkpoint_path)
        model.load_state_dict(torch.load(folder_dir + '/' + name))
        start_epoch = checkpoint_epoch
        
        train_loss_csv = folder_dir + '/training_validation_losses.csv'
        copy_file(train_loss_csv, model_dir, False)
        train_losses = []  
        val_losses = []  
        
        if os.path.exists(train_loss_csv):
            with open(train_loss_csv, 'r') as file:
                reader = csv.reader(file)
                for row in reader:
                    train_losses.append(float(row[1]))
                    val_losses.append(float(row[2]))
        
        logging.info("start_epoch=%s", start_epoch)
    
    return model, start_epoch, train_losses, val_losses



# 定义一个函数来更新checkpoint列表和删除旧的checkpoint
def update_checkpoint_list(new_checkpoint_path, current_checkpoints):
    # 设定最大保存的checkpoint数量
    max_checkpoints = 4
    current_checkpoints.append(new_checkpoint_path)
    # 如果当前保存的checkpoint数量超过了限制，则删除最旧的文件
    while len(current_checkpoints) > max_checkpoints:
        oldest_checkpoint = current_checkpoints.pop(0)  # 移除并返回列表中的第一个元素
        if os.path.exists(oldest_checkpoint):
            os.remove(oldest_checkpoint)
            print(f"Removed old checkpoint: {oldest_checkpoint}")
            


def find_latest_checkpoint(output_folder_name):
    # 使用正则表达式匹配文件名中的epoch数字
    pattern = re.compile(r'checkpoint_epoch_(\d+)\.pth\.tar$')

    # 获取所有匹配的checkpoint文件及其对应的epoch
    checkpoints = []
    for filename in os.listdir(output_folder_name):
        match = pattern.match(filename)
        if match:
            epoch = int(match.group(1))
            checkpoints.append((epoch, filename))

    # 如果没有找到任何checkpoint文件，返回None
    if not checkpoints:
        return None

    # 找到epoch数最大的那个文件
    _, latest_checkpoint_filename = max(checkpoints, key=lambda x: x[0])
    # 构建最新checkpoint文件的完整路径
    latest_checkpoint_path = os.path.join(output_folder_name, latest_checkpoint_filename)
    return latest_checkpoint_path


def KANs_find_latest_checkpoint(output_folder_name):
    import re  # Add this line
    pattern = re.compile(r'KANs_ckpt_(\d+)_state$')

    # Get all matching checkpoint files and their corresponding steps
    checkpoints = []
    for filename in os.listdir(output_folder_name):
        match = pattern.match(filename)
        if match:
            step = int(match.group(1))
            checkpoints.append((step, filename))

    # If no checkpoint files are found, return None and 0
    if not checkpoints:
        return None, 0

    # Find the file with the largest step number
    max_step, latest_checkpoint_filename = max(checkpoints, key=lambda x: x[0])
    # Construct the full path of the latest checkpoint file
    latest_checkpoint_path = os.path.join(output_folder_name, latest_checkpoint_filename)
    return latest_checkpoint_path, max_step


def copy_file(src, dest_folder, overwrite=False):
    """
    Copies a file from src to dest_folder, optionally overwriting the existing file.

    Parameters:
    - src (str): The path of the source file.
    - dest_folder (str): The path of the destination folder.
    - overwrite (bool): If True, the destination file will be overwritten if it exists.

    Returns:
    - None
    """
    # Ensure the destination folder exists
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
        print(f"Destination folder {dest_folder} created.")

    # Extract the file name from the source path
    file_name = os.path.basename(src)
    
    # Construct the full destination path
    dest = os.path.join(dest_folder, file_name)
    
    # Check if the destination file exists
    if os.path.exists(dest):
        if overwrite:
            try:
                os.remove(dest)
                print(f"Existing file {dest} removed.")
            except Exception as e:
                print(f"Failed to remove existing file {dest}: {e}")
                return
        else:
            print(f"File {dest} already exists. Set overwrite=True to overwrite it.")
            return
    
    # Attempt to copy the file
    try:
        shutil.copy(src, dest)
        print(f"File successfully copied from {src} to {dest}.")
    except Exception as e:
        print(f"An error occurred while copying from {src} to {dest}: {e}")

def boxcox_transform_and_visualize(specie_names, data, start_index, small_value=1e-15):
    """
    对数据执行Box-Cox变换，并进行可视化对比。
    
    :param data: 待变换的数据矩阵
    :param start_index: 开始变换的索引
    :param small_value: 添加的小值以确保所有值都是正值
    """
    transformed_data = np.copy(data)
    lambdas = {}
    for i in range(start_index, data.shape[0]):  # 根据start_index决定是否跳过某些变量
        # 添加一个很小的常数以确保所有值都是正值
        data_with_small_value = data[i] + small_value
        # 执行Box-Cox变换
        transformed, lambda_val = stats.boxcox(data_with_small_value)
        transformed_data[i] = transformed  # 更新变换后的值
        lambdas[specie_names[i-start_index]] = lambda_val
        # 可视化变换前后的分布对比
        from plotSampleCantera import visualize_distribution_comparison
        visualize_distribution_comparison(data[i], transformed, specie_names[i-start_index], start_index)
    return transformed_data, lambdas


def reverse_boxcox_transform(transformed_data, lambdas):
    """
    Reverse the Box-Cox transformation on the data.

    Parameters:
    - transformed_data: np.ndarray, the Box-Cox transformed data.
    - lambdas: dict, the lambda values used for Box-Cox transformation for each feature.

    Returns:
    - original_data: np.ndarray, the data in its original scale.
    """
    original_data = np.copy(transformed_data)
    for i, (feature_name, lambda_val) in enumerate(lambdas.items()):
        if lambda_val is not None:  # Check if lambda_val is not None
            original_data[:,i] = inv_boxcox(transformed_data[:, i], lambda_val)
    return original_data

def load_lambdas(lambda_path):
    """
    从指定路径读取Box-Cox变换的λ值，并将其存储在一个字典中。
    """
    lambda_df = pd.read_csv(lambda_path)
    lambda_dict = lambda_df.set_index('Feature')['Lambda'].to_dict()
    print("lambda_dict=",lambda_dict)
    return lambda_dict


def reverse_transform_output(transformed_data, scaler, lambdas=None, small_value=1e-15, bool_temperature=True, bool_output_boxcox=True):
    """
    Reverses the StandardScaler and Box-Cox transformations on the output data.

    Parameters:
    - transformed_data: np.ndarray, the transformed data.
    - scaler: sklearn.preprocessing.StandardScaler, the scaler used for output data.
    - lambdas: dict, the lambda values used for Box-Cox transformation for each feature.

    Returns:
    - original_data: np.ndarray, the data in its original scale.
    """
    # Ensure transformed_data is 2D
    transformed_data = transformed_data.reshape(1, -1) if transformed_data.ndim == 1 else transformed_data
    unscaled_data = scaler.inverse_transform(transformed_data)
    if not bool_output_boxcox:
        if unscaled_data.shape[0] == 1:
            unscaled_data = unscaled_data.flatten()
        return unscaled_data
    # Inverse Box-Cox transformation
    original_data = np.copy(unscaled_data)
    if bool_temperature:
        for i, lambda_val in enumerate(lambdas.values()):
            if transformed_data.ndim == 1:
                original_data[i+1] = inv_boxcox(unscaled_data[i+1], lambda_val) - small_value
            else:
                original_data[:, i+1] = inv_boxcox(unscaled_data[:, i+1], lambda_val) - small_value
    else:
        for i, lambda_val in enumerate(lambdas.values()):
            if lambda_val != 0:
                if transformed_data.ndim == 1:
                    original_data[i] = inv_boxcox(unscaled_data[i], lambda_val) - small_value
                else:
                    original_data[:, i] = inv_boxcox(unscaled_data[:, i], lambda_val) - small_value
    if original_data.shape[0] == 1:
        original_data = original_data.flatten()
    return original_data


def reverse_output(transformed_data, scaler):
    """
    Reverses the StandardScaler on the output data.
    Parameters:
    - transformed_data: np.ndarray, the transformed data.
    - scaler: sklearn.preprocessing.StandardScaler, the scaler used for output data.
    Returns:
    - original_data: np.ndarray, the data in its original scale.
    """
    transformed_data = transformed_data.reshape(1, -1) if transformed_data.ndim == 1 else transformed_data
    original_data = scaler.inverse_transform(transformed_data)
    if original_data.shape[0] == 1:
        original_data = original_data.flatten()

    return original_data


def load_scaler_from_csv(scaler_path, scaler_type):
    from sklearn.preprocessing import StandardScaler
    from sklearn.preprocessing import MinMaxScaler
    scaler_df = pd.read_csv(scaler_path)
    if scaler_type == "StandardScaler":
        scaler = StandardScaler()
        scaler.mean_ = scaler_df['mean'].values
        scaler.var_ = scaler_df['var'].values
        scaler.scale_ = np.sqrt(scaler_df['var'].values)
    elif scaler_type == "MinMaxScaler":
        scaler = MinMaxScaler()
        scaler.min_ = scaler_df['min'].values
        scaler.data_range_ = scaler_df['data_range'].values
        scaler.scale_ = 1.0 / scaler_df['data_range'].values
    else:
        raise ValueError(f"Unsupported scaler type: {scaler_type}")

    return scaler

def print_data_info(inputs_train, labels_train, inputs_test, labels_test):
    """
    Prints the shape and first few rows of the training and test data.

    Parameters:
    - inputs_train: ndarray, training inputs.
    - labels_train: ndarray, training labels.
    - inputs_test: ndarray, test inputs.
    - labels_test: ndarray, test labels.
    """
    print(f"Training input shape: {inputs_train.shape}")
    print(f"Training label shape: {labels_train.shape}")
    print(f"Test input shape: {inputs_test.shape}")
    print(f"Test label shape: {labels_test.shape}")

    print("First few rows of training data:")
    print(inputs_train[:3])
    print(labels_train[:3])

    print("First few rows of test data:")
    print(inputs_test[:3])
    print(labels_test[:3])

def parse_and_prepare():
    parser = argparse.ArgumentParser(description='Train and evaluate models.')
    parser.add_argument('json_path', type=str, help='Path to the JSON configuration file')
    args = parser.parse_args()
    
    with open(args.json_path, 'r') as file:
        props = json.load(file)
    
    output_folder = os.path.join("models", props['output_folder'])
    output_dir = os.path.join(output_folder, props['output_path'])
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    props['output_dir'] = output_dir
    props['json_path'] = args.json_path
    print("output_dir=",props['output_dir'])    
    print("json_path=", props['json_path'])
    # copy_file(args.json_path, output_dir, overwrite=True)

    return props

class TorchScaler:
    def __init__(self, scaler):
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.mean_ = torch.tensor(scaler.mean_, dtype=torch.float32, device=device)
        self.scale_ = torch.tensor(scaler.scale_, dtype=torch.float32, device=device)
        self.var_ = torch.tensor(scaler.var_, dtype=torch.float32, device=device)
    
    def inverse_transform(self, X):
        return X * self.scale_ + self.mean_


def update_json_grid_size(json_path, new_grid_size):
    with open(json_path, 'r') as f:
        data = json.load(f)
    data['grid_size'] = int(new_grid_size)
    
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=4)

def update_json_LR(json_path, new_LR):
    with open(json_path, 'r') as f:
        data = json.load(f)
    data['learning_rate'] = float(new_LR) 
    
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=4)
        
        
def get_nonpremixd_counterflow_noCEM(gas, props):
    """ Get 1D non-premixed counterflow flame
    """
    f = ct.CounterflowDiffusionFlame(gas, width=props['width'])

    # Define the operating pressure and boundary conditions
    f.P = props['P']
    f.fuel_inlet.mdot = 0.02884 /2 # kg/m^2/s
    f.fuel_inlet.X = props['fuel']
    f.fuel_inlet.T = props['Tfuel']
    f.oxidizer_inlet.mdot = 0.025 /2 # kg/m^2/s
    f.oxidizer_inlet.X = props['oxid']
    f.oxidizer_inlet.T = props['Toxid']

    # strain rate is du/dx, and a = max(du/dt) at the cold side

    f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
    f.set_initial_guess()
    
    temperature_limit_extinction = props['Text'] #[K]
    def interrupt_extinction(t):
        if np.max(f.T) < temperature_limit_extinction:
            raise FlameExtinguished('Flame extinguished')
        return 0.
    f.set_interrupt(interrupt_extinction)

    # Initialize and solve
    print('Creating the initial solution')
    f.solve(loglevel=0, auto=True)

    # STRAIN RATE LOOP
    strain_factor = 1.2
    

    # Exponents for the initial solution variation with changes in strain rate
    exp_d_a = - 1. / 2.
    exp_u_a = 1. / 2.
    exp_V_a = 1.
    exp_lam_a = 2.
    exp_mdot_a = 1. / 2.

    # Restore initial solution
    file_name = os.path.join(props['path'], "flame_strain_rate_iter=%03d.csv"%0)
    f.write_csv(file_name, quiet=False)

    # Do the strain rate loop
    n = 0
    while np.max(f.T) > temperature_limit_extinction:
        n += 1
        print('strain rate iteration', n)
        # Create an initial guess based on the previous solution and Update grid
        f.flame.grid *= strain_factor ** exp_d_a
        normalized_grid = f.grid / (f.grid[-1] - f.grid[0])
        
        f.fuel_inlet.mdot *= strain_factor ** exp_mdot_a
        f.oxidizer_inlet.mdot *= strain_factor ** exp_mdot_a
        f.set_profile('velocity', normalized_grid, f.velocity * strain_factor ** exp_u_a)
        f.set_profile('spread_rate', normalized_grid, f.spread_rate * strain_factor ** exp_V_a)
        f.set_profile('lambda', normalized_grid, f.L * strain_factor ** exp_lam_a) # Update pressure curvature
        try:
            # Try solving the flame
            f.solve(loglevel=0)#, auto=True)
            file_name =  os.path.join(props['path'], "flame_strain_rate_iter=%03d.csv"%n)
            f.write_csv(file_name, quiet=False)
        except FlameExtinguished:
            print('Flame extinguished')
            break
        except ct.CanteraError as e:
            print('Error occurred while solving:', e)
            break
        
        


class FlameExtinguished(Exception):
    pass


def analyze_flame(gas, props, data, figs_path=None, doplot=False):
    """ Analyze and plot flame resutls
        data: x, u, spread_rate, lambda, T, rho, Y1, Y2, ..., Yn
    """
    xpos = data[:,0] # grid x [m]
    xvel = data[:,1] # velocity u [m/s]
    T = data[:,4]    # temperature T [K]
    MF = np.zeros_like(T)

    # convert X to Y, and get MF
    X = data[:, 6:-5]
    for i,x in enumerate(xpos):
        gas.TPX = T[i], props['P'], X[i,:]
        MF[i] = gas.mixture_fraction(props['fuel'], props['oxyd'])

    # compute progress variable (not normalized)
    C = data[:, 6+gas.species_index("CO")] + data[:, 6+gas.species_index("CO2")]

    # please refer to cantera's jupyter notebook example for the compuation of a
    # https://cantera.org/examples/jupyter/flames/twin_premixed_flame_axisymmetric.ipynb.html
    pos_maxdudx, a = diffMax(xpos, xvel)
    pos_minu = xvel[:pos_maxdudx].argmin()
    pos_a, a = curvMax(xpos[:pos_minu], xvel[:pos_minu])

    # plot anlysis results
    # fig, axs = plt.subplots(3, 1, figsize=c2i(6,9))
    # fig.suptitle(name)
    
    # axs[0].plot(C, T)
    # axs[0].set_xlim([0,1])
    # axs[0].set_xlabel(r"Progress Variable C")
    # axs[0].set_ylabel("Temperature [K]")
    
    # axs[1].plot(C, data[:, 6+gas.species_index("CO")])
    # axs[1].set_xlim([0,1])
    # axs[1].set_xlabel(r"Progress Variable C")
    # axs[1].set_ylabel(r"$Y_{CO}$")
    
    # axs[2].plot(C, data[:, 6+gas.species_index("CO2")])
    # axs[2].set_xlim([0,1])
    # axs[2].set_xlabel(r"Progress Variable C")
    # axs[2].set_ylabel(r"$Y_{CO2}$ ")

    # fig.subplots_adjust(top=0.93, bottom=0.10, left=0.19, right=0.97, hspace=0.40, wspace=0.20)
    # fig.savefig(figs_path+"post_"+name+".png")
    
    # if doplot == True:
    #     plt.show()
    # plt.close()
    
    return MF, C, T, a, data

def analyze_premixed_flame(gas, data, P, fuel_component, oxid_component, name="", figs_path=None, doplot=False):
    """ Analyze and plot flame resutls
        data: x, u, spread_rate, lambda, T, rho, Y1, Y2, ..., Yn
    """
    xpos = data[:,0] # grid x [m]
    xvel = data[:,1] # velocity u [m/s]
    T = data[:,4]    # temperature T [K]
    MF = np.zeros_like(T)

    # convert X to Y, and get MF
    Y = data[:, 6:-3]
    for i,x in enumerate(xpos):
        gas.TPY = T[i], P, Y[i,:]
        MF[i] = gas.mixture_fraction(fuel_component, oxid_component)

    # compute progress variable (not normalized)
    # C = data[:, 6+gas.species_index("CO")] + data[:, 6+gas.species_index("CO2")]

    # please refer to cantera's jupyter notebook example for the compuation of a
    # https://cantera.org/examples/jupyter/flames/twin_premixed_flame_axisymmetric.ipynb.html
    pos_maxdudx, a = diffMax(xpos, xvel)
    pos_minu = xvel[:pos_maxdudx].argmin()
    pos_a, a = curvMax(xpos[:pos_minu], xvel[:pos_minu])

    # plot anlysis results
    # fig, axs = plt.subplots(3, 1, figsize=c2i(6,9))
    # fig.suptitle(name)
    
    # axs[0].plot(C, T)
    # axs[0].set_xlim([0,1])
    # axs[0].set_xlabel(r"Progress Variable C")
    # axs[0].set_ylabel("Temperature [K]")
    
    # axs[1].plot(C, data[:, 6+gas.species_index("CO")])
    # axs[1].set_xlim([0,1])
    # axs[1].set_xlabel(r"Progress Variable C")
    # axs[1].set_ylabel(r"$Y_{CO}$")
    
    # axs[2].plot(C, data[:, 6+gas.species_index("CO2")])
    # axs[2].set_xlim([0,1])
    # axs[2].set_xlabel(r"Progress Variable C")
    # axs[2].set_ylabel(r"$Y_{CO2}$ ")

    # fig.subplots_adjust(top=0.93, bottom=0.10, left=0.19, right=0.97, hspace=0.40, wspace=0.20)
    # fig.savefig(figs_path+"post_"+name+".png")
    
    # if doplot == True:
    #     plt.show()
    # plt.close()
    
    return MF, T, a, data


def computeStrainRatesTwinPremixedFlame(premixFlame, rhou):
    uu = premixFlame.reactants.mdot / rhou
            # estimate strain rate
    zz = premixFlame.flame.grid
    dz = zz[-1] - zz[0]
    a = 2 * uu / dz
    L = - rhou * a**2

    return a
# def computeStrainRates(oppFlame):
#     # Compute the derivative of axial velocity to obtain normal strain rate
#     strainRates = derivative(oppFlame.grid, oppFlame.velocity)

#     # Obtain the location of the max. strain rate upstream of the pre-heat zone.
#     # This is the characteristic strain rate
#     maxStrLocation = abs(strainRates).argmax()
#     minVelocityPoint = oppFlame.velocity[:maxStrLocation].argmin()

#     # Characteristic Strain Rate = K
#     strainRatePoint = abs(strainRates[:minVelocityPoint]).argmax()
#     # strainRatePoint= None
#     K = abs(strainRates[strainRatePoint])

#     return strainRates, strainRatePoint, K

def computeConsumptionSpeed(oppFlame):

    Tb = max(oppFlame.T)
    Tu = min(oppFlame.T)
    rho_u = max(oppFlame.density)

    integrand = oppFlame.heat_release_rate/oppFlame.cp

    total_heat_release = np.trapz(integrand, oppFlame.grid)
    try:
        # Sc = total_heat_release/(Tb - Tu)/rho_u
        Sc = total_heat_release/rho_u
        
        if math.isinf(Sc):
            Sc = 99999
    except ZeroDivisionError:
        Sc = 99999

    return Sc


def derivative(x, y):
    dydx = np.zeros(y.shape, y.dtype.type)

    dx = np.diff(x)
    dy = np.diff(y)
    dydx[0:-1] = dy/dx

    dydx[-1] = (y[-1] - y[-2])/(x[-1] - x[-2])

    return dydx

