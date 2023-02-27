##==============================================================================
## BHPTNRSur2dq1e3 : arXiv:XXXX.YYYYY
## Description : loads surrogate fits data
## Author : Katie Rink, Dec 2022 [krink@utexas.edu]
##==============================================================================

import numpy as np
import h5py
import os
from os import path
import hashlib
import copy

#----------------------------------------------------------------------------------------------------
def md5(fname,h5_data_dir):
    """ Compute has from file. code taken from 
    https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file"""
    
    # download file if not already there
    if path.isfile('%s/%s'%(h5_data_dir,fname))==False:
        print('BHPTNRSur1dq1e4.h5 file is not found in the directory - PATH-TO/BHPTNRSurrogate/data/')
        print('... downloading h5 file from zenodo')
        print('... this might take some time')
        #KR - change: os.system('wget https://zenodo.org/record/7125742/files/BHPTNRSur1dq1e4.h5 -P %s'%h5_data_dir)
        #print('... downloaded')
    
    hash_md5 = hashlib.md5()
    with open('%s/%s'%(h5_data_dir,fname), "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def chars_to_string(chars):
    """ Function converts array of  unicode characters to string.
            - Copied from gwsurrogate: surrogatIO.py
            - Needed for reading in hdf5 file data.
    """
    return "".join(chr(cc) for cc in chars)

#----------------------------------------------------------------------------------------------------
def load_surrogate(h5_data_dir):
    
    """ 
    Loads all interpolation data for the following modes
    modes=[(2,2),(2,1),(3,1),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5)]
    #KR - change: ,(6,4),(6,5),(6,6),(7,5),(7,6),(7,7),(8,6),(8,7),(8,8),(9,7),(9,8),(9,9),(10,8),(10,9)]
    #KR - change: Assumes the file BHPTNRSur1dq1e4.h5 is located in the same directory as this file.
    """
    
    file_hash = md5('EMRISur2dq1e4_TEST_02-08.h5',h5_data_dir)
    #KR - final: file_hash = md5('BHPTNRSur2dq1e3.h5',h5_data_dir)
    #KR - change: zenodo_current_hash = "58a3a75e8fd18786ecc88cf98f694d4a"

    #KR - final: if file_hash != zenodo_current_hash:
        #KR - change: raise AttributeError("EMRISur1dq1e4.h5 out of date.\n Please download new version from https://zenodo.org/record/7125742")

    with h5py.File('%s/EMRISur2dq1e4_TEST_02-08.h5'%h5_data_dir, 'r') as f:
    #KR - final: with h5py.File('%s/BHPTNRSur2dq1e3.h5'%h5_data_dir, 'r') as f:

        modes = [(2,2),(2,1),
                   (3,1),(3,2),(3,3),
                   (4,2),(4,3),(4,4)]
                   #(5,3), (5,4), (5,5),
                   #(6,4),(6,5),(6,6),
                   #(7,5),(7,6),(7,7),
                   #(8,6),(8,7),(8,8),
                   #(9,7),(9,8),(9,9),
                   #(10,8),(10,9)]

        # fit info

        '''
        h_eim_re_spline_dict = {}
        h_eim_im_spline_dict = {}
        B_re_dict = {}
        B_im_dict = {}
        eim_indicies_im_dict = {}
        eim_indicies_re_dict = {}

        for mode in modes:

            lmode,mmode=mode

            if mode==(2,2):

                eim_indicies_amp_dataset=f['l%s_m%s/eim_indices'%(lmode,mmode)]
                eim_indicies_amp=eim_indicies_amp_dataset[:]
                eim_indicies_ph_dataset=f['l%s_m%s/eim_indices_phase'%(lmode,mmode)]
                eim_indicies_ph=eim_indicies_ph_dataset[:]
                B_ph_dataset=f['l%s_m%s/B_phase'%(lmode,mmode)]
                B_ph=np.transpose(B_ph_dataset[:])
                B_amp_dataset=f['l%s_m%s/B'%(lmode,mmode)]
                B_amp=np.transpose(B_amp_dataset[:])
                time_dataset=f['l%s_m%s/times'%(lmode,mmode)]
                time=time_dataset[:]
                degree_dataset=f['l%s_m%s/degree'%(lmode,mmode)]
                degree=degree_dataset[:]
                knots_dataset_amp=f['l%s_m%s/spline_knots_amp'%(lmode,mmode)]
                knots_amp=knots_dataset_amp[:]
                knots_dataset_ph=f['l%s_m%s/spline_knots_phase'%(lmode,mmode)]
                knots_ph=knots_dataset_ph[:]
                h_spline_amp_dataset=f['l%s_m%s/fitparams_amp'%(lmode,mmode)]
                h_spline_amp=h_spline_amp_dataset[:]
                h_spline_ph_dataset=f['l%s_m%s/fitparams_phase'%(lmode,mmode)]
                h_spline_ph=h_spline_ph_dataset[:]

                h_eim_amp_spline=[(knots_amp[flag], h_spline_amp[flag],int(degree)) for flag in range(len(eim_indicies_amp))]
                h_eim_ph_spline=[(knots_ph[flag], h_spline_ph[flag],int(degree)) for flag in range(len(eim_indicies_ph))]

            else:

                eim_indicies_re_dataset=f['l%s_m%s/eim_indices'%(lmode,mmode)]
                eim_indicies_re_dict[(mode)]=eim_indicies_re_dataset[:]
                eim_indicies_im_dataset=f['l%s_m%s/eim_indices_im'%(lmode,mmode)]
                eim_indicies_im_dict[(mode)]=eim_indicies_im_dataset[:]
                B_im_dataset=f['l%s_m%s/B_im'%(lmode,mmode)]
                B_im_dict[(mode)]=np.transpose(B_im_dataset[:])
                B_re_dataset=f['l%s_m%s/B'%(lmode,mmode)]
                B_re_dict[(mode)]=np.transpose(B_re_dataset[:])
                time_dataset=f['l%s_m%s/times'%(lmode,mmode)]
                time=time_dataset[:]
                degree_dataset=f['l%s_m%s/degree'%(lmode,mmode)]
                degree=degree_dataset[:]
                knots_dataset_re=f['l%s_m%s/spline_knots_re'%(lmode,mmode)]
                knots_re=knots_dataset_re[:]
                knots_dataset_im=f['l%s_m%s/spline_knots_im'%(lmode,mmode)]
                knots_im=knots_dataset_im[:]
                h_spline_re_dataset=f['l%s_m%s/fitparams_re'%(lmode,mmode)]
                h_spline_re=h_spline_re_dataset[:]
                h_spline_im_dataset=f['l%s_m%s/fitparams_im'%(lmode,mmode)]
                h_spline_im=h_spline_im_dataset[:]

                h_eim_re_spline_dict[(mode)]=[(knots_re[flag], h_spline_re[flag],int(degree)) for flag in range(len(eim_indicies_re_dict[(mode)]))]
                h_eim_im_spline_dict[(mode)]=[(knots_im[flag], h_spline_im[flag],int(degree)) for flag in range(len(eim_indicies_im_dict[(mode)]))]
           '''
# dicts to copy .h5 file data into (needed for surrogate generation)
        h_eim_amp_gpr, h_eim_ph_gpr = {}, {}
        B_amp, B_ph = {}, {}
        eim_indicies_amp, eim_indicies_ph = {}, {}
        time = {}
            
        for spin_sign in ['negative_spin', 'positive_spin']:
            time[spin_sign] = copy.deepcopy(f[spin_sign]["l2_m2"]["times"][()]) # same times for all modes
                
            for mode in modes: # modes defined at top of script
                    
                # Copy data groups we need to access from hdf5 file into output dicts.
                    
                f_mode = f[spin_sign]['l%s_m%s'%(mode[0], mode[1])]
                h_amp_file, h_ph_file = dict(f_mode['gpr_amp']), dict(f_mode['gpr_phase'])
                h_eim_amp_gpr[mode], h_eim_ph_gpr[mode] = {}, {} 

                B_amp[mode] = copy.deepcopy(f_mode["B_amp"][()])
                B_ph[mode] = copy.deepcopy(f_mode["B_phase"][()])
                eim_indicies_amp[mode] = copy.deepcopy(f_mode["eim_indicies_amp"][()])
                eim_indicies_ph[mode] = copy.deepcopy(f_mode["eim_indicies_phase"][()])
                    
                # iterate through GPR fit params once per data piece (amp, ph)
                    
                h_params_file = [h_amp_file, h_ph_file]
                h_params_gpr = [h_eim_amp_gpr[mode], h_eim_ph_gpr[mode]]
                eim_indicies = [eim_indicies_amp[mode], eim_indicies_ph[mode]]
                    
                for i in range(2): # run through once for amp data, then again for phase data
                    h_file = h_params_file[i] # hdf5 source data
                    h_gpr = h_params_gpr[i] # empty dictionary
                        
                    # iterate through time interpolation nodes (EIM indicies)
                    for eim_indx in range(len(eim_indicies[i])):

                        h_file_node = h_file['node%s'%eim_indx]
                        h_file_node_params = h_file_node['GPR_params']

                        h_gpr['node%s'%eim_indx] = {}
                        h_gpr_node = h_gpr['node%s'%eim_indx]
                            
                        h_gpr_node['GPR_params'] = {}
                        h_gpr_node_params = h_gpr_node['GPR_params']

                        h_gpr_node['lin_reg_params'] = {}
                        h_gpr_node_params['kernel_'] = {}
                        h_gpr_node_params['kernel_']['k1__k1'] = {}
                        h_gpr_node_params['kernel_']['k2'] = {}
                        h_gpr_node_params['kernel_']['k1__k2'] = {}
                        h_gpr_node_params['kernel_']['k1'] = {}
                        h_gpr_node_params['kernel_']['k1']['k1'] = {}
                        h_gpr_node_params['kernel_']['k1']['k2'] = {}

                        h_gpr_node['fitType'] = chars_to_string(h_file['fitType'][()])
                        h_gpr_node['data_mean'] = copy.deepcopy(h_file_node['data_mean'][()])
                        h_gpr_node['data_std'] = copy.deepcopy(h_file_node['data_std'][()])

                        h_gpr_node['lin_reg_params']['coef_'] = copy.deepcopy(h_file_node['lin_reg_params']['coef_'][()])
                        h_gpr_node['lin_reg_params']['intercept_'] = copy.deepcopy(h_file_node['lin_reg_params']['intercept_'][()])

                        h_gpr_node_params['X_train_'] = copy.deepcopy(h_file_node_params['X_train_'][()])
                        h_gpr_node_params['alpha_'] = copy.deepcopy(h_file_node_params['alpha_'][()])
                        h_gpr_node_params['_y_train_mean'] = copy.deepcopy(h_file_node_params['_y_train_mean'][()])
                        h_gpr_node_params['L_'] = copy.deepcopy(h_file_node_params['L_'][()])

                        h_gpr_node_params['kernel_']['name'] = chars_to_string(h_file_node_params['kernel_']['name'])
                        h_gpr_node_params['kernel_']['k2__noise_level'] = copy.deepcopy(h_file_node_params['kernel_']['k2__noise_level'][()])
                        h_gpr_node_params['kernel_']['k2__noise_level_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k2__noise_level_bounds'][()])
                        h_gpr_node_params['kernel_']['k1__k2__length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2__length_scale'][()])
                        h_gpr_node_params['kernel_']['k1__k2__length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2__length_scale_bounds'][()])
                        h_gpr_node_params['kernel_']['k1__k1__constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1__constant_value'][()])
                        h_gpr_node_params['kernel_']['k1__k1__constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1__constant_value_bounds'][()])

                        h_gpr_node_params['kernel_']['k1']['name'] = chars_to_string(h_file_node_params['kernel_']['k1']['name'])
                        h_gpr_node_params['kernel_']['k1']['k1__constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1__constant_value'][()])
                        h_gpr_node_params['kernel_']['k1']['k1__constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1__constant_value_bounds'][()])
                        h_gpr_node_params['kernel_']['k1']['k2__length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2__length_scale'][()])
                        h_gpr_node_params['kernel_']['k1']['k2__length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2__length_scale_bounds'][()])

                        h_gpr_node_params['kernel_']['k1']['k1']['name'] = chars_to_string(h_file_node_params['kernel_']['k1']['k1']['name'])
                        h_gpr_node_params['kernel_']['k1']['k1']['constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1']['constant_value'][()])
                        h_gpr_node_params['kernel_']['k1']['k1']['constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k1']['constant_value_bounds'][()])                        
                        h_gpr_node_params['kernel_']['k1']['k2']['name'] = chars_to_string(h_file_node_params['kernel_']['k1']['k2']['name'])
                        h_gpr_node_params['kernel_']['k1']['k2']['length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2']['length_scale'][()])
                        h_gpr_node_params['kernel_']['k1']['k2']['length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1']['k2']['length_scale_bounds'][()])     

                        h_gpr_node_params['kernel_']['k1__k1']['name'] = chars_to_string(h_file_node_params['kernel_']['k1__k1']['name'])
                        h_gpr_node_params['kernel_']['k1__k1']['constant_value'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1']['constant_value'][()])
                        h_gpr_node_params['kernel_']['k1__k1']['constant_value_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k1']['constant_value_bounds'][()])

                        h_gpr_node_params['kernel_']['k1__k2']['name'] = chars_to_string(h_file_node_params['kernel_']['k1__k2']['name'])
                        h_gpr_node_params['kernel_']['k1__k2']['length_scale'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2']['length_scale'][()])
                        h_gpr_node_params['kernel_']['k1__k2']['length_scale_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k1__k2']['length_scale_bounds'][()])

                        h_gpr_node_params['kernel_']['k2']['name'] = chars_to_string(h_file_node_params['kernel_']['k2']['name'])
                        h_gpr_node_params['kernel_']['k2']['noise_level'] = copy.deepcopy(h_file_node_params['kernel_']['k2']['noise_level'][()])
                        h_gpr_node_params['kernel_']['k2']['noise_level_bounds'] = copy.deepcopy(h_file_node_params['kernel_']['k2']['noise_level_bounds'][()])

        # nr calibration info

        alpha_coeffs = {}
        for mode in [(2,2),(3,3),(4,4),(5,5)]:
            #KR - final: alpha_coeffs[mode] = f["nr_calib_params/(%d,%d)"%(mode[0],mode[1])]['alpha'][:]
            alpha_coeffs[mode] = 1
        #KR - final: beta_coeffs = f["nr_calib_params/(2,2)"]['beta'][:]
        beta_coeffs = 1

    return time, eim_indicies_amp, eim_indicies_ph, B_amp, B_ph, h_eim_amp_gpr, h_eim_ph_gpr, alpha_coeffs, beta_coeffs \
#        eim_indicies_re_dict, eim_indicies_im_dict, B_re_dict, B_im_dict, h_eim_re_spline_dict, h_eim_im_spline_dict,\ alpha_coeffs, beta_coeffs
    
