# -*- coding: utf-8 -*-
"""
Created on Tue 4 June 2019

@author: Sandor Zoletnik and Mate Lampert

This is the general flap module for MDSPlus archives.
"""
import numpy as np
import copy
import configparser
import io
import pickle
import os
import math
import MDSplus as mds

import flap
from flap.tools import *

class FlapEFITObject(dict):
    """
    Wrapper object for all the EFIT nodes present in the configuration file.
    """
    def __getattr__(self,key):
        return self.get(key)
 
    def __setattr__(self,key,value):
        self[key] = value    
        
    def get_data(self, data_source, exp_id=None, options=None):
        default_options = {'Verbose': False}
        _options = flap.config.merge_options(default_options,options,data_source='NSTX_GPI')
        try:
            efit_dictionary=flap.config.get_all_section('EFIT')
        except:
            raise ValueError("Couldn't read the EFIT section of the configuration file.")
        
        for index in efit_dictionary:
            if not (index == None):
                try:
                    self[index]=flap.get_data(data_source,
                                                   name=efit_dictionary[index],
                                                   exp_id=exp_id,
                                                   object_name=index)
                except:
                    self[index]=None
                    if (_options['Verbose']):
                        print("Couldn't read "+efit_dictionary[index]+ " for shot "+str(exp_id))
    def load(self, filename=None):
        raise NotImplemented("Not implemented yet.")
        
    def save(self, filename=None):
        raise NotImplemented("Not implemented yet.")

def mds_virtual_names(data_name, exp_id, channel_config_file):

    """
    Translates virtual data names to MDSPlus entries.
    Returns a list of virtual names and translated full names.
    The list is read from the configuration file. The file is a standard
    configuration file with one section [Virtual names]. Each entry looks
    like one of the following templates.
        Simple translation:
            <name> = <mdsname>
        Translation in an expID range, including ends:
            <name>(startID-endID) = <mdsname>
          One of startID or endID may be omitted.
        <mdsname> can be either a single MDS name or complex(mds1, mds2) which
        will construct a complex signal from the two MDS+ entries.
    Wildcards are allowed in data names. Data names can be a single string or
        list of strings.
    Return value:
        virt_names, mds_str, mds_names
        virt names is a list of the interpreted names.
        mds_str is the MDSPlus description from the virtual name translation file
        mds_names is a list of the same length as virt_names. Elements can be the following:
        Elements can be the following:
            string: a single MDS+ name
            list: [<type>, <mds1>, <mds2>, ...]
                type can be
                  'complex': two MDS entries are expected (real, imag), complex
                  signal will be created from them.

    """
    if ((type(exp_id) is not int) and (type(exp_id) is not str)):
        raise TypeError("exp_id should be integer or string.")

    if (type(exp_id) is str):
        try:
            exp_id_num = int(exp_id[:8]+exp_id[9:])
        except Exception as e:
            raise ValueError("Invalid exp_id value: '{:s}".format(exp_id))
    if (type(exp_id) is int):
        exp_id_num=exp_id
  
    config = configparser.ConfigParser()
    config.optionxform = str
    read_ok = config.read(channel_config_file)
    if (read_ok == []):
        raise ValueError("Invalid MDSPlus virtual name file "+channel_config_file) 
    try:
        entries = config.items('Virtual names')
    except Exception as e:
        raise ValueError("Invalid mds virtual name file "+channel_config_file)
    entry_descr = [e[0] for e in entries]
    entry_values = [e[1] for e in entries]
    mds_names = []
    entry_names = []
    # Finding etries with matching exp_id
    for e,ev in zip(entry_descr, entry_values):
        start_brace = e.find('(')
        stop_brace = e.find(')')
        if (start_brace * stop_brace < 0):
            raise ValueError("Invalid entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
        if (start_brace > 0):
            # ExpID range given
            exp_id_range = e[start_brace+1:stop_brace].split('-')
            if (len(exp_id_range) != 2):
                raise ValueError("Invalid exp_id range in entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
            if (exp_id_range[0].strip() == ''):
                exp_id_start = None
            else:
                try:
                    exp_id_start = int(exp_id_range[0][:8]+exp_id_range[0][9:])
                except ValueError:
                    raise ValueError("Invalid exp_id start in entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
            if (exp_id_range[1].strip() == ''):
                exp_id_stop = None
            else:
                try:
                    exp_id_stop = int(exp_id_range[1][:8]+exp_id_range[1][9:])
                except ValueError:
                    raise ValueError("Invalid exp_id stop in entry '{:s}' in virtual name file: {:s}".format(e, channel_config_file))
            if ((exp_id_start is not None) and (exp_id_num < exp_id_start) or \
                (exp_id_stop is not None) and (exp_id_num > exp_id_stop)) :
                continue
            entry_names.append(e[:start_brace])
        else:
            entry_names.append(e)
        mds_names.append(ev)
    if (type(data_name ) is not list):
        _data_name = [data_name]
    else:
        _data_name = data_name
    select_list = []
    select_mds_list = []
    for i,dn in enumerate(_data_name):    
        try:
            sl, si = flap.select_signals(entry_names, dn)
            select_list += sl
            select_mds_list += [mds_names[i] for i in si]
        except ValueError as e:
            select_list.append(None)
            select_mds_list.append(dn)
            
    mds_descr = []
    for descr in select_mds_list:
        start_brace = descr.find('(')
        stop_brace = descr.find(')')
        if (start_brace * stop_brace < 0):
            raise ValueError("Invalid value '{:s}' in virtual name file: {:s}".format(descr, channel_config_file))
        if (start_brace > 0):
            mds_list = descr[start_brace+1:stop_brace].split(',')
            mds_type = descr[:start_brace]
            mds_descr.append([mds_type] + mds_list)
        else:
            mds_descr.append(descr)    
    return select_list, select_mds_list, mds_descr        

def mdsplus_get_data(exp_id=None, data_name=None, no_data=False, options=None, 
                     coordinates=None, data_source=None):
    """ Data read function for the MDSplus database
    exp_id: exp_id number, integer or YYYYMMDD.xxx
    data_name: Channel names [\]tree::node
               or virtual names:
                   CR-x: Correlation reflectometry, x is antenna A,B,C,D
    coordinates: List of flap.Coordinate() or a single flap.Coordinate
                 Defines read ranges. The following coordinates are interpreted:
                     'Sample': The read samples
                     'Time': The read times
                     Only a single equidistant range is interpreted in c_range.
    options: Dictionary. Defaults will be read from <data_source> section in configuration file.
            'Protocol': For ssh connection use 'ssh://'.
            'Server': Server name (default: mds-trm-1.ipp-hgw.mpg.de)
            'User': User name for access. Password-free access should be set up for this user.
            'Virtual name file': A file name to translate virtual names to MDS+ entries. For 
                                 format see mds_virtual_names()
            'Verbose': (bool) Write progress information during data read.
            'Cache data': (bool) Cache data to options['Cache directory'] and read it from there
            'Cache directory': (str) Name of the cache directory
    """
    if (exp_id is None):
        raise ValueError('exp_id should be set for reading mdsplus data.')

    default_options = {'Server': None,
                       'Protocol': None,
                       'User': None,
                       'Virtual name file': None,
                       'Verbose': True,
                       'Cache data': False,
                       'Cache directory': None
                       }
    
    _options = flap.config.merge_options(default_options,options,data_source=data_source)
    if (data_source is None):
        data_source = 'MDSPlus'
        
    if (exp_id is None):
        raise ValueError('exp_id should be set for MDSPlus.')
    if (type(exp_id) == int):
        exp_id_mds = exp_id
    if (type(exp_id) is str):       #This is for w7X, could be changed to find and option in the config file
        exp_id_split = exp_id.split('.')
        if ((len(exp_id_split) is not 2) or (len(exp_id_split[0]) != 8) or (len(exp_id_split[1]) != 3)):
            raise ValueError("exp_id format error: must be a string YYYYMMDD.nnn")
            exp_id_mds = int(exp_id_split[0][2:] + exp_id_split[1])
    elif (type(exp_id) is not int):
        raise TypeError("exp_id must be an integer or string")

    if (_options['Server'] is None):
        raise ValueError("Option 'Server' should be set for using MDSPlus.")                    
    #if no username and protocol then open a server (e.g. NSTX or KSTAR)
    if ((_options['Protocol'] is None) and (_options['User'] is None)):
        connection_name = _options['Server']
        
    #if protocol and username exists then use the following syntax (e.g. W7-X)
    if ((_options['Protocol'] is not None) and (_options['User'] is not None)):
        connection_name = _options['Protocol'] + _options['User'] + '@' + _options['Server']
        
    #Error handing if one of the parameters is not set.
    if (((_options['Protocol'] is not None) and (_options['User'] is None)) or
       ((_options['Protocol'] is None) and (_options['User'] is not None))):
        raise ValueError("If Protocol is set then Username must be set, as well.")
    
    if ((type(data_name) is not str) and (type(data_name) is not list)):
        raise ValueError("data_name should be a string or list of strings.")
        
    if (_options['Virtual name file'] is not None):
        try:
            virt_names, virt_mds_txt, virt_mds = mds_virtual_names(data_name, exp_id, _options['Virtual name file'])
        except Exception as e:
            raise e
    else:
        print(_options)
        raise ValueError('No virtual names available.')
        
    read_range = None
    read_samplerange = None
    if (coordinates is not None):
        if (type(coordinates) is not list):
            _coordinates = [coordinates]
        else:
            _coordinates = coordinates
        for coord in _coordinates:
            if (type(coord) is not flap.Coordinate):
                raise TypeError("Coordinate description should be flap.Coordinate.")
            if (coord.unit.name is 'Time'):
                if (coord.mode.equidistant):
                    read_range = [float(coord.c_range[0]),float(coord.c_range[1])]
                else:
                    raise NotImplementedError("Non-equidistant Time axis is not implemented yet.")
                break
                if (coord.unit.unit == 'Millisecond'):
                    read_range = [read_range[0]*1e-3, read_range[1]*1e-3]
                elif (coord.unit.unit == 'Microsecond'):
                    read_range = [read_range[0]*1e-6, read_range[1]*1e-6]
                elif (coord.unit.unit == 'Nanosecond'):
                    read_range = [read_range[0]*1e-9, read_range[1]*1e-9]
                elif (coord.uni.unit == 'Second'):
                    pass
                else:
                    raise ValueError("Unknown time unit '"+coord.unit.unit+"'. Valid: Second, Millisecond, Microsecond, Nanosecond.")

    signal_list = []
    data_list = []
    common_time = None    
    for name, mds_descr in zip(virt_names,virt_mds):
        # Assembling a list of MDS nodes needed for this data
        mds_request_list = []
        if (name is None):
            # This was not recognized as virtual name
            mds_request_list = [mds_descr]
            signal_list.append(mds_descr)
            readtype = 0
        else:
            # This was recongnized as a virtual signal
            signal_list.append(name)
            if (type(mds_descr) is not list):
                # This was recognized as a single MDS node
                mds_request_list = [mds_descr]
                readtype = 0
            else:
                # This is a composite virtual signal
                mds_request_list = mds_descr[1:]
                readtype = 1 
        # Reading the required nodes
        this_data_list = []
        for mds_name in mds_request_list:
            mds_name_split = mds_name.split('::')
            if (len(mds_name_split) != 2):
                raise ValueError("Invalid mds name '{:s}', missing tree name? Data name is tree::node".format(mds_name))
            tree_name = mds_name_split[0].strip()
            if (tree_name[0] == '\\'):
                tree_name = tree_name[1:]
            node_name = mds_name_split[1]
            if ((_options['Cache data']) and (_options['Cache directory'] != None)):
                filename = str(exp_id)+'_'+mds_name
                for c in ['\\',':']:
                   filename = filename.replace(c,'_')
                #filename = os.path.join(_options['Cache directory'],filename+'.pickle')
                directory=os.path.join(_options['Cache directory'],str(exp_id))
                if not (os.path.exists(directory)):
                    try:
                        os.mkdir(directory)
                    except:
                        raise SystemError("The shot folder cannot be created. Cache directory might not be present.")
                                      
                filename = os.path.join(directory,filename+'.pickle')
                try:
                    f = io.open(filename,'rb')
                    mdsdata_pickle = pickle.load(f)
                    f.close()
                    try:
                        if (mdsdata_pickle['MDSdata cache']):
                            mdsdata = mdsdata_pickle['Data']
                            mdsdata_unit = mdsdata_pickle['Data unit']
                            mdsdata_spat = mdsdata_pickle['Data dimension']
                            mdsdata_spat_unit = mdsdata_pickle['Data dimension unit']
                            mdsdata_time = mdsdata_pickle['Time']
                            mdsdata_time_unit = mdsdata_pickle['Time unit']
                            mdsdata_time_step = mdsdata_pickle['Time step']
                            del mdsdata_pickle
                            data_cached = True
                    except:
                        data_cached = False
                except:
                    data_cached = False
            else:
                data_cached = False
            if (not data_cached):
                try:
                    conn
                except NameError:
                    try:
                        if (_options['Verbose']):
                            print("Connecting to "+connection_name)
                        conn = mds.Connection(connection_name)
                    except Exception as e:
                        raise e
                    try:
                        conn.openTree(tree_name,exp_id_mds)
                    except mds.MdsException as e:
                        raise RuntimeError("Error connecting to tree {:s}, experiment {:s}".format(tree_name,str(exp_id))) 
                if (_options['Verbose']):
                    print("Reading "+node_name)
                try:
                    mdsdata = conn.get(node_name).data() # NEEDS TO BE CHECKED IF IT WORKS FOR W7-X
                    mdsdata_unit = conn.get('units('+node_name+')').data()
                    mdsdata_spat=[]
                    mdsdata_spat_unit=[]
                    for dim_ind in range(1,len(mdsdata.shape)):
                        mdsdata_spat.append(conn.get('dim_of('+node_name+','+str(dim_ind)+')').data())
                        mdsdata_spat_unit.append(conn.get('units(dim_of('+node_name+','+str(dim_ind)+'))').data())
                        #print(conn.get('dim_of('+node_name+','+str(dim_ind)+')').data().shape)
                    
                    mdsdata_time = conn.get('dim_of('+node_name+',0)').data()
                    mdsdata_time_unit = conn.get('units(dim_of('+node_name+'))').data()
                    #print(np.asarray(mdsdata).shape,
                    #      np.asarray(mdsdata_spat).shape,
                    #      np.asarray(mdsdata_time).shape,
                    #      mdsdata_spat)
                    if not (len(mdsdata_time) < 2):
                        mdsdata_time_step=mdsdata_time[1]-mdsdata_time[0]
                    else:
                        mdsdata_time_step=0.
#                    t = mdsdata_time*mdsdata_time_unit
#                    ind = np.nonzero(np.logical_and(t > 4.4,t < 4.8))[0]
#                    try:
#                        figcount
#                        figcount += 1
#                    except NameError:
#                        figcount = 1
#                    plt.figure(figcount)
#                    plt.plot(t[ind],mdsdata[ind])
#                    plt.title(mds_name)
                #except mds.MDSplusException as e:
                except:
                    raise RuntimeError("Cannot read MDS node: {:s}".format(node_name))
            if (not data_cached and (_options['Cache data']) and (_options['Cache directory'] is not None)):
                while True:
                    try:
           
                        f = io.open(filename,"wb")
                    except:
                        print("Warning: Cannot open cache file: "+filename)
                        break
                    mdsdata_pickle = {}
                    mdsdata_pickle['MDSdata cache'] = True
                    mdsdata_pickle['Data'] = copy.deepcopy(mdsdata)
                    mdsdata_pickle['Data unit'] = mdsdata_unit
                    mdsdata_pickle['Data dimension'] = mdsdata_spat
                    mdsdata_pickle['Data dimension unit'] = mdsdata_spat_unit
                    mdsdata_pickle['Time'] = mdsdata_time
                    mdsdata_pickle['Time unit'] = mdsdata_time_unit
                    mdsdata_pickle['Time step'] = mdsdata_time_step
                    try:
                        pickle.dump(mdsdata_pickle,f)
                        del mdsdata_pickle
                    except Exception as e:
                        print("Warning: Cannot write cache file: "+filename)
                        break
                    try:
                        f.close()
                    except Exception as e:
                        print("Warning: Cannot write cache file: "+filename)
                    break
                
            if (read_range is not None):
                read_ind = np.nonzero(np.logical_and(mdsdata_time * mdsdata_time_step >= read_range[0],
                                                     mdsdata_time * mdsdata_time_step <= read_range[1]
                                                     )
                                      )[0]
                mdsdata_time = mdsdata_time[read_ind]
                mdsdata = mdsdata[read_ind]
            else:
                read_ind = [0, len(mdsdata)]
            
            mdsdata_time_unit_int = time_unit_translation(mdsdata_time_unit)
            if (common_time is not None):
                if ((len(common_time) != len(mdsdata_time) or \
                    (math.fabs(common_time_unit - mdsdata_time_step)) / common_time_unit > 0.001) or \
                    (np.nonzero(np.abs(common_time - mdsdata_time) \
                        > math.fabs(common_time[1] - common_time[0]) * 0.1)[0].size != 0)):
                    raise ValueError("Different timescales for signals. Not possible to return in one flap.DataObject.")
            else:
                common_time = mdsdata_time
                common_time_unit = mdsdata_time_unit_int
            
            this_data_list.append(mdsdata)
            del mdsdata
        if (readtype == 0):
            data_list.append(this_data_list[0])
        elif (readtype == 1):
            data_list.append(this_data_list[0] + 1j * this_data_list[1])
    # Determining data type
    dtype = int
    for i in range(len(data_list)):
        if (dtype is not complex) and (data_list[i].dtype.kind == 'f'):
                dtype = float
        if (data_list[i].dtype.kind == 'c'):
            dtype = complex
    if (len(data_list) == 1):
        data = data_list[0]
        signal_dim = []
    else:    
        data = np.empty((len(data_list[0]),len(data_list)),dtype=dtype)
        for i in range(len(data_list)):
            data[:,i] = data_list[i].astype(dtype)
        signal_dim = [1]
        
            
    dt = (common_time[1:] - common_time[:-1])
    if (np.nonzero(np.abs(dt - dt[0]) / dt[0] > 0.01)[0].size != 0):
        coord_type = flap.CoordinateMode(equidistant=False)
    else:
        coord_type = flap.CoordinateMode(equidistant=True)
    
    coord = []
    if (coord_type.equidistant):
        coord.append(copy.deepcopy(flap.Coordinate(name='Time',
                                                   unit=mdsdata_time_unit,
                                                   mode=coord_type,
                                                   shape = [],
                                                   start=common_time[0] * common_time_unit,
                                                   step=(common_time[1] - common_time[0]) * common_time_unit,
                                                   dimension_list=[0])
                                    ))
    else:
        coord.append(copy.deepcopy(flap.Coordinate(name='Time',
                                                   unit=mdsdata_time_unit,
                                                   mode=coord_type,
                                                   values=common_time * common_time_unit,
                                                   shape = common_time.shape,
                                                   dimension_list=[0])
                                    ))
    if ((len(data.shape) > 1) and (len(np.asarray(mdsdata_spat).shape) > 2)):
        for dim_ind in range(0,len(data.shape)-1):
            #now only the zeroth time is gotten from the vector
            spatial_data=np.asarray(mdsdata_spat[dim_ind][0][:])
            if (((mdsdata_spat_unit[dim_ind] == 'm') or
                (mdsdata_spat_unit[dim_ind] == 'mm')) and
                (dim_ind < 3)):
                names=['Device R','Device z']
        
                coord.append(copy.deepcopy(flap.Coordinate(name=names[dim_ind],
                                                           unit=mdsdata_spat_unit[dim_ind],
                                                           mode=flap.CoordinateMode(equidistant=False),
                                                           values=spatial_data,
                                                           shape=spatial_data.shape,
                                                           dimension_list=[dim_ind+1])))                                                           
            else:
                coord.append(copy.deepcopy(flap.Coordinate(name='Dimension '+str(dim_ind+1),
                                                           unit=mdsdata_spat_unit[dim_ind],
                                                           mode=flap.CoordinateMode(equidistant=False),
                                                           values=spatial_data,
                                                           shape=spatial_data.shape,
                                                           dimension_list=[dim_ind+1])))
    coord.append(copy.deepcopy(flap.Coordinate(name='Sample',
                                               unit=mdsdata_unit,
                                               mode=flap.CoordinateMode(equidistant=True),
                                               start=read_ind[0],
                                               step=1,
                                               dimension_list=[0])
                                               #dimension_list=[0])
                               ))
    coord.append(copy.deepcopy(flap.Coordinate(name='Signal name',
                                               shape=tuple([len(signal_list)]),
                                               unit='',
                                               mode=flap.CoordinateMode(equidistant=False),
                                               values=signal_list,
                                               dimension_list=signal_dim)
                                 ))
    coord.append(copy.deepcopy(flap.Coordinate(name='Channel',
                                               shape=tuple([len(signal_list)]),
                                               unit='',
                                               mode=flap.CoordinateMode(equidistant=False),
                                               values=np.arange(len(signal_list)) + 1,
                                               dimension_list=signal_dim)))
    
    coord.append(copy.deepcopy(flap.Coordinate(name='mds description',
                                               shape=tuple([len(virt_mds_txt)]),
                                               unit='',
                                               mode=flap.CoordinateMode(equidistant=False),
                                               values=virt_mds_txt,
                                               dimension_list=signal_dim)))
    
    data_title = data_source.split('_')[0]+' '+data_source.split('_')[1]+' data'
    
    d = flap.DataObject(data_array=data,
                        data_unit=flap.Unit(name=''.join(signal_list),unit=mdsdata_unit),
                        coordinates=coord,
                        exp_id=exp_id,
                        data_title=data_title,
                        data_source=data_source)
    return d


def add_coordinate(data_object, coordinates, options=None):

    raise NotImplementedError("Not implemented.")

def register(data_source):
    flap.register_data_source(data_source, get_data_func=mdsplus_get_data,
                              add_coord_func=add_coordinate)
