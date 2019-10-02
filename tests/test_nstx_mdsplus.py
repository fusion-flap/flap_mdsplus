# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:37:34 2019

@author: Mate Lampert  (mlampert@pppl.gov)

Test program for NSTX MDSPlus.

Before usinig this program the following has to be set in test_w7x_mdsplus.cfg in the same directory 
where this program resides:
[Module NSTX_MDSPlus]
 Server = 'skylark.pppl.gov'
 User = 'mlampert'
 Virtual name file = 'Virtual_names.cfg'
 Cache data =  True
 Cache directory = ''
 
 'Server' is the MDSPlus server name, 'User' is the user name. Password-free access should be set up for this user..
 'Virtual name file' points to a file which translated virrtual signal names to MDSPlus entries. An example 
 is included. It should be in the working directory or the path should be added to the entry.
 'Cache data' can be set to True, in this case loaded MDSplus data will be stored locally in directory 
 'Cache directory'. When the same signal read next time it will be loaded automatically from this cache.
"""

import matplotlib.pyplot as plt
import os

import flap
import flap_mdsplus

flap_mdsplus.register('NSTX_MDSPlus')

def test_mdsplus(): 
    plt.close('all')
    print("**** Reading an explicit MDS signal.")
    flap.delete_data_object('*')
    try:
        # Explicit MDSPlus reference
       #d=flap.get_data('NSTX_MDSPlus',
       #                 name='IP',
       #                 exp_id=141398,
       #                 object_name='TEST_MDS'
       #                 )
       d=flap.get_data('NSTX_MDSPlus',
                        name='\EFIT01::\PSIRZ',
                        exp_id=141399,
                        object_name='TEST_MDS'
                        )
       #print(d.coordinate('Time')[0].shape)
       #print(d.coordinate('Device R')[0].shape)
       #print(d.coordinate('Device Z')[0].shape)
       #print(d.coordinate('Dimension 1')[0])
       print(d.data)
       #print(d.data.shape)
    except Exception as e:
        raise e
    flap.plot('TEST_MDS',plot_type='animation',axes=['Device R','Device z','Time'],options={'Z range':[0,0.05],'Wait':0.0,'Clear':False})
    #flap.plot('TEST_MDS',plot_type='contour',axes=['Time','Dimension 1'],options={})
    #flap.plot('TEST_MDS')
    flap.list_data_objects()

def test_mdsplus_efit():

    plt.close('all')
    print("**** Reading all EFIT MDS signals.")
    flap.delete_data_object('*')
    try:
       d=flap_mdsplus.FlapEFITObject()
       d.get_data('NSTX_MDSPlus',exp_id=141399)
       print(d)
    except Exception as e:
        raise e
        
#    flap.plot('TEST_MDS',plot_type='animation',axes=['Device R','Device Z','Time'],options={'Z range':[0,512],'Wait':0.0,'Clear':False})
    flap.list_data_objects()
   
# Reading configuration file in the test directory
thisdir = os.path.dirname(os.path.realpath(__file__))
fn = os.path.join(thisdir,"test_nstx_mdsplus.cfg")
flap.config.read(file_name=fn)

test_mdsplus()
#test_mdsplus_efit()
