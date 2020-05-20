#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json, os
import numpy as np
import pandas as pd
import threading
import logging
from simulation import depthAnalysis  

def getDepthInfo(Patient, pts_list, obl_mat_orig, areolar_points, isleft): ## function call to get depth info based on breast 
    obl_mat = obl_mat_orig
    depth=0 ## Initializing depth, serves as a return flag
    
    xlocation=-1 ## Initializing xposition, serves as a return flag
    zlocation=-1 ## Initializing xposition, serves as a return flag
    posloc = {'x':xlocation,'y':zlocation}
    if not os.path.exists(Patient): ## If visit id directory doesnt exist
        print("This when Patient directory doesnt exist", Patient)
        os.mkdir(Patient)
        np.savetxt(Patient + '/prevSegPts.txt', pts_list) ## save passed left seg points for comparison
        logging.info("Main    : before creating thread")
        x = threading.Thread(target=depthAnalysis, args=(Patient, Template, pts_list, obl_mat, areolar_points, isleft, ))  ## create a thread that calls depth analysis 
        logging.info("Main    : before running thread")
        x.start()
        l = threading.enumerate()
        print(threading.enumerate())
        logging.info("Main    : wait for the thread to finish")
        return {'msg': 'Simulation Initiated','depth':depth,'position':posloc, 'response': 0} ## immediate return status immaterial of the thread
  
## global definitions 
dataPath = "/home/shailesh/57_areolar/" 
Template= "/home/shailesh/alphaAnalysis/57_areolar/template/"  ##Template folder for constants

inputs = pd.read_csv("Input.csv")  ## import csv file as a dataframe
Patient = np.array(inputs['Patient'],dtype=np.str) ## reading patient column
    

for subject in range(len(Patient)):   ## each row represents one patient call
    PatientID = Patient[subject].split('/') ## split Niramai ID into MU/Patient/Visit
    GPatientL = dataPath + PatientID[0] + '_' + PatientID[1] + '_' + PatientID[2] + '_' +'left'
    GPatientR = dataPath + PatientID[0] + '_' + PatientID[1] + '_' + PatientID[2] + '_' + 'right'
    pts_list_left = json.loads(inputs['pts_list_left'][subject]) ## load a json array of seg points
    pts_list_right = json.loads(inputs['pts_list_right'][subject]) ## load a json array of seg points
    PTemplate = Template + PatientID[1] + '/'
    obl_mat_left =  PTemplate + inputs['obl_mat_left'][subject] 
    obl_mat_right = PTemplate + inputs['obl_mat_right'][subject]
    areolar_points_left = json.loads(inputs['areolar_points_left'][subject]) ## load a json array of seg points
    areolar_points_right = json.loads(inputs['areolar_points_right'][subject]) ## load a json array of seg points

    msgl=getDepthInfo(GPatientL, pts_list_left, obl_mat_left, areolar_points_left, isleft=1)
    #msgr=getDepthInfo(GPatientR, pts_list_right, obl_mat_right, areolar_points_right, isleft=0)   

