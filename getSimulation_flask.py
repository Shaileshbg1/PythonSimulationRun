# flask code which receives inputs via an API,invokes the simulations and returns the depth value
from flask import Flask, request

from flask_restful import Resource, Api
import requests, json, os, shutil 
import numpy as np
import threading
import logging
from io import BytesIO
import base64

from dataGeneration import depthAnalysis  ## main function

app = Flask(__name__)
api = Api(app)

## global definitions 
dataPath = "/simulation_build/solver_data/" ## path as per docker containers
Template= "/simulation_build/template/"  ##Template folder for constants

def getDepthInfo(Patient, pts_list, obl_mat_orig): ## function call to get depth info based on breast 
    obl_mat = BytesIO(base64.b64decode(obl_mat_orig))
    depth=0 ## Initializing depth, serves as a return flag
    
    xlocation=-1 ## Initializing xposition, serves as a return flag
    zlocation=-1 ## Initializing xposition, serves as a return flag
    posloc = {'x':xlocation,'y':zlocation}
  
    if not os.path.exists(Patient): ## If visit id directory doesnt exist
        print("This when Patient directory doesnt exist", Patient)
        os.mkdir(Patient)
        np.savetxt(Patient + '/prevSegPts.txt', pts_list) ## save passed left seg points for comparison
        logging.info("Main    : before creating thread")
        x = threading.Thread(target=depthAnalysis, args=(Patient, Template, pts_list, obl_mat,)) ## create a thread that calls depth analysis 
        logging.info("Main    : before running thread")
        x.start()
        logging.info("Main    : wait for the thread to finish")
        return {'msg': 'Simulation Initiated','depth':depth,'position':posloc, 'response': 0} ## immediate return status immaterial of the thread
        ## upon first call generates simulation folders and success/error msg
        ## on second call skips generation checks if texts exists and chooses appropriate actions    
        ## seg points are saved for every first run, and checked for changes for every future run
           
    elif os.path.exists(Patient + '/error.txt'):  ## if error txt exists - remove and re-run the simulation
        print("This is when error exists", Patient)

        ptsListLoaded  = np.loadtxt(Patient + '/prevSegPts.txt')
        if (np.sum(np.abs((pts_list) - ptsListLoaded)))== 0:
            return {'msg': 'Error in simulation','depth':depth,'position':posloc,'response': 2}
        else:
            shutil.rmtree(Patient) ## remove the patient directory
            os.mkdir(Patient) ## create a new directory and respawn the thread
            np.savetxt(Patient + '/prevSegPts.txt', pts_list) ## save the seg points passed as old directory deleted
            logging.info("Main    : before creating thread")
            x = threading.Thread(target=depthAnalysis, args=(Patient, Template, pts_list, obl_mat,)) ## new thread for depth analysis
            logging.info("Main    : before running thread")
            x.start()
            logging.info("Main    : wait for the thread to finish")
            return {'msg': 'Segmentation points changed , simulation respun','depth':depth,'position':posloc,'response': 0}
    
        
    ## If the simulation is running without any error
    elif not os.path.exists(Patient + '/success.txt') and not os.path.exists(Patient + '/depth.txt') :
        print("This is when success and depth both do not exist", Patient)
        return {'msg': 'Simulation running','depth':depth,'position':posloc,'response': 0 } ## if success and depth txt doesnt exist
    
    ## If patient directory exists and success text exists, generation is done call analyzer
    elif os.path.exists(Patient + '/success.txt') and not os.path.exists(Patient + '/depth.txt') : ## only loop that passes to analyzer after simulation generation
        print("This is first depth call, removes VTK", Patient)
        ## check if seg points have been changed before passing to analysis
        ptsListLoaded  = np.loadtxt(Patient + '/prevSegPts.txt') ## load segmentation points for comparison         
        if np.sum(np.abs(pts_list - ptsListLoaded))== 0: ## sum of absolute difference
            ## no change in seg points call analysis
            data = {'Patient':Patient, 'Template':Template, 'pts_list':pts_list, 'obl_mat':obl_mat_orig}
            r = requests.post('http://analyzer:9090/analyzer', json=data) ## use requests to send http request to analyzer docker, json-convert points list into json array  
            print(r.status_code, Patient)
            
            try:
                shutil.rmtree(Patient + '/simulation', ignore_errors=True) ## remove the simulation directory
                shutil.rmtree(Patient + '/VTK', ignore_errors=True) ## remove VTK directory to save space
                print("Removed the simulation and VTK directories", Patient)
            except:
                print("Could not remove the simulation and VTK directories", Patient)
                
            depth, xlocation, zlocation,depth_msg = np.loadtxt(Patient + '/depth.txt',dtype=np.str,delimiter='/n') ## save depth value and x and z coords of hottest point
            if np.float32(depth) < 0:
                
                return {'msg': 'Error in simulation','depth':depth,'position':posloc,'response': 2}
            else:
                return {'msg': depth_msg,'depth':depth,'position':posloc, 'response': 1  }
            
            ## remove simulation/VTKfiles

    
        else:
            print("This is when new set of seg points are passed", Patient)
            ## new set of seg points passed recalculate
            shutil.rmtree(Patient) ## remove the patient directory
            os.mkdir(Patient) ## create a new directory and respawn the thread
            np.savetxt(Patient + '/prevSegPts.txt', pts_list) ## save the newly passed seg points for future checks
            logging.info("Main    : before creating thread")
            x = threading.Thread(target=depthAnalysis, args=(Patient, Template, pts_list, obl_mat,)) ## new thread as seg points have been changed
            logging.info("Main    : before running thread")
            x.start()
            logging.info("Main    : wait for the thread to finish")   
            return {'msg': 'Segmentation points changed , simulation respun','depth':depth,'position':posloc,'response': 0}

    elif os.path.exists(Patient + '/depth.txt'): ## full forward flow has run before
        print("this is when depth exists", Patient)
        ptsListLoaded  = np.loadtxt(Patient + '/prevSegPts.txt') ## load seg points for checking 
        print("Checking pts list done", Patient)
        ## check if there are any changes in seg points passed
        if np.sum(np.abs(pts_list - ptsListLoaded))== 0: ##sum of absolute difference
           ## no change in seg points return previously calculated values
            print("No change in points list", Patient)
            depth, xlocation, zlocation,depth_msg = np.loadtxt(Patient + '/depth.txt',dtype=np.str,delimiter='/n')
            print("loading text done", Patient)
            if np.float32(depth) < 0:
                print("this is when depth = -20", depth, Patient)
                print("This is when there is an error in simulation", Patient)
                return {'msg': 'Error in simulation','depth':depth,'position':posloc,'response': 2}
            else:
                print("This is when there is no error in simulation", Patient)
                return {'msg': depth_msg,'depth':depth,'position':posloc, 'response': 1  }

            
        else:
            ## new set of seg points passed recalculate
            shutil.rmtree(Patient) ## remove the patient directory
            os.mkdir(Patient) ## create a new directory and respawn the thread
            np.savetxt(Patient + '/prevSegPts.txt', pts_list) ## save the newly passed seg points for future checks
            logging.info("Main    : before creating thread")
            x = threading.Thread(target=depthAnalysis, args=(Patient, Template, pts_list, obl_mat,)) ## new thread as seg points have been changed
            logging.info("Main    : before running thread")
            x.start()
            logging.info("Main    : wait for the thread to finish")
            return {'msg': 'Segmentation points changed , simulation respun','depth':depth,'position':posloc,'response': 0}
            
class HelloWorld(Resource):   ## Get request for checking flask/apache/docker setup
    def get(self):
        return {'about': 'Hello World!'}
            
class Solver(Resource):    ## Post request for solver endpoint
    def post(self):   
        req_data = request.get_json(force=True)
        
        try:
            GPatientID = req_data['NiramaiID']   ## Requesting Niramai ID as saving separately for future reference if need be
            pts_list_left = json.loads(req_data['pts_list_left']) ## Requesting left obl breast segmentation points, json-convert points list into json array
            pts_list_right = json.loads(req_data['pts_list_right']) ## Requesting right obl breast segmentation points, json-convert points list into json array
            
              
            PatientID = GPatientID.split('/') ## split Niramai ID into MU/Patient/Visit
            GPatientL = dataPath + PatientID[0] + '_' + PatientID[1] + '_' + PatientID[2] + '_' +'left'
            GPatientR = dataPath + PatientID[0] + '_' + PatientID[1] + '_' + PatientID[2] + '_' + 'right'
            
            msgl=getDepthInfo(GPatientL, pts_list_left, req_data['obl_mat_left'])
            msgr=getDepthInfo(GPatientR, pts_list_right, req_data['obl_mat_right'])   
        except:
            msgl= {'msg': 'Error in simulation','depth':-20,'position':{'x':-1,'y':-1},'response': 2}
            msgr= {'msg': 'Error in simulation','depth':-20,'position':{'x':-1,'y':-1},'response': 2}
        return {'left':msgl, 'right': msgr}  

api.add_resource(HelloWorld, '/')
api.add_resource(Solver, '/solver')

if __name__ == '__main__':
    app.run(host="0.0.0.0", port=8080, debug=True)
    
