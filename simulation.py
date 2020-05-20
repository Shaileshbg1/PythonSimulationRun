## code to impose actual skin temperature as a boundary, iteratively run simulations for diffferent positions of the tumor and generate VTKs for analysis
import os, shutil,traceback
import numpy as np
from PyFoam.Basics.TemplateFile import TemplateFile ## for template conversions
from PyFoam.Execution.BasicRunner import BasicRunner ## run openfoam commands
from PyFoam.Error import error ## output error msg
from PyFoam.Applications.ClearCase import ClearCase ## clear case
from PyFoam.Applications.CloneCase import CloneCase ## make a copy 

## Constants throughout the simulation- defined as global variables
## default geometry and aspect ratio in each direction
xmax = 0.072
xmin = -xmax
zmax = 0.072
zmin = -zmax
m = n = 166
y_var = [0.030, 0.032, 0.034, 0.036, 0.038, 0.040]
xtr = 0
ztr = 0

## Inputs: Tumor radius, breast radius and state of malignancy
radius = 0.005
gr = 0.072
state = "malignant"

## function to run openFOAM model
def FOAM_model(xtr, y, ztr, Patient, Template): 
    errorCode = True ## True when simulation has succesfully run
    ## Specify simulation folder
    Simulation = Patient + "/simulation"

    ## Clear case
    ClearCase(args=["--clear-history", Simulation])
    print("Complete cleaning the case done")

    if not os.path.exists(Simulation):  ## if simulation directory doesnt exist
        ## Clone template onto simulation folder
        CloneCase(args=[Template , Simulation])
        print("Copied generic case to patient specific folder")
        
    ## copy betavSolid and actual skin temperature data onto the gland 0 folder    
    shutil.copyfile(Template +"/0/gland/betavSolid", Simulation +"/0/gland/betavSolid")
    shutil.copyfile(Patient +"/actualSkinData", Simulation +"/0/gland/actualSkinData")

    ## define different cell zones using topoSetDict
    bmName = os.path.join(Simulation,'system', "topoSetDict")
    template = TemplateFile(bmName + ".template", expressionDelimiter="$")
    template.writeToFile(bmName,
                         {'x': xtr, 'y': y, 'z': ztr, 'r': radius, 'gr': gr})

    print("Setting template file for topoSet done")

    ## Run topoSet
    topoSetRun = BasicRunner(argv=["topoSet", "-case" , Simulation], silent=True, server=False, logname='log.topoSet')
    topoSetRun.start()
    if not topoSetRun.runOK():
        error("There was a problem with topoSet")
    print("topoSet done")
    print(xtr,y,ztr)

    ## Split mesh regions based on toposet
    splitMeshRegionsRun = BasicRunner(argv=["splitMeshRegions -cellZones -overwrite", "-case" , Simulation], silent=True, server=False, logname='log.splitMeshRegions')
    splitMeshRegionsRun.start()
    if not splitMeshRegionsRun.runOK():
        error("There was a problem with split mesh regions")
    print("split mesh regions done")

    ## Run change dictionary for gland region
    changeDictionaryGlandRun = BasicRunner(argv=[" changeDictionary -region gland", "-case" , Simulation], silent=True, server=False, logname='log.changeDictionaryGland')
    changeDictionaryGlandRun.start()
    if not changeDictionaryGlandRun.runOK():
        error("There was a problem with change dictionary for gland")
    print("change dictionary gland done")

    ## Run change dictionary for tumor region
    changeDictionaryTumorRun = BasicRunner(argv=[" changeDictionary -region tumor", "-case" , Simulation], silent=True, server=False, logname='log.changeDictionaryTumor')
    changeDictionaryTumorRun.start()
    if not changeDictionaryTumorRun.runOK():
        error("There was a problem with change dictionary for tumor")
    print("change dictionary tumor done")

    ## Run setFields for gland region
    setFieldsGlandRun = BasicRunner(argv=["setFields -region gland", "-case" , Simulation], silent=True, server=False, logname='log.setFieldsGland')
    setFieldsGlandRun.start()
    if not setFieldsGlandRun.runOK():
        error("There was a problem with setFields for gland")
    print("set fields for gland done")

    ## define gland anisotropic thermal conductivity
    bmName = os.path.join(Simulation,'constant','gland', "thermophysicalProperties")
    template = TemplateFile(bmName + ".template", expressionDelimiter="$")
    template.writeToFile(bmName,
                         {'x': xtr, 'y': y, 'z': ztr})

    print("Setting anisotropic thermal conductivity for gland done")

    ## define tumor anisotropic thermal conductivity
    bmName = os.path.join(Simulation,'constant','tumor', "thermophysicalProperties")
    template = TemplateFile(bmName + ".template", expressionDelimiter="$")
    template.writeToFile(bmName,
                         {'x': xtr, 'y': y, 'z': ztr})

    print("Setting anisotropic thermal conductivity for tumor done")

    ## removing fvoptions if benign tumor
    if state == 'benign':
        if not os.path.exists(Simulation + "/constant/tumor/fvOptions"):
            print("Removing heat sources for benign tumors done")
        else:
            os.remove(Simulation + "/constant/tumor/fvOptions")
            print("Removing heat sources for benign tumors done")

    ## multi region simple foam with two heat sources specified for tumor region
    print ("Running")
    theRun = BasicRunner(argv=["chtMultiRegionSimpleFoam", "-case",  Simulation], silent=True, server=False, logname='log.solver')
    #"-postProcess", "-func", "surfaces"
    theRun.start()
    errorCode = theRun.endSeen

    if not theRun.runOK():
        error("There was a problem while running the solver")
    print("Solver run done")
    
    ## converting latest simulation step to VTK- gland
    print ("Converting gland region to VTK")
    VTKGlandRun = BasicRunner(argv=["foamToVTK -fields '(T)' -latestTime -ascii -region gland", "-case", Simulation], silent=True, server=False, logname='log.VTKGland')
    VTKGlandRun.start()
    if not VTKGlandRun.runOK():
        error("There was a problem while converting the gland region to VTK for post-processing")
    print("Conversion of Gland region to VTK done")
    
    ## converting latest simulation step to VTK- tumor
    print ("Converting tumor region to VTK")
    VTKTumorRun = BasicRunner(argv=["foamToVTK -fields '(T)' -latestTime -ascii -region tumor", "-case", Simulation], silent=True, server=False, logname='log.VTKTumor')
    VTKTumorRun.start()
    if not VTKTumorRun.runOK():
        error("There was a problem while converting the tumor region to VTK for post-processing")
    print("Conversion of Tumor region to VTK done")
    
    ## Moving VTK for post processing by rounding off to two decimal places 
    if ((y*100)%1)==0:
        y_str=str(round(y*100)/100)+'0'
    else:
        y_str=str(y)
    shutil.move(Simulation + "/VTK", Patient + "/VTK" + "/VTK" + y_str)
    return errorCode

##Function to define the y-variation and subsequent calling of simulation for each y
def depthAnalysis(Patient,Template, pts_list, obl_mat): ## main function that is called from flaskAPI
    ## Try except for error handling
    #print("Check if thread is spawned")
    try:
      
        ## for every y in the tumor loop
        errorCode = True
        for idx_y,y in enumerate(y_var): 
            ## Function3: OpenFOAM simulation for temeratures calculation         
            errorCode = FOAM_model(xtr, y, ztr, Patient, Template)
            if not errorCode:
                break
        if errorCode:
            print(str(Patient) +"Data generated ready for VTK post-processing")
            ## save success text for successful simulation
            np.savetxt(Patient + '/success.txt', np.array([0]))
        else:
            shutil.rmtree(Patient)
            
        return {'msg': str(Patient) + 'Simulation completed'}
    ## If run executed successfully create a success text, return msg
    except:
        ## save error text for erreneous simulation
        traceback.print_exc(file=open(Patient + "/error.txt", "a"))
        
        return {'msg': str(Patient) + 'Simulation error'}
    ## If any error encountered create an error text and return msg
