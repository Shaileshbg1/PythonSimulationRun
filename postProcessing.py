import sys
import numpy as np
from cv2 import fillPoly

from scipy.io import loadmat
sys.path.append("/postProcess_app/paraview/build/lib/python3.5/site-packages/")
sys.path.append("/postProcess_app/paraview/build/lib/python3.5/site-packages/vtkmodules/")
sys.path.append("/postProcess_app/paraview/build/VTK/ThirdParty/vtkm/vtk-m/lib/")
#### import the simple module from the paraview
from paraview.simple import *

# default geometry and number of cells in each direction
xmax = 0.072
xmin = -xmax
zmax = 0.072
zmin = -zmax
m = n = 55

def getPolyMask(pts,height,width):
    im = np.zeros([height,width],dtype=np.uint8)
    pts = np.reshape(pts,(1,pts.shape[0],pts.shape[1]))
    fillPoly(im, pts, 255)
    return im

def getmappedAreolar(mat_path, pts_list, x_areolar_mapped, z_areolar_mapped):
    t = loadmat(mat_path)
    orig_img = t['mappedTemperatureImage']
    pts = np.transpose([np.array(pts_list[0::2]),np.array(pts_list[1::2])])
    ROI=orig_img*(getPolyMask(pts,orig_img.shape[0],orig_img.shape[1])>0)

    # On the mask crop the image as per min max of non zero values
    c1=np.min(np.argwhere(np.sum(ROI,0)>0))
    c2=np.max(np.argwhere(np.sum(ROI,0)>0))
    r1=np.min(np.argwhere(np.sum(ROI,1)>0))
    r2=np.max(np.argwhere(np.sum(ROI,1)>0))
    cropped = ROI[r1:r2, c1:c2]

    # areolar points

    #ArIm=np.zeros(cropped.shape)
    #ArIm[(z_areolar_mapped-r1-1):(z_areolar_mapped-r1+1),(x_areolar_mapped-c1-1):(x_areolar_mapped-c1+1)]=1
    #ArIm = np.flip(ArIm, axis=0)
    #ArIm = np.flip(ArIm, axis=1)
    #xar = round(55*np.mean(np.argwhere(ArIm)[:,1]))/(c2-c1)
    #zar = round(55*np.mean(np.argwhere(ArIm)[:,0]))/(r2-r1)

    xar=np.float32(x_areolar_mapped)-c1  ## changing co-ordinates to fit in the bounding box
    zar=np.float32(z_areolar_mapped)-r1

    xar = n*xar/(c2-c1)
    zar = m*zar/(r2-r1)  ## modifying length to fit to the changed domain

    xar = -((xmax-xmin)/(n)) * (xar) + xmax
    zar = -((zmax-zmin)/(m)) * (zar) + zmax ## linear transformation as per two point formula

    print(xar,zar)

    return xar,zar


def call_paraview(Simulation, x, y, z, r, gr, mat_path, pts_list, x_areolar_mapped, z_areolar_mapped):
    xar, zar = getmappedAreolar(mat_path, pts_list, x_areolar_mapped, z_areolar_mapped)
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'Legacy VTK Reader'
    simulation_ = LegacyVTKReader(FileNames=[Simulation + '/VTK/gland/simulation_0.vtk', Simulation + '/VTK/gland/simulation_10.vtk', Simulation + '/VTK/gland/simulation_20.vtk', Simulation + '/VTK/gland/simulation_30.vtk', Simulation + '/VTK/gland/simulation_40.vtk', Simulation + '/VTK/gland/simulation_50.vtk', Simulation + '/VTK/gland/simulation_60.vtk', Simulation + '/VTK/gland/simulation_70.vtk', Simulation + '/VTK/gland/simulation_80.vtk', Simulation + '/VTK/gland/simulation_90.vtk', Simulation + '/VTK/gland/simulation_100.vtk'])
    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # rename source object
    RenameSource('gland', simulation_)

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [933, 545]

    # show data in view
    simulation_Display = Show(simulation_, renderView1)

    # trace defaults for the display properties.
    simulation_Display.Representation = 'Surface'
    simulation_Display.ColorArrayName = [None, '']
    simulation_Display.OSPRayScaleArray = 'T'
    simulation_Display.OSPRayScaleFunction = 'PiecewiseFunction'
    simulation_Display.SelectOrientationVectors = 'None'
    simulation_Display.ScaleFactor = 0.01440073996782303
    simulation_Display.SelectScaleArray = 'None'
    simulation_Display.GlyphType = 'Arrow'
    simulation_Display.GlyphTableIndexArray = 'None'
    simulation_Display.GaussianRadius = 0.0007200369983911514
    simulation_Display.SetScaleArray = ['POINTS', 'T']
    simulation_Display.ScaleTransferFunction = 'PiecewiseFunction'
    simulation_Display.OpacityArray = ['POINTS', 'T']
    simulation_Display.OpacityTransferFunction = 'PiecewiseFunction'
    simulation_Display.DataAxesGrid = 'GridAxesRepresentation'
    simulation_Display.PolarAxes = 'PolarAxesRepresentation'
    simulation_Display.ScalarOpacityUnitDistance = 0.002408635511542687

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    simulation_Display.ScaleTransferFunction.Points = [29.589799880981445, 0.0, 0.5, 0.0, 38.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    simulation_Display.OpacityTransferFunction.Points = [29.589799880981445, 0.0, 0.5, 0.0, 38.0, 1.0, 0.5, 0.0]

    # reset view to fit data
    renderView1.ResetCamera()

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(simulation_Display, ('POINTS', 'T'))

    # rescale color and/or opacity maps used to include current data range
    simulation_Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    simulation_Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')

    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')

    animationScene1.GoToLast()

    # create a new 'Legacy VTK Reader'
    simulation__1 = LegacyVTKReader(FileNames=[Simulation + '/VTK/tumor/simulation_0.vtk', Simulation + '/VTK/tumor/simulation_10.vtk', Simulation + '/VTK/tumor/simulation_20.vtk', Simulation + '/VTK/tumor/simulation_30.vtk', Simulation + '/VTK/tumor/simulation_40.vtk', Simulation + '/VTK/tumor/simulation_50.vtk', Simulation + '/VTK/tumor/simulation_60.vtk', Simulation + '/VTK/tumor/simulation_70.vtk', Simulation + '/VTK/tumor/simulation_80.vtk', Simulation + '/VTK/tumor/simulation_90.vtk', Simulation + '/VTK/tumor/simulation_100.vtk'])
    # rename source object
    RenameSource('tumor', simulation__1)

    # show data in view
    simulation__1Display = Show(simulation__1, renderView1)

    # trace defaults for the display properties.
    simulation__1Display.Representation = 'Surface'
    simulation__1Display.ColorArrayName = [None, '']
    simulation__1Display.OSPRayScaleArray = 'T'
    simulation__1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    simulation__1Display.SelectOrientationVectors = 'None'
    simulation__1Display.ScaleFactor = 0.0044626730494201185
    simulation__1Display.SelectScaleArray = 'None'
    simulation__1Display.GlyphType = 'Arrow'
    simulation__1Display.GlyphTableIndexArray = 'None'
    simulation__1Display.GaussianRadius = 0.00022313365247100591
    simulation__1Display.SetScaleArray = ['POINTS', 'T']
    simulation__1Display.ScaleTransferFunction = 'PiecewiseFunction'
    simulation__1Display.OpacityArray = ['POINTS', 'T']
    simulation__1Display.OpacityTransferFunction = 'PiecewiseFunction'
    simulation__1Display.DataAxesGrid = 'GridAxesRepresentation'
    simulation__1Display.PolarAxes = 'PolarAxesRepresentation'
    simulation__1Display.ScalarOpacityUnitDistance = 0.002475065553975941

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    simulation__1Display.ScaleTransferFunction.Points = [31.96739959716797, 0.0, 0.5, 0.0, 35.922000885009766, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    simulation__1Display.OpacityTransferFunction.Points = [31.96739959716797, 0.0, 0.5, 0.0, 35.922000885009766, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # set scalar coloring
    ColorBy(simulation__1Display, ('POINTS', 'T'))

    # rescale color and/or opacity maps used to include current data range
    simulation__1Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    simulation__1Display.SetScalarBarVisibility(renderView1, True)

    #### saving camera placements for all active views

    # current camera placement for renderView1
    renderView1.CameraPosition = [0.0, 0.035999998450279236, 0.41729901512298256]
    renderView1.CameraFocalPoint = [0.0, 0.035999998450279236, 0.0]
    renderView1.CameraParallelScale = 0.1080049326163527

    # create a new 'Sphere'
    sphere1 = Sphere()

    # rename source object
    RenameSource('areolar', sphere1)

    # Properties modified on sphere1
    sphere1.Center = [np.int32(xar), 0.07, np.int32(zar)]
    sphere1.Radius = 0.004

    # show data in view
    sphere1Display = Show(sphere1, renderView1)

    # trace defaults for the display properties.
    sphere1Display.Representation = 'Surface'
    sphere1Display.ColorArrayName = [None, '']
    sphere1Display.OSPRayScaleArray = 'Normals'
    sphere1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    sphere1Display.SelectOrientationVectors = 'None'
    sphere1Display.ScaleFactor = 0.000800000037997961
    sphere1Display.SelectScaleArray = 'None'
    sphere1Display.GlyphType = 'Arrow'
    sphere1Display.GlyphTableIndexArray = 'None'
    sphere1Display.GaussianRadius = 4.0000001899898055e-05
    sphere1Display.SetScaleArray = ['POINTS', 'Normals']
    sphere1Display.ScaleTransferFunction = 'PiecewiseFunction'
    sphere1Display.OpacityArray = ['POINTS', 'Normals']
    sphere1Display.OpacityTransferFunction = 'PiecewiseFunction'
    sphere1Display.DataAxesGrid = 'GridAxesRepresentation'
    sphere1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    sphere1Display.ScaleTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    sphere1Display.OpacityTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # Rescale transfer function
    tLUT.RescaleTransferFunction(29.589799880981445, 37.0)

    # Rescale transfer function
    tPWF.RescaleTransferFunction(29.589799880981445, 37.0)

    # change solid color
    sphere1Display.DiffuseColor = [1.0, 0.0, 1.0]

    # set active source
    SetActiveSource(simulation_)

    # create a new 'Slice'
    slice1 = Slice(Input=simulation_)
    slice1.SliceType = 'Plane'
    slice1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1.SliceType.Origin = [0.0, 0.035999998450279236, 0.0]

    # rename source object
    RenameSource('slice@tumorCentre', slice1)

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=slice1.SliceType)

    # Properties modified on slice1.SliceType
    slice1.SliceType.Origin = [np.float32(x), 0.035999998450279236, np.float32(z)]
    slice1.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice1Display = Show(slice1, renderView1)

    # trace defaults for the display properties.
    slice1Display.Representation = 'Surface'
    slice1Display.ColorArrayName = ['POINTS', 'T']
    slice1Display.LookupTable = tLUT
    slice1Display.OSPRayScaleArray = 'T'
    slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1Display.SelectOrientationVectors = 'None'
    slice1Display.ScaleFactor = 0.01412540227174759
    slice1Display.SelectScaleArray = 'None'
    slice1Display.GlyphType = 'Arrow'
    slice1Display.GlyphTableIndexArray = 'None'
    slice1Display.GaussianRadius = 0.0007062701135873795
    slice1Display.SetScaleArray = ['POINTS', 'T']
    slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1Display.OpacityArray = ['POINTS', 'T']
    slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    slice1Display.ScaleTransferFunction.Points = [30.138792037963867, 0.0, 0.5, 0.0, 37.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice1Display.OpacityTransferFunction.Points = [30.138792037963867, 0.0, 0.5, 0.0, 37.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(simulation_, renderView1)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # hide data in view
    Hide(slice1, renderView1)

    # set active source
    SetActiveSource(simulation_)

    # create a new 'Slice'
    slice1_1 = Slice(Input=simulation_)
    slice1_1.SliceType = 'Plane'
    slice1_1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    slice1_1.SliceType.Origin = [0.0, 0.035999998450279236, 0.0]

    # rename source object
    RenameSource('slice@areolarCentre', slice1_1)

    # Properties modified on slice1_1.SliceType
    slice1_1.SliceType.Origin = [np.int32(xar), 0.035999998450279236, np.int32(zar)]
    slice1_1.SliceType.Normal = [0.0, 0.0, 1.0]

    # show data in view
    slice1_1Display = Show(slice1_1, renderView1)

    # trace defaults for the display properties.
    slice1_1Display.Representation = 'Surface'
    slice1_1Display.ColorArrayName = ['POINTS', 'T']
    slice1_1Display.LookupTable = tLUT
    slice1_1Display.OSPRayScaleArray = 'T'
    slice1_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    slice1_1Display.SelectOrientationVectors = 'None'
    slice1_1Display.ScaleFactor = 0.01437838226556778
    slice1_1Display.SelectScaleArray = 'None'
    slice1_1Display.GlyphType = 'Arrow'
    slice1_1Display.GlyphTableIndexArray = 'None'
    slice1_1Display.GaussianRadius = 0.000718919113278389
    slice1_1Display.SetScaleArray = ['POINTS', 'T']
    slice1_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    slice1_1Display.OpacityArray = ['POINTS', 'T']
    slice1_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    slice1_1Display.DataAxesGrid = 'GridAxesRepresentation'
    slice1_1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    slice1_1Display.ScaleTransferFunction.Points = [30.791973114013672, 0.0, 0.5, 0.0, 37.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    slice1_1Display.OpacityTransferFunction.Points = [30.791973114013672, 0.0, 0.5, 0.0, 37.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(simulation_, renderView1)

    # show color bar/color legend
    slice1_1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(simulation_)

    # show data in view
    simulation_Display = Show(simulation_, renderView1)

    # show color bar/color legend
    simulation_Display.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    Hide(slice1_1, renderView1)

    # set active source
    SetActiveSource(slice1_1)

    # show data in view
    slice1_1Display = Show(slice1_1, renderView1)

    # show color bar/color legend
    slice1_1Display.SetScalarBarVisibility(renderView1, True)

    # toggle 3D widget visibility (only when running from the GUI)
    Hide3DWidgets(proxy=slice1_1.SliceType)

    # set active source
    SetActiveSource(slice1)

    # show data in view
    slice1Display = Show(slice1, renderView1)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(simulation_)

    # create a new 'Threshold'
    threshold1 = Threshold(Input=simulation_)
    threshold1.Scalars = ['POINTS', 'T']
    threshold1.ThresholdRange = [29.589799880981445, 37.0]

    # rename source object
    RenameSource('Range_temperature', threshold1)

    # Properties modified on threshold1
    threshold1.ThresholdRange = [32.0, 37.0]

    # show data in view
    threshold1Display = Show(threshold1, renderView1)

    # trace defaults for the display properties.
    threshold1Display.Representation = 'Surface'
    threshold1Display.ColorArrayName = ['POINTS', 'T']
    threshold1Display.LookupTable = tLUT
    threshold1Display.OSPRayScaleArray = 'T'
    threshold1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    threshold1Display.SelectOrientationVectors = 'None'
    threshold1Display.ScaleFactor = 0.014282599836587907
    threshold1Display.SelectScaleArray = 'None'
    threshold1Display.GlyphType = 'Arrow'
    threshold1Display.GlyphTableIndexArray = 'None'
    threshold1Display.GaussianRadius = 0.0007141299918293953
    threshold1Display.SetScaleArray = ['POINTS', 'T']
    threshold1Display.ScaleTransferFunction = 'PiecewiseFunction'
    threshold1Display.OpacityArray = ['POINTS', 'T']
    threshold1Display.OpacityTransferFunction = 'PiecewiseFunction'
    threshold1Display.DataAxesGrid = 'GridAxesRepresentation'
    threshold1Display.PolarAxes = 'PolarAxesRepresentation'
    threshold1Display.ScalarOpacityFunction = tPWF
    threshold1Display.ScalarOpacityUnitDistance = 0.0031277133422768504

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    threshold1Display.ScaleTransferFunction.Points = [32.0, 0.0, 0.5, 0.0, 37.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    threshold1Display.OpacityTransferFunction.Points = [32.0, 0.0, 0.5, 0.0, 37.0, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(simulation_, renderView1)

    # show color bar/color legend
    threshold1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(simulation_)

    # show data in view
    simulation_Display = Show(simulation_, renderView1)

    # show color bar/color legend
    simulation_Display.SetScalarBarVisibility(renderView1, True)

    # create a new 'Contour'
    contour1 = Contour(Input=simulation_)
    contour1.ContourBy = ['POINTS', 'T']
    contour1.Isosurfaces = [33.29489994049072]
    contour1.PointMergeMethod = 'Uniform Binning'

    # rename source object
    RenameSource('Isosurfaces', contour1)

    # show data in view
    contour1Display = Show(contour1, renderView1)

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = ['POINTS', 'T']
    contour1Display.LookupTable = tLUT
    contour1Display.OSPRayScaleArray = 'T'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 0.014400428533554077
    contour1Display.SelectScaleArray = 'T'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'T'
    contour1Display.GaussianRadius = 0.0007200214266777039
    contour1Display.SetScaleArray = ['POINTS', 'T']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'T']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [33.294898986816406, 0.0, 0.5, 0.0, 33.302711486816406, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [33.294898986816406, 0.0, 0.5, 0.0, 33.302711486816406, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(simulation_, renderView1)

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # hide data in view
    Hide(slice1, renderView1)

    # hide data in view
    Hide(slice1_1, renderView1)

    # hide data in view
    Hide(threshold1, renderView1)

    # set active source
    SetActiveSource(threshold1)

    # show data in view
    threshold1Display = Show(threshold1, renderView1)

    # show color bar/color legend
    threshold1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(slice1_1)

    # show data in view
    slice1_1Display = Show(slice1_1, renderView1)

    # show color bar/color legend
    slice1_1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(slice1)

    # show data in view
    slice1Display = Show(slice1, renderView1)

    # show color bar/color legend
    slice1Display.SetScalarBarVisibility(renderView1, True)

    # set active source
    SetActiveSource(simulation_)

    # show data in view
    simulation_Display = Show(simulation_, renderView1)

    # show color bar/color legend
    simulation_Display.SetScalarBarVisibility(renderView1, True)

    # set active view
    SetActiveView(renderView1)

    #### saving camera placements for all active views

    # current camera placement for renderView1
    renderView1.CameraPosition = [0.2590206542859738, 0.23971438808021892, -0.2560219055736826]
    renderView1.CameraFocalPoint = [0.0, 0.035999998450279236, 0.0]
    renderView1.CameraViewUp = [-0.4006009340893404, 0.8701246519142241, 0.28705745372294134]
    renderView1.CameraParallelScale = 0.1080049326163527

    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    SaveState(Simulation + 'stateFile.pvsm')
    print('Saved state')
