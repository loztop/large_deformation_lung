try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

a1447_C_0_5_0_ = LegacyVTKReader( FileNames=['/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_0.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_1.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_2.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_3.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_4.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_5.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_6.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_7.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_8.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_9.vtk', '/auto/users/lorenzb/mount_point/plot_data/1447_C_0.5_0_10.vtk'] )

AnimationScene1 = GetAnimationScene()
AnimationScene1.EndTime = 10.0
AnimationScene1.PlayMode = 'Snap To TimeSteps'

RenderView1 = GetRenderView()
a1_Pressure_PVLookupTable = GetLookupTableForArray( "Pressure", 1 )

DataRepresentation1 = Show()
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation1.SelectionPointFieldDataArrayName = 'Pressure'
DataRepresentation1.SelectionCellFieldDataArrayName = 'Radius'
DataRepresentation1.ColorArrayName = 'Pressure'
DataRepresentation1.LookupTable = a1_Pressure_PVLookupTable
DataRepresentation1.ScaleFactor = 18.267810440063478

RenderView1.CenterOfRotation = [-58.20110422372818, -135.76050186157227, -105.43095207214355]

CellDatatoPointData1 = CellDatatoPointData()

RenderView1.CameraPosition = [-395.3635057096588, 94.48759614155034, 227.64027359409067]
RenderView1.CameraFocalPoint = [-58.20110422372818, -135.76050186157227, -105.43095207214355]
RenderView1.CameraClippingRange = [261.68556021234485, 862.0330601247457]
RenderView1.CameraParallelScale = 136.37310641898483

a1_Radius_PVLookupTable = GetLookupTableForArray( "Radius", 1 )

DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.SelectionPointFieldDataArrayName = 'Radius'
DataRepresentation2.SelectionCellFieldDataArrayName = 'Radius'
DataRepresentation2.ColorArrayName = 'Radius'
DataRepresentation2.LookupTable = a1_Radius_PVLookupTable
DataRepresentation2.ScaleFactor = 18.267810440063478

DataRepresentation1.Visibility = 0

Tube1 = Tube()

Tube1.Scalars = ['POINTS', 'Radius']
Tube1.Vectors = ['POINTS', '']
Tube1.Radius = 1.8267810440063477

DataRepresentation3 = Show()
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation3.SelectionPointFieldDataArrayName = 'Radius'
DataRepresentation3.SelectionCellFieldDataArrayName = 'Radius'
DataRepresentation3.ColorArrayName = 'Radius'
DataRepresentation3.LookupTable = a1_Radius_PVLookupTable
DataRepresentation3.ScaleFactor = 18.559183502197268

DataRepresentation2.Visibility = 0

RenderView1.CameraClippingRange = [256.5265193355142, 868.5721837199143]

Tube1.VaryRadius = 'By Absolute Scalar'

WriteImage('/auto/users/lorenzb/mount_point/lung_pics/tree_pic.png')


RenderView1.CameraClippingRange = [258.70074863117185, 870.9894236564842]
RenderView1.Background = [0.319996948195621, 0.34000152590219, 0.429999237048905]
RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]

DataRepresentation1.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation1.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation2.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation2.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation3.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation3.AmbientColor = [1.0, 1.0, 1.0]

Render()
