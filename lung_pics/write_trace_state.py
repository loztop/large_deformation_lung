try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

AnimationScene1 = GetAnimationScene()
RenderView1 = CreateRenderView()
RenderView1.CameraViewUp = [0.06918088592473422, 0.41765841130432196, 0.9059665868504306]
RenderView1.CameraPosition = [-485.2393502939059, -189.51658139104825, -47.79687360738539]
RenderView1.CameraClippingRange = [262.2335556733704, 657.6852084190525]
RenderView1.UseCache = 0
RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]

a5036_C_1_0_00_es_000 = ExodusIIReader()

Slice1 = Slice()

a5036_C_0_25_0_00_es_000 = ExodusIIReader()

a5036_C_0_25_0_ = LegacyVTKReader()

CellDatatoPointData1 = CellDatatoPointData()

Tube1 = Tube()

Sphere1 = Sphere()

Sphere2 = Sphere()

a1_s_p_PiecewiseFunction = CreatePiecewiseFunction()

a1_s_p_PVLookupTable = GetLookupTableForArray( "s_p", 1 )

ScalarBarWidgetRepresentation1 = CreateScalarBar( TitleColor=[0.0, 0.0, 0.0], LabelColor=[0.0, 0.0, 0.0] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

a1_Pressure_PiecewiseFunction = CreatePiecewiseFunction()

a1_Pressure_PVLookupTable = GetLookupTableForArray( "Pressure", 1 )

a1_Radius_PiecewiseFunction = CreatePiecewiseFunction()

a1_Radius_PVLookupTable = GetLookupTableForArray( "Radius", 1, RGBPoints=[0.05999999865889549, 0.23, 0.299, 0.754, 7.62936019897461, 0.706, 0.016, 0.15] )

TimeAnimationCue1 = GetTimeTrack()
TimeAnimationCue1.KeyFrames = []

DataRepresentation2 = Show()
DataRepresentation2.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5]

DataRepresentation4 = Show()
DataRepresentation4.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5]
DataRepresentation4.AmbientColor = [0.0, 0.0, 0.0]

DataRepresentation6 = Show()
DataRepresentation6.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation6.EdgeColor = [0.0, 0.0, 0.5]
DataRepresentation6.AmbientColor = [0.0, 0.0, 0.0]

DataRepresentation7 = Show()
DataRepresentation7.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation7.EdgeColor = [0.0, 0.0, 0.5]
DataRepresentation7.AmbientColor = [0.0, 0.0, 0.0]

DataRepresentation8 = Show()
DataRepresentation8.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation8.EdgeColor = [0.0, 0.0, 0.5]
DataRepresentation8.AmbientColor = [0.0, 0.0, 0.0]

DataRepresentation9 = Show()
DataRepresentation9.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation9.EdgeColor = [0.0, 0.0, 0.5]
DataRepresentation9.AmbientColor = [0.0, 0.0, 0.0]

DataRepresentation10 = Show()
DataRepresentation10.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation10.EdgeColor = [0.0, 0.0, 0.5]
DataRepresentation10.AmbientColor = [0.0, 0.0, 0.0]

DataRepresentation11 = Show()
DataRepresentation11.CubeAxesColor = [0.0, 0.0, 0.0]
DataRepresentation11.EdgeColor = [0.0, 0.0, 0.5]

AnimationScene1.AnimationTime = 9.0
AnimationScene1.ViewModules = RenderView1
AnimationScene1.TimeKeeper = []
AnimationScene1.PlayMode = 'Snap To TimeSteps'
AnimationScene1.EndTime = 10.0

WriteImage('/auto/users/lorenzb/mount_point/lung_pics/load_pic.png')


DataRepresentation2.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.500007629510948]

DataRepresentation4.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation4.AmbientColor = [1.0, 1.0, 1.0]

ScalarBarWidgetRepresentation1.TitleColor = [1.0, 1.0, 1.0]
ScalarBarWidgetRepresentation1.LabelColor = [1.0, 1.0, 1.0]

DataRepresentation6.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation6.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation6.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation7.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation7.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation7.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation8.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation8.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation8.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation9.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation9.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation9.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation10.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation10.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation10.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation11.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation11.EdgeColor = [0.0, 0.0, 0.500007629510948]

RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]

Render()
