try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

AnimationScene1 = GetAnimationScene()
RenderView1 = CreateRenderView()
RenderView1.CameraViewUp = [0.13874578729053702, 0.018580777807876803, 0.9901537058483311]
RenderView1.OrientationAxesLabelColor = [0.0, 0.0, 0.0]
RenderView1.CameraClippingRange = [312.5118357397997, 804.7612595820765]
RenderView1.Background = [1.0, 1.0, 1.0]
RenderView1.UseCache = 0

LegacyVTKReader2 = LegacyVTKReader()

CellDatatoPointData2 = CellDatatoPointData()

Tube2 = Tube()

TimeAnimationCue1 = GetTimeTrack()
TimeAnimationCue1.AnimatedDomainName = ''
TimeAnimationCue1.KeyFrames = []

PiecewiseFunction3 = CreatePiecewiseFunction()

a1_Pressure_PVLookupTable = GetLookupTableForArray( "Pressure", 1 )

PiecewiseFunction4 = CreatePiecewiseFunction()

a1_Radius_PVLookupTable = GetLookupTableForArray( "Radius", 1, RGBPoints=[0.024000000208616257, 0.23, 0.299, 0.754, 7.62936019897461, 0.706, 0.016, 0.15] )

a1_Pressure_PiecewiseFunction = CreatePiecewiseFunction()

a1_Radius_PiecewiseFunction = CreatePiecewiseFunction()

PiecewiseFunction1 = CreatePiecewiseFunction()

PiecewiseFunction2 = CreatePiecewiseFunction()

my_representation3 = Show()

my_representation4 = Show()

my_representation5 = Show()

AnimationScene1.PlayMode = 'Snap To TimeSteps'
AnimationScene1.ViewModules = RenderView1
AnimationScene1.EndTime = 10.0
AnimationScene1.TimeKeeper = []

WriteImage('/auto/users/lorenzb/mount_point/lung_pics/tree.png')


RenderView1.Background = [0.319996948195621, 0.34000152590219, 0.429999237048905]
RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]

Render()
