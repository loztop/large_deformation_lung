try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

a1447_W_0_0_1_00_es_000 = ExodusIIReader( FileName=['/auto/users/lorenzb/mount_point/plot_data/1447_W_0_0.1_00.e-s.000'] )

AnimationScene1 = GetAnimationScene()
a1447_W_0_0_1_00_es_000.NodeMapArrayStatus = []
a1447_W_0_0_1_00_es_000.FaceVariables = []
a1447_W_0_0_1_00_es_000.ElementVariables = []
a1447_W_0_0_1_00_es_000.XMLFileName = 'Invalid result'
a1447_W_0_0_1_00_es_000.FaceSetResultArrayStatus = []
a1447_W_0_0_1_00_es_000.PointVariables = []
a1447_W_0_0_1_00_es_000.FaceSetArrayStatus = []
a1447_W_0_0_1_00_es_000.FaceMapArrayStatus = []
a1447_W_0_0_1_00_es_000.FileRange = [0, 0]
a1447_W_0_0_1_00_es_000.SideSetResultArrayStatus = []
a1447_W_0_0_1_00_es_000.ElementSetArrayStatus = []
a1447_W_0_0_1_00_es_000.EdgeVariables = []
a1447_W_0_0_1_00_es_000.FilePrefix = '/auto/users/lorenzb/mount_point/plot_data/1447_W_0_0.1_00.e-s.000'
a1447_W_0_0_1_00_es_000.FilePattern = '%s'
a1447_W_0_0_1_00_es_000.EdgeSetArrayStatus = []
a1447_W_0_0_1_00_es_000.SideSetArrayStatus = []
a1447_W_0_0_1_00_es_000.GlobalVariables = []
a1447_W_0_0_1_00_es_000.NodeSetArrayStatus = []
a1447_W_0_0_1_00_es_000.NodeSetResultArrayStatus = []
a1447_W_0_0_1_00_es_000.ElementMapArrayStatus = []
a1447_W_0_0_1_00_es_000.EdgeSetResultArrayStatus = []
a1447_W_0_0_1_00_es_000.ModeShape = 0
a1447_W_0_0_1_00_es_000.EdgeMapArrayStatus = []
a1447_W_0_0_1_00_es_000.ElementSetResultArrayStatus = []

AnimationScene1.EndTime = 10.0
AnimationScene1.PlayMode = 'Snap To TimeSteps'

RenderView1 = GetRenderView()
DataRepresentation10 = Show()
DataRepresentation10.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation10.SelectionPointFieldDataArrayName = 'const'
DataRepresentation10.SelectionCellFieldDataArrayName = 'GlobalElementId'
DataRepresentation10.ScalarOpacityUnitDistance = 15.481535761542942
DataRepresentation10.ExtractedBlockIndex = 2
DataRepresentation10.ScaleFactor = 18.225154876708984

a1447_W_0_0_1_00_es_000.FaceBlocks = []
a1447_W_0_0_1_00_es_000.ElementBlocks = ['Unnamed block ID: 0 Type: TETRA4']
a1447_W_0_0_1_00_es_000.EdgeBlocks = []
a1447_W_0_0_1_00_es_000.ElementVariables = ['s_p', 'p_nu', 'vol_ref', 'J']
a1447_W_0_0_1_00_es_000.PointVariables = ['s_u', 's_v', 's_w', 's_p', 'mono_f_vel_u', 'mono_f_vel_v', 'mono_f_vel_w', 'const', 'u_nu', 'v_nu', 'w_nu', 'p_nu', 'x_nu', 'y_nu', 'z_nu', 'const_', 'u_ref', 'v_ref', 'w_ref', 'vol_ref', 'x_ref', 'y_ref', 'z_ref', 'const__', 'I1', 'I2', 'I3', 'J', 'p1', 'p2', 'p3']

ddd4 = ddd()

RenderView1.CameraClippingRange = [258.673344198543, 886.5910779446125]

DataRepresentation11 = Show()
DataRepresentation11.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation11.SelectionPointFieldDataArrayName = 'const'
DataRepresentation11.SelectionCellFieldDataArrayName = 'GlobalElementId'
DataRepresentation11.ScalarOpacityUnitDistance = 15.481535761542942
DataRepresentation11.ExtractedBlockIndex = 2
DataRepresentation11.ScaleFactor = 18.225154876708984

DataRepresentation10.Visibility = 0

WriteImage('/auto/users/lorenzb/mount_point/lung_pics/test2.png')


AnimationScene1.AnimationTime = 5.0

RenderView1.CacheKey = 5.0
RenderView1.CameraClippingRange = [252.70317025090668, 898.3694421260597]
RenderView1.ViewTime = 5.0
RenderView1.UseCache = 0
RenderView1.Background = [0.319996948195621, 0.34000152590219, 0.429999237048905]
RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]

DataRepresentation10.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation10.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation10.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation11.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation11.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation11.AmbientColor = [1.0, 1.0, 1.0]

Render()
