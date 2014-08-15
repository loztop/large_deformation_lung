try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

a5036_C_0_01_0_00_es_000 = ExodusIIReader( FileName=['/auto/users/lorenzb/mount_point/plot_data/5036_C_0.01_0_00.e-s.000'] )

AnimationScene1 = GetAnimationScene()
a5036_C_0_01_0_00_es_000.NodeMapArrayStatus = []
a5036_C_0_01_0_00_es_000.FaceVariables = []
a5036_C_0_01_0_00_es_000.ElementVariables = []
a5036_C_0_01_0_00_es_000.XMLFileName = 'Invalid result'
a5036_C_0_01_0_00_es_000.FaceSetResultArrayStatus = []
a5036_C_0_01_0_00_es_000.PointVariables = []
a5036_C_0_01_0_00_es_000.FaceSetArrayStatus = []
a5036_C_0_01_0_00_es_000.FaceMapArrayStatus = []
a5036_C_0_01_0_00_es_000.FileRange = [0, 0]
a5036_C_0_01_0_00_es_000.SideSetResultArrayStatus = []
a5036_C_0_01_0_00_es_000.ElementSetArrayStatus = []
a5036_C_0_01_0_00_es_000.EdgeVariables = []
a5036_C_0_01_0_00_es_000.FilePrefix = '/auto/users/lorenzb/mount_point/plot_data/5036_C_0.01_0_00.e-s.000'
a5036_C_0_01_0_00_es_000.FilePattern = '%s'
a5036_C_0_01_0_00_es_000.EdgeSetArrayStatus = []
a5036_C_0_01_0_00_es_000.SideSetArrayStatus = []
a5036_C_0_01_0_00_es_000.GlobalVariables = []
a5036_C_0_01_0_00_es_000.NodeSetArrayStatus = []
a5036_C_0_01_0_00_es_000.NodeSetResultArrayStatus = []
a5036_C_0_01_0_00_es_000.ElementMapArrayStatus = []
a5036_C_0_01_0_00_es_000.EdgeSetResultArrayStatus = []
a5036_C_0_01_0_00_es_000.ModeShape = 0
a5036_C_0_01_0_00_es_000.EdgeMapArrayStatus = []
a5036_C_0_01_0_00_es_000.ElementSetResultArrayStatus = []

AnimationScene1.EndTime = 10.0
AnimationScene1.PlayMode = 'Snap To TimeSteps'

RenderView1 = GetRenderView()
DataRepresentation5 = Show()
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation5.SelectionPointFieldDataArrayName = 'const'
DataRepresentation5.SelectionCellFieldDataArrayName = 'GlobalElementId'
DataRepresentation5.ScalarOpacityUnitDistance = 9.762185810774751
DataRepresentation5.ExtractedBlockIndex = 2
DataRepresentation5.ScaleFactor = 19.204430389404298

a5036_C_0_01_0_00_es_000.FaceBlocks = []
a5036_C_0_01_0_00_es_000.ElementBlocks = ['Unnamed block ID: 0 Type: TETRA4']
a5036_C_0_01_0_00_es_000.EdgeBlocks = []
a5036_C_0_01_0_00_es_000.ElementVariables = ['s_p', 'p_nu', 'vol_ref', 'J']
a5036_C_0_01_0_00_es_000.PointVariables = ['s_u', 's_v', 's_w', 's_p', 'mono_f_vel_u', 'mono_f_vel_v', 'mono_f_vel_w', 'const', 'u_nu', 'v_nu', 'w_nu', 'p_nu', 'x_nu', 'y_nu', 'z_nu', 'const_', 'u_ref', 'v_ref', 'w_ref', 'vol_ref', 'x_ref', 'y_ref', 'z_ref', 'const__', 'I1', 'I2', 'I3', 'J', 'p1', 'p2', 'p3']

ddd2 = ddd()

RenderView1.CameraClippingRange = [270.2277957596091, 923.9191626608174]

DataRepresentation6 = Show()
DataRepresentation6.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation6.SelectionPointFieldDataArrayName = 'const'
DataRepresentation6.SelectionCellFieldDataArrayName = 'GlobalElementId'
DataRepresentation6.ScalarOpacityUnitDistance = 9.762185810774751
DataRepresentation6.ExtractedBlockIndex = 2
DataRepresentation6.ScaleFactor = 19.204430389404298

DataRepresentation5.Visibility = 0

AnimationScene1.AnimationTime = 9.0

RenderView1.ViewTime = 9.0
RenderView1.CacheKey = 9.0
RenderView1.UseCache = 0
RenderView1.CameraClippingRange = [251.44671979784474, 955.7124270950221]

Render()
