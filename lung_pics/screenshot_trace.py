try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

a5036_W_0_1_00_es_000 = GetActiveSource()
DataRepresentation14 = GetDisplayProperties(a5036_W_0_1_00_es_000)
Slice2 = FindSource("Slice2")
DataRepresentation15 = GetDisplayProperties(Slice2)
a5036_W_0_0_1_00_es_000 = FindSource("5036_W_0_0.1_00.e-s.000")
DataRepresentation17 = GetDisplayProperties(a5036_W_0_0_1_00_es_000)
Slice3 = FindSource("Slice3")
DataRepresentation18 = GetDisplayProperties(Slice3)
a5036_C_0_1_0_00_es_000 = FindSource("5036_C_0.1_0_00.e-s.000")
DataRepresentation20 = GetDisplayProperties(a5036_C_0_1_0_00_es_000)
Slice4 = FindSource("Slice4")
DataRepresentation21 = GetDisplayProperties(Slice4)
ddd2 = FindSource("ddd2")
DataRepresentation24 = GetDisplayProperties(ddd2)
Transform4 = FindSource("Transform4")
DataRepresentation25 = GetDisplayProperties(Transform4)
ddd3 = FindSource("ddd3")
DataRepresentation26 = GetDisplayProperties(ddd3)
Transform5 = FindSource("Transform5")
DataRepresentation27 = GetDisplayProperties(Transform5)
ddd4 = FindSource("ddd4")
DataRepresentation28 = GetDisplayProperties(ddd4)
Transform6 = FindSource("Transform6")
DataRepresentation29 = GetDisplayProperties(Transform6)
a5036_W_0_1_ = FindSource("5036_W_0_1_*")
DataRepresentation31 = GetDisplayProperties(a5036_W_0_1_)
ddd5 = FindSource("ddd5")
DataRepresentation32 = GetDisplayProperties(ddd5)
RenderView1 = GetRenderView()
WriteImage('/auto/users/lorenzb/mount_point/lung_pics/lung_bla.pdf')


DataRepresentation14.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation14.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation14.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation15.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation15.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation15.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation17.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation17.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation17.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation18.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation18.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation18.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation20.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation20.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation20.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation21.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation21.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation21.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation24.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation24.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation24.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation25.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation25.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation25.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation26.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation26.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation26.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation27.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation27.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation27.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation28.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation28.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation28.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation29.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation29.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation29.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation31.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation31.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation31.AmbientColor = [1.0, 1.0, 1.0]

DataRepresentation32.CubeAxesColor = [1.0, 1.0, 1.0]
DataRepresentation32.EdgeColor = [0.0, 0.0, 0.500007629510948]
DataRepresentation32.AmbientColor = [1.0, 1.0, 1.0]

RenderView1.OrientationAxesLabelColor = [1.0, 1.0, 1.0]

Render()
