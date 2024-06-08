import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from steps.API_2.geom import TetMesh

def PlotTriangles(ax, tris, color):
    ax.add_collection(Poly3DCollection(
        [tri.verts for tri in tris],
        facecolor=color,
        edgecolors='black',
        linewidth=0.1
    ))

def PlotMesh(path):
    mesh = TetMesh.LoadGmsh(path, 1e-6)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')
    PlotTriangles(ax, mesh.surface, (0.5, 0.5, 0.5, 0.1))
    ax.set_xlim(mesh.bbox.min.x, mesh.bbox.max.x)
    ax.set_ylim(mesh.bbox.min.y, mesh.bbox.max.y)
    ax.set_zlim(mesh.bbox.min.z, mesh.bbox.max.z)
    ax.set_xlabel('x position [m]')
    ax.set_ylabel('y position [m]')
    ax.set_zlabel('z position [m]')
    plt.gca().set_aspect('equal')
    plt.show()
