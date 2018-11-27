import vtk
import time

# Visualization
def visualize(voxel_db, side_exposed, nx, ny, nz):
    # Create a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.SetWindowName("Tetrahedron")
    render_window.AddRenderer(renderer)
    render_window_interactor = vtk.vtkWin32RenderWindowInteractor()
    render_window_interactor.SetRenderWindow(render_window)
    colors = vtk.vtkNamedColors()

    points = vtk.vtkPoints()

    m = 0
    k = 0
    for c in range(nz):
        for b in range(ny):
            for a in range(nx):
                #                if(c > cz and b > cy and a < cx):
                if side_exposed[a, b, c] > 0:

                    for i in range(2, 26):
                        if voxel_db[k, i].mat != 0:
                            # I made changes to this next bit
                            v_db = [voxel_db[k, i].p1.x, voxel_db[k, i].p1.y, voxel_db[k, i].p1.z]
                            points.InsertNextPoint(v_db)
                            v_db = [voxel_db[k, i].p2.x, voxel_db[k, i].p2.y, voxel_db[k, i].p2.z]
                            points.InsertNextPoint(v_db)
                            v_db = [voxel_db[k, i].p3.x, voxel_db[k, i].p3.y, voxel_db[k, i].p3.z]
                            points.InsertNextPoint(v_db)
                            v_db = [voxel_db[k, i].p4.x, voxel_db[k, i].p4.y, voxel_db[k, i].p4.z]
                            points.InsertNextPoint(v_db)

                            unstructured_grid = vtk.vtkUnstructuredGrid()
                            unstructured_grid.SetPoints(points)

                            tetra = vtk.vtkTetra()

                            tetra.GetPointIds().SetId(0, m)
                            tetra.GetPointIds().SetId(1, m + 1)
                            tetra.GetPointIds().SetId(2, m + 2)
                            tetra.GetPointIds().SetId(3, m + 3)
                            m = m + 4

                            cell_array = vtk.vtkCellArray()
                            cell_array.InsertNextCell(tetra)
                            unstructured_grid.SetCells(vtk.VTK_TETRA, cell_array)

                            # Create a mapper and actor
                            mapper = vtk.vtkDataSetMapper()
                            if vtk.VTK_MAJOR_VERSION <= 5:
                                mapper.SetInputConnection(unstructured_grid.GetProducerPort())
                            else:
                                mapper.SetInputData(unstructured_grid)

                            actor = vtk.vtkActor()
                            actor.SetMapper(mapper)
                            actor.GetProperty().SetColor(colors.GetColor3d("pink"))

                            # Add the actor to the scene
                            renderer.AddActor(actor)
                k = k + 1

    renderer.SetBackground(colors.GetColor3d("white"))
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(-10)
    renderer.GetActiveCamera().Elevation(-20)

    # Render and interact
    render_window.Render()
    time.sleep(10)
    render_window_interactor.Start()
# if __name__ == '__main__':
#    main()
