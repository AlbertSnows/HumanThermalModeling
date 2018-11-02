import vtk

# Visualization
            
def visualize(voxel_db,side_exposed,nx,ny,nz,cx,cy,cz):
    
    # Create a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetWindowName("Tetrahedron")
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    
    colors = vtk.vtkNamedColors() 
    
    points = vtk.vtkPoints()
    
    m = 0
    k = 0
    for c in range(nz):
        for b in range(ny):
            for a in range(nx):
#                if(c > cz and b > cy and a < cx):
                if(side_exposed[a,b,c] > 0):

                    for i in range(2,26):
                        if(voxel_db[k,i].mat != 0):
                            points.InsertNextPoint(voxel_db[k,i].p1.x,voxel_db[k,i].p1.y,voxel_db[k,i].p1.z)
                            points.InsertNextPoint(voxel_db[k,i].p2.x,voxel_db[k,i].p2.y,voxel_db[k,i].p2.z)
                            points.InsertNextPoint(voxel_db[k,i].p3.x,voxel_db[k,i].p3.y,voxel_db[k,i].p3.z)
                            points.InsertNextPoint(voxel_db[k,i].p4.x,voxel_db[k,i].p4.y,voxel_db[k,i].p4.z)
                                
                            unstructuredGrid = vtk.vtkUnstructuredGrid()
                            unstructuredGrid.SetPoints(points)
                        
                            tetra = vtk.vtkTetra()        
                           
                            tetra.GetPointIds().SetId(0, m)
                            tetra.GetPointIds().SetId(1, m+1)
                            tetra.GetPointIds().SetId(2, m+2)
                            tetra.GetPointIds().SetId(3, m+3)
                            m = m + 4
                        
                            cellArray = vtk.vtkCellArray()
                            cellArray.InsertNextCell(tetra)
                            unstructuredGrid.SetCells(vtk.VTK_TETRA, cellArray)
                        
                            # Create a mapper and actor
                            mapper = vtk.vtkDataSetMapper()
                            if vtk.VTK_MAJOR_VERSION <= 5:
                                mapper.SetInputConnection(unstructuredGrid.GetProducerPort())
                            else:
                                mapper.SetInputData(unstructuredGrid)
                            
                                            
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
    renderWindow.Render()
    renderWindowInteractor.Start()


#if __name__ == '__main__':
#    main() 