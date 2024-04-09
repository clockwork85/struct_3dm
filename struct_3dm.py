#!/usr/bin/env python 

import sys, os, subprocess


class CreateStructuredMesh:
    ''' Common base class for creating structured 3dms with variable meshing'''

    def  __init__(self, sim_name='temp', dz=3.5, top_soil=1.5, res=0.01,
            width=0.02, xy_ratio=2, max_aspect=10):
        self.sim_name = sim_name
        self.dz = dz
        self.top_soil = top_soil
        self.res = res
        self.width = width
        self.xy_ratio = xy_ratio
        self.max_aspect = max_aspect

    def determine_max_aspect_ratio(self):
        ''' This algorithm determines the spacing for the adaptive mesh that gives
        you a number that is very close to the maximum aspect ratio given at the 
        bottom of the mesh - There is probably a much better way to do this'''
        height = self.dz - self.top_soil
        final_ar = 1000000 
        ntry = 1.5  # This is just arbitrarily large
        while True:
            ratios=[]
            nres = self.res
            ch = height
            lh = ch
            while True:
                nres *= ntry
                ar = lh/self.res
                lh = ch
                ch -= nres
                # This means we have overshot and are under the mesh
                if ch < 0: 
                    ratios.append(ch)
                    final_ar = ar  
                    ntry -= 0.0001
                    break
                else:
                    ratios.append(ch)
            if ntry <= 0:
                print("ERROR: final ar is:", final_ar)
            elif ratios[-3]/self.res < self.max_aspect: 
                break
            else:
                pass
        ratios = ratios[:-2]
        ratios.append(0)
        return ratios, final_ar

    def create_struct_mesh(self):
            #sim_name, dz, top_soil, res, width, xy_ratio, max_aspect):
        ''' This takes a list of parameters and makes a structured mesh from those
        parameters.  You must have a working version of tetgen in your path. ''' 

        start_of_soil = self.dz - self.top_soil 
        layers_hr_z = int(self.top_soil/self.res)
        #layers_hr_z = int(dz/res) # for single resolution
        
#=============================================================================
#     Creating the top nodes to be extruded
#=============================================================================    
        nnodesx = int(self.width/self.res + 1)

        nodes = []
        x_count = 0 
        for i in range(nnodesx):
            nodes.append((x_count, 0, self.dz))
            x_count += self.res*self.xy_ratio

        top_nodes = []
        y_count = 0  
        for i in range(nnodesx):
            top_nodes += [(x,y+y_count,z) for x,y,z in nodes]
            y_count += self.res*self.xy_ratio
      
#=============================================================================
# Putting together the node numbers of the top boundary nodes for bcs file
#=============================================================================
        top_bnd_nodes = [x+1 for x in range(len(top_nodes))]

#=============================================================================
#     Extruding the nodes in the -z direction
#=============================================================================
        nodes = list(top_nodes)
        z_count = self.res
        for i in range(layers_hr_z):
            nodes += [(x,y,z-z_count) for x,y,z in top_nodes]
            z_count += self.res 
       
        if self.dz == self.top_soil: 
            print("Creating static resolution mesh")
            points = 0
            final_ar = 0.0
            len_points = 0
        else:
            points, final_ar = self.determine_max_aspect_ratio()
            len_points = len(points)

                
#=============================================================================
# Printing out mesh statistics
#=============================================================================
        print("Total height of mesh: ", self.dz)
        print("Start of soil layer: ", start_of_soil)
        print("Resolution in Z: ", self.res)
        print("XY aspect ratio: ", self.xy_ratio)
        print("Max aspect ratio: ", round(final_ar,4)) 
        print(layers_hr_z, "layers of high resolution")
        print(f"{len_points} layers of variable resolution")
        print(f"Points: {points}")

#=============================================================================
#     Creating the nodes that are at the top of the soil layer, then using the
#     points generated to create the bottom adaptive mesh
#=============================================================================
        bottom_nodes = [(x,y,z-self.top_soil) for x,y,z in top_nodes]
        for i in range(len_points):
            nodes += [(x,y,z-(start_of_soil-points[i])) for x,y,z in bottom_nodes]

        boundary_nodes = [x for x in range(len(nodes)-nnodesx**2+1,len(nodes)+1)]
#=============================================================================
#     Printing out the nodes to a particular .node file for tetgen to
#     tetrahedralize
#=============================================================================
        cwd = os.getcwd()
        os.chdir(cwd)
        f1 = open(self.sim_name+'.node','w')
        f1.write('%i 3 0 1\n'% len(nodes))
        cnt = 1
        for x,y,z in nodes:
            if z != self.dz:
                f1.write('%i %f %f %f\n' % (cnt,x,y,z))
                cnt += 1
            else:
                f1.write('%i %f %f %f 1\n' % (cnt,x,y,z))
                cnt += 1
        f1.close()

#=============================================================================
# Externally calling tetgen to create the element file to create the 3dm
#=============================================================================
        print("Running tetgen")
        cmd = 'tetgen '+self.sim_name+'.node'
        failure, output = subprocess.getstatusoutput(cmd)

#=============================================================================
# Creating the 3dm from the .node and .ele files that are produced from tetgen, 
# I changed this from calling an external program.
#=============================================================================
        print("Creating 3dm")
        #cmd = 'ctet3dm.py ' + sim_name
        outfile = open(self.sim_name+'.3dm','w')
        outfile.write('MESH3D\n')
        with open(self.sim_name+'.1.ele') as f1:
            f1.readline()
            for line in f1.readlines():
                try:
                    int(line.split()[0])    
                    line = line.split(); line[-1] = line[-1].rstrip('\n') 
                    line = [int(x) for x in line]
                    nn,n1,n2,n3,n4 = line
                    outfile.write('E4T %d %d %d %d %d 1\n' % (nn,n1,n2,n3,n4))
                except:
                    pass
        with open(self.sim_name+'.1.node') as f1:
            f1.readline()
            for line in f1.readlines():
                try:
                    int(line.split()[0])
                    nn,x,y,z,j = line.split()
                    nn=int(nn); x=float(x); y=float(y); z=float(z)
                    outfile.write('ND %d %.5g %.5g %.5g\n' % (nn,x,y,z)) 
                except:
                    pass
        outfile.close()

#=============================================================================
# Finding and printing out boundary conditions
#=============================================================================
        print("Creating boundary nodes and facets")
        outfile = open(self.sim_name+'.bcs', 'w')
        i = 0 
        with open(self.sim_name+'.3dm') as f1:
            for line in f1:
                if line.split()[0] == "E4T":
                    element = line.split()[2:6]
                    element = [int(e) for e in element]
                    not_in_face = [e for e in element if e not in top_bnd_nodes]
                    if len(not_in_face) == 1:
                        not_in_face = int(not_in_face[0])
                        index = element.index(not_in_face)
#---------------------------------------------------------                         
#  Different cases for FCS cards                                                   
#  Case 1: Nodes are 2, 3, 4                                                       
#  Case 2: Nodes are 1, 4, 3                                                       
#  Case 3: Nodes are 1, 2, 4                                                       
#  Case 4: Nodes are 1, 3, 2                                                       
#---------------------------------------------------------                         
                        if index is 0:                                                     
                            vert_not_in_tet = 1                                            
                        elif index is 1:                                                   
                            vert_not_in_tet = 2                                            
                        elif index is 2:                                                   
                            vert_not_in_tet = 3                                            
                        elif index is 3:                                                   
                            vert_not_in_tet = 4                                            
                        else:                                                              
                            print("There is some problem with your index\n")                
                            sys.exit()                                                     
                        line = "FCS %s %d 2\n" % (line.split()[1], vert_not_in_tet)                  
                        outfile.write(line)                                    
                        i += 1
                elif line.split()[0]=="ND":
                    break 
        for i,n in enumerate(boundary_nodes):
            line = "NDS %d 1\n" % n
            outfile.write(line)
        outfile.close()

#=============================================================================
# Cleanup
#=============================================================================
        tetfiles = ['.1.node','.1.ele','.1.face','.node']
        [os.remove(self.sim_name+x) for x in tetfiles]
        print("--DONE--")

if __name__ == '__main__':
    import argparse
    
#=============================================================================
#     Grabbing arguments dz, top_soil, width, res, xyratio, and max_aspect all 
#     of which have a default value
#=============================================================================    
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', '-n', action='store', type=str, 
            dest='sim_name', default='temp', help='The name of your simulation.')
    parser.add_argument('--height','-z', action='store', type=float, 
            default=3.5, dest='dz', help='The total height of the mesh.')
    parser.add_argument('--top-soil','-s', action='store', type=float, 
            dest='top_soil', default=1.5, 
            help='The height at which the soil layer begins.')
    parser.add_argument('--width','-w', action='store', type=float, 
            dest='width', default=0.02,
            help='The width of the x and y of the mesh.')
    parser.add_argument('--res', '-r', action='store', type=float, dest='res',
            default=0.01, help='The resolution in the z direction of the\
                    non-variable mesh section')
    parser.add_argument('--xyratio','-x', action='store', type=int, default=2,
            dest='xy_ratio', help='The ratio of x and y to resolution in z.')
    parser.add_argument('--max-aspect', '-m', action='store', type=int,
            dest='max_aspect', default=10, 
            help='The maximum aspect ratio between x and y and z in the\
                    variable section of the mesh.')
    
    results = parser.parse_args()

#=============================================================================
#     Parsing out the results for easy passing to the function
#=============================================================================
    sim_name = results.sim_name
    dz = results.dz
    top_soil = results.top_soil 
    res = results.res
    width = results.width
    xy_ratio = results.xy_ratio
    max_aspect = results.max_aspect
#=============================================================================
#    Running the program and creating the mesh.
#=============================================================================
    structmesh = CreateStructuredMesh(sim_name,dz,top_soil, res, width, 
            xy_ratio, max_aspect)
    structmesh.create_struct_mesh()

    



