# -*- coding: utf-8 -*-
"""
FTES Cable structures analysis code

Code created as part of master thesis for
FTES - Máster en Fachadas tecnológicas y envolventes sostenibles
UPM - Universidad Politécnica de Madrid
2022-2023

@author: Gumersindo Freire
gumerfreire@gmail.com
Last update: may 2023
"""

import copy
import math
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.colors
# from glob import glob
from numpy import genfromtxt
import ezdxf
import pandas as pd
from datetime import datetime


class node:
    """
    Structure data to store information of a node in a structure
    
    This class provides structure to store information about coordinates,
    restraints and results of structure analysis of a node in a structure
    composed by nodes and members.
    The structure is internally used by class FTES_Structure.
    """

    def __init__(self,number,x,y,z):
        """
        Initialization of node object
        
        Parameters
        ----------
        number : Integer
            Reference number of the node
        x : Float
            X coord of the node in initial configuration
        y : Float
            Y coord of the node in initial configuration
        z : Float
            Z coord of the node in initial configuration

        Returns
        -------
        None.

        """

        self.number = number
        
        self.x = x
        self.y = y
        self.z = z

        self.restraint_x = False
        self.restraint_y = False
        self.restraint_z = False

        self.force_x = 0
        self.force_y = 0
        self.force_z = 0

        # Data to be filled after structure analysis
        self.ux = np.nan
        self.uy = np.nan
        self.uz = np.nan

        self.disp_xy = np.nan
        self.disp_xz = np.nan
        self.disp_yz = np.nan
        self.disp_total = np.nan

        self.displaced_x = np.nan
        self.displaced_y = np.nan
        self.displaced_z = np.nan

        self.reaction_x = np.nan
        self.reaction_y = np.nan
        self.reaction_z = np.nan



    def restraint(self):
        """
        Retunns bool stating if a node has restraints.

        Returns
        -------
        bool
            True if the node has been locked in X,Y or Z.

        """
        if self.restraint_x == True or self.restraint_y == True or self.restraint_z == True:
            return True
        else:
            return False

    def restraint_vector(self):
        '''
        funcion para devolver el vector de restricciones para componer el vector completo
        '''
        pass

    def set_U(self,ux,uy,uz):
        """
        Sets the displacement data of the node after structure analysis

        Parameters
        ----------
        ux : Float
            X displacement of the node
        uy : TYPE
            Y displacement of the node
        uz : TYPE
            Z displacement of the node

        Returns
        -------
        None.

        """

        self.ux = ux
        self.uy = uy
        self.uz = uz

        self.disp_xy = np.sqrt(ux**2 + uy**2)
        self.disp_xz = np.sqrt(ux**2 + uz**2)
        self.disp_yz = np.sqrt(uy**2 + uz**2)
        self.disp_total = np.sqrt(ux**2 + uy**2 + uz**2)

        self.displaced_x = self.x + ux
        self.displaced_y = self.y + uy
        self.displaced_z = self.z + uz


class member:
    """
    Structure data to store information of a member in a structure
    
    This class provides structure to store information about geometrical, 
    mechanical and other properties of a member in a structure.
    The structure is internally used by class FTES_Structure.
    """

    def __init__(self,number,i,j,E,A,sw,P,cable):
        """
        Initialization of member object

        Parameters
        ----------
        number : Integer
            Reference number fot the member
        i : Integer
            Node i of the member
        j : Integer
            Node j of the member
        E : Float
            Elasticity modulus of the material considered for the member
        A : Float
            Cross section area of the member
        sw : Float
            Self weight. Weight per linear meter of the member
        P : Float
            Prestress of the member
        cable : Bool
            Bool to store if the member is a cable or not (i.e. can resist compression)

        Returns
        -------
        None.

        """

        self.number = number
        self.node_i = i.number
        self.node_j = j.number
        self.ix = i.x
        self.iy = i.y
        self.iz = i.z
        self.jx = j.x
        self.jy = j.y
        self.jz = j.z
        self.E = E
        self.A = A
        self.initial_A = A #Initial Area. This variable saves the initial area en case of slack member overrides value of A=0.
        self.sw = sw
        self.P = P
        self.cable = cable
        self.length = np.sqrt((self.jx-self.ix)**2+(self.jy-self.iy)**2+(self.jz-self.iz)**2)

        # Data to be filled after structure analysis
        self.displaced_ix = np.nan
        self.displaced_iy = np.nan
        self.displaced_iz = np.nan
        self.displaced_jx = np.nan
        self.displaced_jy = np.nan
        self.displaced_jz = np.nan

        self.mbrforce = np.nan
        self.displaced_length = np.nan
        self.tension = np.nan

    def set_results(self,nodes, mbrforce):
        """
        Set results of structura alalysis and perform postprocessing of data

        Parameters
        ----------
        nodes : node object
            Nodes of the structure
        mbrforce : FLoat
            Force resulting in the member after analysis

        Returns
        -------
        None.

        """

        self.mbrforce = mbrforce

        if self.A == 0:
            self.tension = 0
        else:
            self.tension = self.mbrforce / self.A

        self.displaced_ix = nodes[self.node_i].displaced_x
        self.displaced_iy = nodes[self.node_i].displaced_y
        self.displaced_iz = nodes[self.node_i].displaced_z
        self.displaced_jx = nodes[self.node_j].displaced_x
        self.displaced_jy = nodes[self.node_j].displaced_y
        self.displaced_jz = nodes[self.node_j].displaced_z

        self.displaced_length = np.sqrt((self.displaced_jx-self.displaced_ix)**2+(self.displaced_jy-self.displaced_iy)**2+(self.displaced_jz-self.displaced_iz)**2)



class FTES_Structure:
    """
    Main class of FTES Structure cable analysis
    
    This class contains the methods to create a structure composed by nodes and
    bars representing cables or rigid bars for the solution of a system of 
    equations to obtain forces and displacements in the structure.
    
    The structure may be prestressed and the method uses a nonlinear approach
    to solve the analysis taking into account big displacements and nonlinear
    behaviour of the rigidity of the structure, which changes under applied 
    forces.
    """

    def __init__(self):
        """
        Initalization of an FTES_Structure objetc.
        
        The initialization focuses mainly on creating the data structures to
        store the information for the analysis.

        Returns
        -------
        None.

        """

        self.distance_unit = 'm'
        self.force_unit = 'N'

        self.default_E = 130*1000*(1000^2)
        self.default_A = 288/(1000^2)
        self.default_sw = 1
        self.default_P = 35000
        self.default_cable = True

        self.solver_forceIncrements = 20
        self.solver_convThreshold = 10
        self.solver_iterationLimit = 5000
        self.solver_checkSlack = 100 #steps before relaxing cables under compression
        self.solver_selfweight = True

        self.solved_structure = False #Flag to indicate if the structure has been analyzed

        self.str_nodes = {}
        self.str_members = {}

        self.Member_filters = {}
        self.Node_filters = {}
        self.Glazing_filters = {}
        # self.glassunits = {}

        # Dataframes to be filled after analysis
        self.nodes_dataFrame = np.nan
        self.members_dataFrame = np.nan

        self.processed_nodes_dataFrame = np.nan
        self.processed_members_dataFrame = np.nan
        
        # Log record dataframe
        self.printloginscreen = True
        self.log_data = pd.DataFrame(columns=['Date', 'Time', 'Log record'])
        
        # Displacement postprocessing
        self.angles_variation = pd.DataFrame(columns=['Node', 'Bar 1', 'Bar 2', 'Initial angle', 'Final angle', 'Variation (ABS)'])
        self.glazing_variation = False

    def FDM(self):
        '''
        Force Density Method implementation
        
        Force Density methos implementation as described by H.J.Scheck (1974).
        This method allows to find a possible equilibrium configuration for a
        set of nodes and members under certain forces.
        
        This implementation uses the data of the structure (in termos of member
        connections and loads) to find a compatible form. The method does not
        know the original coordinates of the nodes, only the map of connections
        between members.
        
        This method exports a DXF file. The function can be further developed to
        allow FTES_Structure to use the found form to perform the analysis.
        '''
        
        m = len(self.str_members) #number of members in structure
        ns = len(self.str_nodes) #number of nodes in structure
        
        C_c = np.zeros([m,ns])
        l = np.zeros(m)
        s = np.zeros((m,1))
        
        px_c = np.zeros((ns,1))
        py_c = np.zeros((ns,1))
        pz_c = np.zeros((ns,1))
        x_c = np.zeros((ns,1))
        y_c = np.zeros((ns,1))
        z_c = np.zeros((ns,1))
        
        for i in self.str_members:
            C_c[i,self.str_members[i].node_i] = 1
            C_c[i,self.str_members[i].node_j] = -1
            l[i] = self.str_members[i].length
            s[i,0] = self.str_members[i].P
        
        for i in self.str_nodes:
            px_c[i,0] = self.str_nodes[i].force_x
            py_c[i,0] = self.str_nodes[i].force_y
            pz_c[i,0] = self.str_nodes[i].force_z
            
            x_c[i,0] = self.str_nodes[i].x
            y_c[i,0] = self.str_nodes[i].y
            z_c[i,0] = self.str_nodes[i].z
        
        freenodes = []
        fixednodes = []
        for i in self.str_nodes:
            
            if self.str_nodes[i].restraint_x == True or self.str_nodes[i].restraint_y == True or self.str_nodes[i].restraint_z == True:
                fixednodes.append(i)
            else:
                freenodes.append(i)
        
        C = C_c[:,freenodes]
        Cf = C_c[:,fixednodes]
        
        px = px_c[freenodes,:]
        py = py_c[freenodes,:]
        pz = pz_c[freenodes,:]
        
        xf = x_c[fixednodes,:]
        yf = y_c[fixednodes,:]
        zf = z_c[fixednodes,:]
        
        L = np.diag(l)
        Q = np.matmul(np.linalg.inv(L),s)
        
        D = np.matrix(np.matmul(C.T, (Q*C)))
        Df = np.matrix(np.matmul(C.T, (Q*Cf)))
        print('D ', D.shape)
        print(D)
        print('Df ', Df.shape)
        print(Df)
        
        x = np.linalg.inv(D) * (px - Df*xf)
        y = np.linalg.inv(D) * (py - Df*yf)
        z = np.linalg.inv(D) * (pz - Df*zf)
        
        #Reorder deformed coords
        x_fdm = x_c
        y_fdm = y_c
        z_fdm = z_c
        
        indice = 0
        for i in freenodes:
            
            x_fdm[i,0] = x[indice]
            y_fdm[i,0] = y[indice]
            z_fdm[i,0] = z[indice]
            indice = indice+1
            
        #DXF export
        doc = ezdxf.new()
        msp = doc.modelspace()

        doc.layers.new(name="FDM structure", dxfattribs={"color": "30"})
        doc.layers.new(name="Original structure", dxfattribs={"color": "3"})

        for i in self.str_members:

            ix = self.str_members[i].ix
            iy = self.str_members[i].iy
            iz = self.str_members[i].iz
            jx = self.str_members[i].jx
            jy = self.str_members[i].jy
            jz = self.str_members[i].jz

            msp.add_line((ix, iy, iz), (jx, jy, jz), dxfattribs={"layer": "Original structure"})

            dix = x_fdm[self.str_members[i].node_i,0]
            diy = y_fdm[self.str_members[i].node_i,0]
            diz = z_fdm[self.str_members[i].node_i,0]
            djx = x_fdm[self.str_members[i].node_j,0]
            djy = y_fdm[self.str_members[i].node_j,0]
            djz = z_fdm[self.str_members[i].node_j,0]
            
            msp.add_line((dix, diy, diz), (djx, djy, djz), dxfattribs={"layer": "FDM structure"})

        doc.saveas('FDM.dxf')

    def log_record(self, inputMessage, *args):
        """
        Method to print messages on screen while recording them into a dataframe
        
        This method works as a substitute of the print function, acting as a
        print function but also recording the messages to a dataframe to allow
        the export of the recorded data as a log into the results file.

        Parameters
        ----------
        inputMessage : String
            String to be shown on screen and recorded to a log.
        *args : various
            Additional arguments to add to the string.

        Returns
        -------
        None.

        """
        
        texto = inputMessage
        for i in args:
            texto = texto+str(i)
            
        print(inputMessage)
        
        now = datetime.now()
        fecha = now.strftime("%d-%m-%Y")
        hora = now.strftime("%H:%M:%S")
        logrow = pd.DataFrame.from_dict({'Date': [fecha], 'Time': [hora], 'Log record': [texto]})
        
        self.log_data = self.log_data.append(logrow, ignore_index=True)
    
    def set_defaultMemberValues(self,E,A,P,sw, cable):
        """
        Set default mechanical values to assign to the structure

        Parameters
        ----------
        E : Float
            Default elasticity modulus
        A : Float
            Default cross section member area
        P : Float
            Default prestress
        sw : Float
            Default self weight
        cable : Bool
            Default definition as cable or rigid bar

        Returns
        -------
        None.

        """
        self.default_E = E
        self.default_A = A
        self.default_P = P
        self.default_sw = sw
        self.default_cable = cable

    def set_solverParameters(self,forceIncrements,convThreshold, iterationLimit):
        """
        Method to set solver parameters.

        Parameters
        ----------
        forceIncrements : Integer
            Number of increments to divide the applied force into incremental steps
        convThreshold : Integer
            Threshold to consider convergence in the iterative method
        iterationLimit : Integer
            Maximum interations before considering non-convergence

        Returns
        -------
        None.

        """
        self.solver_forceIncrements = forceIncrements
        self.solver_convThreshold = convThreshold
        self.solver_iterationLimit = iterationLimit
        
    def import_CsvStructure(self,vertices,edges):
        """
        Method to import structure data from CSV files
        
        This method allows to introduce structure data from tho CSV files. The
        vertices file expects a CSV file containing a table with three columns.
        Each row corresponds to a node and the columns correspond to its X,Y,Z
        coordinates.
        The CSV file definig the members expects a file containing a table with
        two columns. Each row corresponds to a member, and the columns store the
        numbers for the initial and final nodes of the member.

        Parameters
        ----------
        vertices : TYPE
            DESCRIPTION.
        edges : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.str_nodes = {}
        self.str_members = {}
        nodelist = []
        memberlist = []

        nodes = genfromtxt(vertices, delimiter=',')
        members = genfromtxt(edges, delimiter=',')

        for i in range(nodes.shape[0]):
            self.str_nodes[i] = node(i,nodes[i,0],nodes[i,1],nodes[i,2])
            nodelist.append(i)

        for i in range(members.shape[0]):
            self.str_members[i] = member(i,self.str_nodes[members[i,0]],self.str_nodes[members[i,1]],self.default_E,self.default_A,self.default_sw,self.default_P,self.default_cable)
            memberlist.append(i)
            
        self.nNodes = len(self.str_nodes) #number of nodes
        self.nMembers = len(self.str_members) #number of members
        self.nDoF = len(self.str_nodes)*3
        
        self.Node_filters['Allnodes'] = nodelist
        self.Member_filters['Allmembers'] = memberlist
        
        self.forceVector = np.zeros((self.nDoF,1))
        self.restrainedIndex = []

        self.log_record(f"Imported structure data with {self.nNodes} nodes and  {self.nMembers} members")

        lengths = []
        for i in self.str_members:
            lengths.append(self.str_members[i].length)
        if min(lengths) < 1e-5:
            self.log_record("MEMBERS DETECTED WITH LEGTH 0. PLEASE CORRECT STRUCTURE DATA")

    def import_DxfStructure(self,file):
        """
        Method to import structure data from a DXF file.
        
        This method allows the introduction of a structure geometry and data
        from a DXF file containing a series of lines in layer 0. The function
        interprets each line with a member and detects the endpoints in contact
        with other members to set the same node for various members.
        
        It is possible to import as well member filters and node filters through
        lines and points in layers starting with MF_ and NF_ respectively. See
        documentation for more information on member, node and glazing filters.

        Parameters
        ----------
        file : String
            String with path and name of a DXF file to import.

        Raises
        ------
        ValueError
            The nodes or members in filter layers must be coincident with elements
            in the structure.

        Returns
        -------
        None.

        """
        
        #precission parameters for DXF import
        #modify these parameters to facilitate DXF import where points have small inaccuracies.
        decimals_precision = 4 #decimals to round coordinates
        dist_precision = 0.1 #max distance from points to nodes to consider the same point
        
        self.str_nodes = {}
        self.str_members = {}

        doc = ezdxf.readfile(file)
        msp = doc.modelspace()
        layers = doc.layers

        MF_layers = [] #member filter layers
        NF_layers = [] #node filter layers
        GF_layers = [] #glazing unit filter layers
        for i in layers:
            if i.dxf.name[0:3] == 'MF_': MF_layers.append(i.dxf.name)
            if i.dxf.name[0:3] == 'NF_': NF_layers.append(i.dxf.name)
            if i.dxf.name[0:3] == 'GF_': GF_layers.append(i.dxf.name)

        lines = msp.query('LINE[layer=="0"]')

        nodes = np.empty([0,3])
        nMembers = 0

        for line in lines:

            nodes = np.append(nodes,np.array([[line.dxf.start.x,line.dxf.start.y,line.dxf.start.z]]),axis = 0) #append start point
            nodes = np.append(nodes,np.array([[line.dxf.end.x,line.dxf.end.y,line.dxf.end.z]]),axis = 0) #append start point

            nMembers +=1

        nodes = np.around(nodes, decimals=decimals_precision)
        nodes = np.unique(nodes, axis = 0)
        nNodes = nodes.shape[0]

        members = np.empty([nMembers,2])

        for m,line in enumerate(lines):

            start_dist = np.zeros([nNodes])
            end_dist = np.zeros([nNodes])

            for i in range(nNodes):
                 start_dist[i] = np.sqrt((nodes[i,0]-line.dxf.start.x)**2 + (nodes[i,1]-line.dxf.start.y)**2 + (nodes[i,2]-line.dxf.start.z)**2)
                 end_dist[i] = np.sqrt((nodes[i,0]-line.dxf.end.x)**2 + (nodes[i,1]-line.dxf.end.y)**2 + (nodes[i,2]-line.dxf.end.z)**2)

            members[m,0] = np.where(start_dist == start_dist.min())[0]
            members[m,1] = np.where(end_dist == end_dist.min())[0]

        # Export data to structure info
        nodelist = []
        memberlist = []
        for i in range(nNodes):
            self.str_nodes[i] = node(i,nodes[i,0],nodes[i,1],nodes[i,2])
            nodelist.append(i)

        for i in range(nMembers):
            self.str_members[i] = member(i,self.str_nodes[members[i,0]],self.str_nodes[members[i,1]],self.default_E,self.default_A,self.default_sw,self.default_P,self.default_cable)
            memberlist.append(i)
            
        self.Node_filters['Allnodes'] = nodelist
        self.Member_filters['Allmembers'] = memberlist


        # Collect member filters
        for capa in MF_layers:

            filterlines = msp.query('LINE[layer=="' + capa + '"]')

            filtermembers = np.empty([len(filterlines),2])
            for m,line in enumerate(filterlines):
                start_dist = np.zeros([nNodes])
                end_dist = np.zeros([nNodes])

                for i in range(nNodes):
                    # print('line start ',line.dxf.start.x,line.dxf.start.y,line.dxf.start.z)
                    # print('line end ',line.dxf.end.x,line.dxf.end.y,line.dxf.end.z)
                    # print('Node ',nodes[i,0],nodes[i,1],nodes[i,2])
                    start_dist[i] = np.sqrt((nodes[i,0]-line.dxf.start.x)**2 + (nodes[i,1]-line.dxf.start.y)**2 + (nodes[i,2]-line.dxf.start.z)**2)
                    end_dist[i] = np.sqrt((nodes[i,0]-line.dxf.end.x)**2 + (nodes[i,1]-line.dxf.end.y)**2 + (nodes[i,2]-line.dxf.end.z)**2)
                    # print('distance ', start_dist[i])
                    # print('distance ', end_dist[i])

                if start_dist.min() > dist_precision: raise ValueError('Members in filter not coincident with any bar in structure')
                if end_dist.min() > dist_precision: raise ValueError('Members in filter not coincident with any bar in structure')

                filtermembers[m,0] = np.where(start_dist == start_dist.min())[0]
                filtermembers[m,1] = np.where(end_dist == end_dist.min())[0]

            sortedMembers = np.zeros([nMembers,2])
            for i in range(nMembers):
                sortedMembers[i] = np.sort(members[i])

            sortedfilterMembers = np.zeros([len(filterlines),2])
            for i in range(len(filterlines)):
                sortedfilterMembers[i] = np.sort(filtermembers[i])

            filterIndex = []

            for i in range(len(filterlines)):

                mbr = sortedfilterMembers[i]

                for j in range(nMembers):
                    if sortedMembers[j,0] == mbr[0] and sortedMembers[j,1] == mbr[1]: filterIndex.append(j)

            self.Member_filters[capa[3:]] = filterIndex

        # Collect node filters
        for capa in NF_layers:
            filterpoints = msp.query('POINT[layer=="' + capa + '"]')

            # filterNodes = np.empty([len(filterpoints)], dtype='int')
            filterNodes = []

            for m,point in enumerate(filterpoints):
                point_dist = np.zeros([nNodes])
                x = point.dxf.location.x
                y = point.dxf.location.y
                z = point.dxf.location.z

                for i in range(nNodes):
                    point_dist[i] = np.sqrt((nodes[i,0]-x)**2 + (nodes[i,1]-y)**2 + (nodes[i,2]-z)**2)

                if point_dist.min() > dist_precision: raise ValueError('Nodes in filter not coincident with any node in structure')

                filterNodes.append(int(np.where(point_dist == point_dist.min())[0]))
            self.Node_filters[capa[3:]] = filterNodes

        # Collect glazing filters
        for capa in GF_layers:
            filterpoints = msp.query('POINT[layer=="' + capa + '"]')

            # filterNodes = np.empty([len(filterpoints)], dtype='int')
            filterNodes = []

            for m,point in enumerate(filterpoints):
                point_dist = np.zeros([nNodes])
                x = point.dxf.location.x
                y = point.dxf.location.y
                z = point.dxf.location.z

                for i in range(nNodes):
                    point_dist[i] = np.sqrt((nodes[i,0]-x)**2 + (nodes[i,1]-y)**2 + (nodes[i,2]-z)**2)

                if point_dist.min() > dist_precision: raise ValueError('Nodes in filter not coincident with any node in structure')

                filterNodes.append(int(np.where(point_dist == point_dist.min())[0]))
            if len(filterNodes)==4:
                self.Glazing_filters[capa[3:]] = filterNodes
            else:
                self.log_record(capa + ' - Glazing defined with vertices number other than 4. Not added to model info.')

        self.nNodes = nNodes #number of nodes
        self.nMembers = nMembers #number of members
        self.nDoF = nNodes*3

        self.forceVector = np.zeros((self.nDoF,1))
        self.restrainedIndex = []

        self.log_record(f"Imported structure data with {self.nNodes} nodes and  {self.nMembers} members")
        if len(MF_layers) > 0: self.log_record("Imported Member filters: ", MF_layers)
        if len(NF_layers) > 0: self.log_record("Imported Node filters: ", NF_layers)
        if len(self.Glazing_filters.keys()) > 0: self.log_record("Imported Glazing filters: ", list(self.Glazing_filters.keys()))

        #Check length 0 members in data
        lengths = []
        for i in self.str_members:
            lengths.append(self.str_members[i].length)
        if min(lengths) < 1e-5:
            self.log_record("MEMBERS DETECTED WITH LEGTH 0. PLEASE CORRECT STRUCTURE DATA")

        #Check overlapping members in data
        sortedMembers = np.zeros([nMembers,2])
        for i in range(nMembers):
            sortedMembers[i] = np.sort(members[i])

        if np.unique(sortedMembers, axis = 0).shape[0] < sortedMembers.shape[0]:
            self.log_record("OVERLAPPED MEMBERS DETECTED. PLEASE CORRECT STRUCTURE DATA")
            
    def export_DxfStructure(self, file):
        """
        Method to generate a DXF file with the original and deformed structure
        
        This methos allows to draw the analyzed structure ina  DXF file containing
        the original configuration of the structure and the deformed state after
        the anslysis.

        Parameters
        ----------
        file : String
            Path and name of the DXF file to be created

        Returns
        -------
        None.

        """

        doc = ezdxf.new()
        msp = doc.modelspace()

        doc.layers.new(name="Deflected structure", dxfattribs={"color": "30"})
        doc.layers.new(name="Original structure", dxfattribs={"color": "3"})

        for i in self.str_members:

            ix = self.str_members[i].ix
            iy = self.str_members[i].iy
            iz = self.str_members[i].iz
            jx = self.str_members[i].jx
            jy = self.str_members[i].jy
            jz = self.str_members[i].jz

            msp.add_line((ix, iy, iz), (jx, jy, jz), dxfattribs={"layer": "Original structure"})

            if self.solved_structure:
                dix = self.str_members[i].displaced_ix
                diy = self.str_members[i].displaced_iy
                diz = self.str_members[i].displaced_iz
                djx = self.str_members[i].displaced_jx
                djy = self.str_members[i].displaced_jy
                djz = self.str_members[i].displaced_jz
                msp.add_line((dix, diy, diz), (djx, djy, djz), dxfattribs={"layer": "Deflected structure"})

        doc.saveas(file)

    def plotStructure(self, plane='auto', **kwargs):
        """
        Method to plot ths structure on screen and save it to a PDF file
        
        This method allows to plot theimported structure on screen showing the
        automatic numeration of nodes and members. It is possible to save the
        plot as a PDF file for further consultation of the node and member numbers.
        
        The function detects the main plane of the structure to produce a 2D
        plot either in XY, XZ or YZ plane.
        
        This function can be further developed for more advanced visualizations
        of the structure and to show more information such as forces, axis or
        other selected information.

        Parameters
        ----------
        plane : String, optional
            with the default option 'auto' the mathod detects the main plane
            to plot. Otherwise a plane can be forced for representation, either
            XY, XZ or YZ.
        **kwargs : Dictionary
            if savePdf kwarg present, a PDF file of the plot is saved.

        Returns
        -------
        None.

        """
        
        #Select plotting plane
        if plane == 'auto':
            
            coordx = []
            coordy = []
            coordz = []
            
            for i,node in enumerate(self.str_nodes):
                # print(self.str_nodes[i].x)
                coordx.append(self.str_nodes[i].x)
                coordy.append(self.str_nodes[i].y)
                coordz.append(self.str_nodes[i].z)
            
            Xdif = abs(max(coordx) - min(coordx))
            Ydif = abs(max(coordy) - min(coordy))
            Zdif = abs(max(coordz) - min(coordz))
            
            if Xdif < Ydif and Xdif < Zdif:
                self.log_record('Plane autodetected for plotting: YZ')
                plotplane = 'YZ'
                plotdim1 = Ydif
                plotdim2 = Zdif
                txtdim = np.max([plotdim1,plotdim2])/3
            elif Ydif < Xdif and Ydif < Zdif:
                self.log_record('Plane autodetected for plotting: XZ')
                plotplane = 'XZ'
                plotdim1 = Xdif
                plotdim2 = Zdif
                txtdim = np.max([plotdim1,plotdim2])/3
            elif Zdif < Xdif and Zdif < Ydif:
                self.log_record('Plane autodetected for plotting: XY')
                plotplane = 'XY'
                plotdim1 = Xdif
                plotdim2 = Zdif
                txtdim = np.max([plotdim1,plotdim2])/3
                
        elif plane == 'XY':
            self.log_record('Plotting structure. Selected plane to plot: XY')
            plotplane = 'XY'
        elif plane == 'XZ':
            self.log_record('Plotting structure. Selected plane to plot: XZ')
            plotplane = 'XZ'
        elif plane == 'YZ':
            self.log_record('Plotting structure. Selected plane to plot: YZ')
            plotplane = 'YZ'
        else:
            print('Plane selected for plotting not valid. Please select XY, XZ, YZ or auto.')
            return

        #Kwargs processing
        NF_show = False
        MF_show = False
        GF_show = False
        NF_indices=[]
        MF_indices=[]
        GF_indices=[]
        
        if 'Nodefilter' in kwargs:
            NF_show = True
            NF_indices = self.Node_filters[kwargs['Nodefilter']]
        
        if 'Memberfilter' in kwargs:
            MF_show = True
            MF_indices = self.Member_filters[kwargs['Memberfilter']]
            
        if 'Glazingfilter' in kwargs:
            GF_show = True
            GF_indices = self.Glazing_filters[kwargs['Glazingfilter']]
        
        #Plotting
        fig, ax = plt.subplots(figsize=(plotdim1,plotdim2))
        text_offset = 0.05

        if plotplane == 'XY':
            for i in self.str_members:       
                if MF_show == True and (i in MF_indices):
                    ax.plot([self.str_members[i].ix,self.str_members[i].jx],[self.str_members[i].iy,self.str_members[i].jy], 'r', linewidth = 0.5)
                    ax.text((self.str_members[i].ix+self.str_members[i].jx)/2+text_offset,(self.str_members[i].iy+self.str_members[i].jy)/2+text_offset, str(i), fontsize=txtdim, color='red')
                else:
                    ax.plot([self.str_members[i].ix,self.str_members[i].jx],[self.str_members[i].iy,self.str_members[i].jy], 'k', linewidth = 0.5)
                    ax.text((self.str_members[i].ix+self.str_members[i].jx)/2+text_offset,(self.str_members[i].iy+self.str_members[i].jy)/2+text_offset, str(i), fontsize=txtdim, color='grey')
    
            for i in self.str_nodes:
                
                if (NF_show == True and (i in NF_indices)) or (GF_show == True and (i in GF_indices)):
                    ax.plot(self.str_nodes[i].x, self.str_nodes[i].y, 'r.')
                    ax.text(self.str_nodes[i].x+text_offset, self.str_nodes[i].y+text_offset, str(i), fontsize=txtdim, color='red')
                else:
                    ax.plot(self.str_nodes[i].x, self.str_nodes[i].y, 'k.')
                    ax.text(self.str_nodes[i].x+text_offset, self.str_nodes[i].y+text_offset, str(i), fontsize=txtdim)
            ax.set_xlabel('X Coord (m)')
            ax.set_ylabel('Y Coord (m)')
            ax.set_title('STRUCTURE PLANE XY')

        elif plotplane == 'XZ':
            for i in self.str_members:       
                if MF_show == True and (i in MF_indices):
                    ax.plot([self.str_members[i].ix,self.str_members[i].jx],[self.str_members[i].iz,self.str_members[i].jz], 'r', linewidth = 0.5)
                    ax.text((self.str_members[i].ix+self.str_members[i].jx)/2+text_offset,(self.str_members[i].iz+self.str_members[i].jz)/2+text_offset, str(i), fontsize=txtdim, color='red')
                else:
                    ax.plot([self.str_members[i].ix,self.str_members[i].jx],[self.str_members[i].iz,self.str_members[i].jz], 'k', linewidth = 0.5)
                    ax.text((self.str_members[i].ix+self.str_members[i].jx)/2+text_offset,(self.str_members[i].iz+self.str_members[i].jz)/2+text_offset, str(i), fontsize=txtdim, color='grey')
    
            for i in self.str_nodes:
                
                if (NF_show == True and (i in NF_indices)) or (GF_show == True and (i in GF_indices)):
                    ax.plot(self.str_nodes[i].x, self.str_nodes[i].z, 'r.')
                    ax.text(self.str_nodes[i].x+text_offset, self.str_nodes[i].z+text_offset, str(i), fontsize=txtdim, color='red')
                else:
                    ax.plot(self.str_nodes[i].x, self.str_nodes[i].z, 'k.')
                    ax.text(self.str_nodes[i].x+text_offset, self.str_nodes[i].z+text_offset, str(i), fontsize=txtdim)
            ax.set_xlabel('X Coord (m)')
            ax.set_ylabel('Z Coord (m)')
            ax.set_title('STRUCTURE PLANE XZ')
                    
        elif plotplane == 'YZ':
            for i in self.str_members:       
                if MF_show == True and (i in MF_indices):
                    ax.plot([self.str_members[i].iy,self.str_members[i].jy],[self.str_members[i].iz,self.str_members[i].jz], 'r', linewidth = 0.5)
                    ax.text((self.str_members[i].iy+self.str_members[i].jy)/2+text_offset,(self.str_members[i].iz+self.str_members[i].jz)/2+text_offset, str(i), fontsize=txtdim, color='red')
                else:
                    ax.plot([self.str_members[i].iy,self.str_members[i].jy],[self.str_members[i].iz,self.str_members[i].jz], 'k', linewidth = 0.5)
                    ax.text((self.str_members[i].iy+self.str_members[i].jy)/2+text_offset,(self.str_members[i].iz+self.str_members[i].jz)/2+text_offset, str(i), fontsize=txtdim, color='grey')
    
            for i in self.str_nodes:
                
                if (NF_show == True and (i in NF_indices)) or (GF_show == True and (i in GF_indices)):
                    ax.plot(self.str_nodes[i].y, self.str_nodes[i].z, 'r.')
                    ax.text(self.str_nodes[i].y+text_offset, self.str_nodes[i].z+text_offset, str(i), fontsize=txtdim, color='red')
                else:
                    ax.plot(self.str_nodes[i].y, self.str_nodes[i].z, 'k.')
                    ax.text(self.str_nodes[i].y+text_offset, self.str_nodes[i].z+text_offset, str(i), fontsize=txtdim)
            ax.set_xlabel('Y Coord (m)')
            ax.set_ylabel('Z Coord (m)')
            ax.set_title('STRUCTURE PLANE YZ')

        ax.set_aspect('equal', adjustable='box')
        
        if 'savePdf' in kwargs:
            fig.savefig(kwargs['savePdf'])




    def structureInfo(self):
        """
        Method to show general information of the structure on screen.
        
        The development of this method was stopped in favour of exporting the
        information of the structure to a XLS file from a dataframe.

        Returns
        -------
        None.

        """
        
        print('\n·············································')
        print('\nFTES. Cable structure analysis')
        print('\n·············································')
        print(f'\nImported structure with {self.nNodes} nodes and {self.nMembers} members\n')

        x=[]
        y=[]
        z=[]
        for i in self.str_nodes:
            x.append(self.str_nodes[i].x)
            y.append(self.str_nodes[i].y)
            z.append(self.str_nodes[i].z)

        dimx = max(x) - min(x)
        dimy = max(y) - min(y)
        dimz = max(z) - min(z)

        print(f'General dimensions: X = {dimx} m, Y = {dimy} m, Z = {dimz} m')
        print(f'X coordinates vary from {min(x)} to {max(x)} m')
        print(f'Y coordinates vary from {min(y)} to {max(y)} m')
        print(f'Z coordinates vary from {min(z)} to {max(z)} m')

        print('\n·············································')

        l=[]
        for i in self.str_members:
            l.append(self.str_members[i].length)

        print(f'Length of members vary from {min(l)} to {max(l)} m')
        print('\n·············································')

        print('\nNODE COORDINATES')
        print('Node','\t','X','\t\t','Y','\t\t','Z')
        for i in self.str_nodes:
            print(str(self.str_nodes[i].number),'\t\t',str(self.str_nodes[i].x),'\t',str(self.str_nodes[i].y),'\t',str(self.str_nodes[i].z))
        print('\n·············································')
        print('\nRESTRAINTS')
        print('Node','\t','X','\t\t','Y','\t\t','Z')
        for i in self.str_nodes:
            print(str(self.str_nodes[i].number),'\t\t',str(self.str_nodes[i].restraint_x),'\t',str(self.str_nodes[i].restraint_y),'\t',str(self.str_nodes[i].restraint_z))
        print('\n·············································')
        print('\nMEMBERS GEOMETRY')
        print('Member','\t','Node i','\t','Node j','\t','Length')
        for i in self.str_members:
            print(str(self.str_members[i].number),'\t\t',str(self.str_members[i].node_i),'\t\t\t',str(self.str_members[i].node_j),'\t\t\t',str(self.str_members[i].length) )
        print('\n·············································')

    def addForce(self, force, axis, **kwargs):
        """
        Method to add forces to selected nodes in the structure.

        Parameters
        ----------
        force : Float
            Magnitude of force to de applied
        axis : String
            Direction (X,Y or Z) in which the force may be applied
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of node numbers to apply force
                Nodefilter - Name of a existing node filter
                Nodelist - List of node numbers to apply the force
                Savefilter - Save the node list as a filter

        Returns
        -------
        int
            0 if nothing added

        """

        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Nodefilter' in kwargs:
            indices = self.Node_filters[kwargs['Nodefilter']]
            
        elif 'Nodelist' in kwargs:
            indices = kwargs['Nodelist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:
            if not kwargs['Savefilter'] in self.Node_filters: self.Node_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')
            
              
        if isinstance(axis, str): axis = [axis]

        if 'x' in axis:
            for i in indices:
                self.forceVector[3*i] = self.forceVector[3*i] + force
                self.str_nodes[i].force_x = self.str_nodes[i].force_x + force

        if 'y' in axis:
            for i in indices:
                self.forceVector[3*i+1] = self.forceVector[3*i+1] + force
                self.str_nodes[i].force_y = self.str_nodes[i].force_y + force

        if 'z' in axis:
            for i in indices:
                self.forceVector[3*i+2] = self.forceVector[3*i+2] + force
                self.str_nodes[i].force_z = self.str_nodes[i].force_z + force

    def addRestraint(self, axis, **kwargs):
        """
        Method to add restraints to selected nodes in the structure.

        Parameters
        ----------
        axis : String
            Direction (X,Y or Z) in which the restraint may be applied
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of node numbers to apply force
                Nodefilter - Name of a existing node filter
                Nodelist - List of node numbers to apply the force
                Savefilter - Save the node list as a filter

        Returns
        -------
        int
            0 if nothing added

        """

        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Nodefilter' in kwargs:
            indices = self.Node_filters[kwargs['Nodefilter']]
            
        elif 'Nodelist' in kwargs:
            indices = kwargs['Nodelist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:  
            if not kwargs['Savefilter'] in self.Node_filters: self.Node_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')
        
        if 'x' in axis:
            for n in indices:
                self.str_nodes[n].restraint_x = True

        if 'y' in axis:
            for n in indices:
                self.str_nodes[n].restraint_y = True

        if 'z' in axis:
            for n in indices:
                self.str_nodes[n].restraint_z = True

        # Update estrainedIndex and freeDoF lists
        self.restrainedIndex = []

        for n in self.str_nodes:
            # print(n)
            if self.str_nodes[n].restraint_x == True: self.restrainedIndex.append(3*self.str_nodes[n].number)
            if self.str_nodes[n].restraint_y == True: self.restrainedIndex.append(3*self.str_nodes[n].number+1)
            if self.str_nodes[n].restraint_z == True: self.restrainedIndex.append(3*self.str_nodes[n].number+2)

    def export_Data(self, file):
        """
        Method to export analysis data to a XLS file
        
        This method gathers the information of the structure after an analysis
        and creates three dataframes to export into a single XLS file as a 
        complete report of the result of analysis.
        The XLS file may contain the following information:
            - LOG. Contains the complete history of messages during the analysis
            - NODES. Contains complete information about nodes in the structure
            - MEMBERS. Contains complete information about the members in the structure
            - ANGLES. Contains iformation about angles variation in cables
            - GLAZING. Contains information about deformation in defined glazing units

        Parameters
        ----------
        file : String
            Path and name of XLS file to be created.

        Returns
        -------
        None.

        """

        if not self.solved_structure:
            print('The structure has not been analyzed. Please perform structural analysis before exporting data.')
        else:

            #Nodes info
            nodeTable = np.empty((0,23))

            for i in self.str_nodes:
                nodeRow = np.empty(23)

                nodeRow[0] = self.str_nodes[i].number
                nodeRow[1] = self.str_nodes[i].x
                nodeRow[2] = self.str_nodes[i].y
                nodeRow[3] = self.str_nodes[i].z

                nodeRow[4] = self.str_nodes[i].restraint_x
                nodeRow[5] = self.str_nodes[i].restraint_y
                nodeRow[6] = self.str_nodes[i].restraint_z

                nodeRow[7] = self.str_nodes[i].ux
                nodeRow[8] = self.str_nodes[i].uy
                nodeRow[9] = self.str_nodes[i].uz

                nodeRow[10] = self.str_nodes[i].disp_xy
                nodeRow[11] = self.str_nodes[i].disp_xz
                nodeRow[12] = self.str_nodes[i].disp_yz
                nodeRow[13] = self.str_nodes[i].disp_total

                nodeRow[14] = self.str_nodes[i].displaced_x
                nodeRow[15] = self.str_nodes[i].displaced_y
                nodeRow[16] = self.str_nodes[i].displaced_z

                nodeRow[17] = self.str_nodes[i].reaction_x
                nodeRow[18] = self.str_nodes[i].reaction_y
                nodeRow[19] = self.str_nodes[i].reaction_z

                nodeRow[20] = self.str_nodes[i].force_x
                nodeRow[21] = self.str_nodes[i].force_y
                nodeRow[22] = self.str_nodes[i].force_z

                nodeTable = np.append(nodeTable, [nodeRow], axis=0)

            nodelabels = ['Node number',
                          'Coord. X',
                          'Coord. Y',
                          'Coord. Z',
                          'Restrained X',
                          'Restrained Y',
                          'Restrained Z',
                          'Ux',
                          'Uy',
                          'Uz',
                          'Displacement in plane XY',
                          'Displacement in plane XZ',
                          'Displacement in plane YZ',
                          'Total displacement',
                          'Displaced Coord. X',
                          'Displaced Coord. Y',
                          'Displaced Coord. Z',
                          'Reaction X',
                          'Reaction Y',
                          'Reaction Z',
                          'Force X',
                          'Force Y',
                          'Force Z']

            nodes_df = pd.DataFrame(nodeTable, columns = nodelabels).set_index(['Node number'])
            self.nodes_dataFrame = nodes_df

            #Postprocessing to make it more readable
            nodes_df['Restrained X'] = nodes_df['Restrained X'].replace([0,1],['False','True'])
            nodes_df['Restrained Y'] = nodes_df['Restrained Y'].replace([0,1],['False','True'])
            nodes_df['Restrained Z'] = nodes_df['Restrained Z'].replace([0,1],['False','True'])

            nodes_df['Ux'] = round(nodes_df['Ux']*1000,4)
            nodes_df['Uy'] = round(nodes_df['Uy']*1000,4)
            nodes_df['Uz'] = round(nodes_df['Uz']*1000,4)

            nodes_df['Displacement in plane XY'] = round(nodes_df['Displacement in plane XY']*1000,4)
            nodes_df['Displacement in plane XZ'] = round(nodes_df['Displacement in plane XZ']*1000,4)
            nodes_df['Displacement in plane YZ'] = round(nodes_df['Displacement in plane YZ']*1000,4)
            nodes_df['Total displacement'] = round(nodes_df['Total displacement']*1000,4)

            nodes_df['Reaction X'] = round(nodes_df['Reaction X']/1000,4)
            nodes_df['Reaction Y'] = round(nodes_df['Reaction Y']/1000,4)
            nodes_df['Reaction Z'] = round(nodes_df['Reaction Z']/1000,4)

            nodes_df['Force X'] = round(nodes_df['Force X']/1000,4)
            nodes_df['Force Y'] = round(nodes_df['Force Y']/1000,4)
            nodes_df['Force Z'] = round(nodes_df['Force Z']/1000,4)

            processed_nodelabels = ['Node number',
                          'Coord. X',
                          'Coord. Y',
                          'Coord. Z',
                          'Restrained X',
                          'Restrained Y',
                          'Restrained Z',
                          'Ux (mm)',
                          'Uy (mm)',
                          'Uz (mm)',
                          'Displacement in plane XY (mm)',
                          'Displacement in plane XZ (mm)',
                          'Displacement in plane YZ (mm)',
                          'Total displacement (mm)',
                          'Displaced Coord. X',
                          'Displaced Coord. Y',
                          'Displaced Coord. Z',
                          'Reaction X (kN)',
                          'Reaction Y (kN)',
                          'Reaction Z (kN)',
                          'Force X (kN)',
                          'Force Y (kN)',
                          'Force Z (kN)']
            dict_labels = {}

            for i,lab in enumerate(nodelabels):
                dict_labels[lab] = processed_nodelabels[i]

            nodes_df = nodes_df.rename(dict_labels, axis='columns')
            self.processed_nodes_dataFrame = nodes_df



            # Members info
            memberTable = np.empty((0,12))

            for i in self.str_members:

                memberRow = np.empty(12)

                memberRow[0] = self.str_members[i].number
                memberRow[1] = self.str_members[i].node_i
                memberRow[2] = self.str_members[i].node_j
                memberRow[3] = self.str_members[i].E
                memberRow[4] = self.str_members[i].initial_A
                memberRow[5] = self.str_members[i].sw
                memberRow[6] = self.str_members[i].P
                memberRow[7] = self.str_members[i].cable
                memberRow[8] = self.str_members[i].length
                memberRow[9] = self.str_members[i].mbrforce
                memberRow[10] = self.str_members[i].tension
                memberRow[11] = self.str_members[i].displaced_length

                memberTable = np.append(memberTable, [memberRow], axis=0)

            memberlabels = ['Member number',
                          'node i',
                          'node j',
                          'E',
                          'A',
                          'self weight (kg/m)',
                          'Prestress',
                          'Cable',
                          'Initial length',
                          'Member force',
                          'Tension',
                          'Final length']

            members_df = pd.DataFrame(memberTable, columns = memberlabels).set_index(['Member number'])
            self.members_dataFrame = members_df

            #Postprocessing
            members_df['E'] = round(members_df['E']/1e6,4)
            members_df['A'] = round(members_df['A']*1e6,4)
            members_df['Prestress'] = round(members_df['Prestress']/1000,4)
            members_df['Initial length'] = round(members_df['Initial length'],4)
            members_df['Cable'] = members_df['Cable'].replace([0,1],['False','True'])
            members_df['Member force'] = round(members_df['Member force']/1000,4)
            members_df['Tension'] = round(members_df['Tension']/1e6,4)
            members_df['Final length'] = round(members_df['Final length'],4)

            processed_memberlabels = ['Member number',
                          'node i',
                          'node j',
                          'E (N/mm2)',
                          'A (mm2)',
                          'self weight (kg/m)',
                          'Prestress (kN)',
                          'Cable',
                          'Initial length',
                          'Member force (kN)',
                          'Tension (N/mm2)',
                          'Final length']

            dict_labels = {}

            for i,lab in enumerate(memberlabels):
                dict_labels[lab] = processed_memberlabels[i]

            members_df = members_df.rename(dict_labels, axis='columns')
            self.processed_members_dataFrame = members_df

            with pd.ExcelWriter(file) as writer:
                self.log_data.to_excel(writer, sheet_name="Analysis")
                nodes_df.to_excel(writer, sheet_name="Nodes info")
                members_df.to_excel(writer, sheet_name="Members info")
                self.angles_variation.to_excel(writer, sheet_name='Angles variation')
                if len(self.Glazing_filters) > 0: self.glazing_variation.to_excel(writer, sheet_name='Glazing deviation')
                
    def set_Cable(self, **kwargs):
        """
        Method to define members of the structure as cable.

        Parameters
        ----------
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of member numbers to set as cable
                Memberfilter - Name of a existing member filter
                Memberlist - List of member numbers to define as cable
                Savefilter - Save the member list as a filter

        Returns
        -------
        int
            0 if nothing added.

        """
        
        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Memberfilter' in kwargs:
            indices = self.Member_filters[kwargs['Memberfilter']]
            
        elif 'Memberlist' in kwargs:
            indices = kwargs['Memberlist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:
            
            if not kwargs['Savefilter'] in self.Member_filters: self.Member_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')
        
        for n in indices:
            self.str_members[n].cable = True

    def set_Bar(self, **kwargs):
        """
        Method to define members of the structure as bar.

        Parameters
        ----------
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of member numbers to set as bar
                Memberfilter - Name of a existing member filter
                Memberlist - List of member numbers to define as bar
                Savefilter - Save the member list as a filter

        Returns
        -------
        int
            0 if nothing added.

        """
        
        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Memberfilter' in kwargs:
            indices = self.Member_filters[kwargs['Memberfilter']]
            
        elif 'Memberlist' in kwargs:
            indices = kwargs['Memberlist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:
            
            if not kwargs['Savefilter'] in self.Member_filters: self.Member_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')
        
        for n in indices:
            self.str_members[n].cable = False


    def set_A(self, Area, **kwargs):
        """
        Method to define cross section area of selected members.

        Parameters
        ----------
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of member numbers to set A
                Memberfilter - Name of a existing member filter
                Memberlist - List of member numbers to set A
                Savefilter - Save the member list as a filter

        Returns
        -------
        int
            0 if nothing added.

        """
        
        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Memberfilter' in kwargs:
            indices = self.Member_filters[kwargs['Memberfilter']]
            
        elif 'Memberlist' in kwargs:
            indices = kwargs['Memberlist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:
            
            if not kwargs['Savefilter'] in self.Member_filters: self.Member_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')
        
        for n in indices:
            self.str_members[n].A = Area
            self.str_members[n].initial_A = Area

    def set_E(self, E, **kwargs):
        """
        Method to define elasticity modulus of selected members.

        Parameters
        ----------
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of member numbers to set E
                Memberfilter - Name of a existing member filter
                Memberlist - List of member numbers to set E
                Savefilter - Save the member list as a filter

        Returns
        -------
        int
            0 if nothing added.

        """
        
        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Memberfilter' in kwargs:
            indices = self.Member_filters[kwargs['Memberfilter']]
            
        elif 'Memberlist' in kwargs:
            indices = kwargs['Memberlist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:
            
            if not kwargs['Savefilter'] in self.Member_filters: self.Member_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')
        
        for n in indices:
            self.str_members[n].E = E


    def set_P(self, P, **kwargs):
        """
        Method to define prestress of selected modules.

        Parameters
        ----------
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of member numbers to set P
                Memberfilter - Name of a existing member filter
                Memberlist - List of member numbers to set P
                Savefilter - Save the member list as a filter

        Returns
        -------
        int
            0 if nothing added.

        """
        
        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Memberfilter' in kwargs:
            indices = self.Member_filters[kwargs['Memberfilter']]
            
        elif 'Memberlist' in kwargs:
            indices = kwargs['Memberlist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:
            
            if not kwargs['Savefilter'] in self.Member_filters: self.Member_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')

        for n in indices:
            self.str_members[n].P = P
            

    def set_sw(self, sw, **kwargs):
        """
        Method to define self weight of selected modules.

        Parameters
        ----------
        **kwargs : TYPE
            Accepted:
                Csvfile - CSV containing list of member numbers to set sw
                Memberfilter - Name of a existing member filter
                Memberlist - List of member numbers to set sw
                Savefilter - Save the member list as a filter

        Returns
        -------
        int
            0 if nothing added.

        """
        
        if 'Csvfile' in kwargs:
            indices = genfromtxt(kwargs['Csvfile'], delimiter=',')
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
        
        elif 'Memberfilter' in kwargs:
            indices = self.Member_filters[kwargs['Memberfilter']]
            
        elif 'Memberlist' in kwargs:
            indices = kwargs['Memberlist']
            indices = np.int_(indices)
            if isinstance(indices, int): indices = [indices]
            
        else:
            print('No data introduced. Please specify a file, list or filter')
            return 0
            
        if 'Savefilter' in kwargs:
            
            if not kwargs['Savefilter'] in self.Member_filters: self.Member_filters[kwargs['Savefilter']] = indices
            else: print('Filter name already present in filters list. Not saved.')
        
        for n in indices:
            self.str_members[n].sw = sw



    def calculateTransMatrices(self,UG):

        """
        Internal method to be called by Calculate_Structure
        
        Calculate transformation matrices for each member based on current deformed shape for the structure.
        Current deformed shape is initial position plus cumulative displacements UG.
        
        """
        #Initialise container to hold transformation matrices
        TransformationMatrices = np.zeros((self.nMembers, 3, 6))

        for n in self.str_members:
            node_i = self.str_members[n].node_i
            node_j = self.str_members[n].node_j

            #Index of DoF for this member
            i_x = 3*node_i #index of x-DoF at node i
            i_y = 3*node_i+1 #index of y-DoF at node i
            i_z = 3*node_i+2 #index of z-DoF at node i
            j_x = 3*node_j #index of x-DoF at node j
            j_y = 3*node_j+1 #index of y-DoF at node j
            j_z = 3*node_j+2 #index of z-DoF at node j

            #New positions = initial pos + cum deflection
            ix = self.str_members[n].ix + UG[i_x,0]
            iy = self.str_members[n].iy + UG[i_y,0]
            iz = self.str_members[n].iz + UG[i_z,0]
            jx = self.str_members[n].jx + UG[j_x,0]
            jy = self.str_members[n].jy + UG[j_y,0]
            jz = self.str_members[n].jz + UG[j_z,0]

            #Local transformation matrix for each element, TM
            TM = np.zeros([3,6])

            dx = jx-ix #x-component of vector along member
            dy = jy-iy #y-component of vector along member
            dz = jz-iz #z-component of vector along member
            length = math.sqrt(dx**2 + dy**2 + dz**2) #Magnitude of vector (length of member)

            if(abs(dx)<0.001 and abs(dy)<0.001):
                #Member is vertical - offset in negative global x to define local x-y plane
                i_offset = np.array([ix-1, iy, iz]) #Offset node i by 1m in negative global x-direction
                j_offset = np.array([jx-1, jy, jz]) #Offset node j by 1m in negative global x-direction
            else:
                #Member is not vertical - offset in positive global z to define local x-y plane
                i_offset = np.array([ix, iy, iz+1]) #Offset node i by 1m in global z-direction
                j_offset = np.array([jx, jy, jz+1]) #Offset node j by 1m in global z-direction
            node_k = i_offset + 0.5*(j_offset-i_offset) #Point in the local x-y plane

            #Local x-vector in global RF running along the member
            local_x_vector = np.array([jx, jy, jz]) - np.array([ix, iy, iz])#Vector along local x-axis
            local_x_unit =  local_x_vector/length #Local unit vector defining local x-axis

            #Local y-vector in global RF using Gram-Schmidt process
            node_i = np.array([ix, iy, iz])
            vector_in_plane = node_k- node_i #Vector in the x-y plane
            local_y_vector = vector_in_plane - np.dot(vector_in_plane,local_x_unit)*local_x_unit #local y-vector in global RF (Gram-Schmidt)
            magY = math.sqrt(local_y_vector[0]**2 + local_y_vector[1]**2+local_y_vector[2]**2) #Length of local y-vector
            local_y_unit = local_y_vector/magY #Local unit vector defining local y-axis

            #Local z-vector in global RF using matrix cross product
            local_z_unit = np.cross(local_x_unit, local_y_unit) #Local unit vector defining local z-axis

            #Combine reference frame into standard rotation matrix for the element vector x,y,z => column 1,2,3
            rotationMatrix = np.array([local_x_unit, local_y_unit, local_z_unit]).T
            TM = np.concatenate((-rotationMatrix.T,rotationMatrix.T),axis=1)

            #Store transformation matrix for current member
            TransformationMatrices[n,:,:] = TM

        return TransformationMatrices

    def initPretension(self,TMs):
        """
        Internal method to be called by Calculate_Structure
        
        P = axial pre-tension specified for each bar
        Calculate the force vector [La] for each bar [La] = [T'][AA'][P]
        Combine into an overal vector representing the internal force system and return
        
        """
        F_pre = np.array([np.zeros(self.nDoF)]).T #Initialse internal force vector

        for n in self.str_members:

            TM = TMs[n,:,:]
            AAp = np.array([[1],[0],[0]])
            F_pre_global = np.matmul(TM.T,AAp)*self.str_members[n].P

            F_pre[3*self.str_members[n].node_i,0] = F_pre[3*self.str_members[n].node_i,0] + F_pre_global[0][0]
            F_pre[3*self.str_members[n].node_i+1,0] = F_pre[3*self.str_members[n].node_i+1,0] + F_pre_global[1][0]
            F_pre[3*self.str_members[n].node_i+2,0] = F_pre[3*self.str_members[n].node_i+2,0] + F_pre_global[2][0]
            F_pre[3*self.str_members[n].node_j,0] = F_pre[3*self.str_members[n].node_j,0] + F_pre_global[3][0]
            F_pre[3*self.str_members[n].node_j+1,0] = F_pre[3*self.str_members[n].node_j+1,0] + F_pre_global[4][0]
            F_pre[3*self.str_members[n].node_j+2,0] = F_pre[3*self.str_members[n].node_j+2,0] + F_pre_global[5][0]

        # print(F_pre)
        return F_pre


    def buildStructureStiffnessMatrix(self, UG, TMs):
        """
        Internal method to be called by Calculate_Structure
        
        Standard construction of Primary and Structure stiffness matrix
        Construction of non-linear element stiffness matrix handled in a child function.
        
        """
        Kp = np.zeros([self.nDoF,self.nDoF]) #Initialise the primary stiffness matrix

        for n in self.str_members:

            # Element stiffness matrix

            #Index of DoF for this member
            i_x = 3*self.str_members[n].node_i
            i_y = 3*self.str_members[n].node_i+1
            i_z = 3*self.str_members[n].node_i+2
            j_x = 3*self.str_members[n].node_j
            j_y = 3*self.str_members[n].node_j+1
            j_z = 3*self.str_members[n].node_j+2

            #Displacements
            d_ix = UG[i_x,0]
            d_iy = UG[i_y,0]
            d_iz = UG[i_z,0]
            d_jx = UG[j_x,0]
            d_jy = UG[j_y,0]
            d_jz = UG[j_z,0]

            #Extract current version of transformation matrix [T]
            TM = TMs[n,:,:]

            localDisp = np.matmul(TM,np.array([[d_ix, d_iy, d_iz, d_jx, d_jy, d_jz]]).T)
            u = localDisp[0].item()
            v = localDisp[1].item()
            w = localDisp[2].item()

            #Calculate extension, e
            Lo = self.str_members[n].length
            e = math.sqrt((Lo+u)**2 + v**2 + w**2)-Lo

            #Calculate matrix [AA]
            a1 = (Lo+u)/(Lo+e)
            a2 = v/(Lo+e)
            a3 = w/(Lo+e)
            AA = np.array([[a1,a2,a3]])

            #Calculate axial load, P
            P = self.str_members[n].P + (self.str_members[n].E*self.str_members[n].A/Lo)*e

            #Calculate matrix [d]
            d11 = P*(v**2 + w**2)
            d12 = -P*v*(Lo+u)
            d13 = -P*w*(Lo+u)
            d21 = -P*v*(Lo+u)
            d22 = P*((Lo+u)**2+w**2)
            d23 = -P*v*w
            d31 = -P*w*(Lo+u)
            d32 = -P*v*w
            d33 = P*((Lo+u)**2+v**2)
            denominator = (Lo+e)**3

            d = (1/denominator)*np.array([[d11,d12,d13],
                                          [d21,d22,d23],
                                          [d31,d32,d33]])

            #Calculate element stiffness matrix
            NL = np.matrix((AA.T*(self.str_members[n].E * self.str_members[n].A/Lo)*AA)+d)
            k = TM.T*NL*TM

            #Return element stiffness matrix in quadrants
            K11 = k[0:3,0:3]
            K12 = k[0:3,3:6]
            K21 = k[3:6,0:3]
            K22 = k[3:6,3:6]

            Kp[i_x:i_z+1,i_x:i_z+1] = Kp[i_x:i_z+1,i_x:i_z+1] + K11
            Kp[i_x:i_z+1,j_x:j_z+1] = Kp[i_x:i_z+1,j_x:j_z+1] + K12
            Kp[j_x:j_z+1,i_x:i_z+1] = Kp[j_x:j_z+1,i_x:i_z+1] + K21
            Kp[j_x:j_z+1,j_x:j_z+1] = Kp[j_x:j_z+1,j_x:j_z+1] + K22

        #Reduce to structure stiffness matrix by deleting rows and columns for restrained DoF
        Ks = np.delete(Kp,self.restrainedIndex,0) #Delete rows
        Ks = np.delete(Ks,self.restrainedIndex,1) #Delete columns
        Ks = np.matrix(Ks) # Convert Ks from numpy.ndarray to numpy.matrix to use build in inverter function

        # print('Ks')
        # print(Ks)

        return Ks

    def solveDisplacements(self, Ks, F_inequilibrium):
        """
        Internal method to be called by Calculate_Structure
        
        Standard solving for structural displacements.
        
        """

        forceVectorRed = copy.copy(F_inequilibrium)# Make a copy of forceVector so the copy can be edited, leaving the original unchanged
        forceVectorRed = np.delete(forceVectorRed,self.restrainedIndex,0) #Delete rows corresponding to restrained DoF
        # print(Ks.shape)
        # print('F_inequilibrium')
        # print(forceVectorRed.shape)
        # print('forcevectorRed')
        # print(forceVectorRed.shape)
        U = Ks.I*forceVectorRed
        # U = np.linalg.solve(Ks,forceVectorRed)
        # print(U)
        #Build the global displacement vector inclusing zeros as restrained degrees of freedom
        UG = np.zeros(self.nDoF) #Initialise an array to hold the global displacement vector
        c=0 #Initialise a counter to track how many restraints have been imposed
        for i in np.arange(self.nDoF):
            if i in self.restrainedIndex:
                #Impose zero displacement
                UG[i] = 0
            else:
                #Assign actual displacement
                UG[i] = U[c]
                c=c+1

        UG = np.array([UG]).T

        return UG

    def updateInternalForceSystem(self, UG, TMs):
        """
        Internal method to be called by Calculate_Structure
        
        Calculate the vector of internal forces associated with the incremental displacements UG
        [Ln] = [T'][AA'][P].
        
        """

        F_int = np.zeros([self.nDoF,1]) #Initialise a container for all internal forces

        #Cycle through each member and calculate nodal forces in global coordinates then add to overal internal force system vector
        for n in self.str_members:

            i_x = 3*self.str_members[n].node_i #index 0 (e.g. node 1)
            i_z = 3*self.str_members[n].node_i+2 #index 2 (e.g. node 1)
            j_x = 3*self.str_members[n].node_j #index 3 (e.g. node 2)
            j_z = 3*self.str_members[n].node_j+2 #index 5 (e.g. node 2)

            #Displacements
            d_ix = UG[i_x,0]
            d_iy = UG[i_x+1,0]
            d_iz = UG[i_z,0]
            d_jx = UG[j_x,0]
            d_jy = UG[j_x+1,0]
            d_jz = UG[j_z,0]

            #Extract current version of transformation matrix [T]
            TM = TMs[n,:,:]

            #Calculate local displacements [u, v] using global cumulative displacements UG
            localDisp = np.matmul(TM,np.array([[d_ix, d_iy, d_iz, d_jx, d_jy, d_jz]]).T)
            u = localDisp[0].item()
            v = localDisp[1].item()
            w = localDisp[2].item()

            #Calculate extension, e
            Lo = self.str_members[n].length
            e = math.sqrt((Lo+u)**2 + v**2 + w**2)-Lo

            #Calculate matrix [AA]
            a1 = (Lo+u)/(Lo+e)
            a2 = v/(Lo+e)
            a3 = w/(Lo+e)
            AA = np.array([[a1,a2,a3]])

            #Calculate axial load, P, due to incremental deflections
            P = (self.str_members[n].E*self.str_members[n].A/Lo)*e

            #Determine axial load in global coords
            F_global = np.matmul(TM.T,AA.T)*P

            #Add member pre-tension to overall record
            F_int[i_x,0] = F_int[i_x,0] + F_global[0][0]
            F_int[i_x+1,0] = F_int[i_x+1,0] + F_global[1][0]
            F_int[i_z,0] = F_int[i_z,0] + F_global[2][0]
            F_int[j_x,0] = F_int[j_x,0] + F_global[3][0]
            F_int[j_x+1,0] = F_int[j_x+1,0] + F_global[4][0]
            F_int[j_z,0] = F_int[j_z,0] + F_global[5][0]

        return F_int

    def testForConvergence(self, it, F_inequilibrium):
        """
        Internal method to be called by Calculate_Structure
        
        Test if structure has converged by comparing the maximum force in the equilibrium
        force vector against a threshold for the simulation.
        
        """
        notConverged = True #Initialise the convergence flag
        maxIneq = 0
        if it>0:
            maxIneq = np.max(abs(F_inequilibrium[self.freeDoF]))
            if(maxIneq<self.solver_convThreshold):
                notConverged = False

        return notConverged, maxIneq

    def calculateMbrForces(self, UG):
        """
        Internal method to be called by Calculate_Structure
        
        Calculates the member forces based on change in length of each member
        Takes in the cumulative global displacement vector as UG.
        
        """

        mbrForces = np.zeros(self.nMembers) #Initialise a container to hold axial forces

        for n in self.str_members:

            i_x = 3*self.str_members[n].node_i #index 0 (e.g. node 1)
            i_z = 3*self.str_members[n].node_i+2 #index 2 (e.g. node 1)
            j_x = 3*self.str_members[n].node_j #index 3 (e.g. node 2)
            j_z = 3*self.str_members[n].node_j+2 #index 5 (e.g. node 2)

            #New positions = initial pos + cum deflection
            ix = self.str_members[n].ix + UG[i_x,0]
            iy = self.str_members[n].iy + UG[i_x+1,0]
            iz = self.str_members[n].iz + UG[i_z,0]
            jx = self.str_members[n].jx + UG[j_x,0]
            jy = self.str_members[n].jy + UG[j_x+1,0]
            jz = self.str_members[n].jz + UG[j_z,0]

            dx = jx-ix #x-component of vector along member
            dy = jy-iy #y-component of vector along member
            dz = jz-iz #z-component of vector along member
            newLength = math.sqrt(dx**2 + dy**2 + dz**2) #Magnitude of vector (length of member)

            deltaL = newLength - self.str_members[n].length #Change in length
            force = self.str_members[n].P + deltaL*self.str_members[n].E*self.str_members[n].A/self.str_members[n].length #Axial force due to change in length and any pre-tension
            mbrForces[n] = force #Store member force

        return mbrForces
    
    def displacement_postprocess(self):
        """
        Internal method to be called by Calculate_Structure
        
        Performs postprocessing calculations of angle variation in bars concurring
        in each node whose initial angle is 180 degrees. In this case it is assumed
        that the bars belong to the same cable and the calculation of the angle after
        the deformation.
        
        """
        
        # NODES
        for i in self.str_nodes:
            membersinnode = []
            print('node', i)
            for j in self.str_members:
                if self.str_members[j].node_i == i or self.str_members[j].node_j == i: membersinnode.append(j)
            print('members concurring in this node ', membersinnode)
            
            if len(membersinnode) > 1:
                members_180degrees = np.empty([0,2])
                
                for i in membersinnode:
                    for j in membersinnode:
                        if i != j: members_180degrees = np.vstack([members_180degrees, [i,j]])
                
                members_180degrees = np.sort(members_180degrees, axis=1) #sort values in rows of array
                members_180degrees = np.unique(members_180degrees, axis=0) #delete duplicate pairs of values
                
                print('pairs of members', members_180degrees)
                
                for i in range(np.shape(members_180degrees)[0]):
                    member1i = self.str_members[members_180degrees[i,0]].node_i #member 1, node i
                    member1j = self.str_members[members_180degrees[i,0]].node_j #member 1, node j
                    
                    member2i = self.str_members[members_180degrees[i,1]].node_i #member 2, node i
                    member2j = self.str_members[members_180degrees[i,1]].node_j #member 2, node j
                    
                    print('MIEMBROS', member1i, member1j, member2i, member2j)
                    
                    #Select which node is the common node (C) and which ones are de different (D)
                    if member1i == member2i:
                        member1C = member1i
                        member2C = member2i
                        member1D = member1j
                        member2D = member2j
                        currentnode = member1i
                    elif member1i == member2j:
                        member1C = member1i
                        member2C = member2j
                        member1D = member1j
                        member2D = member2i
                        currentnode = member1i
                    elif member1j == member2j:
                        member1C = member1j
                        member2C = member2j
                        member1D = member1i
                        member2D = member2i
                        currentnode = member1j
                    elif member1j == member2i:
                        member1C = member1j
                        member2C = member2i
                        member1D = member1i
                        member2D = member2j
                        currentnode = member1j
                    
                    #Angle between two bars in undeformed state:
                    v1 = np.array([self.str_nodes[member1D].x,self.str_nodes[member1D].y,self.str_nodes[member1D].z]) - np.array([self.str_nodes[member1C].x,self.str_nodes[member1C].y,self.str_nodes[member1C].z])
                    v2 = np.array([self.str_nodes[member2D].x,self.str_nodes[member2D].y,self.str_nodes[member2D].z]) - np.array([self.str_nodes[member2C].x,self.str_nodes[member2C].y,self.str_nodes[member2C].z])
                    #norms of vectors
                    mag_v1 = np.linalg.norm(v1)
                    mag_v2 = np.linalg.norm(v2)
                    dot_product = np.dot(v1, v2)
                    initial_angle = np.arccos(dot_product / (mag_v1 * mag_v2))
                    
                    print('BARS ', members_180degrees[i,0],members_180degrees[i,1])
                    print('Angle between bars',np.degrees(initial_angle))
                    
                    if abs(180 - np.degrees(initial_angle)) < 0.01:
                        
                        v1d = np.array([self.str_nodes[member1D].displaced_x,self.str_nodes[member1D].displaced_y,self.str_nodes[member1D].displaced_z]) - np.array([self.str_nodes[member1C].displaced_x,self.str_nodes[member1C].displaced_y,self.str_nodes[member1C].displaced_z])
                        v2d = np.array([self.str_nodes[member2D].displaced_x,self.str_nodes[member2D].displaced_y,self.str_nodes[member2D].displaced_z]) - np.array([self.str_nodes[member2C].displaced_x,self.str_nodes[member2C].displaced_y,self.str_nodes[member2C].displaced_z])
                        
                        mag_v1d = np.linalg.norm(v1d)
                        mag_v2d = np.linalg.norm(v2d)
                        dot_productd = np.dot(v1d, v2d)
                        final_angle = np.arccos(dot_productd / (mag_v1d * mag_v2d))
                        print('BARS ', members_180degrees[i,0],members_180degrees[i,1])
                        print('Final angle',np.degrees(final_angle))
                        
                        initial_angle = np.degrees(initial_angle)
                        final_angle = np.degrees(final_angle)
                        
                        angle_variation = abs(initial_angle-final_angle)
                        
                        datarow = pd.DataFrame.from_dict({'Node':[int(currentnode)], 'Bar 1':[int(members_180degrees[i,0])], 'Bar 2':[int(members_180degrees[i,1])], 'Initial angle':[round(initial_angle,2)], 'Final angle':[round(final_angle,2)], 'Variation (ABS)':[round(angle_variation,2)]})
                        self.angles_variation = self.angles_variation.append(datarow, ignore_index=True)
                    
        # GLAZING UNITS
        if len(self.Glazing_filters) > 0:
            
            self.glazing_variation = pd.DataFrame(columns=['Glas unit', 'Node 1', 'Node 2', 'Node 3', 'Node 4', 'Deviation from plane', 'Deviation ratio'])
            
            for glass in self.Glazing_filters:                       
                nodes = self.Glazing_filters[glass]
                
                #Point coordinates for initial diagonal distance
                p1x = self.str_nodes[nodes[0]].x
                p1y = self.str_nodes[nodes[0]].y
                p1z = self.str_nodes[nodes[0]].z
                
                p2x = self.str_nodes[nodes[1]].x
                p2y = self.str_nodes[nodes[1]].y
                p2z = self.str_nodes[nodes[1]].z
                
                p3x = self.str_nodes[nodes[2]].x
                p3y = self.str_nodes[nodes[2]].y
                p3z = self.str_nodes[nodes[2]].z
                
                p4x = self.str_nodes[nodes[3]].x
                p4y = self.str_nodes[nodes[3]].y
                p4z = self.str_nodes[nodes[3]].z
                
                #Distancies between points to get the diagonal
                
                dist_p4_to_p1 = np.sqrt((p4x-p1x)**2 + (p4y-p1y)**2 + (p4z-p1z)**2)
                dist_p4_to_p2 = np.sqrt((p4x-p2x)**2 + (p4y-p2y)**2 + (p4z-p2z)**2)
                dist_p4_to_p3 = np.sqrt((p4x-p3x)**2 + (p4y-p3y)**2 + (p4z-p3z)**2)
                
                diagonal = np.array([dist_p4_to_p1,dist_p4_to_p2,dist_p4_to_p3]).max()
                
                print('GLASS', glass)
                print('point 1 ',p1x,p1y,p1z)
                print('point 2 ',p2x,p2y,p2z)
                print('point 3 ',p3x,p3y,p3z)
                print('point 4 ',p4x,p4y,p4z)
                
                print('Distance point 4 to 1 ', dist_p4_to_p1)
                print('Distance point 4 to 2 ', dist_p4_to_p2)
                print('Distance point 4 to 3 ', dist_p4_to_p3)
                
                print('diagonal ',diagonal)
                
                #Displaced coordinates of three points to form the hypotetical plane of glass
                p1x_d = self.str_nodes[nodes[0]].displaced_x
                p1y_d = self.str_nodes[nodes[0]].displaced_y
                p1z_d = self.str_nodes[nodes[0]].displaced_z
                point1 = np.array([p1x_d, p1y_d, p1z_d])
                
                p2x_d = self.str_nodes[nodes[1]].displaced_x
                p2y_d = self.str_nodes[nodes[1]].displaced_y
                p2z_d = self.str_nodes[nodes[1]].displaced_z
                point2 = np.array([p2x_d, p2y_d, p2z_d])
                
                p3x_d = self.str_nodes[nodes[2]].displaced_x
                p3y_d = self.str_nodes[nodes[2]].displaced_y
                p3z_d = self.str_nodes[nodes[2]].displaced_z
                point3 = np.array([p3x_d, p3y_d, p3z_d])
                
                #Fourth point to calculate the distance to the plane generated by the other three points
                p4x_d = self.str_nodes[nodes[3]].displaced_x
                p4y_d = self.str_nodes[nodes[3]].displaced_y
                p4z_d = self.str_nodes[nodes[3]].displaced_z
                point4 = np.array([p4x_d, p4y_d, p4z_d])
                
                #Calculate normal vector to plane
                v1 = point2 - point1
                v2 = point3 - point1
                normal = np.cross(v1, v2)
                
                #Calculate distance as dot product between the normal and any vector from the plane to the point4
                distance = abs(np.dot(normal, point4 - point1)) / np.linalg.norm(normal)
                
                print('Deviation of node 4 to the plane of glass ', distance)
                
                #Ratio
                metric = distance/diagonal
                print('metric ', metric)
                metric_angle = np.degrees(np.arctan(metric))
                print('metric_angle ',metric_angle)
                
                datarow = pd.DataFrame.from_dict({'Glas unit':[glass], 'Node 1':[nodes[0]], 'Node 2':[nodes[1]], 'Node 3':[nodes[2]], 'Node 4':[nodes[3]], 'Deviation from plane':[round(distance,4)], 'Deviation ratio':[round(metric_angle,4)]})
                self.glazing_variation = self.glazing_variation.append(datarow, ignore_index=True)
                
    def calculate_structure(self):
        """
        Main method of FTES_Structure. Performs structural analysis
        
        This method launches the analysis of the structure executing the iterative
        method of solving the equation f = Kd increasing the load in steps according
        to the configured parameters and iterating to increase the rigidity of the structure
        with the addition of loads to achieve equilibrium.
        
        Results are stored in the data structures of the class, ready to be exported
        with the provided functions.

        Returns
        -------
        None.

        """

        self.log_record("STRUCTURE ANALYSIS STARTED")
        #Add self weight to nodes
        if(self.solver_selfweight):
            SW_at_supports = np.empty((0,2))
            for n in self.str_members:

                F_node = self.str_members[n].length * self.str_members[n].sw * 9.81 / 2 #(N) Self-weight distributed into each node

                self.forceVector[self.str_members[n].node_i*3 + 2] = self.forceVector[self.str_members[n].node_i*3 + 2] - F_node #
                self.forceVector[self.str_members[n].node_j*3 + 2] = self.forceVector[self.str_members[n].node_j*3 + 2]-F_node #

               #Check if SW needs to be directly added to supports (if bars connect to supports)
                if(self.str_members[n].node_i * 3 + 2 in self.restrainedIndex):
                    supportSW = np.array([self.str_members[n].node_i * 3 + 2, F_node])
                    SW_at_supports = np.append(SW_at_supports, [supportSW], axis=0) #Store z-DoF at support and force to be added
                if(self.str_members[n].node_j * 3 + 2 in self.restrainedIndex):
                    supportSW = np.array([self.str_members[n].node_j * 3 + 2, F_node])
                    SW_at_supports = np.append(SW_at_supports, [supportSW], axis=0) #Store z-DoF at support and force to be added

        self.freeDoF = np.delete( np.arange(0,self.nDoF),self.restrainedIndex)
        self.forceIncrement = self.forceVector/self.solver_forceIncrements
        self.maxForce= self.forceVector
        self.forceVector = self.forceIncrement

        #Initialise a container to hold the set of global displacements for each external load increment
        UG_FINAL = np.empty([self.nDoF,0])

        #Initialise a container to hold the set of internal forces for each external load increment
        FI_FINAL = np.empty([self.nDoF,0])

        #Initialise a container to hold the set of axial forces for each external load increment
        EXTFORCES = np.empty([self.nDoF,0])

        #Initialise a container to hold the set of axial forces for each external load increment
        MBRFORCES = np.empty([self.nMembers,0])

        #Initialise global disp vector
        UG = np.zeros([self.nDoF,1]) #Initialise global displacement vector to zero (undeformed state)

        #Calculate initial transformation matrices for all members based on undeformed position
        TMs = self.calculateTransMatrices(UG)

        #Calculate internal force system based on any pre-tension in members [Lo]
        F_pre = self.initPretension(TMs)

        #Initialise a container to store incremental displacements calculated for each iteration [Xa], [Xb] etc.
        UG_inc = np.empty([self.nDoF,0])
        UG_inc = np.append(UG_inc, UG, axis=1) #Add the initial (zero) displacement record

        #Initialise a container to store incremental internal forces calculated for each iteration [La], [Lb] etc.
        F_inc = np.empty([self.nDoF,0])
         
        F_inc = np.append(F_inc, F_pre, axis=1) #Add the initial pre-tension force record (external forces less any pre-tension)

        i=0 #Initialise an iteration counter (zeros out for each load increment)
        inc=0 #Initialise load increment counter
        notConverged = True #Initialise convergence flag

        while notConverged and i<self.solver_iterationLimit:

            #Calculate the cumulative internal forces L_total = La + Lb + Lc + ...
            Fi_total = np.matrix(np.sum(F_inc,axis=1)).T #Sum across columns of F_inc

            #Calculate the cumulative incremental displacements X_total = Xa + Xb + Xc + ...
            UG_total = np.matrix(np.sum(UG_inc,axis=1)).T #Sum across columns of UG_inc

            #Inequilibrium force vector used in this iteration [L]-[Ln] or externalForces - (cumulative) InternalForceSystem
            F_inequilibrium = self.forceVector - Fi_total

            #Build the structure stiffness matrix based on current position (using cumulative displacements)
            Ks = self.buildStructureStiffnessMatrix(UG_total, TMs)

            #Solve for global (incremental) displacement vector [Xn] for this iteration
            UG = self.solveDisplacements(Ks, F_inequilibrium)

            #Calculate a new transformation matrix for each member based on cum disp up to previous iteration
            TMs = self.calculateTransMatrices(UG_total)

            #Calculate the internal force system based on new incremental displacements, [Ln]
            F_int = self.updateInternalForceSystem(UG, TMs)

            #Save incremental displacements and internal forces for this iteration
            UG_inc = np.append(UG_inc, UG, axis=1)
            F_inc = np.append(F_inc, F_int, axis=1)

            #Test for convergence
            notConverged, maxIneq = self.testForConvergence(i, F_inequilibrium)

            i+=1

            #If system has converged, test for slack elements, save converged displacements, forces and increment external loading
            if not notConverged:

                hasSlackElements = False #Initialise a flag to indicate presence of new slack elements
                mbrForces = self.calculateMbrForces(UG_total) #Calculate member forces based on current set of displacements

                #Test for compression in cable elements if designated number of converged increments reached
                if inc > self.solver_checkSlack:
                    for n in self.str_members:
                        if self.str_members[n].cable == True and mbrForces[n]<0:
                            self.log_record(f'Compression in cable element from from nodes {self.str_members[n].node_i} to {self.str_members[n].node_j}')
                            hasSlackElements = True #Switch slack elements flag
                            self.str_members[n].A = 0 #Eliminate member stiffness by seting its cross-sectional area to zero

                if not hasSlackElements:
                    self.log_record(f'System has converged for load increment {inc} after {i-1} iterations')
                    UG_FINAL = np.append(UG_FINAL, UG_total, axis=1) #Add the converged displacement record
                    UG_inc = np.empty([self.nDoF,0]) #Zero out the record of incremental displacements for the next load increment
                    UG_inc = np.array(np.append(UG_inc, UG_total, axis=1)) #Add the initial displacement record for next load increment (manually cast as ndarray instead of matrix)

                    FI_FINAL = np.append(FI_FINAL, Fi_total, axis=1) #Add the converged force record
                    F_inc = np.empty([self.nDoF,0]) #Zero out the record of incremental forces for the next load increment
                    F_inc = np.array(np.append(F_inc, Fi_total, axis=1)) #Add the initial force record for next load increment (manually cast as ndarray instead of matrix)

                    mbrForces = self.calculateMbrForces(UG_FINAL[:,-1]) #Calculate the member forces based on change in mbr length
                    MBRFORCES = np.append(MBRFORCES, np.matrix(mbrForces).T, axis=1) #Add the converged axial forces record

                    EXTFORCES = np.append(EXTFORCES, self.forceVector, axis=1) #Add the external force vector for this load increment

                    #Test if all external loading has been applied
                    if abs(sum(self.forceVector).item()) < abs(sum(self.maxForce).item()):
                        i=0 #Reset counter for next load increment
                        inc +=1
                        self.forceVector = self.forceVector + self.forceIncrement #Increment the applied load
                        notConverged = True #Reset notConverged flag for next load increment
                else:
                    self.log_record('Rerun load increment with slack element stiffness set to zero')
                    F_inc = np.empty([self.nDoF,0]) #Zero out the record of incremental forces
                    F_inc = np.array(np.append(F_inc, FI_FINAL[:,-1], axis=1)) #Reset F_inc to values from previous converged load increment

                    UG_inc = np.empty([self.nDoF,0]) #Zero out the record of incremental displacements for the next load increment
                    UG_inc = np.array(np.append(UG_inc, UG_FINAL[:,-1], axis=1)) #Reset UG_inc to values from previous converged load increment

                    #Reset the external forceVector to last converged value
                    self.forceVector = np.array([EXTFORCES[:,-1]]).T

                    i=0 #Reset iteration counter
                    notConverged = True #Reset notConverged flag for next load increment

        # POST-PROCESSING, PLOTTING AND DATA EXPORT
        self.solved_structure = True
        self.log_record("STRUCTURE ANALYSIS FINISHED")
        
        reactionsFlag = False #Initialise a flag so we can plot a message re. reactions when necessary
        if(self.solver_selfweight):
            if SW_at_supports.size>0:
                reactionsFlag = True
                for SW in SW_at_supports:
                    index = int(SW[0]) #Index of the global force vector 'FG' to update
                    FI_FINAL[index,:] = FI_FINAL[index,:] + SW[1] #Add nodal SW force directly to FG

        #Generate output statements
        self.log_record("OUTSTANDING FORCE IMBALANCE")
        maxInequality = round(max(abs(F_inequilibrium[self.freeDoF,0])).item()/1000,3)
        self.log_record(f"(max = {maxInequality} kN)")

        for i in np.arange(0,self.nDoF):
            if i not in self.restrainedIndex:
                self.log_record(f'Remaining force imbalance at DoF {i} is {round(F_inequilibrium[i,0]/1000,3)} kN')

        print("")
        self.log_record("REACTIONS")
        f_int = FI_FINAL[:,-1]
        for i in np.arange(0,len(self.restrainedIndex)):
            index = self.restrainedIndex[i]
            self.log_record(f"Reaction at DoF {index}: {round(f_int[index].item()/1000,2)} kN")

        for i in range(self.nNodes):

            dof_x = 3*i
            if dof_x in self.restrainedIndex:
                self.str_nodes[i].reaction_x = f_int[dof_x].item()
            dof_y = 3*i+1
            if dof_y in self.restrainedIndex:
                self.str_nodes[i].reaction_y = f_int[dof_y].item()
            dof_z = 3*i+2
            if dof_z in self.restrainedIndex:
                self.str_nodes[i].reaction_z = f_int[dof_z].item()

        print("")
        self.log_record("NODAL DISPLACEMENTS")
        ug = UG_FINAL[:,-1]
        for n in self.str_nodes:
            i_x = 3*n #x-DoF index for this node
            i_y = 3*n+1 #y-DoF index for this node
            i_z = 3*n+2 #z-DoF index for this node
            self.str_nodes[n].set_U(ug[i_x,0],ug[i_y,0],ug[i_z,0])

            ux = round(ug[i_x,0],5) #x-direction nodal displacement
            uy = round(ug[i_y,0],5) #y-direction nodal displacement
            uz = round(ug[i_z,0],5) #z-direction nodal displacement
            self.log_record(f"Node {n+1}: Ux = {ux} m, Uy = {uy} m, Uz = {uz}")

        print("")
        self.log_record("MEMBER FORCES (incl. any pre-tension)")
        for n in self.str_members:
            self.log_record(f"Force in member {n} (nodes {self.str_members[n].node_i} to {self.str_members[n].node_j}) is {round(mbrForces[n]/1000,2)} kN")
            self.str_members[n].set_results(self.str_nodes,mbrForces[n])
        
        self.displacement_postprocess() #postprocess angle variation and deformation of glazing panels